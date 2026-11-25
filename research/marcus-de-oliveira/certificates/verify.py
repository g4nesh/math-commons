#!/usr/bin/env python3
"""Exact finite-instance MDO verifier over Q(i). No floating arithmetic is used."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path
import re
import sys

MAX_DIMENSION = 6
MAX_INPUT_BYTES = 1_048_576
RATIONAL = re.compile(r"^[+-]?[0-9]{1,64}(?:/[1-9][0-9]{0,63})?$")


class InvalidWitness(ValueError):
    """The input does not define a supported exact normal-pair witness."""


@dataclass(frozen=True)
class QComplex:
    real: Fraction = Fraction(0)
    imag: Fraction = Fraction(0)

    def __post_init__(self):
        for name in ("real", "imag"):
            value = getattr(self, name)
            if not isinstance(value, (int, Fraction)) or isinstance(value, bool):
                raise TypeError("Gaussian-rational components must be exact integers or Fractions")
            object.__setattr__(self, name, Fraction(value))

    def __add__(self, other):
        other = coerce(other)
        return QComplex(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return QComplex(-self.real, -self.imag)

    def __sub__(self, other):
        return self + -coerce(other)

    def __rsub__(self, other):
        return coerce(other) - self

    def __mul__(self, other):
        other = coerce(other)
        return QComplex(self.real * other.real - self.imag * other.imag,
                        self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = coerce(other)
        norm = other.real * other.real + other.imag * other.imag
        if norm == 0:
            raise ZeroDivisionError("Division by zero Gaussian rational")
        numerator = self * other.conjugate()
        return QComplex(numerator.real / norm, numerator.imag / norm)

    def conjugate(self):
        return QComplex(self.real, -self.imag)

    def pair(self):
        return (self.real, self.imag)

    def encode(self):
        return [str(self.real), str(self.imag)]


def coerce(value):
    if isinstance(value, QComplex):
        return value
    if isinstance(value, (int, Fraction)) and not isinstance(value, bool):
        return QComplex(Fraction(value))
    raise TypeError("Arithmetic requires exact Gaussian rationals")


ZERO, ONE = QComplex(), QComplex(Fraction(1))


def parse_complex(value, label="complex value"):
    if not isinstance(value, list) or len(value) != 2:
        raise InvalidWitness(f"{label} must be [real-rational-string, imag-rational-string]")
    for item in value:
        if not isinstance(item, str) or not RATIONAL.fullmatch(item):
            raise InvalidWitness(f"{label} requires integer or p/q strings with positive denominators and at most 64 digits per integer")
    return QComplex(Fraction(value[0]), Fraction(value[1]))


def reject_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InvalidWitness(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load_witness(path):
    raw = Path(path).read_bytes()
    if len(raw) > MAX_INPUT_BYTES:
        raise InvalidWitness("Witness exceeds the 1 MiB input limit")
    try:
        return json.loads(raw, object_pairs_hook=reject_duplicates,
                          parse_constant=lambda value: (_ for _ in ()).throw(InvalidWitness(f"Non-JSON constant: {value}")))
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise InvalidWitness(f"Invalid UTF-8 JSON: {error}") from error


def parse_witness(witness):
    allowed = {"schemaVersion", "name", "source", "notes", "a", "b", "U", "expectedDeterminant"}
    if not isinstance(witness, dict) or set(witness) - allowed:
        raise InvalidWitness("Witness must be an object containing only the documented fields")
    if type(witness.get("schemaVersion")) is not int or witness["schemaVersion"] != 1:
        raise InvalidWitness("schemaVersion must be integer 1")
    for key in ("name", "source", "notes"):
        if key in witness and not isinstance(witness[key], str):
            raise InvalidWitness(f"{key} must be text")
    if not isinstance(witness.get("a"), list) or not 1 <= len(witness["a"]) <= MAX_DIMENSION:
        raise InvalidWitness(f"Dimension must be between 1 and {MAX_DIMENSION}")
    n = len(witness["a"])
    if not isinstance(witness.get("b"), list) or len(witness["b"]) != n:
        raise InvalidWitness("a and b must have the same length")
    if not isinstance(witness.get("U"), list) or len(witness["U"]) != n or any(not isinstance(row, list) or len(row) != n for row in witness["U"]):
        raise InvalidWitness("U must be an n-by-n matrix matching the spectra")
    a = [parse_complex(value, f"a[{index}]") for index, value in enumerate(witness["a"])]
    b = [parse_complex(value, f"b[{index}]") for index, value in enumerate(witness["b"])]
    unitary = [[parse_complex(value, f"U[{row}][{col}]") for col, value in enumerate(values)] for row, values in enumerate(witness["U"])]
    expected = parse_complex(witness["expectedDeterminant"], "expectedDeterminant") if "expectedDeterminant" in witness else None
    return a, b, unitary, expected


def adjoint(matrix):
    return [[matrix[j][i].conjugate() for j in range(len(matrix))] for i in range(len(matrix))]


def multiply(left, right):
    n = len(left)
    return [[sum((left[i][k] * right[k][j] for k in range(n)), ZERO) for j in range(n)] for i in range(n)]


def diagonal(values):
    return [[values[i] if i == j else ZERO for j in range(len(values))] for i in range(len(values))]


def identity(n):
    return diagonal([ONE] * n)


def determinant(matrix):
    """Gaussian elimination over Q(i), with exact pivoting and singular handling."""
    values = [[coerce(value) for value in row] for row in matrix]
    if any(len(row) != len(values) for row in values):
        raise ValueError("Determinant requires a square matrix")
    n, result = len(values), ONE
    for col in range(n):
        pivot = next((row for row in range(col, n) if values[row][col] != ZERO), None)
        if pivot is None:
            return ZERO
        if pivot != col:
            values[col], values[pivot] = values[pivot], values[col]
            result = -result
        divisor = values[col][col]
        result = result * divisor
        for row in range(col + 1, n):
            factor = values[row][col] / divisor
            values[row][col] = ZERO
            for index in range(col + 1, n):
                values[row][index] = values[row][index] - factor * values[col][index]
    return result


def permutation_points(a, b):
    result = []
    for permutation in permutations(range(len(a))):
        value = ONE
        for index, partner in enumerate(permutation):
            value *= a[index] + b[partner]
        result.append((tuple(index + 1 for index in permutation), value))
    return result


def cross(origin, first, second):
    return ((first[0] - origin[0]) * (second[1] - origin[1])
            - (first[1] - origin[1]) * (second[0] - origin[0]))


def dot(normal, point):
    return normal[0] * point[0] + normal[1] * point[1]


def convex_hull(points):
    """Counterclockwise extreme vertices, with point and segment degeneracies."""
    points = sorted(set(points))
    if len(points) <= 1:
        return points
    lower, upper = [], []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def exact_scalar(value, label="coordinate"):
    if not isinstance(value, (int, Fraction)) or isinstance(value, bool):
        raise TypeError(f"{label} must be an exact integer or Fraction")
    return Fraction(value)


def exact_point(point, label="point"):
    if not isinstance(point, (tuple, list)) or len(point) != 2:
        raise TypeError(f"{label} must have two exact coordinates")
    return tuple(exact_scalar(value, label) for value in point)


def certify_hull(points, target):
    """Return exact sparse weights or a strict separating linear functional.

    Result weights refer to unique coordinate pairs, not permutation indices.
    A separate direct check below verifies the certificate against ALL points.
    """
    unique = sorted({exact_point(point) for point in points})
    target = exact_point(target)
    if not unique:
        raise ValueError("Cannot certify an empty hull")
    hull = convex_hull(unique)
    if len(hull) == 1:
        if target == hull[0]:
            certificate = {"membership": "inside", "weights": [(hull[0], Fraction(1))]}
        else:
            normal = (target[0] - hull[0][0], target[1] - hull[0][1])
            certificate = {"membership": "outside", "normal": normal}
    elif len(hull) == 2:
        first, second = hull
        displacement = (second[0] - first[0], second[1] - first[1])
        orientation = cross(first, second, target)
        if orientation != 0:
            normal = (-displacement[1], displacement[0])
            if orientation < 0:
                normal = (-normal[0], -normal[1])
            certificate = {"membership": "outside", "normal": normal}
        else:
            coordinate = 0 if displacement[0] else 1
            weight = (target[coordinate] - first[coordinate]) / displacement[coordinate]
            if 0 <= weight <= 1:
                certificate = {"membership": "inside", "weights": [(first, 1 - weight), (second, weight)]}
            else:
                normal = displacement if weight > 1 else (-displacement[0], -displacement[1])
                certificate = {"membership": "outside", "normal": normal}
    else:
        certificate = None
        for index, first in enumerate(hull):
            second = hull[(index + 1) % len(hull)]
            if cross(first, second, target) < 0:
                certificate = {"membership": "outside", "normal": (second[1] - first[1], first[0] - second[0])}
                break
        if certificate is None:
            for index in range(1, len(hull) - 1):
                first, second, third = hull[0], hull[index], hull[index + 1]
                area = cross(first, second, third)
                weights = (cross(target, second, third) / area,
                           cross(first, target, third) / area,
                           cross(first, second, target) / area)
                if all(value >= 0 for value in weights):
                    certificate = {"membership": "inside", "weights": list(zip((first, second, third), weights))}
                    break
            if certificate is None:
                raise ArithmeticError("Hull triangulation failed; no result certified")
    certificate["hull"] = hull
    if certificate["membership"] == "inside":
        certificate["weights"] = [(point, weight) for point, weight in certificate["weights"] if weight != 0]
    else:
        normal = certificate["normal"]
        certificate["bound"] = max(dot(normal, point) for point in unique)
        certificate["targetValue"] = dot(normal, target)
        certificate["gap"] = certificate["targetValue"] - certificate["bound"]
    check_hull_certificate(unique, target, certificate)
    return certificate


def check_hull_certificate(points, target, certificate):
    """Verify directly over exact scalar inputs; never accept float/bool values.

    This helper is callable separately from the generator. Every scalar used in
    its arithmetic is validated here, so approximate inputs cannot smuggle in a
    floating-point equality through a supplied certificate.
    """
    points = [exact_point(point, "hull point") for point in points]
    target = exact_point(target, "target")
    if not points or not isinstance(certificate, dict):
        raise ArithmeticError("A certificate needs a nonempty point set")
    if certificate.get("membership") == "inside":
        supplied_weights = certificate.get("weights")
        if not isinstance(supplied_weights, (list, tuple)) or not 1 <= len(supplied_weights) <= 3:
            raise ArithmeticError("Invalid convex-combination support")
        weights = []
        for entry in supplied_weights:
            if not isinstance(entry, (tuple, list)) or len(entry) != 2:
                raise ArithmeticError("Invalid convex-combination entry")
            weights.append((exact_point(entry[0], "support point"), exact_scalar(entry[1], "weight")))
        if any(point not in points or weight < 0 for point, weight in weights):
            raise ArithmeticError("Invalid convex-combination support")
        if sum(weight for _, weight in weights) != 1:
            raise ArithmeticError("Convex weights do not sum to one")
        reconstructed = tuple(sum(point[coordinate] * weight for point, weight in weights) for coordinate in (0, 1))
        if reconstructed != target:
            raise ArithmeticError("Convex combination does not equal the determinant")
    elif certificate.get("membership") == "outside":
        normal = exact_point(certificate.get("normal"), "separating normal")
        bound = exact_scalar(certificate.get("bound"), "separating bound")
        gap = exact_scalar(certificate.get("gap"), "separation gap")
        target_value = dot(normal, target)
        if "targetValue" in certificate and exact_scalar(certificate["targetValue"], "target value") != target_value:
            raise ArithmeticError("Claimed target functional value is incorrect")
        if normal == (0, 0) or any(dot(normal, point) > bound for point in points):
            raise ArithmeticError("Invalid separating bound")
        if target_value <= bound or gap != target_value - bound:
            raise ArithmeticError("Separation gap is not strictly positive and exact")
    else:
        raise ArithmeticError("Unknown membership certificate")
    return True


def verify_witness(witness):
    a, b, unitary, expected = parse_witness(witness)
    n = len(a)
    star = adjoint(unitary)
    if multiply(unitary, star) != identity(n) or multiply(star, unitary) != identity(n):
        raise InvalidWitness("U is not exactly unitary; approximate inputs are not rounded or repaired")
    matrix_a = diagonal(a)
    matrix_b = multiply(multiply(unitary, diagonal(b)), star)
    if multiply(matrix_a, adjoint(matrix_a)) != multiply(adjoint(matrix_a), matrix_a):
        raise ArithmeticError("Constructed A failed exact normality")
    if multiply(matrix_b, adjoint(matrix_b)) != multiply(adjoint(matrix_b), matrix_b):
        raise ArithmeticError("Constructed B failed exact normality")
    matrix_sum = [[matrix_a[row][col] + matrix_b[row][col] for col in range(n)] for row in range(n)]
    value = determinant(matrix_sum)
    if expected is not None and value != expected:
        raise InvalidWitness(f"Expected determinant mismatch: recomputed {value.encode()}, supplied {expected.encode()}")
    sigma_points = permutation_points(a, b)
    point_to_permutation = {}
    for permutation, point in sigma_points:
        point_to_permutation.setdefault(point.pair(), permutation)
    certificate = certify_hull(list(point_to_permutation), value.pair())
    result = {
        "schemaVersion": 1,
        "name": witness.get("name", "unnamed finite witness"),
        "dimension": n,
        "arithmetic": "exact Gaussian rationals (fractions.Fraction); no floating-point tolerance",
        "unitary": True,
        "normalA": True,
        "normalB": True,
        "determinant": value.encode(),
        "singularSum": value == ZERO,
        "permutationCount": len(sigma_points),
        "uniquePointCount": len(point_to_permutation),
        "permutationPoints": [{"permutation": list(permutation), "point": point.encode()} for permutation, point in sigma_points],
        "hull": [[str(coordinate) for coordinate in point] for point in certificate["hull"]],
        "hullDimension": min(2, len(certificate["hull"]) - 1),
        "membership": certificate["membership"],
        "AplusB": [[value.encode() for value in row] for row in matrix_sum],
        "scope": "One finite exact input over Q(i), dimension at most 6. This is not a proof of the universal conjecture or independent peer review.",
    }
    if certificate["membership"] == "inside":
        result["convexCombination"] = [{"permutation": list(point_to_permutation[point]), "point": [str(coordinate) for coordinate in point], "weight": str(weight)} for point, weight in certificate["weights"]]
    else:
        result["separator"] = {"normal": [str(coordinate) for coordinate in certificate["normal"]], "bound": str(certificate["bound"]), "targetValue": str(certificate["targetValue"]), "gap": str(certificate["gap"]), "interpretation": "normal[0]*Re(z) + normal[1]*Im(z) <= bound for every permutation point; determinant exceeds bound"}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("witness", help="JSON witness file using exact rational strings")
    parser.add_argument("--output", help="Write the result JSON to this file instead of stdout")
    args = parser.parse_args(argv)
    try:
        result = verify_witness(load_witness(args.witness))
    except (InvalidWitness, OSError, ArithmeticError) as error:
        print(json.dumps({"error": str(error), "certified": False}), file=sys.stderr)
        return 2
    text = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf8")
    else:
        print(text, end="")
    # Outside is a successful computation requiring review, not an execution failure.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
