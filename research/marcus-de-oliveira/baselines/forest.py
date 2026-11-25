#!/usr/bin/env python3
"""Fixed exact Q(i) checks for the forest-gate theorem; no random search.

Reuses the adjacent certificate verifier's exact scalar and matrix arithmetic.
This is an author-side regression, not an independent or formal proof checker.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations, permutations
import json
from pathlib import Path
import platform
import resource
import sys
import time

VERIFY = Path(__file__).resolve().parents[1] / "certificates" / "verify.py"
spec = importlib.util.spec_from_file_location("baseline_exact_arithmetic", VERIFY)
exact = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = exact
spec.loader.exec_module(exact)
Q, Z, O = exact.QComplex, exact.ZERO, exact.ONE
mul, adj, eye, diag, det = (exact.multiply, exact.adjoint, exact.identity,
                          exact.diagonal, exact.determinant)


def norm2(z):
    return z.real * z.real + z.imag * z.imag


def product(values):
    out = O
    for value in values:
        out *= value
    return out


def determinant_leibniz(matrix):
    """Second determinant algorithm, sharing only exact scalar arithmetic."""
    result, n = Z, len(matrix)
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        result += sign * product(matrix[i][p[i]] for i in range(n))
    return result


def validate_gate(n, gate):
    i, j, block = gate
    if type(n) is not int or not 1 <= n <= 6:
        raise ValueError("Fixture dimension must be an integer in 1..6")
    if type(i) is not int or type(j) is not int or not 0 <= i < j < n:
        raise ValueError("Gate endpoints must be distinct sorted indices")
    if len(block) != 2 or any(len(row) != 2 for row in block):
        raise ValueError("Gate must contain a 2-by-2 block")
    if any(not isinstance(z, Q) for row in block for z in row):
        raise TypeError("Gate entries must be exact Gaussian rationals")
    if mul(block, adj(block)) != eye(2) or mul(adj(block), block) != eye(2):
        raise ValueError("Gate block is not exactly unitary")


def validate_forest(n, gates):
    if type(n) is not int or not 1 <= n <= 6:
        raise ValueError("Fixture dimension must be an integer in 1..6")
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for gate in gates:
        validate_gate(n, gate)
        i, j, _ = gate
        a, b = root(i), root(j)
        if a == b:
            raise ValueError("Repeated edge or cycle is outside the forest theorem")
        parent[a] = b


def build_unitary(n, gates, require_forest=True):
    if require_forest:
        validate_forest(n, gates)
    result = eye(n)
    for gate in gates:
        validate_gate(n, gate)
        i, j, block = gate
        embedded = eye(n)
        for row, actual_row in enumerate((i, j)):
            for col, actual_col in enumerate((i, j)):
                embedded[actual_row][actual_col] = block[row][col]
        result = mul(embedded, result)
    return result


def permutation_weights(n, gates, require_forest=True):
    """Weights on maps p: input column -> output row, G_m ... G_1 order."""
    if require_forest:
        validate_forest(n, gates)
    weights = {tuple(range(n)): F(1)}
    for gate in gates:
        validate_gate(n, gate)
        i, j, block = gate
        q = norm2(block[0][1])
        following = {}
        for p, weight in weights.items():
            swapped = tuple(j if v == i else i if v == j else v for v in p)
            for candidate, factor in ((p, 1 - q), (swapped, q)):
                following[candidate] = following.get(candidate, F(0)) + weight * factor
        weights = {p: w for p, w in following.items() if w}
    return weights


def vertex(a, b, p):
    # p(j) is the row receiving b_j, so sigma=p^{-1}.
    return product(a[p[j]] + b[j] for j in range(len(a)))


def phased_block(c=F(3, 5), s=F(4, 5), variant=0):
    """L [[c,s],[-s,c]] R with exact independent unit phases."""
    phases = [O, Q(0, 1), Q(F(3, 5), F(4, 5)), Q(F(5, 13), F(-12, 13))]
    left = diag([phases[variant % 4], phases[(variant + 1) % 4]])
    right = diag([phases[(variant + 2) % 4], phases[(variant + 3) % 4]])
    return mul(mul(left, [[Q(c), Q(s)], [Q(-s), Q(c)]]), right)


def monomial(p, phases):
    result = [[Z for _ in p] for _ in p]
    for j, i in enumerate(p):
        result[i][j] = phases[i]
    return result


def fixtures():
    specifications = [
        ("singleton", 1, [], False),
        ("complex-two-coordinate", 2, [(0, 1)], False),
        ("disjoint-pairs", 6, [(0, 1), (2, 3), (4, 5)], False),
        ("overlapping-path", 4, [(1, 2), (0, 1), (2, 3)], False),
        ("overlapping-star", 5, [(0, 3), (0, 1), (0, 4), (0, 2)], False),
        ("two-components", 6, [(0, 2), (3, 5), (1, 2), (4, 5)], False),
        ("permuted-phased-tree", 6, [(1, 4), (0, 2), (2, 5), (2, 4), (3, 4)], True),
        ("endpoint-gates", 4, [(0, 1), (1, 2), (2, 3)], False),
        ("repeated-zero-spectra", 4, [(0, 1), (1, 2), (1, 3)], False),
        ("singular-sum-point-hull", 4, [(0, 1), (1, 2), (2, 3)], False),
        ("real-segment-hull", 3, [(1, 2), (0, 1)], False),
    ]
    for name, n, edges, wrapped in specifications:
        gates = [(i, j, phased_block(variant=k)) for k, (i, j) in enumerate(edges)]
        if name == "endpoint-gates":
            gates[0] = (0, 1, phased_block(F(1), F(0)))
            gates[1] = (1, 2, phased_block(F(0), F(1), 1))
        a = [Q(i - 2, (i % 3) - 1) for i in range(n)]
        b = [Q(3 - i, (2 * i % 5) - 2) for i in range(n)]
        if name == "repeated-zero-spectra":
            a, b = [Z, Z, Q(2, 1), Q(2, 1)], [Z, Q(1, -1), Z, Q(1, -1)]
        if name == "singular-sum-point-hull":
            a, b = [Q(2, 1)] * n, [Q(-2, -1)] * n
        if name == "real-segment-hull":
            a, b = [Q(0), Q(2), Q(5)], [Q(1), Q(3), Q(7)]
        p = tuple(reversed(range(n))) if wrapped else tuple(range(n))
        r = tuple((i + 1) % n for i in range(n)) if wrapped else tuple(range(n))
        yield {"name": name, "n": n, "gates": gates, "a": a, "b": b,
               "leftPermutation": p, "rightPermutation": r, "wrapped": wrapped}


def check_fixture(fixture):
    n, gates = fixture["n"], fixture["gates"]
    unitary, weights = build_unitary(n, gates), permutation_weights(n, gates)
    left, right = fixture["leftPermutation"], fixture["rightPermutation"]
    if fixture["wrapped"]:
        phases = [Q(0, 1) if i % 2 else Q(F(3, 5), F(4, 5)) for i in range(n)]
        unitary = mul(mul(monomial(left, phases), unitary), monomial(right, list(reversed(phases))))
        weights = {tuple(left[p[right[j]]] for j in range(n)): w for p, w in weights.items()}
    if mul(unitary, adj(unitary)) != eye(n) or mul(adj(unitary), unitary) != eye(n):
        raise AssertionError("Unitary identity failed")
    if sum(weights.values(), F(0)) != 1 or any(w < 0 for w in weights.values()):
        raise AssertionError("Weights are not a probability distribution")
    checked = 0
    for k in range(n + 1):
        subsets = list(combinations(range(n), k))
        for rows in subsets:
            for columns in subsets:
                minor = det([[unitary[i][j] for j in columns] for i in rows])
                predicted = sum((w for p, w in weights.items()
                                 if tuple(sorted(p[j] for j in columns)) == rows), F(0))
                if norm2(minor) != predicted:
                    raise AssertionError((fixture["name"], rows, columns, norm2(minor), predicted))
                checked += 1
    a, b = fixture["a"], fixture["b"]
    orbit = mul(mul(unitary, diag(b)), adj(unitary))
    summed = [[orbit[i][j] + (a[i] if i == j else Z) for j in range(n)] for i in range(n)]
    actual = det(summed)
    expected = sum((w * vertex(a, b, p) for p, w in weights.items()), Z)
    if actual != expected or actual != determinant_leibniz(summed):
        raise AssertionError("Determinant convex combination or second algorithm failed")
    return {
        "name": fixture["name"], "dimension": n,
        "gatesInApplicationOrder": [{"edgeZeroBased": [i, j], "block": [[z.encode() for z in row] for row in block]} for i, j, block in gates],
        "leftPermutationColumnToRowZeroBased": list(left),
        "rightPermutationColumnToRowZeroBased": list(right),
        "a": [z.encode() for z in a], "b": [z.encode() for z in b],
        "U": [[z.encode() for z in row] for row in unitary],
        "weights": [{"columnToRowZeroBased": list(p), "weight": str(w)} for p, w in sorted(weights.items())],
        "compoundEntriesChecked": checked,
        "determinant": actual.encode(), "allChecksPassed": True,
    }


def cycle_obstruction():
    block = [[Q(F(3, 5)), Q(F(4, 5))], [Q(F(-4, 5)), Q(F(3, 5))]]
    gates = [(0, 1, block), (1, 2, block), (0, 2, block)]
    unitary = build_unitary(3, gates, require_forest=False)
    weights = permutation_weights(3, gates, require_forest=False)
    actual = norm2(unitary[0][0])
    predicted = sum((w for p, w in weights.items() if p[0] == 0), F(0))
    assert unitary[0][0] == Q(F(109, 125))
    assert actual == F(11881, 15625) and predicted == F(6121, 15625)
    # Isolates this compound coefficient in d(U): a=(0,1,1), b=(1,0,0).
    a, b = [Z, O, O], [O, Z, Z]
    orbit = mul(mul(unitary, diag(b)), adj(unitary))
    summed = [[orbit[i][j] + (a[i] if i == j else Z) for j in range(3)] for i in range(3)]
    assert det(summed) == Q(actual)
    assert sum((w * vertex(a, b, p) for p, w in weights.items()), Z) == Q(predicted)
    return {"scope": "These product weights fail; this is not a counterexample to MOC.",
            "gateEdgesInApplicationOrderZeroBased": [[0, 1], [1, 2], [0, 2]],
            "rotationC": "3/5", "rotationS": "4/5", "U11": unitary[0][0].encode(),
            "actualSquaredModulusAndDeterminant": str(actual), "predictedProductWeight": str(predicted),
            "difference": str(actual - predicted),
            "spectra": {"a": [z.encode() for z in a], "b": [z.encode() for z in b]},
            "originalHull": ["0", "1"], "originalHullMembership": "inside"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    # POSIX CPU guard applies to this finite suite, never to external model work.
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    started = time.process_time()
    cases = [check_fixture(fixture) for fixture in fixtures()]
    obstruction = cycle_obstruction()
    report = {"schemaVersion": 1, "arithmetic": "exact Gaussian rationals; no floating comparisons",
              "evidenceLevel": "same-owner finite regression checks, not formal verification",
              "pythonVersion": platform.python_version(), "seed": None,
              "bounds": {"fixedForestFixtures": len(cases), "cycleFixtures": 1,
                         "maximumDimension": 6, "maximumForestGates": 5, "cpuLimitSeconds": 60},
              "processCpuSeconds": round(time.process_time() - started, 6),
              "dependenciesSha256": {"forest.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                     "../certificates/verify.py": hashlib.sha256(VERIFY.read_bytes()).hexdigest()},
              "compoundEntriesChecked": sum(case["compoundEntriesChecked"] for case in cases),
              "fixtures": cases, "cycleObstruction": obstruction}
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
        print(json.dumps({"output": str(args.output), "fixtures": len(cases),
                          "compoundEntriesChecked": report["compoundEntriesChecked"],
                          "processCpuSeconds": report["processCpuSeconds"], "allChecksPassed": True}))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
