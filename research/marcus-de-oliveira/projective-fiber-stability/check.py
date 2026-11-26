"""Positive laws near four simultaneous spectral collisions. Standard library only."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "four-cycle-spectral"))
import spectral as sp

Q, ZERO, ONE = sp.Q, sp.Z, sp.O
PERMS = tuple(permutations(range(4)))
POLES = tuple(-b for b in sp.B_SPECTRUM)
FLOOR = F(1, 2000)


@lru_cache(None)
def initial_gaps():
    """Coefficient of each square-free spectral monomial in vertex minus determinant."""
    orbit = sp.matrix_sum([ZERO] * 4)
    columns = [dict() for _ in PERMS]
    for size in range(5):
        for subset in combinations(range(4), size):
            subset = frozenset(subset)
            outside = sorted(set(range(4)) - subset)
            minor = sp.leibniz([[orbit[i][j] for j in outside] for i in outside])
            for col, p in zip(columns, PERMS):
                col[subset] = sp.product(b for b, i in zip(sp.B_SPECTRUM, p) if i not in subset) - minor
    return columns


def add_term(poly, key, coefficient):
    poly[key] = poly.get(key, ZERO) + coefficient


def translated(column, subset):
    """Shift every spectral entry by delta before introducing free deviations."""
    subset = frozenset(subset)
    out = {}
    for larger, value in column.items():
        if subset <= larger and value != ZERO:
            add_term(out, (len(larger) - len(subset), 0), value)
    return out


def pencil(column, subset, pole):
    out = translated(column, subset)
    for (d, e), value in translated(column, (*subset, 1)).items():
        add_term(out, (d, e), pole * value)
        add_term(out, (d, e + 1), value)
        add_term(out, (d + 1, e), -value)
    return out


def powers(value, exponent):
    return sp.product([value] * exponent)


def evaluate(poly, delta, eta):
    return sum((v * powers(delta, d) * powers(eta, e) for (d, e), v in poly.items()), ZERO)


@lru_cache(None)
def system(pole):
    # A row is a real or imaginary part of a Gaussian-rational polynomial.
    # The fourth assigned column marginal is implied by the other three and mass one.
    yes = [{(0, 0): ONE} for _ in PERMS]
    rows = [("real", yes)]
    for j in range(3):
        rows.append(("real", [{(0, 0): Q(int(p[j] == 3)) - sp.U[3][j] * sp.U[3][j]}
                              for p in PERMS]))
    for subset in ((1,), (0,), (0, 2)):
        entries = [translated(col, subset) if subset == (1,) else pencil(col, subset, pole)
                   for col in initial_gaps()]
        rows.extend((part, entries) for part in ("real", "imag"))
    return rows


def matrix(rows, delta=ZERO, eta=ZERO):
    return [[getattr(evaluate(poly, delta, eta), part) for poly in entries] for part, entries in rows]


def solve_weights(basis, mat):
    if len(basis) != len(mat) or len(set(basis)) != len(basis) or any(type(i) is not int or not 0 <= i < 24 for i in basis):
        raise ValueError("Invalid square basis")
    fixed = sorted(set(range(24)) - set(basis))
    square = [[row[j] for j in basis] for row in mat]
    rhs = [int(i == 0) - FLOOR * sum(row[j] for j in fixed) for i, row in enumerate(mat)]
    result = [FLOOR] * 24
    for j, mass in zip(basis, sp.exact_solve(square, rhs)):
        result[j] = mass
    if any(sum(a * b for a, b in zip(row, result)) != int(i == 0) for i, row in enumerate(mat)):
        raise ArithmeticError("A probability equation failed")
    return result


def real_rank(mat):
    rows = [row.copy() for row in mat]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank + 1, len(rows)):
            factor = rows[i][col] / rows[rank][col]
            rows[i] = [x - factor * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def regularized_rank(pole):
    rows = [("real", [{(0, 0): ONE} for _ in PERMS])]
    for subset in ((1,), (0,), (2,), (0, 2)):
        entries = [translated(col, subset) if subset == (1,) else pencil(col, subset, pole)
                   for col in initial_gaps()]
        rows.extend((part, entries) for part in ("real", "imag"))
    return real_rank(matrix(rows))


def l1(z):
    return abs(z.real) + abs(z.imag)


def neighborhood_bounds(pole, basis, radius):
    rows = system(pole)
    mat = matrix(rows)
    masses = solve_weights(basis, mat)
    if min(masses) <= 0:
        raise ValueError("Center law is not strictly positive")
    n = len(rows)
    block = [[row[j] for j in basis] for row in mat]
    inverse_columns = [sp.exact_solve(block, [int(i == k) for i in range(n)]) for k in range(n)]
    inverse = list(zip(*inverse_columns))
    variation = []
    for _, entries in rows:
        variation.append([sum(l1(v) * (2 * radius) ** (d + e) for (d, e), v in poly.items() if d + e)
                          for poly in entries])
    rho = nu = F(0)
    for index, inverse_row in enumerate(inverse):
        bounds = [sum(abs(q) * variation[k][j] for k, q in enumerate(inverse_row)) / masses[basis[index]]
                  for j in range(24)]
        rho = max(rho, sum(bounds[j] * masses[j] for j in basis))
        nu = max(nu, sum(b * m for b, m in zip(bounds, masses)))
    return rho, nu, min(masses)


def weights_at(pole_index, delta, eta, data):
    radius = F(data["radius"])
    if max(abs(q) for z in (delta, eta) for q in (z.real, z.imag)) > radius:
        raise ValueError("Outside certified square")
    return solve_weights(data["bases"][pole_index], matrix(system(POLES[pole_index]), delta, eta))


def verify(data):
    radius = F(data["radius"])
    if radius <= 0 or len(data["bases"]) != len(POLES):
        raise ValueError("Expected a positive radius and four bases")
    if sp.e.multiply(sp.U, sp.e.adjoint(sp.U)) != sp.e.identity(4):
        raise ArithmeticError("Fixture is not unitary")
    for i, (pole, basis) in enumerate(zip(POLES, data["bases"])):
        if regularized_rank(pole) != 8:
            raise ArithmeticError("Unexpected collision rank")
        rho, nu, minimum = neighborhood_bounds(pole, basis, radius)
        if rho >= 1 or nu / (1 - rho) >= F(1, 2):
            raise ArithmeticError("Perturbation bound does not ensure positive masses")
        delta, eta = Q(radius, -radius), Q(-radius, radius)
        masses = weights_at(i, delta, eta, data)
        if min(masses) <= 0:
            raise ArithmeticError("Corner law is not positive")
        for x, y in ((ZERO, ONE), (Q(-2, 3), Q(5, -1)), (pole + eta, pole + eta)):
            a = [x, pole + eta, y, delta]
            predicted = sum((m * sp.vertex(a, p) for m, p in zip(masses, PERMS)), ZERO)
            if predicted != sp.e.determinant(sp.matrix_sum(a)):
                raise ArithmeticError("Direct determinant failed")
        print(f"pole {i + 1}: rho < 1, relative displacement < 1/2, center floor {minimum}")
    print(f"All four joint neighborhoods verified at radius {radius}.")


if __name__ == "__main__":
    verify(json.loads((HERE / "bases.json").read_text()))
