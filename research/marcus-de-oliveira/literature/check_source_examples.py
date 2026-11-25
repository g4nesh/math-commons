#!/usr/bin/env python3
"""Independent exact arithmetic checks of two sourced finite examples and a formula.

This verifies neither the full n=4 preprint nor the original MOC.
Only Python's standard library is used. No downloaded code is executed.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
import json


@dataclass(frozen=True)
class QComplex:
    re: F = F(0)
    im: F = F(0)

    def __add__(self, other):
        other = as_q(other)
        return QComplex(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return QComplex(-self.re, -self.im)

    def __sub__(self, other):
        return self + -as_q(other)

    def __mul__(self, other):
        other = as_q(other)
        return QComplex(self.re * other.re - self.im * other.im,
                        self.re * other.im + self.im * other.re)

    __rmul__ = __mul__

    def conjugate(self):
        return QComplex(self.re, -self.im)


def as_q(value):
    return value if isinstance(value, QComplex) else QComplex(F(value))


def determinant(matrix):
    total = as_q(0)
    for p in permutations(range(len(matrix))):
        inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
        term = as_q((-1) ** inversions)
        for i, j in enumerate(p):
            term *= matrix[i][j]
        total += term
    return total


def gram(matrix):
    return [[sum(as_q(x).conjugate() * y for x, y in zip(row_i, row_j))
             for row_j in matrix] for row_i in matrix]


def assert_orthogonal(matrix):
    n = len(matrix)
    assert gram(matrix) == [[as_q(int(i == j)) for j in range(n)] for i in range(n)]


def minor_squared(matrix, rows, columns):
    d = determinant([[matrix[i][j] for j in columns] for i in rows])
    return d.re**2 + d.im**2


def psi(matrix):
    # Indices here are zero-based; the source uses 12,13,24,34.
    return (sum(F(matrix[i][i])**2 for i in range(4))
            - minor_squared(matrix, (0, 1), (0, 2))
            - minor_squared(matrix, (0, 1), (1, 3))
            + minor_squared(matrix, (0, 2), (0, 1))
            + minor_squared(matrix, (0, 2), (2, 3))) / 2


def common_weights_obstruction():
    matrix = [[F(x, 7) for x in row] for row in (
        (1, 4, 4, 4), (-4, 1, 4, -4), (-4, -4, 1, 4), (-4, 4, -4, 1))]
    assert_orthogonal(matrix)
    assert determinant(matrix) == as_q(1)
    values = []
    for p in permutations(range(4)):
        permutation_matrix = [[int(j == p[i]) for j in range(4)] for i in range(4)]
        values.append(psi(permutation_matrix))
    histogram = {str(value): values.count(value) for value in sorted(set(values))}
    assert histogram == {'0': 13, '1': 10, '2': 1}
    assert psi(matrix) == F(-158, 2401)
    return {'orthogonal_exactly': True, 'determinant': '1',
            'psi': str(psi(matrix)), 'permutation_psi_histogram': histogram,
            'scope': 'Obstruction to common permutation weights for all squared minors; not an MOC counterexample.'}


def goyal_hull_membership():
    a = [as_q(0), as_q(1), QComplex(F(0), F(1))]
    b = [QComplex(F(-1), F(-1)), as_q(1), QComplex(F(-1), F(1))]
    u = [[F(x, 3) for x in row] for row in ((-1, 2, 2), (2, -1, 2), (2, 2, -1))]
    assert_orthogonal(u)
    mixed = [[sum(u[i][k] * b[k] * u[j][k] for k in range(3))
              + (a[i] if i == j else 0) for j in range(3)] for i in range(3)]
    d = determinant(mixed)
    points = []
    for p in permutations(range(3)):
        z = as_q(1)
        for i in range(3):
            z *= a[i] + b[p[i]]
        points.append(z)
    expected = [(6, -2), (2, 0), (2, 1), (0, -1), (0, 2), (2, -2)]
    assert points == [QComplex(F(x), F(y)) for x, y in expected]
    assert d == QComplex(F(2, 3), F(2, 9))
    weights = [F(0), F(0), F(1, 3), F(13, 27), F(5, 27), F(0)]
    assert all(w >= 0 for w in weights) and sum(weights) == 1
    assert sum(w * z for w, z in zip(weights, points)) == d
    return {'unitary_exactly': True,
            'determinant': {'re': str(d.re), 'im': str(d.im)},
            'sigma_points_lexicographic': expected,
            'convex_weights_lexicographic': [str(w) for w in weights],
            'scope': 'Certifies original convex-hull membership for this witness; does not check the smaller polygon topology.'}


def permuted_prefactor_check():
    # Corollary 3, p2611 v2 p.4: identity Givens, i=1, j=2, pi=(3,2,1).
    # The outside product must also permute b. All indices here are zero-based.
    a, b, p = (1, 2, 3), (4, 5, 6), (2, 1, 0)
    actual = determinant([[a[i] + b[p[i]] if i == j else 0
                           for j in range(3)] for i in range(3)])
    inside = (a[0] + b[p[0]]) * (a[1] + b[p[1]])
    printed = (a[2] + b[2]) * inside
    corrected = (a[2] + b[p[2]]) * inside
    assert actual == as_q(corrected) == as_q(343)
    assert printed == 441 and printed != corrected
    return {'actual_determinant': '343', 'printed_unpermuted_prefactor': '441',
            'corrected_permuted_prefactor': '343',
            'scope': 'Local displayed-formula discrepancy; not a disproof of the segment result or MOC.'}


def main():
    print(json.dumps({'arithmetic': 'fractions.Fraction; exhaustive Leibniz determinants',
                      'common_weights_obstruction': common_weights_obstruction(),
                      'goyal_witness': goyal_hull_membership(),
                      'permuted_prefactor_check': permuted_prefactor_check(),
                      'general_conjecture_resolved': False}, indent=2))


if __name__ == '__main__':
    main()
