"""Exact six-dimensional checks of the pole identity and lifted permutation masses."""
import itertools
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'certificates'))
import verify as v

Q = v.QComplex


def inverse(matrix):
    n = len(matrix)
    rows = [left[:] + right[:] for left, right in zip(matrix, v.identity(n))]
    for i in range(n):
        pivot = next(k for k in range(i, n) if rows[k][i] != v.ZERO)
        rows[i], rows[pivot] = rows[pivot], rows[i]
        scale = rows[i][i]
        rows[i] = [x / scale for x in rows[i]]
        for k in range(n):
            if k != i:
                factor = rows[k][i]
                rows[k] = [x-factor*y for x, y in zip(rows[k], rows[i])]
    return [r[n:] for r in rows]


def add_diag(matrix, values):
    return [[x + (values[i] if i == j else v.ZERO) for j, x in enumerate(row)]
            for i, row in enumerate(matrix)]


def prod(xs):
    return math.prod(xs, start=v.ONE)


def fixture(circular):
    # Shifted inverse with poles on rows 0,1 and four retained eigenvalues 2,4,6,8.
    g0 = v.diagonal(list(map(Q, [2, 4, 2, 4, 6, 8])))
    g0[0][2] = g0[2][0] = v.ONE
    g0[1][3] = g0[3][1] = Q(F(1, 2))
    scatter = v.identity(6)
    for i, j in itertools.product(range(2, 6), repeat=2):
        scatter[i][j] = Q(int(i == j) - F(1, 2))
    assert v.multiply(scatter, v.adjoint(scatter)) == v.identity(6)
    g = v.multiply(v.multiply(scatter, g0), v.adjoint(scatter))
    shifts = list(map(Q, [1, 3, F(7, 2), F(9, 2), 6, 8]))
    pole = Q(0, 1) if circular else v.ZERO
    if circular:
        g = add_diag(g, [Q(0, 1)] * 6)
        shifts = [x + Q(0, 1) for x in shifts]
    a = add_diag(inverse(g), [-pole] * 6)
    assert v.multiply(a, v.adjoint(a)) == v.multiply(v.adjoint(a), a)
    spectrum = [v.ONE / x - pole for x in shifts]
    return a, g, pole, spectrum


def check(circular):
    a, g, pole, spectrum = fixture(circular)
    retained = [r[2:] for r in g[2:]]
    hs = [Q(x, int(circular)) for x in (2, 4, 6, 8)]
    # Check the advertised eigenvalues without numerical diagonalization.
    for root in hs:
        assert v.determinant(add_diag(retained, [-root]*4)) == v.ZERO
    for root in spectrum:
        assert v.determinant(add_diag(a, [-root]*6)) == v.ZERO
    assert (hs[1]-hs[0])*(hs[3]-hs[2])/((hs[2]-hs[0])*(hs[3]-hs[1])) == Q(F(1, 4))

    exceptional = [Q(0, 2), Q(1, 2), Q(-2, 1), Q(3, -1)]
    b = [pole, pole] + exceptional
    beta = [v.ONE/(bj-pole) for bj in exceptional]
    multiplier = v.determinant(add_diag(a, [pole]*6)) * prod(bj-pole for bj in exceptional)
    direct = v.determinant(add_diag(a, b))
    assert direct == multiplier * v.determinant(add_diag(retained, beta))

    # The reduced matrix is the order-four flat Householder applied to diag(hs).
    # Its uniform eight oriented three-cycles lift by two binary injections.
    weighted = v.ZERO
    count = 0
    for perm in itertools.permutations(range(4)):
        if sum(i != j for i, j in enumerate(perm)) != 3:
            continue
        for u, w in itertools.product((0, 1), (2, 3)):
            embedding = [u, w, 4, 5]
            chosen = [embedding[i] for i in perm]
            unused = sorted(set(range(6)) - set(chosen))
            assert len(unused) == 2
            weighted += prod(spectrum[j] + bj for j, bj in zip(unused + chosen, b)) * F(1, 32)
            count += 1
    assert count == 32 and weighted == direct
    return {'circular': circular, 'dimension': 6, 'nonrealExceptions': 4,
            'liftedTerms': count, 'crossRatio': '1/4', 'determinant': direct.encode()}


if __name__ == '__main__':
    print(check(False))
    print(check(True))
