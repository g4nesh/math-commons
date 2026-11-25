#!/usr/bin/env python3
"""Independent exact re-check of the finite examples in literature/README.md.

Standard library only. Shares no code with literature/check_source_examples.py:
determinants use fraction-exact Gaussian elimination over Q(i) instead of the
Leibniz expansion, and hull membership of the Goyal witness is re-derived by
searching all sigma-point triangles for exact barycentric coordinates before
the author's weights are checked.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
import json
import os

ZERO, ONE = (F(0), F(0)), (F(1), F(0))


def q(x, y=0):
    return (F(x), F(y))


def add(u, v):
    return (u[0] + v[0], u[1] + v[1])


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def mul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def div(u, v):
    n = v[0] ** 2 + v[1] ** 2
    return ((u[0] * v[0] + u[1] * v[1]) / n, (u[1] * v[0] - u[0] * v[1]) / n)


def conj(u):
    return (u[0], -u[1])


def det(m):
    """Exact determinant over Q(i) by Gaussian elimination with pivot search."""
    a = [row[:] for row in m]
    n, d = len(a), ONE
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != ZERO), None)
        if p is None:
            return ZERO
        if p != c:
            a[c], a[p] = a[p], a[c]
            d = (-d[0], -d[1])
        d = mul(d, a[c][c])
        for r in range(c + 1, n):
            f = div(a[r][c], a[c][c])
            a[r] = [sub(x, mul(f, y)) for x, y in zip(a[r], a[c])]
    return d


def is_unitary(u):
    n = len(u)
    for i in range(n):
        for j in range(n):
            s = ZERO
            for k in range(n):
                s = add(s, mul(conj(u[k][i]), u[k][j]))
            if s != (ONE if i == j else ZERO):
                return False
    return True


def sq_minor(v, rows, cols):
    d = det([[v[i][j] for j in cols] for i in rows])
    return d[0] ** 2 + d[1] ** 2


def psi(v):
    # 2 psi(V) = sum |V_ii|^2 - |V[12,13]|^2 - |V[12,24]|^2 + |V[13,12]|^2 + |V[13,34]|^2
    s = sum(v[i][i][0] ** 2 + v[i][i][1] ** 2 for i in range(4))
    s += -sq_minor(v, (0, 1), (0, 2)) - sq_minor(v, (0, 1), (1, 3))
    s += sq_minor(v, (0, 2), (0, 1)) + sq_minor(v, (0, 2), (2, 3))
    return s / 2


def common_weight_obstruction():
    o = [[q(F(x, 7)) for x in row] for row in
         ((1, 4, 4, 4), (-4, 1, 4, -4), (-4, -4, 1, 4), (-4, 4, -4, 1))]
    assert is_unitary(o) and all(x[1] == 0 for r in o for x in r)
    assert det(o) == ONE
    hist = {}
    for p in permutations(range(4)):
        val = psi([[ONE if j == p[i] else ZERO for j in range(4)] for i in range(4)])
        hist[str(val)] = hist.get(str(val), 0) + 1
    assert hist == {'0': 13, '1': 10, '2': 1}
    assert psi(o) == F(-158, 2401)
    return {'orthogonal': True, 'det': '1', 'psi': str(psi(o)), 'perm_psi_hist': hist}


def bary(z, a, b, c):
    """Exact barycentric coordinates of z in triangle abc (None if degenerate)."""
    d = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
    if d == 0:
        return None
    l1 = ((z[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (z[1] - a[1])) / d
    l2 = ((b[0] - a[0]) * (z[1] - a[1]) - (z[0] - a[0]) * (b[1] - a[1])) / d
    return (1 - l1 - l2, l1, l2)


def goyal_witness():
    a = [q(0), q(1), q(0, 1)]
    b = [q(-1, -1), q(1), q(-1, 1)]
    u = [[q(F(x, 3)) for x in row] for row in ((-1, 2, 2), (2, -1, 2), (2, 2, -1))]
    assert is_unitary(u)
    m = []
    for i in range(3):
        row = []
        for j in range(3):
            s = a[i] if i == j else ZERO
            for k in range(3):
                s = add(s, mul(mul(u[i][k], b[k]), conj(u[j][k])))
            row.append(s)
        m.append(row)
    d = det(m)
    assert d == (F(2, 3), F(2, 9))
    pts = []
    for p in permutations(range(3)):  # lexicographic: 123,132,213,231,312,321
        z = ONE
        for i in range(3):
            z = mul(z, add(a[i], b[p[i]]))
        pts.append(z)
    assert pts == [q(6, -2), q(2), q(2, 1), q(0, -1), q(0, 2), q(2, -2)]
    found = [t for t in combinations(range(6), 3)
             if (w := bary(d, *(pts[k] for k in t))) and min(w) >= 0]
    assert found, 'determinant not in any sigma-point triangle'
    w = [F(0), F(0), F(1, 3), F(13, 27), F(5, 27), F(0)]
    combo = ZERO
    for wk, z in zip(w, pts):
        combo = add(combo, (wk * z[0], wk * z[1]))
    assert min(w) >= 0 and sum(w) == 1 and combo == d
    return {'unitary': True, 'det': ['2/3', '2/9'],
            'containing_triangles': [[k + 1 for k in t] for t in found],
            'author_weights_ok': True}


def corollary3_outside_factor():
    # G = I: the Givens bracket reduces to (a_i+b_pi(i))(a_j+b_pi(j)).
    a, b, p = (1, 2, 3), (4, 5, 6), (2, 1, 0)
    actual = det([[q(a[i] + b[p[i]]) if i == j else ZERO for j in range(3)] for i in range(3)])
    inner = (a[0] + b[p[0]]) * (a[1] + b[p[1]])
    assert actual == q(343) and (a[2] + b[p[2]]) * inner == 343 and (a[2] + b[2]) * inner == 441
    # Non-identity rotation c=3/5, s=4/5 in plane (1,2): my own expansion gives
    # bracket = (a1+B1)(a2+B2) + (a1-a2)(B1-B2) s^2 with B = b o pi.
    c, s = F(3, 5), F(4, 5)
    g = [[q(c), q(-s), ZERO], [q(s), q(c), ZERO], [ZERO, ZERO, ONE]]
    beta = [b[p[k]] for k in range(3)]
    m = [[q(a[i] if i == j else 0) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                m[i][j] = add(m[i][j], mul(mul(g[i][k], q(beta[k])), g[j][k]))
    bracket = (a[0] + beta[0]) * (a[1] + beta[1]) + (a[0] - a[1]) * (beta[0] - beta[1]) * s * s
    assert det(m) == q((a[2] + beta[2]) * bracket) != q((a[2] + b[2]) * bracket)
    return {'G=I': {'actual': 343, 'permuted_outside': 343, 'unpermuted_outside': 441},
            'givens_3_4_5': {'actual': str(det(m)[0]), 'unpermuted_outside': str((a[2] + b[2]) * bracket)}}


def main():
    out = {'common_weight_obstruction': common_weight_obstruction(),
           'goyal_witness': goyal_witness(),
           'corollary3_outside_factor': corollary3_outside_factor()}
    rec = 'research/marcus-de-oliveira/literature/source-example-results.json'
    if os.path.exists(rec):
        r = json.load(open(rec))
        assert r['common_weights_obstruction']['psi'] == out['common_weight_obstruction']['psi']
        assert r['goyal_witness']['determinant'] == {'re': '2/3', 'im': '2/9'}
        assert r['permuted_prefactor_check']['printed_unpermuted_prefactor'] == '441'
        out['recorded_results_agree'] = True
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
