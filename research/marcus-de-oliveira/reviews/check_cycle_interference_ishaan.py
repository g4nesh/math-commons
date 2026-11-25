#!/usr/bin/env python3
"""Independent exact re-check of the single-cycle study (task mdo-cycle-interference).

Standard library only. Imports no repository module; reads the two published
JSON files only to compare stored data with fresh recomputation. Permutations
are column-to-row tuples p (P_p e_j = e_{p(j)}), as in the reviewed README.
Run from the repository root; exits nonzero on the first failed assertion.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / 'single-cycle'


class G:  # exact Gaussian rational
    def __init__(self, re, im=0):
        self.re, self.im = F(re), F(im)
    def __add__(s, o): o = o if isinstance(o, G) else G(o); return G(s.re + o.re, s.im + o.im)
    __radd__ = __add__
    def __neg__(s): return G(-s.re, -s.im)
    def __sub__(s, o): return s + (-(o if isinstance(o, G) else G(o)))
    def __mul__(s, o):
        o = o if isinstance(o, G) else G(o)
        return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    __rmul__ = __mul__
    def conj(s): return G(s.re, -s.im)
    def abs2(s): return s.re * s.re + s.im * s.im
    def __truediv__(s, o): n = o.abs2(); t = s * o.conj(); return G(t.re / n, t.im / n)
    def __eq__(s, o): o = o if isinstance(o, G) else G(o); return s.re == o.re and s.im == o.im
    def __bool__(s): return bool(s.re or s.im)


def eye(n): return [[G(int(i == j)) for j in range(n)] for i in range(n)]
def mul(a, b): return [[sum((a[i][k] * b[k][j] for k in range(len(b))), G(0)) for j in range(len(b[0]))] for i in range(len(a))]
def adj(a): return [[a[j][i].conj() for j in range(len(a))] for i in range(len(a[0]))]
def sub(u, s, t): return [[u[i][j] for j in t] for i in s]


def det(m):
    m, d = [r[:] for r in m], G(1)
    for k in range(len(m)):
        p = next((r for r in range(k, len(m)) if m[r][k]), None)
        if p is None: return G(0)
        if p != k: m[k], m[p], d = m[p], m[k], -d
        d = d * m[k][k]
        for r in range(k + 1, len(m)):
            f = m[r][k] / m[k][k]
            m[r] = [m[r][c] - f * m[k][c] for c in range(len(m))]
    return d


def inverse(m):
    n = len(m); a = [m[i][:] + eye(n)[i] for i in range(n)]
    for k in range(n):
        p = next(r for r in range(k, n) if a[r][k]); a[k], a[p] = a[p], a[k]
        a[k] = [v / a[k][k] for v in a[k]]
        for r in range(n):
            if r != k and a[r][k]: f = a[r][k]; a[r] = [a[r][c] - f * a[k][c] for c in range(2 * n)]
    return [row[n:] for row in a]


def cayley(h):  # (I - iH)(I + iH)^{-1}: exact unitary for Hermitian rational H
    n = len(h); ih = [[G(0, 1) * h[i][j] for j in range(n)] for i in range(n)]
    u = mul([[eye(n)[i][j] - ih[i][j] for j in range(n)] for i in range(n)],
            inverse([[eye(n)[i][j] + ih[i][j] for j in range(n)] for i in range(n)]))
    assert mul(u, adj(u)) == eye(n) and mul(adj(u), u) == eye(n)
    return u


def gate(n, i, j, c, s):
    g = eye(n); g[i][i], g[i][j], g[j][i], g[j][j] = G(c), G(s), G(-s), G(c); return g


def product(n, gates):  # gates in application order: U = G_m ... G_1
    u = eye(n)
    for i, j, c, s in gates: u = mul(gate(n, i, j, c, s), u)
    return u


def hits(p, s, t): return tuple(sorted(p[j] for j in t)) == s


def cert_check(u, w):  # all compound layers, all sizes
    n = len(u); assert sum(w.values()) == 1 and min(w.values()) >= 0
    count = 0
    for k in range(n + 1):
        subs = list(combinations(range(n), k))
        for s in subs:
            for t in subs:
                if det(sub(u, s, t)).abs2() != sum((v for p, v in w.items() if hits(p, s, t)), F(0)):
                    return count, (s, t)
                count += 1
    return count, None


def tri(d, x=None):  # the README's six-weight table, table order = permutations(range(3))
    lo = max(F(0), d[2][2] - d[0][1], d[1][1] - d[0][2]); hi = min(d[0][0], d[1][1], d[2][2])
    x = lo if x is None else x
    vals = (x, d[0][0] - x, d[2][2] - x, d[0][2] - d[1][1] + x, d[0][1] - d[2][2] + x, d[1][1] - x)
    return dict(zip(permutations(range(3)), vals)), lo, hi


def two(u): return {(0, 1): u[0][0].abs2(), (1, 0): u[0][1].abs2()}
def sq(u): return [[z.abs2() for z in row] for row in u]


def embed(n, sup, block, w):
    u = eye(n)
    for a, x in enumerate(sup):
        for b, y in enumerate(sup): u[x][y] = block[a][b]
    out = {}
    for p, v in w.items():
        q = list(range(n))
        for j, x in enumerate(sup): q[x] = sup[p[j]]
        out[tuple(q)] = v
    return u, out


def conv(x, y):  # certificate of X*Y: independent p ~ w_X, q ~ w_Y, use p o q
    (u, wu), (v, wv) = x, y; out = {}
    for p, a in wu.items():
        for q, b in wv.items():
            r = tuple(p[q[j]] for j in range(len(q))); out[r] = out.get(r, F(0)) + a * b
    return mul(u, v), out


def det_identity(u, w, a, b):  # eq. (4): det(D_a + U D_b U*) = sum_p w_p prod_j (a_p(j) + b_j)
    n = len(u); m = mul(mul(u, [[b[i] if i == j else G(0) for j in range(n)] for i in range(n)]), adj(u))
    lhs = det([[m[i][j] + (a[i] if i == j else G(0)) for j in range(n)] for i in range(n)])
    rhs = G(0)
    for p, v in w.items():
        t = G(1)
        for j in range(n): t = t * (a[p[j]] + b[j])
        rhs = rhs + v * t
    return lhs, rhs


def rank(rows):
    rows, r = [list(map(F, x)) for x in rows], 0
    for c in range(len(rows[0])):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]: f = rows[i][c] / rows[r][c]; rows[i] = [rows[i][k] - f * rows[r][k] for k in range(len(rows[0]))]
        r += 1
    return r


out, H5 = {}, (F(3, 5), F(4, 5))
S3 = list(permutations(range(3)))
# 1. Triangle table: marginal map S3 -> 3x3 has rank 5, kernel = sign vector = x-coefficients.
M = [[int(p[j] == i) for p in S3] for i in range(3) for j in range(3)]
sign = [1, -1, -1, 1, 1, -1]
assert rank(M) == 5 and all(sum(r[k] * sign[k] for k in range(6)) == 0 for r in M)
H = [[[G(1), G(1, 1), G(0)], [G(1, -1), G(0), G(2)], [G(0), G(2), G(-1)]],
     [[G(2), G(0, 3), G(1, -1)], [G(0, -3), G(-1), G(F(1, 2))], [G(1, 1), G(F(1, 2)), G(0)]],
     [[G(0), G(5, 2), G(-3)], [G(5, -2), G(1), G(0, 1)], [G(-3), G(0, -1), G(3)]]]
tri_checks = 0
for h in H:
    u = cayley(h)
    w0, lo, hi = tri(sq(u)); assert lo <= hi
    for x in (lo, (lo + hi) / 2, hi):
        w = tri(sq(u), x)[0]
        n_, bad = cert_check(u, w); assert bad is None and n_ == 20; tri_checks += 1
        for a, b in (([G(1, 2), G(-1), G(0, 3)], [G(2, -1), G(1, 1), G(-3)]), ([G(0)] * 2 + [G(1)], [G(1)] + [G(0)] * 2)):
            lhs, rhs = det_identity(u, w, a, b); assert lhs == rhs
    w_any = {p: F(k + 1, 21) for k, p in enumerate(S3)}  # converse: a distribution is recovered from x = w(id)
    marg = [[sum(v for p, v in w_any.items() if p[j] == i) for j in range(3)] for i in range(3)]
    assert tri(marg, w_any[(0, 1, 2)])[0] == w_any
out['triangleCayleyChecks'] = tri_checks
# 2. Baseline triangle U = R13 R23 R12 (zero-based gates 01, 12, 02 in application order).
U3 = product(3, [(0, 1) + H5, (1, 2) + H5, (0, 2) + H5])
w3, lo, hi = tri(sq(U3))
assert (lo, hi) == (F(1881, 15625), F(2025, 15625))
assert [w3[p] for p in S3] == [F(k, 15625) for k in (1881, 10000, 144, 3456, 0, 144)]
assert cert_check(U3, w3) == (20, None) and U3[0][0].abs2() == F(11881, 15625) == w3[(0, 1, 2)] + w3[(0, 2, 1)]
old = {(0, 1, 2): F(1)}  # independent stay/swap rule of the forest baseline
for i, j in ((0, 1), (1, 2), (0, 2)):
    nxt = {}
    for p, v in old.items():
        sw = tuple(j if k == i else i if k == j else k for k in p)
        nxt[p] = nxt.get(p, F(0)) + v * F(9, 25); nxt[sw] = nxt.get(sw, F(0)) + v * F(16, 25)
    old = nxt
assert sum(v for p, v in old.items() if p[0] == 0) == F(6121, 15625)
# 3. Gluing: n=6 triangle core + pendant blocks on {2,3} (left), {0,4} (right), {3,5} (left); n=5 two triangles on pivot 2.
B2 = [cayley([[G(1), G(2, -1)], [G(2, 1), G(-2)]]), cayley([[G(0), G(1, 3)], [G(1, -3), G(1)]]), cayley([[G(3), G(0, 1)], [G(0, -1), G(0)]])]
T = [cayley(h) for h in H]
core = embed(6, [0, 1, 2], T[0], tri(sq(T[0]))[0])
for sup, blk, left in (([2, 3], B2[0], True), ([0, 4], B2[1], False), ([3, 5], B2[2], True)):
    g = embed(6, sup, blk, two(blk)); core = conv(g, core) if left else conv(core, g)
five = conv(embed(5, [0, 1, 2], T[1], tri(sq(T[1]))[0]), embed(5, [2, 3, 4], T[2], tri(sq(T[2]))[0]))
glue_counts = []
for u, w in (core, five):
    n_, bad = cert_check(u, w); assert bad is None; glue_counts.append(n_)
    lhs, rhs = det_identity(u, w, [G(k - 2, k % 3 - 1) for k in range(len(u))], [G(3 - k, 2 * k % 5 - 2) for k in range(len(u))]); assert lhs == rhs
assert glue_counts == [924, 252]
bad2 = conv(embed(4, [0, 1, 2], T[0], tri(sq(T[0]))[0]), embed(4, [1, 2, 3], T[1], tri(sq(T[1]))[0]))
assert cert_check(*bad2)[1] is not None  # negative control: overlap of two coordinates
out['gluingEntries'] = glue_counts
# 4. Four-cycle U = R14 R34 R23 R12 and its 13-pair cover.
U4 = product(4, [(0, 1) + H5, (1, 2) + H5, (2, 3) + H5, (0, 3) + H5])
assert [[625 * z.re for z in r] for r in U4] == [[-31, 492, -240, 300], [-300, 225, 500, 0], [240, -180, 225, 500], [-492, -256, -180, 225]]
assert all(z.im == 0 for r in U4 for z in r) and mul(U4, adj(U4)) == eye(4)
COVER = [((0,), (0,)), ((1,), (3,)), ((2,), (1,)), ((3,), (2,)), ((0, 1), (0, 2)), ((0, 1), (0, 3)), ((0, 1), (1, 3)),
         ((0, 2), (0, 1)), ((0, 2), (0, 3)), ((0, 2), (1, 2)), ((0, 3), (0, 1)), ((0, 3), (0, 2)), ((0, 3), (2, 3))]
S4 = list(permutations(range(4)))
h = [sum(hits(p, s, t) for s, t in COVER) for p in S4]
assert h == [1, 4, 4, 7, 4, 4, 1, 1, 4, 4, 1, 1, 1, 4, 1, 4, 1, 1, 1, 4, 1, 1, 4, 1]
val = sum(det(sub(U4, s, t)).abs2() for s, t in COVER); assert val == F(14953, 15625)
stored = json.loads((HERE / 'four-cycle-obstruction.json').read_text())
assert stored['U'] == [[str(z.re) for z in r] for r in U4]
assert [(tuple(c['rows']), tuple(c['columns']), F(c['squaredMinor'])) for c in stored['cover']] == [(s, t, det(sub(U4, s, t)).abs2()) for s, t in COVER]
r = [[z.abs2() for z in row] for row in U4]; greedy = {}  # lexicographic exact Birkhoff decomposition
while any(v for row in r for v in row):
    p = next(p for p in S4 if all(r[p[j]][j] > 0 for j in range(4))); m = min(r[p[j]][j] for j in range(4))
    greedy[p] = greedy.get(p, F(0)) + m
    for j in range(4): r[p[j]][j] -= m
assert all(sum(v for p, v in greedy.items() if p[j] == i) == U4[i][j].abs2() for i in range(4) for j in range(4))
assert det(sub(U4, (0, 1), (0, 1))).abs2() == F(81, 625) and sum(v for p, v in greedy.items() if hits(p, (0, 1), (0, 1))) == F(3321, 15625)
a = [G(-2, -1), G(-1), G(0, 1), G(1, -1)]; b = [G(3, -2), G(2), G(1, 2), G(0, -1)]
actual, weighted = det_identity(U4, greedy, a, b)
assert actual == G(F(269802, 390625), F(-5707764, 390625)) and weighted == G(F(1040842, 390625), F(-5724244, 390625))
res = json.loads((HERE / 'results.json').read_text())['fourCycleGreedyRuleFailure']['determinant']
assert [G(*map(F, z)) for z in res['a']] == a and [G(*map(F, z)) for z in res['b']] == b
hull, tot = G(0), F(0)
for e in res['originalHullConvexWeights']:
    p, wt = tuple(e['columnToRow']), F(e['weight']); pt = G(1)
    for j in range(4): pt = pt * (a[p[j]] + b[j])
    assert wt >= 0 and pt == G(*map(F, e['point'])); hull, tot = hull + wt * pt, tot + wt
assert tot == 1 and hull == actual
out['fourCycle'] = {'coverValue': str(val), 'gap': str(1 - val), 'minHits': min(h), 'greedyMiddleEntry': '3321/15625'}
# 5. Interleaved schedule on the triangle+pendant graph: 01,12,23,02 has no certificate; 01,12,02,23 has one.
H5r = (F(4, 5), F(3, 5))
Ui = product(4, [(0, 1) + H5, (1, 2) + H5, (2, 3) + H5, (0, 2) + H5r])
assert mul(Ui, adj(Ui)) == eye(4)
CI = [((1,), (3,)), ((0, 2), (2, 3)), ((1, 2), (0, 2)), ((1, 2), (0, 3)), ((1, 2), (1, 3)), ((1, 3), (0, 1)), ((1, 3), (1, 3)),
      ((1, 3), (2, 3)), ((2, 3), (0, 1)), ((2, 3), (0, 3)), ((2, 3), (1, 2)), ((0, 1, 2), (0, 2, 3)), ((0, 1, 3), (1, 2, 3)),
      ((0, 2, 3), (0, 1, 2)), ((1, 2, 3), (0, 1, 3))]
hi_ = [sum(hits(p, s, t) for s, t in CI) for p in S4]; vi = sum(det(sub(Ui, s, t)).abs2() for s, t in CI)
assert min(hi_) >= 1 and vi == F(14356, 15625)
Ut = product(3, [(0, 1) + H5, (1, 2) + H5, (0, 2) + H5r])
rec = conv(embed(4, [2, 3], gate(2, 0, 1, *H5), two(gate(2, 0, 1, *H5))), embed(4, [0, 1, 2], Ut, tri(sq(Ut))[0]))
assert rec[0] == product(4, [(0, 1) + H5, (1, 2) + H5, (0, 2) + H5r, (2, 3) + H5]) and cert_check(*rec) == (70, None)
out['interleaved'] = {'schedule': '01,12,23,02 (last gate c=4/5,s=3/5)', 'coverPairs': len(CI), 'minHits': min(hi_),
                      'coverValue': str(vi), 'gap': str(1 - vi), 'recursiveSchedule01_12_02_23': 'certified, 70 entries'}
print(json.dumps(out))
