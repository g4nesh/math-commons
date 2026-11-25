#!/usr/bin/env python3
"""Independent exact re-check of the four-cycle spectral certificates
(task mdo-four-cycle-spectral-weights).

Standard library only. Imports no repository module; reads the published JSON
files only to compare stored data with a fresh recomputation. With U and b
fixed, det(D_a + U D_b U^*) and every pairing value z_p(a) = prod_j (a_{p(j)}
+ b_j) are expanded symbolically (Leibniz) as multi-affine polynomials in the
four eigenvalues a0..a3 of A; every claim is then a coefficient comparison.
Permutations p are column-to-row tuples. Run from the repository root; exits
nonzero on the first failed assertion.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
SPEC = BASE / 'four-cycle-spectral'


class G:  # exact Gaussian rational
    def __init__(s, re, im=0): s.re, s.im = F(re), F(im)
    def __add__(s, o): o = o if isinstance(o, G) else G(o); return G(s.re + o.re, s.im + o.im)
    __radd__ = __add__
    def __neg__(s): return G(-s.re, -s.im)
    def __sub__(s, o): return s + (-(o if isinstance(o, G) else G(o)))
    def __mul__(s, o):
        o = o if isinstance(o, G) else G(o)
        return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    __rmul__ = __mul__
    def __eq__(s, o): o = o if isinstance(o, G) else G(o); return s.re == o.re and s.im == o.im
    def conj(s): return G(s.re, -s.im)


def g(pair): return G(F(pair[0]), F(pair[1]))
def cross(u, v): return u.re * v.im - u.im * v.re


R = [(-31, 492, -240, 300), (-300, 225, 500, 0), (240, -180, 225, 500), (-492, -256, -180, 225)]
U = [[G(F(x, 625)) for x in row] for row in R]
Ustar = [[U[j][i].conj() for j in range(4)] for i in range(4)]
b = [G(3, -2), G(2), G(1, 2), G(0, -1)]
a = [G(-2, -1), G(-1), G(0, 1), G(1, -1)]
I = G(0, 1)
mm = lambda X, Y: [[sum((X[i][k] * Y[k][j] for k in range(4)), G(0)) for j in range(4)] for i in range(4)]
Id = [[G(int(i == j)) for j in range(4)] for i in range(4)]
assert mm(U, Ustar) == Id and mm(Ustar, U) == Id, 'U is not unitary'
B = mm(mm(U, [[b[i] if i == j else G(0) for j in range(4)] for i in range(4)]), Ustar)
assert mm(B, [[B[j][i].conj() for j in range(4)] for i in range(4)]) == mm([[B[j][i].conj() for j in range(4)] for i in range(4)], B)
assert B[0][1] != G(0)  # [D_a,B]_{01} = (a0-a1) B01, and Re(a0(t)-a1) = -1 for every real t
src = json.loads((BASE / 'single-cycle' / 'results.json').read_text())['fourCycleGreedyRuleFailure']
assert [[g(x) for x in row] for row in src['U']] == U
assert [g(x) for x in src['determinant']['a']] == a and [g(x) for x in src['determinant']['b']] == b

# ---- multi-affine polynomials: {frozenset(vars): G}
def pmul(p, q):
    out = {}
    for s, x in p.items():
        for t, y in q.items():
            assert not s & t
            out[s | t] = out.get(s | t, G(0)) + x * y
    return out
def padd(p, q, c=1):
    out = dict(p)
    for k, v in q.items(): out[k] = out.get(k, G(0)) + c * v
    return {k: v for k, v in out.items() if v != G(0)}
def subst(p, fixed):  # substitute fixed values {var: G}
    out, fk = {}, frozenset(fixed)
    for s, x in p.items():
        for v in s & fk: x = x * fixed[v]
        k = s - fk; out[k] = out.get(k, G(0)) + x
    return {k: v for k, v in out.items() if v != G(0)}
def sign(p): return (-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))

DET = {}
for s in permutations(range(4)):
    term = {frozenset(): G(sign(s))}
    for i in range(4):
        entry = {frozenset(): B[i][s[i]]}
        if s[i] == i: entry[frozenset([i])] = G(1)
        term = pmul(term, entry)
    DET = padd(DET, term)
PERMS = list(permutations(range(4)))
Z = {}
for p in PERMS:
    z = {frozenset(): G(1)}
    for j in range(4): z = pmul(z, {frozenset([p[j]]): G(1), frozenset(): b[j]})
    Z[p] = z
def affine_a0(p):  # fix a1,a2,a3 at the fixture; return (constant, coefficient of a0)
    q = subst(p, {1: a[1], 2: a[2], 3: a[3]})
    return q.get(frozenset(), G(0)), q.get(frozenset([0]), G(0))

# ---- section 1-2: interval family a0(t) = (-2-i) + i t
def in_t(p):  # affine-in-t (constant, slope)
    c0, c1 = affine_a0(p); return (c0 + c1 * a[0], c1 * I)
d = in_t(DET)
assert d == (G(F(269802, 390625), F(-5707764, 390625)), G(F(-1235313, 390625), F(1850391, 390625))), 'formula (1)'
TRI = [(1, 3, 2, 0), (0, 2, 1, 3), (1, 2, 0, 3)]
assert [list(p) for p in TRI] == [w['columnToRow'] for w in src['determinant']['originalHullConvexWeights']]
zs = [in_t(Z[p]) for p in TRI]
assert zs == [(G(-48, -64), G(-4, 28)), (G(30, -10), G(-8, 6)), (G(12, 16), G(14, 2))], 'vertex table'
def across(u, v):  # [u, v] for affine-in-t complex u, v -> ascending real quadratic
    return [cross(u[0], v[0]), cross(u[0], v[1]) + cross(u[1], v[0]), cross(u[1], v[1])]
sub = lambda u, v: (u[0] - v[0], u[1] - v[1])
area = across(sub(zs[1], zs[0]), sub(zs[2], zs[0]))
assert area == [3000, -2000, 500]
N = [[15625 * c for c in across(sub(zs[(j + 1) % 3], d), sub(zs[(j + 2) % 3], d))] for j in range(3)]
K = [15625 * c for c in area]
assert N == [[13203996, -5027220, -131806], [14562000, -10814944, 6202444], [19109004, -15407836, 1741862]]
assert K == [46875000, -31250000, 7812500]
assert [sum(col) for col in zip(*N)] == K
def qtimes(n, z):  # real quadratic times complex affine -> complex cubic, ascending in t
    out = [G(0)] * 4
    for k in range(3):
        for m in range(2): out[k + m] = out[k + m] + F(n[k]) * z[m]
    return out
lhs = [G(0)] * 4
for j in range(3): lhs = [x + y for x, y in zip(lhs, qtimes(N[j], zs[j]))]
assert lhs == qtimes(K, d)  # sum_j N_j z_j = K d identically in t
ev = lambda c, t: c[0] + c[1] * t + c[2] * t * t
def qmin(c, l, r):  # exact minimum of a quadratic on [l, r]
    pts = [l, r] + ([F(-c[1], 2 * c[2])] if c[2] and l < F(-c[1], 2 * c[2]) < r else [])
    return min(ev(c, F(x)) for x in pts)
mins = [qmin(n, -40, 1) for n in N]
assert min(mins) > 0 and ev(K, 2) == 15625000 and K[2] > 0
assert N[1][1] ** 2 - 4 * N[1][0] * N[1][2] == -244316944388864  # N2 > 0 on all of R
for (j, l, r) in [(0, -41, -40), (0, 2, 3), (2, 1, F(3, 2)), (2, 7, 8)]:
    assert ev(N[j], F(l)) * ev(N[j], F(r)) < 0  # one root per bracket; brackets are disjoint and ordered
def root(c, sgn): disc = c[1] ** 2 - 4 * c[0] * c[2]; return (-c[1] + sgn * disc ** 0.5) / (2 * c[2])
alpha, gamma = root(N[0], 1), root(N[2], -1)
w = lambda t: [F(ev(n, F(t)), ev(K, F(t))) for n in N]
assert w(0) == [F(x['weight']) for x in src['determinant']['originalHullConvexWeights']]
assert w(0)[0] == F(1100333, 3906250) and w(1)[0] == F(804497, 2343750)

# ---- section 3: five fixed masses for every complex a0
PLANE = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3), (1, 0, 2, 3), (3, 1, 2, 0)]
MP = [F(22689, 390625), F(8848, 15625), F(1296, 15625), F(144, 15625), F(110736, 390625)]
assert sum(MP) == 1 and min(MP) > 0
dz = affine_a0(DET)
assert dz == (G(F(2735271, 390625), F(-1386747, 390625)), G(F(1850391, 390625), F(1235313, 390625)))
for k in range(2):
    assert sum((m * affine_a0(Z[p])[k] for m, p in zip(MP, PLANE)), G(0)) == dz[k]

# ---- stored certificate.json against this recomputation
cert = json.loads((SPEC / 'certificate.json').read_text())
assert [g(x) for x in cert['determinantAffineConstantSlope']] == list(d)
assert [[g(x) for x in v] for v in cert['verticesAffineConstantSlope']] == [list(z) for z in zs]
assert [F(x) for x in cert['orientedAreaCoefficients']] == area
assert [[F(x) for x in n] for n in cert['weightNumeratorsAscending']] == N
assert [F(x) for x in cert['weightCommonDenominatorAscending']] == K
assert [F(x) for x in cert['complexPlaneExtension']['weights']] == MP
for chk in cert['exactParameterChecks']:
    t = F(chk['t']); assert -40 <= t <= 1 and [F(x) for x in chk['weights']] == w(t)
    assert g(chk['determinant']) == d[0] + d[1] * t

# ---- section 4: twelve fixed masses on C^3 at a3 = 1-i; separator on C^4
fib = json.loads((SPEC / 'fiber_certificate.json').read_text())
M12 = {tuple(r['columnToRow']): F(r['mass']) for r in fib['weightsAtFixedFourthEigenvalue']}
assert len(M12) == 12 and sum(M12.values()) == 1 and min(M12.values()) > 0
mix = {}
for p, m in M12.items(): mix = padd(mix, subst(Z[p], {3: a[3]}), m)
assert padd(mix, subst(DET, {3: a[3]}), -1) == {}, 'all 8 coefficients in a0,a1,a2'
pt = lambda p, x: subst(p, {i: G(x[i]) for i in range(4)}).get(frozenset(), G(0))
def inside(P, Q, S, X):  # X in the closed triangle PQS (nondegenerate)
    A = cross(Q - P, S - P)
    return A != 0 and min(cross(Q - X, S - X) / A, cross(S - X, P - X) / A, cross(P - X, Q - X) / A) >= 0
def in_hull(pts, X):  # Caratheodory in the plane: a point, a segment or a triangle
    on_seg = lambda P, Q: cross(Q - P, X - P) == 0 and (X - P).re * (Q - X).re + (X - P).im * (Q - X).im >= 0
    return (any(P == X for P in pts) or any(on_seg(P, Q) for P, Q in combinations(pts, 2))
            or any(inside(*tri, X) for tri in combinations(pts, 3)))
scores, Ld = {p: F(0) for p in PERMS}, F(0)
for term in fib['separatingFunctional']:
    x, c, e_ = term['diagonal'], F(term['realCoefficient']), F(term['imagCoefficient'])
    assert len(x) == 4 and set(x) <= {0, 1}
    dv = pt(DET, x); Ld += c * dv.re + e_ * dv.im
    pts = [pt(Z[p], x) for p in PERMS]
    for p, v in zip(PERMS, pts): scores[p] += c * v.re + e_ * v.im
    # pointwise: d(x) lies in the hull of its own 24 values (masses may differ per diagonal)
    assert in_hull(pts, dv), f'diagonal {x} outside its hull'
assert max(scores.values()) == 0 and Ld == F(11201152, 15625)

print(json.dumps({
    'unitaryAndFixtureIdentity': True, 'determinantFormula1': True, 'vertexTable': True,
    'orientedArea': [str(x) for x in area], 'numeratorsMatch': True, 'polynomialIdentities': True,
    'exactMinimaOnInterval': [str(m) for m in mins],
    'nonnegativeRangeApprox': [round(alpha, 6), round(gamma, 6)],
    'fivePointPlaneCertificate': True, 'storedCertificateMatches': True,
    'twelveMassCoefficientsChecked': 8, 'separatorMaxPairingScore': '0',
    'separatorDeterminantScore': str(Ld), 'pairingsAttainingZero': sum(v == 0 for v in scores.values()),
    'separatorDiagonalsInsideHull': len(fib['separatingFunctional'])}, indent=1))
