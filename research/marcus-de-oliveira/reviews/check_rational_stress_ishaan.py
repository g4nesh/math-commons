#!/usr/bin/env python3
"""Reviewer's exact re-check of rational-stress/ (task mdo-rational-stress).

Standard library only. The one component shared with the artifact is the
input generator screen.make_witness (needed because report.json stores seeds,
not matrices); every saved witness file is also compared with it. Everything
else is separate: Gaussian rationals as Fraction pairs, a Leibniz determinant
over Gaussian integers after clearing denominators, gift-wrapping hulls over
integer-scaled points, full-matrix commutators and normality, own psi and
own (A,B) distinctness hash. Run from the repository root:
    python3 -I -B research/marcus-de-oliveira/reviews/check_rational_stress_ishaan.py
"""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations
from math import lcm
from pathlib import Path
import hashlib, json, sys

DIR = Path(__file__).resolve().parent.parent / 'rational-stress'
sys.path.insert(0, str(DIR))
import screen  # generator only

def cm(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def ca(x, y): return (x[0]+y[0], x[1]+y[1])
def cj(x): return (x[0], -x[1])
def pc(v): return (F(v[0]), F(v[1]))
def mat(X, Y):
    n = len(X); out = []
    for i in range(n):
        row = []
        for j in range(n):
            s = (F(0), F(0))
            for k in range(n): s = ca(s, cm(X[i][k], Y[k][j]))
            row.append(s)
        out.append(row)
    return out
def adj(X): return [[cj(X[j][i]) for j in range(len(X))] for i in range(len(X))]
def sign(p):
    s, seen = 1, [False]*len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            if L % 2 == 0: s = -s
    return s
def leibniz(M):
    n = len(M); D = lcm(*[x.denominator for row in M for z in row for x in z])
    Z = [[(int(z[0]*D), int(z[1]*D)) for z in row] for row in M]
    tot = (0, 0)
    for p in permutations(range(n)):
        t = (sign(p), 0)
        for i in range(n): t = cm(t, Z[i][p[i]])
        tot = ca(tot, t)
    return (F(tot[0], D**n), F(tot[1], D**n))
def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
def gift_wrap(pts):
    """Extreme vertices, counterclockwise; integer points; collinear points skipped."""
    pts = sorted(set(pts))
    if len(pts) == 1: return pts
    if all(cross(pts[0], pts[-1], p) == 0 for p in pts): return [pts[0], pts[-1]]
    hull, cur = [], pts[0]
    while True:
        hull.append(cur); cand = pts[0] if pts[0] != cur else pts[1]
        for p in pts:
            c = cross(cur, cand, p)
            d2 = lambda q: (q[0]-cur[0])**2+(q[1]-cur[1])**2
            if c < 0 or (c == 0 and d2(p) > d2(cand)): cand = p
        cur = cand
        if cur == hull[0]: return hull
def psi(U):
    sq = lambda z: z[0]*z[0]+z[1]*z[1]
    m = lambda r, c: sq(ca(cm(U[r[0]][c[0]], U[r[1]][c[1]]), cm((-1, 0), cm(U[r[0]][c[1]], U[r[1]][c[0]]))))
    return (sum(sq(U[i][i]) for i in range(4)) - m((0,1),(0,2)) - m((0,1),(1,3)) + m((0,2),(0,1)) + m((0,2),(2,3))) / 2

def check(w, rec, out):
    a, b = [pc(v) for v in w['a']], [pc(v) for v in w['b']]
    U = [[pc(v) for v in row] for row in w['U']]; n = len(a)
    I = [[(F(int(i == j)), F(0)) for j in range(n)] for i in range(n)]
    assert mat(U, adj(U)) == I and mat(adj(U), U) == I, 'unitarity'
    B = mat(mat(U, [[b[i] if i == j else (F(0), F(0)) for j in range(n)] for i in range(n)]), adj(U))
    assert mat(B, adj(B)) == mat(adj(B), B), 'normality of B'
    A = [[a[i] if i == j else (F(0), F(0)) for j in range(n)] for i in range(n)]
    C = [[ca(x, cm((-1, 0), y)) for x, y in zip(r1, r2)] for r1, r2 in zip(mat(A, B), mat(B, A))]
    norm = sum(z[0]**2+z[1]**2 for row in C for z in row)
    assert norm > 0 and norm == F(rec['commutatorNormSquared']), 'commutator'
    assert any(z[1] for row in U for z in row) and rec['complexUnitary'], 'complex U'
    d = leibniz([[ca(A[i][j], B[i][j]) for j in range(n)] for i in range(n)])
    assert d == pc(rec['determinant']), 'determinant'
    pts = {}
    for p in permutations(range(n)):
        z = (F(1), F(0))
        for i in range(n): z = cm(z, ca(a[i], b[p[i]]))
        pts[tuple(j+1 for j in p)] = z
    assert len(pts) == rec['permutationCount']; out['points'] += len(pts)
    cc = rec['convexCombination']; ws = [F(x['weight']) for x in cc]
    assert 1 <= len(cc) <= 3 and min(ws) >= 0 and sum(ws) == 1, 'weights'
    rc = (F(0), F(0))
    for x, wt in zip(cc, ws):
        assert pts[tuple(x['permutation'])] == pc(x['point']), 'support point'
        rc = ca(rc, cm(pts[tuple(x['permutation'])], (wt, F(0))))
    assert rc == d, 'convex combination'
    D = lcm(*[c.denominator for z in pts.values() for c in z])
    ip = [(int(z[0]*D), int(z[1]*D)) for z in pts.values()]
    hull = gift_wrap(ip); dd = (d[0]*D, d[1]*D)
    assert sorted((F(x, D), F(y, D)) for x, y in hull) == sorted(pc(h) for h in rec['hull']), 'hull'
    if len(hull) == 2:
        assert cross(hull[0], hull[1], dd) == 0 and all(min(h[k] for h in hull) <= dd[k] <= max(h[k] for h in hull) for k in (0, 1))
        cls, m = 'degenerate-inside', None
    else:
        L = max(max(p[k] for p in hull)-min(p[k] for p in hull) for k in (0, 1))
        m = min(cross(p, q, dd) / ((abs(q[0]-p[0])+abs(q[1]-p[1]))*L) for p, q in zip(hull, hull[1:]+hull[:1]))
        assert m >= 0 and m == F(rec['relativeLinfSupportMargin']), 'margin'
        cls = 'boundary' if m == 0 else 'interior'
    assert cls == rec['classification'], 'classification'
    if n == 4:
        v = psi(U); assert v == F(rec['exoticPsi']), 'psi'; out['psi<0'] += v < 0
        out['psi<0 general'] += v < 0 and rec['spec']['family'] in ('generic', 'thin-complex', 'clustered')
    key = json.dumps([[str(x) for z in a for x in z], [[str(x) for x in z] for row in B for z in row]])
    out['pairs'].add(hashlib.sha256(key.encode()).hexdigest())
    s = rec['spec']; out['cls'][(s['n'], cls)] += 1; out['ctrl' if s['family'] in ('singular', 'two-value', 'hermitian') else 'general'] += 1
    if s['family'] in ('generic', 'thin-complex'):
        out['re'][s['family']].add(tuple(z[0] for z in a) + tuple(sorted(z[0] for z in b)))
    if s['family'] == 'thin-complex' and s['mixing'] == 'near-vertex' and m:
        t = F(1, 2 ** (4, 8, 10)[s['repetition'] % 3]); eps = F(1, 2 ** (8 + 2*(s['repetition'] % 6)))
        out['ratio'].append(float(m / (t*t*eps)))

def main():
    O = [[(F(x, 7), F(0)) for x in r] for r in ((1,4,4,4), (-4,1,4,-4), (-4,-4,1,4), (-4,4,-4,1))]
    assert psi(O) == F(-158, 2401)
    assert all(psi([[(F(int(p[i] == j)), F(0)) for j in range(4)] for i in range(4)]) >= 0 for p in permutations(range(4)))
    out = dict(points=0, pairs=set(), cls=Counter(), general=0, ctrl=0, re={'generic': set(), 'thin-complex': set()}, ratio=[])
    out['psi<0'] = out['psi<0 general'] = 0
    saved = 0
    for rel, wdir in (('report.json', 'witnesses'), ('holdout/report.json', 'holdout/witnesses')):
        rep = json.loads((DIR / rel).read_text())
        for rec in rep['records']:
            w = screen.make_witness(rec['spec'])
            canon = json.dumps(w, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()
            assert hashlib.sha256(canon).hexdigest() == rec['witnessSha256'], 'input hash'
            if rec['spec']['id'] in rep['selectedWitnesses']:
                f = DIR / wdir / (rec['spec']['id'] + '.json')
                assert json.loads(f.read_text()) == w, 'saved witness'; saved += 1
            check(w, rec, out)
    print(json.dumps({'cases': sum(out['cls'].values()), 'byDimensionAndClass': {'%d %s' % k: v for k, v in sorted(out['cls'].items())},
                      'generalComplex': out['general'], 'controls': out['ctrl'], 'permutationPoints': out['points'],
                      'distinctExactPairs': len(out['pairs']), 'savedWitnessFilesMatched': saved,
                      'n4NegativePsi': out['psi<0'], 'n4NegativePsiGeneral': out['psi<0 general'],
                      'distinctRealPartConfigs': {k: len(v) for k, v in out['re'].items()},
                      'thinNearVertexMarginOverT2Eps': [round(min(out['ratio']), 3), round(max(out['ratio']), 3)],
                      'result': 'all exact checks passed'}, indent=1))

if __name__ == '__main__':
    main()
