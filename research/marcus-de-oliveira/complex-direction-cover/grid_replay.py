"""Recompute the full theorem on an exact unisolvent grid, independent of check.py.

The weighted gap has degrees at most (d+1, R+1, 1, 1). Evaluating on
(d+2)(R+2) times four distinct grid points proves polynomial identity.
"""
from fractions import Fraction as F
from itertools import product, permutations
from math import comb
from pathlib import Path
import hashlib
import json
SOURCE = Path(__file__).with_name('patches.json')
ZERO = (F(0), F(0))
ONE = (F(1), F(0))
U = [[F(x, 625) for x in row] for row in ((-31, 492, -240, 300), (-300, 225, 500, 0), (240, -180, 225, 500), (-492, -256, -180, 225))]
B_SPEC = [(F(3), F(-2)), (F(2), F(0)), (F(1), F(2)), (F(0), F(-1))]
PERMS = tuple(permutations(range(4)))

def add(z, w):
    return (z[0] + w[0], z[1] + w[1])

def neg(z):
    return (-z[0], -z[1])

def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

def scale(x, z):
    return (x * z[0], x * z[1])

def div(z, w):
    d = w[0] ** 2 + w[1] ** 2
    assert d
    return ((z[0] * w[0] + z[1] * w[1]) / d, (z[1] * w[0] - z[0] * w[1]) / d)

def prod(values):
    ans = ONE
    for v in values:
        ans = mul(ans, v)
    return ans

def det_dense(matrix):
    a = [list(r) for r in matrix]
    ans = ONE
    for j in range(4):
        pivot = next((i for i in range(j, 4) if a[i][j] != ZERO), None)
        if pivot is None:
            return ZERO
        if pivot != j:
            (a[pivot], a[j]) = (a[j], a[pivot])
            ans = neg(ans)
        ans = mul(ans, a[j][j])
        for i in range(j + 1, 4):
            ratio = div(a[i][j], a[j][j])
            for k in range(j + 1, 4):
                a[i][k] = add(a[i][k], neg(mul(ratio, a[j][k])))
            a[i][j] = ZERO
    return ans
for (i, j) in product(range(4), repeat=2):
    assert sum((U[i][k] * U[j][k] for k in range(4))) == int(i == j)
B = [[ZERO for _ in range(4)] for _ in range(4)]
for (i, j, k) in product(range(4), repeat=3):
    B[i][j] = add(B[i][j], scale(U[i][k] * U[j][k], B_SPEC[k]))

def direct_values(x, middle, y):
    a = [(F(x), F(0)), middle, (F(y), F(0)), ZERO]
    matrix = [[add(B[i][j], a[i] if i == j else ZERO) for j in range(4)] for i in range(4)]
    vertices = [prod((add(a[p[j]], B_SPEC[j]) for j in range(4))) for p in PERMS]
    return (det_dense(matrix), vertices)
patches = json.loads(SOURCE.read_text())
reports = []
for (index, patch) in enumerate(patches):
    chart = patch['chart']
    (lo, hi) = map(F, patch['interval'])
    d = patch['degree']
    assert chart in ('h', 'v') and -1 <= lo < hi <= 1 and isinstance(d, int) and (d >= 0)
    decoded = []
    roots = {j: None for j in range(d + 1)}
    masses = {j: F(0) for j in roots}
    seen = set()
    for (j, word, origins, mass) in patch['terms']:
        p = tuple((int(c) - 1 for c in word))
        centers = tuple(map(F, origins))
        mass = F(mass)
        assert p in PERMS and j in roots and (mass > 0) and ((j, p, centers) not in seen)
        seen.add((j, p, centers))
        decoded.append((j, PERMS.index(p), centers, mass))
        roots[j] = set(centers) if roots[j] is None else roots[j].intersection(centers)
        masses[j] += mass
    assert all((r == set() for r in roots.values()))
    assert all((m == 1 for m in masses.values()))
    R = 2 * max((len(c) for (_, _, c, _) in decoded))
    slopes = [lo + (hi - lo) * F(j, d + 1) for j in range(d + 2)]
    radii = list(range(-(R + 2) // 2, -(R + 2) // 2 + R + 2))
    count = 0
    for (slope, radius) in product(slopes, radii):
        alpha = (slope - lo) / (hi - lo)
        numerators = [F(0)] * 24
        for (j, p, centers, mass) in decoded:
            value = mass * comb(d, j) * alpha ** j * (1 - alpha) ** (d - j)
            for c in centers:
                value *= (F(radius) - c) ** 2
            numerators[p] += value
        denominator = sum(numerators)
        assert denominator > 0 and min(numerators) >= 0
        middle = (F(radius), slope * radius) if chart == 'h' else (slope * radius, F(radius))
        for (x, y) in product((0, 1), repeat=2):
            (target, vertices) = direct_values(x, middle, y)
            mixture = ZERO
            for (mass, vertex) in zip(numerators, vertices):
                mixture = add(mixture, scale(mass, vertex))
            assert mixture == scale(denominator, target), (index, slope, radius, x, y)
            count += 1
    reports.append({'patch': index, 'chart': chart, 'interval': [str(lo), str(hi)], 'angularDegreeOfGap': d + 1, 'radialDegreeOfGap': R + 1, 'slopeNodes': len(slopes), 'radialNodes': len(radii), 'freeEigenvalueCorners': 4, 'denseDeterminants': count, 'vertexEvaluations': 24 * count, 'positiveTerms': len(decoded), 'globalControlRootIntersectionsEmpty': True})
    print('patch', index, 'passed', count, 'original determinant evaluations', flush=True)
for chart in ('h', 'v'):
    intervals = sorted(((F(p['interval'][0]), F(p['interval'][1])) for p in patches if p['chart'] == chart))
    assert intervals[0][0] == -1 and intervals[-1][1] == 1
    assert all((b == c for ((_, b), (c, _)) in zip(intervals, intervals[1:])))
out = {'patchesSha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(), 'method': 'Independent rational Gaussian elimination on the original dense matrices; tensor interpolation grid', 'patchCount': len(reports), 'positiveSquareProductTerms': sum((p['positiveTerms'] for p in reports)), 'originalDenseDeterminants': sum((p['denseDeterminants'] for p in reports)), 'allVertexEvaluations': sum((p['vertexEvaluations'] for p in reports)), 'patches': reports}
assert out['positiveSquareProductTerms'] == 878 and out['originalDenseDeterminants'] == 576 and (out['allVertexEvaluations'] == 13824)
print(json.dumps({k: v for (k, v) in out.items() if k != 'patches'}, indent=2))
