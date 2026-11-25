#!/usr/bin/env python3
"""Exact constructive triangle repair and one-vertex gluing regressions."""
from fractions import Fraction as F
from itertools import combinations, permutations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import resource
import sys
import time

BASE = Path(__file__).resolve().parents[1] / 'baselines' / 'forest.py'
spec = importlib.util.spec_from_file_location('cycle_forest', BASE)
f = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = f
spec.loader.exec_module(f)
Q, Z, O, mul, adj, eye, diag, det = f.Q, f.Z, f.O, f.mul, f.adj, f.eye, f.diag, f.det


def stochastic(d):
    n = len(d)
    if not 1 <= n <= 6 or any(len(row) != n for row in d):
        raise ValueError('Expected a square matrix, dimension 1..6')
    if any(isinstance(v, bool) or not isinstance(v, (int, F)) or v < 0 for row in d for v in row):
        raise ValueError('Expected exact nonnegative rational entries')
    if any(sum(row) != 1 for row in d) or any(sum(d[i][j] for i in range(n)) != 1 for j in range(n)):
        raise ValueError('Matrix is not doubly stochastic')


def triangle_weights(d, x=None):
    stochastic(d)
    if len(d) != 3:
        raise ValueError('Triangle formula requires dimension three')
    lo = max(F(0), d[2][2]-d[0][1], d[1][1]-d[0][2])
    hi = min(d[0][0], d[1][1], d[2][2])
    if x is None:
        x = lo
    if isinstance(x, bool) or not isinstance(x, (int, F)) or not lo <= x <= hi:
        raise ValueError('Parameter must be exact and in the feasible interval')
    values = (x, d[0][0]-x, d[2][2]-x, d[0][2]-d[1][1]+x,
              d[0][1]-d[2][2]+x, d[1][1]-x)
    weights = {p: F(w) for p, w in zip(permutations(range(3)), values) if w}
    assert sum(weights.values()) == 1 and all(w >= 0 for w in weights.values())
    return weights, (lo, hi)


def greedy_birkhoff(d):
    """Lexicographic exact decomposition; used only to expose the n=4 limit."""
    stochastic(d)
    r = [list(map(F, row)) for row in d]
    n, weights = len(r), {}
    while any(v for row in r for v in row):
        p = next(p for p in permutations(range(n)) if all(r[p[j]][j] > 0 for j in range(n)))
        amount = min(r[p[j]][j] for j in range(n))
        weights[p] = weights.get(p, F(0)) + amount
        for j in range(n):
            r[p[j]][j] -= amount
    return weights


def certificate(u):
    n = len(u)
    if mul(u, adj(u)) != eye(n) or mul(adj(u), u) != eye(n):
        raise ValueError('Matrix is not exactly unitary')
    d = [[f.norm2(v) for v in row] for row in u]
    if n == 3:
        return triangle_weights(d)[0]
    if n <= 2:
        return greedy_birkhoff(d)
    raise ValueError('No general certificate constructor above dimension three')


def embed(n, support, block, weights):
    if len(set(support)) != len(support) or len(block) != len(support):
        raise ValueError('Invalid support')
    if any(type(i) is not int or not 0 <= i < n for i in support):
        raise ValueError('Invalid support coordinates')
    u = eye(n)
    for i, a in enumerate(support):
        for j, b in enumerate(support):
            u[a][b] = block[i][j]
    expanded = {}
    for p, w in weights.items():
        q = list(range(n))
        for j, a in enumerate(support):
            q[a] = support[p[j]]
        expanded[tuple(q)] = w
    return u, expanded, frozenset(support)


def glue(left, right):
    x, wx, a = left
    y, wy, b = right
    if len(a & b) > 1:
        raise ValueError('One-vertex gluing requires support overlap at most one')
    if len(x) != len(y):
        raise ValueError('Dimensions differ')
    for matrix, weights, support in (left, right):
        n = len(matrix)
        if any(matrix[i][j] != (O if i == j else Z)
               for i in range(n) for j in range(n) if i not in support or j not in support):
            raise ValueError('Matrix acts outside its declared coordinate support')
        if any(any(p[j] != j for j in range(n) if j not in support) for p in weights):
            raise ValueError('A permutation acts outside its declared support')
    weights = {}
    for p, wp in wx.items():
        for q, wq in wy.items():
            r = tuple(p[q[j]] for j in range(len(x)))
            weights[r] = weights.get(r, F(0)) + wp*wq
    return mul(x, y), weights, a | b


def compound_check(u, weights):
    n, count, first = len(u), 0, None
    assert sum(weights.values()) == 1 and all(w >= 0 for w in weights.values())
    for k in range(n+1):
        subsets = list(combinations(range(n), k))
        for s in subsets:
            for t in subsets:
                actual = f.norm2(det([[u[i][j] for j in t] for i in s]))
                predicted = sum((w for p, w in weights.items() if tuple(sorted(p[j] for j in t)) == s), F(0))
                if actual != predicted and first is None:
                    first = {'rows': list(s), 'columns': list(t), 'actual': str(actual), 'predicted': str(predicted)}
                count += 1
    return count, first


def determinant_check(u, weights, a=None, b=None):
    n = len(u)
    a = a if a is not None else [Q(i-2, (i % 3)-1) for i in range(n)]
    b = b if b is not None else [Q(3-i, (2*i % 5)-2) for i in range(n)]
    orbit = mul(mul(u, diag(b)), adj(u))
    matrix = [[orbit[i][j]+(a[i] if i == j else Z) for j in range(n)] for i in range(n)]
    actual = det(matrix)
    assert actual == f.determinant_leibniz(matrix)
    predicted = sum((w*f.vertex(a, b, p) for p, w in weights.items()), Z)
    point_permutation = {f.vertex(a, b, p).pair(): p for p in permutations(range(n))}
    points = list(point_permutation)
    hull = f.exact.certify_hull(points, actual.pair())
    assert f.exact.check_hull_certificate(points, actual.pair(), hull)
    return {'a': [z.encode() for z in a], 'b': [z.encode() for z in b],
            'actual': actual.encode(), 'weighted': predicted.encode(), 'equal': actual == predicted,
            'originalHullMembership': hull['membership'],
            'originalHullConvexWeights': [{'columnToRow': list(point_permutation[point]),
                                          'point': list(map(str, point)), 'weight': str(weight)}
                                         for point, weight in hull.get('weights', [])]}


def block3(phased=False, order=(0,1,2)):
    base = [[Q(F(3,5)), Q(F(4,5))], [Q(F(-4,5)), Q(F(3,5))]]
    edges = [(0,1), (1,2), (0,2)]
    gates = [(edges[k][0], edges[k][1], f.phased_block(variant=k) if phased else base) for k in order]
    return f.build_unitary(3, gates, require_forest=False), gates


def fixtures():
    for phased in (False, True):
        for order in permutations(range(3)):
            u, _ = block3(phased, order)
            yield 'triangle-'+('complex' if phased else 'real')+'-'+''.join(map(str,order)), u, certificate(u)
    u, _ = block3(True)
    core = embed(6, [0,1,2], u, certificate(u))
    for support, left in (([2,3], True), ([0,4], False), ([3,5], True)):
        block = f.phased_block(variant=sum(support))
        gate = embed(6, support, block, certificate(block))
        core = glue(gate, core) if left else glue(core, gate)
    yield 'complex-triangle-three-pendant-gates-both-ends', core[0], core[1]
    u2, _ = block3(True, (2,0,1))
    two = glue(embed(5, [0,1,2], u, certificate(u)), embed(5, [2,3,4], u2, certificate(u2)))
    yield 'two-complex-triangle-blocks-one-pivot', two[0], two[1]


def report():
    cases = []
    for name, u, w in fixtures():
        count, failure = compound_check(u, w)
        assert failure is None
        d = determinant_check(u, w)
        assert d['equal'] and d['originalHullMembership'] == 'inside'
        cases.append({'name': name, 'dimension': len(u), 'U': [[v.encode() for v in row] for row in u],
                      'weights': [{'columnToRow': list(p), 'weight': str(v)} for p, v in sorted(w.items())],
                      'compoundEntriesChecked': count, 'determinant': d})
    triangle, gates = block3()
    independent = f.permutation_weights(3, gates, require_forest=False)
    corrected, interval = triangle_weights([[f.norm2(z) for z in row] for row in triangle])
    old = sum((w for p,w in independent.items() if p[0] == 0), F(0))
    new = sum((w for p,w in corrected.items() if p[0] == 0), F(0))
    assert old == F(6121,15625) and new == F(11881,15625)
    c = [[Q(F(3,5)), Q(F(4,5))], [Q(F(-4,5)), Q(F(3,5))]]
    four = f.build_unitary(4, [(0,1,c),(1,2,c),(2,3,c),(0,3,c)], require_forest=False)
    four_w = greedy_birkhoff([[f.norm2(z) for z in row] for row in four])
    _, mismatch = compound_check(four, four_w)
    assert mismatch is not None and len(mismatch['rows']) == 2
    result = {'schemaVersion': 1, 'evidenceLevel': 'exact same-owner regressions plus separate written proof; no independent review',
              'pythonVersion': platform.python_version(), 'arithmetic': 'Gaussian rationals, no tolerance',
              'bounds': {'fixedFixtures':len(cases), 'maximumDimension':6, 'cpuLimitSeconds':60},
              'compoundEntriesChecked':sum(c['compoundEntriesChecked'] for c in cases), 'cases':cases,
              'triangleRepair': {'oldProductReturnProbability':str(old), 'correctedReturnProbability':str(new),
                                 'feasibleXInterval':list(map(str,interval)), 'chosenX':str(interval[0])},
              'fourCycleGreedyRuleFailure': {'scope':'This particular singleton Birkhoff decomposition fails; no impossibility of all weights or MOC claimed',
                                           'U':[[v.encode() for v in row] for row in four], 'firstCompoundMismatch':mismatch,
                                           'greedyWeights': [{'columnToRow': list(p), 'weight': str(w)} for p,w in sorted(four_w.items())],
                                           'determinant':determinant_check(four, four_w)},
              'dependencySha256': {str(p.relative_to(BASE.parents[1])):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE, f.VERIFY]}}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (60,60))
    started = time.process_time()
    result = report()
    result['processCpuSeconds'] = round(time.process_time()-started,6)
    encoded = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
        print(json.dumps({k:result[k] for k in ('bounds','compoundEntriesChecked','triangleRepair','processCpuSeconds')}))
    else:
        print(encoded,end='')

if __name__ == '__main__':
    main()
