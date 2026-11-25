#!/usr/bin/env python3
"""Bounded exact triangle LPs, with rational primal/dual certificates.

No external solver or repository imports. The verifier below checks each witness
independently of the simplex tableau that found it. No MOC counterexample claim.
"""
from fractions import Fraction as F
from itertools import combinations
import json


def packing(p):
    """Maximize total triangle mass under edge capacities p_i*p_j."""
    edges = list(combinations(range(len(p)), 2))
    triangles = list(combinations(range(len(p)), 3))
    m, v = len(edges), len(triangles)
    rows = [
        [F(set(edge) <= set(triangle)) for triangle in triangles]
        + [F(i == j) for j in range(m)]
        + [p[edge[0]] * p[edge[1]]]
        for i, edge in enumerate(edges)
    ]
    objective = [F(-1)] * v + [F(0)] * (m + 1)
    basis = list(range(v, v + m))
    steps = 0
    while any(value < 0 for value in objective[:-1]):
        # Bland's rule gives deterministic anti-cycling choices.
        column = next(j for j, value in enumerate(objective[:-1]) if value < 0)
        eligible = [i for i in range(m) if rows[i][column] > 0]
        assert eligible, 'Unexpected unbounded triangle packing'
        pivot = min(eligible, key=lambda i: (rows[i][-1] / rows[i][column], basis[i]))
        divisor = rows[pivot][column]
        rows[pivot] = [value / divisor for value in rows[pivot]]
        for i in range(m):
            if i != pivot:
                factor = rows[i][column]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[pivot])]
        factor = objective[column]
        objective = [x - factor * y for x, y in zip(objective, rows[pivot])]
        basis[pivot] = column
        steps += 1
        assert steps < 10000, 'Bounded regression pivot limit exceeded'
    primal = [F(0)] * v
    for i, j in enumerate(basis):
        if j < v:
            primal[j] = rows[i][-1]
    return primal, objective[v:v + m], steps


def verify(p, primal, dual):
    """Feasibility plus weak-duality equality certifies optimality exactly."""
    assert p and all(value > 0 for value in p) and sum(p) == 1
    edges = list(combinations(range(len(p)), 2))
    triangles = list(combinations(range(len(p)), 3))
    assert len(primal) == len(triangles) and len(dual) == len(edges)
    assert all(value >= 0 for value in primal + dual)
    for (i, j) in edges:
        used = sum(value for triangle, value in zip(triangles, primal) if i in triangle and j in triangle)
        assert used <= p[i] * p[j]
    for i, j, k in triangles:
        covered = sum(value for edge, value in zip(edges, dual) if set(edge) <= {i, j, k})
        assert covered >= 1
    lower = sum(primal)
    upper = sum(p[i] * p[j] * value for (i, j), value in zip(edges, dual))
    assert lower == upper
    return lower


def case_report(name, p):
    primal, dual, steps = packing(p)
    optimum = verify(p, primal, dual)
    edges = list(combinations(range(len(p)), 2))
    triangles = list(combinations(range(len(p)), 3))
    e2 = sum(p[i] * p[j] for i, j in edges)
    cut = max(
        mass * (1 - mass)
        for mask in range(1 << (len(p) - 1))
        for mass in [sum(p[i] for i in range(len(p)) if mask >> i & 1)]
    )
    beta = e2 - optimum
    return {
        'name': name,
        'p': list(map(str, p)),
        'triangle_packing_optimum': str(optimum),
        'e2': str(e2),
        'fractional_triangle_free_optimum': str(beta),
        'formal_c_threshold': str(1 / beta),
        'physical_c_threshold': str(min(F(4), 1 / beta)),
        'constant_edge_lower_bound': str(2 * e2 / 3),
        'cut_lower_bound': str(cut),
        'shortcut_is_exact': beta == max(2 * e2 / 3, cut),
        'half_turn_has_common_minor_certificate': 4 * beta <= 1,
        'primal_triangles': [{'vertices_zero_based': list(t), 'weight': str(w)} for t, w in zip(triangles, primal) if w],
        'dual_edges': [{'vertices_zero_based': list(e), 'cover_weight': str(w)} for e, w in zip(edges, dual)],
        'simplex_pivots': steps,
    }


def main():
    reports = [case_report(f'flat-{n}', [F(1, n)] * n) for n in range(3, 8)]
    for report in reports:
        n = len(report['p'])
        assert F(report['formal_c_threshold']) == F(3 * n, n - 1)
    nonflat = case_report('trace-and-cut-shortcut-fails', [F(i, 43) for i in [4, 2, 6, 1, 19, 8, 3]])
    assert F(nonflat['triangle_packing_optimum']) == F(649, 5547)
    assert F(nonflat['fractional_triangle_free_optimum']) == F(1388, 5547)
    assert not nonflat['shortcut_is_exact']
    assert not nonflat['half_turn_has_common_minor_certificate']
    assert 4 * F(nonflat['cut_lower_bound']) < 1
    assert 4 * F(nonflat['constant_edge_lower_bound']) < 1
    reports.append(nonflat)
    print(json.dumps({
        'status': 'passed',
        'arithmetic': 'exact rational, Python standard library',
        'scope': 'common-minor certificates only; no original-MOC counterexample',
        'cases': reports,
    }, indent=2))


if __name__ == '__main__':
    main()
