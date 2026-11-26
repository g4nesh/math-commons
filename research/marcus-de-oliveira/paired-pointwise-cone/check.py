"""Exact, solver-free obstruction to the 45-generator spectral support claim.

The seven orbit representatives are attributed to PR81. This checks the
specific obstruction and repairs; it does not reprove PR81's completeness.
"""
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations_with_replacement(range(6), 2))
EDGES = tuple(product(range(3), range(3, 6)))
EDGE_PAIRS = tuple(combinations_with_replacement(EDGES, 2))
SIGNS = (1, 1, 1, -1, -1, -1)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parity(p):
    return (-1) ** sum(x > y for i, x in enumerate(p) for y in p[i + 1:])


PERMS = tuple(sorted(permutations(range(3)), key=lambda p: -parity(p)))


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def rank(rows):
    """Fraction-free row elimination, with integer gcd reduction."""
    pivots = {}
    for source in rows:
        row = list(source)
        for j, pivot in pivots.items():
            if row[j]:
                a, b = pivot[j], row[j]
                row = [a * x - b * y for x, y in zip(row, pivot)]
                divisor = gcd(*row)
                if divisor > 1:
                    row = [x // divisor for x in row]
        j = next((i for i, x in enumerate(row) if x), None)
        if j is not None:
            pivots[j] = row
    return len(pivots)


def generators():
    return [tuple(int(i in e and j in f) + int(i in f and j in e)
                  for i, j in PAIRS) for e, f in EDGE_PAIRS]


def monomial_rows():
    """All coefficients of the 21 products of two pairing polynomials."""
    result = {}
    for column, (i, j) in enumerate(PAIRS):
        assignment = PERMS[i] + PERMS[j]
        for take_b in product((False, True), repeat=6):
            a = tuple(sum(not take_b[k] and k % 3 == r for k in range(6))
                      for r in range(3))
            b = tuple(sum(take_b[k] and assignment[k] == r for k in range(6))
                      for r in range(3))
            result.setdefault(a + b, [0] * 21)[column] += 1
    return result


def primitive(row):
    divisor = gcd(*row)
    return tuple(x // divisor for x in row)


def orbit_rows(representatives):
    actions = [p for p in permutations(range(6)) if all(
        SIGNS[i] * SIGNS[j] == SIGNS[p[i]] * SIGNS[p[j]] for i, j in PAIRS)]
    require(len(actions) == 72, "parity-preserving symmetry group")
    result = []
    for sparse in representatives:
        row = [0] * 21
        for i, j, value in sparse:
            row[PAIRS.index((i, j))] = value
        result.append({tuple(row[PAIRS.index(tuple(sorted((p[i], p[j]))))]
                             for i, j in PAIRS) for p in actions})
    return result


def verify(data):
    require(data['schema'] == 'mdo.paired-pointwise-cone.obstruction.v1', 'schema')
    require(tuple(map(tuple, data['permutations'])) == PERMS, 'permutation order')
    q = [F(x, data['dual_denominator']) for x in data['dual_numerators']]
    require(len(q) == 21 and sum(q) == 1, 'normalized dual')
    cols = generators()
    values = [dot(q, col) for col in cols]
    require(min(values) >= 0, 'dual must satisfy all 45 inequalities')
    active = [col for col, value in zip(cols, values) if not value]
    require(rank(cols) == 15 and rank(active) == 14, 'normalized dual vertex ranks')
    require(len(active) == 20, 'active generators')
    require(rank(active + [[1] * 21]) == 15, 'independent normalization')

    rows = monomial_rows()
    require(len(rows) == 141 and rank(rows.values()) == 15, 'grouped coefficient space')
    require(min(dot(q, row) for row in rows.values()) >= 0, 'all grouped coefficients')
    kernel = [[int(i == j == v) if i == j else SIGNS[i] * SIGNS[j] * int(v in (i, j))
               for i, j in PAIRS] for v in range(6)]
    require(rank(kernel) == 6, 'six independent polynomial relations')
    require(all(dot(k, row) == 0 for k in kernel for row in rows.values()), 'polynomial relations')

    p = []
    for perm in PERMS:
        value = (1, 0)
        for i, j in enumerate(perm):
            value = mul(value, add(data['a'][i], data['b'][j]))
        p.append(value)
    products = [mul(p[i], p[j]) for i, j in PAIRS]
    scores = [mul(data['theta'], value)[0] for value in products]
    minimum = min(scores)
    gap = [score - minimum for score in scores]
    require(minimum == -850560 and dot(q, gap) == -20832, 'strict exact spectral obstruction')
    require(all(dot(k, gap) == 0 for k in kernel), 'gap lies in the 15-dimensional space')
    target = tuple(sum(q[k] * value[part] for k, value in enumerate(products)) for part in (0, 1))
    require(target == (120416, 375488), 'abstract signed target')

    orbits = orbit_rows(data['pr81_orbit_representatives'])
    require(list(map(len, orbits)) == [9, 18, 18, 36, 36, 36, 36], 'PR81 orbit sizes')
    all_rays = set().union(*orbits)
    require(len(all_rays) == 189, 'disjoint PR81 orbits')
    require({primitive(c) for c in cols} == set().union(*orbits[:3]), '45 generator rays')
    require(all(dot(k, ray) == 0 for k in kernel for ray in all_rays), 'ray relations')
    negative = {ray for ray in all_rays if dot(q, ray) < 0}
    counts = [sum(dot(q, ray) < 0 for ray in orbit) for orbit in orbits]
    require(counts == [0, 0, 0, 2, 1, 2, 2], 'seven violated rays in four missing orbits')
    repaired = set()
    repair_report = []
    for record in data['repairs']:
        extra = tuple(record['extra_ray'])
        require(extra in negative and extra not in repaired, 'repair ray')
        require(extra in orbits[record['orbit'] - 1], 'repair orbit')
        require(dot(q, extra) == F(record['dual_value']), 'ray dual value')
        total = [F(0)] * 21
        seen = set()
        for index, coefficient in record['positive_terms']:
            weight = F(coefficient)
            require(0 <= index <= 45 and index not in seen and weight > 0, 'positive repair term')
            seen.add(index)
            column = extra if index == 45 else cols[index]
            total = [x + weight * y for x, y in zip(total, column)]
        require(45 in seen and len(seen) == 14 and total == gap, 'exact repaired gap identity')
        repaired.add(extra)
        repair_report.append({'orbit': record['orbit'], 'dual_value': str(dot(q, extra)), 'positive_terms': len(seen)})
    require(repaired == negative, 'all seven single-ray repairs')
    return {'permutations': [list(x) for x in PERMS], 'pairing_values': [list(x) for x in p],
            'minimum_vertex_score': minimum, 'minimum_pairs': [list(PAIRS[k]) for k, x in enumerate(scores) if x == minimum],
            'abstract_target': [int(x) for x in target], 'abstract_target_score': int(dot(q, scores)),
            'separation': int(dot(q, gap)), 'generator_rank': 15, 'active_count': 20, 'active_rank': 14,
            'generator_dual_value_counts': dict(sorted(Counter(map(str, values)).items())),
            'grouped_coefficient_count': len(rows), 'grouped_coefficient_minimum': str(min(dot(q, row) for row in rows.values())),
            'violations_per_pr81_orbit': counts, 'single_ray_repairs': repair_report}


def main():
    import sys
    report = verify(json.loads((HERE / 'certificate.json').read_text()))
    if '--write-report' in sys.argv:
        (HERE / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    else:
        require(report == json.loads((HERE / 'report.json').read_text()), 'archived report')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
