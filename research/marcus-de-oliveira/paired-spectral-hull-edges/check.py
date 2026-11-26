"""Exact edge certificates and open-neighborhood rank checks; no numerical solver."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import gcd
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
PERMS = ((0, 1, 2), (1, 2, 0), (2, 0, 1), (0, 2, 1), (1, 0, 2), (2, 1, 0))
PAIRS = tuple(combinations_with_replacement(range(6), 2))
EDGES = tuple(product(range(3), range(3, 6)))
EDGE_PAIRS = tuple(combinations_with_replacement(EDGES, 2))
SIGNS = (1, 1, 1, -1, -1, -1)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def rank(rows):
    """Exact integer elimination, with primitive reduction to bound integers."""
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


def hull(points):
    """Strict CCW hull, starting at the lexicographically smallest point."""
    unique = {}
    for i, p in enumerate(points):
        unique.setdefault(tuple(p), i)
    ids = [unique[p] for p in sorted(unique)]

    def half(indices):
        out = []
        for i in indices:
            while len(out) > 1 and cross(points[out[-2]], points[out[-1]], points[i]) <= 0:
                out.pop()
            out.append(i)
        return out

    if len(ids) <= 1:
        return ids
    return half(ids)[:-1] + half(ids[::-1])[:-1]


def products(a, b):
    answer = []
    for pi in PERMS:
        value = (1, 0)
        for i, j in enumerate(pi):
            value = mul(value, (a[i][0] + b[j][0], a[i][1] + b[j][1]))
        answer.append(value)
    return answer


def generators():
    return [tuple(int(i in e and j in f) + int(i in f and j in e)
                  for i, j in PAIRS) for e, f in EDGE_PAIRS]


def monomial_rows():
    rows = {}
    for column, (i, j) in enumerate(PAIRS):
        assignment = PERMS[i] + PERMS[j]
        for take in product((False, True), repeat=6):
            a = tuple(sum(not take[k] and k % 3 == r for k in range(6)) for r in range(3))
            b = tuple(sum(take[k] and assignment[k] == r for k in range(6)) for r in range(3))
            rows.setdefault(a + b, [0] * 21)[column] += 1
    return rows


def polynomial_space():
    cols = generators()
    rows = monomial_rows()
    primitive = lambda row: tuple(x // gcd(*row) for x in row)
    positive_rows = {primitive(row) for row in rows.values()}
    require(len(rows) == 141, '141 grouped Cauchy-Binet coefficient rows')
    require(len(set(cols)) == 45 and all(primitive(c) in positive_rows for c in cols),
            'each of 45 generators is a positive multiple of a grouped coefficient row')
    kernel = [[int(i == j == v) if i == j else SIGNS[i] * SIGNS[j] * int(v in (i, j))
               for i, j in PAIRS] for v in range(6)]
    require(rank(kernel) == 6 and rank(cols) == 15, 'polynomial space dimensions')
    require(all(dot(k, c) == 0 for k in kernel for c in cols), 'generator kernel')
    require(all(dot(k, row) == 0 for k in kernel for row in rows.values()), 'all polynomial identities')
    require(all(sum(k) == 0 for k in kernel), 'constant row lies in coefficient space')
    tables = {tuple(sum(PERMS[i][r] == c for i in pair) for r in range(3) for c in range(3))
              for pair in PAIRS}
    actual = {tuple(sum(i // 2 == r and p[i] // 2 == c for i in range(6))
                    for r in range(3) for c in range(3)) for p in permutations(range(6))}
    require(len(tables) == 21 and tables == actual, 'all physical pairing contingencies')
    return cols, kernel


def spectrum(value):
    require(isinstance(value, list) and len(value) == 3, 'three spectral values')
    require(all(isinstance(z, list) and len(z) == 2 and
                all(isinstance(x, int) and not isinstance(x, bool) for x in z) for z in value),
            'Gaussian-integer spectra')
    require(cross(*value) != 0, 'noncollinear simple triple')
    return value


def verify(data):
    require(data.get('schema') == 'mdo.paired-spectral-hull-edges.v1', 'schema')
    a, b = spectrum(data['a']), spectrum(data['b'])
    p = products(a, b)
    require(tuple(sum(SIGNS[i] * p[i][j] for i in range(6)) for j in range(2)) == (0, 0),
            'balanced half-products')
    hp = hull(p)
    require(hp == data['half_hull'] and len(hp) >= 4, 'non-triangular half-product hull')
    half_margin = min(cross(p[l], p[r], p[k]) for l, r in zip(hp, hp[1:] + hp[:1])
                      for k in range(6) if k not in (l, r))
    require(half_margin > 0, 'strict half-product hull combinatorics')
    v = [mul(p[i], p[j]) for i, j in PAIRS]
    hv = hull(v)
    require(hv == data['product_hull'] and len(hv) >= 3, 'strict product hull')
    cols, kernel = polynomial_space()
    require(len(data['edges']) == len(hv), 'every polygon edge covered exactly once')
    reports = []
    for record, (l, r) in zip(data['edges'], zip(hv, hv[1:] + hv[:1])):
        require(record['ends'] == [l, r], 'edge order and orientation')
        gap = [cross(v[l], v[r], x) for x in v]
        require([k for k, x in enumerate(gap) if x == 0] == sorted((l, r)), 'exactly two exposed coordinates')
        require(all(x > 0 for k, x in enumerate(gap) if k not in (l, r)), 'strict inward edge gaps')
        allowed = [k for k, c in enumerate(cols) if c[l] == c[r] == 0]
        values = record['weights']
        require(isinstance(values, list) and len(values) == 45 and all(isinstance(x, str) for x in values),
                '45 exact rational weight strings')
        weights = [F(x) for x in values]
        require(all((w > 0 if k in allowed else w == 0) for k, w in enumerate(weights)),
                'strictly positive face representation')
        require([sum(w * c[k] for w, c in zip(weights, cols)) for k in range(21)] == gap,
                'exact edge certificate identity')
        face_rank = rank([cols[k] for k in allowed])
        coordinate_rows = [[int(i == j) for i in range(21)] for j in (l, r)]
        intersection_dimension = 21 - rank(kernel + coordinate_rows)
        require(face_rank == intersection_dimension == 13, 'full edge face span')
        require(all(dot(k, gap) == 0 for k in kernel), 'edge gap belongs to polynomial space')
        # One fixed column basis gives the continuous right inverse used in the proof.
        basis = []
        for k in allowed:
            if rank([cols[j] for j in basis + [k]]) > len(basis):
                basis.append(k)
        require(len(basis) == 13, 'fixed face basis')
        reports.append({'ends': [l, r], 'pair_labels': [list(PAIRS[l]), list(PAIRS[r])],
                        'allowed_generators': len(allowed), 'face_rank': face_rank,
                        'intersection_dimension': intersection_dimension,
                        'minimum_positive_gap': min(x for x in gap if x),
                        'minimum_positive_weight': str(min(w for w in weights if w)),
                        'fixed_basis': basis})
    return {'status': 'PASS', 'half_products': [list(z) for z in p],
            'half_hull': hp, 'half_hull_minimum_strict_gap': half_margin,
            'product_hull': hv, 'product_vertices': [list(v[k]) for k in hv],
            'all_relative_unitaries': True, 'open_spectral_neighborhood': True,
            'explicit_neighborhood_radius': None, 'edges': reports}


def main():
    import sys
    report = verify(json.loads((HERE / 'certificate.json').read_text()))
    if '--write-report' in sys.argv:
        (HERE / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    else:
        require(report == json.loads((HERE / 'report.json').read_text()), 'archived exact report')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
