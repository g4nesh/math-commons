"""Exact split-pair zonotope criterion and a strict physical open subclass.

Python 3.9+ stdlib. Gaussian arithmetic is the repository's exact verifier.
No LP solver is needed to replay or evaluate the certified neighborhood.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd, prod
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('split_pair_exact', HERE.parent / 'certificates' / 'verify.py')
exact = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = exact
spec.loader.exec_module(exact)
Q = exact.QComplex
PERMS = tuple(permutations(range(6)))
ROW = (0, 0, 1, 1, 2, 3)
COL = (0, 0, 1, 1, 2, 2)
SIGNS = (1, 1, 1, -1, -1, -1)
P3 = ((0, 1, 2), (1, 2, 0), (2, 0, 1), (0, 2, 1), (1, 0, 2), (2, 1, 0))
PAIRS = tuple(combinations_with_replacement(range(6), 2))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def table(p):
    return tuple(sum(ROW[i] == r and COL[p[i]] == s for i in range(6))
                 for r in range(4) for s in range(3))


def independent_rows(rows):
    pivots, selected = {}, []
    for i, source in enumerate(rows):
        row = list(map(F, source))
        for j, pivot in pivots.items():
            if row[j]:
                factor = row[j]
                row = [x - factor * y for x, y in zip(row, pivot)]
        j = next((j for j, x in enumerate(row) if x), None)
        if j is not None:
            factor = row[j]
            pivots[j] = [x / factor for x in row]
            selected.append(i)
    return selected


def determinant_integer(source):
    a = [list(row) for row in source]
    sign, previous, n = 1, 1, len(a)
    for k in range(n - 1):
        if not a[k][k]:
            j = next((j for j in range(k + 1, n) if a[j][k]), None)
            if j is None:
                return 0
            a[k], a[j] = a[j], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, 'exact Bareiss division')
                a[i][j] = numerator // previous
        previous = pivot
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[-1][-1]


def facet_normals(b):
    result = set()
    for subset in combinations(range(12), 5):
        rows = [[row[j] for row in b] for j in subset]
        normal = [(-1) ** j * determinant_integer([row[:j] + row[j + 1:] for row in rows])
                  for j in range(6)]
        divisor = gcd(*normal)
        if not divisor:
            continue
        normal = [x // divisor for x in normal]
        if next(x for x in normal if x) < 0:
            normal = [-x for x in normal]
        result.add(tuple(normal))
    return sorted(result)


@lru_cache(None)
def structure():
    tables = sorted({table(p) for p in PERMS})
    lookup = {t: i for i, t in enumerate(tables)}
    representatives = {}
    counts = Counter()
    for p in PERMS:
        representatives.setdefault(table(p), p)
        counts[lookup[table(p)]] += 1
    rows = {}
    for column, t in enumerate(tables):
        p = representatives[t]
        for take_b in product((0, 1), repeat=6):
            key = tuple(sum(not take_b[i] and ROW[i] == r for i in range(6)) for r in range(4))
            key += tuple(sum(take_b[i] and COL[p[i]] == s for i in range(6)) for s in range(3))
            rows.setdefault(key, [0] * 33)[column] += 1
    keys = sorted(rows)
    a = [rows[key] for key in keys]
    coarse = {tuple(int(P3[i][r] == s) + int(P3[j][r] == s) for r in range(3) for s in range(3)): k
              for k, (i, j) in enumerate(PAIRS)}
    fibers = [[] for _ in PAIRS]
    for i, t in enumerate(tables):
        fibers[coarse[t[:6] + tuple(t[6 + s] + t[9 + s] for s in range(3))]].append(i)
    doubles = [(parent, f) for parent, f in enumerate(fibers) if len(f) == 2]
    odd = [[row[f[0]] - row[f[1]] for _, f in doubles] for row in a]
    selected = independent_rows(odd)
    b = [odd[i] for i in selected]
    kernel = [[int(i == j == v) if i == j else SIGNS[i] * SIGNS[j] * int(v in (i, j))
               for i, j in PAIRS] for v in range(6)]
    return dict(tables=tables, lookup=lookup, keys=keys, matrix=a, counts=counts,
                fibers=fibers, doubles=doubles, odd=odd, odd_selected=selected,
                b=b, kernel=kernel, normals=facet_normals(b))


def parent_data(q, repairs):
    s = structure()
    require(len(q) == 33 and len(repairs) == 6, '33 coefficients and six repairs')
    coarse = [sum(q[j] for j in f) for f in s['fibers']]
    masses = [coarse[p] + sum(repairs[v] * s['kernel'][v][p] for v in range(6)) for p in range(21)]
    old_odd = [q[f[0]] - q[f[1]] for _, f in s['doubles']]
    target = [dot(row, old_odd) for row in s['b']]
    return coarse, masses, target


def repair_inequalities(q):
    """Return 83 exact inequalities (coefficients, bound) for six repair variables."""
    s = structure()
    coarse, _, target = parent_data(q, [F(0)] * 6)
    result = [([-s['kernel'][v][p] for v in range(6)], coarse[p]) for p in range(21)]
    for h in s['normals']:
        magnitudes = [abs(sum(h[i] * s['b'][i][j] for i in range(6))) for j in range(12)]
        constant = sum(m * coarse[p] for m, (p, _) in zip(magnitudes, s['doubles']))
        coefficients = [-sum(m * s['kernel'][v][p] for m, (p, _) in zip(magnitudes, s['doubles']))
                        for v in range(6)]
        result.extend((coefficients, constant + sign * dot(h, target)) for sign in (-1, 1))
    return result


def lift(q, repairs, imbalances):
    s = structure()
    require(len(imbalances) == 12, 'twelve imbalances')
    _, masses, target = parent_data(q, repairs)
    require(min(masses) >= 0, 'nonnegative parent masses')
    require(all(abs(x) <= masses[p] for x, (p, _) in zip(imbalances, s['doubles'])), 'split-pair bounds')
    require([dot(row, imbalances) for row in s['b']] == target, 'six odd moments')
    result = [F(0)] * 33
    for p, f in enumerate(s['fibers']):
        if len(f) == 1:
            result[f[0]] = masses[p]
    for imbalance, (p, f) in zip(imbalances, s['doubles']):
        result[f[0]], result[f[1]] = (masses[p] + imbalance) / 2, (masses[p] - imbalance) / 2
    return result


def fixture(data):
    u = exact.identity(6)
    for record in data['vectors']:
        vector = [Q(F(a), F(b)) for a, b in record]
        norm = sum(x.real * x.real + x.imag * x.imag for x in vector)
        h = [[Q(int(i == j)) - 2 * vector[i] * vector[j].conjugate() / Q(norm)
              for j in range(6)] for i in range(6)]
        u = exact.multiply(h, u)
    return u


def check_unitary(u):
    require(len(u) == 6 and all(len(row) == 6 and all(isinstance(x, Q) for x in row) for row in u), 'exact Gaussian matrix')
    require(exact.multiply(exact.adjoint(u), u) == exact.identity(6), 'exact unitarity')


def canonical(u):
    check_unitary(u)
    s = structure()
    det = exact.determinant(u)
    result = [Q(0)] * 33
    for p in PERMS:
        result[s['lookup'][table(p)]] += parity(p) * prod(u[i][j] for i, j in enumerate(p)) * det.conjugate()
    return [x.real for x in result]


def load():
    return json.loads((HERE / 'certificate.json').read_text())


def weights_at(u, data=None):
    data = load() if data is None else data
    verify(data, minors=False)
    check_unitary(u)
    base = fixture(data)
    radius = F(data['radius'])
    require(all((x - y).real ** 2 + (x - y).imag ** 2 <= radius ** 2
                for row, old in zip(u, base) for x, y in zip(row, old)), 'outside certified ball')
    correction = [F(w) - F(q) for w, q in zip(data['weights'], data['canonical'])]
    result = [q + c for q, c in zip(canonical(u), correction)]
    require(min(result) >= F(data['mass_floor_in_ball']) and sum(result) == 1, 'certified weights')
    return result


def verify(data=None, minors=True):
    data = load() if data is None else data
    s = structure()
    require(data['schema'] == 'mdo.split-pair.zonotope.v1', 'schema')
    require(len(s['tables']) == 33 and len(s['keys']) == 192 and len(independent_rows(s['matrix'])) == 21, 'coefficient dimensions')
    require(Counter(map(len, s['fibers'])) == {1: 9, 2: 12}, 'split fibers')
    require(len(s['odd_selected']) == 6 and len(s['normals']) == 31, 'odd dimensions and normals')
    require(data['tables'] == [list(x) for x in s['tables']] and data['keys'] == [list(x) for x in s['keys']], 'coordinate convention')
    require(data['odd_selected_rows'] == s['odd_selected'] and data['odd_matrix'] == s['b'], 'odd rows')
    require(data['normals'] == [list(x) for x in s['normals']], 'all facet normals')
    require(len(independent_rows(s['kernel'])) == 6, 'six even repairs')
    even_matrix = [[sum(F(row[j], len(f)) for j in f) for f in s['fibers']] for row in s['matrix']]
    require(len(independent_rows(even_matrix)) == 15, 'complete even rank')
    for k in s['kernel']:
        even = [F(0)] * 33
        for p, f in enumerate(s['fibers']):
            for j in f:
                even[j] = F(k[p], len(f))
        require(all(dot(row, even) == 0 for row in s['matrix']), 'polarized even kernel')
    u = fixture(data)
    require(data['U'] == [[x.encode() for x in row] for row in u], 'fixture entries')
    q = canonical(u)
    require(q == list(map(F, data['canonical'])), 'canonical coefficients')
    weights = list(map(F, data['weights']))
    floor = F(data['floor'])
    require(min(weights) == floor == F(1, 1000) and sum(weights) == sum(q) == 1, 'strict law floor')
    repairs = list(map(F, data['even_repairs']))
    imbalances = list(map(F, data['imbalances']))
    require(lift(q, repairs, imbalances) == weights, 'lift recovers the exact law')
    require(all(dot(row, [x - y for x, y in zip(weights, q)]) == 0 for row in s['matrix']), 'all 192 polynomial coefficients')
    inequalities = repair_inequalities(q)
    require(len(inequalities) == 83 and all(dot(row, repairs) <= bound for row, bound in inequalities), 'finite zonotope criterion')
    checked = 0
    if minors:
        physical = {key: F(0) for key in s['keys']}
        for size in range(7):
            for rows in combinations(range(6), size):
                for cols in combinations(range(6), size):
                    value = exact.determinant([[u[i][j] for j in cols] for i in rows])
                    key = tuple(sum(i not in rows and ROW[i] == r for i in range(6)) for r in range(4))
                    key += tuple(sum(COL[j] == c for j in cols) for c in range(3))
                    physical[key] += value.real ** 2 + value.imag ** 2
                    checked += 1
        require(checked == 924 and all(dot(row, weights) == physical[key] for key, row in zip(s['keys'], s['matrix'])), 'independent squared-minor coefficients')
    maximum = max(s['counts'].values())
    lipschitz = 21 * maximum
    radius = F(data['radius'])
    require(maximum == 32 and lipschitz == 672 and radius == floor / (2 * lipschitz), 'uniform radius')
    require(F(data['mass_floor_in_ball']) == floor - lipschitz * radius == F(1, 2000), 'neighborhood floor')
    return {'tables': 33, 'monomials': 192, 'coefficient_rank': 21, 'even_rank': 15, 'odd_rank': 6,
            'zonotope_normal_directions': 31, 'repair_inequalities': 83,
            'canonical_negative_masses': sum(x < 0 for x in q),
            'zero_unitary_entries': sum(x == Q(0) for row in u for x in row),
            'base_mass_floor': str(floor), 'squared_minors_checked': checked,
            'maximum_permutations_per_table': maximum, 'entrywise_lipschitz': lipschitz,
            'closed_entrywise_unitary_ball_radius': str(radius),
            'mass_floor_throughout_ball': data['mass_floor_in_ball']}


if __name__ == '__main__':
    report = verify()
    if '--write-report' in sys.argv:
        (HERE / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    else:
        require(report == json.loads((HERE / 'report.json').read_text()), 'archived report')
    print(json.dumps(report, indent=2))
