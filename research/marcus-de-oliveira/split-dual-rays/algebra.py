"""Exact reconstruction, ray tests and physical relabelings for the split cone."""
from fractions import Fraction as F
from itertools import permutations, product
from math import gcd, lcm

ROW = (0, 0, 1, 1, 2, 3)
COL = (0, 0, 1, 1, 2, 2)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rref(source):
    matrix = [list(map(F, row)) for row in source]
    pivots, r = [], 0
    for c in range(len(matrix[0])):
        k = next((k for k in range(r, len(matrix)) if matrix[k][c]), None)
        if k is None:
            continue
        matrix[r], matrix[k] = matrix[k], matrix[r]
        scale = matrix[r][c]
        matrix[r] = [x / scale for x in matrix[r]]
        for k in range(len(matrix)):
            if k != r and matrix[k][c]:
                scale = matrix[k][c]
                matrix[k] = [x - scale * y for x, y in zip(matrix[k], matrix[r])]
        pivots.append(c)
        r += 1
        if r == len(matrix):
            break
    return matrix, pivots


def primitive(values):
    values = list(map(F, values))
    scale = lcm(*(x.denominator for x in values))
    ints = [int(x * scale) for x in values]
    divisor = gcd(*ints)
    require(divisor != 0, 'nonzero ray')
    return tuple(x // divisor for x in ints)


def structure():
    representatives = {}
    for p in permutations(range(6)):
        table = tuple(sum(ROW[i] == r and COL[p[i]] == s for i in range(6))
                      for r in range(4) for s in range(3))
        representatives.setdefault(table, p)
    tables = sorted(representatives)
    require(len(tables) == 33, '33 physical tables')
    rows = {}
    for j, table in enumerate(tables):
        p = representatives[table]
        for take in product((0, 1), repeat=6):
            key = tuple(sum(not take[i] and ROW[i] == r for i in range(6)) for r in range(4))
            key += tuple(sum(take[i] and COL[p[i]] == s for i in range(6)) for s in range(3))
            rows.setdefault(key, [0] * 33)[j] += 1
    keys = sorted(rows)
    matrix = [rows[key] for key in keys]
    reduced, pivots = rref(matrix)
    require(len(keys) == 192 and len(pivots) == 21, 'coefficient matrix shape and rank')
    kernel = []
    for c in range(33):
        if c in pivots:
            continue
        v = [F(0)] * 33
        v[c] = 1
        for i, j in enumerate(pivots):
            v[j] = -reduced[i][c]
        kernel.append(primitive(v))
    require(len(kernel) == 12, 'kernel dimension')
    require(all(sum(x * y for x, y in zip(k, row)) == 0 for k in kernel for row in matrix),
            'orthogonality of the exact kernel')
    return tables, keys, matrix, kernel


def actions(tables):
    lookup = {t: i for i, t in enumerate(tables)}
    result = []
    for pair_swap, single_swap in product((0, 1), repeat=2):
        row = (1, 0, 2, 3) if pair_swap else (0, 1, 2, 3)
        if single_swap:
            row = row[:2] + (3, 2)
        for col in permutations(range(3)):
            action = []
            for table in tables:
                changed = [0] * 12
                for r in range(4):
                    for s in range(3):
                        changed[3 * row[r] + col[s]] = table[3 * r + s]
                action.append(lookup[tuple(changed)])
            result.append(tuple(action))
    require(len(set(result)) == 24, 'physical group order')
    return result


def act(vector, permutation):
    out = [0] * len(vector)
    for i, j in enumerate(permutation):
        out[j] = vector[i]
    return tuple(out)


def verify_ray(ray, kernel):
    require(len(ray) == 33 and all(isinstance(x, int) and not isinstance(x, bool) and x >= 0 for x in ray),
            '33 nonnegative integer coordinates')
    require(any(ray) and gcd(*ray) == 1, 'primitive nonzero ray')
    require(all(sum(x * y for x, y in zip(k, ray)) == 0 for k in kernel), 'row-space membership')
    support = [i for i, x in enumerate(ray) if x]
    restricted = [[row[i] for i in support] for row in kernel]
    require(len(support) - len(rref(restricted)[1]) == 1, 'one-dimensional support kernel')


def invariant(ray, tables, group):
    lookup = {t: i for i, t in enumerate(tables)}
    swap = [lookup[t[:6] + t[9:12] + t[6:9]] for t in tables]
    orbit = {act(ray, g) for g in group}
    return {'support': sum(bool(x) for x in ray), 'mass': sum(ray), 'maximum': max(ray),
            'singleton_equal_mass': sum(x for x, t in zip(ray, tables) if t[6:9] == t[9:12]),
            'singleton_imbalance_l1': sum(abs(ray[i] - ray[swap[i]]) for i in range(33)),
            'orbit_size': len(orbit), 'stabilizer_size': 24 // len(orbit)}


def read_lrs(path):
    rays = set()
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 34 or fields[0] != '1':
            continue
        values = tuple(map(F, fields[1:]))
        require(sum(values) == 1, 'normalized lrs vertex')
        rays.add(primitive(values))
    return rays
