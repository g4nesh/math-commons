"""Exact four-facet proof for paired six-dimensional spectra. Python standard library."""
import hashlib
import json
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(4), 2))
ROW_PAIRS = tuple(combinations(range(6), 2))
PERMS = ((0, 1, 2), (1, 2, 0), (2, 0, 1), (0, 2, 1), (1, 0, 2), (2, 1, 0))
# PR81, orbit numbers four through seven, in its parity order.
RAYS = (
    ((1,2,1),(1,5,1),(2,4,1),(4,4,1),(5,5,1)),
    ((1,2,2),(1,5,2),(2,3,1),(2,4,1),(3,4,1),(5,5,2)),
    ((0,2,1),(0,5,1),(1,1,1),(1,4,1),(2,3,1),(3,4,1),(5,5,1)),
    ((0,2,1),(0,5,1),(1,2,1),(1,4,1),(2,3,2),(3,3,2),(4,5,1)),
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


class Polynomial:
    """Sparse integer polynomials. A key holds five base-64 exponent digits.

    X and Y are affine. Each matrix entry has total degree at most
    2+2+2+4 = 10. No determinant intermediate exceeds degree 60. Keys
    can therefore be added without a carry between exponent digits.
    """
    def __init__(self, terms):
        self.terms = {key: value for key, value in terms.items() if value}

    def __add__(self, other):
        terms = self.terms.copy()
        pairs = other.terms.items() if isinstance(other, Polynomial) else [(0, other)]
        for key, value in pairs:
            terms[key] = terms.get(key, 0) + value
        return Polynomial(terms)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            return Polynomial({key: value * other for key, value in self.terms.items()})
        terms = {}
        for left, a in self.terms.items():
            for right, b in other.terms.items():
                key = left + right
                terms[key] = terms.get(key, 0) + a*b
        return Polynomial(terms)

    __rmul__ = __mul__


def parity(indices):
    return (-1) ** sum(a > b for i, a in enumerate(indices) for b in indices[i+1:])


def transpose(matrix):
    return list(zip(*matrix))


def multiply(left, right):
    columns = transpose(right)
    return [[sum(x*y for x,y in zip(row, column)) for column in columns] for row in left]


def minor(matrix, rows, cols):
    i,j = rows
    k,l = cols
    return matrix[i][k]*matrix[j][l] - matrix[i][l]*matrix[j][k]


def compound(matrix, rows, cols):
    return [[minor(matrix, r, c) for c in cols] for r in rows]


def frames(a, b, c, d, e):
    first = [[1,0],[0,1],[a,0],[0,b],[c,d],[0,e]]
    complement = [[-a,0,-c,0],[0,-b,-d,-e],[1,0,0,0],
                  [0,1,0,0],[0,0,1,0],[0,0,0,1]]
    return first, complement


def table_of_permutations(i, j):
    return tuple(int(PERMS[i][r]==s)+int(PERMS[j][r]==s)
                 for r in range(3) for s in range(3))


def table_of_row_sets(sets):
    return tuple(sum(t//2==r for t in indices) for r in range(3) for indices in sets)


def facet_matrix(first, complement, ray):
    """Twice the real symmetric matrix in the numerator, before division by det(S)."""
    scores = {table_of_permutations(i,j): v for i,j,v in ray}
    gram = multiply(transpose(complement), complement)
    star = [[parity(t+s) if len(set(t+s))==4 else 0 for t in PAIRS] for s in PAIRS]
    y2 = compound(complement, ROW_PAIRS, PAIRS)
    tail = multiply(multiply(y2, star), compound(gram, PAIRS, PAIRS))
    first_minors = [minor(first, i, (0,1)) for i in ROW_PAIRS]
    terms = []
    for i, left in enumerate(ROW_PAIRS):
        for j, middle in enumerate(ROW_PAIRS):
            if set(left).intersection(middle):
                continue
            right = tuple(k for k in range(6) if k not in left+middle)
            weight = scores.get(table_of_row_sets((left,middle,right)), 0)
            if weight:
                scale = weight * parity(left+middle+right) * first_minors[i]
                terms.append((scale, y2[j], tail[ROW_PAIRS.index(right)]))
    raw = [[sum(scale * v[i] * w[j] for scale,v,w in terms)
            for j in range(6)] for i in range(6)]
    return [[raw[i][j]+raw[j][i] for j in range(6)] for i in range(6)]


def positivity_certificate(matrix):
    @lru_cache(None)
    def determinant(rows, cols):
        if not rows:
            return 1
        return sum((-1)**j * matrix[rows[0]][col] *
                   determinant(rows[1:], cols[:j]+cols[j+1:])
                   for j,col in enumerate(cols))

    result = []
    for size in range(1,7):
        indices = tuple(range(size))
        poly = determinant(indices, indices)
        rows = []
        for key, coefficient in poly.terms.items():
            powers = [(key >> (6*k)) % 64 for k in range(5)]
            require(key < 64**5 and sum(powers)<=10*size, 'degree bound')
            require(all(p%2==0 for p in powers) and coefficient>0, 'positive even monomial')
            rows.append(powers+[coefficient])
        require(rows, 'zero leading minor')
        payload = json.dumps(sorted(rows), separators=(',', ':')).encode('ascii')
        result.append({'order':size, 'monomials':len(rows),
                       'minimum_coefficient':min(r[-1] for r in rows),
                       'sha256':hashlib.sha256(payload).hexdigest()})
    return result


def compute():
    source = json.loads((HERE.parent/'paired-pointwise-cone'/'certificate.json').read_text())
    require(source['pr81_orbit_representatives'][3:] ==
            [[list(term) for term in ray] for ray in RAYS], 'orbit representatives')
    variables = [Polynomial({64**i:1}) for i in range(5)]
    first, complement = frames(*variables)
    return {str(i):positivity_certificate(facet_matrix(first, complement, ray))
            for i,ray in enumerate(RAYS, 4)}


if __name__ == '__main__':
    result = compute()
    require(result==json.loads((HERE/'grassmann-results.json').read_text()),
            'coefficient fingerprints differ from the archive')
    for i,data in result.items():
        print('Orbit',i,'positive even monomial counts:',[x['monomials'] for x in data])
