"""Replay the two pointwise pencil certificates and the signed global determinant."""
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'certificates'))
import verify as v

Q = v.QComplex
N = 9
PIVOT = 4
A = [Q(2, 1), Q(-1, 3), Q(3, -2), Q(0, 1), Q(1, -1),
     Q(1, 2), Q(-2, -1), Q(1), Q(3, 1)]
B = [Q(1, -2), Q(-3, 1), Q(1, 3), Q(2, -1), Q(-1),
     Q(-1, 2), Q(2, 2), Q(3), Q(-1, -2)]


def reflection(n):
    return [[Q(int(i == j) - Fraction(2, n)) for j in range(n)] for i in range(n)]


def embed(block, start, n):
    result = v.identity(n)
    for i, row in enumerate(block):
        for j, entry in enumerate(row):
            result[start + i][start + j] = entry
    return result


def value(a, b, u):
    conjugated = v.multiply(v.multiply(u, v.diagonal(b)), v.adjoint(u))
    return v.determinant([[entry + (a[i] if i == j else v.ZERO)
                           for j, entry in enumerate(row)] for i, row in enumerate(conjugated)])


def pencil(left, a, b, h):
    def at(t):
        return value(a, b[:-1] + [t], v.adjoint(h)) if left else value([t] + a[1:], b, h)
    constant = at(v.ZERO)
    return constant, at(v.ONE) - constant


def vertices(left, a, b):
    # pi maps a local A index to the local B index paired with it.
    for pi in itertools.permutations(range(len(a))):
        i = pi.index(len(b) - 1) if left else 0
        slope = math.prod((a[k] + b[j] for k, j in enumerate(pi) if k != i), start=v.ONE)
        leftover = a[i] if left else b[pi[i]]
        yield pi, (leftover * slope, slope)


def check_certificate(left, a, b, h, weights):
    target = pencil(left, a, b, h)
    options = dict(vertices(left, a, b))
    masses = {tuple(row['permutation']): Fraction(row['mass']) for row in weights}
    assert len(masses) == len(weights) <= 5
    assert min(masses.values()) > 0 and sum(masses.values()) == 1
    for k in range(2):
        assert sum((mass * options[pi][k] for pi, mass in masses.items()), v.ZERO) == target[k]
    return masses


def global_pair(p, q):
    out = [None] * N
    for i, j in enumerate(p):
        out[i] = PIVOT + q[0] if j == PIVOT else j
    for i, j in enumerate(q[1:], 1):
        out[PIVOT + i] = PIVOT + j
    assert sorted(out) == list(range(N))
    return out


def check():
    h = reflection(5)
    assert v.multiply(h, v.adjoint(h)) == v.identity(5)
    u = v.multiply(embed(h, 0, N), embed(h, PIVOT, N))
    assert v.multiply(u, v.adjoint(u)) == v.identity(N)
    weights = json.loads((HERE / 'pencil-weights.json').read_text())
    left = check_certificate(True, A[:5], B[:5], h, weights['left'])
    right = check_certificate(False, A[4:], B[4:], h, weights['right'])
    total = v.ZERO
    for p, mass in left.items():
        for q, weight in right.items():
            pi = global_pair(p, q)
            total += mass * weight * math.prod((a + B[j] for a, j in zip(A, pi)), start=v.ONE)
    assert total == value(A, B, u)

    # Zero slopes on the left, right, and both. Check the arrow formula without division.
    for force_left, force_right in itertools.product((False, True), repeat=2):
        a, b = A.copy(), B.copy()
        if force_left:
            a[:5], b[0] = [Q(2)] * 5, Q(-2)
        if force_right:
            b[4:], a[5] = [Q(3)] * 5, Q(-3)
        l0, l1 = pencil(True, a[:5], b[:5], h)
        r0, r1 = pencil(False, a[4:], b[4:], h)
        assert (l1 == v.ZERO) == force_left and (r1 == v.ZERO) == force_right
        assert l0 * r1 + l1 * r0 == value(a, b, u)
    return {'dimension': N, 'leftTerms': len(left), 'rightTerms': len(right),
            'pairings': len(left) * len(right), 'determinant': total.encode(), 'arrowCases': 4}


if __name__ == '__main__':
    print(json.dumps(check()))
