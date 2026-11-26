"""Exact deletion, moment and permutation identities for one nonconcyclic pair."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations, product
import json
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "certificates"))
import verify as v

Q = v.QComplex


def reflector(n, indices):
    matrix = v.identity(n)
    for i, j in product(indices, repeat=2):
        matrix[i][j] -= F(2, len(indices))
    return matrix


def conjugate(unitary, values):
    return v.multiply(v.multiply(unitary, v.diagonal(values)), v.adjoint(unitary))


def shift(matrix, values):
    return [[x + (values[i] if i == j else 0) for j, x in enumerate(row)]
            for i, row in enumerate(matrix)]


def cross_ratio(values):
    a, b, c, d = values
    return (a-b)*(c-d)/((a-c)*(b-d))


def sum_square(x):
    return (x*x.conjugate()).real


def compression(unitary, eigenvalues):
    """Compare the two squared commutator eigenvalues to the centered moments."""
    g = conjugate(unitary, eigenvalues)
    c = [row[1:] for row in g[1:]]
    cc = v.multiply(c, v.adjoint(c))
    cs = v.multiply(v.adjoint(c), c)
    comm = [[a-b for a, b in zip(row, other)] for row, other in zip(cc, cs)]
    p = [sum_square(x) for x in unitary[0]]
    mu = sum((w*x for w, x in zip(p, eigenvalues)), v.ZERO)
    var = sum(w*sum_square(x-mu) for w, x in zip(p, eigenvalues))
    pseudo = sum((w*(x-mu)*(x-mu) for w, x in zip(p, eigenvalues)), v.ZERO)
    defect_squared = var*var - sum_square(pseudo)
    assert 2*defect_squared == sum(sum_square(x) for row in comm for x in row)
    assert sum((comm[i][i] for i in range(len(c))), v.ZERO) == v.ZERO
    assert v.determinant(comm) == v.ZERO  # three-dimensional compression
    return c, defect_squared


def projector_weights(c, roots):
    """Rows are physical coordinates, columns label eigenvalues of the compression."""
    n = len(c)
    rows = [[] for _ in c]
    for root in roots:
        projection = v.identity(n)
        for other in roots:
            if other != root:
                factor = shift(c, [-other]*n)
                projection = v.multiply(projection, [[x/(root-other) for x in row] for row in factor])
        assert v.adjoint(projection) == projection == v.multiply(projection, projection)
        for i, row in enumerate(rows):
            assert projection[i][i].imag == 0
            row.append(projection[i][i].real)
    assert all(sum(row) == 1 for row in rows)
    assert all(sum(column) == 1 for column in zip(*rows))
    return rows


def birkhoff(rows):
    weights = {}
    while sum(rows[0]):
        perm = next(p for p in permutations(range(len(rows))) if all(row[j] > 0 for row, j in zip(rows, p)))
        mass = min(row[j] for row, j in zip(rows, perm))
        weights[perm] = mass
        for row, j in zip(rows, perm):
            row[j] -= mass
    assert sum(weights.values()) == 1 and all(x == 0 for row in rows for x in row)
    return weights


def interpolation(values, low, high, target):
    weight = (values[high]-target)/(values[high]-values[low])
    assert 0 <= weight <= 1
    return ((low, weight), (high, 1-weight))


def lift(reduced):
    real = [F(1), F(3), F(8)]
    weights = defaultdict(F)
    for large, outer in interpolation(real, 1, 2, F(17, 3)):
        remaining = sorted(set(range(3))-{large})
        for small, inner in interpolation(real, *remaining, F(4, 3)):
            embed = [small, large, 3]
            left, = set(range(4))-set(embed)
            for perm, mass in reduced.items():
                weights[(left, *(embed[j] for j in perm))] += inner*outer*mass
    assert all(w > 0 for w in weights.values()) and sum(weights.values()) == 1
    return weights


def run():
    gamma = [Q(1), Q(3), Q(8), Q(4, 1)]
    a = [v.ONE/x for x in gamma]
    u = v.multiply(reflector(4, (1, 2, 3)), reflector(4, (0, 1, 2)))
    assert v.multiply(u, v.adjoint(u)) == v.identity(4)
    c, good = compression(u, gamma)
    _, bad = compression(reflector(4, tuple(range(4))), gamma)
    assert good == 0 and bad > 0
    roots = [Q(F(4, 3)), Q(F(17, 3)), Q(4, 1)]
    assert all(v.determinant(shift(c, [-r]*3)) == v.ZERO for r in roots)
    weights = lift(birkhoff(projector_weights(c, roots)))
    matrix = conjugate(u, a)

    def check(b):
        direct = v.determinant(shift(matrix, b))
        mixture = sum((mass*prod((a[j]+y for j, y in zip(p, b)), start=v.ONE)
                       for p, mass in weights.items()), v.ZERO)
        assert mixture == direct
        return direct

    # Multi-affinity promotes the eight exact grid checks to an identity in
    # all three remaining complex eigenvalues, including zeros and collisions.
    for choices in product((v.ZERO, v.ONE), repeat=3):
        check([v.ZERO, *choices])
    b = [v.ZERO, Q(1, 1), Q(2, -1), Q(-1, 3)]
    assert cross_ratio(a).imag and cross_ratio(b).imag
    assert all(matrix[0][j] != v.ZERO for j in range(1, 4))
    return dict(squaredDefect=str(bad), determinant=check(b).encode(),
                crossRatioA=cross_ratio(a).encode(), crossRatioB=cross_ratio(b).encode(),
                weights={"".join(str(i+1) for i in p): str(w) for p, w in sorted(weights.items())})


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
