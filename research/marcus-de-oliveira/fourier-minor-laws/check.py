"""Replay a positive common minor law at F5 and at a dense unitary perturbation.

All arithmetic uses integer polynomials modulo the fifth cyclotomic polynomial,
with a common integer denominator. Run with Python 3 and no dependencies.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


class Cyclo:
    """Q[z]/(1+z+z^2+z^3+z^4), with the embedding z=exp(2 pi i / 5)."""

    def __init__(self, coeff=(0,), den=1):
        c = [0] * 5
        for k, a in enumerate(coeff):
            c[k % 5] += a
        c = [x - c[4] for x in c[:4]]
        g = gcd(den, *c)
        if den == 0:
            raise ZeroDivisionError
        if den < 0:
            g = -g
        self.c = tuple(a // g for a in c)
        self.den = den // g

    @staticmethod
    def scalar(x):
        if isinstance(x, Cyclo):
            return x
        r = Fraction(x)
        return Cyclo((r.numerator,), r.denominator)

    def __add__(self, other):
        r = self.scalar(other)
        return Cyclo([a*r.den + b*self.den for a, b in zip(self.c, r.c)], self.den*r.den)

    __radd__ = __add__

    def __neg__(self):
        return Cyclo([-a for a in self.c], self.den)

    def __sub__(self, other):
        return self + -self.scalar(other)

    def __rsub__(self, other):
        return self.scalar(other) + -self

    def __mul__(self, other):
        r = self.scalar(other)
        c = [0]*7
        for i, a in enumerate(self.c):
            for j, b in enumerate(r.c):
                c[i+j] += a*b
        return Cyclo(c, self.den*r.den)

    __rmul__ = __mul__

    def conjugate(self):
        a, b, c, d = self.c
        return Cyclo([a, 0, d, c, b], self.den)

    def real_part(self):
        return (self + self.conjugate()) * Fraction(1, 2)

    def __eq__(self, other):
        r = self.scalar(other)
        return self.c == r.c and self.den == r.den

    def positive(self):
        """Exact sign of a real element, using z^2+z^3=-(1+sqrt(5))/2."""
        a, b, c, d = self.c
        assert b == 0 and c == d
        h, k = 2*a-c, -c  # numerator of h + k sqrt(5)
        if h >= 0 and k >= 0:
            return h+k > 0
        if h <= 0 and k <= 0:
            return False
        if h > 0:
            return h*h > 5*k*k
        return 5*k*k > h*h


def product(seq):
    out = Cyclo.scalar(1)
    for x in seq:
        out *= x
    return out


def parity(p):
    return (-1)**sum(x > y for i, x in enumerate(p) for y in p[i+1:])


N = range(5)
PERMS = list(permutations(N))
SIGNS = {p: parity(p) for p in PERMS}
AFFINE = {tuple((h*j+c) % 5 for j in N) for c in N for h in (2, -2)}
Z = [Cyclo([0]*j + [1]) for j in N]
Q = (1 + 2*(Z[1]+Z[4])) * Fraction(1, 5)
F = [[Q*Z[(i*j) % 5] for j in N] for i in N]


def minor(matrix, rows, cols):
    return sum((parity(p) * product(matrix[rows[p[j]]][col] for j, col in enumerate(cols))
                for p in permutations(range(len(rows)))), Cyclo())


def square(z):
    return z*z.conjugate()


def base_law():
    return {p: (1-Q)*Fraction(1, 120) + (Q*Fraction(1, 10) if p in AFFINE else 0)
            for p in PERMS}


def signed_law(matrix):
    factor = minor(matrix, N, N).conjugate()
    return {p: (factor * SIGNS[p] * product(matrix[p[j]][j] for j in N)).real_part()
            for p in PERMS}


def check(matrix, weights):
    assert sum(weights.values()) == 1
    assert all((w - Fraction(1, 6144)).positive() for w in weights.values())
    assert all(sum(matrix[i][j]*matrix[k][j].conjugate() for j in N) == int(i == k)
               for i in N for k in N)
    count = 0
    for k in range(6):
        subsets = list(combinations(N, k))
        for cols in subsets:
            masses = {rows: Cyclo() for rows in subsets}
            for p, w in weights.items():
                rows = tuple(sorted(p[j] for j in cols))
                masses[rows] += w
            for rows in subsets:
                assert masses[rows] == square(minor(matrix, rows, cols)), (rows, cols)
                count += 1

    for a, b in [([0, 1, Z[1], Z[2], 2+Z[3]], [1+Z[1], 0, Z[3], -2, 3-Z[2]]),
                 ([0, 0, 0, 1, Z[1]], [2, 2, Z[1], Z[1], Z[4]])]:
        target = [[(a[i] if i == k else 0) +
                   sum(matrix[i][j]*b[j]*matrix[k][j].conjugate() for j in N)
                   for k in N] for i in N]
        assert minor(target, N, N) == sum(w*product(a[p[j]] + b[j] for j in N)
                                         for p, w in weights.items())
    assert count == 252
    return count


def main():
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    assert Q.positive() and Q*Q*5 == 1 and len(AFFINE) == 10
    base = base_law()
    t = Fraction(1, 256)
    c, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    moved = [row[:] for row in F]
    moved[0] = [c*x - s*y for x, y in zip(F[0], F[1])]
    moved[1] = [s*x + c*y for x, y in zip(F[0], F[1])]
    assert any(moved[0][j] != F[0][j] for j in N)
    assert all((Fraction(1, 128**2) - square(moved[i][j] - F[i][j])).positive()
               for i in N for j in N)
    before, after = signed_law(F), signed_law(moved)
    moved_law = {p: w + after[p] - before[p] for p, w in base.items()}
    for name, u, weights in [('Fourier F5', F, base), ('rotated F5', moved, moved_law)]:
        print(f'{name}: {check(u, weights)} exact minor identities, 120 positive weights, '
              '2 exact determinant identities')


if __name__ == '__main__':
    main()
