"""Exact common permutation law at F7. Requires only Python 3 standard library."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb

N = 7
LABELS = range(N)
# The integer triple encodes 240 times a mass in the basis 1, alpha, alpha^2.
ORBITS = [
    ((6, 5, 4, 3, 2, 1, 0), (6, 0, 0)),
    ((6, 5, 4, 1, 0, 3, 2), (685, -240, -240)),
    ((6, 5, 2, 0, 4, 1, 3), (-2309, 720, 960)),
    ((6, 4, 1, 5, 3, 0, 2), (1843, -480, -720)),
]
SUBSETS = [list(combinations(LABELS, k)) for k in range(N+1)]
SIGNS = [[(p, (-1)**sum(a > b for i, a in enumerate(p) for b in p[i+1:]))
          for p in permutations(range(k))] for k in range(N+1)]


def square_minor(rows, cols):
    """Integer determinant monomials in zeta_7, then multiply by the conjugate."""
    k = len(cols)
    counts = [0]*N
    for p, sign in SIGNS[k]:
        exponent = sum(rows[p[i]]*col for i, col in enumerate(cols)) % N
        counts[exponent] += sign
    # a[d] is the coefficient of zeta^d in the product with its conjugate.
    a = [sum(v*counts[(i-d) % N] for i, v in enumerate(counts)) for d in LABELS]
    assert a[1:4] == a[6:3:-1]
    # zeta^d+zeta^-d for d=1,2,3 is alpha, alpha^2-2, 1-alpha-alpha^2.
    return tuple(Fraction(c, N**k) for c in
                 (a[0] - 2*a[2] + a[3], a[1] - a[3], a[2] - a[3]))


def orbit_counts(seed):
    """Count p' = L p R^-1 with L(i)=h*i+b, R(j)=e*h^-1*j+c."""
    perms = Counter()
    for h, e, b, c in product(range(1, N), (1, -1), LABELS, LABELS):
        inverse_r = (e*h) % N
        p = tuple((h*seed[(inverse_r*(j-c)) % N]+b) % N for j in LABELS)
        perms[p] += 1
    return perms


def check():
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    # alpha is the unique positive root of x^3+x^2-2x-1.
    lo, hi = Fraction(12469, 10000), Fraction(12470, 10000)
    f = lambda x: x**3+x*x-2*x-1
    assert f(lo) < 0 < f(hi)
    masses = {}
    group_size = 6*2*7*7
    for seed, coeff in ORBITS:
        lower = sum(min(a*lo**i, a*hi**i) for i, a in enumerate(coeff))
        assert lower > 0
        orbit = orbit_counts(seed)
        assert sum(orbit.values()) == group_size
        for p, count in orbit.items():
            row = masses.setdefault(p, [0, 0, 0])
            for j, c in enumerate(coeff):
                row[j] += c*count
    assert tuple(sum(c[j] for _, c in ORBITS) for j in range(3)) == (225, 0, 0)
    checks = 0
    for k, subsets in enumerate(SUBSETS):
        background = Fraction(1, 16*comb(N, k))
        for cols in subsets:
            sums = {rows: [0, 0, 0] for rows in subsets}
            for p, coeff in masses.items():
                dst = sums[tuple(sorted(p[j] for j in cols))]
                for m, value in enumerate(coeff):
                    dst[m] += value
            for rows, value in sums.items():
                observed = tuple(Fraction(a, 240*group_size) + (background if m == 0 else 0)
                                 for m, a in enumerate(value))
                assert observed == square_minor(rows, cols), (rows, cols, observed)
                checks += 1
    assert checks == comb(2*N, N) == 3432
    # The signed-mass Lipschitz bound and half-reserve at the closed ball's edge.
    assert (Fraction(7, 64) + Fraction(1, 49)) / 32768 < Fraction(1, 161280)
    print(f'{checks} exact Fourier squared minors verified at orders 0 through 7; '
          'four positive orbit masses plus 1/16 uniform; perturbation bound verified.')


if __name__ == '__main__':
    check()
