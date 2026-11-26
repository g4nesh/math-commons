"""Exact common laws for F4, F6, F8 and an order-two obstruction for F9."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, gcd

if not __debug__:
    raise SystemExit(
        "The Fourier certificate checker requires assertions. From the repository root run "
        "python3 -E -B research/marcus-de-oliveira/fourier-composite-laws/check.py"
    )

# Coefficients are (r, s), representing mass (r + s*sqrt(2))/denominator.
LAWS = {
    4: (1, [((3, 2, 0, 1), (1, 0))]),
    6: (8, [((5, 4, 1, 0, 3, 2), (2, 0)), ((5, 4, 0, 2, 1, 3), (1, 0)),
            ((5, 3, 1, 4, 0, 2), (5, 0))]),
    8: (4, [((7, 6, 5, 1, 0, 3, 2, 4), (3, -2)),
            ((7, 6, 4, 0, 2, 5, 1, 3), (-4, 4)),
            ((7, 6, 4, 0, 2, 1, 5, 3), (3, -2)),
            ((7, 5, 2, 0, 4, 6, 1, 3), (2, 0))]),
}


def averaged_weights(n, seeds):
    units = [h for h in range(n) if gcd(h, n) == 1]
    labels = range(n)
    group_size = len(units)*2*n*n
    weights = {}
    for seed, (r, s) in seeds:
        assert sorted(seed) == list(labels)
        assert min(r+s*Fraction(1414, 1000), r+s*Fraction(1415, 1000)) > 0
        for h, e, b, c in product(units, (-1, 1), labels, labels):
            # R(j)=e/h*j+c, L(i)=h*i+b; sample L p R^-1, counting repetitions.
            p = tuple((h*seed[(e*h*(j-c)) % n]+b) % n for j in labels)
            a, z = weights.get(p, (0, 0))
            weights[p] = (a+r, z+s)
    return weights, group_size


def minor_numerator(n, rows, cols, swaps):
    # Compute det of the unnormalized minor as a sum of powers of zeta_n.
    count = [0]*n
    for p, sign in swaps:
        power = sum(rows[p[i]]*j for i, j in enumerate(cols))
        count[power % n] += sign
    # Multiply by conjugate; split the symmetric terms into 1 and sqrt(2).
    a = [sum(x*count[(i-d) % n] for i, x in enumerate(count)) for d in range(n)]
    assert all(a[d] == a[-d] for d in range(1, n))
    rational = a[0]-a[n//2]
    if n == 4:
        return rational, 0
    if n == 6:
        return rational+a[1]-a[2], 0
    assert n == 8
    return rational, a[1]-a[3]


def check_law(n):
    denominator, seeds = LAWS[n]
    assert tuple(sum(c[j] for _, c in seeds) for j in (0, 1)) == (denominator, 0)
    weights, group_size = averaged_weights(n, seeds)
    checks = 0
    for k in range(n+1):
        subs = list(combinations(range(n), k))
        swaps = [(p, (-1)**sum(a > b for i, a in enumerate(p) for b in p[i+1:]))
                 for p in permutations(range(k))]
        for cols in subs:
            totals = {rows: [0, 0] for rows in subs}
            for p, (r, s) in weights.items():
                hit = totals[tuple(sorted(p[j] for j in cols))]
                hit[0] += r
                hit[1] += s
            for rows in subs:
                numerator = minor_numerator(n, rows, cols, swaps)
                assert all(totals[rows][i]*n**k == numerator[i]*group_size*denominator
                           for i in (0, 1)), (n, rows, cols)
                checks += 1
    assert checks == comb(2*n, n)
    print(f'F{n}: {checks} exact minor probabilities, {len(seeds)} seed orbit(s).')


def check_nine():
    # Output translations preserve all paired differences and send p(0) to 0.
    maxima = {}
    admissible = 0
    pairs = list(combinations(range(9), 2))
    for tail in permutations(range(1, 9)):
        p = (0,) + tail
        counts = Counter(((j-k)*(p[j]-p[k])) % 9 for j, k in pairs)
        t0 = counts[0]
        difference = counts[4]+counts[5]-counts[2]-counts[7]
        maxima[t0] = max(maxima.get(t0, -36), difference)
        admissible += t0 == 0
    assert admissible == 5184
    assert maxima == {0: 6, 2: 9, 3: 10, 5: 8, 9: 9}
    assert all(2*d-3*z <= 12 for z, d in maxima.items())
    # alpha=2cos(2pi/9) solves x^3-3x+1=0 on the increasing branch 1<x<2.
    x = Fraction(61, 40)
    assert x**3-3*x+1 < 0 and 3*(2*x-3)*(x+2) > Fraction(1, 2)
    assert Fraction(1215, 4096) < Fraction(1, 2)
    print('F9: all 40320 translation representatives checked; '
          'T4-T2-(3/2)T0 <= 6, contradicted throughout the 1/4096 ball.')


if __name__ == '__main__':
    for dimension in LAWS:
        check_law(dimension)
    check_nine()
