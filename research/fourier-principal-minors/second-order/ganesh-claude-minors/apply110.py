#!/usr/bin/env python3
"""The second-order criterion on a full-cover flag for base 10 above 11 (order 110).

Above 11 the field Q(zeta_10) has four charts (11 = 1 mod 10).  The chain
    {5,8,9}  <  {5,7,8,9}  <  {3,5,6,7,8,9}  <  {2,3,5,6,7,8,9}
has one member singular in each chart, so every index set K in Z_10 x Z_11
whose level sets contain it is first-order obstructed (frontier report; this is
the situation of Wang-Zhang, Proposition 6.2).

1. Maximal flag: occupancies 10, 9, ..., 1 along the order 5,8,9,7,3,6,2,0,1,4.
   All occupancies are distinct, so by Corollary 4 the second-order term must
   vanish in every chart.  Checked on random K.
2. A flag with ties: the same chain as level sets, but equal occupancies inside
   the shells next to the singular sets.  The second-order term is then usually
   nonzero in some chart, which proves those minors nonzero.

Usage: python3 apply110.py [SAMPLES]      (default 1000)
Standard library only.
"""
import random
import sys

import so_lib as so

M, Q = 10, 11
CHAIN = [[5, 8, 9], [5, 7, 8, 9], [3, 5, 6, 7, 8, 9], [2, 3, 5, 6, 7, 8, 9]]
ORDER = [5, 8, 9, 7, 3, 6, 2, 0, 1, 4]
TIED = {5: 9, 8: 9, 9: 9, 7: 8, 3: 6, 6: 6, 2: 5, 0: 3, 1: 3, 4: 3}


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    rng = random.Random(110)
    charts = so.charts(M, Q)
    assert len(charts) == 4
    where = [[i for i, ch in enumerate(charts) if ch.minor(S, S) < 0] for S in CHAIN]
    assert sorted(i for w in where for i in w) == [0, 1, 2, 3]
    print("chain members and the chart in which each is singular: " + "; ".join("%s -> chart %d" % (S, w[0]) for S, w in zip(CHAIN, where)))

    for _ in range(50):
        fib = {a: sorted(rng.sample(range(Q), 10 - j)) for j, a in enumerate(ORDER)}
        levels = so.level_sets(fib, Q)
        assert all(S in levels for S in CHAIN)
        for ch in charts:
            e, pre, A0, Lam = so.second_order(ch, fib)
            assert A0 < 0 and Lam < 0
    print("1. maximal flag (occupancies 10, 9, ..., 1): A0 = 0 and Lambda = 0 in all four charts for 50 random K, as Corollary 4 requires.")

    per_chart, resolved = [0, 0, 0, 0], 0
    for _ in range(n):
        fib = {a: sorted(rng.sample(range(Q), c)) for a, c in TIED.items()}
        levels = so.level_sets(fib, Q)
        assert all(S in levels for S in CHAIN)
        out = [so.second_order(ch, fib) for ch in charts]
        assert all(o[2] < 0 for o in out)
        nz = [o[3] >= 0 for o in out]
        for i, z in enumerate(nz):
            per_chart[i] += z
        resolved += any(nz)
    print("2. flag with ties (occupancies %s): A0 = 0 in all four charts for every sample; Lambda != 0 in charts 0..3 for %s of %d samples; "
          "nonzero in at least one chart, hence minor nonzero, for %d of %d." % (TIED, per_chart, n, resolved, n))


if __name__ == "__main__":
    main()
