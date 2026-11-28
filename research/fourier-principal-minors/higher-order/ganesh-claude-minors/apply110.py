#!/usr/bin/env python3
"""The window polynomial on maximal flags for base 10 above 11 (order 110).

Above 11 there are four charts, and the chain
    {5,8,9}  <  {5,7,8,9}  <  {3,5,6,7,8,9}  <  {2,3,5,6,7,8,9}
has one member singular in each.  A maximal flag has occupancies 10, 9, ..., 1
along the order 5,8,9,7,3,6,2,0,1,4, so its level sets contain the chain and all
occupancies are distinct.  The second-order report showed that the second-order
term vanishes on every maximal flag in every chart.

1. Random maximal flags: least order j <= 9 of a nonzero window coefficient, per
   chart and over the four charts.  By Theorem A any such j proves the minor
   nonzero.
2. Nested maximal flags (Ferrers shape): all window coefficients vanish, and the
   minor is nonzero by Proposition F; checked modulo a prime l = 1 (mod 110).

Usage: python3 apply110.py PATH_TO_window_BINARY [SAMPLES]      (default 1000)
Standard library only.
"""
import random
import sys

import ho_lib as ho

so = ho.so
M, Q = 10, 11
CHAIN = [[5, 8, 9], [5, 7, 8, 9], [3, 5, 6, 7, 8, 9], [2, 3, 5, 6, 7, 8, 9]]
ORDER = [5, 8, 9, 7, 3, 6, 2, 0, 1, 4]


def show(xs):
    keys = sorted(set(xs), key=lambda k: (k is None, k))
    return ", ".join("%s: %d" % ("none" if k is None else k, xs.count(k)) for k in keys)


def main():
    binary = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    rng = random.Random(1100)
    charts = so.charts(M, Q)
    assert len(charts) == 4
    where = [[i for i, ch in enumerate(charts) if ch.minor(S, S) < 0] for S in CHAIN]
    assert sorted(i for w in where for i in w) == [0, 1, 2, 3]
    L = Q - 1

    fibs = [{a: sorted(rng.sample(range(Q), 10 - j)) for j, a in enumerate(ORDER)} for _ in range(n)]
    for fib in fibs[:20]:
        levels = so.level_sets(fib, Q)
        assert all(S in levels for S in CHAIN)
    per = [[next((j for j, x in enumerate(d) if x >= 0), None) for d in ho.window_c(binary, ch, fibs, L)] for ch in charts]
    assert all(o is None or o >= 2 for p in per for o in p)
    best = [min((v for v in t if v is not None), default=None) for t in zip(*per)]
    print("1. %d random maximal flags.  Least order j <= 9 of a nonzero window coefficient:" % n)
    for i, p in enumerate(per):
        print("   chart %d: %s" % (i, show(p)))
    print("   over the four charts: %s.  Minor proved nonzero by Theorem A for %d of %d." % (show(best), sum(b is not None for b in best), n))

    mod = ho.ModN(M, Q)
    fibs = []
    for _ in range(100):
        o = rng.sample(range(Q), Q)
        fibs.append({a: sorted(o[:10 - j]) for j, a in enumerate(ORDER)})
    per = [[next((j for j, x in enumerate(d) if x >= 0), None) for d in ho.window_c(binary, ch, fibs, L)] for ch in charts]
    blind = sum(all(p[i] is None for p in per) for i in range(len(fibs)))
    nonzero = sum(mod.minor([(a, b) for a in sorted(f) for b in f[a]]) != 0 for f in fibs)
    print("2. 100 nested maximal flags: all window coefficients of order <= 9 vanish in all four charts for %d; minor nonzero modulo %d for %d (Proposition F)." % (blind, mod.l, nonzero))


if __name__ == "__main__":
    main()
