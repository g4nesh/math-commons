#!/usr/bin/env python3
"""Validate the leading and the second-order term against direct determinants.

For each configuration the principal minor is computed directly in
F[t]/(t^(q-1)) (zeta_q = 1 + t) and compared with the prediction of
so_lib.second_order:
  * coefficients of t^0 .. t^(e-1) vanish,
  * the coefficient of t^e is pre * A0,
  * if A0 = 0, the coefficient of t^(e+1) is pre * Lam.
It also checks the row-sum identity of Corollary 4 on every shell.
Only configurations with e + 1 <= q - 2 are used, so that both coefficients are
visible in the truncated ring.

Two kinds of configurations:
  random    random small index sets;
  targeted  a carrier S that is singular in the chart is placed as a level set
            of multiplicity one, at the bottom level (upper shell only) or at
            the second level (upper and lower shell).
  tall      for (11, 23): a singular carrier of rank 3 at the third level, so that
            the levels c = 2 and c = 3 of Theorem 3 occur.

Usage: python3 validate.py
Standard library only; about half a minute.
"""
import itertools
import random
from math import comb

import so_lib as so


def check_row_sums(ch, fib, stats):
    """Corollary 4: for every level c and every a in the shell T_(c-1) minus T_c,
    sum over a2 in the shell of Cof(a, a2) * E(a, a2) equals d(T_(c-1)) * d(T_c)."""
    K, q = ch.K, ch.q
    T = so.level_sets(fib, q)
    for c in range(1, q):
        shell = [a for a in T[c - 1] if a not in T[c]]
        target = K.mul(ch.minor(T[c - 1], T[c - 1]), ch.minor(T[c], T[c]))
        for a in shell:
            acc = -1
            for a2 in shell:
                acc = K.add(acc, K.mul(ch.cofactor(T[c - 1], a, a2), ch.bordered(T[c], a, a2)))
            assert acc == target, "row-sum identity fails"
            stats["row sums"] += 1


def compare(ch, fib, stats):
    K, q = ch.K, ch.q
    check_row_sums(ch, fib, stats)
    e, pre, A0, Lam = so.second_order(ch, fib)
    assert e + 1 <= q - 2
    ser = so.direct_series(ch, fib)
    assert all(ser[k] < 0 for k in range(e)), "lower coefficients do not vanish"
    assert ser[e] == K.mul(pre, A0), "leading coefficient mismatch"
    stats["configurations"] += 1
    if A0 < 0:
        assert ser[e + 1] == K.mul(pre, Lam), "second-order coefficient mismatch"
        stats["obstructed"] += 1
        stats["Lam nonzero"] += Lam >= 0


def random_fibres(M, q, rng, max_points):
    while True:
        supp = rng.sample(range(M), rng.randint(2, min(M, 6)))
        fib = {a: sorted(rng.sample(range(q), rng.choice([1, 1, 2, 2, 3]))) for a in supp}
        e = sum(comb(len(v), 2) for v in fib.values())
        if e + 1 <= q - 2 and sum(len(v) for v in fib.values()) <= max_points:
            return fib


def singular_sets(ch, max_rank):
    out = []
    for r in range(2, max_rank + 1):
        for S in itertools.combinations(range(ch.M), r):
            if ch.minor(list(S), list(S)) < 0:
                out.append(list(S))
    return out


def targeted(ch, S, rng, max_points):
    """S singular; returns fibres with S as a level set of multiplicity one, or None."""
    q, M = ch.q, ch.M
    outside = [a for a in range(M) if a not in S]
    for _ in range(200):
        fib = {}
        if rng.random() < 0.5:
            # bottom level: sizes 1 (upper shell U) or 2 on S, nothing else
            U = set(rng.sample(S, rng.randint(1, len(S))))
            for a in S:
                fib[a] = sorted(rng.sample(range(q), 1 if a in U else 2))
        else:
            # second level: sizes 2 (upper shell) or 3 on S, size 1 on a lower shell L
            U = set(rng.sample(S, rng.randint(1, len(S))))
            for a in S:
                fib[a] = sorted(rng.sample(range(q), 2 if a in U else 3))
            for a in rng.sample(outside, rng.randint(1, min(3, len(outside)))):
                fib[a] = sorted(rng.sample(range(q), 1))
        e = sum(comb(len(v), 2) for v in fib.values())
        if e + 1 <= q - 2 and sum(len(v) for v in fib.values()) <= max_points:
            return fib
    return None


def tall(ch, S, rng):
    """S singular of rank 3, placed at the third level: sizes 3 on S (upper shell),
    size 2 on one fibre outside S (lower shell).  Eleven points, e = 10."""
    q, M = ch.q, ch.M
    fib = {a: sorted(rng.sample(range(q), 3)) for a in S}
    a = rng.choice([x for x in range(M) if x not in S])
    fib[a] = sorted(rng.sample(range(q), 2))
    return fib


def main():
    rng = random.Random(20261127)
    plan = [  # (M, q, random configurations per chart, targeted per chart, largest rank searched for singular sets)
        (6, 7, 30, 0, 0),
        (10, 3, 30, 12, 6),
        (10, 11, 40, 40, 4),
        (14, 13, 20, 20, 5),
        (14, 29, 20, 30, 5),
        (15, 7, 20, 30, 6),
        (21, 5, 20, 0, 0),
    ]
    for M, q, n_random, n_targeted, max_rank in plan:
        stats = {"configurations": 0, "obstructed": 0, "Lam nonzero": 0, "row sums": 0}
        for ch in so.charts(M, q):
            for _ in range(n_random):
                compare(ch, random_fibres(M, q, rng, 9), stats)
            if n_targeted:
                sing = singular_sets(ch, max_rank)
                done = 0
                while sing and done < n_targeted:
                    fib = targeted(ch, rng.choice(sing), rng, 10)
                    if fib is None:
                        break
                    compare(ch, fib, stats)
                    done += 1
        print("(M, q) = (%d, %d): %d configurations agree; %d of them obstructed (A0 = 0), "
              "second-order term nonzero in %d; %d row-sum identities hold" % (M, q, stats["configurations"], stats["obstructed"], stats["Lam nonzero"], stats["row sums"]), flush=True)
    # higher levels: (11, 23) has singular carriers of rank 3 and q - 2 = 21 leaves room for e = 10
    M, q = 11, 23
    stats = {"configurations": 0, "obstructed": 0, "Lam nonzero": 0, "row sums": 0}
    ch = so.charts(M, q)[0]
    sing = [S for S in singular_sets(ch, 3) if len(S) == 3]
    for _ in range(6):
        compare(ch, tall(ch, rng.choice(sing), rng), stats)
    print("(M, q) = (%d, %d), singular carrier at the third level: %d configurations agree; %d of them obstructed (A0 = 0), "
          "second-order term nonzero in %d; %d row-sum identities hold" % (M, q, stats["configurations"], stats["obstructed"], stats["Lam nonzero"], stats["row sums"]), flush=True)
    print("all comparisons agree")


if __name__ == "__main__":
    main()
