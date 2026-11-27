"""Hyp(14, 11) by exact norms: for every nonempty S in Z_14, 11 does not divide Nm(D_S), D_S = det(w^(ij))_{S},
w a primitive 14th root of unity. Since Nm is the product over all Galois conjugates, this says D_S is nonzero
modulo every prime of Z[w] above 11, i.e. for every primitive 14th root of unity in F_1331. Empty S: D = 1.

usage: hyp14.py
"""
import math
import random
from fractions import Fraction

import numpy as np

from normcrt import all_subsets, det_mod, norms, split_primes


def exact_det(M):
    M = [[Fraction(x) for x in row] for row in M]
    n, d = len(M), Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return d


def self_tests():
    rng = random.Random(1)
    q = split_primes(14, 1)[0]
    ok = True
    for m in range(1, 8):
        for _ in range(30):
            M = [[rng.randrange(q) if rng.random() > 0.3 else 0 for _ in range(m)] for _ in range(m)]
            if rng.random() < 0.2 and m > 1:
                M[-1] = list(M[0])  # singular
            ok &= int(det_mod(np.array([M], dtype=np.int64), q)[0]) == int(exact_det(M)) % q
    print("det mod q agrees with exact rational determinant on 210 random matrices (incl. singular):", ok)
    checks = [
        ("N=4, S={0,2} (non-square-free, vanishes)", norms(4, [[0, 2]])[0] == 0),
        ("N=9, S={0,3,6} (vanishes)", norms(9, [[0, 3, 6]])[0] == 0),
        ("N=7, S={0,1,3}: |Nm| = 2^3 7^3 (CLMP Prop 4.5)", abs(norms(7, [[0, 1, 3]])[0]) == 2 ** 3 * 7 ** 3),
        ("N=6, S={0,1,3}: |Nm| = 12 (CLMP Prop 5.4)", abs(norms(6, [[0, 1, 3]])[0]) == 12),
        ("N=14, singletons: Nm = 1", all(x == 1 for x in norms(14, [[i] for i in range(14)]))),
        ("N=14, full set: |Nm| = 14^(14*6/2)", abs(norms(14, [list(range(14))])[0]) == 14 ** 42),
    ]
    for name, good in checks:
        print(("PASS " if good else "FAIL ") + name)
        ok &= good
    return ok


def main():
    ok = self_tests()
    print("\nsize  subsets  Nm=0  11|Nm   smallest |Nm|   largest |Nm|")
    total = zero = div11 = 0
    for m in range(1, 15):
        nm = norms(14, all_subsets(14, m))
        z = sum(1 for x in nm if x == 0)
        d = sum(1 for x in nm if x != 0 and x % 11 == 0)
        total += len(nm)
        zero += z
        div11 += d
        print(f"{m:4d} {len(nm):8d} {z:5d} {d:6d}   {min(abs(x) for x in nm):>14d}  {max(abs(x) for x in nm):.3e}", flush=True)
    print(f"\nnonempty subsets {total} (+ the empty set, D = 1); Nm = 0: {zero}; 11 | Nm: {div11}")
    ok &= zero == 0 and div11 == 0 and total == 2 ** 14 - 1
    print("Hyp(14, 11) HOLDS at every prime above 11" if ok else "CHECK FAILED")


if __name__ == "__main__":
    main()
