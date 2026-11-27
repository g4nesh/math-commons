#!/usr/bin/env python3
"""Sanity check, not a proof: draw random members K of the family that P(105)
reduces to (some level set of K equals the punctured cross C) and test the
principal minor of F_105 on K exactly modulo the prime ell of lifts105_exact.py.
A nonzero value modulo ell proves that minor is nonzero in C.

Usage: python3 sample105.py . COUNT
"""
import random, sys
sys.path.insert(0, sys.argv[1])
import lifts105_exact as L
ell = 997500421
w = None
for g in range(2, 1000):
    c = pow(g, (ell - 1) // 105, ell)
    if all(pow(c, 105 // r, ell) != 1 for r in (3, 5, 7)):
        w = c; break
wp = [pow(w, e, ell) for e in range(105)]
crt = {(k % 15, k % 7): k for k in range(105)}
C = {3, 5, 6, 9, 10, 12}
rng = random.Random(105)
n = int(sys.argv[2]); zero = 0; sizes = []
for _ in range(n):
    t = rng.randint(1, 7)
    ks = []
    for a in range(15):
        c = rng.randint(t, 7) if a in C else rng.randint(0, t - 1)
        for b in rng.sample(range(7), c):
            ks.append(crt[(a, b)])
    sizes.append(len(ks))
    if L.det_mod([[wp[(k * l) % 105] for l in ks] for k in ks], ell) == 0:
        zero += 1
        print("ZERO mod ell:", sorted(ks))
print("random members of the residual family for 105: %d sampled, sizes %d..%d, vanishing modulo ell: %d" % (n, min(sizes), max(sizes), zero))
