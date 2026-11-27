#!/usr/bin/env python3
"""Exact check for the smallest slice of the family that P(105) reduces to.

Write Z_105 = Z_15 x Z_7 and C = {3, 5, 6, 9, 10, 12}, the punctured cross.  For
every b in (Z_7)^6 the set K_b = {(a, b_a) : a in C} has C as its only level set,
so both first-order criteria are silent about it.  This script shows that the
principal minor of F_105 on K_b is nonzero for all 7^6 = 117649 choices of b.

Method: ell is a prime = 1 mod 105 and w an element of order 105 in F_ell.  The
map zeta_105 -> w is a ring homomorphism Z[zeta_105] -> F_ell, so a minor that
is nonzero in F_ell is nonzero in C.

Usage: python3 lifts105_exact.py
Standard library only; about five seconds.
"""
import itertools


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def det_mod(A, ell):
    n, d = len(A), 1
    A = [r[:] for r in A]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d = d * A[c][c] % ell
        inv = pow(A[c][c], ell - 2, ell)
        for r in range(c + 1, n):
            if A[r][c]:
                fac = A[r][c] * inv % ell
                A[r] = [(x - fac * y) % ell for x, y in zip(A[r], A[c])]
    return d % ell


def main():
    ell = next(n for n in range(105 * 9500000 + 1, 10 ** 10, 105) if is_prime(n))
    w = None
    for g in range(2, 1000):
        c = pow(g, (ell - 1) // 105, ell)
        if all(pow(c, 105 // r, ell) != 1 for r in (3, 5, 7)):
            w = c
            break
    assert w is not None and pow(w, 105, ell) == 1
    wp = [pow(w, e, ell) for e in range(105)]
    cross = [3, 5, 6, 9, 10, 12]
    crt = {(k % 15, k % 7): k for k in range(105)}
    zero = 0
    for b in itertools.product(range(7), repeat=6):
        ks = [crt[(a, bb)] for a, bb in zip(cross, b)]
        if det_mod([[wp[(k * l) % 105] for l in ks] for k in ks], ell) == 0:
            zero += 1
    print("ell = %d, lifts checked: %d, minors vanishing modulo ell: %d" % (ell, 7 ** 6, zero))
    assert zero == 0


if __name__ == "__main__":
    main()
