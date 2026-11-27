"""Exact integer checks for the substitutions in README.md (standard library only).

Run: python3 arithmetic_checks.py
"""
from math import gcd


def order(a, n):
    k, x = 1, a % n
    while x != 1:
        x, k = x * a % n, k + 1
    return k


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def mul_mod_x7(f, g):
    """Multiply two integer polynomials modulo x^7 - 1."""
    h = [0] * 7
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            h[(i + j) % 7] += a * b
    return h


print("154 = prime * cofactor:", [(p, 154 // p) for p in (2, 7, 11)])
for cofactor, p in ((14, 11), (22, 7), (77, 2)):
    f, ph = order(p, cofactor), phi(cofactor)
    print(f"cofactor {cofactor}, prime {p}: ord = {f}, phi = {ph}, "
          f"primes above {p}: {ph // f}, inert: {f == ph}")

# Caragea-Lee-Malikiosis-Pfander Theorem 1.5 at (p1, p2, p3) = (2, 7, 11).
print("Thm 1.5, j = 1: the bound (2/2)^(2*phi(2)/4) is 1, and 7 > 1 holds")
bound = (14 // 2) ** (14 * phi(14) // 4)
print("Thm 1.5, j = 2: 11 >", bound, "is", 11 > bound)

# Norm of the principal minor {0, 1, 3} of F_7 (their Proposition 4.5).
# det [[1,1,1],[1,x,x^3],[1,x^3,x^9]] = 3x^3 - x - x^2 - x^6 with x^7 = 1.
D = [0, -1, -1, 3, 0, 0, -1]
prod = [1, 0, 0, 0, 0, 0, 0]
for k in range(1, 7):
    conj = [0] * 7
    for e, c in enumerate(D):
        conj[e * k % 7] += c
    prod = mul_mod_x7(prod, conj)
# A Galois-invariant element is c0 + c1 (x + ... + x^6) = c0 - c1.
assert len(set(prod[1:])) == 1
norm = prod[0] - prod[1]
print("norm of the {0,1,3} minor of F_7:", norm, "= 2^3 * 7^3:", norm == 2**3 * 7**3)
print("so the 2-adic valuation of its norm from Q(zeta_77) is", 3 * (phi(77) // phi(7)),
      "and beta(77, 2) >=", 2 * 3 * (phi(77) // phi(7)), "= phi(77) =", phi(77))
