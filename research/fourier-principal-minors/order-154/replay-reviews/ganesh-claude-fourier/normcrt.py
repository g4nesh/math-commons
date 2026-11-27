"""Exact norms of principal minors of Fourier matrices, by CRT over split primes.

For a subset S of Z_N let D_S = det(w^(ij))_{i,j in S} with w = exp(2 pi i/N), an algebraic integer in Z[w]. Its
norm Nm(D_S) = prod_{c in (Z/N)^*} D_S(w^c) is a rational integer. For a prime q = 1 mod N and an element u of
order N in F_q, reducing modulo a prime of Z[w] above q gives Nm(D_S) = prod_c det(u^(c ij))_{S} (mod q).
Hadamard's inequality gives |D_S(w^c)| <= m^(m/2) with m = |S|, so |Nm(D_S)| <= m^(m phi(N)/2). Enough primes are
used for the product to exceed twice that bound, plus one further prime on which the result is re-checked, so the
integers returned are exact.

All arithmetic is integer: numpy int64 with entries below 2^31, so products stay below 2^62.
"""
import itertools
import math

import numpy as np

MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)  # deterministic for n < 3.3e24


def is_prime(n):
    if n < 2:
        return False
    for p in MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_factors(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    return out + ([n] if n > 1 else [])


def split_primes(N, count, below=2 ** 31):
    """The `count` largest primes q < below with q = 1 (mod N)."""
    out, q = [], below - 1 - ((below - 2) % N)
    while len(out) < count:
        if is_prime(q):
            out.append(q)
        q -= N
    return out


def root_of_unity(N, q):
    """An element of order exactly N in F_q (requires q = 1 mod N)."""
    for g in range(2, q):
        u = pow(g, (q - 1) // N, q)
        if all(pow(u, N // r, q) != 1 for r in prime_factors(N)):
            return u
    raise ValueError


def modinv(a, q):
    """Elementwise inverse mod q of a nonzero int64 array, by Fermat."""
    r = np.ones_like(a)
    e, b = q - 2, a % q
    while e:
        if e & 1:
            r = r * b % q
        b = b * b % q
        e >>= 1
    return r


def det_mod(A, q):
    """Determinants mod q of a batch A (B x m x m, entries in [0, q)), by Gaussian elimination with pivoting."""
    A = A.copy()
    B, m, _ = A.shape
    det = np.ones(B, dtype=np.int64)
    rows = np.arange(B)
    for c in range(m):
        nz = A[:, c:, c] != 0
        has = nz.any(axis=1)
        piv = np.argmax(nz, axis=1) + c
        top, other = A[rows, c, :].copy(), A[rows, piv, :].copy()
        A[rows, c, :], A[rows, piv, :] = other, top
        det = np.where(piv != c, (q - det) % q, det)
        p = A[:, c, c]
        det = np.where(has, det * p % q, 0)
        inv = modinv(np.where(has, p, 1), q)
        f = A[:, c + 1:, c] * inv[:, None] % q
        A[:, c + 1:, :] = (A[:, c + 1:, :] - f[:, :, None] * A[:, None, c, :] % q) % q
    return det


def norms(N, subsets):
    """Exact Nm(D_S) for a list of subsets of Z_N, all of the same size m >= 1."""
    m = len(subsets[0])
    units = [c for c in range(1, N) if math.gcd(c, N) == 1]
    bound = m ** (m * len(units) / 2)
    need, primes = 0, []
    while True:
        primes = split_primes(N, need)
        if math.prod(primes) > 2 * bound:
            break
        need += 1
    check = split_primes(N, need + 1)[-1]
    idx = np.array(subsets, dtype=np.int64)
    E = (idx[:, :, None] * idx[:, None, :]) % N  # exponents ij mod N
    residues = []
    for q in primes + [check]:
        u = root_of_unity(N, q)
        powers = np.array([pow(u, k, q) for k in range(N)], dtype=np.int64)
        r = np.ones(len(subsets), dtype=np.int64)
        for c in units:
            r = r * det_mod(powers[(c * E) % N], q) % q
        residues.append(r)
    M = math.prod(primes)
    out = []
    for k in range(len(subsets)):
        x = 0
        for q, res in zip(primes, residues):
            Mq = M // q
            x = (x + int(res[k]) * Mq * pow(Mq, -1, q)) % M
        if x > M // 2:
            x -= M
        assert x % check == int(residues[-1][k]), "CRT result fails the check prime"
        out.append(x)
    return out


def valuation(n, p):
    if n == 0:
        return math.inf
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def all_subsets(N, m):
    return [list(S) for S in itertools.combinations(range(N), m)]
