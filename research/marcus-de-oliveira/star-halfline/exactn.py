"""Exact helpers shared by interval_exact.py and verify_interval.py, for degree 3 or 4.
Rational unitaries (Gaussian rationals) come from ../star-exact/exact_space.py."""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd, comb
import numpy as np
sys.path.insert(0, '../star-exact')
from exact_space import rational_unitary, is_unitary, det2, CQ
def space(deg):
    MONS = list(itertools.combinations_with_replacement(range(16), deg)); blocks = {}
    for i, m in enumerate(MONS):
        r = [0] * 4; c = [0] * 4
        for v in m: r[v // 4] += 1; c[v % 4] += 1
        blocks.setdefault((tuple(r), tuple(c)), []).append(i)
    return MONS, list(blocks.values())
def F_exact(U, s):
    """F_s = (1-s) D[12,12] + s(1-s) g0 + s^2 D[13,14] at a rational unitary, exactly"""
    B = lambda p, j: U[p][j].abs2(); D = lambda S, T: det2(U, S, T).abs2()
    g0 = B(0, 0) + B(2, 2) + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
    return (1 - s) * D((0, 1), (0, 1)) + s * (1 - s) * g0 + s * s * D((0, 2), (0, 3))
def bernstein(d, k, tau): return comb(d, k) * tau**k * (1 - tau)**(d - k)
def lcm(a, b): return a * b // gcd(a, b)
def integer_form(U):
    L = 1
    for row in U:
        for z in row: L = lcm(L, lcm(z.re.denominator, z.im.denominator))
    return [int(z.re * L) for row in U for z in row], [int(z.im * L) for row in U for z in row], L
def monomials_int(X, Y, MONS):
    re = []; im = []
    for m in MONS:
        a, b = X[m[0]], Y[m[0]]
        for v in m[1:]: a, b = a * X[v] - b * Y[v], a * Y[v] + b * X[v]
        re.append(a); im.append(b)
    return re, im
def unknown_index(A0):
    return [(bi, i, j) for bi, a in enumerate(A0) for i in range(len(a)) for j in range(i, len(a))]
# ---------- modular evaluation (vectorised) ----------
def entries_mod(Us, p):
    K = len(Us); R = np.zeros((K, 16), dtype=np.int64); I = np.zeros((K, 16), dtype=np.int64)
    for k, U in enumerate(Us):
        for v in range(16):
            z = U[v // 4][v % 4]
            R[k, v] = z.re.numerator % p * pow(z.re.denominator, -1, p) % p; I[k, v] = z.im.numerator % p * pow(z.im.denominator, -1, p) % p
    return R, I
def monomials_mod(R, I, p, MONS):
    MARR = np.array(MONS); a = R[:, MARR[:, 0]].copy(); b = I[:, MARR[:, 0]].copy()
    for k in range(1, MARR.shape[1]):
        c = R[:, MARR[:, k]]; d = I[:, MARR[:, k]]
        a, b = (a * c - b * d) % p, (a * d + b * c) % p
    return a, b
def gram_columns_mod(a, b, A0, A1, s, blocks, p):
    """K x N matrix of the coefficients of the unknowns T_b[i,j] (i <= j) in sum_b q^* T q, q = (A0 + s A1) m, modulo p"""
    sm = s.numerator % p * pow(s.denominator, -1, p) % p; cols = []
    for blk, a0, a1 in zip(blocks, A0, A1):
        if not a0: continue
        Am = (np.array(a0, dtype=np.int64) + sm * (np.array(a1, dtype=np.int64) % p)) % p
        qa = (a[:, blk] @ Am.T) % p; qb = (b[:, blk] @ Am.T) % p; r = len(a0); iu = np.triu_indices(r)
        v = (qa[:, iu[0]] * qa[:, iu[1]] + qb[:, iu[0]] * qb[:, iu[1]]) % p
        v[:, iu[0] != iu[1]] = (2 * v[:, iu[0] != iu[1]]) % p
        cols.append(v)
    return np.hstack(cols)
