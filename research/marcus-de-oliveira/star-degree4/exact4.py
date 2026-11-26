"""Exact helpers at degree 4: rational unitaries in integer form, exact and modular evaluation of amplitudes."""
import sys, json, random
from fractions import Fraction as Fr
from math import gcd
import numpy as np
sys.path.insert(0, '../star-exact')
from exact_space import rational_unitary, is_unitary, det2, CQ
from space4 import MONS, MARR, BLOCK_LIST
def lcm(a, b): return a * b // gcd(a, b)
def integer_form(U):
    """entries u = (X + iY) / L with integers X, Y (lists of 16) and a common positive integer L"""
    L = 1
    for row in U:
        for z in row: L = lcm(L, lcm(z.re.denominator, z.im.denominator))
    X = [int(z.re * L) for row in U for z in row]; Y = [int(z.im * L) for row in U for z in row]
    return X, Y, L
def monomials_int(X, Y):
    """degree-4 monomials as Gaussian integers (numerators over L^4)"""
    re = []; im = []
    for m in MONS:
        a, b = X[m[0]], Y[m[0]]
        for v in m[1:]: a, b = a * X[v] - b * Y[v], a * Y[v] + b * X[v]
        re.append(a); im.append(b)
    return re, im
def load_amplitudes(path):
    d = json.load(open(path)); amp = []
    for a in d['amplitudes']:
        rows = [[Fr(q) for q in row] for row in a]
        assert all(q.denominator == 1 for row in rows for q in row), "amplitude rows must be integer"
        amp.append([[int(q) for q in row] for row in rows])
    return Fr(d['t']), amp
def amplitudes_int(re, im, amp):
    out = []
    for b, a in zip(BLOCK_LIST, amp):
        qs = []
        for row in a:
            x = y = 0
            for c, mi in zip(row, b):
                if c: x += c * re[mi]; y += c * im[mi]
            qs.append((x, y))
        out.append(qs)
    return out
def f_exact(U, tq):
    B = lambda s, j: U[s][j].abs2(); D = lambda S, T: det2(U, S, T).abs2()
    g0 = B(0, 0) + B(2, 2) + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
    return tq * g0 + (1 + tq) * D((0, 1), (0, 1)) + tq * tq * D((0, 2), (0, 3))
def unknown_index(amp):
    idx = []; n = 0
    for bi, a in enumerate(amp):
        r = len(a)
        for i in range(r):
            for j in range(i, r): idx.append((bi, i, j)); n += 1
    return idx
# ---------- modular evaluation (vectorised) ----------
def entries_mod(Us, p):
    """arrays (K,16) of real and imaginary parts of the entries modulo p"""
    K = len(Us); R = np.zeros((K, 16), dtype=np.int64); I = np.zeros((K, 16), dtype=np.int64)
    for k, U in enumerate(Us):
        for v in range(16):
            z = U[v // 4][v % 4]
            R[k, v] = z.re.numerator % p * pow(z.re.denominator, -1, p) % p; I[k, v] = z.im.numerator % p * pow(z.im.denominator, -1, p) % p
    return R, I
def monomials_mod(R, I, p):
    a = R[:, MARR[:, 0]].copy(); b = I[:, MARR[:, 0]].copy()
    for k in range(1, 4):
        c = R[:, MARR[:, k]]; d = I[:, MARR[:, k]]
        a, b = (a * c - b * d) % p, (a * d + b * c) % p
    return a, b
def amplitude_columns_mod(a, b, amp, p):
    """matrix (K x N) of the coefficients of the unknowns T_b[i,j] (i <= j) in sum_b q^* T q, modulo p"""
    cols = []
    for blk, A in zip(BLOCK_LIST, amp):
        if not A: continue
        Am = np.array(A, dtype=np.int64) % p
        qa = (a[:, blk] @ Am.T) % p; qb = (b[:, blk] @ Am.T) % p; r = len(A); iu = np.triu_indices(r)
        v = (qa[:, iu[0]] * qa[:, iu[1]] + qb[:, iu[0]] * qb[:, iu[1]]) % p
        v[:, iu[0] != iu[1]] = (2 * v[:, iu[0] != iu[1]]) % p
        cols.append(v)
    return np.hstack(cols)
