"""Shared exact machinery for the sum-of-squares certificate  beta - Phi(O) >= 0  on O(5).

Variables x_0..x_24, x_{5i+j} = O[i][j].  A monomial is a sorted tuple of variable indices; a polynomial is a dict
{monomial: Fraction}.  Phi is the functional of phi.py written as a polynomial in the entries of a REAL matrix.

Sign characters.  A monomial has a row parity vector and a column parity vector.  Cubic monomials are grouped into 225 blocks
by this pair; a product of two monomials of the same block has trivial character (all row and column degrees even).

Symmetry group GAMMA (40 elements): x_{ij} -> x_{g(i) g(j)} for the 20 affine maps g of Z_5, with or without transposition.
Phi is invariant under GAMMA."""
import itertools
from fractions import Fraction as Fr
N = 5
def var(i, j): return 5 * (i % 5) + (j % 5)
def char(m):
    r = [0] * 5; c = [0] * 5
    for v in m: r[v // 5] ^= 1; c[v % 5] ^= 1
    return tuple(r), tuple(c)
TRIVIAL = ((0,) * 5, (0,) * 5)
AFF = [[(a * i + b) % 5 for i in range(5)] for a in (1, 2, 3, 4) for b in range(5)]
GAMMA = [[5 * g[v // 5] + g[v % 5] for v in range(25)] for g in AFF] + [[5 * g[v % 5] + g[v // 5] for v in range(25)] for g in AFF]
def act(g, m): return tuple(sorted(g[v] for v in m))
_canon = {}
def canon(m):
    """canonical representative of the GAMMA-orbit of a monomial, and the orbit size"""
    r = _canon.get(m)
    if r is None:
        imgs = {act(g, m) for g in GAMMA}; rep = min(imgs); r = (rep, len(imgs))
        for im in imgs: _canon[im] = r
    return r
def padd(p, q, c=1):
    for m, v in q.items():
        nv = p.get(m, 0) + c * v
        if nv: p[m] = nv
        else: p.pop(m, None)
    return p
def pmul(p, q):
    out = {}
    for m1, v1 in p.items():
        for m2, v2 in q.items():
            m = tuple(sorted(m1 + m2)); nv = out.get(m, 0) + v1 * v2
            if nv: out[m] = nv
            else: out.pop(m, None)
    return out
def average(p):
    """orbit coordinates of the GAMMA-average of p: {orbit representative: coefficient of each monomial of the orbit}"""
    out = {}
    for m, v in p.items():
        rep, size = canon(m); nv = out.get(rep, 0) + Fr(v) / size
        if nv: out[rep] = nv
        else: out.pop(rep, None)
    return out
def minor_sq(S, T):
    """(x_ij x_kl - x_il x_kj)^2 as a polynomial, S = (i, k), T = (j, l); or x_ij^2 for singletons"""
    if len(S) == 1: return {tuple(sorted((var(S[0], T[0]),) * 2)): Fr(1)}
    i, k = S; j, l = T; d = {(min(var(i, j), var(k, l)), max(var(i, j), var(k, l))): Fr(1)}; padd(d, {(min(var(i, l), var(k, j)), max(var(i, l), var(k, j))): Fr(1)}, -1)
    return pmul(d, d)
def phi_terms():
    out = [((v,), (v,), -1) for v in range(5)]
    for v in range(5):
        a = tuple(sorted(((v + 1) % 5, (v - 1) % 5))); b = tuple(sorted(((v + 2) % 5, (v - 2) % 5))); out += [(a, b, 1), (b, a, 1)]
        for d in (1, 2):
            s = tuple(sorted((v, (v + d) % 5))); t = tuple(sorted((v, (v - d) % 5))); out += [(s, t, -1), (t, s, -1)]
    return out
def phi_poly():
    p = {}
    for S, T, c in phi_terms(): padd(p, minor_sq(S, T), c)
    return p
def relations():
    """the 30 orthogonality relations: rows i <= k: sum_j x_ij x_kj - delta_ik ; columns j <= l: sum_i x_ij x_il - delta_jl"""
    out = []
    for i in range(5):
        for k in range(i, 5):
            r = {}
            for j in range(5): padd(r, {tuple(sorted((var(i, j), var(k, j)))): Fr(1)})
            if i == k: r[()] = Fr(-1)
            out.append(r)
    for j in range(5):
        for l in range(j, 5):
            r = {}
            for i in range(5): padd(r, {tuple(sorted((var(i, j), var(i, l)))): Fr(1)})
            if j == l: r[()] = Fr(-1)
            out.append(r)
    return out
CUBICS = list(itertools.combinations_with_replacement(range(25), 3))
def cubic_blocks():
    """dict character -> list of cubic monomials"""
    b = {}
    for m in CUBICS: b.setdefault(char(m), []).append(m)
    return b
def pair_orbits(blocks):
    """orbits under GAMMA of unordered pairs of cubic monomials in a common block.  Returns dict pair -> orbit key, where a
    pair is (m, m') with m <= m' and the key is the smallest image of the pair."""
    key = {}
    for ms in blocks.values():
        for a in range(len(ms)):
            for b in range(a, len(ms)):
                pr = (ms[a], ms[b])
                if pr in key: continue
                imgs = set()
                for g in GAMMA:
                    u, w = act(g, pr[0]), act(g, pr[1]); imgs.add((u, w) if u <= w else (w, u))
                k = min(imgs)
                for im in imgs: key[im] = k
    return key
def ldl_positive(M):
    """exact test of positive definiteness of a symmetric rational matrix by Gaussian elimination: all pivots > 0"""
    n = len(M); A = [row[:] for row in M]
    for k in range(n):
        if A[k][k] <= 0: return False
        for i in range(k + 1, n):
            if A[i][k]:
                f = A[i][k] / A[k][k]
                for j in range(k, n): A[i][j] -= f * A[k][j]
    return True
