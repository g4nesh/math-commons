"""Build a rational sum-of-squares certificate for  c - Phi(x) >= 0  on O(5).

Phi is the functional of task mdo-orthogonal-hull-five (found by the agent
ganesh-claude-minors).  The certificate is an identity of polynomials in the
25 entries x_ij of a real 5x5 matrix:

    c - Phi(x) = sum_b v_b(x)^T G_b v_b(x) + sum_k h_k(x) g_k(x) + R(x)

  * v_b: all monomials of degree <= 3 with a fixed parity pattern of row and
    column degrees (block b), G_b a rational positive definite matrix;
  * g_k: the 30 orthogonality relations (rows and columns), h_k rational
    polynomials of degree <= 4;
  * R: whatever is left after rounding.  It is computed exactly by the checker,
    and on O(5) every monomial has absolute value <= 1, so Phi <= c + |R|_1.

Usage:  python generate_certificate.py [t0]
Writes certificate.json next to this file.  Needs numpy, scipy, cvxpy, clarabel.
"""
import itertools, json, os, sys, time
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import cvxpy as cp

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-5   # eigenvalue margin
DEN = 10 ** 9                                              # rounding denominator
start = time.time()
def log(*a):
    print(f"[{time.time() - start:6.1f}s]", *a, flush=True)

N = 5
def var(i, j):
    return N * i + j

# ---- polynomials: dict {sorted tuple of variable indices: coefficient} ----
def mul(p, q):
    r = {}
    for a, ca in p.items():
        for b, cb in q.items():
            m = tuple(sorted(a + b))
            r[m] = r.get(m, 0) + ca * cb
    return {m: c for m, c in r.items() if c != 0}

def add(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + s * c
    return {m: c for m, c in r.items() if c != 0}

def minor2(I, J):
    (a, b), (c, d) = sorted(I), sorted(J)
    return {tuple(sorted((var(a, c), var(b, d)))): 1, tuple(sorted((var(a, d), var(b, c)))): -1}

def phi_poly():
    P = {}
    for v in range(N):
        P = add(P, {(var(v, v), var(v, v)): -1})
        A = {(v + 1) % N, (v - 1) % N}
        B = {(v + 2) % N, (v - 2) % N}
        for I, J in ((A, B), (B, A)):
            m = minor2(I, J)
            P = add(P, mul(m, m))
        for d in (1, 2):
            for I, J in (({v, (v + d) % N}, {v, (v - d) % N}), ({v, (v - d) % N}, {v, (v + d) % N})):
                m = minor2(I, J)
                P = add(P, mul(m, m), -1)
    return P

def parity(m):
    r = 0
    for x in m:
        i, j = divmod(x, N)
        r ^= (1 << i) | (1 << (N + j))
    return r

def monomials(maxdeg):
    out = []
    for d in range(maxdeg + 1):
        out.extend(itertools.combinations_with_replacement(range(N * N), d))
    return out

# ---- relations ----
relations = []   # (label, poly)
for i in range(N):
    for k in range(i, N):
        g = {}
        for j in range(N):
            g = add(g, {tuple(sorted((var(i, j), var(k, j)))): 1})
        if i == k:
            g = add(g, {(): -1})
        relations.append((["row", i, k], g))
for j in range(N):
    for l in range(j, N):
        g = {}
        for i in range(N):
            g = add(g, {tuple(sorted((var(i, j), var(i, l)))): 1})
        if j == l:
            g = add(g, {(): -1})
        relations.append((["col", j, l], g))

PHI = phi_poly()
log("Phi terms", len(PHI))

# ---- blocks ----
blocks = {}
for m in monomials(3):
    blocks.setdefault(parity(m), []).append(m)
blocks = [blocks[k] for k in sorted(blocks)]
log("blocks", len(blocks), "sizes", sorted({len(b) for b in blocks}))

by_par4 = {}
for m in monomials(4):
    by_par4.setdefault(parity(m), []).append(m)

# ---- linear map: unknowns -> coefficients ----
mon_index = {}
def idx(m):
    if m not in mon_index:
        mon_index[m] = len(mon_index)
    return mon_index[m]

rows, cols, vals = [], [], []
ncol = 0
block_offsets = []
for b in blocks:
    n = len(b)
    block_offsets.append(ncol)
    for a in range(n):
        for c in range(n):
            rows.append(idx(tuple(sorted(b[a] + b[c])))); cols.append(ncol + a * n + c); vals.append(1.0)
    ncol += n * n
nG = ncol
mult_cols = []   # (relation number, monomial)
for k, (lab, g) in enumerate(relations):
    pg = parity(next(iter(g)))
    for mu in by_par4.get(pg, []):
        for m, c in g.items():
            rows.append(idx(tuple(sorted(mu + m)))); cols.append(ncol); vals.append(float(c))
        mult_cols.append((k, mu))
        ncol += 1
nL = ncol - nG
for m in PHI:
    idx(m)
idx(())
nrow = len(mon_index)
A = sp.csr_matrix((vals, (rows, cols)), shape=(nrow, ncol))
rhs = np.zeros(nrow)
for m, c in PHI.items():
    rhs[mon_index[m]] -= c          # c - Phi: the constant c is a variable
log("unknowns: Gram entries", nG, "multipliers", nL, "coefficient equations", nrow)

# ---- SDP: minimise c ----
Gs = [cp.Variable((len(b), len(b)), symmetric=True) for b in blocks]
lam = cp.Variable(nL)
c = cp.Variable()
xg = cp.hstack([cp.vec(G, order="C") for G in Gs])
e0 = np.zeros(nrow); e0[mon_index[()]] = 1.0
cons = [A[:, :nG] @ xg + A[:, nG:] @ lam - c * e0 == rhs]
cons += [G - T0 * np.eye(G.shape[0]) >> 0 for G in Gs]
prob = cp.Problem(cp.Minimize(c), cons)
prob.solve(solver=cp.CLARABEL, verbose=False)
log("SDP status", prob.status, "c* =", c.value)

# ---- round ----
def rnd(x):
    return int(round(x * DEN))          # numerator over the common denominator DEN

cstar = Fraction(int(np.ceil(c.value * 10 ** 6)), 10 ** 6)
cert_blocks = []
for b, G in zip(blocks, Gs):
    M = 0.5 * (G.value + G.value.T)
    n = len(b)
    gram = [[rnd(M[a, cc]) for cc in range(a, n)] for a in range(n)]   # upper triangle rows
    cert_blocks.append({"monomials": [[list(divmod(x, N)) for x in m] for m in b], "gram_upper": gram})
mults = {}
lv = lam.value
for (k, mu), val in zip(mult_cols, lv):
    q = rnd(val)
    if q != 0:
        mults.setdefault(k, []).append([[list(divmod(x, N)) for x in mu], q])
cert = {
    "statement": "c - Phi(x) = sum_b v_b^T G_b v_b + sum_k h_k g_k + R(x), with R computed by the checker; on O(5), Phi <= c + sum |coefficients of R|",
    "functional": "Phi from task mdo-orthogonal-hull-five (agent ganesh-claude-minors)",
    "variables": "x[i][j] = entry (i, j) of a real 5x5 matrix, i, j in 0..4; a monomial is a list of [i, j] pairs with repetition",
    "relations": "row i k: sum_j x[i][j] x[k][j] - delta(i,k);  col j l: sum_i x[i][j] x[i][l] - delta(j,l)",
    "c": str(cstar),
    "denominator": DEN,
    "note": "every Gram entry and multiplier coefficient is an integer numerator over 'denominator'",
    "eigenvalue_margin_requested": T0,
    "blocks": cert_blocks,
    "multipliers": [{"relation": relations[k][0], "terms": t} for k, t in sorted(mults.items())],
}
with open(os.path.join(HERE, "certificate.json"), "w") as f:
    json.dump(cert, f, separators=(",", ":"))
log("wrote certificate.json; c =", cstar, "; multiplier terms", sum(len(t) for t in mults.values()))
