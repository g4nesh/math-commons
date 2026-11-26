"""Can PMSH(n-1) imply PMSH(n) through linear constraints?  Test at n = 4 against all qhull facets.

D = all squared minors |det U[S,T]|^2 (1 <= |S| = |T| <= 3); m(U D_c U*) = M D + const.
Relaxation R_triv: D >= 0 on the affine hull of {D(U) : U in U(4)}.
Relaxation R_inh: R_triv plus, for every position t and every shift gamma in {0} + extra,
  (a) proj_{T not containing t}(m)            lies in proj(P(c))
  (b) proj_{T containing t}(L_gamma m)        lies in proj(L_gamma P(c)),  L_gamma m := m(H + gamma I).
Every constraint in (a), (b) follows from PMSH(3): (a) by the interlacing coupling (Lemma C of
../hermitian-lifting/README.md) applied to the compression that deletes t; (b) by Jacobi's identity
det H[T,T] = det H * det(H^-1)[T^c,T^c] and (a) for (H + gamma I)^-1.
The script reports, for the facets of P(c) exceeded by R_triv ("hard" facets), how many are still
exceeded by R_inh, and the largest excess (in qhull's unit-normal scaling).

--control replaces the vertex set of every block (a) by 12 of the 24 permutations; the excess must
then drop, which shows the constraints are really imposed.

usage: inherited_n4.py c1,c2,c3,c4 gammas|none [--control]
"""
import itertools
import sys

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

from facet_search import PM

n = 4
P = PM(n)
rng = np.random.default_rng(0)
cs = np.array([float(x) for x in sys.argv[1].split(",")])
gammas = [] if sys.argv[2] == "none" else [float(x) for x in sys.argv[2].split(",")]
control = "--control" in sys.argv
V = P.verts(cs)
g = V.mean(0)
N, d = V.shape
SUBS = P.SUBS
pairs = [(S, T) for k in range(1, n) for S in itertools.combinations(range(n), k) for T in itertools.combinations(range(n), k)]
si = {T: i for i, T in enumerate(SUBS)}
M = np.zeros((d, len(pairs)))
for col, (S, T) in enumerate(pairs):
    M[si[S], col] = np.prod(cs[list(T)])
const = np.zeros(d)
const[si[tuple(range(n))]] = np.prod(cs)


def Dvec(U):
    return np.array([abs(np.linalg.det(U[np.ix_(S, T)])) ** 2 for S, T in pairs])


samp = np.array([np.r_[Dvec(P.unitary(rng.standard_normal(2 * n * n))), 1.0] for _ in range(2 * len(pairs))])
_, sv, vt = np.linalg.svd(samp, full_matrices=True)
REL = vt[int((sv > 1e-9 * sv[0]).sum()):]  # rows r: r[:-1].D + r[-1] = 0 for every U
nD = len(pairs)


def shift(gm):
    """Matrix and offset of m -> m(H + gm I) on nonempty subsets (m_empty = 1)."""
    L = np.zeros((d, d))
    c0 = np.zeros(d)
    for i, T in enumerate(SUBS):
        c0[i] = gm ** len(T)
        for k in range(1, len(T) + 1):
            for R in itertools.combinations(T, k):
                L[i, si[R]] += gm ** (len(T) - k)
    return L, c0


blocks = []  # (A, off, W): constraint A (M D + const) + off = W^T lambda, lambda a probability vector
for t in range(n):
    sel_a = [i for i, T in enumerate(SUBS) if t not in T]
    blocks.append((np.eye(d)[sel_a], np.zeros(len(sel_a)), (V[:12] if control else V)[:, sel_a]))
    sel_b = [i for i, T in enumerate(SUBS) if t in T]
    for gm in [0.0] + gammas:
        L, c0 = shift(gm)
        blocks.append((L[sel_b], c0[sel_b], (V @ L.T + c0)[:, sel_b]))
nv = nD + sum(W.shape[0] for _, _, W in blocks)
rows, rhs = [], []
for r in REL:
    row = np.zeros(nv)
    row[:nD] = r[:-1]
    rows.append(row)
    rhs.append(-r[-1])
base = nD
for A, off, W in blocks:
    k = W.shape[0]
    for j in range(A.shape[0]):
        row = np.zeros(nv)
        row[:nD] = A[j] @ M
        row[base:base + k] = -W[:, j]
        rows.append(row)
        rhs.append(-(A[j] @ const + off[j]))
    row = np.zeros(nv)
    row[base:base + k] = 1
    rows.append(row)
    rhs.append(1.0)
    base += k
Aeq, beq = np.array(rows), np.array(rhs)

X = V - g
_, s, vt2 = np.linalg.svd(X, full_matrices=False)
Bas = vt2[: int((s > 1e-9 * s[0]).sum())]
hull = ConvexHull(X @ Bas.T)
eqr = np.round(hull.equations / np.abs(hull.equations[:, :-1]).max(1, keepdims=True), 8)
_, idx = np.unique(eqr, axis=0, return_index=True)
n_hard = n_still = 0
worst_triv = worst_inh = 0.0
for f in hull.equations[idx]:
    w = f[:-1] @ Bas  # facet: w.(x - g) + f[-1] <= 0 on P(c)
    off = w @ (const - g) + f[-1]
    r1 = linprog(-(M.T @ w), A_eq=REL[:, :-1], b_eq=-REL[:, -1], bounds=(0, None), method="highs")
    v1 = -r1.fun + off
    if v1 <= 1e-9:
        continue
    n_hard += 1
    worst_triv = max(worst_triv, v1)
    obj = np.zeros(nv)
    obj[:nD] = -(M.T @ w)
    r2 = linprog(obj, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    v2 = -r2.fun + off if r2.status == 0 else np.nan
    if not v2 <= 1e-9:
        n_still += 1
    worst_inh = max(worst_inh, v2)
print(f"c={cs.tolist()} gammas={gammas} control={control}: facets {len(idx)}; exceeded by R_triv {n_hard} "
      f"(max excess {worst_triv:.4e}); still exceeded by R_inh {n_still} (max excess {worst_inh:.4e})")
