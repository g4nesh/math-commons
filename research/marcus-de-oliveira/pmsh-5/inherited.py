"""Are the hard inequalities of P(c) linear consequences of PMSH(n-1)?

Input: hard_n{n}_seed{seed}.npz from facet_classify.py (spectrum c and hard normals l, scaled so that
l.(x - g) <= 1 on P(c), computed with the shifted spectrum c + 1).
For each hard inequality this reports its maximum over
  R_triv: D >= 0 on aff{D(U)}  (D = all squared minors, m = M D + const), and
  R_inh:  R_triv plus, for every position t and shift gamma in {0} + gammas,
          (a) proj_{T not containing t}(m) in proj(P(c)),
          (b) proj_{T containing t}(m(H + gamma I)) in proj(P(c + gamma)).
(a), (b) follow from PMSH(n-1) (see inherited_n4.py and README N3). If the R_inh maximum is <= 1, the
hard inequality is a linear consequence of PMSH(n-1) and nonnegativity of squared minors.

usage: inherited.py hard_n5_seed31.npz gammas|none
"""
import itertools
import sys

import numpy as np
from scipy.optimize import linprog

from facet_search import PM

data = np.load(sys.argv[1])
gammas = [] if sys.argv[2] == "none" else [float(x) for x in sys.argv[2].split(",")]
c = data["c"]
n = len(c)
cs = c + 1.0
P = PM(n)
rng = np.random.default_rng(0)
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
REL = vt[int((sv > 1e-9 * sv[0]).sum()):]
nD = len(pairs)


def shift(gm):
    L = np.zeros((d, d))
    c0 = np.zeros(d)
    for i, T in enumerate(SUBS):
        c0[i] = gm ** len(T)
        for k in range(1, len(T) + 1):
            for R in itertools.combinations(T, k):
                L[i, si[R]] += gm ** (len(T) - k)
    return L, c0


blocks = []
for t in range(n):
    sel_a = [i for i, T in enumerate(SUBS) if t not in T]
    Wa = np.unique(np.round(V[:, sel_a], 12), axis=0)  # projections of distinct permutations coincide
    blocks.append((np.eye(d)[sel_a], np.zeros(len(sel_a)), Wa))
    sel_b = [i for i, T in enumerate(SUBS) if t in T]
    for gm in [0.0] + gammas:
        L, c0 = shift(gm)
        Wb = np.unique(np.round((V @ L.T + c0)[:, sel_b], 12), axis=0)
        blocks.append((L[sel_b], c0[sel_b], Wb))
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

out = []
for l, isf in zip(data["normals"], data["facet"]):
    off = l @ (const - g)
    r1 = linprog(-(M.T @ l), A_eq=REL[:, :-1], b_eq=-REL[:, -1], bounds=(0, None), method="highs")
    obj = np.zeros(nv)
    obj[:nD] = -(M.T @ l)
    r2 = linprog(obj, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    v1 = -r1.fun + off
    v2 = -r2.fun + off if r2.status == 0 else np.nan
    out.append((bool(isf), v1, v2))
still = [o for o in out if not o[2] <= 1 + 1e-9]
print(f"{sys.argv[1]}: c={np.round(c, 6).tolist()} gammas={gammas}: hard inequalities {len(out)} "
      f"(true facets {sum(o[0] for o in out)}); max over R_triv {max((o[1] for o in out), default=float('nan')):.6f}; "
      f"still above 1 over R_inh: {len(still)} (true facets {sum(o[0] for o in still)}), "
      f"max over R_inh {max((o[2] for o in out), default=float('nan')):.6f}")
