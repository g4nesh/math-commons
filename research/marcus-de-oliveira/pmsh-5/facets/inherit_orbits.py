"""Inheritance test of ../inherited.py on one representative per orbit, with a reversal-closed shift set.

The constraints inherited from PMSH(4) (README N3 of ../README.md) over ALL shifts gamma are invariant under
the symmetry group of orbits.py: positions permute t, and the reversal H -> 3 - H (spectrum c + 1) maps
(b) at gamma to (b) at -3 - gamma. So whether a hard facet follows linearly from PMSH(4) is a property of
its orbit, and a finite shift set should be closed under gamma -> -3 - gamma. This script uses the 40
log-spaced shifts of README N3 in [-40, 40], their images -3 - gamma, and 0, -3 (82 shifts), with sparse
LP matrices (the dense ones of ../inherited.py would need about 3 GB). With GRID=union in the environment
it also adds the 15 shifts of the first N3 run (which include -1 and -1.5, inside the range of c + 1)
and their images (112 shifts).

Reads orbits.json; writes one line per representative: orbit, max over R_triv, max over R_inh
(scaled so that the facet is "<= 1").

usage: [GRID=union] inherit_orbits.py [first:last]
"""
import itertools
import json
import os
import sys

import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from facet_search import PM  # noqa: E402

n = 5
P = PM(n)
cs = np.array([1.0, 1.25, 1.5, 1.75, 2.0])
V = P.verts(cs)
g = V.mean(0)
N, d = V.shape
SUBS = P.SUBS
si = {T: i for i, T in enumerate(SUBS)}
pairs = [(S, T) for k in range(1, n) for S in itertools.combinations(range(n), k) for T in itertools.combinations(range(n), k)]
nD = len(pairs)
M = np.zeros((d, nD))
for col, (S, T) in enumerate(pairs):
    M[si[S], col] = np.prod(cs[list(T)])
const = np.zeros(d)
const[si[tuple(range(n))]] = np.prod(cs)
rng = np.random.default_rng(0)


def Dvec(U):
    return np.array([abs(np.linalg.det(U[np.ix_(S, T)])) ** 2 for S, T in pairs])


samp = np.array([np.r_[Dvec(P.unitary(rng.standard_normal(2 * n * n))), 1.0] for _ in range(2 * nD)])
_, sv, vt = np.linalg.svd(samp, full_matrices=True)
REL = vt[int((sv > 1e-9 * sv[0]).sum()):]

G40 = [float(x) for x in np.r_[-np.geomspace(40, 0.05, 20), np.geomspace(0.05, 40, 20)]]
G15 = [-10.0, -5.0, -2.5, -1.5, -1.0, -0.75, -0.5, -0.25, 0.25, 0.5, 1.0, 1.5, 2.5, 5.0, 10.0]  # README N3 grid
BASE = G40 + (G15 if os.environ.get("GRID") == "union" else [])
GAMMAS = sorted(set(round(x, 10) for x in BASE + [-3 - x for x in BASE] + [0.0, -3.0]))


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
    blocks.append((np.eye(d)[sel_a], np.zeros(len(sel_a)), np.unique(np.round(V[:, sel_a], 12), axis=0)))
    sel_b = [i for i, T in enumerate(SUBS) if t in T]
    for gm in GAMMAS:
        L, c0 = shift(gm)
        blocks.append((L[sel_b], c0[sel_b], np.unique(np.round((V @ L.T + c0)[:, sel_b], 12), axis=0)))
nv = nD + sum(W.shape[0] for _, _, W in blocks)
ri, ci, vals, rhs = [], [], [], []


def put(r, cols, v):
    nz = np.flatnonzero(v)
    ri.extend([r] * len(nz))
    ci.extend((cols[0] + nz).tolist())
    vals.extend(v[nz].tolist())


r = 0
for row in REL:
    put(r, (0,), row[:-1])
    rhs.append(-row[-1])
    r += 1
base = nD
for A, off, W in blocks:
    k = W.shape[0]
    for j in range(A.shape[0]):
        sc = max(np.abs(W[:, j]).max(), 1e-12)
        put(r, (0,), A[j] @ M / sc)
        put(r, (base,), -W[:, j] / sc)
        rhs.append(-(A[j] @ const + off[j]) / sc)
        r += 1
    put(r, (base,), np.ones(k))
    rhs.append(1.0)
    r += 1
    base += k
Aeq = sps.csr_matrix((vals, (ri, ci)), shape=(r, nv))
beq = np.array(rhs)


def main():
    orbits = json.load(open(os.path.join(HERE, "orbits.json")))["orbits"]
    lo, hi = (int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "0:1000").split(":"))
    scale = np.array([4.0 ** len(T) for T in SUBS])  # orbits.json uses the spectrum 4 (c + 1)
    for o in orbits[lo:hi]:
        w = np.array(o["l"], float) * scale
        lh = w / (o["h"] - w @ g)
        off = lh @ (const - g)
        r1 = linprog(-(M.T @ lh), A_eq=REL[:, :-1], b_eq=-REL[:, -1], bounds=(0, None), method="highs")
        obj = np.zeros(nv)
        obj[:nD] = -(M.T @ lh)
        r2 = linprog(obj, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
        v2 = -r2.fun + off if r2.status == 0 else np.nan
        print(o["orbit"], f"{-r1.fun + off:.9f}", f"{v2:.9f}", flush=True)


if __name__ == "__main__":
    main()
