"""Facet-driven numerical test of PMSH(n) at one spectrum.

PMSH(n) at a spectrum c says: for every unitary U the principal-minor vector
m(U D_c U*) = (det H[T,T])_T lies in P(c) = conv{v_pi(c)}, v_pi(c)_T = prod_{t in T} c_pi(t).

A plain ascent of "l.m(U) - max_pi l.v_pi" stalls at permutation matrices, where the
gap is 0. This script avoids that:
  1. harvest facets of P(c): shoot a ray from the centroid g toward m(U); the dual of
     the ray LP is the facet the ray leaves through, scaled so that l.(x - g) <= 1 on P(c);
  2. for every distinct facet, maximize l.(m(U) - g) over U(n) from several starts.
PMSH(n) at c fails iff some facet value exceeds 1.
Conditioning: PMSH is invariant under c -> c + 1 (Proposition H2(b) of ../minors/pmsh.md), so
all computations use the shifted spectrum, whose products lie in [1, 2^n]. A harvested
"facet" is kept only if it is valid on every vertex (l.(v - g) <= 1 + 1e-7); the count of
discarded ones is reported. A candidate violation is re-checked with an independent LP
(L1 distance from m(U) to P(c)) and reported only if that distance is positive.
Each facet maximizer is classified by the smallest row maximum of |U_ij|^2 ("perm" if it
is a permutation matrix up to phases) and by its block structure (connected components
of the graph i ~ j when |U_ij| > 1e-4).

usage: facet_search.py n seed nrays nstarts [c1,c2,...]   (writes one JSON line per run)
"""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import linprog, minimize


class PM:
    def __init__(self, n):
        self.n = n
        self.SUBS = [T for k in range(1, n + 1) for T in itertools.combinations(range(n), k)]
        self.PA = np.array(list(itertools.permutations(range(n))))
        self.MASK = np.array([[t in T for t in range(n)] for T in self.SUBS])
        byk = {}
        for i, T in enumerate(self.SUBS):
            byk.setdefault(len(T), []).append((i, T))
        self.IDX = {k: (np.array([i for i, _ in v]), np.array([T for _, T in v])) for k, v in byk.items()}

    def unitary(self, x):
        n = self.n
        K = x[: n * n].reshape(n, n) + 1j * x[n * n:2 * n * n].reshape(n, n)
        Q, R = np.linalg.qr(K)
        return Q * (np.diag(R) / abs(np.diag(R)))

    def minors(self, H):
        out = np.empty(len(self.SUBS))
        for _, (ii, TT) in self.IDX.items():
            out[ii] = np.linalg.det(H[TT[:, :, None], TT[:, None, :]]).real
        return out

    def mU(self, U, c):
        return self.minors(U @ np.diag(c) @ U.conj().T)

    def verts(self, c):
        CP = np.asarray(c)[self.PA]
        return np.prod(np.where(self.MASK[None], CP[:, None, :], 1.0), axis=2)


def blocks(U, tol=1e-4):
    n = len(U)
    adj = np.abs(U) > tol
    seen, sizes = set(), []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], set()
        while stack:
            i = stack.pop()
            if i in comp:
                continue
            comp.add(i)
            for j in range(n):
                # rows i and j are linked if they share a nonzero column
                if j not in comp and np.any(adj[i] & adj[j]):
                    stack.append(j)
        seen |= comp
        sizes.append(len(comp))
    return sorted(sizes, reverse=True)


def main():
    n, seed, nrays, nstarts = map(int, sys.argv[1:5])
    P = PM(n)
    rng = np.random.default_rng(seed)
    if len(sys.argv) > 5:
        c = np.array([float(x) for x in sys.argv[5].split(",")])
    else:
        c = np.r_[0.0, np.sort(rng.random(n - 2)), 1.0]
    cs = c + 1.0  # shifted spectrum used for all computations
    V = P.verts(cs)
    g = V.mean(0)
    N, d = V.shape
    A = np.zeros((d + 1, N + 1))
    A[:d, :N] = V.T
    A[d, :N] = 1
    cost = np.r_[np.zeros(N), -1.0]
    bounds = [(0, None)] * N + [(None, None)]

    def ray(m):
        A[:d, N] = -(m - g)
        r = linprog(cost, A_eq=A, b_eq=np.r_[g, 1.0], bounds=bounds, method="highs")
        if r.status != 0:
            return np.full(d, np.nan)  # failed LP: the caller discards non-finite normals
        y = r.eqlin.marginals
        l, h = y[:d], -y[d]
        return l / (h - l @ g)

    def rand_x(kind):
        x = rng.standard_normal(2 * n * n)
        if kind == "near-perm":
            p = rng.permutation(n)
            M = np.zeros((n, n))
            M[np.arange(n), p] = 8.0
            x[: n * n] += M.ravel()
        return x

    def l1dist(m):
        Nv = len(V)
        Aq = np.zeros((d + 1, Nv + 2 * d))
        Aq[:d, :Nv] = V.T
        Aq[:d, Nv:Nv + d] = np.eye(d)
        Aq[:d, Nv + d:] = -np.eye(d)
        Aq[d, :Nv] = 1
        r = linprog(np.r_[np.zeros(Nv), np.ones(2 * d)], A_eq=Aq, b_eq=np.r_[m, 1.0], bounds=(0, None), method="highs")
        return r.fun if r.status == 0 else np.nan

    facets, keys, invalid = [], set(), 0
    for i in range(nrays):
        U = P.unitary(rand_x("near-perm" if i % 3 == 2 else "haar"))
        lh = ray(P.mU(U, cs))
        if not np.all(np.isfinite(lh)) or ((V - g) @ lh).max() > 1 + 1e-7:
            invalid += 1
            continue
        key = tuple(np.round(lh / np.abs(lh).max(), 6))
        if key not in keys:
            keys.add(key)
            facets.append(lh)

    out = []
    for j, lh in enumerate(facets):
        f = lambda x: -(lh @ (P.mU(P.unitary(x), cs) - g))
        best, bx = -np.inf, None
        for _ in range(nstarts):
            r = minimize(f, rand_x("haar"), method="L-BFGS-B", options=dict(maxiter=3000))
            if -r.fun > best:
                best, bx = -r.fun, r.x
        U = P.unitary(bx)
        B = np.abs(U) ** 2
        dist = l1dist(P.mU(U, cs)) if best > 1 + 1e-9 else 0.0
        out.append(dict(value=best, rowmax_min=float(B.max(1).min()), blocks=blocks(U), l1=dist))
        if best > 1 + 1e-9 and not dist <= 1e-9:
            np.save(f"violation_n{n}_seed{seed}_facet{j}.npy", np.r_[c, bx, lh])
    vals = np.array([o["value"] for o in out])
    nonperm = [o for o in out if o["rowmax_min"] < 1 - 1e-6 and o["value"] > 1 - 1e-7]
    shapes = {}
    for o in nonperm:
        k = "+".join(map(str, o["blocks"]))
        shapes[k] = shapes.get(k, 0) + 1
    flagged = [o for o in out if o["value"] > 1 + 1e-9]
    confirmed = [o for o in flagged if not o["l1"] <= 1e-9]
    print(json.dumps(dict(n=n, seed=seed, c=[round(float(x), 6) for x in c], rays=nrays, facets=len(facets),
                          invalid_discarded=invalid, starts=nstarts, max_value=float(vals.max()), min_value=float(vals.min()),
                          above_1_then_rechecked=len(flagged), confirmed_outside=len(confirmed),
                          touched_off_permutations=len(nonperm), touching_block_shapes=shapes)), flush=True)


if __name__ == "__main__":
    main()
