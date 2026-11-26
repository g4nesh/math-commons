"""Harvest facets of P(c) by ray shooting, keep the HARD ones, and maximize those over U(n).

Harvest as in facet_search.py (shifted spectrum c + 1, validity check, facets keyed by their vertex
sets), with half of the rays aimed at slightly perturbed permuted block-diagonal unitaries: at n = 4
the hard facets are touched at block-reducible U, and plain Haar rays rarely reach them.
Rays that leave through a lower-dimensional face return a supporting hyperplane that is not a facet;
these are kept (any valid inequality is a test of PMSH) and reported separately from true facets
(vertex set of affine dimension dim P(c) - 1). Classification: with D = all squared minors |det U[S,T]|^2 (1 <= |S| = |T| <= n-1) and
m(U D_c U*) = M D + const, a facet l.(x - g) <= 1 is EASY iff its maximum over the linear relaxation
Q_triv = {D >= 0} on aff{D(U)} is <= 1 (then it follows from nonnegativity of squared minors), and
HARD otherwise. Only hard facets can be violated; each is maximized over U(n) from nstarts starts.
PMSH(n) at c fails iff such a maximum exceeds 1 (re-checked by an independent LP).
Contact shapes: block structures (support components, threshold 1e-3) of maximizers with value
above 1 - 1e-7.

The hard inequalities are saved to hard_n{n}_seed{seed}.npz (for inherited.py). nstarts = 0 skips
the maximization over U(n).

usage: facet_classify.py n seed nrays nstarts [c1,c2,...]
"""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import linprog, minimize

from facet_search import PM, blocks


def main():
    n, seed, nrays, nstarts = map(int, sys.argv[1:5])
    P = PM(n)
    rng = np.random.default_rng(seed)
    if len(sys.argv) > 5:
        c = np.array([float(x) for x in sys.argv[5].split(",")])
    else:
        c = np.r_[0.0, np.sort(rng.random(n - 2)), 1.0]
    cs = c + 1.0
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

    pairs = [(S, T) for k in range(1, n) for S in itertools.combinations(range(n), k) for T in itertools.combinations(range(n), k)]
    si = {T: i for i, T in enumerate(P.SUBS)}
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

    def relax_value(l):
        r = linprog(-(M.T @ l), A_eq=REL[:, :-1], b_eq=-REL[:, -1], bounds=(0, None), method="highs")
        return -r.fun + l @ (const - g) if r.status == 0 else np.nan

    def l1dist(m):
        Aq = np.zeros((d + 1, N + 2 * d))
        Aq[:d, :N] = V.T
        Aq[:d, N:N + d] = np.eye(d)
        Aq[:d, N + d:] = -np.eye(d)
        Aq[d, :N] = 1
        r = linprog(np.r_[np.zeros(N), np.ones(2 * d)], A_eq=Aq, b_eq=np.r_[m, 1.0], bounds=(0, None), method="highs")
        return r.fun if r.status == 0 else np.nan

    def rand_x(kind):
        x = rng.standard_normal(2 * n * n)
        if kind == "near-perm":
            p = rng.permutation(n)
            Mp = np.zeros((n, n))
            Mp[np.arange(n), p] = 8.0
            x[: n * n] += Mp.ravel()
        return x

    def block_unitary():
        """Permuted block-diagonal unitary with random Haar blocks, times a small random rotation."""
        sizes, left = [], n
        while left:
            k = int(rng.integers(1, left + 1)) if len(sizes) else int(rng.integers(1, n))
            sizes.append(k)
            left -= k
        B = np.zeros((n, n), complex)
        o = 0
        for k in sizes:
            Q, R = np.linalg.qr(rng.standard_normal((k, k)) + 1j * rng.standard_normal((k, k)))
            B[o:o + k, o:o + k] = Q * (np.diag(R) / abs(np.diag(R)))
            o += k
        B = B[rng.permutation(n)][:, rng.permutation(n)]
        K = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        K = (K - K.conj().T) * (10.0 ** rng.uniform(-4, -1))
        w, Q = np.linalg.eigh(-1j * K)
        return B @ (Q @ np.diag(np.exp(1j * w)) @ Q.conj().T)

    dimP = np.linalg.matrix_rank(V - g, tol=1e-9)

    def is_facet(key):
        W = V[list(key)]
        return np.linalg.matrix_rank(W[1:] - W[0], tol=1e-9) == dimP - 1

    facets, keys = [], set()
    for i in range(nrays):
        U = block_unitary() if i % 2 else P.unitary(rand_x("near-perm" if i % 4 == 1 else "haar"))
        lh = ray(P.mU(U, cs))
        if not np.all(np.isfinite(lh)) or ((V - g) @ lh).max() > 1 + 1e-7:
            continue
        key = tuple(np.flatnonzero(np.abs((V - g) @ lh - 1) < 1e-7))
        if key not in keys:
            keys.add(key)
            facets.append((lh, len(key), is_facet(key)))

    hard = []
    for lh, nv, isf in facets:
        rv = relax_value(lh)
        if rv > 1 + 1e-9:
            hard.append((lh, nv, rv, isf))
    np.savez(f"hard_n{n}_seed{seed}.npz", c=c, normals=np.array([h[0] for h in hard]).reshape(len(hard), d),
             vertices=np.array([h[1] for h in hard]), relaxation=np.array([h[2] for h in hard]),
             facet=np.array([h[3] for h in hard]))
    res = []
    for lh, nv, rv, isf in (hard if nstarts > 0 else []):
        f = lambda x: -(lh @ (P.mU(P.unitary(x), cs) - g))
        best, bx, shapes = -np.inf, None, []
        for _ in range(nstarts):
            r = minimize(f, rand_x("haar"), method="L-BFGS-B", options=dict(maxiter=4000))
            if -r.fun > 1 - 1e-7:
                shapes.append("+".join(map(str, blocks(P.unitary(r.x), 1e-3))))
            if -r.fun > best:
                best, bx = -r.fun, r.x
        dist = l1dist(P.mU(P.unitary(bx), cs)) if best > 1 + 1e-9 else 0.0
        if best > 1 + 1e-9 and not dist <= 1e-9:
            np.save(f"violation_n{n}_seed{seed}_hard{len(res)}.npy", np.r_[c, bx, lh])
        res.append(dict(vertices=nv, facet=isf, relaxation=rv, value=best, shapes=sorted(set(shapes)), l1=dist))
    tally = {}
    for r in res:
        k = ",".join(r["shapes"]) or "not reached"
        tally[k] = tally.get(k, 0) + 1
    print(json.dumps(dict(n=n, seed=seed, c=[round(float(x), 6) for x in c], rays=nrays, dim=int(dimP),
                          supporting=len(facets), facets=sum(1 for f in facets if f[2]),
                          hard=len(hard), hard_facets=sum(1 for h in hard if h[3]), hard_vertex_counts=sorted(set(h[1] for h in hard)),
                          max_relaxation_excess=max([h[2] - 1 for h in hard], default=None),
                          starts=nstarts, max_over_U=max([r["value"] for r in res], default=None),
                          min_over_U=min([r["value"] for r in res], default=None),
                          confirmed_outside=sum(1 for r in res if r["value"] > 1 + 1e-9 and not r["l1"] <= 1e-9),
                          contact_shapes=tally)), flush=True)


if __name__ == "__main__":
    main()
