"""Try to harvest the HARD facets of P(c) from random vertices of the linear relaxation.

D(U) = all squared minors |det U[S,T]|^2 (1 <= |S| = |T| <= n-1); m(U D_c U*) = M_c D(U) + const.
Q_triv = {D >= 0} intersected with the affine hull of {D(U) : U in U(n)}: the relaxation that only
knows nonnegativity of squared minors. A facet of P(c) is EASY iff it holds on M_c(Q_triv)
(Farkas; see classify.py); otherwise it is HARD.
Harvest: maximize random linear functionals on the D-space over Q_triv (LP) to get generic
vertices D* (functionals pulled back from the m-space almost always peak at a permutation
vertex, so they never see the excess); whenever
M_c D* + const lies outside P(c), the ray from the centroid toward it exits through a facet
violated by a point of Q_triv, i.e. a hard facet. Facets are identified by their vertex sets.
Then each hard facet's functional is maximized over U(n) from several starts.

usage: relaxation_vertices.py n seed ndirs nstarts [c1,...]
At n = 4 this harvest finds NO hard facet (logs/relaxation_vertices_n4.txt), although 96 exist:
the excess of the relaxation over P(c) is too thin for random vertices. Kept as a negative result.
"""
import itertools
import json
import sys

import numpy as np
from scipy.optimize import linprog, minimize

from facet_search import PM, blocks


def main():
    n, seed, ndirs, nstarts = map(int, sys.argv[1:5])
    P = PM(n)
    rng = np.random.default_rng(seed)
    c = np.array([float(x) for x in sys.argv[5].split(",")]) if len(sys.argv) > 5 else np.r_[0.0, np.sort(rng.random(n - 2)), 1.0]
    cs = c + 1.0
    V = P.verts(cs)
    g = V.mean(0)
    N, d = V.shape

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
    rank = int((sv > 1e-9 * sv[0]).sum())
    REL = vt[rank:]
    nD = len(pairs)

    A = np.zeros((d + 1, N + 1))
    A[:d, :N] = V.T
    A[d, :N] = 1
    bounds_ray = [(0, None)] * N + [(None, None)]

    def ray(m):
        A[:d, N] = -(m - g)
        r = linprog(np.r_[np.zeros(N), -1.0], A_eq=A, b_eq=np.r_[g, 1.0], bounds=bounds_ray, method="highs")
        y = r.eqlin.marginals
        l, h = y[:d], -y[d]
        return -r.fun, l / (h - l @ g)

    hard, keys, outside = [], set(), 0
    for i in range(ndirs):
        phi = rng.standard_normal(nD)
        r = linprog(-phi, A_eq=REL[:, :-1], b_eq=-REL[:, -1], bounds=(0, None), method="highs")
        if r.status != 0:
            continue
        x = M @ r.x + const
        depth, lh = ray(x)
        if depth >= 1 - 1e-9:
            continue
        outside += 1
        if not np.all(np.isfinite(lh)) or ((V - g) @ lh).max() > 1 + 1e-7:
            continue
        key = tuple(np.flatnonzero(np.abs((V - g) @ lh - 1) < 1e-7))
        if key not in keys:
            keys.add(key)
            hard.append((lh, len(key)))

    res = []
    for lh, nv in hard:
        f = lambda x: -(lh @ (P.mU(P.unitary(x), cs) - g))
        best, shapes = -np.inf, []
        for _ in range(nstarts):
            r = minimize(f, rng.standard_normal(2 * n * n), method="L-BFGS-B", options=dict(maxiter=4000))
            if -r.fun > 1 - 1e-7:
                shapes.append("+".join(map(str, blocks(P.unitary(r.x), 1e-3))))
            best = max(best, -r.fun)
        res.append(dict(vertices=nv, value=best, shapes=sorted(set(shapes))))
        if best > 1 + 1e-9:
            np.save(f"hard_violation_n{n}_seed{seed}_{len(res)}.npy", np.r_[c, lh])
    tally = {}
    for r in res:
        k = ",".join(r["shapes"]) or "none"
        tally[k] = tally.get(k, 0) + 1
    print(json.dumps(dict(n=n, seed=seed, c=[round(float(x), 6) for x in c], relaxation_dirs=ndirs,
                          relaxation_points_outside=outside, hard_facets=len(res),
                          vertex_counts=sorted(set(r["vertices"] for r in res)),
                          max_value=max([r["value"] for r in res], default=None),
                          min_value=min([r["value"] for r in res], default=None),
                          contact_shapes=tally)), flush=True)


if __name__ == "__main__":
    main()
