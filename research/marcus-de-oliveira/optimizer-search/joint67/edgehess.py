"""Second-order test on transposition edges of H (n = 6, 7).

If z_sigma and z_pi are adjacent hull vertices with pi = sigma composed with the
transposition of rows i, j, the curve U(t) = P_sigma with rows i, j rotated by
angle t gives d(U(t)) = cos^2 t z_sigma + sin^2 t z_pi, which runs along that
edge of H. With theta the outward normal of the edge, f(U) = Re(e^{-i theta} d(U))
- h_H(theta) vanishes on the curve. MOC forces f <= 0 on U(n), so at every point
of the curve the gradient of f must vanish and its Hessian must be negative
semidefinite. This script computes both (analytic gradient of the toolkit
objective in the Cayley chart U = U(t) cayley(K(X)), Hessian by central
differences of that gradient) for random spectra, and reports the largest
gradient norm and the largest Hessian eigenvalue, in units of diam(H).
A clearly positive eigenvalue would give a violation by stepping along its
eigenvector (checked with margin and hp_margin if it ever happens).

Usage: python edgehess.py N SAMPLES SEED OUT.json [--time T]
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import z_all, hull_indices, perms, margin, diameter  # noqa
from moc_param import random_spectra, canonical_spectra, perm_matrix, cayley, skew_from_real  # noqa
from moc_opt import Problem, pack  # noqa


def edge_curve(n, sg, i, j, t):
    U = perm_matrix(sg).astype(complex)
    c, s = math.cos(t), math.sin(t)
    U[i, sg[i]], U[i, sg[j]], U[j, sg[i]], U[j, sg[j]] = c, s, -s, c
    return U


def gradX(prob, a, b, X, th):
    n = prob.n
    _, g = prob.fg(pack(a, b, X, th))
    return -g[4 * n:4 * n + n * n]           # gradient of +F w.r.t. X (F = num / s)


def test_point(a, b, U, th, s_over_diam, h=1e-5):
    n = len(a)
    prob = Problem(n, U0=U, tau=0.01)
    X0 = np.zeros((n, n))
    g0 = gradX(prob, a, b, X0, th) * s_over_diam
    N = n * n
    H = np.zeros((N, N))
    for p in range(N):
        E = np.zeros(N)
        E[p] = h
        gp = gradX(prob, a, b, E.reshape(n, n), th)
        gm = gradX(prob, a, b, -E.reshape(n, n), th)
        H[:, p] = (gp - gm) / (2 * h)
    H = 0.5 * (H + H.T) * s_over_diam
    w, V = np.linalg.eigh(H)
    return float(np.abs(g0).max()), float(w[-1]), V[:, -1], H


def main():
    args = sys.argv[1:]
    n, S, seed, out = int(args[0]), int(args[1]), int(args[2]), args[3]
    tlim = float(args[args.index("--time") + 1]) if "--time" in args else 480.0
    P = perms(n)
    t0 = time.time()
    recs = []
    best = None
    n_edges = 0
    for r in range(S):
        if time.time() - t0 > tlim:
            break
        rng = np.random.default_rng([seed, r])
        a, b = canonical_spectra(*random_spectra(n, rng, "mixed"))
        pts = np.concatenate([a, -b])
        M = np.column_stack([np.abs(pts) ** 2, pts.real, pts.imag, np.ones(2 * n)])
        sv = np.linalg.svd(M, compute_uv=False)
        if sv[-1] / sv[0] < 1e-8:
            continue    # a and -b on one circle: collinear hull, f is identically 0 (known case)
        z = z_all(a, b)
        hv = hull_indices(z)
        diam = diameter(z[hv])
        zc = z - z.mean()
        s = float(np.sqrt(np.mean(np.abs(zc) ** 2)))
        for k in range(len(hv)):
            p1, p2 = hv[k], hv[(k + 1) % len(hv)]
            diff = np.nonzero(P[p1] != P[p2])[0]
            if len(diff) != 2:
                continue
            i, j = int(diff[0]), int(diff[1])
            e = z[p2] - z[p1]
            th = float(np.angle(-1j * e))          # outward normal for CCW order
            for t in (math.pi / 8, math.pi / 4, 3 * math.pi / 8):
                U = edge_curve(n, P[p1], i, j, t)
                gmax, lmax, v, H = test_point(a, b, U, th, s / diam)
                n_edges += 1
                rec = dict(sample=r, edge=[int(p1), int(p2)], rows=[i, j], t=t, grad_max=gmax, lam_max=lmax)
                recs.append(rec)
                if best is None or lmax > best["lam_max"]:
                    best = dict(rec, a=[[x.real, x.imag] for x in a], b=[[x.real, x.imag] for x in b],
                                theta=th, eigvec=v.tolist())
    lam = np.array([x["lam_max"] for x in recs])
    gr = np.array([x["grad_max"] for x in recs])
    summary = dict(n=n, seed=seed, samples_done=r + 1 if recs else 0, edge_points_tested=n_edges,
                   max_lam=float(lam.max()) if len(lam) else None, max_grad=float(gr.max()) if len(gr) else None,
                   lam_quantiles=[float(q) for q in np.quantile(lam, [0, 0.5, 0.9, 1])] if len(lam) else [],
                   n_lam_gt_1e6=int(np.sum(lam > 1e-6)), seconds=time.time() - t0,
                   command="python edgehess.py " + " ".join(args))
    json.dump(dict(summary=summary, best=best, points=recs), open(out, "w"))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
