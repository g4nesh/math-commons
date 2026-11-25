"""Joint optimization of the MOC margin at n = 5 over (a, b, U).

Usage (from this folder):
  python joint5.py --family F --restarts R --seed S --out runs/x.json
                   [--sep s] [--dperm p] [--dblock q] [--dv v] [--time T] [--no-polish]

Families of starting points (restart r uses numpy default_rng([seed, r])):
  haar     random spectra (gaussian or circle, chosen per spectrum) and Haar U
  ostar    U = (O* (+) 1) with O* = the rational 4x4 orthogonal matrix of the
           facet example, multiplied by expm(eps K) with eps in [0.05, 0.8] and
           Haar-random K, plus random diagonal phases; gaussian spectra
  circle   a, b on circles with random centres and radii
  line     a on one line, b on another line (different directions)
  cluster  each spectrum in two or three tight clusters
  roots    a = r1 * 5th roots of unity, b = c + r2 * rotated 5th roots
  herm     nearly real spectra (real parts gaussian, imaginary parts 0.05-0.3)
  mixed    each spectrum drawn independently from the kinds above
The optimizer is the toolkit's smooth_maximize (annealed L-BFGS, analytic
gradient) followed by the toolkit's polish. --dv adds one more quadratic
penalty to the toolkit objective: |d - z_sigma| / s >= dv for every sigma
(s = rms spread of the z_sigma), which keeps d away from the points z_sigma
and so away from the vertex fibres. Polish is kept only if every requested
bound still holds at 0.9 of its value.
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import perms, margin, z_all, dist_genperm, dist_blockdiag, min_separation  # noqa: E402
from moc_param import (cayley, skew_from_real, random_unitary, canonical_spectra)      # noqa: E402
from moc_opt import Problem, pack, unpack, polish, describe, _hol_to_real, _gradX      # noqa: E402

N = 5
O_STAR = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7.0
KINDS = ["haar", "ostar", "circle", "line", "cluster", "roots", "herm"]


def cgauss(rng, n):
    return rng.standard_normal(n) + 1j * rng.standard_normal(n)


def spectrum(kind, rng, n=N):
    if kind in ("haar", "ostar"):
        if rng.random() < 0.5:
            return cgauss(rng, n)
        return np.exp(2j * np.pi * rng.random(n))
    if kind == "circle":
        return cgauss(rng, 1)[0] + rng.uniform(0.3, 2) * np.exp(2j * np.pi * rng.random(n))
    if kind == "line":
        return cgauss(rng, 1)[0] + np.exp(1j * np.pi * rng.random()) * rng.standard_normal(n)
    if kind == "cluster":
        k = int(rng.integers(2, 4))
        cs = cgauss(rng, k)
        return cs[rng.integers(0, k, n)] + rng.uniform(0.01, 0.2) * cgauss(rng, n)
    if kind == "roots":
        w = np.exp(2j * np.pi * np.arange(n) / n)
        return cgauss(rng, 1)[0] * 0.5 + rng.uniform(0.3, 2) * w * np.exp(2j * np.pi * rng.random())
    if kind == "herm":
        return rng.standard_normal(n) + 1j * rng.uniform(0.05, 0.3) * rng.standard_normal(n)
    raise ValueError(kind)


def start(family, rng):
    if family == "mixed":
        ka, kb = rng.choice(KINDS[2:]), rng.choice(KINDS[2:])
        a, b = spectrum(ka, rng), spectrum(kb, rng)
        U = random_unitary(N, rng)
        return a, b, U, "mixed:%s/%s" % (ka, kb)
    a, b = spectrum(family, rng), spectrum(family, rng)
    if family == "ostar":
        W = np.eye(N, dtype=complex)
        W[:4, :4] = O_STAR
        H = cgauss(rng, N * N).reshape(N, N)
        K = (H - H.conj().T) / 2
        K /= np.linalg.norm(K, 2)
        eps = rng.uniform(0.05, 0.8)
        U = W @ expm(eps * K)
        U = U * np.exp(2j * np.pi * rng.random(N))[None, :]
        p = rng.permutation(N)
        U = U[p][:, rng.permutation(N)]
        return a, b, U, "ostar:eps=%.3f" % eps
    return a, b, random_unitary(N, rng), family


class ProblemV(Problem):
    """Toolkit objective plus the vertex-distance penalty (|d - z_sigma|/s >= dv)."""

    def __init__(self, n, dv=0.0, **kw):
        super().__init__(n, **kw)
        self.dv = dv

    def fg(self, x):
        obj, grad = super().fg(x)
        if self.dv <= 0:
            return obj, grad
        n, P = self.n, self.P
        a, b, X, th = unpack(x, n)
        I = np.eye(n)
        K = skew_from_real(X)
        W = np.linalg.inv(I + K)
        V = (I - K) @ W
        U = self.U0 @ V
        Db = np.diag(b)
        M = np.diag(a) + U @ Db @ U.conj().T
        adj = np.linalg.det(M) * np.linalg.inv(M) if abs(np.linalg.det(M)) > 1e-300 else None
        if adj is None:
            return obj, grad
        d = np.linalg.det(M)
        Fm = a[None, :] + b[P]
        z = Fm.prod(1)
        zbar = z.mean()
        Nn = len(P)
        s2 = float(np.mean(np.abs(z - zbar) ** 2))
        q = np.abs(d - z) ** 2
        hit = np.nonzero(q / s2 < self.dv ** 2)[0]
        if len(hit) == 0:
            return obj, grad
        pre = np.ones((Nn, n), complex)
        suf = np.ones((Nn, n), complex)
        pre[:, 1:] = np.cumprod(Fm[:, :-1], axis=1)
        suf[:, :-1] = np.cumprod(Fm[:, ::-1][:, :-1], axis=1)[:, ::-1]
        Zk = pre * suf
        cz = np.conj(z - zbar)
        Ha = cz @ Zk
        Hb = np.zeros(n, complex)
        for j in range(n):
            Hb[j] = np.sum(cz[:, None] * Zk * (P == j))
        gs2 = np.concatenate([*_hol_to_real((2.0 / Nn) * Ha), *_hol_to_real((2.0 / Nn) * Hb)])
        ga_d = np.diag(adj)
        gb_d = np.diag(U.conj().T @ adj @ U)
        adjV = self.U0.conj().T @ adj @ self.U0
        C = -W @ (Db @ V.conj().T @ adjV) @ (I + V) + (I + V.conj().T) @ (adjV @ V @ Db) @ W.conj().T
        for k in hit:
            r = math.sqrt(q[k] / s2)
            h = self.dv - r
            obj += self.mu * h * h
            c = np.conj(d - z[k])
            gq = np.zeros(4 * n)
            ra, ia = _hol_to_real(2 * c * ga_d)
            rb, ib = _hol_to_real(2 * c * gb_d)
            ra2, ia2 = _hol_to_real(2 * c * Zk[k])
            rb2 = np.zeros(n)
            ib2 = np.zeros(n)
            rb2[P[k]], ib2[P[k]] = ra2, ia2
            gq[0:n], gq[n:2 * n] = ra - ra2, ia - ia2
            gq[2 * n:3 * n], gq[3 * n:4 * n] = rb - rb2, ib - ib2
            gX = _gradX(2 * c * C).ravel()
            # f = q / s2 ; r = sqrt(f); d r = d f / (2 r)
            dfs = gq / s2 - q[k] / s2 ** 2 * gs2
            coef = -2 * self.mu * h / (2 * max(r, 1e-300))
            grad[:4 * n] += coef * dfs
            grad[4 * n:4 * n + n * n] += coef * gX / s2
        return obj, grad


def smooth_maximize_v(a, b, U, taus=(0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4), maxiter=400, **opts):
    n = len(a)
    nev = 0
    for tau in taus:
        a, b = canonical_spectra(a, b)
        theta = margin(a, b, U)["theta"]
        prob = ProblemV(n, U0=U, tau=tau, **opts)
        x0 = pack(a, b, np.zeros((n, n)), theta)
        res = minimize(prob.fg, x0, jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, gtol=1e-12, ftol=1e-15))
        a, b, U, theta = prob.config(res.x)
        nev += prob.nev
    return a, b, U, theta, nev


def vertex_info(a, b, U):
    """Distance from d to the nearest z_sigma (over all sigma) and to the nearest
    hull vertex, both divided by diam; and the hull-edge position of d."""
    r = margin(a, b, U)
    z = z_all(a, b)
    dd = np.abs(r["d"] - z) / r["diam"]
    k = int(np.argmin(dd))
    return float(dd[k]), k


def classify(rec):
    """Trivial-contact labels for a configuration that ended with m ~ 0."""
    if not rec["m"] > -1e-9:
        return "interior"
    lab = []
    if rec["dist_genperm"] < 1e-4:
        lab.append("genperm")
    if rec["dist_blockdiag"] < 1e-4:
        lab.append("block")
    if min(rec["sep"]) < 1e-3:
        lab.append("repeated-eig")
    if rec["vertex_dist"] < 1e-5:
        lab.append("vertex")
    if rec.get("degenerate"):
        lab.append("degenerate-hull")
    return "+".join(lab) if lab else "NONTRIVIAL"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default="haar")
    ap.add_argument("--restarts", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sep", type=float, default=0.0)
    ap.add_argument("--dperm", type=float, default=0.0)
    ap.add_argument("--dblock", type=float, default=0.0)
    ap.add_argument("--dv", type=float, default=0.0)
    ap.add_argument("--time", type=float, default=540)
    ap.add_argument("--keep", type=int, default=40)
    ap.add_argument("--no-polish", action="store_true")
    A = ap.parse_args()
    opts = dict(sep=A.sep, dperm=A.dperm, dblock=A.dblock, dv=A.dv)
    t0 = time.time()
    rows, best = [], []
    nev_tot = 0
    for r in range(A.first, A.first + A.restarts):
        if time.time() - t0 > A.time:
            break
        rng = np.random.default_rng([A.seed, r])
        a, b, U, tag = start(A.family, rng)
        try:
            a, b, U, th, nev = smooth_maximize_v(a, b, U, **opts)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            continue
        nev_tot += nev
        if not A.no_polish:
            try:
                pa, pb, pU, pth, _ = polish(a, b, U, th)
                ok = (dist_genperm(pU)[0] >= 0.9 * A.dperm and dist_blockdiag(pU)[0] >= 0.9 * A.dblock
                      and min(min_separation(pa, pb)) >= 0.9 * A.sep
                      and vertex_info(pa, pb, pU)[0] * margin(pa, pb, pU)["diam"]
                      >= 0.9 * A.dv * float(np.sqrt(np.mean(np.abs(z_all(pa, pb) - z_all(pa, pb).mean()) ** 2))))
                if ok and margin(pa, pb, pU)["m"] > margin(a, b, U)["m"]:
                    a, b, U, th = pa, pb, pU, pth
            except (np.linalg.LinAlgError, ValueError):
                pass
        a, b = canonical_spectra(a, b)
        rec = describe(a, b, U)
        vd, vk = vertex_info(a, b, U)
        rec.update(seed=A.seed, restart=r, start=tag, vertex_dist=vd, nearest_perm=vk)
        rec["label"] = classify(rec)
        rows.append(dict(r=r, m=rec["m"], ok=rec["ok"], dp=rec["dist_genperm"], db=rec["dist_blockdiag"],
                         sep=min(rec["sep"]), vd=vd, nact=rec["n_active"], label=rec["label"], start=tag))
        best.append(rec)
        best.sort(key=lambda t: -t["m"] if t["m"] == t["m"] else 1e300)
        del best[A.keep:]
    ms = np.array([x["m"] for x in rows])
    labs = {}
    for x in rows:
        labs[x["label"]] = labs.get(x["label"], 0) + 1
    summary = dict(n=N, family=A.family, seed=A.seed, first=A.first, restarts_done=len(rows), opts=opts,
                   polish=not A.no_polish, max_m=float(np.nanmax(ms)) if len(ms) else None,
                   n_m_gt_1e7=int(np.sum(ms > 1e-7)), n_m_gt_1e10=int(np.sum(ms > 1e-10)),
                   quantiles=[float(v) for v in np.nanquantile(ms, [0, 0.5, 0.9, 0.99, 1])] if len(ms) else [],
                   labels_all=labs, objective_evals=nev_tot, seconds=time.time() - t0,
                   command="python joint5.py " + " ".join(sys.argv[1:]))
    with open(A.out, "w") as f:
        json.dump(dict(summary=summary, rows=rows, best=best), f, indent=1)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
