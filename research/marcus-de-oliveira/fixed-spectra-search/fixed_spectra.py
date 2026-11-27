"""Fixed nondegenerate spectra: maximize the support gap over U(n) one hull edge at a time.

Usage (from the repository root):
    python -I -B research/marcus-de-oliveira/fixed-spectra-search/fixed_spectra.py N SEED SPECTRA TLIM OUT.json [DELTA]

With DELTA > 0 a quadratic penalty keeps U at distance >= DELTA from generalized permutations and
block-diagonal unitaries (moc_core.dist_genperm / dist_blockdiag), as in the n = 4 lane; only optima
with both distances >= 0.9 DELTA are recorded.

For spectra index r (rng = default_rng([SEED, N, r])): random spectra (kind cycles gaussian/circle/
mixed), canonicalised, rejected unless min|a_i + b_j| >= 0.1, all eigenvalue gaps >= 0.1 and
area/diam^2 >= 0.01 (the n = 4 lane's filter, optimizer-search/lifted/n4_contacts.py). For each hull
edge with outward unit normal nu and offset c, maximize

    g_e(U) = (Re(conj(nu) d(U)) - c) / diam,   d(U) = det(D_a + U D_b U*),

which is smooth in U. U = U0 cayley(K(X)) with an exact analytic gradient; U0 is recentred after each
of 4 L-BFGS rounds. Starts per edge: two Haar points and one perturbation of each endpoint
permutation matrix. The recorded value is the robust margin m (max over all edges, optimizer-search/
lifted/lift_lib.robust_margin), positive iff d is outside the hull. Every configuration with
m > -1e-9 is re-evaluated at 50 digits with an exactly unitary Cayley rebuild (moc_hp.hp_margin).
A near-zero float m is not a certificate of anything.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OPT = os.path.join(HERE, "..", "optimizer-search")
sys.path.insert(0, os.path.join(OPT, "lib"))
sys.path.insert(0, os.path.join(OPT, "lifted"))
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from moc_core import z_all, det_d, dist_genperm, dist_blockdiag, perms
from moc_param import random_unitary, random_spectra, canonical_spectra, skew_from_real, cayley, perm_matrix
from moc_opt import _adjugate, _gradX
from lift_lib import hull_indices, robust_margin
from moc_hp import hp_margin


def hull_edges(z):
    idx = hull_indices(z)
    V = z[idx]
    k = len(V)
    diam = max(abs(V[i] - V[j]) for i in range(k) for j in range(k))
    area = 0.5 * abs(sum((V[i].conj() * V[(i + 1) % k]).imag for i in range(k)))
    E = []
    for i in range(k):
        e = V[(i + 1) % k] - V[i]
        nrm = -1j * e / abs(e)
        E.append((np.conj(nrm), float((V[i] * np.conj(nrm)).real), int(idx[i]), int(idx[(i + 1) % k])))
    return E, diam, area


def edge_fg(x, U0, a, Db, e, c, diam):
    """(-g_e, -grad) for minimization; gradient formula as in moc_opt.Problem.fg (U part)."""
    n = len(a)
    I = np.eye(n)
    K = skew_from_real(x.reshape(n, n))
    W = np.linalg.inv(I + K)
    V = (I - K) @ W
    U = U0 @ V
    M = np.diag(a) + U @ Db @ U.conj().T
    adj, d = _adjugate(M)
    val = ((e * d).real - c) / diam
    adjV = U0.conj().T @ adj @ U0
    Xm = Db @ V.conj().T @ adjV
    Ym = adjV @ V @ Db
    C = -W @ Xm @ (I + V) + (I + V.conj().T) @ Ym @ W.conj().T
    return -val, -(_gradX(e * C) / diam).ravel()


def self_check(rng):
    """Analytic gradient against central differences on a random instance (assert)."""
    n = 5
    a, b = canonical_spectra(*random_spectra(n, rng, "gaussian"))
    U0 = random_unitary(n, rng)
    e = np.exp(-1j * 0.7)
    x = 0.1 * rng.standard_normal(n * n)
    f0, g = edge_fg(x, U0, a, np.diag(b), e, 0.0, 1.0)
    h = 1e-6
    for k in rng.choice(n * n, 6, replace=False):
        dx = np.zeros(n * n); dx[k] = h
        fd = (edge_fg(x + dx, U0, a, np.diag(b), e, 0.0, 1.0)[0] - edge_fg(x - dx, U0, a, np.diag(b), e, 0.0, 1.0)[0]) / (2 * h)
        assert abs(fd - g[k]) < 1e-6 * max(1.0, abs(fd)), (k, fd, g[k])


def penalty(x, U0, delta, mu=100.0):
    n = U0.shape[0]
    U = U0 @ cayley(skew_from_real(x.reshape(n, n)))
    return mu * (max(0.0, delta - dist_genperm(U)[0]) ** 2 + max(0.0, delta - dist_blockdiag(U)[0]) ** 2)


def penalised_fg(x, U0, a, Db, e, c, diam, delta):
    """edge_fg plus the n = 4 lane's quadratic penalty keeping U delta away from generalized permutations
    and block-diagonal unitaries; the penalty gradient is by central differences, only when active."""
    f, g = edge_fg(x, U0, a, Db, e, c, diam)
    p = penalty(x, U0, delta)
    if p > 0:
        h = 1e-7
        gp = np.empty_like(x)
        for k in range(len(x)):
            dx = np.zeros_like(x); dx[k] = h
            gp[k] = (penalty(x + dx, U0, delta) - penalty(x - dx, U0, delta)) / (2 * h)
        f, g = f + p, g + gp
    return f, g


def climb(U, a, Db, e, c, diam, rounds=4, maxiter=300, delta=0.0):
    n = len(a)
    fun, args = (edge_fg, (U, a, Db, e, c, diam)) if delta <= 0 else (penalised_fg, (U, a, Db, e, c, diam, delta))
    for _ in range(rounds):
        args = (U,) + args[1:]
        res = minimize(fun, np.zeros(n * n), args=args, jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, gtol=1e-11, ftol=1e-15))
        U = U @ cayley(skew_from_real(res.x.reshape(n, n)))
    return U


def contacts(d, z, diam, tol=1e-6):
    """Permutations within tol*diam of the best supporting line through the hull point nearest d."""
    E, _, _ = hull_edges(z)
    best = max(E, key=lambda t: (t[0] * d).real - t[1])
    e, c = best[0], best[1]
    return int(np.sum(np.abs((e * z).real - c) <= tol * diam))


def rec_of(a, b, U, m, edge, start, z, diam):
    return dict(m=float(m), edge=edge, start=start, dist_genperm=float(dist_genperm(U)[0]),
                dist_blockdiag=float(dist_blockdiag(U)[0]), contacts=contacts(det_d(a, b, U), z, diam),
                a=[[v.real, v.imag] for v in a], b=[[v.real, v.imag] for v in b],
                U_re=U.real.tolist(), U_im=U.imag.tolist())


def main(n, seed, nspec, tlim, out, delta=0.0):
    t0 = time.time()
    self_check(np.random.default_rng(12345))
    P = perms(n)
    spectra, top = [], []
    r = tried = 0
    while len(spectra) < nspec and time.time() - t0 < tlim:
        rng = np.random.default_rng([seed, n, r])
        kind = ["gaussian", "circleA", "mixed"][r % 3]
        r += 1
        if kind == "circleA":
            # moc_param's 'circle' kind puts a and -b on one circle, whose hull is a segment; use a on
            # the unit circle against a Gaussian b instead
            a0 = np.exp(2j * np.pi * rng.random(n))
            b0 = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
            a, b = canonical_spectra(a0, b0)
        else:
            a, b = canonical_spectra(*random_spectra(n, rng, kind))
        sab = float(np.abs(a[:, None] + b[None, :]).min())
        sa = min(abs(a[i] - a[j]) for i in range(n) for j in range(i))
        sb = min(abs(b[i] - b[j]) for i in range(n) for j in range(i))
        z = z_all(a, b)
        E, diam, area = hull_edges(z)
        tried += 1
        if sab < 0.1 or sa < 0.1 or sb < 0.1 or area / diam ** 2 < 0.01:
            continue
        Db = np.diag(b)
        best, nonperm_best, per_edge = None, None, []
        for ei, (e, c, i0, i1) in enumerate(E):
            ham = int(np.sum(P[i0] != P[i1]))   # 2 iff the edge's endpoints differ by one transposition
            g_best = -np.inf
            starts = [("haar", random_unitary(n, rng)), ("haar", random_unitary(n, rng))]
            for s in (i0, i1):
                Kp = skew_from_real(0.3 * rng.standard_normal((n, n)))
                starts.append(("perm%d" % s, perm_matrix(P[s]) @ expm(Kp)))
            for name, U in starts:
                U = climb(U, a, Db, e, c, diam, delta=delta)
                if delta > 0 and (dist_genperm(U)[0] < 0.9 * delta or dist_blockdiag(U)[0] < 0.9 * delta):
                    continue   # constrained lane: keep only points that respect the distance limits
                d = det_d(a, b, U)
                g_best = max(g_best, float(((e * d).real - c) / diam))
                m = robust_margin(d, z)
                rec = rec_of(a, b, U, m, ei, name, z, diam)
                if best is None or m > best["m"]:
                    best = rec
                if rec["dist_genperm"] > 0.05 and (nonperm_best is None or m > nonperm_best["m"]):
                    nonperm_best = rec
            per_edge.append(dict(edge=ei, hamming=ham, best_edge_gap=g_best))
            if time.time() - t0 > tlim:
                break
        if best is None:
            continue   # constrained lane: no start stayed away from the trivial sets
        far = [t["best_edge_gap"] for t in per_edge if t["hamming"] > 2 and t["best_edge_gap"] > -np.inf]
        spectra.append(dict(index=r - 1, kind=kind, hull_vertices=len(E), min_abs_a_plus_b=sab,
                            min_gap=float(min(sa, sb)), fatness=float(area / diam ** 2),
                            best=best, best_away_from_permutations=nonperm_best, per_edge=per_edge,
                            best_gap_on_non_transposition_edges=max(far) if far else None))
        top.append(best)
        print(n, r - 1, kind, len(E), "%.3e" % best["m"], best["contacts"],
              "%.3e" % (nonperm_best["m"] if nonperm_best else float("nan")),
              "far_edges=%d" % len(far), "far_best=%.3e" % (max(far) if far else float("nan")), flush=True)
    top.sort(key=lambda t: -t["m"])
    hp = []
    for t in top[:5]:
        if t["m"] > -1e-9:
            a = np.array([complex(*v) for v in t["a"]]); b = np.array([complex(*v) for v in t["b"]])
            U = np.array(t["U_re"]) + 1j * np.array(t["U_im"])
            h = hp_margin(a, b, U, dps=50)
            hp.append(dict(float_m=t["m"], hp_m=str(h.get("m_str", h.get("m")))))
    ms = np.array([s["best"]["m"] for s in spectra]) if spectra else np.array([])
    mn = np.array([s["best_away_from_permutations"]["m"] for s in spectra if s["best_away_from_permutations"]])
    summary = dict(n=n, seed=seed, delta=delta, spectra_tried=tried, spectra_used=len(spectra), seconds=round(time.time() - t0, 1),
                   max_m=float(ms.max()) if len(ms) else None,
                   max_m_away_from_permutations=float(mn.max()) if len(mn) else None,
                   quantiles_m=np.quantile(ms, [0, .5, .9, 1]).tolist() if len(ms) else [],
                   quantiles_m_away=np.quantile(mn, [0, .5, .9, 1]).tolist() if len(mn) else [],
                   non_transposition_edges=int(sum(1 for s in spectra for t in s["per_edge"] if t["hamming"] > 2)),
                   max_gap_on_non_transposition_edges=max((s["best_gap_on_non_transposition_edges"] for s in spectra
                                                           if s["best_gap_on_non_transposition_edges"] is not None), default=None),
                   hp_rechecks=hp, numpy=np.__version__)
    with open(out, "w") as f:
        json.dump(dict(summary=summary, spectra=spectra), f, indent=1)
    print(json.dumps(summary))


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5],
         float(sys.argv[6]) if len(sys.argv) > 6 else 0.0)
