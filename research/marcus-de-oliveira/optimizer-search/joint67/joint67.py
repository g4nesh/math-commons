"""Joint (a, b, U, theta) maximization of the MOC margin at n = 6 and n = 7.

Usage (from this folder):
  python joint67.py N FAMILIES RESTARTS SEED OUT.json [--dperm x] [--dblock x] [--sep x]
                    [--time T] [--nopolish]
FAMILIES is a comma list drawn from: generic, ostar, dft, tensor, dsum, edge.
Restart r uses rng = default_rng([SEED, r]) and family FAMILIES[r % len].
Uses the shared toolkit in ../lib (smooth_maximize: annealed L-BFGS with
log-sum-exp over all n! permutations, vectorized; polish: active-set SLSQP).
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import margin, z_all, hull_indices, perms, dist_genperm, dist_blockdiag, min_separation  # noqa
from moc_param import (cayley, skew_from_real, random_unitary, random_spectra,  # noqa
                       canonical_spectra, perm_matrix)
from moc_opt import smooth_maximize, polish, describe  # noqa

O_STAR = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]) / 7.0


def rot(n, rng, eps):
    """Small random unitary cayley(eps * K), K from a standard normal real matrix."""
    return cayley(skew_from_real(eps * rng.standard_normal((n, n))))


def rand_perm_matrix(n, rng):
    return perm_matrix(rng.permutation(n))


def rand_diag_phase(n, rng):
    return np.diag(np.exp(2j * np.pi * rng.random(n)))


def start(n, family, rng):
    """Returns (a, b, U, info) for one restart."""
    a, b = random_spectra(n, rng, "mixed")
    info = dict(family=family)
    if family == "generic":
        U = random_unitary(n, rng)
    elif family == "ostar":
        # O* (4 x 4, the common-weights counterexample) plus a Haar block, rotated out of block form
        W = np.eye(n, dtype=complex)
        W[:4, :4] = O_STAR
        W[4:, 4:] = random_unitary(n - 4, rng)
        eps = float(rng.choice([0.0, 0.1, 0.3]))
        U = rand_perm_matrix(n, rng) @ W @ rand_perm_matrix(n, rng) @ rot(n, rng, eps)
        info["eps"] = eps
    elif family == "dft":
        # Fourier matrix: for prime n every minor is nonzero (Chebotarev); n = 6 is not prime
        F = np.exp(2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n) / math.sqrt(n)
        eps = float(rng.choice([0.0, 0.05, 0.2]))
        U = rand_diag_phase(n, rng) @ rand_perm_matrix(n, rng) @ F @ rand_diag_phase(n, rng) @ rot(n, rng, eps)
        info["eps"] = eps
    elif family == "tensor":
        T = np.kron(random_unitary(2, rng), random_unitary(3, rng))
        W = np.eye(n, dtype=complex)
        W[:6, :6] = T
        eps = float(rng.choice([0.0, 0.1, 0.3]))
        U = rand_perm_matrix(n, rng) @ W @ rot(n, rng, eps)
        info["eps"] = eps
    elif family == "dsum":
        k = int(rng.integers(2, n - 1))
        W = np.zeros((n, n), complex)
        W[:k, :k] = random_unitary(k, rng)
        W[k:, k:] = random_unitary(n - k, rng)
        eps = float(rng.choice([0.1, 0.3, 0.6]))
        U = W @ rot(n, rng, eps)
        info.update(eps=eps, k=k)
    elif family == "edge":
        # start on a hull edge whose endpoints differ by a transposition:
        # U = P_sigma with rows i, j rotated, so d(U) = c^2 z_sigma + s^2 z_pi exactly
        z = z_all(a, b)
        P = perms(n)
        hv = hull_indices(z)
        cands = []
        for t in range(len(hv)):
            s1, s2 = P[hv[t]], P[hv[(t + 1) % len(hv)]]
            diff = np.nonzero(s1 != s2)[0]
            if len(diff) == 2:
                cands.append((hv[t], int(diff[0]), int(diff[1])))
        if cands:
            k, i, j = cands[int(rng.integers(len(cands)))]
            sg = P[k]
            ang = float(rng.uniform(0.2, math.pi / 2 - 0.2))
            U = perm_matrix(sg).astype(complex)
            c, s = math.cos(ang), math.sin(ang)
            U[i, sg[i]], U[i, sg[j]], U[j, sg[i]], U[j, sg[j]] = c, s, -s, c
            eps = float(rng.choice([0.02, 0.1]))
            U = U @ rot(n, rng, eps)
            info.update(edge=[int(k), i, j], angle=ang, eps=eps, n_transp_edges=len(cands), n_hull=len(hv))
        else:
            U = random_unitary(n, rng)
            info.update(edge=None, n_hull=len(hv))
    else:
        raise ValueError(family)
    return a, b, U, info


def main():
    args = sys.argv[1:]
    n, fams, R, seed, out = int(args[0]), args[1].split(","), int(args[2]), int(args[3]), args[4]
    opts = dict(dperm=0.0, dblock=0.0, sep=0.0)
    tlim, do_polish = 540.0, True
    k = 5
    while k < len(args):
        if args[k] == "--nopolish":
            do_polish = False
            k += 1
            continue
        key, val = args[k][2:], float(args[k + 1])
        if key == "time":
            tlim = val
        else:
            opts[key] = val
        k += 2
    t0 = time.time()
    rows, best, contacts = [], [], []
    nev_total = 0
    for r in range(R):
        if time.time() - t0 > tlim:
            break
        rng = np.random.default_rng([seed, r])
        fam = fams[r % len(fams)]
        a, b, U, info = start(n, fam, rng)
        m_start = margin(a, b, U)["m"]
        try:
            a, b, U, th, sinfo = smooth_maximize(a, b, U, **opts)
            nev_total += sinfo["nev"]
            if do_polish:
                try:
                    pa, pb, pU, pth, pinfo = polish(a, b, U, th)
                    okc = (dist_genperm(pU)[0] >= 0.9 * opts["dperm"]
                           and dist_blockdiag(pU)[0] >= 0.9 * opts["dblock"]
                           and min(min_separation(pa, pb)) >= 0.9 * opts["sep"])
                    if okc and margin(pa, pb, pU)["m"] >= margin(a, b, U)["m"]:
                        a, b, U, th = pa, pb, pU, pth
                except (np.linalg.LinAlgError, ValueError):
                    pass
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as e:
            rows.append(dict(r=r, family=fam, error=str(e)))
            continue
        rec = describe(a, b, U)
        rec.update(seed=seed, restart=r, start=info, m_start=m_start)
        rows.append(dict(r=r, family=fam, m=rec["m"], m_start=m_start, dperm=rec["dist_genperm"],
                         dblock=rec["dist_blockdiag"], sep=min(rec["sep"]), n_active=rec["n_active"],
                         ok=rec["ok"]))
        best.append(rec)
        best.sort(key=lambda t: -t["m"] if t["m"] == t["m"] else 1e300)
        del best[30:]
        if (rec["m"] > -1e-9 and rec["dist_genperm"] > 0.1 and rec["dist_blockdiag"] > 0.1
                and min(rec["sep"]) > 0.05 and len(contacts) < 60):
            contacts.append(rec)
    ms = np.array([x["m"] for x in rows if "m" in x])
    fam_stats = {}
    for f in fams:
        mf = np.array([x["m"] for x in rows if x.get("family") == f and "m" in x])
        if len(mf):
            fam_stats[f] = dict(count=int(len(mf)), max_m=float(mf.max()), median_m=float(np.median(mf)),
                                n_gt_m1e9=int(np.sum(mf > -1e-9)), n_gt_1e7=int(np.sum(mf > 1e-7)))
    summary = dict(n=n, families=fams, seed=seed, restarts_done=len(rows), opts=opts, polish=do_polish,
                   smooth_evaluations=nev_total, seconds=time.time() - t0,
                   max_m=float(ms.max()) if len(ms) else None,
                   n_m_gt_1e7=int(np.sum(ms > 1e-7)), n_m_gt_m1e9=int(np.sum(ms > -1e-9)),
                   quantiles=[float(q) for q in np.quantile(ms, [0, 0.1, 0.5, 0.9, 1])] if len(ms) else [],
                   family_stats=fam_stats, n_nontrivial_near_contacts=len(contacts),
                   command="python joint67.py " + " ".join(args))
    json.dump(dict(summary=summary, rows=rows, best=best, near_contacts=contacts), open(out, "w"))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
