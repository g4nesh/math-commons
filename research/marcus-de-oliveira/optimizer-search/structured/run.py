"""Driver: for each structured unitary U (family, n, instance), maximize the
margin m over spectra from many random starts, in two modes:
  free  : no penalties (drifts toward the trivial zero-margin sets)
  guard : sep = 0.1 (relative eigenvalue gap), fat = 0.2 (no collinear a or b)
          and zfat = 0.05 (no nearly collinear z cloud, i.e. no segment hull)
Each run is followed by the fixed-U active-set polish (with the same guards).
Every record carries hull_width (minimum width / diameter of H); a record is
non-degenerate when hull_width >= 1e-3. Free-mode starts run every FREE_EVERY-th restart.
Usage: python run.py OUT.json SEED TIME_LIMIT_S spec [spec ...]
  spec = family:n:instances:restarts, e.g. fourier:7:1:20
Spectra start r of instance i uses default_rng([SEED, family_id, n, i, r]).
"""
import json
import math
import sys
import time

import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fixedu import maximize_spectra, polish_fixed, hull_width_ratio  # noqa: E402
from moc_core import margin, dist_genperm, dist_blockdiag, z_all, hull_indices, hull_lean

FREE_EVERY = int(os.environ.get("MOC_FREE_EVERY", "4"))  # 0 disables free mode
GUARD = dict(sep=0.1, fat=0.2, zfat=0.05)


def width_of(a, b):
    n = len(a)
    if n <= 8:
        z = z_all(a, b)
        V = z[hull_indices(z)]
    else:
        V, _ = hull_lean(a, b)
    return hull_width_ratio(V)
from moc_param import random_spectra
import unitaries as UN

TAUS2 = (0.01, 1e-3, 1e-4)
FAMILIES = ["fourier", "circulant", "twisted", "sylvester8", "ostar_I", "ostar_F", "ostar_O",
            "ostar_cplx", "ostar_kron", "cycle_big", "cycle_small", "cycle_short", "fixed_file"]


def make_U(fam, n, rng, arg=None):
    if fam == "fourier":
        return UN.fourier(n)
    if fam == "circulant":
        return UN.circulant(n, rng)
    if fam == "twisted":
        p = {6: 2, 8: 2, 9: 3, 4: 2}[n]
        return UN.twisted_hadamard(p, n // p, rng)
    if fam == "sylvester8":
        return UN.sylvester8().astype(complex)
    if fam == "ostar_I":
        return UN.ostar_sum(n, rng, mix=0.6)
    if fam == "ostar_F":
        return UN.ostar_sum(n, rng, mix=0.6, second="fourier")
    if fam == "ostar_O":
        return UN.ostar_sum(n, rng, mix=0.6, second="ostar")
    if fam == "ostar_cplx":
        return UN.ostar_sum(n, rng, mix=0.6, complex_phase=True, second="fourier")
    if fam == "ostar_kron":
        return UN.ostar_kron(rng, t=rng.uniform(0.3, np.pi / 2 - 0.3))
    if fam == "cycle_big":
        return UN.cycle_rotations(n, rng, big=True)
    if fam == "cycle_small":
        return UN.cycle_rotations(n, rng, big=False, reps=2)
    if fam == "cycle_short":
        return UN.cycle_rotations(n, rng, length=min(n, 4), big=True, reps=2)
    raise ValueError(fam)


def fatness(v):
    c = v - v.mean()
    X = np.stack([c.real, c.imag])
    ev = np.linalg.eigvalsh(X @ X.T)
    return math.sqrt(max(ev[0], 0) / ev[1]) if ev[1] > 0 else 0.0


def rel_sep(a, b):
    n = len(a)
    rms = math.sqrt((np.sum(np.abs(a - a.mean()) ** 2) + np.sum(np.abs(b - b.mean()) ** 2)) / (2 * n))
    iu = np.triu_indices(n, 1)
    return [float(np.abs(v[:, None] - v[None, :])[iu].min() / rms) for v in (a, b)]


def cplx_list(v):
    return [[float(x.real), float(x.imag)] for x in np.asarray(v, complex).ravel()]


def main():
    out, seed, tlim = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    specs = sys.argv[4:]
    t0 = time.time()
    recs = []
    summary = []
    stop = False
    for spec in specs:
        parts = spec.split(":")
        fam, n, inst, R = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
        fid = FAMILIES.index(fam)
        for i in range(inst):
            Urng = np.random.default_rng([seed, fid, n, i, 999999])
            U = make_U(fam, n, Urng)
            ud = float(np.abs(U @ U.conj().T - np.eye(n)).max())
            best = {"free": None, "guard": None}
            cnt = {"free": 0, "guard": 0}
            for r in range(R):
                for mode in (("free", "guard") if FREE_EVERY > 0 and r % FREE_EVERY == 0 else ("guard",)):
                    if time.time() - t0 > tlim:
                        stop = True
                        break
                    rng = np.random.default_rng([seed, fid, n, i, r, mode == "guard"])
                    a, b = random_spectra(n, rng, "mixed")
                    kw = dict(GUARD) if mode == "guard" else {}
                    try:
                        a, b, th, nev = maximize_spectra(U, a, b, rng, maxiter=100, taus=TAUS2, **kw)
                        m_smooth = margin(a, b, U)["m"]
                        pa, pb, pth, pinfo = polish_fixed(U, a, b, th, rng=rng, **kw)
                        if margin(pa, pb, U)["m"] > m_smooth:
                            a, b, th = pa, pb, pth
                    except (np.linalg.LinAlgError, ValueError, FloatingPointError, ZeroDivisionError):
                        continue
                    r_ = margin(a, b, U)
                    cnt[mode] += 1
                    rec = dict(family=fam, n=n, instance=i, restart=r, mode=mode, seed=seed,
                               m=float(r_["m"]), m_smooth=float(m_smooth), hull_width=width_of(a, b),
                               active=[int(x) if not isinstance(x, tuple) else list(x) for x in r_["active"][:12]],
                               n_active=len(r_["active"]), ok=bool(r_["ok"]), n_vertices=int(r_["n_vertices"]),
                               theta=float(r_["theta"]), d=[r_["d"].real, r_["d"].imag],
                               diam=float(r_["diam"]), rel_sep=rel_sep(a, b),
                               fat=[fatness(a), fatness(b)], nev=nev,
                               offset_over_rms=float(abs(a.mean() + b.mean()) / np.sqrt((np.sum(np.abs(a - a.mean()) ** 2) + np.sum(np.abs(b - b.mean()) ** 2)) / (2 * n))),
                               a=cplx_list(a), b=cplx_list(b))
                    key = mode if (mode == "free" or rec["hull_width"] >= 1e-3) else "free"
                    if best[key] is None or rec["m"] > best[key]["m"]:
                        best[key] = rec
                    recs.append(rec)
                if stop:
                    break
            dg = dist_genperm(U)[0]
            db = dist_blockdiag(U)[0]
            summary.append(dict(family=fam, n=n, instance=i, unitarity_defect=ud,
                                dist_genperm=float(dg), dist_blockdiag=float(db),
                                restarts_free=cnt["free"], restarts_guard=cnt["guard"],
                                best_free=None if best["free"] is None else best["free"]["m"],
                                best_guard=None if best["guard"] is None else best["guard"]["m"],
                                U_re=U.real.tolist(), U_im=U.imag.tolist()))
            print(fam, n, i, "best(any/degenerate) %.3e best(nondeg guard) %.3e" % (
                best["free"]["m"] if best["free"] else float("nan"),
                best["guard"]["m"] if best["guard"] else float("nan")),
                "dperm %.3f dblock %.3f" % (dg, db), "t=%.0f" % (time.time() - t0), flush=True)
            if stop:
                break
        if stop:
            break
    recs.sort(key=lambda t: -t["m"])
    nd = [r for r in recs if r["hull_width"] >= 1e-3 and r["mode"] == "guard"]
    keep = recs[:10] + sorted(nd, key=lambda t: -t["m"])[:25]
    allm = [r["m"] for r in recs]
    ndm = [r["m"] for r in nd]
    json.dump(dict(argv=sys.argv[1:], seconds=time.time() - t0, n_runs=len(recs),
                   max_m=max(allm) if allm else None,
                   n_gt_1e7=int(sum(m > 1e-7 for m in allm)),
                   max_m_nondegenerate_guard=max(ndm) if ndm else None, n_nondegenerate_guard=len(ndm),
                   n_gt_1e7_nondegenerate=int(sum(m > 1e-7 for m in ndm)),
                   units=summary, best=keep), open(out, "w"), indent=1)
    print("done", len(recs), "runs, max m", max(allm) if allm else None,
          "max m nondeg guard", max(ndm) if ndm else None, "seconds", time.time() - t0)


if __name__ == "__main__":
    main()
