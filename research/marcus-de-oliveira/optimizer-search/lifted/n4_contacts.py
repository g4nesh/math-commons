"""Search for NONTRIVIAL n = 4 contacts with the spectra held FIXED and non-degenerate.
Usage: python -I -B n4_contacts.py SEED R DELTA
For restart r (rng = default_rng([SEED, r])): random spectra (kind cycles gaussian/circle/mixed),
canonicalized, rejected unless min|a_i + b_j| >= 0.1, min |a_i - a_j|, |b_i - b_j| >= 0.1 and the
hull is fat (area / diam^2 >= 0.01). Then maximize m(U) over U(4) (U = U_s expm(K(X)), L-BFGS-B with
numerical gradient, 3 Haar starts) with quadratic penalties keeping dist_genperm and dist_blockdiag
>= DELTA. The hull is the robust one (lift_lib.hull_indices). Writes out/n4fixed_sSEED.json."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from lift_lib import hull_indices
from moc_core import z_all, det_d, dist_genperm, dist_blockdiag
from moc_param import random_unitary, random_spectra, canonical_spectra, skew_from_real

T0 = time.time(); TLIM = 540.0


def edges(z):
    idx = hull_indices(z)
    V = z[idx]; k = len(V)
    diam = max(abs(V[i] - V[j]) for i in range(k) for j in range(k))
    N, C = [], []
    for i in range(k):
        e = V[(i + 1) % k] - V[i]
        nrm = -1j * e / abs(e)
        N.append(np.conj(nrm)); C.append((V[i] * np.conj(nrm)).real)
    area = 0.5 * abs(sum((V[i].conj() * V[(i + 1) % k]).imag for i in range(k)))
    return np.array(N), np.array(C), diam, area, idx


def main(seed, R, delta):
    out, allbest = [], []
    for r in range(R):
        if time.time() - T0 > TLIM:
            break
        rng = np.random.default_rng([seed, r])
        kind = ["gaussian", "circle", "mixed"][r % 3]
        a, b = canonical_spectra(*random_spectra(4, rng, kind))
        sab = np.abs(a[:, None] + b[None, :]).min()
        sa = min(abs(a[i] - a[j]) for i in range(4) for j in range(i))
        sb = min(abs(b[i] - b[j]) for i in range(4) for j in range(i))
        z = z_all(a, b)
        N, C, diam, area, idx = edges(z)
        if sab < 0.1 or sa < 0.1 or sb < 0.1 or area / diam ** 2 < 0.01:
            continue
        best = None
        for s in range(3):
            Us = random_unitary(4, rng)

            def U_of(x):
                return Us @ expm(skew_from_real(x.reshape(4, 4)))

            def f(x):
                U = U_of(x)
                d = det_d(a, b, U)
                m = np.max((N * d).real - C) / diam
                pen = max(0.0, delta - dist_genperm(U)[0]) ** 2 + max(0.0, delta - dist_blockdiag(U)[0]) ** 2
                return -m + 100.0 * pen
            res = minimize(f, np.zeros(16), method="L-BFGS-B", options=dict(maxiter=3000, maxfun=60000))
            U = U_of(res.x)
            d = det_d(a, b, U)
            m = float(np.max((N * d).real - C) / diam)
            dg, db = dist_genperm(U)[0], dist_blockdiag(U)[0]
            feas = dg >= 0.9 * delta and db >= 0.9 * delta
            rec = dict(restart=r, start=s, kind=kind, m=m, dist_genperm=dg, dist_blockdiag=db, feasible=bool(feas),
                       nfev=int(res.nfev), a=[[v.real, v.imag] for v in a], b=[[v.real, v.imag] for v in b],
                       U_re=U.real.tolist(), U_im=U.imag.tolist(), min_abs_a_plus_b=float(sab), fatness=float(area / diam ** 2))
            if feas and (best is None or m > best["m"]):
                best = rec
        if best:
            allbest.append(best)
            print(r, kind, "%.3e" % best["m"], "%.3f %.3f" % (best["dist_genperm"], best["dist_blockdiag"]), flush=True)
    allbest.sort(key=lambda t: -t["m"])
    ms = np.array([t["m"] for t in allbest])
    summ = dict(seed=seed, delta=delta, spectra_used=len(allbest), max_m=float(ms.max()) if len(ms) else None,
                quantiles=np.quantile(ms, [0, .5, .9, 1]).tolist() if len(ms) else [], seconds=time.time() - T0)
    json.dump(dict(summary=summ, best=allbest[:30]), open(os.path.join(HERE, "out", "n4fixed_s%d.json" % seed), "w"), indent=1)
    print(json.dumps(summ))


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]))
