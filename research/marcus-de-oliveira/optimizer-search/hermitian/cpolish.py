"""Constrained polish (spectra fixed): maximize t s.t. Re(e^{-i th}(d(U) - z_s))/s0 >= t for all s,
dist_genperm(U) >= L, dist_blockdiag(U) >= L. U = U0 cayley(K(X)). SLSQP, finite-difference Jacobians.
Usage: python cpolish.py runs/cont_<tag>.json idx L out.json"""
import json, os, sys
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import z_all, det_d, dist_genperm, dist_blockdiag, min_separation
from moc_param import cayley, skew_from_real
from robust import robust_margin, robust_margin_hp
from moc_hp import hp_margin
import mpmath as mp
f, idx, L, outp = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
r = json.load(open(f))["records"][idx]
a = np.array([complex(*v) for v in r["a"]]); b = np.array([complex(*v) for v in r["b"]])
U0 = np.array(r["U_re"]) + 1j * np.array(r["U_im"]); n = len(a)
z = z_all(a, b); s0 = float(np.sqrt(np.mean(np.abs(z - z.mean()) ** 2)))
rm = robust_margin(a, b, U0); th0 = rm["theta"]
def U_of(w): return U0 @ cayley(skew_from_real(w[:n * n].reshape(n, n)))
def gaps(w):
    U = U_of(w); d = det_d(a, b, U); e = np.exp(-1j * w[n * n])
    return (e * (d - z)).real / s0 - w[-1]
proj = (np.exp(-1j * th0) * z).real
work = np.argsort(-proj)[:12]
cons = [dict(type="ineq", fun=lambda w: gaps(w)[work]),
        dict(type="ineq", fun=lambda w: np.array([dist_genperm(U_of(w))[0] - L, dist_blockdiag(U_of(w))[0] - L]))]
w0 = np.concatenate([np.zeros(n * n), [th0], [rm["m"] * rm["diam"] / s0]])
for rnd in range(4):
    res = minimize(lambda w: -w[-1], w0, jac=lambda w: -np.eye(len(w))[-1], method="SLSQP",
                   constraints=cons, options=dict(maxiter=500, ftol=1e-16))
    g = gaps(res.x); new = set(np.nonzero(g < -1e-13)[0]) - set(work)
    w0 = res.x
    if not new: break
    work = np.array(sorted(set(work) | new))
U = U_of(res.x); rr = robust_margin(a, b, U)
h = hp_margin(a, b, U, dps=60); lb = robust_margin_hp(a, b, U, [rr["theta"]], dps=60)
dz = float(np.abs(det_d(a, b, U) - z).min() / rr["diam"])
proj = (np.exp(-1j * rr["theta"]) * z).real; act = np.nonzero(proj >= proj.max() - 1e-8 * rr["diam"])[0]
out = dict(source=os.path.basename(f), idx=idx, L=L, float_m_robust=rr["m"], hp_m=h["m_str"], hp_lower_bound=mp.nstr(lb, 20),
           dist_genperm=dist_genperm(U)[0], dist_blockdiag=dist_blockdiag(U)[0], sep=min_separation(a, b),
           d_to_nearest_z=dz, active=act.tolist(), a=[[v.real, v.imag] for v in a], b=[[v.real, v.imag] for v in b],
           U_re=U.real.tolist(), U_im=U.imag.tolist(), theta=rr["theta"], absU2=(np.abs(U) ** 2).tolist())
json.dump(out, open(outp, "w"), indent=1)
print(json.dumps({k: out[k] for k in ("source", "idx", "float_m_robust", "hp_m", "hp_lower_bound", "dist_genperm", "dist_blockdiag", "sep", "d_to_nearest_z", "active")}))
print(np.round(np.abs(U) ** 2, 3))
