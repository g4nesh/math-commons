"""Controls for the facial-reduction test.
(A) Positive control: F = sum_b |C_b m_b|^2 is a sum of squares on the complement spaces by construction, so the reduced
    SDP must report lambda > 0 for F.
(B) Mixture: F_eps = f + eps * F.  The smallest eps for which the reduced SDP becomes feasible measures how far f is from
    the reduced SOS cone.
(C) Independent zeros: the complement polynomials must vanish at freshly sampled and polished zeros."""
import sys, time, json, numpy as np, cvxpy as cp
from reduced_haar import *
from polish import polish
from zeros import U_of, f_t
from scipy.optimize import minimize
t = 1.0; comp = load_comp(t); E, pcs = setup(); K = E.shape[0]
rng = np.random.default_rng(0); Us = [haar(rng) for _ in range(K)]; M = np.array([mono(U) for U in Us])
F = np.zeros(K)
for b, c in zip(BLOCK_LIST, comp):
    if c: F += (np.abs(M[:, b] @ np.array([[float(q) for q in row] for row in c]).T)**2).sum(1)
f = t * pcs[:, 0] + (1 + t) * pcs[:, 1] + t * t * pcs[:, 2]
u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-10 * sv[0]).sum()); P = u[:, :rk].T; PE = P @ E; As = build(comp, PE)
A0 = np.hstack([A for A in As if A is not None]); uu, ss, _ = np.linalg.svd(A0, full_matrices=False); Qr = uu[:, ss > 1e-9 * ss[0]].T
def best_lambda(target):
    Pt = P @ target; out = np.linalg.norm(Pt - Qr.T @ (Qr @ Pt)); cons = []; expr = 0; lam = cp.Variable()
    for A, c in zip(As, comp):
        if A is None: continue
        r = len(c); S = cp.Variable((r, r), symmetric=True); cons.append(S - lam * np.eye(r) >> 0); expr = expr + (Qr @ A) @ cp.vec(S, order='C')
    prob = cp.Problem(cp.Maximize(lam), cons + [expr == Qr @ Pt]); prob.solve(solver='CLARABEL'); return lam.value, out
l, o = best_lambda(F); print("(A) positive control F: best lambda = %+.4e (expected > 0) ; F outside the range: %.1e" % (l, o), flush=True)
for eps in (1e-6, 1e-5, 1e-4, 1e-3):
    l, o = best_lambda(f + eps * F); print("(B) f + %.0e * F: best lambda = %+.4e" % (eps, l), flush=True)
rng2 = np.random.default_rng(777); res = []
for r in range(60):
    rr = minimize(lambda v: f_t(U_of(v), t), rng2.uniform(-np.pi, np.pi, 16), method='BFGS', options={'gtol': 1e-13})
    if rr.fun < 1e-10:
        U = polish(U_of(rr.x), t); m = mono(U)
        res.append(max(abs(np.array([[float(q) for q in row] for row in c]) @ m[b]).max() for b, c in zip(BLOCK_LIST, comp) if c))
print("(C) %d fresh polished zeros: max |c . m(U)| over all complement polynomials = %.1e" % (len(res), max(res)))
