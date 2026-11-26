"""Facial-reduction test at degree 3 for several t, on the complement spaces of the common zero set
{g0 = D[12,12] = D[13,14] = 0} (complement_t1.json).  lambda < 0 means: no exact degree-3 certificate at that t."""
import sys, time, numpy as np, cvxpy as cp
from reduced_haar import *
comp = load_comp(1.0); E, pcs = setup()
u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-10 * sv[0]).sum()); P = u[:, :rk].T; PE = P @ E; As = build(comp, PE)
A0 = np.hstack([A for A in As if A is not None]); uu, ss, _ = np.linalg.svd(A0, full_matrices=False); Qr = uu[:, ss > 1e-9 * ss[0]].T
for t in [float(x) for x in sys.argv[1].split(',')]:
    f = t * pcs[:, 0] + (1 + t) * pcs[:, 1] + t * t * pcs[:, 2]; Pf = P @ f; out = np.linalg.norm(Pf - Qr.T @ (Qr @ Pf)) / np.linalg.norm(Pf)
    cons = []; expr = 0; lam = cp.Variable()
    for A, c in zip(As, comp):
        if A is None: continue
        r = len(c); S = cp.Variable((r, r), symmetric=True); cons.append(S - lam * np.eye(r) >> 0); expr = expr + (Qr @ A) @ cp.vec(S, order='C')
    prob = cp.Problem(cp.Maximize(lam), cons + [expr == Qr @ Pf]); prob.solve(solver='CLARABEL')
    print("t = %-5g relative part of f_t outside the reduced range %.1e ; best lambda = %+.4e (%s)" % (t, out, lam.value, prob.status), flush=True)
