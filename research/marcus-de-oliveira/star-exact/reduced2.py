"""Reduced SDP on the amplitude spaces of amplitudes_t<t>.json (after the second reduction), Haar sample.
Maximise lambda with T_b >= lambda I.  lambda > 0 means a strictly feasible point exists, which can be rounded."""
import sys, time, json, pickle, numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from reduced_haar import setup, build
from numeric import BLOCK_SLICES
t = float(sys.argv[1]); t0 = time.time()
amp = [[[Fr(q) for q in row] for row in a] for a in json.load(open('amplitudes_t%g.json' % t))['amplitudes']]
E, pcs = setup(); f = t * pcs[:, 0] + (1 + t) * pcs[:, 1] + t * t * pcs[:, 2]
u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-10 * sv[0]).sum()); P = u[:, :rk].T; PE = P @ E; Pf = P @ f; As = build(amp, PE)
A0 = np.hstack([A for A in As if A is not None]); uu, ss, _ = np.linalg.svd(A0, full_matrices=False); Qr = uu[:, ss > 1e-9 * ss[0]].T
print("unknowns %d ; independent conditions %d ; relative part of f outside the range %.1e" % (sum(len(a) * (len(a) + 1) // 2 for a in amp), Qr.shape[0], np.linalg.norm(Pf - Qr.T @ (Qr @ Pf)) / np.linalg.norm(Pf)))
cons = []; expr = 0; lam = cp.Variable(); Tv = []
for A, a in zip(As, amp):
    if A is None: Tv.append(None); continue
    r = len(a); T = cp.Variable((r, r), symmetric=True); Tv.append(T); cons.append(T - lam * np.eye(r) >> 0); expr = expr + (Qr @ A) @ cp.vec(T, order='C')
prob = cp.Problem(cp.Maximize(lam), cons + [expr == Qr @ Pf]); prob.solve(solver='CLARABEL')
print("t = %g: best lambda = %+.4e (%s)  [%.0fs]" % (t, lam.value, prob.status, time.time() - t0))
pickle.dump({'t': t, 'T': [None if T is None else T.value for T in Tv], 'lam': float(lam.value)}, open('reduced2_t%g.pkl' % t, 'wb'))
