"""Relative interior of the reduced feasible set at a given t: solve the reduced identity with zero objective (the
interior-point solution approximates the analytic centre, i.e. a maximal-rank solution) and list the eigenvalues of the
reduced blocks.  Eigenvalues at noise level are further forced kernel directions."""
import sys, time, pickle, numpy as np, cvxpy as cp
from reduced_haar import *
t = float(sys.argv[1]); comp = load_comp(1.0); E, pcs = setup()
u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-10 * sv[0]).sum()); P = u[:, :rk].T; PE = P @ E; As = build(comp, PE)
A0 = np.hstack([A for A in As if A is not None]); uu, ss, _ = np.linalg.svd(A0, full_matrices=False); Qr = uu[:, ss > 1e-9 * ss[0]].T
f = t * pcs[:, 0] + (1 + t) * pcs[:, 1] + t * t * pcs[:, 2]; Pf = P @ f
cons = []; expr = 0; Svars = []; tr = 0
for A, c in zip(As, comp):
    if A is None: Svars.append(None); continue
    r = len(c); S = cp.Variable((r, r), symmetric=True); Svars.append(S); cons.append(S >> 0); expr = expr + (Qr @ A) @ cp.vec(S, order='C'); tr = tr + cp.trace(S)
prob = cp.Problem(cp.Minimize(0), cons + [expr == Qr @ Pf]); prob.solve(solver='CLARABEL')
print("t = %g: status %s" % (t, prob.status))
eig = np.concatenate([np.linalg.eigvalsh(S.value) for S in Svars if S is not None]); lg = np.log10(np.maximum(np.abs(eig), 1e-300))
print("eigenvalues of the reduced blocks (%d): min %.2e" % (len(eig), eig.min()))
for lo in range(-14, 2):
    c = ((lg >= lo) & (lg < lo + 1)).sum()
    if c: print("   1e%d .. 1e%d : %d" % (lo, lo + 1, c))
pickle.dump({'t': t, 'S': [None if S is None else S.value for S in Svars]}, open('center_t%g.pkl' % t, 'wb'))
