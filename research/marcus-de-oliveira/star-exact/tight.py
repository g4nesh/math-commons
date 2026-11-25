"""Is the degree-3 SOS bound for f_t exactly 0, or slightly negative?  Solve  max gamma : f_t - gamma is SOS  with tight
solver tolerances, for f_t scaled so that max f_t over the samples is 1."""
import sys, time, numpy as np, cvxpy as cp, scipy.sparse as sps
from numeric import *
Us, rows, pcs = sample(); E = np.array([[float(x) for x in row] for row in rows])
u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-9 * sv[0]).sum()); P = u[:, :rk].T; PE = P @ E; ones = P @ np.ones(E.shape[0])
def bound(f, tol):
    s = cp.Variable(NW); gam = cp.Variable(); cons = []
    for idx, n in BLOCK_SLICES:
        if n == 1: cons.append(s[idx[0]] >= 0); continue
        S = sps.lil_matrix((n * n, len(idx))); k = 0
        for i in range(n):
            for j in range(i, n):
                if i == j: S[i * n + i, k] = 1
                else: S[i * n + j, k] = 0.5; S[j * n + i, k] = 0.5
                k += 1
        cons.append(cp.reshape(sps.csr_matrix(S) @ s[idx], (n, n), order='C') >> 0)
    cons.append(PE @ s == P @ f - gam * ones)
    prob = cp.Problem(cp.Maximize(gam), cons)
    prob.solve(solver='CLARABEL', tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol, max_iter=500)
    return gam.value, prob.status
for t in [float(x) for x in sys.argv[1].split(',')]:
    f = np.array([float(t * g + (1 + t) * a + t * t * b) for g, a, b in pcs]); sc = f.max(); f = f / sc
    for tol in (1e-8, 1e-10, 1e-12):
        t0 = time.time(); g, st = bound(f, tol)
        print("t = %-4g tolerance %.0e : best gamma (f scaled to max 1) = %+.4e  (%s)  [%.0fs]" % (t, tol, g, st, time.time() - t0), flush=True)
