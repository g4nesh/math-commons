"""One certificate for a whole half-line of the parameter:
   f0 + tau f1 + tau^2 f2 >= 0 for all tau >= 0   via   G(tau) = [I, tau I] M0 [I; tau I] + tau M1,  M0, M1 PSD blockwise."""
import sys, time, numpy as np, cvxpy as cp, scipy.sparse as sps
from sos2 import *
def solve_param(L, f0, f1, f2, slack=True, solver='CLARABEL'):
    gs = [cp.Variable(L.nvar) for _ in range(3)]; cons = []; gamma = cp.Variable()
    for (lo, hi, n) in L.slices:
        iu = np.triu_indices(n); m = hi - lo
        M0 = cp.Variable((2 * n, 2 * n), PSD=True); M1 = cp.Variable((n, n), PSD=True)
        r, c = iu
        cons.append(gs[0][lo:hi] == M0[r, c])
        cons.append(gs[1][lo:hi] == M0[r, n + c] + M0[c, n + r] + M1[r, c])
        cons.append(gs[2][lo:hi] == M0[n + r, n + c])
    # f0 - gamma (gamma <= 0 measures failure; 0 means certified)
    cons.append(L.Ar @ gs[0] == L.P @ f0 - gamma * L.ones)
    cons.append(L.Ar @ gs[1] == L.P @ f1)
    cons.append(L.Ar @ gs[2] == L.P @ f2)
    prob = cp.Problem(cp.Maximize(gamma), cons + [gamma <= 1]); prob.solve(solver=solver)
    return prob.status, gamma.value
if __name__ == '__main__':
    d, K = int(sys.argv[1]), int(sys.argv[2]); t0 = time.time(); L = Level(d, K)
    for label, bv in [("D[13,14]", b_), ("D[13,23]", bw)]:
        # t = 1 + tau:  (1+tau) g0 + (2+tau) a + (1+tau)^2 b
        f0 = L.X @ (g0 + 2 * a_ + bv); f1 = L.X @ (g0 + a_ + 2 * bv); f2 = L.X @ bv
        st, gam = solve_param(L, f0, f1, f2)
        print("degree %d, b = %s: one certificate for ALL t >= 1: slack gamma = %+.3e (%s)  [%.0fs]" % (d, label, gam, st, time.time() - t0), flush=True)
