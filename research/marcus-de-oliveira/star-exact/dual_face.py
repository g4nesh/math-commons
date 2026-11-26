"""Forced kernels from the dual side.  A linear functional L on the 261-dimensional function space with
   L(f) = 0,  L(1) = 1,  moment blocks  M_b(L) = ( L(Re conj(m_p) m_q) )_{p,q}  PSD
forces  G_b M_b(L) = 0  for every exact certificate f = sum_b m_b^* G_b m_b.  A central (maximal-rank) such L therefore
gives the forced kernels  N_b = range M_b(L)  and the usable spaces  C_b = ker M_b(L).
L is parametrised as a combination of evaluations at the rational unitaries."""
import sys, time, pickle, numpy as np, cvxpy as cp
from numeric import *
def moment_blocks(y, E):
    v = E.T @ y if not isinstance(y, cp.Expression) else E.T @ y
    return v
if __name__ == '__main__':
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0; t0 = time.time()
    Us, rows, pcs = sample()
    E = np.array([[float(x) for x in row] for row in rows])
    f = np.array([float(t * g + (1 + t) * a + t * t * b) for g, a, b in pcs])
    u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-9 * sv[0]).sum())
    # functional L in coordinates z on the row space of E:  values on W-basis  v = Vt[:rk]^T z
    z = cp.Variable(rk); v = vt[:rk].T @ z
    # L(1) and L(f): express 1 and f in W via least squares on samples (both are in the span restricted to U(4))
    c1 = np.linalg.lstsq(E, np.ones(len(f)), rcond=None)[0]; cf = np.linalg.lstsq(E, f, rcond=None)[0]
    cons = [c1 @ v == 1, cf @ v == 0]
    for idx, n in BLOCK_SLICES:
        if n == 1: cons.append(v[idx[0]] >= 0); continue
        rowsM = [[None] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i, n):
                rowsM[i][j] = v[idx[k]]; rowsM[j][i] = v[idx[k]]; k += 1
        cons.append(cp.bmat(rowsM) >> 0)
    prob = cp.Problem(cp.Minimize(0), cons); prob.solve(solver='CLARABEL')
    vv = vt[:rk].T @ z.value
    print("t = %g: status %s  [%.0fs]" % (t, prob.status, time.time() - t0))
    eigs = []; ranks = []; sizes = []
    for idx, n in BLOCK_SLICES:
        M = np.zeros((n, n)); k = 0
        for i in range(n):
            for j in range(i, n): M[i, j] = M[j, i] = vv[idx[k]]; k += 1
        w = np.linalg.eigvalsh(M); eigs.extend(w.tolist()); sizes.append(n)
    eigs = np.sort(np.abs(np.array(eigs)))
    print("moment eigenvalues (%d): histogram by decade" % len(eigs))
    for lo, hi in [(0, 1e-11), (1e-11, 1e-9), (1e-9, 1e-7), (1e-7, 1e-5), (1e-5, 1e-3), (1e-3, 1e-1), (1e-1, 10)]:
        print("   [%g, %g): %d" % (lo, hi, ((eigs >= lo) & (eigs < hi)).sum()))
    pickle.dump({'t': t, 'v': vv}, open('dual_t%g.pkl' % t, 'wb'))
