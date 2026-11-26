"""Numerical side of the exact-certificate pipeline, in the same bookkeeping as exact_space.py.
target(U)      : the function f_t = t*g0 + (1+t)*a + t^2*b evaluated exactly at a rational unitary (pieces returned).
solve_numeric  : max gamma s.t. f - gamma = sum_b m_b^* G_b m_b on U(4), G_b PSD; the identity is imposed at rational
                 unitaries converted to floats (261 independent ones suffice by Lemma 1)."""
import random, sys, time, pickle
import numpy as np, cvxpy as cp, scipy.sparse as sps
from fractions import Fraction as Fr
from exact_space import *

def pieces(U):
    """exact values of g0, a = D[12,12], b = D[13,14] at a rational unitary"""
    B = lambda s, j: U[s][j].abs2()
    D = lambda S, T: det2(U, S, T).abs2()
    g0 = B(0, 0) + B(2, 2) + D((0, 3), (2, 3)) - D((0, 1), (0, 3))
    return g0, D((0, 1), (0, 1)), D((0, 2), (0, 3))

def sample(nU=300, seed=20261125):
    rng = random.Random(seed); Us = [rational_unitary(rng) for _ in range(nU)]
    rows = [evaluation_row(U) for U in Us]; pcs = [pieces(U) for U in Us]
    return Us, rows, pcs

BLOCK_SLICES = []        # for each block: (list of pair indices in triu order, n)
for b in BLOCK_LIST:
    n = len(b); idx = [PIDX[(b[i], b[j])] for i in range(n) for j in range(i, n)]
    BLOCK_SLICES.append((idx, n))

def gram_from_s(s, n, idx):
    G = np.zeros((n, n)); k = 0
    for i in range(n):
        for j in range(i, n):
            G[i, j] = G[j, i] = s[idx[k]] if i == j else s[idx[k]] / 2; k += 1
    return G

def solve_numeric(E, f, solver='CLARABEL'):
    """E: float (nU x NW), f: float vector.  Returns gamma, s (coefficients in the e_{pq} basis)."""
    q, r, piv = None, None, None
    u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-9 * sv[0]).sum())
    P = u[:, :rk].T
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
    cons.append((P @ E) @ s == P @ f - gam * (P @ np.ones(len(f))))
    prob = cp.Problem(cp.Maximize(gam), cons); prob.solve(solver=solver)
    return gam.value, s.value, rk

if __name__ == '__main__':
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0; t0 = time.time()
    Us, rows, pcs = sample()
    E = np.array([[float(x) for x in row] for row in rows])
    f = np.array([float(t * g + (1 + t) * a + t * t * b) for g, a, b in pcs])
    gam, s, rk = solve_numeric(E, f)
    print("t = %g: constraint rank %d ; best gamma = %+.3e  [%.0fs]" % (t, rk, gam, time.time() - t0))
    eigs = []; kdim = 0; tot = 0; info = []
    for bi, (idx, n) in enumerate(BLOCK_SLICES):
        w = np.linalg.eigvalsh(gram_from_s(s, n, idx)); eigs.extend(w.tolist()); info.append(w)
    eigs = np.sort(np.array(eigs))
    print("eigenvalues of all blocks: %d total" % len(eigs))
    for lo, hi in [(-1, 1e-9), (1e-9, 1e-7), (1e-7, 1e-5), (1e-5, 1e-3), (1e-3, 1e-1), (1e-1, 1e3)]:
        print("   in [%g, %g): %d" % (lo, hi, ((eigs >= lo) & (eigs < hi)).sum()))
    pickle.dump({'t': t, 's': s, 'gamma': gam}, open('numeric_t%g.pkl' % t, 'wb'))
