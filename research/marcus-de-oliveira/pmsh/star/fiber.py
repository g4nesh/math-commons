"""Cross-check: PMSH certificate LP (weights on S_4)  vs  fiber LP (move N_U along eta q^T inside the compound polytope Q).
Projection used: pi_p(x) = (B, N p), p = eigen-matching sums of the spectrum c; q is orthogonal to p and to (1,1,1)."""
import sys, numpy as np
from scipy.optimize import linprog
from base import *
EX = list(exotic_facets().values())
def p_of(c): return np.array([c[T[0]] * c[T[1]] + c[Tc[0]] * c[Tc[1]] for (T, Tc) in MATCH])
def q_of(p):
    q = np.cross(p, np.ones(3)); return q / np.abs(q).max()
def proj(x, p):
    Nm, _ = n_matrix(x); return np.concatenate([x[:16], Nm @ p])
def cert_lp(xU, p):
    """min t s.t. |sum lam_pi proj(x_pi) - proj(xU)|_inf <= t"""
    V = np.array([proj(v.astype(float), p) for v in VX]); t = proj(xU, p); d = V.shape[1]
    A = np.vstack([np.hstack([V.T, -np.ones((d, 1))]), np.hstack([-V.T, -np.ones((d, 1))])])
    r = linprog(np.r_[np.zeros(24), 1.0], A_ub=A, b_ub=np.r_[t, -t], A_eq=np.r_[np.ones(24), 0][None], b_eq=[1.0], bounds=[(0, None)] * 25, method='highs')
    return r.fun
def fiber_lp(xU, q):
    """max s s.t. x = xU + shift(eta, q) satisfies D >= s and g_F >= s for all exotic F (B >= 0 unchanged)."""
    cols = [shift(np.eye(3)[k] - np.eye(3)[2], q) for k in range(2)]      # eta = (e1, e2, -e1-e2)
    rows, rhs = [], []
    for k in range(16, NX):                                                  # D_k >= s
        rows.append([-cols[0][k], -cols[1][k], 1.0]); rhs.append(xU[k])
    for w in EX:
        rows.append([-(w @ cols[0]), -(w @ cols[1]), 1.0]); rhs.append(w @ xU)
    r = linprog([0, 0, -1.0], A_ub=np.array(rows), b_ub=np.array(rhs), bounds=[(None, None)] * 3, method='highs')
    return -r.fun
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    mism = 0; ncap = 0; tot = 0; worst_margin = np.inf
    for trial in range(4000):
        U = haar(rng); xU = x_U(U)
        gmin = min(w @ xU for w in EX)
        if gmin >= 0 and trial % 40: continue                               # mostly test cap unitaries
        tot += 1; ncap += gmin < 0
        c = np.sort(rng.standard_normal(4) * 2); p = p_of(c); q = q_of(p)
        t = cert_lp(xU, p); s = fiber_lp(xU, q)
        ok_cert = t < 1e-9; ok_fib = s > -1e-9
        mism += ok_cert != ok_fib; worst_margin = min(worst_margin, s)
    print("tested %d unitaries (%d in exotic caps): certificate-LP vs fiber-LP disagreements: %d; min fiber margin %.3e" % (tot, ncap, mism, worst_margin))
