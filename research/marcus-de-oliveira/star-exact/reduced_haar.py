"""Reduced SDP with the identity imposed on Haar-random unitaries (better conditioned than the rational sample).
Maximise lambda with S_b >= lambda I on the complement spaces from complement_t*.json."""
import sys, time, json, pickle, numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from exact_space import MONS, BLOCK_LIST, PAIRS, PIDX
from numeric import BLOCK_SLICES
from zeros import mono, parts
from reduced_sdp import load_comp, build
def haar(rng):
    z = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / np.sqrt(2); Q, R = np.linalg.qr(z); d = np.diag(R); return Q * (d / np.abs(d))
def setup(K=2500, seed=0):
    rng = np.random.default_rng(seed); Us = [haar(rng) for _ in range(K)]
    M = np.array([mono(U) for U in Us]); p = np.array([a for a, b in PAIRS]); q = np.array([b for a, b in PAIRS])
    E = np.real(np.conj(M[:, p]) * M[:, q]); pcs = np.array([parts(U) for U in Us])
    return E, pcs
if __name__ == '__main__':
    t = float(sys.argv[1]); t0 = time.time(); comp = load_comp(t); E, pcs = setup()
    f = t * pcs[:, 0] + (1 + t) * pcs[:, 1] + t * t * pcs[:, 2]
    u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-10 * sv[0]).sum())
    print("Haar sample: rank %d ; singular values of the evaluation matrix: max %.2e, smallest kept %.2e" % (rk, sv[0], sv[rk - 1]))
    P = u[:, :rk].T; PE = P @ E; Pf = P @ f; As = build(comp, PE)
    A0 = np.hstack([A for A in As if A is not None]); uu, ss, _ = np.linalg.svd(A0, full_matrices=False); keep = ss > 1e-9 * ss[0]
    Qr = uu[:, keep].T; print("independent reduced conditions %d (singular values %.2e .. %.2e) ; f outside their range: %.1e" % (keep.sum(), ss[0], ss[keep][-1], np.linalg.norm(Pf - Qr.T @ (Qr @ Pf))))
    Svars = []; cons = []; expr = 0; lam = cp.Variable()
    for A, c in zip(As, comp):
        if A is None: Svars.append(None); continue
        r = len(c); S = cp.Variable((r, r), symmetric=True); Svars.append(S); cons.append(S - lam * np.eye(r) >> 0); expr = expr + (Qr @ A) @ cp.vec(S, order='C')
    prob = cp.Problem(cp.Maximize(lam), cons + [expr == Qr @ Pf]); prob.solve(solver='CLARABEL')
    print("t = %g: best lambda = %+.4e (%s)  [%.0fs]" % (t, lam.value, prob.status, time.time() - t0))
    vals = [None if S is None else S.value for S in Svars]
    pickle.dump({'t': t, 'S': vals, 'lam': float(lam.value)}, open('reduced_haar_t%g.pkl' % t, 'wb'))
