"""Reduced SDP on the complement spaces saved by reduce.py.
mode 'lsq'    : minimise the residual of the identity subject to S_b PSD (is the reduced identity feasible at all?)
mode 'center' : identity exact, maximise lambda with S_b >= lambda I (is it strictly feasible?)"""
import sys, time, json, pickle, numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from numeric import *
def load_comp(t): return [[[Fr(q) for q in row] for row in c] for c in json.load(open('complement_t%g.json' % t))['complement']]
def build(comp, PE):
    """for each block with r > 0: matrix A_b (rk x r*r) such that the block's contribution to P E s is A_b vec(S_b)"""
    out = []
    for (idx, n), c in zip(BLOCK_SLICES, comp):
        r = len(c)
        if r == 0: out.append(None); continue
        Cf = np.array([[float(q) for q in row] for row in c]); A = np.zeros((PE.shape[0], r, r)); k = 0
        for i in range(n):
            for j in range(i, n):
                w = np.outer(Cf[:, i], Cf[:, j]); w = w + w.T if i != j else w
                A += PE[:, idx[k]][:, None, None] * w[None]; k += 1
        out.append(A.reshape(PE.shape[0], r * r))
    return out
if __name__ == '__main__':
    t = float(sys.argv[1]); mode = sys.argv[2]; t0 = time.time()
    comp = load_comp(t); Us, rows, pcs = sample(); E = np.array([[float(x) for x in row] for row in rows])
    f = np.array([float(t * g + (1 + t) * a + t * t * b) for g, a, b in pcs])
    u, sv, vt = np.linalg.svd(E, full_matrices=False); rk = int((sv > 1e-9 * sv[0]).sum()); P = u[:, :rk].T; PE = P @ E; Pf = P @ f
    As = build(comp, PE)
    # keep only independent combinations of the identity conditions on the reduced unknowns
    Afull0 = np.hstack([A for A in As if A is not None]); uu, ss, _ = np.linalg.svd(Afull0, full_matrices=False)
    Qr = uu[:, ss > 1e-9 * ss[0]].T; print("independent reduced conditions: %d ; component of f outside their range: %.2e" % (Qr.shape[0], np.linalg.norm(Pf - Qr.T @ (Qr @ Pf))))
    As = [None if A is None else Qr @ A for A in As]; Pf = Qr @ Pf
    Svars = []; cons = []; expr = 0; lam = cp.Variable()
    for A, c in zip(As, comp):
        if A is None: Svars.append(None); continue
        r = len(c); S = cp.Variable((r, r), symmetric=True); Svars.append(S)
        cons.append(S - (lam if mode == 'center' else 0) * np.eye(r) >> 0)
        expr = expr + A @ cp.vec(S, order='C')
    nun = sum(len(c) * (len(c) + 1) // 2 for c in comp)
    print("unknowns %d" % nun)
    if mode == 'lsq':
        prob = cp.Problem(cp.Minimize(cp.norm(expr - Pf, 2)), cons); prob.solve(solver='CLARABEL')
        print("t = %g: least-squares residual of the reduced identity = %.3e (%s)  [%.0fs]" % (t, prob.value, prob.status, time.time() - t0))
    else:
        prob = cp.Problem(cp.Maximize(lam), cons + [expr == Pf]); prob.solve(solver='CLARABEL')
        print("t = %g: best lambda = %s (%s)  [%.0fs]" % (t, None if lam.value is None else "%+.4e" % lam.value, prob.status, time.time() - t0))
    if Svars and Svars[[i for i, s in enumerate(Svars) if s is not None][0]].value is not None:
        vals = [None if S is None else S.value for S in Svars]
        mins = [np.linalg.eigvalsh(v).min() for v in vals if v is not None]
        print("   min eigenvalue over reduced blocks %.3e ; blocks with min eigenvalue < 1e-6: %d of %d" % (min(mins), sum(m < 1e-6 for m in mins), len(mins)))
        pickle.dump({'t': t, 'S': vals}, open('reduced_%s_t%g.pkl' % (mode, t), 'wb'))
