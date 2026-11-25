"""Reduced semidefinite problem at degree 4.
Any exact degree-4 certificate of  f_t = t*g0 + (1+t)*D[12,12] + t^2*D[13,14]  has the form
    f_t = sum_b (C_b m_b)^* S_b (C_b m_b),   S_b PSD,
with C_b the quartics of complement4.json (those vanishing on the zero set).  The identity is imposed on a Haar sample.
We maximise lambda subject to S_b >= lambda I.   lambda < 0: no exact certificate.  lambda > 0: strictly feasible."""
import sys, json, time, pickle, numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from space4 import *
def load_comp(path='complement4.json'):
    return [np.array([[float(Fr(q)) for q in row] for row in c]) if c else None for c in json.load(open(path))['complement']]
class Reduced:
    def __init__(self, comp, K=3000, seed=0):
        self.comp = comp; rng = np.random.default_rng(seed); self.Us = haar(rng, K); M = mono_batch(self.Us)
        cols = []; self.slices = []; pos = 0
        for b, C in zip(BLOCK_LIST, comp):
            if C is None: self.slices.append(None); continue
            Q = M[:, b] @ C.T; r = C.shape[0]; iu = np.triu_indices(r)
            v = np.real(np.conj(Q[:, iu[0]]) * Q[:, iu[1]]); v[:, iu[0] != iu[1]] *= 2
            cols.append(v); self.slices.append((pos, pos + len(iu[0]), r)); pos += len(iu[0])
        self.A = np.hstack(cols); self.N = pos
        u, s, vt = np.linalg.svd(self.A, full_matrices=False); keep = s > 1e-9 * s[0]
        self.Qr = u[:, keep].T; self.rank = int(keep.sum()); self.QA = self.Qr @ self.A; self.smin = s[keep][-1]; self.smax = s[0]
        self.g0, self.a, self.b = parts_batch(self.Us)
        self.F = sum((np.abs(M[:, b] @ C.T)**2).sum(1) for b, C in zip(BLOCK_LIST, comp) if C is not None)     # positive control
    def f(self, t): return t * self.g0 + (1 + t) * self.a + t * t * self.b
    def solve(self, target, objective='lambda'):
        x = cp.Variable(self.N); lam = cp.Variable(); cons = []; ones = [s[0] for s in self.slices if s is not None and s[2] == 1]
        cons.append(x[ones] >= lam)
        for s in self.slices:
            if s is None or s[2] == 1: continue
            lo, hi, r = s; iu = np.triu_indices(r); rows = [[None] * r for _ in range(r)]
            import scipy.sparse as sps
            S = sps.lil_matrix((r * r, hi - lo))
            for k, (i, j) in enumerate(zip(*iu)): S[i * r + j, k] = 1; S[j * r + i, k] = 1
            cons.append(cp.reshape(sps.csr_matrix(S) @ x[lo:hi], (r, r), order='C') - lam * np.eye(r) >> 0)
        rhs = self.Qr @ target; out = np.linalg.norm(target - self.Qr.T @ rhs) / np.linalg.norm(target)
        cons.append(self.QA @ x == rhs)
        if objective == 'center': cons.append(lam == 0)
        prob = cp.Problem(cp.Maximize(lam) if objective == 'lambda' else cp.Minimize(0), cons); prob.solve(solver='CLARABEL')
        return (None if lam.value is None else float(lam.value)), prob.status, out, x.value
if __name__ == '__main__':
    t0 = time.time(); R = Reduced(load_comp())
    print("reduced unknowns %d ; independent conditions %d (singular values %.2e .. %.2e)  [%.0fs]" % (R.N, R.rank, R.smax, R.smin, time.time() - t0), flush=True)
    lam, st, out, _ = R.solve(R.F); print("positive control F = sum |C m|^2 : lambda = %+.4e (%s) ; outside the range %.1e  [%.0fs]" % (lam, st, out, time.time() - t0), flush=True)
    for t in [float(v) for v in sys.argv[1].split(',')]:
        lam, st, out, x = R.solve(R.f(t))
        print("t = %-5g relative part of f_t outside the reduced range %.1e ; best lambda = %s (%s)  [%.0fs]" % (t, out, None if lam is None else "%+.4e" % lam, st, time.time() - t0), flush=True)
        if x is not None: pickle.dump({'t': t, 'x': x, 'lam': lam}, open('reduced4_t%g.pkl' % t, 'wb'))
