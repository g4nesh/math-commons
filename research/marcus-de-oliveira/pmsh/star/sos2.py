"""Faster sampled torus-symmetric SOS on U(4) (same model as sos_torus.py, vectorised for degree 4)."""
import itertools, sys, time, numpy as np, cvxpy as cp, scipy.sparse as sps
from base import *
class Level:
    def __init__(self, d, K, seed=0):
        self.d = d; mons = list(itertools.combinations_with_replacement(range(16), d)); blocks = {}
        for i, m in enumerate(mons):
            key = (tuple(np.bincount([v // 4 for v in m], minlength=4)), tuple(np.bincount([v % 4 for v in m], minlength=4)))
            blocks.setdefault(key, []).append(i)
        self.blocks = list(blocks.values()); self.keys = list(blocks.keys()); self.mons = mons
        rng = np.random.default_rng(seed); self.U = [haar(rng) for _ in range(K)]
        flat = np.array([U.reshape(-1) for U in self.U]); M = np.ones((K, len(mons)), complex)
        for i, m in enumerate(mons):
            for v in m: M[:, i] *= flat[:, v]
        self.slices = []; cols = []; pos = 0
        for b in self.blocks:
            n = len(b); iu = np.triu_indices(n)
            vals = np.real(np.conj(M[:, [b[p] for p in iu[0]]]) * M[:, [b[q] for q in iu[1]]]); vals[:, iu[0] != iu[1]] *= 2
            cols.append(vals); self.slices.append((pos, pos + len(iu[0]), n)); pos += len(iu[0])
        A = np.hstack(cols); self.nvar = pos
        w, V = np.linalg.eigh(A @ A.T)                       # K x K
        keep = w > 1e-10 * w[-1]; self.rank = int(keep.sum())
        self.P = V[:, keep].T; self.Ar = self.P @ A; self.X = np.array([x_U(U) for U in self.U]); self.ones = self.P @ np.ones(K)
    def solve(self, fvals, maximize_gamma=True, solver='CLARABEL', verbose=False):
        g = cp.Variable(self.nvar); cons = []; ones = [lo for (lo, hi, n) in self.slices if n == 1]
        cons.append(g[ones] >= 0)
        for (lo, hi, n) in self.slices:
            if n == 1: continue
            iu = np.triu_indices(n); S = sps.lil_matrix((n * n, hi - lo))
            for k, (p, q) in enumerate(zip(*iu)): S[p * n + q, k] = 1; S[q * n + p, k] = 1
            cons.append(cp.reshape(sps.csr_matrix(S) @ g[lo:hi], (n, n), order='C') >> 0)
        gamma = cp.Variable()
        cons.append(self.Ar @ g == self.P @ fvals - (gamma * self.ones if maximize_gamma else 0))
        if not maximize_gamma: cons.append(gamma == 0)
        prob = cp.Problem(cp.Maximize(gamma), cons); prob.solve(solver=solver, verbose=verbose)
        self.g = g.value
        return prob.status, gamma.value
W = lambda terms: sum(c * np.eye(NX)[iB(a, b) if k == 'B' else iD(a, b)] for c, k, a, b in terms)
g0 = W([(1,'B',0,0),(1,'B',2,2),(1,'D',(0,3),(2,3)),(-1,'D',(0,1),(0,3))])
a_ = W([(1,'D',(0,1),(0,1))]); b_ = W([(1,'D',(0,2),(0,3))]); bw = W([(1,'D',(0,2),(1,2))])
if __name__ == '__main__':
    d, K = int(sys.argv[1]), int(sys.argv[2]); ts = [float(v) for v in sys.argv[3].split(',')]; t0 = time.time()
    L = Level(d, K)
    print("degree %d: %d monomials, %d blocks, %d Gram unknowns, constraint rank %d from %d samples  [%.0fs]" % (d, len(L.mons), len(L.blocks), L.nvar, L.rank, K, time.time() - t0), flush=True)
    for t in ts:
        st, gam = L.solve(L.X @ (t * g0 + (1 + t) * a_ + t * t * b_))
        print("  t = %-5g certified lower bound / t = %+.3e  (%s)  [%.0fs]" % (t, gam / t, st, time.time() - t0), flush=True)
