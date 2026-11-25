"""Sampled, torus-symmetric sums-of-squares certificates on U(4).
A function f on U(4) (invariant under row and column phases) is certified nonnegative at degree d if
    f(U) = sum over weight blocks  m_b(U)^* G_b m_b(U)   for all U,   G_b real symmetric PSD,
where m_b lists the degree-d monomials in the entries of U with a given (row multiset, column multiset).
The identity is imposed on random unitaries (more than the rank of the constraint system)."""
import itertools, sys, time, numpy as np, cvxpy as cp
from base import *
class Level:
    def __init__(self, d, K, seed=0):
        self.d = d
        mons = list(itertools.combinations_with_replacement(range(16), d))
        blocks = {}
        for i, m in enumerate(mons):
            key = (tuple(np.bincount([v // 4 for v in m], minlength=4)), tuple(np.bincount([v % 4 for v in m], minlength=4)))
            blocks.setdefault(key, []).append(i)
        self.blocks = list(blocks.values()); self.mons = mons
        rng = np.random.default_rng(seed)
        self.U = [haar(rng) for _ in range(K)]
        M = np.ones((K, len(mons)), complex)
        flat = np.array([U.reshape(-1) for U in self.U])
        for i, m in enumerate(mons):
            for v in m: M[:, i] *= flat[:, v]
        # variables: upper-triangular entries of each block
        self.slices = []; cols = []; pos = 0
        for b in self.blocks:
            n = len(b); iu = np.triu_indices(n)
            vals = np.real(np.conj(M[:, [b[p] for p in iu[0]]]) * M[:, [b[q] for q in iu[1]]])
            vals[:, iu[0] != iu[1]] *= 2
            cols.append(vals); self.slices.append((pos, pos + len(iu[0]), n)); pos += len(iu[0])
        self.A = np.hstack(cols); self.nvar = pos
        # independent rows via QR with pivoting on A^T
        u, s, vt = np.linalg.svd(self.A, full_matrices=False)
        self.rank = int((s > 1e-9 * s[0]).sum())
        self.P = u[:, :self.rank].T                         # projects sample equations onto independent combinations
        self.Ar = self.P @ self.A
        self.X = np.array([x_U(U) for U in self.U])
    def gram_blocks(self, g):
        out = []
        for (lo, hi, n) in self.slices:
            iu = np.triu_indices(n); G = cp.Variable((n, n), symmetric=True) if False else None
            out.append((lo, hi, n))
        return out
    def solve(self, fvals_list, maximize_gamma=False, solver='CLARABEL', verbose=False):
        """Half-line ansatz for f_0 + t f_1 + t^2 f_2 >= 0 (t >= 0) if three functions are given; plain SOS if one.
        With maximize_gamma: certify f - gamma >= 0 for the largest gamma (one function)."""
        cons = []; T = len(fvals_list)
        gvars = [cp.Variable(self.nvar) for _ in range(T)]
        for (lo, hi, n) in self.slices:
            iu = np.triu_indices(n)
            def mat(expr_of):   # symmetric n x n from a triu vector expression
                Mx = [[None] * n for _ in range(n)]
                for k, (p, q) in enumerate(zip(*iu)):
                    Mx[p][q] = expr_of(k); Mx[q][p] = expr_of(k)
                return cp.bmat(Mx)
            if T == 1:
                cons.append(mat(lambda k: gvars[0][lo + k]) >> 0)
            else:
                M0 = cp.Variable((2 * n, 2 * n), PSD=True); M1 = cp.Variable((n, n), PSD=True)
                for k, (p, q) in enumerate(zip(*iu)):
                    cons.append(gvars[0][lo + k] == M0[p, q])
                    cons.append(gvars[1][lo + k] == M0[p, n + q] + M0[q, n + p] + M1[p, q])
                    cons.append(gvars[2][lo + k] == M0[n + p, n + q])
        gamma = cp.Variable() if maximize_gamma else None
        for t in range(T):
            rhs = self.P @ fvals_list[t]
            if maximize_gamma: cons.append(self.Ar @ gvars[t] == rhs - gamma * (self.P @ np.ones(len(self.U))))
            else: cons.append(self.Ar @ gvars[t] == rhs)
        prob = cp.Problem(cp.Maximize(gamma) if maximize_gamma else cp.Minimize(0), cons)
        prob.solve(solver=solver, verbose=verbose)
        return prob.status, (gamma.value if maximize_gamma else None)
W = lambda terms: sum(c * np.eye(NX)[iB(a, b) if k == 'B' else iD(a, b)] for c, k, a, b in terms)
g0 = W([(1,'B',0,0),(1,'B',2,2),(1,'D',(0,3),(2,3)),(-1,'D',(0,1),(0,3))])
a_ = W([(1,'D',(0,1),(0,1))]); b_ = W([(1,'D',(0,2),(0,3))])
if __name__ == '__main__':
    d = int(sys.argv[1]); K = int(sys.argv[2]); t0 = time.time()
    L = Level(d, K)
    print("degree %d: %d monomials, %d blocks (max size %d), %d Gram unknowns, constraint rank %d from %d samples  [%.0fs]" %
          (d, len(L.mons), len(L.blocks), max(len(b) for b in L.blocks), L.nvar, L.rank, K, time.time() - t0), flush=True)
    st, gam = L.solve([L.X @ g0], maximize_gamma=True)
    print("  best certified lower bound for g0 at this degree: %s (status %s); true minimum -0.067442  [%.0fs]" % (None if gam is None else "%.6f" % gam, st, time.time() - t0), flush=True)
    for t in [0.3, 1.0, 3.0]:
        st, gam = L.solve([L.X @ (t * g0 + (1 + t) * a_ + t * t * b_)], maximize_gamma=True)
        print("  t=%g: best certified lower bound for t*g0+(1+t)*D[12,12]+t^2*D[13,14]: %s (%s)  [%.0fs]" % (t, None if gam is None else "%.3e" % gam, st, time.time() - t0), flush=True)
