"""Adversarial ray-shooting test for the m-colour POVM determinantal-colouring property."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'moc')
from povm_lp import count_stats, colouring_matrix

def povm_from(Q, ranks):
    E = []; s = 0
    for r in ranks:
        B = Q[s:s+r, :]; s += r; E.append(B.conj().T @ B)
    return E

class RS:
    def __init__(self, n, m, ranks, rng):
        self.n, self.m, self.ranks = n, m, ranks
        Q = np.linalg.qr(rng.standard_normal((sum(ranks), n)) + 0j)[0]
        st = count_stats(povm_from(Q, ranks), n, m); self.keys = sorted(st.keys())
        cols, V = colouring_matrix(n, m, self.keys)
        self.cbar = V.mean(0); X = V - self.cbar
        u, s, vt = np.linalg.svd(X, full_matrices=False); r = (s > 1e-9 * s[0]).sum()
        self.B = vt[:r]; self.Vr = X @ self.B.T; self.dim = r
    def s_star(self, Q):
        st = count_stats(povm_from(Q, self.ranks), self.n, self.m)
        p = np.array([st[k] for k in self.keys]).real
        q = p - self.cbar
        off = np.linalg.norm(q - (q @ self.B.T) @ self.B)
        qr = q @ self.B.T; M = self.Vr.shape[0]
        c = np.zeros(M + 1); c[-1] = -1
        A_eq = np.vstack([np.hstack([self.Vr.T, -qr[:, None]]), np.hstack([np.ones((1, M)), [[0]]])])
        r = linprog(c, A_eq=A_eq, b_eq=np.concatenate([np.zeros(self.dim), [1.0]]), bounds=[(0, None)] * (M + 1), method='highs')
        return (r.x[-1] if r.status == 0 else np.nan), off

def search(n, m, ranks, restarts, iters, seed):
    rng = np.random.default_rng(seed); rs = RS(n, m, ranks, rng); D = sum(ranks)
    best = np.inf; offmax = 0
    for k in range(restarts):
        Z = rng.standard_normal((D, n)) + 1j * rng.standard_normal((D, n)); Q = np.linalg.qr(Z)[0]
        s, off = rs.s_star(Q); step = 0.3; offmax = max(offmax, off)
        for it in range(iters):
            ok = False
            for _ in range(6):
                Z2 = Q + step * (rng.standard_normal((D, n)) + 1j * rng.standard_normal((D, n))) / np.sqrt(D * n)
                Q2 = np.linalg.qr(Z2)[0]
                s2, off = rs.s_star(Q2); offmax = max(offmax, off)
                if s2 < s: Q, s, ok = Q2, s2, True; break
            step = step * 1.3 if ok else step * 0.6
            if step < 1e-5: break
        best = min(best, s)
    return best, rs.dim, offmax

if __name__ == '__main__':
    for spec in sys.argv[1:]:
        n, m, rk = spec.split(':'); n, m = int(n), int(m); ranks = list(map(int, rk.split(',')))
        b, dim, off = search(n, m, ranks, 12, 250, 1)
        print("n=%d colours=%d ranks=%s affdim=%d : min s* = %.6f   (max off-span residual %.1e)" % (n, m, ranks, dim, b, off), flush=True)
