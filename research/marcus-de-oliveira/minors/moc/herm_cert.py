"""Real-a ('Hermitian A') certificates.  For real a, d = sum_T b^T phi_T(p_U) with phi_T(x) = sum_{|S|=|T|} a^{S^c} x_{S,T}
real.  If Phi_a(p_U) is in conv{Phi_a(p_sigma)}, MOC holds at (U, a) for EVERY complex b (weights independent of b)."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
sys.path.insert(0, 'moc')
from moclib import haar_unitary

def all_minor_sq(U):
    n = U.shape[0]; out = {}
    for k in range(0, n + 1):
        for S in itertools.combinations(range(n), k):
            for T in itertools.combinations(range(n), k):
                out[(S, T)] = 1.0 if k == 0 else abs(np.linalg.det(U[np.ix_(S, T)]))**2
    return out

def Phi(xd, a, n):
    v = []
    for k in range(1, n):
        for T in itertools.combinations(range(n), k):
            v.append(sum(np.prod([a[i] for i in range(n) if i not in S]) * xd[(S, T)] for S in itertools.combinations(range(n), k)))
    return np.array(v)

class HermCert:
    def __init__(self, n):
        self.n = n; self.perms = list(itertools.permutations(range(n)))
        self.Pd = []
        for s in self.perms:
            P = np.zeros((n, n)); P[range(n), s] = 1; self.Pd.append(all_minor_sq(P))
    def depth(self, U, a):
        """ray-shooting depth s* from the vertex centroid in Phi_a space (s*<1 <=> no certificate)."""
        V = np.array([Phi(pd, a, self.n) for pd in self.Pd]); p = Phi(all_minor_sq(U), a, self.n)
        c = V.mean(0); X = V - c
        u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); B = vt[:r]
        q = p - c; off = np.linalg.norm(q - (q @ B.T) @ B) / (1 + np.linalg.norm(q))
        Vr = X @ B.T; qr = q @ B.T; m = len(V)
        res = linprog(np.r_[np.zeros(m), -1.0], A_eq=np.vstack([np.hstack([Vr.T, -qr[:, None]]), np.r_[np.ones(m), 0][None]]),
                      b_eq=np.r_[np.zeros(r), 1.0], bounds=[(0, None)] * (m + 1), method='highs')
        return (res.x[-1] if res.status == 0 else np.nan), r, off

if __name__ == '__main__':
    hc = HermCert(4); rng = np.random.default_rng(0)
    Ud = np.load('moc/runs/deepest4.npy'); Uo = np.load('moc/runs/orth_best.npy')
    for name, U in [('deepest unitary violator', Ud), ('deepest orthogonal violator', Uo)]:
        ds = []
        for _ in range(300):
            a = rng.standard_normal(4) * rng.choice([0.3, 1, 3]); s, r, off = hc.depth(U, a); ds.append(s)
        ds = np.array(ds); print("%s: rank(Phi_a on aff Q) = %d, off-span residual %.1e; depth s* over 300 random real a: min %.4f, median %.4f (s*<1 => fails)" % (name, r, off, ds.min(), np.median(ds)))
    bad = 0; mins = []
    for _ in range(300):
        U = haar_unitary(4, rng); a = rng.standard_normal(4); s, r, off = hc.depth(U, a); mins.append(s)
    print("300 Haar U x random real a: min depth %.4f, fraction with s*<1: %.3f" % (min(mins), np.mean(np.array(mins) < 1 - 1e-9)))
