"""Adversarial test of the pattern compound-Birkhoff property.
s*(U) = max{s : cbar + s (p_U - cbar) in conv(class vectors)}  (ray shooting from centroid).
s* < 1  <=>  p_U outside the polytope (no probabilistic certificate).  We minimise s* over U."""
import sys, itertools, numpy as np, time
from scipy.optimize import linprog
from scipy.linalg import expm
sys.path.insert(0, 'moc')
from moclib import haar_unitary

class Pattern:
    def __init__(self, pattern):
        self.pattern = pattern; n = self.n = sum(pattern)
        lab = np.empty(n, int); s = 0
        for j, p in enumerate(pattern): lab[s:s+p] = j; s += p
        self.lab = lab; m = len(pattern)
        # index structures for all k-minors
        self.idx = []
        keys = {}
        for k in range(1, n):
            Ss = list(itertools.combinations(range(n), k))
            Ts = Ss
            Sa = np.array(Ss); Ta = np.array(Ts)
            kap = [tuple(np.bincount(lab[list(T)], minlength=m)) for T in Ts]
            col = np.empty((len(Ss), len(Ts)), int)
            for i, S in enumerate(Ss):
                for j, T in enumerate(Ts):
                    key = (k, S, kap[j])
                    if key not in keys: keys[key] = len(keys)
                    col[i, j] = keys[key]
            self.idx.append((k, Sa, Ta, col))
        self.nkeys = len(keys)
        # class representatives
        classes = {}
        for s_ in itertools.permutations(range(n)):
            key = tuple(lab[list(s_)])
            if key not in classes: classes[key] = s_
        V = []
        for s_ in classes.values():
            P = np.zeros((n, n)); P[range(n), s_] = 1
            V.append(self.stats(P))
        V = np.array(V)
        self.cbar = V.mean(0)
        X = V - self.cbar
        u, sv, vt = np.linalg.svd(X, full_matrices=False)
        r = (sv > 1e-9 * sv[0]).sum()
        self.B = vt[:r]          # basis of affine span directions
        self.Vr = X @ self.B.T   # reduced vertex coords
        self.dim = r; self.ncls = len(V)

    def stats(self, U):
        out = np.zeros(self.nkeys)
        for k, Sa, Ta, col in self.idx:
            sub = U[Sa[:, None, :, None], Ta[None, :, None, :]]   # (|S|,|T|,k,k)
            dets = np.abs(np.linalg.det(sub))**2
            np.add.at(out, col.ravel(), dets.ravel())
        return out

    def ray(self, U):
        p = self.stats(U)
        q = (p - self.cbar)
        res_off = np.linalg.norm(q - (q @ self.B.T) @ self.B)
        qr = q @ self.B.T
        m = self.ncls
        # vars lambda (m), s ; maximise s : Vr^T lambda - s qr = 0, sum lambda = 1
        c = np.zeros(m + 1); c[-1] = -1
        A_eq = np.vstack([np.hstack([self.Vr.T, -qr[:, None]]), np.hstack([np.ones((1, m)), [[0]]])])
        b_eq = np.concatenate([np.zeros(self.dim), [1.0]])
        r = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * m + [(0, None)], method='highs')
        if r.status != 0:
            return np.nan, None, res_off
        s = r.x[-1]
        y = r.eqlin.marginals[:self.dim]   # dual for direction constraints
        return s, y, res_off

def skew(x, n):
    H = np.zeros((n, n), complex); iu = np.triu_indices(n, 1); m = len(iu[0])
    H[np.diag_indices(n)] = 1j * x[:n]
    H[iu] = x[n:n+m] + 1j * x[n+m:]
    return H - np.conj(H.T) * (np.triu(np.ones((n, n)), 1).T > 0) - np.triu(H, 1).conj().T * 0

def minimise(pat, restarts, iters, seed=0, verbose=True):
    rng = np.random.default_rng(seed)
    n = pat.n; dimU = n * n
    best = (np.inf, None)
    for rs in range(restarts):
        U = haar_unitary(n, rng)
        s, y, off = pat.ray(U)
        step = 0.3
        for it in range(iters):
            # random-direction descent with adaptive step (robust to LP kinks)
            improved = False
            for _ in range(6):
                x = rng.standard_normal(dimU)
                H = np.zeros((n, n), complex); iu = np.triu_indices(n, 1); m = len(iu[0])
                H[np.diag_indices(n)] = 1j * x[:n]; H[iu] = x[n:n+m] + 1j * x[n+m:]
                H = H - np.triu(H, 1).conj().T
                U2 = U @ expm(step * H / np.linalg.norm(x))
                s2, y2, off2 = pat.ray(U2)
                if s2 < s:
                    U, s = U2, s2; improved = True; break
            step = step * 1.3 if improved else step * 0.6
            if step < 1e-5: break
        if s < best[0]: best = (s, U)
        if verbose: print("  restart %d: s*=%.6f (best %.6f)" % (rs, s, best[0]), flush=True)
    return best

if __name__ == '__main__':
    pat = Pattern(tuple(map(int, sys.argv[1].split(','))))
    R, I = int(sys.argv[2]), int(sys.argv[3])
    print("pattern", pat.pattern, "classes", pat.ncls, "affine dim", pat.dim, flush=True)
    t = time.time()
    best = minimise(pat, R, I, seed=int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    print("RESULT pattern %s: min s* = %.6f  (s*<1 means certificate fails)  [%.0fs]" % (pat.pattern, best[0], time.time() - t))
    np.save('moc/runs/adv_%s.npy' % '-'.join(map(str, pat.pattern)), best[1])
