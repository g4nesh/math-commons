"""Walk along the zero set: from each known zero, perturb by a random unitary close to the identity and re-minimise.
Repeat for several rounds and track the dimension of the monomial spans per weight block until it stabilises."""
import sys, time, pickle, numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm, logm
from exact_space import MONS, BLOCK_LIST
from zeros import U_of, f_t, mono
def span_dims(Zs):
    M = np.array([mono(U) for U in Zs]); dims = []
    for b in BLOCK_LIST:
        X = np.vstack([M[:, b].real, M[:, b].imag]); sv = np.linalg.svd(X, compute_uv=False)
        dims.append(int((sv > 1e-6 * sv[0]).sum()) if sv[0] > 1e-8 else 0)
    return np.array(dims)
def skew(rng):
    H = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4)); return (H - H.conj().T) / 2
def reminimise(U0, t, rng, eps):
    V = U0 @ expm(eps * skew(rng))
    # parametrise near V:  U = V expm(K(v))
    iu = np.triu_indices(4, 1)
    def K(v):
        H = np.zeros((4, 4), complex); H[iu] = v[:6] + 1j * v[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j * v[12:16]; return H
    res = minimize(lambda v: f_t(V @ expm(K(v)), t), np.zeros(16), method='BFGS', options={'gtol': 1e-13, 'maxiter': 2000})
    return (V @ expm(K(res.x))) if res.fun < 1e-10 else None
if __name__ == '__main__':
    t = float(sys.argv[1]); rounds = int(sys.argv[2]); per = int(sys.argv[3]); rng = np.random.default_rng(5); t0 = time.time()
    Z = pickle.load(open('zeros_t%g.pkl' % t, 'rb'))['Z']
    sizes = np.array([len(b) for b in BLOCK_LIST]); d = span_dims(Z)
    print("start: %d zeros, kernel span %d, complement %d" % (len(Z), d.sum(), (sizes - d).sum()), flush=True)
    front = Z[::max(1, len(Z) // per)][:per]
    for r in range(rounds):
        new = []
        for U0 in front:
            for eps in (0.15, 0.5):
                U1 = reminimise(U0, t, rng, eps)
                if U1 is not None: new.append(U1)
        Z = Z + new; d = span_dims(Z)
        print("round %d: +%d zeros (total %d), kernel span %d, complement %d  [%.0fs]" % (r + 1, len(new), len(Z), d.sum(), (sizes - d).sum(), time.time() - t0), flush=True)
        idx = rng.permutation(len(new))[:per]; front = [new[i] for i in idx]
    pickle.dump({'t': t, 'Z': Z}, open('zeros_walk_t%g.pkl' % t, 'wb'))
    for n in sorted(set(sizes)):
        dd = d[sizes == n]; print("   blocks of size %d: %d ; kernel-dimension histogram %s" % (n, (sizes == n).sum(), dict(zip(*np.unique(dd, return_counts=True)))))
