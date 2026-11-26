"""Sample the zero set Z = {U in U(4): f_t(U) = 0} of  f_t = t*g0 + (1+t)*a + t^2*b  by local minimisation, and compute,
for every weight block, the real span N_b of {Re m_b(U), Im m_b(U) : U in Z}.  Any exact certificate must have
N_b inside the kernel of the Gram block G_b."""
import sys, time, pickle, numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
from exact_space import MONS, BLOCK_LIST
iu = np.triu_indices(4, 1)
def U_of(v):
    H = np.zeros((4, 4), complex); H[iu] = v[:6] + 1j * v[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j * v[12:16]; return expm(H)
def parts(U):
    B = np.abs(U)**2; d = lambda S, T: abs(np.linalg.det(U[np.ix_(S, T)]))**2
    return B[0, 0] + B[2, 2] + d([0, 3], [2, 3]) - d([0, 1], [0, 3]), d([0, 1], [0, 1]), d([0, 2], [0, 3])
def f_t(U, t):
    g, a, b = parts(U); return t * g + (1 + t) * a + t * t * b
def mono(U):
    flat = U.reshape(-1); return np.array([np.prod(flat[list(m)]) for m in MONS])
if __name__ == '__main__':
    t = float(sys.argv[1]); R = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0); t0 = time.time()
    Z = []; vals = []
    for r in range(R):
        res = minimize(lambda v: f_t(U_of(v), t), rng.uniform(-np.pi, np.pi, 16), method='BFGS', options={'gtol': 1e-13, 'maxiter': 3000})
        vals.append(res.fun)
        if res.fun < 1e-10: Z.append(U_of(res.x))
    vals = np.array(vals)
    print("t = %g: %d local minimisations, %d reached f < 1e-10 ; min f %.2e ; largest local-min value %.4f  [%.0fs]" % (t, R, len(Z), vals.min(), vals.max(), time.time() - t0))
    P = np.array([parts(U) for U in Z])
    print("at the zeros: max |g0| = %.1e, max a = %.1e, max b = %.1e" % (np.abs(P[:, 0]).max(), P[:, 1].max(), P[:, 2].max()))
    M = np.array([mono(U) for U in Z])
    dims = []
    for b in BLOCK_LIST:
        X = np.vstack([M[:, b].real, M[:, b].imag]); sv = np.linalg.svd(X, compute_uv=False)
        dims.append(int((sv > 1e-6 * max(sv[0], 1e-300)).sum()) if sv[0] > 1e-8 else 0)
    dims = np.array(dims); sizes = np.array([len(b) for b in BLOCK_LIST])
    print("sum over blocks: block size %d, kernel span %d, complement %d" % (sizes.sum(), dims.sum(), (sizes - dims).sum()))
    for n in sorted(set(sizes)):
        d = dims[sizes == n]; print("   blocks of size %d: %d ; kernel-dimension histogram %s" % (n, (sizes == n).sum(), dict(zip(*np.unique(d, return_counts=True)))))
    pickle.dump({'t': t, 'Z': Z}, open('zeros_t%g.pkl' % t, 'wb'))
