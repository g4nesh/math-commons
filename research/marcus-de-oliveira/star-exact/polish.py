"""Polish sampled zeros of f_t to near machine precision (Newton on the exact gradient, pseudo-inverse of the Hessian),
then recompute the monomial spans per weight block with a clean singular-value gap."""
import sys, time, pickle, numpy as np
from scipy.linalg import expm
from exact_space import MONS, BLOCK_LIST
from zeros import mono, parts
TERMS = lambda t: [(t, 'B', 0, 0), (t, 'B', 2, 2), (t, 'D', (0, 3), (2, 3)), (-t, 'D', (0, 1), (0, 3)), (1 + t, 'D', (0, 1), (0, 1)), (t * t, 'D', (0, 2), (0, 3))]
BASIS = []
for j in range(4):
    K = np.zeros((4, 4), complex); K[j, j] = 1j; BASIS.append(K)
for j in range(4):
    for k in range(j + 1, 4):
        K = np.zeros((4, 4), complex); K[j, k] = 1; K[k, j] = -1; BASIS.append(K)
        K = np.zeros((4, 4), complex); K[j, k] = 1j; K[k, j] = 1j; BASIS.append(K)
BASIS = np.array(BASIS)
def wirtinger(U, t):
    """G[i,j] = d f / d conj(u_ij)"""
    G = np.zeros((4, 4), complex)
    for c, kind, a, b in TERMS(t):
        if kind == 'B': G[a, b] += c * U[a, b]
        else:
            (s1, s2), (t1, t2) = a, b; d = U[s1, t1] * U[s2, t2] - U[s1, t2] * U[s2, t1]
            G[s1, t1] += c * d * np.conj(U[s2, t2]); G[s2, t2] += c * d * np.conj(U[s1, t1])
            G[s1, t2] -= c * d * np.conj(U[s2, t1]); G[s2, t1] -= c * d * np.conj(U[s1, t2])
    return G
def grad(U, t):
    G = wirtinger(U, t); X = G.conj().T @ U
    return np.array([2 * np.real(np.trace(X @ K)) for K in BASIS])
def polish(U, t, iters=12):
    for _ in range(iters):
        g = grad(U, t)
        if np.linalg.norm(g) < 1e-15: break
        H = np.zeros((16, 16)); h = 1e-5
        for k in range(16):
            H[:, k] = (grad(U @ expm(h * BASIS[k]), t) - grad(U @ expm(-h * BASIS[k]), t)) / (2 * h)
        H = (H + H.T) / 2; w, V = np.linalg.eigh(H); keep = w > 1e-6 * w.max()
        step = -V[:, keep] @ ((V[:, keep].T @ g) / w[keep])
        U = U @ expm(np.einsum('k,kij->ij', step, BASIS))
    u, s, vh = np.linalg.svd(U); return u @ vh
if __name__ == '__main__':
    t = float(sys.argv[1]); src = sys.argv[2]; t0 = time.time()
    Z = pickle.load(open(src, 'rb'))['Z']
    Zp = [polish(U, t) for U in Z]
    res = np.array([[abs(x) for x in parts(U)] for U in Zp]); gn = np.array([np.linalg.norm(grad(U, t)) for U in Zp])
    ok = (res.max(1) < 1e-13)
    print("polished %d zeros: max gradient norm %.1e ; max of |g0|, a, b: %.1e ; accepted (all three < 1e-13): %d  [%.0fs]" % (len(Z), gn.max(), res.max(), ok.sum(), time.time() - t0))
    Zp = [U for U, o in zip(Zp, ok) if o]
    M = np.array([mono(U) for U in Zp]); sizes = np.array([len(b) for b in BLOCK_LIST]); dims = []; gaps = []
    for b in BLOCK_LIST:
        X = np.vstack([M[:, b].real, M[:, b].imag]); sv = np.linalg.svd(X, compute_uv=False)
        sv = np.concatenate([sv, np.zeros(max(0, len(b) - len(sv)))])[:len(b)]
        r = int((sv > 1e-8 * max(sv[0], 1e-300)).sum()) if sv[0] > 1e-10 else 0
        dims.append(r)
        if 0 < r < len(b): gaps.append(sv[r] / sv[r - 1])
    dims = np.array(dims)
    print("kernel span %d, complement %d ; largest ratio (first discarded / last kept singular value) = %.1e" % (dims.sum(), (sizes - dims).sum(), max(gaps) if gaps else 0))
    for n in sorted(set(sizes)):
        d = dims[sizes == n]; print("   blocks of size %d: %d ; kernel-dimension histogram %s" % (n, (sizes == n).sum(), dict(zip(*np.unique(d, return_counts=True)))))
    pickle.dump({'t': t, 'Z': Zp, 'dims': dims}, open('zeros_polished_t%g.pkl' % t, 'wb'))
