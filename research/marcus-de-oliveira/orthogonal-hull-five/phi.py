"""The functional Phi on squared minors of 5x5 matrices attached to the cyclic one-factorisation of K_6 (indices mod 5):

Phi(U) = - sum_v |u_vv|^2 + sum_v ( |det U[{v+1,v-1},{v+2,v-2}]|^2 + |det U[{v+2,v-2},{v+1,v-1}]|^2 )
         - sum_v sum_{d=1,2} ( |det U[{v,v+d},{v,v-d}]|^2 + |det U[{v,v-d},{v,v+d}]|^2 ).

Run as a script it reports, in floating point: the values on the 120 permutation matrices; local maxima over O(5) and over
U(5) from Haar-random starts (BFGS on the group); the number of zero Hessian eigenvalues at the orthogonal maxima.
usage: phi.py [orthogonal starts] [unitary starts]"""
import itertools, sys, numpy as np
n = 5
def terms():
    out = [((v,), (v,), -1) for v in range(n)]
    for v in range(n):
        a = tuple(sorted(((v + 1) % n, (v - 1) % n))); b = tuple(sorted(((v + 2) % n, (v - 2) % n))); out += [(a, b, 1), (b, a, 1)]
        for d in (1, 2):
            s = tuple(sorted((v, (v + d) % n))); t = tuple(sorted((v, (v - d) % n))); out += [(s, t, -1), (t, s, -1)]
    return out
TERMS = terms()
def phi(M): return sum(c * abs(np.linalg.det(M[np.ix_(S, T)]))**2 for S, T, c in TERMS)
if __name__ == '__main__':
    from scipy.optimize import minimize
    from scipy.linalg import expm
    rng = np.random.default_rng(0); NO = int(sys.argv[1]) if len(sys.argv) > 1 else 300; NUs = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    vals = [round(phi(np.eye(n)[list(p)])) for p in itertools.permutations(range(n))]; print("values on the 120 permutation matrices (value: count):", {int(k): int(c) for k, c in zip(*np.unique(vals, return_counts=True))})
    IU = np.triu_indices(n, 1)
    def skew(x): K = np.zeros((n, n)); K[IU] = x; return K - K.T
    def herm(x):
        H = np.zeros((n, n), complex); H[IU] = x[:10] + 1j * x[10:20]; H = H + H.conj().T; H[np.diag_indices(n)] = x[20:]; return H
    def haar_o():
        q, r = np.linalg.qr(rng.standard_normal((n, n))); return q * np.sign(np.diag(r))
    def haar_u():
        z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
    def hessian_zeros(O, h=1e-4):
        f = lambda x: phi(O @ expm(skew(x))); H = np.zeros((10, 10))
        for i in range(10):
            for j in range(i, 10):
                ei = np.zeros(10); ej = np.zeros(10); ei[i] = h; ej[j] = h; H[i, j] = H[j, i] = (f(ei + ej) - f(ei - ej) - f(-ei + ej) + f(-ei - ej)) / (4 * h * h)
        ev = np.linalg.eigvalsh(H); return int((np.abs(ev) < 1e-4).sum()), ev.max()
    ro = []; nz = {}; notperm = 0
    for k in range(NO):
        O0 = haar_o(); r = minimize(lambda x: -phi(O0 @ expm(skew(x))), np.zeros(10), method='BFGS', options={'gtol': 1e-11}); ro.append(-r.fun)
        if k < 40:
            O1 = O0 @ expm(skew(r.x)); z, top = hessian_zeros(O1); nz[z] = nz.get(z, 0) + 1; notperm += np.abs(np.abs(O1) - np.round(np.abs(O1))).max() > 1e-3
    print("O(5): %d local searches: largest value %.10f, smallest %.10f" % (NO, max(ro), min(ro)))
    print("      at the first 40 maxima: number of zero Hessian eigenvalues (count): %s ; maxima that are not signed permutation matrices: %d" % (nz, notperm))
    ru = []
    for _ in range(NUs):
        U0 = haar_u(); r = minimize(lambda x: -phi(U0 @ expm(1j * herm(x))), np.zeros(25), method='BFGS', options={'gtol': 1e-10}); ru.append(-r.fun)
    print("U(5): %d local searches: largest value %.10f, smallest %.10f ; 37/25 = %.2f" % (NUs, max(ru), min(ru), 37 / 25))
    s = np.array([phi(haar_o()) for _ in range(100000)]); print("100,000 Haar orthogonal samples: largest value %.5f" % s.max())
