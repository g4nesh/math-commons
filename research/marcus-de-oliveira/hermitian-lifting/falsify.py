"""Falsification attempt for Corollary L4 of README.md: A Hermitian n x n, B normal with exactly four nonreal eigenvalues.

Signed depth of det(A + U B U^*) in the hull of the n! products: the largest s such that the point is an affine combination
of the products with all weights >= s (negative = outside).  Random triples (a, b, U), then Nelder-Mead minimisation of the
depth over spectra and U.  Floating point; a test, not a proof.
usage: falsify.py n [random triples] [adversarial runs] [seed]"""
import itertools, sys, numpy as np
from scipy.optimize import linprog, minimize
from scipy.linalg import expm
n = int(sys.argv[1]); NR = int(sys.argv[2]); NA = int(sys.argv[3]); rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
PERMS = np.array(list(itertools.permutations(range(n))))
def haar():
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
def depth(a, b, U):
    z = np.linalg.det(np.diag(a) + (U * b) @ U.conj().T); pts = np.prod(a[None, :] + b[PERMS], axis=1); sc = np.abs(pts).max(); m = len(pts)
    P = np.column_stack([pts.real, pts.imag]) / sc; zz = np.array([z.real, z.imag]) / sc
    res = linprog(np.r_[np.zeros(m), -1.0], A_ub=np.hstack([-np.eye(m), np.ones((m, 1))]), b_ub=np.zeros(m), A_eq=np.vstack([np.hstack([P.T, np.zeros((2, 1))]), np.r_[np.ones(m), 0][None]]), b_eq=np.r_[zz, 1.0], bounds=[(None, None)] * (m + 1), method='highs')
    return -res.fun if res.status == 0 else np.nan
def spectra():
    a = rng.standard_normal(n); b = rng.standard_normal(n).astype(complex); b[:4] += 1j * rng.standard_normal(4); return a, b
d = [depth(*spectra(), haar()) for _ in range(NR)]; print("n=%d: %d random triples (A Hermitian, B with four nonreal eigenvalues): smallest depth %.3e, fraction below 1e-9: %.3f" % (n, NR, np.nanmin(d), np.mean(np.array(d) < 1e-9)), flush=True)
IU = np.triu_indices(n, 1); m2 = len(IU[0])
def herm(x):
    H = np.zeros((n, n), complex); H[IU] = x[:m2] + 1j * x[m2:2 * m2]; H = H + H.conj().T; H[np.diag_indices(n)] = x[2 * m2:]; return H
worst = np.inf
for _ in range(NA):
    U0 = haar(); a0, b0 = spectra()
    def f(x):
        a = x[:n]; b = x[n:2 * n].astype(complex); b[:4] += 1j * x[2 * n:2 * n + 4]; return depth(a, b, U0 @ expm(1j * herm(x[2 * n + 4:])))
    x0 = np.r_[a0, b0.real, b0.imag[:4], np.zeros(n * n)]; r = minimize(f, x0, method='Nelder-Mead', options={'maxiter': 2500, 'xatol': 1e-8, 'fatol': 1e-12, 'adaptive': True}); worst = min(worst, r.fun)
print("n=%d: %d adversarial runs over spectra and U (Nelder-Mead on the depth): smallest depth reached %.3e" % (n, NA, worst), flush=True)
