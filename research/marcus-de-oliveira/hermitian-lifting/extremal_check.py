"""Numerical check of the first-order identities of README.md (Steps 2-4) at maximisers of l(det(A + U B U^*)) over U(n).

A = diag(a) Hermitian, B normal with one real eigenvalue beta and n-1 nonreal ones, A + beta positive definite,
l(z) = Re(conj(omega) z) for a random direction omega.  A maximiser is found by BFGS from 6 Haar starts.  Checked there:
N = conj(omega) [B', adj M] is Hermitian;  (A + beta) w = 2 l(z0) v;  l(D1) = l(z0) v^*(A + beta)^{-1} v;  and
D1 = det of the compressed pair, v^*(A + beta)^{-1} v = prod(alpha_i + beta)/prod(a_i + beta).  Floating point.
usage: extremal_check.py n [seed] [number of instances]"""
import sys, numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
n = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
IU = np.triu_indices(n, 1); m2 = len(IU[0])
def herm(x):
    H = np.zeros((n, n), complex); H[IU] = x[:m2] + 1j * x[m2:2 * m2]; H = H + H.conj().T; H[np.diag_indices(n)] = x[2 * m2:]; return H
def haar():
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
def adj(M): return np.linalg.det(M) * np.linalg.inv(M)
worst = np.zeros(4)
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 6):
    a = np.sort(rng.standard_normal(n)); beta = -a[0] + rng.uniform(0.2, 1.5); b = np.r_[rng.standard_normal(n - 1) + 1j * rng.standard_normal(n - 1), beta]
    A = np.diag(a).astype(complex); om = np.exp(1j * rng.uniform(0, 2 * np.pi)); ell = lambda z: (np.conj(om) * z).real
    best = None
    for s in range(6):
        U0 = haar(); f = lambda x: -ell(np.linalg.det(A + (U0 @ expm(1j * herm(x)) * b) @ (U0 @ expm(1j * herm(x))).conj().T))
        r = minimize(f, np.zeros(n * n), method='BFGS', options={'gtol': 1e-11, 'maxiter': 3000})
        if best is None or r.fun < best[0]: best = (r.fun, U0 @ expm(1j * herm(r.x)))
    U = best[1]; Bp = (U * b) @ U.conj().T; M = A + Bp; z0 = np.linalg.det(M); v = U[:, n - 1]
    N = np.conj(om) * (Bp @ adj(M) - adj(M) @ Bp); e0 = np.abs(N - N.conj().T).max() / max(1, np.abs(N).max())
    x = adj(M) @ v; y = adj(M).conj().T @ v; w = np.conj(om) * x + om * y
    e1 = np.abs((A + beta * np.eye(n)) @ w - 2 * ell(z0) * v).max() / max(1, abs(z0))
    D1 = (v.conj() @ adj(M) @ v); g = (v.conj() @ np.linalg.inv(A + beta * np.eye(n)) @ v).real
    e2 = abs(ell(D1) - ell(z0) * g) / max(1, abs(z0))
    # compression spectrum and the ratio q/p
    Q = np.linalg.qr(np.column_stack([v, rng.standard_normal((n, n - 1)) + 1j * rng.standard_normal((n, n - 1))]))[0][:, 1:]
    alpha = np.linalg.eigvalsh(Q.conj().T @ A @ Q); e3 = abs(g - np.prod(alpha + beta) / np.prod(a + beta)); D1b = np.linalg.det(Q.conj().T @ M @ Q); e3 = max(e3, abs(D1 - D1b) / max(1, abs(D1)))
    worst = np.maximum(worst, [e0, e1, e2, e3])
print("n=%d: largest errors at the maximisers: N Hermitian %.1e ; (A+beta)w = 2 l(z0) v %.1e ; l(D1) = l(z0) v*(A+beta)^-1 v %.1e ; compression identities %.1e" % (n, *worst))
