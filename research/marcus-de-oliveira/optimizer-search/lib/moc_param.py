"""Parametrizations: unitaries from skew-Hermitian K, normalized spectra, embeddings."""
import numpy as np
from scipy.linalg import expm, schur


def skew_from_real(X):
    """Real n x n X -> skew-Hermitian K = (X - X^T)/2 + i (X + X^T)/2 (a bijection)."""
    X = np.asarray(X, float)
    return 0.5 * (X - X.T) + 0.5j * (X + X.T)


def real_from_skew(K):
    K = np.asarray(K, complex)
    return K.real + K.imag


def cayley(K):
    """U = (I - K)(I + K)^{-1}; unitary for skew-Hermitian K; never has eigenvalue -1."""
    I = np.eye(len(K))
    return np.linalg.solve(I + K, I - K)


def inverse_cayley(U):
    """K = (I - U)(I + U)^{-1} (requires -1 not an eigenvalue of U)."""
    I = np.eye(len(U))
    return np.linalg.solve(I + U, I - U)


def expm_unitary(K):
    return expm(np.asarray(K, complex))


def log_unitary(U):
    """Skew-Hermitian K with expm(K) = U (principal branch)."""
    T, Z = schur(np.asarray(U, complex), output="complex")
    return Z @ np.diag(1j * np.angle(np.diag(T))) @ Z.conj().T


def cayley_coords(U, tries=64, rng=None):
    """Find diagonal phases R and K with U @ diag(R) = cayley(K), choosing R so
    that the eigenvalues of U diag(R) stay far from -1 (well-conditioned K).
    d(U diag(R)) = d(U), so (a, b, U diag(R)) is an equivalent configuration."""
    rng = np.random.default_rng(0) if rng is None else rng
    n = len(U)
    best = None
    for t in range(tries):
        R = np.ones(n) if t == 0 else np.exp(2j * np.pi * rng.random(n))
        V = np.asarray(U) * R[None, :]
        gap = np.abs(np.linalg.eigvals(V) + 1).min()
        if best is None or gap > best[0]:
            best = (gap, R, V)
        if gap > 0.7:
            break
    gap, R, V = best
    return inverse_cayley(V), R, float(gap)


def polar_unitary(U):
    """Nearest unitary (polar factor)."""
    W, _, Vh = np.linalg.svd(np.asarray(U, complex))
    return W @ Vh


def random_unitary(n, rng):
    """Haar-random unitary."""
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))[None, :]


def random_spectra(n, rng, kind="gaussian"):
    """kind: 'gaussian' (complex normal), 'circle' (unit circle), 'mixed'
    (independent random choice of the two per spectrum)."""
    def one(k):
        if k == "gaussian":
            return rng.standard_normal(n) + 1j * rng.standard_normal(n)
        if k == "circle":
            return np.exp(2j * np.pi * rng.random(n))
        raise ValueError(k)
    if kind == "mixed":
        return one(rng.choice(["gaussian", "circle"])), one(rng.choice(["gaussian", "circle"]))
    return one(kind), one(kind)


def canonical_spectra(a, b):
    """Fix the shift and scale invariances: (a, b) -> (c(a + t), c(b - t)) with
    t = mean(b) (so mean b = 0), |c| making the rms of (a - mean a, b) equal
    to 1, and arg c making a[0] - a[1] real positive (when nonzero).
    d and every z_sigma are multiplied by the same c^n, so m is unchanged."""
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    t = b.mean()
    a2, b2 = a + t, b - t
    ca = a2 - a2.mean()
    rms = np.sqrt((np.sum(np.abs(ca) ** 2) + np.sum(np.abs(b2) ** 2)) / (2 * len(a)))
    c = 1.0 / rms if rms > 0 else 1.0
    if len(a) > 1 and abs(a2[0] - a2[1]) > 0:
        c = c * np.exp(-1j * np.angle(a2[0] - a2[1]))
    return c * a2, c * b2


def embed(a, b, U, a_extra, b_extra, U_extra=None):
    """Direct sum: (a + a_extra, b + b_extra, U (+) U_extra), U_extra = I by default."""
    a_extra = np.atleast_1d(np.asarray(a_extra, complex))
    b_extra = np.atleast_1d(np.asarray(b_extra, complex))
    k, e = len(U), len(a_extra)
    W = np.eye(k + e, dtype=complex)
    W[:k, :k] = U
    if U_extra is not None:
        W[k:, k:] = U_extra
    return (np.concatenate([np.asarray(a, complex), a_extra]),
            np.concatenate([np.asarray(b, complex), b_extra]), W)


def perm_matrix(sigma):
    """P with P[i, sigma(i)] = 1, so det(D_a + P D_b P*) = prod a_i + b_sigma(i)."""
    n = len(sigma)
    P = np.zeros((n, n))
    P[np.arange(n), list(sigma)] = 1
    return P
