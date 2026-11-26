"""Squared minors of O in SO(4) and det(D_a + O D_b O^T) in rotation-pair coordinates (README.md, remark on the general case).

Lambda^2 O = A (+) B on self-dual / anti-self-dual 2-forms, A, B in SO(3) (basis as in ../star-orthogonal/README.md).
With x = A o A + B o B and y = A o B (entrywise products), numerically:
  (1) O_ij^2 = (1 + sum_{M,N} eps[i,M] eps[j,N] y[M,N]) / 4 ;
  (2) det O[S_M,S_N]^2 = det O[S_M^c,S_N^c]^2 = (x + 2y)[M,N] / 4 ;  det O[S_M,S_N^c]^2 = det O[S_M^c,S_N]^2 = (x - 2y)[M,N] / 4 ;
  (3) det(D_a + O D_b O^T) = e4(a) + e4(b) + (e3(a) e1(b) + e1(a) e3(b))/4 + p^T x p'/4 + sum beta[M,N] y[M,N]."""
import itertools, numpy as np
rng = np.random.default_rng(1)
P2 = list(itertools.combinations(range(4), 2)); IDX = {S: k for k, S in enumerate(P2)}
SM = [(0, 1), (0, 2), (0, 3)]; SMc = [(2, 3), (1, 3), (1, 2)]; sg = [1, -1, 1]; r2 = np.sqrt(2)
B6 = np.zeros((6, 6))          # columns: sigma_1..3 = (e_S + sg e_S^c)/sqrt2, tau_1..3 = (e_S - sg e_S^c)/sqrt2
for M in range(3):
    B6[IDX[SM[M]], M] = 1 / r2; B6[IDX[SMc[M]], M] = sg[M] / r2; B6[IDX[SM[M]], 3 + M] = 1 / r2; B6[IDX[SMc[M]], 3 + M] = -sg[M] / r2
eps = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
def so4():
    q, r = np.linalg.qr(rng.standard_normal((4, 4))); q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0: q[:, 0] *= -1
    return q
def compound(O): return np.array([[np.linalg.det(O[np.ix_(S, T)]) for T in P2] for S in P2])
def e(a, k): return sum(np.prod([a[i] for i in S]) for S in itertools.combinations(range(4), k))
def coeffs(a):
    p = np.array([a[S[0]] * a[S[1]] + a[T[0]] * a[T[1]] for S, T in zip(SM, SMc)]); q = np.array([a[S[0]] * a[S[1]] - a[T[0]] * a[T[1]] for S, T in zip(SM, SMc)])
    return p, q, eps.T @ np.array([np.prod(np.delete(a, i)) for i in range(4)]), eps.T @ a
w = np.zeros(4)
for _ in range(500):
    O = so4(); C = compound(O); G = B6.T @ C @ B6; A, B = G[:3, :3], G[3:, 3:]; x = A * A + B * B; y = A * B
    w[0] = max(w[0], np.abs(G[:3, 3:]).max(), np.abs(G[3:, :3]).max(), abs(np.linalg.det(A) - 1), abs(np.linalg.det(B) - 1))
    w[1] = max(w[1], np.abs(O**2 - (1 + eps @ y @ eps.T) / 4).max())
    for M in range(3):
        for N in range(3):
            w[2] = max(w[2], abs(C[IDX[SM[M]], IDX[SM[N]]]**2 - (x + 2 * y)[M, N] / 4), abs(C[IDX[SMc[M]], IDX[SMc[N]]]**2 - (x + 2 * y)[M, N] / 4),
                       abs(C[IDX[SM[M]], IDX[SMc[N]]]**2 - (x - 2 * y)[M, N] / 4), abs(C[IDX[SMc[M]], IDX[SM[N]]]**2 - (x - 2 * y)[M, N] / 4))
    a = rng.standard_normal(4) + 1j * rng.standard_normal(4); b = rng.standard_normal(4) + 1j * rng.standard_normal(4)
    p, q, s, t = coeffs(a); p2, q2, s2, t2 = coeffs(b); beta = (np.outer(s, t2) + np.outer(t, s2)) / 4 - np.outer(q, q2) / 2
    val = e(a, 4) + e(b, 4) + (e(a, 3) * e(b, 1) + e(a, 1) * e(b, 3)) / 4 + (p @ x @ p2) / 4 + (beta * y).sum()
    w[3] = max(w[3], abs(val - np.linalg.det(np.diag(a) + O @ np.diag(b) @ O.T)))
print("500 random O in SO(4), complex a, b. Largest errors:")
print("   Lambda^2 O = A (+) B with A, B in SO(3): %.1e ; (1) squared entries: %.1e ; (2) squared 2x2 minors: %.1e ; (3) determinant: %.1e" % tuple(w))
print("ALL CHECKS PASSED" if w.max() < 1e-10 else "A CHECK FAILED")
