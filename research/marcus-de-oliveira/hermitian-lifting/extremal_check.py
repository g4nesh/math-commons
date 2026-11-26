"""Numerical check of Lemma E and of Steps 2-5 of README.md at critical points that are NOT permutation configurations.

Why not maximisers.  If the inclusion holds, a generic linear functional is maximised on the determinantal range at one of
the n! products, where B' commutes with A and every identity of Lemma E is trivial.  So this script looks for arbitrary
critical points of U -> l(det(A + U B U^*)): it solves the criticality equation "N = conj(omega) [B', adj M] is Hermitian"
by a Levenberg-Marquardt root search from Haar-random starts, and keeps the solutions at which [B', adj M] is not zero.

A = diag(a), B normal with one real eigenvalue beta and n-1 nonreal ones, l(z) = Re(conj(omega) z).
Two modes: A + beta positive definite, and A + beta indefinite.  At each critical point kept it reports
  (i)   |N - N^*| ;
  (ii)  |(A + beta) w - 2 l(z0) v|,  w = conj(omega) adj(M) v + omega adj(M)^* v ;
  (iii) |p l(D1) - q l(z0)|,  p = det(A + beta), q = det(A_1 + beta), D1 = det of the compressed pair ;
  (iv)  in the definite mode: the distance of (p/q) D1 from the hull K of the n! products (linear programme), and
        l(z0) - max_K l.  In the indefinite mode (p/q) D1 may leave K; the script counts how often.
Floating point.  usage: extremal_check.py n [seed] [instances] [starts per instance]"""
import itertools, sys, numpy as np
from scipy.optimize import least_squares, linprog
from scipy.linalg import expm
n = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
NI = int(sys.argv[3]) if len(sys.argv) > 3 else 6; NS = int(sys.argv[4]) if len(sys.argv) > 4 else 12
IU = np.triu_indices(n, 1); m2 = len(IU[0]); PERMS = np.array(list(itertools.permutations(range(n))))
def herm(x):
    H = np.zeros((n, n), complex); H[IU] = x[:m2] + 1j * x[m2:2 * m2]; H = H + H.conj().T; H[np.diag_indices(n)] = x[2 * m2:]; return H
def haar():
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
def adj(M):
    u, s, vh = np.linalg.svd(M); cof = np.array([np.prod(np.delete(s, i)) for i in range(n)])          # works for singular M too
    return np.linalg.det(u) * np.linalg.det(vh) * (vh.conj().T * cof) @ u.conj().T
def hull_distance(pts, z):
    sc = np.abs(pts).max(); m = len(pts); P = np.column_stack([pts.real, pts.imag]) / sc; zz = np.array([z.real, z.imag]) / sc
    A_ub = np.vstack([np.hstack([P.T, -np.ones((2, 1))]), np.hstack([-P.T, -np.ones((2, 1))])])
    r = linprog(np.r_[np.zeros(m), 1.0], A_ub=A_ub, b_ub=np.r_[zz, -zz], A_eq=np.r_[np.ones(m), 0][None], b_eq=[1.0], bounds=[(0, None)] * (m + 1), method='highs'); return r.fun
for mode in ('definite', 'indefinite'):
    err = np.zeros(4); kept = 0; found = 0; outside = 0; worst_dist = 0.0; worst_gap = -np.inf
    for inst in range(NI):
        a = np.sort(rng.standard_normal(n)); beta = (-a[0] + rng.uniform(0.2, 1.5)) if mode == 'definite' else -rng.uniform(a[0] + 0.1 * (a[-1] - a[0]), a[-1] - 0.1 * (a[-1] - a[0]))
        b = np.r_[rng.standard_normal(n - 1) + 1j * rng.standard_normal(n - 1), beta]; A = np.diag(a).astype(complex); om = np.exp(1j * rng.uniform(0, 2 * np.pi))
        ell = lambda z: (np.conj(om) * z).real; pts = np.prod(a[None, :] + b[PERMS], axis=1); kmax = ell(pts).max(); scale = np.abs(pts).max()
        for s in range(NS):
            U0 = haar()
            def resid(x):
                U = U0 @ expm(1j * herm(x)); Bp = (U * b) @ U.conj().T; M = A + Bp; N = np.conj(om) * (Bp @ adj(M) - adj(M) @ Bp); K = (N - N.conj().T) / scale
                return np.r_[K[IU].real, K[IU].imag, np.diag(K).imag]
            sol = least_squares(resid, np.zeros(n * n), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
            if np.abs(sol.fun).max() > 1e-11: continue
            found += 1; U = U0 @ expm(1j * herm(sol.x)); Bp = (U * b) @ U.conj().T; M = A + Bp; aM = adj(M)
            if np.abs(Bp @ aM - aM @ Bp).max() / scale < 1e-6: continue                              # a permutation-type point: trivial
            kept += 1; z0 = np.linalg.det(M); v = U[:, n - 1]; N = np.conj(om) * (Bp @ aM - aM @ Bp)
            w = np.conj(om) * (aM @ v) + om * (aM.conj().T @ v)
            Q = np.linalg.qr(np.column_stack([v, rng.standard_normal((n, n - 1)) + 1j * rng.standard_normal((n, n - 1))]))[0][:, 1:]
            A1 = Q.conj().T @ A @ Q; D1 = np.linalg.det(Q.conj().T @ M @ Q); p = np.prod(a + beta); q = np.linalg.det(A1 + beta * np.eye(n - 1)).real
            err = np.maximum(err, [np.abs(N - N.conj().T).max() / scale, np.abs((A + beta * np.eye(n)) @ w - 2 * ell(z0) * v).max() / scale,
                                   abs(p * ell(D1) - q * ell(z0)) / (scale * max(1, abs(p))), abs(D1 - v.conj() @ aM @ v) / scale])
            if abs(q) > 1e-9:
                dist = hull_distance(pts, (p / q) * D1)
                if mode == 'definite': worst_dist = max(worst_dist, dist); worst_gap = max(worst_gap, (ell(z0) - kmax) / scale)
                else: outside += dist > 1e-9
    print("n=%d, A + beta %s: %d critical points found, %d not of permutation type" % (n, mode, found, kept))
    print("   largest errors there: (i) N Hermitian %.1e ; (ii) %.1e ; (iii) %.1e ; cofactor identity %.1e" % tuple(err))
    if mode == 'definite': print("   (iv) largest distance of (p/q) D1 from K: %.1e ; largest l(z0) - max_K l: %.1e (relative)" % (worst_dist, worst_gap))
    else: print("   (p/q) D1 lies outside K at %d of these points (Lemma I does not apply here; this is why Lemma M is needed)" % outside)
