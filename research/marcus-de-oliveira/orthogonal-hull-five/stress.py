"""Stress test (floating point): the extremal unitary V of exact_ustar.py as the relative eigenbasis in dimension five.

 (1) Marcus-de Oliveira: is det(D_a + V D_b V^*) in the hull of the 120 products, for random complex spectra a, b, and for
     spectra chosen adversarially by local optimisation of the distance to the hull?
 (2) PMSH(5): are the principal minors of V D_c V^* in the hull of the 120 permuted product vectors, for real spectra c?
The numbers printed are distances from the hull relative to the size of the data (0 = inside).  This is a test, not a proof."""
import itertools, sys, numpy as np
from scipy.optimize import linprog, minimize
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); NT = int(sys.argv[2]) if len(sys.argv) > 2 else 400
s6 = np.sqrt(6); t = 1j * np.sqrt(15)
V = np.array([[4, 4 * s6, 4 * s6, 4 * s6, 4 * s6],
              [4 * s6, -1 - t, -6 - 2 * t, -6 + 2 * t, 9 + t],
              [4 * s6, -6 - 2 * t, -1 + t, 9 - t, -6 + 2 * t],
              [4 * s6, -6 + 2 * t, 9 - t, -1 + t, -6 - 2 * t],
              [4 * s6, 9 + t, -6 + 2 * t, -6 - 2 * t, -1 - t]]) / 20
assert np.abs(V.conj().T @ V - np.eye(5)).max() < 1e-14
PERMS = np.array(list(itertools.permutations(range(5)))); SUBS = [T for k in range(1, 5) for T in itertools.combinations(range(5), k)]
def lp_dist(points, z):
    """smallest max-norm distance from z to the convex hull of the rows of `points` (real vectors), relative to scale"""
    sc = max(np.abs(points).max(), 1e-300); d = points.shape[1]; m = points.shape[0]
    A = np.vstack([np.hstack([points.T / sc, -np.ones((d, 1))]), np.hstack([-points.T / sc, -np.ones((d, 1))])])
    r = linprog(np.r_[np.zeros(m), 1.0], A_ub=A, b_ub=np.r_[z / sc, -z / sc], A_eq=np.r_[np.ones(m), 0][None], b_eq=[1.0], bounds=[(0, None)] * (m + 1), method='highs'); return r.fun
def moc_dist(a, b, U=V):
    z = np.linalg.det(np.diag(a) + (U * b) @ U.conj().T); pts = np.prod(a[None, :] + b[PERMS], axis=1)
    return lp_dist(np.column_stack([pts.real, pts.imag]), np.array([z.real, z.imag]))
def pmsh_dist(c, U=V):
    H = (U * c) @ U.conj().T; m = np.array([np.linalg.det(H[np.ix_(T, T)]).real for T in SUBS])
    return lp_dist(np.array([[np.prod(c[p[list(T)]]) for T in SUBS] for p in PERMS]), m)
def cplx(): return rng.standard_normal(5) + 1j * rng.standard_normal(5)
w1 = max(moc_dist(cplx(), cplx()) for _ in range(NT)); print("(1) Marcus-de Oliveira at V, %d random complex spectra: largest distance from the hull %.1e" % (NT, w1))
worst = 0
for _ in range(max(NT // 20, 5)):
    x0 = rng.standard_normal(20); f = lambda x: -moc_dist(x[:5] + 1j * x[5:10], x[10:15] + 1j * x[15:20])
    r = minimize(f, x0, method='Nelder-Mead', options={'maxiter': 1500, 'xatol': 1e-7, 'fatol': 1e-12}); worst = max(worst, -r.fun)
print("    adversarial spectra (%d Nelder-Mead runs): largest distance %.1e" % (max(NT // 20, 5), worst))
w2 = max(pmsh_dist(np.sort(rng.standard_normal(5))) for _ in range(NT // 4)); print("(2) PMSH(5) at V, %d random real spectra: largest distance from the hull %.1e" % (NT // 4, w2))
worst = 0
for _ in range(max(NT // 40, 5)):
    r = minimize(lambda x: -pmsh_dist(x), rng.standard_normal(5), method='Nelder-Mead', options={'maxiter': 600, 'xatol': 1e-7, 'fatol': 1e-12}); worst = max(worst, -r.fun)
print("    adversarial real spectra (%d Nelder-Mead runs): largest distance %.1e" % (max(NT // 40, 5), worst))
