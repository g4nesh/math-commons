"""Stress test (floating point): the extremal unitary V of exact_ustar.py as the relative eigenbasis in dimension five.

 (1) Marcus-de Oliveira: det(D_a + V D_b V^*) against the hull of the 120 products, for complex spectra a, b.
 (2) PMSH(5): the principal minors of V D_c V^* against the hull of the 120 permuted product vectors, for real spectra c.
Each test reports a signed depth: the largest s such that the point is an affine combination of the hull points with all
weights >= s (data scaled to size 1).  Depth >= 0 means inside the hull, depth < 0 outside.  Spectra are drawn at random and
then chosen adversarially by Nelder-Mead minimisation of the depth.  This is a test, not a proof."""
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
def depth(points, z):
    """max s with z = sum w_i p_i, sum w_i = 1, w_i >= s; minus the distance to the affine hull if z is not in it"""
    sc = max(np.abs(points).max(), 1e-300); m, d = points.shape
    u, sv, vt = np.linalg.svd((points - points.mean(0)) / sc, full_matrices=False); r = int((sv > 1e-9 * sv[0]).sum()); B = vt[:r]       # coordinates in the affine hull
    y = ((points - points.mean(0)) / sc) @ B.T; yz = ((z - points.mean(0)) / sc) @ B.T; off = np.abs(((z - points.mean(0)) / sc) - yz @ B).max()
    if off > 1e-7: return -off
    res = linprog(np.r_[np.zeros(m), -1.0], A_ub=np.hstack([-np.eye(m), np.ones((m, 1))]), b_ub=np.zeros(m), A_eq=np.vstack([np.hstack([y.T, np.zeros((r, 1))]), np.r_[np.ones(m), 0][None]]),
                  b_eq=np.r_[yz, 1.0], bounds=[(None, None)] * (m + 1), method='highs')
    return -res.fun if res.status == 0 else -1.0
def moc_depth(a, b, U=V):
    z = np.linalg.det(np.diag(a) + (U * b) @ U.conj().T); pts = np.prod(a[None, :] + b[PERMS], axis=1)
    return depth(np.column_stack([pts.real, pts.imag]), np.array([z.real, z.imag]))
def pmsh_depth(c, U=V):
    H = (U * c) @ U.conj().T; mvec = np.array([np.linalg.det(H[np.ix_(T, T)]).real for T in SUBS])
    return depth(np.array([[np.prod(c[p[list(T)]]) for T in SUBS] for p in PERMS]), mvec)
def cplx(): return rng.standard_normal(5) + 1j * rng.standard_normal(5)
d1 = min(moc_depth(cplx(), cplx()) for _ in range(NT)); print("(1) Marcus-de Oliveira at V: %d random complex spectra: smallest depth %.2e" % (NT, d1))
runs = max(NT // 10, 5); worst = np.inf
for _ in range(runs):
    r = minimize(lambda x: moc_depth(x[:5] + 1j * x[5:10], x[10:15] + 1j * x[15:20]), rng.standard_normal(20), method='Nelder-Mead', options={'maxiter': 3000, 'xatol': 1e-9, 'fatol': 1e-13}); worst = min(worst, r.fun)
print("    %d Nelder-Mead minimisations of the depth over (a, b): smallest depth reached %.2e" % (runs, worst))
d2 = min(pmsh_depth(np.sort(rng.standard_normal(5))) for _ in range(NT // 4)); print("(2) PMSH(5) at V: %d random real spectra: smallest depth %.2e" % (NT // 4, d2))
runs = max(NT // 20, 5); worst = np.inf
for _ in range(runs):
    r = minimize(pmsh_depth, rng.standard_normal(5), method='Nelder-Mead', options={'maxiter': 1500, 'xatol': 1e-9, 'fatol': 1e-13}); worst = min(worst, r.fun)
print("    %d Nelder-Mead minimisations of the depth over c: smallest depth reached %.2e" % (runs, worst))
