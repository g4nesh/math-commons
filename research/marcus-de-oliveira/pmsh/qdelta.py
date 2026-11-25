"""Q_delta = Q_triv cut by g_F >= -delta for all 48 exotic facets F (g in 'vertex units': values 0,1,2 on vertices).
Find delta_crit(a) = sup{delta : Phi_a(Q_delta) = Phi_a(Q)}.  If delta_crit(a) > 0.0674422 for all a, PMSH(4) reduces to
the single bound  min_{U in U(4)} g(U) >= -0.0674422... (Drury's constant)."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
from qtriv import E, e, VX, NX, Phi_matrix, xvec_U, perms
# exotic facets of Q in x-coordinates
cen = VX.mean(0); X = VX - cen
u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); Bq = vt[:r]
h = ConvexHull(X @ Bq.T)
eqr = np.round(h.equations / np.abs(h.equations[:, :-1]).max(1, keepdims=True), 8)
_, idx = np.unique(eqr, axis=0, return_index=True); FQ = h.equations[idx]
G = []   # exotic facet functionals as (w, w0): g(x) = w.x + w0, with g = 0 on facet, values 1,2 on other vertices
for f in FQ:
    vals = (X @ Bq.T) @ f[:-1] + f[-1]
    on = np.abs(vals) < 1e-9
    if on.sum() == 13:
        unit = np.abs(vals[~on]).min()
        w = -(f[:-1] @ Bq) / unit; w0 = -(f[-1] - (f[:-1] @ Bq) @ cen) / unit
        gv = VX @ w + w0; assert np.allclose(np.sort(np.unique(np.round(gv, 6))), [0, 1, 2])
        G.append((w, w0))
assert len(G) == 48
GW = np.array([g[0] for g in G]); G0 = np.array([g[1] for g in G])
# sanity: deepest unitary
def excess(a, delta):
    M = Phi_matrix(a); V = VX @ M.T; c0 = V.mean(0); Y = V - c0
    u, sv, vt = np.linalg.svd(Y, full_matrices=False); rr = (sv > 1e-9 * sv[0]).sum(); Bas = vt[:rr]
    hh = ConvexHull(Y @ Bas.T)
    eq2 = np.round(hh.equations / np.abs(hh.equations[:, :-1]).max(1, keepdims=True), 8)
    _, id2 = np.unique(eq2, axis=0, return_index=True)
    worst = 0.0
    A_ub = -GW; b_ub = G0 + delta          # g(x) >= -delta  <=>  -w.x <= w0 + delta
    for f in hh.equations[id2]:
        w = (f[:-1] @ Bas) @ M; const = f[-1] - (f[:-1] @ Bas) @ c0
        res = linprog(-w, A_ub=A_ub, b_ub=b_ub, A_eq=E, b_eq=e, bounds=[(0, None)] * NX, method='highs')
        worst = max(worst, -res.fun + const)
    return worst
def delta_crit(a, hi=1.0):
    lo = 0.0
    if excess(a, hi) <= 1e-9: return hi
    for _ in range(30):
        mid = (lo + hi) / 2
        if excess(a, mid) <= 1e-9: lo = mid
        else: hi = mid
    return lo
if __name__ == '__main__':
    import json
    Ud = np.load(sys.argv[1]) if len(sys.argv) > 1 else None
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    out = []
    tests = [np.array([0.3, 1.1, 2.0, 3.7]), np.array([1, 2, 3, 4.]), np.array([-2, 0.5, 1, 9.]), np.array([1, 1.01, 5, 9.])] + [np.sort(rng.standard_normal(4) * 3) for _ in range(8)]
    for a in tests:
        d = delta_crit(a); out.append((a.tolist(), d))
        print("a =", np.round(a, 3), " delta_crit = %.6f   (Drury depth 0.0674422; ratio %.3f)" % (d, d / 0.0674422489), flush=True)
