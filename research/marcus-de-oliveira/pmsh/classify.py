"""For each facet of Phi_a(Q): minimal total exotic weight in a Farkas representation over Q's facets, and the
maximum of its functional over (i) Q_triv and (ii) actual unitaries (numerical)."""
import sys, itertools, numpy as np
from scipy.optimize import linprog, minimize
from scipy.spatial import ConvexHull
from scipy.linalg import expm
from qtriv import E, e, VX, NX, Phi_matrix, xvec_U
from qdelta import FQ, X as XQ, Bq, cen as cenQ
# all 82 facets of Q as functionals F_k(x) = w.x + w0 <= 0 on Q (normalised: vertex values in {0,-1,-2,..})
Fs = []
for f in FQ:
    vals = (XQ @ Bq.T) @ f[:-1] + f[-1]; on = np.abs(vals) < 1e-9; unit = np.abs(vals[~on]).min()
    Fs.append(((f[:-1] @ Bq) / unit, (f[-1] - (f[:-1] @ Bq) @ cenQ) / unit, on.sum()))
FW = np.array([f[0] for f in Fs]); F0 = np.array([f[1] for f in Fs]); FK = np.array([f[2] for f in Fs])
exotic = FK == 13
iu = np.triu_indices(4, 1)
def U_of(x):
    H = np.zeros((4, 4), complex); H[iu] = x[:6] + 1j * x[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j * x[12:16]; return expm(H)
if __name__ == '__main__':
    a = np.array([float(t) for t in sys.argv[1].split(',')]) if len(sys.argv) > 1 else np.array([1, 2, 3, 4.])
    M = Phi_matrix(a); V = VX @ M.T; c0 = V.mean(0); Y = V - c0
    u, sv, vt = np.linalg.svd(Y, full_matrices=False); rr = (sv > 1e-9 * sv[0]).sum(); Bas = vt[:rr]
    hh = ConvexHull(Y @ Bas.T)
    eq2 = np.round(hh.equations / np.abs(hh.equations[:, :-1]).max(1, keepdims=True), 8)
    _, id2 = np.unique(eq2, axis=0, return_index=True)
    rng = np.random.default_rng(1)
    rows = []
    for f in hh.equations[id2]:
        w = (f[:-1] @ Bas) @ M; const = f[-1] - (f[:-1] @ Bas) @ c0
        nver = int((np.abs(Y @ Bas.T @ f[:-1] + f[-1]) < 1e-9).sum())
        # Farkas: w = sum_k alpha_k FW_k + E^T mu ; const = sum alpha_k F0_k + ... ; minimise exotic weight
        # (w.x + const <= 0 on aff Q) as nonneg combination of facets plus affine equations: solve LP in (alpha>=0, mu free, t)
        nF = len(Fs); nE = E.shape[0]
        A_eq = np.hstack([FW.T, E.T, np.zeros((NX, 1))])          # sum alpha FW + E^T mu = w
        b_eq = w
        A_eq2 = np.hstack([F0[None, :], -e[None, :], -np.ones((1, 1))])  # sum alpha F0 - e.mu - s = const, s >= 0 slack
        cobj = np.r_[exotic.astype(float), np.zeros(nE), 0.0]
        res = linprog(cobj, A_eq=np.vstack([A_eq, A_eq2]), b_eq=np.r_[b_eq, const],
                      bounds=[(0, None)] * nF + [(None, None)] * nE + [(0, None)], method='highs')
        exw = res.fun if res.status == 0 else np.nan
        triv = linprog(-w, A_eq=E, b_eq=e, bounds=[(0, None)] * NX, method='highs'); exc_triv = -triv.fun + const
        rows.append((nver, exw, exc_triv, w, const))
    rows.sort(key=lambda t: -t[1])
    nontriv = [r for r in rows if r[1] > 1e-7]
    print("a =", a, ": projected facets", len(rows), "; needing exotic facets:", len(nontriv),
          "; vertex counts of those:", sorted(set(r[0] for r in nontriv)))
    # for nontrivial facets: numerical max over unitaries of the facet functional (should be <= 0 if PMSH true)
    worst = -np.inf
    for k, (nver, exw, exc, w, const) in enumerate(nontriv):
        best = -np.inf
        for _ in range(25):
            res = minimize(lambda x: -(w @ xvec_U(U_of(x)) + const), rng.uniform(-np.pi, np.pi, 16), method='BFGS')
            best = max(best, -res.fun)
        worst = max(worst, best)
        if k < 8 or best > -1e-6:
            print("  facet %2d: vertices %d, min exotic weight %.4f, excess over Q_triv %.4f, max over U(4) %.3e" % (k, nver, exw, exc, best))
    print("max over nontrivial facets of (max over U(4) of functional): %.3e   (<= 0 means PMSH holds at this a numerically)" % worst)
