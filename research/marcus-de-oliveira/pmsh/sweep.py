"""Across many spectra: number of facets of Phi_a(Q), which need exotic facets (Farkas LP), their vertex counts,
and (optionally) a numerical max over U(4) of each nontrivial facet functional."""
import sys, json, numpy as np
from scipy.optimize import linprog, minimize
from scipy.spatial import ConvexHull
from scipy.linalg import expm
from qtriv import E, e, VX, NX, Phi_matrix, xvec_U
from classify import FW, F0, exotic, U_of
def analyse(a, uopt_restarts=0, rng=None):
    M = Phi_matrix(a); V = VX @ M.T; c0 = V.mean(0); Y = V - c0
    u, sv, vt = np.linalg.svd(Y, full_matrices=False); rr = (sv > 1e-9 * sv[0]).sum(); Bas = vt[:rr]
    hh = ConvexHull(Y @ Bas.T)
    eq2 = np.round(hh.equations / np.abs(hh.equations[:, :-1]).max(1, keepdims=True), 8)
    _, id2 = np.unique(eq2, axis=0, return_index=True)
    nF = FW.shape[0]; nE = E.shape[0]; stats = []; worstU = -np.inf
    for f in hh.equations[id2]:
        w = (f[:-1] @ Bas) @ M; const = f[-1] - (f[:-1] @ Bas) @ c0
        nver = int((np.abs(Y @ Bas.T @ f[:-1] + f[-1]) < 1e-9).sum())
        res = linprog(exotic.astype(float).tolist() + [0] * nE + [0], A_eq=np.vstack([np.hstack([FW.T, E.T, np.zeros((NX, 1))]), np.hstack([F0[None, :], -e[None, :], -np.ones((1, 1))])]),
                      b_eq=np.r_[w, const], bounds=[(0, None)] * nF + [(None, None)] * nE + [(0, None)], method='highs')
        exw = res.fun
        stats.append((nver, exw > 1e-7))
        if exw > 1e-7 and uopt_restarts:
            best = -np.inf
            for _ in range(uopt_restarts):
                r2 = minimize(lambda x: -(w @ xvec_U(U_of(x)) + const), rng.uniform(-np.pi, np.pi, 16), method='BFGS')
                best = max(best, -r2.fun)
            worstU = max(worstU, best)
    nontriv = [s[0] for s in stats if s[1]]
    return len(stats), len(nontriv), sorted(set(nontriv)), sorted(set(s[0] for s in stats if not s[1])), worstU
if __name__ == '__main__':
    seed, k, ur = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = np.random.default_rng(seed); out = []
    for i in range(k):
        kind = i % 3
        if kind == 0: a = np.sort(rng.standard_normal(4) * rng.choice([0.5, 2, 5]))
        elif kind == 1: a = np.sort(np.r_[rng.uniform(1, 2, 2), rng.uniform(1, 2, 2) * rng.choice([3, 30])])   # clustered pairs
        else: a = np.sort(rng.uniform(0.1, 1, 4)) ** rng.choice([1, 3, 6])                                      # skewed ratios
        nf, nt, vt_nt, vt_tr, wU = analyse(a, ur, rng)
        out.append({'a': a.tolist(), 'facets': nf, 'nontrivial': nt, 'nontrivial_vertex_counts': vt_nt, 'trivial_vertex_counts': vt_tr, 'max_over_U': wU})
        print("a=%s facets=%d nontrivial=%d (vertex counts %s; trivial ones %s) max_U(nontrivial)=%.2e" % (np.round(a, 3), nf, nt, vt_nt, vt_tr, wU), flush=True)
    json.dump(out, open('sweep_s%d.json' % seed, 'w'))
