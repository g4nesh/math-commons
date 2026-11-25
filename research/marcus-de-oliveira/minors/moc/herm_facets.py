"""Facets of the projected polytope Phi_a(Q) for n=4 and a fixed real a; express each facet as a set function y
on column sets T (functional  sum_T y_T phi_T ), and check whether supp(y) is a chain."""
import sys, itertools, numpy as np
from scipy.spatial import ConvexHull
sys.path.insert(0, 'moc')
from herm_cert import HermCert, Phi, all_minor_sq
n = 4; hc = HermCert(n)
Ts = [T for k in range(1, n) for T in itertools.combinations(range(n), k)]
def analyse(a, verbose=True):
    V = np.array([Phi(pd, a, n) for pd in hc.Pd]); c = V.mean(0); X = V - c
    u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); B = vt[:r]
    h = ConvexHull(X @ B.T)
    eq = h.equations; eqr = np.round(eq / np.abs(eq[:, :-1]).max(1, keepdims=True), 7)
    _, idx = np.unique(eqr, axis=0, return_index=True); F = eq[idx]
    res = []
    for f in F:
        y = f[:-1] @ B           # functional on Phi-coordinates (one representative, modulo span relations)
        on = np.abs(X @ B.T @ f[:-1] + f[-1]) < 1e-8
        res.append((y, on))
    return r, res
def chain_representable(y_rel, tol=1e-9):
    # does there exist a representative y (mod relations killing all of aff Q) supported on a chain?  (LP-free heuristic:
    # report support of the minimum-L1 representative via linear programming)
    return None
a = np.array([0.3, 1.1, 2.0, 3.7])
r, res = analyse(a)
print("real a =", a, ": projected dim", r, ", facets", len(res))
cnt = {}
for y, on in res: cnt[on.sum()] = cnt.get(on.sum(), 0) + 1
print("vertices-per-facet histogram:", dict(sorted(cnt.items())))
np.save('moc/runs/herm_facets_example.npy', np.array([np.r_[y, on] for y, on in res]))
