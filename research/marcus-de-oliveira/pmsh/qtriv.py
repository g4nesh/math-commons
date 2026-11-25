"""PMSH(4) via trivial facets.  x = (B: 16 entries b_{s,j}; D: 36 entries D_{S,T}, |S|=|T|=2).
Q_triv := {x : identities (i)-(iv) hold, B >= 0, D >= 0}.  Q := conv of permutation points.
Test: does Phi_a(Q_triv) = Phi_a(Q)  (max over Q_triv of each facet functional of Phi_a(Q) equals its max over Q)?"""
import itertools, numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
n = 4
P2 = list(itertools.combinations(range(n), 2))
iB = lambda s, j: s * 4 + j
iD = lambda S, T: 16 + P2.index(S) * 6 + P2.index(T)
NX = 52
def comp(S): return tuple(sorted(set(range(n)) - set(S)))
def equations():
    rows, rhs = [], []
    def add(coeffs, r):
        v = np.zeros(NX)
        for k, c in coeffs: v[k] += c
        rows.append(v); rhs.append(r)
    for s in range(n): add([(iB(s, j), 1) for j in range(n)], 1)           # (i)
    for j in range(n): add([(iB(s, j), 1) for s in range(n)], 1)
    for S in P2: add([(iD(S, T), 1) for T in P2], 1)                         # (ii)
    for T in P2: add([(iD(S, T), 1) for S in P2], 1)
    for S in P2:                                                             # (iii)
        for j in range(n): add([(iD(S, T), 1) for T in P2 if j in T] + [(iB(s, j), -1) for s in S], 0)
    for T in P2:
        for i in range(n): add([(iD(S, T), 1) for S in P2 if i in S] + [(iB(i, t), -1) for t in T], 0)
    for S in P2:                                                             # (iv) Jacobi complement
        for T in P2: add([(iD(S, T), 1), (iD(comp(S), comp(T)), -1)], 0)
    return np.array(rows), np.array(rhs, float)
E, e = equations()
perms = list(itertools.permutations(range(n)))
def xvec_perm(p):
    x = np.zeros(NX)
    for s in range(n): x[iB(s, p[s])] = 1
    for S in P2: x[iD(S, tuple(sorted(p[s] for s in S)))] = 1
    return x
def xvec_U(U):
    x = np.zeros(NX)
    for s in range(n):
        for j in range(n): x[iB(s, j)] = abs(U[s, j])**2
    for S in P2:
        for T in P2: x[iD(S, T)] = abs(np.linalg.det(U[np.ix_(S, T)]))**2
    return x
VX = np.array([xvec_perm(p) for p in perms])
assert np.abs(VX @ E.T - e).max() < 1e-12
# dimension of the solution space of (i)-(iv)
rank = np.linalg.matrix_rank(E); aff_dim_eq = NX - rank
aff_dim_Q = np.linalg.matrix_rank(VX - VX[0])
def Phi_matrix(a):
    Ts = [T for k in (1, 2, 3) for T in itertools.combinations(range(n), k)]
    M = np.zeros((len(Ts), NX))
    for r, T in enumerate(Ts):
        if len(T) == 1:
            t = T[0]
            for s in range(n): M[r, iB(s, t)] = np.prod([a[i] for i in range(n) if i != s])
        elif len(T) == 2:
            for S in P2: M[r, iD(S, T)] = np.prod([a[i] for i in comp(S)])
        else:
            t = comp(T)[0]
            for s in range(n): M[r, iB(s, t)] = a[s]
    return M
def test(a, verbose=False):
    M = Phi_matrix(a); V = VX @ M.T; cen = V.mean(0); X = V - cen
    u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); Bas = vt[:r]
    h = ConvexHull(X @ Bas.T); worst = 0.0; nfac = 0
    eqr = np.round(h.equations / np.abs(h.equations[:, :-1]).max(1, keepdims=True), 8)
    _, idx = np.unique(eqr, axis=0, return_index=True)
    for f in h.equations[idx]:
        nfac += 1
        w = (f[:-1] @ Bas) @ M                     # functional on x
        const = f[-1] - (f[:-1] @ Bas) @ cen       # facet: w.x + const <= 0 on Q
        res = linprog(-w, A_eq=E, b_eq=e, bounds=[(0, None)] * NX, method='highs')
        assert res.status == 0
        val = -res.fun + const                       # max over Q_triv of (w.x + const)
        worst = max(worst, val)
    return worst, nfac, r
if __name__ == '__main__':
    print("rank of identities:", rank, "-> affine dim of solution set:", aff_dim_eq, "; affine dim of Q:", aff_dim_Q)
    rng = np.random.default_rng(0)
    for a in [np.array([0.3, 1.1, 2.0, 3.7]), np.array([1, 2, 3, 4.]), np.array([1, 1.01, 5, 9.]), np.array([-2, 0.5, 1, 9.])] + [np.sort(rng.standard_normal(4) * 3) for _ in range(6)]:
        w, nf, r = test(a)
        print("a =", np.round(a, 3), "proj dim", r, "facets", nf, "| max excess of Q_triv over Phi_a(Q) facets: %.3e" % w)
