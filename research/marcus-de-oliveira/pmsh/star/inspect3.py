"""Look for a simple degree-3 certificate of  F1 = g0 + 2 D[12,12] + D[13,14] >= 0  (the t = 1 member):
minimise the trace of the Gram matrix subject to the exact identity, then list the active blocks."""
import numpy as np, cvxpy as cp, scipy.sparse as sps
from sos2 import *
L = Level(3, 3000)
f = L.X @ (g0 + 2 * a_ + b_)
g = cp.Variable(L.nvar); cons = []; tr = 0
for (lo, hi, n) in L.slices:
    iu = np.triu_indices(n)
    if n == 1: cons.append(g[lo] >= 0); tr = tr + g[lo]; continue
    S = sps.lil_matrix((n * n, hi - lo))
    for k, (p, q) in enumerate(zip(*iu)): S[p * n + q, k] = 1; S[q * n + p, k] = 1
    G = cp.reshape(sps.csr_matrix(S) @ g[lo:hi], (n, n), order='C'); cons.append(G >> 0); tr = tr + cp.trace(G)
cons.append(L.Ar @ g == L.P @ f + 1e-6 * L.ones)
prob = cp.Problem(cp.Minimize(tr), cons); prob.solve(solver='CLARABEL')
print("status", prob.status, " min trace %.6f" % prob.value)
gv = g.value; nm = lambda v: "u%d%d" % (v // 4 + 1, v % 4 + 1)
act = 0
for (lo, hi, n), b in zip(L.slices, L.blocks):
    G = np.zeros((n, n)); iu = np.triu_indices(n); G[iu] = gv[lo:hi]; G = G + G.T - np.diag(np.diag(G))
    w, V = np.linalg.eigh(G)
    if w.max() < 1e-6: continue
    act += 1
    if act <= 40:
        mons = ["*".join(nm(v) for v in L.mons[i]) for i in b]
        terms = []
        for lam, vec in zip(w, V.T):
            if lam > 1e-6: terms.append("%.4f*|%s|^2" % (lam, " ".join("%+.3f*%s" % (c, m) for c, m in zip(vec, mons) if abs(c) > 1e-4)))
        print("  ", " + ".join(terms))
print("active blocks:", act, "of", len(L.blocks))
