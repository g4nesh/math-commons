"""Exact (integer) check of two rewritings of the base exotic facet on aff(Q), and numerical probes.
  g0 = B11 + B33 + D[14,34] - D[12,14]
     = D[12,12] + D[14,12] + D[14,14] + D[14,34] - D[13,24]                       (form 1)
     = D[13,13] + D[14,14] + D[13,14] + D[13,23] - D[12,34]                       (form 2)
So  g0 + D[12,12] + 2 sqrt(D[12,12] D[13,14]) >= 0   is
    D[12,34] <= D[13,13] + D[14,14] + D[13,23] + (sqrt D[12,12] + sqrt D[13,14])^2 ."""
import numpy as np, sympy as sp
from scipy.optimize import minimize
from scipy.linalg import expm
from base import *
def W(terms):
    w = np.zeros(NX, dtype=int)
    for c, kind, a, b in terms: w[iB(a, b) if kind == 'B' else iD(a, b)] += c
    return w
g0 = W([(1,'B',0,0),(1,'B',2,2),(1,'D',(0,3),(2,3)),(-1,'D',(0,1),(0,3))])
f1 = W([(1,'D',(0,1),(0,1)),(1,'D',(0,3),(0,1)),(1,'D',(0,3),(0,3)),(1,'D',(0,3),(2,3)),(-1,'D',(0,2),(1,3))])
f2 = W([(1,'D',(0,2),(0,2)),(1,'D',(0,3),(0,3)),(1,'D',(0,2),(0,3)),(1,'D',(0,2),(1,2)),(-1,'D',(0,1),(2,3))])
# two affine functionals agree on aff(Q) iff they agree on the 24 vertices (which affinely span it)
print("form 1 equals g0 on all 24 vertices:", bool((VX @ g0 == VX @ f1).all()), "; form 2:", bool((VX @ g0 == VX @ f2).all()))
iu = np.triu_indices(4, 1)
def U_of(v):
    H = np.zeros((4,4), complex); H[iu] = v[:6] + 1j*v[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j*v[12:16]; return expm(H)
rng = np.random.default_rng(5); err = 0
for _ in range(200):
    x = x_U(haar(rng)); err = max(err, abs(g0 @ x - f1 @ x), abs(g0 @ x - f2 @ x))
print("max difference on 200 Haar unitaries: %.1e" % err)
def probe(fun, label, R=40):
    best = np.inf
    for r in range(R):
        res = minimize(lambda v: fun(x_U(U_of(v))), rng.uniform(-np.pi, np.pi, 16), method='BFGS', options={'gtol': 1e-10})
        best = min(best, res.fun)
    print("min over U(4) of %-46s = %+.6f" % (label, best), flush=True)
probe(lambda x: g0 @ x, "g0")
probe(lambda x: g0 @ x + x[iD((0,1),(0,1))], "g0 + D[12,12]")
probe(lambda x: g0 @ x + x[iD((0,1),(0,1))] + x[iD((0,2),(0,3))], "g0 + D[12,12] + D[13,14]")
probe(lambda x: g0 @ x + 2 * x[iD((0,1),(0,1))], "g0 + 2 D[12,12]")
