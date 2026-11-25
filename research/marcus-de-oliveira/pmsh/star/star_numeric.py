"""Numerical check of the six spectrum-free inequalities attached to g0 = B11 + B33 + D[14,34] - D[12,14]:
   g0 + D_a + 2 sqrt(D_a D_b) >= 0   for (a,b) in the list below (1-based row sets , column sets)."""
import sys, numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
from base import *
PAIRS = [(((0,1),(0,1)), ((0,2),(0,3))), (((0,1),(0,1)), ((0,2),(1,2))),      # arc A
         (((0,3),(0,3)), ((0,1),(0,2))), (((0,3),(0,3)), ((0,1),(1,3))),      # arc B
         (((0,2),(0,2)), ((0,3),(0,1))), (((0,2),(0,2)), ((0,3),(2,3)))]      # arc C
g0 = np.zeros(NX); g0[iB(0,0)] += 1; g0[iB(2,2)] += 1; g0[iD((0,3),(2,3))] += 1; g0[iD((0,1),(0,3))] -= 1
iu = np.triu_indices(4, 1)
def U_of(v):
    H = np.zeros((4,4), complex); H[iu] = v[:6] + 1j*v[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j*v[12:16]; return expm(H)
def star(x, a, b): return g0 @ x + x[iD(*a)] + 2*np.sqrt(max(x[iD(*a)], 0) * max(x[iD(*b)], 0))
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); R = int(sys.argv[2]) if len(sys.argv) > 2 else 60
nm = lambda S, T: "D[%s,%s]" % ("".join(str(i+1) for i in S), "".join(str(i+1) for i in T))
for a, b in PAIRS:
    best = np.inf; bx = None
    for r in range(R):
        res = minimize(lambda v: star(x_U(U_of(v)), a, b), rng.uniform(-np.pi, np.pi, 16), method='Nelder-Mead', options={'maxiter': 6000, 'xatol': 1e-9, 'fatol': 1e-13})
        if res.fun < best: best, bx = res.fun, x_U(U_of(res.x))
    print("min over U(4) of g0 + %s + 2 sqrt(%s * %s) = %+.3e   (at the minimiser: g0 = %+.4f, %s = %.4f, %s = %.4f)" %
          (nm(*a), nm(*a), nm(*b), best, g0 @ bx, nm(*a), bx[iD(*a)], nm(*b), bx[iD(*b)]), flush=True)
