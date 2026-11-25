"""Parameter-free strengthening of (*):  is
     L(U) := D[12,34] - D[13,13] - D[14,14] - D[13,23]  <=  | alpha * U[12,12] + beta * U[13,14] |^2
for some fixed polynomial weights alpha, beta with |alpha|, |beta| <= 1 on U(4)?  (Then (*) follows by the triangle inequality.)
Torus weights force  wt(alpha) - wt(beta) = (e3 - e2 ; e4 - e2).  We test entry weights and both signs."""
import sys, numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
iu = np.triu_indices(4, 1)
def U_of(v):
    H = np.zeros((4,4), complex); H[iu] = v[:6] + 1j*v[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j*v[12:16]
    U = expm(H); return U / np.linalg.det(U) ** 0.25
m = lambda U, S, T: np.linalg.det(U[np.ix_([s-1 for s in S], [t-1 for t in T])])
def L(U): return abs(m(U,(1,2),(3,4)))**2 - abs(m(U,(1,3),(1,3)))**2 - abs(m(U,(1,4),(1,4)))**2 - abs(m(U,(1,3),(2,3)))**2
cands = {
 "u34*X + u22*Y": lambda U: U[2,3]*m(U,(1,2),(1,2)) + U[1,1]*m(U,(1,3),(1,4)),
 "u34*X - u22*Y": lambda U: U[2,3]*m(U,(1,2),(1,2)) - U[1,1]*m(U,(1,3),(1,4)),
 "conj(u22)*X + conj(u34)*Y": lambda U: np.conj(U[1,1])*m(U,(1,2),(1,2)) + np.conj(U[2,3])*m(U,(1,3),(1,4)),
 "conj(u22)*X - conj(u34)*Y": lambda U: np.conj(U[1,1])*m(U,(1,2),(1,2)) - np.conj(U[2,3])*m(U,(1,3),(1,4)),
}
rng = np.random.default_rng(1)
for name, f in cands.items():
    best = np.inf
    for r in range(30):
        res = minimize(lambda v: abs(f(U_of(v)))**2 - L(U_of(v)), rng.uniform(-np.pi, np.pi, 16), method='BFGS')
        best = min(best, res.fun)
    print("min over U(4) of |%s|^2 - L = %+.6f" % (name, best), flush=True)
# reference: the true statement
best = np.inf
for r in range(30):
    res = minimize(lambda v: (abs(m(U_of(v),(1,2),(1,2))) + abs(m(U_of(v),(1,3),(1,4))))**2 - L(U_of(v)), rng.uniform(-np.pi, np.pi, 16), method='Nelder-Mead', options={'maxiter': 4000, 'fatol': 1e-13})
    best = min(best, res.fun)
print("reference: min of (|X|+|Y|)^2 - L = %+.3e" % best)
