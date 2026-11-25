"""Deepest violation of the exotic facet g over U(4) versus over O(4).
g(U) = |u11|^2 + |u33|^2 + |det U[{1,4},{3,4}]|^2 - |det U[{1,2},{1,4}]|^2   (1-based)."""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
def g(U):
    m = lambda I, J: abs(np.linalg.det(U[np.ix_(I, J)]))**2
    return abs(U[0, 0])**2 + abs(U[2, 2])**2 + m([0, 3], [2, 3]) - m([0, 1], [0, 3])
iu = np.triu_indices(4, 1)
def O_of(x):
    K = np.zeros((4, 4)); K[iu] = x; K = K - K.T; return expm(K)
def U_of(x):
    H = np.zeros((4, 4), complex); H[iu] = x[:6] + 1j * x[6:12]; H = H - H.conj().T; H[np.diag_indices(4)] = 1j * x[12:16]; return expm(H)
rng = np.random.default_rng(2026)
def best(f, dim, R):
    vals = []
    for r in range(R):
        x0 = rng.uniform(-np.pi, np.pi, dim)
        res = minimize(f, x0, method='BFGS', options={'gtol': 1e-12, 'maxiter': 5000})
        res = minimize(f, res.x, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 20000})
        vals.append((res.fun, res.x))
    vals.sort(key=lambda t: t[0]); return vals
vo = best(lambda x: g(O_of(x)), 6, 400)
vu = best(lambda x: g(U_of(x)), 16, 200)
print("min over O(4) (400 restarts): %.10f   (next distinct: %s)" % (vo[0][0], sorted(set(np.round([v[0] for v in vo[:50]], 8)))[:4]))
print("min over U(4) (200 restarts): %.10f   (next distinct: %s)" % (vu[0][0], sorted(set(np.round([v[0] for v in vu[:50]], 8)))[:4]))
print("Zheng's rational orthogonal example O*: -158/2401 = %.10f" % (-158 / 2401))
Ob = O_of(vo[0][1]); Ub = U_of(vu[0][1])
np.set_printoptions(precision=6, suppress=True)
print("best O:\n", Ob, "\n|best U|^2:\n", np.abs(Ub)**2)
np.save('moc/runs/orth_best.npy', Ob); np.save('moc/runs/unit_best.npy', Ub)
