"""Numerical sanity checks for README.md (floating point; not part of the proof).

 (a) minimum of E(A, B) = M22^2 + M33^2 + N23^2 + (|M11| + |M23|)^2 - N11^2 over SO(3) x SO(3), many BFGS starts;
 (b) (star) on random matrices of O(4) with both signs of the determinant, and E = 4 x margin;
 (c) the PMSH certificate LP for random real symmetric 4x4 matrices: spectra with cross-ratio below 1/10 (not covered by
     the interval certificates of ../star-halfline) and orthogonal matrices from Drury's extremal family;
 (d) the same LP for complex unitaries with small cross-ratio (covered only conditionally, by Theorem H)."""
import itertools, numpy as np
from scipy.optimize import minimize, linprog
rng = np.random.default_rng(11)
def rot(p):
    a, b, c, d = p
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)], [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
def E_of(A, B):
    M = A + B; N = A - B; return M[1, 1]**2 + M[2, 2]**2 + N[1, 2]**2 + (abs(M[0, 0]) + abs(M[1, 2]))**2 - N[0, 0]**2
def Eq(x): return E_of(rot(x[:4] / np.linalg.norm(x[:4])), rot(x[4:] / np.linalg.norm(x[4:])))
vals = np.array([minimize(Eq, rng.standard_normal(8), method='BFGS').fun for _ in range(2000)])
print("(a) min of E over SO(3) x SO(3), 2000 BFGS starts: %.2e ; starts ending below 1e-8: %d" % (vals.min(), (vals < 1e-8).sum()))
def haar_o(n):
    q, r = np.linalg.qr(rng.standard_normal((n, n))); return q * np.sign(np.diag(r))
def margin(O):
    D = lambda S, T: np.linalg.det(O[np.ix_(S, T)])**2
    return D((0, 2), (0, 2)) + D((0, 3), (0, 3)) + D((0, 2), (1, 2)) + (np.sqrt(D((0, 1), (0, 1))) + np.sqrt(D((0, 2), (0, 3))))**2 - D((0, 1), (2, 3))
mins = [np.inf, np.inf]
for _ in range(200000):
    O = haar_o(4); mins[int(np.linalg.det(O) < 0)] = min(mins[int(np.linalg.det(O) < 0)], margin(O))
print("(b) smallest margin of (star) on 200,000 Haar matrices of O(4): det +1: %.3e ; det -1: %.3e" % tuple(mins))
PERMS = list(itertools.permutations(range(4))); SUBS = [T for k in range(1, 4) for T in itertools.combinations(range(4), k)]
def gap(U, c):
    H = (U * c) @ U.conj().T; m = np.array([np.linalg.det(H[np.ix_(T, T)]).real for T in SUBS])
    V = np.array([[np.prod([c[p[t]] for t in T]) for T in SUBS] for p in PERMS]); sc = np.abs(V).max(); d = V.shape[1]
    A = np.vstack([np.hstack([V.T / sc, -np.ones((d, 1))]), np.hstack([-V.T / sc, -np.ones((d, 1))])])
    return linprog(np.r_[np.zeros(24), 1.0], A_ub=A, b_ub=np.r_[m / sc, -m / sc], A_eq=np.r_[np.ones(24), 0][None], b_eq=[1.0], bounds=[(0, None)] * 25, method='highs').fun
def spectrum(lo, hi):
    while True:
        c = np.sort(rng.standard_normal(4)); x, y, z = np.diff(c); lam = x * z / ((x + y) * (y + z))
        if lo <= lam <= hi: return 1 + (c - c[0]) / (c[3] - c[0])
S = np.array([[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]], float)
fam = [y * np.eye(4) + np.sqrt((1 - y * y) / 3) * S for y in (0.168231, 1 / 7, 0.1, 0.25)]
fam = [P @ F @ Q for F in fam for P, Q in [(np.eye(4)[list(p)], np.eye(4)[list(q)]) for p, q in zip(PERMS[::4], PERMS[::5])]]
worst = 0; n = 0
for _ in range(40):
    c = spectrum(0.0005, 0.1)
    for O in [haar_o(4) for _ in range(20)] + fam[::3]: worst = max(worst, gap(O, c)); n += 1
print("(c) PMSH certificate LP, real symmetric X = O D_c O^T, 40 spectra with cross-ratio in [0.0005, 0.1], %d pairs incl. Drury's family: largest infeasibility %.2e" % (n, worst))
def haar_u(n):
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
worst = 0; n = 0
for _ in range(40):
    c = spectrum(0.0005, 0.1)
    for U in [haar_u(4) for _ in range(20)]: worst = max(worst, gap(U, c)); n += 1
print("(d) the same for complex unitaries (covered by Theorem H only, i.e. conditionally), %d pairs: largest infeasibility %.2e" % (n, worst))
