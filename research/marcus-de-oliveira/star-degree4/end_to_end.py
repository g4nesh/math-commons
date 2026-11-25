"""End-to-end numerical sanity check of the spectrum theorems (not part of the proof).
For the spectra c of the theorems (default: harmonic; with argument `more`: cross-ratio classes 1/3, 1/4, 2/5, 1/6) and many unitaries U, including ones deep in the exotic caps, solve the PMSH certificate LP directly
(weights on S_4 matching all principal minors of H = U D_c U^*) and report the largest infeasibility.
Also checks the Marcus-de Oliveira inclusion det(D_a + H) in conv{prod(a_i + c_sigma(i))} for random complex a."""
import sys, itertools, numpy as np
from scipy.optimize import linprog
from scipy.linalg import expm
sys.path.insert(0, '../pmsh/star')
from base import haar
PERMS = list(itertools.permutations(range(4))); SUBS = [T for k in range(1, 4) for T in itertools.combinations(range(4), k)]
def pm(H): return np.array([np.linalg.det(H[np.ix_(T, T)]).real for T in SUBS])
def gap(U, c):
    H = (U * c) @ U.conj().T; m = pm(H); V = np.array([[np.prod([c[p[t]] for t in T]) for T in SUBS] for p in PERMS]); sc = np.abs(V).max()
    d = V.shape[1]; A = np.vstack([np.hstack([V.T / sc, -np.ones((d, 1))]), np.hstack([-V.T / sc, -np.ones((d, 1))])])
    r = linprog(np.r_[np.zeros(24), 1.0], A_ub=A, b_ub=np.r_[m / sc, -m / sc], A_eq=np.r_[np.ones(24), 0][None], b_eq=[1.0], bounds=[(0, None)] * 25, method='highs')
    return r.fun, r.x[:24]
def cap_unitary(y):
    S = np.array([[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]], float); return (y * np.eye(4) + np.sqrt((1 - y * y) / 3) * S).astype(complex)
rng = np.random.default_rng(0); worst = 0; worst_moc = -np.inf; n = 0
MORE = len(sys.argv) > 1 and sys.argv[1] == 'more'
SPECTRA = [[0, 2, 3, 4.], [0, 1, 2, 4.], [0, 1, 2, 3.], [-3, -1, 1, 3.], [0, 1, 2, 6.], [0, 4, 5, 6.], [0, 1, 3, 5.], [0, 2, 4, 5.]] if MORE else [[0, 2, 3, 6.], [1, 5, 7, 13.], [-3, -1, 0, 3.]]
for c in map(np.array, SPECTRA):
    Us = [haar(rng) for _ in range(200)] + [cap_unitary(y) for y in (0.168231, 1 / 7, 0.1, 0.25)]
    Us += [P @ cap_unitary(0.168231) @ Qm for P, Qm in [(np.eye(4)[list(p)], np.eye(4)[list(q)]) for p, q in zip(PERMS[::5], PERMS[::7])]]
    for U in Us:
        g, lam = gap(U, c); worst = max(worst, g); n += 1
        a = rng.standard_normal(4) + 1j * rng.standard_normal(4); H = (U * c) @ U.conj().T
        dval = np.linalg.det(np.diag(a) + H); z = np.array([np.prod([a[i] + c[p[i]] for i in range(4)]) for p in PERMS])
        worst_moc = max(worst_moc, abs(dval - lam @ z) / (1 + abs(dval)))
print(("spectra of classes 1/3, 1/4, 2/5, 1/6" if MORE else "harmonic spectra") + ", %d (U, c) pairs incl. deep cap unitaries: largest LP infeasibility %.2e ; largest |det(D_a+H) - sum lam z| (relative) using the SAME weights for random complex a: %.2e" % (n, worst, worst_moc))
c = np.array([0, 1, 2, 3.]); g, _ = gap(cap_unitary(0.168231), c)
if not MORE: print("for comparison, non-harmonic spectrum (0,1,2,3) at the deepest cap unitary: LP infeasibility %.2e (PMSH there is not covered by this theorem)" % g)
