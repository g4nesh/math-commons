"""End-to-end numerical sanity check of the open-set theorem (not part of the proof).
Random spectra c1 < c2 < c3 < c4 with cross-ratio lambda = (c2-c1)(c4-c3)/((c3-c1)(c4-c2)) in [a, 1-a]; for each, unitaries
including ones deep in the exotic caps and their permuted copies; solve the PMSH certificate LP (weights on S_4 matching all
principal minors of H = U D_c U^*) and report the largest infeasibility.  Usage: python3 end_to_end_interval.py <a> [n]"""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
PERMS = list(itertools.permutations(range(4))); SUBS = [T for k in range(1, 4) for T in itertools.combinations(range(4), k)]
def haar(rng):
    z = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / np.sqrt(2); Q, R = np.linalg.qr(z); d = np.diag(R); return Q * (d / np.abs(d))
def gap(U, c):
    H = (U * c) @ U.conj().T; m = np.array([np.linalg.det(H[np.ix_(T, T)]).real for T in SUBS])
    V = np.array([[np.prod([c[p[t]] for t in T]) for T in SUBS] for p in PERMS]); sc = np.abs(V).max(); d = V.shape[1]
    A = np.vstack([np.hstack([V.T / sc, -np.ones((d, 1))]), np.hstack([-V.T / sc, -np.ones((d, 1))])])
    r = linprog(np.r_[np.zeros(24), 1.0], A_ub=A, b_ub=np.r_[m / sc, -m / sc], A_eq=np.r_[np.ones(24), 0][None], b_eq=[1.0], bounds=[(0, None)] * 25, method='highs')
    return r.fun
def cap_unitary(y):
    S = np.array([[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]], float); return (y * np.eye(4) + np.sqrt((1 - y * y) / 3) * S).astype(complex)
a = float(Fr(sys.argv[1])); n = int(sys.argv[2]) if len(sys.argv) > 2 else 40; rng = np.random.default_rng(5); worst = 0; pairs = 0; lams = []
caps = [P @ cap_unitary(y) @ Q for y in (0.168231, 1 / 7, 0.25) for P, Q in [(np.eye(4)[list(p)], np.eye(4)[list(q)]) for p, q in zip(PERMS[::4], PERMS[::5])]]
while len(lams) < n:
    c = np.sort(rng.standard_normal(4)) * 3; x, y, z = np.diff(c); lam = x * z / ((x + y) * (y + z))
    if not (a <= lam <= 1 - a): continue
    lams.append(lam)
    for U in [haar(rng) for _ in range(20)] + caps: worst = max(worst, gap(U, c)); pairs += 1
print("%d random spectra with cross-ratio in [%.4f, %.4f] (observed range %.3f .. %.3f), %d (U, c) pairs incl. %d cap unitaries each: largest LP infeasibility %.2e" % (n, a, 1 - a, min(lams), max(lams), pairs, len(caps), worst))
