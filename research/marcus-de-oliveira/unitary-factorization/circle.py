"""Numerical sanity check of Corollaries P4' and P4'' (floating point; not part of the proof).

Random normal 4x4 pairs A = D_a, B = U D_b U^* with b on a random line or circle and a arbitrary complex: is det(A + B) in
the convex hull of the 24 products prod(a_i + b_sigma(i))?  The test is a linear programme in the plane; the number printed
is the largest distance from the hull found (0 = inside).  A control with unrestricted complex b is included: the
conjecture predicts 0 there as well, but nothing in this directory proves it."""
import itertools, sys, numpy as np
from scipy.optimize import linprog
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); NT = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
PERMS = list(itertools.permutations(range(4)))
def haar():
    z = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / np.sqrt(2); q, r = np.linalg.qr(z); d = np.diag(r); return q * (d / np.abs(d))
def outside(a, b, U):
    z = np.linalg.det(np.diag(a) + (U * b) @ U.conj().T); pts = np.array([np.prod(a + b[list(p)]) for p in PERMS]); sc = np.abs(pts).max()
    A = np.array([np.r_[pts.real / sc, 0], np.r_[pts.imag / sc, 0]]); rhs = np.array([z.real, z.imag]) / sc
    A_ub = np.vstack([np.hstack([A[:, :24], -np.ones((2, 1))]), np.hstack([-A[:, :24], -np.ones((2, 1))])])
    r = linprog(np.r_[np.zeros(24), 1.0], A_ub=A_ub, b_ub=np.r_[rhs, -rhs], A_eq=np.r_[np.ones(24), 0][None], b_eq=[1.0], bounds=[(0, None)] * 25, method='highs'); return r.fun
def cplx(k=4): return rng.standard_normal(k) + 1j * rng.standard_normal(k)
w = {'line': 0, 'circle': 0, 'unitary B': 0, 'repeated eigenvalue': 0, 'control: general b': 0}
for _ in range(NT):
    a = cplx(); U = haar(); c0 = cplx(1)[0]; dr = np.exp(1j * rng.uniform(0, 2 * np.pi)); rad = rng.uniform(0.2, 3)
    w['line'] = max(w['line'], outside(a, c0 + dr * rng.standard_normal(4), U))
    w['circle'] = max(w['circle'], outside(a, c0 + rad * np.exp(1j * rng.uniform(0, 2 * np.pi, 4)), U))
    w['unitary B'] = max(w['unitary B'], outside(a, np.exp(1j * rng.uniform(0, 2 * np.pi, 4)), U))
    b = cplx(); b[3] = b[2]; w['repeated eigenvalue'] = max(w['repeated eigenvalue'], outside(a, b, U))
    w['control: general b'] = max(w['control: general b'], outside(a, cplx(), U))
print("%d random pairs per class; largest relative distance from the hull:" % NT)
for k, v in w.items(): print("   %-22s %.1e" % (k, v))
