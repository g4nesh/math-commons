"""Orientation determinants of all 12-subsets of the 24 projected vertices, as integer quadratics in rho.
Projected vertex (affine coords): B[0:3,0:3] (9 entries), n_11 + rho n_12, n_21 + rho n_22   [p = (1, rho, 0)].
det(rho) = c0 + c1 rho + c2 rho^2, obtained from exact-integer-valued determinants at rho = 0, 1, -1 (computed in
floating point and rounded; every entry is a small integer).  Output: the set of real roots = the only values of rho
at which the oriented matroid, hence the facet structure, of the projected polytope can change."""
import itertools, time, numpy as np
from fractions import Fraction as Fr
from base import *
t0 = time.time()
Bc = np.array([[v[iB(s, j)] for s in range(3) for j in range(3)] for v in VX])
Nn = np.array([[[n_matrix(v.astype(float))[0][a, b] for b in range(3)] for a in range(3)] for v in VX]).astype(int)
def pts(rho): return np.hstack([np.ones((24, 1), dtype=int), Bc, (Nn[:, 0, 0] + rho * Nn[:, 0, 1])[:, None], (Nn[:, 1, 0] + rho * Nn[:, 1, 1])[:, None]])
P0, P1, Pm = pts(0).astype(float), pts(1).astype(float), pts(-1).astype(float)
subs = np.array(list(itertools.combinations(range(24), 12)), dtype=np.int8)
print("12-subsets:", len(subs), flush=True)
C = np.zeros((len(subs), 3), dtype=np.int64); maxdev = 0
for lo in range(0, len(subs), 150000):
    idx = subs[lo:lo + 150000].astype(int)
    d = [np.linalg.det(P[idx]) for P in (P0, P1, Pm)]
    r = [np.rint(x) for x in d]; maxdev = max(maxdev, max(np.abs(x - y).max() for x, y in zip(d, r)))
    d0, d1, dm = [x.astype(np.int64) for x in r]
    assert ((d1 - dm) % 2 == 0).all()
    C[lo:lo + 150000, 0] = d0; C[lo:lo + 150000, 1] = (d1 - dm) // 2; C[lo:lo + 150000, 2] = (d1 + dm) // 2 - d0
print("max |float det - nearest integer| = %.2e ; max |coefficient| = %d  [%.0fs]" % (maxdev, np.abs(C).max(), time.time() - t0), flush=True)
np.save('chiro_coeffs.npy', C)
nz = np.abs(C).sum(1) > 0
print("identically zero determinants: %d ; nonzero: %d" % ((~nz).sum(), nz.sum()))
# normalise each nonzero quadratic to a primitive tuple with positive leading nonzero coefficient; collect distinct ones
U_, cnt = np.unique(C[nz], axis=0, return_counts=True)
prim = {}
for row, k in zip(U_, cnt):
    g = np.gcd.reduce(np.abs(row)); r = tuple(int(v) // int(g) for v in row)
    lead = next(v for v in r[::-1] if v != 0)
    if lead < 0: r = tuple(-v for v in r)
    prim[r] = prim.get(r, 0) + int(k)
print("distinct primitive quadratics (c0, c1, c2) and how many subsets carry them:")
import sympy as sp
x = sp.symbols('rho'); roots = set()
for r, k in sorted(prim.items(), key=lambda kv: -kv[1]):
    poly = r[0] + r[1] * x + r[2] * x**2
    rr = sp.real_roots(sp.Poly(poly, x)) if sp.degree(poly, x) > 0 else []
    roots |= set(rr)
    print("   %-16s x %-8d factor: %-28s real roots: %s" % (r, k, sp.factor(poly), sorted(set(rr))))
print("ALL critical values of rho:", sorted(roots))
