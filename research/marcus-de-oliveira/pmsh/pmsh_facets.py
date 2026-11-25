"""Facets of the PMSH polytope P_c = conv{(prod_{t in T} c_pi(t))_T : pi in S_4}, as exact rational inequalities
sum_T y_T m_T <= h  (m_T = det H[T,T]).  Classify which are implied by trivial positivity (via the compound polytope),
and give one representative per S_4 orbit of the nontrivial ones."""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.spatial import ConvexHull
import sympy as sp
n = 4
Ts = [T for k in (1, 2, 3) for T in itertools.combinations(range(n), k)]
perms = list(itertools.permutations(range(n)))
def vert(c, p): return [np.prod([c[p[t]] for t in T]) for T in Ts]
def facets(c):
    V = np.array([vert(c, p) for p in perms], float); cen = V.mean(0); X = V - cen
    u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); B = vt[:r]
    h = ConvexHull(X @ B.T)
    eqr = np.round(h.equations / np.abs(h.equations[:, :-1]).max(1, keepdims=True), 8)
    _, idx = np.unique(eqr, axis=0, return_index=True)
    out = []
    Vex = [[Fr(x) for x in vert([Fr(ci) for ci in c], p)] for p in perms]
    for f in h.equations[idx]:
        on = np.where(np.abs((X @ B.T) @ f[:-1] + f[-1]) < 1e-9)[0]
        # exact facet normal: y with y.(v_i - v_j) = 0 for on-facet vertices, normalised; solve nullspace exactly
        A = sp.Matrix([[Vex[i][k] - Vex[on[0]][k] for k in range(len(Ts))] for i in on[1:]])
        # also enforce level-constant gauge: y orthogonal to level indicator vectors (remove the 3 trivial directions)
        lev = [[1 if len(T) == k else 0 for T in Ts] for k in (1, 2, 3)]
        ns = (A.col_join(sp.Matrix(lev))).nullspace()
        assert len(ns) == 1
        y = ns[0]; den = sp.ilcm(*[sp.fraction(v)[1] for v in y]); y = y * den
        g = sp.igcd(*[int(v) for v in y if v != 0]); y = y / g
        vals = [sum(y[k] * Vex[i][k] for k in range(len(Ts))) for i in range(24)]
        hval = max(vals)
        if vals[on[0]] != hval: y = -y; vals = [-v for v in vals]; hval = max(vals)
        out.append((tuple(int(v) for v in y), hval, tuple(sorted(on))))
    return out
def orbit_key(y):
    best = None
    for s in perms:
        Tidx = {T: k for k, T in enumerate(Ts)}
        z = tuple(y[Tidx[tuple(sorted(s[t] for t in T))]] for T in Ts)
        best = z if best is None or z < best else best
    return best
if __name__ == '__main__':
    c = [int(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3, 4, 6, 12]
    F = facets(c)
    orbs = {}
    for y, h, on in F: orbs.setdefault(orbit_key(y), []).append((y, h, on))
    print("c =", c, ": facets", len(F), "; S_4-orbits:", len(orbs))
    lab = ["m_" + "".join(str(t + 1) for t in T) for T in Ts]
    for key, mem in sorted(orbs.items(), key=lambda kv: len(kv[1][0][2])):
        y, h, on = mem[0]
        terms = " ".join("%+d*%s" % (v, l) for v, l in zip(y, lab) if v != 0)
        print("  orbit size %2d, %2d vertices on facet:  %s  <=  %s" % (len(mem), len(on), terms, h))
