"""Hard facets of the projected polytope in closed form:  g_F + alpha*D_a + alpha'*D_b >= 0.
For each exotic facet F and each pair of squared 2x2 minors (D_a, D_b), solve for alpha, alpha' >= 0 such that the
functional is constant along every shift eta q^T (so it is a function of the projected coordinates).  Then compare
the tight-vertex sets with the facets of the projected polytope computed independently by qhull."""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.spatial import ConvexHull
from base import *
from fiber import p_of, q_of, proj
EXd = exotic_facets(); EX = list(EXd.values())
def hull_facets(p):
    V = np.array([proj(v.astype(float), p) for v in VX]); cen = V.mean(0); X = V - cen
    u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); Bs = vt[:r]
    h = ConvexHull(X @ Bs.T); out = set()
    for f in h.equations:
        on = frozenset(np.where(np.abs((X @ Bs.T) @ f[:-1] + f[-1]) < 1e-8)[0]); out.add(on)
    return r, out
def candidates(q):
    sh = [shift(np.eye(3)[k] - np.eye(3)[2], q) for k in range(2)]
    Dk = list(range(16, NX)); out = []
    for fi, w in enumerate(EX):
        rhs = -np.array([w @ sh[0], w @ sh[1]])
        for a, b in itertools.combinations(Dk, 2):
            A = np.array([[sh[0][a], sh[0][b]], [sh[1][a], sh[1][b]]])
            if abs(np.linalg.det(A)) < 1e-12: continue
            al = np.linalg.solve(A, rhs)
            if al.min() < 1e-12: continue                    # need both strictly positive (else not 'hard')
            vals = VX @ w + al[0] * VX[:, a] + al[1] * VX[:, b]
            assert vals.min() > -1e-12
            out.append((fi, a, b, al, frozenset(np.where(np.abs(vals) < 1e-9)[0])))
    return out
def name(k):
    k -= 16; S, T = P2[k // 6], P2[k % 6]
    return "D[%s,%s]" % ("".join(str(s + 1) for s in S), "".join(str(t + 1) for t in T))
if __name__ == '__main__':
    c = np.array([float(Fr(v)) for v in (sys.argv[1].split(',') if len(sys.argv) > 1 else ['3', '4', '6', '12'])])
    p = p_of(c); q = q_of(p)
    print("c =", c, " p =", p, " q =", np.round(q, 6))
    r, HF = hull_facets(p)
    hard = {f for f in HF if len(f) == 11}
    print("projected polytope: dim %d, %d facets, %d with 11 vertices" % (r, len(HF), len(hard)))
    cand = candidates(q)
    cand_fac = [cnd for cnd in cand if len(cnd[4]) >= 11]
    sets = {cnd[4] for cnd in cand_fac}
    print("closed-form candidates: %d total, %d with >= 11 tight vertices, %d distinct tight sets" % (len(cand), len(cand_fac), len(sets)))
    print("hard facets reproduced by a closed form: %d of %d ; closed forms that are not facets: %d" % (len(hard & sets), len(hard), len(sets - HF)))
    per = {}
    for cnd in cand_fac:
        if cnd[4] in hard: per.setdefault(cnd[0], []).append(cnd)
    print("exotic facets contributing hard facets: %d ; hard facets per exotic facet: %s" % (len(per), sorted({len({x[4] for x in v}) for v in per.values()})))
    # print the forms for the base cap g0
    g0key = None
    for k, (key, w) in enumerate(EXd.items()):
        if w[iB(0, 0)] == 1 and w[iB(2, 2)] == 1 and w[iD((0, 3), (2, 3))] == 1 and w[iD((0, 1), (0, 3))] == -1 and np.abs(w).sum() == 4: g0key = k
    print("base cap g0 = B11 + B33 + D[14,34] - D[12,14]; its hard facets at this spectrum:")
    seen = set()
    for cnd in per.get(g0key, []):
        if cnd[4] in seen: continue
        seen.add(cnd[4]); print("   g0 + %.6f*%s + %.6f*%s >= 0" % (cnd[3][0], name(cnd[1]), cnd[3][1], name(cnd[2])))
