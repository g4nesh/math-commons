"""Audit replay for final_referee_moc4_T11_uniform_20261022 (primary proof of paper case T_12).

Own code only; nothing from the supplement is read or executed.  Uses the shared
library ../lib/moc4lib.py.  Prints one line per check; exits 1 if an EXACT check fails.
Run:  python -I -B final_referee_moc4_T11_uniform_20261022__check.py [--numeric]
"""
import os, sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from moc4lib import (PERMS, P, F_VALUES, F_IS_MINOR, CONTACT_SETS, AUX7, generators, find_g,
                     orbit_rep, nullspace, values, leibniz_weights_np, so4_from_params)
import sympy as sp

FAIL = []
def check(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name)
    if not ok and level == "EXACT":
        FAIL.append(name)

T = frozenset({0, 2, 3, 5, 9, 12})
FV = lambda j: F_VALUES[j]

# --- A. contact set is paper T_12 ------------------------------------------------
check("contact set {0,2,3,5,9,12} equals paper T_12 (identity element)", T == CONTACT_SETS[12])

# --- source facets ---------------------------------------------------------------
gens = generators(T)
check("generators(T) = {0,1,2,5,7,9,12,39,40,41,43}", gens == [0, 1, 2, 5, 7, 9, 12, 39, 40, 41, 43])
check("non-minor generators are exactly f40, f43", [j for j in gens if not F_IS_MINOR[j]] == [40, 43])

# --- linear identities on W_T = {y in im P : y_T = 0} -----------------------------
WT = [values(c) for c in nullspace([P[k] for k in sorted(T)])]
check("dim W_T = 8", len(WT) == 8)
def on_WT(f):
    return all(f(y) == 0 for y in WT)
def expansion_ok(terms):
    return on_WT(lambda y: int(any(sum(cf(y) * FV(j)[s] for cf, j in terms) != y[s] for s in range(24))))
g = lambda k: (lambda y: y[k])
A = lambda y: y[16] + y[6] - y[8]
B = lambda y: y[13] + y[18] - y[7]
d = lambda y: -A(y)
e = lambda y: -B(y)
check("(1) y = y23 f0 + y13 f2 + y11 f5 + y16 f7 + y6 f12 + y17 f39 + y7 f41 + d f40 on W_T",
      expansion_ok([(g(23), 0), (g(13), 2), (g(11), 5), (g(16), 7), (g(6), 12), (g(17), 39), (g(7), 41), (d, 40)]))
check("(2) y = y23 f1 + y21 f5 + y14 f7 + y13 f9 + y19 f12 + y15 f39 + y18 f41 + e f43 on W_T",
      expansion_ok([(g(23), 1), (g(21), 5), (g(14), 7), (g(13), 9), (g(19), 12), (g(15), 39), (g(18), 41), (e, 43)]))
check("(3a) y1 = y17 + y7 + d on W_T", on_WT(lambda y: y[1] - y[17] - y[7] - d(y)))
check("(3b) y1 = y15 + y18 + e on W_T", on_WT(lambda y: y[1] - y[15] - y[18] - e(y)))
for L, R in [((40, 1), (0, 7, 39)), ((40, 9), (2, 7, 41)), ((43, 0), (1, 5, 12, 41)),
             ((43, 2), (5, 9, 12, 39)), ((40, 43), (5, 7, 12, 39, 41))]:
    check("replacement f%s = f%s" % ("+f".join(map(str, L)), "+f".join(map(str, R))),
          all(sum(FV(j)[s] for j in L) == sum(FV(j)[s] for j in R) for s in range(24)))
AB = {j: (A(FV(j)), B(FV(j))) for j in gens}
print("INFO (A,B) on generators:", AB)
check("A+B >= 0 on every generator (hence on N_T by Lemma 4.2)", all(a + b >= 0 for a, b in AB.values()))
check("A,B >= 0 => minor decomposition: A(f)<=0 off {1,9,43}, B(f)<=0 off {0,2,40}, A(f1)=A(f9)=B(f0)=B(f2)=1, A(f40)=B(f43)=-1",
      all(AB[j][0] <= 0 for j in gens if j not in (1, 9, 43)) and all(AB[j][1] <= 0 for j in gens if j not in (0, 2, 40))
      and AB[1][0] == AB[9][0] == AB[0][1] == AB[2][1] == 1 and AB[40][0] == AB[43][1] == -1)
check("g4 = g1 on W_T", on_WT(lambda y: y[4] - y[1]))

# --- seventh-contact witnesses, with the paper's Table 6 stages --------------------
STAGE = {}
for s_, l_ in {1: [1, 2, 3, 4, 5, 6, 7, 15, 18, 19, 20, 21, 22, 23, 25, 27, 28, 29, 31, 34, 35, 36, 37, 39, 40, 42, 43, 44, 72, 88],
               2: [10, 26, 46, 90], 3: [9, 11, 12, 47], 4: [8, 13, 17, 49], 5: [14, 16, 24, 33, 38, 70]}.items():
    for k in l_:
        STAGE[k] = s_
REP = {orbit_rep(S): k for k, S in CONTACT_SETS.items()}
def allminor(S):
    return all(F_IS_MINOR[j] for j in generators(S))
OWN_STAGE = 3
all_ok = True
for j in range(24):
    if j in T:
        continue
    U = set(T) | {j}
    early = []
    if allminor(U):
        early.append("all generators of T+j minor")
    for S6 in itertools.combinations(sorted(U), 6):
        k = REP.get(orbit_rep(S6))
        if k is not None and STAGE.get(k, 99) < OWN_STAGE:
            early.append("T_%d (stage %d) via subset %s" % (k, STAGE[k], S6))
    print("INFO j=%2d witnesses before stage %d: %s" % (j, OWN_STAGE, early[:2]))
    all_ok &= bool(early)
check("every j outside T has a witness: T+j all-minor or contains an image of a paper case of an earlier stage", all_ok)

# --- chart identities (symbolic) -----------------------------------------------
x, t, r, u, z, w, v = sp.symbols("x t r u z w v")
Asp = (z, 1, x, 0); Cl = [w, v, None, t * v]
Pp = [sp.expand(-sp.prod([Asp[i] - Cl[p[i]] for i in range(4) if p[i] != 2])) for p in PERMS]
Rr = [sp.cancel((Pp[k] - Pp[0]) / (v * (z - w))) for k in range(24)]
for k, val in {2: t * (x - 1), 3: x - t, 5: 1 - t, 1: x * (v - 1) / v}.items():
    check("R%d = %s" % (k, val), sp.simplify(Rr[k] - val) == 0)
check("(P2-P0)/(P3-P5) = t", sp.simplify((Pp[2] - Pp[0]) / (Pp[3] - Pp[5]) - t) == 0)
a_ = t * (v - 1) - u; D = t * (x - 1) - u; M = t + u
E1 = w * (v * (r + t - x) + z * (x - t * v)) + z * v * (-r + t * v - t)
E2 = a_ * z - D * w - t * (v - x)
check("(6) P9-P0-r v(z-w) = E1", sp.expand(Pp[9] - Pp[0] - r * v * (z - w) - E1) == 0)
check("(7) P12-P0-u v(z-w) = v E2", sp.expand(Pp[12] - Pp[0] - u * v * (z - w) - v * E2) == 0)
zsol = (D * w + t * (v - x)) / a_
quad = D * (x - t * v) * w**2 + (a_ * v * (r + t - x) + t * (v - x) * (x - t * v) + D * v * (-r + t * v - t)) * w + t * (v - x) * v * (-r + t * v - t)
check("(9) a(v) E1|z=(8) = quadratic", sp.expand(sp.cancel(a_ * E1.subs(z, zsol)) - quad) == 0)
check("z = w: E2 = t(v-x)(z-1) up to sign", sp.expand(E2.subs(z, w) + t * (x - v) * (w - 1)) == 0 or sp.expand(E2.subs(z, w) - t * (v - x) * (w - 1)) == 0)
check("z = w: E1 = z(v-z)(tv-x) up to sign", sp.expand(E1.subs(w, z) - z * (v - z) * (x - t * v)) == 0 or sp.expand(E1.subs(w, z) - z * (v - z) * (t * v - x)) == 0)
q0 = sp.factor(quad.subs(v, 0))
check("quadratic at v=0 is x w (D w - t x)", sp.expand(q0 - x * w * (D * w - t * x)) == 0)
dq = sp.diff(quad, w).subs(v, 0)
check("both v=0 roots simple (d/dw = -t x^2, t x^2)", sp.simplify(dq.subs(w, 0) + t * x**2) == 0 and sp.simplify(dq.subs(w, t * x / D) - t * x**2) == 0)
for wb, expect_set in [(t * x / D, {1, 4, 15, 17}), (0, {1, 4, 7, 10, 18, 20})]:
    zb = sp.cancel(zsol.subs({v: 0, w: wb}))
    res = [sp.cancel(((Pp[k] - Pp[0]) / (z - w)).subs({v: 0, z: zb, w: wb})) for k in range(24)]
    nz = {k for k in range(24) if res[k] != 0}
    check("branch (z,w)=(%s,%s): z-w != 0 and residues = -x on %s exactly" % (zb, wb, sorted(expect_set)),
          sp.simplify(zb - wb) != 0 and nz == expect_set and all(sp.simplify(res[k] + x) == 0 for k in nz))
    check("  residues of dR = R8-R16-R6 and eR = R7-R13-R18 vanish",
          sp.simplify(res[8] - res[16] - res[6]) == 0 and sp.simplify(res[7] - res[13] - res[18]) == 0)

# --- numeric falsification (optional) --------------------------------------------
if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(11)
    def minphi(y, restarts=6):
        y = np.asarray(y, float); best = np.inf
        for _ in range(restarts):
            rr = minimize(lambda th: float(y @ leibniz_weights_np(so4_from_params(th))), rng.uniform(-np.pi, np.pi, 6), method="BFGS")
            best = min(best, rr.fun)
        return best
    found = tried = 0; worst = np.inf
    while tried < 300000 and found < 40:
        tried += 1
        xx, tt, rr_, uu = rng.normal(size=4) * 2; mu = rng.choice([-1, 1]); vv = complex(rng.normal() * 2, rng.normal() * 2)
        DD = tt * (xx - 1) - uu; aa = tt * (vv - 1) - uu
        coeffs = [DD * (xx - tt * vv), aa * vv * (rr_ + tt - xx) + tt * (vv - xx) * (xx - tt * vv) + DD * vv * (-rr_ + tt * vv - tt), tt * (vv - xx) * vv * (-rr_ + tt * vv - tt)]
        for ww in np.roots(coeffs):
            zz = (DD * ww + tt * (vv - xx)) / aa
            A4 = (zz, 1, xx, 0); C4 = [ww, vv, None, tt * vv]
            Pv = np.array([-np.prod([A4[i] - C4[p[i]] for i in range(4) if p[i] != 2]) for p in PERMS])
            gg = -mu * ((Pv - Pv[0]) / (vv * (zz - ww))).imag
            gg = gg / abs(gg).max()
            if gg.min() < -1e-12:
                continue
            found += 1; worst = min(worst, minphi(gg))
    print("NUMERIC chart sampling (seed 11): %d draws, %d supported points, worst min_O Phi = %.3e" % (tried, found, worst))

print("SUMMARY: %d EXACT failures" % len(FAIL))
sys.exit(1 if FAIL else 0)
