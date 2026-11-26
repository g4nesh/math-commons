"""Audit replay for final_referee_moc4_T12_uniform_20261022 (primary proof of paper case T_13).

Own code only; nothing from the supplement is read or executed.  Uses ../lib/moc4lib.py.
Prints one line per check; exits 1 if an EXACT check fails.
Run:  python -I -B final_referee_moc4_T12_uniform_20261022__check.py [--numeric]
"""
import os, sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from moc4lib import (PERMS, F_VALUES, F_IS_MINOR, CONTACT_SETS, generators, orbit_rep,
                     leibniz_weights_np, so4_from_params)
import sympy as sp

FAIL = []
def check(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name)
    if not ok and level == "EXACT":
        FAIL.append(name)

T = frozenset({1, 2, 3, 4, 8, 12})
FV = lambda j: F_VALUES[j]
check("contact set {1,2,3,4,8,12} equals paper T_13 (identity element)", T == CONTACT_SETS[13])
gens = generators(T)
check("generators(T) = {0,1,2,5,6,9,10,53,54,55,56,57}", gens == [0, 1, 2, 5, 6, 9, 10, 53, 54, 55, 56, 57])
check("non-minor generators are exactly f54, f56, f57", [j for j in gens if not F_IS_MINOR[j]] == [54, 56, 57])
check("value at index 0: 1 for f53..f57, 0 for the other generators (so the sum of non-minor coefficients <= g0)",
      all(FV(j)[0] == (1 if j >= 53 else 0) for j in gens))
check("f54, f56, f57 are >= 1 at each of 0, 5, 21, 23", all(FV(j)[s] >= 1 for j in (54, 56, 57) for s in (0, 5, 21, 23)))

STAGE = {}
for s_, l_ in {1: [1, 2, 3, 4, 5, 6, 7, 15, 18, 19, 20, 21, 22, 23, 25, 27, 28, 29, 31, 34, 35, 36, 37, 39, 40, 42, 43, 44, 72, 88],
               2: [10, 26, 46, 90], 3: [9, 11, 12, 47], 4: [8, 13, 17, 49], 5: [14, 16, 24, 33, 38, 70]}.items():
    for k in l_:
        STAGE[k] = s_
REP = {orbit_rep(S): k for k, S in CONTACT_SETS.items()}
def allminor(S):
    return all(F_IS_MINOR[j] for j in generators(S))
OWN_STAGE = 4
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
check("T+j all-minor exactly for j in {0,5,21,23}", {j for j in range(24) if j not in T and allminor(set(T) | {j})} == {0, 5, 21, 23})

x, t, r, u, z, w, v = sp.symbols("x t r u z w v")
Asp = (z, 1, 0, x); Cl = [w, v, None, t * v]
Pp = [sp.expand(-sp.prod([Asp[i] - Cl[p[i]] for i in range(4) if p[i] != 2])) for p in PERMS]
Rr = [sp.cancel((Pp[k] - Pp[1]) / (v * (z - w))) for k in range(24)]
for k, val in {2: x - t, 3: t * (x - 1), 4: 1 - t, 0: x * (v - 1) / v}.items():
    check("R%d = %s" % (k, val), sp.simplify(Rr[k] - val) == 0)
check("(P3-P1)/(P2-P4) = t", sp.simplify((Pp[3] - Pp[1]) / (Pp[2] - Pp[4]) - t) == 0)
a_ = t * (v - 1) - u; D = x - t - u; M = t + u
E1 = w * (v * (r + t - x) + z * (x - t * v)) + z * v * (-r + t * v - t)
E2 = a_ * z - D * w - t * v + x
check("(3) P8-P1-r v(z-w) = E1", sp.expand(Pp[8] - Pp[1] - r * v * (z - w) - E1) == 0)
check("(4) P12-P1-u v(z-w) = v E2", sp.expand(Pp[12] - Pp[1] - u * v * (z - w) - v * E2) == 0)
zsol = (D * w + t * v - x) / a_
quad = D * (x - t * v) * w**2 + (a_ * v * (r + t - x) + (t * v - x) * (x - t * v) + D * v * (-r + t * v - t)) * w + (t * v - x) * v * (-r + t * v - t)
check("(6) a(v) E1|z=(5) = quadratic", sp.expand(sp.cancel(a_ * E1.subs(z, zsol)) - quad) == 0)
check("z = w: E2 = (tv-x)(z-1)", sp.expand(E2.subs(z, w) - (t * v - x) * (w - 1)) == 0)
check("z = w: E1 = z(v-z)(tv-x)", sp.expand(E1.subs(w, z) - z * (v - z) * (t * v - x)) == 0)
check("quadratic at v=0 is x w (D w - x)", sp.expand(quad.subs(v, 0) - x * w * (D * w - x)) == 0)
dq = sp.diff(quad, w).subs(v, 0)
check("d/dw at w=0, x/D are -x^2, x^2", sp.simplify(dq.subs(w, 0) + x**2) == 0 and sp.simplify(dq.subs(w, x / D) - x**2) == 0)
for wb, m in [(x / D, 53), (0, 55)]:
    zb = sp.cancel(zsol.subs({v: 0, w: wb}))
    res = [sp.cancel(((Pp[k] - Pp[1]) / (z - w)).subs({v: 0, z: zb, w: wb})) for k in range(24)]
    check("branch (z,w)=(%s,%s): z-w != 0 and residue vector = -x f%d" % (zb, wb, m),
          sp.simplify(zb - wb) != 0 and all(sp.simplify(res[k] + x * FV(m)[k]) == 0 for k in range(24)))
check("f53, f55 are minor vectors", F_IS_MINOR[53] and F_IS_MINOR[55])

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
    found = tried = 0; worst = np.inf; worstbar = np.inf
    while tried < 300000 and found < 40:
        tried += 1
        xx, tt, rr_, uu = rng.normal(size=4) * 2; mu = rng.choice([-1, 1]); vv = complex(rng.normal() * 2, rng.normal() * 2)
        DD = xx - tt - uu; aa = tt * (vv - 1) - uu
        coeffs = [DD * (xx - tt * vv), aa * vv * (rr_ + tt - xx) + (tt * vv - xx) * (xx - tt * vv) + DD * vv * (-rr_ + tt * vv - tt), (tt * vv - xx) * vv * (-rr_ + tt * vv - tt)]
        for ww in np.roots(coeffs):
            zz = (DD * ww + tt * vv - xx) / aa
            A4 = (zz, 1, 0, xx); C4 = [ww, vv, None, tt * vv]
            Pv = np.array([-np.prod([A4[i] - C4[p[i]] for i in range(4) if p[i] != 2]) for p in PERMS])
            gg = -mu * ((Pv - Pv[1]) / (vv * (zz - ww))).imag
            gg = gg / abs(gg).max()
            if gg.min() < -1e-12:
                continue
            found += 1; worst = min(worst, minphi(gg))
            for _ in range(50):
                W = leibniz_weights_np(so4_from_params(rng.uniform(-np.pi, np.pi, 6)))
                CO = max(0, -(np.array(FV(54)) @ W), -(np.array(FV(56)) @ W), -(np.array(FV(57)) @ W))
                worstbar = min(worstbar, gg @ W + CO * gg[0])
    print("NUMERIC chart sampling (seed 11): %d draws, %d supported points, worst min_O Phi = %.3e, worst L_O + C_O g0 over 50 random O each = %.3e" % (tried, found, worst, worstbar))

print("SUMMARY: %d EXACT failures" % len(FAIL))
sys.exit(1 if FAIL else 0)
