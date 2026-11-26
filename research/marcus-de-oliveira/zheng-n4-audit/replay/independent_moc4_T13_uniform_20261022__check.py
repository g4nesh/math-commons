"""Exact replay and numeric tests for the supplement note
"The contact family T13: two vanishing barriers" (independent_moc4_T13_uniform_20261022),
the primary proof of paper case T_14 = {0, 2, 3, 5, 11, 12}.

Own code only.  Nothing from the supplement is read or executed.  Uses the shared audit
library ../lib/moc4lib.py (generators, F_VALUES, group action, P).
Every printed line is one check.  Exit status is nonzero iff an EXACT check fails.
"""
import sys
import itertools
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from moc4lib import (P, PERMS, F_VALUES, F_IS_MINOR, generators, find_g, orbit_rep,
                     CONTACT_SETS, leibniz_weights_np, so4_from_params)
import sympy as sp

FAIL = []


def check(name, ok, level="EXACT", info=""):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}" + (f"   {info}" if info else ""), flush=True)
    if not ok and level == "EXACT":
        FAIL.append(name)


def is_zero(e):
    return sp.cancel(sp.together(e)) == 0


T = [0, 2, 3, 5, 11, 12]
FV = F_VALUES

# ---------------------------------------------------------------- A. contact set, generators
check("T is paper T_14 (identity group element)", find_g(CONTACT_SETS[14], set(T)) is not None)
gens = generators(T)
check("generators(T) = {0,1,2,7,8,9,10,12,39,40,41,44}", gens == [0, 1, 2, 7, 8, 9, 10, 12, 39, 40, 41, 44], info=str(gens))
check("non-minor generators are exactly f40, f44", [j for j in gens if not F_IS_MINOR[j]] == [40, 44])

# ---------------------------------------------------------------- linear face
cs = sp.symbols("c0:14")
y = sp.Matrix(P) * sp.Matrix(cs)
sol = sp.solve([y[j] for j in T], cs, dict=True)[0]
yW = [sp.expand(y[j].subs(sol)) for j in range(24)]
free = sorted(set().union(*[e.free_symbols for e in yW]), key=str)
check("face {Pc : (Pc)_T = 0} has dimension 8", len(free) == 8)
piv = [1, 6, 7, 8, 9, 13, 16, 18]
q = sp.symbols("q1 q6 q7 q8 q9 q13 q16 q18")
s2 = sp.solve([yW[p] - qq for p, qq in zip(piv, q)], free, dict=True)
check("q = (g1,g6,g7,g8,g9,g13,g16,g18) are coordinates on the face", len(s2) == 1)
g = [sp.expand(e.subs(s2[0])) for e in yW]
B, C, K = g[6], g[16], g[23]
al = g[7] - g[16] - g[18]
be = g[8] - g[6] - g[9] - g[16]
check("(1) alpha = g1 - g17 - g20", is_zero(al - (g[1] - g[17] - g[20])))
check("(1) beta = g1 - g7 - g17", is_zero(be - (g[1] - g[7] - g[17])))
check("(1) alpha + beta = -C - K", is_zero(al + be + C + K))
check("g4 = g1 on the face", is_zero(g[4] - g[1]))


def decomposition_ok(coefs):
    return all(is_zero(sum(cf * FV[j][s] for j, cf in coefs) - g[s]) for s in range(24))


check("(2) g = g21 f2 + g14 f8 + C f9 + K f10 + g22 f12 + g17 f39 + g20 f41 + alpha f44 (24 coords)",
      decomposition_ok([(2, g[21]), (8, g[14]), (9, C), (10, K), (12, g[22]), (39, g[17]), (41, g[20]), (44, al)]))
check("(3) g = K f0 + g13 f2 + C f7 + g9 f8 + B f12 + g17 f39 + g7 f41 + beta f40 (24 coords)",
      decomposition_ok([(0, K), (2, g[13]), (7, C), (8, g[9]), (12, B), (39, g[17]), (41, g[7]), (40, be)]))

# minor cone dual (16 extreme rays)
minors = [0, 1, 2, 7, 8, 9, 10, 12, 39, 41]
Gq = [sp.Matrix([FV[j][p] for p in piv]) for j in minors]
check("minor cone has dimension 8", sp.Matrix.hstack(*Gq).rank() == 8)
facets = set()
for sub in itertools.combinations(range(10), 7):
    M = sp.Matrix.hstack(*[Gq[i] for i in sub]).T
    if M.rank() != 7:
        continue
    nv = M.nullspace()[0]
    vals = [(nv.T * Gq[i])[0] for i in range(10)]
    if all(v >= 0 for v in vals) or all(v <= 0 for v in vals):
        if all(v <= 0 for v in vals):
            nv = -nv
        facets.add(tuple(nv / max(abs(t) for t in nv)))
claimed = [g[j] for j in (1, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 20, 21, 23)] + [-al, -be]


def nrm(e):
    v = sp.Matrix([sp.expand(e).coeff(qq) for qq in q])
    return tuple(v / max(abs(t) for t in v))


check("minor-cone dual = the 14 coordinate rays and -alpha, -beta (16 rays)",
      len(facets) == 16 and set(nrm(e) for e in claimed) == facets, info=f"{len(facets)} facets")

# ---------------------------------------------------------------- seventh-contact inventory
reps = {orbit_rep(CONTACT_SETS[k]): k for k in CONTACT_SETS}
claimed_family = {1: 9, 4: None, 6: 2, 7: 4, 8: 3, 9: 1, 10: 5, 13: 1, 14: 6, 15: 10, 16: 7,
                  17: 11, 18: 5, 19: 3, 20: 9, 21: 6, 22: 8, 23: 7}
# supplement historical labels 0..10 -> paper T_1, T_2, T_6, T_3, T_4, T_5, T_7, T_8, T_9, T_10, T_11
STAGE_T14 = 5
STAGE = {1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 4, 9: 3, 10: 2, 11: 3}
for j, k in claimed_family.items():
    S = set(T) | {j}
    if k is None:
        check(f"new contact j=4: on the face g4=g1, so it reduces to j=1", is_zero(g[4] - g[1]))
        continue
    hit = any(find_g(CONTACT_SETS[k], set(sub)) is not None for sub in itertools.combinations(sorted(S), 6))
    check(f"new contact j={j}: T u {{j}} contains an image of paper T_{k} (stage {STAGE[k]} < {STAGE_T14})",
          hit and STAGE[k] < STAGE_T14)

# ---------------------------------------------------------------- spectral chart
z, x, c, v, w, s, t, r, u, X, Z = sp.symbols("z x c v w s t r u X Z")
Asp, Csp = [z, 0, 1, x], [c, v, w, None]


def Fp(j):
    p = PERMS[j]
    return sp.prod([Asp[i] - Csp[p[i]] for i in range(4) if p[i] != 3])


N = (z - c) * (x - v)
check("F5 - F3 = (z-c)(x-v) (affine chart, factor at infinity omitted)", is_zero(Fp(5) - Fp(3) - N))
H = {j: (Fp(j) - Fp(0)) / N for j in range(24)}
sdef, tdef = (v - w) / (x - v), w * (1 - x) / (x - v)
rdef, udef = z * (w - c) * (1 - v) / N, (z - x) * (v - c) * (1 - w) / N
check("H0..H5 contacts: H2 = s, H3 = s+t, H5 = 1+s+t", is_zero(H[2] - sdef) and is_zero(H[3] - sdef - tdef) and is_zero(H[5] - 1 - sdef - tdef))
check("H12 = s + r, H11 = A + u", is_zero(H[12] - sdef - rdef) and is_zero(H[11] - 1 - sdef - tdef - udef))
D, E, A, Q, k = 1 + t + s * X, t + (s + 1) * X, 1 + s + t, (s + 1) * X + t * Z, (s + 1) / s
sx = {x: 1 - X}
vv, ww, zz = x * (s * X + t) / E, t * x / E, x / (1 - Z)
check("v = x(sX+t)/E, w = tx/E invert the definitions of s, t",
      is_zero(sdef.subs({v: vv, w: ww}).subs(sx) - s) and is_zero(tdef.subs({v: vv, w: ww}).subs(sx) - t))
csol = sp.solve(sp.Eq(rdef.subs({v: vv, w: ww, z: zz}), r), c)
cspin = x * (t * D - r * E) / (E * (D - r * (1 - Z)))
check("r determines c uniquely and c equals the first C spinor ratio", len(csol) == 1 and is_zero((csol[0] - cspin).subs(sx)))
S = {v: vv, w: ww, z: zz, c: cspin}


def onchart(e):
    return sp.cancel(sp.together(e.subs(S).subs(sx)))


ufun = onchart(udef)
G6 = u * D * Q - A * Z * (X * (s * D + r) + r * (s * X + t) * Z)
check("(6) the contact equation u = u(X,Z) is the cubic G6 = 0", is_zero(G6.subs(u, ufun)))
Hp = {j: onchart(H[j]) for j in range(24)}
Zu = (s + 1) * u / (s * u + A * (s + r))
check("slice X = -t/s: u -> Z(u) is a rational inverse", is_zero(ufun.subs(X, -t / s).subs(Z, Zu) - u))


def real_const(name, expr):
    d = sp.together(expr.subs(u, ufun))
    phi = sp.factor(sp.cancel(d.subs(X, -t / s).subs(Z, Zu)))
    ok = is_zero(d - phi.subs(u, ufun))
    check(name + " is a constant on every level set of u, rational in (s,t,r,u) with rational coefficients", ok, info=f"const = {phi}")


U, W = u / Z, r * (Z - X) / D
table = {1: s * X, 6: U, 7: (s + 1) / A * (1 - X) * U + s * X, 8: U - (s + t / X) * W, 9: -t * u / (A * X),
         13: r * (1 - X) / D, 16: W, 18: s * X - k * D * U / A}
for j, tj in table.items():
    real_const(f"(7) H{j} - table entry", Hp[j] - tj)
real_const("(8) K - (C - k(1+t)B/A)", Hp[23] - (Hp[16] - k * (1 + t) * Hp[6] / A))
real_const("(8) alpha - (kB - C)", (Hp[7] - Hp[16] - Hp[18]) - (k * Hp[6] - Hp[16]))
real_const("(8) beta - (-kB - K)", (Hp[8] - Hp[6] - Hp[9] - Hp[16]) - (-k * Hp[6] - Hp[23]))
real_const("H23 - (sX + rZ - kDU/A)", Hp[23] - (s * X + r * Z - k * D * U / A))
# spinors and normalization (10)
spA = [(x, 1 - Z), (0, 1), (1, 1), (x, 1)]
spC = [(x * (t * D - r * E), E * (D - r * (1 - Z))), (x * (s * X + t), E), (t * x, E), (1, 0)]
br = lambda a, b: a[0] * b[1] - a[1] * b[0]
pp = lambda j: sp.prod([br(spA[i], spC[PERMS[j][i]]) for i in range(4)])
Delta = sp.expand((pp(5) - pp(3)).subs(sx))
check("(10) Delta = F5 - F3 = -x^2 X D E ((s+1)X + tZ)", is_zero((Delta + x**2 * X * D * E * Q).subs(sx)))
check("spinor products reproduce (F_j - F_0)/Delta = H_j for all 24 j on (6)",
      all(is_zero((((pp(j) - pp(0)) / Delta).subs(sx) - Hp[j]).subs(u, ufun)) for j in range(24)))
gq = sp.factor(G6.subs(Z, -(s + 1) * X / t))
check("Q = 0 on (6) forces rE - tD = 0 (given X, s, s+1, A, t nonzero)",
      is_zero(gq + X**2 * s * (s + 1) * A * (r * E - t * D) / t**2), info=str(gq))
check("D/E has determinant st - (s+1)(1+t) = -A", is_zero(s * t - (s + 1) * (1 + t) + A))
check("Z = 0 on (6) forces u(s+1) D X = 0", is_zero(G6.subs(Z, 0) - X * u * (s + 1) * D))
check("restriction to X = 0 is tZ(u(1+t) - ArZ), not identically 0", is_zero(G6.subs(X, 0) - t * Z * (u * (1 + t) - A * r * Z)))
check("restriction to D = 0 is -ArZ(X - Z) with X = -(1+t)/s", is_zero(G6.subs(X, -(1 + t) / s) + A * r * Z * (-(1 + t) / s - Z)))
Tt = sp.symbols("T")
Gh = sp.expand(Tt**3 * G6.subs({X: X / Tt, Z: Z / Tt}))
check("leading form = -A s X Z (sX + rZ)", is_zero(Gh.subs(Tt, 0) + A * s * X * Z * (s * X + r * Z)))
g1c = sp.expand(Gh.subs(X, 1))
dZ, dT = sp.diff(g1c, Z).subs({Z: 0, Tt: 0}), sp.diff(g1c, Tt).subs({Z: 0, Tt: 0})
check("[1:0:0] smooth, Z/X = u(s+1)/(As) T + O(T^2), so U -> As/(s+1)",
      is_zero(dZ + s**2 * A) and is_zero(-dT / dZ - u * (s + 1) / (A * s)))
g2c = sp.expand(Gh.subs(Z, 1))
grad0 = (sp.diff(g2c, X).subs({X: 0, Tt: 0}), sp.diff(g2c, Tt).subs({X: 0, Tt: 0}))
gradr = sp.diff(g2c, X).subs({X: -r / s, Tt: 0})
check("[0:1:0] and [r:-s:0] smooth (gradient nonzero when s,r,A != 0)", is_zero(grad0[0] + r * s * A) and is_zero(gradr - r * s * A))

# ---------------------------------------------------------------- NUMERIC falsification
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(2026)


def minphi(gv, restarts=12):
    f = lambda th: float(gv @ leibniz_weights_np(so4_from_params(th)))
    return min(minimize(f, rng.uniform(-np.pi, np.pi, 6), method="BFGS").fun for _ in range(restarts))


def prods(a, b):
    return np.array([np.prod([a[i][0] * b[p[i]][1] - a[i][1] * b[p[i]][0] for i in range(4)]) for p in PERMS])


acc = tried = 0
worst = np.inf
worst9 = np.inf
while acc < 40 and tried < 200000:
    tried += 1
    sv, tv, rv, uv = rng.normal(size=4) * rng.choice([0.3, 1, 3], size=4)
    Xv = complex(*rng.normal(size=2))
    Dv, Ev, Av, xv = 1 + tv + sv * Xv, tv + (sv + 1) * Xv, 1 + sv + tv, 1 - Xv
    for Zv in np.roots([Av * rv * (sv * Xv + tv), Av * Xv * (sv * Dv + rv) - uv * Dv * tv, -uv * Dv * (sv + 1) * Xv]):
        a = [(xv, 1 - Zv), (0, 1), (1, 1), (xv, 1)]
        b = [(xv * (tv * Dv - rv * Ev), Ev * (Dv - rv * (1 - Zv))), (xv * (sv * Xv + tv), Ev), (tv * xv, Ev), (1, 0)]
        F = prods(a, b)
        Hv = (F - F[0]) / (F[5] - F[3])
        if max(abs(Hv[T].imag)) > 1e-8 * max(1, abs(Hv).max()):
            continue
        for mu in (1, -1):
            gv = -mu * Hv.imag
            if gv.min() < -1e-12:
                continue
            acc += 1
            worst = min(worst, minphi(gv / np.abs(gv).max()))
check(f"main statement on {acc} sampled realizable T14 sources (seed 2026, 12 BFGS restarts each): min Phi >= -1e-9",
      worst > -1e-9, level="NUMERIC", info=f"worst = {worst:.3e}, tried {tried}")

print("EXACT failures:", len(FAIL))
sys.exit(1 if FAIL else 0)
