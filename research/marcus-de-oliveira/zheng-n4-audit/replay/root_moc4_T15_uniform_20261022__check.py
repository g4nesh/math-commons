"""Exact replay and numeric tests for the supplement note
"A six-contact support inequality from a quadratic spectral curve"
(root_moc4_T15_uniform_20261022), the primary proof of paper case
T_16 = {0, 2, 5, 6, 9, 12}; the note works with its transpose T = {0, 2, 5, 6, 8, 18}.

Own code only.  Nothing from the supplement is read or executed.  Uses ../lib/moc4lib.py.
Every printed line is one check.  Exit status is nonzero iff an EXACT check fails.
"""
import sys
import itertools
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from moc4lib import (P, PERMS, F_VALUES, F_IS_MINOR, generators, find_g, CONTACT_SETS, perm_index, inverse,
                     leibniz_weights_np, so4_from_params)
import sympy as sp

FAIL = []


def check(name, ok, level="EXACT", info=""):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}" + (f"   {info}" if info else ""), flush=True)
    if not ok and level == "EXACT":
        FAIL.append(name)


def is_zero(e):
    return sp.cancel(sp.together(e)) == 0


T = [0, 2, 5, 6, 8, 18]
FV = F_VALUES
check("T = {sigma^-1 : sigma in paper T_16} (transpose), and some g in G maps T_16 onto T",
      sorted(perm_index(inverse(PERMS[j])) for j in CONTACT_SETS[16]) == T and find_g(CONTACT_SETS[16], set(T)) is not None)
gens = generators(T)
check("generators(T) = {1,2,4,5,6,9,39,45,46,47}", gens == [1, 2, 4, 5, 6, 9, 39, 45, 46, 47], info=str(gens))
check("only f46 is non-minor", [j for j in gens if not F_IS_MINOR[j]] == [46])
add = lambda js: [sum(FV[j][s] for j in js) for s in range(24)]
check("f46 + f1 = f4 + f6 + f45", add([46, 1]) == add([4, 6, 45]))
check("f46 + f9 = f4 + f5 + f47", add([46, 9]) == add([4, 5, 47]))
check("d = g12 - g16 is 1 on f46, -1 on f1 and f9, 0 on the other generators",
      all(FV[j][12] - FV[j][16] == {46: 1, 1: -1, 9: -1}.get(j, 0) for j in gens))
cs = sp.symbols("c0:14")
y = sp.Matrix(P) * sp.Matrix(cs)
sol = sp.solve([y[j] for j in T], cs, dict=True)[0]
g = [sp.expand(y[j].subs(sol)) for j in range(24)]
d, a, b, c = g[12] - g[16], g[13] - g[7] - g[12], g[1] - g[3], g[7]
coefs = [(2, a), (4, g[16]), (5, g[10]), (6, g[9] - g[7] - d), (39, b), (45, g[3] - g[7] - d), (46, d), (47, g[7])]
check("(1) g = a f2 + g16 f4 + g10 f5 + (g9-g7-d) f6 + b f39 + (g3-g7-d) f45 + d f46 + g7 f47 (24 coords)",
      all(is_zero(sum(cf * FV[j][s] for j, cf in coefs) - g[s]) for s in range(24)))
lin = {2: lambda V: V[13] - V[7] - V[12], 4: lambda V: V[16], 5: lambda V: V[10],
       6: lambda V: V[9] - V[7] - (V[12] - V[16]), 39: lambda V: V[1] - V[3],
       45: lambda V: V[3] - V[7] - (V[12] - V[16]), 47: lambda V: V[7]}
check("the seven coefficients other than d are >= 0 on every generator (hence on N_T)",
      all(fn(FV[j]) >= 0 for fn in lin.values() for j in gens))
check("(2) hence 0 < d <= g3 - g7 <= g3 when d > 0 (coefficient of f45 and g7 >= 0)", True and all(lin[45](FV[j]) >= 0 for j in gens))
zs = sorted(s for s in range(24) if FV[46][s] == 0 and s not in T)
check("f46 vanishes off T exactly at {4,7,10,16,19,22,23}", zs == [4, 7, 10, 16, 19, 22, 23], info=str(zs))
check("f46 > 0 at every other index j not in T", all(FV[46][s] > 0 for s in range(24) if s not in T and s not in zs))
# witnesses: supplement labels T1,T7,T3,T14,T1,T3,T7 = paper T_2, T_8, T_3, T_15, T_2, T_3, T_8; T_16 is stage 5
STAGE = {2: 1, 3: 1, 8: 4, 15: 1}
for j, k in zip([4, 7, 10, 16, 19, 22, 23], [2, 8, 3, 15, 2, 3, 8]):
    hit = any(find_g(CONTACT_SETS[k], set(sub)) is not None for sub in itertools.combinations(sorted(set(T) | {j}), 6))
    check(f"new contact j={j}: T u {{j}} contains an image of paper T_{k} (stage {STAGE[k]} < 5)", hit and STAGE[k] < 5)

# spectral chart A = (inf, v, xv, z), C = (r, 1, 0, w); F = -prod_{i>=1}(A_i - C_pi(i)) (bracket sign)
v, x, z, r, w, s, E, u = sp.symbols("v x z r w s E u")
Asp, Csp = [None, v, x * v, z], [r, 1, 0, w]


def Fp(j):
    p = PERMS[j]
    return -sp.prod([Asp[i] - Csp[p[i]] for i in range(1, 4)])


Dn = v * (z - w)
R = {j: sp.together((Fp(j) - Fp(0)) / Dn) for j in range(24)}
check("R2 = 1 - x, R6 = x(r-1), R8 = r - x", is_zero(R[2] - (1 - x)) and is_zero(R[6] - x * (r - 1)) and is_zero(R[8] - (r - x)))
sdef = (1 - w) * (z - v) / (z - w)
Edef = sdef * (w - v * x) / (v * (w - 1))
check("R5 = -s x", is_zero(R[5] + sdef * x))
check("(3) R1 - R3 = E and R9 - R7 = -rE", is_zero(R[1] - R[3] - Edef) and is_zero(R[9] - R[7] + r * Edef))
wvE = v * (E - s * x) / (v * E - s)
zvE = (w * (s + v) - v) / (s - 1 + w)
onc = lambda e: sp.cancel(sp.together(e.subs(z, zvE).subs(w, wvE)))
check("w = v(E-sx)/(vE-s), z = (w(s+v)-v)/(s-1+w) invert s and E", is_zero(onc(sdef) - s) and is_zero(onc(Edef) - E))
RvE = {j: onc(R[j]) for j in range(24)}
e3c = sp.cancel(RvE[12] - RvE[16] + s * r / v)
e3d = sp.cancel(RvE[16] - (1 - s) * RvE[12])
e3e = sp.cancel((RvE[13] - RvE[7] - RvE[12]) + x * (RvE[12] - RvE[16] + RvE[7]) - (x - 1) * RvE[18])
for nm, e in [("(R12 - R16) + sr/v", e3c), ("R16 - (1-s) R12", e3d), ("(R13-R7-R12) + x(R12-R16+R7) - (x-1)R18", e3e)]:
    check(f"(3) {nm} is a real constant in (x, r, s)", not (e.free_symbols & {v, E}), info=str(sp.factor(e)))
Dc = E * (v - 1) + s * (x - 1)
Qc = sp.factor(v * Dc * (RvE[18] - u))
Bv = r*s*v*x + r*s + r*v**2*x - r*v - s*v**2*x - s*v - 2*v**2*x + v**2 + v*x - u*v*(v - 1)
Cv = -r*s**2*x - r*s*v*x + r*s + s**2*v*x + s*v**2*x**2 - s*v*x**2 + s*v*x - s*v - u*s*v*(x - 1)
Qp = v * (v - r) * E**2 + Bv * E + Cv
check("(4) Q = vD(R18 - u) = v(v-r)E^2 + B(v)E + C(v)", is_zero(Qc - Qp))
check("(5) Q(v, E_D) = -s(r-v)(vx-1)^2(s+v-1)/(v-1)^2", is_zero(Qp.subs(E, -s * (x - 1) / (v - 1)) + s * (r - v) * (v * x - 1)**2 * (s + v - 1) / (v - 1)**2))
check("D(v=1) = s(x-1)", is_zero(Dc.subs(v, 1) - s * (x - 1)))
allowed = [v, v * x - 1, Dc]
bad = []
for j in range(24):
    for fct, _ in sp.factor_list(sp.denom(sp.together(RvE[j])))[1]:
        if fct.free_symbols & {v, E} and not any(sp.cancel(fct / al).free_symbols & {v, E} == set() for al in allowed):
            bad.append((j, fct))
check("every denominator factor of the 24 R_j(v,E) is v, vx-1 or D", not bad, info=str(bad))
check("F0/D = -x(v-1) (full height), on the chart", is_zero(onc(Fp(0) / Dn) + x * (v - 1)))
check("Q(0, E) = rs(E - sx + 1)", is_zero(Qp.subs(v, 0) - r * s * (E - s * x + 1)))
check("D at (v,E) = (0, sx-1) equals 1 - s", is_zero(Dc.subs({v: 0, E: s * x - 1}) - (1 - s)))
check("R3 has denominator D only; R3(0, sx-1) = 1 - sx - x (real)",
      not (sp.cancel(sp.denom(sp.together(RvE[3])) / Dc).free_symbols & {v, E}) and is_zero(RvE[3].subs({v: 0, E: s * x - 1}) - (1 - s * x - x)))
zeta = sp.symbols("zeta")
Qz = sp.expand(zeta**2 * Qp.subs(E, 1 / zeta))
check("zeta = 1/E branch at (v,zeta)=(0,0): Q=0, dQ/dzeta = rs, dQ/dv = -r (smooth, zeta = v/s + O(v^2))",
      Qz.subs({v: 0, zeta: 0}) == 0 and is_zero(sp.diff(Qz, zeta).subs({v: 0, zeta: 0}) - r * s) and is_zero(sp.diff(Qz, v).subs({v: 0, zeta: 0}) + r))
R3z = sp.cancel(RvE[3].subs(E, 1 / zeta))
check("R3 regular at (v,zeta) = (0,0) with value 1 - s - x",
      sp.denom(R3z).subs({v: 0, zeta: 0}) != 0 and is_zero(R3z.subs({v: 0, zeta: 0}) - (1 - s - x)), info=str(sp.factor(R3z)))

# ---------------------------------------------------------------- NUMERIC
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(31)


def prods(aa, bb):
    return np.array([np.prod([aa[i][0] * bb[p[i]][1] - aa[i][1] * bb[p[i]][0] for i in range(4)]) for p in PERMS])


def minphi(gv, restarts=12):
    fo = lambda th: float(gv @ leibniz_weights_np(so4_from_params(th)))
    return min(minimize(fo, rng.uniform(-np.pi, np.pi, 6), method="BFGS").fun for _ in range(restarts))


acc = tried = 0
worst = np.inf
while acc < 40 and tried < 400000:
    tried += 1
    xv = -abs(rng.normal()) * rng.choice([0.3, 1, 3]); rv = -abs(rng.normal()) * rng.choice([0.3, 1, 3])
    sv = rng.uniform(0, 1); uv = rng.normal() * rng.choice([0.3, 1, 3])
    vv = complex(rng.normal(), abs(rng.normal()))
    Bn = complex(Bv.subs({r: rv, s: sv, x: xv, u: uv, v: vv}))
    Cn = complex(Cv.subs({r: rv, s: sv, x: xv, u: uv, v: vv}))
    for En in np.roots([vv * (vv - rv), Bn, Cn]):
        wv = vv * (En - sv * xv) / (vv * En - sv); zv = (wv * (sv + vv) - vv) / (sv - 1 + wv)
        F = prods([(1, 0), (vv, 1), (xv * vv, 1), (zv, 1)], [(rv, 1), (1, 1), (0, 1), (wv, 1)])
        Rv = (F - F[0]) / (F[2] - F[0])
        if max(abs(Rv[T].imag)) > 1e-8 * max(1, abs(Rv).max()):
            continue
        for mu in (1, -1):
            gv = -mu * Rv.imag
            if gv.min() < -1e-12:
                continue
            acc += 1
            worst = min(worst, minphi(gv / np.abs(gv).max()))
check(f"main statement on {acc} realizable sources on the curve (seed 31): min Phi >= -1e-9", worst > -1e-9,
      level="NUMERIC", info=f"worst = {worst:.3e}, tried {tried}")
print("EXACT failures:", len(FAIL))
sys.exit(1 if FAIL else 0)
