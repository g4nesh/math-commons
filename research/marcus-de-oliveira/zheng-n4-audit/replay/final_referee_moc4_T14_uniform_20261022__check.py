"""Exact replay and numeric tests for the supplement note
"A six-contact inequality from a reciprocal half-plane obstruction"
(final_referee_moc4_T14_uniform_20261022), the primary proof of paper case
T_15 = {1, 2, 3, 4, 6, 7}.

Own code only.  Nothing from the supplement is read or executed.  Uses ../lib/moc4lib.py.
Every printed line is one check.  Exit status is nonzero iff an EXACT check fails.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from moc4lib import (P, PERMS, F_VALUES, F_IS_MINOR, F_TABLE, generators, find_g, CONTACT_SETS,
                     leibniz_weights_np, so4_from_params)
import sympy as sp

FAIL = []


def check(name, ok, level="EXACT", info=""):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}" + (f"   {info}" if info else ""), flush=True)
    if not ok and level == "EXACT":
        FAIL.append(name)


def is_zero(e):
    return sp.cancel(sp.together(e)) == 0


T = [1, 2, 3, 4, 6, 7]
FV = F_VALUES
check("T is paper T_15 (identity group element)", find_g(CONTACT_SETS[15], set(T)) is not None)
gens = generators(T)
check("generators(T) = {0,...,8,53,54}", gens == [0, 1, 2, 3, 4, 5, 6, 7, 8, 53, 54], info=str(gens))
check("only f54 is non-minor", [j for j in gens if not F_IS_MINOR[j]] == [54])
check("g16 - g0 is 1 on f1,f4,f7, -1 on f54, 0 on the other generators",
      all(FV[j][16] - FV[j][0] == {1: 1, 4: 1, 7: 1, 54: -1}.get(j, 0) for j in gens))
add = lambda js: [sum(FV[j][s] for j in js) for s in range(24)]
check("f54 + f1 = f0 + f6 + f53", add([54, 1]) == add([0, 6, 53]))
check("f54 + f4 = f2 + f3 + f6 + f53", add([54, 4]) == add([2, 3, 6, 53]))
check("f54 + f7 = f0 + f5 + f8 + f53", add([54, 7]) == add([0, 5, 8, 53]))
cs = sp.symbols("c0:14")
y = sp.Matrix(P) * sp.Matrix(cs)
sol = sp.solve([y[j] for j in T], cs, dict=True)[0]
g = [sp.expand(y[j].subs(sol)) for j in range(24)]
d = g[0] - g[16]
coefs = [(0, g[22]), (2, g[13]), (3, g[12]), (5, g[10]), (6, g[17]), (8, g[8]), (53, g[16]), (54, d)]
check("(2) g = g22 f0 + g13 f2 + g12 f3 + g10 f5 + g17 f6 + g8 f8 + g16 f53 + d f54 (24 coords)",
      all(is_zero(sum(cf * FV[j][s] for j, cf in coefs) - g[s]) for s in range(24)))
check("(3) g11 = b + r + d", is_zero(g[11] - (g[10] + g[17] + d)))
check("(3) g19 = a + s + d", is_zero(g[19] - (g[13] + g[22] + d)))

# spectral chart A = (z,1,0,x), C = (w,v,inf,tv), P_pi = -prod_{pi(i) != 2}(A_i - C_pi(i))
z, x, w, v, t = sp.symbols("z x w v t")
Asp, Csp = [z, 1, 0, x], [w, v, None, t * v]


def Pp(j):
    p = PERMS[j]
    return -sp.prod([Asp[i] - Csp[p[i]] for i in range(4) if p[i] != 2])


R = {j: sp.together((Pp(j) - Pp(1)) / (v * (z - w))) for j in range(24)}
check("(4) R2 = x - t, R3 = t(x-1), R4 = 1 - t",
      is_zero(R[2] - (x - t)) and is_zero(R[3] - t * (x - 1)) and is_zero(R[4] - (1 - t)))
check("(P3 - P1)/(P2 - P4) = t", is_zero((Pp(3) - Pp(1)) / (Pp(2) - Pp(4)) - t))
J = x - t - R[6]
check("J = (x - tv)[z(v+1-w) - v]/(v(z-w))", is_zero(J - (x - t * v) * (z * (v + 1 - w) - v) / (v * (z - w))))
check("R7 + t = tv/(x - tv) J", is_zero(R[7] + t - t * v / (x - t * v) * J))
zz = v / (v + 1 - w)
Rz = {j: sp.cancel(sp.together(R[j].subs(z, zz))) for j in range(24)}
c11 = sp.cancel(Rz[11] - (1 - x) * Rz[10])
c19 = sp.cancel(Rz[19] - (t - 1) / t * Rz[13])
check("(6) R11 - (1-x) R10 is a real constant (gives g11 = (1-x) b)", not (c11.free_symbols & {v, w}), info=str(c11))
check("(6) R19 - (t-1)/t R13 is a real constant (gives g19 = (t-1)/t a)", not (c19.free_symbols & {v, w}), info=str(c19))
p, qq = v / (w - 1), 1 / v
Kx = p * qq - p * t - t + 1
L = ((x - 1) * p - x) / (p * qq - p + 1)
check("K = R10 = pq - pt - t + 1", is_zero(Rz[10] - Kx))
check("z = p/(p-1) and z - w = (pq-p+1)/(qp(p-1))", is_zero(zz - p / (p - 1)) and is_zero(zz - w - (p * qq - p + 1) / (qq * p * (p - 1))))
check("R13 = -t(L+1)", is_zero(Rz[13] + t * (L + 1)))
TT, XX = -t, -x
Hh, Cc, kk = 1 + TT, 1 + XX, TT * XX - 1
U, V = Hh - Kx, XX - L
check("U = -p(T+q), V = p(Xq+1)/(pq-p+1)", is_zero(U + p * (TT + qq)) and is_zero(V - p * (XX * qq + 1) / (p * qq - p + 1)))
D1, D2 = Cc * U + kk, kk - Hh * V
Zz = (TT - XX) * U * V - Hh * XX * U - TT * Cc * V
check("(8) R12 = (T-X) - Z/D1", is_zero(Rz[12] - ((TT - XX) - Zz / D1)))
check("(8) R8 = (T-X) + Z/D2", is_zero(Rz[8] - ((TT - XX) + Zz / D2)))
Tq, Xq, u1, v1, Bq, Aq = sp.symbols("T X u1 v1 B A", real=True)
Uf, Vf = u1 + sp.I * Bq, v1 + sp.I * Aq
Hq, Cq, kq = 1 + Tq, 1 + Xq, Tq * Xq - 1
D1f, D2f = Cq * Uf + kq, kq - Hq * Vf
Zf = (Tq - Xq) * Uf * Vf - Hq * Xq * Uf - Tq * Cq * Vf
check("(9) Z = (X-T)/(CH) D1 D2 - TC/H D1 + XH/C D2",
      is_zero(Zf - ((Xq - Tq) / (Cq * Hq) * D1f * D2f - Tq * Cq / Hq * D1f + Xq * Hq / Cq * D2f)))
m = sp.im(sp.expand(Zf * sp.conjugate(D1f)))      # |D1|^2 Im(Z/D1)
n = -sp.im(sp.expand(Zf * sp.conjugate(D2f)))     # -|D2|^2 Im(Z/D2)
a2 = lambda e: sp.expand(e * sp.conjugate(e))
check("positive-weight identity T C^2 |D1|^2 m + X H^2 |D2|^2 n = (T-X)(T C A |D1|^2 + X H B |D2|^2)",
      sp.expand(Tq * Cq**2 * m + Xq * Hq**2 * n - (Tq - Xq) * (Tq * Cq * Aq * a2(D1f) + Xq * Hq * Bq * a2(D2f))) == 0)
check("Im D1 = C B and Im D2 = -H A", sp.expand(sp.im(D1f) - Cq * Bq) == 0 and sp.expand(sp.im(D2f) + Hq * Aq) == 0)

# ---------------------------------------------------------------- NUMERIC
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(7)
print("f2, f5, f54 =", F_TABLE[2], F_TABLE[5], F_TABLE[54])
Vv = np.array(FV, dtype=float)
ph = lambda j, O: float(Vv[j] @ leibniz_weights_np(O))
m1 = lambda O: np.linalg.det(O[np.ix_([0, 2], [2, 3])])
m2 = lambda O: np.linalg.det(O[np.ix_([0, 1], [1, 3])])
f = lambda th: (lambda O: ph(54, O) + 2 * abs(m1(O) * m2(O)))(so4_from_params(th))
best = min(minimize(f, rng.uniform(-np.pi, np.pi, 6), method="Nelder-Mead",
                    options={"maxiter": 4000, "xatol": 1e-10, "fatol": 1e-13}).fun for _ in range(80))
check("cited SOC (1): Phi_f54 >= -2 |det O[13,34] det O[12,24]| (80 Nelder-Mead starts, seed 7); "
      "with Phi_f2 = det O[13,34]^2, Phi_f5 = det O[12,24]^2 this gives (1)", best > -1e-9, level="NUMERIC", info=f"min = {best:.3e}")
for sgn in (1, -1):
    fs = lambda th: (lambda O: ph(54, O) + sgn * 2 * m1(O) * m2(O))(so4_from_params(th))
    bs = min(minimize(fs, rng.uniform(-np.pi, np.pi, 6), method="BFGS").fun for _ in range(40))
    print(f"[INFO] min Phi_f54 + {sgn}*2 m1 m2 = {bs:.4f} (neither signed bound alone holds)")


def prods(a, b):
    return np.array([np.prod([a[i][0] * b[p[i]][1] - a[i][1] * b[p[i]][0] for i in range(4)]) for p in PERMS])


def minphi(gv, restarts=12):
    fo = lambda th: float(gv @ leibniz_weights_np(so4_from_params(th)))
    return min(minimize(fo, rng.uniform(-np.pi, np.pi, 6), method="BFGS").fun for _ in range(restarts))


acc = tried = nd = 0
worst = worst_ab = worst_TX = np.inf
while acc < 40 and tried < 400000:
    tried += 1
    xv = -abs(rng.normal()) * rng.choice([0.3, 1, 3]); tv = -abs(rng.normal()) * rng.choice([0.3, 1, 3])
    vv = complex(*rng.normal(size=2)); wv = complex(*rng.normal(size=2)); zv = vv / (vv + 1 - wv)
    F = prods([(zv, 1), (1, 1), (0, 1), (xv, 1)], [(wv, 1), (vv, 1), (1, 0), (tv * vv, 1)])
    Rv = (F - F[1]) / (F[4] - F[1])
    if max(abs(Rv[T].imag)) > 1e-8 * max(1, abs(Rv).max()):
        continue
    for mu in (1, -1):
        gv = -mu * Rv.imag
        if gv.min() < -1e-12:
            continue
        dv = gv[0] - gv[16]
        if dv > 0:
            nd += 1
            worst_ab = min(worst_ab, (gv[13] * gv[10] - dv * dv) / dv**2)
            worst_TX = min(worst_TX, -tv - (-xv))
        acc += 1
        worst = min(worst, minphi(gv / np.abs(gv).max()))
check(f"main statement on {acc} realizable sources from the chart (seed 7): min Phi >= -1e-9", worst > -1e-9,
      level="NUMERIC", info=f"worst = {worst:.3e}, tried {tried}")
check(f"ab >= d^2 and T >= X on the {nd} sampled sources with d > 0", nd == 0 or (worst_ab >= 0 and worst_TX >= 0),
      level="NUMERIC", info=f"min (ab-d^2)/d^2 = {worst_ab:.3g}, min T-X = {worst_TX:.3g}")
# the cited six-point-circle stratum (v real): sample it directly, J real chosen, z solved
acc = tried = 0
worst = np.inf
worst_ab_c = np.inf
while acc < 40 and tried < 400000:
    tried += 1
    xv, tv, vv, j0 = rng.normal(size=4) * rng.choice([0.3, 1, 3], size=4)
    wv = complex(*rng.normal(size=2))
    zv = (vv - j0 * vv * wv / (xv - tv * vv)) / ((vv + 1 - wv) - j0 * vv / (xv - tv * vv))
    F = prods([(zv, 1), (1, 1), (0, 1), (xv, 1)], [(wv, 1), (vv, 1), (1, 0), (tv * vv, 1)])
    Rv = (F - F[1]) / (F[4] - F[1])
    if max(abs(Rv[T].imag)) > 1e-8 * max(1, abs(Rv).max()):
        continue
    for mu in (1, -1):
        gv = -mu * Rv.imag
        if gv.min() < -1e-12:
            continue
        acc += 1
        dv = gv[0] - gv[16]
        if dv > 0:
            worst_ab_c = min(worst_ab_c, (gv[13] * gv[10] - dv * dv) / dv**2)
        worst = min(worst, minphi(gv / np.abs(gv).max()))
check(f"six-circle stratum (v real) of T_15: {acc} sampled realizable sources (seed 7 stream), min Phi >= -1e-9",
      worst > -1e-9, level="NUMERIC", info=f"worst = {worst:.3e}, tried {tried}, min (ab-d^2)/d^2 over d>0 = {worst_ab_c:.3g}")
print("EXACT failures:", len(FAIL))
sys.exit(1 if FAIL else 0)
