"""Exact replay for root_shared_coreO_uniform_20261023 (seven-element contact set
O = {0,1,2,7,14,16,23} = Table 5 set no. 15).

Everything is rebuilt from the paper's definitions through the shared library (P, f_j,
the cubic C of (A.2)).  Nothing from the supplement is executed, imported or read.
Polynomial work uses python-flint fmpz_mpoly; determinants by fraction-free Bareiss
(full pivoting, every division exact or the run aborts).
Run: python -I -B <this file>.  Exits 1 if an EXACT check fails.  About 30 s.
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint, sympy

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""), flush=True)
    if not ok:
        FAIL.append(name)

# ---------------------------------------------------------------- helpers
def third_tensor():
    ms, Cexpr = L.cubic_sympy()
    T = {}
    for i in range(14):
        for j in range(i, 14):
            dij = sympy.diff(Cexpr, ms[i], ms[j])
            for k in range(14):
                v = sympy.diff(dij, ms[k])
                if v != 0:
                    T[(i, j, k)] = int(v)
                    T[(j, i, k)] = int(v)
    return T

def hessian(T, m, zero):
    H = [[zero for _ in range(len(m))] for _ in range(14)]
    for (i, j, k), v in T.items():
        H[i][j] = H[i][j] + v * m[k]
    return H

def bareiss_det(M):
    n = len(M); A = [r[:] for r in M]; sign = 1; prev = None; ndiv = 0
    for k in range(n - 1):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                if not A[i][j].is_zero():
                    s = len(A[i][j])
                    if best is None or s < best[0]:
                        best = (s, i, j)
        if best is None:
            return A[0][0] * 0, ndiv
        _, pi, pj = best
        if pi != k:
            A[k], A[pi] = A[pi], A[k]; sign = -sign
        if pj != k:
            for r in A:
                r[k], r[pj] = r[pj], r[k]
            sign = -sign
        p = A[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = p * A[i][j] - A[i][k] * A[k][j]
                if prev is not None:
                    num = num / prev          # DomainError unless exact
                    ndiv += 1
                A[i][j] = num
        prev = p
    d = A[n - 1][n - 1]
    return (d if sign == 1 else -d), ndiv

def zcoeffs(D, zi, nv):
    """Split a polynomial into coefficients of powers of variable zi."""
    out = {}
    for mon, c in D.to_dict().items():
        k = mon[zi]
        mon2 = tuple(0 if t == zi else mon[t] for t in range(nv))
        out.setdefault(k, {})
        out[k][mon2] = out[k].get(mon2, 0) + int(c)
    return out

# ---------------------------------------------------------------- A. combinatorics
OC = frozenset({0, 1, 2, 7, 14, 16, 23})
check("O is Table 5 set no. 15", L.AUX7[15] == OC)
check("generators(O) = {f2,f5,f11,f12,f15,f20,f23,f26}", L.generators(OC) == [2, 5, 11, 12, 15, 20, 23, 26], L.generators(OC))
gens = [2, 5, 11, 12, 15, 23, 26]
check("f2,f5,f11,f12,f15,f23 minor; f20, f26 non-minor",
      all(L.F_IS_MINOR[j] for j in gens[:6]) and not L.F_IS_MINOR[20] and not L.F_IS_MINOR[26],
      [L.f_name(j) for j in gens + [20]])
check("P has rank 14; rows at O have rank 7", L.rank(L.P) == 14 and L.rank([L.P[s] for s in sorted(OC)]) == 7)
check("[f2,f5,f11,f12,f15,f23,f26] independent (span the 7-dim kernel)", L.rank([L.F_VECS[j] for j in gens]) == 7)
r1 = [15, 10, 17, 8, 4, 9]
check("rows [15,10,17,8,4,9] = first six rows of I_7",
      [[L.F_VALUES[j][r] for j in gens] for r in r1] == [[1 if i == k else 0 for k in range(7)] for i in range(6)])
r2 = [13, 20, 12, 22, 18, 3]
check("rows [13,20,12,22,18,3] read a+d, b+d, q+d, c+d, e+d, v+d",
      [[L.F_VALUES[j][r] for j in gens] for r in r2] == [[1 if (k == i or k == 6) else 0 for k in range(7)] for i in range(6)])
lhs = [L.F_VECS[20][k] + L.F_VECS[26][k] for k in range(14)]
rhs = [sum(L.F_VECS[j][k] for j in gens[:6]) for k in range(14)]
check("circuit (2): f20 + f26 = f2 + f5 + f11 + f12 + f15 + f23", lhs == rhs)
# group facts
def gmap_f(g):
    """induced permutation of the 82 generators (index -> index)."""
    out = {}
    for j in range(82):
        img = tuple(L.g_act_values(g, list(L.F_VALUES[j])))
        out[j] = L.F_VALUES.index(img)
    return out
idp = (0, 1, 2, 3)
gT = (idp, idp, 1)
mT = gmap_f(gT)
check("transpose (rho = rho' = id, eps = 1) stabilizes O", L.g_act_set(gT, OC) == OC)
check("transpose exchanges f20, f26 and acts 2<->5, 11<->12, 15<->23",
      (mT[20], mT[26], mT[2], mT[5], mT[11], mT[12], mT[15], mT[23]) == (26, 20, 5, 2, 12, 11, 23, 15))
stab = [g for g in L.G_ELEMENTS if L.g_act_set(g, OC) == OC]
perms_ace = {}
for g in stab:
    mg = gmap_f(g)
    if mg[26] == 26 and set(mg[j] for j in gens[:6]) == set(gens[:6]) and set(mg[j] for j in (2, 12, 15)) == {2, 12, 15}:
        perms_ace[(mg[2], mg[12], mg[15])] = g
print("INFO |Stab(O)| = %d" % len(stab))
check("elements fixing O and f26 induce all 6 permutations of (f2,f12,f15) and permute the 6 minors",
      len(perms_ace) == 6, sorted(perms_ace))
for (i, j) in ((2, 12), (2, 15), (12, 15)):
    hit = None
    for g in L.G_ELEMENTS:
        mg = gmap_f(g)
        if mg[22] == 26 and {mg[1], mg[9]} == {i, j}:
            hit = g; break
    check("some g maps the SOC reference triple (f1,f9;f22) to (f%d,f%d;f26)" % (i, j), hit is not None, hit)

# ---------------------------------------------------------------- B. cubic along the cone
names = ('a', 'b', 'q', 'c', 'e', 'v', 'd', 'z')
ctx = flint.fmpz_mpoly_ctx.get(names, 'lex')
a, b, q, c, e, v, d, z = ctx.gens()
zero = ctx.from_dict({})
coef = {2: a, 5: b, 11: q, 12: c, 15: e, 23: v, 26: d}
w = [zero for _ in range(14)]
for j, s in coef.items():
    for k in range(14):
        assert L.F_VECS[j][k].denominator == 1
        w[k] = w[k] + int(L.F_VECS[j][k]) * s
m = w[:]; m[13] = m[13] + z
zc = zcoeffs(L.cubic_C(m), 7, 8)
check("F(w + z kappa) has z-degree 1", set(zc) <= {0, 1}, sorted(zc))
Cp = ctx.from_dict(zc[0]); Qp = ctx.from_dict(zc[1])
check("C = ace + (b+d)(q+d)(v+d)  (eq. (2) of the cubic display)", Cp == a * c * e + (b + d) * (q + d) * (v + d))
Qdoc = -a*b + a*c + a*e + b*q + b*v + b*d - q*c + q*v + q*d + c*e - e*v + v*d + d**2
check("Q = -ab+ac+ae+bq+bv+bd-qc+qv+qd+ce-ev+vd+d^2  (eq. (6))", Qp == Qdoc)

# ---------------------------------------------------------------- C. Hessian determinant (7)
TT = third_tensor()
H = hessian(TT, m, zero)
t0 = time.time()
D, nd = bareiss_det(H)
print("INFO det H(w + z kappa): %d terms, %d exact divisions, %.1f s" % (len(D), nd, time.time() - t0))
Dk = zcoeffs(D, 7, 8)
check("det H(w + z kappa) has z-degree <= 4", max(Dk) <= 4, sorted(Dk))
LHS = zero
for k, v_ in Dk.items():
    LHS = LHS + ctx.from_dict(v_) * (-Cp) ** k * Qp ** (4 - k)
G = None
try:
    S = LHS / (-4 * d ** 2 * Qp ** 2)
    G = S.sqrt()
    check("identity (7): Q^4 det H(g - C kappa/Q) = -4 d^2 Q^2 G^2", G * G == S)
except Exception as ex:
    check("identity (7): Q^4 det H(g - C kappa/Q) = -4 d^2 Q^2 G^2", False, ex)
if G is not None:
    gd = {tuple(int(u) for u in mn): int(cf) for mn, cf in G.to_dict().items()}
    d8 = (0,) * 6 + (8, 0)
    if gd.get(d8, 0) < 0:
        G = -G; gd = {k_: -v_ for k_, v_ in gd.items()}
    check("[d^8]G = 1", gd.get(d8, 0) == 1)
    check("G has 415 monomials, homogeneous of degree 8", len(gd) == 415 and all(sum(mn) == 8 for mn in gd), len(gd))
    check("individual degrees (2,2,2,2,2,2,8)", tuple(max(mn[i] for mn in gd) for i in range(7)) == (2, 2, 2, 2, 2, 2, 8),
          tuple(max(mn[i] for mn in gd) for i in range(7)))
    # (10)
    G00 = G.subs({'c': 0, 'e': 0})
    check("(10): G(a,b,q,0,0,v,d) = (b+d)(a^2 b+abd+ad^2+bd^2+d^3)(v+d)^2(q+d)^2",
          G00 == (b + d) * (a**2 * b + a*b*d + a*d**2 + b*d**2 + d**3) * (v + d)**2 * (q + d)**2)
    # charts, d = 1; variables X, Y, Z, b, q, v
    cx = flint.fmpz_mpoly_ctx.get(('X', 'Y', 'Z', 'b', 'q', 'v'), 'lex')
    X, Y, Z, B, Qv, V = cx.gens()
    one_ = cx.from_dict({(0,)*6: 1})
    iA, iC, iE = 0, 3, 4                                  # positions of a, c, e
    def balanced(k_ace):
        Dn = one_ + X + Y + Z
        numA, numC, numE = Y + Z, X + Z, X + Y
        bal = cx.from_dict({})
        for mn, cf in gd.items():
            i, j, l = mn[iA], mn[iC], mn[iE]
            bal += cf * powr(numA, i) * powr(numC, j) * powr(numE, l) * powr(Dn, k_ace - i - j - l) \
                   * powr(B, mn[1]) * powr(Qv, mn[2]) * powr(V, mn[5])
        return {k_: int(v_) for k_, v_ in bal.to_dict().items()}
    def dominant(kL, kY, kZ):
        Lx = one_ + X
        dom = cx.from_dict({})
        for mn, cf in gd.items():
            i, j, l = mn[iA], mn[iC], mn[iE]
            dom += cf * powr(one_ + Y * Lx ** 2, i) * powr(Z, j) * powr(Lx, kL - i - j - l) * powr(one_ + Y, kY - i) \
                   * powr(one_ + Z, kZ - j - l) * powr(B, mn[1]) * powr(Qv, mn[2]) * powr(V, mn[5])
        return {k_: int(v_) for k_, v_ in dom.to_dict().items()}
    def powr(p_, k_):
        assert k_ >= 0
        return one_ if k_ == 0 else p_ ** k_
    kmin = max(mn[iA] + mn[iC] + mn[iE] for mn in gd)
    for kk in sorted({kmin, 6}):
        bd = balanced(kk)
        check("balanced chart (8): D^%d G(a,b,q,c,e,v,1) has %d coefficients, all positive, constant 1 (doc: 2259)" % (kk, len(bd)),
              all(x > 0 for x in bd.values()) and bd.get((0,)*6, 0) == 1, len(bd))
    kYm = max(mn[iA] for mn in gd); kZm = max(mn[iC] + mn[iE] for mn in gd)
    for (kL, kY, kZ) in sorted({(kmin, kYm, kZm), (6, 2, 4)}):
        dd = dominant(kL, kY, kZ)
        check("dominant chart (9): L^%d (1+Y)^%d (1+Z)^%d G has %d coefficients, all positive (doc: 3096)" % (kL, kY, kZ, len(dd)),
              all(x > 0 for x in dd.values()), len(dd))
        check("  ... contains the monomial 8X; constant term 0", dd.get((1, 0, 0, 0, 0, 0), 0) == 8 and dd.get((0,)*6, 0) == 0,
              (dd.get((1, 0, 0, 0, 0, 0)), dd.get((0,)*6, 0)))

# chart inverses (sympy)
Xs, Ys, Zs, as_, cs_, es_ = sympy.symbols('X Y Z a c e', positive=True)
Ds = 1 + Xs + Ys + Zs
A8, C8, E8 = (Ys + Zs) / Ds, (Xs + Zs) / Ds, (Xs + Ys) / Ds
inv8 = [(C8 + E8 - A8) / (2 - A8 - C8 - E8), (A8 + E8 - C8) / (2 - A8 - C8 - E8), (A8 + C8 - E8) / (2 - A8 - C8 - E8)]
check("chart (8): stated inverse recovers (X,Y,Z)", [sympy.simplify(u - t_) for u, t_ in zip(inv8, (Xs, Ys, Zs))] == [0, 0, 0])
check("chart (8): a+c+e = 2(X+Y+Z)/D < 2 and c+e-a = 2X/D etc.",
      sympy.simplify(A8 + C8 + E8 - 2 * (Xs + Ys + Zs) / Ds) == 0 and sympy.simplify(C8 + E8 - A8 - 2 * Xs / Ds) == 0)
Ls = 1 + Xs
A9, C9, E9 = (1 + Ys * Ls**2) / (Ls * (1 + Ys)), Zs / (Ls * (1 + Zs)), 1 / (Ls * (1 + Zs))
S9 = C9 + E9
check("chart (9): S = 1/L, X = 1/S - 1, Y = (aS - S^2)/(1 - aS), Z = c/e",
      sympy.simplify(S9 - 1 / Ls) == 0 and sympy.simplify(1 / S9 - 1 - Xs) == 0
      and sympy.simplify((A9 * S9 - S9**2) / (1 - A9 * S9) - Ys) == 0 and sympy.simplify(C9 / E9 - Zs) == 0)
# forward direction of (9): for a >= S > 0, aS < 1, e > 0 the parameters are X > 0, Y >= 0, Z >= 0
check("chart (9) domain: a >= S, 0 < S, aS < 1 imply S < 1 (so X > 0) and Y >= 0",
      True, "PROOF-READ: S^2 <= aS < 1; aS - S^2 = S(a - S) >= 0 and 1 - aS > 0")

# ---------------------------------------------------------------- D. numeric tests
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(20261125)
Pn = L.P_numpy()
def phi_of(j, th):
    return float(L.Phi_np(L.F_VECS[j], L.so4_from_params(th)))
worst = {}
for (i, j) in ((2, 12), (2, 15), (12, 15), (1, 9)):
    k = 22 if (i, j) == (1, 9) else 26
    fobj = lambda th: phi_of(k, th) + 2 * np.sqrt(max(phi_of(i, th), 0) * max(phi_of(j, th), 0))
    best = np.inf
    for r in range(60):
        th0 = rng.uniform(-np.pi, np.pi, 6)
        res = minimize(fobj, th0, method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 4000})
        best = min(best, res.fun)
    worst[(i, j, k)] = best
print("NUMERIC SOC atoms: min over SO(4) of Phi_fk + 2 sqrt(Phi_fi Phi_fj), 60 Nelder-Mead restarts each, seed 20261125:",
      {str(kk): float(round(vv, 12)) for kk, vv in worst.items()})
# H3 bound (5): a f2 + c f12 + e f15 + H3(a,c,e) f26
wH = np.inf
for t in range(25):
    acev = rng.exponential(1.0, 3)
    if t % 5 == 0:
        acev[0] = acev[1] + acev[2] + rng.exponential(1.0)
    s_, m_ = acev.sum(), acev.max()
    h3 = s_ / 2 if 2 * m_ <= s_ else np.sqrt(m_ * (s_ - m_))
    cvec = acev[0] * np.array([float(x) for x in L.F_VECS[2]]) + acev[1] * np.array([float(x) for x in L.F_VECS[12]]) \
        + acev[2] * np.array([float(x) for x in L.F_VECS[15]]) + h3 * np.array([float(x) for x in L.F_VECS[26]])
    val, _ = L.min_phi_over_SO4(list(cvec), restarts=12, seed=int(t))
    wH = min(wH, val / (s_ + h3))
print("NUMERIC H3 bound (5): 25 random (a,c,e) (5 dominant), 12 BFGS restarts each: worst normalized min Phi = %.3e" % wH)


# stratum test: realizable points with contacts containing O must have d <= 0 after the sign normalization, or H3 >= d
resO = []
for seed in range(20):
    s_ = L.random_point_on_contact_stratum(sorted(OC), seed=seed, tries=8, tol=1e-10)
    if s_ is None: continue
    gg = np.real(s_['gaps']); sc = np.mean(gg[[k for k in range(24) if k not in OC]])
    a_, b_, q_, c_, e_, v_ = [gg[r] for r in (15, 10, 17, 8, 4, 9)]
    d_ = gg[13] - a_
    if d_ < 0:   # transposed form (3): transpose sends f5->f2, f11->f12, f23->f15, f20->f26
        a_, c_, e_, d_ = b_ + d_, q_ + d_, v_ + d_, -d_
    acev = np.array([a_, c_, e_]); s3, m3 = acev.sum(), acev.max()
    h3 = s3 / 2 if 2 * m3 <= s3 else np.sqrt(max(m3 * (s3 - m3), 0))
    resO.append(((h3 - d_) / sc, d_ / sc))
print("NUMERIC realizable points on the O stratum (seeds 0..19): %d found; d/mean-gap values %s; min (H3 - d)/mean-gap = %.3e"
      % (len(resO), [round(float(r[1]), 4) for r in resO], min(r[0] for r in resO) if resO else float('nan')))

if FAIL:
    print("FAILED:", FAIL); sys.exit(1)
print("ALL EXACT CHECKS PASSED")
