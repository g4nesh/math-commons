"""Exact replay for blind_shared_coreN_uniform_20261023 (seven-element contact set
N = {0,3,7,11,12,16,20} = Table 5 set no. 24).

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
N = frozenset({0, 3, 7, 11, 12, 16, 20})
check("N is Table 5 set no. 24", L.AUX7[24] == N)
gens = [2, 8, 12, 15, 29, 39, 30]
check("generators(N) = {f2,f8,f12,f15,f29,f39,f30}", sorted(L.generators(N)) == sorted(gens), L.generators(N))
check("f2,f8,f12,f15,f29,f39 minor; f30 non-minor",
      all(L.F_IS_MINOR[j] for j in gens[:6]) and not L.F_IS_MINOR[30], [L.f_name(j) for j in gens])
check("P has rank 14", L.rank(L.P) == 14)
check("rows of P at N have rank 7", L.rank([L.P[s] for s in sorted(N)]) == 7)
check("the 7 generator vectors are independent (so span the 7-dim kernel)",
      L.rank([L.F_VECS[j] for j in gens]) == 7)
rows = [13, 14, 6, 5, 10, 1, 23]
Mrows = [[L.F_VALUES[j][r] for j in gens] for r in rows]
check("rows [13,14,6,5,10,1,23] of the generators = diag(1,1,1,1,1,1,2)",
      Mrows == [[(2 if (i == j == 6) else (1 if i == j else 0)) for j in range(7)] for i in range(7)], Mrows)

# ---------------------------------------------------------------- B. cubic along the cone
names = ('a', 'b', 'c', 'e', 'f', 'g', 'd', 'z')
ctx = flint.fmpz_mpoly_ctx.get(names, 'lex')
a, b, c, e, f, g, d, z = ctx.gens()
zero = ctx.from_dict({})
coef = {2: a, 8: b, 12: c, 15: e, 29: f, 39: g, 30: d}
w = [zero for _ in range(14)]
for j, s in coef.items():
    for k in range(14):
        assert L.F_VECS[j][k].denominator == 1
        w[k] = w[k] + int(L.F_VECS[j][k]) * s
m = w[:]; m[13] = m[13] + z
Fz = L.cubic_C(m)
zc = zcoeffs(Fz, 7, 8)
check("F(w + z kappa) has z-degree 1", set(zc) <= {0, 1}, sorted(zc))
Cp = ctx.from_dict({k: v for k, v in zc[0].items()})
Qp = ctx.from_dict({k: v for k, v in zc[1].items()})
check("C has 38 monomials, all coefficients positive", len(Cp) == 38 and all(int(x) > 0 for x in Cp.coeffs()), len(Cp))
check("[d^3]C = 1", Cp.to_dict().get((0,)*6 + (3, 0), 0) == 1)
check("Q has 19 monomials, all coefficients positive", len(Qp) == 19 and all(int(x) > 0 for x in Qp.coeffs()), len(Qp))
check("[d^2]Q = 1", (0,)*6 + (2, 0) in Qp.monoms() and Qp.to_dict()[(0,)*6 + (2, 0)] == 1)

# ---------------------------------------------------------------- C. the 14x14 Hessian determinant
TT = third_tensor()
H = hessian(TT, m, zero)
t0 = time.time()
D, nd = bareiss_det(H)
print("INFO det H(w + z kappa): %d terms, %d exact divisions, %.1f s" % (len(D), nd, time.time() - t0))
Dk = zcoeffs(D, 7, 8)
check("det H(w + z kappa) has z-degree <= 4", max(Dk) <= 4, sorted(Dk))
Dk = {k: ctx.from_dict(v) for k, v in Dk.items()}
Cn = -Cp
LHS = zero
for k, Pk in Dk.items():
    LHS = LHS + Pk * Cn ** k * Qp ** (4 - k)
S = None
try:
    S = LHS / (-4 * Qp ** 2)
    check("Q^4 det A is divisible by -4 Q^2", True)
except Exception as ex:
    check("Q^4 det A is divisible by -4 Q^2", False, ex)
G = None
if S is not None:
    try:
        G = S.sqrt()
        check("Q^4 det A / (-4 Q^2) is a perfect square G^2 (identity (4))", G * G == S)
    except Exception as ex:
        check("Q^4 det A / (-4 Q^2) is a perfect square G^2 (identity (4))", False, ex)
if G is not None:
    gd = G.to_dict()
    d9 = (0,)*6 + (9, 0)
    if gd.get(d9, 0) < 0:
        G = -G; gd = G.to_dict()
    check("[d^9]G = 1", gd.get(d9, 0) == 1, gd.get(d9, 0))
    check("G has 498 terms", len(G) == 498, len(G))
    check("G homogeneous of degree 9", all(sum(mn) == 9 for mn in gd))
    check("individual degrees of G in (a,b,c,e,f,g,d) = (2,2,2,2,2,2,9)",
          tuple(max(mn[i] for mn in gd) for i in range(7)) == (2, 2, 2, 2, 2, 2, 9),
          tuple(max(mn[i] for mn in gd) for i in range(7)))
    # coefficient signs by power of d
    bypow = {}
    for mn, cf in gd.items():
        bypow.setdefault(mn[6], []).append(int(cf))
    check("coefficients of d^k, k >= 4: all positive", all(all(x > 0 for x in bypow[k]) for k in bypow if k >= 4))
    check("coefficients of d^k, k <= 2: all negative", all(all(x < 0 for x in bypow[k]) for k in bypow if k <= 2))
    print("INFO d^3 coefficient signs: +%d -%d" % (sum(x > 0 for x in bypow.get(3, [])), sum(x < 0 for x in bypow.get(3, []))))
    # R = d G_d - 3G
    R = d * G.derivative('d') - 3 * G
    rd = R.to_dict()
    check("R = d G_d - 3 G: 427 monomials, all positive", len(R) == 427 and all(int(x) > 0 for x in rd.values()), len(R))
    check("R contains 6 d^9", rd.get(d9, 0) == 6)
    check("every monomial of R has d-degree >= 4 or d-degree <= 2 (d^3 cancels)", all(mn[6] != 3 for mn in rd))
    # G(x,0) = -16 abcefg C(x,0)
    G0 = G.subs({'d': 0})
    C0 = Cp.subs({'d': 0})
    check("G(x,0) = -16 a b c e f g C(x,0)", G0 == -16 * a * b * c * e * f * g * C0)
    check("C(x,0) has 16 monomials, all positive cubic", len(C0) == 16 and all(int(x) > 0 for x in C0.coeffs())
          and all(sum(mn) == 3 for mn in C0.monoms()), len(C0))
    # sample point (10)
    one = {nm: 1 for nm in names[:7]}
    Cv, Qv, Gv = Cp.subs(one), Qp.subs(one), G.subs(one)
    check("at a=..=d=1: C = 38, Q = 19, G = 0", (int(Cv.coeffs()[0]) if not Cv.is_zero() else 0, int(Qv.coeffs()[0]), Gv.is_zero()) == (38, 19, True))
    # numeric matrix A at that point, exact Fractions
    zval = Fr(-38, 19)
    mv = [Fr(sum(int(L.F_VECS[j][k]) for j in gens)) for k in range(14)]
    mv[13] += zval
    Hv = [[Fr(0)] * 14 for _ in range(14)]
    for (i, j, k), v in TT.items():
        Hv[i][j] += v * mv[k]
    I = list(range(2, 14))
    detAI = L.det([[Hv[i][j] for j in I] for i in I])
    check("at that point z = -2 and det A_I = -49152 (I = {2..13})", zval == -2 and detAI == -49152, detAI)
    check("at that point det A = 0 and rank A = 12", L.det(Hv) == 0 and L.rank(Hv) == 12, L.rank(Hv))
    # three border 13x13 determinants
    for (rs, cs, lab) in (([1] + I, [1] + I, "({1}+I,{1}+I)"), ([0] + I, [0] + I, "({0}+I,{0}+I)"),
                          ([1] + I, [0] + I, "({1}+I,{0}+I)")):
        sub = [[H[i][j] for j in cs] for i in rs]
        Db, ndb = bareiss_det(sub)
        bk = zcoeffs(Db, 7, 8)
        okdeg = max(bk) <= 2 if bk else True
        num = zero
        for k, v in bk.items():
            num = num + ctx.from_dict(v) * Cn ** k * Qp ** (2 - k)
        try:
            qq = num / G
            ok = True
        except Exception:
            ok = False
        check("border %s: z-degree <= 2 and Q^2 det(..)|_{z=-C/Q} divisible by G" % lab, okdeg and ok,
              "z-deg %s, %d terms, %d divisions" % (max(bk) if bk else '-', len(num), ndb))
    # identity (11): second d-derivative of -4G^2/Q^2 at a root is -8 G_d^2/Q^2 (sympy check of the calculus)
    t, g1, g2, q0, q1, q2 = sympy.symbols('t g1 g2 q0 q1 q2')
    Gt, Qt = g1 * t + g2 * t ** 2, q0 + q1 * t + q2 * t ** 2      # G(d0) = 0 at t = 0
    f2 = sympy.diff(-4 * Gt ** 2 / Qt ** 2, t, 2).subs(t, 0)
    check("d^2/dd^2 (-4G^2/Q^2) at a root of G equals -8 G_d^2 / Q^2", sympy.simplify(f2 + 8 * g1 ** 2 / q0 ** 2) == 0)

# ---------------------------------------------------------------- D. numeric sanity: positive root and inertia
import numpy as np
rng = np.random.default_rng(20261125)
if G is not None:
    Gd = {tuple(int(u) for u in mn): int(cf) for mn, cf in G.to_dict().items()}
    worst_minor = None
    CpD = {tuple(int(u) for u in mn): int(cf) for mn, cf in Cp.to_dict().items()}
    QpD = {tuple(int(u) for u in mn): int(cf) for mn, cf in Qp.to_dict().items()}
    bad = 0
    for trial in range(40):
        x = rng.exponential(1.0, 6)
        if trial % 4 == 0:
            x[rng.integers(0, 6)] = 0.0     # a zero wall
        # univariate polynomial in d
        cd = np.zeros(10)
        for mn, cf in Gd.items():
            cd[mn[6]] += float(cf) * np.prod([x[i] ** mn[i] for i in range(6)])
        roots = [r.real for r in np.roots(cd[::-1]) if abs(r.imag) < 1e-9 and r.real > 1e-12]
        if len(roots) > 1:
            bad += 1
        for dr in roots:
            vals = dict(zip(names[:7], list(x) + [dr]))
            Cf = sum(float(cf) * np.prod([vals[names[i]] ** mn[i] for i in range(7)]) for mn, cf in CpD.items())
            Qf = sum(float(cf) * np.prod([vals[names[i]] ** mn[i] for i in range(7)]) for mn, cf in QpD.items())
            mvf = np.zeros(14)
            for j, s in zip(gens, list(x) + [dr]):
                mvf += s * np.array([float(t) for t in L.F_VECS[j]])
            mvf[13] -= Cf / Qf
            Hf = np.zeros((14, 14))
            for (i, j, k), v in TT.items():
                Hf[i, j] += v * mvf[k]
            ev = np.linalg.eigvalsh(Hf)
            sc = max(abs(ev))
            npos, nneg = int((ev > 1e-9 * sc).sum()), int((ev < -1e-9 * sc).sum())
            if (npos, nneg) not in ((7, 5), (5, 7)):
                worst_minor = (list(x), dr, npos, nneg)
    print("NUMERIC 40 random x >= 0 (seed 20261125, every 4th with one zero wall): >1 positive root in %d cases; "
          "inertia at roots not (7,5,2)/(5,7,2) in: %s" % (bad, worst_minor))


# stratum test: a realizable point with contacts containing N must have d = g_23/2 = 0
resN = []
for seed in range(20):
    s_ = L.random_point_on_contact_stratum(sorted(N), seed=seed, tries=8, tol=1e-10)
    if s_ is None: continue
    gg = np.real(s_['gaps'])
    resN.append(gg[23] / np.mean(gg[[k for k in range(24) if k not in N]]))
print("NUMERIC realizable points on the N stratum (lib least squares, seeds 0..19): %d found; max |g_23|/mean off-N gap = %.2e"
      % (len(resN), max(abs(x) for x in resN) if resN else float('nan')))

if FAIL:
    print("FAILED:", FAIL); sys.exit(1)
print("ALL EXACT CHECKS PASSED")
