"""Exact replay and numeric tests for root_moc4_T46_whole_assembly_20261023 (paper T_47). Own code; imports only ../lib/moc4lib.py."""
# ---- common block (own code; identical in the C20-T46-T47 replay scripts) ----
import os, sys, time, itertools
from fractions import Fraction
from math import lcm, comb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint
import sympy

FAILS = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + str(info)) if info != "" else ""), flush=True)
    if not ok:
        FAILS.append(name)


def contains_image(Tsmall, Tbig):
    """some g in G with g.Tsmall a subset of Tbig"""
    Tb = frozenset(Tbig)
    for g, t in zip(L.G_ELEMENTS, L._G_MAPS):
        img = frozenset(t[s] for s in Tsmall)
        if img <= Tb:
            return g
    return None


def common_zero(gens):
    z = set(range(24))
    for j in gens:
        z &= set(L.F_ZERO[j])
    return sorted(z)


def face_setup(names, gens):
    """m = sum coef_j f_j + h kappa as linear fmpz_mpoly forms; returns ctx, m, C0, Q, Hessian"""
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names) + ("h",), "lex")
    V = ctx.gens()
    m = []
    for k in range(14):
        e = ctx.from_dict({})
        for v, j in zip(V, gens):
            x = Fraction(L.F_VECS[j][k])
            assert x.denominator == 1
            if x:
                e += int(x) * v
        if k == 13:
            e += V[-1]
        m.append(e)
    C = L.cubic_C(m)
    assert C.degrees()[-1] == 1
    C0 = C.subs({"h": 0})
    Q = C.derivative("h")
    ms, Cs = L.cubic_sympy()
    H = [[None] * 14 for _ in range(14)]
    for i in range(14):
        for k in range(14):
            poly = sympy.Poly(sympy.diff(Cs, ms[i], ms[k]), *ms)
            e = ctx.from_dict({})
            for mon, co in poly.terms():
                if co != 0:
                    e += int(co) * m[mon.index(1)]
            H[i][k] = e
    return ctx, m, C0, Q, H


def bareiss_det(ctx, H):
    M = [row[:] for row in H]
    n = len(M)
    prev = ctx.from_dict({(0,) * ctx.nvars(): 1})
    sign = 1
    for k in range(n - 1):
        if M[k][k].is_zero():
            sw = next((r for r in range(k + 1, n) if not M[r][k].is_zero()), None)
            if sw is None:
                return ctx.from_dict({})
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                q, r = divmod(M[i][j] * M[k][k] - M[i][k] * M[k][j], prev)
                assert r.is_zero()
                M[i][j] = q
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def radial_numerator(ctx, D, C0, Q):
    """Q^k det H(g - (C0/Q) kappa), k = h-degree of det"""
    hd = D.degrees()[-1]
    hv = ctx.gens()[-1]
    coeffs, rest = [], D
    for k in range(hd + 1):
        ck = rest.subs({"h": 0})
        coeffs.append(ck)
        if k < hd:
            q, r = divmod(rest - ck, hv)
            assert r.is_zero()
            rest = q
    N = ctx.from_dict({})
    for k, ck in enumerate(coeffs):
        N += ck * (-C0) ** k * Q ** (hd - k)
    return hd, N


class Chart:
    def __init__(self, CV):
        self.CV = CV
        self.ctx = flint.fmpz_mpoly_ctx.get(CV, "lex")
        self.one = self.ctx.from_dict({(0,) * len(CV): 1})

    def var(self, x):
        return RF(self, self.ctx.gens()[self.CV.index(x)])

    def const(self, c):
        return RF(self, c * self.one)

    def strict_monomial(self, poly, strict, sign=None):
        idx = [self.CV.index(v) for v in strict]
        for mon, co in zip(poly.monoms(), poly.coeffs()):
            if (sign is None or (co > 0) == (sign > 0)) and all(mon[i] == 0 for i in range(len(self.CV)) if i not in idx):
                return mon, int(co)
        return None

    def pos_poly(self, poly, strict):
        return all(c > 0 for c in poly.coeffs()) and self.strict_monomial(poly, strict) is not None

    def substitute(self, G, names, sub, strict):
        """G(sub) * prod_v D_v^{deg_v G}; every D_v is checked to be a positive polynomial"""
        DEG = dict(zip(names, G.degrees()))
        for v, rf in sub.items():
            if not self.pos_poly(rf.d, strict):
                raise AssertionError(f"denominator of {v} not certified positive")
        pw = {v: ([sub[v].n ** k for k in range(DEG[v] + 1)], [sub[v].d ** k for k in range(DEG[v] + 1)]) for v in sub}
        tot = self.ctx.from_dict({})
        for mon, co in zip(G.monoms(), G.coeffs()):
            term = int(co) * self.one
            for i, v in enumerate(names):
                if v not in sub:
                    assert mon[i] == 0
                    continue
                term = term * pw[v][0][mon[i]] * pw[v][1][DEG[v] - mon[i]]
            tot += term
        return tot


class RF:
    def __init__(s, ch, n, dd=None):
        s.ch = ch
        if dd is None:
            dd = ch.one
        g = n.gcd(dd)
        if not g.is_one():
            n = n // g; dd = dd // g
        if dd.leading_coefficient() < 0:
            n = -n; dd = -dd
        s.n, s.d = n, dd
    def _c(s, o): return o if isinstance(o, RF) else s.ch.const(o)
    def __add__(s, o): o = s._c(o); return RF(s.ch, s.n * o.d + o.n * s.d, s.d * o.d)
    __radd__ = __add__
    def __neg__(s): return RF(s.ch, -s.n, s.d)
    def __sub__(s, o): return s + (-s._c(o))
    def __rsub__(s, o): return s._c(o) - s
    def __mul__(s, o): o = s._c(o); return RF(s.ch, s.n * o.n, s.d * o.d)
    __rmul__ = __mul__
    def __truediv__(s, o): o = s._c(o); return RF(s.ch, s.n * o.d, s.d * o.n)
    def __rtruediv__(s, o): return s._c(o) / s
    def __pow__(s, k): return RF(s.ch, s.n ** k, s.d ** k)


def sign_report(name, poly, ch, strict):
    cs = [int(c) for c in poly.coeffs()]
    neg = sum(1 for c in cs if c < 0); pos = len(cs) - neg
    sgn = -1 if pos == 0 else (1 if neg == 0 else 0)
    sm = ch.strict_monomial(poly, strict) if sgn else None
    check(name + " : single strict sign with a strict monomial", sgn != 0 and sm is not None and len(cs) > 0,
          f"terms {len(cs)}, negative {neg}, positive {pos}, strict monomial (exponents {ch.CV}) {sm}")
    return len(cs), neg, pos


def numeric_minphi(vals, rng, restarts=20):
    import numpy as np
    from scipy.optimize import minimize
    f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
    return min(minimize(f, rng.uniform(-np.pi, np.pi, 6), method="BFGS").fun for _ in range(restarts))
# ---- end common block ----

# ---- root_moc4_T46_whole_assembly_20261023 (paper T_47) ----
import numpy as np
NAMES = ["q", "p", "j", "c", "e", "d", "b", "a"]
GENS = [27, 1, 59, 5, 2, 28, 6, 4]
T = [1, 3, 5, 6, 7, 8]
FV = L.F_VALUES
check("contact set equals paper T_47 literally", sorted(L.CONTACT_SETS[47]) == T, sorted(L.CONTACT_SETS[47]))
check("(1) generators(T) = {f0,f1,f2,f3,f4,f5,f6,f27,f28,f59}; only f28 non-minor",
      L.generators(frozenset(T)) == [0, 1, 2, 3, 4, 5, 6, 27, 28, 59] and [j for j in L.generators(frozenset(T)) if not L.F_IS_MINOR[j]] == [28])
check("(2a) f28 + f0 = f1 + f2 + f5 + f27", all(FV[28][s] + FV[0][s] == FV[1][s] + FV[2][s] + FV[5][s] + FV[27][s] for s in range(24)))
check("(2b) f28 + f3 = f4 + f5 + f27", all(FV[28][s] + FV[3][s] == FV[4][s] + FV[5][s] + FV[27][s] for s in range(24)))
# absorption bookkeeping, symbolically: x0 f0 + x3 f3 + sum y_k f_k with u <= x0, v <= x3, u + v <= d
sy = sympy.symbols("x0 x3 q p j c e d b a u v")
x0, x3, q_, p_, j_, c_, e_, d_, b_, a_, u, v = sy
lhs = [x0 * FV[0][s] + x3 * FV[3][s] + sum(cf * FV[g][s] for cf, g in zip((q_, p_, j_, c_, e_, d_, b_, a_), GENS)) for s in range(24)]
new = {27: q_ + u + v, 1: p_ + u, 59: j_, 5: c_ + u + v, 2: e_ + u, 28: d_ - u - v, 6: b_, 4: a_ + v}
rhs = [(x0 - u) * FV[0][s] + (x3 - v) * FV[3][s] + sum(new[g] * FV[g][s] for g in GENS) for s in range(24)]
check("absorption: consuming u of f0 and v of f3 raises q,c by u+v, p,e by u, a by v and lowers d by u+v (identity on all 24 values)",
      all(sympy.expand(l - r) == 0 for l, r in zip(lhs, rhs)))
check("canonical form (3) has common zero set T, and its 8 rays are linearly independent (so they span the 8-dim contact nullspace)", common_zero(GENS) == T and L.rank([L.F_VECS[g] for g in GENS]) == 8)
for lab, k in (("q", 29), ("p", 36), ("j", 1), ("c", 39), ("e", 42), ("b", 7)):
    drop = GENS[NAMES.index(lab)]
    Z = common_zero([g for g in GENS if g != drop])
    gg = contains_image(L.CONTACT_SETS[k], Z)
    check(f"(4) zero {lab}=0: contact set {Z} contains a G-image of paper T_{k} (with the a resource present)", gg is not None, f"g={gg}")
Z = common_zero([g for g in GENS if g != 4])
check("a=0 face is exactly T u {14} = Table 5 set T^(7)_26", Z == sorted(T + [14]) and sorted(L.AUX7[26]) == Z, Z)
ctx, m, C0, Q, H = face_setup(NAMES, GENS)
vv = dict(zip(ctx.names(), ctx.gens()))
q, p, j, c, e, d, b, a = (vv[x] for x in NAMES)
M_ = (c + 1) * (a + 1 + e); L_ = M_ + b * (a + 1)
check("C0 (d=1) = -(q+1)L - jM, M=(c+1)(a+1+e), L=M+b(a+1); Q (d=1) = jp - j(c+e+1) - (c+1)(e+1)",
      C0.subs({"d": 1}) == -(q + 1) * L_ - j * M_ and Q.subs({"d": 1}) == j * p - j * (c + e + 1) - (c + 1) * (e + 1), f"C0 = {C0}; Q = {Q}")
check("C0 (d=1) has only negative coefficients and constant -1", all(x < 0 for x in C0.subs({"d": 1}).coeffs()) and C0.subs({"d": 1}).to_dict().get((0,) * 9, 0) == -1)
t0 = time.time()
D = bareiss_det(ctx, H)
hd, N = radial_numerator(ctx, D, C0, Q)
fac = N.factor()
G9 = max(fac[1], key=lambda t_: len(t_[0]))[0]
check("radial identity: Q^4 det Hess C(g - (C0/Q) kappa) = -4 Q^2 G9^2, i.e. det = -4 G9^2/Q^2",
      N * Q ** (4 - hd) == -4 * Q ** 2 * G9 ** 2, f"h-degree {hd}; factors {[(len(f), e_) for f, e_ in fac[1]]}; content {fac[0]}; {time.time()-t0:.0f}s")
check("G9 is irreducible over Q with 1556 terms", len(G9) == 1556 and all(len(f) in (len(Q), 1556) for f, _ in fac[1]), f"{len(G9)} terms, degrees {G9.degrees()}")
# (6) partition of positive (q, b)
check("(6) sectors {q,b<=1 or q=b}, {qb<=1, max>1}, {qb>1, q!=b} cover (0,oo)^2", True, "PROOF-READ: if qb<=1 then max(q,b)<=1 or >1; if qb>1 then q=b (sector 1) or q!=b (sector 3)")
# NUMERIC: necessity sanity (G9 vanishes on matching vectors whose gaps vanish on T) and a brief main-statement test
from scipy.optimize import least_squares
mons = np.array([mm[:-1] for mm in G9.monoms()], float); cos = np.array([float(x) for x in G9.coeffs()])
PV = np.array([[float(x) for x in L.values(L.F_VECS[g])] for g in GENS]).T
rng = np.random.default_rng(7); worst, ctrl, cnt = 0.0, [], 0
for it in range(40):
    def res(x):
        aa = x[0:4] + 1j * x[4:8]; bb = x[8:12] + 1j * x[12:16]; al = np.exp(1j * x[16])
        F = L.F_values_np(aa, bb, al).real
        return np.array([F[s] - F[T[0]] for s in T[1:]])
    sol = least_squares(res, rng.normal(size=17), xtol=1e-15, ftol=1e-15, gtol=1e-15)
    if np.max(np.abs(res(sol.x))) > 1e-11:
        continue
    x = sol.x; aa = x[0:4] + 1j * x[4:8]; bb = x[8:12] + 1j * x[12:16]; al = np.exp(1j * x[16])
    F = L.F_values_np(aa, bb, al).real; gap = F - F[T[0]]
    coef = np.linalg.lstsq(PV, gap, rcond=None)[0]
    if np.max(np.abs(PV @ coef - gap)) > 1e-8 * np.max(np.abs(gap)):
        continue
    tv = cos * np.prod(np.power(coef[None, :], mons), axis=1)
    worst = max(worst, abs(tv.sum()) / np.abs(tv).sum()); cnt += 1
    tv = cos * np.prod(np.power(rng.uniform(0.1, 2, 8)[None, :], mons), axis=1)
    ctrl.append(abs(tv.sum()) / np.abs(tv).sum())
    if cnt >= 10:
        break
print(f"NUMERIC necessity: G9 at {cnt} matching vectors with gaps vanishing on T (seed 7; 8-ray span = contact nullspace): max |G9|/sum|terms| = {worst:.2e}; control at random coefficients: min {min(ctrl):.2e}")
worst = np.inf; n = 0
for seed in range(4):
    pt = L.random_point_on_contact_stratum(L.CONTACT_SETS[47], seed=seed, tries=30)
    if pt is None:
        continue
    vmin = numeric_minphi(np.asarray(pt["gaps"]), np.random.default_rng(seed), 15)
    worst = min(worst, vmin / np.max(np.abs(pt["gaps"]))); n += 1
print(f"NUMERIC main statement: {n} points on the T_47 stratum (lib, seeds 0-3), 15 BFGS restarts each: worst min Phi / max gap = {worst:.3e}")
print("SUMMARY:", "all exact checks passed" if not FAILS else f"FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
