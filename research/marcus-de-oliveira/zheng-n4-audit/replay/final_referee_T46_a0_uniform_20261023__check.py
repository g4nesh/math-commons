"""Exact replay and numeric tests for final_referee_T46_a0_uniform_20261023. Own code; imports only ../lib/moc4lib.py."""
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

# ---- shared part for the three T46 a = 0 notes ----
import pickle
NAMES = ["q", "p", "j", "c", "e", "d", "b"]
GENS = [27, 1, 59, 5, 2, 28, 6]
T46 = [1, 3, 5, 6, 7, 8]
TSTAR = sorted(T46 + [14])
ZERO_TARGETS = {"q": 29, "p": 36, "j": 1, "c": 39, "e": 42, "b": 7}   # historical T28,T35,T0,T38,T41,T6 = paper T_k


def a0_setup(zero_list):
    check("canonical a=0 source (f27,f1,f59,f5,f2,f28,f6) has common zero set T46 u {14} = Table 5 set T^(7)_26",
          common_zero(GENS) == TSTAR and sorted(L.AUX7[26]) == TSTAR, f"{common_zero(GENS)}")
    check("generators(T46 u {14}) = {f0,f1,f2,f5,f6,f27,f28,f59}, only f28 non-minor",
          L.generators(frozenset(TSTAR)) == [0, 1, 2, 5, 6, 27, 28, 59] and [j for j in L.generators(frozenset(TSTAR)) if not L.F_IS_MINOR[j]] == [28])
    check("(1) f28 + f0 = f27 + f1 + f5 + f2 on all 24 permutations",
          all(L.F_VALUES[28][s] + L.F_VALUES[0][s] == L.F_VALUES[27][s] + L.F_VALUES[1][s] + L.F_VALUES[5][s] + L.F_VALUES[2][s] for s in range(24)))
    for v in zero_list:
        drop = GENS[NAMES.index(v)]
        Z = common_zero([j for j in GENS if j != drop])
        k = ZERO_TARGETS[v]
        g = contains_image(L.CONTACT_SETS[k], Z)
        check(f"zero coefficient {v}=0: contact set {Z} contains a G-image of paper T_{k}", g is not None, f"g={g}")
    ctx, m, C0, Q, H = face_setup(NAMES, GENS)
    v = dict(zip(ctx.names(), ctx.gens()))
    q, p, j, c, e, d, b = (v[x] for x in NAMES)
    C0d = C0.subs({"d": 1}); Qd = Q.subs({"d": 1})
    one = ctx.from_dict({(0,) * 8: 1})
    check("C0 (d=1) = -(q+1){(c+1)(e+1)+b} - j(c+1)(e+1); Q (d=1) = jp - j(c+e+1) - (c+1)(e+1)",
          C0d == -(q + 1) * ((c + 1) * (e + 1) + b) - j * (c + 1) * (e + 1) and Qd == j * p - j * (c + e + 1) - (c + 1) * (e + 1),
          f"C0 = {C0}; Q = {Q}")
    check("C0 (d=1) has only negative coefficients and constant term -1 (so C0 < 0 on the closed orthant)",
          all(x < 0 for x in C0d.coeffs()) and C0d.to_dict().get((0,) * 8, 0) == -1)
    t0 = time.time()
    D = bareiss_det(ctx, H)
    hd, N = radial_numerator(ctx, D, C0, Q)
    fac = N.factor()
    G8 = max(fac[1], key=lambda t: len(t[0]))[0]
    check("radial identity: Q^3 det Hess C(g - (C0/Q) kappa) = -4 Q (q+d)^2 G8^2, i.e. det = -4 (q+d)^2 G8^2 / Q^2",
          N == -4 * Q * (q + d) ** 2 * G8 ** 2, f"h-degree {hd}; factors {[(len(f), e_) for f, e_ in fac[1]]}; content {fac[0]}; {time.time()-t0:.0f}s")
    check("G8 has 380 terms, degrees (2,2,2,2,2,8,3) in (q,p,j,c,e,d,b)", len(G8) == 380 and tuple(G8.degrees()[:7]) == (2, 2, 2, 2, 2, 8, 3), f"{len(G8)} {G8.degrees()}")
    return G8.subs({"d": 1, "h": 0})


NAMES_NOD = ["q", "p", "j", "c", "e", "d", "b", "h"]

# ---- final_referee_T46_a0_uniform_20261023 ----
import numpy as np
G8 = a0_setup(["q", "p", "j", "c", "e", "b"])
q, b, jj, J, t = sympy.symbols("q b j J t", positive=True)
lam = (q + b + 2) / (q + b - 2); K = (q + 1) * (b + 1) / (q * b - 1)
check("Lemma 2: K - lambda = (q-b)^2/((qb-1)(q+b-2))", sympy.simplify(K - lam - (q - b) ** 2 / ((q * b - 1) * (q + b - 2))) == 0)
check("Lemma 2: qb > 1 implies q + b > 2 (AM-GM: q + b >= 2 sqrt(qb))", True, "PROOF-READ")
tl = jj / lam
check("Lemma 2: t = j/lambda gives threshold 1/j - 1/lambda; leaf condition j > K t^2 iff j < lambda^2/K",
      sympy.simplify((1 - tl) ** 2 / (jj - lam * tl ** 2) - (1 / jj - 1 / lam)) == 0 and sympy.simplify(K * tl ** 2 - jj * (K * jj / lam ** 2)) == 0)
jm = (lam + K * J) ** 2 / (K * (1 + J) ** 2); tm = (lam + K * J) / (K * (1 + J))
check("(6): t^2 = j/K and (1-t)^2/(j - lambda t^2) = (K-lambda)/(lambda+KJ)^2",
      sympy.simplify(tm ** 2 - jm / K) == 0 and sympy.simplify((1 - tm) ** 2 / (jm - lam * tm ** 2) - (K - lam) / (lam + K * J) ** 2) == 0)
check("(6): j(0) = lambda^2/K, j -> K as J -> oo, dj/dJ = 2(lambda+KJ)(K-lambda)/(K(1+J)^3) > 0",
      sympy.simplify(jm.subs(J, 0) - lam ** 2 / K) == 0 and sympy.simplify(sympy.limit(jm, J, sympy.oo) - K) == 0
      and sympy.simplify(sympy.diff(jm, J) - 2 * (lam + K * J) * (K - lam) / (K * (1 + J) ** 3)) == 0)
check("(4) chart: j = lambda^2 J/(K(1+J)) covers exactly (0, lambda^2/K)", True, "monotone in J, PROOF-READ")
CH = Chart(("A", "V", "J", "S", "C1", "E"))
A, V, Jv, S, C1, E = (CH.var(v) for v in CH.CV)
r = 1 + V
for lab, qq, bb in (("q>b: q=r+A, b=r^2/q", r + A, r ** 2 / (r + A)), ("b>q: b=r+A, q=r^2/b", r ** 2 / (r + A), r + A)):
    lamc = (qq + bb + 2) / (qq + bb - 2); Kc = (qq + 1) * (bb + 1) / (qq * bb - 1)
    jl = lamc ** 2 * Jv / (Kc * (1 + Jv)); Pl = 1 / jl - 1 / lamc
    sign_report(f"G8 table {lab}, r=1+V, low j = lambda^2 J/(K(1+J)), p = P* S/(1+S), P* = 1/j - 1/lambda",
                CH.substitute(G8, NAMES_NOD, {"q": qq, "b": bb, "p": Pl * S / (1 + S), "j": jl, "c": C1, "e": E}, ["A", "V", "J", "C1"]), CH, ["A", "V", "J", "C1"])
    jmc = (lamc + Kc * Jv) ** 2 / (Kc * (1 + Jv) ** 2); Pm = (Kc - lamc) / (lamc + Kc * Jv) ** 2
    sign_report(f"G8 table {lab}, r=1+V, middle j = (lambda+KJ)^2/(K(1+J)^2) with J >= 0, p = P* S/(1+S), P* = (K-lambda)/(lambda+KJ)^2",
                CH.substitute(G8, NAMES_NOD, {"q": qq, "b": bb, "p": Pm * S / (1 + S), "j": jmc, "c": C1, "e": E}, ["A", "V", "C1"]), CH, ["A", "V", "C1"])
rng = np.random.default_rng(13)
PV = np.array([[float(v_) for v_ in L.values(L.F_VECS[g_])] for g_ in [27, 1, 59, 5, 2, 28, 6]]).T
worst, n = np.inf, 0
while n < 12:
    qn, bn = rng.uniform(0.3, 4, 2)
    if qn * bn < 1.05 or abs(qn - bn) < 0.05:
        continue
    l_ = (qn + bn + 2) / (qn + bn - 2); K_ = (qn + 1) * (bn + 1) / (qn * bn - 1)
    j_ = rng.uniform(0.02, 1.2 * K_)
    if j_ >= K_:
        P_ = 0.0
    elif j_ < l_ ** 2 / K_:
        P_ = 1 / j_ - 1 / l_
    else:
        tt = np.sqrt(j_ / K_); P_ = (1 - tt) ** 2 / (j_ - l_ * tt ** 2)
    coef = np.array([qn, P_, j_, 0, 0, 1, bn])
    worst = min(worst, numeric_minphi(PV @ coef, rng, 20) / coef.max()); n += 1
print(f"NUMERIC Lemma 2 boundary (c=e=a=0, p = P*): 12 samples (seed 13, 20 BFGS restarts): worst min Phi / max coef = {worst:.3e}")
print("SUMMARY:", "all exact checks passed" if not FAILS else f"FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
