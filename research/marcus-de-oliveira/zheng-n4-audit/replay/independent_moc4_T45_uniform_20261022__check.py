"""Exact replay and numeric tests for independent_moc4_T45_uniform_20261022 (paper T_46). Own code; imports only ../lib/moc4lib.py."""
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

# ---- independent_moc4_T45_uniform_20261022 (paper T_46) ----
# usage: python -I -B <this> [all|setup|direct LO HI|square|numeric]
# optional env MOC4_CACHE: directory where the computed G9 is cached between partial runs.
import pickle
stage = sys.argv[1] if len(sys.argv) > 1 else "all"
NAMES = ["p", "r", "a", "b", "c", "q", "d", "s"]
GENS = [0, 1, 2, 12, 59, 67, 70, 72]
T = [3, 4, 5, 9, 10, 12]


def get_G9():
    cache = os.environ.get("MOC4_CACHE")
    fn = os.path.join(cache, "G9_T45face.pkl") if cache else None
    if fn and os.path.exists(fn):
        dct = pickle.load(open(fn, "rb"))
        ctx = flint.fmpz_mpoly_ctx.get(tuple(NAMES) + ("h",), "lex")
        print("INFO G9 loaded from cache (it was computed by this script in an earlier run)")
        return ctx.from_dict(dct["G9"]), ctx.from_dict(dct["C0"]), ctx.from_dict(dct["Q"]), ctx
    t0 = time.time()
    ctx, m, C0, Q, H = face_setup(NAMES, GENS)
    D = bareiss_det(ctx, H)
    hd, N = radial_numerator(ctx, D, C0, Q)
    fac = N.factor()
    d_ = dict(zip(ctx.names(), ctx.gens()))
    Qdoc = (d_["b"] - d_["c"]) * d_["a"] + d_["c"] * (d_["r"] + d_["s"]) - (d_["b"] + d_["d"]) * (d_["d"] + d_["r"] + d_["s"])
    G9 = max(fac[1], key=lambda t: len(t[0]))[0]
    # (9): Q^4 det = -4 Q^2 G9^2, i.e. Q^hd det = -4 Q^(hd-2) G9^2
    ok = (N * Q ** (4 - hd) == -4 * Q ** 2 * G9 ** 2)
    check("(9) Q^4 det Hess C(g - (C/Q) kappa) = -4 Q^2 G9^2 (full 14x14 Hessian of the paper's cubic (A.2), Bareiss determinant)", ok,
          f"h-degree of det {hd}; factorization content {fac[0]}, factors (terms, exponent) {[(len(f), e) for f, e in fac[1]]}; {time.time()-t0:.0f}s")
    check("G9 has 1158 terms, total degree 9", len(G9) == 1158 and G9.total_degree() == 9, f"terms {len(G9)}, degrees {G9.degrees()}")
    if fn:
        pickle.dump({"G9": G9.to_dict(), "C0": C0.to_dict(), "Q": Q.to_dict()}, open(fn, "wb"))
    return G9, C0, Q, ctx


if stage in ("all", "setup"):
    check("contact set equals paper T_46 literally", sorted(L.CONTACT_SETS[46]) == T, sorted(L.CONTACT_SETS[46]))
    check("generators(T) = {f0,f1,f2,f12,f59,f67,f70,f72}", L.generators(frozenset(T)) == GENS, L.generators(frozenset(T)))
    check("only f70 is non-minor", [j for j in GENS if not L.F_IS_MINOR[j]] == [70])
    ev = [[L.F_VALUES[j][s] for s in (18, 16, 13, 8, 2, 14, 11, 7)] for j in GENS]
    check("(1) evaluations at (18,16,13,8,2,14,11,7) form the identity matrix", ev == [[int(i == k) for k in range(8)] for i in range(8)])
    check("nullspace of the six contact equations has dimension 8", L.rank([L.P[s] for s in T]) == 6)
    # zero faces: a=0 (f2), c=0 (f59), r=0 (f1)
    for lab, drop, newc, k, img in (("a", 2, 13, 34, [3, 5, 9, 10, 12, 13]), ("c", 59, 2, 7, [2, 3, 4, 5, 9, 10]), ("r", 1, 16, 15, [3, 4, 5, 9, 10, 16])):
        Z = common_zero([j for j in GENS if j != drop])
        g = contains_image(L.CONTACT_SETS[k], Z)
        check(f"zero face {lab}=0: contact set {Z} = T u {{{newc}}} and contains a G-image of paper T_{k} (historical label as printed)", Z == sorted(T + [newc]) and g is not None, f"g={g}")
        check(f"zero face {lab}=0: printed image {img} is a G-image of paper T_{k}", L.find_g(L.CONTACT_SETS[k], img) is not None)
    G9, C0, Q, ctx = get_G9()
    v = dict(zip(ctx.names(), ctx.gens()))
    p, r, a, b, c, q, d, s = (v[x] for x in NAMES)
    Cdoc = c * r * s - (d + q + s) * ((b + d) * (d + r) + d * p)
    Qdoc = (b - c) * a + c * (r + s) - (b + d) * (d + r + s)
    check("(7) C(g + h kappa) = C + h Q with the printed C and Q (paper cubic (A.2))", C0 == Cdoc and Q == Qdoc, f"C0 = {C0}; Q = {Q}")
    R = r + s; T0 = (b + d) * (d + r) + d * p
    Pp = d ** 3 * R + d ** 2 * (r ** 2 + r * s + s ** 2) + R * (d + s) * (b * (d + r) + d * p) + R * q * T0
    check("(8) R C = (cR - dR - d^2) r s - P exactly", R * Cdoc == (c * R - d * R - d ** 2) * r * s - Pp)
    check("(8) P has positive coefficients and contains the monomial d^3 r (so P > 0 when d, r > 0)",
          all(x > 0 for x in Pp.coeffs()) and Pp.to_dict().get((0, 1, 0, 0, 0, 0, 3, 0, 0), 0) > 0)

CH = Chart(("T", "U", "V", "A", "M", "P", "Q0"))


def chart_sub(sx, sy, sz, use14, pmode):
    Tv, Uv, Vv, Av, Mv, Pv, Qv = (CH.var(x) for x in CH.CV)
    h = Tv / (1 + Tv)
    x = {-1: Uv / (1 + Uv), 0: CH.const(1), 1: 1 + Uv}[sx]
    y = {-1: Vv / (1 + Vv), 0: CH.const(1), 1: 1 + Vv}[sy]
    z = {-1: Av / (1 + Av), 0: CH.const(1), 1: 1 + Av}[sz]
    c = (y + z) * Mv / ((y + z - 1) * (1 + Mv)) if use14 else Mv
    a = h / c; b = x * c / h; r = y / c; s = z / c
    p = {"PB": (b / (x - 1)) * Pv / (1 + Pv) if sx == 1 else None, "PR": (r / (y - 1)) * Pv / (1 + Pv) if sy == 1 else None, "free": Pv}[pmode]
    q = (s / (z - 1)) * Qv / (1 + Qv) if sz == 1 else Qv
    return {"p": p, "r": r, "a": a, "b": b, "c": c, "q": q, "s": s}


def strict_of(sx, sy, sz):
    st = ["T", "M", "V"] + (["U"] if sx == 1 else []) + (["A"] if sz == 1 else [])
    return st


if stage in ("all", "direct"):
    G9, C0, Q, ctx = get_G9()
    G9d = G9.subs({"d": 1, "h": 0})
    triples = [t for t in itertools.product((-1, 0, 1), repeat=3) if t != (1, 1, -1)]
    lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if stage == "direct" else (0, 26)
    for (sx, sy, sz) in triples[lo:hi]:
        t0 = time.time()
        use14 = sy >= 0 and sz >= 0
        pmode = "PB" if sx == 1 else ("PR" if sy == 1 else "free")
        st = strict_of(sx, sy, sz)
        poly = CH.substitute(G9d, NAMES, chart_sub(sx, sy, sz, use14, pmode), st)
        cs = [int(x) for x in poly.coeffs()]
        neg = sum(1 for x in cs if x < 0)
        sm = CH.strict_monomial(poly, st, sign=1)
        check(f"Lemma 1 domain signs(x-1,y-1,z-1)={(sx,sy,sz)} p-bound {pmode} (14) {'used' if use14 else 'not used'}: nonnegative coefficients and a positive monomial in strictly positive variables",
              neg == 0 and sm is not None, f"terms {len(cs)}, negative {neg}, strict monomial (T,U,V,A,M,P,Q0) {sm}, {time.time()-t0:.1f}s")

if stage in ("all", "square"):
    G9, C0, Q, ctx = get_G9()
    G9d = G9.subs({"d": 1, "h": 0})
    st = ["T", "U", "V", "M"]
    Tv, Uv, Vv = CH.var("T"), CH.var("U"), CH.var("V")
    Lrf = Tv * Uv * (1 + Vv) / ((1 + Tv) * (1 + Uv) * Vv)
    sub_pb, sub_pr = chart_sub(1, 1, -1, False, "PB"), chart_sub(1, 1, -1, False, "PR")
    # (15): b/U <= r/V  <=>  M^2 <= L
    lhs = sub_pb["b"] / Uv - sub_pb["r"] / Vv
    Mv = CH.var("M")
    ratio = lhs / (Mv ** 2 - Lrf)
    check("(15) b/U - r/V = (M^2 - L) * (positive rational function)", CH.pos_poly(ratio.n, st) and CH.pos_poly(ratio.d, st), f"ratio = ({ratio.n})/({ratio.d})")
    iM = CH.CV.index("M")
    for side, sub in (("PB", sub_pb), ("PR", sub_pr)):
        poly = CH.substitute(G9d, NAMES, sub, st)
        mins = min(mm[iM] for mm in poly.monoms()); degM = max(mm[iM] for mm in poly.monoms())
        coefM = {}
        for mon, co in zip(poly.monoms(), poly.coeffs()):
            mm = list(mon); k = mm[iM] - mins; mm[iM] = 0
            coefM.setdefault(k, {})[tuple(mm)] = int(co)
        n = degM - mins
        cM = [CH.ctx.from_dict(coefM.get(k, {})) for k in range(n + 1)]
        if side == "PR":
            cM = cM[::-1]; Ln, Ld = Lrf.d, Lrf.n   # w = 1/M, w^2 in [0, 1/L]
        else:
            Ln, Ld = Lrf.n, Lrf.d                  # M^2 in [0, L]
        for part, co in (("even", cM[0::2]), ("odd", cM[1::2])):
            mdeg = len(co) - 1
            bern = []
            for k in range(mdeg + 1):
                ck = CH.ctx.from_dict({})
                for i in range(k + 1):
                    if not co[i].is_zero():
                        ck += co[i] * Ln ** i * Ld ** (mdeg - i) * comb(mdeg - i, k - i)
                bern.append(ck)
            negs = [sum(1 for x in bb.coeffs() if x < 0) for bb in bern]
            sm = CH.strict_monomial(bern[-1], ["T", "U", "V"], sign=1)
            check(f"3.1 square interval, {side} side (polynomial degree {n} in {'M' if side=='PB' else 'w=1/M'} after removing M^{mins}), {part} part of degree {mdeg} in the squared variable: all Bernstein-type coefficients nonnegative, last one has a positive monomial in T,U,V",
                  all(x == 0 for x in negs) and sm is not None, f"term counts {[len(bb) for bb in bern]}, negative counts {negs}, strict monomial {sm}")

if stage in ("all", "numeric"):
    import numpy as np
    from scipy.optimize import least_squares
    G9, C0, Q, ctx = get_G9()
    mons = np.array([m_[:-1] for m_ in G9.monoms()], float); cos = np.array([float(x) for x in G9.coeffs()])
    PV = np.array([[float(x) for x in L.values(L.F_VECS[j])] for j in GENS]).T
    rng = np.random.default_rng(7)
    worst, ctrl, cnt, skipped = 0.0, [], 0, 0
    check("G9 vanishes identically on r = d = 0 (so the numeric necessity test skips points near it)", G9.subs({'r': 0, 'd': 0}).is_zero())
    for it in range(40):
        def res(x):
            a_ = x[0:4] + 1j * x[4:8]; b_ = x[8:12] + 1j * x[12:16]; al = np.exp(1j * x[16])
            F = L.F_values_np(a_, b_, al).real
            return np.array([F[s_] - F[T[0]] for s_ in T[1:]])
        sol = least_squares(res, rng.normal(size=17), xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.max(np.abs(res(sol.x))) > 1e-11:
            continue
        x = sol.x; a_ = x[0:4] + 1j * x[4:8]; b_ = x[8:12] + 1j * x[12:16]; al = np.exp(1j * x[16])
        F = L.F_values_np(a_, b_, al).real; gap = F - F[T[0]]
        coef = np.linalg.lstsq(PV, gap, rcond=None)[0]
        if min(abs(coef[1]), abs(coef[6])) < 1e-3 * np.max(np.abs(coef)):
            skipped += 1   # G9 vanishes identically on {r = d = 0}; the relative test is meaningless there
            continue
        tv = cos * np.prod(np.power(coef[None, :], mons), axis=1)
        worst = max(worst, abs(tv.sum()) / np.abs(tv).sum()); cnt += 1
        tv = cos * np.prod(np.power(rng.uniform(0.1, 2, 8)[None, :], mons), axis=1)
        ctrl.append(abs(tv.sum()) / np.abs(tv).sum())
        if cnt >= 10:
            break
    print(f"NUMERIC necessity: G9 at {cnt} matching vectors whose gaps vanish on T (seed 7, not sign-constrained): max |G9|/sum|terms| = {worst:.2e}; control at random coefficients: min {min(ctrl):.2e}; skipped {skipped} points with r or d below 1e-3 relative (found with a_i + b_j = 0)")
    worst = np.inf; n = 0
    for seed in range(4):
        pt = L.random_point_on_contact_stratum(L.CONTACT_SETS[46], seed=seed, tries=30)
        if pt is None:
            continue
        v = numeric_minphi(np.asarray(pt["gaps"]), np.random.default_rng(seed), 15)
        worst = min(worst, v / np.max(np.abs(pt["gaps"]))); n += 1
    print(f"NUMERIC main statement: {n} points on the T_46 stratum (lib, seeds 0-3), 15 BFGS restarts each: worst min Phi / max gap = {worst:.3e}")

print("SUMMARY:", "all exact checks passed" if not FAILS else f"FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
