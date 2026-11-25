# ---------------------------------------------------------------------------------
# Common block (identical in the four K1-cores-BCD replay scripts).
# Own code only. The supplement is read as DATA (json.load) and never executed.
# Usage: python -I -B <this script> [--skip-heavy]
#   env MOC4_SUPP = root of the supplement copy (needed only for the JSON data of G).
# Prints one line per check; exits 1 if an EXACT check fails.
# ---------------------------------------------------------------------------------
import sys, os, json, time, math, random, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint

FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""), flush=True)
    if not ok:
        FAILS.append(name)

def info(name, msg):
    print("INFO " + name + " : " + str(msg), flush=True)

SUPP = os.environ.get("MOC4_SUPP")
SKIP_HEAVY = "--skip-heavy" in sys.argv
def supp_json(name):
    if not SUPP:
        return None
    p = os.path.join(SUPP, "proofs", "manin_transfer_20261021", name)
    return json.load(open(p)) if os.path.exists(p) else None

def make_ctx(names):
    return flint.fmpq_mpoly_ctx.get(tuple(names), 'lex')

def zero(ctx):
    return ctx.from_dict({})

def const(ctx, v):
    if isinstance(v, Fr):
        v = flint.fmpq(v.numerator, v.denominator)
    return ctx.from_dict({(0,) * ctx.nvars(): v})

def jpoly(ctx, lst, keep=None):
    """JSON term list [[exps, "coef"], ...] -> fmpq_mpoly; keep = indices of exponent slots
    mapped to ctx variables (others must be 0 or are evaluated at 1 when listed in ones)."""
    d = {}
    for e, c in lst:
        e2 = tuple(e[i] for i in keep) if keep is not None else tuple(e)
        d[e2] = d.get(e2, 0) + Fr(c)
    return ctx.from_dict({k: flint.fmpq(v.numerator, v.denominator) for k, v in d.items() if v != 0})

def lift(ctx, coefs, src, z):
    """w = sum_k coefs[k] f_{src[k]} + z kappa  (14 coordinates)."""
    w = [zero(ctx) for _ in range(14)]
    for k, j in enumerate(src):
        for i in range(14):
            v = L.F_VECS[j][i]
            if v != 0:
                w[i] = w[i] + coefs[k] * flint.fmpq(v.numerator, v.denominator)
    w[13] = w[13] + z
    return w

def split_z(ctx, P, zi):
    """coefficients of P as a polynomial in variable index zi."""
    d = P.to_dict()
    deg = max((k[zi] for k in d), default=0)
    out = []
    for e in range(deg + 1):
        out.append(ctx.from_dict({k[:zi] + (0,) + k[zi + 1:]: v for k, v in d.items() if k[zi] == e}))
    return out

_MCTX = make_ctx(["m%d" % i for i in range(14)])
_CM = L.cubic_C(_MCTX.gens())

def hessian_at(ctx, w):
    H = [[zero(ctx) for _ in range(14)] for _ in range(14)]
    for i in range(14):
        di = _CM.derivative(i)
        for j in range(14):
            dij = di.derivative(j)
            if not dij.is_zero():
                H[i][j] = dij.compose(*w, ctx=ctx)
    return H

def det_dp(ctx, H, rows=None, cols=None):
    """Exact determinant by row-by-row Laplace expansion over column subsets."""
    rows = list(range(len(H))) if rows is None else rows
    cols = list(range(len(H))) if cols is None else cols
    n = len(rows)
    dp = {0: const(ctx, 1)}
    for r in rows:
        nd = {}
        for mask, val in dp.items():
            for ci, c in enumerate(cols):
                if mask >> ci & 1 or H[r][c].is_zero():
                    continue
                term = val * H[r][c]
                if bin(mask >> (ci + 1)).count("1") % 2:
                    term = -term
                nm = mask | (1 << ci)
                nd[nm] = nd[nm] + term if nm in nd else term
        dp = {k: v for k, v in nd.items() if not v.is_zero()}
    return dp.get((1 << n) - 1, zero(ctx))

def subs_clear(P, ctx_out, subs):
    """P in ctx_in variables v_i; subs[i] = (N_i, D_i) fmpq_mpoly in ctx_out (or None to keep
    nothing: variable must be absent).  Returns P(N/D) * prod D_i^{deg_i P}."""
    d = P.to_dict()
    nv = len(subs)
    degs = [max((k[i] for k in d), default=0) for i in range(nv)]
    powN = [[const(ctx_out, 1)] for _ in range(nv)]
    powD = [[const(ctx_out, 1)] for _ in range(nv)]
    for i in range(nv):
        if degs[i] == 0:
            continue
        N, D = subs[i]
        for _ in range(degs[i]):
            powN[i].append(powN[i][-1] * N)
            powD[i].append(powD[i][-1] * D)
    out = zero(ctx_out)
    for k, c in d.items():
        t = const(ctx_out, c)
        for i in range(nv):
            if degs[i]:
                t = t * powN[i][k[i]] * powD[i][degs[i] - k[i]]
        out = out + t
    return out

def sign_profile(P, strict_vars):
    """(all coefficients >= 0, all <= 0, exists nonzero monomial in strict vars only)."""
    d = P.to_dict()
    pos = all(v >= 0 for v in d.values())
    neg = all(v <= 0 for v in d.values())
    strict = any(all(e == 0 or i in strict_vars for i, e in enumerate(k)) and v != 0 for k, v in d.items())
    return pos, neg, strict, len(d)

def rat(ctx, num, den=None):
    return (num, den if den is not None else const(ctx, 1))

# =================================================================================
# blind_moc4_coreB_uniform_20261023 : seven-contact set T = {0,1,2,7,9,10,16}
# =================================================================================
T = {0, 1, 2, 7, 9, 10, 16}
SRC = [0, 2, 3, 11, 12, 15, 26]          # q a t u b c d
EXCL = [23, 15, 14, 17, 8, 4, 3]
t0 = time.time()
check("T equals Table 5 seven-set AUX7[8]", frozenset(T) == L.AUX7[8])
check("generators(T) == source list", L.generators(T) == SRC, L.generators(T))
check("rank of P rows at T is 7", L.rank([L.P[s] for s in sorted(T)]) == 7)
check("exclusive rows form I7", [[L.F_VALUES[j][e] for e in EXCL] for j in SRC] == [[int(i == k) for k in range(7)] for i in range(7)])
check("first six generators minor, f26 non-minor", [L.F_IS_MINOR[j] for j in SRC] == [True] * 6 + [False])

names = ["q", "a", "t", "u", "b", "c", "z"]
ctx = make_ctx(names)
q, a, t, u, b, c, z = ctx.gens()
one = const(ctx, 1)
w = lift(ctx, [q, a, t, u, b, c, one], SRC, z)
Cw = L.cubic_C(w)
parts = split_z(ctx, Cw, 6)
check("C(w(z)) is affine in z", len(parts) == 2)
C0, Q1 = parts[0], parts[1]
check("(8) C = (1+u)(1+q+t)+abc  (d=1)", C0 == (1 + u) * (1 + q + t) + a * b * c)
check("(8) Q = (1+u)(1+t)+ab+ac+bc-ub-tc  (d=1)", Q1 == (1 + u) * (1 + t) + a * b + a * c + b * c - u * b - t * c)

J = supp_json("blind_moc4_coreB_uniform_exact_20261023.json")
if J is None:
    info("G", "MOC4_SUPP not set: G-dependent checks skipped")
else:
    G = jpoly(ctx, J["G"], keep=[0, 1, 2, 3, 4, 5, 7])     # d = 1, height slot -> z (always 0)
    check("JSON G has 397 terms (homogeneous)", len(J["G"]) == 397, len(J["G"]))
    check("G has no z", G.degrees()[6] == 0)
    if not SKIP_HEAVY:
        H = hessian_at(ctx, w)
        detH = det_dp(ctx, H)
        Dk = split_z(ctx, detH, 6)
        info("det H(w(z))", "terms %d, z-degree %d, %.1fs" % (len(detH), len(Dk) - 1, time.time() - t0))
        lhs = zero(ctx)
        for e, Dke in enumerate(Dk):
            lhs = lhs + Dke * (-C0) ** e * Q1 ** (4 - e)
        check("(9) Q^4 det H(w(-C/Q)) = -4 Q^2 G^2 (own Hessian, d=1)", lhs == -4 * Q1 ** 2 * G ** 2)
    gu = split_z(ctx, G, 3)
    check("(10) G quadratic in u", len(gu) == 3)
    gam, bet, alp = gu
    check("(10) beta = alpha + gamma + a(a+1) b^2 c^2 t^2", bet == alp + gam + a * (a + 1) * b ** 2 * c ** 2 * t ** 2)
    # c = 0 factor, homogeneous: use d as a variable
    ctxh = make_ctx(["q", "a", "t", "u", "b", "c", "d"])
    qh, ah, th, uh, bh, ch, dh = ctxh.gens()
    Gh = jpoly(ctxh, J["G"], keep=[0, 1, 2, 3, 4, 5, 6])
    Gc0 = Gh.subs({"c": 0}) if False else Gh.compose(qh, ah, th, uh, bh, zero(ctxh), dh, ctx=ctxh)
    fac = dh ** 2 * (uh + dh) * (ah + dh) * (qh + th + dh) ** 2 * (qh * (uh + dh) * (bh + dh) + uh * (bh ** 2 + bh * dh + dh ** 2 - ah * bh) + (bh + dh) * (dh ** 2 - ah * bh))
    check("c=0 factorization of G", Gc0 == fac)
    check("G homogeneous of degree 9", all(sum(k) == 9 for k in Gh.to_dict()))
    # the bracket is positive when ab < d^2: rewrite with s = d^2 - ab >= 0 ... check: bracket =
    # q(u+d)(b+d) + u(b^2+bd+s) + (b+d) s with s = d^2-ab
    info("c=0 bracket", "q(u+d)(b+d) + u(b^2+bd+(d^2-ab)) + (b+d)(d^2-ab): each term >= 0 and (b+d)(d^2-ab) > 0 when ab<d^2, d>0 (PROOF-READ)")

    # ---------------- the thirty chart certificates ----------------
    cx = make_ctx(["Q", "A", "T", "B", "C"])
    Qv, Av, Tv, Bv, Cv = cx.gens()
    o = const(cx, 1)
    mains = {
        "R0": [(Av, 1 + Av), (Bv, 1 + Bv), (Cv, 1 + Cv)],
        "R1": [(1 + Av, o), (Bv, (1 + Av) * (1 + Bv + Cv)), (Cv, (1 + Av) * (1 + Bv + Cv))],
        "R2": [(Av, (1 + Bv) * (1 + Av + Cv)), (1 + Bv, o), (Cv, (1 + Bv) * (1 + Av + Cv))],
        "R3": [(Av, (1 + Cv) * (1 + Av + Bv)), (Bv, (1 + Cv) * (1 + Av + Bv)), (1 + Cv, o)],
        "c1": [(Av, 1 + Av + Bv), (Bv, 1 + Av + Bv), (o, o)],
    }
    # strict variable indices: Q=0 A=1 T=2 B=3 C=4
    def qt_charts(mname, cN, cD):
        # c = cN/cD
        if mname in ("R0", "R1", "R2"):
            return [("cq<1", (Qv * cD, cN * (1 + Qv)), (Tv, o), {4}),
                    ("cq>1", ((1 + Qv) * cD, cN), ((1 + Qv) * Tv * cD, Qv * cN * (1 + Tv)), {0, 4}),
                    ("cq=1", (cD, cN), (Tv, o), {4})]
        if mname == "R3":
            return [("q<1", (Qv, 1 + Qv), (Tv, Cv * (1 + Tv)), {4}),
                    ("q>1", (1 + Qv, o), ((1 + Qv) * Tv, ((1 + Cv) * (1 + Qv) - 1) * (1 + Tv)), {0, 4}),
                    ("q=1", (o, o), (Tv, Cv * (1 + Tv)), {4})]
        return [("q<1", (Qv, 1 + Qv), (Tv, o), set()),
                ("q>1", (1 + Qv, o), ((1 + Qv) * Tv, Qv * (1 + Tv)), {0}),
                ("q=1", (o, o), (Tv, o), set())]
    # alpha, gamma as polynomials in (q,a,t,b,c) : drop u (absent) and z
    ctx5 = make_ctx(["q", "a", "t", "b", "c"])
    def to5(P):
        return ctx5.from_dict({(k[0], k[1], k[2], k[4], k[5]): v for k, v in P.to_dict().items()})
    A5, G5 = to5(alp), to5(gam)
    ntab = 0
    allok = True
    for mname, (ma, mb, mc) in mains.items():
        for (cname, qs, ts, strict) in qt_charts(mname, mc[0], mc[1]):
            subs = [qs, ma, ts, mb, mc]
            for pname, P in (("alpha", A5), ("gamma", G5)):
                for (Nn, Dd) in subs:
                    dp_, dn_, ds_, _ = sign_profile(Dd, strict)
                    allok &= dp_ and ds_
                R = subs_clear(P, cx, subs)
                pos, neg, st, nt = sign_profile(R, strict)
                ntab += 1
                ok = pos and st
                allok &= ok
                check("(15) %s %s %s: coefficients >= 0 and strict monomial present" % (mname, cname, pname), ok, "terms %d" % nt)
    check("(15) thirty tables", ntab == 30)
    check("all chart denominators have nonnegative coefficients and a strict monomial", allok)
    # positivity of every denominator used (all are products of 1+..., Q, C, c-terms on strict sets)
    info("denominators", "all chart denominators are products of factors 1+X (X>=0), Q (>0 where used), C (>0 where used), (1+C)(1+Q)-1 = C+Q+CQ (>0 when C>0)")

# ---------------- chart inverse maps (exact symbolic) ----------------
import sympy as sp
A_, B_, C_, Q_, T_ = sp.symbols("A B C Q T", nonnegative=True)
a_, b_, c_, q_, t_ = sp.symbols("a b c q t", positive=True)
# R1 inverse: A=a-1, B=ab/(1-a(b+c)), C=ac/(1-a(b+c))
th = 1 - a_ * (b_ + c_)
sub = {A_: a_ - 1, B_: a_ * b_ / th, C_: a_ * c_ / th}
okR1 = all(sp.simplify(e.subs(sub) - v) == 0 for e, v in [(1 + A_, a_), (B_ / ((1 + A_) * (1 + B_ + C_)), b_), (C_ / ((1 + A_) * (1 + B_ + C_)), c_)])
check("(11) R1 inverse formula", okR1)
# chain chart cq>1 inverse: T = (cq-1)t/[q-(cq-1)t]
Tinv = (c_ * q_ - 1) * t_ / (q_ - (c_ * q_ - 1) * t_)
Qinv = c_ * q_ - 1
check("(12) cq>1 chart inverse", sp.simplify(((1 + Q_) * T_ / (Q_ * c_ * (1 + T_))).subs({Q_: Qinv, T_: Tinv}) - t_) == 0 and sp.simplify(((1 + Q_) / c_).subs(Q_, Qinv) - q_) == 0)

# ---------------- region logic, exact on a rational grid ----------------
def H3(aa, bb, cc):
    s = aa + bb + cc; M = max(aa, bb, cc)
    return s / 2 if 2 * M <= s else math.sqrt(M * (s - M))
grid = [Fr(k, 8) for k in range(0, 25)]
viol = 0
for aa in grid:
    for bb in grid:
        for cc in grid:
            if H3(float(aa), float(bb), float(cc)) >= 1:
                continue
            for (x, y) in ((aa, bb + cc), (bb, aa + cc), (cc, aa + bb)):
                if x >= 1 and x * y >= 1:
                    viol += 1
check("(7) on a 25^3 rational grid: H3<1 and a main coeff >=1 imply product with the others' sum <1", viol == 0)

# ---------------------------------------------------------------------------------
# NUMERIC part (common): capacity inequalities and realizable points on the stratum.
# ---------------------------------------------------------------------------------
def numeric_part(T, SRC, EXCL, CAP_TESTS, Gjson, gvars_slots, seed=2026, starts=40, npoints=4):
    if "--no-numeric" in sys.argv:
        info("numeric", "skipped"); return
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    Pn = L.P_numpy()
    vals = {j: Pn @ np.array([float(x) for x in L.F_VECS[j]]) for j in range(82)}
    def phis(th, js):
        wv = L.leibniz_weights_np(L.so4_from_params(th))
        return [float(vals[j] @ wv) for j in js]
    def sq(x): return math.sqrt(max(x, 0.0))
    for test in CAP_TESTS:
        kind = test[0]
        if kind == "SOC":
            _, i, j = test
            f = lambda th: (lambda A, B, N: N + 2 * sq(A) * sq(B))(*phis(th, [i, j, 26]))
        elif kind == "oldH":
            _, i, j = test        # t f_i + c f_j + d f26 with t (c - d) >= d^2
            f = lambda th: (lambda A, B, N: (B + N) + 2 * sq(A) * sq(B))(*phis(th, [i, j, 26]))
        else:
            _, k0, i, j = test    # q f_k0 + t f_i + c f_j + d f26 with c >= d^2/q + d^2/t
            f = lambda th: (lambda A0, A, B, N: N + 2 * sq(B) * (sq(A0) + sq(A)))(*phis(th, [k0, i, j, 26]))
        best = np.inf
        for _ in range(starts):
            r = minimize(f, rng.uniform(-np.pi, np.pi, 6), method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
            best = min(best, r.fun)
        print("NUMERIC capacity %s : min over SO(4) of the sufficient-condition margin = %.3e (%d Nelder-Mead starts, seed %d)%s"
              % (str(test), best, starts, seed, "" if best > -1e-9 else "  <-- NEGATIVE"), flush=True)
    # realizable points: coefficients, G value, min Phi
    if Gjson is None:
        return
    for k in range(npoints):
        pt = L.random_point_on_contact_stratum(T, seed=seed + k, tries=15, margin=1e-3)
        if pt is None:
            info("stratum point %d" % k, "none found"); continue
        g = pt["gaps"]
        if min(g) < -1e-8 or not np.all(np.isfinite(g)) or abs(g[EXCL[-1]]) < 1e-6:
            info("stratum point %d" % k, "rejected (degenerate solver output)"); continue
        coefs = [g[e] for e in EXCL]
        # reconstruct c and check it equals sum coefs f_j
        cvec = np.zeros(14)
        for cf, j in zip(coefs, SRC):
            cvec += cf * np.array([float(x) for x in L.F_VECS[j]])
        recon = float(np.max(np.abs(Pn @ cvec - g)))
        dd = coefs[gvars_slots.index(6)] if 6 in gvars_slots else None
        sc = coefs[-1]
        x = [cf / sc for cf in coefs]          # normalize d = 1 (homogeneous)
        Gv = 0.0; Gabs = 0.0
        for e, cc in Gjson:
            term = float(Fr(cc))
            for slot, ex in enumerate(e[:7]):
                if ex:
                    term *= x[slot] ** ex
            Gv += term; Gabs += abs(term)
        mv, _ = L.min_phi_over_SO4(cvec, restarts=20, seed=seed + k)
        print("NUMERIC stratum point %d: contact %s, coefficients(d=1) %s, |Pc-g| %.1e, G/sum|terms| = %.2e, min Phi (20 BFGS) = %.4e"
              % (k, pt["contact"], [round(v, 4) for v in x], recon, Gv / max(Gabs, 1e-300), mv), flush=True)

numeric_part(T, SRC, EXCL, [("SOC", 2, 12), ("SOC", 2, 15), ("SOC", 12, 15), ("oldH", 3, 15), ("chain", 0, 3, 15)],
             J["G"] if J is not None else None, list(range(7)))

print("SUMMARY: %d EXACT failures" % len(FAILS), FAILS)
sys.exit(1 if FAILS else 0)
