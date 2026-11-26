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
# blind_moc4_coreC_uniform_20261023 : seven-contact set T = {0,1,2,7,9,16,17}
# =================================================================================
T = {0, 1, 2, 7, 9, 16, 17}
SRC = [0, 2, 3, 5, 12, 15, 26]          # q a t u b c d
EXCL = [23, 15, 14, 10, 8, 4, 3]
t0 = time.time()
check("T equals Table 5 seven-set AUX7[11]", frozenset(T) == L.AUX7[11])
check("generators(T) == source list", L.generators(T) == SRC, L.generators(T))
check("rank of P rows at T is 7", L.rank([L.P[s] for s in sorted(T)]) == 7)
check("exclusive rows form I7", [[L.F_VALUES[j][e] for e in EXCL] for j in SRC] == [[int(i == k) for k in range(7)] for i in range(7)])
check("first six generators minor, f26 non-minor", [L.F_IS_MINOR[j] for j in SRC] == [True] * 6 + [False])

# (2) ordering involution: some g in G fixes T and maps f_src[k] -> f_img[k]
IMG = [0, 15, 5, 3, 12, 2, 26]
found = []
for g in L.G_ELEMENTS:
    if L.g_act_set(g, T) != frozenset(T):
        continue
    if all(L.g_act_c(g, L.F_VECS[SRC[k]]) == L.F_VECS[IMG[k]] for k in range(7)):
        found.append(g)
check("(2) a group element fixes T and sends the source list to [0,15,5,3,12,2,26]", len(found) > 0, "%d elements" % len(found))

names = ["q", "a", "t", "u", "b", "c", "z"]
ctx = make_ctx(names)
q, a, t, u, b, c, z = ctx.gens()
one = const(ctx, 1)
w = lift(ctx, [q, a, t, u, b, c, one], SRC, z)
parts = split_z(ctx, L.cubic_C(w), 6)
check("C(w(z)) is affine in z", len(parts) == 2)
C0, Q1 = parts
check("(8) C0 = q + (t+1)(u+1) + abc (d=1)", C0 == q + (t + 1) * (u + 1) + a * b * c)
check("(8) Q0 = ab+ac+bc-au-tc+tu+t+u+1 (d=1)", Q1 == a * b + a * c + b * c - a * u - t * c + t * u + t + u + 1)

J = supp_json("blind_moc4_coreC_uniform_exact_20261023.json")
if J is None:
    info("G8", "MOC4_SUPP not set: G-dependent checks skipped")
else:
    G = jpoly(ctx, J["G"], keep=[0, 1, 2, 3, 4, 5, 7])
    check("JSON G8 has 312 terms", len(J["G"]) == 312, len(J["G"]))
    ctxh = make_ctx(["q", "a", "t", "u", "b", "c", "d"])
    qh, ah, th, uh, bh, ch, dh = ctxh.gens()
    Gh = jpoly(ctxh, J["G"], keep=[0, 1, 2, 3, 4, 5, 6])
    check("G8 homogeneous of degree 8", all(sum(k) == 8 for k in Gh.to_dict()))
    fl = Gh.factor()
    info("factor(G8)", "content %s, %d factors, degrees %s" % (fl[0], len(fl[1]), [f.total_degree() for f, e in fl[1]]))
    check("G8 irreducible over Q", len(fl[1]) == 1 and fl[1][0][1] == 1)
    if not SKIP_HEAVY:
        H = hessian_at(ctx, w)
        detH = det_dp(ctx, H)
        Dk = split_z(ctx, detH, 6)
        info("det H(w(z))", "terms %d, z-degree %d, %.1fs" % (len(detH), len(Dk) - 1, time.time() - t0))
        lhs = zero(ctx)
        for e, Dke in enumerate(Dk):
            lhs = lhs + Dke * (-C0) ** e * Q1 ** (4 - e)
        check("(9) Q0^4 det H(w(-C0/Q0)) = -4 d^2 Q0^2 G8^2 (own Hessian, d=1)", lhs == -4 * Q1 ** 2 * G ** 2)
    # (14) a = c = 0
    Gac0 = Gh.compose(qh, zero(ctxh), th, uh, bh, zero(ctxh), dh, ctx=ctxh)
    W4 = jpoly(ctxh, J["a_c_zero_factor"], keep=[0, 1, 2, 3, 4, 5, 6])
    check("(14) G8|a=c=0 = d^2 {qd+(t+d)(u+d)} W4", Gac0 == dh ** 2 * (qh * dh + (th + dh) * (uh + dh)) * W4)
    dW = W4.to_dict()
    check("(14) W4 has only positive coefficients and contains d^4 with coefficient 1", all(v > 0 for v in dW.values()) and dW.get((0, 0, 0, 0, 0, 0, 4)) == 1, "terms %d" % len(dW))

    # ---------------- the 32 chart certificates ----------------
    cx = make_ctx(["Q", "A", "T", "U", "B", "C"])
    Qv, Av, Tv, Uv, Bv, Cv = cx.gens()
    o = const(cx, 1)
    ctx6 = make_ctx(["q", "a", "t", "u", "b", "c"])
    G6 = ctx6.from_dict({k[:6]: v for k, v in G.to_dict().items()})
    def mul(x, y): return (x[0] * y[0], x[1] * y[1])
    def addc(x, k): return (x[0] + k * x[1], x[1])    # x + k
    def inv(x): return (x[1], x[0])
    def chain_t(qq, cc, T_):
        # t = q T / ((q c - 1)(1 + T)) with q = qq, c = cc as (num, den)
        qc1 = (qq[0] * cc[0] - qq[1] * cc[1], qq[1] * cc[1])
        return (qq[0] * qc1[1] * T_, qq[1] * qc1[0] * (1 + T_))
    charts = []   # (name, [q,a,t,u,b,c] as (num,den), strict set)
    S = {"Q": 0, "A": 1, "T": 2, "U": 3, "B": 4, "C": 5}
    def st(*ns): return {S[n] for n in ns}
    mainsA = {
        "R0": [(Av, 1 + Av + Cv), (Bv, 1 + Bv), (Av + Cv, 1 + Av + Cv)],
        "R1": [(Av, (1 + Bv) * (1 + 2 * Av + Cv)), (1 + Bv, o), (Av + Cv, (1 + Bv) * (1 + 2 * Av + Cv))],
    }
    for mn, (ma, mb, mc) in mainsA.items():
        qa_ = (ma[1], ma[0]); qc_ = (mc[1], mc[0])               # 1/a, 1/c
        q0 = (Qv * mc[1], mc[0] * (1 + Qv))                         # Q/(c(1+Q))
        q1 = (ma[0] * mc[1] + Qv * mc[0] * ma[1], ma[0] * mc[0] * (1 + Qv))   # (a+Qc)/(ac(1+Q)) -- uses a,c as fractions
        # careful: (a + Q c)/(a c (1+Q)) with a = an/ad, c = cn/cd: numerator (an cd + Q cn ad)/(ad cd), denominator an cn (1+Q)/(ad cd)
        q2 = (ma[1] * (1 + Qv), ma[0])                              # (1+Q)/a
        charts += [(mn + "_cell0", [q0, ma, (Tv, o), (Uv, o), mb, mc], st("A")),
                   (mn + "_cell1", [q1, ma, chain_t(q1, mc, Tv), (Uv, o), mb, mc], st("Q", "A", "C")),
                   (mn + "_cell2", [q2, ma, chain_t(q2, mc, Tv), chain_t(q2, ma, Uv), mb, mc], st("Q", "A")),
                   (mn + "_qc1", [qc_, ma, (Tv, o), (Uv, o), mb, mc], st("A")),
                   (mn + "_qa1", [qa_, ma, chain_t(qa_, mc, Tv), (Uv, o), mb, mc], st("A", "C"))]
    # R2 and c = 1 (a > 0)
    def r2like(prefix, ma, mb, mc, tlow, strict_base, low_q=True):
        qa_ = (ma[1], ma[0])
        qlow = (Qv, 1 + Qv)
        qmid = (ma[0] + Qv * ma[1], ma[0] * (1 + Qv))               # (a+Q)/(a(1+Q))
        qhi = (ma[1] * (1 + Qv), ma[0])
        return [(prefix + "_qlow", [qlow, ma, tlow, (Uv, o), mb, mc], strict_base),
                (prefix + "_q1", [(o, o), ma, tlow, (Uv, o), mb, mc], strict_base),
                (prefix + "_qmid", [qmid, ma, chain_t(qmid, mc, Tv), (Uv, o), mb, mc], strict_base | st("Q")),
                (prefix + "_qhi", [qhi, ma, chain_t(qhi, mc, Tv), chain_t(qhi, ma, Uv), mb, mc], strict_base | st("Q")),
                (prefix + "_qa1", [qa_, ma, chain_t(qa_, mc, Tv), (Uv, o), mb, mc], strict_base)]
    charts += r2like("R2", (Av, (1 + Cv) * (1 + Av + Bv)), (Bv, (1 + Cv) * (1 + Av + Bv)), (1 + Cv, o), (Tv, Cv * (1 + Tv)), st("A", "C"))
    charts += r2like("c1", (Av, 1 + Av + Bv), (Bv, 1 + Av + Bv), (o, o), (Tv, o), st("A"))
    # a = 0 wall, c > 0
    z0 = (zero(cx), o)
    for rn, (mb, mc) in {"r0": ((Bv, 1 + Bv), (Cv, 1 + Cv)), "r1": ((1 + Bv, o), (Cv, (1 + Bv) * (1 + Cv)))}.items():
        ql = (Qv * mc[1], mc[0] * (1 + Qv)); qh_ = (mc[1] * (1 + Qv), mc[0]); qw = (mc[1], mc[0])
        charts += [("a0_" + rn + "_low", [ql, z0, (Tv, o), (Uv, o), mb, mc], st("C")),
                   ("a0_" + rn + "_hi", [qh_, z0, chain_t(qh_, mc, Tv), (Uv, o), mb, mc], st("C", "Q")),
                   ("a0_" + rn + "_wall", [qw, z0, (Tv, o), (Uv, o), mb, mc], st("C"))]
    mb, mc = (Bv, (1 + Cv) * (1 + Bv)), (1 + Cv, o)
    charts += [("a0_r2_low", [(Qv, 1 + Qv), z0, (Tv, Cv * (1 + Tv)), (Uv, o), mb, mc], st("C")),
               ("a0_r2_hi", [(1 + Qv, o), z0, chain_t((1 + Qv, o), mc, Tv), (Uv, o), mb, mc], st("C", "Q")),
               ("a0_r2_wall", [(o, o), z0, (Tv, Cv * (1 + Tv)), (Uv, o), mb, mc], st("C"))]
    mb, mc = (Bv, 1 + Bv), (o, o)
    charts += [("a0_c1_low", [(Qv, 1 + Qv), z0, (Tv, o), (Uv, o), mb, mc], set()),
               ("a0_c1_hi", [(1 + Qv, o), z0, chain_t((1 + Qv, o), mc, Tv), (Uv, o), mb, mc], st("Q")),
               ("a0_c1_wall", [(o, o), z0, (Tv, o), (Uv, o), mb, mc], set())]
    check("thirty-two charts built", len(charts) == 32, len(charts))
    for name, subs, strict in charts:
        t1 = time.time()
        dens_ok = all(sign_profile(Dd, strict)[0] and sign_profile(Dd, strict)[2] for (Nn, Dd) in subs)
        check("chart %s: every coordinate denominator has nonnegative coefficients and a strict monomial" % name, dens_ok)
        if SKIP_HEAVY and name.split("_")[1] in ("cell1", "cell2", "qhi", "qmid") and not name.startswith("a0"):
            info("table " + name, "skipped (--skip-heavy)"); continue
        R = subs_clear(G6, cx, subs)
        pos, neg, stt, nt = sign_profile(R, strict)
        if neg and not pos:
            pos, stt = True, sign_profile(-R, strict)[2]       # sign convention of G8 is arbitrary
        check("(table) %s: one sign and strict monomial present" % name, pos and stt, "terms %d (%.1fs)" % (nt, time.time() - t1))

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

numeric_part(T, SRC, EXCL, [("SOC", 2, 12), ("SOC", 2, 15), ("SOC", 12, 15), ("oldH", 3, 15), ("oldH", 5, 2), ("chain", 0, 3, 15), ("chain", 0, 5, 2)],
             J["G"] if J is not None else None, list(range(7)))

print("SUMMARY: %d EXACT failures" % len(FAILS), FAILS)
sys.exit(1 if FAILS else 0)
