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
DOC = "low_pair"

# =================================================================================
# coreD (low pair and uniform) : seven-contact set T = {0,1,2,7,9,16,23}
# =================================================================================
T = {0, 1, 2, 7, 9, 16, 23}
SRC = [2, 3, 5, 11, 12, 15, 26]          # a t u v b c d
EXCL = [15, 14, 10, 17, 8, 4, 3]
t0 = time.time()
check("T equals Table 5 seven-set AUX7[12]", frozenset(T) == L.AUX7[12])
check("generators(T) == source list", L.generators(T) == SRC, L.generators(T))
check("rank of P rows at T is 7", L.rank([L.P[s] for s in sorted(T)]) == 7)
check("exclusive rows form I7", [[L.F_VALUES[j][e] for e in EXCL] for j in SRC] == [[int(i == k) for k in range(7)] for i in range(7)])
check("first six generators minor, f26 non-minor", [L.F_IS_MINOR[j] for j in SRC] == [True] * 6 + [False])
# source S3: g fixes T, f26, and permutes the pairs (f2,f5), (f12,f11), (f15,f3) as pairs
PAIRS = [(2, 5), (12, 11), (15, 3)]
fidx = {tuple(L.F_VECS[j]): j for j in range(82)}
perms_seen = set()
for g in L.G_ELEMENTS:
    if L.g_act_set(g, T) != frozenset(T):
        continue
    im = {j: fidx.get(tuple(L.g_act_c(g, L.F_VECS[j]))) for j in SRC}
    if im[26] != 26:
        continue
    pi = []
    for (m, x) in PAIRS:
        tgt = [k for k, (m2, x2) in enumerate(PAIRS) if im[m] == m2 and im[x] == x2]
        pi.append(tgt[0] if tgt else None)
    if None not in pi:
        perms_seen.add(tuple(pi))
check("source S3: group elements realise all 6 permutations of the (main, opposite auxiliary) pairs", len(perms_seen) == 6, sorted(perms_seen))

names = ["a", "t", "u", "v", "b", "c", "z"]
ctx = make_ctx(names)
a, t, u, v, b, c, z = ctx.gens()
one = const(ctx, 1)
w = lift(ctx, [a, t, u, v, b, c, one], SRC, z)
parts = split_z(ctx, L.cubic_C(w), 6)
check("C(w(v,z)) is affine in z", len(parts) == 2)
C0, Q1 = parts
Nn = a * b * c + t * u + t + u + 1
Rr = t + u + 1
check("low-pair (2): F = N + R v + Q z with N = abc+tu+t+u+1, R = t+u+1", C0 == Nn + Rr * v)
check("uniform (2): C = abc+tu+tv+uv+t+u+v+1", C0 == a * b * c + t * u + t * v + u * v + t + u + v + 1)
check("(2) Q = ab+ac+bc-au-tc-vb+tu+tv+uv+t+u+v+1", Q1 == a * b + a * c + b * c - a * u - t * c - v * b + t * u + t * v + u * v + t + u + v + 1)

J = supp_json("blind_moc4_coreD_low_pair_exact_20261023.json")
KEEP = [0, 1, 2, 3, 4, 5, 7]
if J is None:
    info("G", "MOC4_SUPP not set: G-dependent checks skipped")
    G = None
else:
    G = jpoly(ctx, J["G"], keep=KEEP)
    check("JSON G has 385 terms", len(J["G"]) == 385, len(J["G"]))
    ctxh = make_ctx(["a", "t", "u", "v", "b", "c", "d"])
    Gh = jpoly(ctxh, J["G"], keep=[0, 1, 2, 3, 4, 5, 6])
    check("G homogeneous of degree 9", all(sum(k) == 9 for k in Gh.to_dict()))
    fl = Gh.factor()
    check("G irreducible over Q", len(fl[1]) == 1 and fl[1][0][1] == 1, "content %s" % fl[0])
    info("G degrees (a,t,u,v,b,c,d)", Gh.degrees())

def pseudo_rem_v(P, G2, vi):
    """pseudo-remainder of P by G2 (quadratic in variable vi): returns R with alpha^k P = q G2 + R, deg_v R < 2."""
    gc = split_z(ctx, G2, vi)
    lc = gc[-1]
    dg = len(gc) - 1
    vv = ctx.gens()[vi]
    while True:
        pc = split_z(ctx, P, vi)
        dp_ = len(pc) - 1
        if P.is_zero() or dp_ < dg:
            return P
        P = lc * P - pc[-1] * vv ** (dp_ - dg) * G2

if G is not None and not SKIP_HEAVY:
    H = hessian_at(ctx, w)
    detH = det_dp(ctx, H)
    Dk = split_z(ctx, detH, 6)
    info("det H(w(v,z))", "terms %d, z-degree %d, %.1fs" % (len(detH), len(Dk) - 1, time.time() - t0))
    lhs = zero(ctx)
    for e, Dke in enumerate(Dk):
        lhs = lhs + Dke * (-C0) ** e * Q1 ** (4 - e)
    check("low-pair (3) / uniform (3): Q^4 det H(w(v,-(N+Rv)/Q)) = -4 Q^2 G^2 (own Hessian, d=1)", lhs == -4 * Q1 ** 2 * G ** 2)

if G is not None:
    gv = split_z(ctx, G, 3)
    check("(3) G quadratic in v", len(gv) == 3)
    gam, bet, alp = gv
    F0 = jpoly(ctx, J["Herons"][0], keep=KEEP)
    F1 = jpoly(ctx, J["Herons"][1], keep=KEEP)
    check("F0, F1 have 77 terms each and are v-free", len(F0) == 77 and len(F1) == 77 and F0.degrees()[3] == 0 and F1.degrees()[3] == 0)
    check("(5) beta^2 - 4 alpha gamma = F+ F-", bet ** 2 - 4 * alp * gam == F0 * F1)
    ctxh = make_ctx(["a", "t", "u", "v", "b", "c", "d"])
    for k, Fk in enumerate(J["Herons"]):
        fl = jpoly(ctxh, Fk, keep=[0, 1, 2, 3, 4, 5, 6]).factor()
        check("F%d irreducible octic" % k, len(fl[1]) == 1 and fl[1][0][1] == 1 and fl[1][0][0].total_degree() == 8)
    check("F0 != +-F1 (so F+F- is not a square up to sign)", F0 != F1 and F0 != -F1)
    if not SKIP_HEAVY and DOC == "low_pair":
        I = [i for i in range(14) if i not in (0, 5)]
        dI = det_dp(ctx, H, I, I)
        check("(4) principal block on I = {0..13}\\{0,5} is height-independent", dI.degrees()[6] == 0)
        P9 = jpoly(ctx, J["principal_P9"], keep=KEEP)
        check("(4) det H_I = (v+b+1)^2 P (P from certificate, 675 terms)", dI == (v + b + 1) ** 2 * P9 and len(P9) == 675)
        pc = split_z(ctx, P9, 3)
        check("(4) P quadratic in v", len(pc) == 3)
        p0, p1, p2 = pc
        U1 = alp * p1 - p2 * bet
        U0 = alp * p0 - p2 * gam
        Wp = jpoly(ctx, J["norm_W"], keep=KEEP)
        check("(5) U0^2 - U1 (beta p0 - p1 gamma) = W^2", U0 ** 2 - U1 * (bet * p0 - p1 * gam) == Wp ** 2)
        check("W, H+- are v-free (so coprime to the irreducible, v-dependent G)", Wp.degrees()[3] == 0)
        Uv = U1 * v + U0
        for tf in J["trace_factors"]:
            sgn = tf["sign"]; Hs = jpoly(ctx, tf["square"], keep=KEEP)
            Fe = F1 if tf["octic_index"] == 1 else F0
            S = 2 * alp * U0 - bet * U1 + sgn * 2 * alp * Wp
            check("(6) S_%+d = F_%d H^2" % (sgn, tf["octic_index"]), S == Fe * Hs ** 2 and Hs.degrees()[3] == 0)
            Rm = P9 * S - (Uv + sgn * Wp) ** 2
            check("(6) P S_%+d = (U %+d W)^2 mod G" % (sgn, sgn), pseudo_rem_v(Rm, G, 3).is_zero())
        # bordered 13-minors, forced height
        for added in (0, 5):
            idx = sorted(I + [added])
            d13 = det_dp(ctx, H, idx, idx)
            zc = split_z(ctx, d13, 6)
            kz = len(zc) - 1
            cleared = zero(ctx)
            for e, ce in enumerate(zc):
                cleared = cleared + ce * (-C0) ** e * Q1 ** (kz - e)
            check("bordered 13-minor on I+{%d}: z-degree %d; Q^%d * minor(z=-(N+Rv)/Q) divisible by G" % (added, kz, kz), pseudo_rem_v(cleared, G, 3).is_zero())
        info("rank argument", "second derivative of -4G^2/Q^2 at a simple root of G is -8 G_v^2/Q^2 (PROOF-READ)")

    # ---------------- sixteen sign tables (8) ----------------
    cx = make_ctx(["A", "T", "U", "B", "C"])
    Av, Tv, Uv_, Bv, Cv = cx.gens()
    o = const(cx, 1)
    ctx6 = make_ctx(["a", "t", "u", "b", "c"])
    def to5(P): return ctx6.from_dict({(k[0], k[1], k[2], k[4], k[5]): val for k, val in P.to_dict().items()})
    polys = {"gamma": (to5(gam), +1), "Rbeta-Nalpha": (to5(Rr * bet - Nn * alp), +1),
             "RF0-Nalpha": (to5(Rr * F0 - Nn * alp), -1), "RF1-Nalpha": (to5(Rr * F1 - Nn * alp), -1)}
    check("table polynomials are v-free", all(P.degrees()[3] == 0 for P in (gam, Rr * bet - Nn * alp, Rr * F0 - Nn * alp, Rr * F1 - Nn * alp)))
    for pn, Dd in (("a+c<1", 1 + 2 * Av + Cv), ("a+c=1", 2 * Av + Cv)):
        aa = (Av, Dd); cc = (Av + Cv, Dd)
        for bn, bb in (("b<1", ((Av + Cv) + Bv * Dd, Dd * (1 + Bv))), ("b=1", (o, o))):
            subs = [aa, (Tv, o), (Uv_, o), bb, cc]
            for name, (Pp, sg) in polys.items():
                Rx = subs_clear(Pp, cx, subs)
                pos, neg, stt, nt = sign_profile(Rx if sg > 0 else -Rx, {0})
                check("(8) %s %s %s has sign %+d with a strict monomial (A>0)" % (pn, bn, name, sg), pos and stt, "terms %d" % nt)
                if pn == "a+c<1":
                    k0 = Rx.to_dict().get((0, 0, 0, 0, 0), 0)
                    check("(8) %s %s %s: nonzero constant coefficient (strict also at A = 0, used by coreD zero-main)" % (pn, bn, name), k0 != 0, k0)

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

numeric_part(T, SRC, EXCL, [], J["G"] if J is not None else None, list(range(7)), npoints=3)

print("SUMMARY: %d EXACT failures" % len(FAILS), FAILS)
sys.exit(1 if FAILS else 0)
