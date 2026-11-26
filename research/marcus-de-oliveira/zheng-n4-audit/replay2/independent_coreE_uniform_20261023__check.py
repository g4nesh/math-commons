"""Second-pass independent replay for the supplement document
independent_coreE_uniform_20261023 ("A seven-contact inequality from a four-coefficient
and chain allocation", 5 pages).

Everything is rebuilt from the paper's definitions through the shared audit library
(../lib/moc4lib.py). Nothing from the supplement is executed, imported or read.

Prints one line per check. Exit status 1 if any EXACT check fails.
Run:  python -I -B independent_coreE_uniform_20261023__check.py [--no-numeric]
"""
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import flint  # noqa: E402
import sympy as sp  # noqa: E402

T0 = time.time()
FAIL = []


def check(name, ok, extra=""):
    print(f"[EXACT] {name}: {'PASS' if ok else 'FAIL'} {extra}".rstrip(), flush=True)
    if not ok:
        FAIL.append(name)


E = [0, 1, 2, 9, 10, 16, 17]
ROLE = ["q", "p", "t", "b", "c", "v", "d"]
GEN = dict(q=0, p=2, t=3, b=12, c=15, v=13, d=26)       # (1) of the document
ROW = dict(q=23, p=15, t=14, b=8, c=4, v=7, d=3)        # rows forming I_7

# ------------------------------------------------------------------ A/(1): the source
check("E equals the paper's Table 5 seven-set no. 17", frozenset(E) == L.AUX7[17])
gens = L.generators(E)
check("generators(E) = {f0,f2,f3,f12,f13,f15,f26}", sorted(gens) == sorted(GEN.values()),
      str([L.f_name(j) for j in gens]))
check("f0,f2,f3,f12,f15,f13 are minor vectors and f26 = h_{4132,-} is not",
      all(L.F_IS_MINOR[GEN[k]] for k in "qptbcv") and not L.F_IS_MINOR[26]
      and L.f_name(26) == "h4132,-")
check("rows [23,15,14,8,4,7,3] of the 7 generators form I_7",
      all(L.F_VALUES[GEN[k]][ROW[k2]] == (1 if k == k2 else 0) for k in ROLE for k2 in ROLE))
check("the 7 generators are linearly independent and the E-rows of P have rank 7",
      L.rank([L.F_VECS[GEN[k]] for k in ROLE]) == 7 and L.rank([L.P[s] for s in E]) == 7)

# ------------------------------------------------------------------ B: group images
fidx = {L.F_VALUES[j]: j for j in range(82)}


def gf(t, j):
    vals = [None] * 24
    for s in range(24):
        vals[t[s]] = L.F_VALUES[j][s]
    return fidx[tuple(vals)]


def images(src, targets):
    return [g for g, t in zip(L.G_ELEMENTS, L._G_MAPS) if tuple(gf(t, j) for j in src) in targets]


for (i, j) in [(2, 12), (2, 15), (12, 15)]:
    r = images((1, 9, 22), {(i, j, 26), (j, i, 26)})
    check(f"(2): some g maps the SOC reference (f1,f9;f22) to (f{i},f{j};f26)", len(r) > 0,
          f"({len(r)} elements, e.g. {r[0] if r else None})")
r = images((0, 6, 10, 54), {(0, 13, 3, 26), (13, 0, 3, 26)})
check("(3): some g maps (f0,f6,f10;f54) to (f0/f13, f13/f0, f3; f26) with the t-role f10 -> f3",
      len(r) > 0, f"({len(r)} elements)")
r = images((2, 0, 12, 70), {(15, 0, 3, 26), (15, 3, 0, 26)})
check("(4): some g maps the chain (f2;f0,f12;f70) to (f15;f0,f3;f26)", len(r) > 0,
      f"({len(r)} elements)")
r = images((2, 0, 12, 70), {(15, 13, 3, 26), (15, 3, 13, 26)})
check("(4'): some g maps the chain to (f15;f13,f3;f26) (q f0 replaced by v f13)", len(r) > 0,
      f"({len(r)} elements)")

# ------------------------------------------------------------------ C: (5), (6), Lemma 1
q, p, t, b, c, v, d, dl, s_, w = sp.symbols("q p t b c v d delta s w", positive=True)
z = lambda e: sp.simplify(sp.together(e)) == 0
# (5): proportional split of two chains
check("(5): split of d, t, c in proportions q/w reproduces the chain condition",
      z((d * q / w) ** 2 / q + (d * q / w) ** 2 / (t * q / w) - (q / w) * (d ** 2 / w + d ** 2 / t)))
# (6): two SOC atoms with shared central coefficient
check("(6): SOC atoms x*y = weight^2 and central sum d^2/s",
      z((d ** 2 * b / s_ ** 2) * b - (d * b / s_) ** 2) and z((d ** 2 * c / s_ ** 2) * c - (d * c / s_) ** 2)
      and z((d ** 2 * b / (b + c) ** 2 + d ** 2 * c / (b + c) ** 2) - d ** 2 / (b + c)))
# Lemma 1
lam = d - dl
A_ = dl + dl ** 2 / v
B_ = lam ** 2 / c
K = A_ + B_
W = q + dl - K
X = W * A_ / K
qF = dl ** 2 / v + X
qC = q - qF
tF = A_ + A_ ** 2 / X
tC = B_ + B_ ** 2 / (W - X)
check("Lemma 1: W - X = W B/K and q_C = B + (W - X)", z(W - X - W * B_ / K) and z(qC - (B_ + W - X)))
check("Lemma 1: q_F v (t_F - delta) = delta^2 (q_F + v + t_F + delta)  [(3) with equality]",
      z(qF * v * (tF - dl) - dl ** 2 * (qF + v + tF + dl)))
check("Lemma 1: lambda^2/q_C + lambda^2/t_C = c  [(4) with equality]", z(lam ** 2 / qC + lam ** 2 / tC - c))
check("Lemma 1: t_F + t_C = K + K^2/W = (q+delta)K/(q+delta-K)",
      z(tF + tC - (K + K ** 2 / W)) and z(K + K ** 2 / W - (q + dl) * K / W))
R_d = d * (v + d) * (q * (c + v) + d * v) / ((c + v) * (q * (c + v) - d ** 2))
sub = {dl: d * v / (c + v)}
check("Lemma 1: delta = dv/(c+v) gives W = q - d^2/(c+v) and total = R_d of (7)",
      z(W.subs(sub) - (q - d ** 2 / (c + v))) and z(((q + dl) * K / W).subs(sub) - R_d))
check("Lemma 1: R_d(q,0,c) = d^2 q/(cq - d^2) = the chain threshold of (4) at v = 0",
      z(R_d.subs(v, 0) - d ** 2 * q / (c * q - d ** 2)))
check("failure of (3) <=> t < T_F when qv > d^2: qv(t-d) - d^2(q+v+t+d) = (qv-d^2)(t-T_F)",
      z(q * v * (t - d) - d ** 2 * (q + v + t + d)
        - (q * v - d ** 2) * (t - d * (q + d) * (v + d) / (q * v - d ** 2))))

# ------------------------------------------------------------------ C: (8), (9)
names = ("q", "p", "t", "b", "c", "v", "d", "z")
ctx = flint.fmpq_mpoly_ctx.get(names, "lex")
gq, gp, gt, gb, gc, gv, gd, gz = ctx.gens()
coef = dict(q=gq, p=gp, t=gt, b=gb, c=gc, v=gv, d=gd)
fq = lambda x: flint.fmpq(x.numerator, x.denominator)
gvec = [ctx.from_dict({}) for _ in range(14)]
for k, sym in coef.items():
    for r_ in range(14):
        if L.F_VECS[GEN[k]][r_] != 0:
            gvec[r_] = gvec[r_] + fq(L.F_VECS[GEN[k]][r_]) * sym
wv = [gvec[r_] + (gz if r_ == 13 else 0) for r_ in range(14)]   # w(z) = g + z kappa
Cw = L.cubic_C(wv).to_dict()
C0 = ctx.from_dict({m: co for m, co in Cw.items() if m[7] == 0})
Qz = ctx.from_dict({tuple(list(m[:7]) + [0]): co for m, co in Cw.items() if m[7] == 1})
check("(8): C(g + z kappa) is affine in z", max(m[7] for m in Cw) == 1)
check("(8): constant term = d{(d+q)(d+v)+dt} + cbp (paper's cubic (A.2), factor 1)",
      C0 == gd * ((gd + gq) * (gd + gv) + gd * gt) + gc * gb * gp)
check("(8): z-coefficient Q = d^2 + dt + cb - ct + cp + bp",
      Qz == gd ** 2 + gd * gt + gc * gb - gc * gt + gc * gp + gb * gp)
ms, Cs = L.cubic_sympy()
H = [[None] * 14 for _ in range(14)]
for i in range(14):
    for j in range(14):
        acc = ctx.from_dict({})
        for mon, co in sp.Poly(sp.diff(Cs, ms[i], ms[j]), *ms).terms():
            ks = [k for k in range(14) if mon[k]]
            co = Fraction(int(co.p), int(co.q))
            acc = acc + (fq(co) if not ks else fq(co) * wv[ks[0]])
        H[i][j] = acc
M = [row[:] for row in H]
prev = ctx.from_dict({(0,) * 8: 1})
sign = 1
rem_ok = True
for k in range(13):
    if M[k][k].is_zero():
        piv = next(r_ for r_ in range(k + 1, 14) if not M[r_][k].is_zero())
        M[k], M[piv] = M[piv], M[k]
        sign = -sign
    for i in range(k + 1, 14):
        for j in range(k + 1, 14):
            qq, rr = divmod(M[i][j] * M[k][k] - M[i][k] * M[k][j], prev)
            rem_ok = rem_ok and rr.is_zero()
            M[i][j] = qq
    prev = M[k][k]
D = (sign * M[13][13]).to_dict()
check("(9): every Bareiss division in det H(g + z kappa) is exact", rem_ok,
      f"(det: {len(D)} terms, z-degree {max(m[7] for m in D)}, {time.time() - T0:.0f}s)")
N = ctx.from_dict({})
for k in range(max(m[7] for m in D) + 1):
    ak = ctx.from_dict({tuple(list(m[:7]) + [0]): co for m, co in D.items() if m[7] == k})
    N = N + ak * (-C0) ** k * Qz ** (4 - k)
S, rr = divmod(N, gd ** 2 * Qz ** 2)
check("(9): Q^4 det H(w(-C/Q)) is divisible by d^2 Q^2", rr.is_zero())
z7 = flint.fmpz_mpoly_ctx.get(names[:7], "lex")
Sd = S.to_dict()
den = 1
for co in Sd.values():
    den = den * int(co.q) // __import__("math").gcd(den, int(co.q))
Si = z7.from_dict({m[:7]: int(co * den) for m, co in Sd.items()})
okG = den == 1 and Si.content() == 4
G = (-(Si / 4)).sqrt() if okG else None
okG = okG and G is not None and -4 * G ** 2 == Si
check("(9): Q^4 det H(w(-C/Q)) = -4 d^2 Q^2 G^2 with G integral", okG)
gq7, gp7, gt7, gb7, gc7, gv7, gd7 = z7.gens()


def at_zero(F, i):
    return z7.from_dict({m: co for m, co in F.to_dict().items() if m[i] == 0})


A0 = gd7 ** 2 + gd7 * (gq7 + gt7 + gv7) + gq7 * gv7
B0 = ((gp7 + gb7 + gd7) * (gq7 + gd7) * (gv7 + gd7) + gp7 * gb7 * (gd7 + gq7 + gv7)
      + gt7 * ((gp7 + gd7) * (gb7 + gd7) - gq7 * gv7))
target = gd7 * A0 * ((gd7 + gq7) * (gd7 + gv7) - gp7 * gb7) * B0
if at_zero(G, 4) == -target:
    G = -G      # fix the sign of the square root
Gd = G.to_dict()
check("(9): G has 408 terms and is homogeneous of degree 8",
      len(Gd) == 408 and all(sum(m) == 8 for m in Gd))


def swap(F, i, j):
    return z7.from_dict({tuple(m[j] if k == i else m[i] if k == j else m[k] for k in range(7)): co
                         for m, co in F.to_dict().items()})


check("(10): G(q,p,t,b,c,v,d) = G(v,p,t,b,c,q,d) = G(q,b,t,p,c,v,d)",
      swap(G, 0, 5) == G and swap(G, 1, 3) == G)
check("c = 0: G = d A {(d+q)(d+v) - pb} B (with the displayed A, B)", at_zero(G, 4) == target)
Bs = (p + b + d) * (q + d) * (v + d) + p * b * (d + q + v) + t * ((p + d) * (b + d) - q * v)
TF = d * (q + d) * (v + d) / (q * v - d ** 2)
check("c = 0: (qv - d^2) B(T_F) = (q+d)(v+d){(p+b)qv + dpb} + pb(d+q+v)(qv-d^2)",
      sp.expand(sp.cancel((q * v - d ** 2) * Bs.subs(t, TF))
                - ((q + d) * (v + d) * ((p + b) * q * v + d * p * b) + p * b * (d + q + v) * (q * v - d ** 2))) == 0)
print(f"[INFO] G degrees in (q,p,t,b,c,v,d): {G.degrees()}", flush=True)

# ------------------------------------------------------------------ C: the 11 tables
P_, B_v, C_, V_, Q_, T_ = sp.symbols("P B C V Q T")
pctx = flint.fmpz_mpoly_ctx.get(("P", "B", "C", "V", "Q", "T"), "lex")
degs = G.degrees()


def R1(qq, vv, cc):
    return (vv + 1) * (qq * (cc + vv) + vv) / ((cc + vv) * (qq * (cc + vv) - 1))


def mainL():
    return dict(p=P_ / (1 + P_), b=B_v / (1 + B_v), c=C_ / (1 + C_))


def mainP():
    dn = (1 + P_) * (1 + B_v + C_)
    return dict(p=1 + P_, b=B_v / dn, c=C_ / dn)


def aux11(m):
    vv = 1 + V_
    qq = vv + Q_
    return dict(v=vv, q=qq, t=R1(qq, vv, m["c"]) * T_ / (1 + T_))


def aux12(m, kind):
    vv = V_ / (1 + V_)
    cc = m["c"]
    if kind == "<":
        return dict(v=vv, q=Q_ / ((cc + vv) * (1 + Q_)), t=T_)
    if kind == ">":
        qq = (1 + Q_) / (cc + vv)
        return dict(v=vv, q=qq, t=R1(qq, vv, cc) * T_ / (1 + T_))
    return dict(v=vv, q=1 / (cc + vv), t=T_)


def domc(kind):
    cc = 1 + C_
    dn = (1 + C_) * (1 + P_ + B_v)
    m = dict(c=cc, p=P_ / dn, b=B_v / dn)
    if kind == "<":
        ww, tt = Q_ / (cc * (1 + Q_)), T_
    elif kind == ">":
        ww = (1 + Q_) / cc
        tt = (ww / (cc * ww - 1)) * T_ / (1 + T_)
    else:
        ww, tt = 1 / cc, T_
    m.update(q=ww / (1 + V_), v=ww * V_ / (1 + V_), t=tt)
    return m


def to_flint(e):
    pe = sp.Poly(sp.expand(e), P_, B_v, C_, V_, Q_, T_)
    return pctx.from_dict({mon: int(co) for mon, co in pe.terms()})


def cleared(mp):
    nd = {}
    for var in ["q", "p", "t", "b", "c", "v"]:
        num, dn = sp.fraction(sp.cancel(sp.together(mp[var])))
        Fn, Fd = to_flint(num), to_flint(dn)
        if all(x < 0 for x in Fd.to_dict().values()):
            Fn, Fd = -Fn, -Fd
        nd[var] = (Fn, Fd)
    one = pctx.from_dict({(0,) * 6: 1})
    cache = {}

    def pw(var, k, which):
        if (var, k, which) not in cache:
            cache[(var, k, which)] = nd[var][which] ** k if k else one
        return cache[(var, k, which)]
    tot = pctx.from_dict({})
    for mon, co in Gd.items():
        term = pctx.from_dict({(0,) * 6: int(co)})
        for i, var in enumerate(["q", "p", "t", "b", "c", "v"]):
            term = term * pw(var, mon[i], 0) * pw(var, degs[i] - mon[i], 1)
        tot = tot + term      # d = 1
    return tot, nd


ROWS15 = [("L", "(11)", 11763, lambda: {**mainL(), **aux11(mainL())}),
          ("P", "(11)", 84446, lambda: {**mainP(), **aux11(mainP())}),
          ("L", "q(c+v)<1", 4181, lambda: {**mainL(), **aux12(mainL(), "<")}),
          ("L", "q(c+v)>1", 11409, lambda: {**mainL(), **aux12(mainL(), ">")}),
          ("L", "q(c+v)=1", 1091, lambda: {**mainL(), **aux12(mainL(), "=")}),
          ("P", "q(c+v)<1", 23837, lambda: {**mainP(), **aux12(mainP(), "<")}),
          ("P", "q(c+v)>1", 93285, lambda: {**mainP(), **aux12(mainP(), ">")}),
          ("P", "q(c+v)=1", 6427, lambda: {**mainP(), **aux12(mainP(), "=")}),
          ("(12)", "cw<1", 27531, lambda: domc("<")),
          ("(12)", "cw>1", 23799, lambda: domc(">")),
          ("(12)", "cw=1", 3955, lambda: domc("="))]
for main, aux, nterms, mk in ROWS15:
    tot, nd = cleared(mk())
    dd = tot.to_dict()
    strict_vars = ([2] if main in ("L", "P") else []) + ([4] if ">" in aux else [])
    only = lambda mon: all(mon[i] == 0 for i in range(6) if i not in strict_vars)
    strict = sorted((sum(m), m, co) for m, co in dd.items() if only(m))
    dens_pos = all(all(x > 0 for x in Fd.to_dict().values()) and any(only(m) for m in Fd.to_dict())
                   for (Fn, Fd) in nd.values())
    mk_str = (f"{strict[0][2]}*" + "".join(f"{n}^{e}" for n, e in zip("PBCVQT", strict[0][1]) if e)) if strict else "none"
    check(f"(15) row {main:4s} {aux:9s}: all coefficients > 0, a monomial in the strictly positive "
          f"parameters only, positive reduced denominators",
          all(co > 0 for co in dd.values()) and bool(strict) and dens_pos,
          f"(terms {len(dd)}, document says {nterms}; lowest strict monomial {mk_str})")

# inverse maps (coverage of the main partition)
pp, bb, cc = sp.symbols("p b c", positive=True)
invP = {P_: pp - 1, B_v: pp * bb / (1 - pp * (bb + cc)), C_: pp * cc / (1 - pp * (bb + cc))}
mP = mainP()
check("row P inverse: P = p-1, B = pb/[1-p(b+c)], C = pc/[1-p(b+c)]",
      all(sp.simplify(mP[k].subs(invP) - val) == 0 for k, val in (("p", pp), ("b", bb), ("c", cc))))
invC = {C_: cc - 1, P_: cc * pp / (1 - cc * (pp + bb)), B_v: cc * bb / (1 - cc * (pp + bb))}
mC = domc("<")
check("dominant-c inverse: C = c-1, P = cp/[1-c(p+b)], B = cb/[1-c(p+b)]",
      all(sp.simplify(mC[k].subs(invC) - val) == 0 for k, val in (("p", pp), ("b", bb), ("c", cc))))
print(f"[INFO] exact part finished in {time.time() - T0:.0f}s", flush=True)

# ------------------------------------------------------------------ D: numerical tests
if "--no-numeric" not in sys.argv:
    import numpy as np
    from scipy.optimize import minimize, least_squares
    FV = {j: np.array(L.F_VALUES[j], float) for j in range(82)}
    ph = lambda j, O: float(FV[j] @ L.leibniz_weights_np(O))
    sq = lambda x: np.sqrt(max(x, 0.0))

    def msmin(fun, restarts, seed):
        r_ = np.random.default_rng(seed)
        best = np.inf
        for _ in range(restarts):
            th = r_.uniform(-np.pi, np.pi, 6)
            res = minimize(lambda x: fun(L.so4_from_params(x)), th, method="Nelder-Mead",
                           options=dict(maxiter=4000, xatol=1e-10, fatol=1e-14))
            res = minimize(lambda x: fun(L.so4_from_params(x)), res.x, method="BFGS")
            best = min(best, res.fun)
        return best
    for (i, j) in [(2, 12), (2, 15), (12, 15)]:
        m = msmin(lambda O: ph(26, O) + 2 * sq(ph(i, O)) * sq(ph(j, O)), 30, 1)
        print(f"[NUMERIC] (2) sharp form min_O Phi26 + 2 sqrt(Phi{i} Phi{j}) (30 starts, seed 1): {m:.2e}", flush=True)
    for qj in (0, 13):
        m = msmin(lambda O: ph(26, O) + 2 * sq(ph(15, O)) * (sq(ph(qj, O)) + sq(ph(3, O))), 30, 2)
        print(f"[NUMERIC] (4) sharp form with f{qj}: min_O Phi26 + 2 sqrt(Phi15)(sqrt(Phi{qj})+sqrt(Phi3)) (30 starts, seed 2): {m:.2e}", flush=True)
    worst = np.inf
    for (qq, vv) in [(1.2, 1.2), (2, 2), (5, 1.5), (10, 10), (3, 0.5), (50, 0.05), (1.01, 1.01), (100, 100)]:
        tt = (qq + 1) * (vv + 1) / (qq * vv - 1)
        worst = min(worst, msmin(lambda O: qq * ph(0, O) + vv * ph(13, O) + tt * ph(3, O) + ph(26, O), 12, 3))
    print(f"[NUMERIC] (3) at the boundary t = (q+1)(v+1)/(qv-1), d = 1, 8 (q,v) pairs x 12 starts: worst min {worst:.2e}", flush=True)

    # realizable points with contact set exactly E and d >= 0.3
    others = [s_ for s_ in range(24) if s_ not in E]

    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])

    def resid(x):
        a, b_, al = unpack(x)
        F = L.F_values_np(a, b_, al).real
        g = F - F[0]
        r_ = [g[s_] for s_ in E[1:]] + [min(0.0, g[s_]) for s_ in others]
        r_.append(min(0.0, g[3] - 0.3))
        r_.append(np.mean(g[others]) - 1.0)
        return np.array(r_)
    Gt = [(tuple(int(e) for e in mm), float(int(co))) for mm, co in Gd.items()]

    def Gval(x):
        terms = [co * np.prod([x[i] ** mm[i] for i in range(7)]) for mm, co in Gt]
        return abs(sum(terms)) / sum(abs(u) for u in terms)
    rng = np.random.default_rng(7)
    found = starts = 0
    wG = 0.0
    wphi = np.inf
    exactE = 0
    while starts < 300 and found < 20:
        starts += 1
        sol = least_squares(resid, rng.normal(size=17), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        if np.max(np.abs(resid(sol.x))) > 1e-10:
            continue
        a, b_, al = unpack(sol.x)
        F = L.F_values_np(a, b_, al).real
        gaps = F - F[0]
        x = [gaps[ROW[k]] for k in ROLE]
        cvec = sum(x[i] * np.array([float(y) for y in L.F_VECS[GEN[k]]]) for i, k in enumerate(ROLE))
        found += 1
        exactE += sorted(int(s_) for s_ in np.where(np.abs(gaps) < 1e-7)[0]) == E
        wG = max(wG, Gval(x))
        val, _ = L.min_phi_over_SO4(cvec, restarts=15, seed=found)
        wphi = min(wphi, val / np.abs(gaps).max())
    print(f"[NUMERIC] realizable points with contacts >= E and d >= 0.3 (seed 7): {found} found in {starts} starts, "
          f"{exactE} with contact set exactly E; worst |G|/sum|terms of G| = {wG:.1e}; "
          f"worst min_SO4 Phi / max gap (15 starts each) = {wphi:.1e}", flush=True)

print(f"[INFO] total time {time.time() - T0:.0f}s; EXACT failures: {len(FAIL)} {FAIL}")
sys.exit(1 if FAIL else 0)
