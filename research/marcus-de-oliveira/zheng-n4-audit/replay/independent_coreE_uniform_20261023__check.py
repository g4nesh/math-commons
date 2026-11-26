"""Exact replay for independent_coreE_uniform_20261023 (seven-contact set E = {0,1,2,9,10,16,17},
Table 5 seven-set number 17 of the paper, literally).

Everything is rebuilt from the paper's definitions through the shared library moc4lib
(P, f_j, the cubic (A.2), the group G).  Nothing of the supplement is executed or imported.
If MOC4_SUPP is set, the document's own G (408 coefficients, read from its JSON as data) is
compared with the G derived here; the proof checks use only the G derived here.

Run: python -I -B <this file> [--numeric].  Exits 1 if an EXACT check fails.
Runtime: about 2-3 minutes (the 14 x 14 symbolic determinant takes about 45 s).
"""
import os, sys, time, json, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
import flint

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""), flush=True)
    if not ok:
        FAIL.append(name)

E = [0, 1, 2, 9, 10, 16, 17]
SRC = [0, 2, 3, 12, 15, 13, 26]           # coefficients q, p, t, b, c, v, d (document eq. (1))
ROWS = [23, 15, 14, 8, 4, 7, 3]
NAMES = ['q', 'p', 't', 'b', 'c', 'v', 'd']

# ---------------------------------------------------------------- A. source and contact cone
check("E equals Table 5 seven-set no. 17 literally", frozenset(E) == frozenset(L.AUX7[17]))
check("generators(E) = {f0,f2,f3,f12,f13,f15,f26}", L.generators(E) == sorted(SRC), L.generators(E))
check("rows of P at E have rank 7 (kernel of dimension 7)", L.rank([L.P[s] for s in E]) == 7)
check("f0,f2,f3,f12,f15,f13 are minor vectors and f26 = h_{4132,-} is not",
      [L.F_IS_MINOR[j] for j in SRC] == [True] * 6 + [False], [L.f_name(j) for j in SRC])
M7 = [[L.F_VALUES[j][r] for j in SRC] for r in ROWS]
check("values of the 7 generators at rows [23,15,14,8,4,7,3] form I_7 (so coefficients = gaps)",
      M7 == [[int(i == k) for k in range(7)] for i in range(7)])
check("the 7 generators are independent in V_R", L.rank([L.F_VECS[j] for j in SRC]) == 7)

# ---------------------------------------------------------------- B. reference actions
def find_ordered(mother, target):
    for g in L.G_ELEMENTS:
        if all(tuple(L.g_act_values(g, L.F_VALUES[m])) == L.F_VALUES[t] for m, t in zip(mother, target)):
            return g
    return None
for mother, target, what in [((1, 9, 22), (2, 12, 26), "SOC (p,b)"), ((1, 9, 22), (2, 15, 26), "SOC (p,c)"),
                             ((1, 9, 22), (12, 15, 26), "SOC (b,c)"),
                             ((0, 6, 10, 54), (0, 13, 3, 26), "four-coefficient (q,v,t)"),
                             ((2, 0, 12, 70), (15, 0, 3, 26), "chain (c;q,t)"),
                             ((2, 0, 12, 70), (15, 13, 3, 26), "chain (c;v,t)")]:
    g = find_ordered(mother, target)
    check(f"group element maps reference {mother} to {target} in order [{what}]", g is not None, g)

# ---------------------------------------------------------------- C. cubic along the lift
ctx8 = flint.fmpz_mpoly_ctx.get(tuple(NAMES + ['z']), 'degrevlex')
g8 = ctx8.gens()
q, p, t, b, c, v, d, z = g8
def lift(ctx, gens):
    m = []
    for k in range(14):
        e = ctx.from_dict({})
        for gi, j in enumerate(SRC):
            x = Fr(L.F_VECS[j][k]); assert x.denominator == 1
            if x: e = e + int(x) * gens[gi]
        if k == 13: e = e + gens[7]
        m.append(e)
    return m
m8 = lift(ctx8, g8)
Cw = L.cubic_C(m8)
Cdoc = d * ((d + q) * (d + v) + d * t) + c * b * p
Qdoc = d**2 + d*t + c*b - c*t + c*p + b*p
check("(8) C(w(z)) = C + Q z with the document's C and Q (paper cubic (A.2), paper coordinates)",
      Cw == Cdoc + Qdoc * z)
check("kappa direction: C(kappa)=0 and grad C(kappa)=0 explain linearity in z", Cw.degrees()[7] == 1)

# ---------------------------------------------------------------- D. Hessian determinant, d = 1
ms, Cs = L.cubic_sympy()
def hess_lin(m, ctx):
    H = []
    for i in range(14):
        row = []
        for j in range(14):
            ex = sp.expand(sp.diff(Cs, ms[i], ms[j]))
            e = ctx.from_dict({})
            for k in range(14):
                co = ex.coeff(ms[k])
                if co != 0: e = e + int(co) * m[k]
            row.append(e)
        H.append(row)
    return H
def bareiss_det(Min, one):
    M = [r[:] for r in Min]; n = len(M); sign = 1; prev = one; ndiv = 0
    for k in range(n - 1):
        if M[k][k].is_zero():
            r = next((i for i in range(k + 1, n) if not M[i][k].is_zero()), None)
            if r is None: return one * 0, ndiv
            M[k], M[r] = M[r], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                qq, rr = divmod(M[k][k] * M[i][j] - M[i][k] * M[k][j], prev)
                if not rr.is_zero(): raise RuntimeError("inexact Bareiss division")
                M[i][j] = qq; ndiv += 1
        prev = M[k][k]
    return sign * M[n - 1][n - 1], ndiv

ctx7 = flint.fmpz_mpoly_ctx.get(('q', 'p', 't', 'b', 'c', 'v', 'z'), 'degrevlex')
g7 = ctx7.gens()
one7 = ctx7.from_dict({(0,) * 7: 1})
gens_d1 = list(g7[:6]) + [one7, g7[6]]
m7 = lift(ctx7, gens_d1)
t0 = time.time()
detH, ndiv = bareiss_det(hess_lin(m7, ctx7), one7)
print(f"INFO det H(w(z)) at d=1: {len(detH)} terms, degrees {detH.degrees()}, {ndiv} exact divisions, {time.time()-t0:.1f}s")
# split in powers of z and substitute z = -C/Q (cleared by Q^4)
q7, p7, t7, b7, c7, v7, z7 = g7
C1 = ((1 + q7) * (1 + v7) + t7) + c7 * b7 * p7
Q1 = 1 + t7 + c7 * b7 - c7 * t7 + c7 * p7 + b7 * p7
hz = {}
for mono, co in detH.terms():
    k = mono[6]; mm = tuple(mono[:6]) + (0,)
    hz[k] = hz.get(k, ctx7.from_dict({})) + ctx7.from_dict({mm: co})
check("det H(w(z)) has z-degree <= 4 (H(kappa) has rank 4)", max(hz) <= 4, max(hz))
N = ctx7.from_dict({})
for k, h in hz.items():
    N = N + h * (-C1) ** k * Q1 ** (4 - k)
qq, rr = divmod(N, -4 * Q1 ** 2)
check("Q^4 det H(w(-C/Q)) is divisible by -4 Q^2", rr.is_zero())
Gd1 = qq.sqrt()        # raises if not a perfect square
check("Q^4 det H(w(-C/Q)) / (-4 Q^2) is a perfect square G^2 (d = 1)", Gd1 * Gd1 == qq)
# re-homogenize G to degree 8 in (q,p,t,b,c,v,d)
ctxG = flint.fmpz_mpoly_ctx.get(tuple(NAMES), 'degrevlex')
gG = ctxG.gens()
terms = {}
for mono, co in Gd1.terms():
    deg = sum(mono[:6]); assert deg <= 8
    terms[tuple(mono[:6]) + (8 - deg,)] = int(co)
G = ctxG.from_dict(terms)
if int(dict(zip(G.monoms(), G.coeffs())).get((0,) * 6 + (8,), 0)) < 0:
    G = -G
print(f"INFO derived G: {len(G)} terms, degrees {G.degrees()}, coefficient of d^8 = "
      f"{dict(zip(G.monoms(), G.coeffs())).get((0,)*6+(8,))}")
check("derived G has 408 terms (document: 408 terms, homogeneous of degree 8)", len(G) == 408, len(G))
qG, pG, tG, bG, cG, vG, dG = gG
check("(10) G(q,p,t,b,c,v,d) = G(v,p,t,b,c,q,d)", G.compose(vG, pG, tG, bG, cG, qG, dG) == G)
check("(10) G(q,p,t,b,c,v,d) = G(q,b,t,p,c,v,d)", G.compose(qG, bG, tG, pG, cG, vG, dG) == G)
supp = os.environ.get("MOC4_SUPP")
if supp:
    fn = os.path.join(supp, "proofs", "manin_transfer_20261021", "independent_coreE_uniform_exact_20261023.json")
    with open(fn) as fh: data = json.load(fh)
    Gdoc = ctxG.from_dict({tuple(mo[:7]): int(co) for mo, co in data["G"]})
    check("derived G equals the document's G (its JSON, read as data)", Gdoc == G or Gdoc == -G)
else:
    print("INFO MOC4_SUPP not set: comparison with the document's G skipped")

# ---------------------------------------------------------------- E. c = 0 factor and threshold
A_ = dG**2 + dG * (qG + tG + vG) + qG * vG
B_ = (pG + bG + dG) * (qG + dG) * (vG + dG) + pG * bG * (dG + qG + vG) + tG * ((pG + dG) * (bG + dG) - qG * vG)
Gc0 = G.subs({'c': 0})
check("c = 0 factorization G = d A {(d+q)(d+v) - pb} B (up to sign)",
      Gc0 == dG * A_ * ((dG + qG) * (dG + vG) - pG * bG) * B_ or Gc0 == -(dG * A_ * ((dG + qG) * (dG + vG) - pG * bG) * B_))
qs, ps, ts, bs, cs, vs, ds = sp.symbols('q p t b c v d')
Bs = (ps + bs + ds) * (qs + ds) * (vs + ds) + ps * bs * (ds + qs + vs) + ts * ((ps + ds) * (bs + ds) - qs * vs)
TF = ds * (qs + ds) * (vs + ds) / (qs * vs - ds**2)
lhs = sp.cancel((qs * vs - ds**2) * Bs.subs(ts, TF))
rhs = (qs + ds) * (vs + ds) * ((ps + bs) * qs * vs + ds * ps * bs) + ps * bs * (ds + qs + vs) * (qs * vs - ds**2)
check("threshold identity (qv-d^2) B(T_F) = (q+d)(v+d){(p+b)qv+dpb} + pb(d+q+v)(qv-d^2)", sp.expand(lhs - rhs) == 0)
check("failure of (3) with qv > d^2 is t < T_F", sp.expand(
    (qs*vs*(ts-ds) - ds**2*(qs+vs+ts+ds)) - ((qs*vs - ds**2)*ts - ds*(qs+ds)*(vs+ds))) == 0)

# ---------------------------------------------------------------- F. Lemma 1 and sufficient tests
dl, ql, vl, cl = sp.symbols('delta q v c', positive=True)
lam = ds - dl; A1 = dl + dl**2 / vl; B1 = lam**2 / cl; K1 = A1 + B1; W1 = ql + dl - K1
X1 = W1 * A1 / K1
qF = dl**2 / vl + X1; qC = ql - qF; tF = A1 + A1**2 / X1; tC = B1 + B1**2 / (W1 - X1)
check("Lemma 1: qF v (tF - delta) = delta^2 (qF + v + tF + delta)",
      sp.cancel(qF * vl * (tF - dl) - dl**2 * (qF + vl + tF + dl)) == 0)
check("Lemma 1: lambda^2/qC + lambda^2/tC = c", sp.cancel(lam**2 / qC + lam**2 / tC - cl) == 0)
check("Lemma 1: tF + tC = K + K^2/W = (q+delta)K/W", sp.cancel(tF + tC - (ql + dl) * K1 / W1) == 0
      and sp.cancel(tF + tC - K1 - K1**2 / W1) == 0)
check("Lemma 1: qC = B + (W - X) and W - X = W B / K", sp.cancel(qC - B1 - (W1 - X1)) == 0 and sp.cancel(W1 - X1 - W1 * B1 / K1) == 0)
sub = {dl: ds * vl / (cl + vl)}
check("Lemma 1: at delta = dv/(c+v), W = q - d^2/(c+v)", sp.cancel(W1.subs(sub) - (ql - ds**2 / (cl + vl))) == 0)
Rd = ds * (vl + ds) * (ql * (cl + vl) + ds * vl) / ((cl + vl) * (ql * (cl + vl) - ds**2))
check("Lemma 1: at delta = dv/(c+v), (q+delta)K/W = R_d(q,v,c) of (7)", sp.cancel(((ql + dl) * K1 / W1).subs(sub) - Rd) == 0)
check("Lemma 1: R_d(q,0,c) = d^2 q/(cq - d^2) = chain threshold", sp.cancel(Rd.subs(vl, 0) - ds**2 * ql / (cl * ql - ds**2)) == 0)
# (6) and (5): additive splitting
pp_, bb_, cc_, ss_ = sp.symbols('p b c s', positive=True)
check("(6) splitting: central sum d^2 b/s^2 + d^2 c/s^2 = d^2/s with s=b+c, SOC atoms tight",
      sp.simplify((ds*bb_/(bb_+cc_))**2 / bb_ - ds**2*bb_/(bb_+cc_)**2) == 0)
ww = sp.symbols('w', positive=True)
check("(5) splitting: chain atoms with weights q/w, v/w sum to d^2/w + d^2/t",
      sp.cancel((ds*ql/(ql+vl))**2/ql + (ds*ql/(ql+vl))**2/(ts*ql/(ql+vl)) + (ds*vl/(ql+vl))**2/vl
                + (ds*vl/(ql+vl))**2/(ts*vl/(ql+vl)) - ds**2/(ql+vl) - ds**2/ts) == 0)

# ---------------------------------------------------------------- G. coverage (inverse maps)
P_, B_s, C_, V_, Q_, T_ = sp.symbols('P B C V Q T', nonnegative=True)
# row P inverse
pmap = {ps: 1 + P_, bs: B_s / ((1 + P_) * (1 + B_s + C_)), cs: C_ / ((1 + P_) * (1 + B_s + C_))}
inv = {P_: ps - 1, B_s: ps * bs / (1 - ps * (bs + cs)), C_: ps * cs / (1 - ps * (bs + cs))}
check("row P is inverted by P=p-1, B=pb/(1-p(b+c)), C=pc/(1-p(b+c))",
      all(sp.cancel(pmap[k].subs(inv, simultaneous=True) - k) == 0 for k in pmap))
dmap = {cs: 1 + C_, ps: P_ / ((1 + C_) * (1 + P_ + B_s)), bs: B_s / ((1 + C_) * (1 + P_ + B_s))}
inv = {C_: cs - 1, P_: cs * ps / (1 - cs * (ps + bs)), B_s: cs * bs / (1 - cs * (ps + bs))}
check("dominant-c map (13) is inverted by C=c-1, P=cp/(1-c(p+b)), B=cb/(1-c(p+b))",
      all(sp.cancel(dmap[k].subs(inv, simultaneous=True) - k) == 0 for k in dmap))
check("(14) middle row: failure of (5) for w,t>0 is t < w/(cw-1) when cw>1",
      sp.expand((cs*ww*ts - ww - ts) - ((cs*ww - 1)*ts - ww)) == 0)

# ---------------------------------------------------------------- H. eleven positive tables
ctxT = flint.fmpz_mpoly_ctx.get(('P', 'B', 'C', 'V', 'Q', 'T'), 'degrevlex')
symT = [P_, B_s, C_, V_, Q_, T_]
def to_flint(ex):
    ex = sp.expand(ex)
    pol = sp.Poly(ex, *symT)
    return ctxT.from_dict({mo: int(co) for mo, co in pol.terms()})
def R1(qx, vx, cx):
    return (vx + 1) * (qx * (cx + vx) + vx) / ((cx + vx) * (qx * (cx + vx) - 1))
MAIN = {
    'L': {ps: P_ / (1 + P_), bs: B_s / (1 + B_s), cs: C_ / (1 + C_)},
    'P': {ps: 1 + P_, bs: B_s / ((1 + P_) * (1 + B_s + C_)), cs: C_ / ((1 + P_) * (1 + B_s + C_))},
}
def aux_maps(cx):
    out = {}
    vx = 1 + V_; qx = vx + Q_
    out['(11)'] = ({qs: qx, vs: vx, ts: R1(qx, vx, cx) * T_ / (1 + T_)}, [C_])
    vx = V_ / (1 + V_)
    out['q(c+v)<1'] = ({vs: vx, qs: Q_ / ((cx + vx) * (1 + Q_)), ts: T_}, [C_])
    qx = (1 + Q_) / (cx + vx)
    out['q(c+v)>1'] = ({vs: vx, qs: qx, ts: R1(qx, vx, cx) * T_ / (1 + T_)}, [C_, Q_])
    out['q(c+v)=1'] = ({vs: vx, qs: 1 / (cx + vx), ts: T_}, [C_])
    return out
charts = []
for row in ('L', 'P'):
    for name, (amap, strict) in aux_maps(MAIN[row][cs]).items():
        mp = dict(MAIN[row]); mp.update(amap); charts.append((row + ' ' + name, mp, strict))
cx = 1 + C_
for name, wx, tx, strict in [('cw<1', Q_ / (cx * (1 + Q_)), T_, []),
                             ('cw>1', (1 + Q_) / cx, (1 + Q_) / (cx * Q_) * T_ / (1 + T_), [Q_]),
                             ('cw=1', 1 / cx, T_, [])]:
    mp = dict(dmap); mp.update({qs: wx / (1 + V_), vs: wx * V_ / (1 + V_), ts: tx})
    charts.append(('dominant-c ' + name, mp, strict))
check("(14) middle row: t = w/(cw-1) T/(1+T) with w=(1+Q)/c, cw-1 = Q",
      sp.cancel((((1 + Q_) / cx) / (cx * ((1 + Q_) / cx) - 1)) - (1 + Q_) / (cx * Q_)) == 0)
degs = G.degrees()[:6]
print("INFO spectral degrees of G (q,p,t,b,c,v):", degs)
EXPECT = {'L (11)': 11763, 'P (11)': 84446, 'L q(c+v)<1': 4181, 'L q(c+v)>1': 11409, 'L q(c+v)=1': 1091,
          'P q(c+v)<1': 23837, 'P q(c+v)>1': 93285, 'P q(c+v)=1': 6427,
          'dominant-c cw<1': 27531, 'dominant-c cw>1': 23799, 'dominant-c cw=1': 3955}
Gd1_terms = [(mo, co) for mo, co in G.subs({'d': 1}).terms()]
order = [qs, ps, ts, bs, cs, vs]
for name, mp, strict in charts:
    t1 = time.time()
    ND = []
    for k, s in enumerate(order):
        fr = sp.cancel(sp.together(mp[s]))
        num, den = sp.fraction(fr)
        Nf, Df = to_flint(num), to_flint(den)
        # make the denominator have positive leading behaviour
        ND.append((Nf, Df))
    # denominators: nonnegative coefficients and a positive monomial in strict variables only
    den_ok = True
    for Nf, Df in ND:
        cs_ = [int(x) for x in Df.coeffs()]
        if min(cs_) < 0: den_ok = False
        strict_idx = [symT.index(sv) for sv in strict]
        if not any(all(e == 0 or i in strict_idx for i, e in enumerate(mo)) for mo in Df.monoms()): den_ok = False
    # also every numerator must be nonnegative on the domain, check coefficients
    num_ok = all(min(int(x) for x in Nf.coeffs()) >= 0 for Nf, Df in ND if len(Nf) > 0)
    pw = []
    for k, (Nf, Df) in enumerate(ND):
        pw.append(([Nf ** e for e in range(degs[k] + 1)], [Df ** e for e in range(degs[k] + 1)]))
    tot = ctxT.from_dict({})
    for mo, co in Gd1_terms:
        term = ctxT.from_dict({(0,) * 6: int(co)})
        for k in range(6):
            e = mo[k]
            term = term * pw[k][0][e] * pw[k][1][degs[k] - e]
        tot = tot + term
    coeffs = [int(x) for x in tot.coeffs()]
    allpos = min(coeffs) > 0
    strict_idx = [symT.index(sv) for sv in strict]
    marker = [(mo, co) for mo, co in tot.terms() if all(e == 0 or i in strict_idx for i, e in enumerate(mo))]
    check(f"table [{name}]: cleared G has only positive coefficients", allpos, f"{len(tot)} terms (document {EXPECT[name]}), {time.time()-t1:.1f}s")
    check(f"table [{name}]: has a monomial in the strictly positive parameters {[str(s) for s in strict]} only",
          len(marker) > 0, marker[:2])
    check(f"table [{name}]: every reduced denominator has nonnegative coefficients and a strictly positive monomial", den_ok)
    print(f"INFO table [{name}]: numerators have nonnegative coefficients: {num_ok}; term count matches document: {len(tot) == EXPECT[name]}")

# ---------------------------------------------------------------- NUMERIC (only with --numeric)
if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    Pn = L.P_numpy()
    FV = {j: np.array([float(x) for x in L.F_VECS[j]]) for j in range(82)}
    rng = np.random.default_rng(2026)
    def minphi(co, restarts=8):
        vals = Pn @ sum(w * FV[j] for j, w in co.items())
        f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
        return min(minimize(f, rng.uniform(-np.pi, np.pi, 6), method='BFGS').fun for _ in range(restarts)) / np.abs(vals).max()
    lu = lambda: float(np.exp(rng.uniform(-2.5, 2.5)))
    def run(name, gen, N, scale=1.0):
        worst = (np.inf, None)
        for _ in range(N):
            co = gen(); co = {k: (w * scale if k == 26 else w) for k, w in co.items()}
            m = minphi(co)
            if m < worst[0]: worst = (m, co)
        print(f"NUMERIC {name}{'' if scale == 1 else f' [d scaled by {scale}: sensitivity test]'}: N={N}, 8 BFGS restarts, "
              f"worst min Phi/max|Pc| = {worst[0]:.2e} at {({k: round(w, 4) for k, w in worst[1].items()})}", flush=True)
    TSET = E
    def soc(i, j):
        def gen():
            x = lu(); return {i: x, j: 1 / x, 26: 1.0}
        return gen
    def four():
        while True:
            q_, v_ = lu(), lu()
            if q_ * v_ > 1: return {0: q_, 13: v_, 3: (q_ + 1) * (v_ + 1) / (q_ * v_ - 1), 26: 1.0}
    def chain(lead):
        def gen():
            q_, t_ = lu(), lu(); return {lead: q_, 3: t_, 15: 1 / q_ + 1 / t_, 26: 1.0}
        return gen
    def lemma1():
        while True:
            q_, v_, c_ = lu(), lu(), lu()
            if q_ * (c_ + v_) > 1: break
        R = (v_ + 1) * (q_ * (c_ + v_) + v_) / ((c_ + v_) * (q_ * (c_ + v_) - 1))
        return {0: q_, 13: v_, 15: c_, 3: R, 2: rng.uniform(0, 0.3), 12: rng.uniform(0, 0.3), 26: 1.0}
    FAMILIES = [("SOC (2) p f2 + b f12 + f26, pb = 1", soc(2, 12), 40, 1.0),
                ("SOC (2) p f2 + c f15 + f26, pc = 1", soc(2, 15), 40, 1.0),
                ("SOC (2) b f12 + c f15 + f26, bc = 1", soc(12, 15), 40, 1.0),
                ("SOC (2) p f2 + b f12 + f26", soc(2, 12), 6, 2.0),
                ("four (3) q f0 + v f13 + t f3 + f26 at equality", four, 40, 1.0),
                ("four (3)", four, 6, 2.0),
                ("chain (4) q f0 + t f3 + c f15 + f26 at equality", chain(0), 40, 1.0),
                ("chain (4) with v f13", chain(13), 40, 1.0),
                ("chain (4)", chain(0), 6, 2.0),
                ("Lemma 1 at t = R_1(q,v,c), p,b in [0,0.3]", lemma1, 40, 1.0)]
    for name, gen, N, sc in FAMILIES:
        run(name, gen, N, sc)
    Gterms = [(tuple(int(e) for e in mo), float(int(co))) for mo, co in G.terms()]
    def Gnum(x):
        ts = [co * np.prod([x[i] ** mo[i] for i in range(7)]) for mo, co in Gterms]
        return sum(ts) / sum(abs(t) for t in ts)
    Fv = np.array([[L.F_VALUES[j][s] for j in SRC] for s in range(24)], float)
    for seed in range(2):
        pt = L.random_point_on_contact_stratum(TSET, seed=seed, tries=40, tol=1e-10, margin=0.01)
        if pt is None:
            print(f"NUMERIC stratum seed {seed}: no point found"); continue
        g = pt['gaps']; coef = np.linalg.lstsq(Fv, g, rcond=None)[0]
        c_ = np.linalg.lstsq(Pn, g, rcond=None)[0]
        vals = Pn @ c_
        f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
        m = min(minimize(f, rng.uniform(-np.pi, np.pi, 6), method='BFGS').fun for _ in range(10)) / np.abs(vals).max()
        print(f"NUMERIC realizable point (stratum search seed {seed}, margin 0.01): contact set {pt['contact']}, "
              f"coefficients {np.round(coef, 4).tolist()}, G/sum|terms| = {Gnum(coef):.1e}, min Phi/max|Pc| = {m:.1e}", flush=True)

print("SUMMARY", "ALL EXACT CHECKS PASS" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
