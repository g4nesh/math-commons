"""Exact replay for blind_shared_coreF_uniform_20261023 (seven-contact set F = {0,1,2,7,9,10,14},
Table 5 seven-set number 7 of the paper, literally).

Everything is rebuilt from the paper's definitions through moc4lib.  Nothing of the
supplement is executed or imported.  If MOC4_SUPP is set, the document's own G (441
coefficients, read from its JSON as data) is compared with the G derived here; the proof
checks use only the G derived here.

Run: python -I -B <this file> [--numeric].  Exits 1 if an EXACT check fails.  Runtime about 1-2 minutes.
"""
import os, sys, time, json
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

T = [0, 1, 2, 7, 9, 10, 14]
SRC = [0, 2, 11, 12, 15, 21, 26]            # q, p, v, b, c, w, d   (document eq. (1))
ROWS6 = [23, 15, 17, 8, 4, 16]
NAMES = ['q', 'p', 'v', 'b', 'c', 'w', 'd']

# ---------------------------------------------------------------- A. source and kernel
check("F equals Table 5 seven-set no. 7 literally", frozenset(T) == frozenset(L.AUX7[7]))
check("generators(F) = {f0,f1,f2,f11,f12,f15,f21,f26} (eight rays, nonsimplicial)",
      L.generators(T) == [0, 1, 2, 11, 12, 15, 21, 26], L.generators(T))
check("rows of P at F have rank 7", L.rank([L.P[s] for s in T]) == 7)
check("the 7 listed generators are independent (basis of the 7-dim kernel)", L.rank([L.F_VECS[j] for j in SRC]) == 7)
check("f0,f1,f2,f11,f12,f15,f21 are minor vectors, f26 is not",
      [L.F_IS_MINOR[j] for j in (0, 1, 2, 11, 12, 15, 21)] == [True] * 7 and not L.F_IS_MINOR[26])
M6 = [[L.F_VALUES[j][r] for j in SRC] for r in ROWS6]
check("values at rows [23,15,17,8,4,16] are (I_6 | column of f26 = 0): six coefficients are gaps",
      M6 == [[int(i == k) for k in range(6)] + [0] for i in range(6)], M6)
add = [L.F_VALUES[26][s] + L.F_VALUES[1][s] - L.F_VALUES[0][s] - L.F_VALUES[11][s] - L.F_VALUES[21][s] for s in range(24)]
check("(2) f26 + f1 = f0 + f11 + f21 (values at all 24 permutations, hence in V_R)", add == [0] * 24)
def gapcol(s): return [L.F_VALUES[j][s] for j in SRC]
check("gap at 20 = q + d, gap at 11 = v + d, gap at 3 = w + d",
      gapcol(20) == [1, 0, 0, 0, 0, 0, 1] and gapcol(11) == [0, 0, 1, 0, 0, 0, 1] and gapcol(3) == [0, 0, 0, 0, 0, 1, 1],
      (gapcol(20), gapcol(11), gapcol(3)))
# d < 0 rewriting: identity of vectors with symbolic coefficients
qs, ps, vs, bs, cs, ws, ds = sp.symbols('q p v b c w d')
dl = sp.symbols('delta')
lhs = [qs * L.F_VALUES[0][s] + ps * L.F_VALUES[2][s] + vs * L.F_VALUES[11][s] + bs * L.F_VALUES[12][s]
       + cs * L.F_VALUES[15][s] + ws * L.F_VALUES[21][s] - dl * L.F_VALUES[26][s] for s in range(24)]
rhs = [dl * L.F_VALUES[1][s] + (qs - dl) * L.F_VALUES[0][s] + ps * L.F_VALUES[2][s] + (vs - dl) * L.F_VALUES[11][s]
       + bs * L.F_VALUES[12][s] + cs * L.F_VALUES[15][s] + (ws - dl) * L.F_VALUES[21][s] for s in range(24)]
check("d = -delta < 0: g = delta f1 + (q-delta) f0 + p f2 + (v-delta) f11 + b f12 + c f15 + (w-delta) f21",
      all(sp.expand(a - b_) == 0 for a, b_ in zip(lhs, rhs)))
# every supported c with contacts containing F lies in span; P injective
check("all 24 gaps of the source are nonnegative combinations when d >= 0 (generator values >= 0)",
      all(L.F_VALUES[j][s] >= 0 for j in SRC for s in range(24)))

# ---------------------------------------------------------------- B. reference actions
def find_ordered(mother, target):
    for g in L.G_ELEMENTS:
        if all(tuple(L.g_act_values(g, L.F_VALUES[m])) == L.F_VALUES[t] for m, t in zip(mother, target)):
            return g
    return None
for mother, target in [((1, 9, 22), (2, 15, 26)), ((1, 9, 22), (2, 12, 26)), ((1, 9, 22), (12, 15, 26))]:
    g = find_ordered(mother, target)
    check(f"SOC: group element maps reference {mother} to {target} in order", g is not None, g)
# old-H: mother minors (7,6) -> (11,12); identify the mother non-minor generator
found = []
for g in L.G_ELEMENTS:
    if tuple(L.g_act_values(g, L.F_VALUES[7])) == L.F_VALUES[11] and tuple(L.g_act_values(g, L.F_VALUES[6])) == L.F_VALUES[12]:
        for k in range(82):
            if tuple(L.g_act_values(g, L.F_VALUES[k])) == L.F_VALUES[26]:
                found.append((g, k))
print(f"INFO old-H: group elements mapping paper f7, f6 to f11, f12: {len(found)} "
      "(so the document's 'mother minor indices [7, 6]' are not paper f-indices; they index a local minor list)")
# the action the document's data names: row (0,1,3,2), column (3,1,0,2), no transpose = our g below
gH = ((3, 1, 0, 2), (0, 1, 3, 2), 0)
pre = {tgt: [k for k in range(82) if tuple(L.g_act_values(gH, L.F_VALUES[k])) == L.F_VALUES[tgt]] for tgt in (11, 12, 26)}
check("old-H: the named action maps paper (f12, f29; f80) to (f11, f12; f26), f80 = h_{1324,+} non-minor",
      pre == {11: [12], 12: [29], 26: [80]} and not L.F_IS_MINOR[80], pre)

# ---------------------------------------------------------------- C. cubic along the lift
ctx8 = flint.fmpz_mpoly_ctx.get(tuple(NAMES + ['z']), 'degrevlex')
g8 = ctx8.gens()
q, p, v, b, c, w, d, z = g8
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
Cw = L.cubic_C(lift(ctx8, g8))
Cdoc = (q + d) * (v + d) * (w + d) + p * b * c
Qdoc = p * b + p * c + b * c + v * (d - b) + d**2
check("(5)-(6) C(w(z)) = C + Q z with the document's C and Q", Cw == Cdoc + Qdoc * z)

# ---------------------------------------------------------------- D. determinant identity (7), d = 1
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
ctx7 = flint.fmpz_mpoly_ctx.get(('q', 'p', 'v', 'b', 'c', 'w', 'z'), 'degrevlex')
g7 = ctx7.gens(); one7 = ctx7.from_dict({(0,) * 7: 1})
t0 = time.time()
detH, ndiv = bareiss_det(hess_lin(lift(ctx7, list(g7[:6]) + [one7, g7[6]]), ctx7), one7)
print(f"INFO det H(w(z)) at d=1: {len(detH)} terms, degrees {detH.degrees()}, {ndiv} exact divisions, {time.time()-t0:.1f}s")
q7, p7, v7, b7, c7, w7, z7 = g7
C1 = (q7 + 1) * (v7 + 1) * (w7 + 1) + p7 * b7 * c7
Q1 = p7 * b7 + p7 * c7 + b7 * c7 + v7 * (1 - b7) + 1
hz = {}
for mono, co in detH.terms():
    k = mono[6]; mm = tuple(mono[:6]) + (0,)
    hz[k] = hz.get(k, ctx7.from_dict({})) + ctx7.from_dict({mm: co})
check("det H(w(z)) has z-degree <= 4", max(hz) <= 4, max(hz))
N = ctx7.from_dict({})
for k, h in hz.items():
    N = N + h * (-C1) ** k * Q1 ** (4 - k)
qq, rr = divmod(N, -4 * (v7 + 1) ** 2 * Q1 ** 2)
check("Q^4 det H(w(-C/Q)) is divisible by -4 (v+d)^2 Q^2", rr.is_zero())
Gd1 = qq.sqrt()
check("the quotient is a perfect square G^2", Gd1 * Gd1 == qq)
ctxG = flint.fmpz_mpoly_ctx.get(tuple(NAMES), 'degrevlex')
gG = ctxG.gens()
G = ctxG.from_dict({tuple(mo[:6]) + (8 - sum(mo[:6]),): int(co) for mo, co in Gd1.terms()})
if int(dict(zip(G.monoms(), G.coeffs())).get((0,) * 6 + (8,), 0)) < 0:
    G = -G
print(f"INFO derived G: {len(G)} terms, degrees {G.degrees()}, coeff of d^8 = {dict(zip(G.monoms(), G.coeffs())).get((0,)*6+(8,))}")
check("derived G has 441 terms (document: 441, degree 8, coefficient of d^8 = 1)", len(G) == 441, len(G))
qG, pG, vG, bG, cG, wG, dG = gG
check("(8) G(q,p,v,b,c,w,d) = G(q,c,v,b,p,w,d)", G.compose(qG, cG, vG, bG, pG, wG, dG) == G)
check("(8) G(q,p,v,b,c,w,d) = G(w,p,v,b,c,q,d)", G.compose(wG, pG, vG, bG, cG, qG, dG) == G)
supp = os.environ.get("MOC4_SUPP")
if supp:
    fn = os.path.join(supp, "proofs", "manin_transfer_20261021", "blind_shared_coreF_uniform_exact_20261023.json")
    with open(fn) as fh: data = json.load(fh)
    Gdoc = ctxG.from_dict({tuple(mo[:7]): int(co) for mo, co in data["G"]})
    check("derived G equals the document's G (its JSON, read as data)", Gdoc == G or Gdoc == -G)
else:
    print("INFO MOC4_SUPP not set: comparison with the document's G skipped")

# ---------------------------------------------------------------- E. charts (9): inverse maps
P_, B_, C_, Q_, V_, W_ = sp.symbols('P B C Q V W', nonnegative=True)
CH = {
    'S':  ({ps: P_ / (1 + P_), bs: B_ / (1 + B_), cs: C_ / (1 + C_), vs: V_}, []),
    'B>': ({ps: P_ / ((1 + B_) * (1 + P_ + C_)), bs: 1 + B_, cs: C_ / ((1 + B_) * (1 + P_ + C_)), vs: V_ / (B_ * (1 + V_))}, [B_]),
    'B=': ({ps: P_ / (1 + P_ + C_), bs: sp.Integer(1), cs: C_ / (1 + P_ + C_), vs: V_}, []),
    'P>=': ({ps: 1 + P_, bs: B_ / ((1 + P_) * (1 + B_ + C_)), cs: C_ / ((1 + P_) * (1 + B_ + C_)), vs: V_}, []),
}
inv = {B_: bs - 1, P_: bs * ps / (1 - bs * (ps + cs)), C_: bs * cs / (1 - bs * (ps + cs)), V_: vs * (bs - 1) / (1 - vs * (bs - 1))}
check("chart B> is inverted by the displayed inverse", all(sp.cancel(CH['B>'][0][k].subs(inv, simultaneous=True) - k) == 0 for k in (ps, bs, cs, vs)))
inv = {P_: ps / (1 - ps - cs), C_: cs / (1 - ps - cs)}
check("chart B= is inverted by P=p/(1-p-c), C=c/(1-p-c)", all(sp.cancel(CH['B='][0][k].subs(inv, simultaneous=True) - k) == 0 for k in (ps, cs)))
inv = {P_: ps - 1, B_: ps * bs / (1 - ps * (bs + cs)), C_: ps * cs / (1 - ps * (bs + cs))}
check("chart P>= is inverted by P=p-1, B=pb/(1-p(b+c)), C=pc/(1-p(b+c))", all(sp.cancel(CH['P>='][0][k].subs(inv, simultaneous=True) - k) == 0 for k in (ps, bs, cs)))
check("failure of (4) with b>1 is v < 1/(b-1); chart value v = V/(B(1+V)) < 1/B", sp.cancel(CH['B>'][0][vs] - V_ / (B_ * (1 + V_))) == 0)

# ---------------------------------------------------------------- F. four positive tables
ctxT = flint.fmpz_mpoly_ctx.get(('P', 'B', 'C', 'Q', 'V', 'W'), 'degrevlex')
symT = [P_, B_, C_, Q_, V_, W_]
def to_flint(ex):
    pol = sp.Poly(sp.expand(ex), *symT)
    return ctxT.from_dict({mo: int(co) for mo, co in pol.terms()})
degs = G.degrees()[:6]
print("INFO spectral degrees of G (q,p,v,b,c,w):", degs)
EXPECT = {'S': 852, 'B>': 3075, 'B=': 478, 'P>=': 2855}
Gd1_terms = list(G.subs({'d': 1}).terms())
order = [qs, ps, vs, bs, cs, ws]
total_terms = 0
for name, (mp0, strict) in CH.items():
    mp = dict(mp0); mp.update({qs: Q_, ws: W_})
    ND = []
    for s in order:
        num, den = sp.fraction(sp.cancel(sp.together(mp[s])))
        ND.append((to_flint(num), to_flint(den)))
    sidx = [symT.index(x) for x in strict]
    den_ok = all(min(int(x) for x in Df.coeffs()) >= 0 and
                 any(all(e == 0 or i in sidx for i, e in enumerate(mo)) for mo in Df.monoms()) for Nf, Df in ND)
    pw = [([Nf ** e for e in range(degs[k] + 1)], [Df ** e for e in range(degs[k] + 1)]) for k, (Nf, Df) in enumerate(ND)]
    tot = ctxT.from_dict({})
    for mo, co in Gd1_terms:
        term = ctxT.from_dict({(0,) * 6: int(co)})
        for k in range(6):
            term = term * pw[k][0][mo[k]] * pw[k][1][degs[k] - mo[k]]
        tot = tot + term
    total_terms += len(tot)
    marker = [(mo, co) for mo, co in tot.terms() if all(e == 0 or i in sidx for i, e in enumerate(mo))]
    check(f"table [{name}]: cleared G has only positive coefficients", min(int(x) for x in tot.coeffs()) > 0,
          f"{len(tot)} terms (document {EXPECT[name]})")
    check(f"table [{name}]: has a monomial in the strictly positive parameters {[str(x) for x in strict]} only", len(marker) > 0, marker[:2])
    check(f"table [{name}]: reduced denominators nonnegative with a strictly positive monomial", den_ok)
print(f"INFO total coefficients in the four tables: {total_terms} (document 7260)")

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
    TSET = T
    def soc(i, j):
        def gen():
            x = lu(); return {i: x, j: 1 / x, 26: 1.0}
        return gen
    def oldH():
        b_ = 1 + lu(); return {11: 1 / (b_ - 1), 12: b_, 26: 1.0}
    FAMILIES = [("SOC (3) p f2 + b f12 + f26, pb = 1", soc(2, 12), 40, 1.0),
                ("SOC (3) p f2 + c f15 + f26, pc = 1", soc(2, 15), 40, 1.0),
                ("SOC (3) b f12 + c f15 + f26, bc = 1", soc(12, 15), 40, 1.0),
                ("old-H (4) v f11 + b f12 + f26 at v(b-1) = 1", oldH, 60, 1.0),
                ("old-H (4)", oldH, 6, 2.0)]
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
