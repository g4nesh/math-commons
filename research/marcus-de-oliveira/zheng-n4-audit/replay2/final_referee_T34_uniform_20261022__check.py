"""Second-pass independent replay for final_referee_T34_uniform_20261022 (paper case T_35).

Everything is rebuilt from the paper's definitions through ../lib/moc4lib.py.  Nothing
from the supplement is executed or imported.  If MOC4_SUPP is set, the supplement's stored
degree-nine polynomial G is loaded as DATA (json.load) only to compare it with ours.

Usage:  python -I -B final_referee_T34_uniform_20261022__check.py [--no-numeric]
Prints one line per check; exits 1 if an EXACT check fails.
"""
import sys, os, json, time, random, importlib.util
from pathlib import Path
import flint
import sympy as S

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("moc4lib", HERE.parent / "lib" / "moc4lib.py")
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)

OK = True
def report(name, cond, extra=""):
    global OK
    OK &= bool(cond)
    print(("PASS " if cond else "FAIL ") + "EXACT " + name + (("  " + extra) if extra else ""), flush=True)

t0 = time.time()
# ------------------------------------------------------------------ 1. contact set, source
T = {2, 3, 4, 6, 9, 12}
report("document contact set == paper T_35 (Appendix E), literally", L.CONTACT_SETS[35] == frozenset(T))
gen = L.generators(T)
report("generators(T_35) == [0,1,2,5,7,9,53,67,68]", gen == [0, 1, 2, 5, 7, 9, 53, 67, 68])
report("only f68 is non-minor", [j for j in gen if not L.F_IS_MINOR[j]] == [68])
report("rank of the nine generators is 8", L.rank([L.F_VECS[j] for j in gen]) == 8)
ns = L.nullspace([[L.F_VECS[j][k] for j in gen] for k in range(14)])
rel = {j: x for j, x in zip(gen, ns[0])}
report("unique relation is f68 + f2 = f0 + f9 + f67",
       len(ns) == 1 and all(rel[j] * rel[68] == {68: 1, 2: 1, 0: -1, 9: -1, 67: -1}.get(j, 0) for j in gen))

# ------------------------------------------------------------------ 2. cubic along g + h kappa
names = ('q', 'p', 'j', 'c', 'e', 'd', 'b', 'a', 'h')
ctx = flint.fmpz_mpoly_ctx.get(names, 'deglex')
q, p, j, c, e, d, b, a, h = ctx.gens()
zero = ctx.from_dict({}); one = ctx.from_dict({(0,) * 9: 1})
coef = {67: q, 5: p, 53: j, 9: c, 1: e, 68: d, 7: b, 0: a}
w = [zero] * 14
for k, var in coef.items():
    for i in range(14):
        x = L.F_VECS[k][i]
        assert x.denominator == 1
        if x:
            w[i] = w[i] + int(x) * var
w[13] = w[13] + h            # kappa = e_13
Cw = L.cubic_C(w)
Cstar = Cw.subs({'h': 0}); hQ = Cw - Cstar
Q = hQ.subs({'h': 1})
report("C(g + h kappa) is affine in h", Cw.degrees()[8] == 1 and hQ == Q * h)
report("C* has only positive coefficients and contains d^3",
       all(int(x) > 0 for x in Cstar.coeffs()) and (0, 0, 0, 0, 0, 3, 0, 0, 0) in Cstar.monoms())
report("Q == c + ce + cj - cp + 1 + e + j + ej (homogenized with d)",
       Q == c * d + c * e + c * j - c * p + d * d + d * e + d * j + e * j)

# ------------------------------------------------------------------ 3. full 14x14 Hessian determinant
mctx = flint.fmpz_mpoly_ctx.get(tuple('m%d' % i for i in range(14)), 'deglex')
Cm = L.cubic_C(list(mctx.gens()))
def ev(lin):
    out = zero
    for mon, cf in lin.to_dict().items():
        out = out + int(cf) * w[mon.index(1)]
    return out
H = [[ev(Cm.derivative(r).derivative(s)) for s in range(14)] for r in range(14)]

def det_dp(M):
    """Exact determinant by dynamic programming over column subsets (Laplace, row by row)."""
    n = len(M)
    order = sorted(range(n), key=lambda r: sum(1 for s in range(n) if not M[r][s].is_zero()))
    cur = {0: one}
    for r in order:
        nxt = {}
        for Sm, Mm in cur.items():
            for s in range(n):
                if Sm >> s & 1 or M[r][s].is_zero():
                    continue
                term = Mm * M[r][s]
                if bin(Sm >> (s + 1)).count('1') % 2:
                    term = -term
                Tm = Sm | (1 << s)
                nxt[Tm] = nxt[Tm] + term if Tm in nxt else term
        cur = {Sm: Mm for Sm, Mm in nxt.items() if not Mm.is_zero()}
    out = cur[(1 << n) - 1]
    inv = sum(1 for i in range(n) for k in range(i + 1, n) if order[i] > order[k])
    return -out if inv % 2 else out

det = det_dp(H)
rng = random.Random(5)
agree = True
for _ in range(3):
    vals = [rng.randint(-5, 5) for _ in range(9)]
    agree &= L.det([[int(H[r][s](*vals)) for s in range(14)] for r in range(14)]) == int(det(*vals))
report("symbolic Hessian determinant agrees with Fraction det at 3 integer points", agree)
Dk = {}
for mon, cf in det.to_dict().items():
    Dk.setdefault(int(mon[8]), {})[tuple(mon[:8]) + (0,)] = int(cf)
R = zero
for k, dct in Dk.items():
    R = R + ctx.from_dict(dct) * (-Cstar) ** k * Q ** (4 - k)
quo, rem = divmod(-R, 4 * Q * Q)
G = quo.sqrt() if rem.is_zero() else None
if G is not None and G(*([1] * 9)) < 0:
    G = -G
report("(2) Q^4 det Hess C(g - (C*/Q) kappa) == -4 Q^2 G^2 for a polynomial G",
       G is not None and G * G == quo)
report("G: 1475 terms, homogeneous of degree 9, degree 2 in p",
       len(G) == 1475 and all(sum(m) == 9 for m in G.monoms()) and G.degrees()[1] == 2)
supp = os.environ.get("MOC4_SUPP")
if supp:
    fn = Path(supp) / "proofs/manin_transfer_20261021/final_referee_T34_full_radial_20261022.json"
    rad = json.load(open(fn))
    terms = next(r['terms'] for r in rad['radial_factors'] if r['degree'] == 9)
    Gs = ctx.from_dict({tuple(ex): int(co) for ex, co in terms})
    print(("PASS " if Gs == G else "INFO-MISMATCH ") + "DATA our G equals the supplement's stored G (same sign)", flush=True)
print("time %.1fs" % (time.time() - t0), flush=True)

# ------------------------------------------------------------------ 4. positivity tables (d = 1)
ORIG = ('q', 'p', 'j', 'c', 'e', 'b', 'a')
CH = ('B', 'U', 'V', 'X', 'Y', 'C', 'E', 'J')
ALL = ORIG + CH
cx = flint.fmpz_mpoly_ctx.get(ALL, 'deglex')
gg = dict(zip(ALL, cx.gens()))
one2 = cx.from_dict({(0,) * len(ALL): 1}); zero2 = cx.from_dict({})
dct = {}
for m, cf in G.to_dict().items():
    qq_, pp_, jj_, cc_, ee_, dd_, bb_, aa_, hh_ = [int(t) for t in m]
    key = (qq_, pp_, jj_, cc_, ee_, bb_, aa_) + (0,) * len(CH)
    dct[key] = dct.get(key, 0) + int(cf)
G7 = cx.from_dict(dct)
q, p, j, c, e, b, a = [gg[n] for n in ORIG]
B, U, V, X, Y, C, E, J = [gg[n] for n in CH]

def split(P, v):
    i = ALL.index(v); out = {}
    for m, cf in P.to_dict().items():
        m2 = list(m); k = int(m2[i]); m2[i] = 0
        out.setdefault(k, {})[tuple(m2)] = cf
    return {k: cx.from_dict(dd) for k, dd in out.items()}

mG = {v: G7.degrees()[ALL.index(v)] for v in ORIG}

def subst(P, v, num, den, m):
    """P(v = num/den) * den^m, with m >= deg_v P (a positive clearing factor)."""
    parts = split(P, v)
    R = zero2
    for k, Pk in parts.items():
        R = R + Pk * num ** k * den ** (m - k)
    return R

def nonneg(P): return all(int(x) >= 0 for x in P.coeffs())
def setv(P, **kw): return P.subs(kw)

kk = split(G7, 'p'); k0, k1, k2 = kk.get(0, zero2), kk.get(1, zero2), kk.get(2, zero2)
G00 = setv(k0, e=0, j=0)
report("G00 = G|p=e=j=0 == (1+b)(c+1+b)^2(a+1)^3(q+1)^3", G00 == (1 + b) * (c + 1 + b) ** 2 * (a + 1) ** 3 * (q + 1) ** 3)
report("k0|b=0 == (a+1)(q+1)(e+1+a)^2(c+1)^2(q+j+1)^2",
       setv(k0, b=0) == (a + 1) * (q + 1) * (e + 1 + a) ** 2 * (c + 1) ** 2 * (q + j + 1) ** 2)
Gc0 = setv(G7, c=0)
report("G|c=b=0 == (qp+q+qa+p+pa+1+a)(e+1+a)^2(q+j+1)^2",
       setv(Gc0, b=0) == (q * p + q + q * a + p + p * a + 1 + a) * (e + 1 + a) ** 2 * (q + j + 1) ** 2)
report("G|c=0,p=e=j=0 == (1+b)^3(a+1)^3(q+1)^3", setv(Gc0, p=0, e=0, j=0) == (1 + b) ** 3 * (a + 1) ** 3 * (q + 1) ** 3)

def cases(bnum, bden):
    out = []
    for ca in ('below', 'above'):
        for cq in ('below', 'above'):
            s = {}
            if ca == 'below':
                s['a'] = (X * bden, bnum * (1 + X)); s['e'] = (E, one2)
            else:
                an, ad = (1 + X) * bden, bnum
                s['a'] = (an, ad); s['e'] = ((an + ad) * (bnum + bden) * E, ad * bden * X * (1 + E))
            if cq == 'below':
                s['q'] = (Y * bden, bnum * (1 + Y)); s['j'] = (J, one2)
            else:
                qn, qd = (1 + Y) * bden, bnum
                s['q'] = (qn, qd); s['j'] = ((qn + qd) * (bnum + bden) * J, qd * bden * Y * (1 + J))
            out.append((ca + '/' + cq, s))
    return out

def apply(P, s, bsub, csub=None):
    R = P
    for v in ('a', 'e', 'q', 'j'):
        R = subst(R, v, *s[v], m=mG[v])
    R = subst(R, 'b', *bsub, m=mG['b'])
    if csub:
        R = subst(R, 'c', *csub, m=mG['c'])
    return R

counts7 = []
for lab, s in cases(B, one2):
    K = apply(k0, s, (B, one2), (C, one2))
    ref = setv(setv(apply(G00, s, (B, one2), (C, one2)), E=0), J=0)
    counts7.append(len(K))
    report("(7) table %s: nonnegative coefficients" % lab, nonneg(K), "terms=%d" % len(K))
    report("(7) table %s: (E,J)-free part equals the cleared G00" % lab, setv(setv(K, E=0), J=0) == ref)
report("(7) term counts == 2227, 3268, 3268, 4945", counts7 == [2227, 3268, 3268, 4945])
counts9 = []
for lab, s in cases(B, one2):
    K = apply(Gc0, s, (B, one2))
    ref = setv(setv(apply(setv(Gc0, p=0, e=0, j=0), s, (B, one2)), E=0), J=0)
    counts9.append(len(K))
    report("(9) table %s: nonnegative coefficients (all p)" % lab, nonneg(K), "terms=%d" % len(K))
    report("(9) table %s: (E,J,p)-free part equals cleared (1+b)^3(a+1)^3(q+1)^3" % lab,
           setv(setv(setv(K, E=0), J=0), p=0) == ref)
report("(9) term counts == 1419, 2110, 2110, 3211", counts9 == [1419, 2110, 2110, 3211])
countsb = []
for lab, s in cases(B, 1 + B):
    Ks = [apply(x, s, (B, 1 + B), (C, one2)) for x in (k0, k1, k2)]
    Nn, Dd = C + 1, C
    mid = 2 * Dd * Ks[0] + Nn * Ks[1]
    end = Dd * Dd * Ks[0] + Nn * Dd * Ks[1] + Nn * Nn * Ks[2]
    countsb.append((len(mid), len(end)))
    report("b<=1 (P*=1+1/c) table %s: middle Bernstein numerator nonnegative" % lab, nonneg(mid), "terms=%d" % len(mid))
    report("b<=1 (P*=1+1/c) table %s: endpoint Bernstein numerator nonnegative" % lab, nonneg(end), "terms=%d" % len(end))
report("b<=1 term counts == (3612,3612),(7230,7230),(7230,7230),(12258,12258)",
       countsb == [(3612, 3612), (7230, 7230), (7230, 7230), (12258, 12258)])
bb = 1 + B; Dc = 1 + U + V
sub = {'a': (U + X * Dc, bb * Dc), 'q': (V + Y * Dc, bb * Dc),
       'e': (U * (bb + 1) * (X * Dc + U * (1 + bb)), bb * Dc * Dc * X),
       'j': (V * (bb + 1) * (Y * Dc + V * (1 + bb)), bb * Dc * Dc * Y)}
def apply2(P):
    R = P
    for v in ('a', 'e', 'q', 'j'):
        R = subst(R, v, *sub[v], m=mG[v])
    R = subst(R, 'b', bb, one2, m=mG['b'])
    return subst(R, 'c', C, one2, m=mG['c'])
Ks = [apply2(x) for x in (k0, k1, k2)]
Nn = bb + C * Dc; Dd = C * bb * Dc * Dc
mid = 2 * Dd * Ks[0] + Nn * Ks[1]
end = Dd * Dd * Ks[0] + Nn * Dd * Ks[1] + Nn * Nn * Ks[2]
report("b>=1 chart (15): base k0 nonnegative", nonneg(Ks[0]), "terms=%d" % len(Ks[0]))
report("b>=1 chart (15): middle Bernstein numerator nonnegative", nonneg(mid), "terms=%d" % len(mid))
report("b>=1 chart (15): endpoint Bernstein numerator nonnegative", nonneg(end), "terms=%d" % len(end))
report("b>=1 term counts == 117110, 209489, 289848", [len(Ks[0]), len(mid), len(end)] == [117110, 209489, 289848])
print("time %.1fs" % (time.time() - t0), flush=True)

# ------------------------------------------------------------------ 5. allocation identities (sympy)
a_, b_, dl, e_, r_, U_, V_, X_, Y_, B_, C_, N_, D_, t_, p_, K0, K1, K2, E_ = S.symbols(
    'a b delta e r U V X Y B C N D t p k0 k1 k2 E')
Ea = dl * (b_ + 1) * (a_ + dl) / (a_ * b_ - dl)
z = lambda ex: S.simplify(ex) == 0
report("(11) derivative of E_a", z(S.diff(Ea, dl) - (b_ + 1) * (a_**2 * b_ + 2 * a_ * b_ * dl - dl**2) / (a_ * b_ - dl)**2))
report("(11) E_a(1) = (a+1)(b+1)/(ab-1)", z(Ea.subs(dl, 1) - (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)))
report("(13) holds identically at E = E_a(delta)", z(a_ * (dl * b_) * (Ea - dl) - dl**2 * (a_ + dl * b_ + Ea + dl)))
report("(13) is (5) for the piece rescaled by 1/delta",
       z((a_ / dl * b_ * (Ea / dl - 1) - (a_ / dl + b_ + Ea / dl + 1)) * dl**2 - (a_ * (dl * b_) * (Ea - dl) - dl**2 * (a_ + dl * b_ + Ea + dl))))
report("(5) fails <=> L_e e < 1", S.expand(a_ * b_ * (e_ - 1) - (a_ + b_ + e_ + 1) - (e_ * (a_ * b_ - 1) - (a_ + 1) * (b_ + 1))) == 0)
report("(14) is (4) for the residual piece rescaled by 1/r", z(r_ * (1 / b_ + r_ / C_) - (r_**2 / C_ + r_ / b_)))
Dch = 1 + U_ + V_; bch = 1 + B_; d1 = U_ / Dch; d2 = V_ / Dch
ach = (U_ + X_ * Dch) / (bch * Dch); qch = (V_ + Y_ * Dch) / (bch * Dch)
report("(15) X = ab - delta1, Y = bq - delta2, r = 1 - delta1 - delta2",
       z(ach * bch - d1 - X_) and z(qch * bch - d2 - Y_) and z(1 - d1 - d2 - 1 / Dch))
report("(15) e = E_a(delta1) in chart form",
       z(Ea.subs({a_: ach, b_: bch, dl: d1}) - U_ * (bch + 1) * (X_ * Dch + U_ * (1 + bch)) / (bch * Dch**2 * X_)))
report("(16) P* = r^2/c + r/b = (b + C D)/(C b D^2)", z(1 / (Dch**2 * C_) + 1 / (Dch * bch) - (bch + C_ * Dch) / (C_ * bch * Dch**2)))
bern = K0 * (1 - t_)**2 + (2 * D_ * K0 + N_ * K1) / D_ * t_ * (1 - t_) + (D_**2 * K0 + N_ * D_ * K1 + N_**2 * K2) / D_**2 * t_**2
report("(10) Bernstein coefficients of k0 + k1 p + k2 p^2 on [0, N/D]", z(bern - (K0 + K1 * p_ + K2 * p_**2).subs(p_, t_ * N_ / D_)))
aab = (1 + X_) / B_; eab = (aab + 1) * (B_ + 1) / X_ * E_ / (1 + E_)
report("(8) above-one chart: ab = 1+X and 1 - L_e e = 1/(1+E)",
       z(aab * B_ - 1 - X_) and z(1 - (aab * B_ - 1) / ((aab + 1) * (B_ + 1)) * eab - 1 / (1 + E_)))

# ------------------------------------------------------------------ 6. group transport of cited pointwise inequalities
FI = {v: i for i, v in enumerate(L.F_VALUES)}
def gmap(gm, jdx):
    new = [0] * 24
    for s_ in range(24):
        new[gm[s_]] = L.F_VALUES[jdx][s_]
    return FI[tuple(new)]
def exists(base, target, sym2):
    for gm in L._G_MAPS:
        im = [gmap(gm, x) for x in base]
        if (set(im[:2]) == set(target[:2]) and im[2:] == target[2:]) if sym2 else im == target:
            return True
    return False
report("group maps four-coefficient base (f0,f6;f10;f54) to ({7,0};f1;f68)", exists([0, 6, 10, 54], [7, 0, 1, 68], True))
report("group maps four-coefficient base (f0,f6;f10;f54) to ({7,67};f53;f68)", exists([0, 6, 10, 54], [7, 67, 53, 68], True))
report("group maps chain base (f2;f0,f12;f70) to (f5;f7,f9;f68)", exists([2, 0, 12, 70], [5, 7, 9, 68], False))
report("group maps old-H base (f3,f5;f54) to (f9,f5;f68) in this order", exists([3, 5, 54], [9, 5, 68], False))
print("time %.1fs" % (time.time() - t0), flush=True)

# ------------------------------------------------------------------ 7. NUMERIC tests
if "--no-numeric" not in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    FV = np.array(L.F_VALUES, float)
    def minphi(vals, restarts, rg):
        f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
        return min(minimize(f, rg.uniform(-np.pi, np.pi, 6), method='BFGS').fun for _ in range(restarts))
    rg = np.random.default_rng(7)
    comb = lambda dd: sum(v * FV[k] for k, v in dd.items())
    for name in ('oldH (3)', 'chain (4)', 'four (5)', 'four (6)'):
        ws, ins = [], []
        for _ in range(12):
            if name.startswith('oldH'):
                cc = np.exp(rg.uniform(-2, 2)); pp = 1 + 1 / cc
                dd, di = {9: cc, 5: pp, 68: 1}, {9: cc, 5: 0.97 * pp, 68: 1}
            elif name.startswith('chain'):
                b1, c1 = np.exp(rg.uniform(-2, 2, 2)); pp = 1 / b1 + 1 / c1
                dd, di = {5: pp, 7: b1, 9: c1, 68: 1}, {5: 0.97 * pp, 7: b1, 9: c1, 68: 1}
            else:
                a1, b1 = np.exp(rg.uniform(-1, 2.5, 2))
                if a1 * b1 <= 1.05:
                    a1 = 1.5 / b1
                ee = (a1 + 1) * (b1 + 1) / (a1 * b1 - 1)
                side, aux = (0, 1) if '(5)' in name else (67, 53)
                dd, di = {7: b1, side: a1, aux: ee, 68: 1}, {7: b1, side: a1, aux: 0.97 * ee, 68: 1}
            ws.append(minphi(comb(dd), 12, rg)); ins.append(minphi(comb(di), 12, rg))
        print("NUMERIC pointwise %s on its boundary: 12 coefficient samples x 12 BFGS restarts (seed 7): worst min Phi = %.2e; at 0.97 x threshold worst = %.2e"
              % (name, min(ws), min(ins)), flush=True)
    Gt = [([int(x) for x in m], int(cf)) for m, cf in G.to_dict().items()]
    def Gev(x): return sum(cf * np.prod([v**k for v, k in zip(x, ex)]) for ex, cf in Gt)
    def Gabs(x): return sum(abs(cf) * np.prod([abs(v)**k for v, k in zip(x, ex)]) for ex, cf in Gt)
    A8 = FV[[67, 5, 53, 9, 1, 68, 7, 0]].T
    npts = int(os.environ.get("T34_NPTS", "6")); found = 0; rels = []; phis = []; branch_ok = True; nd = 0
    for sd in range(200):
        if found >= npts:
            break
        pt = L.random_point_on_contact_stratum(sorted(T), seed=7000 + sd, tries=8, margin=0.01)
        if pt is None:
            continue
        found += 1
        gaps = np.array(pt['gaps'], float)
        co = np.linalg.lstsq(A8, gaps, rcond=None)[0]
        mx = np.max(np.abs(co))
        phis.append(minphi(gaps, 10, rg))
        if co[5] > 1e-6 * mx and co.min() > -1e-7 * mx:
            nd += 1
            qv, pv, jv, cv, ev_, dv, bv, av = co / co[5]
            x = [qv, pv, jv, cv, ev_, 1.0, bv, av, 0.0]
            rels.append(abs(Gev(x)) / Gabs(x))
            # which certificate of the document covers this realizable cost?
            if av * bv * (ev_ - 1) >= av + bv + ev_ + 1 or bv * qv * (jv - 1) >= bv + qv + jv + 1:
                br = 'four'
            elif bv <= 1:
                br = 'oldH' if pv >= 1 + 1 / cv else 'UNCOVERED'
            else:
                from scipy.optimize import brentq
                def dsol(aa, ee):
                    if aa <= 0 or ee <= 0:
                        return 0.0
                    f = lambda t: t * (bv + 1) * (aa + t) / (aa * bv - t) - ee
                    return brentq(f, 0, min(1.0, aa * bv) * (1 - 1e-15))
                d1, d2 = dsol(av, ev_), dsol(qv, jv)
                rr = 1 - d1 - d2
                br = 'double-four' if rr <= 0 else ('chain' if pv >= rr**2 / cv + rr / bv else 'UNCOVERED')
            branch_ok &= br != 'UNCOVERED'
            print("NUMERIC realizable point (d>0): coefficients/d = %s, relative |G| = %.1e, covering certificate = %s"
                  % (np.round(co / co[5], 4).tolist(), rels[-1], br), flush=True)
    print("NUMERIC stratum T_35 (lib random_point_on_contact_stratum, margin 0.01, seeds 7000+): %d points, %d with d>0; max relative |G| = %s; every d>0 point covered by a certificate: %s; worst min Phi = %s"
          % (found, nd, ("%.1e" % max(rels)) if rels else None, branch_ok, ("%.2e" % min(phis)) if phis else None), flush=True)

print("time %.1fs" % (time.time() - t0))
print("ALL EXACT CHECKS PASS" if OK else "SOME EXACT CHECK FAILED")
sys.exit(0 if OK else 1)
