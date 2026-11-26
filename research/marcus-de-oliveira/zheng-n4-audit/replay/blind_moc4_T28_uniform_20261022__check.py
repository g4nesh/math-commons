# ---- common header (own code; nothing from the supplement is executed or imported) ----
# Usage: MOC4_SUPP=<supplement root> python -I -B <this file> [--numeric]
# The supplement's JSON files are read as data only (json.load).
import sys, os, json, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
from sympy.polys.rings import ring
from sympy import ZZ
SUPP = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else None)
DATA = os.path.join(SUPP, "proofs", "manin_transfer_20261021") if SUPP else None
NUMERIC = "--numeric" in sys.argv
FAILED = []
def check(label, cond):
    print(("PASS " if cond else "FAIL ") + label, flush=True)
    if not cond:
        FAILED.append(label)
    return cond
def info(*a):
    print("INFO", *a, flush=True)
V = L.F_VALUES

def radial_identity(jsonname, names, idx, npts=10, seed=12345):
    """Own cubic (lib cubic_C, paper (A.2)) on g + z*kappa, g = sum of the named generators.
    Compares C, Q with the JSON's C_terms/Q_terms, and checks
    Q^4 det Hess(g - C/Q kappa) = -4 Q^2 G^2 at npts random rational points (randomized identity test),
    with G the degree-nine factor read from the JSON as data.  Returns the JSON.
    After multiplying by Q^10 both sides are polynomials of degree <= 50 in the 8 coefficients, so by
    Schwartz-Zippel a false identity survives one uniform integer point in [1,1e6]^7 with probability
    <= 5e-5, and all 10 independent points with probability < 1e-42."""
    J = json.load(open(os.path.join(DATA, jsonname)))
    check(f"{jsonname}: variable order {names}", J["variables"][:len(names)] == names)
    syms = sp.symbols(names); z = sp.Symbol("z"); allv = list(syms) + [z]
    def frompoly(terms):
        return sp.Add(*[sp.Rational(c) * sp.Mul(*[v**e for v, e in zip(allv, ex)]) for ex, c in terms])
    m = [sum(s * L.F_VECS[i][k] for s, i in zip(syms, idx)) for k in range(14)]
    m[13] = m[13] + z
    ms, Cexpr = L.cubic_sympy()
    Fm = sp.Poly(sp.expand(Cexpr.subs(dict(zip(ms, m)), simultaneous=True)), z)
    check(f"{jsonname}: F(g + z kappa) is affine in z (F(kappa)=0, grad F(kappa)=0)", Fm.degree() == 1)
    Cm = Fm.coeff_monomial(1); Qm = Fm.coeff_monomial(z)
    check(f"{jsonname}: own C equals the JSON C", sp.expand(Cm - frompoly(J["C_terms"])) == 0)
    check(f"{jsonname}: own Q equals the JSON Q", sp.expand(Qm - frompoly(J["Q_terms"])) == 0)
    Qf = frompoly(J["radial_factors"][0]["terms"]); G = frompoly(J["radial_factors"][1]["terms"])
    check(f"{jsonname}: first radial factor = +-Q", sp.expand(Qf - Qm) == 0 or sp.expand(Qf + Qm) == 0)
    H = sp.hessian(Cexpr, ms); Hf = sp.lambdify(ms, H, "sympy")
    rng = random.Random(seed); ratios = set()
    for _ in range(npts):
        pt = [sp.Integer(rng.randint(1, 10**6)) for _ in syms]   # Schwartz-Zippel: integer points from a set of size 1e6
        pt[names.index("d")] = sp.Integer(1)
        sb = dict(zip(syms, pt)); Qv = Qm.subs(sb)
        if Qv == 0: continue
        zv = -Cm.subs(sb) / Qv
        mv = [sp.Rational(x.subs(sb).subs(z, zv)) for x in m]   # exact; do NOT use nsimplify (it rounds huge rationals)
        import flint
        Hm = Hf(*mv)
        dv = sp.Rational(str(flint.fmpq_mat([[flint.fmpq(int(sp.numer(e)), int(sp.denom(e))) for e in row] for row in Hm.tolist()]).det()))
        Gv = G.subs(sb)
        ratios.add(Qv**4 * dv / (Qv**2 * Gv**2) if Gv != 0 else None)
    info(f"{jsonname}: set of values of Q^4 det H / (Q^2 G^2) over the points: {ratios}")
    check(f"{jsonname}: Q^4 det H(g - C/Q kappa) = -4 Q^2 G^2 at {npts} random integer points", ratios == {-4})
    return J, Cm, Qm, syms

class Eng:
    """Exact rational-substitution engine on a sympy ZZ ring."""
    def __init__(self, R): self.R = R
    def degs(self, P):
        dg = [0] * self.R.ngens
        for mon in P.monoms():
            for i, k in enumerate(mon): dg[i] = max(dg[i], k)
        return dg
    def coeffs_in(self, P, vi):
        out = {}
        for mon, co in P.terms():
            k = mon[vi]; mm = list(mon); mm[vi] = 0
            out[k] = out.get(k, self.R(0)) + self.R({tuple(mm): co})
        return out
    def subst(self, P, sub, dg):
        """numerator of P(sub) times prod den_i^dg_i (sub: index -> (num, den))."""
        R = self.R; cache = {}
        def pw(i, k):
            if (i, k) not in cache:
                n, d = sub[i]; cache[(i, k)] = n**k * d**(dg[i] - k)
            return cache[(i, k)]
        groups = {}
        for mon, co in P.terms():
            key = tuple(mon[i] for i in sorted(sub)); rest = list(mon)
            for i in sub: rest[i] = 0
            groups[key] = groups.get(key, R(0)) + R({tuple(rest): co})
        out = R(0)
        for key, poly in groups.items():
            t = poly
            for i, k in zip(sorted(sub), key): t = t * pw(i, k)
            out += t
        return out
def nonneg(P): return all(c >= 0 for c in P.coeffs())
def load_ring_poly(R, terms, names, vmap):
    P = R(0)
    for ex, co in terms:
        mm = R(int(co))
        for nm, k in zip(names, ex):
            if k: mm *= vmap[nm]**k
        P += mm
    return P
_X = sp.Symbol("x")
def univ_nonneg(coeffs):
    """exact: polynomial sum coeffs[k] x^k is >= 0 on [0, oo) (sign changes only at real roots;
    squarefree factorization, real-root isolation of odd-multiplicity factors)."""
    if all(v >= 0 for v in coeffs.values()): return True
    P = sp.Poly({(k,): v for k, v in coeffs.items()}, _X)
    if P.is_zero: return True
    if P.LC() < 0: return False
    for f, mult in P.sqf_list()[1]:
        if mult % 2 == 1 and any(r > 0 for r in sp.real_roots(f)): return False
    return True
def grouped_nonneg(P, var):
    """group P by all monomials except variable index var; each group must be >= 0 for var >= 0."""
    gl = {}
    for mon, co in P.terms():
        key = tuple(m if i != var else 0 for i, m in enumerate(mon))
        gl.setdefault(key, {}); gl[key][mon[var]] = gl[key].get(mon[var], 0) + co
    bad = [g for g in gl.values() if any(v < 0 for v in g.values())]
    return all(univ_nonneg(g) for g in bad), len(bad)
def fm(x, y): return (x[0] * y[0], x[1] * y[1])
def fa(x, y): return (x[0] * y[1] + y[0] * x[1], x[1] * y[1])
# ---- end header ----

# ===== blind_moc4_T28_uniform_20261022 : paper case T_29 = {0,4,5,6,7,8} =====
T = L.CONTACT_SETS[29]
check("doc contact set {0,4,5,6,7,8} is the paper's T_29 literally", sorted(T) == [0, 4, 5, 6, 7, 8])
gens = L.generators(T)
check("eleven generators 0,1,2,3,4,5,6,33,34,45,46 = lib generators(T_29)", gens == [0, 1, 2, 3, 4, 5, 6, 33, 34, 45, 46])
check("non-minor generators are exactly 34 and 46", [g for g in gens if not L.F_IS_MINOR[g]] == [34, 46])
def vec(combo): return [sum(cc * V[jj][s] for jj, cc in combo) for s in range(24)]
for l, r_ in (([(34, 1), (46, 1)], [(3, 1), (6, 2), (33, 1), (45, 1)]), ([(2, 1), (34, 1)], [(0, 1), (6, 1), (33, 1)]),
              ([(4, 1), (34, 1)], [(1, 1), (3, 1), (6, 1), (33, 1)]), ([(0, 1), (46, 1)], [(2, 1), (3, 1), (6, 1), (45, 1)]),
              ([(1, 1), (46, 1)], [(4, 1), (6, 1), (45, 1)])):
    check(f"identity {l} = {r_} on all 24 permutation values", vec(l) == vec(r_))
b1 = [3, 6, 1, 33, 45, 34, 0, 5]; b2 = [4, 3, 2, 33, 5, 6, 46, 45]
check("basis (1) = generators minus {2,4,46}, rank 8", sorted(b1) == sorted(set(gens) - {2, 4, 46}) and L.rank([list(V[x]) for x in b1]) == 8)
check("basis (2) = generators minus {0,1,34}, rank 8", sorted(b2) == sorted(set(gens) - {0, 1, 34}) and L.rank([list(V[x]) for x in b2]) == 8)

# ---- radial necessity, both bases
J34, C34, Q34, s34 = radial_identity("blind_moc4_T28_bad34_full_radial_20261022.json", ["q", "p", "j", "c", "e", "d", "b", "a"], b1)
pC = sp.Poly(C34, *s34)
check("(1): every coefficient of C is positive and C contains d^3", all(x > 0 for x in pC.coeffs()) and pC.coeff_monomial(s34[5]**3) == 1)
q_, p_, j_, c_, e_, d_, b_, a_ = s34
check("(1): Q = ad+ae+aq+cd+ce+cq+d^2+de+dq-ej+eq", sp.expand(Q34 - (a_*d_ + a_*e_ + a_*q_ + c_*d_ + c_*e_ + c_*q_ + d_**2 + d_*e_ + d_*q_ - e_*j_ + e_*q_)) == 0)
J46, C46, Q46, s46 = radial_identity("blind_moc4_T28_bad46_full_radial_20261022.json", ["p", "r", "a", "b", "c", "q", "d", "sa"], b2)
p_, r_, a_, b_, c_, q_, d_, s_ = s46
M = b_ + c_; Uu = d_ + r_ + s_
check("(2): Q = -M a + U M + d U + r s", sp.expand(Q46 - (-M * a_ + Uu * M + d_ * Uu + r_ * s_)) == 0)
MC = sp.Poly(sp.expand(M * C46 - b_ * c_ * Q46), *s46)
check(f"(2): M C - b c Q has nonnegative coefficients ({len(MC.terms())} terms) and contains b d^3 + c d^3",
      all(x > 0 for x in MC.coeffs()) and MC.coeff_monomial(b_ * d_**3) > 0 and MC.coeff_monomial(c_ * d_**3) > 0)

# ---- nonminor34 branch certificates
R, q, p, j, c, e, b, a, B, Tt, A = ring("q,p,j,c,e,b,a,B,T,A", ZZ)
IQ, IP, IJ, IC, IE, IB, IA = range(7)
eng = Eng(R)
G = load_ring_poly(R, J34["radial_factors"][1]["terms"], J34["variables"],
                   {"q": q, "p": p, "j": j, "c": c, "e": e, "d": R(1), "b": b, "a": a, "z": R(1)})
kk = eng.coeffs_in(G, IJ)
check("G34 has 1376 terms, degree 2 in j", len(G.terms()) == 1376 and max(kk) == 2)
k0, k1, k2 = kk[0], kk[1], kk[2]
F8, rem = k0.div(p + 1)
check("k0 = (p+1) F8 exactly", rem == 0)
Fa = eng.coeffs_in(F8, IA); F0, F1, F2 = Fa[0], Fa.get(1, R(0)), Fa.get(2, R(0))
check("F8 is quadratic in a", max(Fa) == 2)
check("F8 a^0 coefficient: 204 terms, all positive, constant 1",
      len(F0.terms()) == 204 and all(x > 0 for x in F0.coeffs()) and F0.compose([(x, R(0)) for x in (q, p, c, e, b)]) == 1)
check("e=0: k2 = 0 and k1 has nonnegative coefficients", k2.compose(e, R(0)) == 0 and nonneg(k1.compose(e, R(0))))
def Lval(fb_, fp_): return (fb_[0] * fp_[0] - fb_[1] * fp_[1], (fb_[0] + fb_[1]) * (fp_[0] + fp_[1]))
def transform(label, Fs, fb_, fp_, regime):
    dg = [max(x, y, z) for x, y, z in zip(eng.degs(Fs[0]), eng.degs(Fs[1]), eng.degs(Fs[2]))]
    Fz = [eng.subst(x, {IB: fb_, IP: fp_}, dg) for x in Fs]
    if regime == "bp<1": cs = Fz
    else:
        n, d = Lval(fb_, fp_)
        check(f"{label} {regime}: L numerator/denominator nonnegative", nonneg(n) and nonneg(d))
        cs = [Fz[0], d * Fz[1] + 2 * n * Fz[0], d * d * Fz[2] + n * d * Fz[1] + n * n * Fz[0]]
    for i, x in enumerate(cs): check(f"{label} {regime}: transformed coefficient {i} of (4) nonnegative ({len(x.terms())} terms)", nonneg(x))
transform("F8", [F0, F1, F2], (B, R(1)), (Tt, B * (1 + Tt)), "bp<1")
transform("F8", [F0, F1, F2], (B, R(1)), (1 + Tt, B), "bp>1")
Gb = k0 + k1 * b + k2 * b**2
Gba = eng.coeffs_in(Gb, IA); check("G34(b) is quadratic in a", max(Gba) == 2)
H0, H1, H2 = Gba[0], Gba[1], Gba[2]
fb01 = (B, 1 + B)
transform("G34(b), b=B/(1+B)", [H0, H1, H2], fb01, (Tt * (1 + B), B * (1 + Tt)), "bp<1")
transform("G34(b), b=B/(1+B)", [H0, H1, H2], fb01, ((1 + Tt) * (1 + B), B), "bp>1")
H0s = eng.subst(H0, {IB: fb01}, eng.degs(H0))
check("G34(b) a^0 coefficient at b=B/(1+B): nonnegative with constant 1", nonneg(H0s) and H0s.compose([(x, R(0)) for x in (q, p, c, e, B)]) == 1)
H01 = H0.compose(b, R(1))
check("G34(b) a^0 coefficient at b=1: nonnegative with constant 16", nonneg(H01) and H01.compose([(x, R(0)) for x in (q, p, c, e)]) == 16)
dgall = [max(x, y, z) for x, y, z in zip(eng.degs(k0), eng.degs(k1), eng.degs(k2))]
def charts_p_a(fb_):
    fp1 = (Tt * fb_[1], fb_[0] * (1 + Tt)); fa1 = (A, R(1))
    fp2 = ((1 + Tt) * fb_[1], fb_[0]); bp = fm(fb_, fp2)
    fa2 = ((fb_[0] + fb_[1]) * (fp2[0] + fp2[1]) * bp[1] * A, fb_[1] * fp2[1] * (bp[0] - bp[1]) * (1 + A))
    return [("bp<1", fp1, fa1), ("bp>1", fp2, fa2)]
for nm, fp_, fa_ in charts_p_a(fb01):
    check(f"b<=1 {nm}: chart numerators/denominators nonnegative", all(nonneg(x) for x in (*fp_, *fa_)))
    K = [eng.subst(x, {IB: fb01, IP: fp_, IA: fa_}, dgall) for x in (k0, k1, k2)]
    bn, bd = fb01
    check(f"b<=1 {nm}: left middle 2k0 + b k1 >= 0", nonneg(2 * bd * K[0] + bn * K[1]))
    check(f"b<=1 {nm}: cut G34(b) >= 0", nonneg(bd * bd * K[0] + bn * bd * K[1] + bn * bn * K[2]))
    check(f"b<=1 {nm}: right middle 2e k0 + (be+e+1) k1 + 2b(e+1) k2 >= 0", nonneg(bd * 2 * e * K[0] + (bn * e + bd * e + bd) * K[1] + 2 * bn * (e + 1) * K[2]))
    check(f"b<=1 {nm}: endpoint e^2 k0 + (e+1) e k1 + (e+1)^2 k2 >= 0", nonneg(e * e * K[0] + (e + 1) * e * K[1] + (e + 1)**2 * K[2]))
fb1 = (1 + B, R(1))
for nm, fp_, fa_ in charts_p_a(fb1):
    K = [eng.subst(x, {IB: fb1, IP: fp_, IA: fa_}, dgall) for x in (k0, k1, k2)]
    N = fb1[0] + e; Dd = e * fb1[0]
    for lab, P in (("middle", 2 * Dd * K[0] + N * K[1]), ("endpoint", Dd * Dd * K[0] + N * Dd * K[1] + N * N * K[2])):
        if nonneg(P): check(f"b>=1 {nm} {lab} (J=1/e+1/b): coefficientwise nonnegative", True)
        else:
            ok, nb_ = grouped_nonneg(P, IE)
            check(f"b>=1 {nm} {lab} (J=1/e+1/b): nonnegative by grouping in e ({nb_} univariate certificates)", ok)

# ---- nonminor46 branch: Lemma 2 charts and certificates
R2, p, r, a, b, c, q, s, Tt, Lq, U, W, Y = ring("p,r,a,b,c,q,s,T,L,U,W,Y", ZZ)
IP, IR, IA, IB, IC, IQ, IS, IT, IL, IU, IW, IY = range(12)
eng2 = Eng(R2)
G = load_ring_poly(R2, J46["radial_factors"][1]["terms"], J46["variables"],
                   {"p": p, "r": r, "a": a, "b": b, "c": c, "q": q, "d": R2(1), "sa": s, "z": R2(1)})
kk = eng2.coeffs_in(G, IA)
check("G46 has 1549 terms, degree 2 in a", len(G.terms()) == 1549 and max(kk) == 2)
k0, k1, k2 = kk[0], kk[1], kk[2]
check("G46 k0: 612 terms, all positive, constant 1", len(k0.terms()) == 612 and all(x > 0 for x in k0.coeffs()) and k0.compose([(x, R2(0)) for x in (p, r, b, c, q, s)]) == 1)
dgall = [max(x, y, z) for x, y, z in zip(eng2.degs(k0), eng2.degs(k1), eng2.degs(k2))]
# symbolic chart identities: J_chart = right side of (6) at a feasible (beta, delta)
Ts, Ls, Us, Ws = sp.symbols("T L U W", positive=True); th = Ts / (1 + Ts)
def rhs6(p_, b_, c_, q_, b1, d1, d2): return 1 - d1 - d2 + b1**2 / b_ + (1 - b1)**2 / c_ + d1**2 / p_ + d2**2 / q_
pp = 2*th*Us/(1+Us); qq = 2*(1-th)*Ws/(1+Ws)
check("HH: J = 1 + L - (p+q)/4 equals (6) at beta1=theta, delta=(p/2,q/2); delta_i < beta_i",
      sp.simplify(rhs6(pp, th/Ls, (1-th)/Ls, qq, th, pp/2, qq/2) - (1 + Ls - (pp+qq)/4)) == 0)
pp = 2*th*Us/(1+Us); qq = 2*(1-th)*(1+Ws); cc = 2*(1-th)*(1+Ws)/(2*Ls*(1+Ws)+Ws)
check("HC (and CH by exchange): J = L + 1/2 + theta/(2(1+U)) equals (6) at delta=(p/2, beta2), beta2 <= q/2",
      sp.simplify(rhs6(pp, th/Ls, cc, qq, th, pp/2, 1-th) - (Ls + sp.Rational(1, 2) + th/(2*(1+Us)))) == 0)
pp = 2*th*(1+Us); qq = 2*(1-th)*(1+Us+Ws); cc = 2*(1-th)*(1+Us)*(1+Us+Ws)/(2*Ls*(1+Us)*(1+Us+Ws)+Ws)
check("CC (and exchange): J = L + 1/(2(1+U)) equals (6) at delta = beta",
      sp.simplify(rhs6(pp, th/Ls, cc, qq, th, th, 1-th) - (Ls + 1/(2*(1+Us)))) == 0)
e1, e2, bH, cH = sp.symbols("e1 e2 bH cH", positive=True)
check("Lemma 2 Cauchy step: e1^2/bH + e2^2/cH - (e1+e2)^2/(bH+cH) = (cH e1 - bH e2)^2/(bH cH (bH+cH))",
      sp.simplify(e1**2/bH + e2**2/cH - (e1+e2)**2/(bH+cH) - (cH*e1 - bH*e2)**2/(bH*cH*(bH+cH))) == 0)
TH = (Tt, 1 + Tt)
charts = {
 "HH": (dict(b=(Tt, (1+Tt)*Lq), c=(R2(1), (1+Tt)*Lq), p=(2*Tt*U, (1+Tt)*(1+U)), q=(2*W, (1+Tt)*(1+W))), "HH"),
 "HC": (dict(b=(Tt, (1+Tt)*Lq), p=(2*Tt*U, (1+Tt)*(1+U)), q=(2*(1+W), 1+Tt), c=(2*(1+W), (1+Tt)*(2*Lq*(1+W)+W))), "HC"),
 "CH": (dict(c=(Tt, (1+Tt)*Lq), q=(2*Tt*U, (1+Tt)*(1+U)), p=(2*(1+W), 1+Tt), b=(2*(1+W), (1+Tt)*(2*Lq*(1+W)+W))), "HC"),
 "CC": (dict(p=(2*Tt*(1+U), 1+Tt), q=(2*(1+U+W), 1+Tt), b=(Tt, (1+Tt)*Lq), c=(2*(1+U)*(1+U+W), (1+Tt)*(2*Lq*(1+U)*(1+U+W)+W))), "CC"),
 "CC'": (dict(q=(2*Tt*(1+U), 1+Tt), p=(2*(1+U+W), 1+Tt), c=(Tt, (1+Tt)*Lq), b=(2*(1+U)*(1+U+W), (1+Tt)*(2*Lq*(1+U)*(1+U+W)+W))), "CC"),
}
def Jchart(kind, ch):
    if kind == "HH":
        pq = fa(ch["p"], ch["q"]); return fa((1 + Lq, R2(1)), (-pq[0], 4 * pq[1]))
    if kind == "HC": return fa((2 * Lq + 1, R2(2)), (Tt, 2 * (1 + Tt) * (1 + U)))
    return fa((Lq, R2(1)), (R2(1), 2 * (1 + U)))
idx = {"p": IP, "b": IB, "c": IC, "q": IQ}
for name, (ch, kind) in charts.items():
    K = [eng2.subst(x, {idx[k]: v for k, v in ch.items()}, dgall) for x in (k0, k1, k2)]
    n, d0 = Jchart(kind, ch)
    check(f"{name}: chart numerators/denominators and J = n/d0 have nonnegative coefficients",
          all(nonneg(x) for v in ch.values() for x in v) and nonneg(n) and nonneg(d0))
    P0 = d0**2 * K[0]; P1 = d0 * (2 * d0 * K[0] + n * K[1]); P2 = d0**2 * K[0] + n * d0 * K[1] + n**2 * K[2]
    NUM = P0 + P1 * Y + P2 * Y**2
    check(f"{name}: r=s=0 part of P0 + P1 Y + P2 Y^2 has nonnegative coefficients", nonneg(NUM.compose([(r, R2(0)), (s, R2(0))])))
    check(f"{name}: P0 nonnegative", nonneg(P0))
    if nonneg(NUM):
        check(f"{name}: whole numerator coefficientwise nonnegative", True); continue
    if kind == "CC":
        ok, nb_ = grouped_nonneg(P2, IL)
        check(f"{name}: P1 nonnegative, P2 nonnegative by grouping in L ({nb_} univariate certificates)", nonneg(P1) and ok)
        continue
    groups = {}
    for mon, co in NUM.terms():
        key = tuple(m if i not in (IU, IY) else 0 for i, m in enumerate(mon))
        groups.setdefault(key, {}); groups[key][(mon[IU], mon[IY])] = groups[key].get((mon[IU], mon[IY]), 0) + co
    bad = [(k, g) for k, g in groups.items() if any(v < 0 for v in g.values())]
    ok = True
    for k, g in bad:
        Ac = {u: v for (u, y), v in g.items() if y == 2}; Bc = {u: v for (u, y), v in g.items() if y == 1}; Cc = {u: v for (u, y), v in g.items() if y == 0}
        if not (univ_nonneg(Ac) and univ_nonneg(Cc)): ok = False; continue
        if univ_nonneg(Bc): continue
        X = sp.Symbol("x"); mk = lambda dct: sp.Poly({(u,): v for u, v in dct.items()} or {(0,): 0}, X)
        disc = 4 * mk(Ac) * mk(Cc) - mk(Bc)**2
        if not univ_nonneg({m[0]: v for m, v in disc.terms()}): ok = False
    check(f"{name}: {len(bad)} groups (all variables except U,Y) with a negative coefficient, each A,C >= 0 and (B >= 0 or 4AC-B^2 >= 0) on U >= 0", ok)
    check(f"{name}: every exceptional group contains a positive power of r or s", all(k[IR] > 0 or k[IS] > 0 for k, g in bad))

if NUMERIC:
    import numpy as np
    rng = np.random.default_rng(7)
    def fvec(combo): return [sum(float(cc) * float(L.F_VECS[jj][k]) for jj, cc in combo) for k in range(14)]
    U_ = lambda: float(np.exp(rng.uniform(-2, 2)))
    tests = [("old-H (45,1,34) j=1+1/e", lambda: (lambda e: [(45, e), (1, 1 + 1/e), (34, 1)])(U_())),
             ("chain (1;0,45;34) j=1/e+1/b", lambda: (lambda x, y: [(0, x), (1, 1/x + 1/y), (45, y), (34, 1)])(U_(), U_())),
             ("four-coefficient (6,0,5,34) bp(a-1)=b+p+a+1", lambda: (lambda pp, bb: [(6, pp), (0, bb), (5, (bb*pp + bb + pp + 1)/(bb*pp - 1)), (34, 1)])(U_() + 1.5, U_() + 1.0)),
             ("old-H (33,2,46)", lambda: (lambda e: [(33, e), (2, 1 + 1/e), (46, 1)])(U_())),
             ("old-H (5,2,46)", lambda: (lambda e: [(5, e), (2, 1 + 1/e), (46, 1)])(U_())),
             ("chain (2;4,33;46)", lambda: (lambda x, y: [(4, x), (2, 1/x + 1/y), (33, y), (46, 1)])(U_(), U_())),
             ("chain (2;5,6;46)", lambda: (lambda x, y: [(5, x), (2, 1/x + 1/y), (6, y), (46, 1)])(U_(), U_()))]
    for name, mk in tests:
        worst = min(L.min_phi_over_SO4(fvec(mk()), restarts=12, seed=sd)[0] for sd in range(6))
        info(f"NUMERIC cited bound at its boundary, {name}: worst min Phi (6 samples x 12 BFGS restarts) = {worst:.3e}")

print("SUMMARY:", "ALL EXACT CHECKS PASSED" if not FAILED else f"{len(FAILED)} FAILED: {FAILED}")
sys.exit(1 if FAILED else 0)
