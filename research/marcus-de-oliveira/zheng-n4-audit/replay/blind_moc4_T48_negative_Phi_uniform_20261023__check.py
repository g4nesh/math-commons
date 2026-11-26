"""Replay for the note 'Excluding a negative discriminant factor on a six-contact cone'
(blind_moc4_T48_negative_Phi_uniform_20261023, supplement).

Own code only; nothing from the supplement is imported or executed.  The degree-nine
polynomial G (1772 terms) is loaded AS DATA from the supplement file
proofs/manin_transfer_20261021/blind_moc4_T48_native_full_radial_20261023.json
(location from env MOC4_SUPP or argument 2) and every claim about it is re-derived here.
Argument 1: path of the audit library moc4lib.  Exit status nonzero iff an EXACT check fails.
"""
import os, sys, json, random, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "lib"))
SUPP = sys.argv[2] if len(sys.argv) > 2 else os.environ.get("MOC4_SUPP", "")
import moc4lib as L
import sympy as sp
import flint
from fractions import Fraction
import numpy as np

FAIL = []
def check(name, ok, level="EXACT"):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}")
    if not ok and level == "EXACT":
        FAIL.append(name)
t0 = time.time()
T = frozenset([0, 3, 5, 6, 9, 12])
FV = L.F_VALUES
G8 = [7, 5, 29, 39, 2, 40, 41, 0]
names = ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a']
def add(*js):
    return tuple(sum(FV[j][s] for j in js) for s in range(24))
check("T = (0,3,5,6,9,12) is paper T_49; generators {0,1,2,5,7,9,29,39,40,41}",
      T == L.CONTACT_SETS[49] and L.generators(T) == [0, 1, 2, 5, 7, 9, 29, 39, 40, 41])
check("f40+f1 = f0+f7+f39 and f40+f9 = f2+f7+f41", add(40, 1) == add(0, 7, 39) and add(40, 9) == add(2, 7, 41))
# remaining-scope transport (27,1,59,5,28,6) -> (7,5,29,39,40,41)
src, dst = (27, 1, 59, 5, 28, 6), (7, 5, 29, 39, 40, 41)
hits = []
for g in L.G_ELEMENTS:
    t = L.g_index_map(g)
    def img(a):
        v = [None]*24
        for s_ in range(24): v[t[s_]] = FV[a][s_]
        return tuple(v)
    if all(img(a) == FV[b] for a, b in zip(src, dst)):
        hits.append(g)
check(f"remaining scope: exactly one g maps the roles (27,1,59,5,28,6) of paper T_47 to (7,5,29,39,40,41) in order "
      f"(g = {[(L.perm_to_oneline(g[0]), L.perm_to_oneline(g[1]), g[2]) for g in hits]})",
      len(hits) == 1 and set(src) <= set(L.generators(L.CONTACT_SETS[47])))

# ---------------- cubic along the height line ---------------------------------------
sy = sp.symbols('q p j c e d b a')
q, p, j, c, e, d, b, a = sy
hh = sp.Symbol('h')
cv = [sum(sp.Rational(L.F_VECS[G8[k]][i].numerator, L.F_VECS[G8[k]][i].denominator)*sy[k] for k in range(8)) for i in range(14)]
w = list(cv); w[13] = w[13] + hh
Fw = sp.Poly(sp.expand(L.cubic_C(w)), hh)
C0 = Fw.coeff_monomial(1); Q1 = Fw.coeff_monomial(hh)
check("C(m + h kappa) is affine in h", Fw.degree() == 1)
Qclaim = (c + d)*(d + e + j) - p*(e + j)
check("Q = (c+d)(d+e+j) - p(e+j) exactly (paper normalization of the cubic)", sp.expand(Q1 - Qclaim) == 0)
GJ = None
if SUPP:
    GJ = json.load(open(os.path.join(SUPP, "proofs", "manin_transfer_20261021", "blind_moc4_T48_native_full_radial_20261023.json")))
    vs = [sp.Symbol(x) for x in GJ['variables'][:8]]
    def poly_from(terms):
        return sp.Add(*[sp.Rational(cf)*sp.Mul(*[vv**ex for vv, ex in zip(vs, mono[:8])]) for mono, cf in terms])
    check("data: variable order is (q,p,j,c,e,d,b,a,z) and z never occurs",
          GJ['variables'] == ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a', 'z']
          and all(m[8] == 0 for m, _ in GJ['radial_factors'][1]['terms']))
    check("our C(m) equals the companion C_terms (same normalization)", sp.expand(C0 - poly_from(GJ['C_terms'])) == 0)
PC = sp.Poly(sp.expand((e + j)*C0 - e*j*Q1), *sy)
coefs = PC.coeffs()
check(f"(e+j)C - ejQ has {len(coefs)} terms, all positive (claim: 42); contains d^3 e and d^3 j",
      len(coefs) == 42 and all(x > 0 for x in coefs) and PC.coeff_monomial(d**3*e) > 0 and PC.coeff_monomial(d**3*j) > 0)
print("  time", round(time.time() - t0, 1))

# ---------------- (2): det H(w) Q^2 = -4 G^2 at random rational points -----------------
ms, Cs = L.cubic_sympy()
Hsym = [[sp.diff(Cs, ms[i], ms[k]) for k in range(14)] for i in range(14)]
Hlin = [[sp.lambdify(ms, Hsym[i][k]) for k in range(14)] for i in range(14)]
def hess_at(wv):
    return [[Fraction(Hlin[i][k](*wv)) if not isinstance(Hsym[i][k], sp.Integer) else Fraction(int(Hsym[i][k]))
             for k in range(14)] for i in range(14)]
cv_f = sp.lambdify(sy, cv); C0_f = sp.lambdify(sy, C0); Q_f = sp.lambdify(sy, Q1)
if GJ is not None:
    Gpoly = poly_from(GJ['radial_factors'][1]['terms'])
    G_f = sp.lambdify(sy, Gpoly)
    rng = random.Random(7)
    oks = []
    for _ in range(12):
        pt = [Fraction(rng.randint(1, 40), rng.randint(1, 9)) for _ in range(8)]
        Qv = Q_f(*pt)
        if Qv == 0: continue
        hv = -C0_f(*pt)/Qv
        wv = list(cv_f(*pt)); wv[13] += hv
        Dv = L.det(hess_at([Fraction(x) for x in wv]))
        oks.append(Dv*Qv**2 == -4*G_f(*pt)**2)
    check(f"(2) det H(w) Q^2 = -4 G^2 at h = -C/Q on {len(oks)} random positive rational points (probabilistic identity test)",
          all(oks) and len(oks) >= 10, "EXACT-SAMPLED")
print("  time", round(time.time() - t0, 1))

# ---------------- flint: gcd, discriminant (3), conic (4) ------------------------------
if GJ is not None:
    ctx = flint.fmpq_mpoly_ctx.get(('q', 'p', 'j', 'c', 'e', 'd', 'b', 'a'), 'lex')
    gens = ctx.gens()
    fq, fp, fj, fc, fe, fd, fb, fa = gens
    Gf = ctx.from_dict({tuple(m[:8]): flint.fmpq(int(cf)) for m, cf in GJ['radial_factors'][1]['terms']})
    check(f"G has {len(Gf.to_dict())} terms and total degree 9", len(Gf.to_dict()) == 1772 and Gf.total_degree() == 9)
    Gj = Gf.derivative('j')
    check("gcd(G, G_j) = 1", Gf.gcd(Gj).total_degree() == 0)
    # coefficients in j
    coeffj = {}
    for mono, cf in Gf.to_dict().items():
        k = mono[2]; mm = list(mono); mm[2] = 0
        coeffj.setdefault(k, {})[tuple(mm)] = cf
    check("G has degree 2 in j", max(coeffj) == 2)
    Aj = ctx.from_dict(coeffj.get(2, {})); Bj = ctx.from_dict(coeffj.get(1, {})); Cj = ctx.from_dict(coeffj.get(0, {}))
    disc = Bj*Bj - 4*Aj*Cj
    V = fa + fd + fe; Tq = fq*V - (fa + fd)*(fd + fe)
    Phi = ((fq + fd)*V*fb + fp*Tq)**2 + 4*fd*fp**2*(fa + fd)*(fq - fe)*V
    Y = fb + fc + fd
    fac = Y**2*(fq + fd)**2*Phi
    H6, rem = divmod(disc, fac)
    check(f"(3) B^2 - 4AC0 = (b+c+d)^2 (q+d)^2 Phi H6 exactly; H6 has {len(H6.to_dict())} terms (claim 159) and no j",
          rem == 0 and len(H6.to_dict()) == 159 and all(m[2] == 0 for m, _ in H6.to_dict().items()))
    # conic structure: A quadratic in c; B = -(q+d) Y (L1 c + L0); C0 = -K Y^2
    def deg_in(P, idx):
        return max((m[idx] for m, _ in P.to_dict().items()), default=0)
    LB, r1 = divmod(-Bj, (fq + fd)*Y)
    Kc, r2 = divmod(-Cj, Y**2)
    ok_struct = deg_in(Aj, 3) <= 2 and r1 == 0 and deg_in(LB, 3) <= 1 and r2 == 0 and deg_in(Kc, 3) == 0
    check("G = (A2 c^2 + A1 c + A0) j^2 - (q+d) Y (L1 c + L0) j - K Y^2 with A_i, L_i, K free of c and j", ok_struct)
    def cpart(P, k):
        return ctx.from_dict({m: cf for m, cf in P.to_dict().items() if m[3] == k})
    def strip_c(P, k):
        return ctx.from_dict({tuple(mm[:3]) + (0,) + tuple(mm[4:]): cf for mm, cf in P.to_dict().items() if mm[3] == k})
    A2, A1, A0 = strip_c(Aj, 2), strip_c(Aj, 1), strip_c(Aj, 0)
    L1, L0 = strip_c(LB, 1), strip_c(LB, 0)
    Mm = [[2*A2, A1, -(fq + fd)*L1], [A1, 2*A0, -(fq + fd)*L0], [-(fq + fd)*L1, -(fq + fd)*L0, -2*Kc]]
    detM = (Mm[0][0]*(Mm[1][1]*Mm[2][2] - Mm[1][2]*Mm[2][1]) - Mm[0][1]*(Mm[1][0]*Mm[2][2] - Mm[1][2]*Mm[2][0])
            + Mm[0][2]*(Mm[1][0]*Mm[2][1] - Mm[1][1]*Mm[2][0]))
    Wf = fb*fq*V - fe*fp*(fa + fd)
    check("(4) det M = -2 (e-b) p W (q+d)^2 Phi^2", detM == -2*(fe - fb)*fp*Wf*(fq + fd)**2*Phi**2)
    # identity modulo W: substitute b = e p (a+d)/(q V)  (W is linear in b with coefficient qV)
    qs = sp.symbols('q p j c e d b a')
    Vs = qs[7] + qs[5] + qs[4]
    Phis = ((qs[0] + qs[5])*Vs*qs[6] + qs[1]*(qs[0]*Vs - (qs[7] + qs[5])*(qs[5] + qs[4])))**2 \
        + 4*qs[5]*qs[1]**2*(qs[7] + qs[5])*(qs[0] - qs[4])*Vs
    brace = (-qs[7]*qs[5]*qs[4] + qs[7]*qs[5]*qs[0] + qs[7]*qs[0]**2 - qs[5]**2*qs[4] + qs[5]**2*qs[0]
             + qs[5]*qs[0]**2 + qs[4]*qs[0]**2)
    lhs = qs[4]**2*(qs[7] + qs[5])**2*Phis - qs[6]**2*Vs**2*brace**2
    bsub = qs[4]*qs[1]*(qs[7] + qs[5])/(qs[0]*Vs)
    check("e^2 (a+d)^2 Phi = b^2 V^2 {...}^2 modulo W", sp.cancel(lhs.subs(qs[6], bsub)) == 0)
print("  time", round(time.time() - t0, 1))

# ---------------- gcd(G, P) = 1 via restriction to a line ------------------------------
keep = [i for i in range(14) if i not in (0, 5)]
check("H_kappa is supported on rows/columns {0,5} only after deleting (0,5): the (0,5)-deleted principal block is height-free",
      all(L.hess_kappa()[i][k] == 0 for i in keep for k in keep))
for seed_line in ((11, 1, 2) if GJ is not None else ()):
    rng = random.Random(seed_line)
    p0 = [Fraction(rng.randint(-9, 9)) for _ in range(8)]
    v0 = [Fraction(rng.randint(1, 9)) for _ in range(8)]
    tt = sp.Symbol('t')
    # top-degree part of G at v0 nonzero  =>  any common factor restricts to a nonconstant polynomial
    Gtop = sum(sp.Rational(int(cf))*sp.Mul(*[sp.Rational(v0[k].numerator, v0[k].denominator)**mono[k] for k in range(8)])
               for mono, cf in GJ['radial_factors'][1]['terms'] if sum(mono[:8]) == 9)
    def P_at(pt):
        wv = list(cv_f(*pt))
        Hm = hess_at([Fraction(x) for x in wv])
        return L.det([[Hm[i][k] for k in keep] for i in keep])
    xs = [Fraction(k) for k in range(13)]
    ys = [P_at([p0[i] + xk*v0[i] for i in range(8)]) for xk in xs]
    Pline = sp.expand(sp.interpolate(list(zip([sp.Rational(x.numerator, x.denominator) for x in xs],
                                                [sp.Rational(y_.numerator, y_.denominator) for y_ in ys])), tt))
    Gline = sp.expand(Gpoly.subs({sy[i]: sp.Rational(p0[i].numerator) + tt*sp.Rational(v0[i].numerator) for i in range(8)}))
    g_ = sp.gcd(Pline, Gline)
    okline = Gtop != 0 and sp.degree(g_, tt) <= 0 and Pline != 0
    print(f"  line seed {seed_line}: G_top(v) = {Gtop}, deg P|line = {sp.degree(Pline, tt)}, gcd of restrictions = {sp.factor(g_)}")
    if okline:
        break
check("gcd(G, P) = 1: on some rational line with G_top(v) != 0 the restrictions of G and P are coprime "
      "(a common factor would restrict to a nonconstant common factor)", okline)
print("  time", round(time.time() - t0, 1))

# ---------------- NUMERIC tests ----------------------------------------------------------
Hnum = [[sp.lambdify(ms, Hsym[i][k], 'numpy') for k in range(14)] for i in range(14)]
def hess_np(wv):
    return np.array([[float(Hnum[i][k](*wv)) for k in range(14)] for i in range(14)])
def minors12(Hm):
    out = []
    for i, k in itertools.combinations(range(14), 2):
        idx = [x for x in range(14) if x not in (i, k)]
        out.append(np.linalg.det(Hm[np.ix_(idx, idx)]))
    return np.array(out)
if GJ is not None:
    import mpmath
    Gs = [sp.lambdify(sy, sp.Poly(Gpoly, j).coeff_monomial(j**k)) for k in (2, 1, 0)]
    Phif = sp.lambdify(sy, Phis.subs(dict(zip(qs, sy))))
    H6f = None
    rng = np.random.default_rng(2026)
    n_pts = n_neg = n_rank12 = n_bad = n_le12 = 0; worst_ratio = 0.0
    for trial in range(3000):
        prm = rng.uniform(0.05, 3.0, 8)
        if trial % 2: prm[3] = 0.0 if trial % 4 == 1 else prm[3]   # include c = 0
        args = list(prm)
        cA, cB, cC = [g(*[args[i] if i != 2 else 0.0 for i in range(8)]) for g in Gs]
        roots = np.roots([cA, cB, cC]) if abs(cA) > 1e-14 else np.array([-cC/cB])
        for jr in roots:
            if abs(jr.imag) > 1e-12: continue
            pt = list(prm); pt[2] = jr.real
            Qv = Q_f(*pt)
            if abs(Qv) < 1e-8: continue
            n_pts += 1
            wv = np.array(cv_f(*pt), float); wv[13] += -C0_f(*pt)/Qv
            Hm = hess_np(wv)
            sv = np.linalg.svd(Hm, compute_uv=False)
            if sv[-2] < 1e-7*sv[0]: n_le12 += 1
            if sv[-2] < 1e-7*sv[0] and sv[-3] > 1e-5*sv[0]: n_rank12 += 1
            if Phif(*pt) < 0:
                n_neg += 1
                mins = minors12(Hm)
                scale = np.max(np.abs(mins))
                if np.max(mins) > 1e-7*scale: n_bad += 1
    print(f"  numeric (seed 2026): {n_pts} real points of G = 0 with Q != 0 (positive orthant except j), "
          f"{n_le12} have numerical rank <= 12, {n_rank12} exactly 12; {n_neg} have Phi < 0, of which {n_bad} have a positive 12-principal minor")
    check("NUMERIC: rank of H(w) is <= 12 at real points of G = 0 (claim)", n_pts > 0 and n_le12 == n_pts, "NUMERIC")
    check("NUMERIC: at real points of G = 0 with Phi < 0 every 12-principal minor is <= 0", n_bad == 0, "NUMERIC")
# cited global theorem: 12-principal minors of Hess C(Re x) >= 0 for x in X
rng = np.random.default_rng(99)
nbad = 0; ntot = 0; worst = 0.0
for trial in range(300):
    av = rng.normal(size=4) + 1j*rng.normal(size=4); bv = rng.normal(size=4) + 1j*rng.normal(size=4)
    al = np.exp(1j*rng.uniform(0, 2*np.pi))
    mvec = L.matching_vector_np(av, bv, al)
    uvec = mvec.real
    Hm = hess_np(uvec)
    mins = minors12(Hm)
    sc = np.max(np.abs(mins)); ntot += 1
    worst = min(worst, np.min(mins)/sc)
    if np.min(mins) < -1e-8*sc: nbad += 1
print(f"  numeric (seed 99): cited Hessian theorem, 300 random complex spectra and phases: {nbad} with a negative "
      f"12-principal minor of Hess C(Re x); worst relative value {worst:.3e}")
check("NUMERIC: cited global theorem (all 12-principal minors of Hess C(Re x) >= 0) not falsified", nbad == 0, "NUMERIC")
print("  total time", round(time.time() - t0, 1))
print("FAILED:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)
