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

# ===== root_moc4_T30_uniform_20261022 : paper case T_31 = {0,1,5,6,9,12} =====
T = L.CONTACT_SETS[31]
check("doc contact set {0,1,5,6,9,12} is the paper's T_31 literally", sorted(T) == [0, 1, 5, 6, 9, 12])
gens = L.generators(T)
b3 = [2, 29, 1, 9, 5, 35, 0, 7]; ev = [15, 4, 17, 7, 11, 3, 18, 14]
check("source list (1) = lib generators(T_31)", sorted(b3) == gens)
check("contact equations have rank 6 (kernel dim 8)", L.rank([L.P[s] for s in T]) == 6)
check("values of the 8 generators at permutations (15,4,17,7,11,3,18,14) form the identity matrix",
      all(V[b3[i]][ev[k]] == (1 if i == k else 0) for i in range(8) for k in range(8)))
check("only f35 is non-minor", [j for j in b3 if not L.F_IS_MINOR[j]] == [35])

# radial necessity
Jr, Cm, Qm, syms = radial_identity("root_moc4_T30_full_radial_20261022.json",
                                   ["q", "p", "j", "c", "e", "d", "b", "a"], b3)
pC = sp.Poly(Cm, *syms)
check("C has only negative coefficients and contains -d^3", all(c < 0 for c in pC.coeffs()) and pC.coeff_monomial(syms[5]**3) == -1)

# Lemma 1 algebra
a_, b_, t_, c_, e_ = sp.symbols("a b t c e", positive=True)
E = t_ * (a_ + 1) * (b_ + t_) / (a_ * b_ - t_)
check("Lemma 1: E'(t) formula", sp.simplify(sp.diff(E, t_) - (a_ + 1) * (a_ * b_**2 + 2 * a_ * b_ * t_ - t_**2) / (a_ * b_ - t_)**2) == 0)
check("Lemma 1: homogeneous four-coefficient bound t a' b (e - t) >= t^2 (a' + b + e + t), a'=ta, solves to e >= E(t)",
      [sp.simplify(s - E) for s in sp.solve(sp.Eq(t_ * a_ * b_ * (e_ - t_), t_**2 * (t_ * a_ + b_ + e_ + t_)), e_)] == [0])
check("Lemma 1: E(1) = (a+1)(b+1)/(ab-1)", sp.simplify(E.subs(t_, 1) - (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)) == 0)

# sign certificates on G9 (read from the JSON as data, identity-tested above)
R, a, b, c, e, q, p, j, A, Tt, Vv = ring("a,b,c,e,q,p,j,A,T,V", ZZ)
IA, IB, IC, IE = 0, 1, 2, 3
eng = Eng(R)
G = load_ring_poly(R, Jr["radial_factors"][1]["terms"], Jr["variables"],
                   {"q": q, "p": p, "j": j, "c": c, "e": e, "d": R(1), "b": b, "a": a, "z": R(1)})
kk = eng.coeffs_in(G, 6)
check("G9 has 1183 terms and degree 2 in j", len(G.terms()) == 1183 and max(kk) == 2)
k0, k1, k2 = kk[0], kk[1], kk[2]
check("k2 has nonnegative coefficients", nonneg(k2))
Bpol = (a + 1)**3 * (b + 1)**3
check("k0 at c=e=q=p=0 equals B=(a+1)^3(b+1)^3", k0.compose([(c, R(0)), (e, R(0)), (q, R(0)), (p, R(0))]) == Bpol)
check("k0 at e=0 minus B has nonnegative coefficients", nonneg(k0.compose(e, R(0)) - Bpol))
check("c=0: k1, k2 nonnegative coefficients", nonneg(k1.compose(c, R(0))) and nonneg(k2.compose(c, R(0))))
check("a=0: k0,k1,k2 nonnegative coefficients and k0 constant 1",
      all(nonneg(x.compose(a, R(0))) for x in (k0, k1, k2)) and k0.compose([(x, R(0)) for x in (a, b, c, e, q, p)]) == 1)
dgall = [max(x, y, z) for x, y, z in zip(eng.degs(k0), eng.degs(k1), eng.degs(k2))]
# (5) four charts
for sw in (0, 1):
    for regime in ("ab<1", "ab>1"):
        Rr = 1 + A
        zz = (Tt, 1 + Tt) if regime == "ab<1" else (1 + Tt, R(1))
        na, da, nb, db = zz[0] * Rr, zz[1], zz[0], zz[1] * Rr
        if sw: na, da, nb, db = nb, db, na, da
        if regime == "ab<1":
            ev_ = (Vv, R(1)); fac = R(1)
        else:
            ev_ = ((na + da) * (nb + db) * Vv, (na * nb - da * db) * (1 + Vv)); fac = 1 + Vv  # 1 - Le*e = 1/(1+V)
        dg = eng.degs(k0)
        K0 = eng.subst(k0, {IA: (na, da), IB: (nb, db), IE: ev_}, dg)
        Dn = da**dg[IA] * db**dg[IB] * ev_[1]**dg[IE]
        expr = K0 * da**3 * db**3 * fac**2 - (na + da)**3 * (nb + db)**3 * Dn
        check(f"(5) k0 >= B(1-Le e)^2, chart swap={sw} {regime}: {len(expr.terms())} terms, nonnegative (denominators positive: {nonneg(Dn)})",
              nonneg(expr) and nonneg(Dn) and nonneg(ev_[1]))
# seven domains
def domain(d):
    one = R(1)
    if d == 1: fa_, fb_ = (A, 1 + A), (A * Tt, (1 + A) * (1 + Tt))
    elif d == 2: fa_, fb_ = (A * Tt, (1 + A) * (1 + Tt)), (A, 1 + A)
    elif d == 3: fa_, fb_ = (1 + A, one), (Tt, (1 + Tt) * (1 + A))
    elif d == 4: fa_, fb_ = (Tt, (1 + Tt) * (1 + A)), (1 + A, one)
    elif d == 5: fa_, fb_ = (A, 1 + A), ((1 + Tt) * (1 + A), A)
    elif d == 6: fa_, fb_ = (1 + A + Tt, 1 + A), ((1 + Tt)**2 * (1 + A), 1 + A + Tt)
    else: fa_, fb_ = ((1 + Tt) * (1 + A), one), (1 + Tt, 1 + A)
    ab = fm(fa_, fb_); Jold = (c + 1, c)
    Em = ((fa_[0] + fa_[1]) * (fb_[0] + fb_[1]) * ab[1], fa_[1] * fb_[1] * (ab[0] - ab[1]))
    if d in (1, 2, 4): fe_, Jf = (Vv, one), Jold
    elif d == 5: fe_, Jf = (Em[0] * Vv, Em[1] * (1 + Vv)), Jold
    elif d == 6: fe_, Jf = (Em[0] * Vv, Em[1] * (1 + Vv)), (fa_[1] * c + fa_[0], fa_[0] * c)
    else:
        ft = (ab[0] * Vv, ab[1] * (1 + Vv)) if d == 3 else (Vv, 1 + Vv)
        a1 = (fa_[0] + fa_[1], fa_[1]); bt = fa(fb_, ft); abmt = fa(ab, (-ft[0], ft[1]))
        num = fm(fm(ft, a1), bt); fe_ = (num[0] * abmt[1], num[1] * abmt[0])
        omt = (ft[1] - ft[0], ft[1])
        Jf = fa((omt[0] * fa_[1], omt[1] * fa_[0]), (omt[0]**2, omt[1]**2 * c))
    return fa_, fb_, fe_, Jf
t0 = time.time()
for d in range(1, 8):
    fa_, fb_, fe_, (N, Dd) = domain(d)
    check(f"domain {d}: chart numerators and denominators of a, b, e, J have nonnegative coefficients",
          all(nonneg(x) for x in (fa_[0], fa_[1], fb_[0], fb_[1], fe_[0], fe_[1], N, Dd)))
    K0, K1, K2 = [eng.subst(k, {IA: fa_, IB: fb_, IE: fe_}, dgall) for k in (k0, k1, k2)]
    if d == 1:
        ln, ld = fa(fa_, (-fb_[0], fb_[1]))
        check("domain 1: cut l = a - b has nonnegative numerator/denominator", nonneg(ln) and nonneg(ld))
        Gl = K0 * ld**2 + K1 * ln * ld + K2 * ln**2
        check("domain 1: [0,l] middle 2k0 + l k1 >= 0", nonneg(2 * K0 * ld + ln * K1))
        Dn = fa_[1]**dgall[IA] * fb_[1]**dgall[IB] * fe_[1]**dgall[IE]
        check("domain 1: (7) G9(l) - B >= 0", nonneg(Gl * fa_[1]**3 * fb_[1]**3 - (fa_[0] + fa_[1])**3 * (fb_[0] + fb_[1])**3 * Dn * ld**2))
        JL = fa((N, Dd), (-ln, ld))
        check("domain 1: J - l has nonnegative numerator/denominator", nonneg(JL[0]) and nonneg(JL[1]))
        mid = 2 * Gl * JL[1] * ld + JL[0] * (K1 * ld + 2 * K2 * ln) * ld**2
        check("domain 1: [l,J] middle 2G(l) + (J-l)G'(l) >= 0", nonneg(mid))
        check("domain 1: endpoint G9(J) >= 0", nonneg(Dd**2 * K0 + N * Dd * K1 + N**2 * K2))
    else:
        check(f"domain {d}: middle 2D k0 + N k1 >= 0", nonneg(2 * Dd * K0 + N * K1))
        EP = Dd**2 * K0 + N * Dd * K1 + N**2 * K2
        if nonneg(EP):
            check(f"domain {d}: endpoint D^2 k0 + N D k1 + N^2 k2 >= 0 (coefficientwise)", True)
        else:
            ok, nb_ = grouped_nonneg(EP, IC)
            info(f"domain {d}: endpoint not coefficientwise nonnegative; {nb_} groups (all variables except c) need a univariate certificate")
            check(f"domain {d}: endpoint >= 0 by grouping in c (exact squarefree + real-root isolation)", ok and d == 3)
info("certificate time", round(time.time() - t0, 1), "s")

# orbit witnesses for the three cited bounds (base profile tuple -> target tuple, by some g in G)
_vidx = {vv: jj for jj, vv in enumerate(V)}
def _gf(g, jj):
    mp = L.g_index_map(g); out = [0] * 24
    for s_ in range(24): out[mp[s_]] = V[jj][s_]
    return _vidx.get(tuple(out))
for base, targets in (((0, 6, 10, 54), [(7, 0, 5, 35), (0, 7, 5, 35)]), ((3, 5, 54), [(9, 1, 35)]), ((2, 0, 12, 70), [(1, 7, 9, 35), (1, 9, 7, 35)])):
    check(f"some g in G maps base profiles {base} to target {targets[0]} (or the stated exchange)",
          any(tuple(_gf(g, jj) for jj in base) in targets for g in L.G_ELEMENTS))

if NUMERIC:
    import numpy as np
    rng = np.random.default_rng(7)
    def vec(combo): return [sum(float(cc) * float(L.F_VECS[jj][k]) for jj, cc in combo) for k in range(14)]
    for name, mk in (("old-H j=1+1/c: c f9 + j f1 + f35", lambda u: [(9, u), (1, 1 + 1 / u), (35, 1)]),
                     ("chain j=1/a+1/c: a f7 + j f1 + c f9 + f35", lambda u: [(7, u), (1, 1 / u + 1 / (2 * u + .3)), (9, 2 * u + .3), (35, 1)]),
                     ("four-coefficient ab(e-1)=a+b+e+1: a f7 + b f0 + e f5 + f35", lambda u: [(7, u), (0, 2 / u + .5), (5, (u * (2 / u + .5) + u + 2 / u + .5 + 1) / (u * (2 / u + .5) - 1)), (35, 1)])):
        worst = min(L.min_phi_over_SO4(vec(mk(float(np.exp(rng.uniform(-2, 2))))), restarts=12, seed=s)[0] for s in range(6))
        info(f"NUMERIC cited bound at its boundary, {name}: worst min Phi over 6 samples x 12 BFGS restarts = {worst:.3e}")

print("SUMMARY:", "ALL EXACT CHECKS PASSED" if not FAILED else f"{len(FAILED)} FAILED: {FAILED}")
sys.exit(1 if FAILED else 0)
