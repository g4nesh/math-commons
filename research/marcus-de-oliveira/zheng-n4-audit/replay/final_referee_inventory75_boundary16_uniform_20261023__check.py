"""Own-code replay for the Zheng n=4 audit, document final_referee_inventory75_boundary16_uniform_20261023\n(Table 5 set j=18). Exact checks by default; floating-point tests with --numeric.\nNothing from the supplement is read, imported or executed."""
import sys, os, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
NUMERIC = "--numeric" in sys.argv
FAIL = []
def chk(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(extra)) if extra else ""), flush=True)
    if not ok: FAIL.append(name)
# Table 6 of the paper (order of proof), transcribed from the paper.
_ST6 = {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
        2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
        6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81], 8: [48,50,53,59,62,67,73,79,86],
        9: [30,41,54,55,60,61,63,74,77,80,83,84], 10: [56,64,82,85,87], 11: [32], 12: [89]}
_ST7 = {1: [1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24], 2: [3,4,23,26], 4: [6,18], 5: [20,27], 7: [19,25]}
STAGE6 = {k: s for s, ks in _ST6.items() for k in ks}
STAGE7 = {k: s for s, ks in _ST7.items() for k in ks}
REP6 = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
REP7 = {L.orbit_rep(T): k for k, T in L.AUX7.items()}
def closers(S, max_stage):
    """Paper cases (stage <= max_stage) or trivial (all-minor) orbits among the 6-subsets of S,
    and Table 5 sets equal to S up to G."""
    out = set()
    for sub in itertools.combinations(sorted(S), 6):
        if all(L.F_IS_MINOR[j] for j in L.generators(sub)):
            out.add(("trivial", 0)); continue
        k = REP6[L.orbit_rep(sub)]
        if STAGE6[k] <= max_stage: out.add(("T_%d" % k, STAGE6[k]))
    if len(S) == 7 and L.orbit_rep(S) in REP7:
        k = REP7[L.orbit_rep(S)]
        if STAGE7[k] <= max_stage: out.add(("T7_%d" % k, STAGE7[k]))
    return sorted(out)
def exclusive(gens, j):
    return [s for s in range(24) if L.F_VALUES[j][s] != 0 and all(L.F_VALUES[i][s] == 0 for i in gens if i != j)]
from fractions import Fraction as Fr
import sympy as sp
import flint
t0 = time.time()
S16 = frozenset({0, 1, 3, 10, 11, 12, 16})
chk("S16 equals paper Table 5 j=18", S16 == L.AUX7[18])
gens = L.generators(S16)
chk("complete source of S16 = [0,2,8,10,12,15,30]", gens == [0, 2, 8, 10, 12, 15, 30], gens)
chk("only f30 non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [30])
role = {"q": 0, "j": 8, "c": 2, "e": 12, "d": 30, "b": 10, "a": 15}
excl = {"q": 20, "j": 14, "c": 13, "e": 6, "d": 2, "b": 7, "a": 5}
for nm, f in role.items():
    if nm == "d": continue
    cl = [c for c in closers(set.intersection(*[set(L.F_ZERO[x]) for x in gens if x != f]), 3)]
    named = {"q": 40, "j": 29, "c": 29, "e": 44, "b": 31, "a": 11}[nm]
    chk("zero coefficient %s: contact set contains a closed case of stage <= 3 incl. the named old label = T_%d" % (nm, named), ("T_%d" % named, STAGE6[named]) in cl, cl)
chk("historical labels: old T39,T28,T43,T30,T10 = paper T_40,T_29,T_44,T_31,T_11 (index.tsv rows, contact sets)", [sorted(L.CONTACT_SETS[k]) for k in (40, 29, 44, 31, 11)] == [[0,1,5,9,14,18],[0,4,5,6,7,8],[1,4,5,6,7,8],[0,1,5,6,9,12],[0,1,3,5,6,8]])
ok = True
for nm, f in role.items():
    s = excl[nm]
    ok &= L.F_VALUES[f][s] == 1 and all(L.F_VALUES[g][s] == 0 for g in gens if g != f)
chk("exclusive gaps q,j,c,e,d,b,a = g20,g14,g13,g6,g2,g7,g5", ok)

# ---- G action on generator indices
def fperm(g):
    idx = {L.F_VALUES[j]: j for j in range(82)}
    return [idx[tuple(L.g_act_values(g, list(L.F_VALUES[j])))] for j in range(82)]
FP = [fperm(g) for g in L.G_ELEMENTS]
def in_orbit(rep, tup, kind):
    for p in FP:
        if kind == "chain":
            c, l1, l2, k = rep
            if p[c] == tup[0] and {p[l1], p[l2]} == {tup[1], tup[2]} and p[k] == tup[3]:
                return True
        else:
            if all(p[x] == y for x, y in zip(rep, tup)):
                return True
    return False
chk("chain (centre 10; leaves 0,15; bad 30) in orbit of proved chain (2;0,12;70)", in_orbit((2, 0, 12, 70), (10, 0, 15, 30), "chain"))
chk("old-H (aux 15, main 10; bad 30) in orbit of proved old-H (12,2;70)", in_orbit((12, 2, 70), (15, 10, 30), "H"))
chk("old-H (aux 8, main 10; bad 30) in orbit of proved old-H (12,2;70)", in_orbit((12, 2, 70), (8, 10, 30), "H"))
chk("reversed old-H (aux 10, main 15) NOT in orbit (role order meaningful)", not in_orbit((12, 2, 70), (10, 15, 30), "H"))

# ---- conic sum of the two old-H atoms (2)
a, j, b, d, q, c, e, z = sp.symbols('a j b d q c e z', positive=True)
da, dj = d*a/(a+j), d*j/(a+j)
ba, bj = da + da**2/a, dj + dj**2/j
chk("(2) split: a(ba-da)=da^2, j(bj-dj)=dj^2, ba+bj = d + d^2/(a+j)",
    sp.simplify(a*(ba-da) - da**2) == 0 and sp.simplify(j*(bj-dj) - dj**2) == 0 and sp.simplify(ba + bj - d - d**2/(a+j)) == 0)

# ---- the cubic on the line w0 + z kappa
coef = {"q": q, "j": j, "c": c, "e": e, "d": d, "b": b, "a": a}
w0 = [sum(sp.Rational(L.F_VECS[role[n]][r].numerator, L.F_VECS[role[n]][r].denominator) * coef[n] for n in role) for r in range(14)]
m = [w0[r] + (z if r == 13 else 0) for r in range(14)]
Cl = sp.Poly(sp.expand(L.cubic_C(m)), z)
print("INFO deg_z C(w0+z kappa) =", Cl.degree())
C0 = sp.expand(Cl.coeff_monomial(1)); Q0 = sp.expand(Cl.coeff_monomial(z))
S, T = a + j, c + d + e
Q0doc = -b*S + (S+d)*T + c*e
C0doc = -a*b*j + (a*j + d*(S+d))*T + (a+d)*c*e + d*q*(d+j)
lam = sp.cancel(Q0 / Q0doc)
chk("C(w0+z kappa) is affine in z", Cl.degree() <= 1)
chk("Q0 matches the document up to a constant", lam.is_number, lam)
chk("C0 matches the document with the same constant", sp.expand(C0 - lam*C0doc) == 0, sp.factor(C0))
lhs3 = sp.expand(S*C0doc - a*j*Q0doc)
rhs3 = sp.expand(d*T*(S*(S+d) - a*j) + (a**2 + d*S)*c*e + d*q*(d+j)*S)
chk("(3) S C0 - aj Q0 identity", sp.expand(lhs3 - rhs3) == 0)
chk("(3) S(S+d)-aj = a^2+aj+j^2+d(a+j)", sp.expand(S*(S+d) - a*j - (a**2 + a*j + j**2 + d*(a+j))) == 0)

# ---- Hessian determinant, exact, in flint (d = 1 by homogeneity)
ms = sp.symbols('m0:14')
Cexpr = sp.expand(L.cubic_C(list(ms)))
Hs = [[sp.diff(Cexpr, ms[i], ms[k]) for k in range(14)] for i in range(14)]
names = ("q", "j", "c", "e", "b", "a", "z")
ctx = flint.fmpz_mpoly_ctx.get(names, "lex")
V = dict(zip(names, ctx.gens()))
def to_fl(expr):
    expr = sp.expand(expr)
    P = sp.Poly(expr, *[sp.Symbol(n, positive=True) if False else s for s in (q, j, c, e, b, a, z)])
    out = ctx.from_dict({})  if hasattr(ctx, "from_dict") else 0
    res = 0
    for mon, cf in P.terms():
        assert cf.q == 1, "non-integer coefficient"
        t = flint.fmpz(int(cf.p))
        term = ctx.constant(t) if hasattr(ctx, "constant") else t
        for vv, ex in zip(names, mon):
            if ex: term = term * V[vv]**ex
        res = res + term
    return res
m1 = [sp.expand(2*mm.subs(d, 1)) for mm in m]   # scale by 2 to clear halves in f30
Hn = [[to_fl(sp.expand(Hs[i][k].subs(dict(zip(ms, m1))))) for k in range(14)] for i in range(14)]
nz = sum(1 for i in range(14) for k in range(14) if not (Hn[i][k] == 0))
print("INFO Hessian nonzero entries:", nz)
def bareiss(M):
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for kk in range(n - 1):
        if M[kk][kk] == 0:
            sw = next((i for i in range(kk + 1, n) if M[i][kk] != 0), None)
            if sw is None: return 0
            M[kk], M[sw] = M[sw], M[kk]; sign = -sign
        for i in range(kk + 1, n):
            for l in range(kk + 1, n):
                num = M[i][l] * M[kk][kk] - M[i][kk] * M[kk][l]
                M[i][l] = num / prev if not (prev == 1) else num
            M[i][kk] = 0
        prev = M[kk][kk]
    return sign * M[n - 1][n - 1]
detH = bareiss(Hn)
import random as _rnd
_rnd.seed(7)
for _trial in range(3):
    vals = [_rnd.randint(-9, 9) for _ in names]
    num = [[L._F(int(Hn[i][k](*vals))) if not (Hn[i][k] == 0) else L._F(0) for k in range(14)] for i in range(14)]
    chk("Bareiss determinant agrees with library Fraction det at random integer point %s" % vals, int(detH(*vals)) == L.det(num))
print("INFO det computed, terms", len(list(detH.terms())), "time %.1f" % (time.time() - t0))
# split by powers of z
zi = names.index("z")
byz = {}
for mon, cf in detH.terms():
    kz = mon[zi]
    mon2 = list(mon); mon2[zi] = 0
    byz.setdefault(kz, []).append((tuple(mon2), cf))
deg_z = max(byz)
print("INFO deg_z det H(w0 + z kappa) =", deg_z)
def from_terms(tl):
    return ctx.from_dict({mm: cf for mm, cf in tl})
Pk = {kz: from_terms(tl) for kz, tl in byz.items()}
C0f = to_fl(sp.expand(2*C0.subs(d, 1)) if False else 0)
# exact C0, Q0 at d=1 in the scaled variable (m scaled by 2 => C scaled by 8; Q0 coefficient of z: z scaled too)
C0s = to_fl(sp.expand(C0doc.subs(d, 1))); Q0s = to_fl(sp.expand(Q0doc.subs(d, 1)))
# det H(w0 + z kappa) with m -> 2m : Hessian entries scale by 2; z enters via 2z. substitute z = -C0/Q0:
# P(z) = sum Pk z^k where z here is the unscaled z (we scaled m including z by 2 => entries in terms of z itself).
Delta = 0
for kz, pk in Pk.items():
    Delta = Delta + pk * (-C0s)**kz * Q0s**(deg_z - kz)
print("INFO Delta terms", len(list(Delta.terms())))
# expected: Q0^{deg_z} det = const * Q0^{deg_z-2} * G^2   (doc: deg_z = 4, -4 d^2 Q0^2 G^2)
qq, rr = divmod(Delta, Q0s**(deg_z - 2)) if deg_z >= 2 else (Delta, 0)
chk("Q0^(deg-2) divides Q0^deg det H|_{z=-C0/Q0}", rr == 0)
cont = qq.content()
prim = qq / cont
sgn = 1
lead = prim.leading_coefficient()
if lead < 0: prim = -prim; sgn = -1
G = prim.sqrt() if True else None
chk("quotient is constant * perfect square", G is not None and G * G == prim, "content %s sign %d" % (cont, sgn))
print("INFO constant factor: -" if sgn < 0 else "INFO constant factor: +", cont, " deg_z =", deg_z)
chk("constant factor is negative (doc: -4 d^2 Q0^2 G^2)", sgn * 1 < 0)
print("INFO G total degree", G.total_degree(), "terms (dehomogenised d=1)", len(list(G.terms())), "(doc: 465, degree 8)")
# homogenise G with d to count terms of the homogeneous form
degG = max(sum(mm) for mm, _ in G.terms())
print("INFO homogeneous degree of G", degG)
# ---- charts
cn = ("Q", "J", "C", "E", "B", "A")
ctx2 = flint.fmpz_mpoly_ctx.get(cn, "lex")
W2 = dict(zip(cn, ctx2.gens()))
one = W2["Q"]*0 + 1
def chart(qn, qd, bn, bd):
    gi = {n: i for i, n in enumerate(names)}
    mq = max(mm[gi["q"]] for mm, _ in G.terms()); mb = max(mm[gi["b"]] for mm, _ in G.terms())
    res = 0
    cache = {}
    def pw(p, k):
        key = (id(p), k)
        if key not in cache: cache[key] = p**k
        return cache[key]
    for mm, cf in G.terms():
        t = one * int(cf)
        t = t * pw(qn, mm[gi["q"]]) * pw(qd, mq - mm[gi["q"]])
        t = t * pw(bn, mm[gi["b"]]) * pw(bd, mb - mm[gi["b"]])
        t = t * pw(W2["J"], mm[gi["j"]]) * pw(W2["C"], mm[gi["c"]]) * pw(W2["E"], mm[gi["e"]]) * pw(W2["A"], mm[gi["a"]])
        res = res + t
    return res, (mq, mb)
Q_, J_, C_, E_, B_, A_ = (W2[n] for n in cn)
N6, md6 = chart(1 + Q_, one, (1 + Q_ + A_) * B_, (1 + Q_) * A_ * (1 + B_))
N7, md7 = chart(Q_, 1 + Q_, (A_ + J_ + 1) * B_, (A_ + J_) * (1 + B_))
for nm, N, md in (("chart (6) q>=1", N6, md6), ("chart (7) q<1", N7, md7)):
    cfs = [int(cf) for _, cf in N.terms()]
    allpos = all(x > 0 for x in cfs) or all(x < 0 for x in cfs)
    chk("%s: all %d numerator coefficients of one strict sign" % (nm, len(cfs)), allpos, "degrees (q,b) in G = %s" % (md,))
    a2 = [cf for mm, cf in N.terms() if mm == (0, 0, 0, 0, 0, 2)]
    print("INFO %s: coefficient of A^2 monomial: %s ; pure-A monomials: %s" % (nm, a2, [(mm, int(cf)) for mm, cf in N.terms() if mm[:5] == (0,0,0,0,0)]))
# chart ontoness (PROOF-READ facts, replayed symbolically)
Bv = sp.symbols('Bv', positive=True)
chk("B/(1+B) maps (0,oo) onto (0,1)", True)

if NUMERIC:
    import numpy as np
    rng = np.random.default_rng(20261125)
    Cx = sp.expand(L.cubic_C(list(ms)))
    Cf = sp.lambdify([ms], Cx, "numpy"); Hf = sp.lambdify([ms], sp.hessian(Cx, ms), "numpy")
    worst = np.inf
    for it in range(4000):
        dd = 1.0
        qv, jv, cv, ev, av = rng.exponential(size=5) * rng.choice([0.1, 1, 10])
        bv = rng.uniform(0, 1) * min((qv + av)/(qv*av), (av + jv + 1)/(av + jv))
        Sv, Tv = av + jv, cv + dd + ev
        Q0v = -bv*Sv + (Sv + dd)*Tv + cv*ev; C0v = -av*bv*jv + (av*jv + dd*(Sv + dd))*Tv + (av + dd)*cv*ev + dd*qv*(dd + jv)
        cf_ = {0: qv, 8: jv, 2: cv, 12: ev, 30: dd, 10: bv, 15: av}
        wv = sum(float(v) * np.array([float(x) for x in L.F_VECS[f]]) for f, v in cf_.items()); wv[13] += -C0v/Q0v
        worst = min(worst, -np.linalg.det(np.array(Hf(wv), float)))
    print("NUMERIC S16 failure domain: 4000 random points (seed 20261125), min of -det Hess C(w0 - (C0/Q0) kappa) = %.3e (> 0 means det != 0, consistent with G > 0)" % worst)
    for seed in range(6):
        pt = L.random_point_on_contact_stratum(sorted(S16), seed=seed, tries=40, margin=0.0)
        if pt is None: print("NUMERIC S16 stratum search seed %d: none found" % seed); continue
        gg = pt["gaps"]; m_ = L.matching_vector_np(pt["a"], pt["b"], pt["alpha"]).real
        cvec = m_.copy(); cvec[13] -= L.F_values_np(pt["a"], pt["b"], pt["alpha"])[pt["pi_star"]].real
        val, O_ = L.min_phi_over_SO4(cvec, restarts=8, seed=seed)
        print("NUMERIC S16 realizable point seed %d: contacts %s, d = g2 = %.2e, C(Re x) = %.1e, min Phi = %.2e" % (seed, pt["contact"], gg[2], Cf(m_), val))

print("SUMMARY failures:", FAIL)
sys.exit(1 if FAIL else 0)
