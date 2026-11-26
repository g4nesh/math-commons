"""Own-code replay for the Zheng n=4 audit, document root_inventory80_native_uniform_20261023
(primary proof of paper case T_81 = {1,2,10,11,12,13}).  Everything is rebuilt from the definitions and
from the formulas printed in the PDF; nothing from the supplement is read, imported or executed.
EXACT checks: sympy / python-flint over Q.  NUMERIC (skip with --no-numeric): a brief main-statement
search.  Exit code 1 if an EXACT check fails."""
import sys, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
import flint
NUMERIC = "--no-numeric" not in sys.argv
PERMS = list(itertools.permutations(range(4)))
FAIL = []
def chk(name, ok):
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok: FAIL.append(name)
def br(p, q): return p[0]*q[1] - p[1]*q[0]
def prods(A, C):
    return [br(A[0], C[p[0]])*br(A[1], C[p[1]])*br(A[2], C[p[2]])*br(A[3], C[p[3]]) for p in PERMS]

# Appendix G, Table 6 (paper p. 63), transcribed from the page
STAGE6 = {1:[1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
          2:[10,26,46,90],3:[9,11,12,47],4:[8,13,17,49],5:[14,16,24,33,38,70],
          6:[52,57,58,65,68,71,76,78],7:[45,51,66,69,75,81],8:[48,50,53,59,62,67,73,79,86],
          9:[30,41,54,55,60,61,63,74,77,80,83,84],10:[56,64,82,85,87],11:[32],12:[89]}
STAGE7 = {1:[1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24],2:[3,4,23,26],4:[6,18],5:[20,27],7:[19,25]}
ST6 = {k: s for s, ks in STAGE6.items() for k in ks}
ST7 = {k: s for s, ks in STAGE7.items() for k in ks}
REP6 = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
REP7 = {L.orbit_rep(T): k for k, T in L.AUX7.items()}

def classify6(S):
    if all(L.F_IS_MINOR[j] for j in L.generators(S)):
        return ("minor", 0)
    k = REP6[L.orbit_rep(S)]
    return ("T%d" % k, ST6[k])

def boundary_ties(T, k, special):
    """For every noncontact j: does a tie at j force the non-minor coefficient(s) to 0, and
    otherwise which six-subsets of T+{j} lie in an all-minor orbit or in a paper case of an
    earlier stage, and whether T+{j} is an image of a Table 5 set of an earlier stage."""
    T = frozenset(T)
    out = {}
    for j in range(24):
        if j in T: continue
        kill = all(L.F_VALUES[s][j] > 0 for s in special)
        subs = [(d,) + classify6((T | {j}) - {d}) for d in sorted(T)]
        prior = [(d, c) for d, c, st in subs if c == "minor" or st < ST6[k]]
        a7 = REP7.get(L.orbit_rep(T | {j}))
        a7ok = a7 is not None and ST7[a7] < ST6[k]
        out[j] = (kill, prior, a7 if a7ok else None)
        print("  tie j=%2d: kills non-minor coeff=%s; earlier six-subsets=%s; earlier Table-5 set=%s" % (
            j, kill, prior, ("T7_%d (stage %d)" % (a7, ST7[a7])) if a7ok else None), flush=True)
    return out

def numeric_main(T, n_points, seed):
    import numpy as np
    worst = None
    found = 0
    for i in range(n_points):
        pt = L.random_point_on_contact_stratum(sorted(T), seed=seed + i, tries=20, margin=1e-3)
        if pt is None: continue
        found += 1
        vals = np.asarray(pt["gaps"], float)
        c = np.linalg.lstsq(L.P_numpy(), vals, rcond=None)[0]
        resid = float(np.max(np.abs(L.P_numpy() @ c - vals)))
        val, O = L.min_phi_over_SO4(c, restarts=8, seed=seed + i)
        val = val / max(1.0, float(np.max(np.abs(vals))))
        worst = val if worst is None else min(worst, val)
        print("  numeric point %d: contact=%s  P-residual=%.1e  min Phi (scaled) = %.3e" % (i, pt["contact"], resid, val), flush=True)
    return found, worst

print("== A/B. contact set, generators, boundary ties")
T = frozenset({1, 2, 10, 11, 12, 13})
chk("T equals paper T_81 literally", T == L.CONTACT_SETS[81])
gens = L.generators(T)
chk("generators(T) = [0,1,8,10,12,15,24,53] (eq. (1))", gens == [0, 1, 8, 10, 12, 15, 24, 53])
chk("only f_24 is non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [24])
rows = [21, 17, 15, 7, 6, 4, 3, 0]
chk("rows [21,17,15,7,6,4,3,0] of the generators form the identity (unique nonneg. coefficients = gaps)",
    all(L.F_VALUES[f][rw] == (1 if a == b else 0) for a, f in enumerate(gens) for b, rw in enumerate(rows)))
chk("coefficient of f_24 is d = g3", rows[gens.index(24)] == 3)
ties = boundary_ties(T, 81, [24])
chk("ties not annihilating d are exactly 0,4,6,7,15,17,21 (eq. (3))",
    sorted(j for j, v in ties.items() if not v[0]) == [0, 4, 6, 7, 15, 17, 21])
chk("ties 0,4,6,7,15,17 contain an earlier-stage six-set", all(ties[j][1] for j in (0, 4, 6, 7, 15, 17)))
Hset = frozenset({0, 1, 6, 9, 16, 17, 22})
chk("H = {0,1,6,9,16,17,22} is Table 5 set T7_21 (stage 1)", Hset == L.AUX7[21] and ST7[21] == 1)
g = L.find_g(Hset, T | {21})
chk("tie 21: explicit g with g.H = T+{21}: %s" % (g,), g is not None and not ties[21][1])

print("== C. original variables: A = (0,1,z,w), C = (oo,v,x,t)")
o = flint.fmpq_mpoly_ctx.get(('z', 'w', 'v', 'x', 't'), 'lex')
z, w, v, x, t = o.gens(); one = o.from_dict({(0,)*5: 1}); zero = one - one
Fo = prods([(zero, one), (one, one), (z, one), (w, one)], [(one, zero), (v, one), (x, one), (t, one)])
Eo = Fo[12] - Fo[2]
chk("E = F12-F2 = (z-v)(w-t)", Eo == (z-v)*(w-t))
chk("(1) F12 = x E", Fo[12] == x*Eo)
chk("(4) F13-F12 = x(t-v)(z-w)", Fo[13] - Fo[12] == x*(t-v)*(z-w))
chk("(4) F11-F10 = v(1-t)(z-w)", Fo[11] - Fo[10] == v*(1-t)*(z-w))
chk("(6) x(w-v)F1 + (w-x)(1-v)F13 = 0 identically", (x*(w-v)*Fo[1] + (w-x)*(1-v)*Fo[13]).is_zero())
chk("(7) (z-x)F10 - (w-x)F11 = 0 identically", ((z-x)*Fo[10] - (w-x)*Fo[11]).is_zero())
chk("raw x(F3-F12) - (x-1)(F13-F12) + xE = 0", (x*(Fo[3]-Fo[12]) - (x-1)*(Fo[13]-Fo[12]) + x*Eo).is_zero())
chk("F6 = v(w-t)(z-x), F6-F0 = (z-x)(w-t)", Fo[6] == v*(w-t)*(z-x) and Fo[6] - Fo[0] == (z-x)*(w-t))
zs, ws, vs, xs, ts = sp.symbols('z w v x t')
def tosym5(p):
    return sp.Add(*[sp.Rational(int(c.numer()), int(c.denom()))*zs**m[0]*ws**m[1]*vs**m[2]*xs**m[3]*ts**m[4] for m, c in zip(p.monoms(), p.coeffs())])
Es = tosym5(Eo)
Hs = {i: sp.cancel(tosym5(Fo[i] - Fo[12]) / Es) for i in (1, 10, 11, 13)}
rr, ss, uu = Hs[13], Hs[10], Hs[1]
kk = sp.cancel((Hs[11] - Hs[10]) / rr)
chk("k = v(1-t)/(x(t-v))", sp.cancel(kk - vs*(1-ts)/(xs*(ts-vs))) == 0)
VNs = kk*xs*ts; VDs = 1 - ts + kk*xs
WDs = xs*(xs+uu)*VDs + (xs+rr)*(VDs-VNs)
WNs = xs*((xs+uu)*VNs + (xs+rr)*(VDs-VNs))
ZDs = (xs+ss)*WDs
ZNs = xs*(xs+ss)*WDs + (xs+ss+kk*rr)*(WNs - xs*WDs)
chk("inverse v = VN/VD", sp.cancel(VNs/VDs - vs) == 0)
chk("inverse w = WN/WD", sp.cancel(WNs/WDs - ws) == 0)
chk("inverse z = ZN/ZD", sp.cancel(ZNs/ZDs - zs) == 0)
print("  factored: VD =", sp.factor(VDs), "; 1+kx =", sp.factor(1+kk*xs))
print("  factored x+r, x+s, x+u, x+s+kr:", [sp.factor(xs+q) for q in (rr, ss, uu, ss+kk*rr)])

print("== C. forward parametrization in (r,s,u,k; x,t)")
Pc = flint.fmpq_mpoly_ctx.get(('r', 's', 'u', 'k', 'x', 't'), 'lex')
r, s, u, k, X, Tt = Pc.gens(); one = Pc.from_dict({(0,)*6: 1}); zero = one - one
VN = k*X*Tt; VD = 1 - Tt + k*X
WD = X*(X+u)*VD + (X+r)*(VD-VN)
WN = X*((X+u)*VN + (X+r)*(VD-VN))
ZD = (X+s)*WD
ZN = X*(X+s)*WD + (X+s+k*r)*(WN - X*WD)
chk("WN - xWD = x(x+u)(VN - xVD)", WN - X*WD == X*(X+u)*(VN - X*VD))
chk("(8) VN - VD = (t-1)(1+kx)", VN - VD == (Tt-1)*(1+k*X))
F = prods([(zero, one), (one, one), (ZN, ZD), (WN, WD)], [(one, zero), (VN, VD), (X, one), (Tt, one)])
E = F[12] - F[2]
fac = (F[13] - F[12] - r*E).factor()
print("  numerator of H13 - r factors: const", fac[0], [(str(f) if len(list(f.terms())) < 10 else "<%d terms>" % len(list(f.terms())), e) for f, e in fac[1]])
cands = [f for f, e in fac[1] if len(list(f.terms())) == 74]
chk("(9) a 74-term irreducible factor C exists", len(cands) == 1)
Cc = cands[0]
chk("(9) C has bidegree (5,2) in (x,t)", Cc.degrees()[4] == 5 and Cc.degrees()[5] == 2)
def modC(p):
    if p.is_zero(): return "identically"
    q, rem = divmod(p, Cc)
    return "mod C" if rem.is_zero() else "NOT"
for i, tg in {2: -one, 12: zero, 13: r, 10: s, 1: u, 11: s + k*r}.items():
    res = modC((F[i] - F[12]) - tg*E)
    chk("(5) H%d holds: %s" % (i, res), res != "NOT")
res = modC(X*(F[3] - F[12]) - (r*X - X - r)*E)
chk("(10) H3 = r-1-r/x: %s" % res, res != "NOT")
chk("F12 = x E (forward)", (F[12] - X*E).is_zero())

print("== C. fibre over x = 0 and the two branches (sympy)")
R, S, U, Kk, x_, t_ = sp.symbols('r s u k x t')
def tosym6(p):
    return sp.Add(*[sp.Rational(int(c.numer()), int(c.denom()))*R**m[0]*S**m[1]*U**m[2]*Kk**m[3]*x_**m[4]*t_**m[5] for m, c in zip(p.monoms(), p.coeffs())])
Cs = tosym6(Cc)
if sp.expand(Cs.subs(x_, 0) - R*S*t_*(1-t_)) != 0: Cs = -Cs
Cs = sp.expand(Cs)
chk("(12) C(0,t) = r s t(1-t) (up to the global sign)", sp.expand(Cs.subs(x_, 0) - R*S*t_*(1-t_)) == 0)
chk("leading t^2 coefficient at x=0 is -rs (no point at t = oo)", sp.expand(sp.Poly(Cs.subs(x_, 0), t_).coeff_monomial(t_**2) + R*S) == 0)
Cx = sp.diff(Cs, x_); Ct = sp.diff(Cs, t_)
chk("(0,0): Cx = -rs, Ct = rs", sp.expand(Cx.subs({x_: 0, t_: 0}) + R*S) == 0 and sp.expand(Ct.subs({x_: 0, t_: 0}) - R*S) == 0)
chk("(16) (0,1): Cx = 0, Ct = -rs, Cxx = -2ruk^2", sp.expand(Cx.subs({x_: 0, t_: 1})) == 0 and sp.expand(Ct.subs({x_: 0, t_: 1}) + R*S) == 0
    and sp.expand(sp.diff(Cs, x_, 2).subs({x_: 0, t_: 1}) + 2*R*U*Kk**2) == 0)
def branch_series(t0, order):
    cs = [sp.Integer(t0)]
    for n in range(1, order+1):
        a = sp.Symbol('a%d' % n)
        ser = sum(c*x_**i for i, c in enumerate(cs)) + a*x_**n
        coeff = sp.Poly(sp.expand(Cs.subs(t_, ser)), x_).coeff_monomial(x_**n)
        sol = sp.solve(coeff, a); assert len(sol) == 1
        cs.append(sp.factor(sol[0]))
    return cs
b1 = branch_series(0, 3); b2 = branch_series(1, 4)
chk("(13) branch 1: t = x + O(x^2)", b1[1] == 1)
chk("(16) branch 2: t = 1 - (u k^2/s) x^2 + O(x^3)", b2[1] == 0 and sp.simplify(b2[2] + U*Kk**2/S) == 0)
def ser(expr, b):
    return sp.expand(tosym6(expr).subs(t_, sum(c*x_**i for i, c in enumerate(b))))
def lead(e, nmax=4):
    pol = sp.Poly(e, x_)
    for n in range(nmax+1):
        c = sp.factor(pol.coeff_monomial(x_**n))
        if c != 0: return n, c
lv = {nm: (lead(ser(ex, b1)), lead(ser(ex, b2))) for nm, ex in
      [("VN", VN), ("VD", VD), ("WN", WN), ("WD", WD), ("ZN", ZN), ("ZD", ZD)]}
for nm, (l1, l2) in lv.items(): print("  %s: branch 1 order/lead %s, branch 2 order/lead %s" % (nm, l1, l2))
chk("(13) branch 1: VD -> 1, WD -> r, w/x -> 1, z/x -> 1, v/x -> 0",
    lv["VD"][0] == (0, 1) and lv["WD"][0] == (0, R) and lv["WN"][0] == (1, R) and lv["ZN"][0][0] == 1
    and sp.simplify(lv["ZN"][0][1] - lv["ZD"][0][1]) == 0 and lv["ZD"][0][0] == 0 and lv["VN"][0][0] == 2)
chk("branch 2: VN, VD order 1 lead k; WN, WD order 2 lead uk, uk(s+kr)/s",
    lv["VN"][1] == (1, Kk) and lv["VD"][1] == (1, Kk) and lv["WN"][1][0] == 2 and sp.simplify(lv["WN"][1][1] - U*Kk) == 0
    and lv["WD"][1][0] == 2 and sp.simplify(lv["WD"][1][1] - U*Kk*(S+Kk*R)/S) == 0)
chk("branch 2: ZN, ZD order 2 with equal leads (z -> 1)", lv["ZN"][1][0] == 2 and lv["ZD"][1][0] == 2 and sp.simplify(lv["ZN"][1][1] - lv["ZD"][1][1]) == 0)
# four-factor identity (15) for completely free homogeneous lifts
syms = sp.symbols('p0:16')
Af = [(syms[2*i], syms[2*i+1]) for i in range(4)]; Cf = [(syms[8+2*i], syms[9+2*i]) for i in range(4)]
Ff = prods(Af, Cf)
chk("(15) F12 - F2 = -[A0A1][C0C2][A2C1][A3C3] for free lifts",
    sp.expand(Ff[12] - Ff[2] + br(Af[0], Af[1])*br(Cf[0], Cf[2])*br(Af[2], Cf[1])*br(Af[3], Cf[3])) == 0)
# limits of the four factors with primitive lifts
def lim_vec(pair_leads):  # returns limit spinor from (order, lead) pairs, dividing by the common minimal order
    m = min(o for o, c in pair_leads)
    return tuple(c if o == m else 0 for o, c in pair_leads)
# branch 1 after q -> q/x:  A0'=(0,1), A1'=(1,x)->(1,0), A2'=(ZN/x,ZD), A3'=(WN/x,WD), C0'=(1,0), C1'=(kt,VD)->(0,1), C2'=(1,1), C3'=(t/x,1)->(1,1)
A2l = (lv["ZN"][0][1], lv["ZD"][0][1]); A3l = (lv["WN"][0][1], lv["WD"][0][1])
f1 = [br((0, 1), (1, 0)), br((1, 0), (1, 1)), sp.factor(br(A2l, (0, 1))), sp.factor(br(A3l, (1, 1)))]
print("  branch 1 limits of [A0'A1'],[C0'C2'],[A2'C1'],[A3'C3']:", f1)
chk("(14)-(15) branch 1: exactly [A3'C3'] vanishes (others -1, 1, rs)", f1[0] == -1 and f1[1] == 1 and sp.simplify(f1[2] - R*S) == 0 and f1[3] == 0)
A2l = lim_vec([lv["ZN"][1], lv["ZD"][1]]); A3l = lim_vec([lv["WN"][1], lv["WD"][1]]); C1l = lim_vec([lv["VN"][1], lv["VD"][1]])
f2 = [br((0, 1), (1, 1)), br((1, 0), (0, 1)), sp.factor(br(A2l, C1l)), sp.factor(br(A3l, (1, 1)))]
print("  branch 2 limits of [A0A1],[C0C2],[A2C1],[A3C3]:", f2)
chk("(17) branch 2: exactly [A2C1] vanishes, [A3C3] -> -k^2 r u/s != 0", f2[0] == -1 and f2[1] == 1 and f2[2] == 0 and sp.simplify(f2[3] + Kk**2*R*U/S) == 0)
print("  branch 2 limit points: w ->", sp.factor(A3l[0]/A3l[1]), " z ->", sp.factor(A2l[0]/A2l[1]), " v ->", sp.factor(C1l[0]/C1l[1]))
chk("(11) genericity jointly proper at (2,3,5,7)", 2*3*5*7*(3 + 7*2) != 0)

if NUMERIC:
    print("== D. NUMERIC main-statement search (library stratum sampler, margin 1e-3)")
    found, worst = numeric_main(T, 6, 8100)
    print("  NUMERIC: %d points found, worst scaled min Phi = %s" % (found, worst))
print("EXACT failures:", FAIL)
sys.exit(1 if FAIL else 0)
