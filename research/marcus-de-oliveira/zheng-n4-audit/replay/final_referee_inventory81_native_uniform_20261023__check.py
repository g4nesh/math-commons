"""Own-code replay for the Zheng n=4 audit, document final_referee_inventory81_native_uniform_20261023
(primary proof of paper case T_82 = {2,5,9,10,12,13}).  Everything is rebuilt from the definitions and
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

print("== A/B. contact set, generators, barriers, boundary ties")
T = frozenset({2, 5, 9, 10, 12, 13})
chk("T equals paper T_82 literally", T == L.CONTACT_SETS[82])
gens = L.generators(T)
chk("generators(T) = [0,1,12,39,45,67,70,72] (eq. (1))", gens == [0, 1, 12, 39, 45, 67, 70, 72])
chk("only f_70 is non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [70])
rows = [18, 16, 8, 4, 3, 14, 11, 7]
chk("coefficients (q,p,j,c,e,b,d,a) = (g18,g16,g8,g4,g3,g14,g11,g7): rows form the identity",
    all(L.F_VALUES[f][rw] == (1 if a == b else 0) for a, f in enumerate(gens) for b, rw in enumerate(rows)))
name = dict(zip(gens, "q p j c e b d a".split()))
for bidx, want in ((6, "j d a"), (15, "c b d"), (19, "q j d")):
    got = sorted(name[f] for f in gens if L.F_VALUES[f][bidx] == 1)
    chk("(2) g%d = %s (all with coefficient 1)" % (bidx, " + ".join(want.split())),
        got == sorted(want.split()) and all(L.F_VALUES[f][bidx] in (0, 1) for f in gens))
ties = boundary_ties(T, 82, [70])
chk("zeros of f_70 outside T are 3,4,7,8,14,16,18", sorted(j for j, v in ties.items() if not v[0]) == [3, 4, 7, 8, 14, 16, 18])
chk("each of them contains an earlier-stage six-set", all(ties[j][1] for j in (3, 4, 7, 8, 14, 16, 18)))
for j, lab in zip([3, 4, 7, 8, 14, 16, 18], [33, 59, 59, 30, 78, 59, 78]):
    T7 = T | {j}
    hit = [d for d in sorted(T7) if L.find_g(L.CONTACT_SETS[lab + 1], T7 - {d}) is not None]
    chk("certificate label T%d read as inventory index = paper T_%d (stage %d) maps into T+{%d}" % (lab, lab + 1, ST6[lab + 1], j),
        bool(hit) and ST6[lab + 1] < ST6[82])

print("== C. original variables: A = (0,1,z,w), C = (oo,v,x,t)")
zs, ws, vs, xs, ts = sp.symbols('z w v x t')
Fo = [sp.expand(f) for f in prods([(0, 1), (1, 1), (zs, 1), (ws, 1)], [(1, 0), (vs, 1), (xs, 1), (ts, 1)])]
Eo = sp.expand(Fo[12] - Fo[2])
chk("E = (z-v)(w-t)", sp.expand(Eo - (zs-vs)*(ws-ts)) == 0)
chk("(3) F12 = x E", sp.expand(Fo[12] - xs*Eo) == 0)
chk("ratio F5/F13", sp.cancel(Fo[5]/Fo[13] + (1-ts)*(zs-xs)/(xs*(zs-ts))) == 0)
chk("ratio F9/F13", sp.cancel(Fo[9]/Fo[13] - vs*(1-xs)/(xs*(ws-vs))) == 0)
chk("ratio F10/F9", sp.cancel(Fo[10]/Fo[9] - (1-ts)*(ws-xs)/((1-xs)*(zs-ts))) == 0)
chk("(12) x F0 = (x-z) F14", sp.expand(xs*Fo[0] - (xs-zs)*Fo[14]) == 0)
Ho = {i: sp.cancel((Fo[i] - Fo[12])/Eo) for i in (5, 9, 10, 13)}
def inv(r, s, u, h, x, z):
    Td = z*(x+r) + x*(s-r); Tn = z*(x**2 + x*(s+1) + r) - x*(x+r)
    Wd = x*(x+s)*(x+u); K = (x+h)*(x-1)*(x+r); Wn = x*Wd + K*(z-x); V = r + x*(1+u-r)
    return Td, Tn, Wd, K, Wn, V
r0, s0, u0, h0 = Ho[13], Ho[5], Ho[9], Ho[10]
Td, Tn, Wd, K, Wn, V = inv(r0, s0, u0, h0, xs, zs)
chk("(4) t = Tn/Td", sp.cancel(Tn/Td - ts) == 0)
chk("(4) w = Wn/Wd", sp.cancel(Wn/Wd - ws) == 0)
chk("(4) v = x(x+u)w/V", sp.cancel(xs*(xs+u0)*ws/V - vs) == 0)
chk("Td - Tn = x(x+s)(1-z)", sp.cancel(Td - Tn - xs*(xs+s0)*(1-zs)) == 0)
print("  factored x+r, x+s, x+u, x+h:", [sp.factor(xs+q) for q in (r0, s0, u0, h0)], " V:", sp.factor(V))
print("  factored r = H13:", sp.factor(r0))

print("== C. forward parametrization in (r,s,u,h; x,z)")
Pc = flint.fmpq_mpoly_ctx.get(('r', 's', 'u', 'h', 'x', 'z'), 'lex')
r, s, u, h, X, Z = Pc.gens(); one = Pc.from_dict({(0,)*6: 1}); zero = one - one
Td, Tn, Wd, K, Wn, V = inv(r, s, u, h, X, Z)
F = prods([(zero, one), (one, one), (Z, one), (Wn, Wd)], [(one, zero), (X*(X+u)*Wn, Wd*V), (X, one), (Tn, Td)])
E = F[12] - F[2]
chk("(3) F12 = x E forward", (F[12] - X*E).is_zero())
fac = (F[13] - F[12] - r*E).factor()
print("  numerator of H13 - r factors: const", fac[0], [(str(f) if len(list(f.terms())) < 10 else "<%d terms>" % len(list(f.terms())), e) for f, e in fac[1]])
cands = [f for f, e in fac[1] if len(list(f.terms())) == 355]
chk("a 355-term irreducible factor C exists", len(cands) == 1)
Cc = cands[0]
chk("C bidegree (6,2) in (x,z)", Cc.degrees()[4] == 6 and Cc.degrees()[5] == 2)
def modC(p):
    if p.is_zero(): return "identically"
    q, rem = divmod(p, Cc)
    return "mod C" if rem.is_zero() else "NOT"
for i, tg in {2: -one, 12: zero, 13: r, 5: s, 9: u, 10: h}.items():
    res = modC((F[i] - F[12]) - tg*E); chk("contact H%d holds: %s" % (i, res), res != "NOT")
res = modC((X-1)*(X+r)*(F[11] - F[12]) - (X*(X+s)*(X+u) - X*(X-1)*(X+r))*E); chk("(5) H11 formula: %s" % res, res != "NOT")
res = modC(X*(Z-1)*(F[19] - F[12]) - (Z*(X*(s+1) + r) - X*r)*E); chk("(9) H19 formula: %s" % res, res != "NOT")
res = modC((X-1)*Wn*(F[15] - F[12]) - (X*(u+1)*Wn - V*Wd)*E); chk("(10) H15 formula: %s" % res, res != "NOT")
vn, vd = X*(X+u)*Wn, Wd*V
res = modC((Z*vd - vn)*(F[6] - F[12]) - (vn*(Z-X) - X*(Z*vd - vn))*E); chk("(11) H6 formula: %s" % res, res != "NOT")

print("== C. sections, local branches, limits (sympy)")
R, S, U, Hh, x, z, tau = sp.symbols('r s u h x z tau')
def tosym(p):
    return sp.Add(*[sp.Rational(int(c.numer()), int(c.denom()))*R**m[0]*S**m[1]*U**m[2]*Hh**m[3]*x**m[4]*z**m[5] for m, c in zip(p.monoms(), p.coeffs())])
Cs = sp.expand(tosym(Cc))
if sp.expand(Cs.subs(x, 1) + (U+1)**2*(S+1)**3*(z-1)) != 0: Cs = -Cs
Cs = sp.expand(Cs)
chk("(7) C(1,z) = -(u+1)^2(s+1)^3(z-1) (up to the global sign)", sp.expand(Cs.subs(x, 1) + (U+1)**2*(S+1)**3*(z-1)) == 0)
chk("(8) C(-r,z) = -r^3(r-u)^2(r-s)^3(z+r)", sp.expand(Cs.subs(x, -R) + R**3*(R-U)**2*(R-S)**3*(z+R)) == 0)
chk("z-degree of C is 2 (so each fibre over x=1, x=-r is one finite point and z = oo)", sp.Poly(Cs, z).degree() == 2)
Ct = sp.expand(tau**2*Cs.subs(z, 1/tau))
for x0 in (1, -R):
    d = sp.factor(sp.diff(Ct, tau).subs({x: x0, tau: 0}))
    print("  d_tau Ct(%s,0) = %s" % (x0, d))
    chk("simple point at (x,z) = (%s, oo)" % x0, sp.expand(d) != 0)
dtau = sp.diff(Ct, tau).subs({x: 1, tau: 0}); dx = sp.diff(Ct, x).subs({x: 1, tau: 0})
chk("d_tau Ct(1,0) = -(s+1)^3(u+1)^2 and d_x Ct(1,0) = (h+1)(r+1)^2(s+1)(u+1)",
    sp.expand(dtau + (S+1)**3*(U+1)**2) == 0 and sp.expand(dx - (Hh+1)*(R+1)**2*(S+1)*(U+1)) == 0)
kap = sp.factor(-dx/dtau)
chk("tau/eps -> kappa = (h+1)(r+1)^2/((s+1)^2(u+1))", sp.simplify(kap - (Hh+1)*(R+1)**2/((S+1)**2*(U+1))) == 0)
_, _, Wd_s, K_s, Wn_s, V_s = inv(R, S, U, Hh, x, z)
wlim = sp.limit(x + K_s*(1/(kap*(x-1)) - x)/Wd_s, x, 1)
chk("at (1,oo): w -> (s+r+2)/(r+1), v/w -> 1, H6 -> v - 1 real",
    sp.simplify(wlim - (S+R+2)/(R+1)) == 0 and sp.simplify((x*(x+U)/V_s).subs(x, 1) - 1) == 0)
cz = sp.diff(Cs, z).subs({x: 1, z: 1})
chk("(1,1) smooth: C_z(1,1) = -(s+1)^3(u+1)^2 != 0", sp.expand(cz + (S+1)**3*(U+1)**2) == 0)
a1 = sp.Symbol('a1'); zz = 1 + a1*(x-1)
wv = x + K_s.subs(z, zz)*(zz - x)/Wd_s.subs(z, zz)
chk("at (1,1): H15 -> u+r+1 for every slope of the branch", sp.simplify(sp.limit((x*(U+1)*wv - V_s)/((x-1)*wv), x, 1) - (U+R+1)) == 0)
H19 = (z*(x*(S+1)+R) - x*R)/(x*(z-1))
chk("at (-r,-r): H19 = r(s+1)/(r+1)", sp.simplify(H19.subs({x: -R, z: -R}) - R*(S+1)/(R+1)) == 0)
chk("at (-r,oo): H19 -> s", sp.simplify(sp.limit(H19.subs(x, -R), z, sp.oo) - S) == 0)
cz2 = sp.diff(Cs, z).subs({x: -R, z: -R})
chk("(-r,-r) smooth: C_z != 0", sp.expand(cz2) != 0)
H11 = x*(x+S)*(x+U)/((x-1)*(x+R)) - x
chk("(5) H11 -> s+u-r+1 at x = oo", sp.simplify(sp.limit(H11, x, sp.oo) - (S+U-R+1)) == 0)
chk("(6) genericity jointly proper at (2,3,5,7)", 2*3*6*4*8*(2-3)*(2-5) != 0)

print("== C. identically-zero specializations")
def coeff(e, i, j): return sp.Poly(e, x, z).coeff_monomial(x**i*z**j)
c = Cs.subs({S: -1, R: -1})
chk("r=s=-1: [x^2 z^0]C = -h(u+1)", sp.expand(coeff(c, 2, 0) + Hh*(U+1)) == 0)
chk("r=s=-1, h=0: [x^6 z]C = -1", sp.expand(coeff(c.subs(Hh, 0), 6, 1) + 1) == 0)
c = Cs.subs({U: -1, R: -1})
chk("r=u=-1: [x^2 z^0]C = h^2(s+1)", sp.expand(coeff(c, 2, 0) - Hh**2*(S+1)) == 0)
chk("r=u=-1, h=0: [x^6 z^2]C - [x z^2]C = 2", sp.expand(coeff(c.subs(Hh, 0), 6, 2) - coeff(c.subs(Hh, 0), 1, 2) - 2) == 0)
for sub, val in (({S: -1, R: 0}, U), ({S: -1, R: U}, 0), ({U: -1, R: 0}, S), ({U: -1, R: S}, 0), ({R: -1, S: -1, U: -1}, 0)):
    chk("H11 on stratum %s is the real constant %s" % (sub, val), sp.simplify(H11.subs(sub) - val) == 0)
zero_any = [sub for sub in ({S: -1, R: 0}, {S: -1, R: U}, {U: -1, R: 0}, {U: -1, R: S}, {R: -1, S: -1, U: -1}) if sp.expand(Cs.subs(sub)) == 0]
print("  INFO: strata on which C is the zero polynomial for all h (none expected):", zero_any)

if NUMERIC:
    print("== D. NUMERIC main-statement search (library stratum sampler, margin 1e-3)")
    found, worst = numeric_main(T, 6, 8200)
    print("  NUMERIC: %d points found, worst scaled min Phi = %s" % (found, worst))
print("EXACT failures:", FAIL)
sys.exit(1 if FAIL else 0)
