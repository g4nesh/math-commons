"""Own-code replay for the Zheng n=4 audit, document final_referee_inventory79_native_uniform_20261023
(primary proof of paper case T_80 = {1,3,4,6,8,12}).  Everything is rebuilt from the definitions and
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
products = None
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
T = frozenset({1, 3, 4, 6, 8, 12})
chk("T equals paper T_80 literally", T == L.CONTACT_SETS[80])
gens = L.generators(T)
chk("generators(T) = [0,1,2,5,6,9,10,53,54,59] (eq. (1))", gens == [0, 1, 2, 5, 6, 9, 10, 53, 54, 59])
chk("only f_54 is non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [54])
chk("f_54(5) = 1 and f_54(15) = 1, so delta <= g5 and delta <= g15 (eq. (2))", L.F_VALUES[54][5] == 1 and L.F_VALUES[54][15] == 1)
ties = boundary_ties(T, 80, [54])
chk("ties forcing delta = 0 are exactly 0,5,9,11,14,15,18,19,20,21,23",
    sorted(j for j, v in ties.items() if v[0]) == [0, 5, 9, 11, 14, 15, 18, 19, 20, 21, 23])
chk("ties 2,7,10,13,16 contain an earlier-stage six-set", all(ties[j][1] for j in (2, 7, 10, 13, 16)))
chk("ties 17,22: no earlier six-set, but T+{j} is an image of Table 5 set 13 (stage 1)",
    all(not ties[j][1] and ties[j][2] == 13 for j in (17, 22)))
M = frozenset({0, 1, 2, 7, 10, 16, 17})
chk("M = {0,1,2,7,10,16,17} is Table 5 set T7_13", M == L.AUX7[13])
for j in (17, 22):
    g = L.find_g(M, T | {j})
    chk("explicit g with g.M = T+{%d}: %s" % (j, g), g is not None)
print("== C. parametrization in original variables (sympy)")
z, x, w, v, t = sp.symbols('z x w v t')
A = [(z, 1), (1, 1), (0, 1), (x, 1)]
Cc = [(w, 1), (v, 1), (1, 0), (t*v, 1)]
F = [sp.expand(f) for f in prods(A, Cc)]
E = sp.expand(F[3] - F[1])
chk("E = F3-F1 = v t (x-1)(z-w)", sp.expand(E - v*t*(x-1)*(z-w)) == 0)
H = [sp.cancel((F[i] - F[1]) / E) for i in range(24)]
r0, h0 = H[4], H[6]
s0 = sp.cancel(H[8] - H[6]); u0 = sp.cancel(H[12] - H[6]); J0 = H[5]
chk("(3) r = (1-t)/(t(x-1))", sp.cancel(r0 - (1-t)/(t*(x-1))) == 0)
chk("(3) s = (z-v)(x-tv)/E", sp.cancel(s0 - (z-v)*(x-t*v)/E) == 0)
chk("(3) u = z(1-w)(x-tv)/E", sp.cancel(u0 - z*(1-w)*(x-t*v)/E) == 0)
al0 = sp.cancel(r0 + 1 - J0)
chk("alpha = (x-tv)/(vt(x-1))", sp.cancel(al0 - (x-t*v)/(v*t*(x-1))) == 0)
Z0 = z
B00 = sp.cancel(u0 - al0*Z0)
chk("(4) v = Z((s-a)Z+u-s)/B0", sp.cancel(v - Z0*((s0-al0)*Z0 + u0 - s0)/B00) == 0)
chk("(4) w = Z(u-a)/B0", sp.cancel(w - Z0*(u0-al0)/B00) == 0)
chk("(5) x = 1+(v-1)/(h+u-uv/Z)", sp.cancel(x - 1 - (v-1)/(h0+u0-u0*v/Z0)) == 0)
chk("(5) t = 1/(1+r(x-1))", sp.cancel(t - 1/(1+r0*(x-1))) == 0)
chk("(6) B0 = a z(1-z)/(z-w)", sp.cancel(B00 - al0*z*(1-z)/(z-w)) == 0)
chk("(6) (h+u-uv/z)(x-1) = v-1", sp.cancel((h0+u0-u0*v/z)*(x-1) - (v-1)) == 0)
print("  factor alpha numerator:", sp.factor(sp.numer(sp.together(al0))))
print("  factor h:", sp.factor(h0))
print("  factor J:", sp.factor(J0))

# forward: free parameters
r, s, u, h, J, Z = sp.symbols('r s u h J Z')
al = r + 1 - J
B0 = u - al*Z
vv = Z*((s-al)*Z + u - s)/B0
ww = Z*(u-al)/B0
xx = sp.cancel(1 + (vv-1)/(h+u-u*vv/Z))
tt = sp.cancel(1/(1+r*(xx-1)))
D = u*(h+s) - Z*(h*al + u*s)
N = (Z-1)*((s-al)*Z + u)
chk("x = (D+N)/D", sp.cancel(xx - (D+N)/D) == 0)
chk("t = D/(D+rN)", sp.cancel(tt - D/(D+r*N)) == 0)
gam = u*(s+h-1)*(r-s-h)
Cp = sp.expand((al-s)*(h*al-r+s*u)*Z*(Z-1) + h*(h-J)*(al-u)*Z + gam*(1-Z))
P = sp.Poly(Cp, J, Z)
nterms = len(sp.Poly(Cp, r, s, u, h, J, Z).terms())
print("  C terms (in r,s,u,h,J,Z):", nterms, " bidegree (J,Z):", P.degree(J), P.degree(Z))
chk("(7) C has 49 terms, bidegree (2,2)", nterms == 49 and P.degree(J) == 2 and P.degree(Z) == 2)
chk("(8) B0 C = D(J B0 - D) + N(r B0 - a D)", sp.expand(B0*Cp - (D*(J*B0 - D) + N*(r*B0 - al*D))) == 0)
print("  factor C:", sp.factor_list(Cp)[1][:3] if len(sp.factor_list(Cp)[1]) > 1 else "irreducible over Q")

H15rhs = h - 1 + (s*(u+1)*Z - u*(h+s-1))/(al*Z)
# leading coefficient and fibres
chk("[J^2]C = h Z^2", sp.expand(P.coeff_monomial(J**2) if False else sp.Poly(Cp, J).coeff_monomial(J**2) - h*Z**2) == 0)
chk("C(J,0) = gamma", sp.expand(Cp.subs(Z, 0) - gam) == 0)
c1 = sp.Poly(Cp, J).coeff_monomial(J)
chk("[J^1 Z^0]C = 0 (no linear term in tangent cone)", sp.expand(c1.subs(Z, 0)) == 0)
beta = -h**2 + h*r - h*s + h*u + h - r + s*u
chk("[J Z]C = beta", sp.expand(sp.Poly(c1, Z).coeff_monomial(Z) - beta) == 0)
T0 = sp.symbols('T0')
Ct = sp.expand(T0**2*Cp.subs(J, 1/T0))
low = sum(c*T0**m[0]*Z**m[1] for m, c in sp.Poly(Ct, T0, Z).terms() if m[0]+m[1] == 2)
chk("(11) tangent cone h Z^2 + beta T0 Z + gamma T0^2", sp.expand(low - (h*Z**2 + beta*T0*Z + gam*T0**2)) == 0
    and min(m[0]+m[1] for m, c in sp.Poly(Ct, T0, Z).terms()) == 2)
chk("[J Z^2]C at h=0 is r - s u", sp.expand(sp.Poly(Cp.subs(h, 0), J, Z).coeff_monomial(J*Z**2) - (r - s*u)) == 0)
chk("(14) C at h=0, r=su", sp.expand(Cp.subs({h: 0}).subs({r: s*u}) - s*u*(s-1)*(u-1)*(1-Z)) == 0)
# full classification of C == 0 identically (as polynomial in J,Z)
coeffs = sp.Poly(Cp, J, Z).coeffs()
# C == 0 identically: h=0 forced by [J^2 Z^2]; then r=su by [J Z^2]; then (14)
c22 = sp.Poly(Cp, J, Z).coeff_monomial(J**2*Z**2)
chk("[J^2 Z^2]C = h", sp.expand(c22 - h) == 0)
# H15 on h=0, r=su strata with fixed Z: limit J->oo
for sub in ({h: 0, r: s*u, s: 1}, {h: 0, r: s*u, u: 1}):
    lim = sp.limit(H15rhs.subs(sub), J, sp.oo)
    chk("H15 -> -1 as J->oo on stratum %s" % sub, sp.simplify(lim + 1) == 0)

print("== C. the 24 products from the inverse, modulo C (python-flint)")
ctx = flint.fmpq_mpoly_ctx.get(('r', 's', 'u', 'h', 'J', 'Z'), 'lex')
rq, sq, uq, hq, Jq, Zq = ctx.gens()
oneq = ctx.from_dict({(0,)*6: 1}); zeroq = oneq - oneq
alq = rq + 1 - Jq
B0q = uq - alq*Zq
Dq = uq*(hq+sq) - Zq*(hq*alq + uq*sq)
Nq = (Zq-1)*((sq-alq)*Zq + uq)
vNq = Zq*((sq-alq)*Zq + uq - sq)
Aq = [(Zq, oneq), (oneq, oneq), (zeroq, oneq), (Dq+Nq, Dq)]
Cq = [(Zq*(uq-alq), B0q), (vNq, B0q), (oneq, zeroq), (Dq*vNq, (Dq+rq*Nq)*B0q)]
Fq = prods(Aq, Cq)
gamq = uq*(sq+hq-1)*(rq-sq-hq)
Cpq = (alq-sq)*(hq*alq-rq+sq*uq)*Zq*(Zq-1) + hq*(hq-Jq)*(alq-uq)*Zq + gamq*(1-Zq)
Eq = Fq[3] - Fq[1]
chk("C irreducible over Q", len(Cpq.factor()[1]) == 1 and Cpq.factor()[1][0][1] == 1)
def divisible(p):
    if p.is_zero(): return "identically"
    q, rem = divmod(p, Cpq)
    return "mod C" if rem.is_zero() else "NOT"
for i, tg in {1: zeroq, 3: oneq, 4: rq, 5: Jq, 6: hq, 8: hq+sq, 12: hq+uq}.items():
    res = divisible((Fq[i] - Fq[1]) - tg*Eq)
    chk("(9) H%d holds: %s" % (i, res), res != "NOT")
res = divisible(alq*Zq*(Fq[15] - Fq[1]) - (alq*Zq*(hq-1) + sq*(uq+1)*Zq - uq*(hq+sq-1))*Eq)
chk("(13) H15 formula: %s" % res, res != "NOT")
q, rem = divmod((Fq[5] - Fq[1]) - Jq*Eq, Cpq)
print("  numerator of H5 - J = C * cofactor, cofactor =", q.factor())

if NUMERIC:
    print("== D. NUMERIC main-statement search (library stratum sampler, margin 1e-3)")
    found, worst = numeric_main(T, 6, 8000)
    print("  NUMERIC: %d points found, worst scaled min Phi = %s" % (found, worst))
print("EXACT failures:", FAIL)
sys.exit(1 if FAIL else 0)
