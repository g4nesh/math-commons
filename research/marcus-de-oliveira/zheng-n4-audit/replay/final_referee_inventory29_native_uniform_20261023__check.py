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

def radial_identity(jsonname, names, idx, npts=8, seed=12345):
    """Own cubic (lib cubic_C, paper (A.2)) on g + z*kappa, g = sum of the named generators.
    Compares C, Q with the JSON's C_terms/Q_terms, and checks
    Q^4 det Hess(g - C/Q kappa) = -4 Q^2 G^2 at npts random rational points (randomized identity test),
    with G the degree-nine factor read from the JSON as data.  Returns the JSON."""
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
        pt = [sp.Rational(rng.randint(0, 40), rng.randint(1, 9)) for _ in syms]
        pt[names.index("d")] = sp.Integer(1)
        sb = dict(zip(syms, pt)); Qv = Qm.subs(sb)
        if Qv == 0: continue
        zv = -Cm.subs(sb) / Qv
        mv = [sp.nsimplify(x.subs(sb).subs(z, zv)) for x in m]
        dv = sp.Matrix(Hf(*mv)).det(method="bareiss"); Gv = G.subs(sb)
        ratios.add(Qv**4 * dv / (Qv**2 * Gv**2) if Gv != 0 else None)
    check(f"{jsonname}: Q^4 det H(g - C/Q kappa) = -4 Q^2 G^2 at {npts} random rational points", ratios == {-4})
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

# ===== final_referee_inventory29_native_uniform_20261023 : paper case T_30 = {2,5,9,12,13,14} =====
T = L.CONTACT_SETS[30]
check("working contacts {2,5,9,12,13,14} are the paper's T_30 literally", sorted(T) == [2, 5, 9, 12, 13, 14])
g0 = L.find_g({0, 1, 3, 7, 11, 12}, T)
check(f"stated contact set {{0,1,3,7,11,12}} is G-equivalent to T_30 (g = {g0})", g0 is not None)
S = [0, 1, 5, 12, 39, 41, 43, 45, 70, 72]
check("S = (0,1,5,12,39,41,43,45,70,72) = lib generators(T_30)", L.generators(T) == S)
check("non-minor generators are exactly 43 and 70", [j for j in S if not L.F_IS_MINOR[j]] == [43, 70])
check("contact kernel has dimension 8", 14 - L.rank([L.P[s] for s in T]) == 8)
for name, pos in (("g11", [11]), ("g6", [6]), ("g15+g16", [15, 16]), ("g19+g16", [19, 16])):
    check(f"(2) Theta <= {name}: every f_i in S has value >= [i in {{43,70}}] there",
          all(sum(V[jj][s_] for s_ in pos) >= (1 if jj in (43, 70) else 0) for jj in S))
# extra contacts: resolve to earlier-stage cases (paper Table 6; T_30 is in stage 9)
STAGE = {}
for st, lst in ((1, [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88]), (2, [10,26,46,90]),
                (3, [9,11,12,47]), (4, [8,13,17,49]), (5, [14,16,24,33,38,70]), (6, [52,57,58,65,68,71,76,78]),
                (7, [45,51,66,69,75,81]), (8, [48,50,53,59,62,67,73,79,86]), (9, [30,41,54,55,60,61,63,74,77,80,83,84]),
                (10, [56,64,82,85,87]), (11, [32]), (12, [89])):
    for nu in lst: STAGE[nu] = st
check("stage table covers T_1..T_90 once", sorted(STAGE) == list(range(1, 91)))
order = sorted(STAGE, key=lambda nu: (STAGE[nu], nu))
for k in range(24):
    if k in T: continue
    if V[43][k] > 0 and V[70][k] > 0:
        check(f"extra contact {k}: f43, f70 both positive there, so both non-minor coefficients vanish", True); continue
    found = None
    for nu in order:
        if STAGE[nu] >= 9: break
        if any(L.g_act_set(g, L.CONTACT_SETS[nu]) <= (T | {k}) for g in L.G_ELEMENTS):
            found = nu; break
    check(f"extra contact {k} (f43={V[43][k]}, f70={V[70][k]}): T_30 + {{{k}}} contains a G-image of T_{found} (stage {STAGE.get(found)} < 9)", found is not None)

# ---- curve identities, from the definitions: minus-bracket products, A=(0,1,z,w), C=(oo,v,x,t)
z, w, v, x, t, r, s, u, h = sp.symbols("z w v x t r s u h")
Aa = [0, 1, z, w]; Cc = [None, v, x, t]
def Fk(k):
    sig = L.PERMS[k]; pr = sp.Integer(-1)          # [a, oo] = -1 in homogeneous brackets
    for i in range(4):
        if sig[i] != 0: pr *= (Aa[i] - Cc[sig[i]])
    return sp.expand(pr)
Fs = [Fk(k) for k in range(24)]
E = Fs[12] - Fs[2]
check("E = F12 - F2 = (z-v)(w-t)", sp.expand(E - (z - v) * (w - t)) == 0)
H = lambda i: (Fs[i] - Fs[12]) / E
check("H2 = -1, H12 = 0, F12/E = x", sp.simplify(H(2) + 1) == 0 and sp.simplify(H(12)) == 0 and sp.simplify(Fs[12] / E - x) == 0)
check("H14 = x(1-z)/(z-v)", sp.simplify(H(14) - x * (1 - z) / (z - v)) == 0)
check("F5/F13 = -(1-t)(z-x)/(x(z-t)) and F9/F13 = v(1-x)/(x(w-v))",
      sp.simplify(Fs[5] / Fs[13] + (1 - t) * (z - x) / (x * (z - t))) == 0 and sp.simplify(Fs[9] / Fs[13] - v * (1 - x) / (x * (w - v))) == 0)
Vx = r + x * (1 + u - r); Td = z * (x + r) + x * (s - r); Tn = z * (x**2 + x * (s + 1) + r) - x * (x + r)
check("Td - Tn = x(x+s)(1-z)", sp.expand(Td - Tn - x * (x + s) * (1 - z)) == 0)
sub = {t: Tn / Td, v: ((x + h) * z - x) / h, w: Vx * ((x + h) * z - x) / (h * x * (x + u))}
Hs = {i: sp.cancel(H(i).subs(sub, simultaneous=True)) for i in (5, 6, 9, 11, 13, 14, 15, 16, 19)}
check("under (4): H14 = h identically", sp.simplify(Hs[14] - h) == 0)
check("under (4): (H5+x)/(H13+x) = (x+s)/(x+r) and (H9+x)/(H13+x) = (x+u)/(x+r) identically",
      sp.simplify((Hs[5] + x) * (x + r) - (x + s) * (Hs[13] + x)) == 0 and sp.simplify((Hs[9] + x) * (x + r) - (x + u) * (Hs[13] + x)) == 0)
fl = sp.factor_list(sp.numer(sp.together(Hs[13] - r)))
Cpoly = [f for f, m_ in fl[1] if sp.Poly(f, x, z).degree(x) == 3 and sp.Poly(f, x, z).degree(z) == 2]
others = [f for f, m_ in fl[1] if not (sp.Poly(f, x, z).degree(x) == 3 and sp.Poly(f, x, z).degree(z) == 2)]
check("numerator of H13 - r = const * (x+r) * C, C of bidegree (3,2) in (x,z) with 61 terms",
      len(Cpoly) == 1 and [sp.expand(f) for f in others] in ([x + r], [r + x]) and len(sp.Poly(Cpoly[0], x, z, r, s, u, h).terms()) == 61)
C = sp.expand(Cpoly[0])
check("C is irreducible over Q(r,s,u,h)[x,z]", len(sp.factor_list(C, x, z)[1]) == 1)
def zero_mod_C(expr):
    n = sp.expand(sp.numer(sp.together(expr)))
    return sp.expand(sp.prem(sp.Poly(n, z), sp.Poly(C, z)).as_expr()) == 0
check("H5 = s, H9 = u, H13 = r hold modulo C", zero_mod_C(Hs[5] - s) and zero_mod_C(Hs[9] - u) and zero_mod_C(Hs[13] - r))
H11c = x * (x + s) * (x + u) / ((x - 1) * (x + r)) - x
check("(5) H11 = x(x+s)(x+u)/((x-1)(x+r)) - x modulo C", zero_mod_C(Hs[11] - H11c))
check("(5) H11 -> s+u-r+1 (real) as x -> oo", sp.limit(H11c, x, sp.oo) == s + u - r + 1)
sgn = sp.Integer(-1)   # our C equals -(the document's C); the sign is checked here
check("(7) C(1,z) = (u+1)(z-1)((r+1)(h+1)(z-1)+s+1) (up to the global sign of C)",
      sp.expand(sgn * C.subs(x, 1) - (u + 1) * (z - 1) * ((r + 1) * (h + 1) * (z - 1) + s + 1)) == 0)
check("(8) C(-r,z) = -r^3 (r-u)(r-s)(z-1) (same global sign)", sp.expand(sgn * C.subs(x, -r) + r**3 * (r - u) * (r - s) * (z - 1)) == 0)
H6c = z * (x * (h + 1) - (x + h) * z) / (x * (z - 1))
H15c = (x * (u + 1) * sub[w] - Vx) / ((x - 1) * sub[w]); H16c = x * (s + z) / (x - z); H19c = (z * (x * (s + 1) + r) - x * r) / (x * (z - 1))
check("(9) H6 formula holds identically under (4)", sp.simplify(Hs[6] - H6c) == 0)
check("(10)-(12) H15, H16, H19 formulas hold modulo C", zero_mod_C(Hs[15] - H15c) and zero_mod_C(Hs[16] - H16c) and zero_mod_C(Hs[19] - H19c))
Cx = sp.factor(sgn * sp.diff(C, x).subs({x: 1, z: 1})); Cz = sp.factor(sgn * sp.diff(C, z).subs({x: 1, z: 1}))
check("C_x(1,1) = -h(s+1)(r+1), C_z(1,1) = (s+1)(u+1)", sp.expand(Cx + h * (s + 1) * (r + 1)) == 0 and sp.expand(Cz - (s + 1) * (u + 1)) == 0)
zp = -Cx / Cz
check("z'(1) = h(r+1)/(u+1)", sp.simplify(zp - h * (r + 1) / (u + 1)) == 0)
lim = sp.limit(H6c.subs(z, 1 + zp * (x - 1)), x, 1)
check("H6 -> (u+1)/(r+1) - h - 1 at x=z=1 along the branch (0/0 of first order)", sp.simplify(lim - ((u + 1) / (r + 1) - h - 1)) == 0)
zs = 1 - (s + 1) / ((r + 1) * (h + 1))
check("H6 at x=1, z=z* is real: -z*(h+1)", sp.simplify(H6c.subs({x: 1, z: zs}) + zs * (h + 1)) == 0)
check("x=-r, z=1: w -> 1, H15 -> r, H16 -> r(s+1)/(r+1)",
      sp.simplify(sub[w].subs({x: -r, z: 1}) - 1) == 0 and sp.simplify(H15c.subs({x: -r, z: 1}) - r) == 0 and sp.simplify(H16c.subs({x: -r, z: 1}) - r * (s + 1) / (r + 1)) == 0)
check("x=-r, z->oo: H16 -> r, H19 -> s", sp.limit(H16c.subs(x, -r), z, sp.oo) == r and sp.limit(H19c.subs(x, -r), z, sp.oo) == s)
check("x F0 = (x - z) F14", sp.expand(x * Fs[0] - (x - z) * Fs[14]) == 0)
Pc = sp.Poly(sgn * C, x, z)
check("walls: [x^2 z^0]C = s+1 at r=u=-1 and [x^3 z^0]C = h(u+1) at r=s=-1",
      sp.expand(Pc.coeff_monomial(x**2).subs({r: -1, u: -1}) - (s + 1)) == 0 and sp.expand(Pc.coeff_monomial(x**3).subs({r: -1, s: -1}) - h * (u + 1)) == 0)
check("walls: H11 is the real constant s, 0, u, 0 at (u,r)=(-1,0), (u=-1,r=s), (s,r)=(-1,0), (s=-1,r=u); 0 at r=s=u=-1",
      [sp.simplify(H11c.subs(dd)) for dd in ({u: -1, r: 0}, {u: -1, r: s}, {s: -1, r: 0}, {s: -1, r: u}, {r: -1, s: -1, u: -1})] == [s, 0, u, 0, 0])

print("SUMMARY:", "ALL EXACT CHECKS PASSED" if not FAILED else f"{len(FAILED)} FAILED: {FAILED}")
sys.exit(1 if FAILED else 0)
