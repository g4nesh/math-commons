"""Own-code replay for final_referee_inventory68_native_harmonic_20261023 (primary proof of paper
case T_69 = {0,2,5,7,9,12}). Audit batch C37-T69-T71.
Uses only the shared audit library (../lib/moc4lib.py) and the audit's own dependency map
(../deps/dependency-map.json, for Appendix G stages). Executes nothing from the supplement.
Every formula is re-derived from the definitions F_sigma = prod_i (A_i + B_sigma(i)), with the
factor for B = infinity equal to 1."""
import sys, os, json, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)
def info(name, s):
    print("INFO " + name + " :: " + str(s))

STAGE = {c["case"]: c["stage_paper"] for c in json.load(open(os.path.join(HERE, "..", "deps", "dependency-map.json")))["cases"]}
MINOR_ORBIT_CACHE = {}

def comb(coeffs):
    out = [Fr(0)] * 14
    for j, c in coeffs.items():
        for k in range(14):
            out[k] += c * L.F_VECS[j][k]
    return out

def case_of(S6):
    """paper case number whose contact set is G-equivalent to S6, or 'minor' if every generator
    of S6 is a minor vector, else None"""
    if all(L.F_IS_MINOR[j] for j in L.generators(S6)):
        return "minor"
    rep = L.orbit_rep(S6)
    for nu, T in L.CONTACT_SETS.items():
        if L.orbit_rep(T) == rep:
            return nu
    return None

def tie_covers(T, sigma):
    """all six-subsets of T+{sigma} containing sigma, resolved to paper cases"""
    S = sorted(set(T) | {sigma})
    res = []
    for drop in T:
        S6 = [x for x in S if x != drop]
        res.append((drop, case_of(S6)))
    return res

def cone_checks(T, src, bad, rel_lhs, rel_rhs, keep, meas, dmeas, stage_self, name):
    """generic cone / decomposition / tie checks"""
    G = L.generators(T)
    check(name + ": generators(T) equal the printed source", sorted(G) == sorted(src), G)
    nonminor = [j for j in G if not L.F_IS_MINOR[j]]
    check(name + ": sole non-minor generator", nonminor == [bad], nonminor)
    lhs = L.values(comb({j: 1 for j in rel_lhs})); rhs = L.values(comb({j: 1 for j in rel_rhs}))
    check(name + ": relation in 24 profiles", lhs == rhs)
    check(name + ": relation in 14 coordinates", comb({j: 1 for j in rel_lhs}) == comb({j: 1 for j in rel_rhs}))
    # rank: the 8 kept profiles independent and span {c : (Pc)_T = 0}
    M = [list(L.F_VECS[j]) for j in keep]
    check(name + ": the eight kept generators are independent", L.rank(M) == 8)
    Tspace = L.nullspace([L.P[p] for p in sorted(T)])
    check(name + ": span(kept) = {c : (Pc)_pi = 0 on T} (dim 8)", len(Tspace) == 8 and L.rank(M + [list(v) for v in Tspace]) == 8)
    # measurement functionals: for every kept generator f_j, coefficient readout
    ok = True
    for j in keep:
        v = L.F_VALUES[j]
        for (i, sig) in meas:
            if v[sig] != (1 if i == j else 0):
                ok = False
        for (plus, minus) in dmeas:
            if v[plus] - v[minus] != (1 if j == bad else 0):
                ok = False
    check(name + ": coefficient readouts C_i = g_sigma and d = g_p - g_q on the kept basis", ok)
    # entire-cost absorption: exhausting the resource gives either all-minor or (2) (structural)
    # added ties
    rows = []
    allok = True
    for sig in range(24):
        if sig in T:
            continue
        if L.F_VALUES[bad][sig] > 0:
            # generators of T+{sig} are the source generators vanishing at sig: all minor?
            Gs = L.generators(set(T) | {sig})
            okm = all(L.F_IS_MINOR[j] for j in Gs)
            rows.append((sig, "f_bad>0", "all-minor" if okm else "NOT all minor"))
            allok &= okm
        else:
            cov = tie_covers(T, sig)
            earlier = [c for (_, c) in cov if c == "minor" or (isinstance(c, int) and STAGE["T_%d" % c] < stage_self)]
            rows.append((sig, "f_bad=0", cov, "earlier:" + str(sorted(set(map(str, earlier))))))
            allok &= bool(earlier)
    for r in rows:
        info(name + ": tie", r)
    check(name + ": every added tie covered (all-minor cone, or an earlier-stage case)", allok)
    return rows

# ---------------------------------------------------------------------------------------------
T = [0, 2, 5, 7, 9, 12]
check("T_69 literal equals Appendix E", frozenset(T) == L.CONTACT_SETS[69])
info("stage of T_69 (Table 6)", STAGE["T_69"])
src = [0, 1, 2, 5, 7, 12, 39, 40, 45]
keep = [0, 2, 5, 7, 12, 39, 40, 45]
meas = [(0, 23), (2, 13), (5, 11), (7, 16), (12, 6), (39, 17), (45, 3)]
rows = cone_checks(T, src, 40, [40, 1], [0, 7, 39], keep, meas, [(14, 16), (4, 17)], STAGE["T_69"], "T69")

# ---------------------------------------------------------------------------------------------
# chart
V, W, r, s, u, h = sp.symbols("V W r s u h")
a_, t_, z_, x_, w_ = sp.symbols("a t z x w")
def Fsig(A, B):
    """F_sigma for all 24 sigma; B[3] = infinity contributes factor 1"""
    out = []
    for p in L.PERMS:
        e = 1
        for i in range(4):
            if p[i] != 3:
                e *= (A[i] + B[p[i]])
        out.append(e)
    return out
F0 = Fsig([a_, z_, x_, w_], [t_, 0, 1, None])
E0 = sp.factor(F0[2] - F0[0])
check("E = F2 - F0 = (a+t)(x-z)", sp.expand(E0 - (a_ + t_) * (x_ - z_)) == 0, E0)
check("F5-F0 = (a+t)(w-z)(x+1)", sp.expand(F0[5] - F0[0] - (a_ + t_) * (w_ - z_) * (x_ + 1)) == 0)
check("F12-F2 = x(a-z)(t-1)", sp.expand(F0[12] - F0[2] - x_ * (a_ - z_) * (t_ - 1)) == 0)
check("F9-F7 = a(t-1)(z-w)", sp.expand(F0[9] - F0[7] - a_ * (t_ - 1) * (z_ - w_)) == 0)

M = u * (W - 1) + s * r * V
xs = W - 1; zs = W * (1 - V) - 1; ws = W - 1 + (r - W) * V
as_ = u * zs * (W - 1) / M
ts = as_ * (r - u * W) / (as_ * r + u * W)
sub = {a_: as_, t_: ts, z_: zs, x_: xs, w_: ws}
F = [sp.cancel(f.subs(sub)) for f in F0]
E = sp.cancel(F[2] - F[0])
H = [sp.factor(sp.cancel((F[k] - F[0]) / E)) for k in range(24)]
check("inverse (2): H5 = r", sp.simplify(H[5] - r) == 0)
check("inverse (2): H12 = 1+s", sp.simplify(H[12] - 1 - s) == 0)
check("inverse (2): H9 - H7 = u", sp.simplify(H[9] - H[7] - u) == 0)
check("inverse (2): V = (x-z)/(x+1), W = x+1", sp.simplify((xs - zs) / (xs + 1) - V) == 0)
# forward-then-back: generic numeric spectra -> (r,s,u,V,W) -> (2) recovers spectra
rng = random.Random(7)
okrt = True
for _ in range(5):
    vals = {a_: sp.Rational(rng.randint(-50, 50), 7) + sp.I * sp.Rational(rng.randint(-50, 50), 11),
            t_: sp.Rational(rng.randint(-50, 50), 3), z_: sp.Rational(rng.randint(-50, 50), 13) + 2 * sp.I,
            x_: sp.Rational(rng.randint(-50, 50), 5) - sp.I, w_: sp.Rational(rng.randint(-50, 50), 17)}
    Fv = [sp.nsimplify(f.subs(vals)) for f in F0]
    Ev = Fv[2] - Fv[0]
    rv, sv, uv = [sp.simplify((Fv[i] - Fv[j]) / Ev) for (i, j) in ((5, 0), (12, 2), (9, 7))]
    Vv = sp.simplify((vals[x_] - vals[z_]) / (vals[x_] + 1)); Wv = vals[x_] + 1
    back = {k: sp.simplify(e.subs({V: Vv, W: Wv, r: rv, s: sv, u: uv})) for k, e in sub.items()}
    okrt &= all(sp.simplify(back[k] - vals[k]) == 0 for k in sub)
check("chart round trip on 5 random complex spectra (r,s,u complex here: identity check)", okrt)

Q = V * W - W + 1
P = V * W**2 * u - V * W * u - V * r * s - W**2 * u + W * u
D = sp.factor(H[14] - H[16])
check("D = H14 - H16 formula (4)", sp.simplify(D - (-V * W * r + V * W * u + V * r * s + W * r - r) / Q) == 0, D)
check("H14 = (Wu+rs)/r", sp.simplify(H[14] - (W * u + r * s) / r) == 0)
check("H4 = r+1-r/W", sp.simplify(H[4] - (r + 1 - r / W)) == 0)
fF = sp.cancel(F[0] / E)
check("f = F0/E = -Q/V", sp.simplify(fF + Q / V) == 0, sp.factor(fF))
check("H4 - H17 = D + r - s - u", sp.simplify(H[4] - H[17] - (D + r - s - u)) == 0)
# denominators of all 24 H divide r W (W-1) Q P (radical check)
allowed = [sp.Poly(e, V, W, r, s, u) for e in (r, W, W - 1, Q, P)]
okden = True; dens = set()
for k in range(24):
    den = sp.denom(sp.cancel(sp.together(H[k])))
    for fac, mult in sp.factor_list(den)[1]:
        fp = sp.Poly(fac, V, W, r, s, u)
        if fp.total_degree() == 0:
            continue
        hit = any(sp.cancel(fac / g.as_expr()).is_number for g in allowed)
        dens.add(str(fac))
        okden &= hit
check("every denominator factor of the 24 H lies in {r, W, W-1, Q, P}", okden, sorted(dens))
# curve
num7, den7 = sp.fraction(sp.cancel(sp.together(H[7] - h)))
Cc = sp.Poly(sp.expand(num7), V, W)
info("curve C: degrees (V,W), number of terms in V,W,r,s,u,h", (Cc.degree(V), Cc.degree(W), len(sp.Poly(sp.expand(num7), V, W, r, s, u, h).terms())))
check("C has bidegree (2,3)", Cc.degree(V) == 2 and Cc.degree(W) == 3)
den7p = sp.Poly(den7, V, W)
check("denominator of H7 is linear in V", den7p.degree(V) == 1, sp.factor(den7))
Cexpr = Cc.as_expr()
Acoef = sp.factor(Cc.as_expr().coeff(V, 2))
W3 = sp.factor(sp.expand(Cexpr).coeff(W, 3))
# the printed normalisation: A(W) = W u (r-W)(s+1-W), [W^3]C = u (V-1)^2 ; allow a common scalar
lam = sp.cancel(Acoef / (W * u * (r - W) * (s + 1 - W)))
check("leading V^2 coefficient A(W) = W u (r-W)(s+1-W) up to a constant", lam.free_symbols <= set() or sp.Poly(lam, V, W).total_degree() == 0, (Acoef, lam))
check("[W^3]C = u (V-1)^2 with the same constant", sp.simplify(W3 - lam * u * (V - 1)**2) == 0, W3)
Cn = sp.expand(Cexpr / lam)
check("C (normalised) has 32 terms", len(sp.Poly(Cn, V, W, r, s, u, h).terms()) == 32, len(sp.Poly(Cn, V, W, r, s, u, h).terms()))
R1 = sp.factor(sp.resultant(Cn, sp.expand(P), V))
R2 = sp.factor(sp.resultant(Cn, sp.expand(Q), V))
info("Res_V(C,P)", R1); info("Res_V(C,Q)", R2)
t1 = sp.cancel(R1 / (r * s * (W - 1) * (r - u * W) * (u * W + r * s) * (s + u - u * W)))
t2 = sp.cancel(R2 / (W * h * (W - 1) * (u * W + r * s)))
check("Res_V(C,P) = const * rs(W-1)(r-uW)(uW+rs)(s+u-uW)", t1.is_number, t1)
check("Res_V(C,Q) = const * W h (W-1)(uW+rs)", t2.is_number, t2)
# V = 0 slice
C0 = sp.factor(Cn.subs(V, 0))
Q2 = u * W**2 + u * (h - r - s - 1) * W + r * (s + u)
check("C(V=0) = const (W-1) Q2(W)", sp.cancel(C0 / ((W - 1) * Q2)).is_number, C0)
# homogeneous inverse and E_hom
an, ad = sp.fraction(sp.cancel(as_)); tn, td = sp.fraction(sp.cancel(ts))
def br(Ai, Bj):  # [A_i + B_j] for spinors A=(num,den), B=(num,den)
    return Ai[0] * Bj[1] + Bj[0] * Ai[1]
Ah = [(an, ad), (zs, 1), (xs, 1), (ws, 1)]; Bh = [(tn, td), (0, 1), (1, 1), (1, 0)]
def Fh(p):
    e = 1
    for i in range(4):
        e *= br(Ah[i], Bh[p[i]])
    return sp.expand(e)
Ehom = sp.factor(Fh(L.PERMS[2]) - Fh(L.PERMS[0]))
lamE = sp.cancel(Ehom / (-r * V * W * (W - 1) * Q * P))
check("E_hom = const * (-r V W (W-1) Q P)", lamE.is_number, (Ehom, lamE))
info("reduced pairs a = an/ad, t = tn/td", (sp.factor(an), sp.factor(ad), sp.factor(tn), sp.factor(td)))
# specialisation (2,3,4,5)
spc = {r: 2, s: 3, u: 4, h: 5}
Cs = sp.Poly(Cn.subs(spc), V, W)
lead = sp.Poly(Cs.as_expr().coeff(V, 2), W); const = sp.Poly(Cs.as_expr().coeff(V, 0), W)
check("(2,3,4,5): leading and constant V coefficients coprime", sp.gcd(lead, const).degree() == 0)
check("(2,3,4,5): V=inf projections 0, r, s+1 distinct", len({0, 2, 4}) == 3 and sp.Poly(Acoef.subs(spc), W).degree() == 3)
q2s = sp.Poly(Q2.subs(spc), W)
check("(2,3,4,5): Q2 = 4W^2 - 4W + 14", sp.expand(q2s.as_expr() - (4 * W**2 - 4 * W + 14)) == 0, q2s.as_expr())
check("(2,3,4,5): Q2 squarefree and coprime to both resultants and to W(W-1)",
      sp.discriminant(q2s) != 0 and sp.gcd(q2s, sp.Poly(R1.subs(spc), W)).degree() == 0 and sp.gcd(q2s, sp.Poly(R2.subs(spc), W)).degree() == 0 and sp.gcd(q2s, sp.Poly(W * (W - 1), W)).degree() == 0)
# kernel of Hess C at kappa = singleton (assignment) span E: replayed in the library (3.11)
Hk = L.hess_kappa()
ker = L.nullspace(Hk)
Eperp_ok = all(v[0] == 0 and v[1] == 0 and v[5] == 0 and v[6] == 0 for v in ker) and len(ker) == 10
check("ker Hess C(kappa) = {m0=m1=m5=m6=0}, dim 10 (lib)", Eperp_ok)

# ---------------------------------------------------------------------------------------------
# NUMERIC falsification on the chart: random real (r,s,u,h), random complex W, solve C=0 for V,
# keep supported points (all 24 gaps >= 0), test the barriers and the Leibniz sum.
import numpy as np
Hf = sp.lambdify((V, W, r, s, u), H, "numpy")
Cf = sp.lambdify((W, r, s, u, h), [Cn.coeff(V, 2), Cn.coeff(V, 1), Cn.coeff(V, 0)], "numpy")
Pn = np.array(L.P, dtype=float)
rng = np.random.default_rng(20261125)
nsup = 0; worst = None; bar_bad = 0; supp_pts = []
for trial in range(40000):
    rv, sv, uv, hv = rng.normal(size=4) * 2
    Wv = complex(rng.normal() * 2, rng.normal() * 2)
    c2, c1, c0 = [complex(x) for x in Cf(Wv, rv, sv, uv, hv)]
    if abs(c2) < 1e-12:
        continue
    for Vv in np.roots([c2, c1, c0]):
        try:
            Hv = np.array(Hf(Vv, Wv, rv, sv, uv), dtype=complex)
        except ZeroDivisionError:
            continue
        if not np.all(np.isfinite(Hv)):
            continue
        for eta in (1.0, -1.0):
            g = -eta * Hv.imag
            if g.min() < -1e-9 * max(1, abs(g).max()):
                continue
            nsup += 1
            d = g[14] - g[16]
            if d > g[14] + 1e-9 or d > g[4] + 1e-9 or abs((g[4] - g[17]) - d) > 1e-7 * max(1, abs(g).max()):
                bar_bad += 1
            supp_pts.append(g)
info("NUMERIC: supported chart points found (40000 trials, seed 20261125)", nsup)
check("NUMERIC: barriers d <= g14, d <= g4 and g4-g17 = g14-g16 on all supported points", bar_bad == 0, bar_bad)
# Leibniz sums for a few supported points
mins = []
for g in supp_pts[:25]:
    c = np.linalg.lstsq(Pn, g, rcond=None)[0]
    val, _ = L.min_phi_over_SO4(c, restarts=8, seed=1)
    mins.append(val / max(1e-12, abs(g).max()))
if mins:
    info("NUMERIC: min over SO(4) of Phi/|g|_inf for first %d supported points (8 restarts, seed 1)" % len(mins), min(mins))
    check("NUMERIC: no violation (> -1e-8)", min(mins) > -1e-8)

print("SUMMARY: %d exact/numeric checks failed" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
