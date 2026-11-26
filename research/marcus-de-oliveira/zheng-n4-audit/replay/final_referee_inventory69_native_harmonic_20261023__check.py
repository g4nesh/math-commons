"""Own-code replay for final_referee_inventory69_native_harmonic_20261023 (primary proof of paper
case T_70 = {1,3,4,8,9,12}). Audit batch C37-T69-T71.
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
T = [1, 3, 4, 8, 9, 12]
check("T_70 literal equals Appendix E", frozenset(T) == L.CONTACT_SETS[70])
info("stage of T_70 (Table 6)", STAGE["T_70"])
src = [0, 1, 2, 5, 9, 53, 55, 56, 59]
keep = [0, 1, 2, 9, 53, 55, 56, 59]
meas = [(0, 20), (1, 17), (2, 15), (9, 10), (53, 14), (55, 11), (59, 2)]
rows = cone_checks(T, src, 56, [56, 5], [0, 9, 55], keep, meas, [(7, 10)], STAGE["T_70"], "T70")
# printed tie table (historical labels; paper index = historical + 1 except T4 -> T_4)
printed = {2: 13, 10: 42, 11: 29, 14: 44, 15: 4, 17: 44, 20: 29}
okp = True
for sig, nu in printed.items():
    S7 = set(T) | {sig}
    okp &= any(L.orbit_rep(S6) == L.orbit_rep(L.CONTACT_SETS[nu]) for S6 in
               [sorted(S7 - {x}) for x in T])
    okp &= STAGE["T_%d" % nu] < STAGE["T_70"]
check("T70: printed tie table (2,10,11,14,15,17,20 -> T_13,T_42,T_29,T_44,T_4,T_44,T_29) holds, all earlier stage", okp)
f56pos = sorted(s for s in range(24) if s not in T and L.F_VALUES[56][s] > 0)
check("T70: ties with f56 = 0 are exactly 2,10,11,14,15,17,20", sorted(set(range(24)) - set(T) - set(f56pos)) == [2, 10, 11, 14, 15, 17, 20])

# ---------------------------------------------------------------------------------------------
V, W, r, s, u, h, tau, Y = sp.symbols("V W r s u h tau Y")
a_, t_, z_, x_, w_ = sp.symbols("a t z x w")
def Fsig(A, B):
    out = []
    for p in L.PERMS:
        e = 1
        for i in range(4):
            if p[i] != 3:
                e *= (A[i] + B[p[i]])
        out.append(e)
    return out
F0 = Fsig([a_, z_, x_, w_], [t_, 0, 1, None])
check("E = F3 - F1 = (a+t)(w-z)", sp.expand(F0[3] - F0[1] - (a_ + t_) * (w_ - z_)) == 0)
check("F4-F1 = (a+t)(w+1)(x-z)", sp.expand(F0[4] - F0[1] - (a_ + t_) * (w_ + 1) * (x_ - z_)) == 0)
check("F9-F3 = t(a-w)(z+1)", sp.expand(F0[9] - F0[3] - t_ * (a_ - w_) * (z_ + 1)) == 0)
check("F9-F8 = a(w-x)(z+1)", sp.expand(F0[9] - F0[8] - a_ * (w_ - x_) * (z_ + 1)) == 0)
Lr = W - r
Pp = (W - 1) * (Lr - u * V) + s * (V - 1) * Lr
zs = W / V - 1; xs = r - 1 + Lr / V; ws = W - 1; as_ = Pp / (Lr - u * V); ts = Pp / (u * V)
sub = {a_: as_, t_: ts, z_: zs, x_: xs, w_: ws}
F = [sp.cancel(f.subs(sub)) for f in F0]
E = sp.cancel(F[3] - F[1])
H = [sp.factor(sp.cancel((F[k] - F[1]) / E)) for k in range(24)]
check("inverse (1): H1=0,H3=1,H4=r,H9=1+s,H8=1+s-u",
      [sp.simplify(H[1]), sp.simplify(H[3] - 1), sp.simplify(H[4] - r), sp.simplify(H[9] - 1 - s), sp.simplify(H[8] - 1 - s + u)] == [0] * 5)
# forward: uV/L = a/(a+t), L = (w+1)(w-x)/(w-z), P = t u V
Ef = F0[3] - F0[1]
rf = (F0[4] - F0[1]) / Ef; sf = (F0[9] - F0[3]) / Ef; uf = (F0[9] - F0[8]) / Ef
Vf = (w_ + 1) / (z_ + 1); Wf = w_ + 1
check("forward: L = W - r = (w+1)(w-x)/(w-z)", sp.simplify(Wf - rf - (w_ + 1) * (w_ - x_) / (w_ - z_)) == 0)
check("forward: uV/L = a/(a+t)", sp.simplify(uf * Vf / (Wf - rf) - a_ / (a_ + t_)) == 0)
check("forward: P = t u V", sp.simplify(Pp.subs({V: Vf, W: Wf, r: rf, s: sf, u: uf}) - t_ * uf * Vf) == 0)
check("round trip: inverse(forward) = identity", all(sp.simplify(e.subs({V: Vf, W: Wf, r: rf, s: sf, u: uf}) - k) == 0 for k, e in sub.items()))
Kw = -W * s + W * u + r * s
D = sp.factor(H[7] - H[10])
check("D = H7 - H10 = ruV/(r-W)", sp.simplify(D - r * u * V / (r - W)) == 0, D)
check("H10 = V(s-u) + W", sp.simplify(H[10] - V * (s - u) - W) == 0)
fF = sp.factor(sp.cancel(F[1] / E))
check("f = F1/E = (W-V)/(V-1)", sp.simplify(fF - (W - V) / (V - 1)) == 0, fF)
check("H6 - H11 = D - u and H18 - H20 = D + r", sp.simplify(H[6] - H[11] - D + u) == 0 and sp.simplify(H[18] - H[20] - D - r) == 0)
dens = set(); okden = True
allowed = [r, V, W, W - r, sp.expand(Pp)]
for k in range(24):
    den = sp.denom(sp.cancel(sp.together(H[k])))
    for fac, m in sp.factor_list(den)[1]:
        if sp.Poly(fac, V, W, r, s, u).total_degree() == 0:
            continue
        dens.add(str(fac))
        okden &= any(sp.cancel(fac / g).is_number for g in allowed)
check("every denominator factor of the 24 H lies in {r, V, W, W-r, P}", okden, sorted(dens))
num12, den12 = sp.fraction(sp.cancel(sp.together(H[12] - h)))
Cc = sp.Poly(sp.expand(num12), V, W)
check("C = num(H12-h) has bidegree (2,4)", (Cc.degree(V), Cc.degree(W)) == (2, 4), (Cc.degree(V), Cc.degree(W)))
Aco = sp.expand(Cc.as_expr()).coeff(V, 2)
lam = sp.cancel(Aco / ((r - 1) * Kw**2))
check("A(W) = (r-1) K(W)^2 up to a constant", lam.is_number, lam)
Cn = sp.expand(Cc.as_expr() / lam)
nterms = len(sp.Poly(Cn, V, W, r, s, u, h).terms())
check("C has 80 integer terms", nterms == 80, nterms)
Bco = Cn.coeff(V, 1); C0co = Cn.coeff(V, 0)
info("deg_W of A,B,C coefficients", (sp.degree(Aco, W), sp.degree(Bco, W), sp.degree(C0co, W)))
# nonconstancy, r = 1
H12r1 = sp.cancel(H[12].subs(r, 1))
dV = sp.factor(sp.diff(H12r1, V))
target = s * u * (W * (s - u + 1) - s - 1) / (W * (W - 1) * (V * (s - u) + W - s - 1)**2)
check("r = 1: dH12/dV printed formula", sp.simplify(dV - target) == 0, dV)
# homogeneous inverse
def br(Ai, Bj):
    return Ai[0] * Bj[1] + Bj[0] * Ai[1]
Ah = [(Pp, Lr - u * V), (W - V, V), ((r - 1) * V + Lr, V), (W - 1, 1)]
Bh = [(Pp, u * V), (0, 1), (1, 1), (1, 0)]
check("homogeneous inverse matches (1)", all(sp.simplify(Ah[i][0] / Ah[i][1] - [as_, zs, xs, ws][i]) == 0 for i in range(4)) and sp.simplify(Bh[0][0] / Bh[0][1] - ts) == 0)
def Fh(Ah, Bh, p):
    e = 1
    for i in range(4):
        e *= br(Ah[i], Bh[p[i]])
    return sp.expand(e)
Ehom = sp.factor(Fh(Ah, Bh, L.PERMS[3]) - Fh(Ah, Bh, L.PERMS[1]))
lamE = sp.cancel(Ehom / (V * W * Lr * (V - 1) * Pp))
check("E_hom = const * V W L (V-1) P", lamE.is_number, (Ehom, lamE))
R6 = sp.factor(sp.resultant(Cn, sp.expand(Pp), V))
info("Res_V(C,P)", R6)
t6 = sp.cancel(R6 / (-s * u * (r - W)**2 * (W - 1) * (r - s + u - 1) * ((-s + u - 1) * W + r * s + r)))
check("Res_V(C,P) printed (6) up to constant", t6.is_number, t6)
Pv = sp.Poly(sp.expand(Pp), V)
info("P = p1(W) V + p0(W): p1", sp.factor(Pv.coeff_monomial(V)))
Cw0 = sp.factor(Cn.subs(W, 0))
check("C(V, W=0) = const r^2 s (Vr-V-r)(Vs-s-1)", sp.cancel(Cw0 / (r**2 * s * (V * r - V - r) * (V * s - s - 1))).is_number, Cw0)
check("A(0) != 0 (as polynomial)", sp.expand(Aco.subs(W, 0)) != 0)
# node at W = r
Cwr = sp.factor(Cn.subs(W, r))
check("C(V, W=r) = const V^2 r^2 u^2 (r-1)", sp.cancel(Cwr / (V**2 * r**2 * u**2 * (r - 1))).is_number, Cwr)
Cnode = sp.expand(Cn.subs({V: tau * Y, W: r + tau}))
pn = sp.Poly(Cnode, tau)
low = min(m[0] for m in pn.monoms())
lowpoly = sp.factor(pn.coeff_monomial(tau**low))
info("node: lowest tau power and its coefficient (quadratic in Y)", (low, lowpoly))
check("node: lowest order is tau^2, a quadratic in Y with nonzero Y^2 and Y^0 coefficients",
      low == 2 and sp.degree(lowpoly, Y) == 2 and sp.expand(lowpoly.subs(Y, 0)) != 0)
# tangent cone: no terms of total order < 2 in (V, tau)
Ct = sp.Poly(sp.expand(Cn.subs(W, r + tau)), V, tau)
check("node: (0,r) is a double point (all terms of total order >= 2 in (V, W-r))", min(i + j for i, j in Ct.monoms()) == 2)
Astar = (r - 1) * (1 - u * Y) - s
Ahn = [(Astar, 1 - u * Y), (r, 0), (1 + (r - 1) * Y, Y), (r - 1, 1)]
Bhn = [(Astar, u * Y), (0, 1), (1, 1), (1, 0)]
Fn = [Fh(Ahn, Bhn, p) for p in L.PERMS]
En = sp.factor(Fn[3] - Fn[1])
check("node limit spinors: reference difference = -r A* Y", sp.expand(En + r * Astar * Y) == 0, En)
oklim = True
for k in range(24):
    Hk = sp.cancel(H[k].subs({V: tau * Y, W: r + tau}))
    nk, dk = sp.fraction(Hk)
    d0 = sp.expand(dk.subs(tau, 0))
    if d0 == 0:
        oklim = False; continue
    lim = sp.cancel(nk.subs(tau, 0) / d0)
    oklim &= sp.simplify(lim - (Fn[k] - Fn[1]) / En) == 0
check("node: all 24 H have finite limits along V = tau Y, W = r + tau, equal to the limit-spinor values (identically in Y)", oklim)
check("node limit is a shared-zero configuration: A1 = (r,0) and B3 = (1,0) give bracket 0", br(Ahn[1], Bhn[3]) == 0)
# infinity W -> inf with V = Y W: Newton polygon
CP = sp.Poly(Cn, V, W)
top = max(i + j for i, j in CP.monoms())
lead = sp.factor(sum(c * Y**i for (i, j), c in zip(CP.monoms(), CP.coeffs()) if i + j == top))
info("W=inf: top weighted degree and leading polynomial in Y", (top, lead))
target = (Y * (s - u) + 1) * ((r - 1) * Y * (s - u) - h + r + s - u)
check("W=inf: leading polynomial = const {Y(s-u)+1}{(r-1)Y(s-u)-h+r+s-u}", sp.cancel(lead / target).is_number)
check("W=inf: Newton polygon: max j at i=0,1,2 equal 4,3,2 (collinear, slope 1 only)",
      [max(j for (i, j) in CP.monoms() if i == ii) for ii in (0, 1, 2)] == [4, 3, 2])
# V = inf with finite W: A(W) = (r-1) K^2, root W* = rs/(s-u)
Ws = r * s / (s - u)
check("K(W*) = 0", sp.simplify(Kw.subs(W, Ws)) == 0)
G = sp.cancel(H[15] + H[7])
Gn, Gd = sp.fraction(G)
info("G* = H15 + H7: degrees in V of numerator/denominator", (sp.degree(Gn, V), sp.degree(Gd, V)))
# check readout: g15 + g7 = C2 + C9 + d with C2 = g15, C9 = g10, d = g7 - g10 (linear, trivial) and on kept basis
okG = all(L.F_VALUES[j][15] + L.F_VALUES[j][7] == (1 if j in (2, 9, 56) else 0) for j in keep)
check("G*: -eta Im G* = C2 + C9 + d on the kept basis", okG)
q, rem = sp.div(sp.Poly(Gn, V), sp.Poly(Gd, V))
astar = sp.factor(q.coeff_monomial(V) if q.degree() >= 1 else 0)
check("G* = a* V + b* + c*/(d1 V + d0) with a* = (W-1)K^2/(W(r-W)(K-u))",
      sp.degree(Gd, V) == 1 and q.degree() == 1 and sp.simplify(astar - (W - 1) * Kw**2 / (W * (r - W) * (Kw - u))) == 0, astar)
d1 = sp.Poly(Gd, V).coeff_monomial(V)
info("d1(W)", sp.factor(d1))
# height V = 1
C1 = sp.factor(Cn.subs(V, 1))
info("C(V=1, W)", C1)
Pv1 = sp.factor(Pp.subs(V, 1))
info("P at V=1 (pole divisor met at V=1)", Pv1)
# specialisation (2,3,5,7)
spc = {r: 2, s: 3, u: 5, h: 7}
check("(2,3,5,7): node equation 10(10Y^2+16Y+3) up to sign/scale", sp.cancel(lowpoly.subs(spc) / (10 * (10 * Y**2 + 16 * Y + 3))).is_number, sp.factor(lowpoly.subs(spc)))
check("(2,3,5,7): W* = -3", Ws.subs(spc) == -3)
check("(2,3,5,7): C(1,W) = -(W-1)(7W^3-50W^2+16W+12) up to scale", sp.cancel(C1.subs(spc) / (-(W - 1) * (7 * W**3 - 50 * W**2 + 16 * W + 12))).is_number, C1.subs(spc))
Cs = sp.Poly(Cn.subs(spc), V, W)
la = sp.Poly(sp.expand(Cs.as_expr()).coeff(V, 2), W); lc = sp.Poly(sp.expand(Cs.as_expr()).coeff(V, 0), W)
check("(2,3,5,7): leading and constant V-coefficients coprime", sp.gcd(la, lc).degree() == 0)
q3 = sp.Poly(sp.cancel(C1.subs(spc) / (W - 1)), W)
alld = sp.Poly(sp.expand(W * (W - r) * (W - 1) * R6 * Pv1).subs(spc), W)
check("(2,3,5,7): height cubic squarefree, disjoint from W(W-r)(W-1), Res_V(C,P) and P(V=1)", sp.discriminant(q3) != 0 and sp.gcd(q3, alld).degree() == 0)
check("(2,3,5,7): node quadratic disjoint from A* = 0 and has distinct nonzero roots",
      sp.discriminant(sp.Poly(lowpoly.subs(spc), Y)) != 0 and sp.gcd(sp.Poly(lowpoly.subs(spc), Y), sp.Poly(Astar.subs(spc) * Y, Y)).degree() == 0)
check("(2,3,5,7): d1(W*) != 0 and K(W*)-u != 0 and W* not in {0, r, 1}", sp.simplify(d1.subs(spc).subs(W, -3)) != 0 and Ws.subs(spc) not in (0, 2, 1))
check("(2,3,5,7): W=inf leading roots simple nonzero", sp.discriminant(sp.Poly(lead.subs(spc), Y)) != 0 and sp.expand(lead.subs(spc).subs(Y, 0)) != 0)

# ---------------------------------------------------------------------------------------------
import numpy as np
Hf = sp.lambdify((V, W, r, s, u), H, "numpy")
Cf = sp.lambdify((W, r, s, u, h), [Cn.coeff(V, 2), Cn.coeff(V, 1), Cn.coeff(V, 0)], "numpy")
Pn = np.array(L.P, dtype=float)
rng = np.random.default_rng(20261125)
nsup = 0; bad = 0; pts = []
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
            sc = max(1, abs(g).max())
            if g.min() < -1e-9 * sc:
                continue
            nsup += 1
            d = g[7] - g[10]
            if d > g[15] + g[7] + 1e-9 * sc:
                bad += 1
            dd = -eta * (rv * uv * Vv / (rv - Wv)).imag
            if abs(dd - d) > 1e-7 * sc:
                bad += 1
            pts.append(g)
info("NUMERIC: supported chart points (40000 trials, seed 20261125)", nsup)
check("NUMERIC: d = -eta Im(ruV/(r-W)) and d <= g15 + g7 on all supported points", bad == 0, bad)
mins = []
for g in pts[:25]:
    c = np.linalg.lstsq(Pn, g, rcond=None)[0]
    val, _ = L.min_phi_over_SO4(c, restarts=8, seed=1)
    mins.append(val / max(1e-12, abs(g).max()))
if mins:
    info("NUMERIC: min over SO(4) of Phi/|g|_inf, first %d supported points (8 restarts, seed 1)" % len(mins), min(mins))
    check("NUMERIC: no violation (> -1e-8)", min(mins) > -1e-8)
print("SUMMARY: %d exact/numeric checks failed" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
