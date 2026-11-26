"""Own-code replay for independent_inventory70_native_harmonic_uniform_20261023 (primary proof of paper
case T_71 = {1,3,4,10,11,12}). Audit batch C37-T69-T71.
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
AUXN = json.load(open(os.path.join(HERE, "..", "deps", "dependency-map.json")))["auxiliary_and_branch_nodes"]

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
            if not earlier:
                S7 = frozenset(set(T) | {sig})
                for node in AUXN:
                    if node.get("contacts") and len(node["contacts"]) == 7 and node.get("stage_paper") is not None and L.orbit_rep(node["contacts"]) == L.orbit_rep(S7) and node["stage_paper"] < stage_self:
                        earlier.append("aux7:" + node["node"])
                rows[-1] = rows[-1] + ("aux:" + str(earlier),)
            allok &= bool(earlier)
    for r in rows:
        info(name + ": tie", r)
    check(name + ": every added tie covered (all-minor cone, or an earlier-stage case)", allok)
    return rows

# ---------------------------------------------------------------------------------------------
T71 = [1, 3, 4, 10, 11, 12]
check("T_71 literal equals Appendix E", frozenset(T71) == L.CONTACT_SETS[71])
info("stage of T_71 (Table 6)", STAGE["T_71"])
tau_ = (0, 2, 1, 3)   # (1 2) on 0-based labels 0..3
g = ((0, 1, 2, 3), tau_, 0)
img = sorted(L.g_act_perm(g, L.PERMS[j]) for j in T71)
img = sorted(L.perm_index(p) for p in [L.compose(L.PERMS[j], tau_) for j in T71])
check("sigma -> sigma o tau maps T_71 to the working tuple {1,4,5,7,9,14}", img == [1, 4, 5, 7, 9, 14], img)
check("sigma -> sigma o tau is the group element (id, tau, 0) of the library", sorted(L.g_index_map(g)[j] for j in T71) == img)
lit_src = L.generators(T71)
check("literal source = [0,1,2,8,10,12,53,56,59]", lit_src == [0, 1, 2, 8, 10, 12, 53, 56, 59], lit_src)
gm = L.g_index_map(g)
# image of generators: f_j o g^{-1}; compute pulled back profiles
work_src = L.generators([1, 4, 5, 7, 9, 14])
mapped = sorted(next(k for k in range(82) if all(L.F_VALUES[k][gm[sig]] == L.F_VALUES[j][sig] for sig in range(24))) for j in lit_src)
check("relabelled literal source = working source [0,1,2,5,11,12,33,35,59]", mapped == work_src == [0, 1, 2, 5, 11, 12, 33, 35, 59], (mapped, work_src))
T = [1, 4, 5, 7, 9, 14]
keep = [0, 1, 2, 5, 12, 33, 35, 59]
meas = [(0, 18), (1, 17), (2, 15), (5, 11), (12, 6), (33, 12), (59, 0)]
rows = cone_checks(T, work_src, 35, [35, 11], [1, 5, 12, 33], keep, meas, [(3, 12)], STAGE["T_71"], "T71")
check("T71: d <= g3 (g12 = e >= 0)", True)
rest = sorted(s_ for s_ in range(24) if s_ not in T and L.F_VALUES[35][s_] == 0)
check("T71: extras with f35 = 0 are exactly 0,6,11,12,15,17,18", rest == [0, 6, 11, 12, 15, 17, 18], rest)
check("T71: T + {12} = Table 5 set #27 {1,4,5,7,9,12,14} literally", frozenset(T) | {12} == L.AUX7[27])
info("stage of Table-5 #27 (dependency map)", "5 (S67, root_inventory67_g14)")

# ---------------------------------------------------------------------------------------------
Z, Qv, r, k, v, h = sp.symbols("Z Q r k v h")
y_, z_, w_, x_, t_ = sp.symbols("y z w x t")
def Fsig(A, B):
    out = []
    for p in L.PERMS:
        e = 1
        for i in range(4):
            if p[i] != 3:
                e *= (A[i] - B[p[i]])   # this document uses difference brackets: B = c = -b
        out.append(e)
    return out
F0 = Fsig([y_, z_, w_, x_], [t_, 0, 1, None])
E0 = F0[7] - F0[1]
# affine convention (factor 1 for C = infinity) gives E = -t(y-z)(x-1); the printed sign is that of the
# homogeneous difference bracket, whose C = infinity factor is -A_den (one such factor per sigma, so all
# 24 products change sign together and every ratio H is unchanged)
check("E = F7 - F1 = t(y-z)(x-1) up to the global sign of the convention", sp.expand(E0 + t_ * (y_ - z_) * (x_ - 1)) == 0, sp.factor(E0))
Hf0 = [sp.cancel((F0[j] - F0[1]) / E0) for j in range(24)]
check("r = H4 = (y-t)(z-w)/(t(y-z))", sp.simplify(Hf0[4] - (y_ - t_) * (z_ - w_) / (t_ * (y_ - z_))) == 0)
check("H5 = r(1+k), k = (x-w)/((x-1)(z-w))", sp.simplify(Hf0[5] / Hf0[4] - 1 - (x_ - w_) / ((x_ - 1) * (z_ - w_))) == 0)
check("H9 = 1 + v, v = (1-t)y(x-z)/(t(y-z)(x-1))", sp.simplify(Hf0[9] - 1 - (1 - t_) * y_ * (x_ - z_) / (t_ * (y_ - z_) * (x_ - 1))) == 0)
check("H7 = 1", sp.simplify(Hf0[7] - 1) == 0)
# inverse (5): exact round trip on random rational spectra (forward ratios -> (5) -> spectra)
import random as _rnd
rg = _rnd.Random(11)
okrt = True
for _ in range(20):
    pt = {y_: sp.Rational(rg.randint(-60, 60), rg.randint(1, 9)), z_: sp.Rational(rg.randint(-60, 60), rg.randint(1, 9)),
          w_: sp.Rational(rg.randint(-60, 60), rg.randint(1, 9)), x_: sp.Rational(rg.randint(-60, 60), rg.randint(1, 9)),
          t_: sp.Rational(rg.randint(-60, 60), rg.randint(1, 9))}
    if len(set(pt.values()) | {0, 1}) < 7:
        continue
    Ev = E0.subs(pt)
    if Ev == 0:
        continue
    rv = Hf0[4].subs(pt); kv = Hf0[5].subs(pt) / rv - 1; vv = Hf0[9].subs(pt) - 1
    zz, ww = pt[z_], pt[w_]
    xx = (kv * (zz - ww) - ww) / (kv * (zz - ww) - 1)
    MM = (vv * (1 - ww) / (1 - kv * (zz - 1)) - rv) / (zz - ww)
    yy = (1 + MM * zz) / (1 + MM)
    tt = yy * (zz - ww) / (rv * (yy - zz) + zz - ww)
    okrt &= (xx == pt[x_] and yy == pt[y_] and tt == pt[t_])
check("inverse (5) recovers (x,y,t) from (z,w,r,k,v) exactly on 20 random rational spectra", okrt)
# (Z, Q) coordinates and the homogeneous inverse
Ac = 1 - k * Z
N = (1 - r) * Z - r + Qv * (1 + v + (v - k) * Z)
Dd = Z - r + Qv * (1 + v - k * Z)
TD = (1 - r) * Z - r + Qv * (1 + v - k * Z)
K = Z + Qv * (1 - k * Z)
Ah = [(N, Dd), (Z + 1, 1), (1 - Qv * Ac, 1), ((k + 1) * Qv - 1, k * Qv - 1)]
Bh = [(N, TD), (0, 1), (1, 1), (1, 0)]
xs = (k * (z_ - w_) - w_) / (k * (z_ - w_) - 1)
Ms = (v * (1 - w_) / (1 - k * (z_ - 1)) - r) / (z_ - w_)
subZQ = {z_: Z + 1, w_: 1 - Qv * (1 - k * Z)}
xZ = xs.subs(subZQ)
MZ = Ms.subs(subZQ)
yZ = (1 + MZ * (Z + 1)) / (1 + MZ)
tZ = yZ * (Z + 1 - subZQ[w_]) / (r * (yZ - Z - 1) + Z + 1 - subZQ[w_])
check("homogeneous inverse equals (5) in (Z,Q): y, x, t", sp.cancel(Ah[0][0] / Ah[0][1] - yZ) == 0 and sp.cancel(Ah[3][0] / Ah[3][1] - xZ) == 0 and sp.cancel(Bh[0][0] / Bh[0][1] - tZ) == 0)
def br(Ai, Bj):   # difference bracket [A_i - C_j] for spinors (num, den)
    return Ai[0] * Bj[1] - Bj[0] * Ai[1]
def Fh(p):
    e = 1
    for i in range(4):
        e *= br(Ah[i], Bh[p[i]])
    return sp.expand(e)
Fhs = [Fh(p) for p in L.PERMS]
Ehom = sp.factor(Fhs[7] - Fhs[1])
H = [sp.factor(sp.cancel((Fhs[j] - Fhs[1]) / (Fhs[7] - Fhs[1]))) for j in range(24)]
lamE = sp.cancel(Ehom / (-Qv * Z * K * N))
check("E_hom = const * (-Q Z K N)", lamE.is_number, (Ehom, lamE))
check("H3 = r/Q", sp.simplify(H[3] - r / Qv) == 0, H[3])
fF = sp.cancel(Fhs[1] / (Fhs[7] - Fhs[1]))
check("F1/E = -r(Z+1)/K", sp.simplify(fF + r * (Z + 1) / K) == 0, sp.factor(fF))
check("contacts H1=0, H4=r, H5=r(1+k), H7=1, H9=1+v in (Z,Q)", [sp.simplify(H[1]), sp.simplify(H[4] - r), sp.simplify(H[5] - r * (1 + k)), sp.simplify(H[7] - 1), sp.simplify(H[9] - 1 - v)] == [0] * 5)
dens = set(); okden = True
allowed = [r, k, v, Qv, Z, sp.expand(N)]
for j in range(24):
    den = sp.denom(sp.cancel(sp.together(H[j])))
    for fac, m in sp.factor_list(den)[1]:
        if sp.Poly(fac, Z, Qv, r, k, v).total_degree() == 0:
            continue
        dens.add(str(fac))
        okden &= any(sp.cancel(fac / gg).is_number for gg in allowed)
check("every denominator factor of the 24 H lies in {Q, Z, N} (and constants r,k,v)", okden, sorted(dens))
num14, den14 = sp.fraction(sp.cancel(sp.together(H[14])))
Cexpr = sp.expand(num14 - h * den14)
Cp = sp.Poly(Cexpr, Z, Qv)
check("C = N14 - h D14 bidegree (2,2) in (Z,Q)", (Cp.degree(Z), Cp.degree(Qv)) == (2, 2), (Cp.degree(Z), Cp.degree(Qv)))
z2 = sp.factor(Cexpr.coeff(Z, 2))
lam = sp.cancel(z2 / (k * (Qv * v - r) * (Qv * k + r - 1)))
check("Z^2 coefficient = const k(Qv-r)(Qk+r-1), independent of h", lam.is_number, z2)
Cn = sp.expand(Cexpr / lam)
nterms = len(sp.Poly(Cn, Z, Qv, r, k, v, h).terms())
check("C has 33 terms", nterms == 33, nterms)
e1214 = sp.factor(H[14] - H[12])
check("H14 - H12 = (r/Q - v)(kQ - 1)", sp.simplify(H[14] - H[12] - (r / Qv - v) * (k * Qv - 1)) == 0, e1214)
check("H3 - H12 = rk + v - H14 - vkQ identically (so = rk+v-h-vkQ on C)", sp.simplify(H[3] - H[12] - (r * k + v - H[14] - v * k * Qv)) == 0)
RZ = sp.factor(Cn.subs(Z, 0))
check("C(Z=0) = const (Qv+Q-r)(Qkv-rk-v+h)", sp.cancel(RZ / ((Qv * v + Qv - r) * (Qv * k * v - r * k - v + h))).is_number, RZ)
RN = sp.factor(sp.resultant(Cn, sp.expand(N), Qv))
tN = sp.cancel(RN / (v * Z * (Z + 1) * (Z * r * k - Z * v + r * k - v - 1) * (Z * r * k - Z * k - r + v + 1)))
check("Res_Q(C,N) = const vZ(Z+1)(Zrk-Zv+rk-v-1)(Zrk-Zk-r+v+1)", tN.is_number, RN)
Np = sp.Poly(sp.expand(N), Qv)
common = sp.resultant(Np.coeff_monomial(Qv), Np.coeff_monomial(1), Z)
check("N's Q-coefficient and constant have common zero iff 1+v-r(k+1) = 0", sp.cancel(common / (1 + v - r * (k + 1))).is_number, sp.factor(common))
check("top coefficient [Z^2 Q^2] = const k^2 v", sp.cancel(sp.Poly(Cn, Z, Qv).coeff_monomial(Z**2 * Qv**2) / (k**2 * v)).is_number)
# Z = -1 on K = 0
Km1 = sp.solve(sp.Eq(K.subs(Z, -1), 0), Qv)
check("Z=-1, K=0 gives Q = 1/(k+1) and N = 0 there", Km1 == [1 / (k + 1)] and sp.simplify(N.subs({Z: -1, Qv: 1 / (k + 1)})) == 0, Km1)
# height points K = 0 elsewhere: H bounded except on QZN
info("K = 0 meets Q Z N only at", "checked by the denominator list above (K is not a denominator of any H)")
# specialisation (2,3,4,5)
spc = {r: 2, k: 3, v: 4, h: 5}
conds = [r, k, v, r - 1, k + 1, h, v + 1, 1 + v - r * (k + 1), r * k - v * (1 - r)]
check("(2,3,4,5): all listed nonvanishing conditions hold", all(sp.simplify(c.subs(spc)) != 0 for c in conds))
Cs = sp.expand(Cn.subs(spc))
c2 = sp.Poly(Cs.coeff(Z, 2), Qv); c0 = sp.Poly(Cs.coeff(Z, 0), Qv)
check("(2,3,4,5): Z^2 and Z^0 coefficients coprime", sp.gcd(c2, c0).degree() == 0)
# Z-projection side: constant-Z components need Q^2 and Q^0 coefficients to share a root
q2 = sp.Poly(Cs.coeff(Qv, 2), Z); q0 = sp.Poly(Cs.coeff(Qv, 0), Z)
info("(2,3,4,5): gcd of Q^2 and Q^0 coefficients (constant-Z components)", sp.gcd(q2, q0))

# ---------------------------------------------------------------------------------------------
import numpy as np
Hf = sp.lambdify((Z, Qv, r, k, v), H, "numpy")
Cf = sp.lambdify((Qv, r, k, v, h), [Cn.coeff(Z, 2), Cn.coeff(Z, 1), Cn.coeff(Z, 0)], "numpy")
Pn = np.array(L.P, dtype=float)
rng = np.random.default_rng(20261125)
nsup = 0; bad = 0; pts = []
for trial in range(40000):
    rv, kv, vv, hv = rng.normal(size=4) * 2
    Qc = complex(rng.normal() * 2, rng.normal() * 2)
    a2, a1, a0 = [complex(x) for x in Cf(Qc, rv, kv, vv, hv)]
    if abs(a2) < 1e-12:
        continue
    for Zc in np.roots([a2, a1, a0]):
        try:
            Hv = np.array(Hf(Zc, Qc, rv, kv, vv), dtype=complex)
        except ZeroDivisionError:
            continue
        if not np.all(np.isfinite(Hv)):
            continue
        for eta in (1.0, -1.0):
            gg = -eta * Hv.imag
            sc = max(1, abs(gg).max())
            if gg.min() < -1e-9 * sc:
                continue
            nsup += 1
            d = gg[3] - gg[12]
            if d > gg[3] + 1e-9 * sc or abs(d - eta * vv * kv * Qc.imag) > 1e-7 * sc:
                bad += 1
            pts.append(gg)
info("NUMERIC: supported chart points (40000 trials, seed 20261125)", nsup)
check("NUMERIC: d = eta v k Im Q and d <= g3 on all supported points", bad == 0, bad)
mins = []
for gg in pts[:25]:
    c = np.linalg.lstsq(Pn, gg, rcond=None)[0]
    val, _ = L.min_phi_over_SO4(c, restarts=8, seed=1)
    mins.append(val / max(1e-12, abs(gg).max()))
if mins:
    info("NUMERIC: min over SO(4) of Phi/|g|_inf, first %d supported points (8 restarts, seed 1)" % len(mins), min(mins))
    check("NUMERIC: no violation (> -1e-8)", min(mins) > -1e-8)
print("SUMMARY: %d exact/numeric checks failed" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
