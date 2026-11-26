"""Second-pass independent replay for the supplement document root_moc4_T16_uniform_20261022
(primary proof of paper case T_17 = {0,1,4,5,9,12}).

Everything is rebuilt from definitions: spectra -> sigma-points -> gaps -> generators.
Nothing from the supplement is executed or imported.  The shared library moc4lib (own code
of this audit) is imported by path.  Prints one line per check; exits 1 if an EXACT check fails.

Usage: python -I -B root_moc4_T16_uniform_20261022__check.py [--numeric]
"""
import os
import sys
import itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)


T = frozenset({0, 1, 4, 5, 9, 12})

# ------------------------------------------------------------------------------------------
# 1. Identification with the paper case and the generator list
# ------------------------------------------------------------------------------------------
check("doc contact set equals paper T_17 literally", L.CONTACT_SETS[17] == T)
gens = L.generators(T)
check("generators(T) == {0,1,2,5,7,9,12,35}", gens == [0, 1, 2, 5, 7, 9, 12, 35], gens)
check("only f35 is non-minor among them",
      [j for j in gens if not L.F_IS_MINOR[j]] == [35], [L.f_name(j) for j in gens])
mat = [list(L.F_VALUES[j]) for j in gens]
check("the eight generators are linearly independent (simplicial cone)", L.rank(mat) == 8)

# dimension of V_T = {y in im P : y_T = 0}
rowsT = [L.P[s] for s in sorted(T)]
ns = L.nullspace(rowsT)              # c in Q^14 with (Pc)_T = 0
check("dim V_T = 8", len(ns) == 8, len(ns))


def coeffs(y):
    g = lambda j: y[j]
    return {0: g(18), 1: g(6) - g(7) - g(8) + g(16), 2: -g(2) - g(7) + g(13),
            5: g(6) - g(7) - g(8) + g(10), 7: -g(2) - g(6) + g(8), 9: g(7), 12: g(6), 35: g(2)}


ok = True
for c in ns:
    y = L.values(c)
    cf = coeffs(y)
    rec = [sum(cf[j] * L.F_VALUES[j][s] for j in gens) for s in range(24)]
    ok = ok and all(Fraction(rec[s]) == Fraction(y[s]) for s in range(24))
check("source expansion g = g18 f0 + ... + g2 f35 holds on all of V_T (basis of 8 vectors)", ok)
# also each generator reproduces itself (coefficient functional = dual basis)
ok = True
for k in gens:
    cf = coeffs(list(L.F_VALUES[k]))
    ok = ok and all(cf[j] == (1 if j == k else 0) for j in gens)
check("coefficient functionals are the dual basis of the generators", ok)

z35 = sorted(L.F_ZERO[35] - T)
check("nontied zero indices of f35 are 6,7,11,14,15,17,18", z35 == [6, 7, 11, 14, 15, 17, 18], z35)
pos35 = sorted(set(range(24)) - L.F_ZERO[35])
ok = all(all(L.F_IS_MINOR[j] for j in L.generators(T | {j0})) for j0 in pos35)
check("for every j with f35(j)>0, generators(T u {j}) are all minor", ok, pos35)

# ------------------------------------------------------------------------------------------
# 2. Extra contacts at the seven zero indices of f35: which 6-subsets of T u {j} are
#    trivial (all-minor) orbits or earlier paper cases.  Stage numbers from Appendix G Table 6.
# ------------------------------------------------------------------------------------------
STAGE = {}
for k in (1, 2, 3, 4, 5, 6, 7, 15, 18, 19, 20, 21, 22, 23, 25, 27, 28, 29, 31, 34, 35, 36, 37,
          39, 40, 42, 43, 44, 72, 88):
    STAGE[k] = 1
for k in (10, 26, 46, 90):
    STAGE[k] = 2
for k in (9, 11, 12, 47):
    STAGE[k] = 3
for k in (8, 13, 17, 49):
    STAGE[k] = 4
for k in (14, 16, 24, 33, 38, 70):
    STAGE[k] = 5
REP2CASE = {L.orbit_rep(Tk): k for k, Tk in L.CONTACT_SETS.items()}
CITED = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16}     # T0..T11, T14, T15 (supplement labels)
all_ok_stage = True
all_ok_cited = True
for j in z35:
    U = T | {j}
    opts = []
    for S in itertools.combinations(sorted(U), 6):
        if j not in S:
            continue
        S = frozenset(S)
        if all(L.F_IS_MINOR[q] for q in L.generators(S)):
            opts.append(("minor-only", sorted(S), None))
        else:
            k = REP2CASE.get(L.orbit_rep(S))
            opts.append(("T_%s" % k, sorted(S), k))
    early = [o for o in opts if o[2] is None or STAGE.get(o[2], 99) < 4]
    cited = [o for o in opts if o[2] is None or o[2] in CITED]
    print("INFO extra contact j=%d: %s" % (j, "; ".join("%s %s" % (o[0], o[1]) for o in opts)))
    if early:
        o = early[0]
        if o[2] is not None:
            Tk = L.CONTACT_SETS[o[2]]
            g = L.find_g(Tk, frozenset(o[1]))
            ok = g is not None and L.g_act_set(g, Tk) == frozenset(o[1]) and frozenset(o[1]) <= U
            check("j=%d: group element maps T_%d into T u {j} (g=%s)" % (j, o[2], g), ok)
        else:
            ok = True
            check("j=%d: 6-subset %s of T u {j} is an all-minor face" % (j, o[1]), ok)
    all_ok_stage = all_ok_stage and bool(early)
    all_ok_cited = all_ok_cited and bool(cited)
check("every extra contact j in the f35 zero set is covered by a case of stage < 4 or a minor-only face",
      all_ok_stage)
check("every extra contact is covered by one of the cited T_1..T_12, T_15, T_16 or a minor-only face",
      all_ok_cited)

# real-v (common circle) case supplied by the auditor: contacts become {0..5} u T
U = T | {0, 1, 2, 3, 4, 5}
check("generators(T u {0..5}) are all minor (real-v stratum closes)",
      all(L.F_IS_MINOR[q] for q in L.generators(U)), L.generators(U))
check("perms 0..5 are exactly those with sigma(1)=1",
      [s for s in range(24) if L.PERMS[s][0] == 0] == [0, 1, 2, 3, 4, 5])

# ------------------------------------------------------------------------------------------
# 3. The chart A = (z,1,x,0), C = (w,v,oo,tv): exact symbolic identities
# ------------------------------------------------------------------------------------------
z, w, v, x, t, h, r, u, mu = sp.symbols("z w v x t h r u mu")
p, q = sp.symbols("p q", real=True)
Apts = [(z, 1), (1, 1), (x, 1), (0, 1)]


def Cpts(tt, tv=None):
    return [(w, 1), (v, 1), (1, 0), ((tt * v if tv is None else tv), 1)]


def br(P1, P2):
    return P1[0] * P2[1] - P1[1] * P2[0]


def Fs(tt, tv=None):
    C = Cpts(tt, tv)
    return [sp.expand(sp.Mul(*[br(Apts[i], C[s[i]]) for i in range(4)])) for s in L.PERMS]


F = Fs(t)
R = [sp.cancel((F[j] - F[0]) / (v * (z - w))) for j in range(24)]
check("R1 = x(v-1)/v", sp.simplify(R[1] - x * (v - 1) / v) == 0)
check("R4 = R1 + (t-1)(x-1)", sp.simplify(R[4] - R[1] - (t - 1) * (x - 1)) == 0)
check("R5 = 1 - t", sp.simplify(R[5] - (1 - t)) == 0)
check("(R4-R1)/R5 = 1-x", sp.simplify((R[4] - R[1]) / R[5] - (1 - x)) == 0)
hh = v * (1 - t) / (v - 1)
check("R1/R5 = x/h with h = v(1-t)/(v-1)", sp.simplify(R[1] / R[5] - x / hh) == 0)

tsub = 1 - h * (v - 1) / v
check("t = 1 - h(v-1)/v inverts h", sp.simplify(hh.subs(t, tsub) - h) == 0)
check("tv = (1-h)v + h", sp.expand(tsub * v - ((1 - h) * v + h)) == 0)
check("tv = 1 at h = 1 (cross collision c4 = a2)", sp.simplify((tsub * v).subs(h, 1) - 1) == 0)
Fh = Fs(None, (1 - h) * v + h)   # c4 = tv = (1-h)v + h, polynomial in v
den = v * (z - w) * (1 - tsub)
check("v(z-w)(1-t) = h(v-1)(z-w)", sp.simplify(den - h * (v - 1) * (z - w)) == 0)
H = [sp.cancel((Fh[j] - Fh[0]) / den) for j in range(24)]
check("H1 = x/h", sp.simplify(H[1] - x / h) == 0)
check("H4 = x/h - x + 1", sp.simplify(H[4] - (x / h - x + 1)) == 0)
check("H5 = 1", sp.simplify(H[5] - 1) == 0)
H2claim = (x - 1) * (v / (h * (v - 1)) - 1)
check("H2 = (x-1){v/(h(v-1)) - 1}", sp.simplify(H[2] - H2claim) == 0)
# d = -mu Im H2 with v = p + i q, x, h, mu real
H2pq = H2claim.subs(v, p + sp.I * q)
xr, hr, mur = sp.symbols("xr hr mur", real=True)
H2pq = H2pq.subs({x: xr, h: hr})
dexpr = -mur * sp.im(sp.expand_complex(H2pq))
dclaim = mur * (xr - 1) * q / (hr * ((p - 1) ** 2 + q ** 2))
check("d = -mu Im H2 = mu(x-1)Im v/(h|v-1|^2)", sp.simplify(dexpr - dclaim) == 0)
F0claim = ((h - 1) * v - h) / h
check("F0/(v(z-w)(1-t)) = ((h-1)v-h)/h", sp.simplify(sp.cancel(Fh[0] / den) - F0claim) == 0)

# E9, E12
a_ = tsub * (v - 1) - u * (1 - tsub)
D_ = tsub * (x - 1) - u * (1 - tsub)
b_ = r * (1 - tsub)
E9 = w * (v * (b_ + tsub - x) + z * (x - tsub * v)) + z * v * (-b_ + tsub * v - tsub)
E12 = a_ * z - D_ * w - tsub * (v - x)
num9 = sp.numer(sp.together(H[9] - r))
num12 = sp.numer(sp.together(H[12] - u))
rat9 = sp.factor(sp.cancel(num9 / E9))
rat12 = sp.factor(sp.cancel(num12 / E12))
check("H9 - r = E9 * nonzero factor", not rat9.free_symbols & {z, w, r, u}, "ratio=" + str(sp.factor(sp.cancel((H[9] - r) / E9))))
check("H12 - u = E12 * nonzero factor", not rat12.free_symbols & {z, w, r, u}, "ratio=" + str(sp.factor(sp.cancel((H[12] - u) / E12))))
check("a = -(v-1)(hu+hv-h-v)/v", sp.simplify(a_ + (v - 1) * (h * u + h * v - h - v) / v) == 0)
check("D = (1-h)(x-1) - hu + h(x-1+u)/v",
      sp.simplify(D_ - ((1 - h) * (x - 1) - h * u + h * (x - 1 + u) / v)) == 0)
check("x - tv = x - h - (1-h)v", sp.expand(x - tsub * v - (x - h - (1 - h) * v)) == 0)
# hu+hv-h-v as linear factor in v: coefficient (h-1) nonzero iff h != 1; zero v = h(1-u)/(h-1) real
# curve (3)
Nz = D_ * w + tsub * (v - x)
Q3 = (D_ * (x - tsub * v) * w ** 2
      + (a_ * v * (b_ + tsub - x) + tsub * (v - x) * (x - tsub * v) + D_ * v * (-b_ + tsub * v - tsub)) * w
      + tsub * (v - x) * v * (-b_ + tsub * v - tsub))
check("a * E9(z = N/a) equals the quadratic (3)", sp.simplify(a_ * E9.subs(z, Nz / a_) - Q3) == 0)
# D identically zero analysis: D = c0 + c1/v with c1 = h(x-1+u)
check("D(u=1-x) = x-1 (so D not identically zero)", sp.simplify(D_.subs(u, 1 - x) - (x - 1)) == 0)
# z = w divisor
check("E12|z=w = t(v-x)(z-1)", sp.simplify(E12.subs(w, z) - tsub * (v - x) * (z - 1)) == 0)
check("E9|z=w = z(v-z)(tv-x)", sp.simplify(E9.subs(w, z) - z * (v - z) * (tsub * v - x)) == 0)
check("t = 0 only at real v = -h/(1-h)", sp.simplify(tsub.subs(v, -h / (1 - h))) == 0)

# fiber at v = 1
Q3v1 = sp.factor(sp.simplify(Q3.subs(v, 1)))
check("curve at v=1 is (x-1)^2 w(w-1)", sp.simplify(Q3v1 - (x - 1) ** 2 * w * (w - 1)) == 0, Q3v1)
check("a at v=1 has derivative a1 = 1 - hu", sp.simplify(sp.diff(a_, v).subs(v, 1) - (1 - h * u)) == 0)
check("D at v=1 equals x-1", sp.simplify(D_.subs(v, 1) - (x - 1)) == 0)

# ------------------------------------------------------------------------------------------
# 4. Puiseux/Taylor branches at v = 1 (s = v - 1) and the residues of all 24 H_j
# ------------------------------------------------------------------------------------------
s = sp.symbols("s")
NORD = 4
Qs = sp.Poly(sp.numer(sp.together(Q3.subs(v, 1 + s))), s, w)
dQ = sp.Poly(sp.diff(Qs.as_expr(), w), s, w)


def trunc(e, n):
    e = sp.expand(e)
    return sum(sp.cancel(e.coeff(s, k)) * s ** k for k in range(n + 1))


def branch(w0):
    ws = sp.Integer(w0)
    d0 = sp.cancel(dQ.as_expr().subs({s: 0, w: w0}))
    print("INFO branch w(1)=%d: dQ/dw at the point = %s (series denominators are powers of it)" % (w0, sp.factor(d0)))
    assert d0 != 0
    for k in range(1, NORD + 1):
        val = sp.expand(Qs.as_expr().subs(w, ws))
        ck = sp.cancel(-val.coeff(s, k) / d0)
        ws = ws + ck * s ** k
    val = sp.expand(Qs.as_expr().subs(w, ws))
    ok = all(sp.cancel(val.coeff(s, k)) == 0 for k in range(NORD + 1))
    return ws, ok


DENOMS = set()


def laurent(expr_num, expr_den, nmax=1):
    """Laurent coefficients of num/den in s, num and den given as expressions analytic in s."""
    num = sp.expand(expr_num)
    dd = sp.expand(expr_den)
    m = 0
    while sp.cancel(dd.coeff(s, m)) == 0:
        m += 1
    # den = s^m * (d_m + d_{m+1} s + ...)
    dc = [sp.cancel(dd.coeff(s, m + k)) for k in range(NORD)]
    DENOMS.add(sp.factor(dc[0]))
    nc = [sp.cancel(num.coeff(s, k)) for k in range(NORD)]
    # series quotient
    qc = []
    for k in range(NORD):
        acc = nc[k] - sum(qc[i] * dc[k - i] for i in range(k))
        qc.append(sp.cancel(acc / dc[0]))
    out = {k - m: qc[k] for k in range(NORD) if k - m <= nmax}
    # drop vanishing leading coefficients (numerator and denominator may share a zero)
    for key in sorted(out):
        if key < 0 and sp.cancel(out[key]) == 0:
            del out[key]
        else:
            break
    return out


vS = 1 + s
tS = tsub.subs(v, vS)
aS = sp.cancel(a_.subs(v, vS))
DS = sp.cancel(D_.subs(v, vS))
FhS = [sp.expand(f.subs(v, vS)) for f in Fh]
results = {}
for w0 in (0, 1):
    ws, ok = branch(w0)
    check("branch w(1)=%d solves the curve to order s^%d" % (w0, NORD), ok)
    if w0 == 1:
        c1 = sp.cancel(sp.expand(ws).coeff(s, 1))
        c2 = sp.cancel(sp.expand(ws).coeff(s, 2))
        check("w=1 branch: w = 1 + (1-hu)/(x-1) s^2 + O(s^3)", c1 == 0 and sp.simplify(c2 - (1 - h * u) / (x - 1)) == 0, (c1, c2))
    # z = Nz/a ; multiply num and den by v^2 to make them polynomial in s (v = 1 + s)
    NzS = sp.cancel((Nz.subs(v, vS)).subs(w, ws) * vS ** 2)
    aSv = sp.cancel(aS * vS ** 2)
    NzS = trunc(sp.expand(sp.series(NzS, s, 0, NORD + 1).removeO()), NORD)
    aSv = trunc(sp.expand(aSv), NORD + 2)
    zl = laurent(NzS, aSv, nmax=2)
    if w0 == 0:
        check("w=0 branch: z has a simple pole at s=0", min(zl) == -1 and sp.cancel(zl[-1]) != 0, zl.get(-1))
    else:
        check("w=1 branch: z = 1 + s + O(s^2)", min(zl) == 0 and sp.cancel(zl[0] - 1) == 0 and sp.cancel(zl[1] - 1) == 0,
              (zl.get(0), zl.get(1)))
    # H_j = (F_j - F_0)/(h s (z - w)).  Write z = Nz/a; then (F_j - F_0) is affine in z:
    #   F_j - F_0 = P_j z + Q_j  =>  H_j = (P_j Nz + Q_j a)/(h s (Nz - w a))
    res = {}
    poles_ok = True
    for j in range(24):
        diff = sp.expand(FhS[j] - FhS[0])
        Pj = diff.coeff(z, 1)
        Qj = diff.coeff(z, 0)
        assert sp.expand(diff - Pj * z - Qj) == 0
        num = sp.expand((Pj.subs(w, ws)) * NzS + (Qj.subs(w, ws)) * aSv)
        num = trunc(num, NORD)
        denm = trunc(sp.expand(h * s * (NzS - ws * aSv)), NORD + 1)
        if j == 0:
            continue
        lj = laurent(num, denm, nmax=0)
        if min(lj) < -1:
            poles_ok = False
        res[j] = sp.factor(lj.get(-1, 0))
    check("w(1)=%d branch: every H_j has at most a simple pole" % w0, poles_ok)
    nz = {j: rj for j, rj in res.items() if rj != 0}
    print("INFO w(1)=%d branch: nonzero residues: %s" % (w0, nz))
    results[w0] = res
print("INFO leading denominator coefficients used in Laurent divisions: %s" % sorted(map(str, DENOMS)))
check("w=0 branch: res H2 = (x-1)/h", sp.simplify(results[0][2] - (x - 1) / h) == 0)
check("w=0 branch: res H7 = -(x-1)/h", sp.simplify(results[0][7] + (x - 1) / h) == 0)
check("w=1 branch: res H2 = (x-1)/h", sp.simplify(results[1][2] - (x - 1) / h) == 0)
check("w=1 branch: res H14 = -(x-1)/h", sp.simplify(results[1][14] + (x - 1) / h) == 0)
check("indices 7 and 14 are not contacts and are nontied zero indices of f35",
      7 not in T and 14 not in T and 7 in z35 and 14 in z35)

# ------------------------------------------------------------------------------------------
# 5. Real-v stratum: F_sigma with sigma(1)=1 equals (z-w) times a polynomial free of z, w
# ------------------------------------------------------------------------------------------
ok = True
for j in range(6):
    quo, rem = sp.div(sp.Poly(Fh[j], z, w), sp.Poly(z - w, z, w))
    ok = ok and rem.is_zero and not (quo.as_expr().free_symbols & {z, w})
check("F_0..F_5 = (z-w) * (polynomial in x, v, h only)", ok)
check("F1 - F0 nonzero factor: (F1-F0)/(z-w) = -x(1-v)... ",
      sp.simplify((Fh[1] - Fh[0]) / (z - w) - x * (v - 1)) == 0, sp.factor((Fh[1] - Fh[0]) / (z - w)))

# ------------------------------------------------------------------------------------------
# 6. NUMERIC falsification (only with --numeric)
# ------------------------------------------------------------------------------------------
if "--numeric" in sys.argv:
    import numpy as np
    Pn = L.P_numpy()
    worst = (np.inf, None)
    found = 0
    seeds_used = []
    for seed in range(int(os.environ.get("NSEEDS", "12"))):
        pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=30, margin=1e-3)
        if pt is None:
            continue
        found += 1
        seeds_used.append(seed)
        g = pt["gaps"]
        c = np.linalg.lstsq(Pn, g, rcond=None)[0]
        val, O = L.min_phi_over_SO4(c, restarts=30, seed=seed)
        scale = float(np.mean(g[[s for s in range(24) if s not in T]]))
        print("NUMERIC stratum seed=%d contact=%s min Phi=%.6e (mean off-T gap %.3f)" % (seed, pt["contact"], val, scale))
        if val < worst[0]:
            worst = (val, seed)
    print("NUMERIC summary: %d supported points with contact set exactly T found (seeds %s); worst min Phi = %.6e (seed %s)"
          % (found, seeds_used, worst[0], worst[1]))
    # chart-based sampling: random real x,h,r,u, nonreal v, mu = +-1; keep supported points
    rng = np.random.default_rng(2026)
    fF = sp.lambdify((z, w, v, x, h), Fh, "numpy")
    fQ = sp.lambdify((v, x, h, r, u), [sp.cancel(Q3.coeff(w, k)) for k in (2, 1, 0)], "numpy")
    fz = sp.lambdify((w, v, x, h, r, u), Nz / a_, "numpy")
    kept = 0
    worst2 = np.inf
    tried = 0
    maxdev = 0.0
    while tried < int(os.environ.get("NCHART", "40000")) and kept < 25:
        tried += 1
        xv, hv, rv, uv = rng.normal(size=4) * 2
        vv = complex(rng.normal() * 2, rng.normal() * 2)
        A2, A1, A0 = fQ(vv, xv, hv, rv, uv)
        for wv in np.roots([A2, A1, A0]):
            zv = fz(wv, vv, xv, hv, rv, uv)
            Fv = np.array(fF(zv, wv, vv, xv, hv), complex)
            tv_ = 1 - hv * (vv - 1) / vv
            Hv = (Fv - Fv[0]) / (vv * (zv - wv) * (1 - tv_))
            maxdev = max(maxdev, abs(Hv[9].imag), abs(Hv[12].imag), abs(Hv[1].imag), abs(Hv[4].imag), abs(Hv[5].imag))
            for sgn_mu in (1.0, -1.0):
                g = -sgn_mu * Hv.imag
                if np.min(g) > -1e-9:
                    kept += 1
                    c = np.linalg.lstsq(Pn, g, rcond=None)[0]
                    val, O = L.min_phi_over_SO4(c, restarts=20, seed=kept)
                    worst2 = min(worst2, val / max(1e-300, np.mean(g)))
    print("NUMERIC chart: max |Im H_T| on chart points = %.2e (contacts hold); %d supported points out of %d draws; worst normalized min Phi = %.6e"
          % (maxdev, kept, tried, worst2))

print("SUMMARY: %d EXACT check(s) failed" % len(FAIL))
sys.exit(1 if FAIL else 0)
