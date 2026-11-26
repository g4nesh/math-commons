"""Second-pass replay for supplement document root_moc4_T20_uniform_20261022 (paper case T_21).

Own code only; uses the shared audit library ../lib/moc4lib.py. Nothing from the supplement
is executed or read. Prints one line per check; exits nonzero if an EXACT check fails.
"""
import os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

fails = 0
def report(name, ok, extra=""):
    global fails
    if not ok:
        fails += 1
    print(f"{'PASS' if ok else 'FAIL'} {name} {extra}".rstrip())

T = frozenset({0, 3, 5, 9, 11, 12})
# E1 contact set = paper T_21 (Appendix E), literally
report("E1 doc contact set equals paper T_21 exactly", L.CONTACT_SETS[21] == T,
       f"T_21={sorted(L.CONTACT_SETS[21])}")
print("INFO one-line notation of T:", [L.perm_to_oneline(L.PERMS[s]) for s in sorted(T)])

# E2 generator list
gens = L.generators(T)
report("E2 generators(T) = {0,1,2,7,9,12,29,39,40,41}", gens == [0, 1, 2, 7, 9, 12, 29, 39, 40, 41],
       f"got {gens}")
nonminor = [j for j in gens if not L.F_IS_MINOR[j]]
report("E2b only f40 is non-minor", nonminor == [40], f"nonminor={nonminor} names={[L.f_name(j) for j in gens]}")

# E3 replacement relations
def vadd(*terms):
    out = [Fraction(0)] * 14
    for coef, v in terms:
        out = [o + coef * x for o, x in zip(out, v)]
    return out
F = L.F_VECS
report("E3a f40+f1 = f0+f7+f39", vadd((1, F[40]), (1, F[1])) == vadd((1, F[0]), (1, F[7]), (1, F[39])))
report("E3b f40+f9 = f2+f7+f41", vadd((1, F[40]), (1, F[9])) == vadd((1, F[2]), (1, F[7]), (1, F[41])))
# span of generators: dimension and relation space
M = [F[j] for j in gens]
rk = L.rank(M)
report("E3c generators span an 8-dim space (10 gens, 2 relations)", rk == 8, f"rank={rk}")

# E4 detector functional delta = v8 - v2 - v6 - v16 on each generator
dval = {j: L.F_VALUES[j][8] - L.F_VALUES[j][2] - L.F_VALUES[j][6] - L.F_VALUES[j][16] for j in gens}
want = {j: 0 for j in gens}; want[40] = 1; want[1] = -1; want[9] = -1
report("E4 delta(f40)=1, delta(f1)=delta(f9)=-1, delta(others)=0", dval == want, f"{dval}")

# E5 sigma-point identities, affine chart F_sigma = prod (A_i - C_sigma(i))
A = sp.symbols("A0:4"); C = sp.symbols("C0:4")
def Fs(k):
    s = L.PERMS[k]
    return sp.Mul(*[A[i] - C[s[i]] for i in range(4)])
report("E5a F3*F11 = F5*F9", sp.expand(Fs(3) * Fs(11) - Fs(5) * Fs(9)) == 0)
report("E5b F3-F5 = (A0-C0)(A1-A2)(A3-C1)(C2-C3)",
       sp.expand(Fs(3) - Fs(5) - (A[0]-C[0])*(A[1]-A[2])*(A[3]-C[1])*(C[2]-C[3])) == 0)
report("E5c F3-F9 = (A0-A3)(A1-C2)(A2-C3)(C0-C1)",
       sp.expand(Fs(3) - Fs(9) - (A[0]-A[3])*(A[1]-C[2])*(A[2]-C[3])*(C[0]-C[1])) == 0)
# homogeneous (projective) version: brackets [x,y] = x1 y2 - x2 y1
X = [sp.symbols(f"x{i} y{i}") for i in range(8)]
def br(i, j):
    return X[i][0]*X[j][1] - X[i][1]*X[j][0]
def Fh(k):
    s = L.PERMS[k]
    return sp.Mul(*[br(i, 4 + s[i]) for i in range(4)])
report("E5d homogeneous: F3F11=F5F9, F3-F5=[0,4][1,2][3,5][6,7], F3-F9=[0,3][1,6][2,7][4,5]",
       sp.expand(Fh(3)*Fh(11) - Fh(5)*Fh(9)) == 0 and
       sp.expand(Fh(3) - Fh(5) - br(0,4)*br(1,2)*br(3,5)*br(6,7)) == 0 and
       sp.expand(Fh(3) - Fh(9) - br(0,3)*br(1,6)*br(2,7)*br(4,5)) == 0)
# lib convention check: (P m)_sigma = prod(a_i + b_sigma(i)) with c = -b, i.e. same products
a = [Fraction(2), Fraction(-3), Fraction(5), Fraction(7)]
b = [Fraction(1), Fraction(11), Fraction(-4), Fraction(6)]
m = L.matching_vector(a, b)
z = L.sigma_points(a, b)
report("E5e lib matching vector reproduces sigma-points (c=-b convention)", L.values(m) == z)

# E6 rectangle algebra: H(q3+q11)=H(q5+q9), q3 q11 = q5 q9 and H != 0 imply (q3-q5)(q3-q9)=0
H, q3, q5, q9, q11 = sp.symbols("H q3 q5 q9 q11", real=True)
prod_rel = sp.expand((H + sp.I*q3)*(H + sp.I*q11) - (H + sp.I*q5)*(H + sp.I*q9))
re_rel, im_rel = sp.re(prod_rel), sp.im(prod_rel)
ok6 = sp.simplify(im_rel - H*(q3 + q11 - q5 - q9)) == 0 and sp.simplify(re_rel - (q5*q9 - q3*q11)) == 0
# substitute q11 = q5 + q9 - q3 into real relation
ok6b = sp.expand((q5*q9 - q3*(q5 + q9 - q3)) - (q3 - q5)*(q3 - q9)) == 0
report("E6 rectangle: Im = H(q3+q11-q5-q9), Re = q5q9-q3q11, and with sums equal Re = (q3-q5)(q3-q9)", ok6 and ok6b)

# E7 polarization / Euler for the paper cubic C (A.2): grad C(y+i d) = 0 => C(y) = 0
ys = sp.symbols("y0:14", real=True); ds = sp.symbols("e0:14", real=True)
Cy = sp.expand(L.cubic_C(ys)); Cd = sp.expand(L.cubic_C(ds))
w = [ys[k] + sp.I*ds[k] for k in range(14)]
Cw = sp.expand(L.cubic_C(w))
ok7 = True
for k in range(14):
    gk = sp.expand(sp.diff(Cw, ds[k]) / sp.I)   # = dC/dm_k evaluated at y+i d
    gy = sp.diff(Cy, ys[k]); gd = sp.diff(Cd, ds[k])
    Hyd = sum(sp.diff(Cy, ys[k], ys[l]) * ds[l] for l in range(14))
    if sp.expand(sp.re(gk) - (gy - gd)) != 0 or sp.expand(sp.im(gk) - Hyd) != 0:
        ok7 = False
report("E7a grad C(y+i d) = [gradC(y)-gradC(d)] + i H_C(y) d (C has real coefficients)", ok7)
# y.gradC(y) = y.gradC(d) = (1/2) d^T H(d) y  (identity d^T H(d) = 2 gradC(d)^T for a cubic)
Hd = sp.Matrix(14, 14, lambda i, j: sp.diff(Cd, ds[i], ds[j]))
gradd = sp.Matrix([sp.diff(Cd, ds[k]) for k in range(14)])
report("E7b Euler for gradient: H_C(d) d = 2 grad C(d)", sp.expand(Hd*sp.Matrix(ds) - 2*gradd) == sp.zeros(14, 1))
report("E7c Euler: y.gradC(y) = 3 C(y)", sp.expand(sum(ys[k]*sp.diff(Cy, ys[k]) for k in range(14)) - 3*Cy) == 0)

# E8 cubic substitution identity (2)
q, aa, bb, cc, e, f, d, h = sp.symbols("q a b c e f d h")
coef = {0: q, 2: aa, 7: bb, 12: cc, 29: e, 39: f, 40: d, 41: h}
y0 = [sp.Integer(0)] * 14
for j, s in coef.items():
    y0 = [y0[k] + s * sp.Rational(F[j][k].numerator, F[j][k].denominator) for k in range(14)]
Cy0 = sp.expand(L.cubic_C(y0))
Mexp = (bb + d)*(q + aa + d) + cc*(aa + d)
rhs = (f + d + h)*Mexp + e*((q + aa + cc + d)*(f + d) + aa*cc)
diff = sp.expand(Cy0 - rhs)
report("E8 C(y0) = (f+d+h)M + e{(q+a+c+d)(f+d)+ac} with M=(b+d)(q+a+d)+c(a+d)", diff == 0,
       "" if diff == 0 else f"C(y0)={Cy0}; C+rhs={sp.expand(Cy0+rhs)}")
# P y0 = gap vector g, so y0 is the unique real lift
vals = [sum(sp.Integer(L.F_VALUES[j][s]) * coef[j] for j in coef) for s in range(24)]
Py0 = [sp.expand(sum(sp.Integer(L.P[s][k]) * y0[k] for k in range(14))) for s in range(24)]
report("E8b P y0 = sum coef_j (P f_j)", all(sp.expand(Py0[s] - vals[s]) == 0 for s in range(24)))
# positivity: every coefficient of rhs, as polynomial in the 8 nonnegative variables, is >= 0,
# and rhs - d^3 also has nonnegative coefficients
P_rhs = sp.Poly(sp.expand(rhs), q, aa, bb, cc, e, f, d, h)
P_rhs3 = sp.Poly(sp.expand(rhs - d**3), q, aa, bb, cc, e, f, d, h)
report("E8c rhs and rhs-d^3 have only nonnegative coefficients",
       all(c_ >= 0 for c_ in P_rhs.coeffs()) and all(c_ >= 0 for c_ in P_rhs3.coeffs()))

# E9 replacement procedure: lambda on 10 generators -> form (1)
lam = {j: sp.Symbol(f"l{j}", nonnegative=True) for j in gens}
lhs = [sum(lam[j]*F[j][k] for j in gens) for k in range(14)]
form = {0: lam[0] + lam[1], 2: lam[2] + lam[9], 7: lam[7] + lam[1] + lam[9], 12: lam[12],
        29: lam[29], 39: lam[39] + lam[1], 40: lam[40] - lam[1] - lam[9], 41: lam[41] + lam[9]}
rhs9 = [sum(form[j]*F[j][k] for j in form) for k in range(14)]
report("E9 sum lambda_j f_j = (l0+l1)f0+(l2+l9)f2+(l7+l1+l9)f7+l12 f12+l29 f29+(l39+l1)f39+(l40-l1-l9)f40+(l41+l9)f41",
       all(sp.expand(lhs[k] - rhs9[k]) == 0 for k in range(14)))

# ---------------- NUMERIC falsification ----------------
import numpy as np
from scipy.optimize import linprog
PN = L.P_numpy()
FV = np.array([[float(x) for x in L.F_VECS[j]] for j in range(82)])
minor_gens = [j for j in gens if L.F_IS_MINOR[j]]
nfound = 0; worstH = 0.0; worst_delta = -np.inf; worst_phi = np.inf; lp_fail = 0; min_sep = np.inf
n_distinct = 0; n_coll = 0; coll_H = 0.0; coll_d35 = np.inf
for seed in range(12):
    pt = L.random_point_on_contact_stratum(T, seed=seed, tries=30, tol=1e-10)
    if pt is None:
        continue
    nfound += 1
    a_, b_, al = pt["a"], pt["b"], pt["alpha"]
    labels = list(a_) + list(-b_)
    sep = min(abs(labels[i] - labels[j]) for i in range(8) for j in range(i))
    min_sep = min(min_sep, sep)
    Fv = L.F_values_np(a_, b_, al)
    scale = np.max(np.abs(Fv))
    relsep = sep / max(abs(x) for x in labels)
    if relsep > 1e-4:
        n_distinct += 1
        worstH = max(worstH, abs(Fv[0].real) / scale)
    else:
        n_coll += 1
        coll_H = max(coll_H, abs(Fv[0].real) / scale)
        coll_d35 = min(coll_d35, abs(Fv[3]-Fv[5])/scale, abs(Fv[3]-Fv[9])/scale)
    g = pt["gaps"]
    delta = g[8] - g[2] - g[6] - g[16]
    worst_delta = max(worst_delta, delta / np.max(np.abs(g)))
    # nonnegative combination of the 9 minor generators?
    Aeq = PN @ FV[minor_gens].T
    res = linprog(np.zeros(len(minor_gens)), A_eq=Aeq, b_eq=g, bounds=[(0, None)]*len(minor_gens), method="highs")
    if res.status != 0:
        lp_fail += 1
    cvec = np.linalg.lstsq(PN, g, rcond=None)[0]
    val, _ = L.min_phi_over_SO4(cvec, restarts=8, seed=seed)
    worst_phi = min(worst_phi, val / np.max(np.abs(g)))
print(f"NUMERIC stratum samples found={nfound}/12 (seeds 0..11, tol 1e-10); min label separation={min_sep:.3e}")
print(f"NUMERIC eight-distinct samples (relative label separation > 1e-4): {n_distinct}; samples with a collision: {n_coll}")
print(f"NUMERIC eight-distinct: max |Re(alpha F_0)|/max|F| = {worstH:.3e}   (document: H = 0)")
print(f"NUMERIC collision samples: max |H|/max|F| = {coll_H:.3e}, min min(|F3-F5|,|F3-F9|)/max|F| = {coll_d35:.3e} (H != 0 needs F3=F5 or F3=F9)")
print(f"NUMERIC max delta/max g = {worst_delta:.3e}   (document: delta <= 0 on eight-distinct points)")
print(f"NUMERIC minor-cone LP infeasible in {lp_fail} samples   (document: 0)")
print(f"NUMERIC min over samples of min_SO4 Phi / max g = {worst_phi:.3e} (8 BFGS restarts each)")

# Own sampler on the height-zero locus: alpha = i (so H = Re(alpha F_0) = 0 is imposed),
# equations Re(alpha F_sigma) = 0 on T, mean off-T gap = 1; positivity weight w = 0 (gaps
# may be negative) or w = 1 (gaps >= 0 imposed). Checks C(Re x) = 0 (the necessary identity)
# and records label separation, sign of the off-T gaps and delta.
from scipy.optimize import least_squares
others = [s_ for s_ in range(24) if s_ not in T]
def unpack(x):
    return x[0:4] + 1j*x[4:8], x[8:12] + 1j*x[12:16]
def resid(x, w):
    a_, b_ = unpack(x)
    g = (1j*L.F_values_np(a_, b_)).real
    return np.array([g[s_] for s_ in sorted(T)] + [w*min(0.0, g[s_]) for s_ in others] + [np.mean(g[others]) - 1.0])
rng = np.random.default_rng(7)
for w in (0.0, 1.0):
    n = 0; n8 = 0; worstC = 0.0; dpos_neg = 0; dpos = 0; n8_supported = 0
    for trial in range(40):
        sol = least_squares(lambda x: resid(x, w), rng.normal(size=16), max_nfev=2000)
        if np.max(np.abs(resid(sol.x, w))) > 1e-9:
            continue
        a_, b_ = unpack(sol.x)
        lab = list(a_) + list(-b_)
        rel = min(abs(lab[i_] - lab[j_]) for i_ in range(8) for j_ in range(i_)) / max(abs(z_) for z_ in lab)
        n += 1
        x = 1j * L.matching_vector_np(a_, b_)
        y = x.real
        worstC = max(worstC, abs(float(L.cubic_C(list(y)))) / np.linalg.norm(y)**3)
        g = (1j*L.F_values_np(a_, b_)).real
        dl = g[8] - g[2] - g[6] - g[16]
        if rel > 1e-3:
            n8 += 1
            if min(g[others]) >= -1e-9:
                n8_supported += 1
            if dl > 1e-9:
                dpos += 1
                if min(g[others]) < -1e-9:
                    dpos_neg += 1
    print(f"NUMERIC height-zero sampler w={w}: {n} solutions (40 starts, seed 7), {n8} with relative label separation > 1e-3, "
          f"{n8_supported} of these with all gaps >= 0; max |C(Re x)|/|Re x|^3 = {worstC:.2e}; "
          f"eight-distinct with delta > 0: {dpos}, of which with a negative gap: {dpos_neg}")
print(f"SUMMARY exact failures: {fails}")
sys.exit(1 if fails else 0)
