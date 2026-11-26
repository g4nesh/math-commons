"""Auditor replay for final_referee_inventory82_native_uniform_20261023 (paper case T_83).

Own code only; nothing from the supplement is imported or executed, and no supplement data is
read. Run from any directory:  python -I -B <this file>   (optionally --no-numeric).
Prints one line per check; exits 1 if any EXACT check fails.
"""
import itertools, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + "EXACT " + name + ((" | " + str(info)) if info else ""))
    if not ok: FAIL.append(name)

# ---------- paper proof order (Appendix G, Table 6) ----------
ST6 = {1:[1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
       2:[10,26,46,90],3:[9,11,12,47],4:[8,13,17,49],5:[14,16,24,33,38,70],6:[52,57,58,65,68,71,76,78],
       7:[45,51,66,69,75,81],8:[48,50,53,59,62,67,73,79,86],9:[30,41,54,55,60,61,63,74,77,80,83,84],
       10:[56,64,82,85,87],11:[32],12:[89]}
STAGE6 = {n: s for s, l in ST6.items() for n in l}
ST7 = {1:[1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24],2:[3,4,23,26],4:[6,18],5:[20,27],7:[19,25]}
STAGE7 = {n: s for s, l in ST7.items() for n in l}

def covering(Tj, before_stage):
    """Earlier-stage paper statements (6-sets T_n or 7-sets T^(7)_n) some g in G maps INTO Tj."""
    out = []
    for name, S, st in ([("T%d" % n, S, STAGE6[n]) for n, S in L.CONTACT_SETS.items()] +
                        [("T7_%d" % n, S, STAGE7[n]) for n, S in L.AUX7.items()]):
        if st >= before_stage: continue
        for g, tm in zip(L.G_ELEMENTS, L._G_MAPS):
            if frozenset(tm[j] for j in S) <= Tj:
                out.append((name, st, g)); break
    return out

T = frozenset({3, 5, 6, 9, 14, 18})
MYSTAGE = 9
# ---------- A. contact set is the paper's T_83 ----------
check("contact set equals paper T_83 (identity element)", T == L.CONTACT_SETS[83])
# ---------- source rays (2) and minor profiles ----------
gens = L.generators(T)
check("eq.(2): generators of T are f_1,f_2,f_5,f_9,f_28,f_29,f_39,f_59", gens == [1, 2, 5, 9, 28, 29, 39, 59], gens)
check("rank of contact rows 6, rank of generators 8",
      L.rank([L.P[j] for j in sorted(T)]) == 6 and L.rank([L.F_VECS[j] for j in gens]) == 8)
check("only f_28 is non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [28])
claimed_minor = {1: "n12,34", 2: "n13,34", 5: "n12,24", 9: "n14,23", 29: "n13,12", 39: "n14,13", 59: "n14,14"}
check("seven minor events as listed (0-based pairs (01,23),(02,23),(01,13),(03,12),(02,01),(03,02),(03,03))",
      all(L.f_name(j) == nm for j, nm in claimed_minor.items()))
# (3): theta_28 <= g_sigma whenever f_28(sigma) >= 1 (all coefficients >= 0, unique rep since rank 8)
check("eq.(3): f_28 equals 1 at sigma = 12, 13, 2", all(L.F_VALUES[28][s] == 1 for s in (12, 13, 2)))
extra = [j for j in range(24) if j not in T and L.F_VALUES[28][j] == 0]
check("extra ties not forcing theta_28 = 0 are 0,1,7,8,19,20,22", extra == [0, 1, 7, 8, 19, 20, 22], extra)
for j in extra:
    cv = covering(T | {j}, MYSTAGE)
    check("boundary +%d covered by an earlier-stage statement" % j, len(cv) > 0,
          [(c[0], "stage %d" % c[1]) for c in cv][:4])
H7 = frozenset({0, 1, 6, 9, 16, 17, 22})
check("eq.(5): H = paper Table 5 set T^(7)_21 (stage 1)", H7 == L.AUX7[21])
g = L.find_g(H7, T | {20})
check("+20: some g maps H onto T u {20}", g is not None, g)

# ---------- exact algebra in the normalization A=(inf,0,1,z), C=(x,y,v,t) ----------
PERMS = list(itertools.permutations(range(4)))
x, y, z, v, t, r, s, u, h, al, k, dl, c2 = sp.symbols('x y z v t r s u h alpha k delta c2')
A = [None, 0, 1, z]; Cc = [x, y, v, t]
def Fp(i):
    p = PERMS[i]; e = sp.Integer(1)
    for a in range(1, 4): e *= (A[a] - Cc[p[a]])   # [inf, C] = 1
    return sp.expand(e)
Fs = [Fp(i) for i in range(24)]
E = sp.expand(Fs[5] - Fs[3])
check("eq.(7): E = F5 - F3 = (v-t)(z-y)", sp.expand(E - (v - t) * (z - y)) == 0)
check("F9 - F3 = v(1-t)(x-y)", sp.expand(Fs[9] - Fs[3] - v * (1 - t) * (x - y)) == 0)
H = {i: (Fs[i] - Fs[3]) / E for i in range(24)}
par = {9: r, 6: s, 14: u, 18: h}
claimed8 = {9: v*(1-t)*(x-y) - r*(v-t)*(z-y),
            6: -x*(1-v)*(z-t) + v*(1-t)*(z-y) - s*(v-t)*(z-y),
            14: -y*(1-x)*(z-t) + v*(1-t)*(z-y) - u*(v-t)*(z-y),
            18: -x*(1-y)*(z-v) + v*(1-t)*(z-y) - h*(v-t)*(z-y)}
rows = []
for i in (9, 6, 14, 18):
    Ci = sp.expand(Fs[i] - Fs[3] - par[i] * E)
    check("eq.(8): contact row C%d as displayed" % i, sp.expand(Ci - claimed8[i]) == 0)
    d = sp.Poly(Ci, v, t).as_dict()
    check("C%d linear in (vt, v, t, 1)" % i, set(d) <= {(1, 1), (1, 0), (0, 1), (0, 0)})
    rows.append([d.get((1, 1), 0), d.get((1, 0), 0), d.get((0, 1), 0), d.get((0, 0), 0)])
M = sp.Matrix(rows); B3 = M[:3, :3]; cst = -M[:3, 3]
Delta = sp.expand(B3.det())
N = []
for j in range(3):
    Bj = B3.copy(); Bj[:, j] = cst; N.append(sp.expand(Bj.det()))
check("eq.(9) Cramer: B3 (N0,N1,N2) = Delta * rhs", sp.expand(B3 * sp.Matrix(N) - Delta * cst) == sp.zeros(3, 1))
Draw = sp.cancel(Delta.subs({r: H[9], s: H[6], u: H[14]}))
check("eq.(10): Delta at raw ratios = -xyz(z-1)(x-1)(x-y)/(v-t)",
      sp.cancel(Draw + x*y*z*(z-1)*(x-1)*(x-y)/(v-t)) == 0)
Q = r*(z-y)**2 + (u*x*(y-1) - h*y*(x-1))*(z-y) - x*y*(x-1)*(y-1)
detM = sp.expand(M.det())
check("eq.(11): det M = -xz(z-1)(x-y) Q", sp.expand(detM + x*z*(z-1)*(x-y)*Q) == 0)
check("eq.(12): Q has 15 terms", len(sp.Poly(Q, x, y, z, r, u, h).terms()) == 15)
alpha = y*(x-1)/(z-y); beta = x*(y-1)/(z-y); kk = (z-y)/(x*y)
check("eq.(13)-(14): Q = (z-y)^2 (r + u beta - h alpha - alpha beta)",
      sp.cancel(Q/(z-y)**2 - (r + u*beta - h*alpha - alpha*beta)) == 0)
check("eq.(15): inverse x=1/(1-alpha k), y=1/(1-beta k), z=y+kxy",
      sp.cancel(1/(1-alpha*kk) - x) == 0 and sp.cancel(1/(1-beta*kk) - y) == 0 and sp.cancel(y + kk*x*y - z) == 0)
check("eq.(16a): H20 = H14 - alpha (raw products)", sp.cancel(H[20] - (H[14] - alpha)) == 0)
check("eq.(16b): H12 = H18 + beta (raw products)", sp.cancel(H[12] - (H[18] + beta)) == 0)
check("eq.(16c): H9 - H18 H14 = (alpha - H14)(beta + H18) (raw products)",
      sp.cancel(H[9] - H[18]*H[14] - (alpha - H[14])*(beta + H[18])) == 0)
S = sp.expand(N[0]*Delta - N[1]*N[2])
fl = sp.factor_list(S)
big = [f for f, m in fl[1] if len(sp.Poly(f, x, y, z, r, s, u, h).terms()) > 10]
small = sp.Mul(*[f**m for f, m in fl[1] if f not in big]) * fl[0]
check("eq.(18): N0 Delta - N1 N2 = z(z-1)y(y-z)(x-1)(x-y) * S64 with S64 of 64 terms",
      len(big) == 1 and len(sp.Poly(big[0], x, y, z, r, s, u, h).terms()) == 64 and
      sp.cancel(small / (z*(z-1)*y*(y-z)*(x-1)*(x-y))) in (1, -1))
S64 = big[0]
R = r - h*u
betaA = (h*al - r)/(u - al)
X = 1/(1-al*k); Y = 1/(1-betaA*k); Z = Y + k*X*Y
sub = {x: X, y: Y, z: Z}
num = sp.expand(sp.numer(sp.cancel(sp.together(S64.subs(sub)))))
fl2 = sp.factor_list(num)
bigC = [f for f, m in fl2[1] if len(sp.Poly(f, al, k, r, s, u, h).terms()) > 20]
rest = sp.Mul(*[f**m for f, m in fl2[1] if f not in bigC])
check("eq.(19): numerator = (const) alpha k^2 * C, C of 179 terms, bidegree (5,2)",
      len(bigC) == 1 and sp.cancel(rest/(al*k**2)) in (1, -1) and
      len(sp.Poly(bigC[0], al, k, r, s, u, h).terms()) == 179 and
      sp.Poly(bigC[0], al, k).degree(al) == 5 and sp.Poly(bigC[0], al, k).degree(k) == 2)
Cak = bigC[0]
if sp.expand(Cak.subs(al, u) + R**3*k*(u*k-1)) != 0: Cak = -Cak
check("eq.(20): C(u,k) = -R^3 k(uk-1) (after sign normalization)", sp.expand(Cak.subs(al, u) + R**3*k*(u*k-1)) == 0)

def Hcr(i):
    d = sp.Poly(sp.expand(Fs[i] - Fs[3]), v, t).as_dict()
    nm = d.get((1, 1), 0)*N[0] + d.get((1, 0), 0)*N[1] + d.get((0, 1), 0)*N[2] + d.get((0, 0), 0)*Delta
    return nm/((N[1] - N[2])*(z - y))
def on_curve(expr):
    n = sp.expand(sp.numer(sp.cancel(sp.together(expr.subs(sub)))))
    if n == 0: return "identically"
    return "mod C" if sp.prem(sp.Poly(n, k), sp.Poly(Cak, k)).is_zero else False
H12f = R/(al - u)
for i, val in [(9, r), (6, s), (14, u), (18, h), (5, sp.Integer(1))]:
    res = on_curve(Hcr(i) - val)
    check("Cramer coordinate function H%d equals its fixed value on C=0" % i, bool(res), res)
res = on_curve(Hcr(12) - H12f); check("eq.(17): H12 = R/(alpha-u) on C=0", bool(res), res)
res = on_curve(Hcr(8) - (u + al*k*(s - H12f))); check("eq.(21): H8 on C=0", bool(res), res)
H13f = k*((s - r - betaA)*(1 - al*k) - 1)/((1 - al*k)*(1 + (1 - al)*k))
res = on_curve(Hcr(13) - H13f); check("eq.(22): H13 on C=0", bool(res), res)
H2f = H12f*(1 - al*k) + k*(r + s*betaA*(1 - al*k))/(1 - betaA*k)
res = on_curve(Hcr(2) - H2f); check("eq.(23): H2 on C=0", bool(res), res)
lead = sp.Poly(Cak, k).LC()
check("fibre alpha=u: leading k^2 coefficient = -u R^3", sp.expand(lead.subs(al, u) + u*R**3) == 0)
Ck = sp.diff(Cak, k); Ca = sp.diff(Cak, al)
check("eq.(24): C_k(u,0) = R^3, C_k(u,1/u) = -R^3",
      sp.expand(Ck.subs({al: u, k: 0}) - R**3) == 0 and sp.cancel(Ck.subs({al: u, k: 1/u}) + R**3) == 0)
k1a = sp.cancel(-Ca.subs({al: u, k: 0})/Ck.subs({al: u, k: 0}))
k1b = sp.cancel(-Ca.subs({al: u, k: 1/u})/Ck.subs({al: u, k: 1/u}))
check("eq.(25): dk/dalpha = 1/R at (u,0) and -(R-u)/(u^2 R) at (u,1/u)",
      sp.cancel(k1a - 1/R) == 0 and sp.cancel(k1b + (R - u)/(u**2*R)) == 0, (k1a, k1b))
lim1 = sp.limit(sp.cancel(H13f.subs({al: u + dl, k: dl/R + c2*dl**2}, simultaneous=True)), dl, 0)
check("eq.(26): H13 -> -1 along the k=0 branch (any second-order term)", sp.simplify(lim1 + 1) == 0, lim1)
lim2 = sp.limit(sp.cancel(H2f.subs({al: u + dl, k: 1/u - (R - u)/(u**2*R)*dl + c2*dl**2}, simultaneous=True)), dl, 0)
check("eq.(27): H2 -> -1 along the k=1/u branch (any second-order term)", sp.simplify(lim2 + 1) == 0, lim2)
check("u=0: section C(alpha=u=0, k) = r^3 k, a nonzero polynomial when R = r != 0 (Rouche step input)",
      sp.expand(Cak.subs(u, 0).subs(al, 0) - r**3*k) == 0)

# ---------- D. numeric falsification ----------
if "--no-numeric" not in sys.argv:
    import numpy as np
    Pn = L.P_numpy(); found = 0; worst = np.inf; t0 = time.time()
    for sd in range(500, 700):
        if found >= 6 or time.time() - t0 > 150: break
        pt = L.random_point_on_contact_stratum(T, seed=sd, tries=8, margin=1e-3)
        if pt is None: continue
        found += 1
        m = L.matching_vector_np(pt['a'], pt['b'], pt['alpha'])
        c = m.real.copy(); c[13] -= (Pn @ m)[pt['pi_star']].real
        vals = Pn @ c; val, _ = L.min_phi_over_SO4(c, restarts=10, seed=sd)
        worst = min(worst, val/np.max(np.abs(vals)))
    print("INFO NUMERIC T_83 stratum: %d points (seeds 500.., margin 1e-3), 10 BFGS restarts each; "
          "worst min Phi / max gap = %.3e" % (found, worst))
print("SUMMARY: %d EXACT failures" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
