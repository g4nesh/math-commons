"""Second-pass replay for supplement document final_referee_inventory82_native_uniform_20261023
(primary proof of paper case T_83 = {3,5,6,9,14,18}).

Own code only; uses the shared audit library ../lib/moc4lib.py. Nothing from the supplement is
executed, imported or read. Every formula is rebuilt from the definitions
F_sigma = prod_i [A_i, C_sigma(i)], A = (inf, 0, 1, z), C = (x, y, v, t)  (doc eq. (6)).
Prints one line per check; exits nonzero if an EXACT check fails. NUMERIC lines are evidence only.
"""
import itertools, os, sys, time
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
    print(f"{'PASS' if ok else 'FAIL'} {name} {extra}".rstrip(), flush=True)

t0 = time.time()
T = frozenset({3, 5, 6, 9, 14, 18})
# ---------------------------------------------------------------- A. statement / combinatorics
report("A1 doc contact set equals paper T_83 literally", L.CONTACT_SETS[83] == T)
gens = L.generators(T)
report("A2 generators(T) = (2) {1,2,5,9,28,29,39,59}", gens == [1, 2, 5, 9, 28, 29, 39, 59],
       str([L.f_name(j) for j in gens]))
report("A3 only f28 = h(3421,-) is non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [28])
minor01 = []
for j in gens:
    if L.F_IS_MINOR[j]:
        _, a, b = L.F_TABLE[j]
        minor01.append("".join(str(int(c) - 1) for c in a) + "," + "".join(str(int(c) - 1) for c in b))
report("A4 seven minor generators are the 0-based (I,J) listed in the doc, in order",
       minor01 == ["01,23", "02,23", "01,13", "03,12", "02,01", "03,02", "03,03"], str(minor01))
report("A5 contact rows of P have rank 6; generators have rank 8",
       L.rank([L.P[s] for s in sorted(T)]) == 6 and L.rank([L.F_VECS[j] for j in gens]) == 8)
f28 = L.F_VALUES[28]
report("A6 (3) f28 = 1 at 12, 13, 2 and every f_j >= 0 everywhere (so theta <= g12, g13, g2)",
       f28[12] == 1 and f28[13] == 1 and f28[2] == 1 and all(min(L.F_VALUES[j]) >= 0 for j in gens))
extra = sorted(L.F_ZERO[28] - T)
report("A7 non-forcing extra ties = zeros of f28 off T = {0,1,7,8,19,20,22}", extra == [0, 1, 7, 8, 19, 20, 22])
# all other extra ties j: f28(j) > 0 so theta = 0 (generators of N_{T+j} all minor)
others = [j for j in range(24) if j not in T and j not in extra]
report("A8 every other extra tie leaves only minor generators",
       all(all(L.F_IS_MINOR[i] for i in L.generators(T | {j})) for j in others))
orbit_of = {L.orbit_rep(S): k for k, S in L.CONTACT_SETS.items()}
STAGE = {}
for st, ks in {1: [1, 2, 3, 4, 5, 6, 7, 15, 18, 19, 20, 21, 22, 23, 25, 27, 28, 29, 31, 34, 35, 36, 37, 39, 40, 42, 43, 44, 72, 88],
               2: [10, 26, 46, 90], 3: [9, 11, 12, 47], 4: [8, 13, 17, 49], 5: [14, 16, 24, 33, 38, 70],
               6: [52, 57, 58, 65, 68, 71, 76, 78], 7: [45, 51, 66, 69, 75, 81], 8: [48, 50, 53, 59, 62, 67, 73, 79, 86],
               9: [30, 41, 54, 55, 60, 61, 63, 74, 77, 80, 83, 84], 10: [56, 64, 82, 85, 87], 11: [32], 12: [89]}.items():
    for kk in ks:
        STAGE[kk] = st
for j in [0, 1, 7, 8, 19, 22]:
    U = T | {j}
    cov = []
    for S in itertools.combinations(sorted(U), 6):
        S = frozenset(S)
        if S == T:
            continue
        kk = orbit_of.get(L.orbit_rep(S))
        if kk is not None and STAGE[kk] < STAGE[83]:
            g = L.find_g(L.CONTACT_SETS[kk], S)
            cov.append((kk, STAGE[kk], sorted(S), g is not None))
    best = min(cov, key=lambda c: c[1]) if cov else None
    report(f"A9 tie +{j}: T+{j} contains an image of an earlier-stage paper case",
           best is not None and best[3], f"e.g. T_{best[0]} (stage {best[1]}) -> {best[2]}; all: {[(c[0], c[1]) for c in cov]}")
H = frozenset({0, 1, 6, 9, 16, 17, 22})
g = L.find_g(H, T | {20})
report("A10 (5) H = Table 5 aux set T^(7)_21 (stage 1) and some g in G maps H onto T+{20}",
       H == L.AUX7[21] and g is not None, f"g={g}")
print(f"INFO time {time.time()-t0:.1f}s", flush=True)

# ---------------------------------------------------------------- C. algebra
x, y, z, v, t, r, s, u, h, al, k, dl, c2 = sp.symbols('x y z v t r s u h alpha k delta c2')
PERMS = L.PERMS
Apts = [None, 0, 1, z]
Cpts = [x, y, v, t]
def Fdoc(i):
    p = PERMS[i]
    out = sp.Integer(1)
    for a in range(1, 4):           # [inf, C] = 1
        out *= (Apts[a] - Cpts[p[a]])
    return sp.expand(out)
Fs = [Fdoc(i) for i in range(24)]
E = sp.expand(Fs[5] - Fs[3])
report("C1 (7) E = F5-F3 = (v-t)(z-y)", sp.expand(E - (v - t) * (z - y)) == 0)
report("C2 (7) F9-F3 = v(1-t)(x-y)", sp.expand(Fs[9] - Fs[3] - v * (1 - t) * (x - y)) == 0)
Hraw = [sp.cancel((Fs[i] - Fs[3]) / E) for i in range(24)]
Cdoc = {9: v*(1-t)*(x-y) - r*(v-t)*(z-y),
        6: -x*(1-v)*(z-t) + v*(1-t)*(z-y) - s*(v-t)*(z-y),
        14: -y*(1-x)*(z-t) + v*(1-t)*(z-y) - u*(v-t)*(z-y),
        18: -x*(1-y)*(z-v) + v*(1-t)*(z-y) - h*(v-t)*(z-y)}
par = {9: r, 6: s, 14: u, 18: h}
report("C3 (8) the four contact equations are (F_i - F_3) - H_i E for i = 9,6,14,18",
       all(sp.expand(Fs[i] - Fs[3] - par[i] * E - Cdoc[i]) == 0 for i in par))
rows = []
lin_ok = True
for i in (9, 6, 14, 18):
    P = sp.Poly(sp.expand(Cdoc[i]), v, t)
    lin_ok &= all(m[0] <= 1 and m[1] <= 1 for m in P.monoms())
    rows.append([P.coeff_monomial(v * t), P.coeff_monomial(v), P.coeff_monomial(t), P.coeff_monomial(1)])
report("C4 each contact equation is linear in (vt, v, t, 1)", lin_ok)
M = sp.Matrix(rows)
D3 = M[:3, :3]
Delta = sp.expand(D3.det())
N = []
for j in range(3):
    B = D3.copy()
    B[:, j] = -M[:3, 3]
    N.append(sp.expand(B.det()))
report("C5 (9) Cramer identities D3 (N0,N1,N2)^t = -Delta * const",
       all(sp.expand(sum(D3[i, j] * N[j] for j in range(3)) + Delta * M[i, 3]) == 0 for i in range(3)))
Dsub = sp.cancel(Delta.subs({r: Hraw[9], s: Hraw[6], u: Hraw[14], h: Hraw[18]}))
report("C6 (10) Delta at the raw ratios = -xyz(z-1)(x-1)(x-y)/(v-t)",
       sp.cancel(Dsub + x*y*z*(z-1)*(x-1)*(x-y)/(v-t)) == 0)
Q = r*(z-y)**2 + (u*x*(y-1) - h*y*(x-1))*(z-y) - x*y*(x-1)*(y-1)
report("C7 (11),(12) det M = -xz(z-1)(x-y)Q, Q has 15 terms",
       sp.expand(M.det() + x*z*(z-1)*(x-y)*Q) == 0 and len(sp.Poly(sp.expand(Q), x, y, z, r, u, h).terms()) == 15)
alpha_e = y*(x-1)/(z-y); beta_e = x*(y-1)/(z-y); k_e = (z-y)/(x*y)
report("C8 (13),(14) Q = (z-y)^2 (r + u beta - h alpha - alpha beta)",
       sp.cancel(Q - (z-y)**2*(r + u*beta_e - h*alpha_e - alpha_e*beta_e)) == 0)
report("C9 (15) inverse positions", sp.cancel(1/(1-alpha_e*k_e) - x) == 0 and sp.cancel(1/(1-beta_e*k_e) - y) == 0
       and sp.cancel(y + k_e*x*y - z) == 0)
report("C10 (16) raw identities H20 = H14 - alpha, H12 = H18 + beta, H9 - H18 H14 = (alpha - H14)(beta + H18)",
       sp.cancel(Hraw[20] - (Hraw[14] - alpha_e)) == 0 and sp.cancel(Hraw[12] - (Hraw[18] + beta_e)) == 0
       and sp.cancel(Hraw[9] - Hraw[18]*Hraw[14] - (alpha_e - Hraw[14])*(beta_e + Hraw[18])) == 0)
SG = sp.expand(N[0]*Delta - N[1]*N[2])
pref = sp.expand(z*(z-1)*y*(y-z)*(x-1)*(x-y))
S64, rem = sp.div(SG, pref, x, y, z, r, s, u, h)
nS = len(sp.Poly(S64, x, y, z, r, s, u, h).terms())
report("C11 (18) N0 Delta - N1 N2 = z(z-1)y(y-z)(x-1)(x-y) S64, S64 has 64 terms", rem == 0 and nS == 64,
       f"S64 variables: {sorted(map(str, S64.free_symbols))} (doc writes S64(x,y,z))")
R = r - h*u
bsub = (h*al - r)/(u - al)
X = 1/(1 - al*k); Y = 1/(1 - bsub*k); Z = Y + k*X*Y
num, den = sp.fraction(sp.cancel(sp.together(S64.subs({x: X, y: Y, z: Z}, simultaneous=True))))
fl = sp.factor_list(sp.expand(num))
bigs = [f for f, e in fl[1] if sp.Poly(f, al, k).total_degree() > 2]
small = sorted(str(f**e) for f, e in fl[1] if sp.Poly(f, al, k).total_degree() <= 2)
Cn = bigs[0] if len(bigs) == 1 else None
ratio = sp.cancel(sp.expand(Cn.subs(al, u)) / (-R**3*k*(u*k-1))) if Cn is not None else None
if ratio is not None and ratio.is_number:
    Cn = sp.expand(Cn / ratio)
PC = sp.Poly(Cn, al, k)
report("C12 (19) numerator = (alpha k^2) x one irreducible factor C with 179 terms, bidegree (5,2)",
       Cn is not None and small == ['alpha', 'k**2'] and len(sp.Poly(Cn, al, k, r, s, u, h).terms()) == 179
       and PC.degree(al) == 5 and PC.degree(k) == 2, f"cleared denominator {sp.factor(den)}; small factors {small}")
report("C13 (20) C(u,k) = -R^3 k(uk-1) exactly", sp.expand(Cn.subs(al, u) + R**3*k*(u*k-1)) == 0)
report("C14 k^2-coefficient of C(u,k) is -u R^3 (no point at k = inf over alpha = u when uR != 0)",
       sp.expand(sp.Poly(Cn.subs(al, u), k).coeff_monomial(k**2) + u*R**3) == 0)
Ck = sp.diff(Cn, k); Ca = sp.diff(Cn, al)
report("C15 (24) C_k(u,0) = R^3, C_k(u,1/u) = -R^3",
       sp.expand(Ck.subs({al: u, k: 0}) - R**3) == 0 and sp.cancel(Ck.subs({al: u, k: 1/u}) + R**3) == 0)
d1 = sp.cancel(-Ca.subs({al: u, k: 0}) / Ck.subs({al: u, k: 0}))
d2 = sp.cancel(-Ca.subs({al: u, k: 1/u}) / Ck.subs({al: u, k: 1/u}))
report("C16 (25) dk/dalpha = 1/R at (u,0) and -(R-u)/(u^2 R) at (u,1/u)",
       sp.cancel(d1 - 1/R) == 0 and sp.cancel(d2 + (R-u)/(u**2*R)) == 0)
# coordinate functions by the linearised Cramer route (valid on the Segre locus N0 Delta = N1 N2)
xyz = {x: X, y: Y, z: Z}
Ns = [sp.together(n.subs(xyz, simultaneous=True)) for n in N]
Ds = sp.together(Delta.subs(xyz, simultaneous=True))
def Hlin(i):
    P = sp.Poly(sp.expand(Fs[i] - Fs[3]), v, t)
    cs = [sp.together(P.coeff_monomial(m).subs(xyz, simultaneous=True)) for m in (v*t, v, t, 1)]
    return (cs[0]*Ns[0] + cs[1]*Ns[1] + cs[2]*Ns[2] + cs[3]*Ds) / ((Ns[1] - Ns[2])*(Z - Y))
def zero_on_curve(expr):
    nm, _ = sp.fraction(sp.cancel(sp.together(expr)))
    nm = sp.expand(nm)
    if nm == 0:
        return True, "identically"
    rm = sp.prem(sp.Poly(nm, k), sp.Poly(Cn, k))
    return sp.expand(rm.as_expr()) == 0, "modulo C"
H12c = R/(al - u)
H13c = k*((s - r - bsub)*(1 - al*k) - 1)/((1 - al*k)*(1 + (1 - al)*k))
H2c = H12c*(1 - al*k) + k*(r + s*bsub*(1 - al*k))/(1 - bsub*k)
for name, i, rhs in [("H5 = 1", 5, sp.Integer(1)), ("H9 = r", 9, r), ("H6 = s", 6, s), ("H14 = u", 14, u),
                     ("H18 = h", 18, h), ("(17) H12 = R/(alpha-u)", 12, H12c), ("(16) H20 = u - alpha", 20, u - al),
                     ("(21) H8", 8, u + al*k*(s - H12c)), ("(22) H13", 13, H13c), ("(23) H2", 2, H2c)]:
    ok, how = zero_on_curve(Hlin(i) - rhs)
    report(f"C17 coordinate function {name} on the curve", ok, how)
# (26), (27): limits along the two branches over alpha = u, k = k0 + k1 delta + c2 delta^2 (c2 free)
b_d = R/dl - h
lim13 = sp.cancel(H13c.subs({al: u + dl, k: dl/R + c2*dl**2}).subs(bsub, b_d))
v13 = sp.cancel(sp.cancel(lim13).subs(dl, 0))
report("C18 (25),(26) beta = R/delta - h exactly and H13 -> -1 at (u,0), independent of the delta^2 term",
       sp.cancel(bsub.subs(al, u + dl) - b_d) == 0 and v13 == -1, f"limit {v13}")
k2 = 1/u - (R - u)/(u**2*R)*dl + c2*dl**2
lim2 = sp.cancel(H2c.subs({al: u + dl, k: k2}))
v2 = sp.cancel(lim2.subs(dl, 0))
report("C19 (27) H2 -> -1 at (u,1/u), independent of the delta^2 term", v2 == -1, f"limit {v2}")
print(f"INFO time {time.time()-t0:.1f}s", flush=True)

# ---------------------------------------------------------------- D. numeric falsification
import numpy as np
Cfun = sp.lambdify((al, k, r, s, u, h), Cn, "numpy")
def mobius_coords(a, b):
    a1, a2, a3, a4 = a
    m = lambda w: (w - a2)*(a3 - a1)/((w - a1)*(a3 - a2))
    return m(a4), [m(-bj) for bj in b]
worst_phi = np.inf; worst_id = 0.0; worst_rep = 0.0; npts = 0
for seed in range(8):
    pt = L.random_point_on_contact_stratum(T, seed=seed, tries=40, margin=1e-2, pi_star=3)
    if pt is None:
        print(f"NUMERIC seed {seed}: no stratum point found")
        continue
    a, b, alph = pt["a"], pt["b"], pt["alpha"]
    if sorted(pt["contact"]) != sorted(T):
        print(f"NUMERIC seed {seed}: contact {pt['contact']} != T, skipped")
        continue
    npts += 1
    Fv = L.F_values_np(a, b, 1.0)
    Hn = (Fv - Fv[3])/(Fv[5] - Fv[3])
    zz, (xx, yy, vv, tt) = mobius_coords(a, b)
    rr, ss, uu, hh = Hn[9], Hn[6], Hn[14], Hn[18]
    A_ = yy*(xx-1)/(zz-yy); B_ = xx*(yy-1)/(zz-yy); K_ = (zz-yy)/(xx*yy)
    Rn = rr - hh*uu
    errs = [abs(Hn[i].imag) for i in (9, 6, 14, 18)]
    errs += [abs(Hn[20] - (uu - A_)), abs(Hn[12] - (hh + B_)), abs(Hn[12] - Rn/(A_ - uu)),
             abs(Hn[8] - (uu + A_*K_*(ss - Hn[12]))),
             abs(Hn[13] - K_*((ss-rr-B_)*(1-A_*K_)-1)/((1-A_*K_)*(1+(1-A_)*K_))),
             abs(Hn[2] - (Hn[12]*(1-A_*K_) + K_*(rr + ss*B_*(1-A_*K_))/(1-B_*K_))),
             abs(Cfun(A_, K_, rr.real, ss.real, uu.real, hh.real))/max(1.0, abs(Cfun(A_ + 0.1, K_, rr.real, ss.real, uu.real, hh.real)))]
    worst_id = max(worst_id, max(errs))
    mvec = L.matching_vector_np(a, b, alph)
    Pm = L.P_numpy() @ mvec
    c = mvec.real.copy(); c[13] -= Pm[3].real
    gaps = L.P_numpy() @ c
    # unique generator coefficients (rank 8); check (3)
    G = np.array([[float(q) for q in L.F_VALUES[j]] for j in gens]).T
    lam, *_ = np.linalg.lstsq(G, gaps, rcond=None)
    th = lam[gens.index(28)]
    worst_rep = max(worst_rep, np.max(np.abs(G @ lam - gaps)), -min(lam.min(), 0),
                    max(0.0, th - min(gaps[12], gaps[13], gaps[2])))
    val, O = L.min_phi_over_SO4(c, restarts=20, seed=seed)
    scale = np.max(np.abs(gaps))
    worst_phi = min(worst_phi, val/scale)
    print(f"NUMERIC seed {seed}: min Phi/max gap = {val/scale:.3e}, theta = {th:.3e}, identity residual {max(errs):.1e}")
print(f"NUMERIC summary: {npts} points with contact set exactly T (margin 1e-2), 20 BFGS restarts each; "
      f"worst normalized min Phi = {worst_phi:.3e}; worst identity residual (14),(16),(17),(21)-(23),(19) = {worst_id:.1e}; "
      f"worst generator-representation / (3) defect = {worst_rep:.1e}")
print(f"INFO total time {time.time()-t0:.1f}s; exact failures: {fails}")
sys.exit(1 if fails else 0)
