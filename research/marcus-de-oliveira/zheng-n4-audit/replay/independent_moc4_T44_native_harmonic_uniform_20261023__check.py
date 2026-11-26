"""Own-code replay for independent_moc4_T44_native_harmonic_uniform_20261023 (paper case T_45). Audit batch C18.
Uses only the shared audit library (../lib/moc4lib.py); executes nothing from the supplement."""
# ---- common helpers (own code; duplicated in each C18 check script) ----
import sys, os, itertools, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import flint
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)

def info(name, s):
    print("INFO " + name + " :: " + str(s))

def vec_comb(coeffs):
    """sum_j coeff_j f_j as a V_R vector; coeffs: dict j -> number/sympy"""
    out = [0] * 14
    for j, c in coeffs.items():
        for k in range(14):
            out[k] = out[k] + c * L.F_VECS[j][k]
    return out

def prof(v):
    return L.values(v)

def lin_identity(lhs, rhs):
    """profiles of sum lhs f_j == sum rhs f_j"""
    a = prof(vec_comb({j: 1 for j in lhs}) if isinstance(lhs, list) else vec_comb(lhs))
    b = prof(vec_comb({j: 1 for j in rhs}) if isinstance(rhs, list) else vec_comb(rhs))
    return all(x == y for x, y in zip(a, b))

# Hessian of the cubic: H(m)[a][b] = sum_k HT[a][b][k] m_k
_ms, _Cs = L.cubic_sympy()
HT = [[[int(sp.diff(_Cs, _ms[a], _ms[b], _ms[k])) for k in range(14)] for b in range(14)] for a in range(14)]

def hess_det(m):
    M = [[sum(HT[a][b][k] * m[k] for k in range(14) if HT[a][b][k]) for b in range(14)] for a in range(14)]
    return flint.fmpq_mat(14, 14, [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for row in M for x in row]).det()

def CQ(g):
    """C(g + z kappa) = C0 + z*Q (checked linear in z by caller)"""
    C0 = L.cubic_C(g)
    gk = list(g); gk[13] = gk[13] + 1
    C1 = L.cubic_C(gk)
    return C0, C1 - C0

def radial_G2(g):
    """-Q^2 det H(g - (C/Q) kappa)/4  (= G^2 up to the calibration constant)"""
    C0, Q = CQ(g)
    w = list(g); w[13] = w[13] - Fr(C0) / Fr(Q)
    D = hess_det(w)
    return -flint.fmpq(Fr(Q).numerator, Fr(Q).denominator) ** 2 * D / 4, C0, Q

def to_fmpq(x):
    x = Fr(x); return flint.fmpq(x.numerator, x.denominator)

def interp(ts, ys):
    """Lagrange interpolation over Q -> fmpq_poly"""
    X = flint.fmpq_poly([0, 1])
    P = flint.fmpq_poly([0])
    for i, ti in enumerate(ts):
        term = flint.fmpq_poly([ys[i]])
        for j, tj in enumerate(ts):
            if j != i:
                term = term * (X - to_fmpq(tj)) / (to_fmpq(ti) - to_fmpq(tj))
        P = P + term
    return P
# ---- end common helpers ----

# ================= independent_moc4_T44_native_harmonic_uniform_20261023 (paper T_45) =================
T = frozenset({1, 2, 3, 7, 11, 12})
check("contact set equals paper T_45 literally", T == L.CONTACT_SETS[45])
gens = L.generators(T)
check("generators(T) = {0,1,2,7,8,12,15,18,53}", gens == [0,1,2,7,8,12,15,18,53], gens)
check("only nonminor generator is 18", [j for j in gens if not L.F_IS_MINOR[j]] == [18])
check("identity (1): f18 + f0 = f2 + f1 + f7 + f15", lin_identity([18, 0], [2, 1, 7, 15]))
NAMES = ["q","p","j","c","e","d","b","a"]; IDX = [2, 1, 53, 12, 8, 18, 7, 15]
co = sp.symbols("q p j c e d b a"); qq, pp, jj, cc, ee, dd, bb, aa = co
B8 = [L.F_VECS[j] for j in IDX]
check("the eight generators of (2) are independent", L.rank(B8) == 8)
gv = L.values(vec_comb({IDX[i]: co[i] for i in range(8)}))
# express each coefficient as a single gap where possible
excl = {}
for i, n in enumerate(NAMES):
    excl[n] = [s for s in range(24) if sp.expand(gv[s] - co[i]) == 0]
info("gaps equal to a single coefficient", excl)
check("seven minor coefficients are 'exclusive gaps'", all(excl[n] for n in NAMES if n != "d"))
check("a+d = g4, d = g4 - g18, q+d = g21 + d - ... : g4 = a+d and g18 = a",
      sp.expand(gv[4] - (aa+dd)) == 0 and sp.expand(gv[18] - aa) == 0)
check("q = g21 (so q + d = g21 + d)", sp.expand(gv[21] - qq) == 0)
for n, i in (("q", 0), ("p", 1), ("b", 6), ("a", 7)):
    w = [s for s in range(24) if sp.expand(gv[s] - (co[i] + dd)) == 0]
    check("%s + d is a gap value (needed for d<=0 absorption)" % n, len(w) > 0, w)
# ---- seventh ties ----
f18z = sorted(L.F_ZERO[18] - T)
check("extra zeros of f18 = {0,6,9,18,20,21,23}", f18z == [0,6,9,18,20,21,23], f18z)
STAGE = {}
for st, ks in {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
               2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
               6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81]}.items():
    for k in ks: STAGE[k] = st
REPS = {L.orbit_rep(L.CONTACT_SETS[k]): k for k in L.CONTACT_SETS}
allmin_reps = {}
def classify6(S):
    rep = L.orbit_rep(S)
    if rep in REPS: return ("T_%d" % REPS[rep], STAGE.get(REPS[rep], 99))
    if all(L.F_IS_MINOR[j] for j in L.generators(S)): return ("all-minor", 0)
    return ("??", 99)
okall = True
for s7 in f18z:
    S7 = sorted(T | {s7})
    cls = [classify6(set(S7) - {x}) for x in S7]
    best = min(cls, key=lambda t: t[1])
    good = best[1] <= 6
    okall = okall and good
    print("   seventh tie %2d: 6-subsets -> %s ; earliest = %s" % (s7, cls, best))
check("every seventh tie contains a six-set closed at a stage < 7 (or all-minor)", okall)

# ---- the contact curve chart ----
Z, W, r, s, u, h = sp.symbols("Z W r s u h")
Y = r*Z + 1 - r; beta = 1 - r
A2s = W*Y - beta; Ls = s*A2s - Z*W; Ks = 1 - W*(1 - Z)
A_sp = [(-s*A2s*(Ks + u*Z*W), Ks*(s*A2s - W) - u*s*A2s*W*(1 - Z)), (Z, 1 - Z), (W*Y - beta, 1), (W - 1, 1)]
B_sp = [(s*A2s, s*A2s - W), (0, 1), (1, 1), (1, 0)]
def br(P_, Q_): return P_[0]*Q_[1] + P_[1]*Q_[0]   # (a_i + b_j): B spinors carry b-values, paired via v_{4+j} = (-b_j, 1)
F = [sp.expand(sp.Mul(*[br(A_sp[i], B_sp[sig[i]]) for i in range(4)])) for sig in L.PERMS]
E = sp.factor(F[7] - F[1])
check("E = F7 - F1 = -s A2 W K L", sp.expand(E + s*A2s*W*Ks*Ls) == 0, E)
Hs = [sp.cancel((F[k] - F[1]) / (F[7] - F[1])) for k in range(24)]
claims = {0: r*u*Z, 1: 0, 2: u*(r-1), 3: -u, 4: -u*W*Y, 6: beta + r*(u+1)*Z + r/W, 7: 1, 12: u*(r-1) + 1/s}
for k, v in claims.items():
    check("H_%d formula" % k, sp.simplify(Hs[k] - v) == 0)
check("D := H4 - H18 = W(Y - 1/s) - beta", sp.simplify(Hs[4] - Hs[18] - (W*(Y - 1/s) - beta)) == 0)
check("H21 - H11 + D = (W-1)/(s A2) - 1 - 1/s", sp.simplify(Hs[21] - Hs[11] + (W*(Y - 1/s) - beta) - ((W-1)/(s*A2s) - 1 - 1/s)) == 0)
check("F1/E = u Z W / K", sp.simplify(F[1]/(F[7]-F[1]) - u*Z*W/Ks) == 0)
# denominators of every H_sigma
dens = set()
for k in range(24):
    dn = sp.factor(sp.denom(sp.together(Hs[k])))
    for fac, _ in sp.factor_list(dn)[1]:
        dens.add(sp.factor(fac))
info("irreducible denominator factors of all H_sigma", dens)
allowed = [sp.factor(x) for x in (W, A2s, Ls, s)]
def divides_allowed(fc):
    return any(sp.simplify(sp.cancel(fc / al)).is_number for al in allowed) or fc.free_symbols <= {s}
check("every denominator factor of every H_sigma divides W*A2*L (or is a constant in s)", all(divides_allowed(f_) for f_ in dens))
# the curve
N = sp.cancel(Hs[11] * Ls)
check("H11 * L is a polynomial N in (Z,W)", sp.denom(sp.together(N)).free_symbols <= {s}, sp.denom(sp.together(N)))
N = sp.expand(N)
Cc = sp.Poly(sp.expand(N - h*Ls), W)
check("curve is quadratic in W", Cc.degree() == 2)
Acoef = sp.factor(Cc.coeff_monomial(W**2)); Bcoef = sp.expand(Cc.coeff_monomial(W)); C0coef = sp.expand(Cc.coeff_monomial(1))
check("A(Z) = -{(u+1)Z - 1} Y (sY - 1)", sp.expand(Acoef + ((u+1)*Z - 1)*Y*(s*Y - 1)) == 0, Acoef)
check("C(Z) = curve at W=0 has degree 1 in Z", sp.Poly(C0coef, Z).degree() == 1, sp.factor(C0coef))
curve = sp.expand(N - h*Ls)
info("curve bidegree (Z, W) and B(Z) degree", (sp.Poly(curve, Z).degree(), sp.Poly(curve, W).degree(), sp.Poly(Bcoef, Z).degree()))
check("basepoint (0,1) lies on the curve for all parameters", sp.expand(curve.subs({Z: 0, W: 1})) == 0)
# resultants
RL = sp.factor(sp.resultant(curve, Ls, W)); RA = sp.factor(sp.resultant(curve, A2s, W)); RK = sp.factor(sp.resultant(curve, Ks, W))
info("Res_W(C,L)", RL); info("Res_W(C,A2)", RA); info("Res_W(C,K)", RK)
cl = (r*s*u - s*u - s + 1)*Z*(Z - 1)*(Z*r*s - Z*r - r*s + s)
ca = Z*Y*(Z*r**2*u - Z*r*u - Z*r - r*u - r*h + r + u + h)
check("Res(C,L) = const * (rsu-su-s+1) Z (Z-1)(Zrs-Zr-rs+s)", sp.simplify(RL / cl).is_number or sp.simplify(RL/cl).free_symbols <= {s, r, u} and sp.Poly(sp.cancel(RL/cl), Z).degree() == 0, sp.cancel(RL/cl))
check("Res(C,A2) = const * Z Y {Zr^2u - Zru - Zr - ru - rh + r + u + h}", sp.Poly(sp.cancel(RA/ca), Z).degree() == 0 and sp.denom(sp.cancel(RA/ca)).free_symbols <= {r, s, u}, sp.cancel(RA/ca))
rk = sp.cancel(RK / Z)
check("Res(C,K) = Z * (quadratic in Z) (up to a Z-free factor)", sp.Poly(sp.numer(sp.together(rk)), Z).degree() == 2)
spec = {r: 2, s: 3, u: 4, h: 5}
Q2t = sp.factor(rk.subs(spec)); info("Res(C,K)/Z at (r,s,u,h)=(2,3,4,5)", Q2t)
check("at (2,3,4,5) Q2 proportional to 16Z^2-12Z+4", sp.cancel(Q2t / (16*Z**2 - 12*Z + 4)).is_number)
# W = infinity: roots of A ; Z = infinity: leading equation (5)
Vv = sp.symbols("V")
lead = sp.expand(sp.cancel(curve.subs(W, Vv/Z) * Z**2))
lp = sp.Poly(sp.expand(lead), Z)
top = sp.factor(lp.coeffs()[0]); info("Z=inf: leading coefficient in Z after W=V/Z, degree", (lp.degree(), top))
check("(5): leading equation is -r{(u+1)V+1}(rsV+rs-s+1) (up to constant)", sp.cancel(top / (-r*((u+1)*Vv + 1)*(r*s*Vv + r*s - s + 1))).is_number, top)
# Newton polygon: all monomials Z^i W^k satisfy i - k <= 1, so no branch with Z->inf other than V finite
mons = sp.Poly(curve, Z, W).monoms()
check("Newton polygon: every monomial Z^i W^k of the curve has i - k <= 1 (no other Z->inf branch)", all(i - k <= 1 for i, k in mons), sorted(mons))
check("monomial Z^3 W^2 has coefficient -r^2 s (u+1)", sp.expand(sp.Poly(curve, Z, W).coeff_monomial(Z**3*W**2) + r**2*s*(u+1)) == 0)
# limits at the W=inf points: (W-1)/A2 -> 1/Y, and at Y=0, WY finite: check via Puiseux at test params
Zr0 = sp.Rational(1, 1)*(r-1)/r
# at Y=0 root of A: curve/W^2 = A + B/W + C/W^2; Z(omega) with omega=1/W; Z'(0) = -B/A' at Z0
dA = sp.diff(Acoef, Z)
lim_WY = sp.simplify((r * (-Bcoef / dA)).subs(Z, Zr0))
info("limit of W*Y at the W=inf, Y=0 point (= r * dZ/domega)", lim_WY)
check("that limit is a real rational function of (r,s,u,h) (finite at generic params)", lim_WY.subs(spec).is_real and lim_WY.subs(spec).is_finite)
# genericity at (2,3,4,5)
Ap = sp.Poly(Acoef.subs(spec), Z)
check("A simple roots at (2,3,4,5)", sp.discriminant(Ap) != 0)
check("(5) simple roots at (2,3,4,5)", sp.discriminant(sp.Poly(top.subs(spec), Vv)) != 0)
check("gcd(A, C) = 1 at (2,3,4,5)", sp.gcd(Acoef.subs(spec), C0coef.subs(spec)) .is_number)
# spinors nonzero off the divisors: a spinor vanishes only where W*A2*K*L = 0
zero_ok = True
for P_ in A_sp + B_sp:
    sol = sp.solve([P_[0], P_[1]], [Z, W], dict=True)
    for so in sol:
        val = (W*A2s*Ks*Ls).subs(so)
        zero_ok = zero_ok and sp.simplify(val) == 0
check("every common zero of a spinor's two components lies on W*A2*K*L = 0", zero_ok)
# numeric: random parameters, puncture realness
import numpy as np
nrng = np.random.default_rng(7); bad = []
for _ in range(20):
    pr = {r: sp.Rational(int(nrng.integers(-9, 10)) or 1, int(nrng.integers(1, 5))), s: sp.Rational(int(nrng.integers(-9, 10)) or 2, int(nrng.integers(1, 5))),
          u: sp.Rational(int(nrng.integers(-9, 10)) or 3, int(nrng.integers(1, 5))), h: sp.Rational(int(nrng.integers(-9, 10)), int(nrng.integers(1, 5)))}
    for fac in (RL.subs(pr), RA.subs(pr), sp.Poly(C0coef.subs(pr), Z).as_expr(), top.subs(pr)):
        if fac == 0: continue
        var = Vv if fac.has(Vv) else Z
        rts = sp.Poly(fac, var).nroots() if sp.Poly(fac, var).degree() > 0 else []
        if any(abs(complex(x).imag) > 1e-9 for x in rts): bad.append((pr, fac))
check("all finite-pole and Z=inf projections are real at 20 random rational parameter points", not bad, bad[:2])
# ---- NUMERIC: the chart covers sampled realizable eight-distinct T_45 points ----
def mob(z0, z1, zinf):
    # M with M(z0)=0, M(z1)=-1, M(zinf)=inf
    k = -(z1 - zinf) / (z1 - z0)
    return lambda x: k * (x - z0) / (x - zinf)
nchk = 0; errs = []; minphi = np.inf
for seed in range(10):
    pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=20, margin=1e-3, pi_star=1)
    if pt is None: continue
    av, bv = pt["a"], pt["b"]
    Fv = L.F_values_np(av, bv, 1.0)
    Hn = (Fv - Fv[1]) / (Fv[7] - Fv[1])
    M_ = mob(-bv[1], -bv[2], -bv[3])
    z_, w_, A2_, A0_ = M_(av[1]), M_(av[3]), M_(av[2]), M_(av[0]); t_ = -M_(-bv[0])
    u_ = -Hn[3]; r_ = 1 + Hn[2] / u_; s_ = 1 / (Hn[12] - Hn[2])
    A2c = w_ + r_ * (z_ - w_) / (z_ + 1); tc = s_ * A2c / (s_ * A2c - w_ - 1)
    A0c = -tc * (z_ - w_ + u_ * z_ * (w_ + 1)) / (z_ - w_ - u_ * tc * (w_ + 1))
    imag_rsu = max(abs(complex(r_).imag), abs(complex(s_).imag), abs(complex(u_).imag), abs(complex(Hn[11]).imag))
    errs.append(max(abs(A2c - A2_), abs(tc - t_), abs(A0c - A0_), imag_rsu))
    nchk += 1
    cvec = np.linalg.solve(np.array([[float(x) for x in row] for row in L.P])[L.P_PIVOT_ROWS], pt["gaps"][L.P_PIVOT_ROWS])
    minphi = min(minphi, L.min_phi_over_SO4(cvec, restarts=6, seed=seed)[0])
check("NUMERIC: chart formulas for A2, t, A0 reproduce %d sampled realizable T_45 points; r,s,u,h real" % nchk, nchk > 0 and max(errs) < 1e-6, "max err %.2e" % max(errs))
info("NUMERIC min Phi over SO(4) on those points (6 restarts each)", "%.3g" % minphi)
print("SUMMARY failures:", FAIL); sys.exit(1 if FAIL else 0)
