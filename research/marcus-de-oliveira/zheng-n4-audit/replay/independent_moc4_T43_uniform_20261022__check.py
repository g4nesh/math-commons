"""Own-code replay for independent_moc4_T43_uniform_20261022 (paper case T_44). Audit batch C18.
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

# ================= independent_moc4_T43_uniform_20261022 (paper T_44) =================
import numpy as np
T = frozenset({1, 4, 5, 6, 7, 8})
check("contact set equals paper T_44 literally", T == L.CONTACT_SETS[44])
gens = L.generators(T)
check("generators(T) = {0,1,2,3,4,5,6,33,34,59}", gens == [0,1,2,3,4,5,6,33,34,59], gens)
check("only nonminor generator is 34", [j for j in gens if not L.F_IS_MINOR[j]] == [34])
NAMES = ["q","p","j","c","e","d","b","a"]; IDX = [3, 33, 1, 59, 5, 34, 6, 0]
GAPF = {"q": lambda v: v[14], "p": lambda v: v[13], "j": lambda v: v[16], "c": lambda v: v[0],
        "e": lambda v: v[10], "d": lambda v: v[3]-v[13], "b": lambda v: v[15], "a": lambda v: v[19]}
# contact nullspace of P
rowsT = [L.P[s] for s in sorted(T)]
NS = L.nullspace(rowsT)
check("contact nullspace has dimension 8", len(NS) == 8)
ok = True
for v in NS:
    vals = L.values(v)
    rep = vec_comb({IDX[i]: GAPF[NAMES[i]](vals) for i in range(8)})
    ok = ok and all(x == y for x, y in zip(L.values(rep), vals))
check("identity (1) with gap-coordinate coefficients holds on the whole contact nullspace", ok)
# gap identities used for d <= 0
co = sp.symbols("q p j c e d b a")
gv = L.values(vec_comb({IDX[i]: co[i] for i in range(8)}))
qq, pp, jj, cc, ee, dd, bb, aa = co
check("g3 = p+d, g9 = b+d, g18 = q+d+a, g22 = j+d+a",
      [sp.expand(gv[3]-(pp+dd)), sp.expand(gv[9]-(bb+dd)), sp.expand(gv[18]-(qq+dd+aa)), sp.expand(gv[22]-(jj+dd+aa))] == [0,0,0,0])
check("identity f2+f34 = f0+f6+f33", lin_identity([2,34],[0,6,33]))
check("identity f4+f34 = f1+f3+f6+f33", lin_identity([4,34],[1,3,6,33]))
x_, y_, m_ = sp.symbols("x y m")
lhs = vec_comb({34: -m_}); rhs = vec_comb({2: x_, 4: y_, 0: -x_, 1: -y_, 3: -y_, 6: -m_, 33: -m_})
check("-m f34 = x f2 + y f4 - x f0 - y f1 - y f3 - m f6 - m f33 when x+y=m",
      all(sp.expand((L.values(lhs)[s]-L.values(rhs)[s]).subs(y_, m_-x_)) == 0 for s in range(24)))
# core action
def maps_f(gi, j):
    t = L.g_index_map(gi); v = [None]*24
    for s in range(24): v[t[s]] = L.F_VALUES[j][s]
    return [k for k in range(82) if list(L.F_VALUES[k]) == v]
src = [59, 53, 55, 67, 54, 2]; dst = [59, 5, 6, 0, 34, 1]
acts = [gi for gi in L.G_ELEMENTS if all(maps_f(gi, s_) == [t_] for s_, t_ in zip(src, dst))]
check("some g in G maps (f59,f53,f55,f67,f54,f2) to (f59,f5,f6,f0,f34,f1)", len(acts) > 0, acts[:2])

# ---- cubic restriction ----
h = sp.symbols("h")
g = vec_comb({IDX[i]: co[i] for i in range(8)})
gh = list(g); gh[13] = gh[13] + h
Ch = sp.Poly(sp.expand(L.cubic_C(gh)), h)
check("C(g + h kappa) is affine in h", Ch.degree() <= 1)
Cm = sp.expand(Ch.coeff_monomial(1))
R = sp.Poly(sp.expand(Cm - pp*(ee*(dd-cc)+dd**2)), *co)
check("(4): C - p{e(d-c)+d^2} has only nonnegative coefficients", all(x >= 0 for x in R.coeffs()), "terms=%d" % len(R.coeffs()))
check("(4): its d^3 coefficient is 1", R.coeff_monomial(dd**3) == 1)

# ---- radial factor G9 ----
def cost(vals):  # vals in NAMES order
    return vec_comb({IDX[i]: Fr(vals[i]) for i in range(8)})
BASE = [0,0,0,0,0,1,0,0]
G2b, C0b, Qb = radial_G2(cost(BASE))
info("base point d=1: -Q^2 detH/4, C, Q", (G2b, C0b, Qb))
def G_line(P0, P1, npts, g0sign_val=None):
    ys = []
    for t in range(npts):
        v = [Fr(P0[i]) + t*(Fr(P1[i]) - Fr(P0[i])) for i in range(8)]
        ys.append(radial_G2(cost(v))[0])
    S = interp(list(range(npts)), ys)
    return S, S.sqrt()
def k_coeffs(pt):
    """exact (k0,k1,k2) of G9 in c at pt (c entry ignored); sign normalised by G9(base)=+1"""
    P0 = list(pt); P0[3] = 0
    S, Rt = G_line(BASE, P0, 21)
    if Rt is None: return None
    sgn = 1 if Rt(0) > 0 else -1
    G0 = sgn * Rt(1) / abs(Rt(0))  # G9 at c=0, normalised
    P1 = list(P0); P1[3] = 1
    S2, R2 = G_line(P0, P1, 7)
    if R2 is None or R2.degree() > 2: return None
    sg2 = 1 if (R2(0) > 0) == (G0 > 0) else -1
    cf = [sg2 * R2.coeffs()[i] / abs(Rt(0)) * 0 + sg2 * R2.coeffs()[i] for i in range(R2.degree()+1)]
    scale = G0 / cf[0]
    return [Fr(int((x*scale).p), int((x*scale).q)) for x in cf] + [Fr(0)]*(3-len(cf))
def k0_formula(q, p, j, e, b, a):
    K = (a*b-1)*e - (a+1)*(b+1)
    D = a*a*b*j + a*a*e*q + a*a*j + a*b*j*q + 2*a*b*j + a*e*j*q + a*e*j + a*e*q + a*j*q + 2*a*j + b*j*q + b*j + e*j*q + e*j + j*q + j
    H4 = (a+q+1)*(-K) + D
    P2 = (a+1)*(b+1) + q*(b+e+1) + e
    return (b+1)*(p+1)**2*P2*H4
rng = random.Random(44)
def rnd(): return Fr(round(2.0 ** rng.uniform(-3, 3) * 8), 8)
okk0 = 0; okq = 0; N0 = 6
for it in range(N0):
    pt = [rnd() for _ in range(8)]; pt[5] = Fr(1)
    if it == 1: pt[2] = Fr(0)
    if it == 2: pt[4] = Fr(0)
    ks = k_coeffs(pt)
    if ks is None: continue
    okq += 1
    if ks[0] == k0_formula(pt[0], pt[1], pt[2], pt[4], pt[6], pt[7]): okk0 += 1
check("G9 is a polynomial of degree <= 2 in c (exact square root along c-lines) at %d points" % N0, okq == N0)
check("k0 = (b+1)(p+1)^2 P2 H4 with H4 = (a+q+1)(-K)+D, exact at %d rational points" % N0, okk0 == N0)
# (10): e = j = 0
ok10 = 0
for it in range(3):
    pt = [rnd() for _ in range(8)]; pt[5] = Fr(1); pt[2] = Fr(0); pt[4] = Fr(0)
    ks = k_coeffs(pt)
    q_, p_, _, _, _, _, b_, a_ = pt
    want = [(b_+1)**2*(p_+1)*(a_+q_+1)**2*(a_+1)*(b_+1)*(p_+1), (b_+1)**2*(p_+1)*(a_+q_+1)**2*(a_+b_+1)*(p_+1), (b_+1)**2*(p_+1)*(a_+q_+1)**2*p_]
    if ks == want: ok10 += 1
    else: info("(10) mismatch", (ks, want))
check("(10) factorization of G9 at e=j=0, exact at 3 points", ok10 == 3)

# ---- Bernstein tables (8) evaluated exactly at sampled chart points ----
def chart_points(row, A, Tt, V, U):
    if row == "I-":  a=A; b=Tt/(A*(1+Tt)); e=V; j=V+U; J=1/j
    elif row == "I+":
        a=A; b=(1+Tt)/A; Em=(a+1)*(b+1)/(a*b-1); e=Em*V/(1+V); j=e+U; J=1/j
    elif row == "II":
        a=A; b=(1+Tt)/A; Em=(a+1)*(b+1)/(a*b-1); e=Em*V/(1+V); j=e*U/(1+U); J=1/e-1/Em
    elif row == "IIIa-":
        b0=Tt/(1+Tt); a=b0*V/(1+V); b=b0; e=A; j=A*U/((1+A)*(1+U)); J=1+1/e
    elif row == "IIIa+":
        b0=Tt/(1+Tt); a=b0*V/(1+V); b=b0; e=A; j=A/(1+A)+A*A*U/((1+A)*(1+U)); J=1/j
    elif row == "IIIb+":
        b1=1+Tt; a=V/(b1*(1+V)); b=b1; e=A; j=A*b1/(A+b1)+A*A*U/((A+b1)*(1+U)); J=1/j
    elif row == "IIIb-":
        b1=1+Tt; a=V*V/(b1*(1+V)**2); b=b1; e=A; j=A*b1*U/((A+b1)*(1+U)); w=V/(1+V); J=1/e+b*(1-w)**2/(b-w)**2
    return a, b, e, j, J
nb = 0; bad = []
rows = ["I-","I+","II","IIIa-","IIIa+","IIIb+","IIIb-"]
for row in rows:
    for orient in (0, 1):
        for it in range(2):
            A, Tt, V, U, q_, p_ = [rnd() for _ in range(6)]
            a_, b_, e_, j_, J = chart_points(row, A, Tt, V, U)
            if orient: a_, b_ = b_, a_
            pt = [q_, p_, j_, 0, e_, 1, b_, a_]
            ks = k_coeffs(pt)
            if ks is None: bad.append((row, "nosqrt")); continue
            k0, k1, k2 = ks
            B = [k0, k0 + J*k1/2, k0 + J*k1 + J*J*k2]
            nb += 1
            if min(B) < 0: bad.append((row, orient, [float(x) for x in B]))
check("Bernstein coefficients (8) nonnegative at %d exact chart points (both leaf orientations; 2 per row/orientation)" % nb, not bad, bad)

# ---- realizable sampling ----
Bmat = np.array([[float(v) for v in L.F_VALUES[j]] for j in IDX]).T
def safe(q, p, j, c, e, b, a):
    if c*j >= 1: return True
    if e > 0 and c >= 1 + 1/e: return True
    if (a*b-1)*e - (a+1)*(b+1) >= 0: return True
    if a*b > 1 and e > 0 and c >= 1/e - (a*b-1)/((a+1)*(b+1)): return True
    for (x, y) in ((a, b), (b, a)):
        if x*y < 1 and y >= x and y >= 1 and e > 0:
            w = np.sqrt(x*y)
            if c >= 1/e + y*(1-w)**2/(y-w)**2: return True
    return False
nfound = 0; nform = 0; nhard = 0; minphi = np.inf
for seed in range(12):
    pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=20, margin=1e-3)
    if pt is None: continue
    nfound += 1
    gv = pt["gaps"]; cof = np.array([GAPF[n](gv) for n in NAMES])
    if cof[5] <= 1e-9: continue
    nform += 1
    cof = cof / cof[5]
    q, p, j, c, e, d, b, a = cof
    if not safe(q, p, j, c, e, b, a): nhard += 1
    mval, _ = L.min_phi_over_SO4(np.linalg.solve(np.array([[float(x) for x in row] for row in L.P])[L.P_PIVOT_ROWS], gv[L.P_PIVOT_ROWS]), restarts=6, seed=seed)
    minphi = min(minphi, mval)
info("NUMERIC realizable sampling (seeds 0-11, margin 1e-3)", "found=%d d>0=%d hard=%d min Phi over SO4 = %.3g" % (nfound, nform, nhard, minphi))
check("NUMERIC: no sampled realizable T_44 cost with d>0 fails every coefficient bound", nhard == 0)
print("SUMMARY failures:", FAIL); sys.exit(1 if FAIL else 0)
