"""Own-code replay for blind_moc4_T44_coherent_chain_capacity_20261022 (branch document for paper case T_45). Audit batch C18.
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

# ================= blind_moc4_T44_coherent_chain_capacity_20261022 (branch for paper T_45) =================
import numpy as np
T = frozenset({1, 2, 3, 7, 11, 12})
NAMES = ["q","p","j","c","e","d","b","a"]; IDX = [2, 1, 53, 12, 8, 18, 7, 15]
check("cost generators are the T_45 generator set minus f0", sorted(IDX) == sorted(set(L.generators(T)) - {0}))
Pf = np.array([[float(x) for x in row] for row in L.P])
def cvec(vals):
    v = np.zeros(14)
    for i, j in enumerate(IDX): v += float(vals[i]) * np.array([float(x) for x in L.F_VECS[j]])
    return v
# exact: derivative claim and t=0 reduction
z_, J_, b_ = sp.symbols("z J b", positive=True)
f = (1 + z_**2/b_) / (J_ + z_)**2
check("d/dz (1+z^2/b)/(sqrt j + z)^2 = 2(z sqrt j - b)/(b (sqrt j + z)^3)", sp.simplify(sp.diff(f, z_) - 2*(z_*J_ - b_)/(b_*(J_ + z_)**3)) == 0)
r, k, v, w = sp.symbols("r k v w", real=True)
x = (v*(2*r - k) - 1)/(r - k); y = k*(1 - r*v)/(w*(r - k))
jj, cc = sp.symbols("j c", positive=True)
sub = {v: 0, r: sp.sqrt(jj), k: -z_, w: sp.sqrt(b_)}
check("Prop 2 from Prop 1 at v=0, r=sqrt j, k=-z, w=sqrt b: p >= (1+z^2/b)/(sqrt j+z)^2", sp.simplify((x**2 + y**2).subs(sub) - (1 + z_**2/b_)/(sp.sqrt(jj) + z_)**2) == 0)
zs = sp.sqrt(cc)
check("Prop 2 branch jc<=b^2 value (b+c)/(b(sqrt j+sqrt c)^2) at z=sqrt c", sp.simplify(f.subs({z_: zs, J_: sp.sqrt(jj)}) - (b_ + cc)/(b_*(sp.sqrt(jj) + zs)**2)) == 0)
zb = b_/sp.sqrt(jj)
check("Prop 2 branch jc>=b^2 value 1/(j+b) at z=b/sqrt j", sp.simplify(f.subs({z_: zb, J_: sp.sqrt(jj)}) - 1/(jj + b_)) == 0)
t_ = sp.symbols("t", positive=True); m_ = sp.symbols("m", positive=True); Jr, Cr = sp.symbols("Jr Cr", positive=True)
sub1 = {v: t_/sp.sqrt(1 - t_**2), r: Jr*sp.sqrt(1 - t_**2), k: -Cr*sp.sqrt(1 - t_**2), w: sp.sqrt(b_)}
lhs1 = sp.simplify((x**2 + y**2).subs(sub1))
rhs1 = (((2*Jr + Cr)*t_ - 1)**2/(1 - t_**2) + Cr**2*(1 - Jr*t_)**2/b_)/(Jr + Cr)**2
check("scalar test (1) equals x^2+y^2 under the stated substitution (with c = C^2)", sp.simplify(lhs1 - rhs1) == 0, sp.simplify(lhs1 - rhs1))
check("under that substitution j = r^2(1+v^2) = J^2 and c = k^2(1+v^2) = C^2", sp.simplify((r**2*(1 + v**2)).subs(sub1) - Jr**2) == 0 and sp.simplify((k**2*(1 + v**2)).subs(sub1) - Cr**2) == 0)
check("and a,e >= v^2 = t^2/(1-t^2) <= m  iff  t <= sqrt(m/(1+m))", sp.simplify((t_**2/(1 - t_**2)).subs(t_, sp.sqrt(m_/(1 + m_))) - m_) == 0)
# NUMERIC: Prop 1 at its thresholds
rng = np.random.default_rng(11); worst = np.inf; n1 = 0
for it in range(25):
    rr, kk, vv = rng.normal(size=3) * 1.2; ww = abs(rng.normal()) + 0.05
    if abs(rr - kk) < 0.1: continue
    xx = (vv*(2*rr - kk) - 1)/(rr - kk); yy = kk*(1 - rr*vv)/(ww*(rr - kk))
    vals = [0, xx**2 + yy**2, rr**2*(1 + vv**2), kk**2*(1 + vv**2), vv**2, 1, ww**2, vv**2]
    mval, _ = L.min_phi_over_SO4(cvec(vals), restarts=8, seed=it)
    scale = sum(vals); worst = min(worst, mval/scale); n1 += 1
check("NUMERIC Prop 1: min Phi over SO(4) >= -1e-9 (relative) at %d random threshold costs (8 restarts each)" % n1, worst > -1e-9, "worst min/sum(coeffs) = %.3e" % worst)
# NUMERIC: Prop 2 at threshold
worst2 = np.inf
for it in range(20):
    bb, jv, cv = np.exp(rng.normal(size=3))
    pth = (bb + cv)/(bb*(np.sqrt(jv) + np.sqrt(cv))**2) if jv*cv <= bb*bb else 1/(jv + bb)
    vals = [0, pth, jv, cv, 0, 1, bb, 0]
    mval, _ = L.min_phi_over_SO4(cvec(vals), restarts=8, seed=100 + it)
    worst2 = min(worst2, mval/sum(vals))
check("NUMERIC Prop 2: min Phi >= -1e-9 (relative) at 20 random threshold costs", worst2 > -1e-9, "worst %.3e" % worst2)
# sharpness probe: 3% below the Prop 2 threshold (information only)
below = []
for it in range(6):
    bb, jv, cv = np.exp(rng.normal(size=3))
    pth = (bb + cv)/(bb*(np.sqrt(jv) + np.sqrt(cv))**2) if jv*cv <= bb*bb else 1/(jv + bb)
    vals = [0, 0.97*pth, jv, cv, 0, 1, bb, 0]
    below.append(round(L.min_phi_over_SO4(cvec(vals), restarts=8, seed=200 + it)[0], 5))
info("NUMERIC Prop 2 at 0.97 x threshold, min Phi (negative = threshold is sharp there)", below)
mval, _ = L.min_phi_over_SO4(cvec([0, 0.5, 1, 1, 0, 1, 1, 0]), restarts=20, seed=3)
check("NUMERIC example (1/2) f1 + f53 + f12 + f7 + f18 >= 0 on SO(4)", mval > -1e-9, "min %.3e" % mval)
print("SUMMARY failures:", FAIL); sys.exit(1 if FAIL else 0)
