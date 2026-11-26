"""Own-code replay for blind_moc4_T44_conic_branch_exclusion_20261022 (branch document for paper case T_45). Audit batch C18.
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

# ================= blind_moc4_T44_conic_branch_exclusion_20261022 (branch for paper T_45) =================
import numpy as np
T = frozenset({1, 2, 3, 7, 11, 12})
NAMES = ["q","p","j","c","e","d","b","a"]; IDX = [2, 1, 53, 12, 8, 18, 7, 15]
co = sp.symbols("q p j c e d b a", positive=True); q, p, j, c, e, d, b, a = co
hh = sp.symbols("h")
g = vec_comb({IDX[i]: co[i] for i in range(8)})
gh = list(g); gh[13] = gh[13] + hh
Ch = sp.Poly(sp.expand(L.cubic_C(gh)), hh)
check("C(g + h kappa) affine in h", Ch.degree() <= 1)
Cm = sp.expand(Ch.coeff_monomial(1))
ex = sp.Poly(sp.expand((a + e)*Cm - a*e*(d*(a + e) + d**2 - j*(a + e))), *co)
check("(a+e)C - ae{d(a+e)+d^2-j(a+e)} has only nonnegative coefficients", all(x >= 0 for x in ex.coeffs()), "terms=%d" % len(ex.coeffs()))
check("... and contains (a+e)d^3", ex.coeff_monomial(a*d**3) == 1 and ex.coeff_monomial(e*d**3) == 1)
# pure algebra of Lemma 2
Wd = c + d + b
R = (d + a)*(j - c)*p + j*c*d + j*(d + b)*(d - a) + c*b*(d + a)
Xp = d*j*(b - a)*Wd*(c*(d + a) - j*a)
Phi = R**2 - 4*Xp
X = p*c - a*Wd; Dd = c*b*(d + a) - a*j*(d + b); Yy = -p*Dd + a*j*d*Wd
check("DR = D^2 + dj(b-a)W{c(d+a)-ja} - (d+a)(j-c)Y", sp.expand(Dd*R - (Dd**2 + Xp - (d + a)*(j - c)*Yy)) == 0)
check("on Y=0 (so DR = D^2 + X'): D^2 Phi = (D^2 + X')^2 - 4 D^2 X' = (D^2 - X')^2", sp.expand((Dd**2 + Xp)**2 - 4*Dd**2*Xp - (Dd**2 - Xp)**2) == 0)
check("a=0: Y = -pcbd", sp.expand(Yy.subs(a, 0) + p*c*b*d) == 0)
Zz = p*j - p*c - d*Wd; E1 = (p - b)*Zz - j*d*(p + Wd); E0 = p*(d + b)*Zz
M = p*c + d*Wd; L1 = M**2 + d*j*Wd**2 + (d + b)*Wd*Zz; L0 = d*L1 + p*c*(d + b)*Zz
T3 = p*Zz - d*(p*c + (j + d)*Wd)
check("L1 + W E1 = p{pjW - (d+b)(pc+dW)}", sp.expand(L1 + Wd*E1 - p*(p*j*Wd - (d + b)*(p*c + d*Wd))) == 0)
# sign claims: pj < d^2 => Z<0, T3<0, p{pjW-(d+b)(pc+dW)}<0 ; then L1<0 needs E1>0
rng = np.random.default_rng(9); badsign = 0
for _ in range(2000):
    vals = dict(zip(["p","j","c","b","d"], np.exp(rng.normal(size=5))))
    if vals["p"]*vals["j"] >= vals["d"]**2: continue
    sb = {p: vals["p"], j: vals["j"], c: vals["c"], b: vals["b"], d: vals["d"]}
    zz = float(Zz.subs(sb)); t3 = float(T3.subs(sb)); lw = float((L1 + Wd*E1).subs(sb))
    if not (zz < 0 and t3 < 0 and lw < 0): badsign += 1
check("NUMERIC: pj<d^2 implies Z<0, T3<0, L1+W E1<0 (random samples)", badsign == 0)
# the substitution identity W^2 E1^2 G|_{a=pc/W, e=-E0/E1} = d(d+b) Z T3 (q L1 + L0)^2, up to a global constant
def costv(vals): return vec_comb({IDX[i]: Fr(vals[i]) for i in range(8)})
ratios = []
rr = random.Random(5)
for _ in range(6):
    pv, jv, cv, bv, dv, qv = [Fr(rr.randint(1, 12), rr.randint(1, 6)) for _ in range(6)]
    if pv*jv >= dv*dv: jv = dv*dv/(2*pv)
    Wv = cv + dv + bv; av = pv*cv/Wv
    sb = {p: pv, j: jv, c: cv, b: bv, d: dv}
    E0v = Fr(str(E0.subs(sb))); E1v = Fr(str(E1.subs(sb)))
    if E1v == 0: continue
    ev = -E0v/E1v
    if ev <= 0: continue   # only positive e is relevant
    G2v, C0v, Qv = radial_G2(costv([qv, pv, jv, cv, ev, dv, bv, av]))
    rhs = Fr(str((d*(d + b)*Zz*T3*(q*L1 + L0)**2).subs({**sb, q: qv}))) / (Wv**2 * E1v**2)
    ratios.append(Fr(int(G2v.p), int(G2v.q)) / rhs**2 if rhs != 0 else None)
info("ratio G_own^2 / RHS^2 at sampled points with a=pc/W, e=-E0/E1 > 0", [str(x) for x in ratios])
check("substitution identity holds up to one global constant at %d exact points" % len(ratios), len(ratios) > 0 and len(set(ratios)) == 1)
# realizable sampling: is any realizable T_45 cost in the chamber with Phi < 0 ?
nf = 0; nin = 0; nneg = 0
for seed in range(12):
    pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=20, margin=1e-3)
    if pt is None: continue
    nf += 1
    gv = pt["gaps"]
    cf = np.array([gv[21], gv[23], gv[0], gv[6], gv[9], gv[4] - gv[18], gv[20], gv[18]])
    if cf[5] <= 1e-9: continue
    cf = cf / cf[5]
    qv, pv, jv, cv, ev, dv, bv, av = cf
    if min(qv, pv, jv, cv, ev, bv) > 0 and pv*jv < 1 and jv < 1 + 1/(av + ev):
        nin += 1
        sb = {p: pv, j: jv, c: cv, b: bv, d: 1.0, a: av}
        if float(Phi.subs(sb)) < 0: nneg += 1
info("NUMERIC realizable T_45 samples", "found=%d in chamber=%d with Phi<0=%d" % (nf, nin, nneg))
check("NUMERIC: no sampled realizable T_45 cost lies in the chamber with Phi < 0", nneg == 0)
print("SUMMARY failures:", FAIL); sys.exit(1 if FAIL else 0)
