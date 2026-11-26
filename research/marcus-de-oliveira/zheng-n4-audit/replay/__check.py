"""Own-code replay for  (branch document for paper case T_45). Audit batch C18.
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
print("SUMMARY failures:", FAIL); sys.exit(1 if FAIL else 0)
