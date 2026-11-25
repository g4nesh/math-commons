"""Replay for "The order-four minor polytope and its negative-mass constant"
(root_complete_negative_mass_20261022).  Own code only; nothing from the supplement is read
or executed (their verifier and JSON certificate are NOT used).  One printed line per check;
exit status nonzero iff an EXACT check fails.

Checks
  1. EXACT: for exact Gaussian-rational unitary U, W_sigma = sgn(sigma) prod_i u_{i,sigma(i)} / det U
     satisfies sum_sigma W_sigma 1{sigma(I)=J} = |det U[I,J]|^2 for all I, J (|I| = 0..4), so the
     real parts are a common signed representation with sum 1
  2. EXACT: the facet value of h_{pi,eps} (paper Sect. 4.1, library h_vec) at the response
     equals (t_pi + eps c_pi)/2 with t, c defined from Q_V, V = U[:, pi], as in the document
  3. EXACT: 32 t^3 + 48 t^2 + 28 t - 1 has exactly one positive root, in the stated interval
  4. NUMERIC: nu(U) (linear program) equals rho(U) on random Haar unitaries and orthogonals,
     including points near the maximizer
  5. NUMERIC: max of rho over O(4) and over U(4) (multistart) versus tau*
The facet completeness (82 facets) is not re-proved here; it is the library's exact
Prop. 4.1 replay (moc4lib test_moc4lib.py).
"""
import itertools
import os
import random
import sys
from fractions import Fraction as Fr

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

FAIL = []


def check(name, ok, level="EXACT", info=""):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({info})" if info else ""))
    if not ok and level == "EXACT":
        FAIL.append(name)


QI = L.QI


def qinv(M):
    n_ = len(M)
    A = [list(r) + [QI(1) if i == j else QI(0) for j in range(n_)] for i, r in enumerate(M)]
    for c in range(n_):
        piv = next(r for r in range(c, n_) if A[r][c] != QI(0))
        A[c], A[piv] = A[piv], A[c]
        inv_p = QI(1) / A[c][c]
        A[c] = [x * inv_p for x in A[c]]
        for r in range(n_):
            if r != c and A[r][c] != QI(0):
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n_:] for r in A]


def qmatmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), QI(0)) for j in range(len(B[0]))] for i in range(len(A))]


def gaussian_unitary(rng, size=3):
    S = [[QI(0)] * 4 for _ in range(4)]
    for i in range(4):
        S[i][i] = QI(0, Fr(rng.randint(-size, size), rng.randint(1, size)))
        for j in range(i + 1, 4):
            v = QI(Fr(rng.randint(-size, size), rng.randint(1, size)), Fr(rng.randint(-size, size), rng.randint(1, size)))
            S[i][j] = v
            S[j][i] = -v.conjugate()
    I_ = [[QI(1) if i == j else QI(0) for j in range(4)] for i in range(4)]
    ImS = [[I_[i][j] - S[i][j] for j in range(4)] for i in range(4)]
    IpS = [[I_[i][j] + S[i][j] for j in range(4)] for i in range(4)]
    return qmatmul(ImS, qinv(IpS))


def qdet(M):
    k = len(M)
    if k == 0:
        return QI(1)
    tot = QI(0)
    for perm in itertools.permutations(range(k)):
        t = QI(L.sgn(perm) if k == 4 else sp_sign(perm))
        for i in range(k):
            t = t * M[i][perm[i]]
        tot = tot + t
    return tot


def sp_sign(perm):
    s = 1
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                s = -s
    return s


def sub(M, I, J):
    return [[M[i][j] for j in J] for i in I]


def absq(z):
    return z.re * z.re + z.im * z.im


rng = random.Random(22)
ok1 = ok_sum = ok_unit = ok2 = True
ALLSETS = [tuple(c) for k in range(5) for c in itertools.combinations(range(4), k)]
for trial in range(4):
    U = gaussian_unitary(rng)
    UU = qmatmul(U, [[U[j][i].conjugate() for j in range(4)] for i in range(4)])
    ok_unit &= all(UU[i][j] == (QI(1) if i == j else QI(0)) for i in range(4) for j in range(4))
    dU = qdet(U)
    W = []
    for s in L.PERMS:
        t = QI(L.sgn(s))
        for i in range(4):
            t = t * U[i][s[i]]
        W.append(t / dU)
    ok_sum &= sum(W, QI(0)) == QI(1)
    Q = {}
    for I in ALLSETS:
        for J in ALLSETS:
            if len(I) != len(J):
                continue
            q_ = absq(qdet(sub(U, I, J)))
            Q[(I, J)] = q_
            rep = sum((W[k] for k, s in enumerate(L.PERMS) if tuple(sorted(s[i] for i in I)) == J), QI(0))
            if not (rep.re == q_ and rep.im == 0):
                ok1 = False
    # 2. facet identity
    w = [x.re for x in W]
    for pi in L.PERMS:
        # V = U[:, pi]: column k of V is column pi(k) of U
        def QV(I1, J1):  # 1-based index strings
            Ii = tuple(int(ch) - 1 for ch in I1)
            Jj = tuple(sorted(pi[int(ch) - 1] for ch in J1))
            return Q[(Ii, Jj)]
        t_pi = sum(QV(str(i + 1), str(i + 1)) for i in range(4))
        c_pi = QV("12", "13") + QV("12", "24") - QV("13", "12") - QV("13", "34")
        for eps in (1, -1):
            hv = L.values(L.h_vec(pi, eps))
            lhs = sum(wk * hk for wk, hk in zip(w, hv))
            if lhs != (t_pi + eps * c_pi) / 2:
                ok2 = False
check("Cayley transforms are exactly unitary (4 Gaussian-rational U)", ok_unit)
check("sum_sigma W_sigma = 1", ok_sum)
check("sum_{sigma(I)=J} W_sigma = |det U[I,J]|^2 (imaginary part 0) for all I, J, |I| = 0..4", ok1)
check("facet value of h_{pi,eps} at the response = (t_pi + eps c_pi)/2 (doc's t, c), all 48", ok2)
# the minor vectors n_{I,J} span V_R, so (v(n_{I,J}))_{I,J} is an injective linear image of V_R^*:
# the doc's polytope in minor-response coordinates is affinely isomorphic to the paper's Q,
# and the paper's facet list (Prop. 4.1, library replay) transfers facet by facet
minor_rows = [list(L.F_VECS[j]) for j in range(82) if L.F_IS_MINOR[j]]
check("the 34 minor vectors span V_R (rank 14)", len(minor_rows) == 34 and L.rank(minor_rows) == 14)

# 3. the cubic
x = sp.symbols("x")
cub = 32 * x**3 + 48 * x**2 + 28 * x - 1
lo, hi = sp.Rational("0.033721124434"), sp.Rational("0.033721124435")
check("cubic changes sign on [0.033721124434, 0.033721124435]", cub.subs(x, lo) < 0 < cub.subs(x, hi))
check("cubic has exactly one positive real root (derivative 96x^2+96x+28 > 0)",
      sp.Poly(cub, x).count_roots(0, sp.oo) == 1)
tau = float(sp.nsolve(cub, x, 0.0337))

# 4./5. NUMERIC
import numpy as np  # noqa: E402
from scipy.optimize import linprog, minimize  # noqa: E402

rows_idx = [(I, J) for I in ALLSETS for J in ALLSETS if len(I) == len(J) and 0 < len(I) <= 2]
Mrep = np.array([[1.0 if tuple(sorted(s[i] for i in I)) == J else 0.0 for s in L.PERMS] for (I, J) in rows_idx])
Mrep = np.vstack([np.ones(24), Mrep])
H_VALS = np.array([[float(v) for v in L.values(L.h_vec(pi, e))] for pi in L.PERMS for e in (1, -1)])


def response(U):
    return np.array([1.0] + [abs(np.linalg.det(U[np.ix_(I, J)])) ** 2 for (I, J) in rows_idx])


def nu(U):
    y = response(U)
    # variables w (24, free), s (24, >= 0); min sum s; -w - s <= 0
    c = np.concatenate([np.zeros(24), np.ones(24)])
    A_ub = np.hstack([-np.eye(24), -np.eye(24)])
    A_eq = np.hstack([Mrep, np.zeros((Mrep.shape[0], 24))])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(24), A_eq=A_eq, b_eq=y,
                  bounds=[(None, None)] * 24 + [(0, None)] * 24, method="highs")
    assert res.status == 0
    return res.fun


def leibniz_w(U):
    w = []
    for s, sg in zip(L.PERMS, L.SGN):
        w.append(sg * np.prod([U[i, s[i]] for i in range(4)]))
    return np.real(np.array(w) / np.linalg.det(U))


def rho(U):
    w = leibniz_w(U)
    return max(0.0, float(np.max(-(H_VALS @ w) / 2)))


def haar_unitary(rs):
    Z = (rs.normal(size=(4, 4)) + 1j * rs.normal(size=(4, 4))) / np.sqrt(2)
    Qm, R = np.linalg.qr(Z)
    return Qm * (np.diag(R) / np.abs(np.diag(R)))


rs = np.random.default_rng(99)
worst = 0.0
npos = 0
NS = 300
for k in range(NS):
    U = haar_unitary(rs) if k % 2 == 0 else L.so4_from_params(rs.normal(size=6) * 2).astype(complex)
    a_, r_ = nu(U), rho(U)
    worst = max(worst, abs(a_ - r_))
    npos += r_ > 1e-9
check(f"NUMERIC nu(U) = rho(U) (LP vs formula) on {NS} samples (seed 99; half Haar U(4), half SO(4))",
      worst < 1e-7, "NUMERIC", info=f"max |nu-rho| = {worst:.2e}; rho > 0 in {npos} samples")


def neg_rho_O(theta, refl):
    O = L.so4_from_params(theta)
    if refl:
        O = O @ np.diag([-1.0, 1, 1, 1])
    return -rho(O.astype(complex))


best = 0.0
bestO = None
for k in range(100):
    th0 = rs.normal(size=6) * 2
    res = minimize(neg_rho_O, th0, args=(k % 2,), method="Nelder-Mead",
                   options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 6000})
    if -res.fun > best:
        best = -res.fun
        bestO = (res.x, k % 2)
check(f"NUMERIC max rho over O(4) (100 Nelder-Mead starts) vs tau* = {tau:.12f}",
      best <= tau + 1e-9, "NUMERIC", info=f"best = {best:.12f}, tau* - best = {tau - best:.2e}")
Ob = L.so4_from_params(bestO[0]) @ (np.diag([-1.0, 1, 1, 1]) if bestO[1] else np.eye(4))
check("NUMERIC nu = rho at the best orthogonal point", abs(nu(Ob.astype(complex)) - rho(Ob.astype(complex))) < 1e-7,
      "NUMERIC", info=f"nu = {nu(Ob.astype(complex)):.12f}")


def neg_rho_U(x):
    Hm = np.zeros((4, 4), complex)
    k = 0
    for i in range(4):
        Hm[i, i] = x[k]; k += 1
    for i in range(4):
        for j in range(i + 1, 4):
            Hm[i, j] = x[k] + 1j * x[k + 1]
            Hm[j, i] = np.conj(Hm[i, j]); k += 2
    from scipy.linalg import expm
    return -rho(expm(1j * Hm))


bestU = 0.0
for k in range(50):
    res = minimize(neg_rho_U, rs.normal(size=16) * 1.5, method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 8000})
    bestU = max(bestU, -res.fun)
check("NUMERIC max rho over U(4) (50 Nelder-Mead starts) does not exceed tau*",
      bestU <= tau + 1e-9, "NUMERIC", info=f"best = {bestU:.12f}")

print("SUMMARY", "ALL EXACT CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
