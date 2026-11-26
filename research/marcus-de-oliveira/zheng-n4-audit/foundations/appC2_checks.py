#!/usr/bin/env python3
"""Appendix C, part 2 (Sections C.4-C.7 and the proof of Theorem C.1) of
Zheng, "The Marcus-de Oliveira conjecture for n = 4": own-code replays.

Everything here is derived from the definitions printed in the paper
(matching basis (A.1), cubic (A.2), d(O) of Section C.4).  Nothing from the
supplement is imported or executed.  Optional argument: a section list, e.g.
    python appC2_checks.py core c9 c10 c8 lemC9 lemC10 c13 c11 crit attack
(default: all).  Deterministic seeds throughout.

Evidence levels printed with each check:
  EXACT   integer / rational / polynomial-identity computation
  NUMERIC floating point (sizes, seeds and tolerances printed)
"""
import itertools, sys, time, math, random
from fractions import Fraction as Fr

import numpy as np
import sympy as sp
import flint
from scipy.optimize import linprog, minimize, least_squares
from scipy.linalg import expm

T0 = time.time()
import os
SMALL = os.environ.get("APPC2_SMALL") == "1"     # smoke-test sizes only
def SZ(n):
    return max(1, n // 10) if SMALL else n
SECTIONS = sys.argv[1:] or ["core", "c9", "c10", "c8", "lemC9", "c8prop", "lemC10",
                            "c13", "c11", "crit", "attack"]
FAILS = []


def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


def check(cond, msg):
    log(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ---------------------------------------------------------------------------
# Core objects, from the paper's definitions
# ---------------------------------------------------------------------------
MATCH = [((1, 2), (3, 4), (5, 6), (7, 8)), ((1, 2), (3, 4), (5, 8), (6, 7)),
         ((1, 2), (3, 6), (4, 5), (7, 8)), ((1, 2), (3, 8), (4, 5), (6, 7)),
         ((1, 2), (3, 8), (4, 7), (5, 6)), ((1, 4), (2, 3), (5, 6), (7, 8)),
         ((1, 4), (2, 3), (5, 8), (6, 7)), ((1, 6), (2, 3), (4, 5), (7, 8)),
         ((1, 6), (2, 5), (3, 4), (7, 8)), ((1, 8), (2, 3), (4, 5), (6, 7)),
         ((1, 8), (2, 3), (4, 7), (5, 6)), ((1, 8), (2, 5), (3, 4), (6, 7)),
         ((1, 8), (2, 7), (3, 4), (5, 6)), ((1, 8), (2, 7), (3, 6), (4, 5))]
PERMS = list(itertools.permutations(range(4)))          # lexicographic, pi_0 = id
SUBSETS = {k: list(itertools.combinations(range(4), k)) for k in range(5)}
MINOR_KEYS = [(I, J) for k in range(5) for I in SUBSETS[k] for J in SUBSETS[k]]  # 70


def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def matching_vec(z):
    """m(z) = (B_0(z),...,B_13(z)), [ij] = z_i - z_j, z indexed 1..8 via z[0..7]."""
    out = []
    for M in MATCH:
        v = 1
        for (i, j) in M:
            v = v * (z[i - 1] - z[j - 1])
        out.append(v)
    return out


# multilinear monomial coefficient matrix of the basis: rows = 4-subsets S of 0..7
ZS = sp.symbols("z1:9")
MONS = list(itertools.combinations(range(8), 4))
MON_INDEX = {S: k for k, S in enumerate(MONS)}


def multilinear_coeffs(expr):
    poly = sp.Poly(sp.expand(expr), *ZS)
    vec = [sp.Integer(0)] * 70
    for mon, c in poly.terms():
        S = tuple(i for i in range(8) if mon[i] == 1)
        assert sum(mon) == 4 and len(S) == 4, "not multilinear of degree 4"
        vec[MON_INDEX[S]] = c
    return vec


BASIS_POLYS = matching_vec(ZS)
MMON = sp.Matrix([multilinear_coeffs(b) for b in BASIS_POLYS]).T   # 70 x 14
MMON_rank = MMON.rank()
PIV_ROWS = None


def express_in_basis(vec70):
    """solve MMON rho = vec70 exactly; returns (rho, consistent?)"""
    global PIV_ROWS
    if PIV_ROWS is None:
        _, piv = MMON.T.rref()
        PIV_ROWS = list(piv)
    A = MMON.extract(PIV_ROWS, list(range(14)))
    rho = A.LUsolve(sp.Matrix([vec70[r] for r in PIV_ROWS]))
    ok = all(sp.simplify(sum(MMON[r, j] * rho[j] for j in range(14)) - vec70[r]) == 0
             for r in range(70))
    return rho, ok


# P : (P x)_sigma = x(p_sigma), p_sigma = prod_i [i, 4+sigma(i)]
P_ROWS = []
for s in PERMS:
    expr = 1
    for i in range(4):
        expr *= (ZS[i] - ZS[4 + s[i]])
    rho, ok = express_in_basis(multilinear_coeffs(expr))
    assert ok
    P_ROWS.append(list(rho))
P = sp.Matrix(P_ROWS)                     # 24 x 14 integer
KAPPA = sp.Matrix([0] * 13 + [1])
MS = sp.symbols("m0:14")
m = MS
CUBIC = (m[13] * (m[0] * m[6] - m[1] * m[5])
         - m[0] * m[9] * (m[0] + m[1] + m[2] + m[3] + m[5] + m[7] + m[9]) - m[1] * m[5] * m[9]
         + m[10] * (-m[0] * m[3] - m[0] * m[9] + m[1] * m[2] + m[1] * m[7])
         + m[11] * (-m[0] * m[7] - m[0] * m[9] + m[2] * m[5] + m[3] * m[5])
         + m[12] * (-m[0] * m[9] + m[1] * m[7] - m[2] * m[6] - m[2] * m[9] + m[3] * m[5] + m[3] * m[7])
         + m[2] * m[10] * m[11] - m[3] * m[8] * m[10] - m[4] * m[7] * m[11] + m[4] * m[8] * m[9])
QQ = m[0] * m[6] - m[1] * m[5]
GRAD = [sp.diff(CUBIC, v) for v in MS]
HESS = sp.hessian(CUBIC, MS)
C_fn = sp.lambdify([MS], CUBIC, "sympy")


def Cval(vec):
    return sp.expand(CUBIC.subs(dict(zip(MS, list(vec))), simultaneous=True))


def qval(vec):
    return sp.expand(vec[0] * vec[6] - vec[1] * vec[5])


def Gval(vec):
    d = dict(zip(MS, list(vec)))
    return [sp.expand(g.subs(d, simultaneous=True)) for g in GRAD]


def Hval(vec):
    d = dict(zip(MS, list(vec)))
    return HESS.subs(d, simultaneous=True).applyfunc(sp.expand)


# minor vectors n_{I,J}: P n = indicator{sigma(I)=J}
ROWSEL = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18]
PSUB = P.extract(ROWSEL, list(range(14)))


def solve_P(target24):
    n = PSUB.LUsolve(sp.Matrix([target24[r] for r in ROWSEL]))
    ok = (P * n - sp.Matrix(target24)).is_zero_matrix
    return n, ok


NVEC = {}
for (I, J) in MINOR_KEYS:
    tgt = [1 if tuple(sorted(s[i] for i in I)) == J else 0 for s in PERMS]
    n, ok = solve_P(tgt)
    assert ok, (I, J)
    NVEC[(I, J)] = n


def minor(M, I, J):
    if len(I) == 0:
        return 1
    return M.extract(list(I), list(J)).det() if isinstance(M, sp.MatrixBase) else \
        np.linalg.det(M[np.ix_(I, J)])


def rho_exact(U):
    """rho(U) for an exact (sympy, real or complex) matrix U, via (C.3) coefficients."""
    vec = [0] * 70
    for (I, J) in MINOR_KEYS:
        Ic = [i for i in range(4) if i not in I]
        S = tuple(sorted(Ic + [4 + j for j in J]))
        mi = minor(U, I, J)
        vec[MON_INDEX[S]] = (-1) ** len(J) * sp.expand(mi * sp.conjugate(mi))
    rho, ok = express_in_basis(vec)
    return rho, ok


def rho_poly_real(O):
    """holomorphic version on (complex) orthogonal: |det|^2 -> det^2."""
    vec = [0] * 70
    for (I, J) in MINOR_KEYS:
        Ic = [i for i in range(4) if i not in I]
        S = tuple(sorted(Ic + [4 + j for j in J]))
        mi = minor(O, I, J)
        vec[MON_INDEX[S]] = (-1) ** len(J) * sp.expand(mi * mi)
    rho, ok = express_in_basis(vec)
    return rho, ok


def cayley(S):
    I4 = sp.eye(S.shape[0])
    return (I4 - S) * (I4 + S).inv()


def rand_skew(rng, n=4, H=5):
    S = sp.zeros(n)
    for i in range(n):
        for j in range(i + 1, n):
            v = sp.Rational(rng.randint(-H, H), rng.randint(1, H))
            S[i, j] = v
            S[j, i] = -v
    return S


def d_of(O):
    """d(O) of Section C.4 (0-indexed rows 1,2,3; columns 0,1,2)."""
    pairs = [(1, 2), (1, 3), (2, 3)]
    return sp.Matrix([[O[r, k] * O[s, k] for k in range(3)] for (r, s) in pairs]).det()


def r_of(R):
    """r(R) = det(1/R_ij) prod R_ij, cleared (sum over perms of the complementary products)."""
    n = R.shape[0]
    tot = 0
    for p in itertools.permutations(range(n)):
        prod = 1
        for i in range(n):
            for j in range(n):
                if j != p[i]:
                    prod *= R[i, j]
        tot += perm_sign(p) * prod
    return sp.expand(tot)


# numeric helpers ------------------------------------------------------------
PERM_ARR = np.array(PERMS)
N_np = None


def sigma_points(a, b):
    return np.array([np.prod([a[i] + b[s[i]] for i in range(4)]) for s in PERMS])


def detfun(a, b, U):
    return np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)


def run_core():
    log("== core: basis, P, kappa, cubic, minor vectors, (C.3), (C.5), (C.6), Moebius")
    check(MMON_rank == 14, "EXACT the 14 matching polynomials (A.1) are linearly independent")
    check(PSUB.det() == 1, "EXACT 14x14 submatrix of P with rows (A.1-list) has det 1 (paper App. A)")
    check(P * KAPPA == sp.ones(24, 1), "EXACT P kappa = 1 for kappa = e13 (3.2)")
    check(sp.expand(CUBIC.subs(m[13], m[13] + sp.Symbol('h')) - CUBIC - sp.Symbol('h') * QQ) == 0,
          "EXACT (C.5) C(u + h kappa) = C(u) + h q(u), q = u0 u6 - u1 u5")
    # supplement's printed expansion of the same cubic (one_circle.pdf, App. A), transcribed
    y = MS
    Csupp = (-y[0]**2*y[9] - y[0]*y[1]*y[9] - y[0]*y[10]*y[3] - y[0]*y[10]*y[9] - y[0]*y[11]*y[7]
             - y[0]*y[11]*y[9] - y[0]*y[12]*y[9] + y[0]*y[13]*y[6] - y[0]*y[2]*y[9] - y[0]*y[3]*y[9]
             - y[0]*y[5]*y[9] - y[0]*y[7]*y[9] - y[0]*y[9]**2 + y[1]*y[10]*y[2] + y[1]*y[10]*y[7]
             + y[1]*y[12]*y[7] - y[1]*y[13]*y[5] - y[1]*y[5]*y[9] + y[10]*y[11]*y[2] - y[10]*y[3]*y[8]
             + y[11]*y[2]*y[5] + y[11]*y[3]*y[5] - y[11]*y[4]*y[7] - y[12]*y[2]*y[6] - y[12]*y[2]*y[9]
             + y[12]*y[3]*y[5] + y[12]*y[3]*y[7] + y[4]*y[8]*y[9])
    check(sp.expand(Csupp - CUBIC) == 0, "EXACT paper (A.2) equals the supplement's expanded 28-term cubic")
    # X subset Sing(C): gradient vanishes on matching vectors (symbolic in z)
    mz = matching_vec(ZS)
    dz = dict(zip(MS, mz))
    check(all(sp.expand(g.subs(dz, simultaneous=True)) == 0 for g in GRAD),
          "EXACT grad C vanishes identically on m(z) (X in Sing C), 8 symbolic parameters")
    # (C.3) for exact orthogonal and exact unitary samples; <rho,kappa>=1; <rho,n_IJ>=|det|^2
    rng = random.Random(1)
    ok_all = True
    for t in range(4):
        O = cayley(rand_skew(rng))
        rho, ok = rho_exact(O)
        ok_all &= ok and rho[13] == 1
        ok_all &= all(sp.simplify(sum(rho[j] * NVEC[(I, J)][j] for j in range(14)) - minor(O, I, J)**2) == 0
                      for (I, J) in MINOR_KEYS)
        # determinant identity at a random rational spectrum
        a = [sp.Rational(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(4)]
        b = [sp.Rational(rng.randint(-9, 9), rng.randint(1, 4)) + sp.I * rng.randint(-5, 5) for _ in range(4)]
        lhs = sp.expand((sp.diag(*a) + O * sp.diag(*b) * O.T).det())
        mv = matching_vec(a + [-x for x in b])
        rhs = sp.expand(sum(rho[j] * mv[j] for j in range(14)))
        ok_all &= sp.simplify(lhs - rhs) == 0
    # one exact complex unitary: Cayley of a skew-Hermitian Gaussian-rational matrix
    S = sp.Matrix([[0, sp.Rational(1, 2) + sp.I / 3, sp.Rational(-1, 3), sp.I],
                   [-sp.Rational(1, 2) + sp.I / 3, sp.I / 2, 2, sp.Rational(1, 5) - sp.I],
                   [sp.Rational(1, 3), -2, 0, sp.Rational(3, 7) * sp.I],
                   [sp.I, -sp.Rational(1, 5) - sp.I, sp.Rational(3, 7) * sp.I, -sp.I]])
    assert (S + S.H).is_zero_matrix
    U = (sp.eye(4) - S) * (sp.eye(4) + S).inv()
    U = U.applyfunc(sp.nsimplify)
    check((U * U.H - sp.eye(4)).applyfunc(sp.simplify).is_zero_matrix, "EXACT test unitary is unitary")
    rho, ok = rho_exact(U)
    ok_all &= ok and sp.simplify(rho[13] - 1) == 0
    a = [0, 1, 3, 7]
    b = [2 + sp.I, -1 + 3 * sp.I, 5 - 2 * sp.I, 4 + 4 * sp.I]
    lhs = sp.expand((sp.diag(*a) + U * sp.diag(*b) * U.H).det())
    mv = matching_vec(a + [-x for x in b])
    ok_all &= sp.simplify(lhs - sum(rho[j] * mv[j] for j in range(14))) == 0
    check(ok_all, "EXACT (C.3) on 4 rational orthogonal + 1 Gaussian-rational unitary: det(Da+U Db U*) "
                  "= <rho(U), m(a,-b)>, <rho,kappa>=1, <rho,n_IJ>=|det U[I,J]|^2 (all 70)")
    # (C.6): q(Re(theta m(a,-b))) = 0 for real a, symbolic b, theta
    ar = sp.symbols("a1:5", real=True)
    br = sp.symbols("br1:5", real=True)
    bi = sp.symbols("bi1:5", real=True)
    tr, ti = sp.symbols("tr ti", real=True)
    z = list(ar) + [-(br[j] + sp.I * bi[j]) for j in range(4)]
    mv = [sp.expand((tr + sp.I * ti) * x) for x in matching_vec(z)]
    yv = [sp.re(x) for x in mv]
    check(sp.expand(qval(yv)) == 0, "EXACT (C.6) q(Re(theta m(a,-b))) = 0 for real a (17 symbolic real params)")
    # Moebius scaling of matching vectors
    al, be, ga, de = sp.symbols("alpha beta gamma delta")
    f = [(al * x + be) / (ga * x + de) for x in ZS]
    fac = (al * de - be * ga) ** 4 / sp.prod([ga * x + de for x in ZS])
    mf = matching_vec(f)
    check(all(sp.simplify(mf[j] - fac * mz[j]) == 0 for j in range(14)),
          "EXACT Moebius: m(f(z)) = (ad-bc)^4/prod(gamma z_j+delta) m(z) (proof of Thm C.1, p.57)")


if "core" in SECTIONS:
    run_core()


# ---------------------------------------------------------------------------
# (C.9): reciprocal cofactor identities, quaternion parametrisation (flint)
# ---------------------------------------------------------------------------
def quat_T(ctx):
    a, b, c, d, e, f = ctx.gens()
    one = ctx.from_dict({(0,) * 6: 1})
    w, x, y, zz = one, a, b, c
    L = [[w, -x, -y, -zz], [x, w, -zz, y], [y, zz, w, -x], [zz, -y, x, w]]
    w, x, y, zz = one, d, e, f
    R = [[w, -x, -y, -zz], [x, w, zz, -y], [y, -zz, w, x], [zz, y, -x, w]]
    T = [[sum((L[i][k] * R[k][j] for k in range(4)), ctx.from_dict({})) for j in range(4)] for i in range(4)]
    n = (one + a * a + b * b + c * c) * (one + d * d + e * e + f * f)
    return T, n


def r_flint(R, ctx):
    n = len(R)
    tot = ctx.from_dict({})
    for p in itertools.permutations(range(n)):
        prod = ctx.from_dict({(0,) * 6: perm_sign(p)})
        for i in range(n):
            for j in range(n):
                if j != p[i]:
                    prod = prod * R[i][j]
        tot = tot + prod
    return tot


def det3_flint(M, ctx):
    tot = ctx.from_dict({})
    for p in itertools.permutations(range(3)):
        tot = tot + perm_sign(p) * (M[0][p[0]] * M[1][p[1]] * M[2][p[2]])
    return tot


def run_c9():
    log("== (C.9): (-1)^(i+j) r(O_ij^) = d(O) and det(1/O) prod O = d(O)^2 on SO(4,C) (quaternion chart)")
    ctx = flint.fmpz_mpoly_ctx.get(('a', 'b', 'c', 'd', 'e', 'f'), 'degrevlex')
    T, n = quat_T(ctx)
    TT_ok = all((sum((T[i][k] * T[j][k] for k in range(4)), ctx.from_dict({}))
                 - (n if i == j else ctx.from_dict({}))).is_zero() for i in range(4) for j in range(4))
    check(TT_ok, "EXACT T = L(1,a,b,c) R(1,d,e,f) satisfies T T^t = n I")
    pairs = [(1, 2), (1, 3), (2, 3)]
    D = det3_flint([[T[r][k] * T[s][k] for k in range(3)] for (r, s) in pairs], ctx)
    log(f"     d(T) has total degree {D.total_degree()} and {len(list(D.terms()))} terms")
    ok = True
    for i in range(4):
        for j in range(4):
            sub = [[T[r][c] for c in range(4) if c != j] for r in range(4) if r != i]
            rv = r_flint(sub, ctx)
            ok &= ((-1) ** (i + j) * rv - D).is_zero()
    check(ok, "EXACT (C.9a) all 16 signed reciprocal 3x3 cofactors equal d (polynomial identity in 6 params)")
    t1 = time.time()
    r4 = r_flint(T, ctx)
    check((r4 - D * D).is_zero(), f"EXACT (C.9b) det(1/T) prod T = d(T)^2 (degree {r4.total_degree()}, "
                                  f"{time.time()-t1:.1f}s)")
    # sign-flip invariance (so (C.9) also holds on O(4)^- = D SO(4))
    Os = sp.Matrix(4, 4, sp.symbols("o0:16"))
    Dflip = sp.diag(1, 1, 1, -1)
    ok2 = True
    for Dm in [sp.diag(-1, 1, 1, 1), sp.diag(1, -1, 1, 1), sp.diag(1, 1, 1, -1)]:
        ok2 &= sp.expand(d_of(Dm * Os) - d_of(Os)) == 0 and sp.expand(d_of(Os * Dm) - d_of(Os)) == 0
    for (i, j) in [(0, 0), (1, 2), (3, 3)]:
        R = Os.copy()
        R.row_del(i)
        R.col_del(j)
        ok2 &= sp.expand(r_of(sp.diag(-1, 1, 1) * R) - r_of(R)) == 0
    check(ok2, "EXACT d and r are invariant under any row/column sign flip, so (C.9) holds on all of O(4,C)")
    # transposition / permutations: d(O^T), d(P O) agree with d(O) up to sign on rational samples
    rng = random.Random(9)
    sgns = set()
    for t in range(6):
        O = cayley(rand_skew(rng))
        d0 = d_of(O)
        for Pm in [sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
                   sp.Matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])]:
            sgns.add(sp.simplify(d_of(Pm * O) / d0))
            sgns.add(sp.simplify(d_of(O * Pm) / d0))
        sgns.add(sp.simplify(d_of(O.T) / d0))
    check(sgns <= {1, -1}, f"EXACT d(PO), d(OP), d(O^T) = +-d(O) on 6 rational O (ratios seen {sorted(sgns)})")


if "c9" in SECTIONS:
    run_c9()


# ---------------------------------------------------------------------------
# (C.10): q(c_W) = 0 and C(c_W) = det(w_ij - w_i0 - w_0j + w_00)
# ---------------------------------------------------------------------------
def cW_vec(W):
    return [sp.expand(sum(W[i, j] * NVEC[((i,), (j,))][k] for i in range(4) for j in range(4))) for k in range(14)]


def run_c10():
    log("== (C.10) in sixteen independent weights")
    W = sp.Matrix(4, 4, sp.symbols("w0:16"))
    c = cW_vec(W)
    check(qval(c) == 0, "EXACT q(c_W) = 0 (16 symbolic weights)")
    Dm = sp.Matrix(3, 3, lambda i, j: W[i + 1, j + 1] - W[i + 1, 0] - W[0, j + 1] + W[0, 0])
    check(sp.expand(Cval(c) - Dm.det()) == 0, "EXACT C(c_W) = det(w_ij - w_i0 - w_0j + w_00)_{i,j=1..3}")
    # and <rho(U), c_W> = sum w_ij |U_ij|^2 on a rational orthogonal sample
    O = cayley(rand_skew(random.Random(10)))
    rho, _ = rho_exact(O)
    check(sp.expand(sum(rho[k] * c[k] for k in range(14)) - sum(W[i, j] * O[i, j]**2 for i in range(4) for j in range(4))) == 0,
          "EXACT <rho(O), c_W> = sum w_ij O_ij^2 (rational O)")


if "c10" in SECTIONS:
    run_c10()


# ---------------------------------------------------------------------------
# Lemma C.8: Cauchy form <-> d = 0 (exact Cauchy samples; numeric divisor samples)
# ---------------------------------------------------------------------------
def cauchy_orth_data(a, c):
    """For real nodes a (rows) and c (columns) return exact u_i^2, v_j^2 (rationals) with
    O_ij = u_i v_j / (a_i - c_j) real orthogonal, or None if signs are inconsistent."""
    u2 = []
    for i in range(4):
        num = sp.prod([a[i] - c[j] for j in range(4)])
        den = sp.prod([a[k] - a[i] for k in range(4) if k != i])
        u2.append(-num / den)
    if all(x > 0 for x in u2):
        pass
    elif all(x < 0 for x in u2):
        u2 = [-x for x in u2]
    else:
        return None
    v2 = [1 / sum(u2[i] / (a[i] - c[j]) ** 2 for i in range(4)) for j in range(4)]
    return u2, v2


def cauchy_sq_minors(a, c, u2, v2):
    """exact squared minors of the Cauchy orthogonal matrix (rational)."""
    out = {}
    for (I, J) in MINOR_KEYS:
        if len(I) == 0:
            out[(I, J)] = sp.Integer(1)
            continue
        Cm = sp.Matrix(len(I), len(J), lambda r, s: 1 / (a[I[r]] - c[J[s]]))
        out[(I, J)] = sp.prod([u2[i] for i in I]) * sp.prod([v2[j] for j in J]) * Cm.det() ** 2
    return out


def rand_interlaced(rng, H=30):
    vals = sorted(rng.sample(range(-H, H), 8))
    vals = [sp.Rational(v, rng.choice([1, 2, 3])) if False else sp.Integer(v) for v in vals]
    a = vals[0::2]
    c = vals[1::2]
    # random relabelling of rows / columns (labels need not be sorted)
    rng.shuffle(a)
    rng.shuffle(c)
    return a, c


def num_d(O):
    pairs = [(1, 2), (1, 3), (2, 3)]
    return np.linalg.det(np.array([[O[r, k] * O[s, k] for k in range(3)] for (r, s) in pairs]))


SKEW_BASIS = []
for (i, j) in itertools.combinations(range(4), 2):
    E = np.zeros((4, 4))
    E[i, j], E[j, i] = 1.0, -1.0
    SKEW_BASIS.append(E)


def rand_SO4(nprng):
    Q, R = np.linalg.qr(nprng.standard_normal((4, 4)))
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def divisor_samples(nprng, n, max_tries=None):
    """real O in SO(4) with d(O) = 0 (to ~1e-15) found by Brent root-finding on random geodesics."""
    from scipy.optimize import brentq
    out = []
    tries = 0
    while len(out) < n:
        tries += 1
        O1 = rand_SO4(nprng)
        Om = sum(w * E for w, E in zip(nprng.standard_normal(6), SKEW_BASIS))
        Om /= np.linalg.norm(Om)
        f = lambda t: num_d(O1 @ expm(t * Om))
        ts = np.linspace(0, 3.2, 65)
        vals = [f(t) for t in ts]
        for k in range(len(ts) - 1):
            if vals[k] * vals[k + 1] < 0:
                t0 = brentq(f, ts[k], ts[k + 1], xtol=1e-15)
                out.append(O1 @ expm(t0 * Om))
                break
    return out


def nodes_from_reciprocal(O):
    R = 1.0 / O
    Usvd, s, Vt = np.linalg.svd(R)
    X = Usvd[:, :2] * s[:2]
    Y = Vt[:2, :].T
    a = X[:, 0] / X[:, 1]
    c = -Y[:, 1] / Y[:, 0]
    return a, c, s


def alternates(a, c):
    pts = sorted([(x, 'a') for x in a] + [(x, 'c') for x in c])
    labs = [p[1] for p in pts]
    return all(labs[k] != labs[(k + 1) % 8] for k in range(8))


def run_c8():
    log("== Lemma C.8 (Cauchy form <-> d = 0)")
    rng = random.Random(8)
    n_ok, n_tot, n_noninter_ok = 0, 0, 0
    for t in range(40):
        a, c = rand_interlaced(rng)
        data = cauchy_orth_data(a, c)
        n_tot += 1
        if data is None:
            continue
        u2, v2 = data
        Cm = lambda r, k: 1 / (a[r] - c[k])
        pairs = [(1, 2), (1, 3), (2, 3)]
        dred = sp.Matrix([[Cm(r, k) * Cm(s, k) for k in range(3)] for (r, s) in pairs]).det()
        # orthogonality of O = D_u C D_v via squared quantities: column norms and orthogonality
        orth = all(sp.simplify(sum(u2[i] / ((a[i] - c[j]) * (a[i] - c[k])) for i in range(4))
                               * (1 if j != k else v2[j]) - (1 if j == k else 0)) == 0
                   for j in range(4) for k in range(4))
        if dred == 0 and orth:
            n_ok += 1
    check(n_ok == n_tot, f"EXACT {n_ok}/{n_tot} interlaced rational node sets give a real orthogonal "
                         f"Cauchy matrix with d(O) = 0 (d computed as (u1u2u3)^2 (v0v1v2)^2 det[C_rk C_sk])")
    # non-interlaced nodes never give a real orthogonal Cauchy matrix (sign obstruction)
    bad = 0
    for t in range(40):
        vals = sorted(rng.sample(range(-30, 30), 8))
        a = [sp.Integer(v) for v in vals[:4]]
        c = [sp.Integer(v) for v in vals[4:]]
        if cauchy_orth_data(a, c) is not None:
            bad += 1
    check(bad == 0, "EXACT 40 non-alternating node sets: no real orthogonal Cauchy matrix (consistent with "
                    "the alternation claim in the proof of Lemma C.10)")
    # one fully explicit exact instance with radicals
    a = [sp.Integer(x) for x in (0, 2, 4, 6)]
    c = [sp.Integer(x) for x in (1, 3, 5, 7)]
    u2, v2 = cauchy_orth_data(a, c)
    O = sp.Matrix(4, 4, lambda i, j: sp.sqrt(u2[i]) * sp.sqrt(v2[j]) / (a[i] - c[j]))
    check((O.T * O - sp.eye(4)).applyfunc(sp.nsimplify).applyfunc(sp.simplify).is_zero_matrix
          and sp.simplify(d_of(O)) == 0,
          "EXACT nodes a=(0,2,4,6), c=(1,3,5,7): O (with radicals) is orthogonal and d(O)=0")
    # forward direction, numerically: d(O)=0 samples have rank-2 reciprocal and alternating nodes
    nprng = np.random.default_rng(80)
    samples = divisor_samples(nprng, 200)
    worst_rank = max(nodes_from_reciprocal(O)[2][2] / nodes_from_reciprocal(O)[2][0] for O in samples)
    alt = sum(alternates(*nodes_from_reciprocal(O)[:2]) for O in samples)
    check(worst_rank < 1e-8, f"NUMERIC 200 real O on d=0 (Brent, seed 80): sigma_3/sigma_1 of (1/O_ij) "
                             f"<= {worst_rank:.1e} (rank two, as Lemma C.8 asserts)")
    check(alt == 200, f"NUMERIC {alt}/200 of those have alternating row/column nodes on RP^1")
    return samples


DIV_SAMPLES = None
if "c8" in SECTIONS:
    DIV_SAMPLES = run_c8()


# ---------------------------------------------------------------------------
# Lemma C.9: (C.11) witness on random rational O and the boundary example
# ---------------------------------------------------------------------------
def entry_critical(W, O):
    """derivatives of sum w_ij O_ij^2 along left skew generators (exact)."""
    out = []
    for (i, j) in itertools.combinations(range(4), 2):
        E = sp.zeros(4)
        E[i, j], E[j, i] = 1, -1
        EO = E * O
        out.append(sp.simplify(sum(2 * W[a, b] * O[a, b] * EO[a, b] for a in range(4) for b in range(4))))
    return out


def run_lemC9():
    log("== Lemma C.9: explicit c in K_O with C(c) = d(O) prod xi eta / prod O_ij")
    rng = random.Random(19)
    n_ok = 0
    N = 12
    done, skipped = 0, 0
    while done < N:
        O = cayley(rand_skew(rng, H=9))
        allnz = all(minor(O, I, J) != 0 for (I, J) in MINOR_KEYS)
        x1, x2 = sp.Rational(rng.randint(1, 9), rng.randint(1, 5)), sp.Rational(-rng.randint(1, 9), rng.randint(1, 5))
        x3 = -(O[1, 0] * x1 + O[2, 0] * x2) / O[3, 0]
        xi = sp.Matrix([0, x1, x2, x3])
        eta = O.T * xi
        if not (allnz and all(v != 0 for v in list(xi[1:]) + list(eta[1:]))):
            skipped += 1   # sample lies on a Pluecker divisor (outside the lemma's first case)
            continue
        done += 1
        W = sp.Matrix(4, 4, lambda i, j: 0 if i == 0 or j == 0 else xi[i] * eta[j] / O[i, j])
        conds = [allnz, eta[0] == 0, all(v != 0 for v in list(xi[1:]) + list(eta[1:]))]
        conds.append((W.multiply_elementwise(O) * O.T - xi * xi.T).is_zero_matrix)
        conds.append(all(v == 0 for v in entry_critical(W, O)))
        h = (xi.T * xi)[0]
        conds.append(sp.simplify(sum(W[i, j] * O[i, j] ** 2 for i in range(4) for j in range(4)) - h) == 0)
        cvec = [ci - (h if k == 13 else 0) for k, ci in enumerate(cW_vec(W))]
        Cc = sp.simplify(Cval(cvec))
        pred = d_of(O) * sp.prod(list(xi[1:])) * sp.prod(list(eta[1:])) / sp.prod([O[i, j] for i in range(1, 4) for j in range(1, 4)])
        conds.append(sp.simplify(Cc - pred) == 0 and Cc != 0)
        n_ok += all(conds)
    log(f"     ({skipped} rational samples skipped because a minor or xi_i/eta_j vanished)")
    check(n_ok == N, f"EXACT (C.11) on {n_ok}/{N} random rational O in SO(4) (all 70 minors nonzero): "
                     f"W o O O^t = xi xi^t, criticality, height, C(c_W - h kappa) = d prod(xi eta)/prod O != 0")
    # the paper's boundary example
    O = sp.Matrix([[0, sp.Rational(2, 3), sp.Rational(2, 3), -sp.Rational(1, 3)],
                   [sp.Rational(2, 3), sp.Rational(5, 9), -sp.Rational(4, 9), sp.Rational(2, 9)],
                   [sp.Rational(2, 3), -sp.Rational(4, 9), sp.Rational(5, 9), sp.Rational(2, 9)],
                   [sp.Rational(1, 3), -sp.Rational(2, 9), -sp.Rational(2, 9), -sp.Rational(8, 9)]])
    W = sp.zeros(4)
    W[1:, 1:] = sp.Matrix([[sp.Rational(9, 5), -sp.Rational(9, 4), 18], [-sp.Rational(9, 4), sp.Rational(9, 5), 18], [18, 18, 18]])
    xi = sp.Matrix([0, 1, 1, -4])
    eta = sp.Matrix([0, 1, 1, 4])
    cvec = [ci - (18 if k == 13 else 0) for k, ci in enumerate(cW_vec(W))]
    conds = {
        "O orthogonal, det 1": (O * O.T - sp.eye(4)).is_zero_matrix and O.det() == 1,
        "O^t xi = eta": O.T * xi == eta,
        "(W o O) O^t = xi xi^t": (W.multiply_elementwise(O) * O.T - xi * xi.T).is_zero_matrix,
        "critical": all(v == 0 for v in entry_critical(W, O)),
        "h = 18": sum(W[i, j] * O[i, j] ** 2 for i in range(4) for j in range(4)) == 18,
        "d(O) = -16/729": d_of(O) == sp.Rational(-16, 729),
        "C(c_W - 18 kappa) = -531441/200": Cval(cvec) == sp.Rational(-531441, 200),
        "only zero entry is O_00 (indecomposable, stable)": sum(1 for v in O if v == 0) == 1,
    }
    for k, v in conds.items():
        check(v, f"EXACT Lemma C.9 example: {k}")
    zero_minors = [(I, J) for (I, J) in MINOR_KEYS if minor(O, I, J) == 0]
    log(f"     example: vanishing minors {zero_minors} (entry O_00 and its complementary 3x3 minor)")


if "lemC9" in SECTIONS:
    run_lemC9()


# ---------------------------------------------------------------------------
# (C.8) tested directly: C restricted to K_O versus det of the phase Hessian
# ---------------------------------------------------------------------------
def lam_matrix():
    """L (70 x 14): <rho(U), c> = sum_{IJ} (L c)_{IJ} |det U[I,J]|^2 (index = MINOR_KEYS order)."""
    A = MMON.extract(PIV_ROWS, list(range(14)))
    Ainv = A.inv()
    L = sp.zeros(70, 14)
    for idx, (I, J) in enumerate(MINOR_KEYS):
        Ic = [i for i in range(4) if i not in I]
        S = tuple(sorted(Ic + [4 + j for j in J]))
        r = MON_INDEX[S]
        if r in PIV_ROWS:
            col = PIV_ROWS.index(r)
            for k in range(14):
                L[idx, k] = (-1) ** len(J) * Ainv[k, col]
    return L


def second_order_sq_minors(O, Xi):
    """t^2-coefficient of |det U(t)[I,J]|^2, U(t) = (1 + i t Xi - t^2 Xi^2/2) O."""
    t = sp.Symbol('t', real=True)
    U = (sp.eye(4) + sp.I * t * Xi - t ** 2 * Xi * Xi / 2) * O
    Uc = (sp.eye(4) - sp.I * t * Xi - t ** 2 * Xi * Xi / 2) * O
    out = []
    for (I, J) in MINOR_KEYS:
        if len(I) == 0:
            out.append(0)
            continue
        D1 = sp.expand(U.extract(list(I), list(J)).det())
        D2 = sp.expand(Uc.extract(list(I), list(J)).det())
        out.append(sp.expand(D1 * D2).coeff(t, 2))
    return out


def first_order_sq_minors(O, E):
    t = sp.Symbol('t')
    Ot = (sp.eye(4) + t * E) * O
    out = []
    for (I, J) in MINOR_KEYS:
        if len(I) == 0:
            out.append(0)
            continue
        D = sp.expand(Ot.extract(list(I), list(J)).det())
        out.append(sp.expand(D * D).coeff(t, 1))
    return out


def alpha_test(O, label):
    L = lam_matrix()
    sq0 = [minor(O, I, J) ** 2 for (I, J) in MINOR_KEYS]
    rows = [[sum(sq0[r] * L[r, k] for r in range(70)) for k in range(14)]]
    for (i, j) in itertools.combinations(range(4), 2):
        E = sp.zeros(4)
        E[i, j], E[j, i] = 1, -1
        s1 = first_order_sq_minors(O, E)
        rows.append([sum(s1[r] * L[r, k] for r in range(70)) for k in range(14)])
    Mk = sp.Matrix(rows)
    basis = Mk.nullspace()
    dimK = len(basis)
    # gauge space {D1 + O D2 O^t} and a phase complement
    gauge = []
    for i in range(4):
        D = sp.zeros(4); D[i, i] = 1; gauge.append(D)
        gauge.append(O * D * O.T)
    gvecs = sp.Matrix([[G[i, j] for i in range(4) for j in range(i, 4)] for G in gauge])
    grank = gvecs.rank()
    cands = []
    for (i, j) in itertools.combinations(range(4), 2):
        X = sp.zeros(4); X[i, j] = X[j, i] = 1; cands.append(X)
    phase = []
    cur = gvecs
    for X in cands:
        trial = cur.col_join(sp.Matrix([[X[i, j] for i in range(4) for j in range(i, 4)]]))
        if trial.rank() > cur.rank():
            phase.append(X); cur = trial
        if len(phase) == 3:
            break
    # quadratic form values via polarisation
    def Q(X):
        s2 = second_order_sq_minors(O, X)
        return [sum(s2[r] * L[r, k] for r in range(70)) for k in range(14)]   # linear in c
    Qd = [Q(X) for X in phase]
    Qp = {(a, b): Q(phase[a] + phase[b]) for (a, b) in [(0, 1), (0, 2), (1, 2)]}
    ts = sp.symbols("t1:%d" % (dimK + 1))
    cvec = [sum(ts[q] * basis[q][k] for q in range(dimK)) for k in range(14)]
    lin = lambda coeffs: sp.expand(sum(coeffs[k] * cvec[k] for k in range(14)))
    S = sp.zeros(3)
    for a in range(3):
        S[a, a] = lin(Qd[a])
    for (a, b), v in Qp.items():
        S[a, b] = S[b, a] = sp.expand((lin(v) - lin(Qd[a]) - lin(Qd[b])) / 2)
    detS = sp.expand(S.det())
    Cc = Cval(cvec)
    # gauge null-check: bilinear term of a phase direction with a gauge direction vanishes on K_O
    Gx = O * sp.diag(1, 2, -1, 3) * O.T
    gnull = all(sp.expand(lin([x - y - z for x, y, z in zip(Q(phase[0] + Gx), Qd[0], Q(Gx))])) == 0 for _ in [0])
    if detS == 0:
        ratio = None
    else:
        ratio = sp.cancel(Cc / detS)
    log(f"     [{label}] dim K_O = {dimK}, gauge rank = {grank}, gauge directions null: {gnull}, "
        f"C|K_O / det S_O = {ratio if ratio is None or ratio.is_number else 'NOT CONSTANT'}")
    return dimK, grank, gnull, ratio, d_of(O)


def run_c8prop():
    log("== (C.8) C|K_O = alpha(O) det S_O, tested on exact orthogonal matrices (part-1 input to Lemma C.9)")
    rng = random.Random(88)
    tests = []
    while len(tests) < 2:
        O = cayley(rand_skew(rng, H=4))
        if all(minor(O, I, J) != 0 for (I, J) in MINOR_KEYS):
            tests.append(("random rational O", O))
    tests.append(("Lemma C.9 boundary example", sp.Matrix([[0, sp.Rational(2, 3), sp.Rational(2, 3), -sp.Rational(1, 3)],
                   [sp.Rational(2, 3), sp.Rational(5, 9), -sp.Rational(4, 9), sp.Rational(2, 9)],
                   [sp.Rational(2, 3), -sp.Rational(4, 9), sp.Rational(5, 9), sp.Rational(2, 9)],
                   [sp.Rational(1, 3), -sp.Rational(2, 9), -sp.Rational(2, 9), -sp.Rational(8, 9)]])))
    # an indecomposable d = 0 point with zero entries: row (3/5, 4/5, 0, 0)
    a, b = sp.Rational(3, 5), sp.Rational(4, 5)
    R3 = cayley(sp.Matrix([[0, sp.Rational(1, 2), sp.Rational(1, 3)], [-sp.Rational(1, 2), 0, 1], [-sp.Rational(1, 3), -1, 0]]))
    x, y, z = R3[:, 0], R3[:, 1], R3[:, 2]
    Ob = sp.Matrix([[a, b, 0, 0]] + [[-b * x[i], a * x[i], y[i], z[i]] for i in range(3)])
    tests.append(("indecomposable d=0 point, row (3/5,4/5,0,0)", Ob))
    results = []
    for label, O in tests:
        assert (O * O.T - sp.eye(4)).is_zero_matrix
        results.append((label,) + alpha_test(O, label))
    ok = all(r[1] == 7 and r[2] == 7 and r[3] and r[4] is not None and r[4].is_number and r[4] != 0
             for r in results[:3])
    check(ok, "EXACT on 2 random rational O and the Lemma C.9 example: dim K_O = 7, gauge space has dim 7, "
              "C|K_O is a NONZERO constant multiple of det S_O (so (C.8) holds there and alpha(O) != 0)")
    r = results[3]
    log(f"     info: at the d=0 indecomposable point d={r[5]}, alpha = {r[4]}")


if "c8prop" in SECTIONS:
    run_c8prop()


# ---------------------------------------------------------------------------
# Lemma C.10: common permutation weights on the real Cauchy divisor
# ---------------------------------------------------------------------------
EQ_KEYS = [(I, J) for (I, J) in MINOR_KEYS if len(I) in (1, 2)]   # sizes 0,3,4 follow (Jacobi)
INCID = np.array([[1.0 if tuple(sorted(s[i] for i in I)) == J else 0.0 for s in PERMS] for (I, J) in EQ_KEYS])


def lp_weights(sq, allowed=None):
    """min L1 residual of A p = sq, p >= 0, sum p = 1 (floats).  Returns (residual, p)."""
    idx = list(range(24)) if allowed is None else list(allowed)
    A = np.vstack([INCID[:, idx], np.ones((1, len(idx)))])
    rhs = np.concatenate([np.array([float(sq[k]) for k in EQ_KEYS]), [1.0]])
    nE, nV = A.shape
    cost = np.concatenate([np.zeros(nV), np.ones(2 * nE)])
    Aeq = np.hstack([A, np.eye(nE), -np.eye(nE)])
    res = linprog(cost, A_eq=Aeq, b_eq=rhs, bounds=[(0, None)] * (nV + 2 * nE), method="highs")
    p = np.zeros(24)
    p[idx] = res.x[:nV]
    return res.fun, p


def exact_weights(sq, p_float, allowed=None):
    """exact rational solution on the numerical support; verifies p >= 0 and all identities."""
    supp = [k for k in range(24) if p_float[k] > 1e-10 and (allowed is None or k in allowed)]
    A = sp.Matrix([[int(INCID[r, k]) for k in supp] for r in range(len(EQ_KEYS))] + [[1] * len(supp)])
    rhs = sp.Matrix([sq[k] for k in EQ_KEYS] + [1])
    try:
        sol, params = A.gauss_jordan_solve(rhs)
    except ValueError:
        return False, None
    if params.shape[0] > 0:
        sol = sol.subs({pp: 0 for pp in params})
    p = [sp.Integer(0)] * 24
    for t, k in enumerate(supp):
        p[k] = sp.simplify(sol[t])
    ok = all(x >= 0 for x in p) and sum(p) == 1
    for (I, J) in MINOR_KEYS:
        tot = sum(p[k] for k, s in enumerate(PERMS) if tuple(sorted(s[i] for i in I)) == J)
        ok &= sp.simplify(tot - sq[(I, J)]) == 0
    return ok, p


def noncrossing_matchings(pts):
    if not pts:
        return [[]]
    out = []
    a = pts[0]
    for k in range(1, len(pts), 2):
        b = pts[k]
        for L in noncrossing_matchings(pts[1:k]):
            for R in noncrossing_matchings(pts[k + 1:]):
                out.append([(a, b)] + L + R)
    return out


def lam_support(a, c):
    """permutations sigma (row i -> column sigma(i)) coming from the 14 noncrossing matchings of
    the 8 nodes in their cyclic order (rows at a_i, columns at c_j)."""
    pts = sorted([(float(x), ('r', i)) for i, x in enumerate(a)] + [(float(x), ('c', j)) for j, x in enumerate(c)])
    labels = [p[1] for p in pts]
    perms = set()
    for M in noncrossing_matchings(list(range(8))):
        s = [None] * 4
        good = True
        for (u, v) in M:
            lu, lv = labels[u], labels[v]
            if lu[0] == lv[0]:
                good = False
                break
            r, cc = (lu, lv) if lu[0] == 'r' else (lv, lu)
            s[r[1]] = cc[1]
        if good:
            perms.add(PERMS.index(tuple(s)))
    return sorted(perms)


def sq_from_numeric(O):
    return {(I, J): (np.linalg.det(O[np.ix_(I, J)]) ** 2 if len(I) else 1.0) for (I, J) in MINOR_KEYS}


def run_lemC10():
    global DIV_SAMPLES
    log("== Lemma C.10: common weights for ALL squared minors on d(O) = 0")
    rng = random.Random(10)
    # (1) exact Cauchy samples, incl. near-coincident nodes
    n_ok = n_tot = 0
    worst24 = 0.0
    minp = None
    for t in range(60):
        if t < 40:
            a, c = rand_interlaced(rng)
        else:   # adversarial: a node pair 1/1000 apart, others spread
            base = sorted(rng.sample(range(-50, 50), 7))
            k = rng.randrange(7)
            vals = sorted([sp.Integer(x) for x in base] + [base[k] + sp.Rational(1, 1000)])
            a, c = vals[0::2], vals[1::2]
        u2, v2 = cauchy_orth_data(a, c)
        sq = cauchy_sq_minors(a, c, u2, v2)
        n_tot += 1
        worst24 = max(worst24, lp_weights(sq)[0])
        lam = lam_support(a, c)
        A = sp.Matrix([[int(INCID[r, k]) for k in lam] for r in range(len(EQ_KEYS))] + [[1] * len(lam)])
        rhs = sp.Matrix([sq[k] for k in EQ_KEYS] + [1])
        ok = len(lam) == 14 and A.rank() == 14
        if ok:   # unique solution from 14 independent rows, then exact residual on all rows
            _, piv = A.T.rref()
            sol = A.extract(list(piv), list(range(14))).LUsolve(rhs.extract(list(piv), [0]))
            ok = (A * sol - rhs).is_zero_matrix
        if ok:
            p = [sp.Integer(0)] * 24
            for q, k in enumerate(lam):
                p[k] = sol[q]
            ok = all(x > 0 for x in sol)
            for (I, J) in MINOR_KEYS:
                tot = sum(p[k] for k, s_ in enumerate(PERMS) if tuple(sorted(s_[i] for i in I)) == J)
                ok = ok and sp.simplify(tot - sq[(I, J)]) == 0
            mp = min(sol)
            minp = mp if minp is None else min(minp, mp)
        n_ok += ok
    check(n_ok == n_tot, f"EXACT {n_ok}/{n_tot} exact Cauchy orthogonal matrices on d=0 (20 with two nodes 1/1000 apart): "
                         f"on the 14 permutations of the noncrossing matchings in alternating order (the Lam support) "
                         f"the system is nonsingular and its unique rational solution is STRICTLY positive and "
                         f"reproduces all 70 squared minors (smallest weight {float(minp):.2e})")
    log(f"     same samples, float LP over all 24 permutations: worst L1 residual {worst24:.1e}")
    # (2) numeric divisor samples (independent of the Cauchy construction)
    nprng = np.random.default_rng(81)
    if DIV_SAMPLES is None:
        DIV_SAMPLES = divisor_samples(np.random.default_rng(80), 200)
    worst = worst_lam = 0.0
    for O in DIV_SAMPLES + divisor_samples(nprng, 300):
        sq = sq_from_numeric(O)
        worst = max(worst, lp_weights(sq)[0])
        a, c, _ = nodes_from_reciprocal(O)
        worst_lam = max(worst_lam, lp_weights(sq, lam_support(a, c))[0])
    check(worst < 1e-9 and worst_lam < 1e-9,
          f"NUMERIC 500 real O with d(O)=0 found by root-finding (seeds 80,81): min L1 residual of the weight LP "
          f"<= {worst:.1e} (all 24 perms), <= {worst_lam:.1e} (14 Lam perms)")
    # (3) zero-entry boundary: row (a,b,0,0), [x y z] in SO(3), exact
    n_ok = n_tot = 0
    for t in range(25):
        s = sp.Rational(rng.randint(1, 9), rng.randint(1, 9))
        a, b = (1 - s * s) / (1 + s * s), 2 * s / (1 + s * s)
        R3 = cayley(rand_skew(rng, 3, 6))
        if t % 5 == 4:   # force a zero coordinate x_1 = 0 (Householder onto a rational unit vector (0,s1,s2))
            q = sp.Rational(rng.randint(1, 7), rng.randint(1, 7))
            target = sp.Matrix([0, (1 - q * q) / (1 + q * q), 2 * q / (1 + q * q)])
            w = R3[:, 0] - target
            if w != sp.zeros(3, 1):
                R3 = (sp.eye(3) - 2 * w * w.T / (w.T * w)[0]) * R3
        if R3.det() < 0:
            R3[:, 2] = -R3[:, 2]
        x, y, z = R3[:, 0], R3[:, 1], R3[:, 2]
        O = sp.Matrix([[a, b, 0, 0]] + [[-b * x[i], a * x[i], y[i], z[i]] for i in range(3)])
        assert (O * O.T - sp.eye(4)).is_zero_matrix
        if d_of(O) != 0:
            continue
        sq = {(I, J): minor(O, I, J) ** 2 for (I, J) in MINOR_KEYS}
        n_tot += 1
        res, p = lp_weights(sq)
        n_ok += exact_weights(sq, p)[0]
    check(n_ok == n_tot and n_tot > 0, f"EXACT {n_ok}/{n_tot} rational O with a row (a,b,0,0) (d=0, zero entries): exact common weights")
    # (4) decomposable blocks (Birkhoff products) exact
    n_ok = n_tot = 0
    for t in range(10):
        R3 = cayley(rand_skew(rng, 3, 6))
        s = sp.Rational(rng.randint(1, 9), rng.randint(1, 9))
        R2 = sp.Matrix([[(1 - s * s) / (1 + s * s), -2 * s / (1 + s * s)], [2 * s / (1 + s * s), (1 - s * s) / (1 + s * s)]])
        s2 = sp.Rational(rng.randint(1, 9), rng.randint(1, 9))
        R2b = sp.Matrix([[(1 - s2 * s2) / (1 + s2 * s2), -2 * s2 / (1 + s2 * s2)], [2 * s2 / (1 + s2 * s2), (1 - s2 * s2) / (1 + s2 * s2)]])
        for O in (sp.diag(R3, 1), sp.diag(R2, R2b)):
            Pm = sp.Matrix(PERMS[rng.randrange(24)]).T
            Pr = sp.zeros(4); Pc = sp.zeros(4)
            pr, pc = PERMS[rng.randrange(24)], PERMS[rng.randrange(24)]
            for i in range(4):
                Pr[i, pr[i]] = 1; Pc[i, pc[i]] = 1
            O = Pr * O * Pc
            sq = {(I, J): minor(O, I, J) ** 2 for (I, J) in MINOR_KEYS}
            n_tot += 1
            res, p = lp_weights(sq)
            n_ok += exact_weights(sq, p)[0]
    check(n_ok == n_tot, f"EXACT {n_ok}/{n_tot} decomposable rational O (3+1 and 2+2 blocks, permuted): exact common weights")
    # (5) sanity: generic O (d != 0) usually has NO common weights (so the LP test has teeth)
    infeas = 0
    for t in range(300):
        O = rand_SO4(nprng)
        if lp_weights(sq_from_numeric(O))[0] > 1e-6:
            infeas += 1
    log(f"     sanity NUMERIC: {infeas}/300 Haar-random O in SO(4) admit no common weights (seed 81)")
    # (6) adversarial search: maximise LP residual along the divisor (random walks on d = 0)
    from scipy.optimize import brentq
    worst_walk = 0.0
    O = DIV_SAMPLES[0]
    steps = 0
    for t in range(400):
        Om = sum(w * E for w, E in zip(nprng.standard_normal(6), SKEW_BASIS))
        Om /= np.linalg.norm(Om)
        Om2 = sum(w * E for w, E in zip(nprng.standard_normal(6), SKEW_BASIS))
        O1 = O @ expm(0.3 * Om)
        f = lambda tt: num_d(O1 @ expm(tt * Om2))
        ts = np.linspace(-1.5, 1.5, 31)
        vals = [f(tt) for tt in ts]
        for k in range(30):
            if vals[k] * vals[k + 1] < 0:
                O = O1 @ expm(brentq(f, ts[k], ts[k + 1], xtol=1e-15) * Om2)
                steps += 1
                worst_walk = max(worst_walk, lp_weights(sq_from_numeric(O))[0])
                break
    check(worst_walk < 1e-9, f"NUMERIC random walk on d=0 ({steps} points): worst weight-LP residual {worst_walk:.1e}")


if "lemC10" in SECTIONS:
    run_lemC10()


# ---------------------------------------------------------------------------
# (C.13), (C.14)
# ---------------------------------------------------------------------------
def rat_unit3(rng):
    u, v = sp.Rational(rng.randint(-7, 7), rng.randint(1, 5)), sp.Rational(rng.randint(-7, 7), rng.randint(1, 5))
    n = 1 + u * u + v * v
    return sp.Matrix([2 * u / n, 2 * v / n, (u * u + v * v - 1) / n])


def run_c13():
    log("== (C.13) and (C.14)")
    rng = random.Random(13)
    signs_plus, signs_minus = set(), set()
    for t in range(12):
        R = cayley(rand_skew(rng))
        r = R[:, 0]
        s3 = rat_unit3(rng)
        s = sp.Matrix([0, s3[0], s3[1], s3[2]])
        w = r - s
        H = sp.eye(4) - 2 * w * w.T / (w.T * w)[0]
        O = H * R                       # det -1, first column s, O_00 = 0
        for Om, bucket in ((O * sp.diag(1, 1, 1, -1), signs_plus), (O, signs_minus)):
            assert Om[0, 0] == 0 and (Om * Om.T - sp.eye(4)).is_zero_matrix
            rhs = Om[0, 1] * Om[0, 2] * Om[0, 3] * Om[1, 0] * Om[2, 0] * Om[3, 0]
            dv = d_of(Om)
            bucket.add(sp.simplify(dv / rhs) if rhs != 0 else ('rhs0', dv))
    check(signs_plus == {1}, f"EXACT (C.13) d(O) = O01 O02 O03 O10 O20 O30 when O00 = 0, det O = +1 (12 rational samples): ratios {signs_plus}")
    log(f"     on det O = -1 the ratio is {signs_minus} (the paper normalises det O = 1 first, so this is consistent)")
    # (C.14)
    sgn_seen = set()
    ok = True
    for t in range(15):
        s = sp.Rational(rng.randint(1, 9), rng.randint(1, 9))
        a, b = (1 - s * s) / (1 + s * s), 2 * s / (1 + s * s)
        R3 = cayley(rand_skew(rng, 3, 6))
        x, y, z = R3[:, 0], R3[:, 1], R3[:, 2]
        O = sp.Matrix([[a, b, 0, 0]] + [[-b * x[i], a * x[i], y[i], z[i]] for i in range(3)])
        assert d_of(O) == 0
        tt = sp.Symbol('tt')
        per = []
        for i in (1, 2, 3):
            E = sp.zeros(4)
            E[0, i], E[i, 0] = 1, -1
            dd = sp.expand(d_of((sp.eye(4) + tt * E) * O)).coeff(tt, 1)
            j, k = [q for q in (1, 2, 3) if q != i]
            pred = a * b * x[j - 1] * x[k - 1] * y[i - 1] * z[i - 1]
            per.append(sp.simplify(dd / pred) if pred != 0 else ('pred0', dd))
        if len(set(per)) != 1:
            ok = False
        sgn_seen |= set(per)
    check(ok and sgn_seen <= {1, -1}, f"EXACT (C.14) dd(E_i O) = s * ab x_j x_k y_i z_i with one common sign s per sample "
                                      f"(15 rational samples, signs seen {sgn_seen})")


if "c13" in SECTIONS:
    run_c13()


# ---------------------------------------------------------------------------
# Lemma C.11
# ---------------------------------------------------------------------------
def run_c11():
    log("== Lemma C.11")
    Hk = HESS.subs(dict(zip(MS, list(KAPPA))))
    check(Hk.rank() == 4, "EXACT rank H(kappa) = 4 (3.10)")
    a = [0, 1, 3, 7]
    b = [2 + sp.I, -1 + 3 * sp.I, 5 - 2 * sp.I, 4 + 4 * sp.I]
    mv = [sp.expand(x) for x in matching_vec(a + [-x for x in b])]
    y = [sp.re(x) for x in mv]
    yp = [-36, 108, 189, -569, 344, -126, 378, 288, 168, -784, 496, -512, 512, 16]
    check(y == yp, "EXACT Lemma C.11 example: Re m(a,-b) equals the printed y")
    check(qval(y) == 0 and Cval(y) == 0 and any(g != 0 for g in Gval(y)),
          "EXACT example: q(y) = C(y) = 0 and G(y) != 0")
    yk = [y[k] + (1 if k == 13 else 0) for k in range(14)]
    dH = Hval(yk).det()
    check(dH == -114105877589169227887017984, f"EXACT example: det H(y + kappa) = {dH}")
    h = sp.Symbol('h')
    poly = sp.Poly(sp.expand((Hval(y) + h * Hk).det()), h)
    check(poly.degree() == 2 and poly.coeffs() == [dH] and len(poly.terms()) == 1,
          f"EXACT example: det(H(y) + h H(kappa)) = {dH} h^2 exactly (consistent with s(y+h kappa) = beta h)")
    # random rational instances of y = Re(theta m(a,-b))
    rng = random.Random(11)
    ok = True
    for t in range(6):
        a = [sp.Integer(v) for v in rng.sample(range(-9, 9), 4)]          # distinct real a
        b = [sp.Integer(rng.randint(-6, 6)) + sp.I * rng.choice([-5, -3, -2, -1, 1, 2, 3, 4]) for _ in range(4)]
        th = sp.Integer(rng.randint(-3, 3)) + sp.I * rng.randint(1, 3)
        y = [sp.re(sp.expand(th * x)) for x in matching_vec(a + [-x for x in b])]
        pl = sp.Poly(sp.expand((Hval(y) + h * Hk).det()), h)
        ok &= qval(y) == 0 and Cval(y) == 0 and (not pl.is_zero) and pl.degree() == 2 and len(pl.terms()) == 1
    check(ok, "EXACT 6 random rational (distinct real a, nonreal b, theta): q(y)=C(y)=0 and det(H(y)+hH(kappa)) = c*h^2, c != 0")
    # general u on {C = 0, q = 0} (not of the form Re(theta m)): det(H(u)+hH(kappa)) has degree <= 2 and is a square
    ok = True
    degs = []
    for t in range(6):
        u = [sp.Rational(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(14)]
        u[0] = u[0] or 1
        u[6] = u[1] * u[5] / u[0]
        u4 = sp.Symbol('u4')
        uu = list(u)
        uu[4] = u4
        sol = sp.solve(Cval(uu), u4)
        if not sol:
            continue
        u[4] = sol[0]
        pl = sp.Poly(sp.expand((Hval(u) + h * Hk).det()), h)
        degs.append(pl.degree())
        disc_ok = pl.degree() < 2 or sp.discriminant(pl.as_expr(), h) == 0
        ok &= pl.degree() <= 2 and disc_ok and (pl.degree() != 1)
    check(ok, f"EXACT 6 random rational u on {{C=0,q=0}}: det(H(u)+hH(kappa)) has degree {degs} (<=2) and is "
              f"a constant times a square in h, as the proof of Lemma C.11 requires")


if "c11" in SECTIONS:
    run_c11()


# ---------------------------------------------------------------------------
# Numeric test of the final assembly's key consequence: on generic (a real, b, theta)
# EVERY critical point O of f(O) = Re(theta det(D_a + O D_b O^t)) on SO(4) that is
# indecomposable with d(O) != 0 must have f(O) = 0 and det(D_a + O D_b O^t) = 0
# (the proof of Thm C.1 uses only first-order criticality to reach h = 0).
# ---------------------------------------------------------------------------
def adj4(M):
    adj = np.zeros((4, 4), dtype=complex)
    for i in range(4):
        for j in range(4):
            sub = np.delete(np.delete(M, i, 0), j, 1)
            adj[j, i] = (-1) ** (i + j) * np.linalg.det(sub)
    return adj


def f_grad(O, a, b, th, gens):
    B = O @ np.diag(b) @ O.conj().T
    M = np.diag(a).astype(complex) + B
    adj = adj4(M)
    dm = np.linalg.det(M)
    g = np.array([np.real(th * np.trace(adj @ (E @ B - B @ E))) for E in gens])
    return np.real(th * dm), g, dm


def indecomposable(O, tol=1e-6):
    nz = np.abs(O) > tol
    seen_r, seen_c = {0}, set()
    frontier = [('r', 0)]
    while frontier:
        kind, k = frontier.pop()
        if kind == 'r':
            for j in range(4):
                if nz[k, j] and j not in seen_c:
                    seen_c.add(j); frontier.append(('c', j))
        else:
            for i in range(4):
                if nz[i, k] and i not in seen_r:
                    seen_r.add(i); frontier.append(('r', i))
    return len(seen_r) == 4 and len(seen_c) == 4


def run_crit():
    log("== Critical points of Re(theta det) on SO(4) (consequence of the Thm C.1 argument)")
    nprng = np.random.default_rng(2026)
    n_spec, n_start = SZ(30), SZ(40)
    stats = {"decomposable": 0, "d~0": 0, "ambiguous d": 0, "value 0 & det 0": 0, "OTHER": 0, "not converged": 0}
    worst_other = None
    dmax_zero = 0.0
    for sidx in range(n_spec):
        a = nprng.standard_normal(4)
        b = nprng.standard_normal(4) + 1j * nprng.standard_normal(4)
        th = np.exp(1j * nprng.uniform(0, 2 * np.pi))
        scale = np.max(np.abs(sigma_points(a, b)))
        for st in range(n_start):
            O0 = rand_SO4(nprng)
            fun = lambda w: f_grad(expm(sum(wi * E for wi, E in zip(w, SKEW_BASIS))) @ O0, a, b, th, SKEW_BASIS)[1] / scale
            sol = least_squares(fun, np.zeros(6), method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
            O = expm(sum(wi * E for wi, E in zip(sol.x, SKEW_BASIS))) @ O0
            v, g, dm = f_grad(O, a, b, th, SKEW_BASIS)
            if np.max(np.abs(g)) / scale > 1e-10:
                stats["not converged"] += 1
                continue
            dv = abs(num_d(O))
            if not indecomposable(O):
                stats["decomposable"] += 1
            elif dv < 1e-9:
                stats["d~0"] += 1
                dmax_zero = max(dmax_zero, dv)
            elif dv < 1e-4:
                stats["ambiguous d"] += 1
            elif abs(v) / scale < 1e-8 and abs(dm) / scale < 1e-8:
                stats["value 0 & det 0"] += 1
            else:
                stats["OTHER"] += 1
                cand = (abs(v) / scale, abs(dm) / scale, abs(num_d(O)), np.min(np.abs(O)))
                if worst_other is None or cand[0] > worst_other[0]:
                    worst_other = cand
    log(f"     {n_spec} spectra x {n_start} starts (seed 2026); critical points by class: {stats}")
    log(f"     class d~0 means |d(O)| < 1e-9 (largest seen {dmax_zero:.1e}); ambiguous means 1e-9 <= |d| < 1e-4")
    if worst_other:
        log(f"     largest |f|/scale among OTHER: {worst_other[0]:.3e} (|det|/scale {worst_other[1]:.2e}, "
            f"|d| {worst_other[2]:.2e}, min|O_ij| {worst_other[3]:.2e})")
    check(stats["OTHER"] == 0 and stats["ambiguous d"] == 0, "NUMERIC every converged indecomposable critical point with d != 0 has f = 0 and det = 0")


if "crit" in SECTIONS:
    run_crit()


# ---------------------------------------------------------------------------
# Numeric attack on Theorem C.1 itself, and the stronger common-weight (PMSH) form
# ---------------------------------------------------------------------------
from scipy.spatial import ConvexHull

HERM_GENS = []
for (i, j) in itertools.combinations(range(4), 2):
    E = np.zeros((4, 4), dtype=complex); E[i, j], E[j, i] = 1, -1; HERM_GENS.append(E)
    F = np.zeros((4, 4), dtype=complex); F[i, j] = F[j, i] = 1j; HERM_GENS.append(F)
for i in range(4):
    F = np.zeros((4, 4), dtype=complex); F[i, i] = 1j; HERM_GENS.append(F)


def hull_dirs(pts, nprng, n_random=4):
    P2 = np.column_stack([pts.real, pts.imag])
    dirs = []
    try:
        hull = ConvexHull(P2)
        for eq in hull.equations:
            dirs.append(eq[0] - 1j * eq[1])
    except Exception:
        v = P2.max(0) - P2.min(0)
        n = np.array([-v[1], v[0]]) / (np.linalg.norm(v) + 1e-300)
        dirs += [n[0] - 1j * n[1], -n[0] + 1j * n[1]]
    for _ in range(n_random):
        dirs.append(np.exp(1j * nprng.uniform(0, 2 * np.pi)))
    return [d / abs(d) for d in dirs]


def max_support(a, b, th, nprng, unitary=False, starts=3):
    gens = HERM_GENS if unitary else SKEW_BASIS
    pts = sigma_points(a, b)
    best = -np.inf
    k_best = int(np.argmax(np.real(th * pts)))
    for st in range(starts):
        if st == 0:
            O0 = np.zeros((4, 4)); 
            for i in range(4):
                O0[i, PERMS[k_best][i]] = 1.0
            O0 = O0 @ expm(0.05 * sum(w * E for w, E in zip(nprng.standard_normal(6), SKEW_BASIS)))
        else:
            O0 = rand_SO4(nprng)
        if unitary:
            O0 = O0 @ expm(0.3 * sum(w * E for w, E in zip(nprng.standard_normal(16), HERM_GENS)))

        def obj(w):
            U = expm(sum(wi * E for wi, E in zip(w, gens))) @ O0
            return -np.real(th * detfun(a, b, U))
        res = minimize(obj, np.zeros(len(gens)), method="BFGS", options={"gtol": 1e-11, "maxiter": 400})
        # polish with a few Riemannian gradient steps at the found point
        U = expm(sum(wi * E for wi, E in zip(res.x, gens))) @ O0
        best = max(best, -res.fun)
    return best


def run_attack():
    log("== Numeric attack on Theorem C.1 (margin = hull support - max over the group)")
    nprng = np.random.default_rng(4)
    families = {
        "real a, random b": lambda: (nprng.standard_normal(4), nprng.standard_normal(4) + 1j * nprng.standard_normal(4)),
        "clustered a": lambda: (1 + 1e-3 * nprng.standard_normal(4), nprng.standard_normal(4) + 1j * nprng.standard_normal(4)),
        "near-collinear b": lambda: (nprng.standard_normal(4), nprng.standard_normal(4) + 1e-3j * nprng.standard_normal(4)),
        "large modulus ratio": lambda: (nprng.standard_normal(4) * np.array([1e-2, 1, 1e1, 1e2]),
                                        (nprng.standard_normal(4) + 1j * nprng.standard_normal(4)) * np.array([1e2, 1e-1, 1, 3e1])),
        "b on a circle-ish, a real clustered pairs": lambda: (np.array([0, 1e-2, 1, 1 + 1e-2]) + 0.1 * nprng.standard_normal(4),
                                                              np.exp(1j * nprng.uniform(0, 2 * np.pi, 4))),
    }
    worst = (np.inf, None)
    count = 0
    per_fam = {}
    for name, gen in families.items():
        fam_worst = np.inf
        for t in range(SZ(18)):
            a, b = gen()
            pts = sigma_points(a, b)
            scale = np.max(np.abs(pts)) + 1e-300
            for th in hull_dirs(pts, nprng, 2):
                sup = np.max(np.real(th * pts))
                mx = max_support(a, b, th, nprng, unitary=False, starts=2)
                marg = (sup - mx) / scale
                count += 1
                fam_worst = min(fam_worst, marg)
                if marg < worst[0]:
                    worst = (marg, (name, a, b, th))
        per_fam[name] = fam_worst
    for k, v in per_fam.items():
        log(f"     SO(4), {k}: worst relative margin {v:.3e}")
    # circle spectrum a and real a, searched over U(4) directly
    for name, gen in {"a on unit circle (U(4))": lambda: (np.exp(1j * nprng.uniform(0, 2 * np.pi, 4)),
                                                         nprng.standard_normal(4) + 1j * nprng.standard_normal(4)),
                      "a real (U(4))": lambda: (nprng.standard_normal(4), nprng.standard_normal(4) + 1j * nprng.standard_normal(4))}.items():
        fam_worst = np.inf
        for t in range(SZ(8)):
            a, b = gen()
            pts = sigma_points(a, b)
            scale = np.max(np.abs(pts))
            for th in hull_dirs(pts, nprng, 1):
                sup = np.max(np.real(th * pts))
                mx = max_support(a, b, th, nprng, unitary=True, starts=2)
                marg = (sup - mx) / scale
                count += 1
                fam_worst = min(fam_worst, marg)
                if marg < worst[0]:
                    worst = (marg, (name, a, b, th))
        log(f"     {name}: worst relative margin {fam_worst:.3e}")
    log(f"     total optimisations: {count}; overall worst relative margin {worst[0]:.3e} ({worst[1][0]})")
    check(worst[0] > -1e-9, "NUMERIC no violation of the support inequality of Theorem C.1 found")
    if worst[0] < -1e-9:
        log(f"     worst instance: a={worst[1][1]}, b={worst[1][2]}, theta={worst[1][3]}")

    # common weights for all b at once (supplement's principal-minor route; implies Thm C.1 for that U)
    log("== Stronger form: one p (depending on U, a) with det(D_a+U D_b U*) = sum p_s prod(a_i+b_s(i)) for ALL b")
    def coeff_lp(a, U):
        sq = {(I, J): (abs(np.linalg.det(U[np.ix_(I, J)])) ** 2 if len(I) else 1.0) for (I, J) in MINOR_KEYS}
        rows, rhs = [], []
        for k in range(5):
            for J in SUBSETS[k]:
                lhs = sum(np.prod([a[i] for i in range(4) if i not in I]) * sq[(I, J)] for I in SUBSETS[k])
                row = [np.prod([a[i] for i in range(4) if s[i] not in J]) for s in PERMS]
                rows.append(row); rhs.append(lhs)
        A = np.array(rows, dtype=complex); r = np.array(rhs, dtype=complex)
        A = np.vstack([A.real, A.imag, np.ones((1, 24))]); r = np.concatenate([r.real, r.imag, [1.0]])
        sc = np.max(np.abs(A), axis=1, keepdims=True) + 1e-300
        A, r = A / sc, r / sc[:, 0]
        nE = A.shape[0]
        res = linprog(np.concatenate([np.zeros(24), np.ones(2 * nE)]),
                      A_eq=np.hstack([A, np.eye(nE), -np.eye(nE)]), b_eq=r,
                      bounds=[(0, None)] * (24 + 2 * nE), method="highs")
        return res.fun
    def haar_U():
        Z = nprng.standard_normal((4, 4)) + 1j * nprng.standard_normal((4, 4))
        Q, R = np.linalg.qr(Z)
        return Q * (np.diag(R) / np.abs(np.diag(R)))
    worst_r = worst_c = worst_x = 0.0
    for t in range(SZ(300)):
        a = nprng.standard_normal(4) if t % 3 else 1 + 1e-2 * nprng.standard_normal(4)
        worst_r = max(worst_r, coeff_lp(a, haar_U()))
    for t in range(SZ(150)):
        a = 0.5 + 2 * np.exp(1j * nprng.uniform(0, 2 * np.pi, 4))      # circle of radius 2 centred at 1/2
        worst_c = max(worst_c, coeff_lp(a, haar_U()))
    noncirc = 0
    for t in range(SZ(60)):
        a = nprng.standard_normal(4) + 1j * nprng.standard_normal(4)  # generic: not concyclic
        if coeff_lp(a, haar_U()) > 1e-6:
            noncirc += 1
    check(worst_r < 1e-9, f"NUMERIC 300 Haar U(4), real a (1/3 clustered, seed 4): common-weight LP residual <= {worst_r:.1e}")
    check(worst_c < 1e-9, f"NUMERIC 150 Haar U(4), a on a circle: common-weight LP residual <= {worst_c:.1e}")
    log(f"     info: for generic non-concyclic a, {noncirc}/{SZ(60)} Haar U fail the common-weight LP")
    O0 = np.array([[1, 4, 4, 4], [-4, 1, -4, 4], [-4, 4, 1, -4], [-4, -4, 4, 1]]) / 7.0
    r1 = coeff_lp(np.array([0, 1, 1j, 2]), O0)
    r2 = coeff_lp(np.array([0, 1, 1j, 2]), O0.T)
    r3 = coeff_lp(np.array([0, 1, -1, 2]), O0)
    check(max(r1, r2) > 1e-6 and r3 < 1e-9,
          f"NUMERIC sanity (LP has teeth): explicit orthogonal O0 with non-concyclic a=(0,1,i,2) gives residual "
          f"{r1:.2e} / {r2:.2e} (U=O0 / O0^t), while real a=(0,1,-1,2) gives {r3:.1e}")


if "attack" in SECTIONS:
    run_attack()

log(f"== done; {len(FAILS)} failing checks" + (": " + "; ".join(FAILS) if FAILS else ""))
