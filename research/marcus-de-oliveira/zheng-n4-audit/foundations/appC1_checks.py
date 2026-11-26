#!/usr/bin/env python3
"""Appendix C audit, part 1 (Theorem C.1 statement .. end of Section C.3).

Own code only.  Nothing from the supplement is imported, executed or read: every
object is rebuilt from the definitions printed in the paper (matching basis (A.1),
cubic (A.2), brackets [ij] = z_i - z_j, sigma-points, rho via Cauchy-Binet).

Usage:  python -I -B appC1_checks.py [SECTION ...]
        sections: A B C D E F (default: all).  Deterministic seeds.

Levels printed with each result:
  EXACT    integer / rational / symbolic identity, or exact rank/inertia at an exact point
  NUMERIC  floating point (sample size, seed and tolerance printed)
"""
import sys, itertools, time, math, cmath
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ, QQ

T0 = time.time()


def say(*a):
    print(*a, flush=True)


def hdr(s):
    say("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

# ----------------------------------------------------------------------------
# 0. Core objects rebuilt from the paper
# ----------------------------------------------------------------------------
# (A.1) noncrossing matchings B0..B13, labels 1..8, [ij] = z_i - z_j (Section C.1)
MATCH = [((1, 2), (3, 4), (5, 6), (7, 8)), ((1, 2), (3, 4), (5, 8), (6, 7)),
         ((1, 2), (3, 6), (4, 5), (7, 8)), ((1, 2), (3, 8), (4, 5), (6, 7)),
         ((1, 2), (3, 8), (4, 7), (5, 6)), ((1, 4), (2, 3), (5, 6), (7, 8)),
         ((1, 4), (2, 3), (5, 8), (6, 7)), ((1, 6), (2, 3), (4, 5), (7, 8)),
         ((1, 6), (2, 5), (3, 4), (7, 8)), ((1, 8), (2, 3), (4, 5), (6, 7)),
         ((1, 8), (2, 3), (4, 7), (5, 6)), ((1, 8), (2, 5), (3, 4), (6, 7)),
         ((1, 8), (2, 7), (3, 4), (5, 6)), ((1, 8), (2, 7), (3, 6), (4, 5))]
SUBSETS = list(itertools.combinations(range(1, 9), 4))  # 70 monomials z_T
TIDX = {T: k for k, T in enumerate(SUBSETS)}
PERMS4 = list(itertools.permutations(range(4)))         # lexicographic: pi_0..pi_23


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


def matching_coeffs(pairs):
    """coefficient of z_T (|T|=4) in prod (z_i - z_j) over the pairs."""
    out = {}
    for ch in itertools.product((0, 1), repeat=len(pairs)):
        T = tuple(sorted(pr[c] for pr, c in zip(pairs, ch)))
        out[T] = out.get(T, 0) + (-1) ** sum(ch)
    return out


def T_to_IJ(T):
    """monomial z_T  <->  (I,J): I = [4] minus (T cap [4]), J = (T cap {5..8}) - 4 (0-based)."""
    I = tuple(i - 1 for i in range(1, 5) if i not in T)
    J = tuple(t - 5 for t in T if t > 4)
    return I, J


Bmat = sp.zeros(70, 14)
for k, pr in enumerate(MATCH):
    for T, c in matching_coeffs(pr).items():
        Bmat[TIDX[T], k] = c
Bnp = np.array(Bmat.tolist(), dtype=float)

# P matrix: p_sigma = prod_i [i, 4+sigma(i)]
P_rows = []
for s in PERMS4:
    v = sp.zeros(70, 1)
    for T, c in matching_coeffs(tuple((i + 1, 5 + s[i]) for i in range(4))).items():
        v[TIDX[T]] = c
    P_rows.append(v)
Pmat = sp.zeros(24, 14)
_Bpinv = (Bmat.T * Bmat).inv() * Bmat.T
for r, v in enumerate(P_rows):
    sol = _Bpinv * v
    assert Bmat * sol == v
    Pmat[r, :] = sol.T
Pnp = np.array(Pmat.tolist(), dtype=float)

# 14 rows R of Bmat with invertible block -> rho = Binv rt[R]
def _pick_rows():
    rows, M = [], sp.zeros(0, 14)
    for r in range(70):
        M2 = M.col_join(Bmat[r, :])
        if M2.rank() > M.rank():
            rows.append(r); M = M2
        if len(rows) == 14:
            break
    return rows
RROWS = _pick_rows()
Binv = Bmat.extract(RROWS, list(range(14))).inv()
Binv_np = np.array(Binv.tolist(), dtype=float)
Binv_fr = [[Fr(int(x.p), int(x.q)) for x in Binv.row(i)] for i in range(14)]

# (A.2) cubic, as printed in the paper
msym = sp.symbols('m0:14')
(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13) = msym
C_paper = (m13 * (m0 * m6 - m1 * m5) - m0 * m9 * (m0 + m1 + m2 + m3 + m5 + m7 + m9) - m1 * m5 * m9
           + m10 * (-m0 * m3 - m0 * m9 + m1 * m2 + m1 * m7)
           + m11 * (-m0 * m7 - m0 * m9 + m2 * m5 + m3 * m5)
           + m12 * (-m0 * m9 + m1 * m7 - m2 * m6 - m2 * m9 + m3 * m5 + m3 * m7)
           + m2 * m10 * m11 - m3 * m8 * m10 - m4 * m7 * m11 + m4 * m8 * m9)
# the same cubic as printed (expanded) in the supplement note "one circle", App. A (typed from its text)
y = msym
C_supp = (-y[0]**2*y[9] - y[0]*y[1]*y[9] - y[0]*y[10]*y[3] - y[0]*y[10]*y[9] - y[0]*y[11]*y[7]
          - y[0]*y[11]*y[9] - y[0]*y[12]*y[9] + y[0]*y[13]*y[6] - y[0]*y[2]*y[9] - y[0]*y[3]*y[9]
          - y[0]*y[5]*y[9] - y[0]*y[7]*y[9] - y[0]*y[9]**2 + y[1]*y[10]*y[2] + y[1]*y[10]*y[7]
          + y[1]*y[12]*y[7] - y[1]*y[13]*y[5] - y[1]*y[5]*y[9] + y[10]*y[11]*y[2] - y[10]*y[3]*y[8]
          + y[11]*y[2]*y[5] + y[11]*y[3]*y[5] - y[11]*y[4]*y[7] - y[12]*y[2]*y[6] - y[12]*y[2]*y[9]
          + y[12]*y[3]*y[5] + y[12]*y[3]*y[7] + y[4]*y[8]*y[9])
Cpoly = sp.Poly(C_paper, *msym)
GRAD = [sp.Poly(sp.diff(C_paper, v), *msym) for v in msym]
# Hessian is linear: H(c) = sum_k c_k HK[k]
HK = np.zeros((14, 14, 14))
HKint = [[[0] * 14 for _ in range(14)] for _ in range(14)]
for a in range(14):
    for b in range(14):
        e = sp.Poly(sp.diff(C_paper, msym[a], msym[b]), *msym)
        for mon, co in e.terms():
            if sum(mon) == 0:
                continue
            k = mon.index(1)
            HK[k, a, b] = float(co); HKint[k][a][b] = int(co)


def Hess(c):
    return np.tensordot(c, HK, axes=(0, 0))


def Cval(c):
    return (c @ Hess(c) @ c) / 6.0


def Gval(c):
    return 0.5 * Hess(c) @ c


def Hess_exact(c):
    return sp.Matrix(14, 14, lambda a, b: sum(sp.nsimplify(c[k]) * HKint[k][a][b] for k in range(14)))


def mvec(z):
    """matching vector m(z), z = 8 affine parameters (complex ok)."""
    return np.array([np.prod([z[i - 1] - z[j - 1] for (i, j) in pr]) for pr in MATCH])


def mvec_h(pts):
    """homogeneous version: pts = 8 vectors in C^2, [ij] = det(v_i, v_j)."""
    return np.array([np.prod([pts[i - 1][0] * pts[j - 1][1] - pts[i - 1][1] * pts[j - 1][0]
                              for (i, j) in pr]) for pr in MATCH])


def det_np(M):
    return np.linalg.det(M) if M.shape[0] else 1.0


def rt_np(V):
    """70-vector rt_T = (-1)^{|J|} det V[I,J] det V^{-1}[J,I] (Cauchy-Binet coefficients)."""
    W = np.linalg.inv(V)
    out = np.zeros(70, dtype=complex)
    for k, T in enumerate(SUBSETS):
        I, J = T_to_IJ(T)
        out[k] = (-1) ** len(J) * det_np(V[np.ix_(I, J)]) * det_np(W[np.ix_(J, I)])
    return out


def rho_np(V):
    return Binv_np @ rt_np(V)[RROWS]


def haar_unitary(rng, n=4):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))


def rand_SO4(rng):
    K = rng.standard_normal((4, 4)); K = K - K.T
    from scipy.linalg import expm
    return expm(K)


def sigma_points(a, b):
    return np.array([np.prod([a[i] + b[s[i]] for i in range(4)]) for s in PERMS4])


# ---------------- jets (value, gradient, Hessian) over any field ----------------
class Jet:
    __slots__ = ("c", "g", "H")

    def __init__(self, c, g, H):
        self.c, self.g, self.H = c, g, H

    def __add__(self, o):
        return Jet(self.c + o.c, self.g + o.g, None if self.H is None else self.H + o.H)

    def __sub__(self, o):
        return Jet(self.c - o.c, self.g - o.g, None if self.H is None else self.H - o.H)

    def __mul__(self, o):
        c = self.c * o.c
        g = self.c * o.g + o.c * self.g
        if self.H is None:
            return Jet(c, g, None)
        H = self.c * o.H + o.c * self.H + np.outer(self.g, o.g) + np.outer(o.g, self.g)
        return Jet(c, g, H)

    def scale(self, s):
        return Jet(s * self.c, s * self.g, None if self.H is None else s * self.H)


def jzero(n, order, exact):
    if exact:
        return Jet(Fr(0), np.array([Fr(0)] * n, dtype=object),
                   np.array([[Fr(0)] * n for _ in range(n)], dtype=object) if order == 2 else None)
    return Jet(0j, np.zeros(n, complex), np.zeros((n, n), complex) if order == 2 else None)


def jone(n, order, exact):
    z = jzero(n, order, exact)
    z.c = Fr(1) if exact else 1 + 0j
    return z


PERM_CACHE = {k: [(perm_sign(p), p) for p in itertools.permutations(range(k))] for k in range(5)}


def jdet(M, n, order, exact):
    k = len(M)
    if k == 0:
        return jone(n, order, exact)
    tot = jzero(n, order, exact)
    for s, p in PERM_CACHE[k]:
        t = M[0][p[0]]
        for i in range(1, k):
            t = t * M[i][p[i]]
        tot = tot + t if s > 0 else tot - t
    return tot


def rho_jets(V0, V0inv, Elist, order, exact):
    """jets of rho(V) at V = (I + X) V0, X = sum x_k E_k; V^{-1} = V0^{-1}(I - X + X^2 - ...)."""
    n = len(Elist)
    Vj = [[None] * 4 for _ in range(4)]
    Wj = [[None] * 4 for _ in range(4)]
    EV = [E @ V0 for E in Elist]
    WE = [V0inv @ E for E in Elist]
    if order == 2:
        WEE = [[V0inv @ (Elist[k] @ Elist[l] + Elist[l] @ Elist[k]) for l in range(n)] for k in range(n)]
    for r in range(4):
        for s in range(4):
            if exact:
                g = np.array([EV[k][r, s] for k in range(n)], dtype=object)
                Vj[r][s] = Jet(V0[r, s], g, np.array([[Fr(0)] * n for _ in range(n)], dtype=object) if order == 2 else None)
                g2 = np.array([-WE[k][r, s] for k in range(n)], dtype=object)
                H2 = np.array([[WEE[k][l][r, s] for l in range(n)] for k in range(n)], dtype=object) if order == 2 else None
                Wj[r][s] = Jet(V0inv[r, s], g2, H2)
            else:
                g = np.array([EV[k][r, s] for k in range(n)], dtype=complex)
                Vj[r][s] = Jet(complex(V0[r, s]), g, np.zeros((n, n), complex) if order == 2 else None)
                g2 = np.array([-WE[k][r, s] for k in range(n)], dtype=complex)
                H2 = np.array([[WEE[k][l][r, s] for l in range(n)] for k in range(n)], dtype=complex) if order == 2 else None
                Wj[r][s] = Jet(complex(V0inv[r, s]), g2, H2)
    rt = {}
    for row in RROWS:
        T = SUBSETS[row]
        I, J = T_to_IJ(T)
        d1 = jdet([[Vj[i][j] for j in J] for i in I], n, order, exact)
        d2 = jdet([[Wj[j][i] for i in I] for j in J], n, order, exact)
        p = d1 * d2
        rt[row] = p if len(J) % 2 == 0 else p.scale(-1)
    out = []
    for i in range(14):
        acc = jzero(n, order, exact)
        for jj, row in enumerate(RROWS):
            co = Binv_fr[i][jj] if exact else complex(Binv_np[i, jj])
            if co != 0:
                acc = acc + rt[row].scale(co)
        out.append(acc)
    return out


def frmat(M):
    return np.array([[Fr(x) for x in row] for row in M], dtype=object)


def so_sym_basis():
    so, sym = [], []
    for i in range(4):
        for j in range(i + 1, 4):
            E = np.zeros((4, 4), dtype=int); E[i, j] = 1; E[j, i] = -1; so.append(E)
    for i in range(4):
        E = np.zeros((4, 4), dtype=int); E[i, i] = 1; sym.append(E)
    for i in range(4):
        for j in range(i + 1, 4):
            E = np.zeros((4, 4), dtype=int); E[i, j] = 1; E[j, i] = 1; sym.append(E)
    return so, sym


SO_B, SYM_B = so_sym_basis()


def sym_coords(S):
    return [S[i, i] for i in range(4)] + [S[i, j] for i in range(4) for j in range(i + 1, 4)]


def d_cauchy(O):
    M = [[O[1, j] * O[2, j] for j in range(3)], [O[1, j] * O[3, j] for j in range(3)], [O[2, j] * O[3, j] for j in range(3)]]
    return np.linalg.det(np.array(M))


def num_rank(M, rtol=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > rtol * s[0]).sum()), s


def inertia(A, rtol=1e-9):
    ev = np.linalg.eigvalsh(A)
    t = rtol * max(1e-300, abs(ev).max())
    return int((ev > t).sum()), int((ev < -t).sum()), int((abs(ev) <= t).sum())


def exact_inertia(A):
    """A: sympy rational symmetric matrix; Descartes' rule on the char. polynomial (all roots real)."""
    from flint import fmpz_mat
    den = sp.ilcm(*[sp.fraction(x)[1] for x in A]) if any(sp.fraction(x)[1] != 1 for x in A) else 1
    Ai = (A * den).applyfunc(int)
    cp = fmpz_mat([[int(Ai[i, j]) for j in range(14)] for i in range(14)]).charpoly()
    co = [int(x) for x in cp.coeffs()]  # low -> high
    co += [0] * (15 - len(co))
    z = next(i for i in range(15) if co[i] != 0)
    nz = co[z:]
    def changes(seq):
        s = [x for x in seq if x != 0]
        return sum(1 for u, v in zip(s, s[1:]) if (u > 0) != (v > 0))
    pos = changes(nz)
    neg = changes([c * (-1) ** (i + z) for i, c in enumerate(nz)])
    return pos, neg, z


# ----------------------------------------------------------------------------
# A. exact replays: P, kappa, cubic, X in Sing(C), (C.5), (C.6), Moebius, (C.3), (C.4)
# ----------------------------------------------------------------------------
def section_A():
    hdr("A. Exact replays of Section C.1 identities (own sympy code from (A.1), (A.2))")
    # A1 P, kappa
    allint = all(x.is_integer for x in Pmat)
    sub = Pmat.extract([0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18], list(range(14)))
    say(f"[A1 EXACT] P is integer: {allint}; rank P = {Pmat.rank()}; det of App. A 14x14 row-submatrix = {sub.det()}")
    say(f"[A1 EXACT] column m13 of P is all ones (P kappa = 1 with m(kappa)=e13): {all(Pmat[r, 13] == 1 for r in range(24))}")
    # A2 cubic identity with supplement printing; X in Sing C
    say(f"[A2 EXACT] paper (A.2) cubic == supplement-printed expanded cubic: {sp.expand(C_paper - C_supp) == 0}")
    R8, *zs = ring("z1:9", ZZ)
    mz = []
    for pr in MATCH:
        t = R8.one
        for (i, j) in pr:
            t *= (zs[i - 1] - zs[j - 1])
        mz.append(t)
    def ev_poly(P_, vec, R):
        acc = R.zero
        for mon, co in P_.terms():
            t = R(int(co))
            for k, e in enumerate(mon):
                if e:
                    t *= vec[k] ** e
            acc += t
        return acc
    ok = all(ev_poly(g, mz, R8) == 0 for g in GRAD) and ev_poly(Cpoly, mz, R8) == 0
    say(f"[A2 EXACT] grad C(m(z)) == 0 and C(m(z)) == 0 identically in z1..z8 (X subset Sing C): {ok}")
    # A3 (C.5)
    h = sp.Symbol('h')
    sub5 = {m13: m13 + h}
    lhs = sp.expand(C_paper.subs(sub5) - C_paper - h * (m0 * m6 - m1 * m5))
    say(f"[A3 EXACT] (C.5) C(u + h kappa) = C(u) + h (u0 u6 - u1 u5): {lhs == 0}")
    # (C.6) a real, b complex, theta complex
    ar = sp.symbols('a1:5', real=True)
    br = sp.symbols('br1:5', real=True); bi = sp.symbols('bi1:5', real=True)
    t1, t2 = sp.symbols('t1 t2', real=True)
    zz = list(ar) + [-(br[j] + sp.I * bi[j]) for j in range(4)]
    def mk(k):
        return sp.prod([zz[i - 1] - zz[j - 1] for (i, j) in MATCH[k]])
    th = t1 + sp.I * t2
    Ry = {k: sp.re(sp.expand(th * mk(k))) for k in (0, 1, 5, 6)}
    q = sp.expand(Ry[0] * Ry[6] - Ry[1] * Ry[5])
    say(f"[A3 EXACT] (C.6) q(Re(theta m(a,-b))) == 0 for real a, complex b, theta: {q == 0}")
    # where the circle hypothesis enters: cross-ratio of a1..a4 real <=> q(Re theta m)=0 for all theta
    rng = np.random.default_rng(11)
    vals = {}
    for name in ("concyclic", "collinear-rotated", "generic"):
        worst = 0
        for _ in range(200):
            if name == "concyclic":
                c0 = rng.normal() + 1j * rng.normal(); r = abs(rng.normal()) + .3
                a = c0 + r * np.exp(1j * rng.uniform(0, 2 * np.pi, 4))
            elif name == "collinear-rotated":
                a = (rng.normal() + 1j * rng.normal()) + np.exp(1j * rng.uniform(0, 7)) * rng.normal(size=4)
            else:
                a = rng.normal(size=4) + 1j * rng.normal(size=4)
            b = rng.normal(size=4) + 1j * rng.normal(size=4); thv = rng.normal() + 1j * rng.normal()
            yv = (thv * mvec(np.r_[a, -b])).real
            worst = max(worst, abs(yv[0] * yv[6] - yv[1] * yv[5]) / max(1e-300, np.linalg.norm(yv) ** 2))
        vals[name] = worst
    say(f"[A3 NUMERIC seed 11, 200 each] max |q(Re theta m)|/|y|^2: concyclic a {vals['concyclic']:.1e}, "
        f"collinear a {vals['collinear-rotated']:.1e}, generic a {vals['generic']:.1e}")
    # A4 Moebius
    al, be, ga, de = sp.symbols('alpha beta gamma delta')
    Z = sp.symbols('Z1:9')
    f = lambda u: (al * u + be) / (ga * u + de)
    okm = True
    for k in range(14):
        lhs = sp.prod([f(Z[i - 1]) - f(Z[j - 1]) for (i, j) in MATCH[k]]) * sp.prod([ga * Zi + de for Zi in Z])
        rhs = (al * de - be * ga) ** 4 * sp.prod([Z[i - 1] - Z[j - 1] for (i, j) in MATCH[k]])
        if sp.expand(sp.cancel(sp.together(lhs - rhs))) != 0:
            okm = False
    say(f"[A4 EXACT] Moebius rule m(f(z)) = (ad-bc)^4 / prod(c z_j + d) * m(z), all 14 coords: {okm}")
    # A5 (C.3): rt(V) lies in span(B) for every V.  Jacobi: det V^{-1}[J,I] = (-1)^{s(I)+s(J)} det V[I^c,J^c]/det V
    rng = np.random.default_rng(5)
    V = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); W = np.linalg.inv(V)
    jac_ok = True
    for T in SUBSETS:
        I, J = T_to_IJ(T)
        Ic = tuple(i for i in range(4) if i not in I); Jc = tuple(j for j in range(4) if j not in J)
        sgn = (-1) ** (sum(I) + sum(J))  # 0-based sums have same parity as 1-based when |I|=|J|
        if abs(det_np(W[np.ix_(J, I)]) - sgn * det_np(V[np.ix_(Ic, Jc)]) / np.linalg.det(V)) > 1e-9:
            jac_ok = False
    Vs = sp.Matrix(4, 4, sp.symbols('v0:16'))
    ann = Bmat.T.nullspace()          # 56 vectors annihilating span(B) in Q^70
    vec = []
    for T in SUBSETS:
        I, J = T_to_IJ(T)
        Ic = [i for i in range(4) if i not in I]; Jc = [j for j in range(4) if j not in J]
        d1 = Vs.extract(list(I), list(J)).det() if I else 1
        d2 = Vs.extract(Ic, Jc).det() if Ic else 1
        vec.append((-1) ** len(J) * (-1) ** (sum(I) + sum(J)) * d1 * d2)
    span_ok = all(sp.expand(sum(a_[t] * vec[t] for t in range(70))) == 0 for a_ in ann)
    say(f"[A5 EXACT] Jacobi sign convention validated numerically: {jac_ok}; "
        f"det(V)*rt(V) satisfies all {len(ann)} linear relations of span(B0..B13) as polynomials in 16 entries: {span_ok}")
    say(f"[A5 EXACT] rank of the 70 coefficient functionals n_IJ on V (= rank Bmat) = {Bmat.rank()} (so rho is a linear function of the 70 minor squares; (C.7) follows linearly from Thm 1.2)")
    # numeric (C.3) for unitary U and <rho,kappa>=1, <rho,n_IJ>=|det U[I,J]|^2
    worst = 0; wk = 0; wn = 0
    for _ in range(200):
        U = haar_unitary(rng); a = rng.normal(size=4) + 1j * rng.normal(size=4); b = rng.normal(size=4) + 1j * rng.normal(size=4)
        r = rho_np(U)
        d = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
        worst = max(worst, abs(d - r @ mvec(np.r_[a, -b])) / max(1, abs(d)))
        wk = max(wk, abs(r[13] - 1)); wn = max(wn, abs(r.imag).max())
        rt = rt_np(U)
        for k, T in enumerate(SUBSETS):
            I, J = T_to_IJ(T)
            wn = max(wn, abs((-1) ** len(J) * rt[k] - abs(det_np(U[np.ix_(I, J)])) ** 2))
    say(f"[A5 NUMERIC seed 5, 200 Haar U] (C.3) max rel err {worst:.1e}; |<rho,kappa>-1| {wk:.1e}; rho real and <rho,n_IJ>=|det U[I,J]|^2 max err {wn:.1e}")
    # A6 (C.4) in 16 variables, all 70 pairings
    t0 = time.time()
    R16, *g16 = ring("z1:9,w1:9", ZZ)
    zv, wv = g16[:8], g16[8:]
    def mring(v):
        out = []
        for pr in MATCH:
            t = R16.one
            for (i, j) in pr:
                t *= (v[i - 1] - v[j - 1])
            out.append(t)
        return out
    S = [a + b for a, b in zip(mring(zv), mring(wv))]
    LHS = [ev_poly(g, S, R16) for g in GRAD]
    A = [[R16.one, zv[i], wv[i], zv[i] * wv[i]] for i in range(8)]  # columns
    def det4(cols):
        tot = R16.zero
        for s, p in PERM_CACHE[4]:
            t = R16.one
            for r in range(4):
                t *= cols[p[r]][r]
            tot = tot + t if s > 0 else tot - t
        return tot
    A0 = [A[i] for i in range(4)]; A1 = [A[4 + i] for i in range(4)]
    ok70 = True; bad = []
    for k, T in enumerate(SUBSETS):
        I, J = T_to_IJ(T)
        c1 = list(A0)
        for pos, j in zip(I, J):
            c1[pos] = A1[j]
        c2 = list(A1)
        for pos, i in zip(J, I):
            c2[pos] = A0[i]
        rhs = det4(c1) * det4(c2)
        if len(J) % 2:
            rhs = -rhs
        lhs = R16.zero
        for kk in range(14):
            if Bmat[k, kk] != 0:
                lhs += int(Bmat[k, kk]) * LHS[kk]
        if lhs != rhs:
            ok70 = False; bad.append(T)
    say(f"[A6 EXACT] (C.4) G(m(z)+m(w)) = det(A0)det(A1) rho(A0^-1 A1): all 70 coefficient pairings agree as "
        f"polynomials in z1..z8,w1..w8: {ok70} (mismatches: {bad[:3]}) [{time.time() - t0:.1f}s]")


# ----------------------------------------------------------------------------
# B. Lemma C.2 (n=3) and Lemma C.3 (zero determinant)
# ----------------------------------------------------------------------------
def hull_margin(pts):
    """signed distance-type margin of 0 w.r.t. conv(pts) in C: >0 iff 0 outside (up to rounding)."""
    from scipy.spatial import ConvexHull
    from scipy.optimize import linprog
    P2 = np.c_[pts.real, pts.imag]
    sc = max(1e-300, np.abs(pts).max())
    try:
        hull = ConvexHull(P2 / sc)
        return float(hull.equations[:, 2].max())
    except Exception:
        # degenerate (collinear) set: minimize the L1 norm of a convex combination (LP); 0 => inside
        n = len(pts)
        c = np.r_[np.zeros(n), 1, 1]
        Aeq = np.zeros((3, n + 2)); beq = np.array([0, 0, 1.0])
        Aub = np.zeros((4, n + 2)); bub = np.zeros(4)
        X = P2 / sc
        Aub[0, :n] = X[:, 0]; Aub[0, n] = -1
        Aub[1, :n] = -X[:, 0]; Aub[1, n] = -1
        Aub[2, :n] = X[:, 1]; Aub[2, n + 1] = -1
        Aub[3, :n] = -X[:, 1]; Aub[3, n + 1] = -1
        Aeq[2, :n] = 1
        r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq[2:], b_eq=beq[2:], bounds=[(0, None)] * (n + 2), method="highs")
        return float(r.fun) if r.success else float("nan")


def section_B():
    hdr("B. Lemma C.2 (n=3) and Lemma C.3 (zero-determinant case)")
    from scipy.optimize import linprog, minimize
    rng = np.random.default_rng(2026)
    # C.2 expansion identity and inclusion
    perms3 = list(itertools.permutations(range(3)))
    werr = 0; wm = -1
    for _ in range(500):
        Z = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
        U, _r = np.linalg.qr(Z)
        a = rng.normal(size=3) + 1j * rng.normal(size=3); b = rng.normal(size=3) + 1j * rng.normal(size=3)
        p = abs(U) ** 2
        d = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
        f = np.prod(a) + np.prod(b) + sum(p[i, j] * (b[j] * np.prod([a[k] for k in range(3) if k != i]) +
                                                     a[i] * np.prod([b[l] for l in range(3) if l != j])) for i in range(3) for j in range(3))
        werr = max(werr, abs(d - f) / max(1, abs(d)))
        zs = np.array([np.prod([a[i] + b[s[i]] for i in range(3)]) for s in perms3])
        wm = max(wm, hull_margin(np.r_[zs - d]))
    say(f"[B1 NUMERIC seed 2026, 500] Lemma C.2 expansion identity max rel err {werr:.1e}; max hull margin of det wrt sigma-points {wm:.1e} (<=0 means inside)")
    # C.3 random tests with constructive replay
    stats = dict(trials=0, maxmargin=-9, maxsupport=0, single=0, three=0, Werr=0, deterr=0, sub_margin=-9)
    fams = ["generic", "hermitianA", "repeatedA", "collinearB", "smallscale", "unitaryA"]
    for fam in fams:
        for _ in range(300):
            a = rng.normal(size=4) + 1j * rng.normal(size=4)
            if fam == "hermitianA":
                a = a.real
            if fam == "repeatedA":
                a[1] = a[0]
            if fam == "unitaryA":
                a = np.exp(1j * rng.uniform(0, 2 * np.pi, 4))
            b = rng.normal(size=4) + 1j * rng.normal(size=4)
            if fam == "collinearB":
                b = (1 + 2j) * rng.normal(size=4) + 0.3j
            if fam == "smallscale":
                a *= 1e-3; b *= 1e-3
            U = haar_unitary(rng)
            def dt(b4):
                bb = b.copy(); bb[3] = b4
                return np.linalg.det(np.diag(a) + U @ np.diag(bb) @ U.conj().T)
            f0, f1 = dt(0), dt(1)
            if abs(f1 - f0) < 1e-12:
                continue
            b[3] = -f0 / (f1 - f0)
            M = np.diag(a) + U @ np.diag(b) @ U.conj().T
            zs = sigma_points(a, b)
            stats["trials"] += 1
            stats["maxmargin"] = max(stats["maxmargin"], hull_margin(zs))
            # constructive proof replay
            _, sv, Vh = np.linalg.svd(M); v = Vh[-1].conj()
            p = abs(v) ** 2
            e = -b; q = abs(U.conj().T @ v) ** 2
            Aeq = np.array([np.r_[np.ones(4), np.zeros(4)], np.r_[np.zeros(4), np.ones(4)],
                            np.r_[a.real, -e.real], np.r_[a.imag, -e.imag], np.r_[abs(a) ** 2, -abs(e) ** 2]])
            beq = np.array([1, 1, 0, 0, 0.0])
            assert np.abs(Aeq @ np.r_[p, q] - beq).max() < 1e-8
            r = linprog(rng.normal(size=8), A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * 8, method="highs-ds")
            if not r.success:
                r = linprog(np.zeros(8), A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * 8, method="highs-ds")
            w = r.x; pp, qq = w[:4], w[4:]
            sup_p = [i for i in range(4) if pp[i] > 1e-10]; sup_q = [j for j in range(4) if qq[j] > 1e-10]
            stats["maxsupport"] = max(stats["maxsupport"], len(sup_p) + len(sup_q))
            if len(sup_p) == 1 or len(sup_q) == 1:
                stats["single"] += 1
                if len(sup_p) == 1:
                    i0 = sup_p[0]; gap = min(abs(e[j] - a[i0]) for j in sup_q)
                else:
                    j0 = sup_q[0]; gap = min(abs(e[j0] - a[i]) for i in sup_p)
                stats["deterr"] = max(stats["deterr"], gap / max(1, abs(a).max()))
                continue
            stats["three"] += 1
            I = (sup_p + [i for i in range(4) if i not in sup_p])[:3]
            J = (sup_q + [j for j in range(4) if j not in sup_q])[:3]
            x = np.sqrt(np.maximum(pp[I], 0)); yv = np.sqrt(np.maximum(qq[J], 0))
            Da = np.diag(a[I]); Dg = np.diag(e[J])
            # unitary W with W y = x, W Dg y = Da x (Gram matrices equal)
            Ym = np.c_[yv, Dg @ yv]; Xm = np.c_[x, Da @ x]
            Qy, Ry = np.linalg.qr(np.c_[Ym, rng.normal(size=(3, 1)) + 1j * rng.normal(size=(3, 1))])
            Qx, Rx = np.linalg.qr(np.c_[Xm, rng.normal(size=(3, 1)) + 1j * rng.normal(size=(3, 1))])
            # align phases so that the first two columns match (R factors equal on 2x2 block)
            ph = np.diag(Rx)[:2] / np.diag(Ry)[:2]
            ph = ph / abs(ph)
            Qy2 = Qy.copy(); Qy2[:, :2] = Qy[:, :2] / ph
            Wm = Qx @ Qy2.conj().T
            stats["Werr"] = max(stats["Werr"], np.abs(Wm @ yv - x).max(), np.abs(Wm @ Dg @ yv - Da @ x).max())
            Mn = Da - Wm @ Dg @ Wm.conj().T
            stats["deterr"] = max(stats["deterr"], abs(np.linalg.det(Mn)) / max(1, np.abs(Mn).max() ** 3))
            z3 = np.array([np.prod([a[I[rr]] - e[J[s[rr]]] for rr in range(3)]) for s in itertools.permutations(range(3))])
            stats["sub_margin"] = max(stats["sub_margin"], hull_margin(z3))
    say(f"[B2 NUMERIC seed 2026] Lemma C.3: {stats['trials']} random normal pairs with det(A+B)=0 (families {fams})")
    say(f"      max hull margin of 0 wrt the 24 sigma-points: {stats['maxmargin']:.2e} (<= 0 means 0 in hull)")
    say(f"      proof replay: extreme point support size max {stats['maxsupport']} (claim <= 5); single-support branch {stats['single']}x, "
        f"3x3 branch {stats['three']}x; max |Wy-x|,|W Dg y - Da x| {stats['Werr']:.1e}; max zero-det/coincidence residual {stats['deterr']:.1e}; "
        f"max 3x3 hull margin {stats['sub_margin']:.1e}")
    # adversarial search: maximize hull margin over det=0 configurations
    def margin_of(par):
        a = par[0:4] + 1j * par[4:8]; b = par[8:12] + 1j * par[12:16]
        K = par[16:32].reshape(4, 4); H = (K + K.T) / 2 + 1j * (K - K.T) / 2
        from scipy.linalg import expm
        U = expm(1j * H)
        def dt(b4):
            bb = b.copy(); bb[3] = b4
            return np.linalg.det(np.diag(a) + U @ np.diag(bb) @ U.conj().T)
        f0, f1 = dt(0), dt(1)
        if abs(f1 - f0) < 1e-9:
            return 0.0
        bb = b.copy(); bb[3] = -f0 / (f1 - f0)
        return hull_margin(sigma_points(a, bb))
    best = -9
    for st in range(25):
        x0 = rng.normal(size=32)
        res = minimize(lambda p: -margin_of(p), x0, method="Nelder-Mead", options=dict(maxiter=1500, xatol=1e-8, fatol=1e-12))
        best = max(best, -res.fun)
    say(f"[B3 NUMERIC seed 2026] adversarial Nelder-Mead (25 starts x 1500 it) maximizing the margin over det(A+B)=0: best margin {best:.2e}")


# ----------------------------------------------------------------------------
# C. Corollary C.5 and Proposition C.6
# ----------------------------------------------------------------------------
def principal12(A):
    idx = list(range(14)); out = []
    for drop in itertools.combinations(idx, 2):
        keep = [i for i in idx if i not in drop]
        out.append(np.linalg.det(A[np.ix_(keep, keep)]))
    return np.array(out)


def section_C():
    hdr("C. Corollary C.5 (rank <= 12 on the secant cone) and Proposition C.6 (inertia (6,6,2))")
    rng = np.random.default_rng(606)
    # C1 rank H on X, on Sigma2, on Re x
    rX = []; rS = []; rR = []
    for _ in range(300):
        z = rng.normal(size=8) + 1j * rng.normal(size=8); w = rng.normal(size=8) + 1j * rng.normal(size=8)
        lam = rng.normal() + 1j * rng.normal()
        x = lam * mvec(z)
        rX.append(num_rank(Hess(x))[0]); rS.append(num_rank(Hess(x + mvec(w)))[0]); rR.append(num_rank(Hess(x.real))[0])
    say(f"[C1 NUMERIC seed 606, 300, rtol 1e-9] rank H(x), x in X: {sorted(set(rX))}; rank H(x+x'), x,x' in X: {sorted(set(rS))}; rank H(Re x): {sorted(set(rR))}")
    # explicit kernel vector s x - t y on c = s x + t y
    z = rng.normal(size=8) + 1j * rng.normal(size=8); w = rng.normal(size=8) + 1j * rng.normal(size=8)
    x, yv = mvec(z), mvec(w); s, t = 0.7 - .2j, 1.3 + .5j
    c = s * x + t * yv
    say(f"[C1 NUMERIC] H(sx+ty)(sx-ty) relative size {np.linalg.norm(Hess(c) @ (s * x - t * yv)) / (np.linalg.norm(Hess(c)) * np.linalg.norm(s * x - t * yv)):.1e} (explicit kernel vector on every secant)")
    # dimensions: Sigma2 cone 12, G(Sigma2) cone 10, rank of d(rho) image 10
    def jac(fun, p, h=1e-6):
        f0 = fun(p); J = np.zeros((len(f0), len(p)), complex)
        for k in range(len(p)):
            e = np.zeros(len(p), complex); e[k] = h
            J[:, k] = (fun(p + e) - fun(p - e)) / (2 * h)
        return J
    p0 = rng.normal(size=16) + 1j * rng.normal(size=16)
    J1 = jac(lambda p: mvec(p[:8]) + mvec(p[8:]), p0)
    J2 = jac(lambda p: Gval(mvec(p[:8]) + mvec(p[8:])), p0)
    say(f"[C1 NUMERIC] dim of secant cone (Jacobian rank of (z,w)->m(z)+m(w)) = {num_rank(J1, 1e-7)[0]} (claim 12); "
        f"rank of d(G o secant) = {num_rank(J2, 1e-7)[0]} (claim 10, => rank H(c) <= 12 on Sigma2)")
    # C2 inertia of H(Re x): random generic + degenerate families
    fams = {}
    def sample(fam):
        z = rng.normal(size=8) + 1j * rng.normal(size=8); lam = rng.normal() + 1j * rng.normal()
        if fam == "a1=a2": z[1] = z[0]
        if fam == "b1=b2": z[5] = z[4]
        if fam == "a1+b1=0": z[4] = z[0]
        if fam == "a real": z[:4] = z[:4].real
        if fam == "a concyclic": z[:4] = 0.3 + 1.7 * np.exp(1j * rng.uniform(0, 6.3, 4))
        if fam == "a 3 distinct": z[3] = z[2]
        if fam == "all real (x real)": z = z.real; lam = lam.real
        if fam == "a1=a2=a3": z[1] = z[0]; z[2] = z[0]
        if fam == "b const": z[5:8] = z[4]
        if fam == "huge spread": z = z * np.array([1e3, 1, 1e-2, 1, 1, 1e2, 1, 1])
        return lam * mvec(z)
    tot = {}
    for fam in ["generic", "a1=a2", "b1=b2", "a1+b1=0", "a real", "a concyclic", "a 3 distinct",
                "a1=a2=a3", "b const", "all real (x real)", "huge spread"]:
        n = 600 if fam == "generic" else 150
        cnt = {}; minpm = np.inf; maxpos = -1; maxneg = -1; rank12_haspos = True
        for _ in range(n):
            x = sample(fam); A = Hess(x.real)
            ine = inertia(A); cnt[ine] = cnt.get(ine, 0) + 1
            maxpos = max(maxpos, ine[0]); maxneg = max(maxneg, ine[1])
            if fam in ("generic", "a1=a2", "a real", "huge spread") and _ < 60:
                pm = principal12(A / np.abs(A).max())
                minpm = min(minpm, pm.min() / max(1e-300, abs(pm).max()))
                if ine[2] == 2 and pm.max() <= 0:
                    rank12_haspos = False
        tot[fam] = (cnt, maxpos, maxneg, minpm, rank12_haspos)
        say(f"[C2 NUMERIC seed 606, n={n}] {fam:18s} inertia counts {dict(sorted(cnt.items()))}; max n+ {maxpos}, max n- {maxneg}"
            + (f"; min normalized 12x12 principal minor {minpm:.1e}; rank-12 => some minor >0: {rank12_haspos}" if minpm < np.inf else ""))
    # C3 exact points: rank H(x)=8, rank [A;B]=14 (ker H(x) cap ker H(xbar)=0), exact inertia of A
    say("[C3 EXACT] Gaussian-rational sample points (exact rational arithmetic):")
    from flint import fmpq_mat, fmpq
    def cm(u, v):  # Gaussian-rational multiply, pairs (re, im) of Fractions
        return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])
    def Hq(c):
        return [[sum(c[k] * HKint[k][a][b] for k in range(14)) for b in range(14)] for a in range(14)]
    def fq(x):
        return fmpq(x.numerator, x.denominator)
    for seed in (1, 2, 3):
        r2 = np.random.default_rng(seed)
        zr = [(Fr(int(r2.integers(-9, 10)), int(r2.integers(1, 5))), Fr(int(r2.integers(-9, 10)), int(r2.integers(1, 5)))) for _ in range(8)]
        lam = (Fr(int(r2.integers(1, 5))), Fr(int(r2.integers(-3, 4))))
        xs = []
        for pr in MATCH:
            t = lam
            for (i, j) in pr:
                t = cm(t, (zr[i - 1][0] - zr[j - 1][0], zr[i - 1][1] - zr[j - 1][1]))
            xs.append(t)
        A = Hq([v[0] for v in xs]); B = Hq([v[1] for v in xs])
        big = [[fq(A[i][j]) for j in range(14)] + [fq(-B[i][j]) for j in range(14)] for i in range(14)] + \
              [[fq(B[i][j]) for j in range(14)] + [fq(A[i][j]) for j in range(14)] for i in range(14)]
        rk_x = fmpq_mat(big).rank() // 2          # complex rank via real embedding
        rkAB = fmpq_mat([[fq(v) for v in row] for row in A + B]).rank()
        ine = exact_inertia(sp.Matrix(14, 14, lambda i, j: sp.Rational(A[i][j].numerator, A[i][j].denominator)))
        say(f"      seed {seed}: rank H(x) = {rk_x} (claim 8); rank [H(Re x); H(Im x)] = {rkAB} (claim 14: transversality); "
            f"exact inertia of H(Re x) (n+, n-, n0) = {ine} (claim (6,6,2))")


# ----------------------------------------------------------------------------
# D. Proposition C.4: conormal spaces of M8' vs the secant variety
# ----------------------------------------------------------------------------
def E_elem():
    out = []
    for r in range(4):
        for s in range(4):
            E = np.zeros((4, 4), dtype=int); E[r, s] = 1; out.append(E)
    return out


def tangent_cone_np(V):
    Vinv = np.linalg.inv(V)
    js = rho_jets(V, Vinv, E_elem(), 1, False)
    rho = np.array([j.c for j in js]); D = np.array([j.g for j in js])  # 14 x 16
    return rho, D


def conormal_np(V):
    rho, D = tangent_cone_np(V)
    M = np.c_[rho, D]                     # 14 x 17, columns span the affine tangent of the cone
    from scipy.linalg import null_space
    L = null_space(M.T, rcond=1e-10)                   # bilinear annihilator: c with sum_j c_j M_jk = 0
    return rho, D, M, L, None


def section_D():
    hdr("D. Proposition C.4: (M8')^* = Sigma2 (conormal spaces of the Gale quotient)")
    from scipy.linalg import null_space
    rng = np.random.default_rng(404)
    mons4 = [(i, j) for i in range(4) for j in range(i, 4)]
    res = dict(rk=[], Cres=0, Hrank=[], qrank=[], memb=0, polar=0)
    t0 = time.time()
    for trial in range(40):
        V = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        rho, D, M, L, s = conormal_np(V)
        res["rk"].append(num_rank(M, 1e-9)[0])
        # C on L, rank H on L, restricted quadrics
        for _ in range(3):
            lam = rng.normal(size=L.shape[1]) + 1j * rng.normal(size=L.shape[1])
            c = L @ lam; c = c / np.linalg.norm(c)
            res["Cres"] = max(res["Cres"], abs(Cval(c)) / np.linalg.norm(Hess(c)))
            res["Hrank"].append(num_rank(Hess(c), 1e-8)[0])
        Q = np.zeros((14, len(mons4)), complex)
        for k in range(14):
            Hk = L.T @ HK[k] @ L  # not used directly
        # quadric G_j(L lam) coefficients via polarization
        Gq = np.zeros((14, len(mons4)), complex)
        for (i, j), col in zip(mons4, range(len(mons4))):
            ei = L[:, i]; ej = L[:, j]
            if i == j:
                Gq[:, col] = Gval(ei)
            else:
                Gq[:, col] = Gval(ei + ej) - Gval(ei) - Gval(ej)
        res["qrank"].append(num_rank(Gq, 1e-9)[0])
        # elliptic quartic: pencil of quadrics through the 8 points [I, V]; rulings -> P^1 octads f(A)
        pts = np.c_[np.eye(4), V].T                       # 8 points in C^4
        Qm = np.array([[p[i] * p[j] * (1 if i == j else 2) for (i, j) in mons4] for p in pts])
        pen = null_space(Qm, rcond=1e-10)                              # 10 x 2
        PL = L @ np.linalg.pinv(L)
        octads = []
        for _ in range(3):
            lam = rng.normal() + 1j * rng.normal()
            qv = pen[:, 0] + lam * pen[:, 1]
            Mq = np.zeros((4, 4), complex)
            for (i, j), v in zip(mons4, qv):
                Mq[i, j] = Mq[j, i] = v
            ev, Ev = np.linalg.eig(Mq)
            Ev = Ev / np.sqrt(np.sum(Ev * Ev, axis=0))   # complex-orthogonal normalization
            Rq = np.diag(np.sqrt(ev)) @ Ev.T              # Mq = Rq^T Rq
            uu = (Rq @ pts.T).T
            yy = np.c_[uu[:, 0] + 1j * uu[:, 1], uu[:, 2] + 1j * uu[:, 3], -(uu[:, 2] - 1j * uu[:, 3]), uu[:, 0] - 1j * uu[:, 1]]
            r1 = [np.array([yy[k, 0], yy[k, 2]]) if abs(yy[k, 0]) + abs(yy[k, 2]) > 1e-8 else np.array([yy[k, 1], yy[k, 3]]) for k in range(8)]
            r2 = [np.array([yy[k, 0], yy[k, 1]]) if abs(yy[k, 0]) + abs(yy[k, 1]) > 1e-8 else np.array([yy[k, 2], yy[k, 3]]) for k in range(8)]
            m1v, m2v = mvec_h(r1), mvec_h(r2)
            for mv in (m1v, m2v):
                res["memb"] = max(res["memb"], np.linalg.norm(mv - PL @ mv) / np.linalg.norm(mv))
            g = Gval(m1v + m2v)
            cosang = abs(np.vdot(g, rho)) / (np.linalg.norm(g) * np.linalg.norm(rho))
            res["polar"] = max(res["polar"], 1 - cosang)
    say(f"[D1 NUMERIC seed 404, 40 random complex V, rtol 1e-9] rank of [rho, d rho] (affine tangent of cone over M8') = {sorted(set(res['rk']))} (claim 10, so conormal = P^3)")
    say(f"      on the conormal space L(V): max |C(c)|/|H(c)| = {res['Cres']:.1e}; rank H(c) = {sorted(set(res['Hrank']))} (<=12); "
        f"the 14 quadrics grad C restricted to L span a space of dim {sorted(set(res['qrank']))} (2 = pencil => X cap P(L) is a quartic curve, consistent with an elliptic normal quartic)")
    say(f"      octads f(A) from both rulings of 3 quadrics of the pencil through [I,V]: max distance of m(f(A)) from L = {res['memb']:.1e}; "
        f"G(m(ruling1)+m(ruling2)) parallel to rho(V): max (1-|cos|) = {res['polar']:.1e}")
    # dimension of the dual variety: (V,lambda) -> L(V) lambda
    def Lnorm(V):
        _, _, _, L, _ = conormal_np(V)
        return L @ np.linalg.inv(L[[0, 3, 7, 11], :])
    V0 = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); l0 = rng.normal(size=4) + 1j * rng.normal(size=4)
    def F(p):
        return Lnorm(p[:16].reshape(4, 4)) @ p[16:]
    p0 = np.r_[V0.ravel(), l0]; h = 1e-6
    J = np.zeros((14, 20), complex)
    for k in range(20):
        e = np.zeros(20, complex); e[k] = h
        J[:, k] = (F(p0 + e) - F(p0 - e)) / (2 * h)
    rk, sv = num_rank(J, 1e-7)
    say(f"[D2 NUMERIC] affine dimension of the union of conormal spaces (Jacobian rank of (V,lambda) -> c) = {rk} "
        f"(claim 12, i.e. dim (M8')^* = 11 = dim Sigma2); singular values tail {sv[10:14]}")
    # exact at a rational V
    t0 = time.time()
    r3 = np.random.default_rng(7)
    while True:
        Vi = r3.integers(-3, 4, size=(4, 4))
        if round(np.linalg.det(Vi)) != 0 and np.all(Vi != 0):
            break
    Vq = sp.Matrix(Vi.tolist()); Vinvq = Vq.inv()
    V0f = frmat(Vi.tolist()); Vinvf = np.array([[Fr(int(x.p), int(x.q)) for x in Vinvq.row(i)] for i in range(4)], dtype=object)
    js = rho_jets(V0f, Vinvf, [frmat(E.tolist()) for E in E_elem()], 1, True)
    Mq = sp.Matrix([[sp.Rational(j.c.numerator, j.c.denominator)] + [sp.Rational(x.numerator, x.denominator) for x in j.g] for j in js])
    rkq = Mq.rank()
    Lq = Mq.T.nullspace()
    lam = sp.symbols('l0:4')
    cvec = sum((lam[i] * Lq[i] for i in range(len(Lq))), sp.zeros(14, 1))
    Csub = sp.expand(C_paper.subs(dict(zip(msym, list(cvec)))))
    Gsub = [sp.Poly(sp.expand(g.as_expr().subs(dict(zip(msym, list(cvec))))), *lam) for g in GRAD]
    monos = sorted({mm for g in Gsub for mm in g.monoms()})
    Gm = sp.Matrix([[g.coeff_monomial(mm) for mm in monos] for g in Gsub])
    cpt = cvec.subs(dict(zip(lam, [1, 2, -3, 5])))
    from flint import fmpq_mat, fmpq
    cq = [sp.Rational(v) for v in cpt]
    rkH = fmpq_mat([[sum((fmpq(int(cq[k].p), int(cq[k].q)) * HKint[k][a][b] for k in range(14)), fmpq(0)) for b in range(14)] for a in range(14)]).rank()
    say(f"[D3 EXACT] rational V = {Vi.tolist()}: rank [rho, d rho] = {rkq}; dim conormal L = {len(Lq)}; C|L == 0: {Csub == 0}; "
        f"rank of the 14 restricted quadrics = {Gm.rank()}; rank H(c) at c = L(1,2,-3,5) = {rkH} [{time.time() - t0:.1f}s]")


# ----------------------------------------------------------------------------
# E. Section C.3: K_O, phase Hessian S_O, Lemma C.7 injectivity, (C.8)
# ----------------------------------------------------------------------------
def torus_dirs(O):
    out = []
    for i in range(4):
        D = np.zeros((4, 4)); D[i, i] = 1; out.append(sym_coords(D))
    for j in range(4):
        D = np.zeros((4, 4)); D[j, j] = 1; out.append(sym_coords(O @ D @ O.T))
    return np.array(out, dtype=complex).T  # 10 x 8


def phase_data_np(O):
    El = SO_B + SYM_B
    js = rho_jets(O, O.T, El, 2, False)
    rho0 = np.array([j.c for j in js]); grad = np.array([j.g for j in js]); Hs = np.array([j.H for j in js])
    return rho0, grad, Hs


def section_E():
    hdr("E. Section C.3: K_O (dim 7), phase Hessian S_O, injectivity of B in Lemma C.7, identity (C.8)")
    from scipy.linalg import null_space, expm
    rng = np.random.default_rng(707)
    out = dict(rk7=[], symgrad=0, torus=0, mixed=0, rkS=[], rkB=[], c8=0, alpha=[])
    for trial in range(30):
        if trial < 20:
            O = rand_SO4(rng)
        else:  # complex orthogonal
            K = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); K = (K - K.T) / 2
            O = expm(K)
        rho0, grad, Hs = phase_data_np(O)
        Ao = grad[:, :6]; Asym = grad[:, 6:]
        out["symgrad"] = max(out["symgrad"], np.abs(Asym).max() / np.abs(Ao).max())
        M7 = np.c_[rho0, Ao]
        out["rk7"].append(num_rank(M7)[0])
        K_O = null_space(M7.T, rcond=1e-10)            # 14 x 7, bilinear annihilator
        Tor = torus_dirs(O)
        Qc = null_space(Tor.T, rcond=1e-10)            # 10 x 3 complement (w.r.t. bilinear form) of torus directions
        Bvec = []                         # B(e_i o e_j) = 1/2 d_i d_j rho along the 3 phase directions
        for i in range(3):
            for j in range(i, 3):
                Bvec.append(0.5 * np.einsum('kab,a,b->k', Hs[:, 6:, 6:], Qc[:, i], Qc[:, j]))
        out["rkB"].append(num_rank(np.c_[rho0, Ao, np.array(Bvec).T])[0])
        def S_of(c):
            Hc = np.einsum('k,kab->ab', c, Hs)
            return Hc
        Smaps = []
        for col in range(7):
            c = K_O[:, col]; Hc = S_of(c)
            out["torus"] = max(out["torus"], np.abs(Hc[6:, 6:] @ Tor).max() / max(1e-300, np.abs(Hc).max()))
            out["mixed"] = max(out["mixed"], np.abs(Hc[:6, 6:]).max() / max(1e-300, np.abs(Hc).max()))
            S3 = Qc.T @ Hc[6:, 6:] @ Qc
            Smaps.append([S3[i, j] for i in range(3) for j in range(i, 3)])
        out["rkS"].append(num_rank(np.array(Smaps))[0])
        Cs = []; Ds = []
        for _ in range(25):
            t = rng.normal(size=7) + 1j * rng.normal(size=7)
            c = K_O @ t
            S3 = Qc.T @ S_of(c)[6:, 6:] @ Qc
            Cs.append(Cval(c)); Ds.append(np.linalg.det(S3))
        Cs = np.array(Cs); Ds = np.array(Ds)
        al = np.vdot(Ds, Cs) / np.vdot(Ds, Ds)
        out["c8"] = max(out["c8"], np.linalg.norm(Cs - al * Ds) / np.linalg.norm(Cs))
        out["alpha"].append(abs(al))
    say(f"[E1 NUMERIC seed 707; 20 real SO(4) + 10 complex SO(4,C), rtol 1e-9]")
    say(f"      rank[rho(O), d_SO rho] = {sorted(set(out['rk7']))} (claim 7 => dim K_O = 7)")
    say(f"      d rho along symmetric directions S O vanishes (Gale evenness): max rel {out['symgrad']:.1e}")
    say(f"      rank[rho0, A, B(Sym^2 C^3)] = {sorted(set(out['rkB']))} (claim 13: Lemma C.7 injectivity of B mod C rho0 + im A)")
    say(f"      for c in K_O: torus directions D1+O D2 O^t in kernel of the Hessian (max rel {out['torus']:.1e}); "
        f"mixed so x sym Hessian block vanishes (max rel {out['mixed']:.1e})")
    say(f"      rank of c -> S_O(c) on K_O = {sorted(set(out['rkS']))} (claim 6: surjective onto Sym^2)")
    say(f"      (C.8) C|K_O = alpha(O) det S_O: max relative residual of the best alpha over 25 points = {out['c8']:.1e}; "
        f"|alpha| range [{min(out['alpha']):.2e}, {max(out['alpha']):.2e}] (phase basis is arbitrary, alpha only up to a unit)")
    # exact version at a rational Cayley orthogonal matrix and at the Lemma C.9 example matrix
    def exact_C8(Oq, label):
        t0 = time.time()
        Of = np.array([[Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in Oq.row(i)] for i in range(4)], dtype=object)
        El = [frmat(E.tolist()) for E in SO_B + SYM_B]
        js = rho_jets(Of, Of.T.copy(), El, 2, True)
        Q = lambda x: sp.Rational(x.numerator, x.denominator)
        rho0 = sp.Matrix([Q(j.c) for j in js]); grad = sp.Matrix([[Q(x) for x in j.g] for j in js])
        Hs = [sp.Matrix(16, 16, lambda a, b: Q(j.H[a, b])) for j in js]
        M7 = rho0.row_join(grad[:, :6])
        symg = grad[:, 6:].is_zero_matrix
        KO = M7.T.nullspace()
        Tor = sp.Matrix([sym_coords(sp.Matrix(4, 4, lambda r, s: 1 if (r == s == i) else 0)) for i in range(4)] +
                        [sym_coords(Oq * sp.Matrix(4, 4, lambda r, s: 1 if (r == s == j) else 0) * Oq.T) for j in range(4)]).T
        Qc = Tor.T.nullspace()
        Qc = sp.Matrix.hstack(*Qc)
        tt = sp.symbols('t0:%d' % len(KO))
        c = sum((tt[i] * KO[i] for i in range(len(KO))), sp.zeros(14, 1))
        Hc = sum((c[k] * Hs[k] for k in range(14)), sp.zeros(16, 16))
        torus_ok = sp.expand(Hc[6:, 6:] * Tor).is_zero_matrix
        S3 = (Qc.T * Hc[6:, 6:] * Qc).applyfunc(sp.expand)
        dS = sp.Poly(sp.expand(S3.det()), *tt)
        Cc = sp.Poly(sp.expand(C_paper.subs(dict(zip(msym, list(c))))), *tt)
        # alpha = ratio of leading coefficients; then compare all coefficients
        mon = dS.monoms()[0]
        alpha = Cc.coeff_monomial(mon) / dS.coeff_monomial(mon)
        ident = (Cc - dS * alpha).is_zero
        say(f"[E2 EXACT] {label}: dim K_O = {len(KO)}; sym-gradient zero: {symg}; torus in Hessian kernel: {torus_ok}; "
            f"C|K_O == alpha * det S_O as cubic forms in 7 variables: {ident}; alpha = {alpha} [{time.time() - t0:.1f}s]")
    r5 = np.random.default_rng(55)
    while True:   # rational Cayley orthogonal matrix with all entries nonzero and d(O) != 0
        k6 = r5.integers(-4, 5, size=6)
        Kq = sp.zeros(4, 4); t_ = 0
        for i in range(4):
            for j in range(i + 1, 4):
                Kq[i, j] = sp.Rational(int(k6[t_]), 3); Kq[j, i] = -Kq[i, j]; t_ += 1
        Oc = (sp.eye(4) - Kq) * (sp.eye(4) + Kq).inv()
        dO = sp.Matrix([[Oc[1, j] * Oc[2, j] for j in range(3)], [Oc[1, j] * Oc[3, j] for j in range(3)], [Oc[2, j] * Oc[3, j] for j in range(3)]]).det()
        if all(x != 0 for x in Oc) and dO != 0:
            break
    assert Oc * Oc.T == sp.eye(4) and Oc.det() == 1
    say(f"      Cayley test matrix: all 16 entries nonzero, d(O) = {dO}")
    exact_C8(Oc, f"Cayley O = (I-K)(I+K)^-1, 3K = {(Kq * 3).tolist()}")
    O9 = sp.Matrix([[0, sp.Rational(2, 3), sp.Rational(2, 3), -sp.Rational(1, 3)],
                    [sp.Rational(2, 3), sp.Rational(5, 9), -sp.Rational(4, 9), sp.Rational(2, 9)],
                    [sp.Rational(2, 3), -sp.Rational(4, 9), sp.Rational(5, 9), sp.Rational(2, 9)],
                    [sp.Rational(1, 3), -sp.Rational(2, 9), -sp.Rational(2, 9), -sp.Rational(8, 9)]])
    say(f"      (Lemma C.9 example matrix: O O^t = I: {O9 * O9.T == sp.eye(4)}, det O = {O9.det()}; computed here only as a cross-check for part 2)")
    exact_C8(O9, "Lemma C.9 example O (zero entry O_00)")


# ----------------------------------------------------------------------------
# F. End-to-end consequence of C.4-C.11: at every critical point O of
#    f(O) = Re(theta det(D_a + O D_b O^t)) on SO(4) with [I,O] stable and d(O) != 0,
#    the proof forces h = f(O) = 0 and det(D_a + O D_b O^t) = 0.
# ----------------------------------------------------------------------------
def section_F():
    hdr("F. End-to-end test of the critical-point mechanism of the proof of Thm C.1 (real a, complex b)")
    from scipy.linalg import expm
    from scipy.optimize import least_squares, minimize
    rng = np.random.default_rng(1001)
    def adj4(M):
        A = np.zeros((4, 4), complex)
        for i in range(4):
            for j in range(4):
                r = [k for k in range(4) if k != i]; c = [k for k in range(4) if k != j]
                A[j, i] = (-1) ** (i + j) * np.linalg.det(M[np.ix_(r, c)])
        return A
    cats = dict(decomposable=0, cauchy=0, h_zero=0, VIOLATION=0)
    viol = []
    nmax = dict(decomposable=0, cauchy=0, h_zero=0, VIOLATION=0)
    t0 = time.time()
    ninst = 0
    while time.time() - t0 < 150 and ninst < 40:
        ninst += 1
        a = rng.normal(size=4); b = rng.normal(size=4) + 1j * rng.normal(size=4); th = np.exp(1j * rng.uniform(0, 2 * np.pi))
        Db = np.diag(b); Da = np.diag(a)
        zsig = th * sigma_points(a, b); scale = np.abs(zsig).max()
        def f(O):
            return (th * np.linalg.det(Da + O @ Db @ O.T)).real
        def grad(O):
            Bp = O @ Db @ O.T; M = Da + Bp; Ad = adj4(M)
            return np.array([(th * np.trace(Ad @ (E @ Bp - Bp @ E))).real for E in SO_B])
        def Ofrom(p, O0):
            K = sum(p[k] * SO_B[k] for k in range(6)); return expm(K) @ O0
        found = []
        for st in range(12):
            O0 = rand_SO4(rng)
            if st % 2 == 0:  # local maxima
                r = minimize(lambda p: -f(Ofrom(p, O0)) / scale, np.zeros(6), method="BFGS", options=dict(gtol=1e-11))
                O1 = Ofrom(r.x, O0)
            else:
                O1 = O0
            r2 = least_squares(lambda p: grad(Ofrom(p, O1)) / scale, np.zeros(6), xtol=1e-15, ftol=1e-15, gtol=1e-15)
            O2 = Ofrom(r2.x, O1)
            if np.abs(grad(O2)).max() / scale > 1e-10:
                continue
            found.append((O2, st % 2 == 0))
        for O2, ismax in found:
            absO = np.abs(O2)
            # connectivity of bipartite support graph
            adj = absO > 1e-6
            seen_r, seen_c, frontier = {0}, set(), [("r", 0)]
            while frontier:
                t, i = frontier.pop()
                if t == "r":
                    for j in range(4):
                        if adj[i, j] and j not in seen_c: seen_c.add(j); frontier.append(("c", j))
                else:
                    for k in range(4):
                        if adj[k, i] and k not in seen_r: seen_r.add(k); frontier.append(("r", k))
            indec = len(seen_r) == 4 and len(seen_c) == 4
            dval = abs(d_cauchy(O2))
            h = f(O2) / scale; dv = abs(np.linalg.det(Da + O2 @ Db @ O2.T)) * abs(th) / scale
            if not indec:
                key = "decomposable"
            elif dval < 1e-6:
                key = "cauchy"
            elif abs(h) < 1e-8 and dv < 1e-8:
                key = "h_zero"
            else:
                key = "VIOLATION"; viol.append((ninst, h, dv, dval, absO.min()))
            cats[key] += 1
            if ismax:
                nmax[key] += 1
    say(f"[F NUMERIC seed 1001] {ninst} random instances (a real, b complex, |theta|=1), 12 starts each "
        f"(BFGS local maxima + Gauss-Newton critical points, |grad|/scale < 1e-10); tolerances: support 1e-6, |d(O)| 1e-6, |h|,|det| 1e-8 relative")
    say(f"      critical points by type: {cats}; of which from maximization: {nmax}")
    if viol:
        say(f"      VIOLATIONS (instance, h, |theta det|, |d|, min|O_ij|): {viol[:8]}")


if __name__ == "__main__":
    secs = sys.argv[1:] or list("ABCDEF")
    for s in secs:
        {"A": section_A, "B": section_B, "C": section_C, "D": section_D, "E": section_E, "F": section_F}[s]()
        say(f"[elapsed {time.time() - T0:.1f}s]")
