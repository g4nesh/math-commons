"""Shared exact objects and tools for the independent replay of the group V cases
T_44 = {1,4,5,6,7,8} and T_46 = {3,4,5,9,10,12} of D. Zheng, "The Marcus-de Oliveira
conjecture for n = 4" (October 2026).

PARTS 1 and 2 are COPIED from this owner's earlier helper `../cases-v-a/replay_T29_common.py`,
unchanged except for the table PAPER_T_V (the cases needed here) and one added entry "SOC" in
the table UNIVERSAL of representatives:

  PART 1 (foundations; there copied from `../cases-i-iv/replay_iv_common.py`): Checker,
  permutations, the group of order 1152, exact linear algebra, bracket polynomials, the
  matrix P, the cubic C of (A.2) and its Hessian, Gaussian rationals, matching vectors, the 82
  vectors f_j of Appendix D, extreme rays of a contact cone, check_foundations,
  check_radial_necessity, check_minor_leibniz.  Those objects were rebuilt from the
  DEFINITIONS in the paper; the copied checks are re-run by every script that uses them, so
  nothing is taken on trust from the earlier replays.

  PART 2: sparse integer polynomials with packed exponent keys, exact composition with
  rational charts, an exact determinant for matrices of linear forms, exact division and
  square roots of polynomials, radial_data (C0, Q, det Hess C on the lift, N = Q^4 det and
  the factor G with N = -4 Q^2 G^2), Sturm certificates, transport of the universal
  inequalities under the group, and floating-point helpers used for falsification attempts
  only.

PART 3 was written for T_44 and T_46: sparse integer polynomials on python-flint
(fmpz_mpoly; only +, -, *, **, exact division, gcd and term access are used) for the large
chart tables, reduced fractions of such polynomials, substitution of a rational chart with
the denominators cleared, and the digest of a polynomial.  The classes PR and RF follow the
design of this owner's `../cases-v-a/replay_T35_common_poly.py` (re-typed here, with the
conversion to and from the packed-key polynomials of Part 2 added).

PART 4 (used by T_46 only) is COPIED from this owner's `../cases-v-a/replay_T25_common.py` and
adapted to the class PR: the "coefficient Gram" sum-of-squares representation on SO(4)
(quaternion model, class GramFace, symmetric elimination, sym_gram_check).

Nothing is imported from, executed from, or copied out of the author's supplement.  When the
environment variable MOC4_SUPP is set, a few JSON files are read as DATA: most are only
COMPARED with the objects rebuilt here; the two files read by GramFace (Part 4) supply a
change of basis and sixteen cubic polynomials, and the identity that makes them a certificate
is re-derived and verified here against Leibniz sums rebuilt from the paper.  The audit
library `zheng-n4-audit/lib/moc4lib.py` is not used to build any object; it is only imported,
when present, to compare finished objects.

Exact arithmetic: Python ints and fractions.Fraction; python-flint integers and integer
polynomials; sympy for a few symbolic identities; numpy int64 inside det_linear_forms only,
under an a-priori bound that excludes overflow.  numpy/scipy floating point occurs only in
functions whose names start with `numeric_`.

Conventions (as in the earlier replay)
  * Permutations of {0,1,2,3} are tuples s with s[i] the image of i; PERMS is in
    lexicographic order, so PERMS[j] is the paper's pi_j.
  * Labels 1..8 carry vectors v_l = (x_l, y_l); [ij] = x_i y_j - x_j y_i.  With
    v_i = (a_i, 1), v_{4+j} = (-b_j, 1) one has [i, 4+j] = a_i + b_j.
  * V_R = Q^14 with coordinates m_0..m_13; values(c) = P c is the vector of the 24 numbers
    (P c)_sigma in PERMS order.  P is injective, so vectors of V_R are handled through values.
"""

# ========================================================================================
# PART 1.  Copied from ../cases-i-iv/replay_iv_common.py (same owner).
# ========================================================================================
from fractions import Fraction
import hashlib
import itertools
import os
import sys
import time

# ----------------------------------------------------------------------------------------
# check registry
# ----------------------------------------------------------------------------------------


class Checker:
    """Prints one line per check and keeps the count of failures."""

    def __init__(self, title):
        self.title = title
        self.n = 0
        self.failed = []
        self.t0 = time.time()
        print("=" * 100)
        print(title)
        print("=" * 100)

    def section(self, name):
        print()
        print("--- " + name)

    def check(self, level, name, ok, detail=""):
        """level in {EXACT, NUMERIC}; ok must be a bool."""
        self.n += 1
        tag = "PASS" if ok else "FAIL"
        line = "[%s] [%s] %s" % (tag, level, name)
        if detail:
            line += "  ::  " + str(detail)
        print(line)
        if not ok:
            self.failed.append(name)
        sys.stdout.flush()
        return ok

    def info(self, text):
        print("[INFO] " + str(text))
        sys.stdout.flush()

    def finish(self):
        print()
        dt = time.time() - self.t0
        if self.failed:
            print("FAILED CHECKS:")
            for f in self.failed:
                print("   " + f)
            print("FINAL: FAIL  (%d of %d checks failed, %.1f s)" % (len(self.failed), self.n, dt))
            sys.exit(1)
        print("FINAL: PASS  (%d checks, %.1f s)" % (self.n, dt))
        sys.exit(0)


# ----------------------------------------------------------------------------------------
# permutations and the group G of order 1152
# ----------------------------------------------------------------------------------------
PERMS = [tuple(p) for p in itertools.permutations(range(4))]
PIDX = {p: i for i, p in enumerate(PERMS)}


def sgn(p):
    inv = sum(1 for i in range(len(p)) for j in range(i + 1, len(p)) if p[i] > p[j])
    return -1 if inv % 2 else 1


def compose(p, q):
    """(p o q)(i) = p(q(i))."""
    return tuple(p[q[i]] for i in range(len(q)))


def inverse(p):
    r = [0] * len(p)
    for i, x in enumerate(p):
        r[x] = i
    return tuple(r)


def oneline(s):
    """'4213' (the paper's 1-based one-line notation) -> (3, 1, 0, 2)."""
    return tuple(int(ch) - 1 for ch in s)


def build_group():
    """All maps sigma -> rho sigma rho' and sigma -> rho sigma^{-1} rho' as index maps
    (tuples t with t[j] = index of the image of pi_j).  Returns the list of distinct maps."""
    maps = set()
    for e in (0, 1):
        for r1 in PERMS:
            for r2 in PERMS:
                t = []
                for s in PERMS:
                    u = inverse(s) if e else s
                    t.append(PIDX[compose(r1, compose(u, r2))])
                maps.add(tuple(t))
    return sorted(maps)


G_MAPS = build_group()


def orbit_rep(T):
    """Lexicographically smallest sorted image of the index set T under G."""
    T = list(T)
    return min(tuple(sorted(t[j] for j in T)) for t in G_MAPS)


# ----------------------------------------------------------------------------------------
# exact linear algebra over Q
# ----------------------------------------------------------------------------------------


def rref(rows):
    A = [[Fraction(x) for x in r] for r in rows]
    if not A:
        return A, []
    nr, nc = len(A), len(A[0])
    piv, r = [], 0
    for c in range(nc):
        k = next((i for i in range(r, nr) if A[i][c] != 0), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        iv = 1 / A[r][c]
        A[r] = [x * iv for x in A[r]]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == nr:
            break
    return A, piv


def rank(rows):
    return len(rref(rows)[1])


def nullspace(rows, ncols=None):
    if not rows:
        return [[Fraction(int(i == j)) for j in range(ncols)] for i in range(ncols)]
    R, piv = rref(rows)
    nc = len(rows[0])
    out = []
    for f in [c for c in range(nc) if c not in piv]:
        x = [Fraction(0)] * nc
        x[f] = Fraction(1)
        for r, p in enumerate(piv):
            x[p] = -R[r][f]
        out.append(x)
    return out


def solve_square(A, b):
    """Unique solution of the square system A x = b over Q (raises if singular)."""
    n = len(A)
    R, piv = rref([list(A[i]) + [b[i]] for i in range(n)])
    if piv[:n] != list(range(n)):
        raise ValueError("singular system")
    return [R[i][n] for i in range(n)]


def int_null_vector(rows, ncols):
    """For an integer matrix with ncols columns: if its rank is ncols-1 return a primitive
    integer vector spanning the kernel, else None.  Fraction-free elimination."""
    A = [list(r) for r in rows]
    nr = len(A)
    piv_cols = []
    r = 0
    for c in range(ncols):
        k = next((i for i in range(r, nr) if A[i][c] != 0), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        pr = A[r]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f, g = A[i][c], pr[c]
                A[i] = [g * a - f * b for a, b in zip(A[i], pr)]
        piv_cols.append(c)
        r += 1
        if r == nr:
            break
    if len(piv_cols) != ncols - 1:
        return None
    free = next(c for c in range(ncols) if c not in piv_cols)
    # x_free = D, x_piv = -A[r][free] * D / A[r][piv]
    from math import gcd
    D = 1
    for i, c in enumerate(piv_cols):
        D = D * abs(A[i][c]) // gcd(D, abs(A[i][c]))
    x = [0] * ncols
    x[free] = D
    for i, c in enumerate(piv_cols):
        x[c] = -A[i][free] * D // A[i][c]
        assert x[c] * A[i][c] == -A[i][free] * D
    g = 0
    for v in x:
        g = gcd(g, abs(v))
    return [v // g for v in x]


# ----------------------------------------------------------------------------------------
# bracket polynomials in sixteen variables, the basis (A.1) and the matrix P
# ----------------------------------------------------------------------------------------


def _xk(l):
    return 8 ** (l - 1)


def _yk(l):
    return 8 ** (7 + l)


def padd(p, q, c=1):
    out = dict(p)
    for k, v in q.items():
        w = out.get(k, 0) + c * v
        if w:
            out[k] = w
        else:
            out.pop(k, None)
    return out


def pmul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            k = a + b
            w = out.get(k, 0) + ca * cb
            if w:
                out[k] = w
            else:
                out.pop(k, None)
    return out


def pscale(p, c):
    return {k: c * v for k, v in p.items()} if c else {}


def bracket(i, j):
    """[ij] = x_i y_j - x_j y_i for labels 1..8."""
    if i == j:
        return {}
    return {_xk(i) + _yk(j): 1, _xk(j) + _yk(i): -1}


def bracket_product(pairs):
    p = {0: 1}
    for i, j in pairs:
        p = pmul(p, bracket(i, j))
    return p


# (A.1), read from the page image of p. 35 of the paper.
MATCHINGS = [
    [(1, 2), (3, 4), (5, 6), (7, 8)],   # B0
    [(1, 2), (3, 4), (5, 8), (6, 7)],   # B1
    [(1, 2), (3, 6), (4, 5), (7, 8)],   # B2
    [(1, 2), (3, 8), (4, 5), (6, 7)],   # B3
    [(1, 2), (3, 8), (4, 7), (5, 6)],   # B4
    [(1, 4), (2, 3), (5, 6), (7, 8)],   # B5
    [(1, 4), (2, 3), (5, 8), (6, 7)],   # B6
    [(1, 6), (2, 3), (4, 5), (7, 8)],   # B7
    [(1, 6), (2, 5), (3, 4), (7, 8)],   # B8
    [(1, 8), (2, 3), (4, 5), (6, 7)],   # B9
    [(1, 8), (2, 3), (4, 7), (5, 6)],   # B10
    [(1, 8), (2, 5), (3, 4), (6, 7)],   # B11
    [(1, 8), (2, 7), (3, 4), (5, 6)],   # B12
    [(1, 8), (2, 7), (3, 6), (4, 5)],   # B13
]
B_POLYS = [bracket_product(m) for m in MATCHINGS]


def p_sigma(s):
    """prod_i [i, 4 + sigma(i)] as a polynomial in sixteen variables."""
    return bracket_product([(i + 1, 5 + s[i]) for i in range(4)])


P_SIGMA = [p_sigma(s) for s in PERMS]


def _build_P():
    monos = sorted({k for b in B_POLYS for k in b})
    rows = [[b.get(k, 0) for k in monos] for b in B_POLYS]          # 14 x #monomials
    _, piv = rref(rows)
    if len(piv) != 14:
        raise AssertionError("B_0..B_13 are not linearly independent")
    sel = [monos[c] for c in piv]
    # coefficients c with sum_a c_a B_a[k] = target[k] on the selected monomials
    At = [[B_POLYS[a].get(k, 0) for a in range(14)] for k in sel]
    P = []
    for s in PERMS:
        tgt = p_sigma(s)
        c = solve_square(At, [tgt.get(k, 0) for k in sel])
        recon = {}
        for a in range(14):
            if c[a]:
                if c[a].denominator != 1:
                    raise AssertionError("non-integer entry of P")
                recon = padd(recon, B_POLYS[a], int(c[a]))
        if recon != tgt:
            raise AssertionError("p_sigma is not in the span of the B_a")
        P.append([int(x) for x in c])
    return P


P = _build_P()          # 24 x 14, rows in PERMS order; every row is a verified identity
KAPPA = [Fraction(0)] * 13 + [Fraction(1)]


def values(c):
    return [sum(P[s][a] * c[a] for a in range(14)) for s in range(24)]


def _pivot_rows():
    _, piv = rref([list(col) for col in zip(*P)])   # row-reduce P^t: pivots = independent rows
    return piv


_PIV = _pivot_rows()
assert len(_PIV) == 14


def solve_P(vals):
    """The unique c in Q^14 with P c = vals; raises unless all 24 equations hold."""
    c = solve_square([P[r] for r in _PIV], [Fraction(vals[r]) for r in _PIV])
    if values(c) != [Fraction(v) for v in vals]:
        raise ValueError("value vector is not in the image of P")
    return c


# ----------------------------------------------------------------------------------------
# the cubic (A.2), read from the page image of p. 36 of the paper
# ----------------------------------------------------------------------------------------


def cubic_C(m):
    (m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13) = m
    return (m13 * (m0 * m6 - m1 * m5)
            - m0 * m9 * (m0 + m1 + m2 + m3 + m5 + m7 + m9) - m1 * m5 * m9
            + m10 * (-m0 * m3 - m0 * m9 + m1 * m2 + m1 * m7)
            + m11 * (-m0 * m7 - m0 * m9 + m2 * m5 + m3 * m5)
            + m12 * (-m0 * m9 + m1 * m7 - m2 * m6 - m2 * m9 + m3 * m5 + m3 * m7)
            + m2 * m10 * m11 - m3 * m8 * m10 - m4 * m7 * m11 + m4 * m8 * m9)


_CUBIC_CACHE = {}


def cubic_data():
    """sympy symbols m0..m13, C, its gradient (list) and Hessian (list of lists)."""
    if not _CUBIC_CACHE:
        import sympy
        ms = sympy.symbols("m0:14")
        C = sympy.expand(cubic_C(ms))
        grad = [sympy.diff(C, v) for v in ms]
        hess = [[sympy.diff(g, v) for v in ms] for g in grad]
        _CUBIC_CACHE.update(ms=ms, C=C, grad=grad, hess=hess)
    d = _CUBIC_CACHE
    return d["ms"], d["C"], d["grad"], d["hess"]


def grad_terms():
    """For each j the list of (coefficient, i, k) with dC/dm_j = sum coefficient m_i m_k."""
    import sympy
    ms, C, grad, _ = cubic_data()
    out = []
    for g in grad:
        terms = []
        for mono, coef in sympy.Poly(g, *ms).terms():
            idx = [i for i, e in enumerate(mono) for _ in range(e)]
            assert len(idx) == 2 and coef.is_Integer
            terms.append((int(coef), idx[0], idx[1]))
        out.append(terms)
    return out


def hess_at(u):
    """Hess C(u) as a 14 x 14 list of exact numbers (entries are linear in u)."""
    import sympy
    ms, _, _, hess = cubic_data()
    if "hess_lin" not in _CUBIC_CACHE:
        table = []
        for i in range(14):
            row = []
            for j in range(14):
                lin = []
                if hess[i][j] != 0:
                    for mono, coef in sympy.Poly(hess[i][j], *ms).terms():
                        k = [t for t, ee in enumerate(mono) if ee]
                        assert len(k) == 1 and mono[k[0]] == 1 and coef.is_Integer
                        lin.append((int(coef), k[0]))
                row.append(lin)
            table.append(row)
        _CUBIC_CACHE["hess_lin"] = table
    table = _CUBIC_CACHE["hess_lin"]
    return [[sum((cf * u[k] for cf, k in table[i][j]), Fraction(0)) for j in range(14)]
            for i in range(14)]


# ----------------------------------------------------------------------------------------
# Gaussian rationals and matching vectors
# ----------------------------------------------------------------------------------------


class GQ:
    """Exact complex number re + i im with Fraction parts."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = Fraction(re), Fraction(im)

    @staticmethod
    def of(o):
        return o if isinstance(o, GQ) else GQ(o, 0)

    def __add__(self, o):
        o = GQ.of(o)
        return GQ(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __sub__(self, o):
        o = GQ.of(o)
        return GQ(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return GQ.of(o) - self

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __mul__(self, o):
        o = GQ.of(o)
        return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __eq__(self, o):
        o = GQ.of(o)
        return self.re == o.re and self.im == o.im

    def __hash__(self):
        return hash((self.re, self.im))

    def __repr__(self):
        return "GQ(%s, %s)" % (self.re, self.im)


def eval_poly16(poly, pts):
    """Evaluate a sixteen-variable polynomial at pts[l-1] = (x_l, y_l) (any ring)."""
    tot = 0
    for key, c in poly.items():
        t = c
        k = key
        for slot in range(16):
            e = k % 8
            k //= 8
            if e:
                base = pts[slot][0] if slot < 8 else pts[slot - 8][1]
                for _ in range(e):
                    t = t * base
        tot = tot + t
    return tot


def matching_vector(a, b, alpha=1):
    """x = alpha Psi(v) with v_i = (a_i, 1), v_{4+j} = (-b_j, 1); (P x)_sigma = alpha z_sigma."""
    pts = [(a[i], 1) for i in range(4)] + [(-b[j], 1) for j in range(4)]
    return [alpha * eval_poly16(B_POLYS[r], pts) for r in range(14)]


def z_points(a, b):
    out = []
    for s in PERMS:
        t = 1
        for i in range(4):
            t = t * (a[i] + b[s[i]])
        out.append(t)
    return out


# ----------------------------------------------------------------------------------------
# the 82 vectors f_j of Appendix D (through their values on the 24 permutations)
# ----------------------------------------------------------------------------------------


def n_values(I, J):
    """Values of n_{I,J}: 1 at sigma with sigma(I) = J, else 0 (I, J 1-based)."""
    I0 = [i - 1 for i in I]
    J0 = set(j - 1 for j in J)
    return [1 if set(s[i] for i in I0) == J0 else 0 for s in PERMS]


def t_values(pi):
    out = [0] * 24
    for i in range(4):
        out = [x + y for x, y in zip(out, n_values((i + 1,), (pi[i] + 1,)))]
    return out


def d_values(pi):
    p = [x + 1 for x in pi]
    terms = [(+1, (1, 2), (p[0], p[2])), (+1, (1, 2), (p[1], p[3])),
             (-1, (1, 3), (p[0], p[1])), (-1, (1, 3), (p[2], p[3]))]
    out = [0] * 24
    for sg, I, J in terms:
        out = [x + sg * y for x, y in zip(out, n_values(I, J))]
    return out


def h_values(pi, eps):
    v = [Fraction(x + eps * y, 2) for x, y in zip(t_values(pi), d_values(pi))]
    assert all(x.denominator == 1 for x in v)
    return [int(x) for x in v]


# Appendix D, read from the page image of p. 59: ("n", I, J) or ("h", one-line pi, eps).
F_TABLE = [
    ("n", "1", "4"), ("n", "12", "34"), ("n", "13", "34"), ("n", "14", "34"),     # 0-3
    ("n", "1", "3"), ("n", "12", "24"), ("n", "4", "1"), ("n", "3", "1"),         # 4-7
    ("n", "12", "23"), ("n", "14", "23"), ("n", "13", "24"), ("n", "13", "23"),   # 8-11
    ("n", "14", "24"), ("n", "2", "1"), ("n", "1", "2"), ("n", "12", "14"),       # 12-15
    ("n", "2", "4"), ("h", "3421", +1), ("h", "3412", -1), ("h", "2413", -1),     # 16-19
    ("h", "2431", +1), ("n", "4", "2"), ("h", "3412", +1), ("n", "14", "12"),     # 20-23
    ("h", "4312", +1), ("h", "3142", -1), ("h", "4132", -1), ("n", "3", "2"),     # 24-27
    ("h", "3421", -1), ("n", "13", "12"), ("h", "4321", -1), ("h", "4123", +1),   # 28-31
    ("h", "3124", +1), ("n", "12", "13"), ("h", "4321", +1), ("h", "4312", -1),   # 32-35
    ("n", "2", "3"), ("h", "2341", +1), ("h", "2314", -1), ("n", "14", "13"),     # 36-39
    ("h", "4213", +1),                                                            # 40
    ("n", "4", "3"), ("h", "4123", -1), ("h", "2413", +1), ("h", "2143", -1),     # 41-44
    ("n", "13", "14"), ("h", "3241", -1), ("n", "3", "4"), ("h", "3142", +1),     # 45-48
    ("h", "2143", +1), ("h", "2341", -1), ("h", "1423", -1), ("h", "1342", +1),   # 49-52
    ("n", "13", "13"), ("h", "4231", -1), ("n", "3", "3"), ("h", "4132", +1),     # 53-56
    ("h", "2431", -1), ("h", "2134", +1), ("n", "14", "14"), ("h", "3214", +1),   # 57-60
    ("n", "4", "4"), ("h", "3124", -1), ("h", "2134", -1), ("h", "2314", +1),     # 61-64
    ("h", "1432", +1), ("h", "1324", -1), ("n", "2", "2"), ("h", "4213", -1),     # 65-68
    ("h", "3241", +1), ("h", "4231", +1), ("h", "3214", -1), ("n", "12", "12"),   # 69-72
    ("h", "1243", +1), ("h", "1432", -1), ("h", "1423", +1), ("h", "1234", -1),   # 73-76
    ("h", "1342", -1), ("h", "1234", +1), ("h", "1243", -1), ("h", "1324", +1),   # 77-80
    ("n", "1", "1"),                                                              # 81
]
assert len(F_TABLE) == 82


def _f_values(spec):
    kind, a, b = spec
    if kind == "n":
        return n_values(tuple(int(ch) for ch in a), tuple(int(ch) for ch in b))
    return h_values(oneline(a), b)


F_VALUES = [tuple(_f_values(s)) for s in F_TABLE]
F_IS_MINOR = [s[0] == "n" for s in F_TABLE]
F_COORD = [solve_P(v) for v in F_VALUES]       # raises unless every f_j is in the image of P


def f_name(j):
    k, a, b = F_TABLE[j]
    return "n_{%s,%s}" % (a, b) if k == "n" else "h_{%s,%s}" % (a, "+" if b > 0 else "-")


def generators(T):
    """Indices j with f_j = 0 on every permutation index of T."""
    return [j for j in range(82) if all(F_VALUES[j][s] == 0 for s in T)]


def combo_values(coeffs):
    """Values of sum_j coeffs[j] f_j, for a dict {j: coefficient} (or list of pairs)."""
    items = coeffs.items() if isinstance(coeffs, dict) else coeffs
    out = [0] * 24
    for j, c in items:
        out = [x + c * y for x, y in zip(out, F_VALUES[j])]
    return out


def combo_coord(coeffs):
    """Coordinates in V_R of sum_j coeffs[j] f_j."""
    items = coeffs.items() if isinstance(coeffs, dict) else coeffs
    out = [0] * 14
    for j, c in items:
        out = [x + c * y for x, y in zip(out, F_COORD[j])]
    return out


# ----------------------------------------------------------------------------------------
# contact cone: direct enumeration of the extreme rays of N_T
# ----------------------------------------------------------------------------------------


def extreme_rays_of_contact_cone(T):
    """Extreme rays of N_T = {c : P c >= 0, (P c)_T = 0}, computed from the definition and
    without the list of 82 facets.  Returns (dim of the linear span of the equalities,
    sorted list of primitive integer value vectors (24-tuples) of the extreme rays).

    Method: N_T is a pointed polyhedral cone in the subspace L = {(Pc)_T = 0} of dimension
    d.  A nonzero c in N_T spans an extreme ray iff the inequalities tight at c have rank
    d - 1 on L.  So every extreme ray is the kernel of some d - 1 of the 24 - |T| remaining
    inequality rows; all such subsets are enumerated."""
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    N = nullspace([P[s] for s in T], 14)                     # basis of L, as rows
    d = len(N)
    # integer basis
    from math import gcd
    Nint = []
    for v in N:
        den = 1
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
        Nint.append([int(x * den) for x in v])
    # A[s][k] = (P N_k)_s for s not in T: the inequality rows in coordinates y on L
    A = {s: [sum(P[s][a] * Nint[k][a] for a in range(14)) for k in range(d)] for s in others}
    assert rank([A[s] for s in others]) == d, "cone is not pointed"
    rays = set()
    for sub in itertools.combinations(others, d - 1):
        y = int_null_vector([A[s] for s in sub], d)
        if y is None:
            continue
        w = [sum(A[s][k] * y[k] for k in range(d)) for s in others]
        if all(x >= 0 for x in w):
            pass
        elif all(x <= 0 for x in w):
            w = [-x for x in w]
        else:
            continue
        full = [0] * 24
        for s, x in zip(others, w):
            full[s] = x
        g = 0
        for x in full:
            g = gcd(g, x)
        rays.add(tuple(x // g for x in full))
    return d, sorted(rays)


def check_foundations(ck):
    """The objects of Sections 3 and 4 of the paper that the three cases rely on."""
    ck.section("Foundations rebuilt from the paper's definitions")
    ck.check("EXACT", "P is a 24 x 14 integer matrix; every row is a polynomial identity in "
             "16 variables (checked at construction); rank 14",
             len(P) == 24 and all(len(r) == 14 for r in P) and rank(P) == 14)
    ck.check("EXACT", "P kappa = (1,...,1) for kappa = e_13",
             all(P[s][13] == 1 for s in range(24)))
    pts = [(0, 1)] * 4 + [(1, 1)] * 4
    ck.check("EXACT", "kappa = Psi(v) at v_1..4 = (0,1), v_5..8 = (1,1), so kappa is in X",
             [eval_poly16(B_POLYS[r], pts) for r in range(14)] == [0] * 13 + [1])
    ck.check("EXACT", "paper's pivot rows 0,1,2,3,4,6,7,8,9,10,12,13,16,18 of P are independent",
             rank([P[r] for r in (0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18)]) == 14)
    # gradient of C vanishes on the image of Psi (sixteen-variable identities)
    gt = grad_terms()
    bad = []
    BB = {}
    for j in range(14):
        tot = {}
        for coef, i, k in gt[j]:
            key = (min(i, k), max(i, k))
            if key not in BB:
                BB[key] = pmul(B_POLYS[i], B_POLYS[k])
            tot = padd(tot, BB[key], coef)
        if tot:
            bad.append(j)
    ck.check("EXACT", "dC/dm_j (Psi(v)) = 0 identically in 16 variables for j = 0..13 "
             "(so X lies in the singular locus of C)", not bad, "nonzero: %s" % bad if bad else "")
    # facet vectors
    prof_ok = True
    for j in range(82):
        v = F_VALUES[j]
        cnt = (v.count(0), v.count(1), v.count(2))
        if F_IS_MINOR[j]:
            prof_ok &= cnt in ((18, 6, 0), (20, 4, 0))
        else:
            prof_ok &= cnt == (13, 10, 1)
    ck.check("EXACT", "82 vectors f_j: all lie in the image of P; 34 minor vectors with 18 or "
             "20 zeros, 48 vectors h with value profile 13/10/1; pairwise distinct",
             prof_ok and len(set(F_VALUES)) == 82 and sum(F_IS_MINOR) == 34)
    ck.check("EXACT", "group G: 1152 distinct maps on the 24 permutation indices",
             len(G_MAPS) == 1152)
    fset = {v: j for j, v in enumerate(F_VALUES)}
    ok = True
    for t in G_MAPS:
        for j in range(82):
            w = [0] * 24
            for s in range(24):
                w[t[s]] = F_VALUES[j][s]
            jj = fset.get(tuple(w))
            if jj is None or F_IS_MINOR[jj] != F_IS_MINOR[j]:
                ok = False
    ck.check("EXACT", "G permutes the 82 vectors f_j and maps minor vectors to minor vectors", ok)


def check_radial_necessity(ck, seed=20261125, samples=12):
    """7.3 of the audit's PROOF.md / the 'real part obeys the cubic' paragraph of the
    documents: x = u + i v in X  =>  C(u) = 0 and Hess C(u) v = 0."""
    import random
    import sympy
    ck.section("Radial necessity: the real part of a point of X lies on the cubic")
    ms, C, grad, hess = cubic_data()
    us = sympy.symbols("u0:14")
    vs = sympy.symbols("v0:14")
    t = sympy.Symbol("t")
    su = dict(zip(ms, us))
    sv = dict(zip(ms, vs))
    suv = {m: u + t * v for m, u, v in zip(ms, us, vs)}
    ok = True
    for j in range(14):
        lhs = sympy.expand(grad[j].xreplace(suv))
        hv = sum(hess[j][k].xreplace(su) * vs[k] for k in range(14))
        rhs = sympy.expand(grad[j].xreplace(su) + t * hv + t ** 2 * grad[j].xreplace(sv))
        ok &= sympy.expand(lhs - rhs) == 0
    ck.check("EXACT", "grad C(u + t v) = grad C(u) + t Hess C(u) v + t^2 grad C(v) identically "
             "(t = i: real part grad C(u) - grad C(v), imaginary part Hess C(u) v)", ok)
    euler = sympy.expand(sum(us[j] * grad[j].xreplace(su) for j in range(14)) - 3 * C.xreplace(su))
    ck.check("EXACT", "Euler: u . grad C(u) = 3 C(u) identically", euler == 0)
    pol = sympy.expand(sum(us[j] * grad[j].xreplace(sv) for j in range(14))
                       - sympy.Rational(1, 2) * sum(vs[j] * hess[j][k].xreplace(su) * vs[k]
                                                    for j in range(14) for k in range(14)))
    ck.check("EXACT", "polarization: u . grad C(v) = (1/2) v^t Hess C(u) v identically", pol == 0)
    ck.info("hence grad C(u+iv) = 0 gives grad C(u) = grad C(v) and Hess C(u) v = 0, and "
            "3 C(u) = u.grad C(u) = u.grad C(v) = (1/2) v^t Hess C(u) v = 0")
    rng = random.Random(seed)
    okc = okh = True
    for _ in range(samples):
        def rq():
            return GQ(Fraction(rng.randint(-9, 9), rng.randint(1, 5)),
                      Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
        a = [rq() for _ in range(4)]
        b = [rq() for _ in range(4)]
        al = rq()
        x = matching_vector(a, b, al)
        x = [GQ.of(c) for c in x]
        u = [c.re for c in x]
        v = [c.im for c in x]
        okc &= cubic_C(u) == 0
        H = hess_at(u)
        okh &= all(sum(H[i][k] * v[k] for k in range(14)) == 0 for i in range(14))
        zs = z_points(a, b)
        okc &= [GQ.of(sum(P[s][r] * x[r] for r in range(14))) for s in range(24)] == \
               [al * z for z in zs]
    ck.check("EXACT", "spot check on %d random Gaussian-rational points x = alpha Psi(a,b) "
             "(seed %d): (P x)_sigma = alpha z_sigma and C(Re x) = 0" % (samples, seed), okc)
    ck.check("EXACT", "same points: Hess C(Re x) Im x = 0", okh)


def check_minor_leibniz(ck, minor_indices):
    """Lemma 6.2 / (4.2) for the listed minor vectors, through two polynomial identities in a
    generic 4 x 4 matrix M:
      (L) sum_{sigma(I)=J} sgn(sigma) prod_i M[i,sigma(i)] = eps det M[I,J] det M[I^c,J^c],
      (J) det(adj M)[J,I] = eps det M[I^c,J^c] det(M)^(|I|-1),      eps = (-1)^(sum I + sum J).
    For O in SO(4), adj O = O^t, so (J) gives eps det O[I^c,J^c] = det O^t[J,I] = det O[I,J]
    and (L) becomes Phi_{n_{I,J}}(O) = det O[I,J]^2 >= 0."""
    import sympy
    M = sympy.Matrix(4, 4, lambda i, j: sympy.Symbol("o%d%d" % (i + 1, j + 1)))
    adj = M.adjugate()
    detM = sympy.expand(M.det(method="berkowitz"))
    ok = True
    for j in minor_indices:
        kind, a, b = F_TABLE[j]
        assert kind == "n"
        I = [int(ch) - 1 for ch in a]
        J = [int(ch) - 1 for ch in b]
        Ic = [i for i in range(4) if i not in I]
        Jc = [i for i in range(4) if i not in J]
        eps = (-1) ** (sum(I) + sum(J))            # same parity for 0- and 1-based labels
        phi = 0
        for s, val in zip(PERMS, F_VALUES[j]):
            if val:
                phi += sgn(s) * M[0, s[0]] * M[1, s[1]] * M[2, s[2]] * M[3, s[3]]
        dIJ = M.extract(I, J).det()
        dc = M.extract(Ic, Jc).det()
        ok &= sympy.expand(phi - eps * dIJ * dc) == 0
        dadj = adj.extract(J, I).det()
        ok &= sympy.expand(dadj - eps * dc * detM ** (len(I) - 1)) == 0
    ck.check("EXACT", "Laplace and Jacobi identities for the %d minor generators %s: "
             "Phi_{f_j}(O) = (minor of O)^2 >= 0 on SO(4)" % (len(minor_indices), list(minor_indices)), ok)



# ========================================================================================
# PART 2.  Copied from ../cases-v-a/replay_T29_common.py (same owner; written there for the
#          group V replays of T_29 and T_31).
# ========================================================================================
import heapq
import random

# Appendix E of the paper (p. 60), read from the text of the PDF: the two cases replayed here
# (44, 46) and the cases their documents cite (7, 15, 25, 34).
PAPER_T_V = {
    7: (0, 1, 2, 4, 7, 8),
    15: (1, 2, 3, 4, 6, 7),
    25: (3, 4, 7, 8, 10, 12),
    34: (2, 3, 5, 9, 10, 12),
    44: (1, 4, 5, 6, 7, 8),
    46: (3, 4, 5, 9, 10, 12),
}

# ----------------------------------------------------------------------------------------
# sparse polynomials with integer coefficients: dict {packed exponent key: int}
# ----------------------------------------------------------------------------------------
SH = 10                       # bits per variable; every polynomial must have total degree < 2**SH
MASK = (1 << SH) - 1
DEGLIM = 1 << SH


class Ring:
    """Variable names for packed-key polynomials.  Variable i occupies bits [SH*i, SH*(i+1))."""

    def __init__(self, names):
        self.names = list(names.split(",")) if isinstance(names, str) else list(names)
        self.n = len(self.names)
        self.pos = {nm: i for i, nm in enumerate(self.names)}
        assert len(self.pos) == self.n

    def var(self, name):
        return {1 << (SH * self.pos[name]): 1}

    def gens(self):
        return [self.var(nm) for nm in self.names]

    def key(self, exps):
        k = 0
        for i, e in enumerate(exps):
            k |= e << (SH * i)
        return k

    def exps(self, key):
        return tuple((key >> (SH * i)) & MASK for i in range(self.n))

    def tdeg(self, p):
        return max((sum(self.exps(k)) for k in p), default=0)

    def deg(self, p, name):
        s = SH * self.pos[name]
        return max(((k >> s) & MASK for k in p), default=0)

    def mul(self, p, q):
        """Checked product (total degree stays below the packing limit)."""
        assert self.tdeg(p) + self.tdeg(q) < DEGLIM
        return p_mul(p, q)

    def pow(self, p, n):
        assert self.tdeg(p) * n < DEGLIM
        out = {0: 1}
        for _ in range(n):
            out = p_mul(out, p)
        return out

    def split(self, p, name):
        """{e: coefficient polynomial of name**e} (the variable removed from the keys)."""
        s = SH * self.pos[name]
        out = {}
        for k, c in p.items():
            e = (k >> s) & MASK
            out.setdefault(e, {})[k - (e << s)] = c
        return out

    def subs_zero(self, p, names):
        """Set the listed variables to 0."""
        m = 0
        for nm in names:
            m |= MASK << (SH * self.pos[nm])
        return {k: c for k, c in p.items() if not (k & m)}

    def subs_one(self, p, name):
        """Set one variable to 1."""
        s = SH * self.pos[name]
        out = {}
        for k, c in p.items():
            kk = k - (((k >> s) & MASK) << s)
            v = out.get(kk, 0) + c
            if v:
                out[kk] = v
            else:
                out.pop(kk, None)
        return out

    def group(self, p, keep):
        """Group by the monomial in the variables NOT in `keep`: {outer key: {inner key: coeff}}
        where the inner polynomial only involves the variables in `keep`."""
        m = 0
        for nm in keep:
            m |= MASK << (SH * self.pos[nm])
        out = {}
        for k, c in p.items():
            out.setdefault(k & ~m, {})[k & m] = c
        return out

    def evaluate(self, p, point):
        """Exact value at a dict {name: number} (Fractions or ints)."""
        vals = [point[nm] for nm in self.names]
        pw = [{0: 1} for _ in range(self.n)]
        tot = 0
        for k, c in p.items():
            t = c
            for i in range(self.n):
                e = (k >> (SH * i)) & MASK
                if e:
                    d = pw[i]
                    if e not in d:
                        d[e] = vals[i] ** e
                    t = t * d[e]
            tot = tot + t
        return tot

    def to_sympy(self, p):
        import sympy
        syms = sympy.symbols(",".join(self.names))
        if self.n == 1:
            syms = (syms,)
        return sympy.Add(*[c * sympy.Mul(*[s ** e for s, e in zip(syms, self.exps(k)) if e])
                           for k, c in p.items()])

    def from_sympy(self, expr):
        import sympy
        syms = sympy.symbols(",".join(self.names))
        if self.n == 1:
            syms = (syms,)
        P_ = sympy.Poly(sympy.expand(expr), *syms)
        out = {}
        for mono, c in P_.terms():
            assert c.is_Integer
            out[self.key(mono)] = int(c)
        return out

    def fmt(self, p, limit=12):
        items = sorted(p.items(), key=lambda kc: tuple(-e for e in self.exps(kc[0])))
        parts = []
        for k, c in items[:limit]:
            mon = "*".join(nm if e == 1 else "%s^%d" % (nm, e) for nm, e in zip(self.names, self.exps(k)) if e)
            parts.append(("%+d" % c) + ("*" + mon if mon else ""))
        if len(items) > limit:
            parts.append("... (%d terms)" % len(items))
        return " ".join(parts) if parts else "0"


def p_add(p, q, c=1):
    out = dict(p)
    for k, v in q.items():
        w = out.get(k, 0) + c * v
        if w:
            out[k] = w
        else:
            out.pop(k, None)
    return out


def p_scale(p, c):
    return {k: c * v for k, v in p.items()} if c else {}


_GUARD = 0
for _i in range(24):
    _GUARD |= 1 << (SH * _i + SH - 1)          # the top bit of every exponent field


def p_mul(p, q):
    """Product.  Rigorous packing guard: every exponent of both factors must be below
    2**(SH-1) (no top bit set in any field), so every exponent of the product is below 2**SH
    and no carry between fields can occur."""
    oa = ob = 0
    for k in p:
        oa |= k
    for k in q:
        ob |= k
    assert not ((oa | ob) & _GUARD) and (oa | ob) < (1 << (SH * 24)), "exponent packing limit"
    if len(p) < len(q):
        p, q = q, p
    out = {}
    get = out.get
    for kb, cb in q.items():
        for ka, ca in p.items():
            k = ka + kb
            out[k] = get(k, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def p_min(p):
    return min(p.values(), default=0)


def p_nonneg(p):
    """All coefficients >= 0."""
    return all(v >= 0 for v in p.values())


def p_const(p):
    return p.get(0, 0)


def p_rename(p, Rsrc, Rdst, mapping=None):
    """Move a polynomial to another ring; variable nm of Rsrc becomes mapping.get(nm, nm) of Rdst."""
    mapping = mapping or {}
    tgt = [Rdst.pos[mapping.get(nm, nm)] if mapping.get(nm, nm) in Rdst.pos else None for nm in Rsrc.names]
    out = {}
    for k, c in p.items():
        kk = 0
        for i in range(Rsrc.n):
            e = (k >> (SH * i)) & MASK
            if e:
                assert tgt[i] is not None, "variable %s not in the target ring" % Rsrc.names[i]
                kk += e << (SH * tgt[i])
        out[kk] = out.get(kk, 0) + c
    return {k: v for k, v in out.items() if v}


def p_compose(F, Rsrc, Rdst, subs, degs=None, rename=None):
    """Substitute x = Num_x / Den_x (polynomials over Rdst) for the variables x in `subs`
    (dict name -> (Num, Den)); every other variable of Rsrc is kept (renamed by `rename`).
    Returns (numerator over Rdst, degs) with
        F(subs) = numerator / prod_x Den_x ** degs[x],
    degs[x] = degree of F in x unless a larger value is supplied (common denominators).
    The total degree of every intermediate product is bounded a priori below the packing limit."""
    rename = rename or {}
    names = [nm for nm in Rsrc.names if nm in subs]
    actual = {nm: Rsrc.deg(F, nm) for nm in names}
    degs = dict(degs or {})
    for nm in names:
        degs[nm] = max(degs.get(nm, 0), actual[nm])
        assert degs[nm] >= actual[nm]
    bound = sum(degs[nm] * max(Rdst.tdeg(subs[nm][0]), Rdst.tdeg(subs[nm][1])) for nm in names)
    bound += Rsrc.tdeg(F)
    assert bound < DEGLIM // 2, "packing limit"
    spos = [Rsrc.pos[nm] for nm in names]
    smask = 0
    for i in spos:
        smask |= MASK << (SH * i)
    # group by the exponents of the substituted variables
    groups = {}
    for k, c in F.items():
        groups.setdefault(k & smask, {})[k & ~smask] = c
    # powers
    powN = {nm: [{0: 1}] for nm in names}
    powD = {nm: [{0: 1}] for nm in names}
    for nm in names:
        for _ in range(degs[nm]):
            powN[nm].append(p_mul(powN[nm][-1], subs[nm][0]))
            powD[nm].append(p_mul(powD[nm][-1], subs[nm][1]))
    fac = {}

    def factor(nm, e):
        key = (nm, e)
        if key not in fac:
            fac[key] = p_mul(powN[nm][e], powD[nm][degs[nm] - e])
        return fac[key]
    prefix = {(): {0: 1}}

    def prod_for(es):
        if es not in prefix:
            head = prod_for(es[:-1])
            prefix[es] = p_mul(head, factor(names[len(es) - 1], es[-1]))
        return prefix[es]
    out = {}
    get = out.get
    for gk in sorted(groups):
        es = tuple((gk >> (SH * i)) & MASK for i in spos)
        X = prod_for(es)
        free = p_rename(groups[gk], Rsrc, Rdst, rename)
        for kb, cb in free.items():
            for ka, ca in X.items():
                k = ka + kb
                out[k] = get(k, 0) + ca * cb
    return {k: v for k, v in out.items() if v}, degs


# ----------------------------------------------------------------------------------------
# exact division and square root of sparse polynomials (lexicographic leading terms)
# ----------------------------------------------------------------------------------------


def p_divexact(p, q):
    """Quotient p / q if q divides p exactly, else None.  Keys are compared as integers, which
    is a monomial order (lex with the last variable highest), so the leading key of a product
    is the sum of the leading keys."""
    lq = max(q)
    cq = q[lq]
    rest = [(k, c) for k, c in q.items() if k != lq]
    rem = dict(p)
    heap = [-k for k in rem]
    heapq.heapify(heap)
    quo = {}
    while heap:
        k = -heapq.heappop(heap)
        c = rem.get(k, 0)
        if c == 0:
            continue
        del rem[k]
        kk = k - lq
        if kk < 0 or c % cq:
            return None
        cc = c // cq
        # the monomial k must be divisible by the leading monomial of q
        ok = True
        t, u = k, lq
        while u:
            if (t & MASK) < (u & MASK):
                ok = False
                break
            t >>= SH
            u >>= SH
        if not ok:
            return None
        quo[kk] = cc
        for k2, c2 in rest:
            k3 = kk + k2
            if k3 in rem:
                v = rem[k3] - cc * c2
                if v:
                    rem[k3] = v
                else:
                    rem[k3] = 0
            else:
                rem[k3] = -cc * c2
                heapq.heappush(heap, -k3)
    return quo


def p_sqrt(p):
    """A polynomial g with g*g == p (leading coefficient > 0), or None.  The result is verified
    by the final exact comparison."""
    if not p:
        return {}
    import math
    rem = dict(p)
    heap = [-k for k in rem]
    heapq.heapify(heap)
    lk = max(rem)
    lc = rem[lk]
    if lc <= 0 or lk % 2:
        return None
    r = math.isqrt(lc)
    if r * r != lc:
        return None
    # check every exponent of the leading monomial is even
    t = lk
    while t:
        if (t & MASK) % 2:
            return None
        t >>= SH
    g0k, g0c = lk >> 1, r
    g = {g0k: g0c}
    order = [(g0k, g0c)]
    # remove g0^2
    rem[lk] = 0
    while heap:
        k = -heapq.heappop(heap)
        c = rem.get(k, 0)
        if c == 0:
            continue
        kk = k - g0k
        if kk < 0 or kk >= g0k or c % (2 * g0c):
            return None
        cc = c // (2 * g0c)
        # rem -= 2*cc*x^kk*g_old + cc^2 x^(2kk)
        for k2, c2 in order:
            k3 = kk + k2
            v = rem.get(k3)
            if v is None:
                rem[k3] = -2 * cc * c2
                heapq.heappush(heap, -k3)
            else:
                rem[k3] = v - 2 * cc * c2
        k3 = 2 * kk
        v = rem.get(k3)
        if v is None:
            rem[k3] = -cc * cc
            heapq.heappush(heap, -k3)
        else:
            rem[k3] = v - cc * cc
        g[kk] = cc
        order.append((kk, cc))
        if len(order) > 200000:
            return None
    if any(v for v in rem.values()):
        return None
    return g


# ----------------------------------------------------------------------------------------
# determinant of a matrix of integer linear forms (dense int64 vectors, rigorous bound)
# ----------------------------------------------------------------------------------------
_MONO = {}
_SHIFT = {}


def _monomials(nv, k):
    """All multisets of size k from range(nv) as sorted tuples (itertools order) and their index."""
    key = (nv, k)
    if key not in _MONO:
        lst = list(itertools.combinations_with_replacement(range(nv), k))
        _MONO[key] = (lst, {m: i for i, m in enumerate(lst)})
    return _MONO[key]


def _shift_index(nv, k, i):
    """numpy index array: position among the degree k+1 monomials of (degree k monomial) * x_i."""
    import bisect
    import numpy as np
    key = (nv, k, i)
    if key not in _SHIFT:
        lst, _ = _monomials(nv, k)
        _, idx1 = _monomials(nv, k + 1)
        out = np.empty(len(lst), dtype=np.int64)
        for n, m in enumerate(lst):
            p = bisect.bisect_left(m, i)
            out[n] = idx1[m[:p] + (i,) + m[p:]]
        _SHIFT[key] = out
    return _SHIFT[key]


def det_linear_forms(M, nv):
    """det of an n x n matrix whose entries are integer linear forms in nv variables
    (M[r][c] = dict {variable index: int}).  Laplace expansion along the rows, one row at a
    time, over all column subsets; a polynomial of degree k is a dense int64 vector indexed by
    the multisets of size k.  numpy is used with INTEGER arithmetic only; the l1 norm of every
    intermediate polynomial is bounded a priori (same recursion on the l1 norms of the entries,
    Python ints) and the function asserts that this bound is below 2^62, so no int64 overflow
    can occur and the result is exact.  Returns ({multiset tuple: int}, l1 bound)."""
    import numpy as np
    n = len(M)
    rows = sorted(range(n), key=lambda r: -sum(1 for c in range(n) if M[r][c]))
    inv = sum(1 for a in range(n) for b in range(a + 1, n) if rows[a] > rows[b])
    rsign = -1 if inv % 2 else 1
    l1 = [[sum(abs(v) for v in M[r][c].values()) for c in range(n)] for r in range(n)]
    # pass 1: the l1 bounds (pure Python integers)
    bound = {0: 1}
    maxbound = 1
    for k in range(1, n + 1):
        r = rows[k - 1]
        nb = {}
        for S in itertools.combinations(range(n), k):
            mask = 0
            for c in S:
                mask |= 1 << c
            b = 0
            for c in S:
                if l1[r][c]:
                    b += l1[r][c] * bound.get(mask ^ (1 << c), 0)
            if b:
                nb[mask] = b
                maxbound = max(maxbound, b)
        bound = nb
    assert maxbound < 2 ** 62, "int64 overflow is not excluded"
    # pass 2: the polynomials
    prev = {0: np.ones(1, dtype=np.int64)}
    for k in range(1, n + 1):
        r = rows[k - 1]
        cur = {}
        size_out = len(_monomials(nv, k)[0])
        for S in itertools.combinations(range(n), k):
            mask = 0
            for c in S:
                mask |= 1 << c
            acc = {}
            for pos, c in enumerate(S):
                if not M[r][c]:
                    continue
                src = prev.get(mask ^ (1 << c))
                if src is None:
                    continue
                sg = -1 if (k - 1 + pos) % 2 else 1
                for var, cf in M[r][c].items():
                    w = sg * cf
                    acc[var] = acc[var] + w * src if var in acc else w * src
            if not acc:
                continue
            out = np.zeros(size_out, dtype=np.int64)
            for var, vec in acc.items():
                out[_shift_index(nv, k - 1, var)] += vec      # injective index map
            if out.any():
                cur[mask] = out
        prev = cur
    res = prev.get((1 << n) - 1)
    lst, _ = _monomials(nv, n)
    out = {}
    if res is not None:
        for i in np.nonzero(res)[0]:
            out[lst[i]] = rsign * int(res[i])
    return out, maxbound


def eval_multiset_poly(poly, x):
    tot = 0
    for m, c in poly.items():
        t = c
        for v in m:
            t *= x[v]
        tot += t
    return tot


# ----------------------------------------------------------------------------------------
# the radial equations for an ordered list of eight generators
# ----------------------------------------------------------------------------------------


def lift_linear_forms(order):
    """u = sum_i x_i f_order[i] + x_8 kappa as 14 integer linear forms {variable: coeff}."""
    u = [{} for _ in range(14)]
    for i, jj in enumerate(order):
        for k in range(14):
            v = F_COORD[jj][k]
            assert v.denominator == 1
            if v:
                u[k][i] = int(v)
    u[13][8] = u[13].get(8, 0) + 1
    return u


def hessian_linear_forms(order):
    """Hess C(u) for u as in lift_linear_forms: 14 x 14 matrix of integer linear forms."""
    import sympy
    ms, _, _, hess = cubic_data()
    u = lift_linear_forms(order)
    H = [[{} for _ in range(14)] for _ in range(14)]
    for i in range(14):
        for k in range(14):
            if hess[i][k] == 0:
                continue
            d = {}
            for mono, coef in sympy.Poly(hess[i][k], *ms).terms():
                idx = [t for t, ee in enumerate(mono) if ee]
                assert len(idx) == 1 and mono[idx[0]] == 1 and coef.is_Integer
                for var, cf in u[idx[0]].items():
                    d[var] = d.get(var, 0) + int(coef) * cf
            H[i][k] = {v: c for v, c in d.items() if c}
    return H


def radial_data(ck, order, names, seed=20261125):
    """names: the eight coefficient names in the order of `order`, then the name of the height.
    Computes from the definitions, exactly:
       C(u) = C0 + h Q                 (u = sum coeff_i f_order[i] + h kappa),
       D(h) = det Hess C(u)            (polynomial in nine variables),
       N    = Q^4 D(-C0/Q)             (polynomial in the eight coefficients),
    and the factorisation N = -4 Q^2 G^2.  Every step is checked; returns a dict."""
    import sympy
    R = Ring(names)
    hname = R.names[8]
    u = lift_linear_forms(order)
    syms = sympy.symbols(names)
    uu = [sum(cf * syms[v] for v, cf in u[k].items()) for k in range(14)]
    Cfull = R.from_sympy(sympy.expand(cubic_C(uu)))
    parts = R.split(Cfull, hname)
    ck.check("EXACT", "C(g + h kappa) is of degree <= 1 in the height h: C = C0 + h Q",
             set(parts) <= {0, 1} and 1 in parts)
    C0, Q = parts.get(0, {}), parts[1]
    # Q is the derivative of C in the direction kappa = e_13: m0 m6 - m1 m5 on the lift
    Qdef = R.from_sympy(sympy.expand(uu[0] * uu[6] - uu[1] * uu[5]))
    ck.check("EXACT", "Q = (m0 m6 - m1 m5)(g) (= dC/dm13; it does not involve h)", Qdef == Q)
    H = hessian_linear_forms(order)
    t0 = time.time()
    det, l1b = det_linear_forms(H, 9)
    D = {}
    for m, c in det.items():
        e = [0] * 9
        for v in m:
            e[v] += 1
        D[R.key(e)] = c
    rng = random.Random(seed)
    ok = True
    try:
        import flint
        for _ in range(6):
            x = [rng.randint(-9, 9) for _ in range(9)]
            Mi = flint.fmpz_mat(14, 14, [sum(c * x[v] for v, c in H[i][k].items())
                                         for i in range(14) for k in range(14)])
            ok &= int(Mi.det()) == eval_multiset_poly(det, x)
        how = "python-flint fmpz_mat.det"
    except ImportError:
        for _ in range(3):
            x = [rng.randint(-9, 9) for _ in range(9)]
            Mi = sympy.Matrix(14, 14, lambda i, k: sum(c * x[v] for v, c in H[i][k].items()))
            ok &= int(Mi.det()) == eval_multiset_poly(det, x)
        how = "sympy integer determinant"
    ck.check("EXACT", "det Hess C(g + h kappa) as a polynomial in nine variables: %d terms (subset "
             "Laplace expansion, int64 with l1 bound 2^%d < 2^62); equal to the integer determinant "
             "(%s) at random integer points (seed %d)" % (len(D), l1b.bit_length(), how, seed), ok,
             "%.1f s" % (time.time() - t0))
    Dk = R.split(D, hname)
    ck.check("EXACT", "det Hess C(g + h kappa) has degree <= 4 in h (rank of Hess C(kappa) is 4)",
             max(Dk) <= 4, "h-degrees present: %s" % sorted(Dk))
    A = Dk.get(0, {})
    mC = p_scale(C0, -1)
    pw = {0: 1}
    assert 14 + 4 * 3 < DEGLIM
    for k in range(1, 5):
        pw = p_mul(pw, mC)
        A = p_add(p_mul(A, Q), p_mul(Dk.get(k, {}), pw))
    N = A
    M1 = p_divexact(N, Q)
    M2 = p_divexact(M1, Q) if M1 is not None else None
    okdiv = M2 is not None and all(c % 4 == 0 for c in M2.values())
    G = p_sqrt({k: -c // 4 for k, c in M2.items()}) if okdiv else None
    ident = G is not None and p_scale(p_mul(p_mul(Q, Q), p_mul(G, G)), -4) == N
    ck.check("EXACT", "N := Q^4 det Hess C(g - (C0/Q) kappa) = sum_k D_k (-C0)^k Q^(4-k) has %d terms and "
             "N = -4 Q^2 G^2 with G of %d terms, homogeneous of degree %d (own exact division and "
             "square root; the identity is re-verified by multiplying out)"
             % (len(N), len(G or {}), R.tdeg(G or {})), ident)
    return dict(R=R, C0=C0, Q=Q, D=D, N=N, G=G)


def compare_radial_json(ck, rd, relpath, alias=None):
    """Comparison only: the supplement's stored C, Q and degree-nine factor (read as JSON data)
    against the polynomials rebuilt here."""
    import json
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: the stored radial polynomials are not compared")
        return
    path = os.path.join(root, relpath)
    if not os.path.isfile(path):
        ck.info("not found under MOC4_SUPP: %s (comparison skipped)" % relpath)
        return
    with open(path) as fh:
        data = json.load(fh)
    R = rd["R"]

    def load(terms):
        return {R.key(e): int(c) for e, c in terms}
    facs = {f["degree"]: load(f["terms"]) for f in data["radial_factors"]}
    okC = load(data["C_terms"]) == rd["C0"]
    okQ = load(data["Q_terms"]) == rd["Q"]
    okG = facs.get(9) in (rd["G"], p_scale(rd["G"], -1))
    okQ2 = facs.get(2) in (rd["Q"], p_scale(rd["Q"], -1))
    ck.check("EXACT", "supplement data %s (read as JSON, comparison only): its variable order has the "
             "same eight coefficient names, its C and Q equal C0 and Q, its quadratic factor is +-Q and its "
             "degree-nine factor is +-G" % os.path.basename(relpath),
             [(alias or {}).get(v, v) for v in data["variables"][:8]] == R.names[:8]
             and okC and okQ and okG and okQ2,
             "stored G = %sG" % ("+" if facs.get(9) == rd["G"] else "-"))


# ----------------------------------------------------------------------------------------
# univariate sign certificates (Fractions; Sturm sequences)
# ----------------------------------------------------------------------------------------


def _trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def _polyrem(p, q):
    p = [Fraction(x) for x in p]
    q = [Fraction(x) for x in q]
    while len(p) >= len(q) and p:
        f = p[-1] / q[-1]
        s = len(p) - len(q)
        for i, c in enumerate(q):
            p[i + s] -= f * c
        p = _trim(p)
    return p


def _polygcd(p, q):
    p, q = _trim(p), _trim(q)
    while q:
        p, q = q, _polyrem(p, q)
    return [Fraction(x) / p[-1] for x in p]


def _polydiv(p, q):
    p = [Fraction(x) for x in p]
    q = [Fraction(x) for x in q]
    out = [Fraction(0)] * (len(p) - len(q) + 1)
    while len(p) >= len(q) and p:
        f = p[-1] / q[-1]
        s = len(p) - len(q)
        out[s] = f
        for i, c in enumerate(q):
            p[i + s] -= f * c
        p = _trim(p)
    assert not p, "inexact division"
    return out


def _variations(vals):
    s = [1 if v > 0 else -1 for v in vals if v != 0]
    return sum(1 for a, b in zip(s, s[1:]) if a != b)


def sturm_positive_roots(coeffs):
    """Number of distinct real roots in (0, +infinity) of the polynomial with coefficient list
    coeffs (low degree first), which must have a nonzero constant term.  Sturm's theorem for the
    chain p, p', -rem(...), ... (valid for a non-squarefree p because 0 is not a root)."""
    p = _trim([Fraction(x) for x in coeffs])
    assert p and p[0] != 0
    if len(p) == 1:
        return 0
    seq = [p, _trim([i * p[i] for i in range(1, len(p))])]
    while True:
        r = _polyrem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-x for x in r])
    return _variations([s[0] for s in seq]) - _variations([s[-1] for s in seq])


def uni_positive(cf):
    """cf = {exponent: int}.  True iff x^(-lowest exponent) * p(x) has a positive constant term and
    no root in (0, infinity), i.e. p(x) > 0 for every x > 0."""
    cf = {e: c for e, c in cf.items() if c}
    if not cf:
        return False
    low = min(cf)
    lst = [cf.get(e + low, 0) for e in range(max(cf) - low + 1)]
    return lst[0] > 0 and sturm_positive_roots(lst) == 0


def uni_nonneg(cf):
    """cf = {exponent: int}.  True iff p(x) >= 0 for every x >= 0.  Method: strip the power of x;
    squarefree decomposition by gcds over Q (Yun); the product r of the squarefree factors of odd
    multiplicity has the sign of p away from finitely many points, and r is squarefree, so
    p >= 0 on (0, inf) iff r(0) > 0 and r has no root in (0, inf) (Sturm)."""
    cf = {e: c for e, c in cf.items() if c}
    if not cf:
        return True
    if all(c >= 0 for c in cf.values()):
        return True
    low = min(cf)
    p = [Fraction(cf.get(e + low, 0)) for e in range(max(cf) - low + 1)]
    lead = p[-1]
    if lead < 0 or p[0] < 0:
        return False
    # Yun's squarefree decomposition p = lead * prod a_i^i
    dp = _trim([i * p[i] for i in range(1, len(p))])
    if not dp:
        return p[0] >= 0
    g = _polygcd(p, dp)
    b = _polydiv(p, g)
    c = _polydiv(dp, g)
    d = _trim([x - y for x, y in itertools.zip_longest(c, _trim([i * b[i] for i in range(1, len(b))]),
                                                         fillvalue=Fraction(0))])
    odd = [Fraction(1)]
    i = 1
    while len(b) > 1:
        a = _polygcd(b, d) if d else b
        a = [x / a[-1] for x in a]
        if i % 2 and len(a) > 1:
            # multiply into the odd part
            new = [Fraction(0)] * (len(odd) + len(a) - 1)
            for s, x in enumerate(odd):
                for t, y in enumerate(a):
                    new[s + t] += x * y
            odd = new
        b = _polydiv(b, a)
        if d:
            c = _polydiv(d, a)
        else:
            c = []
        d = _trim([x - y for x, y in itertools.zip_longest(c, _trim([k * b[k] for k in range(1, len(b))]),
                                                             fillvalue=Fraction(0))])
        i += 1
        if i > 200:
            raise AssertionError("squarefree decomposition did not terminate")
    if len(odd) == 1:
        return lead > 0
    # odd is monic and squarefree with odd(0) != 0; sign of p on (0, inf) = sign(lead * odd)
    if odd[0] <= 0:
        return False
    return sturm_positive_roots(odd) == 0


def uni_sympy_nonneg(cf):
    """Independent cross-check with sympy: p >= 0 on [0, inf) via factor_list and count_roots."""
    import sympy
    x = sympy.Symbol("x")
    P_ = sympy.Poly.from_dict({(e,): c for e, c in cf.items() if c}, x)
    if P_.is_zero:
        return True
    unit, facs = sympy.factor_list(P_)
    sign = 1 if unit > 0 else -1
    for f, m in facs:
        if m % 2:
            if f.eval(0) == 0:
                # the factor x
                continue
            if f.count_roots(0, sympy.oo) != 0:
                return False
            if f.eval(0) < 0:
                sign = -sign
    return sign > 0


def uni_mul(p, q):
    out = {}
    for a, x in p.items():
        for b, y in q.items():
            out[a + b] = out.get(a + b, 0) + x * y
    return {k: v for k, v in out.items() if v}


def uni_sub(p, q):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) - v
    return {k: v for k, v in out.items() if v}


def grouped_univariate(R, poly, var):
    """Group `poly` by the monomial in all variables except `var`.  Returns (number of groups,
    list of the univariate coefficient polynomials {exponent: int} that have a negative
    coefficient, with their outer keys)."""
    groups = R.group(poly, [var])
    s = SH * R.pos[var]
    bad = []
    for outer, inner in groups.items():
        if min(inner.values()) < 0:
            bad.append((outer, {(k >> s) & MASK: c for k, c in inner.items()}))
    return len(groups), bad


# ----------------------------------------------------------------------------------------
# the group with explicit elements; transport of the universal inequalities
# ----------------------------------------------------------------------------------------
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_index_map(g):
    """g = (rho, rho', eps): sigma -> rho o sigma o rho' (eps = 0), rho o sigma^-1 o rho' (eps = 1)."""
    r1, r2, e = g
    return tuple(PIDX[compose(r1, compose(inverse(s) if e else s, r2))] for s in PERMS)


G_EXPLICIT = [(g, g_index_map(g)) for g in G_ELEMENTS]
_FIDX = None


def g_facet(gmap, j):
    """Index of g.f_j, where (g.c) has the value of c at sigma on the permutation g.sigma."""
    global _FIDX
    if _FIDX is None:
        _FIDX = {v: i for i, v in enumerate(F_VALUES)}
    w = [0] * 24
    for s in range(24):
        w[gmap[s]] = F_VALUES[j][s]
    return _FIDX[tuple(w)]


def find_transport(rep, target, swap=None):
    """Group elements g with g.f_rep[m] = f_target[m] for all m; with swap = (i, j) the roles i, j
    of the representative may be exchanged.  Returns list of (g, roles used)."""
    reps = [tuple(rep)]
    if swap is not None:
        sw = list(rep)
        sw[swap[0]], sw[swap[1]] = sw[swap[1]], sw[swap[0]]
        reps.append(tuple(sw))
    out = []
    for rp in reps:
        for g, gm in G_EXPLICIT:
            if all(g_facet(gm, a) == b for a, b in zip(rp, target)):
                out.append((g, rp))
    return out


def g_matrix(g, O):
    """(g.O)[i][j] = O[rho'(i)][rho^-1(j)], after transposing O when eps = 1."""
    r1, r2, e = g
    if e:
        O = [[O[j][i] for j in range(4)] for i in range(4)]
    r1i = inverse(r1)
    return [[O[r2[i]][r1i[j]] for j in range(4)] for i in range(4)]


def leibniz(vals, O):
    return sum(v * sgn(s) * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
               for v, s in zip(vals, PERMS) if v)


def det4(O):
    return leibniz([1] * 24, O)


def random_rational_SO4(rng, size=4):
    """Cayley transform (I - S)(I + S)^-1 of a random rational skew matrix: exactly in SO(4)."""
    S = [[Fraction(0)] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            v = Fraction(rng.randint(-size, size), rng.randint(1, size))
            S[i][j], S[j][i] = v, -v
    ImS = [[(1 if i == j else 0) - S[i][j] for j in range(4)] for i in range(4)]
    IpS_t = [[(1 if i == j else 0) + S[j][i] for j in range(4)] for i in range(4)]
    return [solve_square(IpS_t, ImS[i]) for i in range(4)]


UNIVERSAL = {
    # family: (representative roles, exchangeable role positions, coefficient condition as text)
    "old-H": ((12, 2, 70), None, "t f_t + c f_c + d f_k has Phi >= 0 if t (c - d) >= d^2"),
    "chain": ((2, 0, 12, 70), (1, 2), "a f_centre + p f_leaf + b f_leaf' + d f_k has Phi >= 0 if a >= d^2/p + d^2/b"),
    "four-coefficient": ((0, 6, 10, 54), (0, 1), "a f_a + b f_b + c f_c + d f_d has Phi >= 0 if a b (c - d) >= d^2 (a + b + c + d)"),
    # added for T_44 and T_46 (the representative of PROOF.md 7.6, first row):
    "SOC": ((1, 9, 22), (0, 1), "a f_i + b f_j + d f_k has Phi >= 0 if a b >= d^2"),
}


def generic_covariance(g):
    """EXACT, for a generic 4 x 4 matrix M (entries as position labels): for every permutation
    sigma the monomial prod_i (g.M)[i][(g.sigma)(i)] equals prod_i M[i][sigma(i)], and
    sgn(g.sigma) = eps sgn(sigma) with one sign eps for all sigma.  Hence
    Phi_{g.c}(g.M) = eps Phi_c(M) for every c and det(g.M) = eps det M; for O in SO(4) the
    matrix O'' = g.O (with one column negated if eps = -1) is in SO(4), Phi_{g.c}(O'') = Phi_c(O),
    and O -> O'' is a bijection of SO(4)."""
    Mlab = [[(i, k) for k in range(4)] for i in range(4)]
    gM = g_matrix(g, Mlab)
    gm = g_index_map(g)
    eps = None
    for si, s in enumerate(PERMS):
        tau = PERMS[gm[si]]
        if sorted(gM[i][tau[i]] for i in range(4)) != sorted((i, s[i]) for i in range(4)):
            return False
        e = sgn(tau) * sgn(s)
        if eps is None:
            eps = e
        elif e != eps:
            return False
    return True


def check_transport(ck, family, target, rng, trials=4):
    """EXACT: `target` (facet indices in the roles of the family's representative) lies in the
    G-orbit of the representative of PROOF.md 7.6, and the covariance Phi_{g.f}(O'') = Phi_f(O)
    holds on random rational O in SO(4), where O'' = g.O with one column negated if det = -1."""
    rep, swap, text = UNIVERSAL[family]
    hits = find_transport(rep, target, swap)
    ok = bool(hits)
    detail = ""
    if ok:
        g, roles = hits[0]
        okc = generic_covariance(g)
        for _ in range(trials):
            O = random_rational_SO4(rng)
            O2 = g_matrix(g, O)
            if det4(O2) == -1:
                O2 = [[-row[0]] + list(row[1:]) for row in O2]
            okc &= det4(O) == 1 and det4(O2) == 1
            okc &= all(sum(O2[i][k] * O2[j][k] for k in range(4)) == (1 if i == j else 0)
                       for i in range(4) for j in range(4))
            for a_, b_ in zip(roles, target):
                okc &= leibniz(F_VALUES[a_], O) == leibniz(F_VALUES[b_], O2)
        ok = okc
        detail = ("g = (rho=%s, rho'=%s, transpose=%d) maps roles %s to %s; %d group elements; covariance: "
                  "monomial identity for a generic 4x4 matrix (all 24 permutations) and exact values on %d "
                  "random rational O" % (g[0], g[1], g[2], roles, tuple(target), len(hits), trials))
    ck.check("EXACT", "%s: (%s) is a member of the orbit of the representative (%s) with the same roles"
             % (family, ", ".join("f%d" % t for t in target), ", ".join("f%d" % t for t in rep)), ok, detail)
    return ok


# ----------------------------------------------------------------------------------------
# NUMERIC helpers (floating point; falsification attempts only)
# ----------------------------------------------------------------------------------------


def numeric_min_phi(vals, seed, n_random=20000, n_starts=12):
    """NUMERIC.  Smallest value found of Phi_c(O) over O in SO(4) for the value vector vals:
    n_random random rotations (pairs of unit quaternions) and n_starts local minimisations
    (BFGS on the eight quaternion coordinates, normalised).  An upper bound for the minimum."""
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    S = np.array(PERMS)
    sg = np.array([sgn(s) for s in PERMS], float)
    v = np.array(vals, float)

    def mats(pq):
        p = pq[..., :4] / np.linalg.norm(pq[..., :4], axis=-1, keepdims=True)
        q = pq[..., 4:] / np.linalg.norm(pq[..., 4:], axis=-1, keepdims=True)
        a, b, c, d = [p[..., i] for i in range(4)]
        e, f, g, h = [q[..., i] for i in range(4)]
        Lm = np.stack([np.stack([a, -b, -c, -d], -1), np.stack([b, a, -d, c], -1),
                       np.stack([c, d, a, -b], -1), np.stack([d, -c, b, a], -1)], -2)
        Rm = np.stack([np.stack([e, -f, -g, -h], -1), np.stack([f, e, h, -g], -1),
                       np.stack([g, -h, e, f], -1), np.stack([h, g, -f, e], -1)], -2)
        return Lm @ Rm

    def phi(pq):
        O = mats(pq)
        w = O[..., 0, S[:, 0]] * O[..., 1, S[:, 1]] * O[..., 2, S[:, 2]] * O[..., 3, S[:, 3]]
        return (w * (sg * v)).sum(-1)

    X = rng.normal(size=(n_random, 8))
    vals_r = phi(X)
    best = float(vals_r.min())
    order = np.argsort(vals_r)[:n_starts]
    for i in order:
        res = minimize(lambda x: float(phi(x)), X[i], method="BFGS", options={"gtol": 1e-12, "maxiter": 500})
        best = min(best, float(res.fun))
    return best


def _numeric_gaps(x, T0):
    import numpy as np
    S = np.array(PERMS)
    a, b, al = x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])
    F = (al * np.prod(a[None, :] + b[S], axis=1)).real
    return F - F[T0], F[T0], a, b


def numeric_stratum_search(T, seed, starts, eps=0.02, tol=1e-10):
    """NUMERIC.  Search for spectra a, b in C^4 and a phase alpha with all gaps
    Re(alpha z_sigma) - Re(alpha z_pi*) >= 0, zero on T, sum of gaps = 1 (pi* = first index of T).
    Two Levenberg-Marquardt stages from each random start (17 standard normal parameters), as in
    this owner's replay_iv_common.numeric_stratum_search: stage 1 adds residuals eps/(s_kl + 1e-3)
    that push the eight labels a_i, -b_j apart, stage 2 drops them and polishes.  Returns a list
    of dicts: gaps (24 floats), sep (smallest normalised distance of two labels), H."""
    import numpy as np
    from scipy.optimize import least_squares
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(8, 1)

    def seps(a, b):
        z = np.concatenate([a, -b])
        D = np.abs(z[:, None] - z[None, :])[iu]
        return D / np.sqrt(np.mean(D ** 2))

    def base(x):
        g, H, a, b = _numeric_gaps(x, T[0])
        return np.array(list(g[T[1:]]) + list(np.minimum(0.0, g[others])) + [np.sum(g) - 1.0]), g, H, a, b

    out = []
    for _ in range(starts):
        x0 = rng.normal(size=17)
        try:
            s1 = least_squares(lambda x: np.concatenate([base(x)[0], eps / (seps(*base(x)[3:]) + 1e-3)]),
                               x0, method="lm", max_nfev=600)
            s2 = least_squares(lambda x: base(x)[0], s1.x, method="lm", xtol=1e-15, ftol=1e-15,
                               gtol=1e-15, max_nfev=2000)
        except Exception:
            continue
        r, g, H, a, b = base(s2.x)
        if not np.all(np.isfinite(r)) or float(np.max(np.abs(r))) >= tol:
            continue
        out.append(dict(gaps=g, sep=float(np.min(seps(a, b))), H=float(H)))
    return out


def numeric_equality_points(T, seed, n, tol=1e-12, sep_min=0.05):
    """NUMERIC.  Points x = alpha Psi(a, b) of X whose real parts satisfy only the contact
    EQUATIONS Re(F_sigma - F_pi*) = 0 for sigma in T (the other gaps may have either sign):
    from a random start (17 standard normal parameters) Gauss-Newton with least-norm steps on the
    five equations, scaled by the root mean square of the 24 gaps; the result is kept if the
    scaled residual is < tol and the eight labels a_i, -b_j stay apart (normalised distance
    > sep_min), and is then rescaled to sum of squared gaps = 1.  Returns list of dicts
    gaps, H = Re F_pi*, sep."""
    import numpy as np
    T = sorted(T)
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(8, 1)

    def res(x):
        g, H, a, b = _numeric_gaps(x, T[0])
        return g[T[1:]] / np.sqrt(np.mean(g * g))

    out = []
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        x = rng.normal(size=17)
        ok = False
        for _ in range(100):
            r = res(x)
            if not np.all(np.isfinite(r)):
                break
            if np.max(np.abs(r)) < tol:
                ok = True
                break
            Jm = np.empty((len(r), 17))
            for k in range(17):
                hstep = 1e-7 * max(1.0, abs(x[k]))
                xp = x.copy()
                xp[k] += hstep
                xm = x.copy()
                xm[k] -= hstep
                Jm[:, k] = (res(xp) - res(xm)) / (2 * hstep)
            step = np.linalg.lstsq(Jm, -r, rcond=None)[0]
            nrm = np.linalg.norm(step)
            if nrm > 0.5:
                step *= 0.5 / nrm
            x = x + step
        if not ok:
            continue
        g, H, a, b = _numeric_gaps(x, T[0])
        z = np.concatenate([a, -b])
        D = np.abs(z[:, None] - z[None, :])[iu]
        sep = float(np.min(D) / np.sqrt(np.mean(D ** 2)))
        if sep < sep_min:
            continue
        sc = float(np.sqrt(np.sum(g * g)))
        out.append(dict(gaps=g / sc, H=float(H) / sc, sep=sep))
    return out


def numeric_radial_test(rd, order, pts):
    """NUMERIC.  For points with real part in the span of f_order[0..7] and kappa (gap vectors
    normalised to Euclidean norm 1): coordinates by least squares, then
      fit   = max |gaps - sum coords f|,
      cubic = max |C0 + H Q| / (sum|monomials of C0| + |H| sum|monomials of Q|),
      hess  = max |Q G| / (sum|monomials of Q| * sum|monomials of G|)      (N = -4 Q^2 G^2 = 0),
      qmin  = min |Q| / sum|monomials of Q|   (if this is not small, hess small means G = 0),
      gmax  = max |G| / sum|monomials of G|.
    Returns (fit, cubic, hess, qmin, gmax)."""
    import numpy as np
    R, C0, Q, G = rd["R"], rd["C0"], rd["Q"], rd["G"]
    A = np.array([F_VALUES[j] for j in order], float).T
    absC = {k: abs(v) for k, v in C0.items()}
    absQ = {k: abs(v) for k, v in Q.items()}
    absG = {k: abs(v) for k, v in G.items()}
    worst = [0.0, 0.0, 0.0, float("inf"), 0.0]
    for pnt in pts:
        co, *_ = np.linalg.lstsq(A, pnt["gaps"], rcond=None)
        worst[0] = max(worst[0], float(np.max(np.abs(A @ co - pnt["gaps"]))))
        val = dict(zip(R.names[:8], [float(v) for v in co]))
        val[R.names[8]] = 0.0
        av = {k: abs(v) for k, v in val.items()}
        c0, qv, gv = R.evaluate(C0, val), R.evaluate(Q, val), R.evaluate(G, val)
        rc, rq, rg = R.evaluate(absC, av), R.evaluate(absQ, av), R.evaluate(absG, av)
        worst[1] = max(worst[1], abs(c0 + pnt["H"] * qv) / (rc + abs(pnt["H"]) * rq))
        worst[2] = max(worst[2], abs(qv * gv) / (rq * rg))
        worst[3] = min(worst[3], abs(qv) / rq)
        worst[4] = max(worst[4], abs(gv) / rg)
    return tuple(worst)


# ----------------------------------------------------------------------------------------
# provenance
# ----------------------------------------------------------------------------------------
DOC_SHA256_V = {}


def check_document(ck, relpath, sha):
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: document hash not compared (%s)" % relpath)
        return
    path = os.path.join(root, relpath)
    if not os.path.isfile(path):
        ck.check("EXACT", "document present under MOC4_SUPP: " + relpath, False, path)
        return
    with open(path, "rb") as fh:
        h = hashlib.sha256(fh.read()).hexdigest()
    ck.check("EXACT", "file under MOC4_SUPP is the one that was read (sha256): %s" % os.path.basename(relpath),
             h == sha, h)


def crosscheck_audit_library(ck, cases):
    """Comparison of the rebuilt objects with the audit library (not used to build anything)."""
    here = os.path.dirname(os.path.abspath(__file__))
    lib = os.path.normpath(os.path.join(here, "..", "..", "zheng-n4-audit", "lib"))
    if not os.path.isfile(os.path.join(lib, "moc4lib.py")):
        ck.info("audit library not found at %s: comparison skipped" % lib)
        return
    ck.section("Comparison with the audit library moc4lib (not used to build any object)")
    sys.path.insert(0, lib)
    try:
        import moc4lib as L
    finally:
        sys.path.pop(0)
    import sympy
    ms, C, _, _ = cubic_data()
    ck.check("EXACT", "PERMS, P, kappa, the 82 value vectors, the minor flags, the coordinates of all f_j "
             "and the cubic C agree with moc4lib",
             L.PERMS == PERMS and L.P == P and list(L.KAPPA) == KAPPA
             and [tuple(v) for v in L.F_VALUES] == F_VALUES and list(L.F_IS_MINOR) == F_IS_MINOR
             and [list(v) for v in L.F_VECS] == F_COORD and sympy.expand(L.cubic_C(ms) - C) == 0)
    ck.check("EXACT", "the group has the same 1152 index maps as moc4lib; CONTACT_SETS and generators agree "
             "for the cases %s" % list(cases),
             sorted(set(L._G_MAPS)) == G_MAPS
             and all(tuple(sorted(L.CONTACT_SETS[k])) == PAPER_T_V[k] for k in cases)
             and all(L.generators(PAPER_T_V[k]) == generators(PAPER_T_V[k]) for k in cases))


def crosscheck_supplement_tables(ck):
    """Comparison only: the supplement's stored P, cubic and facet table (JSON data)."""
    import json
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: the supplement's stored P, cubic and facet table are not compared")
        return
    f1 = os.path.join(root, "proofs", "manin_transfer_20261021", "root_scout_moc_spinor_20261021.json")
    f2 = os.path.join(root, "proofs", "principal_minor_result_20261022",
                      "root_permutation_polytope_complete_20261022.json")
    if not (os.path.isfile(f1) and os.path.isfile(f2)):
        ck.info("supplement tables not found under MOC4_SUPP: comparison skipped")
        return
    import sympy
    with open(f1) as fh:
        meta = json.load(fh)
    with open(f2) as fh:
        fac = json.load(fh)
    ms, C, _, _ = cubic_data()
    C2 = sum(sympy.Rational(str(co)) * sympy.prod([ms[int(j)] for j in ix])
             for ix, co in meta["cubic_terms"])
    ck.check("EXACT", "supplement data (JSON, comparison only): its 24 x 14 matrix equals P, its cubic "
             "equals the paper's C with the same sign, and its 82 facet value vectors equal F_VALUES in "
             "the same order (the documents' f_j are the paper's f_j)",
             meta["permutation_to_basis"] == P and sympy.expand(C2 - C) == 0
             and [tuple(z["values"]) for z in fac["facets"]] == F_VALUES)


# ========================================================================================
# PART 3.  python-flint polynomials for the large chart tables (written for T_44 and T_46).
# ========================================================================================


class PR:
    """Z[names].  Elements are flint.fmpz_mpoly objects in len(names) variables.  Only
    +, -, * (with a polynomial or an int on the right), **, exact division //, gcd and term
    access of that class are used."""

    def __init__(self, names):
        from flint import fmpz_mpoly
        self.names = list(names.split(",")) if isinstance(names, str) else list(names)
        self.n = len(self.names)
        assert self.n >= 1 and len(set(self.names)) == self.n
        self.gens = tuple(fmpz_mpoly.gens(self.n))
        self.zero = self.gens[0] * 0
        self.one = self.gens[0] ** 0
        self.idx = {nm: i for i, nm in enumerate(self.names)}

    def __getitem__(self, name):
        return self.gens[self.idx[name]]

    def const(self, c):
        from flint import fmpz
        return self.one * fmpz(int(c))

    @staticmethod
    def terms(p):
        """{exponent tuple: int coefficient} of a ring element."""
        return {tuple(int(e) for e in p.exponent_tuple(i)): int(p.coefficient(i)) for i in range(len(p))}

    def compose(self, d, nums, dens=None, degs=None):
        """For a term dict d = {exponent tuple: int} in len(nums) source variables:
              sum_terms  c * prod_v nums[v]^e_v * dens[v]^(degs[v] - e_v),
        i.e. the numerator of d(x_v = nums[v]/dens[v]) over prod_v dens[v]^degs[v].  dens[v] may
        be None (no denominator); degs[v] defaults to the degree of d in variable v.  nums and
        dens are elements of this ring.  Recursion over the variables with cached powers."""
        nv = len(nums)
        dens = list(dens) if dens is not None else [None] * nv
        actual = [max((ex[v] for ex in d), default=0) for v in range(nv)]
        degs = list(degs) if degs is not None else actual
        assert all(len(ex) == nv for ex in d) and all(degs[v] >= actual[v] for v in range(nv))
        cache = {}

        def factor(v, e):
            key = (v, e)
            if key not in cache:
                f = nums[v] ** e if e else self.one
                if dens[v] is not None and degs[v] - e > 0:
                    f = f * dens[v] ** (degs[v] - e)
                cache[key] = f
            return cache[key]

        def rec(items, v):
            if v == nv:
                tot = 0
                for _, c in items:
                    tot += c
                return self.const(tot)
            groups = {}
            for ex, c in items:
                groups.setdefault(ex[v], []).append((ex, c))
            out = self.zero
            for e in sorted(groups):
                inner = rec(groups[e], v + 1)
                if len(inner) == 0:
                    continue
                if e or (dens[v] is not None and degs[v] > 0):
                    inner = inner * factor(v, e)
                out = out + inner
            return out

        return rec(list(d.items()), 0)

    def from_terms(self, d):
        """Ring element with the given {exponent tuple: coefficient}."""
        return self.compose(d, list(self.gens))

    def evaluate(self, p, point):
        """Exact value of a ring element at a list of Fractions/ints (in the order of names)."""
        return eval_terms(PR.terms(p), point)


def eval_terms(d, point):
    """Exact value of a term dict at a point (list of Fractions or ints)."""
    tot = 0
    for ex, c in d.items():
        t = c
        for x, e in zip(point, ex):
            if e:
                t = t * x ** e
        tot += t
    return tot


def packed_terms(Rsrc, p, names):
    """A packed-key polynomial p of the Part 2 ring Rsrc as a term dict {exponent tuple: int}
    in the variable order `names`; every variable of Rsrc that occurs in p must be listed."""
    pos = [Rsrc.pos[nm] for nm in names]
    other = [i for i in range(Rsrc.n) if i not in pos]
    out = {}
    for k, c in p.items():
        ex = Rsrc.exps(k)
        assert all(ex[i] == 0 for i in other), "a variable outside `names` occurs"
        out[tuple(ex[i] for i in pos)] = c
    return out


def split_terms(d, v):
    """Split a term dict by the exponent of variable number v: {e: term dict without v}."""
    out = {}
    for ex, c in d.items():
        out.setdefault(ex[v], {})[ex[:v] + ex[v + 1:]] = c
    return out


def fl_stats(p):
    """(number of terms, number of negative coefficients, smallest coefficient) of an fmpz_mpoly."""
    n = len(p)
    neg = 0
    mn = None
    for i in range(n):
        c = p.coefficient(i)
        if c < 0:
            neg += 1
        if mn is None or c < mn:
            mn = c
    return n, neg, (int(mn) if mn is not None else 0)


def fl_nonneg(p):
    return all(p.coefficient(i) >= 0 for i in range(len(p)))


def fl_positive_on(R, p, positive_vars):
    """Sufficient exact test that p > 0 whenever all variables are >= 0 and the variables named
    in positive_vars are > 0: all coefficients are >= 0 and some monomial involves only
    variables of positive_vars."""
    pv = set(R.idx[nm] for nm in positive_vars)
    found = False
    for i in range(len(p)):
        if p.coefficient(i) < 0:
            return False
        if not found and all((int(e) == 0) or (v in pv) for v, e in enumerate(p.exponent_tuple(i))):
            found = True
    return found


def fl_digest(p):
    """sha256 of the sorted term list 'e1,...,en:coefficient\\n'.  This is the canonical form in
    which the supplement's data files record their tables (learnt by reading a verifier as text);
    equal digests mean equal polynomials in the same variable order."""
    hh = hashlib.sha256()
    for ex, co in sorted(PR.terms(p).items()):
        hh.update((",".join(map(str, ex)) + ":" + str(co) + "\n").encode())
    return hh.hexdigest()


class RF:
    """num/den with num, den in a PR ring, reduced by their gcd, den with positive leading
    coefficient.  Used for the small rational expressions of the charts (exact identities
    between rational functions are tested by cross-multiplication)."""
    __slots__ = ("R", "num", "den")

    def __init__(self, R, num, den=None, reduce=True):
        self.R = R
        if den is None:
            den = R.one
        if isinstance(num, int):
            num = R.const(num)
        if isinstance(den, int):
            den = R.const(den)
        assert len(den) > 0, "zero denominator"
        if len(num) == 0:
            den = R.one
        elif reduce:
            g = num.gcd(den)
            if not g.is_one():
                num, den = num // g, den // g
        if den.leading_coefficient() < 0:
            num, den = -num, -den
        self.num, self.den = num, den

    def _c(self, o):
        return o if isinstance(o, RF) else RF(self.R, o)

    def __add__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den + o.num * self.den, self.den * o.den)
    __radd__ = __add__

    def __neg__(self):
        return RF(self.R, -self.num, self.den, reduce=False)

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.num, self.den * o.den)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den, self.den * o.num)

    def __rtruediv__(self, o):
        return self._c(o) / self

    def __pow__(self, k):
        return RF(self.R, self.num ** k, self.den ** k, reduce=False)

    def __eq__(self, o):
        o = self._c(o)
        return self.num * o.den == o.num * self.den

    def __ne__(self, o):
        return not self.__eq__(o)

    def is_zero(self):
        return len(self.num) == 0

    def nonneg_form(self):
        """True if num has only nonnegative and den only nonnegative coefficients (den != 0):
        then the value is >= 0 wherever the variables are >= 0 and den does not vanish."""
        return fl_nonneg(self.num) and fl_nonneg(self.den)

    def value(self, point):
        return Fraction(eval_terms(PR.terms(self.num), point)) / Fraction(eval_terms(PR.terms(self.den), point))


def eval_terms_common(dicts, point):
    """Exact values of several term dicts (in the same variables) at one rational point, computed
    with integers over the common denominator prod_v den_v^deg_v (deg_v = the largest exponent
    of variable v in the dicts).  Returns a list of Fractions; same values as eval_terms."""
    nv = len(point)
    pt = [Fraction(x) for x in point]
    degs = [max((ex[v] for d in dicts for ex in d), default=0) for v in range(nv)]
    tab = []
    den = 1
    for v in range(nv):
        n_, d_ = pt[v].numerator, pt[v].denominator
        tab.append([n_ ** e * d_ ** (degs[v] - e) for e in range(degs[v] + 1)])
        den *= d_ ** degs[v]
    out = []
    for d in dicts:
        tot = 0
        for ex, c in d.items():
            t = c
            for v in range(nv):
                t *= tab[v][ex[v]]
            tot += t
        out.append(Fraction(tot, den))
    return out


# ========================================================================================
# PART 4.  The "coefficient Gram" sum-of-squares representation on SO(4) (used by T_46).
#          COPIED from this owner's ../cases-v-a/replay_T25_common.py (functions quat_mul,
#          quat_matrix, parse_poly, class GramFace, ldl_certificate, sym_gram_check) and
#          adapted to the ring class PR of Part 3; the mathematics is unchanged.
# ========================================================================================
import json
import re

QNAMES = ["p0", "p1", "p2", "p3", "q0", "q1", "q2", "q3"]


def mat_inv(A):
    """Inverse of a square matrix of Fractions (raises if singular)."""
    n = len(A)
    cols = [solve_square(A, [Fraction(int(i == k)) for i in range(n)]) for k in range(n)]
    return [[cols[k][i] for k in range(n)] for i in range(n)]


def leibniz_poly(vals, O, zero):
    """Phi(O) = sum_sigma sgn(sigma) vals[sigma] prod_i O[i][sigma(i)] for a matrix of flint
    polynomials (the integer factor is multiplied on the right)."""
    tot = zero
    for s, v in zip(PERMS, vals):
        if v:
            tot = tot + O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]] * int(sgn(s) * v)
    return tot


def quat_mul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def quat_matrix(p, q, one, zero):
    """Matrix of x -> p x q in the basis 1, i, j, k: column k is the image of the k-th basis
    quaternion.  (x -> p x q with unit p, q runs over all of SO(4): classical, CITED.)"""
    cols = []
    for k in range(4):
        e = tuple(one if m == k else zero for m in range(4))
        cols.append(quat_mul(quat_mul(p, e), q))
    return [[cols[j][i] for j in range(4)] for i in range(4)]


def parse_poly(text, ring, scale=1):
    """Parse a polynomial written with + - * ** and division by integers (for example
    'p0*p1*q2/2 - p1**2*q3/2') into scale * polynomial; no eval.  Raises unless the scaled
    coefficients are integers and every name is a generator of the PR ring `ring`."""
    s = text.replace(" ", "")
    tot = ring.zero
    for t in re.findall(r"[+-]?[^+-]+", s):
        num = Fraction(-1 if t[0] == "-" else 1)
        t = t.lstrip("+-")
        toks = re.split(r"(\*\*|\*|/)", t)
        mono = ring.one
        k, op = 0, "*"
        while k < len(toks):
            tok = toks[k]
            if tok in ("*", "/"):
                op = tok
                k += 1
                continue
            exp = 1
            if k + 2 < len(toks) and toks[k + 1] == "**":
                exp = int(toks[k + 2])
                k += 2
            if re.fullmatch(r"\d+", tok):
                val = Fraction(int(tok)) ** exp
                num = num * val if op == "*" else num / val
            else:
                if op != "*" or tok not in ring.idx:
                    raise ValueError("cannot parse %r" % text)
                mono = mono * ring[tok] ** exp
            k += 1
        c = num * scale
        if c.denominator != 1:
            raise ValueError("non-integer coefficient in %r" % text)
        tot = tot + mono * int(c)
    return tot


class GramFace:
    """Sum-of-squares representation of the Leibniz sums of one face.

    Data read from the supplement (as text/JSON, never executed): two blocks of eight cubic
    polynomials F_blk in (p, q) ('features'), an invertible rational 8 x 8 matrix T, four
    signed coordinate permutations S_0..S_3 of (p, q), and the list `order` of eight facet
    indices.  With phi_blk = T^{-1} F_blk (the document's basis) define for a symmetric 8 x 8
    matrix Gm

        Sigma(Gm) = sum_{j=0..3} sum_{blk=0,1} phi_blk(S_j y)^t Gm phi_blk(S_j y),   y = (p, q).

    Rebuilt here from the paper's definitions: M(p, q) = matrix of x -> p x q, and
    Phi_{f_m}(M(p,q)) for the eight facets.  Derived here by exact linear algebra on normal
    forms modulo N = |p|^2 - |q|^2 (p0^2 eliminated): the set of pairs (Gm, w) with

        (|p|^2 + |q|^2) * Sigma(Gm)  =  4 * Phi_{sum_m w_m f_{order_m}}(M(p, q))   mod N.   (*)

    On unit quaternions (*) reads Sigma(Gm) = 2 Phi_w(M(p,q)); so Gm >= 0 gives Phi_w >= 0 on
    all of SO(4)."""

    def __init__(self, supp_dir, prefix):
        with open(os.path.join(supp_dir, prefix + "_budget_PSD8_20261022.json")) as fh:
            bud = json.load(fh)
        with open(os.path.join(supp_dir, prefix + "_cubic_reduction_20261022.json")) as fh:
            cub = json.load(fh)
        self.order = [int(x) for x in bud["source_order"]]
        self.stored_R = [[Fraction(x) for x in row] for row in bud["R"]]
        self.keys = [(i, k) for i in range(8) for k in range(i, 8)]
        assert [list(k) for k in self.keys] == bud["entries"]
        Tm = [[Fraction(x) for x in row] for row in bud["T"]]
        Ti = mat_inv(Tm)
        self.Y = Y = PR(QNAMES)
        yv = Y.gens
        self.pp = yv[0] * yv[0] + yv[1] * yv[1] + yv[2] * yv[2] + yv[3] * yv[3]
        self.qq = yv[4] * yv[4] + yv[5] * yv[5] + yv[6] * yv[6] + yv[7] * yv[7]
        feats = [[parse_poly(s, Y, Fraction(2)) for s in blk] for blk in cub["features"]]   # 2 F
        assert len(feats) == 2 and all(len(b) == 8 for b in feats)
        from math import gcd
        den = 1
        for row in Ti:
            for x in row:
                den = den * x.denominator // gcd(den, x.denominator)
        self.scale = (2 * den) ** 2          # (stored phis) = 2 den * phi
        phis = []
        for blk in feats:
            ph = []
            for i in range(8):
                t = Y.zero
                for k in range(8):
                    cf = Ti[i][k] * den
                    if cf:
                        t = t + blk[k] * int(cf)
                ph.append(t)
            phis.append(ph)
        self.signed_perms = [[(int(i), int(s)) for i, s in S] for S in cub["left_four"]]
        ok = True
        for S in self.signed_perms:
            ok = ok and sorted(i for i, _ in S[:4]) == [0, 1, 2, 3] and sorted(i for i, _ in S[4:]) == [4, 5, 6, 7] \
                and all(s in (1, -1) for _, s in S)
        self.signed_perms_ok = ok            # each S_j preserves |p|^2 and |q|^2
        self.allphi = []
        for S in self.signed_perms:
            img = [yv[idx] * sg for idx, sg in S]
            for ph in phis:
                self.allphi.append([Y.compose(PR.terms(f), img) for f in ph])
        nn = self.pp + self.qq
        self.E = {}
        for (i, k) in self.keys:
            tot = Y.zero
            for ph in self.allphi:
                tot = tot + ph[i] * ph[k]
            self.E[(i, k)] = self.red(tot * nn)
        Mq = quat_matrix(yv[:4], yv[4:], Y.one, Y.zero)
        self.M = Mq
        n2 = self.pp * self.qq
        okm = True
        for i in range(4):
            for j in range(4):
                e = Y.zero
                for k in range(4):
                    e = e + Mq[k][i] * Mq[k][j]
                okm = okm and e == (n2 if i == j else Y.zero)
        self.model_ok = okm and leibniz_poly([1] * 24, Mq, Y.zero) == n2 * n2   # M^t M = |p|^2|q|^2 I, det M = (|p|^2|q|^2)^2
        self.cost = [self.red(leibniz_poly(F_VALUES[m], Mq, Y.zero) * (4 * self.scale)) for m in self.order]
        self._solve()

    def red(self, f):
        """normal form modulo N = |p|^2 - |q|^2:  p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2."""
        Y = self.Y
        yv = Y.gens
        sub = self.qq - yv[1] * yv[1] - yv[2] * yv[2] - yv[3] * yv[3]
        groups = {}
        for ex, co in PR.terms(f).items():
            groups.setdefault(ex[0] // 2, {})[(ex[0] % 2,) + ex[1:]] = co
        tot = Y.zero
        spow = Y.one
        for k in range(max(groups) + 1 if groups else 0):
            if k in groups:
                tot = tot + Y.from_terms(groups[k]) * spow
            spow = spow * sub
        return tot

    def _solve(self):
        """All (Gm, w) satisfying (*): rows = monomials, columns = 36 entries of Gm (upper
        triangle, off-diagonal entries counted twice) then the 8 weights (with sign -)."""
        cols = []
        for key in self.keys:
            mult = 1 if key[0] == key[1] else 2
            cols.append({ex: co * mult for ex, co in PR.terms(self.E[key]).items()})
        for cst in self.cost:
            cols.append({ex: -co for ex, co in PR.terms(cst).items()})
        monos = sorted(set().union(*[set(c) for c in cols]))
        rows = [[c.get(m, 0) for c in cols] for m in monos]
        self.n_equations = len(rows)
        perm = list(range(36, 44)) + list(range(36))          # weights first
        R2, piv = rref([[r[c] for c in perm] for r in rows])
        R2 = R2[:len(piv)]
        self.rank = len(piv)
        self.weight_rows = {}          # m -> coefficients over the 36 entries: w_m = sum coeff * G_entry
        self.relations = []            # rows over the 36 entries: sum coeff * G_entry = 0
        for r in R2:
            lead = next(i for i, x in enumerate(r) if x != 0)
            if lead < 8:
                self.weight_rows[lead] = [-x for x in r[8:]]
                assert all(r[i] == 0 for i in range(8) if i != lead)
            else:
                self.relations.append(r[8:])

    def direct_identity(self, Gm_int, w_int):
        """(*) checked directly for an INTEGER symmetric matrix and integer weights, without
        the linear algebra of _solve: returns True iff
        (|p|^2+|q|^2) * sum phi^t Gm phi == 4 * scale * Phi_w(M)  modulo N."""
        Y = self.Y
        nn = self.pp + self.qq
        tot = Y.zero
        for ph in self.allphi:
            for i in range(8):
                for k in range(8):
                    if Gm_int[i][k]:
                        tot = tot + ph[i] * ph[k] * int(Gm_int[i][k])
        lhs = self.red(tot * nn)
        vals = combo_values([(m, int(w)) for m, w in zip(self.order, w_int)])
        rhs = self.red(leibniz_poly(vals, self.M, Y.zero) * (4 * self.scale))
        return lhs == rhs


def ldl_certificate(Gs, order):
    """Symmetric elimination of the sympy matrix Gs along `order`.  Returns (pivots, residual):
    Gs = sum_i col_i col_i^t / pivot_i + residual, the sum over the nonzero pivots.  So Gs is
    positive semidefinite whenever all pivots are >= 0 and the residual is the zero matrix."""
    import sympy as sp
    A = sp.Matrix(Gs)
    piv = []
    for i in order:
        d = sp.cancel(A[i, i])
        if d == 0:
            piv.append(sp.Integer(0))
            continue
        col = A[:, i]
        A = (A - col * col.T / d).applyfunc(sp.cancel)
        piv.append(d)
    return piv, A


def sym_gram_check(ck, gf, label, entries, diag, weights, order, expected_pivots, instance, note=""):
    """One Gram construction of a document, checked exactly:
      entries  {(i,k): sympy expr}  off-diagonal entries (others 0),  diag  list of 8 sympy exprs;
      weights  the eight coefficients the document claims for this Gram matrix;
      order, expected_pivots  an elimination order and the claimed pivots (see ldl_certificate);
      instance {symbol: Fraction}  a rational parameter point for the direct polynomial check.
    Checks: the re-derived relations, the weights, the pivots and zero residual, and the
    identity (*) of GramFace directly at the instance."""
    import sympy as sp
    Gs = sp.zeros(8, 8)
    for i in range(8):
        Gs[i, i] = diag[i]
    for (i, k), v in entries.items():
        Gs[i, k] = v
        Gs[k, i] = v
    rows = [[Gs[i, k] for k in range(8)] for i in range(8)]
    ent = [rows[i][k] for (i, k) in gf.keys]

    def comb(coeffs):
        return sum(sp.Rational(c.numerator, c.denominator) * e for c, e in zip(coeffs, ent) if c != 0)
    rel_ok = all(sp.cancel(comb(r)) == 0 for r in gf.relations)
    w_calc = [comb(gf.weight_rows[m]) for m in range(8)]
    w_ok = all(sp.cancel(wc - we) == 0 for wc, we in zip(w_calc, weights))
    piv, resid = ldl_certificate(Gs, order)
    piv_ok = all(sp.cancel(p - e) == 0 for p, e in zip(piv, expected_pivots)) and resid == sp.zeros(8, 8)
    sub = {s: sp.Rational(v.numerator, v.denominator) for s, v in instance.items()}
    Gnum = Gs.subs(sub)
    wnum = [sp.nsimplify(sp.sympify(w).subs(sub)) for w in weights]
    den = 1
    for x in list(Gnum) + wnum:
        x = sp.nsimplify(x)
        assert x.is_Rational, x
        den = sp.ilcm(den, x.q)
    Gint = [[int(Gnum[i, k] * den) for k in range(8)] for i in range(8)]
    wint = [int(w * den) for w in wnum]
    direct = gf.direct_identity(Gint, wint)
    ck.check("EXACT", "%s: the Gram matrix satisfies the re-derived relations; its weights are %s; elimination in the "
             "order %s has pivots %s and zero residual; identity (*) holds directly at the rational instance %s%s"
             % (label, [str(w) for w in weights], list(order), [str(sp.factor(e)) for e in expected_pivots],
                {str(s): str(v) for s, v in instance.items()}, note),
             rel_ok and w_ok and piv_ok and direct,
             "relations %s, weights %s, pivots %s, direct %s" % (rel_ok, w_ok, piv_ok, direct))
    return rel_ok and w_ok and piv_ok and direct
