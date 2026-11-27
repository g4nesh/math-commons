"""Shared exact objects and tools for the independent replay of the group V cases
T_72 = {1,3,4,7,8,12} and T_88 = {1,5,6,10,13,14} of D. Zheng, "The Marcus-de Oliveira
conjecture for n = 4" (October 2026).

PART 1 (foundations) and PART 2 (polynomial tools, radial data, transport of the universal
inequalities, numeric helpers) are COPIED from this owner's earlier replay
`../cases-v-a/replay_T29_common.py` (Part 1 there was copied from
`../cases-i-iv/replay_iv_common.py`).  Those objects were rebuilt from the DEFINITIONS in the
paper (tables read from the page images of pp. 35, 36, 59 of the paper).  Changes made in the
copy: the Sturm-sequence section and the reader of another document's radial JSON are left out,
the family "SOC" is added to the table UNIVERSAL, and the table of paper contact sets is the one
of Part 3.  The copied checks are re-run by every script that uses them, so nothing is taken on
trust from the earlier replays.

PART 3 was written for these two cases: the action of the group on labels, on V and on X;
polynomial tables in chart variables with python-flint; a second (fraction-free) computation
of the Hessian determinant; the replay of the identities of the SOC representative (copied
from this owner's `../cases-v-a/replay_T27_common.py`); the three-pair SOC capacity H3.

Nothing is imported from, executed from, or copied out of the author's supplement.  When the
environment variable MOC4_SUPP is set, a few JSON files are read as DATA and COMPARED with the
objects rebuilt here; no step depends on them.  The audit library
`zheng-n4-audit/lib/moc4lib.py` is not used to build any object; it is only imported, when
present, to compare finished objects.

Exact arithmetic: Python ints and fractions.Fraction; sympy for symbolic identities in few
variables; python-flint fmpz_mpoly (exact integer polynomials) for the large chart tables and
for the fraction-free elimination; numpy int64 inside det_linear_forms only, under an a-priori
bound that excludes overflow.  numpy/scipy floating point occurs only in functions whose names
start with `numeric_`.

Conventions (as in the earlier replays)
  * Permutations of {0,1,2,3} are tuples s with s[i] the image of i; PERMS is in
    lexicographic order, so PERMS[j] is the paper's pi_j.
  * Labels 1..8 carry vectors v_l = (x_l, y_l); [ij] = x_i y_j - x_j y_i.  With
    v_i = (a_i, 1), v_{4+j} = (-b_j, 1) one has [i, 4+j] = a_i + b_j.
  * V_R = Q^14 with coordinates m_0..m_13; values(c) = P c is the vector of the 24 numbers
    (P c)_sigma in PERMS order.  P is injective, so vectors of V_R are handled through values.
  * A group element is g = (rho, chi, eps): sigma -> rho o sigma o chi (eps = 0) or
    rho o sigma^-1 o chi (eps = 1); (g.c) has the value of c at sigma on the permutation g.sigma.
"""


# ========================================================================================
# PART 1.  Copied from ../cases-v-a/replay_T29_common.py (same owner; there copied from
#          ../cases-i-iv/replay_iv_common.py).
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
# PART 2.  Copied from ../cases-v-a/replay_T29_common.py (same owner), with the changes listed
#          in the module docstring.
# ========================================================================================
import heapq
import random


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
    "SOC": ((1, 9, 22), (0, 1), "a f_i + b f_j + d f_k has Phi >= 0 if a b >= d^2"),
    "old-H": ((12, 2, 70), None, "t f_t + c f_c + d f_k has Phi >= 0 if t (c - d) >= d^2"),
    "chain": ((2, 0, 12, 70), (1, 2), "a f_centre + p f_leaf + b f_leaf' + d f_k has Phi >= 0 if a >= d^2/p + d^2/b"),
    "four-coefficient": ((0, 6, 10, 54), (0, 1), "a f_a + b f_b + c f_c + d f_d has Phi >= 0 if a b (c - d) >= d^2 (a + b + c + d)"),
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
        ck.info("audit library not found at ../../zheng-n4-audit/lib: comparison skipped")
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
             and all(tuple(sorted(L.CONTACT_SETS[k])) == PAPER_T[k] for k in cases)
             and all(L.generators(PAPER_T[k]) == generators(PAPER_T[k]) for k in cases))


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
# PART 3.  Written for the replays of T_72 and T_88 (new, except where a function says that
#          it is copied from another helper module of this owner).
# ========================================================================================

# Appendix E of the paper (p. 60), entries 72 and 88, read from the text of the PDF (the table
# lists T_1..T_90 in order; entry 72 is {1,3,4,7,8,12} and entry 88 is {1,5,6,10,13,14}).
PAPER_T = {72: (1, 3, 4, 7, 8, 12), 88: (1, 5, 6, 10, 13, 14)}

# ----------------------------------------------------------------------------------------
# the action of the group on the labels, on V and on X
# ----------------------------------------------------------------------------------------
G_BY_MAP = {gm: g for g, gm in G_EXPLICIT}
assert len(G_BY_MAP) == 1152


def g_inverse(g):
    """The group element h with h.(g.sigma) = sigma for every sigma."""
    gm = g_index_map(g)
    inv = [0] * 24
    for s in range(24):
        inv[gm[s]] = s
    return G_BY_MAP[tuple(inv)]


def g_relabelling(g):
    """lam = (lam_1, ..., lam_8), a permutation of the labels 1..8, with
           p_sigma(v o lam) = p_{g.sigma}(v)   for every sigma,   (v o lam)_l = v_{lam_l}.
    For g: sigma -> rho sigma chi put lam_j = chi^-1(j) on the first four labels and
    lam_{4+k} = 4 + rho(k) on the last four.  For g: sigma -> rho sigma^-1 chi put
    lam_k = 4 + rho(k) and lam_{4+m} = chi^-1(m); the four brackets are then reversed, which
    changes the sign four times.  The identity is CHECKED, not assumed, in check_group_on_X."""
    r1, r2, e = g
    r2i = inverse(r2)
    if not e:
        return tuple([r2i[j] + 1 for j in range(4)] + [5 + r1[k] for k in range(4)])
    return tuple([5 + r1[k] for k in range(4)] + [r2i[m] + 1 for m in range(4)])


def relabel_poly16(poly, lam):
    """poly(v o lam) for a sixteen-variable polynomial: x_l -> x_{lam_l}, y_l -> y_{lam_l}."""
    out = {}
    for key, c in poly.items():
        k, new = key, 0
        for slot in range(16):
            ex = k % 8
            k //= 8
            if ex:
                if slot < 8:
                    new += ex * _xk(lam[slot])
                else:
                    new += ex * _yk(lam[slot - 8])
        out[new] = out.get(new, 0) + c
    return {k: v for k, v in out.items() if v}


def g_action_matrix(g):
    """The 14 x 14 matrix M_g (list of rows) with (P M_g x)_{g.sigma} = (P x)_sigma for all x.
    solve_P raises unless the permuted value vector of every basis vector is in the image of P."""
    gm = g_index_map(g)
    cols = []
    for a in range(14):
        vals = [0] * 24
        for s in range(24):
            vals[gm[s]] = P[s][a]
        cols.append(solve_P(vals))
    return [[cols[a][i] for a in range(14)] for i in range(14)]


def check_group_on_X(ck, used):
    """6.6 of the audit's PROOF.md, for all 1152 elements (relabelling) and, for the elements in
    `used` (list of (name, g)), the matrix on V, kappa, the facet vectors and the Leibniz sums."""
    bad = 0
    for g, gm in G_EXPLICIT:
        lam = g_relabelling(g)
        if sorted(lam) != list(range(1, 9)) or \
                any(relabel_poly16(P_SIGMA[s], lam) != P_SIGMA[gm[s]] for s in range(24)):
            bad += 1
    ck.check("EXACT", "for each of the 1152 group elements g there is a permutation lam_g of the eight labels with "
             "p_sigma(v o lam_g) = p_{g.sigma}(v) identically in the sixteen variables (24 identities each); "
             "hence (P Psi(v o lam_h))_{g.sigma} = (P Psi(v))_sigma for h = g^-1, i.e. g maps the image of Psi, "
             "and so X, onto itself", bad == 0, "%d elements fail" % bad if bad else "")
    for name, g in used:
        gm = g_index_map(g)
        Mg = g_action_matrix(g)
        Mh = g_action_matrix(g_inverse(g))
        prod = [[sum(Mg[i][k] * Mh[k][j] for k in range(14)) for j in range(14)] for i in range(14)]
        ok = all(x.denominator == 1 for row in Mg for x in row) and \
            all(x.denominator == 1 for row in Mh for x in row) and \
            prod == [[int(i == j) for j in range(14)] for i in range(14)]
        ok &= [row[13] for row in Mg] == KAPPA
        ok &= all([sum(Mg[i][k] * F_COORD[j][k] for k in range(14)) for i in range(14)] == F_COORD[g_facet(gm, j)]
                  for j in range(82))
        ok &= generic_covariance(g)
        ck.check("EXACT", "%s = (rho=%s, chi=%s, transpose=%d): its matrix M_g on V (defined by "
                 "(P M_g x)_{g.sigma} = (P x)_sigma) is an integer matrix with integer inverse, M_g kappa = kappa, "
                 "M_g f_j = f_{g(j)} for all 82 vectors, and Phi_{g.c}(g.M) = +-Phi_c(M) for a generic 4x4 matrix M "
                 "with the sign of det(g.M)/det(M) (so Phi_{g.c}(O'') = Phi_c(O) with O -> O'' a bijection of SO(4))"
                 % (name, g[0], g[1], g[2]), ok)


# ----------------------------------------------------------------------------------------
# polynomial tables in chart variables with python-flint (fmpz_mpoly: exact integers)
# ----------------------------------------------------------------------------------------


class ChartRing:
    """Z[names] with python-flint's fmpz_mpoly; the same names as sympy symbols.  Only +, -, *, **,
    len, coefficient(i), exponent_tuple(i) and divmod of fmpz_mpoly are used."""

    def __init__(self, names):
        import sympy
        from flint import fmpz_mpoly
        self.names = tuple(names.split(",")) if isinstance(names, str) else tuple(names)
        self.n = len(self.names)
        self.pos = {nm: i for i, nm in enumerate(self.names)}
        self.gens = tuple(fmpz_mpoly.gens(self.n))
        self.zero = self.gens[0] * 0
        self.one = self.gens[0] ** 0
        self.syms = sympy.symbols(self.names)
        self.sym = dict(zip(self.names, self.syms))

    def const(self, c):
        from flint import fmpz
        return self.one * fmpz(int(c))

    def from_terms(self, d):
        out = self.zero
        for ex, c in d.items():
            t = self.const(c)
            for g, e in zip(self.gens, ex):
                if e:
                    t = t * g ** int(e)
            out = out + t
        return out

    def poly(self, expr):
        """A sympy expression that is a polynomial with integer coefficients in the names."""
        import sympy
        P_ = sympy.Poly(sympy.expand(expr), *self.syms)
        d = {}
        for mono, co in P_.terms():
            assert co.is_Integer
            d[tuple(int(e) for e in mono)] = int(co)
        return self.from_terms(d)

    @staticmethod
    def terms(p):
        """{exponent tuple: int coefficient}."""
        return {tuple(int(e) for e in p.exponent_tuple(i)): int(p.coefficient(i)) for i in range(len(p))}

    def reduced(self, expr):
        """(N, D): expanded sympy polynomials with integer coefficients, expr = N / D, gcd(N, D) = 1 over Z
        (so in particular no common integer factor), and the coefficients of D not all negative.  This
        fixes N and D uniquely when D has positive coefficients; every property is asserted."""
        import sympy
        expr = sympy.sympify(expr)
        n, d = sympy.fraction(sympy.cancel(expr))
        n, d = sympy.expand(n), sympy.expand(d)
        pn, pd = sympy.Poly(n, *self.syms), sympy.Poly(d, *self.syms)
        assert all(c.is_Integer for c in pn.coeffs()) and all(c.is_Integer for c in pd.coeffs())
        if all(c < 0 for c in pd.coeffs()):
            n, d, pn, pd = -n, -d, -pn, -pd
        assert sympy.cancel(n / d - expr) == 0
        assert sympy.gcd(pn, pd).is_one
        return n, d

    def compose_cleared(self, gterms, pairs, degs, order=None):
        """sum_nu g_nu prod_i N_i^nu_i D_i^(m_i - nu_i)   for gterms = {nu: g_nu}, pairs[i] = (N_i, D_i)
        (fmpz_mpoly) and degs[i] = m_i >= max nu_i: the numerator of the polynomial with term dict gterms at
        x_i = N_i / D_i over the common denominator prod_i D_i^m_i.  Recursion over the variables in `order`."""
        nv = len(pairs)
        order = list(order) if order is not None else list(range(nv))
        assert sorted(order) == list(range(nv))
        for ex in gterms:
            assert len(ex) == nv and all(0 <= ex[i] <= degs[i] for i in range(nv))
        fac = []
        for (num, den), m in zip(pairs, degs):
            pn, pd = [self.one], [self.one]
            for _ in range(m):
                pn.append(pn[-1] * num)
                pd.append(pd[-1] * den)
            fac.append([pn[e] * pd[m - e] for e in range(m + 1)])

        def rec(items, pos):
            if pos == nv:
                return self.const(sum(c for _, c in items))
            v = order[pos]
            groups = {}
            for ex, c in items:
                groups.setdefault(ex[v], []).append((ex, c))
            out = self.zero
            for e in sorted(groups):
                out = out + rec(groups[e], pos + 1) * fac[v][e]
            return out
        return rec(list(gterms.items()), 0)


def eval_terms(terms, point):
    """Exact value of {exponent tuple: coeff} at a tuple of ints or Fractions."""
    n = len(point)
    pw = [{0: 1} for _ in range(n)]
    tot = 0
    for ex, c in terms.items():
        t = c
        for i in range(n):
            e = ex[i]
            if e:
                d = pw[i]
                if e not in d:
                    d[e] = point[i] ** e
                t = t * d[e]
        tot = tot + t
    return tot


def support_monomial(terms, names, allowed):
    """True iff some monomial with a positive coefficient involves only variables in `allowed`.  If all
    coefficients are >= 0 this means: the polynomial is > 0 at every point of the closed orthant at which
    the variables in `allowed` are > 0."""
    idx = [i for i, nm in enumerate(names) if nm not in allowed]
    return any(c > 0 and all(ex[i] == 0 for i in idx) for ex, c in terms.items())


def table_digest(terms):
    """sha256 of the serialisation used by the supplement's stored table hashes (rows [exponents, coefficient
    as a string], sorted by exponent tuple, compact JSON).  Used for COMPARISON with stored hashes only."""
    import json
    rows = [[list(ex), str(c)] for ex, c in sorted(terms.items())]
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def packed_to_terms(R, p, drop=()):
    """Packed polynomial over the Ring R -> {exponent tuple: coeff}, with the listed variables set to 1
    (their exponents dropped)."""
    keep = [i for i, nm in enumerate(R.names) if nm not in drop]
    out = {}
    for k, c in p.items():
        e = R.exps(k)
        key = tuple(e[i] for i in keep)
        out[key] = out.get(key, 0) + c
    return {k: v for k, v in out.items() if v}


# ----------------------------------------------------------------------------------------
# fraction-free elimination (second, independent computation of the Hessian determinant)
# ----------------------------------------------------------------------------------------


def bareiss_det(M, zero, one):
    """Determinant of a square matrix of fmpz_mpoly by fraction-free (Bareiss) elimination with the
    sparsest available pivot.  Copied from this owner's ../cases-v-a/replay_T27_common.bareiss_det, with two
    counters added.  Returns (det, updates, divisions): `updates` counts every entry update
    M[i][j] <- (pivot * M[i][j] - M[i][k] * M[k][j]) / previous pivot (previous pivot = 1 in the first step),
    `divisions` those with a previous pivot different from 1; every division is by divmod and a nonzero
    remainder raises."""
    n = len(M)
    M = [list(row) for row in M]
    sign = 1
    prev = one
    updates = divisions = 0
    for k in range(n - 1):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                L = len(M[i][j])
                if L and (best is None or L < best[0]):
                    best = (L, i, j)
        if best is None:
            return zero, updates, divisions
        _, pi, pj = best
        if pi != k:
            M[k], M[pi] = M[pi], M[k]
            sign = -sign
        if pj != k:
            for row in M:
                row[k], row[pj] = row[pj], row[k]
            sign = -sign
        piv = M[k][k]
        for i in range(k + 1, n):
            mik = M[i][k]
            for j in range(k + 1, n):
                a = piv * M[i][j]
                if len(mik) and len(M[k][j]):
                    a = a - mik * M[k][j]
                updates += 1
                if k:
                    q, r = divmod(a, prev)
                    if len(r):
                        raise ArithmeticError("Bareiss division is not exact")
                    a = q
                    divisions += 1
                M[i][j] = a
            M[i][k] = zero
        prev = piv
    out = M[n - 1][n - 1]
    return (out if sign > 0 else -out), updates, divisions


def check_bareiss(ck, rd, order):
    """The 'fraction-free elimination with 819 exact divisions' of the documents, with own code, on the
    full homogeneous Hessian (nine variables); the result is compared with the Laplace-expansion
    determinant of radial_data."""
    from flint import fmpz, fmpz_mpoly
    R = rd["R"]
    gens = fmpz_mpoly.gens(9)
    zero, one = gens[0] * 0, gens[0] ** 0
    H = hessian_linear_forms(order)
    Mf = [[sum((gens[v] * fmpz(int(c)) for v, c in H[i][k].items()), zero) for k in range(14)] for i in range(14)]
    t0 = time.time()
    det, updates, divisions = bareiss_det(Mf, zero, one)
    D2 = {R.key(tuple(int(e) for e in det.exponent_tuple(i))): int(det.coefficient(i)) for i in range(len(det))}
    ck.check("EXACT", "second computation of det Hess C(g + h kappa): fraction-free (Bareiss) elimination of the "
             "14 x 14 matrix of linear forms in nine variables with python-flint polynomials; %d entry updates "
             "= sum_{k=1}^{13} k^2 = 819, every division exact (%d divisions by a pivot, the other %d updates "
             "divide by 1); the result has %d terms and equals the Laplace-expansion determinant"
             % (updates, divisions, updates - divisions, len(D2)),
             updates == 819 == sum(k * k for k in range(1, 14)) and D2 == rd["D"], "%.1f s" % (time.time() - t0))
    return D2


# ----------------------------------------------------------------------------------------
# the SOC representative (copied from this owner's ../cases-v-a/replay_T27_common.py, where it was
# copied from ../cases-i-iv/replay_T3_common.py; `phi` there is `leibniz` here)
# ----------------------------------------------------------------------------------------


def minor_spec(j):
    """(I, J) as 0-based tuples for a minor vector f_j."""
    k, a, b = F_TABLE[j]
    assert k == "n"
    return tuple(int(c) - 1 for c in a), tuple(int(c) - 1 for c in b)


def minor_det(O, I, J):
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def quat_mul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def quat_matrix(p, q, one, zero, conj_q=False):
    """Matrix of x -> p x q (or p x conj(q)) in the basis 1, i, j, k: column k is the image
    of the k-th basis quaternion.  `one`, `zero` are the ring's 1 and 0."""
    if conj_q:
        q = (q[0], -q[1], -q[2], -q[3])
    cols = []
    for k in range(4):
        e = tuple(one if m == k else zero for m in range(4))
        cols.append(quat_mul(quat_mul(p, e), q))
    return [[cols[j][i] for j in range(4)] for i in range(4)]


def replay_SOC_representative(ck):
    """SOC representative (f1, f9; f22):  a f1 + b f9 + d f22 has Phi >= 0 on SO(4) when
    a, b, d >= 0 and ab >= d^2.  Statement and proof: audit, foundations/so4-atoms.md 3.1.
    Here: every polynomial identity that proof uses, exactly, in the model x -> p x conj(q)."""
    from sympy.polys.rings import ring
    from sympy.polys.domains import ZZ
    ck.section("SOC representative (f1, f9; f22): identities of the audit's proof (so4-atoms.md 3.1), own code")
    R, p0, p1, p2, p3, q0, q1, q2, q3 = ring("p0,p1,p2,p3,q0,q1,q2,q3", ZZ)
    Pq, Qq = (p0, p1, p2, p3), (q0, q1, q2, q3)
    pp, qq = sum(x * x for x in Pq), sum(x * x for x in Qq)
    M = quat_matrix(Pq, Qq, R(1), R(0), conj_q=True)
    MT = [[M[j][i] for j in range(4)] for i in range(4)]
    MM = matmul(MT, M)
    ck.check("EXACT", "model x -> p x conj(q): M^t M = |p|^2|q|^2 I and det M = (|p|^2|q|^2)^2 "
             "(so M is in SO(4) for unit p, q)",
             all(MM[i][j] == (pp * qq if i == j else 0) for i in range(4) for j in range(4))
             and det4(M) == (pp * qq) ** 2)
    a = [x * x for x in Pq]
    A = [x * x for x in Qq]
    vh = [a[i] * qq - A[i] * pp for i in range(4)]          # v_i = a_i - A_i, homogenised
    x1 = minor_det(M, *minor_spec(1))
    x9 = minor_det(M, *minor_spec(9))
    ck.check("EXACT", "Phi_{f1} = x1^2, Phi_{f9} = x9^2 with x1 = det M[12,34], x9 = det M[14,23]; "
             "x1 = v0 + v1 and x9 = v0 + v3 (homogenised)",
             leibniz(F_VALUES[1], M) == x1 ** 2 and leibniz(F_VALUES[9], M) == x9 ** 2
             and x1 == vh[0] + vh[1] and x9 == vh[0] + vh[3])
    rhs = (2 * (a[0] * A[2] + a[2] * A[0] + a[1] * A[3] + a[3] * A[1]) * pp * qq
           - 4 * p0 * p1 * p2 * p3 * qq ** 2 - 4 * q0 * q1 * q2 * q3 * pp ** 2)
    ck.check("EXACT", "Phi_{f22} = 2(a0A2 + a2A0 + a1A3 + a3A1) - 4 p0p1p2p3 - 4 q0q1q2q3 on S^3 x S^3 (homogenised)",
             leibniz(F_VALUES[22], M) == rhs)
    S, a0, a1, a2, a3, b0, b1, b2, b3 = ring("al0,al1,al2,al3,be0,be1,be2,be3", ZZ)
    al, be = (a0, a1, a2, a3), (b0, b1, b2, b3)
    aa, AA = [x * x for x in al], [x * x for x in be]
    v = [aa[i] - AA[i] for i in range(4)]
    Rr = (2 * (aa[0] * AA[2] + aa[1] * AA[3] + aa[2] * AA[0] + aa[3] * AA[1])
          - 4 * a0 * a1 * a2 * a3 - 4 * b0 * b1 * b2 * b3)
    de = a0 * a2 - a1 * a3
    De = b0 * b2 - b1 * b3
    s_ = v[2] + v[3]
    t_ = v[1] + v[2]
    summ = v[0] + v[1] + v[2] + v[3]
    ck.check("EXACT", "scalar step: R = 2 delta^2 + 2 Delta^2 - 2(v0v2 + v1v3); st = v1v3 - v0v2 + v2 sum(v); "
             "(v0+v1)(v0+v3) = st + (v0-v2) sum(v)",
             Rr == 2 * de ** 2 + 2 * De ** 2 - 2 * (v[0] * v[2] + v[1] * v[3])
             and s_ * t_ == v[1] * v[3] - v[0] * v[2] + v[2] * summ
             and (v[0] + v[1]) * (v[0] + v[3]) == s_ * t_ + (v[0] - v[2]) * summ)
    U, r_, s2, u_, w_ = ring("r,s,u,w", ZZ)
    ck.check("EXACT", "increment lemma: (r^2+u^2)(s^2+w^2) - (rs + uw)^2 = (rw - us)^2",
             (r_ ** 2 + u_ ** 2) * (s2 ** 2 + w_ ** 2) - (r_ * s2 + u_ * w_) ** 2 == (r_ * w_ - u_ * s2) ** 2)
    ck.info("the inequality itself (Phi_22 + 2|x1 x9| >= 0, hence a f1 + b f9 + d f22 >= 0 for ab >= d^2) is CITED: "
            "audit, foundations/so4-atoms.md 3.1; its polynomial identities are replayed here and its scalar "
            "argument was read line by line (see the report)")


# ----------------------------------------------------------------------------------------
# the three-pair SOC capacity H3 (the same statement in both documents)
# ----------------------------------------------------------------------------------------


def H3_below_one(a, b, c):
    """Exact test of H3(a, b, c) < 1 for rationals a, b, c >= 0, where H3 = s/2 if 2m <= s and
    sqrt(m (s - m)) if 2m >= s (s = a + b + c, m = max)."""
    s, m = a + b + c, max(a, b, c)
    if 2 * m <= s:
        return s < 2
    return m * (s - m) < 1


def check_H3_allocation(ck, fa, fb, fc, fk):
    """EXACT identities behind 'a f_a + b f_b + c f_c + d f_k has Phi >= 0 if d <= H3(a, b, c)', given the
    SOC inequality for the three pairs.  (fa, fb, fc, fk are only used in the text of the checks.)"""
    import sympy as sp
    a, b, c, S, r, lam = sp.symbols("a b c S r lambda", positive=True)
    dab, dac, dbc = (a + b - c) / 2, (a + c - b) / 2, (b + c - a) / 2
    ck.check("EXACT", "H3, balanced case (2 max <= a+b+c): the atoms (f%d, f%d; f%d), (f%d, f%d; f%d), (f%d, f%d; f%d) "
             "with both minor coefficients and the non-minor coefficient equal to (a+b-c)/2, (a+c-b)/2, (b+c-a)/2 "
             "(each >= 0 exactly in the balanced case, and x y = delta^2 for each atom) add up to "
             "a f%d + b f%d + c f%d + ((a+b+c)/2) f%d"
             % (fa, fb, fk, fa, fc, fk, fb, fc, fk, fa, fb, fc, fk),
             sp.expand(dab + dac - a) == 0 and sp.expand(dab + dbc - b) == 0 and sp.expand(dac + dbc - c) == 0
             and sp.expand(dab + dac + dbc - (a + b + c) / 2) == 0)
    # dominant coefficient m with leaves l1, l2, S = l1 + l2 > 0; r stands for sqrt(m / S), i.e. m = r^2 S
    l1, l2 = sp.symbols("l1 l2", positive=True)
    m = r ** 2 * (l1 + l2)
    x1, x2 = m * l1 / (l1 + l2), m * l2 / (l1 + l2)
    d1, d2 = l1 * r, l2 * r
    ck.check("EXACT", "H3, dominant case (m >= S = l1 + l2 > 0, r := sqrt(m/S)): the atoms with minor coefficients "
             "(m l1/S, l1), (m l2/S, l2) and non-minor coefficients l1 r, l2 r satisfy x y = delta^2, use "
             "m l1/S + m l2/S = m, and their non-minor coefficients add up to S r with (S r)^2 = m S",
             sp.cancel(x1 * l1 - d1 ** 2) == 0 and sp.cancel(x2 * l2 - d2 ** 2) == 0
             and sp.cancel(x1 + x2 - m) == 0 and sp.cancel((d1 + d2) ** 2 - m * (l1 + l2)) == 0)
    ck.check("EXACT", "H3: scaling the non-minor coefficient of an SOC atom by lambda in [0,1] keeps x y >= delta^2 "
             "(x y - (lambda delta)^2 = (x y - delta^2) + (1 - lambda^2) delta^2); the two formulas for H3 agree "
             "at 2m = s (both give m)",
             sp.expand((a * b - (lam * c) ** 2) - ((a * b - c ** 2) + (1 - lam ** 2) * c ** 2)) == 0
             and sp.simplify(sp.sqrt(a * (2 * a - a)) - a) == 0)


def numeric_H3(a, b, c):
    import math
    s, m = a + b + c, max(a, b, c)
    return s / 2 if 2 * m <= s else math.sqrt(max(m * (s - m), 0.0))
