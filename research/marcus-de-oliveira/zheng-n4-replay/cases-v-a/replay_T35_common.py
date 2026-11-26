"""Shared exact objects for the independent replay of the group V cases T_35, T_36, T_37 of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (October 2026).

PART A of this module (everything up to the line "PART B") is COPIED from this owner's earlier
helper `../cases-i-iv/replay_iv_common.py` (functions and tables: Checker, permutations, the
group of order 1152, exact linear algebra over Q, bracket polynomials, the basis (A.1), the
matrix P, the cubic (A.2) with gradient and Hessian, Gaussian rationals, matching vectors, the
82 vectors f_j of Appendix D, generators of a contact cone, direct enumeration of the extreme
rays of N_T, and the checks check_foundations, check_radial_necessity, check_minor_leibniz,
check_linear_identities, eliminate, check_elimination; the last two are at the end of the file).
Those objects are rebuilt there from the DEFINITIONS in the paper (Sections 3.1-3.2, 4.1-4.2,
Appendices A, D, E); every row of P is a verified polynomial identity in sixteen variables.
They are re-checked here each time a replay runs (check_foundations) and compared, as finished
objects only, with the audit library `zheng-n4-audit/lib/moc4lib.py`.

PART B is new: group elements with their (rho, rho', transpose) labels and their action on
matrices, the representatives of the universal inequalities on SO(4) of the audit's PROOF.md
7.6 with exact orbit tests, and floating-point helpers (NUMERIC only).

Nothing is imported from, executed from, or read out of the author's supplement, except that
`crosscheck_with_supplement_data` and `check_document_identity` read data files / hash a PDF
when the environment variable MOC4_SUPP is set (comparison only).

Conventions
  * Permutations of {0,1,2,3} are tuples s with s[i] the image of i; PERMS is in
    lexicographic order, so PERMS[j] is the paper's pi_j.
  * Labels 1..8 carry vectors v_l = (x_l, y_l); [ij] = x_i y_j - x_j y_i.  With
    v_i = (a_i, 1), v_{4+j} = (-b_j, 1) one has [i, 4+j] = a_i + b_j.
  * A polynomial in the sixteen variables x_1..x_8, y_1..y_8 is a dict {key: int} where
    key = sum of 8**slot over the factors (slot l-1 for x_l, slot 7+l for y_l).
  * V_R = Q^14 with coordinates m_0..m_13; values(c) = P c is the vector of the 24 numbers
    (P c)_sigma in PERMS order.  P is injective, so two vectors of V_R are equal iff their
    24 values agree.

Python packages: sympy (identities in few variables); numpy, scipy only inside functions whose
names start with `numeric_` (NUMERIC evidence, labelled as such).
"""
from fractions import Fraction
import hashlib
import itertools
import os
import sys
import time

# ========================================================================================
# PART A  (copied from ../cases-i-iv/replay_iv_common.py)
# ========================================================================================

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


def check_linear_identities(ck, idents, label):
    """idents: list of (left index list, right index list) with claimed sum f_left = sum f_right
    (indices may repeat)."""
    bad = []
    for left, right in idents:
        lv = combo_values([(j, 1) for j in left])
        rv = combo_values([(j, 1) for j in right])
        if lv != rv:
            bad.append((left, right))
    ck.check("EXACT", "%s: %d identities sum f_left = sum f_right hold on all 24 values"
             % (label, len(idents)), not bad, "failing: %s" % bad if bad else "")


# ========================================================================================
# PART B  (new for the group V replays)
# ========================================================================================

# ----------------------------------------------------------------------------------------
# group elements with labels; action on the f_j and on matrices
# (g_matrix / transfer_to_SO4 / random_rational_SO4 follow ../cases-i-iv/replay_T3_common.py)
# ----------------------------------------------------------------------------------------
# g = (rho, rho', eps):  g.sigma = rho o sigma o rho' (eps = 0),  rho o sigma^{-1} o rho' (eps = 1)
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_perm(g, s):
    r1, r2, e = g
    s = inverse(s) if e else tuple(s)
    return compose(r1, compose(s, r2))


G_ELEMENT_MAPS = [tuple(PIDX[g_perm(g, s)] for s in PERMS) for g in G_ELEMENTS]
FIDX = {v: j for j, v in enumerate(F_VALUES)}


def g_values(gmap, vals):
    """(g.c) defined by v_{g.sigma}(g.c) = v_sigma(c)."""
    out = [None] * 24
    for j in range(24):
        out[gmap[j]] = vals[j]
    return tuple(out)


# G_FACET[k][j] = index of g_k . f_j   (KeyError here would mean that G does not permute the f_j)
G_FACET = [[FIDX[g_values(gm, F_VALUES[j])] for j in range(82)] for gm in G_ELEMENT_MAPS]


def find_group_elements(rep, target):
    """All indices k with g_k . f_{rep[m]} = f_{target[m]} for every position m."""
    return [k for k in range(len(G_ELEMENTS))
            if all(G_FACET[k][a] == b for a, b in zip(rep, target))]


SGN = [sgn(p) for p in PERMS]


def leibniz_weights(O):
    return [s * O[0][p[0]] * O[1][p[1]] * O[2][p[2]] * O[3][p[3]] for p, s in zip(PERMS, SGN)]


def phi(vals, O):
    """Phi(O) = sum_sigma sgn(sigma) vals[sigma] prod_i O[i][sigma(i)]."""
    return sum(v * w for v, w in zip(vals, leibniz_weights(O)) if v != 0)


def det4(O):
    return sum(leibniz_weights(O))


def g_matrix(g, O):
    """g.O with Phi_{g.c}(g.O) = (+-1) Phi_c(O): for eps = 0, (g.O)[i][j] = O[rho'(i)][rho^{-1}(j)];
    for eps = 1 the matrix is transposed first."""
    r1, r2, e = g
    if e:
        O = [[O[j][i] for j in range(4)] for i in range(4)]
    r1i = inverse(r1)
    return [[O[r2[i]][r1i[j]] for j in range(4)] for i in range(4)]


def transfer_to_SO4(g, O):
    """For O in SO(4): O'' in SO(4) with Phi_{g.c}(O'') = Phi_c(O) for all c (first column of g.O
    negated when det(g.O) = -1)."""
    Mx = g_matrix(g, O)
    if det4(Mx) == -1:
        Mx = [[-row[0]] + list(row[1:]) for row in Mx]
    return Mx


def random_rational_SO4(rng, size=4):
    """(I - S)(I + S)^{-1} for a random rational skew S: exactly special orthogonal."""
    while True:
        Sk = [[Fraction(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i + 1, 4):
                v = Fraction(rng.randint(-size, size), rng.randint(1, size))
                Sk[i][j], Sk[j][i] = v, -v
        ImS = [[(1 if i == j else 0) - Sk[i][j] for j in range(4)] for i in range(4)]
        IpS_t = [[(1 if i == j else 0) + Sk[j][i] for j in range(4)] for i in range(4)]
        O = [solve_square(IpS_t, ImS[i]) for i in range(4)]
        if all(x != 0 for row in O for x in row):
            return O


def is_SO4(O):
    return all(sum(O[i][k] * O[j][k] for k in range(4)) == (1 if i == j else 0)
               for i in range(4) for j in range(4)) and det4(O) == 1


# ----------------------------------------------------------------------------------------
# the universal inequalities on SO(4) of the audit's PROOF.md 7.6 (CITED results)
# ----------------------------------------------------------------------------------------
# role order of the tuples:   SOC (i, j; k)      a f_i + b f_j + d f_k        if ab >= d^2
#                             old-H (t, c; k)    t f_t + c f_c + d f_k        if t(c - d) >= d^2
#                             chain (c; l, m; k) a f_c + p f_l + b f_m + d f_k if a >= d^2/p + d^2/b
#                             four (a, b; c; d)  a f_a + b f_b + c f_c + d f_d if ab(c-d) >= d^2(a+b+c+d)
# pointwise forms, x_f = the minor of O with Phi_f = x_f^2:
#   SOC   Phi_k + 2|x_i x_j| >= 0
#   old-H Phi_c + Phi_k + 2|x_t x_c| >= 0
#   chain Phi_k + 2|x_c|(|x_l| + |x_m|) >= 0
#   four  Phi_d + Phi_c + 2(|x_a x_b| + |x_a x_c| + |x_b x_c|) >= 0
INEQ_REP = {
    "SOC": (1, 9, 22),
    "old-H": (12, 2, 70),
    "chain": (2, 0, 12, 70),
    "four": (0, 6, 10, 54),
}
INEQ_SYM = {"SOC": (0, 1), "old-H": None, "chain": (1, 2), "four": (0, 1)}
INEQ_ORBIT = {"SOC": 144, "old-H": 288, "chain": 576, "four": 288}


def role_group_element(family, target):
    """A group element index k and the source tuple with g_k . source = target position by
    position, source being the representative or its symmetric variant; None if not in the orbit."""
    rep = INEQ_REP[family]
    variants = [rep]
    if INEQ_SYM[family]:
        i, j = INEQ_SYM[family]
        sw = list(rep)
        sw[i], sw[j] = sw[j], sw[i]
        variants.append(tuple(sw))
    for src in variants:
        ks = find_group_elements(src, target)
        if ks:
            return ks[0], src, len(ks)
    return None


def pointwise_value(family, tup, O, sqrt=None):
    """Left side of the pointwise form for the role tuple `tup` at the matrix O (floats or
    exact numbers with a suitable sqrt)."""
    import math
    sq = sqrt or (lambda x: math.sqrt(max(x, 0.0)))
    ph = lambda j: phi(F_VALUES[j], O)
    ab = lambda j: sq(ph(j))
    if family == "SOC":
        i, j, k = tup
        return ph(k) + 2 * ab(i) * ab(j)
    if family == "old-H":
        t, c, k = tup
        return ph(c) + ph(k) + 2 * ab(t) * ab(c)
    if family == "chain":
        c, l, m, k = tup
        return ph(k) + 2 * ab(c) * (ab(l) + ab(m))
    if family == "four":
        a, b, c, d = tup
        return ph(d) + ph(c) + 2 * (ab(a) * ab(b) + ab(a) * ab(c) + ab(b) * ab(c))
    raise ValueError(family)


def check_roles(ck, roles, seed=20261201, trials=5):
    """roles: list of (family, target tuple, description).  EXACT: the target tuple lies in the
    G-orbit of the representative with the same roles; the transfer O -> O'' of PROOF.md 6.6
    for the group element found is checked on random rational matrices of SO(4)."""
    import random
    rng = random.Random(seed)
    found = {}
    for family, target, desc in roles:
        hit = role_group_element(family, target)
        okmin = all(F_IS_MINOR[j] for j in target[:-1]) and not F_IS_MINOR[target[-1]]
        detail = ""
        ok = hit is not None and okmin
        if hit is not None:
            k, src, cnt = hit
            g = G_ELEMENTS[k]
            cov = True
            for _ in range(trials):
                O = random_rational_SO4(rng)
                O2 = transfer_to_SO4(g, O)
                cov &= is_SO4(O) and is_SO4(O2)
                for a, b in zip(src, target):
                    cov &= phi(F_VALUES[a], O) == phi(F_VALUES[b], O2)
            ok &= cov
            detail = ("g = (rho=%s, rho'=%s, inverse=%d) maps %s to %s position by position (%d group "
                      "elements do); Phi_{g.f}(O'') = Phi_f(O) exactly on %d random rational O in SO(4)"
                      % (g[0], g[1], g[2], src, tuple(target), cnt, trials))
            found[(family, tuple(target))] = g
        ck.check("EXACT", "%s roles %s [%s] lie in the orbit of the representative %s (CITED inequality, "
                 "orbit size %d)" % (family, tuple(target), desc, INEQ_REP[family], INEQ_ORBIT[family]),
                 ok, detail)
    return found


def check_coefficient_forms(ck):
    """The algebra that turns the four pointwise forms into the coefficient forms used by the
    documents (polynomial identities; X, Y, Z >= 0 stand for |x_f|)."""
    import sympy
    ck.section("From the pointwise forms to the coefficient forms (algebra only)")
    a, b, c, d, t, p, X, Y, Z = sympy.symbols("a b c d t p X Y Z")
    e1 = sympy.expand(a * (a * X ** 2 + b * Y ** 2 - 2 * d * X * Y) - ((a * X - d * Y) ** 2 + (a * b - d ** 2) * Y ** 2))
    ck.check("EXACT", "SOC: a(aX^2 + bY^2 - 2dXY) = (aX - dY)^2 + (ab - d^2)Y^2; so ab >= d^2, a > 0 give "
             "a Phi_i + b Phi_j + d Phi_k >= d(Phi_k + 2|x_i x_j|) >= 0  (a = 0 forces d = 0)", e1 == 0)
    e2 = sympy.expand(t * (t * X ** 2 + (c - d) * Y ** 2 - 2 * d * X * Y)
                      - ((t * X - d * Y) ** 2 + (t * (c - d) - d ** 2) * Y ** 2))
    ck.check("EXACT", "old-H: t(tX^2 + (c-d)Y^2 - 2dXY) = (tX - dY)^2 + (t(c-d) - d^2)Y^2; so t(c-d) >= d^2 "
             "gives t Phi_t + c Phi_c + d Phi_k >= d(Phi_c + Phi_k + 2|x_t x_c|) >= 0", e2 == 0)
    e3 = sympy.expand(p * b * (p * X ** 2 + b * Y ** 2 + (d ** 2 / p + d ** 2 / b) * Z ** 2 - 2 * d * Z * (X + Y))
                      - (b * (p * X - d * Z) ** 2 + p * (b * Y - d * Z) ** 2))
    ck.check("EXACT", "chain: pb[pX^2 + bY^2 + (d^2/p + d^2/b)Z^2 - 2dZ(X+Y)] = b(pX - dZ)^2 + p(bY - dZ)^2; so "
             "a >= d^2/p + d^2/b gives p Phi_l + b Phi_m + a Phi_c + d Phi_k >= d(Phi_k + 2|x_c|(|x_l|+|x_m|)) >= 0",
             sympy.simplify(e3) == 0)
    Mx = sympy.Matrix([[a, -d, -d], [-d, b, -d], [-d, -d, c - d]])
    cap = a * b * (c - d) - d ** 2 * (a + b + c + d)
    e4 = sympy.expand(Mx.det() - cap)
    e5 = sympy.expand((a * b - d ** 2) * (c - d) - cap - d ** 2 * (a + b + 2 * d))
    e6 = sympy.expand(a * (c - d) * b - cap - d ** 2 * (a + b + c + d))
    ck.check("EXACT", "four: the form aX^2 + bY^2 + (c-d)Z^2 - 2d(XY + XZ + YZ) has matrix of determinant "
             "ab(c-d) - d^2(a+b+c+d) =: cap, and (ab - d^2)(c-d) = cap + d^2(a+b+2d); so cap >= 0, d > 0 force "
             "c > d, ab > d^2, a > 0: all leading minors a, ab - d^2, cap are >= 0 with the first two > 0, the "
             "form is positive semidefinite, and a Phi_a + b Phi_b + c Phi_c + d Phi_d >= d(pointwise form) >= 0",
             e4 == 0 and e5 == 0 and e6 == 0)


# ----------------------------------------------------------------------------------------
# provenance and comparisons
# ----------------------------------------------------------------------------------------

# Appendix E of the paper (p. 60), read from the text extract of the paper and compared with the
# audit library in crosscheck_with_audit_library
PAPER_T = {35: (2, 3, 4, 6, 9, 12), 36: (2, 3, 4, 6, 11, 12), 37: (0, 5, 6, 7, 9, 12)}

DOC_DIR = os.path.join("proofs", "manin_transfer_20261021")
DOC_SHA256 = {
    "final_referee_T34_uniform_20261022.pdf":
        "de4600291a448d367c13fcdbcca8d0f9fe89cff9dc647ece9264d844bda9df95",
    "final_referee_T35_uniform_20261022.pdf":
        "ee867aee728a58491d9b7c12468b1f85f8a138eb432565774c556ff549390c62",
    "independent_moc4_T36_uniform_20261022.pdf":
        "4b9e286cf7e7742c9ffc96c78388996b7524b83fea1a7f93a8aa3694c666e8dc",
}


def supplement_path(*parts):
    root = os.environ.get("MOC4_SUPP")
    if not root:
        return None
    return os.path.join(root, *parts)


def check_document_identity(ck, name):
    """Provenance only: hash of the proof document that was read (never parsed or executed)."""
    path = supplement_path(DOC_DIR, name)
    if path is None:
        ck.info("MOC4_SUPP not set: the document hash is not compared")
        return
    if not os.path.isfile(path):
        ck.check("EXACT", "document present under MOC4_SUPP: " + name, False, path)
        return
    with open(path, "rb") as fh:
        h = hashlib.sha256(fh.read()).hexdigest()
    ck.check("EXACT", "document under MOC4_SUPP is the one that was read (sha256): " + name,
             h == DOC_SHA256[name], h)


def print_versions(ck):
    import sympy
    import flint
    line = "python %s, sympy %s, python-flint %s" % (sys.version.split()[0], sympy.__version__, flint.__version__)
    try:
        import numpy
        import scipy
        line += ", numpy %s, scipy %s (NUMERIC part only)" % (numpy.__version__, scipy.__version__)
    except Exception:
        pass
    ck.info(line)


def crosscheck_with_supplement_data(ck):
    """If MOC4_SUPP is set: compare the supplement's stored P, cubic and facet table (JSON, read as
    data) with the objects rebuilt here.  Shows only that the documents talk about the same objects."""
    import json
    f1 = supplement_path(DOC_DIR, "root_scout_moc_spinor_20261021.json")
    f2 = supplement_path("proofs", "principal_minor_result_20261022",
                         "root_permutation_polytope_complete_20261022.json")
    if f1 is None:
        ck.info("MOC4_SUPP not set: the supplement's stored P, cubic and facet table are not compared")
        return
    ck.section("Supplement data files (read as JSON, untrusted, comparison only)")
    if not (os.path.isfile(f1) and os.path.isfile(f2)):
        ck.info("data files not found under MOC4_SUPP: comparison skipped")
        return
    import sympy
    with open(f1) as fh:
        meta = json.load(fh)
    with open(f2) as fh:
        fac = json.load(fh)
    ms, C, _, _ = cubic_data()
    C2 = sum(sympy.Rational(str(co)) * sympy.prod([ms[int(j)] for j in ix])
             for ix, co in meta["cubic_terms"])
    ck.check("EXACT", "the supplement's stored 24 x 14 matrix equals P, and its stored cubic equals the "
             "paper's C of (A.2) with the same sign (so the documents' 'F' is C)",
             meta["permutation_to_basis"] == P and sympy.expand(C2 - C) == 0)
    ck.check("EXACT", "the supplement's stored 82 facet value vectors equal F_VALUES in the same order "
             "(the documents' f_j are the paper's f_j of Appendix D)",
             [tuple(z["values"]) for z in fac["facets"]] == F_VALUES)


def crosscheck_with_audit_library(ck):
    here = os.path.dirname(os.path.abspath(__file__))
    lib = os.path.normpath(os.path.join(here, "..", "..", "zheng-n4-audit", "lib"))
    if not os.path.isfile(os.path.join(lib, "moc4lib.py")):
        ck.info("audit library not found at %s: comparison skipped" % lib)
        return
    ck.section("Comparison of the rebuilt objects with the audit library moc4lib (not used to build anything)")
    sys.path.insert(0, lib)
    try:
        import moc4lib as L
    finally:
        sys.path.pop(0)
    ck.check("EXACT", "PERMS, P, kappa agree with moc4lib",
             L.PERMS == PERMS and L.P == P and list(L.KAPPA) == KAPPA)
    ck.check("EXACT", "the 82 value vectors F_VALUES, the minor flags and the coordinates of all f_j agree "
             "with moc4lib", [tuple(v) for v in L.F_VALUES] == F_VALUES and list(L.F_IS_MINOR) == F_IS_MINOR
             and [list(v) for v in L.F_VECS] == F_COORD)
    import sympy
    ms, C, _, _ = cubic_data()
    ck.check("EXACT", "the cubic C agrees with moc4lib.cubic_C", sympy.expand(L.cubic_C(ms) - C) == 0)
    ck.check("EXACT", "T_35, T_36, T_37 agree with moc4lib.CONTACT_SETS; the group has the same 1152 maps",
             all(tuple(sorted(L.CONTACT_SETS[k])) == PAPER_T[k] for k in PAPER_T)
             and sorted(set(L._G_MAPS)) == G_MAPS and sorted(set(G_ELEMENT_MAPS)) == G_MAPS)


def check_contact_set(ck, nu, doc_T, doc_label):
    ck.section("Contact set and numbering")
    ck.check("EXACT", "the document's set %s (its label '%s') equals the paper's T_%d of Appendix E as a set "
             "of lexicographic permutation indices (group element: identity)" % (sorted(doc_T), doc_label, nu),
             tuple(sorted(doc_T)) == PAPER_T[nu])
    reps = {k: orbit_rep(v) for k, v in PAPER_T.items()}
    ck.check("EXACT", "orbit representative under G (1152 maps) of the document's set equals that of T_%d: %s; "
             "T_35, T_36, T_37 lie in three different orbits" % (nu, reps[nu]),
             orbit_rep(doc_T) == reps[nu] and len(set(reps.values())) == 3)


def check_generators(ck, T, doc_sources, doc_nonminor):
    ck.section("Generators of the contact cone")
    gens = generators(T)
    ck.check("EXACT", "F_T = {j : f_j = 0 on T} equals the document's list %s" % sorted(doc_sources),
             gens == sorted(doc_sources), "computed %s" % gens)
    nm = [j for j in gens if not F_IS_MINOR[j]]
    ck.check("EXACT", "non-minor generators equal the document's list", nm == sorted(doc_nonminor),
             "computed %s = %s" % (nm, [f_name(j) for j in nm]))
    t0 = time.time()
    d, rays = extreme_rays_of_contact_cone(T)
    want = sorted(F_VALUES[j] for j in gens)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; "
             "does not use the 82-facet theorem) gives exactly the %d vectors f_j, j in F_T"
             % (d, len(gens)), rays == want, "%d rays, %.1f s" % (len(rays), time.time() - t0))
    return gens


def check_hessian_necessity(ck):
    """The second half of the radial necessity (PROOF.md 7.3): det Hess C(Re x) = 0 for x in X."""
    import sympy
    ms, C, grad, hess = cubic_data()
    xs = sympy.symbols("x0:14")
    ys = sympy.symbols("y0:14")
    sx = dict(zip(ms, xs))
    sy = dict(zip(ms, ys))
    sp = {m: x + y for m, x, y in zip(ms, xs, ys)}
    ok = True
    for j in range(14):
        lhs = sum(hess[j][k].xreplace(sp) * (xs[k] - ys[k]) for k in range(14))
        ok &= sympy.expand(lhs - 2 * (grad[j].xreplace(sx) - grad[j].xreplace(sy))) == 0
    ck.check("EXACT", "Hess C(x + y)(x - y) = 2(grad C(x) - grad C(y)) identically (28 variables): for x in X and "
             "y = conj(x), x + y = 2 Re x and x - y = 2i Im x, so Hess C(Re x) Im x = 0", ok)
    eul = True
    for j in range(14):
        eul &= sympy.expand(sum(hess[j][k] * ms[k] for k in range(14)) - 2 * grad[j]) == 0
    ck.check("EXACT", "Hess C(u) u = 2 grad C(u) identically: if Im x = 0 then u = x in X and Hess C(u) u = 0. "
             "Hence det Hess C(Re x) = 0 for every x in X (Im x != 0, or Re x != 0 = Im x, or x = 0)", eul)
    ck.check("EXACT", "C(kappa) = 0, grad C(kappa) = 0 and Hess C(kappa) = e0 e6^t + e6 e0^t - e1 e5^t - e5 e1^t, "
             "so C(g + h kappa) = C(g) + h dC/dm_13(g) is linear in h and the Hessian is affine in h",
             cubic_C(KAPPA) == 0
             and all(g.xreplace(dict(zip(ms, KAPPA))) == 0 for g in grad)
             and all(hess_at(KAPPA)[i][j] == {(0, 6): 1, (6, 0): 1, (1, 5): -1, (5, 1): -1}.get((i, j), 0)
                     for i in range(14) for j in range(14)))


# ----------------------------------------------------------------------------------------
# NUMERIC helpers (floating point; evidence only)
# ----------------------------------------------------------------------------------------


def numeric_so4(x):
    """SO(4) matrix from 8 real parameters (two quaternions, normalised): x -> p x conj(q)."""
    import numpy as np
    p = x[:4] / np.linalg.norm(x[:4])
    q = x[4:] / np.linalg.norm(x[4:])
    a, b, c, d = p
    Lm = np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]])
    a, b, c, d = q
    Rm = np.array([[a, b, c, d], [-b, a, -d, c], [-c, d, a, -b], [-d, -c, b, a]])
    return Lm @ Rm


def numeric_phi(vals, O):
    import numpy as np
    tot = 0.0
    for v, p, s in zip(vals, PERMS, SGN):
        if v:
            tot += v * s * O[0, p[0]] * O[1, p[1]] * O[2, p[2]] * O[3, p[3]]
    return tot


def numeric_min_phi(vals, seed, starts=6):
    """NUMERIC upper bound for min over SO(4) of Phi with the given 24 values (multistart BFGS)."""
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    vals = [float(v) for v in vals]
    best = None
    for _ in range(starts):
        x0 = rng.normal(size=8)
        r = minimize(lambda x: numeric_phi(vals, numeric_so4(x)), x0, method="BFGS",
                     options=dict(gtol=1e-12, maxiter=400))
        if best is None or r.fun < best:
            best = float(r.fun)
    return best


def numeric_pointwise_attack(ck, roles, seed, samples=4000, starts=4):
    """NUMERIC: try to make the pointwise forms negative for the role tuples used (random matrices
    of SO(4) and a local minimisation from the worst ones)."""
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    worst_all = 0.0
    for family, target, desc in roles:
        f = lambda x: pointwise_value(family, target, numeric_so4(x))
        xs = rng.normal(size=(samples, 8))
        vals = np.array([f(x) for x in xs])
        best = float(vals.min())
        for i in np.argsort(vals)[:starts]:
            r = minimize(f, xs[i], method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-15, maxiter=4000))
            best = min(best, float(r.fun))
        worst_all = min(worst_all, best)
        ck.check("NUMERIC", "pointwise form of %s for roles %s [%s]: no value below -1e-9 on %d random matrices "
                 "of SO(4) plus %d local minimisations (seed %d)" % (family, tuple(target), desc, samples, starts, seed),
                 best > -1e-9, "smallest value found %.3e" % best)
    return worst_all


def numeric_stratum_search(T, seed, starts, extra_residual=None, tol=1e-10, x0s=None):
    """NUMERIC.  Search for spectra a, b in C^4 and a phase alpha such that all gaps
    Re(alpha z_sigma) - Re(alpha z_pi*) are >= 0, vanish on T, and sum to 1 (pi* = first index of
    T).  extra_residual(gaps) may add equations (e.g. to force a coefficient to be positive).
    Levenberg-Marquardt from `starts` random starts (17 standard normal parameters), or from the
    given starting vectors x0s; a point is returned when the maximal residual is below `tol`.
    Adapted from ../cases-i-iv/replay_iv_common.py."""
    import numpy as np
    from scipy.optimize import least_squares
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    Sx = np.array(PERMS)
    rng = np.random.default_rng(seed)

    def base(x):
        a, b, al = x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])
        F = (al * np.prod(a[None, :] + b[Sx], axis=1)).real
        g = F - F[T[0]]
        r = list(g[T[1:]]) + list(np.minimum(0.0, g[others])) + [np.sum(g) - 1.0]
        if extra_residual is not None:
            r += list(extra_residual(g))
        return np.array(r), F, g, a, b

    out = []
    if x0s is None:
        x0s = [rng.normal(size=17) for _ in range(starts)]
    for x0 in x0s:
        try:
            s1 = least_squares(lambda x: base(x)[0], x0, method="lm", max_nfev=800)
            s2 = least_squares(lambda x: base(x)[0], s1.x, method="lm", xtol=1e-15, ftol=1e-15,
                               gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        r, F, g, a, b = base(s2.x)
        if not np.all(np.isfinite(r)) or float(np.max(np.abs(r))) >= tol:
            continue
        z = np.concatenate([a, -b])
        iu = np.triu_indices(8, 1)
        Dm = np.abs(z[:, None] - z[None, :])[iu]
        out.append(dict(H=float(F[T[0]]), gaps=g, a=a, b=b, x=s2.x,
                        sep=float(np.min(Dm) / np.sqrt(np.mean(Dm ** 2))),
                        contacts=tuple(int(s) for s in range(24) if abs(g[s]) < 1e-7)))
    return out


def numeric_hard_points(T, order, d_index, seed, targets=(0.004, 0.01), starts=150, n_perturb=300,
                        sigma=0.08, floor=5e-4):
    """NUMERIC.  Points of the stratum of T whose coefficient on the non-minor generator
    order[d_index] (in the basis `order` of the contact nullspace) is > floor, in the
    normalisation sum of gaps = 1.  Phase 1: searches with the extra hinge equation
    min(0, d - target) = 0 for each target.  Phase 2: every point found is perturbed n_perturb
    times in total (Gaussian noise of size sigma on the 17 parameters) and polished again with
    the hinge at `floor`.  Returns a list of (coefficient vector, point)."""
    import numpy as np
    pinv = np.linalg.pinv(np.array([F_VALUES[j] for j in order], float).T)
    dfun = lambda g: pinv[d_index] @ g
    rng = np.random.default_rng(seed)
    found = []
    for k, target in enumerate(targets):
        found += numeric_stratum_search(T, seed + 1 + k, starts,
                                        extra_residual=lambda g, t=target: [min(0.0, dfun(g) - t)])
    if found and n_perturb:
        x0s = [found[int(rng.integers(len(found)))]["x"] + sigma * rng.normal(size=17) for _ in range(n_perturb)]
        found += numeric_stratum_search(T, seed, 0, extra_residual=lambda g: [min(0.0, dfun(g) - floor)], x0s=x0s)
    out = []
    for pt in found:
        co, res = numeric_coefficients(order, pt["gaps"])
        if res < 1e-8 and co[d_index] > floor * 0.999:
            out.append((co, pt))
    return out


def numeric_coefficients(order, gaps):
    """Least-squares coefficients of the gap vector in the generators f_j, j in order
    (returns coefficients and the residual norm)."""
    import numpy as np
    A = np.array([F_VALUES[j] for j in order], float).T
    co, res, rk, sv = np.linalg.lstsq(A, np.asarray(gaps, float), rcond=None)
    return co, float(np.linalg.norm(A @ co - gaps))


# ----------------------------------------------------------------------------------------
# elimination of non-minor generators (eliminate, check_elimination: COPIED from
# ../cases-i-iv/replay_iv_common.py)
# ----------------------------------------------------------------------------------------


def eliminate(coeffs, nonminor, pair_repl, absorb, rng=None):
    """The document's finite elimination, on a dict {j: Fraction >= 0}.
    pair_repl: {(k, l): right list}; absorb: {k: [(resource, right list), ...]}.
    With rng the order of the pair cancellations and of the absorptions is random (the
    document's argument allows any order).  Returns the new dict (zero entries removed)."""
    c = {j: Fraction(v) for j, v in coeffs.items() if v != 0}
    while True:
        live = [k for k in nonminor if c.get(k, 0) > 0]
        if len(live) < 2:
            break
        if rng is not None:
            rng.shuffle(live)
        k, l = live[0], live[1]
        right = pair_repl[(k, l)] if (k, l) in pair_repl else pair_repl[(l, k)]
        s = min(c[k], c[l])
        c[k] -= s
        c[l] -= s
        for j in right:
            c[j] = c.get(j, 0) + s
        c = {j: v for j, v in c.items() if v != 0}
    live = [k for k in nonminor if c.get(k, 0) > 0]
    if live:
        k = live[0]
        steps = list(absorb[k])
        if rng is not None:
            rng.shuffle(steps)
        for res, right in steps:
            s = min(c.get(k, 0), c.get(res, 0))
            if s == 0:
                continue
            c[k] -= s
            c[res] -= s
            for j in right:
                c[j] = c.get(j, 0) + s
            c = {j: v for j, v in c.items() if v != 0}
    return c


def check_elimination(ck, gens, nonminor, pair_repl, absorb, normal_forms, seed, samples):
    """Structure of the elimination argument plus a randomised exact run of it."""
    import random
    nm = sorted(nonminor)
    pairs_needed = set(itertools.combinations(nm, 2))
    have = set(tuple(sorted(k)) for k in pair_repl)
    ck.check("EXACT", "a replacement is listed for every pair of non-minor generators (%d pairs) and "
             "its right side has only minor generators of the cone" % len(pairs_needed),
             have == pairs_needed and all(F_IS_MINOR[j] and j in gens for r in pair_repl.values() for j in r))
    ok = True
    for k in nm:
        res = set(r for r, _ in absorb[k])
        for r, right in absorb[k]:
            ok &= F_IS_MINOR[r] and r in gens
            ok &= all(F_IS_MINOR[j] and j in gens and j not in res for j in right)
    ck.check("EXACT", "absorptions: every resource is a minor generator, and no right side contains a "
             "non-minor generator or a resource of the same non-minor generator", ok)
    ok = True
    for k in nm:
        res = set(r for r, _ in absorb[k])
        support = set(gens) - (set(nm) - {k}) - res
        ok &= support == set(normal_forms[k]) and len(normal_forms[k]) == 8
    ck.check("EXACT", "after exhausting the resources of the surviving non-minor generator k, the "
             "support is F_T minus the other non-minor generators minus the resources of k, "
             "which is the document's eight-element list for each k", ok)
    rng = random.Random(seed)
    ok = True
    n_minor = n_normal = 0
    for _ in range(samples):
        mode = rng.random()
        c = {}
        for j in gens:
            if rng.random() < mode:
                c[j] = Fraction(rng.randint(0, 50), rng.randint(1, 7))
        before = combo_values(c)
        out = eliminate(c, nm, pair_repl, absorb, rng=rng)
        ok &= combo_values(out) == before and all(v > 0 for v in out.values())
        live = [k for k in nm if k in out]
        if not live:
            n_minor += 1
        else:
            ok &= len(live) == 1 and set(out) <= set(normal_forms[live[0]])
            n_normal += 1
    ck.check("EXACT", "randomised exact run of the elimination, in random order, on %d nonnegative rational "
             "coefficient vectors (seed %d): the vector is unchanged, all coefficients stay >= 0, and the result is "
             "a minor decomposition (%d times) or one of the listed normal forms with d > 0 (%d times)"
             % (samples, seed, n_minor, n_normal), ok and n_minor > 0 and n_normal > 0)
