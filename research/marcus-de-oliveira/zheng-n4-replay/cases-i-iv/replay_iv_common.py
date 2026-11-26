"""Shared exact objects for the independent replay of the group IV cases T_21, T_22, T_23
of D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (October 2026).

Everything here is rebuilt from the DEFINITIONS in the paper (Sections 3.1-3.2, 4.1-4.2 and
Appendices A, D, E; tables read from the page images of pp. 35, 36, 59, 60).  Nothing is
imported from, executed from, or read out of the author's supplement, and the audit library
`zheng-n4-audit/lib/moc4lib.py` is NOT used to build any object: it is only imported, when
present, by `crosscheck_with_audit_library` to compare finished objects.

Exact arithmetic: Python ints, fractions.Fraction, and sympy for polynomial identities in
few variables.  numpy/scipy are imported only inside the functions whose names start with
`numeric_`; their output is NUMERIC evidence and is labelled as such.

Conventions
  * Permutations of {0,1,2,3} are tuples s with s[i] the image of i; PERMS is in
    lexicographic order, so PERMS[j] is the paper's pi_j.
  * Labels 1..8 carry vectors v_l = (x_l, y_l); [ij] = x_i y_j - x_j y_i.  With
    v_i = (a_i, 1), v_{4+j} = (-b_j, 1) one has [i, 4+j] = a_i + b_j.
  * A polynomial in the sixteen variables x_1..x_8, y_1..y_8 is a dict {key: int} where
    key = sum of 8**slot over the factors (slot l-1 for x_l, slot 7+l for y_l); adding
    keys multiplies monomials (all exponents stay below 8).
  * V_R = Q^14 with coordinates m_0..m_13; values(c) = P c is the vector of the 24 numbers
    (P c)_sigma in PERMS order.  P is injective, so two vectors of V_R are equal iff their
    24 values agree; the f_j are handled through their values.
"""
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


# ----------------------------------------------------------------------------------------
# common checks
# ----------------------------------------------------------------------------------------

# sha256 of the documents that were read for this replay (data only; nothing is executed)
DOC_SHA256 = {
    "proofs/manin_transfer_20261021/root_moc4_T20_uniform_20261022.pdf":
        "1ce411c62da4736d675914358ec1eaad99750b69b6eee3f8ebbfb991ad221d2c",
    "proofs/manin_transfer_20261021/root_moc4_T21_T22_uniform_20261022.pdf":
        "016d9480cfb547dc5d4e96658caa8ee8d57d2ac8c85138a5d87271a12fa3c989",
}

# Appendix E, read from the page image of p. 60 of the paper
PAPER_T = {21: (0, 3, 5, 9, 11, 12), 22: (2, 3, 4, 8, 10, 12), 23: (2, 3, 5, 9, 11, 12)}


def check_document_identity(ck, relpath):
    """Provenance only: if MOC4_SUPP is set, confirm that the proof document there is the
    file that was read for this replay.  The file is hashed, never executed or parsed."""
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: the document hash is not compared (no data from the "
                "supplement is needed by this replay)")
        return
    path = os.path.join(root, relpath)
    if not os.path.isfile(path):
        ck.check("EXACT", "document present under MOC4_SUPP: " + relpath, False, path)
        return
    with open(path, "rb") as fh:
        h = hashlib.sha256(fh.read()).hexdigest()
    ck.check("EXACT", "document under MOC4_SUPP is the one that was read (sha256)",
             h == DOC_SHA256[relpath], h)


def print_versions(ck, numeric=True):
    import sympy
    line = "python %s, sympy %s" % (sys.version.split()[0], sympy.__version__)
    if numeric:
        import numpy
        import scipy
        line += ", numpy %s, scipy %s (NUMERIC part only)" % (numpy.__version__, scipy.__version__)
    ck.info(line)


def crosscheck_with_supplement_data(ck):
    """If MOC4_SUPP is set: read two JSON data files of the supplement (as data; nothing is
    executed) and compare the objects stored there, which the documents' companion programs
    use, with the objects rebuilt here from the paper.  This only shows that the documents talk
    about the same P, the same cubic and the same numbering f_0..f_81; it is not evidence for
    any step."""
    import json
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: the supplement's stored P, cubic and facet table are not compared")
        return
    ck.section("Supplement data files (read as JSON, untrusted, comparison only)")
    f1 = os.path.join(root, "proofs", "manin_transfer_20261021", "root_scout_moc_spinor_20261021.json")
    f2 = os.path.join(root, "proofs", "principal_minor_result_20261022",
                      "root_permutation_polytope_complete_20261022.json")
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
    ck.check("EXACT", "the supplement's stored 24 x 14 matrix equals P, and its stored cubic equals "
             "the paper's C of (A.2) (same sign)",
             meta["permutation_to_basis"] == P and sympy.expand(C2 - C) == 0)
    ck.check("EXACT", "the supplement's stored 82 facet value vectors equal F_VALUES in the same order "
             "(the documents' f_j are the paper's f_j of Appendix D)",
             [tuple(z["values"]) for z in fac["facets"]] == F_VALUES)


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


def check_real_part_lemma(ck):
    """Four numbers H + i q_k with product relation F_a F_d = F_b F_c and H != 0 force
    F_a = F_b or F_a = F_c."""
    import sympy
    ck.section("Common real part of a multiplicative rectangle")
    H, qa, qb, qc, qd = sympy.symbols("H q_a q_b q_c q_d", real=True)
    rel = sympy.expand((H + sympy.I * qa) * (H + sympy.I * qd) - (H + sympy.I * qb) * (H + sympy.I * qc))
    re_part = sympy.expand(sympy.re(rel))
    im_part = sympy.expand(sympy.im(rel))
    S = qa + qd - qb - qc
    R = qa * qd - qb * qc
    ck.check("EXACT", "Im(F_a F_d - F_b F_c) = H (q_a + q_d - q_b - q_c)", sympy.expand(im_part - H * S) == 0)
    ck.check("EXACT", "Re(F_a F_d - F_b F_c) = -(q_a q_d - q_b q_c)", sympy.expand(re_part + R) == 0)
    ck.check("EXACT", "(q_a - q_b)(q_a - q_c) = q_a (q_a + q_d - q_b - q_c) - (q_a q_d - q_b q_c), "
             "so H != 0 and the relation force F_a = F_b or F_a = F_c",
             sympy.expand((qa - qb) * (qa - qc) - qa * S + R) == 0)
    # equivalent closed form, eliminating F_d = F_b F_c / F_a
    Fd = (H + sympy.I * qb) * (H + sympy.I * qc) * (H - sympy.I * qa)
    lhs = sympy.expand(sympy.re(sympy.expand(Fd)) - H * (H ** 2 + qa ** 2))
    ck.check("EXACT", "closed form: (H^2 + q_a^2)(Re(F_b F_c / F_a) - H) = -H (q_a - q_b)(q_a - q_c)",
             sympy.expand(lhs + H * (qa - qb) * (qa - qc)) == 0)


def check_rectangle(ck, quad, swaps):
    """quad = (a, b, c, d) with claimed p_a p_d = p_b p_c; swaps = list of
    (i, j, four bracket pairs): claimed p_i - p_j = +- product of the four brackets."""
    a, b, c, d = quad
    ck.check("EXACT", "p_%d p_%d = p_%d p_%d identically in 16 variables" % (a, d, b, c),
             pmul(P_SIGMA[a], P_SIGMA[d]) == pmul(P_SIGMA[b], P_SIGMA[c]))
    for i, j, pairs in swaps:
        diff = padd(P_SIGMA[i], P_SIGMA[j], -1)
        prod = bracket_product(pairs)
        sign = 1 if diff == prod else (-1 if diff == pscale(prod, -1) else 0)
        distinct_labels = all(x != y for x, y in pairs)
        within = sum(1 for x, y in pairs if (x <= 4) == (y <= 4))
        ck.check("EXACT", "p_%d - p_%d = %s%s identically: a product of four brackets, each joining two "
                 "different labels (%d within one spectrum, %d across)"
                 % (i, j, "+" if sign > 0 else "-", "".join("[%d%d]" % pr for pr in pairs),
                    within, 4 - within),
                 sign != 0 and distinct_labels and len(diff) > 0 and within == 2)


def check_collision_delegation(ck, T, swaps, seed=20261126, samples=8):
    """The layer with coincident labels is not proved in the documents; it is delegated to the
    paper's Theorem B.1 (some a_i + b_j = 0) and Lemma B.2 (a repeated eigenvalue), which give
    det(D_a + U D_b U^*) = sum_sigma mu_sigma z_sigma with mu >= 0, sum mu = 1.  Checked here:
    (i) which coincidences the eight-distinct argument really excludes; (ii) the identity (3.5)
    of the paper, Phi_{c(x)}(O) = Re(alpha det(D_a + O D_b O^t) - F_pi*(x)), on exact data, which
    turns a common-weight representation into Phi_{c(x)}(O) = sum mu_sigma gap_sigma >= 0."""
    import random
    ck.section("Coincident labels: what is delegated, and the identity that makes the delegation work")
    need = sorted(set(tuple(sorted(pr)) for _, _, pairs in swaps for pr in pairs))

    def describe(pr):
        k, l = pr
        if k <= 4 and l <= 4:
            return "a_%d = a_%d" % (k, l)
        if k > 4 and l > 4:
            return "b_%d = b_%d" % (k - 4, l - 4)
        return "a_%d + b_%d = 0" % (k, l - 4)
    cross = [pr for pr in need if pr[0] <= 4 < pr[1]]
    ck.check("EXACT", "the argument for H = 0 uses only the non-vanishing of the %d brackets of the two "
             "factorisations: %d of type a_i + b_j (Theorem B.1 if one vanishes) and %d inside one "
             "spectrum (Lemma B.2 if one vanishes)" % (len(need), len(cross), len(need) - len(cross)),
             len(need) == 8 and len(cross) == 4, "; ".join(describe(pr) for pr in need) + "  (1-based)")
    rng = random.Random(seed)
    ok = True
    for _ in range(samples):
        S = [[Fraction(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i + 1, 4):
                v = Fraction(rng.randint(-4, 4), rng.randint(1, 4))
                S[i][j], S[j][i] = v, -v
        ImS = [[(1 if i == j else 0) - S[i][j] for j in range(4)] for i in range(4)]
        IpS = [[(1 if i == j else 0) + S[i][j] for j in range(4)] for i in range(4)]
        # O = (I - S)(I + S)^{-1}: solve X (I + S) = (I - S) row by row
        IpS_t = [list(r) for r in zip(*IpS)]
        O = [solve_square(IpS_t, ImS[i]) for i in range(4)]
        orth = all(sum(O[i][k] * O[j][k] for k in range(4)) == (1 if i == j else 0)
                   for i in range(4) for j in range(4))
        detO = sum(sgn(s) * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]] for s in PERMS)
        ok &= orth and detO == 1

        def rq():
            return GQ(Fraction(rng.randint(-6, 6), rng.randint(1, 4)),
                      Fraction(rng.randint(-6, 6), rng.randint(1, 4)))
        a = [rq() for _ in range(4)]
        b = [rq() for _ in range(4)]
        al = rq()
        Fz = [al * z for z in z_points(a, b)]
        ref = rng.randrange(24)
        gaps = [f.re - Fz[ref].re for f in Fz]
        phi = sum(sgn(s) * gaps[k] * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
                  for k, s in enumerate(PERMS))
        M = [[(a[i] if i == j else GQ(0)) + sum((GQ(O[i][k] * O[j][k]) * b[k] for k in range(4)), GQ(0))
              for j in range(4)] for i in range(4)]
        det = GQ(0)
        for s in PERMS:
            t = GQ(sgn(s))
            for i in range(4):
                t = t * M[i][s[i]]
            det = det + t
        ok &= phi == (al * det).re - Fz[ref].re
    ck.check("EXACT", "identity (3.5): Phi_{c(x)}(O) = Re(alpha det(D_a + O D_b O^t)) - Re F_pi*(x) on %d "
             "random rational O in SO(4) (Cayley transforms) with Gaussian-rational spectra and alpha "
             "(seed %d)" % (samples, seed), ok)
    ck.info("with det = sum mu_sigma z_sigma this gives Phi_{c(x)}(O) = sum mu_sigma gap_sigma >= 0; "
            "Theorem B.1 and Lemma B.2 themselves are CITED, not replayed here")


def all_rectangles():
    """The multiplicative relations p_a p_d = p_b p_c coming from two disjoint transpositions of
    positions: list of (a, b, c, d), each verified as an identity in 16 variables."""
    seen = set()
    out = []
    for s in PERMS:
        for (i, j), (k, l) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
            def sw(p, x, y):
                p = list(p)
                p[x], p[y] = p[y], p[x]
                return tuple(p)
            a, b, c = PIDX[s], PIDX[sw(s, i, j)], PIDX[sw(s, k, l)]
            d = PIDX[sw(sw(s, i, j), k, l)]
            key = frozenset((frozenset((a, d)), frozenset((b, c))))
            if key in seen:
                continue
            seen.add(key)
            assert pmul(P_SIGMA[a], P_SIGMA[d]) == pmul(P_SIGMA[b], P_SIGMA[c])
            out.append((a, b, c, d))
    return out


def rectangle_closure(T):
    """Smallest superset of T containing the fourth corner of every relation of all_rectangles()
    with three corners in the set.  (If H = 0 and all products are nonzero, three purely imaginary
    corners force the fourth corner to be purely imaginary, i.e. to be a contact.)"""
    T = set(T)
    rects = all_rectangles()
    changed = True
    while changed:
        changed = False
        for r in rects:
            if len(T & set(r)) == 3:
                T |= set(r)
                changed = True
    return sorted(T)


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


def check_contact_set(ck, nu, doc_T, doc_name):
    ck.section("Contact set and scope")
    ck.check("EXACT", "the document's set %s (its label '%s') equals the paper's T_%d of Appendix E "
             "as a set of lexicographic indices (group element: identity)" % (sorted(doc_T), doc_name, nu),
             tuple(sorted(doc_T)) == PAPER_T[nu])
    reps = {k: orbit_rep(v) for k, v in PAPER_T.items()}
    ck.check("EXACT", "orbit representative under G of the document's set equals that of T_%d: %s; "
             "T_21, T_22, T_23 lie in three different orbits" % (nu, reps[nu]),
             orbit_rep(doc_T) == reps[nu] and len(set(reps.values())) == 3)


def check_generators(ck, T, doc_sources, doc_nonminor):
    ck.section("Generators of the contact cone")
    gens = generators(T)
    ck.check("EXACT", "F_T = {j : f_j = 0 on T} equals the document's list", gens == sorted(doc_sources),
             "computed %s" % gens)
    nm = [j for j in gens if not F_IS_MINOR[j]]
    ck.check("EXACT", "non-minor generators equal the document's list", nm == sorted(doc_nonminor),
             "computed %s = %s" % (nm, [f_name(j) for j in nm]))
    t0 = time.time()
    d, rays = extreme_rays_of_contact_cone(T)
    want = sorted(F_VALUES[j] for j in gens)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} "
             "(dimension %d, no use of the 82-facet list) gives exactly the %d vectors f_j, j in F_T"
             % (d, len(gens)), rays == want, "%d rays, %.1f s" % (len(rays), time.time() - t0))
    return gens


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


def cubic_on_normal_form(order):
    """sympy symbols r0..r7 and the expanded C(sum_i r_i f_{order[i]})."""
    import sympy
    r = sympy.symbols("r0:8")
    m = [sum(sympy.Rational(F_COORD[j][k].numerator, F_COORD[j][k].denominator) * r[i]
             for i, j in enumerate(order)) for k in range(14)]
    return r, sympy.expand(cubic_C(m))


def check_cubic_identity(ck, name, order, sigma, claimed, d_index):
    """claimed: function of (r0..r7) returning the document's polynomial."""
    import sympy
    r, C = cubic_on_normal_form(order)
    target = sympy.expand(claimed(*r))
    ck.check("EXACT", "%s: %+d * C(sum r_i f_lambda_i) = claimed polynomial, lambda = %s"
             % (name, sigma, list(order)), sympy.expand(sigma * C - target) == 0)
    poly = sympy.Poly(target, *r)
    coeffs = poly.coeffs()
    d = r[d_index]
    rest = sympy.Poly(target - d ** 3, *r)
    ck.check("EXACT", "%s: all %d coefficients of the polynomial are positive integers, the "
             "coefficient of d^3 (d = r%d, the coefficient of f_%d) is 1, and polynomial - d^3 has "
             "only nonnegative coefficients: so it is >= d^3 > 0 on {r >= 0, d > 0}, boundary included"
             % (name, len(coeffs), d_index, order[d_index]),
             all(c > 0 for c in coeffs) and poly.coeff_monomial(d ** 3) == 1
             and all(c >= 0 for c in rest.coeffs()) and not F_IS_MINOR[order[d_index]]
             and all(F_IS_MINOR[j] for i, j in enumerate(order) if i != d_index))
    return sympy.expand(sigma * C)


def exact_sampling_of_normal_form(ck, name, order, sigma, d_index, seed, samples):
    """Randomised EXACT test of sigma * C(g) >= d^3 on the closed region r >= 0, d > 0, with many
    coefficients set to 0 and with very small and very large d (rational arithmetic; the cubic is
    evaluated directly in the 14 matching coordinates, not through the expanded polynomial)."""
    import random
    rng = random.Random(seed)
    worst = None
    ok = True
    for n in range(samples):
        r = []
        mode = rng.randint(0, 3)
        for i in range(8):
            if i == d_index:
                r.append(Fraction(rng.randint(1, 10 ** 6), 10 ** rng.randint(0, 9)))
            elif mode == 0 or rng.random() < (0.0, 0.3, 0.6, 0.9)[mode]:
                r.append(Fraction(0))
            else:
                r.append(Fraction(rng.randint(0, 10 ** 6), 10 ** rng.randint(0, 9)))
        m = combo_coord([(order[i], r[i]) for i in range(8)])
        val = sigma * cubic_C(m)
        d3 = r[d_index] ** 3
        if val < d3:
            ok = False
        ratio = val / d3
        worst = ratio if worst is None or ratio < worst else worst
    ck.check("EXACT", "%s: %d random rational points of {r >= 0, d > 0} incl. faces (seed %d): "
             "sigma*C(g) >= d^3 always" % (name, samples, seed), ok,
             "min of sigma*C(g)/d^3 = %s" % (float(worst),))


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


# ----------------------------------------------------------------------------------------
# optional comparison with the audit library (finished objects only)
# ----------------------------------------------------------------------------------------


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
    ck.check("EXACT", "the 82 value vectors F_VALUES and the minor flags agree with moc4lib",
             [tuple(v) for v in L.F_VALUES] == F_VALUES and list(L.F_IS_MINOR) == F_IS_MINOR)
    ck.check("EXACT", "coordinates of all f_j in V_R agree with moc4lib",
             [list(v) for v in L.F_VECS] == F_COORD)
    import sympy
    ms, C, _, _ = cubic_data()
    ck.check("EXACT", "the cubic C agrees with moc4lib.cubic_C", sympy.expand(L.cubic_C(ms) - C) == 0)
    ck.check("EXACT", "T_21, T_22, T_23 agree with moc4lib.CONTACT_SETS; the group has the same 1152 maps",
             all(tuple(sorted(L.CONTACT_SETS[k])) == PAPER_T[k] for k in PAPER_T)
             and sorted(set(L._G_MAPS)) == G_MAPS)


# ----------------------------------------------------------------------------------------
# NUMERIC falsification helpers (floating point; evidence only)
# ----------------------------------------------------------------------------------------


def numeric_stratum_search(T, gens, seed, starts, h_target=None, eps=0.02, tol=1e-10):
    """NUMERIC.  Search for spectra a, b in C^4 and a phase alpha such that all gaps
    Re(alpha z_sigma) - Re(alpha z_pi*) are >= 0, vanish on T, and sum to 1 (pi* = first index
    of T).  With h_target the extra equation Re(alpha z_pi*) = h_target is imposed.

    Two Levenberg-Marquardt stages from each random start (17 standard normal parameters):
    stage 1 adds the residuals eps / (s_kl + 1e-3) that push the eight labels a_i, -b_j apart
    (s_kl = distance of labels k, l divided by the root mean square of all 28 distances);
    stage 2 drops them and polishes, so a returned point solves the stated equations to `tol`
    in the maximum norm.  Returns a list of dicts with
      H            Re(alpha z_pi*) in the normalisation sum of gaps = 1,
      sep          min s_kl (0 means two labels collide),
      minor_resid  Euclidean distance (NNLS) of the gap vector from the cone spanned by the
                   minor generators of N_T,
      gaps, contacts (indices with |gap| < 1e-7)."""
    import numpy as np
    from scipy.optimize import least_squares, nnls
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    S = np.array(PERMS)
    rng = np.random.default_rng(seed)
    minors = [j for j in gens if F_IS_MINOR[j]]
    A_min = np.array([F_VALUES[j] for j in minors], float).T
    iu = np.triu_indices(8, 1)

    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])

    def seps(a, b):
        z = np.concatenate([a, -b])
        D = np.abs(z[:, None] - z[None, :])[iu]
        return D / np.sqrt(np.mean(D ** 2))

    def base(x):
        a, b, al = unpack(x)
        F = (al * np.prod(a[None, :] + b[S], axis=1)).real
        g = F - F[T[0]]
        r = list(g[T[1:]]) + list(np.minimum(0.0, g[others])) + [np.sum(g) - 1.0]
        if h_target is not None:
            r.append(F[T[0]] - h_target)
        return np.array(r), F, g, a, b

    def res1(x):
        r, F, g, a, b = base(x)
        return np.concatenate([r, eps / (seps(a, b) + 1e-3)])

    def res2(x):
        return base(x)[0]

    out = []
    for _ in range(starts):
        x0 = rng.normal(size=17)
        try:
            s1 = least_squares(res1, x0, method="lm", max_nfev=600)
            s2 = least_squares(res2, s1.x, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                               max_nfev=2000)
        except Exception:
            continue
        r, F, g, a, b = base(s2.x)
        if not np.all(np.isfinite(r)) or float(np.max(np.abs(r))) >= tol:
            continue
        _, rmin = nnls(A_min, np.maximum(g, 0.0))
        z = np.concatenate([a, -b])
        rms = float(np.sqrt(np.mean(np.abs(z[:, None] - z[None, :])[iu] ** 2)))
        out.append(dict(H=float(F[T[0]]), sep=float(np.min(seps(a, b))), minor_resid=float(rmin),
                        gaps=g, contacts=tuple(int(s) for s in range(24) if abs(g[s]) < 1e-7),
                        labels=z / rms))
    return out


def numeric_falsification(ck, nu, T, gens, seed, named_pairs, starts_free=300, starts_forced=80,
                          forced=(0.3, -0.3, 1.0, -1.0), sep_min=1e-3, eps_list=(0.02,), extra=None):
    """NUMERIC attempts to falsify the intermediate claims of the document on the points of the
    stratum to which its argument applies.  named_pairs are the eight label pairs whose brackets
    occur in the two factorisations.  A converged point is called
      covered     if none of the eight named distances is below sep_min (the argument applies),
      8-distinct  if moreover all 28 label distances are above sep_min,
      delegated   if a named distance is below 1e-4 (a named coincidence: Appendix B of the paper).
    Claims tested:  (a) H = 0 at covered points;  (b) the gap vector of a covered point lies in the
    cone of the minor generators;  (a') and (a'') a point with H != 0 is delegated.
    `extra(points, covered_points)` may add case-specific lines.  The free search is repeated for
    every repulsion strength in eps_list (seeds seed, seed + 1, ...), starts_free starts each."""
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")

    def named_sep(p):
        return min(abs(p["labels"][k - 1] - p["labels"][l - 1]) for k, l in named_pairs)

    pts = []
    for k, eps in enumerate(eps_list):
        part = numeric_stratum_search(T, gens, seed + k, starts_free, eps=eps)
        ck.info("free search, repulsion eps = %g: %d starts (seed %d), %d converged (max residual < 1e-10): "
                "%d covered, of which %d 8-distinct; %d delegated"
                % (eps, starts_free, seed + k, len(part),
                   sum(1 for p in part if named_sep(p) > sep_min),
                   sum(1 for p in part if p["sep"] > sep_min),
                   sum(1 for p in part if named_sep(p) < 1e-4)))
        pts += part
    cov = [p for p in pts if named_sep(p) > sep_min]
    if cov:
        from collections import Counter
        cs = Counter(p["contacts"] for p in cov).most_common(3)
        ck.info("covered points: max |H| = %.3e, max distance from the minor cone = %.3e, smallest named "
                "distance in [%.3g, %.3g]; most frequent contact sets: %s"
                % (max(abs(p["H"]) for p in cov), max(p["minor_resid"] for p in cov),
                   min(named_sep(p) for p in cov), max(named_sep(p) for p in cov),
                   "; ".join("%s x%d" % (list(c), n) for c, n in cs)))
        ck.check("NUMERIC", "T_%d (a): no covered point has |H| > 1e-7 (normalisation: sum of gaps = 1)" % nu,
                 all(abs(p["H"]) <= 1e-7 for p in cov), "%d points tested" % len(cov))
        ck.check("NUMERIC", "T_%d (b): no covered point is farther than 1e-7 from the cone of the minor "
                 "generators" % nu, all(p["minor_resid"] <= 1e-7 for p in cov), "%d points tested" % len(cov))
    else:
        ck.info("no covered point was found: claims (a), (b) are NOT tested by this search")
    if pts:
        ck.info("all converged points (delegated ones included; the document claims nothing there): "
                "max |H| = %.3e, max distance from the minor cone = %.3e"
                % (max(abs(p["H"]) for p in pts), max(p["minor_resid"] for p in pts)))
    with_h = [p for p in pts if abs(p["H"]) > 1e-6]
    n_free_h = len(with_h)
    for k, h in enumerate(forced):
        fp = numeric_stratum_search(T, gens, seed + 100 + k, starts_forced, h_target=h)
        with_h += fp
        ck.info("forced H = %+.1f: %d starts (seed %d), %d converged; largest value of the smallest named "
                "distance %.3e" % (h, starts_forced, seed + 100 + k, len(fp),
                                   max([named_sep(p) for p in fp] or [0.0])))
    worst = max([named_sep(p) for p in with_h] or [0.0])
    ck.check("NUMERIC", "T_%d (a'), (a''): every converged point with |H| > 1e-6 (%d from the free search, "
             "%d from the searches with H = h imposed, h in %s) is delegated: one of the eight named "
             "brackets vanishes (normalised distance < 1e-4)"
             % (nu, n_free_h, len(with_h) - n_free_h, list(forced)), worst < 1e-4 and len(with_h) > 0,
             "largest value of the smallest named distance %.3e" % worst)
    if extra is not None:
        extra(pts, cov)


def numeric_hinge_search(T, seed, starts, s0, tol=1e-10):
    """NUMERIC.  One-stage Levenberg-Marquardt search for a point of the stratum of T (gaps >= 0,
    gaps = 0 on T, sum of gaps = 1) whose eight labels keep normalised separation >= s0: the
    residuals min(0, s_kl - s0) vanish exactly on such points, so they are zeros of the system.
    Returns (number of starts that reached max residual < tol, smallest max residual seen,
    largest |H| among the zeros or None)."""
    import numpy as np
    from scipy.optimize import least_squares
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    S = np.array(PERMS)
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(8, 1)

    def parts(x):
        a = x[0:4] + 1j * x[4:8]
        b = x[8:12] + 1j * x[12:16]
        F = (np.exp(1j * x[16]) * np.prod(a[None, :] + b[S], axis=1)).real
        z = np.concatenate([a, -b])
        D = np.abs(z[:, None] - z[None, :])[iu]
        return F, D / np.sqrt(np.mean(D ** 2))

    def res(x):
        F, sp = parts(x)
        g = F - F[T[0]]
        r = list(g[T[1:]]) + list(np.minimum(0.0, g[others])) + [np.sum(g) - 1.0]
        return np.concatenate([np.array(r), np.minimum(0.0, sp - s0)])

    n, best, hmax = 0, float("inf"), None
    for _ in range(starts):
        try:
            s = least_squares(res, rng.normal(size=17), method="lm", xtol=1e-15, ftol=1e-15,
                              gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        rr = float(np.max(np.abs(res(s.x))))
        if not np.isfinite(rr):
            continue
        best = min(best, rr)
        if rr < tol:
            n += 1
            h = abs(float(parts(s.x)[0][T[0]]))
            hmax = h if hmax is None else max(hmax, h)
    return n, best, hmax


def numeric_common_weights(ck, seed, samples=60):
    """NUMERIC evidence for the two CITED results of the paper's Appendix B (they are not part of
    the documents replayed here, and this is no proof of them).

    By the minor expansion, det(D_a + U D_b U^*) = sum_{I,J} prod_{i not in I} a_i prod_{j in J} b_j
    |det U[I,J]|^2, and z_sigma is the same sum with |det U[I,J]|^2 replaced by [sigma(I) = J].
    Common weights mu >= 0 (sum 1) valid for all spectra with one coincidence exist iff a linear
    programme is feasible: after the substitution expressing the coincidence, the coefficient of
    every monomial must agree.  The programme is solved (HiGHS) for random U; coincidences tested:
    a_1 + b_1 = 0 (Theorem B.1), b_3 = b_4 and a_3 = a_4 (Lemma B.2); other index choices are
    equivalent by relabelling rows and columns.  Negative control: with no coincidence the
    programme asks for r(U) in the hull of the r(P_sigma), which fails for generic U."""
    import numpy as np
    from scipy.optimize import linprog
    ck.section("NUMERIC evidence for the CITED Theorem B.1 and Lemma B.2 (common weights; not a proof)")
    rng = np.random.default_rng(seed)
    subsets = [c for k in range(5) for c in itertools.combinations(range(4), k)]
    pairs = [(I, J) for I in subsets for J in subsets if len(I) == len(J)]       # 70 pairs

    def minor2(U, I, J):
        if not I:
            return 1.0
        return abs(np.linalg.det(U[np.ix_(I, J)])) ** 2

    def key_of(mode, I, J):
        A = tuple(i for i in range(4) if i not in I)       # a-variables present
        sign = 1
        if mode == "none":
            return (A, J), sign
        if mode == "B1":                                    # b_0 = -a_0
            e = (0 in A) + (0 in J)
            if 0 in J:
                sign = -1
            return (e, tuple(i for i in A if i), tuple(j for j in J if j)), sign
        if mode == "B2b":                                   # b_3 = b_2
            return (A, tuple(j for j in J if j < 2), sum(1 for j in J if j >= 2)), sign
        if mode == "B2a":                                   # a_3 = a_2
            return (tuple(i for i in A if i < 2), sum(1 for i in A if i >= 2), J), sign
        raise ValueError(mode)

    def feasible(U, mode):
        rows = {}
        for I, J in pairs:
            k, sg = key_of(mode, I, J)
            row = rows.setdefault(k, [np.zeros(24), 0.0])
            for s_idx, s in enumerate(PERMS):
                if set(s[i] for i in I) == set(J):
                    row[0][s_idx] += sg
            row[1] += sg * minor2(U, list(I), list(J))
        Aeq = np.array([r[0] for r in rows.values()])
        beq = np.array([r[1] for r in rows.values()])
        res = linprog(np.zeros(24), A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * 24, method="highs",
                      options={"primal_feasibility_tolerance": 1e-10})
        if res.status != 0:
            return False, None
        err = float(np.max(np.abs(Aeq @ res.x - beq)))
        worst_lp[0] = max(worst_lp[0], err)
        # independent of the monomial bookkeeping: the weights must reproduce the determinant at
        # random spectra having the coincidence
        S = np.array(PERMS)
        for _ in range(3):
            a = rng.normal(size=4) + 1j * rng.normal(size=4)
            b = rng.normal(size=4) + 1j * rng.normal(size=4)
            if mode == "B1":
                b[0] = -a[0]
            elif mode == "B2b":
                b[3] = b[2]
            elif mode == "B2a":
                a[3] = a[2]
            det = np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)
            z = np.prod(a[None, :] + b[S], axis=1)
            err = max(err, float(abs(det - res.x @ z)) / (1.0 + float(np.max(np.abs(z)))))
        return True, err

    worst_lp = [0.0]

    def haar(real):
        Z = rng.normal(size=(4, 4)) if real else rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        Q, Rm = np.linalg.qr(Z)
        Q = Q * (np.diag(Rm) / np.abs(np.diag(Rm)))
        if real and np.linalg.det(Q) < 0:
            Q[:, 0] = -Q[:, 0]
        return Q

    o_star = np.array([[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]], float) / 7.0
    mats = [o_star] + [haar(True) for _ in range(samples)] + [haar(False) for _ in range(samples)]
    for mode, name in (("B1", "Theorem B.1, a_1 + b_1 = 0"), ("B2b", "Lemma B.2, b_3 = b_4"),
                       ("B2a", "Lemma B.2 with the spectra interchanged, a_3 = a_4")):
        nok, worst = 0, 0.0
        worst_lp[0] = 0.0
        for U in mats:
            ok, r = feasible(U, mode)
            if ok:
                nok += 1
                worst = max(worst, r)
        ck.check("NUMERIC", "%s: common weights found by linear programming for all %d matrices (the "
                 "paper's O_* of Example 1.3, %d random in SO(4), %d Haar unitary; seed %d), and they "
                 "reproduce the determinant at 3 random spectra with the coincidence (tolerance 1e-7: "
                 "the weights come from a floating-point LP solver)"
                 % (name, len(mats), samples, samples, seed),
                 nok == len(mats) and worst < 1e-7,
                 "%d feasible; max LP equation residual %.2e; max relative determinant mismatch %.2e"
                 % (nok, worst_lp[0], worst))
    nfeas = sum(1 for U in mats if feasible(U, "none")[0])
    ck.info("negative control without a coincidence: the programme is feasible for %d of the same %d "
            "matrices, and it is %s at O_*" % (nfeas, len(mats),
                                               "feasible" if feasible(o_star, "none")[0] else "infeasible"))
