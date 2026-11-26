"""Shared exact objects for the independent replay of the group V cases T_27 and T_28 of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (October 2026).

Everything is rebuilt from the DEFINITIONS in the paper (Sections 3.1-3.2, 4.1-4.2 and
Appendices A, D, E; the tables (A.1), (A.2), Appendix D and Appendix E were read again from
the page images of pp. 35, 36, 59, 60 for this replay).  Nothing is imported from, executed
from, or read out of the author's supplement to build an object, and the audit library
`zheng-n4-audit/lib/moc4lib.py` is NOT used to build any object: it is only imported, when
present, by `crosscheck_with_audit_library` to compare finished objects.

Provenance of the code.  Sections 1-6 below (permutations, the group of order 1152, exact
linear algebra, bracket polynomials, the matrix P, the cubic, the table of the 82 vectors,
the enumeration of extreme rays, Leibniz sums and rational SO(4)) are COPIED, with small
edits, from this owner's earlier replays `../cases-i-iv/replay_iv_common.py` and
`../cases-i-iv/replay_T3_common.py`; the functions `replay_SOC_representative` and
`replay_four_coefficient_note` are copied from `replay_T3_common.py`.  What the two cases
rely on is re-checked here by `check_foundations` (every row of P is a verified polynomial
identity; grad C vanishes identically on the image of Psi; value profiles; the group
permutes the 82 vectors).  Sections 7-9 (polynomial arithmetic with python-flint, the
Hessian determinant, the radial polynomial, exact square roots, Sturm sequences, orbit and
covariance checks) are new.

Exact arithmetic: Python ints, fractions.Fraction, python-flint (fmpz_mpoly, fmpz_mat,
fmpq_poly) and sympy rings for two copied replays.  numpy/scipy only inside functions
whose names start with `numeric_`; their output is NUMERIC evidence and is labelled so.

Conventions
  * Permutations of {0,1,2,3} are tuples s with s[i] the image of i; PERMS is in
    lexicographic order, so PERMS[j] is the paper's pi_j.
  * Labels 1..8 carry vectors v_l = (x_l, y_l); [ij] = x_i y_j - x_j y_i.  With
    v_i = (a_i, 1), v_{4+j} = (-b_j, 1) one has [i, 4+j] = a_i + b_j.
  * V_R = Q^14 with coordinates m_0..m_13; values(c) = P c is the vector of the 24 numbers
    (P c)_sigma in PERMS order.
  * Group elements are triples g = (r1, r2, e): g.sigma = r1 o sigma o r2 for e = 0 and
    r1 o sigma^{-1} o r2 for e = 1; (g.c) is defined by v_{g.sigma}(g.c) = v_sigma(c).
"""
from fractions import Fraction
import hashlib
import itertools
import os
import sys
import time

# ----------------------------------------------------------------------------------------
# 0. check registry
# ----------------------------------------------------------------------------------------


class Checker:
    """Prints one line per check and keeps the count of failures."""

    def __init__(self, title):
        self.title = title
        self.n = 0
        self.failed = []
        self.counts = {}
        self.t0 = time.time()
        print("=" * 100)
        print(title)
        print("=" * 100)

    def section(self, name):
        print()
        print("--- " + name)
        sys.stdout.flush()

    def check(self, level, name, ok, detail=""):
        """level in {EXACT, NUMERIC, DATA}; ok must be a bool."""
        self.n += 1
        ok = bool(ok)
        self.counts[level] = self.counts.get(level, 0) + 1
        tag = "PASS" if ok else "FAIL"
        line = "[%s] [%s] %s" % (tag, level, name)
        if detail != "":
            line += "  ::  " + str(detail)
        print(line)
        sys.stdout.flush()
        if not ok:
            self.failed.append(name)
        return ok

    def info(self, text):
        print("       " + str(text))
        sys.stdout.flush()

    def finish(self):
        print()
        summary = ", ".join("%d %s" % (v, k) for k, v in sorted(self.counts.items()))
        if self.failed:
            print("FAIL: %d of %d checks failed (%s) in %.0f s: %s"
                  % (len(self.failed), self.n, summary, time.time() - self.t0, self.title))
            for f in self.failed:
                print("   failed: " + f)
            sys.exit(1)
        print("PASS: all %d checks passed (%s) in %.0f s: %s"
              % (self.n, summary, time.time() - self.t0, self.title))
        sys.exit(0)


# ----------------------------------------------------------------------------------------
# 1. permutations and the group G of order 1152
# ----------------------------------------------------------------------------------------
PERMS = [tuple(p) for p in itertools.permutations(range(4))]
PIDX = {p: i for i, p in enumerate(PERMS)}


def sgn(p):
    inv = sum(1 for i in range(len(p)) for j in range(i + 1, len(p)) if p[i] > p[j])
    return -1 if inv % 2 else 1


SGN = [sgn(p) for p in PERMS]


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


G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_perm(g, s):
    r1, r2, e = g
    s = inverse(s) if e else tuple(s)
    return compose(r1, compose(s, r2))


def g_index_map(g):
    """t with t[j] = index of g.pi_j."""
    return tuple(PIDX[g_perm(g, s)] for s in PERMS)


G_MAPS = [g_index_map(g) for g in G_ELEMENTS]


def g_values(gmap, vals):
    """Values of g.c, defined by v_{g.sigma}(g.c) = v_sigma(c)."""
    out = [None] * 24
    for j in range(24):
        out[gmap[j]] = vals[j]
    return tuple(out)


def g_set(gmap, T):
    return frozenset(gmap[j] for j in T)


def orbit_rep(T):
    """Lexicographically smallest sorted image of the index set T under G."""
    T = list(T)
    return min(tuple(sorted(t[j] for j in T)) for t in G_MAPS)


def g_matrix(g, O):
    """g.O with Phi_{g.c}(g.O) = (+-1) Phi_c(O): for e = 0, (g.O)[i][j] = O[r2(i)][r1^{-1}(j)];
    for e = 1 the matrix is transposed first."""
    r1, r2, e = g
    if e:
        O = [[O[j][i] for j in range(4)] for i in range(4)]
    r1i = inverse(r1)
    return [[O[r2[i]][r1i[j]] for j in range(4)] for i in range(4)]


# ----------------------------------------------------------------------------------------
# 2. exact linear algebra over Q
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
    from math import gcd
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
# 3. bracket polynomials in sixteen variables, the basis (A.1) and the matrix P
# ----------------------------------------------------------------------------------------
# A polynomial in x_1..x_8, y_1..y_8 is a dict {key: int}, key = sum of 8**slot over the
# factors (slot l-1 for x_l, slot 7+l for y_l); adding keys multiplies monomials.


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


# (A.1), page image of p. 35 of the paper.
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


def _build_P():
    monos = sorted({k for b in B_POLYS for k in b})
    rows = [[b.get(k, 0) for k in monos] for b in B_POLYS]          # 14 x #monomials
    _, piv = rref(rows)
    if len(piv) != 14:
        raise AssertionError("B_0..B_13 are not linearly independent")
    sel = [monos[c] for c in piv]
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
KAPPA = [0] * 13 + [1]


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
# 4. the cubic (A.2), page image of p. 36 of the paper
# ----------------------------------------------------------------------------------------


def cubic_C(m):
    (m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13) = m
    return (m13 * (m0 * m6 - m1 * m5)
            - m0 * m9 * (m0 + m1 + m2 + m3 + m5 + m7 + m9) - m1 * m5 * m9
            + m10 * (-(m0 * m3) - m0 * m9 + m1 * m2 + m1 * m7)
            + m11 * (-(m0 * m7) - m0 * m9 + m2 * m5 + m3 * m5)
            + m12 * (-(m0 * m9) + m1 * m7 - m2 * m6 - m2 * m9 + m3 * m5 + m3 * m7)
            + m2 * m10 * m11 - m3 * m8 * m10 - m4 * m7 * m11 + m4 * m8 * m9)


class _Mono:
    """Tiny exact polynomial type in m_0..m_13 (dict: sorted index tuple -> int), used only
    to expand the cubic into monomials so that gradient and Hessian are read off by rule."""

    def __init__(self, d):
        self.d = {k: v for k, v in d.items() if v}

    @staticmethod
    def var(i):
        return _Mono({(i,): 1})

    def __add__(self, o):
        out = dict(self.d)
        for k, v in o.d.items():
            out[k] = out.get(k, 0) + v
        return _Mono(out)

    def __neg__(self):
        return _Mono({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        out = {}
        for a, ca in self.d.items():
            for b, cb in o.d.items():
                k = tuple(sorted(a + b))
                out[k] = out.get(k, 0) + ca * cb
        return _Mono(out)


CUBIC_TERMS = sorted(cubic_C([_Mono.var(i) for i in range(14)]).d.items())   # ((i,j,k), coef)
assert all(len(k) == 3 for k, _ in CUBIC_TERMS)


def _hess_table():
    """table[a][b] = list of (coef, k): Hess C(m)[a][b] = sum coef * m_k."""
    tab = [[{} for _ in range(14)] for _ in range(14)]
    for idx, co in CUBIC_TERMS:
        e = [0] * 14
        for i in idx:
            e[i] += 1
        for a in range(14):
            if not e[a]:
                continue
            e1 = list(e)
            c1 = co * e1[a]
            e1[a] -= 1
            for b in range(14):
                if not e1[b]:
                    continue
                e2 = list(e1)
                c2 = c1 * e2[b]
                e2[b] -= 1
                k = e2.index(1)
                assert sum(e2) == 1
                tab[a][b][k] = tab[a][b].get(k, 0) + c2
    return [[sorted((v, k) for k, v in cell.items() if v) for cell in row] for row in tab]


HESS_TABLE = _hess_table()


def _grad_table():
    """table[a] = list of (coef, i, k): dC/dm_a = sum coef * m_i * m_k (i <= k)."""
    tab = [{} for _ in range(14)]
    for idx, co in CUBIC_TERMS:
        for a in set(idx):
            rest = list(idx)
            mult = rest.count(a)
            rest.remove(a)
            key = tuple(rest)
            tab[a][key] = tab[a].get(key, 0) + co * mult
    return [sorted((v, k[0], k[1]) for k, v in cell.items() if v) for cell in tab]


GRAD_TABLE = _grad_table()


def _acc(terms):
    """Sum of a list of ring elements (0 for the empty list); avoids int + element."""
    tot = None
    for t in terms:
        tot = t if tot is None else tot + t
    return 0 if tot is None else tot


def hess_at(u):
    """Hess C(u) as a 14 x 14 list (entries are linear in u); works over any ring
    (ints, Fractions, python-flint polynomials; scalars are multiplied from the right)."""
    return [[_acc([u[k] * cf for cf, k in HESS_TABLE[i][j]]) for j in range(14)] for i in range(14)]


def grad_at(u):
    return [_acc([u[i] * u[k] * cf for cf, i, k in GRAD_TABLE[a]]) for a in range(14)]


# ----------------------------------------------------------------------------------------
# 5. Gaussian rationals, matching vectors, the 82 vectors f_j
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


# Appendix D, page image of p. 59: ("n", I, J) or ("h", one-line pi, eps).
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
assert all(x.denominator == 1 for v in F_COORD for x in v)
F_COORD = [[int(x) for x in v] for v in F_COORD]
FIDX = {v: j for j, v in enumerate(F_VALUES)}
assert len(FIDX) == 82


def f_name(j):
    k, a, b = F_TABLE[j]
    return "n_{%s,%s}" % (a, b) if k == "n" else "h_{%s,%s}" % (a, "+" if b > 0 else "-")


def minor_spec(j):
    """(I, J) as 0-based tuples for a minor vector f_j."""
    k, a, b = F_TABLE[j]
    assert k == "n"
    return tuple(int(c) - 1 for c in a), tuple(int(c) - 1 for c in b)


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


def g_facet(gmap, j):
    """Index of the vector g.f_j."""
    return FIDX[g_values(gmap, F_VALUES[j])]


def find_group_elements(rep, target):
    """All indices k of G_ELEMENTS with g_k . f_{rep[m]} = f_{target[m]} for every position m."""
    out = []
    for k, gm in enumerate(G_MAPS):
        if all(g_facet(gm, a) == b for a, b in zip(rep, target)):
            out.append(k)
    return out


# ----------------------------------------------------------------------------------------
# 6. contact cone, Leibniz sums, rational SO(4)
# ----------------------------------------------------------------------------------------


def extreme_rays_of_contact_cone(T):
    """Extreme rays of N_T = {c : P c >= 0, (P c)_T = 0}, computed from the definition and
    without the list of 82 facets.  Returns (dimension of {(Pc)_T = 0}, sorted list of
    primitive integer value vectors of the extreme rays).  A nonzero c in N_T spans an
    extreme ray iff the inequalities tight at c have rank d - 1 on L = {(Pc)_T = 0}; all
    (d-1)-subsets of the remaining inequality rows are enumerated."""
    from math import gcd
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    N = nullspace([P[s] for s in T], 14)                     # basis of L, as rows
    d = len(N)
    Nint = []
    for v in N:
        den = 1
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
        Nint.append([int(x * den) for x in v])
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


def leibniz_weights(O):
    return [O[0][p[0]] * O[1][p[1]] * O[2][p[2]] * O[3][p[3]] * s for p, s in zip(PERMS, SGN)]


def phi(vals, O):
    """Phi(O) = sum_sigma sgn(sigma) vals[sigma] prod_i O[i][sigma(i)]."""
    return sum(v * w for v, w in zip(vals, leibniz_weights(O)) if v != 0)


def minor_det(O, I, J):
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


def det4(O):
    return sum(leibniz_weights(O))


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def random_rational_SO4(rng, size=4):
    """(I - S)(I + S)^{-1} for a random rational skew S: exactly special orthogonal."""
    while True:
        S = [[Fraction(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i + 1, 4):
                val = Fraction(rng.randint(-size, size), rng.randint(1, size))
                S[i][j], S[j][i] = val, -val
        A = [[Fraction(int(i == j)) - S[i][j] for j in range(4)] for i in range(4)]
        B = [[Fraction(int(i == j)) + S[i][j] for j in range(4)] for i in range(4)]
        Bt = [list(r) for r in zip(*B)]
        O = [solve_square(Bt, A[i]) for i in range(4)]       # X B = A row by row
        if all(x != 0 for row in O for x in row):
            return O


def is_SO4(O):
    T = [[O[j][i] for j in range(4)] for i in range(4)]
    return matmul(O, T) == [[Fraction(int(i == j)) for j in range(4)] for i in range(4)] and det4(O) == 1


def transfer_to_SO4(g, O):
    """For O in SO(4) return O'' in SO(4) with Phi_{g.c}(O'') = Phi_c(O) for all c:
    g.O, with the first column negated when det(g.O) = -1."""
    M = g_matrix(g, O)
    if det4(M) == -1:
        M = [[-row[0]] + list(row[1:]) for row in M]
    return M


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


# ----------------------------------------------------------------------------------------
# 7. sparse multivariate polynomials over Z with python-flint (fmpz_mpoly)
# ----------------------------------------------------------------------------------------
# python-flint 0.6 exposes fmpz_mpoly with +, -, *, **, //, divmod, gcd, len, coefficient(i),
# exponent_tuple(i) and integer evaluation; contexts are identified by the number of
# variables only.  All polynomials of one computation are created through one Ring.


class Ring:
    def __init__(self, names):
        from flint import fmpz_mpoly
        self._cls = fmpz_mpoly
        self.names = tuple(names)
        self.n = len(self.names)
        self.gens = tuple(fmpz_mpoly.gens(self.n))
        self.zero = fmpz_mpoly(0, self.n)
        self.one = fmpz_mpoly(1, self.n)
        self.var = dict(zip(self.names, self.gens))
        self._pow = {}

    def const(self, c):
        return self._cls(int(c), self.n)

    def idx(self, name):
        return self.names.index(name)

    def gpow(self, i, k):
        key = (i, k)
        if key not in self._pow:
            self._pow[key] = self.gens[i] ** k
        return self._pow[key]

    def monomial(self, exps, coef=1):
        m = self.const(coef)
        for i, k in enumerate(exps):
            if k:
                m = m * self.gpow(i, k)
        return m

    def to_dict(self, p):
        return {tuple(int(e) for e in p.exponent_tuple(i)): int(p.coefficient(i)) for i in range(len(p))}

    def from_dict(self, d):
        items = [self.monomial(e, c) for e, c in d.items() if c]
        return self.sum(items)

    def sum(self, items):
        """Balanced sum of a list of polynomials."""
        items = list(items)
        if not items:
            return self.zero
        while len(items) > 1:
            nxt = []
            for i in range(0, len(items) - 1, 2):
                nxt.append(items[i] + items[i + 1])
            if len(items) % 2:
                nxt.append(items[-1])
            items = nxt
        return items[0]

    def coeffs(self, p):
        return [int(p.coefficient(i)) for i in range(len(p))]

    def degrees(self, p):
        """Tuple of the degrees of p in each variable (0 for the zero polynomial)."""
        out = [0] * self.n
        for i in range(len(p)):
            e = p.exponent_tuple(i)
            for k in range(self.n):
                if e[k] > out[k]:
                    out[k] = int(e[k])
        return tuple(out)

    def total_degrees(self, p):
        """Set of total degrees of the monomials of p."""
        return {sum(int(x) for x in p.exponent_tuple(i)) for i in range(len(p))}

    def split(self, p, var):
        """dict k -> coefficient polynomial of gens[var]^k (not involving gens[var])."""
        parts = {}
        for i in range(len(p)):
            e = tuple(int(x) for x in p.exponent_tuple(i))
            k = e[var]
            parts.setdefault(k, {})[e[:var] + (0,) + e[var + 1:]] = int(p.coefficient(i))
        return {k: self.from_dict(d) for k, d in parts.items()}

    def evaluate_int(self, p, point):
        """Exact value of p at an integer point (python-flint's evaluation), as a Python int."""
        return int(p(*[int(x) for x in point]))

    def evaluate(self, p, point):
        """Exact value of p at a point given as a list of Fractions/ints (own term-by-term sum)."""
        tot = Fraction(0)
        pw = [{0: Fraction(1)} for _ in range(self.n)]
        for i in range(len(p)):
            e = p.exponent_tuple(i)
            t = Fraction(int(p.coefficient(i)))
            for k in range(self.n):
                ek = int(e[k])
                if ek:
                    cache = pw[k]
                    if ek not in cache:
                        cache[ek] = Fraction(point[k]) ** ek
                    t *= cache[ek]
            tot += t
        return tot

    def digest(self, p):
        """sha256 of a canonical text form (sorted exponent tuples with coefficients)."""
        d = self.to_dict(p)
        raw = "\n".join(",".join(map(str, e)) + ":" + str(c) for e, c in sorted(d.items()))
        return hashlib.sha256(raw.encode()).hexdigest()


class RF:
    """Rational function num/den over a Ring, without cancellation.  Equality is decided by
    cross-multiplication, so `x == y` is an exact polynomial identity check (denominators are
    products of the denominators that were entered and are never zero polynomials)."""
    __slots__ = ("R", "num", "den")

    def __init__(self, R, num, den=None):
        self.R = R
        self.num = R.const(num) if isinstance(num, int) else num
        self.den = R.one if den is None else (R.const(den) if isinstance(den, int) else den)
        assert len(self.den) > 0

    def _c(self, o):
        return o if isinstance(o, RF) else RF(self.R, o)

    def __add__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den + o.num * self.den, self.den * o.den)

    def __sub__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den - o.num * self.den, self.den * o.den)

    def __mul__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.num, self.den * o.den)

    def __truediv__(self, o):
        o = self._c(o)
        return RF(self.R, self.num * o.den, self.den * o.num)

    def __neg__(self):
        return RF(self.R, -self.num, self.den)

    def __eq__(self, o):
        o = self._c(o)
        return self.num * o.den == o.num * self.den

    def __hash__(self):
        raise TypeError("RF is not hashable")


def convert(src_dict, dst, mapping, fixed=None):
    """Polynomial in `dst` obtained from a term dict by sending source variable i to the
    destination variable index mapping[i]; mapping[i] = None requires exponent 0, unless
    `fixed` gives an integer value for it."""
    out = {}
    for e, c in src_dict.items():
        ne = [0] * dst.n
        cc = c
        for i, k in enumerate(e):
            if not k:
                continue
            j = mapping[i]
            if j is None:
                if fixed is None or i not in fixed:
                    raise ValueError("variable %d occurs but has no image" % i)
                cc = cc * fixed[i] ** k
            else:
                ne[j] += k
        ne = tuple(ne)
        out[ne] = out.get(ne, 0) + cc
    return dst.from_dict(out)


def subst_cleared(R, terms, nums, dens, degs):
    """sum over terms (e, c) of  c * prod_i nums[i]^e_i * dens[i]^(degs[i] - e_i)  in R:
    the polynomial obtained by substituting x_i = nums[i]/dens[i] and multiplying by
    prod_i dens[i]^degs[i].  `terms` is a dict {exponent tuple: int}; degs[i] >= max e_i."""
    nv = len(nums)
    cache = {}

    def factor(i, k):
        key = (i, k)
        if key not in cache:
            cache[key] = (nums[i] ** k) * (dens[i] ** (degs[i] - k))
        return cache[key]

    def rec(items, i):
        if i == nv:
            return R.const(sum(c for _, c in items))
        groups = {}
        for e, c in items:
            groups.setdefault(e[i], []).append((e, c))
        parts = []
        for k, sub in groups.items():
            assert 0 <= k <= degs[i]
            parts.append(rec(sub, i + 1) * factor(i, k))
        return R.sum(parts)

    return rec(list(terms.items()), 0)


def bareiss_det(M, R, log=None):
    """Determinant of a square matrix of fmpz_mpoly by fraction-free (Bareiss) elimination
    with the sparsest available pivot.  Every division is checked to be exact."""
    n = len(M)
    M = [list(row) for row in M]
    sign = 1
    prev = R.one
    t0 = time.time()
    for k in range(n - 1):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                L = len(M[i][j])
                if L and (best is None or L < best[0]):
                    best = (L, i, j)
        if best is None:
            return R.zero
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
                if k:
                    q, r = divmod(a, prev)
                    if len(r):
                        raise ArithmeticError("Bareiss division is not exact")
                    a = q
                M[i][j] = a
            M[i][k] = R.zero
        prev = piv
        if log is not None:
            log("elimination step %d: largest remaining entry has %d terms (%.1f s)"
                % (k, max(len(M[i][j]) for i in range(k + 1, n) for j in range(k + 1, n)), time.time() - t0))
    out = M[n - 1][n - 1]
    return out if sign > 0 else -out


def poly_sqrt(R, F, var=0):
    """Exact square root in Z[x] of F (python-flint polynomial of the Ring R), by descending
    coefficient comparison in gens[var], recursively; returns None if F is not a square.
    The result is verified by squaring."""
    if len(F) == 0:
        return R.zero
    if var == R.n:
        from math import isqrt
        c = int(F.coefficient(0))
        assert len(F) == 1
        if c < 0:
            return None
        s = isqrt(c)
        return R.const(s) if s * s == c else None
    co = R.split(F, var)
    D = max(co)
    if D % 2:
        return None
    m = D // 2
    g = [None] * (m + 1)
    g[m] = poly_sqrt(R, co[D], var + 1)
    if g[m] is None:
        return None
    two = g[m] * 2
    for i in range(m - 1, -1, -1):
        s = co.get(m + i, R.zero)
        for a in range(i + 1, m):
            b = m + i - a
            if i < b < m:
                s = s - g[a] * g[b]
        q, r = divmod(s, two)
        if len(r):
            return None
        g[i] = q
    G = R.sum([g[i] * R.gpow(var, i) if i else g[i] for i in range(m + 1)])
    return G if G * G == F else None


# ----------------------------------------------------------------------------------------
# 8. radial necessity on a normal form  g = sum_i r_i f_{order[i]}
# ----------------------------------------------------------------------------------------


def lift_coordinates(R, order, znames="z"):
    """The 14 coordinates of y = sum_i r_i f_{order[i]} + z kappa as linear forms in R
    (variables r_0.., then z)."""
    zi = R.idx(znames)
    y = []
    for k in range(14):
        terms = [R.gens[i] * F_COORD[j][k] for i, j in enumerate(order) if F_COORD[j][k]]
        if KAPPA[k]:
            terms.append(R.gens[zi] * KAPPA[k])
        y.append(R.sum(terms))
    return y


def radial_data(R, order, log=None):
    """For the normal form with generators `order` (variables R.names[:-1], last variable z):
    returns dict with
      C0, Q   : C(y0 + z kappa) = C0 + z Q  (checked to have degree <= 1 in z),
      det     : det Hess C(y0 + z kappa) as a polynomial in the coefficients and z,
      E       : sum_k det_k (-C0)^k Q^(4-k), i.e. Q^4 det Hess C(y0 - (C0/Q) kappa),
    all python-flint polynomials in R."""
    zi = R.n - 1
    y = lift_coordinates(R, order, R.names[zi])
    Cy = cubic_C(y)
    parts = R.split(Cy, zi)
    assert set(parts) <= {0, 1}, "the cubic is not affine in z on this lift"
    C0 = parts.get(0, R.zero)
    Q = parts.get(1, R.zero)
    H = hess_at(y)
    H = [[h if not isinstance(h, int) else R.const(h) for h in row] for row in H]
    det = bareiss_det(H, R, log=log)
    dz = R.split(det, zi)
    assert max(dz) <= 4
    E = R.sum([dz[k] * ((-C0) ** k) * (Q ** (4 - k)) for k in dz])
    return {"y": y, "C0": C0, "Q": Q, "H": H, "det": det, "det_by_z": dz, "E": E}


def det_spot_check(R, order, det, rng, trials=6, size=9):
    """Independent exact check of the symbolic determinant: at random integer points the value
    of `det` equals the determinant of the integer matrix Hess C(y) computed by python-flint's
    fmpz_mat (a different algorithm from the Bareiss elimination above)."""
    from flint import fmpz_mat
    ok = True
    for _ in range(trials):
        pt = [rng.randint(-size, size) for _ in range(R.n)]
        y = combo_coord([(j, pt[i]) for i, j in enumerate(order)])
        y[13] += pt[R.n - 1]
        Hn = hess_at(y)
        dn = int(fmpz_mat(Hn).det())
        ok = ok and Fraction(dn) == R.evaluate(det, pt)
    return ok


# ----------------------------------------------------------------------------------------
# 9. Sturm sequences (exact, rational)
# ----------------------------------------------------------------------------------------


def sturm_positive_roots(coeffs):
    """Number of distinct real roots in (0, +infinity) of the polynomial sum coeffs[i] x^i
    (integers, coeffs[0] != 0, coeffs[-1] != 0), by Sturm's theorem with the chain
    f, f', -rem(f, f'), ... computed with python-flint fmpq_poly division.  Returns
    (number of roots, chain length)."""
    from flint import fmpq_poly
    assert coeffs[0] != 0 and coeffs[-1] != 0
    f = fmpq_poly([int(c) for c in coeffs])
    chain = [f, f.derivative()]
    while chain[-1].degree() > 0:
        r = chain[-2] % chain[-1]
        if r == 0:
            break
        chain.append(-r)
    if chain[-1] == 0:
        chain.pop()

    def var(signs):
        s = [x for x in signs if x != 0]
        return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))

    def sg(x):
        x = Fraction(int(x.p), int(x.q))
        return (x > 0) - (x < 0)

    at0 = [sg(p[0]) for p in chain]
    atinf = [sg(p[p.degree()]) for p in chain]
    return var(at0) - var(atinf), len(chain)


def descartes_bound(coeffs):
    """Number of sign changes of the coefficient sequence (upper bound for positive roots)."""
    s = [c for c in coeffs if c != 0]
    return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))


# ----------------------------------------------------------------------------------------
# 10. common checks
# ----------------------------------------------------------------------------------------

# sha256 of the files that were read for this replay (data only; nothing is executed)
DOC_SHA256 = {
    "proofs/manin_transfer_20261021/final_referee_T26_uniform_20261022.pdf":
        "c0fd13f3a956e90f50c9aa89a228045bcd6f080ad72862320963c77bc857062a",
    "proofs/manin_transfer_20261021/final_referee_T26_boundary10_uniform_20261022.pdf":
        "ec3451eff2e7d8cb4f791b81a451cc5e1b7d10410a5c86f05094ac21cb99f953",
    "proofs/manin_transfer_20261021/final_referee_T26_boundary16_uniform_20261022.pdf":
        "0d8ee9453494f82c49027107c3cc806c47c26fe1cbbe9eb8969290ca3ca24cf4",
    "proofs/manin_transfer_20261021/blind_moc4_T27_uniform_20261022.pdf":
        "5004369d7570871904590e46cd07b7f9948d469881ddf931eec1273fcd8df722",
    "proofs/manin_transfer_20261021/final_referee_T26_full_radial_20261022.json":
        "a4ba38ff4e4a5783668ccb9c5689db272f82c2bbf824c8a03f8565f214ece7ab",
    "proofs/manin_transfer_20261021/blind_moc4_T27_full_radial_20261022.json":
        "2be64c16d2156f272dcf3f6e95b254ba98e216413b0574976af195a331f47894",
}

# Appendix E, page image of p. 60 of the paper (the four neighbours fix the shift of numbering)
PAPER_T = {25: (3, 4, 7, 8, 10, 12), 26: (2, 3, 5, 7, 8, 12), 27: (2, 3, 4, 7, 8, 12),
           28: (2, 3, 5, 6, 9, 12), 29: (0, 4, 5, 6, 7, 8)}


def supp_path(relpath):
    root = os.environ.get("MOC4_SUPP")
    if not root:
        return None
    path = os.path.join(root, relpath)
    return path if os.path.isfile(path) else None


def check_document_identity(ck, relpath):
    """Provenance only: if MOC4_SUPP is set, confirm that the file there is the one that was
    read for this replay.  The file is hashed, never executed."""
    if not os.environ.get("MOC4_SUPP"):
        ck.info("MOC4_SUPP not set: sha256 of %s not compared" % relpath)
        return
    path = supp_path(relpath)
    if path is None:
        ck.check("DATA", "file present under MOC4_SUPP: " + relpath, False)
        return
    with open(path, "rb") as fh:
        h = hashlib.sha256(fh.read()).hexdigest()
    ck.check("DATA", "file under MOC4_SUPP is the one that was read (sha256): " + os.path.basename(relpath),
             h == DOC_SHA256[relpath], h)


def print_versions(ck, numeric=True):
    import flint
    import sympy
    line = "python %s, python-flint %s, sympy %s" % (sys.version.split()[0], flint.__version__, sympy.__version__)
    if numeric:
        import numpy
        import scipy
        line += ", numpy %s, scipy %s (NUMERIC part only)" % (numpy.__version__, scipy.__version__)
    ck.info(line)


def check_foundations(ck):
    """The objects of Sections 3 and 4 of the paper that the cases rely on."""
    ck.section("Foundations rebuilt from the paper's definitions")
    ck.check("EXACT", "P is a 24 x 14 integer matrix; every row is a polynomial identity in "
             "16 variables (checked at construction); rank 14; P kappa = (1,...,1) for kappa = e_13",
             len(P) == 24 and all(len(r) == 14 for r in P) and rank(P) == 14
             and all(P[s][13] == 1 for s in range(24)))
    pts = [(0, 1)] * 4 + [(1, 1)] * 4
    ck.check("EXACT", "kappa = Psi(v) at v_1..4 = (0,1), v_5..8 = (1,1), so kappa is in X",
             [eval_poly16(B_POLYS[r], pts) for r in range(14)] == [0] * 13 + [1])
    bad = []
    BB = {}
    for j in range(14):
        tot = {}
        for coef, i, k in GRAD_TABLE[j]:
            key = (i, k)
            if key not in BB:
                BB[key] = pmul(B_POLYS[i], B_POLYS[k])
            tot = padd(tot, BB[key], coef)
        if tot:
            bad.append(j)
    ck.check("EXACT", "the cubic (A.2) has %d monomials; dC/dm_j (Psi(v)) = 0 identically in 16 variables for "
             "j = 0..13, so X lies in the singular locus of C (only this inclusion is used)" % len(CUBIC_TERMS),
             not bad and len(CUBIC_TERMS) == 28, "nonzero: %s" % bad if bad else "")
    prof_ok = True
    for j in range(82):
        v = F_VALUES[j]
        cnt = (v.count(0), v.count(1), v.count(2))
        if F_IS_MINOR[j]:
            prof_ok &= cnt in ((18, 6, 0), (20, 4, 0))
        else:
            prof_ok &= cnt == (13, 10, 1)
    ck.check("EXACT", "82 vectors f_j of Appendix D: all lie in the image of P with integer coordinates; 34 minor "
             "vectors with 18 or 20 zeros, 48 vectors h with value profile 13/10/1; pairwise distinct",
             prof_ok and len(set(F_VALUES)) == 82 and sum(F_IS_MINOR) == 34)
    ok = len(set(G_MAPS)) == 1152
    for t in G_MAPS:
        for j in range(82):
            jj = FIDX.get(g_values(t, F_VALUES[j]))
            if jj is None or F_IS_MINOR[jj] != F_IS_MINOR[j]:
                ok = False
    ck.check("EXACT", "group G: 1152 distinct maps on the 24 permutation indices; G permutes the 82 vectors f_j "
             "and maps minor vectors to minor vectors", ok)


def check_radial_necessity(ck, seed=20261127, samples=10):
    """PROOF.md 7.3 / the 'forced radial equation' paragraph of the documents:
    x = u + i v in X  =>  C(u) = 0, Hess C(u) v = 0, det Hess C(u) = 0."""
    import random
    from flint import fmpq_mat, fmpq
    ck.section("Radial necessity: the real part u of a point of X satisfies C(u) = 0 and det Hess C(u) = 0")
    R = Ring(["u%d" % i for i in range(14)] + ["v%d" % i for i in range(14)] + ["t"])
    u = list(R.gens[:14])
    v = list(R.gens[14:28])
    t = R.gens[28]
    uv = [u[i] + t * v[i] for i in range(14)]
    gu, gv, guv = grad_at(u), grad_at(v), grad_at(uv)
    Hu = hess_at(u)
    ok = True
    for j in range(14):
        hv = R.sum([Hu[j][k] * v[k] for k in range(14) if not isinstance(Hu[j][k], int)])
        ok &= (guv[j] - gu[j] - t * hv - t * t * gv[j]) == 0
    ck.check("EXACT", "grad C(u + t v) = grad C(u) + t Hess C(u) v + t^2 grad C(v) identically in 29 variables "
             "(t = i: real part grad C(u) - grad C(v), imaginary part Hess C(u) v)", ok)
    Cu = cubic_C(u)
    euler = R.sum([u[j] * gu[j] for j in range(14)]) - Cu * 3
    ck.check("EXACT", "Euler: u . grad C(u) = 3 C(u) identically", euler == 0)
    pol = R.sum([u[j] * gv[j] for j in range(14)]) * 2 - R.sum(
        [v[j] * Hu[j][k] * v[k] for j in range(14) for k in range(14) if not isinstance(Hu[j][k], int)])
    ck.check("EXACT", "polarization: 2 u . grad C(v) = v^t Hess C(u) v identically", pol == 0)
    hu = R.sum([Hu[0][k] * u[k] for k in range(14) if not isinstance(Hu[0][k], int)]) - gu[0] * 2
    ck.check("EXACT", "Hess C(u) u = 2 grad C(u) (first component, identically; used when Im x = 0)", hu == 0)
    ck.info("so grad C(u+iv) = 0 gives grad C(u) = grad C(v), Hess C(u) v = 0, 3 C(u) = u.grad C(v) = "
            "(1/2) v^t Hess C(u) v = 0; and Hess C(u) is singular (kernel vector v, or u if v = 0)")
    rng = random.Random(seed)
    okc = okh = okd = True
    for _ in range(samples):
        def rq():
            return GQ(Fraction(rng.randint(-9, 9), rng.randint(1, 5)),
                      Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
        a = [rq() for _ in range(4)]
        b = [rq() for _ in range(4)]
        al = rq()
        x = [GQ.of(c) for c in matching_vector(a, b, al)]
        uu = [c.re for c in x]
        vv = [c.im for c in x]
        okc &= cubic_C(uu) == 0
        H = hess_at(uu)
        okh &= all(sum(H[i][k] * vv[k] for k in range(14)) == 0 for i in range(14))
        Hq = fmpq_mat([[fmpq(Fraction(h).numerator, Fraction(h).denominator) for h in row] for row in H])
        okd &= Hq.det() == 0
        zs = z_points(a, b)
        okc &= [GQ.of(sum(P[s][r] * x[r] for r in range(14))) for s in range(24)] == [al * z for z in zs]
    ck.check("EXACT", "spot check on %d random Gaussian-rational points x = alpha Psi(a,b) (seed %d): "
             "(P x)_sigma = alpha z_sigma, C(Re x) = 0, Hess C(Re x) Im x = 0 and det Hess C(Re x) = 0"
             % (samples, seed), okc and okh and okd)


def check_contact_set(ck, nu, doc_T, doc_name):
    ck.section("Contact set and numbering")
    ck.check("EXACT", "the document's set %s (its label '%s') equals the paper's T_%d of Appendix E "
             "as a set of lexicographic indices (group element: identity)" % (sorted(doc_T), doc_name, nu),
             tuple(sorted(doc_T)) == PAPER_T[nu])
    reps = {k: orbit_rep(v) for k, v in PAPER_T.items()}
    ck.check("EXACT", "orbit representative under G of the document's set equals that of T_%d: %s; "
             "T_25 .. T_29 lie in five different orbits, so no other neighbouring case is meant"
             % (nu, reps[nu]), orbit_rep(doc_T) == reps[nu] and len(set(reps.values())) == len(reps))


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
    """idents: list of (left index list, right index list), claimed sum f_left = sum f_right."""
    bad = []
    for left, right in idents:
        if combo_values([(j, 1) for j in left]) != combo_values([(j, 1) for j in right]):
            bad.append((left, right))
    ck.check("EXACT", "%s: %s hold on all 24 values" % (label, "; ".join(
        " + ".join("f%d" % j for j in l) + " = " + " + ".join("f%d" % j for j in r) for l, r in idents)),
        not bad, "failing: %s" % bad if bad else "")


def check_symmetry(ck, T, gens, claimed_action, doc_element, seed):
    """The symmetry of the contact set used to exchange the two non-minor generators.
    claimed_action: dict j -> image index for all generators.  doc_element: (r1, r2, e) in
    this module's convention.  Returns the list of all group elements with that action."""
    import random
    hits = []
    for k, gm in enumerate(G_MAPS):
        if g_set(gm, T) != frozenset(T):
            continue
        if all(g_facet(gm, j) == claimed_action[j] for j in gens):
            hits.append(G_ELEMENTS[k])
    ck.check("EXACT", "some element of G preserves T and acts on the generators by %s"
             % ", ".join("%d->%d" % (j, claimed_action[j]) for j in gens), len(hits) > 0,
             "%d elements: %s" % (len(hits), hits))
    ck.check("EXACT", "the element printed in the document, read as sigma -> rho o sigma^{-1} o rho' with "
             "rho = %s, rho' = %s, is one of them" % (doc_element[0], doc_element[1]), doc_element in hits)
    if not hits:
        return hits
    g = doc_element if doc_element in hits else hits[0]
    rng = random.Random(seed)
    ok = True
    for _ in range(4):
        O = random_rational_SO4(rng)
        O2 = transfer_to_SO4(g, O)
        ok = ok and is_SO4(O) and is_SO4(O2)
        for j in range(82):
            ok = ok and phi(F_VALUES[j], O) == phi(F_VALUES[g_facet(g_index_map(g), j)], O2)
    ck.check("EXACT", "covariance for this element: Phi_{g.f_j}(O'') = Phi_{f_j}(O) for all 82 j on 4 random "
             "rational O in SO(4), with O'' = g.O (first column negated if det = -1) in SO(4)", ok)
    # g maps matching vectors to matching vectors of relabelled spectra
    r1, r2, e = g
    gm = g_index_map(g)
    ok = True
    for _ in range(4):
        def rq():
            return GQ(Fraction(rng.randint(-7, 7), rng.randint(1, 4)), Fraction(rng.randint(-7, 7), rng.randint(1, 4)))
        a = [rq() for _ in range(4)]
        b = [rq() for _ in range(4)]
        z = z_points(a, b)
        r1i = inverse(r1)
        if e == 0:
            a2 = [a[r2[i]] for i in range(4)]
            b2 = [b[r1i[j]] for j in range(4)]
        else:
            a2 = [b[r2[i]] for i in range(4)]
            b2 = [a[r1i[j]] for j in range(4)]
        z2 = z_points(a2, b2)
        ok = ok and all(z2[gm[s]] == z[s] for s in range(24))
    ck.check("EXACT", "the element maps X to X: the permuted products are the products of relabelled spectra "
             "(a' = a o rho', b' = b o rho^{-1}, with a and b exchanged first under transposition), exact on 4 "
             "random Gaussian-rational spectra", ok)
    return hits


def orbit_check(ck, family, rep, target, sym_positions=None, expected_orbit=None):
    """`target` (tuple of facet indices with roles) lies in the G-orbit of `rep`; the roles in
    `sym_positions` may be exchanged (the inequality is symmetric in them).  Returns a group
    element (r1, r2, e) with g.rep = target (possibly after the exchange) or None."""
    reps = [tuple(rep)]
    if sym_positions is not None:
        i, j = sym_positions
        sw = list(rep)
        sw[i], sw[j] = sw[j], sw[i]
        reps.append(tuple(sw))
    hit = None
    for rp in reps:
        ks = find_group_elements(rp, target)
        if ks:
            hit = (G_ELEMENTS[ks[0]], rp, len(ks))
            break
    ck.check("EXACT", "%s: roles %s lie in the G-orbit of the representative %s"
             % (family, tuple(target), tuple(rep)), hit is not None,
             "" if hit is None else "g = (rho=%s, rho'=%s, transpose=%d), %d group elements, source roles %s"
             % (hit[0][0], hit[0][1], hit[0][2], hit[2], hit[1]))
    if expected_orbit is not None:
        orb = {tuple(g_facet(gm, a) for a in rep) for gm in G_MAPS}
        if sym_positions is not None:
            i, j = sym_positions
            orb = {tuple(sorted((t[i], t[j]))) + tuple(t[m] for m in range(len(t)) if m not in (i, j))
                   for t in orb}
        ck.check("EXACT", "%s: orbit size = %d" % (family, expected_orbit), len(orb) == expected_orbit, len(orb))
    return None if hit is None else (hit[0], hit[1])


def covariance_check(ck, label, g, pairs, rng, trials=6):
    """Exact: for random rational O in SO(4), O'' = transfer_to_SO4(g, O) is in SO(4) and
    Phi_{f_target}(O'') = Phi_{f_rep}(O) for each (rep, target) in pairs."""
    ok = True
    for _ in range(trials):
        O = random_rational_SO4(rng)
        O2 = transfer_to_SO4(g, O)
        ok = ok and is_SO4(O) and is_SO4(O2)
        for a, b in pairs:
            ok = ok and phi(F_VALUES[a], O) == phi(F_VALUES[b], O2)
    ck.check("EXACT", "%s: Phi_{g.f}(O'') = Phi_f(O) with O'' in SO(4), exact on %d random rational O"
             % (label, trials), ok)


def check_minor_leibniz(ck, minor_indices, seed=5):
    """Lemma 6.2 / (4.2) for the listed minor vectors, through two polynomial identities in a
    generic 4 x 4 matrix M (sixteen indeterminates):
      (L) sum_{sigma(I)=J} sgn(sigma) prod_i M[i,sigma(i)] = eps det M[I,J] det M[I^c,J^c],
      (J) det(adj M)[J,I] = eps det M[I^c,J^c] det(M)^(|I|-1),      eps = (-1)^(sum I + sum J).
    For O in SO(4), adj O = O^t, so (J) gives eps det O[I^c,J^c] = det O^t[J,I] = det O[I,J]
    and (L) becomes Phi_{n_{I,J}}(O) = det O[I,J]^2 >= 0.  No parametrisation of SO(4) is used."""
    import random
    R = Ring(["o%d%d" % (i + 1, j + 1) for i in range(4) for j in range(4)])
    M = [[R.gens[4 * i + j] for j in range(4)] for i in range(4)]

    def det(rows, cols):
        n = len(rows)
        if n == 1:
            return M[rows[0]][cols[0]]
        tot = None
        for k in range(n):
            t = M[rows[0]][cols[k]] * det(rows[1:], cols[:k] + cols[k + 1:])
            t = t if k % 2 == 0 else -t
            tot = t if tot is None else tot + t
        return tot

    full = (0, 1, 2, 3)
    detM = det(full, full)
    adj = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            cof = det(tuple(x for x in full if x != i), tuple(x for x in full if x != j))
            adj[j][i] = cof if (i + j) % 2 == 0 else -cof
    ok = True
    for jx in minor_indices:
        I, J = minor_spec(jx)
        Ic = tuple(x for x in full if x not in I)
        Jc = tuple(x for x in full if x not in J)
        eps = (-1) ** (sum(I) + sum(J))
        lhs = _acc([w * v for v, w in zip(F_VALUES[jx], leibniz_weights(M)) if v])
        dIJ = det(I, tuple(sorted(J)))
        dc = det(Ic, Jc)
        ok = ok and lhs == dIJ * dc * eps
        Js = tuple(sorted(J))
        if len(I) == 1:
            dadj = adj[Js[0]][I[0]]
            ok = ok and dadj == dc * eps
        else:
            dadj = adj[Js[0]][I[0]] * adj[Js[1]][I[1]] - adj[Js[0]][I[1]] * adj[Js[1]][I[0]]
            ok = ok and dadj == dc * detM * eps
    rng = random.Random(seed)
    for _ in range(3):
        O = random_rational_SO4(rng)
        for jx in minor_indices:
            I, J = minor_spec(jx)
            ok = ok and phi(F_VALUES[jx], O) == minor_det(O, I, tuple(sorted(J))) ** 2
    ck.check("EXACT", "minor generators %s: Laplace and Jacobi identities in a generic 4 x 4 matrix, which give "
             "Phi_{f_j}(O) = (minor of O)^2 >= 0 on SO(4) (adj O = O^t); also checked on 3 random rational O"
             % list(minor_indices), ok)


# ----------------------------------------------------------------------------------------
# 11. replays of the SOC and four-coefficient representatives (copied from replay_T3_common)
# ----------------------------------------------------------------------------------------


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
             phi(F_VALUES[1], M) == x1 ** 2 and phi(F_VALUES[9], M) == x9 ** 2
             and x1 == vh[0] + vh[1] and x9 == vh[0] + vh[3])
    rhs = (2 * (a[0] * A[2] + a[2] * A[0] + a[1] * A[3] + a[3] * A[1]) * pp * qq
           - 4 * p0 * p1 * p2 * p3 * qq ** 2 - 4 * q0 * q1 * q2 * q3 * pp ** 2)
    ck.check("EXACT", "Phi_{f22} = 2(a0A2 + a2A0 + a1A3 + a3A1) - 4 p0p1p2p3 - 4 q0q1q2q3 on S^3 x S^3 (homogenised)",
             phi(F_VALUES[22], M) == rhs)
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
            "audit, foundations/so4-atoms.md 3.1; only its polynomial identities are replayed here")


def replay_four_coefficient_note(ck):
    """The four-coefficient inequality for its representative (f0, f6; f10; f54): the bridge to
    Leibniz sums in the model x -> p x q, the certificate identity modulo N = |p|^2 - |q|^2,
    and the determinant identities behind 'Z is PSD iff capacity'."""
    from sympy.polys.rings import ring
    from sympy.polys.domains import ZZ
    ck.section("Four-coefficient representative (f0, f6; f10; f54): identities of the certificate, own code")
    R, p0, p1, p2, p3, q0, q1, q2, q3, a, b, c, d = ring("p0,p1,p2,p3,q0,q1,q2,q3,a,b,c,d", ZZ)
    Pq, Qq = (p0, p1, p2, p3), (q0, q1, q2, q3)
    pp, qq = sum(x * x for x in Pq), sum(x * x for x in Qq)
    Ssub = qq - p1 ** 2 - p2 ** 2 - p3 ** 2          # p0^2 = Ssub  modulo N

    def red(f):
        """Normal form modulo N = |p|^2 - |q|^2 (N is monic in p0): degree in p0 <= 1."""
        out = R(0)
        for mon, cf in f.terms():
            e0 = mon[0]
            out += R({(e0 % 2,) + tuple(mon[1:]): cf}) * Ssub ** (e0 // 2)
        return out

    def quads(y):
        y0, y1, y2, y3, z0, z1, z2, z3 = y
        A = y0 * z3 + y3 * z0
        B = y1 * z2 - y2 * z1
        D = y0 ** 2 + y2 ** 2 - z0 ** 2 - z2 ** 2
        E = y1 ** 2 - y2 ** 2 - z0 ** 2 + z3 ** 2
        G = y0 ** 2 - y3 ** 2 - z1 ** 2 + z2 ** 2
        return A, B, D, E, G

    y = (p0, p1, p2, p3, q0, q1, q2, q3)
    A, B, D, E, G = quads(y)
    H = 2 * (A ** 2 + B ** 2) - E * G
    n = pp + qq
    M = quat_matrix(Pq, Qq, R(1), R(0), conj_q=False)
    quart = {0: (A - B) ** 2, 6: (A + B) ** 2, 10: D ** 2, 54: H}
    ck.check("EXACT", "bridge: 4 Phi_{f_j}(M(p,q)) = n^2 Q_j mod N for j = 0, 6, 10, 54 with "
             "Q = (A-B)^2, (A+B)^2, D^2, H and n = |p|^2 + |q|^2",
             all(red(4 * phi(F_VALUES[j], M)) == red(n * n * quart[j]) for j in (0, 6, 10, 54)))
    Z = [[a + b + 2 * d, a - b, R(0), R(0)],
         [a - b, a + b + 2 * d, R(0), 2 * d],
         [R(0), R(0), c, d],
         [R(0), 2 * d, d, d]]

    def lm(yy):
        y0, y1, y2, y3, z0, z1, z2, z3 = yy
        A_, B_, D_, E_, G_ = quads(yy)
        l = [-y0 * B_, y0 * A_, -z3 * D_, -z3 * G_]
        m = [z0 * B_, -z0 * A_, -y3 * D_, -y3 * E_]
        return l, m

    def qform(vec):
        return sum(vec[i] * Z[i][j] * vec[j] for i in range(4) for j in range(4))

    Smaps = [y,
             (-p1, p0, -p3, p2, q1, -q0, -q3, q2),
             (-p2, p3, p0, -p1, q2, q3, -q0, -q1),
             (-p3, -p2, p1, p0, q3, -q2, q1, -q0)]
    rhs = R(0)
    for Sy in Smaps:
        l, m = lm(Sy)
        rhs += qform(l) + qform(m)
    lhs = n * (a * (A - B) ** 2 + b * (A + B) ** 2 + c * D ** 2 + d * H)
    ck.check("EXACT", "certificate: n[a(A-B)^2 + b(A+B)^2 + cD^2 + dH] = sum_j l(S_j y)^t Z l(S_j y) + "
             "m(S_j y)^t Z m(S_j y) mod N, with each S_j preserving |p|^2 and |q|^2",
             red(lhs) == red(rhs)
             and all(sum(t * t for t in Sy[:4]) == pp and sum(t * t for t in Sy[4:]) == qq for Sy in Smaps))
    cap = a * b * (c - d) - d * d * (a + b + c + d)
    S3 = [[a + b + 2 * d, a - b, R(0)], [a - b, a + b - 2 * d, -2 * d], [R(0), -2 * d, c - d]]
    det2 = (a + b + 2 * d) * ((a + b - 2 * d) * (c - d) - 4 * d * d) - (a - b) ** 2 * (c - d)
    ck.check("EXACT", "Schur complements: d (Z[:3,:3] - col col^t / d) = d S3, and (c-d) det of the 2x2 "
             "complement = 4 [ab(c-d) - d^2(a+b+c+d)]",
             all(d * Z[i][j] - Z[i][3] * Z[3][j] == d * S3[i][j] for i in range(3) for j in range(3))
             and det2 == 4 * cap)
    ck.info("the inequality itself (a f0 + b f6 + c f10 + d f54 >= 0 on SO(4) when ab(c-d) >= d^2(a+b+c+d)) is "
            "CITED: supplement note final_referee_T5_four_budget_capacity and audit, foundations/so4-atoms.md 3.4")


# ----------------------------------------------------------------------------------------
# 12. comparisons with the audit library and with the supplement's stored objects
# ----------------------------------------------------------------------------------------


def crosscheck_with_audit_library(ck, cases):
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
             L.PERMS == PERMS and L.P == P and [int(x) for x in L.KAPPA] == KAPPA)
    ck.check("EXACT", "the 82 value vectors, the minor flags and the coordinates of all f_j agree with moc4lib",
             [tuple(v) for v in L.F_VALUES] == F_VALUES and list(L.F_IS_MINOR) == F_IS_MINOR
             and [[int(x) for x in v] for v in L.F_VECS] == F_COORD)
    import random
    rng = random.Random(77)
    okc = True
    for _ in range(20):
        m = [Fraction(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(14)]
        okc = okc and L.cubic_C(m) == cubic_C(m)
    ck.check("EXACT", "the cubic agrees with moc4lib.cubic_C at 20 random rational points (both are transcriptions "
             "of (A.2); a cubic in 14 variables is not determined by 20 values, this is a consistency test)", okc)
    ck.check("EXACT", "T_%s agree with moc4lib.CONTACT_SETS, with the same generators; the group has the same "
             "1152 maps" % ", T_".join(str(k) for k in cases),
             all(tuple(sorted(L.CONTACT_SETS[k])) == PAPER_T[k]
                 and list(L.generators(PAPER_T[k])) == generators(PAPER_T[k]) for k in cases)
             and set(L._G_MAPS) == set(G_MAPS))


def crosscheck_with_supplement_data(ck):
    """If MOC4_SUPP is set: read two JSON data files (as data; nothing is executed) and compare
    the stored P, cubic and facet table with the objects rebuilt here.  This only shows that the
    documents talk about the same objects; it is not evidence for any step."""
    import json
    if not os.environ.get("MOC4_SUPP"):
        ck.info("MOC4_SUPP not set: the supplement's stored P, cubic and facet table are not compared")
        return
    ck.section("Supplement data files (read as JSON, untrusted, comparison only)")
    f1 = supp_path("proofs/manin_transfer_20261021/root_scout_moc_spinor_20261021.json")
    f2 = supp_path("proofs/principal_minor_result_20261022/root_permutation_polytope_complete_20261022.json")
    if f1 is None or f2 is None:
        ck.info("data files not found under MOC4_SUPP: comparison skipped")
        return
    with open(f1) as fh:
        meta = json.load(fh)
    with open(f2) as fh:
        fac = json.load(fh)
    mine = {k: v for k, v in CUBIC_TERMS}
    theirs = {}
    for ix, co in meta["cubic_terms"]:
        key = tuple(sorted(int(j) for j in ix))
        theirs[key] = theirs.get(key, 0) + Fraction(str(co))
    ck.check("DATA", "the supplement's stored 24 x 14 matrix equals P, and its stored cubic equals the paper's C "
             "of (A.2) (same sign): the documents' 'matching cubic F' is C",
             meta["permutation_to_basis"] == P and {k: Fraction(v) for k, v in mine.items()} == theirs)
    ck.check("DATA", "the supplement's stored 82 facet value vectors equal F_VALUES in the same order "
             "(the documents' f_j are the paper's f_j of Appendix D)",
             [tuple(z["values"]) for z in fac["facets"]] == F_VALUES)


def load_radial_json(relpath):
    """The stored radial factorisation of the supplement, as data: returns dict with the
    variable names, C_terms, Q_terms (dicts) and the list of (degree, multiplicity, dict)."""
    import json
    path = supp_path(relpath)
    if path is None:
        return None
    with open(path) as fh:
        d = json.load(fh)
    td = lambda terms: {tuple(int(x) for x in e): int(c) for e, c in terms}
    return {"variables": list(d["variables"]), "scale": d["scale"], "constant": d["constant"],
            "C": td(d["C_terms"]), "Q": td(d["Q_terms"]),
            "factors": [(int(f["degree"]), int(f["multiplicity"]), td(f["terms"])) for f in d["radial_factors"]]}


# ----------------------------------------------------------------------------------------
# 13. NUMERIC helpers (floating point; evidence only)
# ----------------------------------------------------------------------------------------


def numeric_stratum_search(T, seed, starts, eps=0.02, tol=1e-10):
    """NUMERIC.  Search for spectra a, b in C^4 and a phase alpha such that all gaps
    Re(alpha z_sigma) - Re(alpha z_pi*) are >= 0, vanish on T, and sum to 1 (pi* = first index
    of T).  Two Levenberg-Marquardt stages from each random start (17 standard normal
    parameters): stage 1 adds residuals eps / (s_kl + 1e-3) that push the eight labels a_i, -b_j
    apart; stage 2 drops them and polishes, so a returned point solves the stated equations to
    `tol` in the maximum norm.  Returns a list of dicts with the gap vector `gaps`, the real
    part `u` of x in the 14 matching coordinates, `h` = Re(alpha z_pi*), and `sep` = the
    smallest normalised label distance."""
    import numpy as np
    from scipy.optimize import least_squares
    T = sorted(T)
    others = [s for s in range(24) if s not in T]
    S = np.array(PERMS)
    rng = np.random.default_rng(seed)
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
        out.append(dict(h=float(F[T[0]]), sep=float(np.min(seps(a, b))), gaps=g, F=F,
                        contacts=tuple(int(s) for s in range(24) if abs(g[s]) < 1e-7)))
    return out


def numeric_cone_coefficients(gaps, gens):
    """NUMERIC.  Nonnegative least squares coefficients of a gap vector on the generators."""
    import numpy as np
    from scipy.optimize import nnls, lsq_linear
    A = np.array([F_VALUES[j] for j in gens], float).T
    g = np.maximum(np.asarray(gaps, float), 0.0)
    try:
        lam, res = nnls(A, g, maxiter=200 * A.shape[1])
    except RuntimeError:                      # scipy's active-set iteration can cycle on degenerate data
        sol = lsq_linear(A, g, bounds=(0.0, np.inf), tol=1e-14)
        lam, res = sol.x, float(np.linalg.norm(A @ sol.x - g))
    return dict(zip(gens, lam)), float(res)


def numeric_min_phi(vals, seed, starts=24):
    """NUMERIC.  Approximate minimum over SO(4) of the Leibniz sum with value vector `vals`,
    by BFGS from random starts in the model x -> p x conj(q) (p, q normalised).  Returns an
    upper bound for the minimum (not certified)."""
    import numpy as np
    from scipy.optimize import minimize
    S = np.array(PERMS)
    w = np.array(SGN, float) * np.asarray(vals, float)
    rng = np.random.default_rng(seed)

    def mat(x):
        p = x[:4] / np.linalg.norm(x[:4])
        q = x[4:] / np.linalg.norm(x[4:])
        M = quat_matrix(tuple(p), tuple(q), 1.0, 0.0, conj_q=True)
        return np.array(M, float)

    def f(x):
        O = mat(x)
        return float(np.sum(w * O[0, S[:, 0]] * O[1, S[:, 1]] * O[2, S[:, 2]] * O[3, S[:, 3]]))

    best = np.inf
    for _ in range(starts):
        r = minimize(f, rng.normal(size=8), method="BFGS", options={"gtol": 1e-12, "maxiter": 400})
        best = min(best, float(r.fun))
    return best
