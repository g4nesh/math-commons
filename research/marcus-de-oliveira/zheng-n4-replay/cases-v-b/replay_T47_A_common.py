"""Shared exact tools for the independent replay of PART A of the paper's case
T_47 = {1,3,5,6,7,8} of D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (October 2026).

Proof documents replayed (author's supplement, read as PDF; never executed).  The documents
call the contact set "T46"; it is the paper's T_47.
  root_moc4_T46_whole_assembly_20261023.pdf                    (top-level assembly)
  final_referee_T46_six_ray_coherent_capacity_20261023.pdf     (six-ray Gram capacity)
  final_referee_T46_core_coherent_capacity_20261023.pdf        (five-ray Gram capacity)
  final_referee_T46_negative_j_Heron_uniform_20261023.pdf      (necessity H6 >= 0)
  final_referee_T46_two_Heron_necessity_20261023.pdf           (necessity H8 >= 0)

Everything is rebuilt from the DEFINITIONS in the paper (Sections 3.1-3.2, 4.1-4.2, Appendices
A, D, E).  Nothing of the supplement is imported or executed.  Supplement files are read only
as DATA, only through the environment variable MOC4_SUPP, and only for
  * sha256 comparison of the files that were read,
  * comparison of stored objects (P, cubic, facet table, the polynomial G, the polynomial W,
    the six stored group actions of the closed faces) with the objects rebuilt here,
  * the sum-of-squares certificate behind the Gram capacities: sixteen cubic polynomials
    ('features'), four signed coordinate permutations and an 8 x 8 change of basis exist only
    in JSON files of the supplement, not in the PDFs.  They are parsed here as text; the
    identity they are claimed to satisfy is verified against Leibniz sums rebuilt from the
    paper's definitions, and the linear relations on the Gram matrix are RE-DERIVED here by
    exact linear algebra (the stored relations are only compared).

Provenance of the code.  This module is a COPY of this owner's earlier helper
../cases-v-a/replay_T25_common.py (sections 0-14: Checker, permutations, the group G of order
1152 with explicit elements, exact linear algebra, the matrix P, the cubic (A.2), Gaussian
rationals, the 82 facet vectors of Appendix D, extreme rays of a contact cone, Leibniz sums,
the quaternion model, flint polynomial ring, the class Radial, the class GramFace, common
checks, numeric searches, symbolic Gram checks), with the file hashes and the table of contact
sets replaced.  The copied checks are re-run by the scripts that use them.  Section 15 is new
(recursive exact square root, the documents' convention for group actions, tangent spaces of
X, exact inertia, a second numeric search copied from ../cases-v-a/replay_T29_common.py).

Exact arithmetic: Python ints / fractions.Fraction, sympy (identities in few variables) and
python-flint 0.6 (fmpz_mpoly: products, exact division, gcd of large integer polynomials).
numpy / scipy are imported only inside functions whose names start with `numeric_`.
"""
from fractions import Fraction
import hashlib
import itertools
import json
import os
import re
import sys
import time

# ----------------------------------------------------------------------------------------
# 0. check registry                                   (copied from replay_iv_common.py)
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
        print("--- " + name + "   [t = %.0f s]" % (time.time() - self.t0))
        sys.stdout.flush()

    def check(self, level, name, ok, detail=""):
        """level in {EXACT, NUMERIC}; ok must be a bool."""
        self.n += 1
        ok = bool(ok)
        tag = "PASS" if ok else "FAIL"
        line = "[%s] [%s] %s" % (tag, level, name)
        if detail != "":
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
# 1. permutations and the group G of order 1152
#    (replay_iv_common.py; the explicit elements and the matrix action from replay_T3_common.py)
# ----------------------------------------------------------------------------------------
PERMS = [tuple(p) for p in itertools.permutations(range(4))]      # lexicographic: PERMS[j] = pi_j
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


# g = (r1, r2, e):  g.sigma = r1 o sigma o r2 (e = 0),  r1 o sigma^{-1} o r2 (e = 1).
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_perm(g, s):
    r1, r2, e = g
    s = inverse(s) if e else tuple(s)
    return compose(r1, compose(s, r2))


def g_index_map(g):
    return tuple(PIDX[g_perm(g, s)] for s in PERMS)


G_MAPS = [g_index_map(g) for g in G_ELEMENTS]


def g_values(gmap, vals):
    """(g.c) defined by v_{g.sigma}(g.c) = v_sigma(c)."""
    out = [None] * 24
    for j in range(24):
        out[gmap[j]] = vals[j]
    return tuple(out)


def g_set(gmap, T):
    return tuple(sorted(gmap[j] for j in T))


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


def g_spectra(g, a, b):
    """Spectra (a', b') whose 24 products are the g-image of those of (a, b):
    z_{g.sigma}(a', b') = z_sigma(a, b).  For e = 0: a' = a o r2, b' = b o r1^{-1};
    for e = 1 the two spectra are exchanged first."""
    r1, r2, e = g
    if e:
        a, b = b, a
    r1i = inverse(r1)
    return [a[r2[i]] for i in range(4)], [b[r1i[j]] for j in range(4)]


# ----------------------------------------------------------------------------------------
# 2. exact linear algebra over Q                      (copied from replay_iv_common.py)
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


def mat_inv(M):
    n = len(M)
    R, piv = rref([list(M[i]) + [Fraction(int(i == j)) for j in range(n)] for i in range(n)])
    if piv[:n] != list(range(n)):
        raise ValueError("singular matrix")
    return [row[n:] for row in R]


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
#                                                     (copied from replay_iv_common.py)
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
KAPPA = [Fraction(0)] * 13 + [Fraction(1)]


def values(c):
    return [sum(P[s][a] * c[a] for a in range(14)) for s in range(24)]


_PIV = rref([list(col) for col in zip(*P)])[1]      # independent rows of P
assert len(_PIV) == 14


def solve_P(vals):
    """The unique c in Q^14 with P c = vals; raises unless all 24 equations hold."""
    c = solve_square([P[r] for r in _PIV], [Fraction(vals[r]) for r in _PIV])
    if values(c) != [Fraction(v) for v in vals]:
        raise ValueError("value vector is not in the image of P")
    return c


# ----------------------------------------------------------------------------------------
# 4. the cubic (A.2), read from the page image of p. 36 of the paper
#                                                     (copied from replay_iv_common.py)
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


def hess_linear_table():
    """table[i][j] = list of (integer coefficient, k) with (Hess C(u))_ij = sum coefficient u_k."""
    import sympy
    if "hess_lin" not in _CUBIC_CACHE:
        ms, _, _, hess = cubic_data()
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
    return _CUBIC_CACHE["hess_lin"]


def hess_at(u):
    """Hess C(u) as a 14 x 14 list of exact numbers (entries are linear in u)."""
    table = hess_linear_table()
    return [[sum((cf * u[k] for cf, k in table[i][j]), Fraction(0)) for j in range(14)]
            for i in range(14)]


# ----------------------------------------------------------------------------------------
# 5. Gaussian rationals and matching vectors          (copied from replay_iv_common.py)
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
# 6. the 82 vectors f_j of Appendix D (through their values on the 24 permutations)
#                                                     (copied from replay_iv_common.py)
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
F_INDEX = {v: j for j, v in enumerate(F_VALUES)}


def f_name(j):
    k, a, b = F_TABLE[j]
    return "n_{%s,%s}" % (a, b) if k == "n" else "h_{%s,%s}" % (a, "+" if b > 0 else "-")


def minor_spec(j):
    """(I, J) as 0-based tuples for a minor facet f_j."""
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
    """Index of the facet g.f_j (KeyError if G did not permute the facets)."""
    return F_INDEX[g_values(gmap, F_VALUES[j])]


def find_group_elements(rep, target):
    """All indices k of G with g_k . f_{rep[m]} = f_{target[m]} for every position m."""
    out = []
    for k, gm in enumerate(G_MAPS):
        if all(g_facet(gm, a) == b for a, b in zip(rep, target)):
            out.append(k)
    return out


def extreme_rays_of_contact_cone(T):
    """Extreme rays of N_T = {c : P c >= 0, (P c)_T = 0}, computed from the definition and
    without the list of 82 facets.  Returns (dim of the linear span of the equalities,
    sorted list of primitive integer value vectors (24-tuples) of the extreme rays).
    Every extreme ray is the kernel of some d - 1 of the 24 - |T| remaining inequality rows
    restricted to L = {(Pc)_T = 0}; all such subsets are enumerated."""
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


# ----------------------------------------------------------------------------------------
# 7. Leibniz sums, minors, exact rational SO(4), quaternion model
#                                                     (from replay_T3_common.py)
# ----------------------------------------------------------------------------------------


def leibniz(vals, O, zero=0):
    """Phi(O) = sum_sigma sgn(sigma) vals[sigma] prod_i O[i][sigma(i)]  (entries first, then the
    integer factor: python-flint 0.6 polynomials only multiply by an int on the right)."""
    tot = zero
    for pm, sg, v in zip(PERMS, SGN, vals):
        if v:
            tot = tot + O[0][pm[0]] * O[1][pm[1]] * O[2][pm[2]] * O[3][pm[3]] * (sg * v)
    return tot


def minor_det(O, I, J):
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


def det4(O):
    return leibniz([1] * 24, O)


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
        O = matmul(A, mat_inv(B))
        if all(x != 0 for row in O for x in row):
            return O


def is_SO4(O):
    Tt = [[O[j][i] for j in range(4)] for i in range(4)]
    return matmul(O, Tt) == [[Fraction(int(i == j)) for j in range(4)] for i in range(4)] and det4(O) == 1


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


def quat_matrix(p, q, one, zero):
    """Matrix of x -> p x q in the basis 1, i, j, k: column k is the image of the k-th basis
    quaternion.  `one`, `zero` are the ring's 1 and 0.  (x -> p x q with unit p, q runs over
    all of SO(4): classical, CITED.)"""
    cols = []
    for k in range(4):
        e = tuple(one if m == k else zero for m in range(4))
        cols.append(quat_mul(quat_mul(p, e), q))
    return [[cols[j][i] for j in range(4)] for i in range(4)]


# ----------------------------------------------------------------------------------------
# 8. integer polynomials in several variables on top of python-flint 0.6 (fmpz_mpoly)
# ----------------------------------------------------------------------------------------
# python-flint 0.6 offers: generators, + - * ** with a polynomial on the left, exact division
# (divmod), gcd, len, coefficient(i), exponent_tuple(i).  Everything else is built here.


class Ring:
    """Z[names] as flint polynomials (lexicographic order)."""

    def __init__(self, names):
        from flint import fmpz_mpoly
        self.names = list(names)
        self.n = len(self.names)
        self.gens = list(fmpz_mpoly.gens(self.n))
        self.g = dict(zip(self.names, self.gens))
        self.zero = fmpz_mpoly(0, self.n)
        self.one = fmpz_mpoly(1, self.n)
        self._pw = [[self.one] for _ in range(self.n)]

    def const(self, c):
        from flint import fmpz_mpoly
        return fmpz_mpoly(int(c), self.n)

    def pw(self, i, e):
        L = self._pw[i]
        while len(L) <= e:
            L.append(L[-1] * self.gens[i])
        return L[e]

    def mono(self, ex):
        m = self.one
        for i, e in enumerate(ex):
            if e:
                m = m * self.pw(i, e)
        return m

    def from_terms(self, terms):
        """terms: iterable of (exponent tuple, int coefficient); equal exponents are added."""
        return tree_sum([self.mono(ex) * int(co) for ex, co in terms if co], self.zero)

    @staticmethod
    def terms(f):
        """list of (exponent tuple, int coefficient)."""
        return [(tuple(int(x) for x in f.exponent_tuple(i)), int(f.coefficient(i)))
                for i in range(len(f))]

    def degrees(self, f):
        d = [0] * self.n
        for i in range(len(f)):
            ex = f.exponent_tuple(i)
            for k in range(self.n):
                if ex[k] > d[k]:
                    d[k] = int(ex[k])
        return d

    def derivative(self, f, i):
        out = []
        for ex, co in self.terms(f):
            if ex[i]:
                e2 = list(ex)
                e2[i] -= 1
                out.append((tuple(e2), co * ex[i]))
        return self.from_terms(out)

    def evaluate(self, f, point):
        """exact value at a rational point (list of Fractions / ints)."""
        tot = Fraction(0)
        pows = [[Fraction(1)] for _ in range(self.n)]
        for ex, co in self.terms(f):
            t = Fraction(co)
            for i, e in enumerate(ex):
                if e:
                    L = pows[i]
                    while len(L) <= e:
                        L.append(L[-1] * point[i])
                    t *= L[e]
            tot += t
        return tot


def tree_sum(items, zero):
    items = list(items)
    if not items:
        return zero
    while len(items) > 1:
        nxt = [items[i] + items[i + 1] for i in range(0, len(items) - 1, 2)]
        if len(items) % 2:
            nxt.append(items[-1])
        items = nxt
    return items[0]


def coeff_stats(f):
    """(number of terms, smallest coefficient, number of negative coefficients)."""
    mn, neg = None, 0
    for i in range(len(f)):
        c = int(f.coefficient(i))
        if mn is None or c < mn:
            mn = c
        if c < 0:
            neg += 1
    return len(f), mn, neg


def poly_sha(f):
    """sha256 of the sorted term list 'e1,e2,...:coefficient' (one term per line)."""
    lines = sorted((tuple(int(x) for x in f.exponent_tuple(i)), int(f.coefficient(i)))
                   for i in range(len(f)))
    txt = "\n".join(",".join(map(str, ex)) + ":" + str(co) for ex, co in lines)
    return hashlib.sha256(txt.encode()).hexdigest()


def subst_clear(terms, nums, dens, ring_out, degs=None):
    """sum over terms of  co * prod_i nums[i]^e_i * dens[i]^(degs[i] - e_i):  the polynomial
    obtained from f = sum co x^e by x_i = nums[i]/dens[i] and multiplication by
    prod_i dens[i]^degs[i].  degs defaults to the degrees of f in each variable.  nums, dens
    are polynomials of ring_out.  Horner-type recursion over the variables."""
    terms = [(tuple(ex), co) for ex, co in terms if co]
    n = len(nums)
    if degs is None:
        degs = [max((ex[i] for ex, _ in terms), default=0) for i in range(n)]
    npow = [[ring_out.one] for _ in range(n)]
    dpow = [[ring_out.one] for _ in range(n)]
    for i in range(n):
        for _ in range(degs[i]):
            npow[i].append(npow[i][-1] * nums[i])
            dpow[i].append(dpow[i][-1] * dens[i])

    def rec(i, ts):
        if i == n:
            tot = 0
            for _, co in ts:
                tot += co
            return ring_out.const(tot)
        groups = {}
        for ex, co in ts:
            groups.setdefault(ex[i], []).append((ex, co))
        tot = ring_out.zero
        for e_, sub in groups.items():
            tot = tot + rec(i + 1, sub) * npow[i][e_] * dpow[i][degs[i] - e_]
        return tot

    if not terms:
        return ring_out.zero
    return rec(0, terms)


def laplace_det(M, zero, row_order=None):
    """Determinant of a square matrix of flint polynomials by Laplace expansion along the rows
    in `row_order`, with memoisation over column subsets.  Returns (det, number of minors)."""
    n = len(M)
    rows = list(range(n)) if row_order is None else list(row_order)
    perm_sign = sgn(rows)
    memo = {}

    def minor(k, cols):
        if k == n:
            return None                      # empty determinant = 1
        if cols in memo:
            return memo[cols]
        r = rows[k]
        tot = zero
        for pos, c in enumerate(cols):
            ent = M[r][c]
            if ent == zero:
                continue
            sub = minor(k + 1, cols[:pos] + cols[pos + 1:])
            if sub is None:
                term = ent
            elif sub == zero:
                continue
            else:
                term = ent * sub
            tot = tot - term if pos % 2 else tot + term
        memo[cols] = tot
        return tot

    d = minor(0, tuple(range(n)))
    return d * perm_sign, len(memo)


def bareiss_det(M, zero, one):
    """Fraction-free (Bareiss) determinant; a second, independent algorithm (slow: minutes)."""
    n = len(M)
    A = [list(r) for r in M]
    sign = 1
    prev = one
    for k in range(n - 1):
        piv = None
        for i in range(k, n):
            if A[i][k] != zero and (piv is None or len(A[i][k]) < len(A[piv][k])):
                piv = i
        if piv is None:
            return zero
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                q, r = divmod(A[i][j] * A[k][k] - A[i][k] * A[k][j], prev)
                assert r == zero
                A[i][j] = q
        prev = A[k][k]
    return A[n - 1][n - 1] * sign


# ----------------------------------------------------------------------------------------
# 9. the radial equations of PROOF.md 7.3 on a canonical representation
# ----------------------------------------------------------------------------------------


class Radial:
    """For a list `order` of eight facet indices with coefficient names `names` (same order):
    the ring Z[names, h], the coordinates u of  g + h kappa,  g = sum names_i f_{order_i},
    the cubic F = C(u) and the determinant of the full 14 x 14 Hessian of C at u."""

    def __init__(self, order, names):
        self.order = list(order)
        self.names = list(names)
        self.R = Ring(self.names + ["h"])
        R = self.R
        self.v = dict(zip(R.names, R.gens))
        for jx in self.order:
            assert all(x.denominator == 1 for x in F_COORD[jx])
        u = []
        for k in range(14):
            t = R.zero
            for w, jx in zip(R.gens[:8], self.order):
                ck = int(F_COORD[jx][k])
                if ck:
                    t = t + w * ck
            u.append(t)
        u[13] = u[13] + R.gens[8]                      # + h kappa, kappa = e_13
        self.u = u
        self.F = cubic_C(u)
        table = hess_linear_table()
        H = [[R.zero for _ in range(14)] for _ in range(14)]
        for i in range(14):
            for k in range(14):
                t = R.zero
                for cf, idx in table[i][k]:
                    t = t + u[idx] * cf
                H[i][k] = t
        self.H = H
        self._det = None

    def det(self):
        if self._det is None:
            R = self.R
            dens = [sum(1 for k in range(14) if self.H[r][k] != R.zero) for r in range(14)]
            ro = sorted(range(14), key=lambda r: -dens[r])
            self._det, self.n_minors = laplace_det(self.H, R.zero, ro)
        return self._det

    def cleared(self, Cst, Q):
        """Q^4 * det Hess C(g + (Cst/Q) kappa) as a polynomial (h eliminated)."""
        R = self.R
        byh = {}
        for ex, co in R.terms(self.det()):
            byh.setdefault(ex[8], []).append((ex[:8] + (0,), co))
        self.h_degree = max(byh)
        assert self.h_degree <= 4
        N = R.zero
        for k, tl in byh.items():
            N = N + R.from_terms(tl) * Cst ** k * Q ** (4 - k)
        return N

    def radial_factor(self, Cst, Q, var):
        """Returns (ok, G, S): S = Q^4 det / (-4 Q^2) (exact division), G = gcd(S, dS/dvar)
        normalised to G^2 = S.  ok is True iff both the division and the square are exact."""
        R = self.R
        N = self.cleared(Cst, Q)
        S, rem = divmod(N, Q * Q * (-4))
        if rem != R.zero:
            return False, None, None
        G = S.gcd(R.derivative(S, R.names.index(var)))
        ok = (G * G == S)
        return ok, G, S


# ----------------------------------------------------------------------------------------
# 10. the "coefficient Gram" sum-of-squares representation on SO(4)
# ----------------------------------------------------------------------------------------
QNAMES = ["p0", "p1", "p2", "p3", "q0", "q1", "q2", "q3"]


def parse_poly(text, ring, scale=1):
    """Parse a polynomial written with + - * ** and division by integers (e.g.
    'p0*p1*q2/2 - p1**2*q3/2') into scale * polynomial; no eval.  Raises unless the scaled
    coefficients are integers and every name is a generator of `ring`."""
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
                if op != "*" or tok not in ring.g:
                    raise ValueError("cannot parse %r" % text)
                mono = mono * ring.g[tok] ** exp
            k += 1
        c = num * scale
        if c.denominator != 1:
            raise ValueError("non-integer coefficient in %r" % text)
        tot = tot + mono * int(c)
    return tot


class GramFace:
    """Sum-of-squares representation of the Leibniz sums of one face.

    Data read from the supplement (as text): two blocks of eight cubic polynomials F_blk in
    (p, q) ('features'), an invertible rational 8 x 8 matrix T, four signed coordinate
    permutations S_0..S_3 of (p, q), and the list `order` of eight facet indices.
    With phi_blk = T^{-1} F_blk (the document's basis) define for a symmetric 8 x 8 matrix Gm

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
        self.Y = Y = Ring(QNAMES)
        yv = Y.gens
        self.pp = tree_sum([x * x for x in yv[:4]], Y.zero)
        self.qq = tree_sum([x * x for x in yv[4:]], Y.zero)
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
                self.allphi.append([subst_clear(Y.terms(f), img, [Y.one] * 8, Y) for f in ph])
        nn = self.pp + self.qq
        self.E = {}
        for (i, k) in self.keys:
            tot = Y.zero
            for ph in self.allphi:
                tot = tot + ph[i] * ph[k]
            self.E[(i, k)] = self.red(tot * nn)
        M = quat_matrix(yv[:4], yv[4:], Y.one, Y.zero)
        self.M = M
        MtM = [[tree_sum([M[k][i] * M[k][j] for k in range(4)], Y.zero) for j in range(4)] for i in range(4)]
        n2 = self.pp * self.qq
        self.model_ok = all(MtM[i][j] == (n2 if i == j else Y.zero) for i in range(4) for j in range(4)) \
            and leibniz([1] * 24, M, Y.zero) == n2 * n2      # M^t M = |p|^2|q|^2 I and det M = (|p|^2|q|^2)^2
        self.cost = [self.red(leibniz(F_VALUES[m], M, Y.zero) * (4 * self.scale)) for m in self.order]
        self._solve()

    def red(self, f):
        """normal form modulo N = |p|^2 - |q|^2:  p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2."""
        Y = self.Y
        yv = Y.gens
        sub = self.qq - yv[1] * yv[1] - yv[2] * yv[2] - yv[3] * yv[3]
        groups = {}
        for ex, co in Y.terms(f):
            groups.setdefault(ex[0] // 2, []).append(((ex[0] % 2,) + ex[1:], co))
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
        Y = self.Y
        cols = []
        for key in self.keys:
            mult = 1 if key[0] == key[1] else 2
            cols.append({ex: co * mult for ex, co in Y.terms(self.E[key])})
        for cst in self.cost:
            cols.append({ex: -co for ex, co in Y.terms(cst)})
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

    def weights_of(self, Gm):
        """the weight vector w forced by (*) for a symmetric matrix Gm (any exact ring)."""
        ent = [Gm[i][k] for (i, k) in self.keys]
        return [sum(c * e for c, e in zip(self.weight_rows[m], ent) if c != 0) for m in range(8)]

    def relation_values(self, Gm):
        ent = [Gm[i][k] for (i, k) in self.keys]
        return [sum(c * e for c, e in zip(r, ent) if c != 0) for r in self.relations]

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
        rhs = self.red(leibniz(vals, self.M, Y.zero) * (4 * self.scale))
        return lhs == rhs


# ----------------------------------------------------------------------------------------
# 11. Bernstein tables
# ----------------------------------------------------------------------------------------


def bernstein_polys(RC, k, N, D):
    """For G = k0 + k1 x + k2 x^2 on [0, N/D] (D > 0): the numerators
    middle = 2 D k0 + N k1, endpoint = D^2 k0 + N D k1 + N^2 k2 of the Bernstein coefficients
    b1 = middle / (2D), b2 = endpoint / D^2  (b0 = k0)."""
    return k[0] * D * 2 + k[1] * N, k[0] * D * D + k[1] * N * D + k[2] * N * N


def bernstein_identity_check():
    """sympy: D^2 G(t N/D) = (1-t)^2 D^2 k0 + t(1-t) D (2 D k0 + N k1) + t^2 (D^2 k0 + N D k1 + N^2 k2)."""
    import sympy
    k0, k1, k2, N, D, t = sympy.symbols("k0 k1 k2 N D t")
    x = t * N / D
    lhs = D ** 2 * (k0 + k1 * x + k2 * x ** 2)
    rhs = (1 - t) ** 2 * D ** 2 * k0 + t * (1 - t) * D * (2 * D * k0 + N * k1) \
        + t ** 2 * (D ** 2 * k0 + N * D * k1 + N ** 2 * k2)
    return sympy.expand(lhs - rhs) == 0


# ----------------------------------------------------------------------------------------
# 12. common checks
# ----------------------------------------------------------------------------------------
SUPP_SUB = os.path.join("proofs", "manin_transfer_20261021")

# sha256 of the supplement files that were read for these replays (data only)
DOC_SHA256 = {
    "root_moc4_T46_whole_assembly_20261023.pdf": "90f2b739c757eda6ca4798f09e946151f27ac0513ebfb70279ec3fac44f41799",
    "final_referee_T46_six_ray_coherent_capacity_20261023.pdf": "f2adf23794db91ff26db546c2b7ec3ef4900de8de6d4246b699192ec7a302a9b",
    "final_referee_T46_core_coherent_capacity_20261023.pdf": "f6a1196fd31f8d2e76781372004bdcddf4846edd0513e2810ca924cdae1eb5d1",
    "final_referee_T46_two_Heron_necessity_20261023.pdf": "0a32065e9eae0f7e2e170136328d12d982393a2350fd8b0c14bd9d19cbfec4ec",
    "final_referee_T46_negative_j_Heron_uniform_20261023.pdf": "ae3a4791889ae704353ce3ded6e1b4a96f2d2788ed4f89d1293381073cc6c0f5",
    "root_moc4_T46_full_radial_20261022.json": "4140ecfc03dfe68448043c4645f3d6a0b8f71a2ee56422147dfb31566e664370",
    "final_referee_T46_budget_PSD8_20261022.json": "7d9f740d5707447ec48d7a021d6d8e7587927de44af79d028e83bfbb2141e7c8",
    "final_referee_T46_cubic_reduction_20261022.json": "649a96d9e9826f2ee74ca22cf2966da839319fff48c9956039ca32893df1d02a",
    "final_referee_T46_j_squareclass_20261022.json": "0d34ebe220ed8575bffacae22e9e89d10d8418eba140e56094ec6ec3363d1b78",
    "final_referee_T46_closed_budget_faces_20261023.json": "6abfc3a4f3e7b7464fc734c5999d097f18e9f9284a0577c8dd144af150aa31b7",
    "independent_moc4_T46_preclosure_allowlist_20261023.json": "aff435ac310b2c24930b68b752d18ca25728f00462f72b6e3f76e7144dd6257c",
    "final_referee_T46_core_T47_transports_20261023.json": "3b1666e48870890ce851a55086de8b1c056b874701af0c8bb93a40f060a4cd82",
    "independent_moc4_T47_budget_PSD8_20261022.json": "1ad79a108e7b1406ba64a8377e2e6f633d2bfa067ee11d6abbde1e95fefba15e",
}

# Appendix E of the paper (p. 60), read from the text of the PDF; compared with the audit library
PAPER_T = {1: (0, 1, 2, 3, 6, 7), 7: (0, 1, 2, 4, 7, 8), 25: (3, 4, 7, 8, 10, 12), 29: (0, 4, 5, 6, 7, 8),
           36: (2, 3, 4, 6, 11, 12), 39: (2, 3, 5, 7, 9, 12), 42: (1, 2, 3, 6, 7, 8),
           47: (1, 3, 5, 6, 7, 8), 48: (0, 2, 5, 9, 12, 13), 65: (0, 4, 5, 6, 9, 12), 67: (2, 4, 5, 7, 9, 12),
           77: (0, 1, 3, 10, 11, 16)}


def supp_dir():
    root = os.environ.get("MOC4_SUPP")
    if not root:
        return None
    d = os.path.join(root, SUPP_SUB)
    return d if os.path.isdir(d) else None


def check_file_identity(ck, names):
    """Provenance only: the files under MOC4_SUPP are the ones read for this replay."""
    d = supp_dir()
    if d is None:
        ck.info("MOC4_SUPP not set (or folder missing): file hashes are not compared")
        return False
    ok_all = True
    for nm in names:
        path = os.path.join(d, nm)
        if not os.path.isfile(path):
            ck.check("EXACT", "file present under MOC4_SUPP: " + nm, False)
            ok_all = False
            continue
        with open(path, "rb") as fh:
            h = hashlib.sha256(fh.read()).hexdigest()
        ok_all &= ck.check("EXACT", "sha256 of %s is the one of the file that was read" % nm,
                           h == DOC_SHA256[nm], h[:16] + "...")
    return ok_all


def print_versions(ck):
    import sympy
    import flint
    line = "python %s, sympy %s, python-flint %s" % (sys.version.split()[0], sympy.__version__, flint.__version__)
    try:
        import numpy
        import scipy
        line += ", numpy %s, scipy %s (NUMERIC parts only)" % (numpy.__version__, scipy.__version__)
    except Exception:
        pass
    ck.info(line)


def check_foundations(ck):
    """The objects of Sections 3 and 4 of the paper that the cases rely on (as in
    replay_iv_common.check_foundations)."""
    ck.section("Foundations rebuilt from the paper's definitions")
    ck.check("EXACT", "P is a 24 x 14 integer matrix; every row is a polynomial identity in 16 variables "
             "(checked at construction); rank 14; P kappa = (1,...,1) for kappa = e_13",
             len(P) == 24 and all(len(r) == 14 for r in P) and rank(P) == 14
             and all(P[s][13] == 1 for s in range(24)))
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
             "(so X lies in the singular locus of the cubic C of (A.2))", not bad)
    prof_ok = True
    for j in range(82):
        v = F_VALUES[j]
        cnt = (v.count(0), v.count(1), v.count(2))
        prof_ok &= (cnt in ((18, 6, 0), (20, 4, 0))) if F_IS_MINOR[j] else (cnt == (13, 10, 1))
    ck.check("EXACT", "82 vectors f_j of Appendix D: all in the image of P; 34 minor vectors with 18 or 20 zeros, "
             "48 vectors h with value profile 13/10/1; pairwise distinct",
             prof_ok and len(set(F_VALUES)) == 82 and sum(F_IS_MINOR) == 34)
    ok = len(set(G_MAPS)) == 1152
    for t in G_MAPS:
        for j in range(82):
            jj = F_INDEX.get(g_values(t, F_VALUES[j]))
            if jj is None or F_IS_MINOR[jj] != F_IS_MINOR[j]:
                ok = False
    ck.check("EXACT", "group G: 1152 distinct maps of the 24 permutation indices; G permutes the 82 vectors f_j "
             "and maps minor vectors to minor vectors", ok)


def check_group_covariance(ck, seed=2501, trials=3):
    """Exact spot checks of PROOF.md 6.6 for ALL 1152 elements: (i) the g-image of the vector of
    24 products of spectra (a, b) is the vector of products of relabelled spectra, so G
    preserves X; (ii) Phi_{g.c}(O'') = Phi_c(O) with O'' in SO(4)."""
    import random
    rng = random.Random(seed)
    ok1 = ok2 = True
    for _ in range(trials):
        a = [GQ(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
             for _ in range(4)]
        b = [GQ(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
             for _ in range(4)]
        z = z_points(a, b)
        O = random_rational_SO4(rng)
        vals = [Fraction(rng.randint(-5, 5)) for _ in range(24)]
        ph = leibniz(vals, O)
        for g, gm in zip(G_ELEMENTS, G_MAPS):
            a2, b2 = g_spectra(g, a, b)
            ok1 &= tuple(z_points(a2, b2)) == g_values(gm, z)
            O2 = transfer_to_SO4(g, O)
            ok2 &= leibniz(g_values(gm, vals), O2) == ph
        ok2 &= is_SO4(O)
    ck.check("EXACT", "for all 1152 g and %d random Gaussian-rational spectra: the products of the relabelled "
             "spectra g_spectra(g, a, b) are the g-image of the products of (a, b)  (G preserves X and kappa)"
             % trials, ok1)
    ck.check("EXACT", "for all 1152 g, %d random rational O in SO(4) and random integer value vectors c: "
             "Phi_{g.c}(O'') = Phi_c(O) with O'' = transfer_to_SO4(g, O) (a row/column permutation of O or O^t, "
             "one column negated if needed)" % trials, ok2)


def check_radial_necessity(ck, seed=20261125, samples=10):
    """PROOF.md 7.3: x = u + i v in X  =>  C(u) = 0 and Hess C(u) v = 0 (so det Hess C(u) = 0).
    (As in replay_iv_common.check_radial_necessity.)"""
    import random
    import sympy
    ck.section("Radial necessity (PROOF.md 7.3): the real part of a point of X")
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
    pol = sympy.expand(sum(us[j] * grad[j].xreplace(sv) for j in range(14))
                       - sympy.Rational(1, 2) * sum(vs[j] * hess[j][k].xreplace(su) * vs[k]
                                                    for j in range(14) for k in range(14)))
    ck.check("EXACT", "Euler u.grad C(u) = 3 C(u) and polarisation u.grad C(v) = (1/2) v^t Hess C(u) v, identically",
             euler == 0 and pol == 0)
    ck.info("hence for x = u + iv with grad C(x) = 0: Hess C(u) v = 0 and 3 C(u) = u.grad C(v) = (1/2) v^t Hess C(u) v = 0; "
            "if v = 0 then Hess C(u) u = 2 grad C(u) = 0; so det Hess C(u) = 0 in every case (u = 0 included)")
    rng = random.Random(seed)
    okc = okh = True
    for _ in range(samples):
        def rq():
            return GQ(Fraction(rng.randint(-9, 9), rng.randint(1, 5)), Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
        a = [rq() for _ in range(4)]
        b = [rq() for _ in range(4)]
        al = rq()
        x = [GQ.of(c) for c in matching_vector(a, b, al)]
        u = [c.re for c in x]
        v = [c.im for c in x]
        okc &= cubic_C(u) == 0
        H = hess_at(u)
        okh &= all(sum(H[i][k] * v[k] for k in range(14)) == 0 for i in range(14))
        okc &= [GQ.of(sum(P[s][r] * x[r] for r in range(14))) for s in range(24)] == [al * z for z in z_points(a, b)]
    ck.check("EXACT", "spot check on %d random Gaussian-rational points x = alpha Psi(a,b) (seed %d): "
             "(P x)_sigma = alpha z_sigma, C(Re x) = 0 and Hess C(Re x) Im x = 0" % (samples, seed), okc and okh)


def check_minor_leibniz(ck, minor_indices):
    """Lemma 6.2 for the listed minor vectors: Laplace and Jacobi identities in a generic 4 x 4
    matrix give Phi_{n_{I,J}}(O) = det O[I,J]^2 on SO(4).  (As in replay_iv_common.)"""
    import sympy
    M = sympy.Matrix(4, 4, lambda i, j: sympy.Symbol("o%d%d" % (i + 1, j + 1)))
    adj = M.adjugate()
    detM = sympy.expand(M.det(method="berkowitz"))
    ok = True
    for j in minor_indices:
        I, J = minor_spec(j)
        I, J = list(I), list(J)
        Ic = [i for i in range(4) if i not in I]
        Jc = [i for i in range(4) if i not in J]
        eps = (-1) ** (sum(I) + sum(J))
        phi = 0
        for s, val in zip(PERMS, F_VALUES[j]):
            if val:
                phi += sgn(s) * M[0, s[0]] * M[1, s[1]] * M[2, s[2]] * M[3, s[3]]
        dIJ = M.extract(I, J).det()
        dc = M.extract(Ic, Jc).det()
        ok &= sympy.expand(phi - eps * dIJ * dc) == 0
        ok &= sympy.expand(adj.extract(J, I).det() - eps * dc * detM ** (len(I) - 1)) == 0
    ck.check("EXACT", "Laplace and Jacobi identities for the %d minor generators %s: "
             "Phi_{f_j}(O) = (minor of O)^2 >= 0 on SO(4)" % (len(minor_indices), list(minor_indices)), ok)


def check_contact_set(ck, nu, doc_T, doc_name, other):
    ck.section("Contact set and numbering")
    ck.check("EXACT", "the document's set %s (its label '%s') equals the paper's T_%d of Appendix E as a set of "
             "lexicographic permutation indices (group element: identity)" % (sorted(doc_T), doc_name, nu),
             tuple(sorted(doc_T)) == PAPER_T[nu])
    r1, r2 = orbit_rep(PAPER_T[nu]), orbit_rep(PAPER_T[other])
    ck.check("EXACT", "orbit representative under G: %s; T_%d (representative %s) lies in a different orbit"
             % (r1, other, r2), orbit_rep(doc_T) == r1 and r1 != r2)


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
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d, no use "
             "of the 82-facet list) gives exactly the %d vectors f_j, j in F_T" % (d, len(gens)),
             rays == want, "%d rays, %.1f s" % (len(rays), time.time() - t0))
    ck.info("generators: " + ", ".join("f%d = %s" % (j, f_name(j)) for j in gens))
    return gens


def crosscheck_with_supplement_data(ck):
    """Read two JSON data files of the supplement (as data) and compare the stored P, cubic and
    facet table with the objects rebuilt here: shows only that the documents talk about the
    same P, the same cubic (same sign) and the same numbering f_0..f_81."""
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: the supplement's stored P, cubic and facet table are not compared")
        return
    f1 = os.path.join(root, SUPP_SUB, "root_scout_moc_spinor_20261021.json")
    f2 = os.path.join(root, "proofs", "principal_minor_result_20261022",
                      "root_permutation_polytope_complete_20261022.json")
    if not (os.path.isfile(f1) and os.path.isfile(f2)):
        ck.info("supplement data files not found: comparison skipped")
        return
    import sympy
    with open(f1) as fh:
        meta = json.load(fh)
    with open(f2) as fh:
        fac = json.load(fh)
    ms, C, _, _ = cubic_data()
    C2 = sum(sympy.Rational(str(co)) * sympy.prod([ms[int(j)] for j in ix]) for ix, co in meta["cubic_terms"])
    ck.check("EXACT", "supplement data (read as JSON): its stored 24 x 14 matrix equals P, its stored cubic equals "
             "the paper's C of (A.2) with the same sign, and its 82 facet value vectors equal mine in the same order",
             meta["permutation_to_basis"] == P and sympy.expand(C2 - C) == 0
             and [tuple(z["values"]) for z in fac["facets"]] == F_VALUES)


def crosscheck_with_audit_library(ck, contact_sets):
    here = os.path.dirname(os.path.abspath(__file__))
    lib = os.path.normpath(os.path.join(here, "..", "..", "zheng-n4-audit", "lib"))
    if not os.path.isfile(os.path.join(lib, "moc4lib.py")):
        ck.info("audit library not found at %s: comparison skipped" % lib)
        return
    sys.path.insert(0, lib)
    try:
        import moc4lib as L
    finally:
        sys.path.pop(0)
    import sympy
    ms, C, _, _ = cubic_data()
    ck.check("EXACT", "audit library moc4lib (imported only for this comparison): PERMS, P, kappa, the 82 value "
             "vectors, the minor flags, the coordinates of all f_j, the cubic C and the 1152 group maps agree with mine",
             L.PERMS == PERMS and L.P == P and list(L.KAPPA) == KAPPA
             and [tuple(v) for v in L.F_VALUES] == F_VALUES and list(L.F_IS_MINOR) == F_IS_MINOR
             and [list(v) for v in L.F_VECS] == F_COORD and sympy.expand(L.cubic_C(ms) - C) == 0
             and set(L._G_MAPS) == set(G_MAPS))
    ck.check("EXACT", "audit library: CONTACT_SETS and generators agree for %s" % sorted(contact_sets),
             all(tuple(sorted(L.CONTACT_SETS[k])) == PAPER_T[k] and L.generators(PAPER_T[k]) == generators(PAPER_T[k])
                 for k in contact_sets))


def orbit_check(ck, family, rep, target, sym_positions=None):
    """`target` (facet indices with roles) lies in the G-orbit of the proved representative
    `rep`; the roles in sym_positions may be exchanged.  Returns a group element or None."""
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
    ck.check("EXACT", "%s: the roles %s lie in the G-orbit of the representative %s"
             % (family, tuple(target), tuple(rep)), hit is not None,
             "" if hit is None else "g = (rho=%s, rho'=%s, transpose=%d); %d group elements; source roles %s"
             % (hit[0][0], hit[0][1], hit[0][2], hit[2], hit[1]))
    return None if hit is None else hit[0]


# ----------------------------------------------------------------------------------------
# 13. NUMERIC helpers (floating point; evidence only)
# ----------------------------------------------------------------------------------------


def numeric_leibniz_weights(Pq, Qq):
    """numpy: the 24 numbers sgn(sigma) prod_i M[i, sigma(i)] for M = matrix of x -> p x q,
    p = Pq/|Pq|, q = Qq/|Qq|."""
    import numpy as np
    p = Pq / np.linalg.norm(Pq)
    q = Qq / np.linalg.norm(Qq)
    cols = []
    for k in range(4):
        e = [0.0, 0.0, 0.0, 0.0]
        e[k] = 1.0
        cols.append(quat_mul(quat_mul(tuple(p), tuple(e)), tuple(q)))
    M = np.array(cols).T
    return np.array([SGN[s] * M[0, PERMS[s][0]] * M[1, PERMS[s][1]] * M[2, PERMS[s][2]] * M[3, PERMS[s][3]]
                     for s in range(24)])


def numeric_min_phi(vals, rng, restarts=12):
    """NUMERIC: smallest value found of Phi_vals over SO(4) (multistart BFGS in quaternion
    coordinates).  An upper bound for the true minimum, not a certified minimum."""
    import numpy as np
    from scipy.optimize import minimize
    vals = np.asarray(vals, float)

    def f(x):
        return float(vals @ numeric_leibniz_weights(x[:4], x[4:]))
    best = None
    for _ in range(restarts):
        r = minimize(f, rng.normal(size=8), method="BFGS", options={"gtol": 1e-12, "maxiter": 400})
        best = r.fun if best is None else min(best, r.fun)
    return best


def numeric_stratum_search(T, seed, starts, eps=0.02, tol=1e-10, extra=None):
    """NUMERIC (as in replay_iv_common.numeric_stratum_search).  Search for spectra a, b in C^4
    and a phase alpha such that all gaps Re(alpha z_sigma) - Re(alpha z_pi*) are >= 0, vanish on
    T, and sum to 1 (pi* = first index of T).  Two Levenberg-Marquardt stages from each random
    start; stage 1 pushes the eight labels a_i, -b_j apart, stage 2 polishes without that term.
    `extra`, if given, maps the gap vector to a list of additional residuals (used to steer the
    search into a part of the stratum).  Returns a list of dicts with H = Re(alpha z_pi*), gaps,
    sep (smallest normalised label distance) and the spectra."""
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
        if extra is not None:
            r += list(extra(g))
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
            s2 = least_squares(res2, s1.x, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
        except Exception:
            continue
        r, F, g, a, b = base(s2.x)
        if not np.all(np.isfinite(r)) or float(np.max(np.abs(r))) >= tol:
            continue
        al = np.exp(1j * s2.x[16])
        out.append(dict(H=float(F[T[0]]), sep=float(np.min(seps(a, b))), gaps=g, a=a, b=b, alpha=al))
    return out


# ----------------------------------------------------------------------------------------
# 14. symbolic checks of Gram constructions (sympy), exact rational determinants
# ----------------------------------------------------------------------------------------


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
    Checks: the ten re-derived relations, the weights, the pivots and zero residual, and the
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
    wnum = [sp.cancel(sp.sympify(w).subs(sub)) for w in weights]
    den = 1
    for x in list(Gnum) + wnum:
        x = sp.cancel(x)
        assert x.is_Rational, x
        den = sp.ilcm(den, x.q)
    Gint = [[int(Gnum[i, k] * den) for k in range(8)] for i in range(8)]
    wint = [int(w * den) for w in wnum]
    direct = gf.direct_identity(Gint, wint)
    ck.check("EXACT", "%s: the Gram matrix satisfies the ten relations; its weights are %s; elimination in the "
             "order %s has pivots %s and zero residual; identity (*) holds directly at the rational instance %s%s"
             % (label, [str(w) for w in weights], list(order), [str(sp.factor(e)) for e in expected_pivots],
                {str(s): str(v) for s, v in instance.items()}, note),
             rel_ok and w_ok and piv_ok and direct,
             "relations %s, weights %s, pivots %s, direct %s" % (rel_ok, w_ok, piv_ok, direct))
    return rel_ok and w_ok and piv_ok and direct


def det_fraction(M):
    """exact determinant of a square matrix of Fractions (Gaussian elimination)."""
    A = [[Fraction(x) for x in r] for r in M]
    n = len(A)
    det = Fraction(1)
    for c in range(n):
        k = next((i for i in range(c, n) if A[i][c] != 0), None)
        if k is None:
            return Fraction(0)
        if k != c:
            A[c], A[k] = A[k], A[c]
            det = -det
        det *= A[c][c]
        iv = 1 / A[c][c]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] * iv
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return det


# ----------------------------------------------------------------------------------------
# 15. additions for T_47 part A (new)
# ----------------------------------------------------------------------------------------


def fsqrt(R, f):
    """Square root of a flint polynomial of the Ring R, or None.  Recursive in the variables:
    for f = sum c_k x^k of degree 2n in its first occurring variable x, w_n = sqrt(c_2n) and
    w_i = (c_{n+i} - sum_{a+b=n+i, i<a,b<n} w_a w_b) / (2 w_n) by exact division.  The result is
    returned only if its square equals f (final exact comparison), so it is a certified root."""
    import math
    if f == R.zero:
        return R.zero
    degs = R.degrees(f)
    v = next((i for i, dg in enumerate(degs) if dg > 0), None)
    if v is None:
        c = int(f.coefficient(0))
        if c < 0:
            return None
        r = math.isqrt(c)
        return R.const(r) if r * r == c else None
    if degs[v] % 2:
        return None
    n = degs[v] // 2
    groups = {}
    for ex, co in R.terms(f):
        groups.setdefault(ex[v], []).append((ex[:v] + (0,) + ex[v + 1:], co))
    cs = {k: R.from_terms(tl) for k, tl in groups.items()}
    wn = fsqrt(R, cs.get(2 * n, R.zero))
    if wn is None or wn == R.zero:
        return None
    w = {n: wn}
    two = wn * 2
    for i in range(n - 1, -1, -1):
        acc = cs.get(n + i, R.zero)
        for a_ in range(i + 1, n):
            b_ = n + i - a_
            if i < b_ < n:
                acc = acc - w[a_] * w[b_]
        quo, rem = divmod(acc, two)
        if rem != R.zero:
            return None
        w[i] = quo
    x = R.gens[v]
    out = R.zero
    xp = R.one
    for i in range(0, n + 1):
        out = out + w[i] * xp
        xp = xp * x
    return out if out * out == f else None


def coeff_in(R, f, var, k):
    """coefficient of (variable var)^k in f, as a polynomial of R not involving var."""
    v = R.names.index(var)
    return R.from_terms([(ex[:v] + (0,) + ex[v + 1:], co) for ex, co in R.terms(f) if ex[v] == k])


def total_degrees(R, f):
    return sorted(set(sum(ex) for ex, _ in R.terms(f)))


def doc_action(row, col, transpose):
    """The supplement's convention for "the action with row r, column c and transpose flag":
    the index map a(sigma) = r o sigma^(+-1) o c of the 24 permutations (this reproduces the six
    stored index maps of the closed faces, see the assembly script), used as a PULLBACK: a
    profile v becomes v o a.  In this module's notation the pullback by a is the action of the
    group element g = a^{-1}.  Returns (index map a, index map of g = a^{-1})."""
    a = g_index_map((tuple(row), tuple(col), 1 if transpose else 0))
    ainv = [0] * 24
    for s in range(24):
        ainv[a[s]] = s
    return a, tuple(ainv)


def doc_action_facets(row, col, transpose, facets):
    """images of the facet indices under the pullback by the document's action."""
    _, ginv = doc_action(row, col, transpose)
    return [g_facet(ginv, j) for j in facets]


def tangent_space_rows(pts):
    """Rows spanning the tangent space of X at Psi(v), v = pts (eight pairs of rationals): the
    16 partial derivatives of (B_0, ..., B_13) (the cone direction is among their span, by
    Euler)."""
    rows = []
    for slot in range(16):
        row = []
        for r in range(14):
            tot = Fraction(0)
            for key, c in B_POLYS[r].items():
                e = (key >> (3 * slot)) & 7
                if not e:
                    continue
                t = Fraction(c * e)
                k = key
                for s2 in range(16):
                    ee = k % 8
                    k //= 8
                    if s2 == slot:
                        ee -= 1
                    if ee:
                        base = pts[s2][0] if s2 < 8 else pts[s2 - 8][1]
                        t *= Fraction(base) ** ee
                tot += t
            row.append(tot)
        rows.append(row)
    return rows


def exact_inertia(M):
    """(n+, n-, n0) of a symmetric matrix of Fractions, by exact symmetric elimination with
    diagonal pivots and, where every remaining diagonal entry vanishes, a 2 x 2 hyperbolic
    step.  Congruence transformations only, so Sylvester's law gives the inertia."""
    A = [[Fraction(x) for x in r] for r in M]
    n = len(A)
    idx = list(range(n))
    npos = nneg = 0
    while idx:
        piv = next((i for i in idx if A[i][i] != 0), None)
        if piv is not None:
            d = A[piv][piv]
            if d > 0:
                npos += 1
            else:
                nneg += 1
            rest = [i for i in idx if i != piv]
            for i in rest:
                f = A[i][piv] / d
                if f:
                    for k in rest:
                        A[i][k] -= f * A[piv][k]
            idx = rest
            continue
        pair = next(((i, k) for i in idx for k in idx if i < k and A[i][k] != 0), None)
        if pair is None:
            break
        i, k = pair
        # replace e_i by e_i + e_k: the new diagonal entry is 2 A[i][k] != 0
        for m in range(n):
            A[i][m] += A[k][m]
        for m in range(n):
            A[m][i] += A[m][k]
    return npos, nneg, n - npos - nneg


def numeric_equality_points(T, seed, n, tol=1e-12, sep_min=0.05):
    """NUMERIC (copied from ../cases-v-a/replay_T29_common.py).  Points x = alpha Psi(a, b) of X
    whose real parts satisfy only the contact EQUATIONS Re(F_sigma - F_pi*) = 0 for sigma in T
    (the other gaps may have either sign): Gauss-Newton with least-norm steps from random
    starts; kept if the scaled residual is < tol and the eight labels stay apart.  Returns a
    list of dicts gaps (sum of squares 1), H = Re F_pi*, sep."""
    import numpy as np
    T = sorted(T)
    S = np.array(PERMS)
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(8, 1)

    def gaps_of(x):
        a, b, al = x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])
        F = (al * np.prod(a[None, :] + b[S], axis=1)).real
        return F - F[T[0]], F[T[0]], a, b

    def res(x):
        g = gaps_of(x)[0]
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
        g, H, a, b = gaps_of(x)
        z = np.concatenate([a, -b])
        Dm = np.abs(z[:, None] - z[None, :])[iu]
        sep = float(np.min(Dm) / np.sqrt(np.mean(Dm ** 2)))
        if sep < sep_min:
            continue
        sc = float(np.sqrt(np.sum(g * g)))
        out.append(dict(gaps=g / sc, H=float(H) / sc, sep=sep))
    return out


def float_eval(terms, x, absolute=False):
    """NUMERIC evaluation of a term list [(exponents, coefficient)] at the float vector x
    (absolute = True: sum of the absolute values of the monomials)."""
    tot = 0.0
    for ex, co in terms:
        t = abs(float(co)) if absolute else float(co)
        for i, e_ in enumerate(ex[:len(x)]):
            if e_:
                t *= (abs(x[i]) if absolute else x[i]) ** e_
        tot += t
    return tot
