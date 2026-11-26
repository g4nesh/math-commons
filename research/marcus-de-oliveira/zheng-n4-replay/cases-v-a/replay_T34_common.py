"""Shared exact objects for the independent replay of case T_34 = {2,3,5,9,10,12} (group V) of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4" (October 2026).

Provenance of the code in this file
  * Sections 1-7 and 9-10 (check registry, permutations, the group of order 1152, exact linear
    algebra, bracket polynomials, the matrix P, the cubic (A.2), the 82 vectors of Appendix D,
    contact cones, Leibniz sums, rational SO(4), quaternion matrices, Gaussian rationals and
    the generic checks) are COPIED, with small adaptations, from this owner's earlier replay
    modules `../cases-i-iv/replay_iv_common.py` and `../cases-i-iv/replay_T3_common.py`
    (same owner, same task family), with the two group conventions merged.  The Hessian
    table (`hess_linear`, `hess_at`) is rewritten so that it works over any ring.  Those modules rebuilt every object from the DEFINITIONS in the paper
    (Sections 3.1-3.2, 4.1-4.2, Appendices A, D, E; tables read from the page images).  The
    objects this replay relies on are re-checked here (`check_foundations`) and compared with
    the audit library `zheng-n4-audit/lib/moc4lib.py` (`crosscheck_with_audit_library`); the
    audit library is not used to build anything.
  * Section 8 (a thin layer over python-flint's sparse integer polynomials) is new.

Nothing here imports or executes anything from the author's supplement.

Exact arithmetic: Python ints, fractions.Fraction, python-flint fmpz_mpoly (integer
polynomials), sympy for a few small symbolic identities.  numpy is used only by functions
whose names start with `numeric_`.

Conventions
  * Permutations of {0,1,2,3} are tuples s with s[i] the image of i; PERMS is in
    lexicographic order, so PERMS[j] is the paper's pi_j.
  * Labels 1..8 carry vectors v_l = (x_l, y_l); [ij] = x_i y_j - x_j y_i.  With
    v_i = (a_i, 1), v_{4+j} = (-b_j, 1) one has [i, 4+j] = a_i + b_j.
  * V_R = Q^14 with coordinates m_0..m_13; values(c) = P c is the vector of the 24 numbers
    (P c)_sigma in PERMS order.  P is injective, so a vector of V_R is determined by its 24
    values; the f_j are handled through their values and through F_COORD = P^{-1}(values).
"""
from fractions import Fraction
import hashlib
import itertools
import os
import sys
import time

# ----------------------------------------------------------------------------------------
# 1. check registry
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
        sys.stdout.flush()

    def check(self, level, name, ok, detail=""):
        """level in {EXACT, NUMERIC}; ok must be a bool."""
        self.n += 1
        ok = bool(ok)
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
# 2. permutations and the group G of order 1152
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


# g = (r1, r2, e):  g.sigma = r1 o sigma o r2 (e = 0),  r1 o sigma^{-1} o r2 (e = 1).
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_perm(g, s):
    r1, r2, e = g
    s = inverse(s) if e else tuple(s)
    return compose(r1, compose(s, r2))


def g_index_map(g):
    return tuple(PIDX[g_perm(g, s)] for s in PERMS)


G_MAPS = [g_index_map(g) for g in G_ELEMENTS]


def orbit_rep(T):
    """Lexicographically smallest sorted image of the index set T under G."""
    T = list(T)
    return min(tuple(sorted(t[j] for j in T)) for t in G_MAPS)


def g_values(gmap, vals):
    """(g.c) defined by v_{g.sigma}(g.c) = v_sigma(c)."""
    out = [None] * 24
    for j in range(24):
        out[gmap[j]] = vals[j]
    return tuple(out)


def g_matrix(g, O):
    """g.O with Phi_{g.c}(g.O) = (+-1) Phi_c(O): for e = 0, (g.O)[i][j] = O[r2(i)][r1^{-1}(j)];
    for e = 1 the matrix is transposed first."""
    r1, r2, e = g
    if e:
        O = [[O[j][i] for j in range(4)] for i in range(4)]
    r1i = inverse(r1)
    return [[O[r2[i]][r1i[j]] for j in range(4)] for i in range(4)]


# ----------------------------------------------------------------------------------------
# 3. exact linear algebra over Q
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
    aug = [list(M[i]) + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    R, piv = rref(aug)
    if piv[:n] != list(range(n)):
        raise ValueError("singular matrix")
    return [row[n:] for row in R]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def ring_matmul(A, B):
    """Matrix product over a ring whose elements cannot be added to a Python int on the left
    (python-flint polynomials): no implicit 0 start value."""
    out = []
    for i in range(len(A)):
        row = []
        for j in range(len(B[0])):
            acc = A[i][0] * B[0][j]
            for k in range(1, len(B)):
                acc = acc + A[i][k] * B[k][j]
            row.append(acc)
        out.append(row)
    return out


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
# 4. bracket polynomials in sixteen variables, the basis (A.1) and the matrix P
#    A polynomial in x_1..x_8, y_1..y_8 is a dict {key: int}, key = sum of 8**slot over the
#    factors (slot l-1 for x_l, slot 7+l for y_l); adding keys multiplies monomials.
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


# ----------------------------------------------------------------------------------------
# 5. the cubic (A.2), read from the page image of p. 36 of the paper
# ----------------------------------------------------------------------------------------


def cubic_C(m):
    (m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13) = m
    return (m13 * (m0 * m6 - m1 * m5)
            - m0 * m9 * (m0 + m1 + m2 + m3 + m5 + m7 + m9) - m1 * m5 * m9
            + m10 * (-m0 * m3 - m0 * m9 + m1 * m2 + m1 * m7)
            + m11 * (-m0 * m7 - m0 * m9 + m2 * m5 + m3 * m5)
            + m12 * (-m0 * m9 + m1 * m7 - m2 * m6 - m2 * m9 + m3 * m5 + m3 * m7)
            + m2 * m10 * m11 - m3 * m8 * m10 - m4 * m7 * m11 + m4 * m8 * m9)


def cubic_terms():
    """C as a dict {(i, j, k) sorted: integer coefficient}: cubic_C is expanded once with sympy
    in fourteen symbols and the result is cached."""
    if not _CUBIC:
        import sympy
        ms = sympy.symbols("m0:14")
        Cs = sympy.Poly(sympy.expand(cubic_C(ms)), *ms)
        for mono, coef in Cs.terms():
            idx = tuple(i for i, e in enumerate(mono) for _ in range(e))
            assert len(idx) == 3 and coef.is_Integer
            _CUBIC[idx] = int(coef)
    return dict(_CUBIC)


_CUBIC = {}


def grad_terms():
    """For each j the dict {(i, k) sorted: coefficient} with dC/dm_j = sum coefficient m_i m_k."""
    out = [dict() for _ in range(14)]
    for idx, co in cubic_terms().items():
        for j in set(idx):
            rest = list(idx)
            mult = rest.count(j)
            rest.remove(j)
            key = tuple(rest)
            out[j][key] = out[j].get(key, 0) + co * mult
    return out


def hess_linear():
    """table[i][j] = list of (coefficient, k) with d^2 C / dm_i dm_j = sum coefficient m_k."""
    gt = grad_terms()
    table = [[[] for _ in range(14)] for _ in range(14)]
    for i in range(14):
        for (a, b), co in gt[i].items():
            # d/dm_j of co m_a m_b
            if a == b:
                table[i][a].append((2 * co, a))
            else:
                table[i][a].append((co, b))
                table[i][b].append((co, a))
    # merge
    for i in range(14):
        for j in range(14):
            acc = {}
            for co, k in table[i][j]:
                acc[k] = acc.get(k, 0) + co
            table[i][j] = sorted((co, k) for k, co in acc.items() if co)
    return table


def hess_at(u, zero=0):
    """Hess C(u) as a 14 x 14 list (entries linear in u; any ring whose elements can be
    multiplied by an int on the right).  `zero` is the ring's zero (used for empty entries)."""
    tab = hess_linear()
    out = []
    for i in range(14):
        row = []
        for j in range(14):
            v = None
            for co, k in tab[i][j]:
                t = u[k] * co
                v = t if v is None else v + t
            row.append(zero if v is None else v)
        out.append(row)
    return out


# ----------------------------------------------------------------------------------------
# 6. the 82 vectors f_j of Appendix D (through their values on the 24 permutations)
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
FIDX = {v: j for j, v in enumerate(F_VALUES)}
assert len(FIDX) == 82


def f_name(j):
    k, a, b = F_TABLE[j]
    return "n_{%s,%s}" % (a, b) if k == "n" else "h_{%s,%s}" % (a, "+" if b > 0 else "-")


def minor_spec(j):
    """(I, J) as 0-based tuples for a minor vector f_j."""
    kind, a, b = F_TABLE[j]
    assert kind == "n"
    return tuple(int(ch) - 1 for ch in a), tuple(int(ch) - 1 for ch in b)


def generators(T):
    """Indices j with f_j = 0 on every permutation index of T."""
    return [j for j in range(82) if all(F_VALUES[j][s] == 0 for s in T)]


def g_facet(gmap, j):
    """Index of the vector g.f_j (KeyError if G did not permute the 82 vectors)."""
    return FIDX[g_values(gmap, F_VALUES[j])]


def find_group_elements(rep, target):
    """All indices k of G with g_k . f_{rep[m]} = f_{target[m]} for every position m."""
    out = []
    for k, gm in enumerate(G_MAPS):
        if all(g_facet(gm, a) == b for a, b in zip(rep, target)):
            out.append(k)
    return out


def orbit_of_tuple(rep):
    return {tuple(g_facet(gm, a) for a in rep) for gm in G_MAPS}


def extreme_rays_of_contact_cone(T):
    """Extreme rays of N_T = {c : P c >= 0, (P c)_T = 0}, computed from the definition and
    without the list of 82 facets.  Returns (dim of the linear span of the equalities,
    sorted list of primitive integer value vectors (24-tuples) of the extreme rays).

    Method: N_T is a pointed polyhedral cone in the subspace L = {(Pc)_T = 0} of dimension
    d.  A nonzero c in N_T spans an extreme ray iff the inequalities tight at c have rank
    d - 1 on L.  So every extreme ray is the kernel of some d - 1 of the 24 - |T| remaining
    inequality rows; all such subsets are enumerated."""
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
# 7. Leibniz sums, minors, exact rational SO(4), quaternion matrices
# ----------------------------------------------------------------------------------------


def leibniz_weights(O):
    out = []
    for p, s in zip(PERMS, SGN):
        out.append(O[0][p[0]] * O[1][p[1]] * O[2][p[2]] * O[3][p[3]] * s)
    return out


def phi(vals, O, weights=None):
    """Phi(O) = sum_sigma sgn(sigma) vals[sigma] prod_i O[i][sigma(i)]."""
    w = leibniz_weights(O) if weights is None else weights
    tot = None
    for v, x in zip(vals, w):
        if v != 0:
            tot = x * v if tot is None else tot + x * v
    return 0 if tot is None else tot


def minor_det(O, I, J):
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


def det4(O):
    w = leibniz_weights(O)
    tot = w[0]
    for x in w[1:]:
        tot = tot + x
    return tot


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
# 8. sparse integer polynomials: a thin layer over python-flint (fmpz_mpoly)
#    python-flint 0.6 exposes only gens, + - * ** (polynomial on the left of an int),
#    len, coefficient(i), exponent_tuple(i), leading_coefficient, and a lex term order.
# ----------------------------------------------------------------------------------------


class PR:
    """Polynomial ring Z[names] (lex order, first name largest)."""

    def __init__(self, names):
        from flint import fmpz_mpoly
        self.names = list(names)
        self.n = len(self.names)
        self.gens = list(fmpz_mpoly.gens(self.n, "lex"))
        self.zero = self.gens[0] - self.gens[0]
        self.one = self.zero + 1
        self.var = dict(zip(self.names, self.gens))

    def const(self, c):
        return self.zero + int(c)

    def monomial(self, exps, coef=1):
        t = self.one
        for g, e in zip(self.gens, exps):
            if e:
                t = t * g ** int(e)
        return t * int(coef)

    def from_terms(self, terms):
        """terms: iterable of (exponent tuple, integer coefficient)."""
        out = self.zero
        chunk = self.zero
        k = 0
        for exps, co in terms:
            chunk = chunk + self.monomial(exps, co)
            k += 1
            if k % 256 == 0:
                out = out + chunk
                chunk = self.zero
        return out + chunk

    @staticmethod
    def terms(f):
        """list of (exponent tuple, python int coefficient)."""
        return [(tuple(int(e) for e in f.exponent_tuple(i)), int(f.coefficient(i)))
                for i in range(len(f))]

    @staticmethod
    def coeffs(f):
        return [int(f.coefficient(i)) for i in range(len(f))]

    @staticmethod
    def to_dict(f):
        return dict(PR.terms(f))

    def degrees(self, f):
        """Tuple of the degrees of f in each variable (zeros for the zero polynomial)."""
        d = [0] * self.n
        for i in range(len(f)):
            e = f.exponent_tuple(i)
            for k in range(self.n):
                if e[k] > d[k]:
                    d[k] = int(e[k])
        return tuple(d)

    def total_degree_set(self, f):
        return sorted({sum(int(x) for x in f.exponent_tuple(i)) for i in range(len(f))})

    def lead(self, f):
        """(exponent tuple, coefficient) of the lex-leading term of a nonzero f."""
        return tuple(int(e) for e in f.exponent_tuple(0)), int(f.coefficient(0))

    def coeff_in(self, f, k):
        """List c with f = sum_j c[j] * gens[k]**j, each c[j] free of gens[k]."""
        buckets = {}
        for exps, co in PR.terms(f):
            e = list(exps)
            j = e[k]
            e[k] = 0
            buckets.setdefault(j, []).append((tuple(e), co))
        top = max(buckets) if buckets else 0
        return [self.from_terms(buckets.get(j, [])) for j in range(top + 1)]

    def exact_div(self, f, g):
        """Exact quotient f / g by leading-term elimination (lex order).  Raises ValueError
        if g does not divide f over Z."""
        if len(g) == 0:
            raise ZeroDivisionError
        ge, gc = self.lead(g)
        rem = f
        qterms = []
        while len(rem) != 0:
            re_, rc = self.lead(rem)
            if any(a < b for a, b in zip(re_, ge)) or rc % gc != 0:
                raise ValueError("not exactly divisible")
            t = self.monomial(tuple(a - b for a, b in zip(re_, ge)), rc // gc)
            qterms.append(t)
            rem = rem - g * t
        out = self.zero
        for t in qterms:
            out = out + t
        return out

    def exact_sqrt(self, f):
        """A polynomial s with s*s = f (leading coefficient positive), by leading-term
        iteration in lex order.  Raises ValueError if f is not a square in Z[...]."""
        from math import isqrt
        if len(f) == 0:
            return self.zero
        fe, fc = self.lead(f)
        if fc <= 0 or isqrt(fc) ** 2 != fc or any(e % 2 for e in fe):
            raise ValueError("leading term is not a square")
        le, lc = tuple(e // 2 for e in fe), isqrt(fc)
        s = self.monomial(le, lc)
        rem = f - s * s
        while len(rem) != 0:
            re_, rc = self.lead(rem)
            if any(a < b for a, b in zip(re_, le)) or rc % (2 * lc) != 0:
                raise ValueError("not a perfect square")
            t = self.monomial(tuple(a - b for a, b in zip(re_, le)), rc // (2 * lc))
            rem = rem - (s * 2 + t) * t
            s = s + t
        return s

    def subs_clear(self, f, nums, dens):
        """prod_i dens[i]**deg_i(f) * f(nums[0]/dens[0], ...), as a polynomial of the ring
        in which nums, dens live.  Exact; the multiplier is a product of powers of the dens."""
        deg = self.degrees(f)
        npow = []
        dpow = []
        for n_, d_, m in zip(nums, dens, deg):
            a = [n_ - n_ + 1]
            b = [d_ - d_ + 1]
            for _ in range(m):
                a.append(a[-1] * n_)
                b.append(b[-1] * d_)
            npow.append(a)
            dpow.append(b)
        zero_t = nums[0] - nums[0]
        out = zero_t

        # nested (Horner-like) accumulation over the variables to keep products small
        def rec(terms, k):
            # terms: list of (exps, co) sharing the exponents of variables < k
            if k == self.n:
                tot = 0
                for _, co in terms:
                    tot += co
                return zero_t + tot
            groups = {}
            for exps, co in terms:
                groups.setdefault(exps[k], []).append((exps, co))
            acc = zero_t
            for j, sub in groups.items():
                acc = acc + rec(sub, k + 1) * (npow[k][j] * dpow[k][deg[k] - j])
            return acc
        return rec(PR.terms(f), 0)

    def digest(self, f):
        """sha256 of the canonical term list 'e1,e2,...:coef' sorted by exponent tuple."""
        lines = "\n".join(",".join(map(str, ex)) + ":" + str(co)
                          for ex, co in sorted(PR.terms(f)))
        return hashlib.sha256(lines.encode()).hexdigest()


def det_laplace(Mx, zero, one):
    """Determinant of a square matrix over a commutative ring by Laplace expansion along the
    rows with memoisation over column subsets (no division).  Rows are processed in order;
    zero entries are skipped."""
    n = len(Mx)
    cache = {0: one}
    # level r uses rows n-r .. n-1 ; subset = bitmask of columns used by those rows
    for r in range(1, n + 1):
        row = Mx[n - r]
        nz = [not ((isinstance(e, int) and e == 0) or (not isinstance(e, int) and len(e) == 0))
              for e in row]
        new = {}
        for cols in itertools.combinations(range(n), r):
            mask = 0
            for c in cols:
                mask |= 1 << c
            tot = None
            for pos, c in enumerate(cols):
                e = row[c]
                if not nz[c]:
                    continue
                sub = cache.get(mask ^ (1 << c))
                if sub is None:
                    continue
                term = sub * e
                if pos % 2:
                    term = -term
                tot = term if tot is None else tot + term
            if tot is not None and len(tot) != 0:
                new[mask] = tot
        cache = new
    return cache.get((1 << n) - 1, zero)


def sha256_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


# ----------------------------------------------------------------------------------------
# 9. Gaussian rationals, matching vectors (copied from ../cases-i-iv/replay_iv_common.py)
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
# 10. common checks (adapted from ../cases-i-iv/replay_iv_common.py)
# ----------------------------------------------------------------------------------------


def check_foundations(ck):
    """The objects of Sections 3 and 4 of the paper that the case relies on."""
    ck.section("Foundations rebuilt from the paper's definitions")
    ck.check("EXACT", "P is a 24 x 14 integer matrix; every row is a polynomial identity in "
             "16 variables (checked at construction); rank 14",
             len(P) == 24 and all(len(r) == 14 for r in P) and rank(P) == 14)
    ck.check("EXACT", "P kappa = (1,...,1) for kappa = e_13",
             all(P[s][13] == 1 for s in range(24)))
    pts = [(0, 1)] * 4 + [(1, 1)] * 4
    ck.check("EXACT", "kappa = Psi(v) at v_1..4 = (0,1), v_5..8 = (1,1), so kappa is in X",
             [eval_poly16(B_POLYS[r], pts) for r in range(14)] == [0] * 13 + [1])
    gt = grad_terms()
    bad = []
    BB = {}
    for j in range(14):
        tot = {}
        for (i, k), coef in gt[j].items():
            key = (min(i, k), max(i, k))
            if key not in BB:
                BB[key] = pmul(B_POLYS[i], B_POLYS[k])
            tot = padd(tot, BB[key], coef)
        if tot:
            bad.append(j)
    ck.check("EXACT", "dC/dm_j (Psi(v)) = 0 identically in 16 variables for j = 0..13 "
             "(so X lies in the singular locus of the cubic C of (A.2))", not bad,
             "nonzero: %s" % bad if bad else "")
    # the Hessian table is consistent with the cubic: m^t Hess(m) m = 6 C(m) on random points
    import random
    rng = random.Random(1)
    okh = True
    for _ in range(5):
        m = [Fraction(rng.randint(-5, 5)) for _ in range(14)]
        H = hess_at(m, Fraction(0))
        okh &= sum(m[i] * H[i][j] * m[j] for i in range(14) for j in range(14)) == 6 * cubic_C(m)
        okh &= all(H[i][j] == H[j][i] for i in range(14) for j in range(14))
    ck.check("EXACT", "own Hessian table of C: symmetric and m^t Hess C(m) m = 6 C(m) at 5 random "
             "integer points (the table is derived from the 37 monomials of (A.2))", okh)
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
    ok = len(set(G_MAPS)) == 1152
    for t in G_MAPS:
        for j in range(82):
            jj = FIDX.get(g_values(t, F_VALUES[j]))
            if jj is None or F_IS_MINOR[jj] != F_IS_MINOR[j]:
                ok = False
    ck.check("EXACT", "group G: 1152 distinct maps on the 24 permutation indices; G permutes the 82 "
             "vectors f_j and maps minor vectors to minor vectors", ok)


def check_radial_necessity(ck, seed=20261125, samples=12):
    """7.3 of the audit's PROOF.md / p. 27 of the paper: x = u + i v in X  =>  C(u) = 0 and
    Hess C(u) v = 0."""
    import random
    import sympy
    ck.section("Radial necessity (PROOF.md 7.3; paper p. 27): the real part of a point of X")
    ms = sympy.symbols("m0:14")
    C = sympy.expand(cubic_C(ms))
    grad = [sympy.diff(C, v) for v in ms]
    hess = [[sympy.diff(g, v) for v in ms] for g in grad]
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
    eul2 = all(sympy.expand(sum(hess[j][k].xreplace(su) * us[k] for k in range(14))
                            - 2 * grad[j].xreplace(su)) == 0 for j in range(14))
    ck.check("EXACT", "Hess C(u) u = 2 grad C(u) identically (used when Im x = 0)", eul2)
    ck.info("hence grad C(u+iv) = 0 gives grad C(u) = grad C(v), Hess C(u) v = 0, and "
            "3 C(u) = u.grad C(u) = u.grad C(v) = (1/2) v^t Hess C(u) v = 0; Hess C(u) is singular "
            "whenever u != 0 (kernel vector v if v != 0, else u)")
    rng = random.Random(seed)
    okc = okh = True
    for _ in range(samples):
        def rq():
            return GQ(Fraction(rng.randint(-9, 9), rng.randint(1, 5)),
                      Fraction(rng.randint(-9, 9), rng.randint(1, 5)))
        a = [rq() for _ in range(4)]
        b = [rq() for _ in range(4)]
        al = rq()
        x = [GQ.of(c) for c in matching_vector(a, b, al)]
        u = [c.re for c in x]
        v = [c.im for c in x]
        okc &= cubic_C(u) == 0
        H = hess_at(u, Fraction(0))
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
        ph = 0
        for s, val in zip(PERMS, F_VALUES[j]):
            if val:
                ph += sgn(s) * M[0, s[0]] * M[1, s[1]] * M[2, s[2]] * M[3, s[3]]
        dIJ = M.extract(I, J).det()
        dc = M.extract(Ic, Jc).det()
        ok &= sympy.expand(ph - eps * dIJ * dc) == 0
        dadj = adj.extract(J, I).det()
        ok &= sympy.expand(dadj - eps * dc * detM ** (len(I) - 1)) == 0
    ck.check("EXACT", "Laplace and Jacobi identities for the %d minor generators %s: "
             "Phi_{f_j}(O) = (minor of O)^2 >= 0 on SO(4)" % (len(minor_indices), list(minor_indices)), ok)


def crosscheck_with_audit_library(ck, nu, T):
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
    import random
    rng = random.Random(7)
    okc = True
    for _ in range(20):
        m = [Fraction(rng.randint(-6, 6), rng.randint(1, 3)) for _ in range(14)]
        okc &= L.cubic_C(m) == cubic_C(m)
    # a cubic form in 14 variables is determined by its values on {0,1,2,3}^14-type grids; the
    # random test is a comparison only, both sides being transcriptions of the same formula (A.2)
    ck.check("EXACT", "the cubic C agrees with moc4lib.cubic_C on 20 random rational points "
             "(comparison of two transcriptions of (A.2))", okc)
    ck.check("EXACT", "T_%d agrees with moc4lib.CONTACT_SETS[%d]; generators agree; the group has "
             "the same 1152 maps" % (nu, nu),
             tuple(sorted(L.CONTACT_SETS[nu])) == tuple(sorted(T))
             and list(L.generators(set(T))) == generators(T)
             and set(L._G_MAPS) == set(G_MAPS))
