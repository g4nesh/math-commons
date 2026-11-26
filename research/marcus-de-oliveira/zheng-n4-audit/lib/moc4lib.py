"""moc4lib: exact objects for the audit of Zheng, "The Marcus-de Oliveira conjecture for n = 4".

Every object is derived from its DEFINITION in the paper (Sections 3.1-3.3, 4.1-4.2,
Appendices A, D, E, G.1); nothing is read from or computed by the supplement.
Exact arithmetic uses Python ints / fractions.Fraction, a small Gaussian-rational class QI,
and sympy for the cubic form.  numpy/scipy appear only in the section marked NUMERIC.

Index conventions (see README.md):
  * Permutations sigma of {0,1,2,3} are tuples with sigma[i] = image of i (0-based).
    PERMS lists all 24 in lexicographic order of one-line notation; PERMS[0] = identity.
    The paper's pi_j (1-based one-line notation) is PERMS[j] with every entry minus 1.
    The paper's P_sigma has a 1 in row i, column sigma(i); a matrix O enters the Leibniz
    sum through prod_i O[i][sigma[i]].
  * compose(p, q) = p o q, i.e. (p o q)(i) = p(q(i)).
  * Bracket labels are 1..8 as in the paper; [ij] = det(v_i, v_j) = x_i y_j - x_j y_i.
  * Vectors c in V_R are lists of 14 coordinates m_0..m_13 (m_j = c(B_j)).
  * Index sets I, J of minors are 1-based tuples, as in the paper.
"""
from fractions import Fraction
import itertools
import random

# =====================================================================================
# 1. Permutations
# =====================================================================================
PERMS = [tuple(p) for p in itertools.permutations(range(4))]   # lexicographic
_PIDX = {p: i for i, p in enumerate(PERMS)}


def perm_index(p):
    """Lexicographic index (0..23) of the 0-based permutation tuple p."""
    return _PIDX[tuple(p)]


def sgn(p):
    inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
    return -1 if inv % 2 else 1


def compose(p, q):
    """p o q : i -> p[q[i]]."""
    return tuple(p[q[i]] for i in range(len(q)))


def inverse(p):
    r = [0] * len(p)
    for i, pi in enumerate(p):
        r[pi] = i
    return tuple(r)


def perm_from_oneline(s):
    """'2143' (paper's 1-based one-line notation) -> (1,0,3,2)."""
    return tuple(int(ch) - 1 for ch in str(s))


def perm_to_oneline(p):
    return "".join(str(x + 1) for x in p)


SGN = [sgn(p) for p in PERMS]
IDENTITY = PERMS[0]

# =====================================================================================
# 2. Exact linear algebra over Fraction (small helpers)
# =====================================================================================


def _F(x):
    """Coerce int -> Fraction; leave Fraction and Gaussian rationals (QI) unchanged."""
    if isinstance(x, Fraction) or type(x).__name__ == "QI":
        return x
    return Fraction(x)


def rref(rows):
    """Reduced row echelon form over Q. Returns (rows, pivot_columns)."""
    A = [[_F(x) for x in r] for r in rows]
    if not A:
        return A, []
    nr, nc = len(A), len(A[0])
    piv, r = [], 0
    for c in range(nc):
        k = next((i for i in range(r, nr) if A[i][c] != 0), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
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


def det(M):
    A = [[_F(x) for x in r] for r in M]
    n, d = len(A), Fraction(1)
    for c in range(n):
        k = next((i for i in range(c, n) if A[i][c] != 0), None)
        if k is None:
            return Fraction(0)
        if k != c:
            A[c], A[k] = A[k], A[c]
            d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return d


def inv(M):
    n = len(M)
    aug = [list(M[i]) + [1 if j == i else 0 for j in range(n)] for i in range(n)]
    R, piv = rref(aug)
    if piv[:n] != list(range(n)):
        raise ValueError("singular matrix")
    return [r[n:] for r in R]


def nullspace(M):
    """Basis of {x : M x = 0} over Q."""
    R, piv = rref(M)
    nc = len(M[0])
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for f in free:
        x = [Fraction(0)] * nc
        x[f] = Fraction(1)
        for r, p in enumerate(piv):
            x[p] = -R[r][f]
        basis.append(x)
    return basis


def matmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), 0) for j in range(len(B[0]))]
            for i in range(len(A))]


def matvec(A, v):
    return [sum((a * x for a, x in zip(row, v)), 0) for row in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def eye(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def det_leibniz4(M):
    """4x4 determinant by the Leibniz formula; works over any commutative ring."""
    tot = 0
    for p, s in zip(PERMS, SGN):
        t = M[0][p[0]] * M[1][p[1]] * M[2][p[2]] * M[3][p[3]]
        tot = tot + t if s > 0 else tot - t
    return tot


# =====================================================================================
# 3. Gaussian rationals (exact complex numbers)
# =====================================================================================
class QI:
    """Exact Gaussian rational re + i*im with Fraction parts."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = _F(re), _F(im)

    @staticmethod
    def _c(o):
        return o if isinstance(o, QI) else QI(o, 0)

    def __add__(self, o):
        o = QI._c(o); return QI(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __sub__(self, o):
        o = QI._c(o); return QI(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return QI._c(o) - self

    def __neg__(self):
        return QI(-self.re, -self.im)

    def __mul__(self, o):
        o = QI._c(o)
        return QI(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = QI._c(o); n = o.re * o.re + o.im * o.im
        return self * QI(o.re / n, -o.im / n)

    def __eq__(self, o):
        try:
            o = QI._c(o)
        except Exception:
            return NotImplemented
        return self.re == o.re and self.im == o.im

    def __hash__(self):
        return hash((self.re, self.im))

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            return NotImplemented
        r = QI(1)
        for _ in range(n):
            r = r * self
        return r

    def conjugate(self):
        return QI(self.re, -self.im)

    def __complex__(self):
        return complex(float(self.re), float(self.im))

    def __repr__(self):
        return f"QI({self.re}, {self.im})"


def re_part(x):
    """Real part of int / Fraction / QI / complex / sympy number (exact where possible)."""
    if isinstance(x, QI):
        return x.re
    if isinstance(x, (int, Fraction)):
        return _F(x)
    if isinstance(x, complex):
        return x.real
    try:
        import sympy
        if isinstance(x, sympy.Basic):
            return sympy.re(x)
    except ImportError:
        pass
    return x.real


# =====================================================================================
# 4. Bracket polynomials, the matching basis B_0..B_13 (A.1) and the matrix P
# =====================================================================================
# An element of R_1 (multi-degree (1,...,1) in eight vectors v_l = (x_l, y_l)) is stored
# as a dict  mask -> integer coefficient, where bit (l-1) of mask set means the factor
# x_l (else y_l) occurs in the monomial.  This represents the polynomial in the sixteen
# variables x_1..x_8, y_1..y_8 EXACTLY (no dehomogenization), so an equality of dicts is
# a polynomial identity in eight symbolic 2-vectors.  Setting y_l = 1, x_l = z_l gives the
# dehomogenized form used in (4.5): [ij] = z_i - z_j; this map is a bijection on
# multi-degree (1,...,1) polynomials.


def bracket(i, j):
    """[ij] = x_i y_j - x_j y_i (labels 1..8)."""
    if i == j:
        return {}
    return {1 << (i - 1): 1, 1 << (j - 1): -1}


def _mul_disjoint(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            assert a & b == 0
            out[a | b] = out.get(a | b, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def matching_poly(pairs):
    """Product of brackets [i j] over the given pairs (each label 1..8 exactly once)."""
    labels = sorted(x for pr in pairs for x in pr)
    assert labels == list(range(1, 9)), pairs
    p = {0: 1}
    for i, j in pairs:
        p = _mul_disjoint(p, bracket(i, j))
    return p


def poly_add(p, q, cq=1):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) + cq * v
    return {k: v for k, v in out.items() if v}


def poly_relabel(p, perm8):
    """Polynomial of the relabelled configuration: v_l replaced by v_{perm8[l]}.
    perm8 is a dict on labels 1..8.  (B o t)(v) = B(t v), where (t v)_l = v_{t(l)}."""
    out = {}
    for mask, c in p.items():
        nm = 0
        for l in range(1, 9):
            if mask >> (l - 1) & 1:
                nm |= 1 << (perm8[l] - 1)
        out[nm] = out.get(nm, 0) + c
    return out


# (A.1), transcribed from the page image (p. 35) and checked against the text extract.
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
B_POLYS = [matching_poly(m) for m in MATCHINGS]


def is_noncrossing(pairs):
    for (a, b), (c, d) in itertools.combinations(pairs, 2):
        a, b = sorted((a, b)); c, d = sorted((c, d))
        if a < c < b < d or c < a < d < b:
            return False
    return True


def _basis_setup():
    masks = sorted({m for p in B_POLYS for m in p})
    rows = [[p.get(m, 0) for m in masks] for p in B_POLYS]   # 14 x |masks|
    # choose 14 monomials (pivot columns) on which the B's are independent
    Rm, pm = rref(rows)
    if len(pm) != 14:
        raise AssertionError("B_0..B_13 not linearly independent")
    sel = [masks[c] for c in pm]
    Msel = [[rows[a][masks.index(m)] for m in sel] for a in range(14)]   # 14x14, [a][k]
    return sel, inv(transpose(Msel))   # coeffs = inv(Msel^T) * target[sel]


_SEL_MASKS, _SEL_INV = _basis_setup()


def express_in_basis(poly):
    """Exact coefficients c_a with poly = sum_a c_a B_a as polynomials in 16 variables.
    Raises ValueError if poly is not in the span (the identity is checked on EVERY
    monomial, so a successful return is a proof of the polynomial identity)."""
    t = [poly.get(m, 0) for m in _SEL_MASKS]
    c = matvec(_SEL_INV, t)
    recon = {}
    for a, ca in enumerate(c):
        for m, v in B_POLYS[a].items():
            recon[m] = recon.get(m, 0) + ca * v
    recon = {k: v for k, v in recon.items() if v != 0}
    target = {k: v for k, v in poly.items() if v != 0}
    if recon != target:
        raise ValueError("polynomial not in span of B_0..B_13")
    return c


def p_sigma_poly(sigma):
    """p_sigma = prod_i [i, 4 + sigma(i)] (paper labels; sigma 0-based tuple)."""
    return matching_poly([(i + 1, 5 + sigma[i]) for i in range(4)])


def _build_P():
    P = []
    for s in PERMS:
        c = express_in_basis(p_sigma_poly(s))
        assert all(x.denominator == 1 for x in c)
        P.append([int(x) for x in c])
    return P


P = _build_P()                                   # 24 x 14 integer matrix, rows = PERMS order
P_PIVOT_ROWS = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18]   # Appendix A
_P_PIV_INV = inv([P[r] for r in P_PIVOT_ROWS])
KAPPA = [Fraction(0)] * 13 + [Fraction(1)]       # (3.2): m(kappa) = e_13


def values(c):
    """(P c)_sigma for all 24 sigma, i.e. v_sigma(c)."""
    return [sum((P[s][a] * c[a] for a in range(14)), 0) for s in range(24)]


def solve_P(vals):
    """Unique c in Q^14 with P c = vals (24 values).  Raises ValueError if vals is not
    in the image of P (all 24 rows are checked)."""
    c = matvec(_P_PIV_INV, [_F(vals[r]) for r in P_PIVOT_ROWS])
    if values(c) != [_F(v) for v in vals]:
        raise ValueError("value vector not in image of P")
    return c


# =====================================================================================
# 5. The cubic C of (A.2) and the transposition matrices T_i
# =====================================================================================


def cubic_C(m):
    """C(m) of (A.2), transcribed from the page image (p. 36). Works on any ring."""
    m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13 = m
    return (m13 * (m0 * m6 - m1 * m5) - m0 * m9 * (m0 + m1 + m2 + m3 + m5 + m7 + m9)
            - m1 * m5 * m9
            + m10 * (-m0 * m3 - m0 * m9 + m1 * m2 + m1 * m7)
            + m11 * (-m0 * m7 - m0 * m9 + m2 * m5 + m3 * m5)
            + m12 * (-m0 * m9 + m1 * m7 - m2 * m6 - m2 * m9 + m3 * m5 + m3 * m7)
            + m2 * m10 * m11 - m3 * m8 * m10 - m4 * m7 * m11 + m4 * m8 * m9)


def cubic_sympy():
    """(symbols m0..m13, sympy expression of C)."""
    import sympy
    ms = sympy.symbols("m0:14")
    return ms, sympy.expand(cubic_C(ms))


# Quadratic monomials listed after (A.2), as index pairs, in the paper's order.
C_PARTIAL_MONOMIALS = [(10, 11), (9, 12), (8, 10), (8, 9), (7, 12), (4, 11), (4, 9), (4, 8),
                       (4, 7), (3, 12), (3, 8), (3, 7), (2, 12), (1, 5)]


def transposition_matrix(i):
    """T_i for the adjacent label transposition t_i = (i i+1), 1 <= i <= 7:
    B_r(t_i v) = sum_s T_i[r][s] B_s(v), so m(t_i v) = T_i m(v) for matching vectors."""
    t = {l: l for l in range(1, 9)}
    t[i], t[i + 1] = i + 1, i
    return [[int(x) for x in express_in_basis(poly_relabel(B_POLYS[r], t))] for r in range(14)]


def hess_kappa():
    """H_kappa = e0 e6^t + e6 e0^t - e1 e5^t - e5 e1^t as claimed in (3.10)."""
    H = [[0] * 14 for _ in range(14)]
    H[0][6] = H[6][0] = 1
    H[1][5] = H[5][1] = -1
    return H


# =====================================================================================
# 6. Entry subspace (Appendix A.2, eq. (3.7))
# =====================================================================================


def E_coords(Q):
    """Coordinates of E(Q) from Appendix A.2; Q is 4x4 (0-based), q_ij = Q[i-1][j-1]."""
    q = lambda i, j: Q[i - 1][j - 1]
    z = q(1, 1) * 0
    return [z, z,
            q(1, 3) - q(1, 4) - q(2, 3) + q(2, 4),          # m2
            q(1, 2) - q(1, 3) - q(2, 2) + q(2, 3),          # m3
            q(1, 1) - q(1, 2) - q(2, 1) + q(2, 2),          # m4
            z, z,
            q(2, 3) - q(2, 4) - q(3, 3) + q(3, 4),          # m7
            q(3, 3) - q(3, 4) - q(4, 3) + q(4, 4),          # m8
            q(2, 2) - q(2, 3) - q(3, 2) + q(3, 3),          # m9
            q(2, 1) - q(2, 2) - q(3, 1) + q(3, 2),          # m10
            q(3, 2) - q(3, 3) - q(4, 2) + q(4, 3),          # m11
            q(3, 1) - q(3, 2) - q(4, 1) + q(4, 2),          # m12
            q(1, 4) + q(2, 3) + q(3, 2) + q(4, 1)]          # m13


# =====================================================================================
# 7. Minor vectors, t, d, h and the facet list f_0..f_81 (Section 4.1, Appendix D)
# =====================================================================================


def n_values(I, J):
    """Value vector v_sigma(n_{I,J}) = 1{sigma(I) = J}; I, J 1-based tuples."""
    I0 = {i - 1 for i in I}; J0 = {j - 1 for j in J}
    return [1 if {s[i] for i in I0} == J0 else 0 for s in PERMS]


def n_vec(I, J):
    """n_{I,J} in Q^14 (eq. (4.1)), solved exactly and checked on all 24 rows."""
    return solve_P(n_values(I, J))


def _comp(I):
    return tuple(sorted(set(range(1, 5)) - set(I)))


# 16 singletons and 18 classes of 2-sets (representative with I in {12,13,14}).
MINOR_KEYS = ([((i,), (j,)) for i in range(1, 5) for j in range(1, 5)] +
              [(I, J) for I in [(1, 2), (1, 3), (1, 4)]
               for J in itertools.combinations(range(1, 5), 2)])


def canonical_minor_key(I, J):
    I, J = tuple(sorted(I)), tuple(sorted(J))
    if len(I) == 2 and 1 not in I:
        I, J = _comp(I), _comp(J)
    return (I, J)


def t_vec(pi):
    out = [Fraction(0)] * 14
    for i in range(4):
        out = [a + b for a, b in zip(out, n_vec((i + 1,), (pi[i] + 1,)))]
    return out


def d_vec(pi):
    p = [x + 1 for x in pi]   # 1-based images pi(1..4)
    terms = [(+1, (1, 2), (p[0], p[2])), (+1, (1, 2), (p[1], p[3])),
             (-1, (1, 3), (p[0], p[1])), (-1, (1, 3), (p[2], p[3]))]
    out = [Fraction(0)] * 14
    for s, I, J in terms:
        out = [a + s * b for a, b in zip(out, n_vec(I, tuple(sorted(J))))]
    return out


def h_vec(pi, eps):
    return [(a + eps * b) / 2 for a, b in zip(t_vec(pi), d_vec(pi))]


# Appendix D, transcribed from the page image (p. 59) and cross-checked with the text.
# ('n', I, J) with digit strings, or ('h', one-line pi, eps).
F_TABLE = [
    ("n", "1", "4"), ("n", "12", "34"), ("n", "13", "34"), ("n", "14", "34"),
    ("n", "1", "3"), ("n", "12", "24"), ("n", "4", "1"), ("n", "3", "1"),
    ("n", "12", "23"), ("n", "14", "23"), ("n", "13", "24"), ("n", "13", "23"),
    ("n", "14", "24"), ("n", "2", "1"), ("n", "1", "2"), ("n", "12", "14"),
    ("n", "2", "4"), ("h", "3421", +1), ("h", "3412", -1), ("h", "2413", -1),
    ("h", "2431", +1), ("n", "4", "2"), ("h", "3412", +1), ("n", "14", "12"),
    ("h", "4312", +1), ("h", "3142", -1), ("h", "4132", -1), ("n", "3", "2"),
    ("h", "3421", -1), ("n", "13", "12"), ("h", "4321", -1), ("h", "4123", +1),
    ("h", "3124", +1), ("n", "12", "13"), ("h", "4321", +1), ("h", "4312", -1),
    ("n", "2", "3"), ("h", "2341", +1), ("h", "2314", -1), ("n", "14", "13"),
    ("h", "4213", +1),                                                   # j = 0..40
    ("n", "4", "3"), ("h", "4123", -1), ("h", "2413", +1), ("h", "2143", -1),
    ("n", "13", "14"), ("h", "3241", -1), ("n", "3", "4"), ("h", "3142", +1),
    ("h", "2143", +1), ("h", "2341", -1), ("h", "1423", -1), ("h", "1342", +1),
    ("n", "13", "13"), ("h", "4231", -1), ("n", "3", "3"), ("h", "4132", +1),
    ("h", "2431", -1), ("h", "2134", +1), ("n", "14", "14"), ("h", "3214", +1),
    ("n", "4", "4"), ("h", "3124", -1), ("h", "2134", -1), ("h", "2314", +1),
    ("h", "1432", +1), ("h", "1324", -1), ("n", "2", "2"), ("h", "4213", -1),
    ("h", "3241", +1), ("h", "4231", +1), ("h", "3214", -1), ("n", "12", "12"),
    ("h", "1243", +1), ("h", "1432", -1), ("h", "1423", +1), ("h", "1234", -1),
    ("h", "1342", -1), ("h", "1234", +1), ("h", "1243", -1), ("h", "1324", +1),
    ("n", "1", "1"),                                                     # j = 41..81
]
assert len(F_TABLE) == 82


def f_spec_vec(spec):
    kind, a, b = spec
    if kind == "n":
        return n_vec(tuple(int(ch) for ch in a), tuple(int(ch) for ch in b))
    return h_vec(perm_from_oneline(a), b)


F_VECS = [f_spec_vec(s) for s in F_TABLE]
F_VALUES = [tuple(int(v) for v in values(f)) for f in F_VECS]      # all values are integers
F_ZERO = [frozenset(s for s in range(24) if F_VALUES[j][s] == 0) for j in range(82)]
F_IS_MINOR = [s[0] == "n" for s in F_TABLE]


def f_name(j):
    k, a, b = F_TABLE[j]
    return f"n{a},{b}" if k == "n" else f"h{a},{'+' if b > 0 else '-'}"


# =====================================================================================
# 8. Leibniz sums
# =====================================================================================


def leibniz_weights(O):
    """w_sigma = sgn(sigma) prod_i O[i][sigma(i)] (Definition 3.2 on SO(4))."""
    out = []
    for p, s in zip(PERMS, SGN):
        t = O[0][p[0]] * O[1][p[1]] * O[2][p[2]] * O[3][p[3]]
        out.append(t if s > 0 else -t)
    return out


def Phi(c, O):
    """Phi_c(O) = sum_sigma sgn(sigma) (Pc)_sigma prod_i O_{i,sigma(i)} (3.4), exact."""
    return sum((v * w for v, w in zip(values(c), leibniz_weights(O))), 0)


def Phi_values(vals, O):
    """Same as Phi but taking the 24 values (P c) directly."""
    return sum((v * w for v, w in zip(vals, leibniz_weights(O))), 0)


def Phi_O4(c, O):
    """Extension to O(4) from Section 4.2: divides the weights by det O."""
    return Phi(c, O) / det_leibniz4(O)


def submatrix(O, I, J):
    return [[O[i - 1][j - 1] for j in J] for i in I]


def minor_det(O, I, J):
    M = submatrix(O, I, J)
    if len(I) == 1:
        return M[0][0]
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


O_STAR = [[Fraction(x, 7) for x in r] for r in
          [[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]]   # Example 1.3


def cayley(S):
    """(I - S)(I + S)^{-1} for a rational skew-symmetric 4x4 S: exactly in SO(4)."""
    I4 = eye(4)
    A = [[I4[i][j] - S[i][j] for j in range(4)] for i in range(4)]
    B = [[I4[i][j] + S[i][j] for j in range(4)] for i in range(4)]
    return matmul(A, inv(B))


def plane_rotation(i, j, p, q):
    """Rational rotation in coordinate plane (i, j) (0-based) by the angle with
    cos = (p^2 - q^2)/(p^2 + q^2), sin = 2pq/(p^2 + q^2)."""
    n = p * p + q * q
    c, s = Fraction(p * p - q * q, n), Fraction(2 * p * q, n)
    R = eye(4)
    R[i][i], R[j][j], R[i][j], R[j][i] = c, c, -s, s
    return R


def random_rational_SO4(rng, size=4):
    """Random exactly-special-orthogonal rational 4x4: a Cayley transform times one or
    two rational plane rotations.  Entries are re-drawn until all 16 are nonzero."""
    while True:
        S = [[Fraction(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i + 1, 4):
                v = Fraction(rng.randint(-size, size), rng.randint(1, size))
                S[i][j], S[j][i] = v, -v
        O = cayley(S)
        for _ in range(rng.randint(0, 2)):
            i, j = rng.sample(range(4), 2)
            O = matmul(O, plane_rotation(i, j, rng.randint(1, 5), rng.randint(1, 5)))
        if all(x != 0 for r in O for x in r):
            return O


def is_special_orthogonal(O):
    return matmul(O, transpose(O)) == eye(4) and det(O) == 1


# =====================================================================================
# 9. Spectral side (Section 3.1-3.2)
# =====================================================================================


def eval_poly(poly, z):
    """Evaluate an R_1 polynomial at v_l = (z_l, 1); z is a list of 8 ring elements
    (index 0 = label 1)."""
    tot = 0
    for mask, c in poly.items():
        t = c
        for l in range(8):
            if mask >> l & 1:
                t = t * z[l]
        tot = tot + t
    return tot


def eval_poly_hom(poly, pts):
    """Evaluate at eight arbitrary 2-vectors pts[l] = (x, y)."""
    tot = 0
    for mask, c in poly.items():
        t = c
        for l in range(8):
            t = t * (pts[l][0] if mask >> l & 1 else pts[l][1])
        tot = tot + t
    return tot


def matching_vector(a, b, alpha=1):
    """m in C^14 with m_r = alpha * B_r(v), v_i = (a_i, 1), v_{4+j} = (-b_j, 1).
    Then (P m)_sigma = alpha * prod_i (a_i + b_{sigma(i)}).  Exact for int/Fraction/QI/sympy."""
    z = list(a) + [-x for x in b]
    return [alpha * eval_poly(B_POLYS[r], z) for r in range(14)]


def matching_vector_from_points(pts, alpha=1):
    return [alpha * eval_poly_hom(B_POLYS[r], pts) for r in range(14)]


def sigma_points(a, b):
    """z_sigma = prod_i (a_i + b_{sigma(i)}) in PERMS order."""
    out = []
    for s in PERMS:
        t = 1
        for i in range(4):
            t = t * (a[i] + b[s[i]])
        out.append(t)
    return out


def gap_vector(m, pi_star):
    """c_{pi*}(m) = Re(m - F_{pi*}(m) kappa) (3.3).  pi_star: index or tuple."""
    k = pi_star if isinstance(pi_star, int) else perm_index(pi_star)
    F = sum((P[k][a] * m[a] for a in range(14)), 0)
    c = [re_part(x) for x in m]
    c[13] = c[13] - re_part(F)
    return c


def contact_set(c, require_nonneg=True):
    """Indices sigma with (Pc)_sigma == 0.  Raises if some value is negative and
    require_nonneg (the contact set is only defined when P c >= 0)."""
    v = values(c)
    if require_nonneg and any(x < 0 for x in v):
        raise ValueError("P c has a negative entry")
    return frozenset(s for s in range(24) if v[s] == 0)


def det_Da_plus_ODbOt(a, b, O):
    """det(D_a + O D_b O^t), exact over the ring of a, b, O."""
    M = [[(a[i] if i == j else 0) + sum((O[i][k] * b[k] * O[j][k] for k in range(4)), 0)
          for j in range(4)] for i in range(4)]
    return det_leibniz4(M)


# =====================================================================================
# 10. Symmetry group G (Section 4.2, Table 3)
# =====================================================================================
# A group element is g = (rho, rho2, eps) with rho, rho2 0-based permutation tuples and
# eps in {0, 1}.  On permutations:  g.sigma = rho o sigma o rho2           (eps = 0)
#                                   g.sigma = rho o sigma^{-1} o rho2      (eps = 1)
# with (p o q)(i) = p(q(i)).  This is the paper's A_{rho,rho',eps} (Section 5.2), and it
# reproduces all seven rows of Table 3 (tested).  On matrices (Section 4.2):
#   eps = 0: g.O = R_{rho2}^{-1} O R_rho^{-1}, i.e. (g.O)[i][j] = O[rho2(i)][rho^{-1}(j)]
#   eps = 1: first O -> O^t, then the eps = 0 formula.
# On V_R: g.c is the unique vector with v_{g.sigma}(g.c) = v_sigma(c).
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_act_perm(g, sigma):
    r1, r2, e = g
    s = inverse(sigma) if e else tuple(sigma)
    return compose(r1, compose(s, r2))


def g_index_map(g):
    """Tuple t with t[j] = index of g.pi_j."""
    return tuple(perm_index(g_act_perm(g, s)) for s in PERMS)


_G_MAPS = [g_index_map(g) for g in G_ELEMENTS]


def g_act_values(g, vals):
    t = g_index_map(g)
    out = [None] * 24
    for j in range(24):
        out[t[j]] = vals[j]
    return out


def g_act_c(g, c):
    return solve_P(g_act_values(g, values(c)))


def g_act_O(g, O):
    r1, r2, e = g
    if e:
        O = transpose(O)
    r1i = inverse(r1)
    return [[O[r2[i]][r1i[j]] for j in range(4)] for i in range(4)]


def g_act_spectra(g, a, b):
    """Spectra (a', b') whose sigma-points satisfy F'_{g.sigma} = F_sigma, so that
    matching_vector(a', b', alpha) = g . matching_vector(a, b, alpha) (tested).
    eps = 0: a' = a o rho2, b' = b o rho^{-1};  eps = 1: swap the roles of a and b first."""
    r1, r2, e = g
    if e:
        a, b = b, a
    r1i = inverse(r1)
    return [a[r2[i]] for i in range(4)], [b[r1i[j]] for j in range(4)]


def g_act_set(g, T):
    t = g_index_map(g) if not isinstance(g, int) else _G_MAPS[g]
    return frozenset(t[j] for j in T)


def orbit_rep(T):
    """Lexicographically smallest sorted tuple among the images of T under G."""
    T = list(T)
    return min(tuple(sorted(t[j] for j in T)) for t in _G_MAPS)


def find_g(T, T2):
    """Some g in G with g.T = T2 (as sets), or None."""
    T2 = frozenset(T2)
    for g, t in zip(G_ELEMENTS, _G_MAPS):
        if frozenset(t[j] for j in T) == T2:
            return g
    return None


# =====================================================================================
# 11. Contact sets (Appendix E) and seven-element auxiliary sets (Table 5)
# =====================================================================================
# Transcribed from the page images (pp. 60, 62-63) and cross-checked with the text.
CONTACT_SETS = {
    1: {0, 1, 2, 3, 6, 7}, 2: {0, 1, 2, 3, 6, 8}, 3: {0, 1, 2, 3, 7, 8},
    4: {0, 1, 4, 5, 7, 8}, 5: {1, 2, 3, 4, 6, 8}, 6: {0, 1, 2, 4, 6, 8},
    7: {0, 1, 2, 4, 7, 8}, 8: {0, 2, 5, 6, 7, 8}, 9: {0, 1, 2, 3, 11, 12},
    10: {0, 1, 3, 5, 7, 8}, 11: {0, 1, 3, 5, 6, 8}, 12: {0, 2, 3, 5, 9, 12},
    13: {1, 2, 3, 4, 8, 12}, 14: {0, 2, 3, 5, 11, 12}, 15: {1, 2, 3, 4, 6, 7},
    16: {0, 2, 5, 6, 9, 12}, 17: {0, 1, 4, 5, 9, 12}, 18: {0, 1, 3, 6, 7, 8},
    19: {0, 1, 4, 6, 7, 8}, 20: {0, 1, 5, 6, 7, 8}, 21: {0, 3, 5, 9, 11, 12},
    22: {2, 3, 4, 8, 10, 12}, 23: {2, 3, 5, 9, 11, 12}, 24: {2, 3, 4, 6, 7, 8},
    25: {3, 4, 7, 8, 10, 12}, 26: {2, 3, 5, 7, 8, 12}, 27: {2, 3, 4, 7, 8, 12},
    28: {2, 3, 5, 6, 9, 12}, 29: {0, 4, 5, 6, 7, 8}, 30: {2, 5, 9, 12, 13, 14},
    31: {0, 1, 5, 6, 9, 12}, 32: {0, 2, 3, 7, 11, 12}, 33: {0, 3, 6, 7, 11, 12},
    34: {2, 3, 5, 9, 10, 12}, 35: {2, 3, 4, 6, 9, 12}, 36: {2, 3, 4, 6, 11, 12},
    37: {0, 5, 6, 7, 9, 12}, 38: {1, 2, 3, 8, 11, 12}, 39: {2, 3, 5, 7, 9, 12},
    40: {0, 1, 5, 9, 14, 18}, 41: {0, 1, 3, 6, 11, 12}, 42: {1, 2, 3, 6, 7, 8},
    43: {1, 2, 4, 6, 7, 8}, 44: {1, 4, 5, 6, 7, 8}, 45: {1, 2, 3, 7, 11, 12},
    46: {3, 4, 5, 9, 10, 12}, 47: {1, 3, 5, 6, 7, 8}, 48: {0, 2, 5, 9, 12, 13},
    49: {0, 3, 5, 6, 9, 12}, 50: {1, 3, 5, 6, 10, 11}, 51: {1, 3, 5, 10, 11, 12},
    52: {1, 3, 5, 10, 11, 14}, 53: {0, 1, 4, 6, 9, 12}, 54: {2, 3, 5, 6, 11, 12},
    55: {1, 2, 3, 6, 11, 12}, 56: {3, 4, 5, 8, 9, 12}, 57: {0, 1, 5, 7, 9, 12},
    58: {0, 3, 5, 7, 9, 12}, 59: {3, 4, 5, 6, 9, 12}, 60: {2, 3, 5, 7, 11, 12},
    61: {1, 2, 3, 9, 11, 12}, 62: {3, 5, 7, 8, 9, 12}, 63: {1, 2, 3, 10, 11, 12},
    64: {0, 3, 5, 9, 14, 18}, 65: {0, 4, 5, 6, 9, 12}, 66: {2, 4, 5, 6, 9, 12},
    67: {2, 4, 5, 7, 9, 12}, 68: {1, 4, 5, 7, 9, 12}, 69: {0, 2, 5, 7, 9, 12},
    70: {1, 3, 4, 8, 9, 12}, 71: {1, 3, 4, 10, 11, 12}, 72: {1, 3, 4, 7, 8, 12},
    73: {0, 4, 5, 7, 9, 12}, 74: {3, 4, 5, 7, 8, 12}, 75: {3, 5, 7, 8, 10, 12},
    76: {0, 1, 3, 10, 11, 12}, 77: {0, 1, 3, 10, 11, 16}, 78: {1, 4, 5, 6, 9, 12},
    79: {1, 3, 4, 8, 11, 12}, 80: {1, 3, 4, 6, 8, 12}, 81: {1, 2, 10, 11, 12, 13},
    82: {2, 5, 9, 10, 12, 13}, 83: {3, 5, 6, 9, 14, 18}, 84: {2, 5, 9, 11, 12, 13},
    85: {0, 5, 7, 9, 14, 18}, 86: {4, 5, 8, 9, 12, 13}, 87: {4, 5, 7, 9, 12, 14},
    88: {1, 5, 6, 10, 13, 14}, 89: {2, 5, 9, 10, 13, 14}, 90: {1, 5, 6, 9, 14, 18},
}
CONTACT_SETS = {k: frozenset(v) for k, v in CONTACT_SETS.items()}

AUX7 = {
    1: {0, 1, 2, 3, 6, 15, 17}, 2: {0, 1, 2, 3, 6, 15, 22}, 3: {0, 1, 2, 3, 6, 17, 19},
    4: {0, 1, 2, 3, 6, 17, 22}, 5: {0, 1, 2, 6, 9, 11, 20}, 6: {0, 1, 2, 6, 9, 11, 23},
    7: {0, 1, 2, 7, 9, 10, 14}, 8: {0, 1, 2, 7, 9, 10, 16}, 9: {0, 1, 2, 7, 9, 10, 17},
    10: {0, 1, 2, 7, 9, 10, 23}, 11: {0, 1, 2, 7, 9, 16, 17}, 12: {0, 1, 2, 7, 9, 16, 23},
    13: {0, 1, 2, 7, 10, 16, 17}, 14: {0, 1, 2, 7, 10, 16, 23},
    15: {0, 1, 2, 7, 14, 16, 23}, 16: {0, 1, 2, 9, 10, 14, 17},
    17: {0, 1, 2, 9, 10, 16, 17}, 18: {0, 1, 3, 10, 11, 12, 16},
    19: {0, 1, 4, 6, 9, 12, 17}, 20: {0, 1, 5, 7, 9, 12, 17},
    21: {0, 1, 6, 9, 16, 17, 22}, 22: {0, 2, 3, 7, 11, 12, 23},
    23: {0, 3, 5, 6, 9, 12, 17}, 24: {0, 3, 7, 11, 12, 16, 20},
    25: {1, 3, 4, 8, 11, 12, 15}, 26: {1, 3, 5, 6, 7, 8, 14}, 27: {1, 4, 5, 7, 9, 12, 14},
}
AUX7 = {k: frozenset(v) for k, v in AUX7.items()}

# Table 3 (p. 31): (j, source case nu, rho, rho', eps, listed image)
TABLE3 = [
    (1, 40, "1234", "1234", 0, {0, 1, 5, 9, 14, 18}),
    (4, 55, "1243", "1234", 0, {0, 4, 5, 7, 9, 18}),
    (6, 40, "1234", "2143", 0, {5, 6, 7, 9, 14, 18}),
    (11, 55, "2143", "1243", 1, {0, 5, 7, 9, 11, 14}),
    (12, 55, "1234", "2134", 1, {0, 5, 7, 12, 14, 18}),
    (15, 55, "2134", "2143", 0, {0, 7, 9, 14, 15, 18}),
    (17, 81, "1324", "2341", 1, {0, 5, 7, 9, 17, 18}),
]


def generators(T):
    """F_T = [j : f_j vanishes at every permutation index in T] (Section 4.2)."""
    T = set(T)
    return [j for j in range(82) if T <= F_ZERO[j]]


def enumerate_six_set_orbits():
    """All G-orbits of 6-subsets of the 24 permutation indices.  Returns a list of
    (lexicographically smallest representative tuple, orbit size)."""
    seen = set()
    orbits = []
    for T in itertools.combinations(range(24), 6):
        key = sum(1 << j for j in T)
        if key in seen:
            continue
        orb = set()
        for t in _G_MAPS:
            orb.add(sum(1 << t[j] for j in T))
        seen |= orb
        orbits.append((T, len(orb)))   # T is lexicographically first in its orbit
    return orbits


# =====================================================================================
# 12. NUMERIC helpers (numpy / scipy; floating point, best effort)
# =====================================================================================
# Everything below is floating-point.  Results are computational evidence only.

def _np():
    import numpy as np
    return np


def P_numpy():
    return _np().array(P, dtype=float)


_SIGMA_ARR = None


def leibniz_weights_np(O):
    """NUMERIC: w_sigma(O) = sgn(sigma) prod_i O[i, sigma(i)] for a float array O."""
    np = _np()
    global _SIGMA_ARR
    if _SIGMA_ARR is None:
        _SIGMA_ARR = np.array(PERMS)
    O = np.asarray(O)
    prod = O[0, _SIGMA_ARR[:, 0]] * O[1, _SIGMA_ARR[:, 1]] * O[2, _SIGMA_ARR[:, 2]] * \
        O[3, _SIGMA_ARR[:, 3]]
    return np.array(SGN) * prod


def Phi_np(c, O):
    """NUMERIC Leibniz sum."""
    np = _np()
    return float(P_numpy() @ np.asarray(c, dtype=float) @ leibniz_weights_np(O))


def skew6(theta):
    np = _np()
    a, b, c, d, e, f = theta
    return np.array([[0, a, b, c], [-a, 0, d, e], [-b, -d, 0, f], [-c, -e, -f, 0]], float)


def so4_from_params(theta):
    """NUMERIC parametrization of SO(4): O = expm(S(theta)), theta in R^6 (surjective)."""
    from scipy.linalg import expm
    return expm(skew6(theta))


def min_phi_over_SO4(c, restarts=40, seed=0):
    """NUMERIC: multistart BFGS minimum of Phi_c over SO(4).  Returns (value, O).
    The value is an upper bound for the true minimum (local search can miss the global
    minimum; it can never report a value below one attained)."""
    np = _np()
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    vals = P_numpy() @ np.asarray([float(x) for x in c])
    f = lambda th: float(vals @ leibniz_weights_np(so4_from_params(th)))
    best = (np.inf, None)
    for _ in range(restarts):
        th0 = rng.uniform(-np.pi, np.pi, 6)
        r = minimize(f, th0, method="BFGS")
        if r.fun < best[0]:
            best = (float(r.fun), so4_from_params(r.x))
    return best


def F_values_np(a, b, alpha=1.0):
    """NUMERIC: alpha * prod_i (a_i + b_sigma(i)) for complex arrays a, b."""
    np = _np()
    a, b = np.asarray(a, complex), np.asarray(b, complex)
    S = np.array(PERMS)
    return alpha * np.prod(a[None, :] + b[S], axis=1)


def matching_vector_np(a, b, alpha=1.0):
    """NUMERIC matching vector (complex array of length 14)."""
    z = [complex(x) for x in a] + [-complex(x) for x in b]
    return _np().array([alpha * eval_poly(B_POLYS[r], z) for r in range(14)], complex)


def random_point_on_contact_stratum(T, seed=0, tries=60, tol=1e-9, pi_star=None, margin=0.0):
    """NUMERIC, best effort: look for spectra a, b in C^4 and alpha = exp(i phi) such
    that all gaps Re(F_sigma - F_{pi*}) are >= 0 and vanish for sigma in T.  Solves the
    real equations by scipy least_squares from random starts, with the normalization
    mean(gaps off T) = 1.  With margin > 0 the gaps off T are pushed to be >= margin,
    so that the contact set is exactly T.  Returns dict(a, b, alpha, pi_star, gaps,
    contact) or None.  A returned point satisfies the conditions to ~tol in floating
    point only; with margin = 0, contacts outside T may also occur (the paper's contact
    conditions are non-strict).  Failure (None) says nothing about whether the stratum
    is empty."""
    np = _np()
    from scipy.optimize import least_squares
    T = sorted(T)
    ps = T[0] if pi_star is None else pi_star
    others = [s for s in range(24) if s not in T]
    rng = np.random.default_rng(seed)

    def unpack(x):
        a = x[0:4] + 1j * x[4:8]
        b = x[8:12] + 1j * x[12:16]
        return a, b, np.exp(1j * x[16])

    def res(x):
        a, b, al = unpack(x)
        F = F_values_np(a, b, al).real
        g = F - F[ps]
        r = [g[s] for s in T if s != ps]
        r += [min(0.0, g[s] - margin) for s in others]
        r.append(np.mean(g[others]) - 1.0)
        return np.array(r)

    for _ in range(tries):
        x0 = rng.normal(size=17)
        try:
            sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception:
            continue
        if np.max(np.abs(res(sol.x))) < tol:
            a, b, al = unpack(sol.x)
            F = F_values_np(a, b, al).real
            g = F - F[ps]
            return dict(a=a, b=b, alpha=al, pi_star=ps, gaps=g,
                        contact=sorted(int(s) for s in range(24) if abs(g[s]) < 1e-7))
    return None
