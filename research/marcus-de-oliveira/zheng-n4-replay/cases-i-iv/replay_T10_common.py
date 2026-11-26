"""Shared exact helpers for the independent replays of T_10, T_11 and T_15 (group III of Zheng's Table 1).

Everything here is rebuilt from the definitions (PROOF.md Sections 3.6, 6.1-6.6; paper Sections 4.1-4.2, Appendix D):
permutations, the 82 facet profiles as integer vectors of length 24 (their values at the 24 permutations), the group of
order 1152, Leibniz sums, and the homogeneous matching products of eight points of P^1.  Nothing from the supplement is
read, imported or executed.  The audit library (zheng-n4-audit/lib/moc4lib.py) is imported ONLY by `crosscheck_library`,
which compares its tables with the ones built here; no result depends on it.

Packages: Python standard library and sympy (exact); numpy only inside the functions whose name starts with `numeric_`.
"""
import itertools
import os
import sys
from fractions import Fraction

import sympy as sp

# ------------------------------------------------------------------------------------------------------------------
# bookkeeping
# ------------------------------------------------------------------------------------------------------------------


class Checker:
    def __init__(self, title):
        self.n = 0
        self.failed = []
        print(title)
        print("=" * len(title))

    def section(self, name):
        print()
        print("--", name)

    def check(self, name, ok, detail=""):
        self.n += 1
        ok = bool(ok)
        if not ok:
            self.failed.append(name)
        print("[%s] %s%s" % ("ok" if ok else "FAIL", name, (" :: " + str(detail)) if detail else ""))
        sys.stdout.flush()
        return ok

    def info(self, text):
        print("       " + str(text))
        sys.stdout.flush()

    def finish(self):
        print()
        if self.failed:
            print("FAIL: %d of %d checks failed: %s" % (len(self.failed), self.n, "; ".join(self.failed)))
            sys.exit(1)
        print("PASS: all %d checks passed" % self.n)
        sys.exit(0)


# ------------------------------------------------------------------------------------------------------------------
# permutations (0-based tuples, lexicographic order of one-line notation; index j is the paper's pi_j)
# ------------------------------------------------------------------------------------------------------------------
PERMS = list(itertools.permutations(range(4)))
PIDX = {p: j for j, p in enumerate(PERMS)}


def sgn(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


SGN = [sgn(p) for p in PERMS]


def compose(p, q):
    """(p o q)(i) = p(q(i))."""
    return tuple(p[q[i]] for i in range(4))


def inverse(p):
    r = [0] * 4
    for i in range(4):
        r[p[i]] = i
    return tuple(r)


# ------------------------------------------------------------------------------------------------------------------
# the 82 facet profiles (values at the 24 permutations), from the definitions of Section 4.1 / PROOF.md 6.1, 6.3
# ------------------------------------------------------------------------------------------------------------------
# Appendix D of the paper (p. 59), own transcription from the text extract: ("n", I, J) or ("h", one-line pi, eps).
APPENDIX_D = [
    ("n", "1", "4"), ("n", "12", "34"), ("n", "13", "34"), ("n", "14", "34"), ("n", "1", "3"), ("n", "12", "24"),
    ("n", "4", "1"), ("n", "3", "1"), ("n", "12", "23"), ("n", "14", "23"), ("n", "13", "24"), ("n", "13", "23"),
    ("n", "14", "24"), ("n", "2", "1"), ("n", "1", "2"), ("n", "12", "14"), ("n", "2", "4"), ("h", "3421", 1),
    ("h", "3412", -1), ("h", "2413", -1), ("h", "2431", 1), ("n", "4", "2"), ("h", "3412", 1), ("n", "14", "12"),
    ("h", "4312", 1), ("h", "3142", -1), ("h", "4132", -1), ("n", "3", "2"), ("h", "3421", -1), ("n", "13", "12"),
    ("h", "4321", -1), ("h", "4123", 1), ("h", "3124", 1), ("n", "12", "13"), ("h", "4321", 1), ("h", "4312", -1),
    ("n", "2", "3"), ("h", "2341", 1), ("h", "2314", -1), ("n", "14", "13"), ("h", "4213", 1),
    ("n", "4", "3"), ("h", "4123", -1), ("h", "2413", 1), ("h", "2143", -1), ("n", "13", "14"), ("h", "3241", -1),
    ("n", "3", "4"), ("h", "3142", 1), ("h", "2143", 1), ("h", "2341", -1), ("h", "1423", -1), ("h", "1342", 1),
    ("n", "13", "13"), ("h", "4231", -1), ("n", "3", "3"), ("h", "4132", 1), ("h", "2431", -1), ("h", "2134", 1),
    ("n", "14", "14"), ("h", "3214", 1), ("n", "4", "4"), ("h", "3124", -1), ("h", "2134", -1), ("h", "2314", 1),
    ("h", "1432", 1), ("h", "1324", -1), ("n", "2", "2"), ("h", "4213", -1), ("h", "3241", 1), ("h", "4231", 1),
    ("h", "3214", -1), ("n", "12", "12"), ("h", "1243", 1), ("h", "1432", -1), ("h", "1423", 1), ("h", "1234", -1),
    ("h", "1342", -1), ("h", "1234", 1), ("h", "1243", -1), ("h", "1324", 1), ("n", "1", "1"),
]
assert len(APPENDIX_D) == 82


def n_profile(I, J):
    """n_{I,J}: 1 at sigma if sigma(I) = J, else 0.  I, J are sets of 1-based labels."""
    I0 = [i - 1 for i in I]
    J0 = set(j - 1 for j in J)
    return tuple(1 if set(s[i] for i in I0) == J0 else 0 for s in PERMS)


def t_profile(pi):
    """t_pi = sum_i n_{i, pi(i)}: the number of agreements with pi (pi is a 0-based tuple)."""
    return tuple(sum(1 for i in range(4) if s[i] == pi[i]) for s in PERMS)


def d_profile(pi):
    """d_pi = n_{12,{pi1,pi3}} + n_{12,{pi2,pi4}} - n_{13,{pi1,pi2}} - n_{13,{pi3,pi4}} (1-based labels)."""
    q = [x + 1 for x in pi]
    a = n_profile((1, 2), (q[0], q[2]))
    b = n_profile((1, 2), (q[1], q[3]))
    c = n_profile((1, 3), (q[0], q[1]))
    e = n_profile((1, 3), (q[2], q[3]))
    return tuple(a[k] + b[k] - c[k] - e[k] for k in range(24))


def h_profile(pi, eps):
    t = t_profile(pi)
    d = d_profile(pi)
    out = []
    for k in range(24):
        v = t[k] + eps * d[k]
        assert v % 2 == 0
        out.append(v // 2)
    return tuple(out)


def _profile(spec):
    kind, a, b = spec
    if kind == "n":
        return n_profile([int(ch) for ch in a], [int(ch) for ch in b])
    return h_profile(tuple(int(ch) - 1 for ch in a), b)


F = [_profile(s) for s in APPENDIX_D]            # F[j][sigma_index] = value of f_j at pi_sigma_index
IS_MINOR = [s[0] == "n" for s in APPENDIX_D]
FIDX = {f: j for j, f in enumerate(F)}


def fname(j):
    k, a, b = APPENDIX_D[j]
    return "n_{%s,%s}" % (a, b) if k == "n" else "h_{%s,%s}" % (a, "+" if b > 0 else "-")


def generators(T):
    """Indices j with f_j = 0 on all of T (the set F_T of the paper's Section 4.2)."""
    return [j for j in range(82) if all(F[j][s] == 0 for s in T)]


def minor_sets(j):
    """For a minor vector f_j = n_{I,J}: the 0-based index lists (I, J)."""
    k, a, b = APPENDIX_D[j]
    assert k == "n"
    return [int(ch) - 1 for ch in a], [int(ch) - 1 for ch in b]


def lin_comb(coeffs):
    """sum_j coeffs[j] * f_j as a 24-vector (coeffs: dict j -> number)."""
    return tuple(sum(c * F[j][s] for j, c in coeffs.items()) for s in range(24))


def functional(expr):
    """A linear functional on gap vectors given as dict sigma_index -> coefficient; returns a function on profiles."""
    return lambda prof: sum(c * prof[s] for s, c in expr.items())


# ------------------------------------------------------------------------------------------------------------------
# the group G = (S4 x S4) x| C2 of order 1152 (PROOF.md 6.6): sigma -> rho o sigma^{+-1} o rho'
# ------------------------------------------------------------------------------------------------------------------
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_perm(g, s):
    r1, r2, e = g
    return compose(r1, compose(inverse(s) if e else s, r2))


def g_map(g):
    """tuple m with m[j] = index of g.pi_j."""
    return tuple(PIDX[g_perm(g, s)] for s in PERMS)


G_MAPS = [g_map(g) for g in G_ELEMENTS]


def g_profile(m, prof):
    """(g.f)(g.sigma) = f(sigma), with m = g_map(g)."""
    out = [None] * 24
    for j in range(24):
        out[m[j]] = prof[j]
    return tuple(out)


def g_set(m, T):
    return frozenset(m[j] for j in T)


def g_facet(m, j):
    """index of the facet vector g.f_j (None if it is not one of the 82)."""
    return FIDX.get(g_profile(m, F[j]))


def find_set_maps(T, T2):
    T2 = frozenset(T2)
    return [k for k, m in enumerate(G_MAPS) if g_set(m, T) == T2]


def find_tuple_maps(src, dst, unordered_prefix=0):
    """All group elements k with g.f_{src[i]} = f_{dst[i]} for all i; the first `unordered_prefix` entries of the
    tuples are compared as sets."""
    out = []
    u = unordered_prefix
    for k, m in enumerate(G_MAPS):
        img = [g_facet(m, j) for j in src]
        if None in img:
            continue
        if set(img[:u]) == set(dst[:u]) and len(set(img[:u])) == len(set(dst[:u])) and img[u:] == list(dst[u:]):
            out.append(k)
    return out


def g_matrix(k, O):
    """The matrix O' with w_{g.sigma}(O') = sgn(rho rho') w_sigma(O) (derived in T10.md, 'Covariance'):
    O'[a][b] = O~[rho'(a)][rho^{-1}(b)], where O~ = O (eps = 0) or O^t (eps = 1)."""
    r1, r2, e = G_ELEMENTS[k]
    if e:
        O = [[O[j][i] for j in range(4)] for i in range(4)]
    r1i = inverse(r1)
    return [[O[r2[a]][r1i[b]] for b in range(4)] for a in range(4)]


def g_sign(k):
    r1, r2, e = G_ELEMENTS[k]
    return sgn(r1) * sgn(r2)


# ------------------------------------------------------------------------------------------------------------------
# Leibniz sums
# ------------------------------------------------------------------------------------------------------------------


def leibniz_weights(O):
    return [SGN[j] * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]] for j, s in enumerate(PERMS)]


def leibniz(prof, O):
    w = leibniz_weights(O)
    return sum(prof[j] * w[j] for j in range(24))


def minor_value(j, O):
    """The entry or 2x2 minor x_j of O with Phi_{f_j} = x_j^2 on SO(4) (PROOF.md Lemma 6.2)."""
    I, J = minor_sets(j)
    if len(I) == 1:
        return O[I[0]][J[0]]
    return O[I[0]][J[0]] * O[I[1]][J[1]] - O[I[0]][J[1]] * O[I[1]][J[0]]


def cayley_so4(skew6):
    """Exact rational SO(4) matrix (I - S)(I + S)^{-1} from the six entries of a skew matrix S."""
    a, b, c, d, e, f = [sp.Rational(x) for x in skew6]
    S = sp.Matrix([[0, a, b, c], [-a, 0, d, e], [-b, -d, 0, f], [-c, -e, -f, 0]])
    O = (sp.eye(4) - S) * (sp.eye(4) + S).inv()
    return [[O[i, j] for j in range(4)] for i in range(4)]


def is_so4(O):
    M = sp.Matrix(O)
    return (M * M.T - sp.eye(4)).is_zero_matrix and M.det() == 1


def rational_so4_samples(n, seed):
    import random
    rng = random.Random(seed)
    out = []
    while len(out) < n:
        sk = [Fraction(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(6)]
        O = cayley_so4(sk)
        # compose with a signed permutation so that the samples are not all near the identity chart
        p = rng.choice(PERMS)
        Pm = [[1 if p[i] == j else 0 for j in range(4)] for i in range(4)]
        if sgn(p) < 0:
            Pm[0] = [-x for x in Pm[0]]
        M = (sp.Matrix(O) * sp.Matrix(Pm))
        out.append([[M[i, j] for j in range(4)] for i in range(4)])
    return out


# ------------------------------------------------------------------------------------------------------------------
# covariance of the Leibniz sums under one group element, as a polynomial identity in 16 free entries
# ------------------------------------------------------------------------------------------------------------------
_OSYM = [[sp.Symbol("o%d%d" % (i, j)) for j in range(4)] for i in range(4)]


def covariance_identity(k, src):
    """For every j in src: Phi_{g.f_j}(O') = sgn(rho rho') Phi_{f_j}(O) identically in the 16 entries of O."""
    m = G_MAPS[k]
    Op = g_matrix(k, _OSYM)
    s = g_sign(k)
    for j in src:
        lhs = leibniz(g_profile(m, F[j]), Op)
        rhs = s * leibniz(F[j], _OSYM)
        if sp.expand(lhs - rhs) != 0:
            return False
    return True


def minor_covariance(k, src):
    """For minor vectors f_j, j in src: the entry/minor x_{g.f_j}(O') equals +- x_{f_j}(O) identically."""
    m = G_MAPS[k]
    Op = g_matrix(k, _OSYM)
    for j in src:
        jj = g_facet(m, j)
        if jj is None or not IS_MINOR[jj] or not IS_MINOR[j]:
            return False
        a = sp.expand(minor_value(jj, Op))
        b = sp.expand(minor_value(j, _OSYM))
        if sp.expand(a - b) != 0 and sp.expand(a + b) != 0:
            # complementary minor: n_{I,J} = n_{I^c,J^c}; on SO(4) the two minors agree up to sign (Jacobi), but
            # they are different polynomials.  Accept the complementary one explicitly.
            I, J = minor_sets(j)
            if len(I) == 2:
                Ic = [i for i in range(4) if i not in I]
                Jc = [i for i in range(4) if i not in J]
                c = sp.expand(_OSYM[Ic[0]][Jc[0]] * _OSYM[Ic[1]][Jc[1]] - _OSYM[Ic[0]][Jc[1]] * _OSYM[Ic[1]][Jc[0]])
                if sp.expand(a - c) == 0 or sp.expand(a + c) == 0:
                    continue
            return False
    return True


# ------------------------------------------------------------------------------------------------------------------
# homogeneous matching products of eight points of P^1 (PROOF.md 3.3 with v_i = (A_i,1), v_{4+j} = (C_j,1), C = -b)
# ------------------------------------------------------------------------------------------------------------------
INF = sp.Symbol("INFINITY_MARKER")


def bracket_AC(Ai, Cj):
    """[v_i, v_{4+j}] = det((A_i,1),(C_j,1)) = A_i - C_j; a point at infinity is the vector (1,0)."""
    if Ai is INF and Cj is INF:
        return sp.Integer(0)
    if Cj is INF:
        return sp.Integer(-1)          # det((A,1),(1,0)) = -1
    if Ai is INF:
        return sp.Integer(1)           # det((1,0),(C,1)) = +1
    return Ai - Cj


def products(A, C):
    """F_sigma = prod_i [v_i, v_{4+sigma(i)}] for the 24 permutations, lexicographic order."""
    return [sp.Mul(*[bracket_AC(A[i], C[s[i]]) for i in range(4)]) for s in PERMS]


# ------------------------------------------------------------------------------------------------------------------
# symbolic helpers.  A nonreal parameter v is represented by two independent symbols (v, vb); every other parameter is
# real.  For a rational function R with real coefficients, conj R is R with v <-> vb, and
#     Im R / Im v = (R - conj R) / (v - vb),
# a rational function with rational coefficients (no imaginary unit needed).
# ------------------------------------------------------------------------------------------------------------------


def swap(expr, v, vb):
    return expr.xreplace({v: vb, vb: v})


def im_ratio(expr, v, vb):
    """Im(expr)/Im(v) as an element of Q(v, vb, real parameters)."""
    return sp.cancel((expr - swap(expr, v, vb)) / (v - vb))


def re_part2(expr, v, vb):
    """2 Re(expr)."""
    return expr + swap(expr, v, vb)


def is_zero(expr):
    """Exact test that a rational function is identically zero."""
    e = sp.together(expr)
    num, den = sp.fraction(e)
    return sp.expand(num) == 0


def same(a, b):
    return is_zero(a - b)


def numer_denom(expr):
    e = sp.cancel(sp.together(expr))
    return sp.fraction(e)


def free_of(expr, syms):
    e = sp.cancel(sp.together(expr))
    return not (e.free_symbols & set(syms))


def poly_coeffs_nonneg(expr, gens):
    """True if expr is a polynomial in gens all of whose coefficients are nonnegative rationals."""
    p = sp.Poly(sp.expand(expr), *gens)
    return all(c >= 0 for c in p.coeffs())


# ------------------------------------------------------------------------------------------------------------------
# cross-check against the audit library (optional; nothing above depends on it)
# ------------------------------------------------------------------------------------------------------------------


def crosscheck_library(ck, cases):
    here = os.path.dirname(os.path.abspath(__file__))
    lib = os.path.join(here, "..", "..", "zheng-n4-audit", "lib")
    if not os.path.isfile(os.path.join(lib, "moc4lib.py")):
        ck.info("audit library not found; cross-check skipped")
        return
    sys.path.insert(0, lib)
    import moc4lib as L
    ck.check("library: PERMS identical to mine", list(L.PERMS) == PERMS)
    ck.check("library: Appendix D table identical to my transcription",
             [tuple(x) for x in L.F_TABLE] == [tuple(x) for x in APPENDIX_D])
    ck.check("library: all 82 profiles P f_j equal my profiles built from the definitions",
             [tuple(v) for v in L.F_VALUES] == F)
    ck.check("library: 1152 index maps of G equal mine", [tuple(m) for m in L._G_MAPS] == G_MAPS)
    for nu, T in cases.items():
        ck.check("library: Appendix E contact set T_%d = %s" % (nu, sorted(T)), L.CONTACT_SETS[nu] == frozenset(T))
        ck.check("library: generators(T_%d) equal mine" % nu, sorted(L.generators(T)) == generators(T))


def basic_selfchecks(ck):
    ck.check("24 permutations, lexicographic; identity is pi_0", len(PERMS) == 24 and PERMS[0] == (0, 1, 2, 3))
    ck.check("82 distinct facet profiles; 34 minor and 48 non-minor",
             len(set(F)) == 82 and sum(IS_MINOR) == 34)
    prof_ok = all(sorted(set(f)) in ([0, 1],) for j, f in enumerate(F) if IS_MINOR[j]) and \
        all((f.count(0), f.count(1), f.count(2)) == (13, 10, 1) for j, f in enumerate(F) if not IS_MINOR[j])
    ck.check("value profiles: minors 0/1; non-minors take 0 thirteen times, 1 ten times, 2 once", prof_ok)
    ck.check("G has 1152 elements acting by 1152 distinct permutations of the 24 indices",
             len(G_MAPS) == 1152 and len(set(G_MAPS)) == 1152)
    ck.check("G permutes the 82 profiles and preserves minor/non-minor",
             all(g_facet(m, j) is not None and IS_MINOR[g_facet(m, j)] == IS_MINOR[j]
                 for m in G_MAPS[::37] for j in range(82)))
    # Lemma 6.2 on exact rational points: Phi_{n_{I,J}}(O) = det(O[I,J])^2
    Os = rational_so4_samples(4, 11)
    ck.check("4 exact rational matrices are in SO(4)", all(is_so4(O) for O in Os))
    ck.check("Phi_{f_j}(O) = x_j(O)^2 for all 34 minor vectors at these matrices (Lemma 6.2)",
             all(leibniz(F[j], O) == minor_value(j, O) ** 2 for O in Os for j in range(82) if IS_MINOR[j]))


# ------------------------------------------------------------------------------------------------------------------
# the universal inequalities of PROOF.md 7.6 that these three cases invoke: representatives and transfer
# ------------------------------------------------------------------------------------------------------------------
SOC_REP = (1, 9, 22)            # Phi_22 + 2|x_1 x_9| >= 0            (so4-atoms.md 3.1), pair unordered
FOUR_REP = (0, 6, 10, 54)       # Phi_54 + Phi_10 + 2(|x_0 x_6| + |x_0 x_10| + |x_6 x_10|) >= 0  (so4-atoms.md 3.4)


def orbit_size(rep, unordered_prefix):
    imgs = set()
    u = unordered_prefix
    for m in G_MAPS:
        img = [g_facet(m, j) for j in rep]
        imgs.add((frozenset(img[:u]),) + tuple(img[u:]))
    return len(imgs)


def check_atom_use(ck, label, rep, used, unordered_prefix):
    """The tuple `used` is g.rep for some g in G (roles respected), and the Leibniz sums transfer along g."""
    ks = find_tuple_maps(rep, used, unordered_prefix)
    ck.check("%s: used tuple %s is in the G-orbit of the representative %s with the same roles (%d group elements)"
             % (label, used, rep, len(ks)), len(ks) > 0)
    if not ks:
        return None
    k = ks[0]
    ck.check("%s: Phi_{g.f}(O') = sgn(g) Phi_f(O) as polynomial identities in the 16 entries of O, for the 3-4 "
             "vectors of the tuple (g = %s)" % (label, (G_ELEMENTS[k],)), covariance_identity(k, rep))
    ck.check("%s: the minors of the used tuple are, up to sign and complementation, the minors of the representative "
             "evaluated at O'" % label, minor_covariance(k, rep[:-1]))
    return k


def numeric_atom_attack(kind, tup, n_random, n_starts, seed):
    """Floating-point falsification attempt of the pointwise form at the orbit member `tup`.
    SOC (i, j, k):        Phi_k + 2 sqrt(Phi_i Phi_j) >= 0
    four (a, b, c, d):    Phi_d + Phi_c + 2(sqrt(Phi_a Phi_b) + sqrt(Phi_a Phi_c) + sqrt(Phi_b Phi_c)) >= 0
    Returns the smallest value found (random Haar-like sampling plus Nelder-Mead from the worst points)."""
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    sig = np.array(PERMS)
    sg = np.array(SGN, dtype=float)
    Fv = np.array(F, dtype=float)

    def so4(pq):
        p = pq[:4] / np.linalg.norm(pq[:4])
        q = pq[4:] / np.linalg.norm(pq[4:])
        a, b, c, d = p
        Lm = np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]])
        a, b, c, d = q
        Rm = np.array([[a, b, c, d], [-b, a, -d, c], [-c, d, a, -b], [-d, -c, b, a]])
        return Lm @ Rm

    def phis(O):
        w = sg * O[0, sig[:, 0]] * O[1, sig[:, 1]] * O[2, sig[:, 2]] * O[3, sig[:, 3]]
        return Fv @ w

    def value(pq):
        ph = phis(so4(pq))
        pos = lambda j: max(ph[j], 0.0)
        if kind == "SOC":
            i, j, k = tup
            return ph[k] + 2 * np.sqrt(pos(i) * pos(j))
        a, b, c, d = tup
        return ph[d] + ph[c] + 2 * (np.sqrt(pos(a) * pos(b)) + np.sqrt(pos(a) * pos(c)) + np.sqrt(pos(b) * pos(c)))

    pts = rng.normal(size=(n_random, 8))
    vals = np.array([value(p) for p in pts])
    best = float(vals.min())
    for idx in np.argsort(vals)[:n_starts]:
        res = minimize(value, pts[idx], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 4000})
        best = min(best, float(res.fun))
    return best


# ------------------------------------------------------------------------------------------------------------------
# exact rational-function arithmetic: sympy fraction fields (sympy.polys.fields.field) over ZZ or QQ.
# Elements are kept as reduced fractions of sparse polynomials; `fzero` tests that a numerator is the zero polynomial.
# ------------------------------------------------------------------------------------------------------------------


def fzero(e):
    """Exact zero test for a fraction-field element (or a number)."""
    if hasattr(e, "numer"):
        return e.numer == 0
    return e == 0


def feq(a, b):
    return fzero(a - b)


def field_products(one, A, Cc):
    """Homogeneous matching products F_sigma for field elements; None marks a point at infinity:
    bracket((A,1),(1,0)) = -1 for C = infinity, bracket((1,0),(C,1)) = +1 for A = infinity."""
    out = []
    for s in PERMS:
        val = one
        for i in range(4):
            a = A[i]
            c = Cc[s[i]]
            if a is None and c is None:
                val = val * 0
            elif c is None:
                val = -val
            elif a is None:
                pass
            else:
                val = val * (a - c)
        out.append(val)
    return out


def fexpr(e):
    """sympy expression of a fraction-field element."""
    return e.numer.as_expr() / e.denom.as_expr()


def denominator_factors(elems):
    """Set of irreducible factors (as strings, sign-normalised by sympy) of the denominators of field elements."""
    out = set()
    for e in elems:
        for f_, _m in sp.factor_list(e.denom.as_expr())[1]:
            out.add(sp.sstr(sp.expand(f_)))
    return out


def factor_strings(exprs):
    """Sign-normalised strings for a list of sympy polynomials, for comparison with `denominator_factors`."""
    out = set()
    for e in exprs:
        for f_, _m in sp.factor_list(sp.expand(e))[1]:
            out.add(sp.sstr(sp.expand(f_)))
    return out


# ------------------------------------------------------------------------------------------------------------------
# floating point (NUMERIC only): vectorised gap vectors of a chart
# ------------------------------------------------------------------------------------------------------------------


def numeric_gaps(A, Cc, kappa, ref=0):
    """A, Cc: lists of 4 numpy arrays (or None for a point at infinity); returns the array (n, 24) of
    Re(kappa (F_sigma - F_ref)), F_sigma the homogeneous matching products."""
    import numpy as np
    cols = []
    for s in PERMS:
        val = 1.0
        for i in range(4):
            a = A[i]
            c = Cc[s[i]]
            if c is None:
                val = -val
            elif a is None:
                pass
            else:
                val = val * (a - c)
        cols.append(val)
    Fm = np.stack([np.broadcast_to(c_, np.shape(kappa)).astype(complex) for c_ in cols], axis=1)
    return (kappa[:, None] * (Fm - Fm[:, [ref]])).real
