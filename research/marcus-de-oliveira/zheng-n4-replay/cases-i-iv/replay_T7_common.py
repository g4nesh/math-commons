"""Shared helper for the independent replay of the group-II cases T_7, T_18, T_19, T_20.

Everything here is re-derived from the definitions in Zheng's paper (Sections 3.1, 4.1, 4.2 and
Appendices D, E).  Nothing from the supplement is executed or imported; the only access to it is
the optional sha256 of a proof PDF in document_hash_check (when MOC4_SUPP is set).  The audit library
`zheng-n4-audit/lib/moc4lib.py` is NOT used to build any object; it is imported (if present) only
by `crosscheck_with_audit_lib()` to compare the 82 value vectors and the four contact sets.

Conventions
-----------
* Permutations of {0,1,2,3}, lexicographic order of one-line notation, index 0..23 (pi_0 = id).
* A vector c of V_R is represented by its 24 values (P c)_sigma.  P is injective, so nothing is
  lost; the minor vectors, the vectors h_{pi,eps} and every gap vector are given in this form.
* Labels: A_i = a_i, C_j = -b_j (points of the projective line).  F_sigma = prod_i [A_i, C_sigma(i)].
* Gaps: g_sigma = Re(kappa (F_sigma - F_0)).

Packages: Python standard library and sympy (polynomial rings / fraction fields over QQ only) for
everything exact.  numpy and scipy are imported lazily, inside the helpers marked NUMERIC.
"""
import itertools
import os
import sys
from fractions import Fraction

from sympy import QQ
from sympy.polys.fields import field as _sym_field

# --------------------------------------------------------------------------------------------
# check harness
# --------------------------------------------------------------------------------------------


class Checks:
    def __init__(self, title):
        self.title = title
        self.n = 0
        self.failed = []
        print("=" * 100)
        print(title)
        print("=" * 100)

    def section(self, text):
        print()
        print("--- " + text)

    def check(self, name, ok, detail=""):
        self.n += 1
        ok = bool(ok)
        tag = "ok  " if ok else "FAIL"
        line = "[%s] %s" % (tag, name)
        if detail:
            line += "   | " + str(detail)
        print(line)
        if not ok:
            self.failed.append(name)
        return ok

    def info(self, text):
        print("       " + str(text))

    def finish(self):
        print()
        if self.failed:
            print("FAIL: %d of %d checks failed: %s" % (len(self.failed), self.n, self.failed))
            sys.exit(1)
        print("PASS: all %d checks passed (%s)" % (self.n, self.title))
        sys.exit(0)


# --------------------------------------------------------------------------------------------
# permutations and the 82 facet vectors (values on the 24 permutations)
# --------------------------------------------------------------------------------------------

PERMS = [tuple(p) for p in itertools.permutations(range(4))]      # lexicographic
PIDX = {p: i for i, p in enumerate(PERMS)}


def n_values(I, J):
    """Values of the minor vector n_{I,J}: 1 if sigma(I) = J, else 0.  I, J are 1-based."""
    I0 = [i - 1 for i in I]
    J0 = set(j - 1 for j in J)
    return tuple(1 if set(s[i] for i in I0) == J0 else 0 for s in PERMS)


def t_values(pi):
    """t_pi = sum_i n_{i, pi(i)}; pi is a 0-based tuple."""
    return tuple(sum(1 for i in range(4) if s[i] == pi[i]) for s in PERMS)


def d_values(pi):
    """d_pi = n_{12,{pi1,pi3}} + n_{12,{pi2,pi4}} - n_{13,{pi1,pi2}} - n_{13,{pi3,pi4}}."""
    p = [x + 1 for x in pi]
    a = n_values((1, 2), (p[0], p[2]))
    b = n_values((1, 2), (p[1], p[3]))
    c = n_values((1, 3), (p[0], p[1]))
    d = n_values((1, 3), (p[2], p[3]))
    return tuple(a[i] + b[i] - c[i] - d[i] for i in range(24))


def h_values(pi, eps):
    t = t_values(pi)
    d = d_values(pi)
    out = []
    for i in range(24):
        v = t[i] + eps * d[i]
        assert v % 2 == 0
        out.append(v // 2)
    return tuple(out)


# Appendix D of the paper (p. 59), own transcription from the text of the PDF.
# "n I J": minor vector n_{I,J} (digits are 1-based indices); "h pi s": h_{pi, s} (pi in one-line
# notation, 1-based; s = + or -).
_APPENDIX_D = """
0 n 1 4 | 1 n 12 34 | 2 n 13 34 | 3 n 14 34 | 4 n 1 3 | 5 n 12 24 | 6 n 4 1 | 7 n 3 1 |
8 n 12 23 | 9 n 14 23 | 10 n 13 24 | 11 n 13 23 | 12 n 14 24 | 13 n 2 1 | 14 n 1 2 |
15 n 12 14 | 16 n 2 4 | 17 h 3421 + | 18 h 3412 - | 19 h 2413 - | 20 h 2431 + | 21 n 4 2 |
22 h 3412 + | 23 n 14 12 | 24 h 4312 + | 25 h 3142 - | 26 h 4132 - | 27 n 3 2 | 28 h 3421 - |
29 n 13 12 | 30 h 4321 - | 31 h 4123 + | 32 h 3124 + | 33 n 12 13 | 34 h 4321 + | 35 h 4312 - |
36 n 2 3 | 37 h 2341 + | 38 h 2314 - | 39 n 14 13 | 40 h 4213 + | 41 n 4 3 | 42 h 4123 - |
43 h 2413 + | 44 h 2143 - | 45 n 13 14 | 46 h 3241 - | 47 n 3 4 | 48 h 3142 + | 49 h 2143 + |
50 h 2341 - | 51 h 1423 - | 52 h 1342 + | 53 n 13 13 | 54 h 4231 - | 55 n 3 3 | 56 h 4132 + |
57 h 2431 - | 58 h 2134 + | 59 n 14 14 | 60 h 3214 + | 61 n 4 4 | 62 h 3124 - | 63 h 2134 - |
64 h 2314 + | 65 h 1432 + | 66 h 1324 - | 67 n 2 2 | 68 h 4213 - | 69 h 3241 + | 70 h 4231 + |
71 h 3214 - | 72 n 12 12 | 73 h 1243 + | 74 h 1432 - | 75 h 1423 + | 76 h 1234 - | 77 h 1342 - |
78 h 1234 + | 79 h 1243 - | 80 h 1324 + | 81 n 1 1
"""


def _parse_appendix_d():
    specs = {}
    for item in _APPENDIX_D.replace("\n", " ").split("|"):
        tok = item.split()
        if not tok:
            continue
        j = int(tok[0])
        if tok[1] == "n":
            specs[j] = ("n", tuple(int(ch) for ch in tok[2]), tuple(int(ch) for ch in tok[3]))
        else:
            specs[j] = ("h", tuple(int(ch) - 1 for ch in tok[2]), 1 if tok[3] == "+" else -1)
    assert sorted(specs) == list(range(82))
    return [specs[j] for j in range(82)]


F_SPEC = _parse_appendix_d()


def _f_values(spec):
    if spec[0] == "n":
        return n_values(spec[1], spec[2])
    return h_values(spec[1], spec[2])


F_VALUES = [_f_values(s) for s in F_SPEC]
F_IS_MINOR = [s[0] == "n" for s in F_SPEC]

# Appendix E of the paper, the four sets of this replay (own transcription).
PAPER_T = {7: (0, 1, 2, 4, 7, 8), 18: (0, 1, 3, 6, 7, 8), 19: (0, 1, 4, 6, 7, 8), 20: (0, 1, 5, 6, 7, 8)}


def f_name(j):
    s = F_SPEC[j]
    if s[0] == "n":
        return "n_{%s,%s}" % ("".join(map(str, s[1])), "".join(map(str, s[2])))
    return "h_{%s,%s}" % ("".join(str(x + 1) for x in s[1]), "+" if s[2] > 0 else "-")


def generators(T):
    """Indices j with f_j = 0 on T (the set F_T of the paper)."""
    return [j for j in range(82) if all(F_VALUES[j][s] == 0 for s in T)]


DOC_SHA256 = {
    "root_moc4_T6_uniform_20261022.pdf": "aa11d6eba3749cc374e62c7d32a8b343ae22f570cbc3113b954978316e1618d7",
    "root_moc4_T17_uniform_20261022.pdf": "95a7f8f4c5ce424844f6fa506e7e95529bb08e6bb9fc76afb56e22cd5fa9671b",
    "root_moc4_T18_uniform_20261022.pdf": "a7e305a493f06c5129eb693a6bf71277a8505862e5d08191c0437b86956979d4",
    "root_moc4_T19_uniform_20261022.pdf": "d045da12e70c3ef1904e15866bc21a6b932b1d986f3b7cee35823e05cb6a7a7c",
}


def document_hash_check(ck, name):
    """Pin the document that was replayed.  The PDF is hashed as bytes, never parsed or executed.
    Skipped (not a failure) when MOC4_SUPP is not set."""
    import hashlib
    root = os.environ.get("MOC4_SUPP")
    if not root:
        ck.info("MOC4_SUPP not set: sha256 comparison of %s skipped" % name)
        return
    path = os.path.join(root, "proofs", "manin_transfer_20261021", name)
    with open(path, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    ck.check("sha256 of %s is the one of the document that was read" % name,
             digest == DOC_SHA256[name], digest)


def facet_table_sanity(ck):
    """Structural facts of Section 4.1 that pin down the 82 vectors (EXACT, complete)."""
    minors = [F_VALUES[j] for j in range(82) if F_IS_MINOR[j]]
    hs = [F_VALUES[j] for j in range(82) if not F_IS_MINOR[j]]
    ck.check("Appendix D: 34 minor vectors and 48 vectors h, all 82 value vectors distinct",
             len(minors) == 34 and len(hs) == 48 and len(set(F_VALUES)) == 82)
    all_minor = set()
    for i in range(1, 5):
        for j in range(1, 5):
            all_minor.add(n_values((i,), (j,)))
    for I in itertools.combinations(range(1, 5), 2):
        for J in itertools.combinations(range(1, 5), 2):
            all_minor.add(n_values(I, J))
    ck.check("the 34 listed minor vectors are exactly the 16 + 18 distinct n_{I,J}",
             len(all_minor) == 34 and set(minors) == all_minor)
    all_h = set(h_values(pi, e) for pi in PERMS for e in (1, -1))
    ck.check("the 48 listed h are exactly all h_{pi,eps}; profile thirteen 0, ten 1, one 2",
             len(all_h) == 48 and set(hs) == all_h
             and all(sorted(v) == [0] * 13 + [1] * 10 + [2] for v in hs))
    ck.check("all 82 vectors are >= 0 on the 24 permutations",
             all(min(v) >= 0 for v in F_VALUES))


def crosscheck_with_audit_lib(ck):
    """Compare own objects with the audit library (cross-check only; not used otherwise)."""
    here = os.path.dirname(os.path.abspath(__file__))
    libdir = os.path.normpath(os.path.join(here, "..", "..", "zheng-n4-audit", "lib"))
    if not os.path.isfile(os.path.join(libdir, "moc4lib.py")):
        ck.info("audit library not found at %s: cross-check skipped" % libdir)
        return
    sys.path.insert(0, libdir)
    try:
        import moc4lib as L
    finally:
        sys.path.pop(0)
    ck.check("cross-check: own 82 value vectors = moc4lib.F_VALUES (same numbering)",
             [tuple(v) for v in L.F_VALUES] == [tuple(v) for v in F_VALUES])
    ck.check("cross-check: own PERMS = moc4lib.PERMS; minor flags agree",
             list(L.PERMS) == PERMS and list(L.F_IS_MINOR) == F_IS_MINOR)
    ck.check("cross-check: T_7, T_18, T_19, T_20 of Appendix E agree with moc4lib.CONTACT_SETS",
             all(set(L.CONTACT_SETS[nu]) == set(PAPER_T[nu]) for nu in PAPER_T))
    # value vectors really are images under P of vectors in V_R (uses the library's P)
    ok = True
    for j in range(82):
        try:
            L.solve_P(list(F_VALUES[j]))
        except Exception:
            ok = False
    ck.check("cross-check: each own value vector lies in the image of moc4lib.P", ok)


# --------------------------------------------------------------------------------------------
# exact linear algebra over Q
# --------------------------------------------------------------------------------------------


def _rref(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    piv = []
    r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return M[:r], piv


def rank(rows):
    return len(_rref(rows)[1]) if rows else 0


def nullspace(rows, ncol):
    """Basis of {x : rows . x = 0} over Q."""
    if not rows:
        return [[Fraction(int(i == j)) for i in range(ncol)] for j in range(ncol)]
    R, piv = _rref(rows)
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for f in free:
        x = [Fraction(0)] * ncol
        x[f] = Fraction(1)
        for i, pc in enumerate(piv):
            x[pc] = -R[i][f]
        basis.append(x)
    return basis


def _primitive(vec):
    """Scale a rational vector by a positive rational to a primitive integer vector."""
    from math import gcd
    den = 1
    for x in vec:
        den = den * x.denominator // gcd(den, x.denominator)
    ints = [int(x * den) for x in vec]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    if g == 0:
        return tuple(ints)
    return tuple(x // g for x in ints)


def coordinate_columns(vectors, r):
    """Greedy choice of r coordinates (permutation indices) on which the vectors have rank r."""
    cols = []
    for c in range(24):
        trial = cols + [c]
        if rank([[v[i] for i in trial] for v in vectors]) == len(trial):
            cols = trial
        if len(cols) == r:
            break
    assert len(cols) == r
    return cols


def cone_facets(gens):
    """All facets of cone(gens) inside its linear span.

    gens: list of integer vectors (here: 24 values).  Let d = dim span.  A facet of a d-dimensional
    polyhedral cone is spanned by d-1 linearly independent generators, and d-1 independent
    generators span a facet iff all generators lie weakly on one side of their hyperplane.  So the
    enumeration over all (d-1)-subsets is complete.  A linear functional on the span is determined
    by its values on the generators; each facet is returned as that value tuple (primitive
    integers, >= 0), together with a representative functional given in d chosen coordinates.

    Returns (d, cols, facets) with facets = dict {values_on_gens: normal_in_coords(cols)}.
    """
    d = rank(gens)
    cols = coordinate_columns(gens, d)
    red = [[Fraction(v[c]) for c in cols] for v in gens]
    facets = {}
    for S in itertools.combinations(range(len(gens)), d - 1):
        rows = [red[i] for i in S]
        if rank(rows) != d - 1:
            continue
        ns = nullspace(rows, d)
        assert len(ns) == 1
        nrm = ns[0]
        vals = [sum(a * b for a, b in zip(nrm, g)) for g in red]
        if all(x >= 0 for x in vals):
            pass
        elif all(x <= 0 for x in vals):
            nrm = [-a for a in nrm]
            vals = [-x for x in vals]
        else:
            continue
        key = _primitive(vals)
        scale = None
        for a, b in zip(vals, key):
            if b != 0:
                scale = Fraction(b) / a
                break
        facets[key] = tuple(a * scale for a in nrm)
    return d, cols, facets


def functional_values(coeffs, vectors):
    """coeffs: dict {permutation index: coefficient}; value of sum coeff * g_index on each vector."""
    return [sum(Fraction(c) * v[i] for i, c in coeffs.items()) for v in vectors]


def proportional_positive(u, w):
    """True iff u = t w for some rational t > 0 (u, w sequences of rationals, w != 0)."""
    t = None
    for a, b in zip(u, w):
        if b != 0:
            t = Fraction(a) / Fraction(b)
            break
    if t is None or t <= 0:
        return False
    return all(Fraction(a) == t * Fraction(b) for a, b in zip(u, w))


def minor_cone_report(ck, T, doc_generators, doc_nonminor, doc_coordinate_rays, doc_extra, doc_nrays):
    """Recompute the contact cone of T and the facet description of its minor sub-cone.

    doc_coordinate_rays: list of permutation indices j with 'g_j >= 0' listed as a dual ray
                         (None if the document does not list them).
    doc_extra: dict name -> {index: coeff}, the two non-coordinate inequalities A, B.
    Returns (gens, minor_idx, nonminor_idx, facets, cols) for further use; facets maps the value
    tuple of a facet functional on the minor generators to its coefficient vector in the
    coordinates g_c, c in cols.
    """
    gens = generators(T)
    ck.check("generators F_T recomputed = document's list", gens == list(doc_generators),
             "F_T = %s" % gens)
    nonminor = [j for j in gens if not F_IS_MINOR[j]]
    minor = [j for j in gens if F_IS_MINOR[j]]
    ck.check("non-minor generators = document's", nonminor == list(doc_nonminor),
             ", ".join("f_%d = %s" % (j, f_name(j)) for j in nonminor))
    allv = [F_VALUES[j] for j in gens]
    minv = [F_VALUES[j] for j in minor]
    ck.check("rank of all generators = rank of the minor generators = 8 "
             "(so every source profile lies in the span of the minor generators)",
             rank(allv) == 8 and rank(minv) == 8, "%d minor generators" % len(minor))
    d, cols, facets = cone_facets(minv)
    ck.check("minor cone has exactly %d facets (dual extreme rays)" % doc_nrays,
             len(facets) == doc_nrays, "found %d; coordinates used %s" % (len(facets), cols))
    # value of every facet functional on the non-minor generators (functional known on the span)
    nmv = {}
    for key, nrm in facets.items():
        nmv[key] = [sum(a * F_VALUES[j][c] for a, c in zip(nrm, cols)) for j in nonminor]
    bad = [key for key in facets if min(nmv[key]) < 0]
    ck.check("exactly two facet functionals are negative on some non-minor generator",
             len(bad) == 2)
    # compare with the document's A and B
    matched = set()
    for name, coeffs in doc_extra.items():
        vals = functional_values(coeffs, minv)
        hit = [key for key in bad if proportional_positive(vals, key)]
        ok = len(hit) == 1
        if ok:
            matched.add(hit[0])
        nm = functional_values(coeffs, [F_VALUES[j] for j in nonminor])
        ck.check("document's %s is one of the two (as a functional on the span)" % name, ok,
                 "values on minor generators %s; on non-minor %s"
                 % ([str(x) for x in vals], [str(x) for x in nm]))
    ck.check("the document's A and B are the two distinct exceptional facets", len(matched) == 2)
    good = [key for key in facets if key not in bad]
    if doc_coordinate_rays is not None:
        m2 = set()
        okc = True
        for j in doc_coordinate_rays:
            vals = functional_values({j: 1}, minv)
            hit = [key for key in good if proportional_positive(vals, key)]
            if len(hit) != 1:
                okc = False
            else:
                m2.add(hit[0])
        ck.check("the other %d facets are exactly the document's coordinate inequalities g_j >= 0, j in %s"
                 % (len(good), list(doc_coordinate_rays)),
                 okc and len(m2) == len(good) == len(doc_coordinate_rays))
    else:
        # every remaining facet is >= 0 on all generators; identify those that are coordinates
        ck.check("the other %d facets are >= 0 on all generators of the contact cone" % len(good),
                 all(min(nmv[key]) >= 0 for key in good))
    # Second, logically independent confirmation of the implication that the proofs use
    # ("all facet inequalities hold  =>  nonnegative combination of minor generators"):
    # enumerate the extreme rays of H = {y in span : l(y) >= 0 for every facet functional l}.
    # A functional l is stored as its values on the minor generators, so the family of the l's is
    # a list of vectors of rank d (then H is pointed and is generated by its extreme rays).  An
    # extreme ray y of H is cut out by d-1 independent l's; cone_facets applied to the l's
    # enumerates exactly these and returns each y as the tuple (l(y))_l.
    keys = list(facets.keys())
    d2, _, rays = cone_facets([list(k) for k in keys])
    back = set(rays.keys())
    prim = set(_primitive([Fraction(k[gi]) for k in keys]) for gi in range(len(minor)))
    ck.check("double-dual check: the facet functionals have rank 8 and every extreme ray of "
             "{all facet inequalities} is a positive multiple of a minor generator",
             d2 == 8 and len(back) > 0 and back <= prim,
             "%d extreme rays; %d minor generators" % (len(back), len(minor)))
    return gens, minor, nonminor, facets, cols


# --------------------------------------------------------------------------------------------
# exact rational-function arithmetic
# --------------------------------------------------------------------------------------------


def make_field(names):
    """Fraction field QQ(names).  Returns (K, dict name -> generator)."""
    res = _sym_field(names, QQ)
    K = res[0]
    gens = res[1:]
    return K, dict(zip(names.replace(" ", "").split(","), gens))


def products(A, C, inf_A=None, inf_C=None):
    """The 24 products F_sigma = prod_i [A_i, C_sigma(i)] for labels given in an affine chart.

    A, C: lists of four field elements; the label with index inf_A (resp. inf_C) is the point at
    infinity and its entry is ignored.  With homogeneous representatives (1,0) for infinity and
    (z,1) otherwise, [A_i, C_j] = A_i - C_j for finite labels, [inf, C_j] = 1 and [A_i, inf] = -1.
    """
    out = []
    for s in PERMS:
        val = None
        for i in range(4):
            j = s[i]
            if inf_A is not None and i == inf_A:
                f = 1
            elif inf_C is not None and j == inf_C:
                f = -1
            else:
                f = A[i] - C[j]
            val = f if val is None else val * f
        out.append(val)
    return out


# --------------------------------------------------------------------------------------------
# exact Gaussian rationals (for exact evaluation at sample points)
# --------------------------------------------------------------------------------------------


class GQ:
    """Element of Q(i) with Fraction coordinates."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = Fraction(re)
        self.im = Fraction(im)

    @staticmethod
    def _c(x):
        return x if isinstance(x, GQ) else GQ(x)

    def __add__(self, o):
        o = GQ._c(o)
        return GQ(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __sub__(self, o):
        o = GQ._c(o)
        return GQ(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return GQ._c(o) - self

    def __mul__(self, o):
        o = GQ._c(o)
        return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    def norm(self):
        return self.re * self.re + self.im * self.im

    def inv(self):
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError("GQ division by zero")
        return GQ(self.re / n, -self.im / n)

    def __truediv__(self, o):
        return self * GQ._c(o).inv()

    def __rtruediv__(self, o):
        return GQ._c(o) * self.inv()

    def __eq__(self, o):
        o = GQ._c(o)
        return self.re == o.re and self.im == o.im

    def __hash__(self):
        return hash((self.re, self.im))

    def __repr__(self):
        return "(%s + %s i)" % (self.re, self.im)


# --------------------------------------------------------------------------------------------
# NUMERIC helpers (floating point; falsification attempts only)
# --------------------------------------------------------------------------------------------


def float_products(a, c):
    """F_sigma = prod_i (a_i - c_sigma(i)) for finite complex labels (numpy arrays)."""
    import numpy as np
    S = np.array(PERMS)
    a = np.asarray(a, complex)
    c = np.asarray(c, complex)
    return np.prod(a[None, :] - c[S], axis=1)


def find_contact_point(T, rng, tries=20, tol=1e-11):
    """NUMERIC, independent of every chart of the documents: finite spectra A, C in C^4 and a unit
    complex kappa such that the contact EQUALITIES Re(kappa (F_s - F_t0)) = 0 hold for s in T
    (t0 = min T).  The other gaps are not constrained in sign (the derivation of the charts uses
    only the equalities and distinctness of the labels); they are normalized to mean square 1.
    Least squares from a random start.  Returns (A, C, kappa, gaps) or None."""
    import numpy as np
    from scipy.optimize import least_squares
    T = sorted(T)
    others = [s for s in range(24) if s not in T]

    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])

    def gaps(x):
        a, c, kap = unpack(x)
        F = float_products(a, c)
        return (kap * (F - F[T[0]])).real

    def res(x):
        g = gaps(x)
        r = [g[s] for s in T[1:]]
        r.append(np.mean(g[others] ** 2) - 1.0)
        return np.array(r)

    for _ in range(tries):
        x0 = rng.normal(size=17)
        try:
            sol = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
        except Exception:
            continue
        if np.max(np.abs(res(sol.x))) < tol:
            a, c, kap = unpack(sol.x)
            return a, c, kap, gaps(sol.x)
    return None


def labels_well_separated(a, c, dmin=1e-2, big=1e3):
    import numpy as np
    lab = np.concatenate([np.asarray(a), np.asarray(c)])
    if np.max(np.abs(lab)) > big:
        return False
    return min(abs(lab[i] - lab[j]) for i in range(8) for j in range(i + 1, 8)) > dmin


def cross_ratio(p, q, r, s):
    """(p - r)(q - s) / ((p - s)(q - r))."""
    return (p - r) * (q - s) / ((p - s) * (q - r))


def moebius_to(p_inf, p_zero, p_one):
    """The Moebius map sending p_inf -> infinity, p_zero -> 0, p_one -> 1 (as a function)."""
    def f(x):
        return (x - p_zero) * (p_one - p_inf) / ((x - p_inf) * (p_one - p_zero))
    return f


def minor_cone_nnls_residual(g, minor_idx):
    """NUMERIC: is the 24-vector g a nonnegative combination of the minor generators?
    Returns the least-squares residual of a nonnegative least squares fit (0 means yes)."""
    import numpy as np
    from scipy.optimize import nnls
    M = np.array([F_VALUES[j] for j in minor_idx], float).T
    x, rn = nnls(M, np.asarray(g, float))
    return rn


# --------------------------------------------------------------------------------------------
# tools shared by the three "contact rectangle" cases T_18, T_19, T_20
# --------------------------------------------------------------------------------------------


def rectangle_checks(ck):
    """Identities for eight general finite labels used by the rectangle argument (EXACT).
    Returns the field and the 24 products for further case-specific identities."""
    K, g = make_field("A0,A1,A2,A3,C0,C1,C2,C3")
    A = [g["A%d" % i] for i in range(4)]
    C = [g["C%d" % i] for i in range(4)]
    F = products(A, C)
    ck.check("rectangle: F0 F7 = F1 F6 for general labels", F[0] * F[7] == F[1] * F[6])
    ck.check("F0 - F1 = (A0-C0)(A1-C1)(A2-A3)(C2-C3)   (nonzero for distinct labels)",
             F[0] - F[1] == (A[0] - C[0]) * (A[1] - C[1]) * (A[2] - A[3]) * (C[2] - C[3]))
    ck.check("F0 - F6 = (A2-C2)(A3-C3)(A0-A1)(C0-C1)   (nonzero for distinct labels)",
             F[0] - F[6] == (A[2] - C[2]) * (A[3] - C[3]) * (A[0] - A[1]) * (C[0] - C[1]))
    ck.check("F6/F0 = (A0-C1)(A1-C0)/((A0-C0)(A1-C1))   (a cross-ratio of A0, A1, C0, C1)",
             F[6] / F[0] == (A[0] - C[1]) * (A[1] - C[0]) / ((A[0] - C[0]) * (A[1] - C[1])))
    ck.check("F8/F6 = (A1-C2)(A2-C0)/((A1-C0)(A2-C2))",
             F[8] / F[6] == (A[1] - C[2]) * (A[2] - C[0]) / ((A[1] - C[0]) * (A[2] - C[2])))
    ck.check("F1/F0 = (A2-C3)(A3-C2)/((A2-C2)(A3-C3))",
             F[1] / F[0] == (A[2] - C[3]) * (A[3] - C[2]) / ((A[2] - C[2]) * (A[3] - C[3])))
    # (H + i q0)(H + i q7) = (H + i q1)(H + i q6): real and imaginary parts
    from sympy.polys.rings import ring as _ring
    Rg, Hh, q0, q1, q6, q7, ii = _ring("H,q0,q1,q6,q7,i", QQ)
    diff = (Hh + ii * q0) * (Hh + ii * q7) - (Hh + ii * q1) * (Hh + ii * q6)
    diff = diff.rem(ii ** 2 + 1)
    ck.check("(H+iq0)(H+iq7) - (H+iq1)(H+iq6) = (q1 q6 - q0 q7) + i H (q0 + q7 - q1 - q6)",
             diff == (q1 * q6 - q0 * q7) + ii * Hh * (q0 + q7 - q1 - q6))
    R4, q0, q1, q6 = _ring("q0,q1,q6", QQ)
    q7 = q1 + q6 - q0
    ck.check("with q7 = q1 + q6 - q0:  q1 q6 - q0 q7 = (q0 - q1)(q0 - q6)   (so H != 0 forces F0 = F1 or F0 = F6)",
             q1 * q6 - q0 * q7 == (q0 - q1) * (q0 - q6))
    return K, A, C, F


def certified_sign(elem):
    """Sign certificate for a rational function in variables that are all > 0: returns +1 (-1) if
    numerator and denominator each have coefficients of one sign and the quotient is positive
    (negative) on the whole open positive orthant; 0 if this test is inconclusive."""
    def one_sign(poly):
        cs = [c for _, c in poly.terms()]
        if not cs:
            return 0
        if all(c > 0 for c in cs):
            return 1
        if all(c < 0 for c in cs):
            return -1
        return 0
    a = one_sign(elem.numer)
    b = one_sign(elem.denom)
    return a * b


def eval_frac(elem, names, point):
    """Evaluate an element of Q(names) at GQ values (dict name -> GQ); exact."""
    def ev(poly):
        tot = GQ(0)
        for mon, coef in poly.terms():
            term = GQ(Fraction(int(coef.numerator), int(coef.denominator)))
            for nm, e in zip(names, mon):
                for _ in range(e):
                    term = term * point[nm]
            tot = tot + term
        return tot
    return ev(elem.numer) / ev(elem.denom)


def gq_conj(x):
    return GQ(x.re, -x.im)


def exact_ratio_gaps(A, C, mu):
    """Gaps g_j = -mu Im(F_j/F_0) for labels A (finite GQ) and C = (infinity, C1, C2, C3)."""
    Fq = []
    for sgm in PERMS:
        val = GQ(1)
        for i in range(4):
            j = sgm[i]
            val = val * (GQ(-1) if j == 0 else (A[i] - C[j]))
        Fq.append(val)
    return [-(mu * (f / Fq[0]).im) for f in Fq]


def float_ratio_gaps(A, C, mu):
    """Float version of exact_ratio_gaps (numpy array)."""
    import numpy as np
    Ff = []
    for sgm in PERMS:
        val = 1.0 + 0j
        for i in range(4):
            j = sgm[i]
            val = val * (-1.0 if j == 0 else (A[i] - C[j]))
        Ff.append(val)
    return np.array([-(mu * (f / Ff[0]).imag) for f in Ff])


def substitute(elem, gen_names, val, Kp):
    """Substitute val (dict generator name -> element of the field Kp) into an element of a
    fraction field with generators gen_names.  Raises KeyError if elem involves a generator that
    is not in val (so the caller learns that the expression was not free of it)."""
    def ev(poly):
        tot = Kp(0)
        for mon, cf in poly.terms():
            term = Kp(Fraction(int(cf.numerator), int(cf.denominator)))
            for nm, ex in zip(gen_names, mon):
                if ex:
                    term = term * val[nm] ** ex
            tot = tot + term
        return tot
    return ev(elem.numer) / ev(elem.denom)


CELL_LABELS = ["<0", "(0,1)", ">1"]


def cell_signs(exprs, gen_names, cell_names, Kp, pos):
    """Certified signs of the expressions on each of the 3^n cells of the real variables
    cell_names (each either < 0, in (0,1) or > 1).  A cell is parametrized bijectively by positive
    variables: t = -p, t = 1/(1+p), t = 1+p; certified_sign then decides the sign on the whole
    cell.  Returns a list of (dict name -> cell label, [signs])."""
    forms = {"<0": lambda p: -p, "(0,1)": lambda p: 1 / (1 + p), ">1": lambda p: 1 + p}
    out = []
    for combo in itertools.product(CELL_LABELS, repeat=len(cell_names)):
        val = {nm: forms[cl](pv) for nm, cl, pv in zip(cell_names, combo, pos)}
        signs = [certified_sign(substitute(e, gen_names, val, Kp)) for e in exprs]
        out.append((dict(zip(cell_names, combo)), signs))
    return out
