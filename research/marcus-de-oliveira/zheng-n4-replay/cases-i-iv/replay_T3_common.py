"""Shared exact tools for the independent replays of T_3, T_4, T_5, T_6 (group III).

Everything here is built from the DEFINITIONS in the paper (Sections 3.1-3.2, 4.1-4.2,
Appendix D), by own code:

  * permutations in lexicographic order, pi_j = PERMS[j] (0-based one-line notation);
  * the 82 facet vectors f_j as their 24 VALUES v_sigma(f_j) (Appendix D transcribed here
    from the paper's text; the matrix P is not needed because P is injective and every
    Leibniz sum is a function of the values);
  * the group G of order 1152 acting on permutation indices, value vectors and matrices;
  * exact linear algebra over Fraction, cone facets by enumeration;
  * Leibniz sums, 2x2 minors, random exact rational SO(4) matrices;
  * the quaternion model x -> p x q of SO(4) as exact polynomials (sympy rings).

Nothing of the supplement is executed or imported.  `supplement_facet_values()` reads one
JSON file of the supplement AS DATA (only if MOC4_SUPP is set) to confirm that the
documents' facet numbering is the paper's Appendix D numbering.

`crosscheck_audit_lib()` compares the objects built here with the audit's lib/moc4lib.py.

`verify_contact_cone()` enumerates a contact cone directly; `GQ` and `exact_paper_gaps()`
compute gap vectors of Gaussian-rational spectra exactly from the paper's definition.

Python packages: sympy (exact rings/fields), python-flint (integer nullspaces in
verify_contact_cone), numpy (falsification attempts only, in the replay scripts).
"""
import itertools
import json
import os
import sys
from fractions import Fraction as Fr

# ---------------------------------------------------------------------------------------
# 0. bookkeeping
# ---------------------------------------------------------------------------------------


class Checker:
    """Prints every check and remembers failures; finish() exits 0 iff all passed."""

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.failed = []
        print("=" * 78)
        print(name)
        print("=" * 78)

    def section(self, title):
        print()
        print("--- " + title)

    def check(self, label, ok, detail=""):
        self.n += 1
        ok = bool(ok)
        tag = "ok  " if ok else "FAIL"
        print("[%s] %s%s" % (tag, label, (" :: " + str(detail)) if detail != "" else ""))
        if not ok:
            self.failed.append(label)
        return ok

    def info(self, text):
        print("       " + str(text))

    def finish(self):
        print()
        if self.failed:
            print("FAIL: %d of %d checks failed in %s" % (len(self.failed), self.n, self.name))
            for f in self.failed:
                print("   failed: " + f)
            sys.exit(1)
        print("PASS: all %d checks passed in %s" % (self.n, self.name))
        sys.exit(0)


# ---------------------------------------------------------------------------------------
# 1. permutations
# ---------------------------------------------------------------------------------------
PERMS = [tuple(p) for p in itertools.permutations(range(4))]      # lexicographic
PIDX = {p: i for i, p in enumerate(PERMS)}


def sgn(p):
    inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
    return -1 if inv % 2 else 1


SGN = [sgn(p) for p in PERMS]


def compose(p, q):
    """(p o q)(i) = p(q(i))."""
    return tuple(p[q[i]] for i in range(4))


def inverse(p):
    r = [0] * 4
    for i, x in enumerate(p):
        r[x] = i
    return tuple(r)


# ---------------------------------------------------------------------------------------
# 2. the 82 facet vectors, as value vectors (paper Section 4.1 and Appendix D)
# ---------------------------------------------------------------------------------------
# Own transcription of Appendix D (paper p. 59): index j -> name.  "nI,J": minor vector with
# value 1{sigma(I) = J} (digits are 1-based).  "hPQRS,e": h_{pi,eps} with pi in one-line
# notation, e = + or -.
APP_D = (
    "n1,4 n12,34 n13,34 n14,34 n1,3 n12,24 n4,1 n3,1 n12,23 n14,23 "          # 0-9
    "n13,24 n13,23 n14,24 n2,1 n1,2 n12,14 n2,4 h3421,+ h3412,- h2413,- "     # 10-19
    "h2431,+ n4,2 h3412,+ n14,12 h4312,+ h3142,- h4132,- n3,2 h3421,- n13,12 "  # 20-29
    "h4321,- h4123,+ h3124,+ n12,13 h4321,+ h4312,- n2,3 h2341,+ h2314,- n14,13 "  # 30-39
    "h4213,+ n4,3 h4123,- h2413,+ h2143,- n13,14 h3241,- n3,4 h3142,+ h2143,+ "  # 40-49
    "h2341,- h1423,- h1342,+ n13,13 h4231,- n3,3 h4132,+ h2431,- h2134,+ n14,14 "  # 50-59
    "h3214,+ n4,4 h3124,- h2134,- h2314,+ h1432,+ h1324,- n2,2 h4213,- h3241,+ "  # 60-69
    "h4231,+ h3214,- n12,12 h1243,+ h1432,- h1423,+ h1234,- h1342,- h1234,+ h1243,- "  # 70-79
    "h1324,+ n1,1"                                                             # 80-81
).split()
assert len(APP_D) == 82


def _ind(I, J):
    """Value vector of 1{sigma(I) = J}; I, J are 0-based index collections."""
    I, J = tuple(I), frozenset(J)
    return [1 if frozenset(s[i] for i in I) == J else 0 for s in PERMS]


def n_values(Istr, Jstr):
    return tuple(_ind([int(c) - 1 for c in Istr], [int(c) - 1 for c in Jstr]))


def h_values(pistr, eps):
    """h_{pi,eps} = (t_pi + eps d_pi)/2 with
    t_pi = sum_i n_{i,pi(i)},
    d_pi = n_{12,{pi1,pi3}} + n_{12,{pi2,pi4}} - n_{13,{pi1,pi2}} - n_{13,{pi3,pi4}}."""
    pi = tuple(int(c) - 1 for c in pistr)
    t = [sum(1 for i in range(4) if s[i] == pi[i]) for s in PERMS]
    a = _ind((0, 1), (pi[0], pi[2]))
    b = _ind((0, 1), (pi[1], pi[3]))
    c = _ind((0, 2), (pi[0], pi[1]))
    d = _ind((0, 2), (pi[2], pi[3]))
    out = []
    for k in range(24):
        num = t[k] + eps * (a[k] + b[k] - c[k] - d[k])
        assert num % 2 == 0
        out.append(num // 2)
    return tuple(out)


def facet_values(name):
    if name[0] == "n":
        I, J = name[1:].split(",")
        return n_values(I, J)
    pi, e = name[1:].split(",")
    return h_values(pi, 1 if e == "+" else -1)


FVAL = [facet_values(nm) for nm in APP_D]            # FVAL[j][s] = v_{pi_s}(f_j)
IS_MINOR = [nm[0] == "n" for nm in APP_D]
FIDX = {v: j for j, v in enumerate(FVAL)}
assert len(FIDX) == 82


def minor_spec(j):
    """(I, J) as 0-based tuples for a minor facet f_j."""
    assert IS_MINOR[j]
    I, J = APP_D[j][1:].split(",")
    return tuple(int(c) - 1 for c in I), tuple(int(c) - 1 for c in J)


def generators(T):
    """Indices j with f_j = 0 at every permutation index of T (the set F_T)."""
    return [j for j in range(82) if all(FVAL[j][s] == 0 for s in T)]


# ---------------------------------------------------------------------------------------
# 3. the group G = (S4 x S4) x| C2 of order 1152
# ---------------------------------------------------------------------------------------
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


def g_facet(gmap, j):
    """Index of the facet g.f_j (KeyError if G did not permute the facets)."""
    return FIDX[g_values(gmap, FVAL[j])]


def g_set(gmap, T):
    return frozenset(gmap[j] for j in T)


def find_group_elements(rep, target):
    """All indices k of G with g_k . f_{rep[m]} = f_{target[m]} for every position m."""
    out = []
    for k, gm in enumerate(G_MAPS):
        if all(g_facet(gm, a) == b for a, b in zip(rep, target)):
            out.append(k)
    return out


def orbit_of_tuple(rep):
    return {tuple(g_facet(gm, a) for a in rep) for gm in G_MAPS}


def g_matrix(g, O):
    """g.O with Phi_{g.c}(g.O) = (+-1) Phi_c(O): for e = 0, (g.O)[i][j] = O[r2(i)][r1^{-1}(j)];
    for e = 1 the matrix is transposed first."""
    r1, r2, e = g
    if e:
        O = [[O[j][i] for j in range(4)] for i in range(4)]
    r1i = inverse(r1)
    return [[O[r2[i]][r1i[j]] for j in range(4)] for i in range(4)]


# ---------------------------------------------------------------------------------------
# 4. exact linear algebra over Q
# ---------------------------------------------------------------------------------------


def rref(rows):
    A = [[Fr(x) for x in r] for r in rows]
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


def nullspace(rows):
    R, piv = rref(rows)
    nc = len(rows[0])
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for f in free:
        x = [Fr(0)] * nc
        x[f] = Fr(1)
        for i, p in enumerate(piv):
            x[p] = -R[i][f]
        basis.append(x)
    return basis


def primitive(vec):
    """Positive rational multiple of vec with coprime integer entries."""
    from math import gcd
    den = 1
    for x in vec:
        den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    ints = [int(Fr(x) * den) for x in vec]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    return tuple(x // g for x in ints) if g else tuple(ints)


def cone_facets(gens):
    """Facet normals of the cone generated by `gens` (vectors in Q^d spanning Q^d):
    all primitive integer l with l.g >= 0 on every generator and with the generators
    in {l = 0} spanning a hyperplane.  Exhaustive over (d-1)-subsets, exact."""
    d = len(gens[0])
    assert rank(gens) == d
    found = set()
    for S in itertools.combinations(range(len(gens)), d - 1):
        ns = nullspace([gens[i] for i in S])
        if len(ns) != 1:
            continue
        l = ns[0]
        vals = [sum(a * b for a, b in zip(l, g)) for g in gens]
        if all(x >= 0 for x in vals):
            found.add(primitive(l))
        elif all(x <= 0 for x in vals):
            found.add(primitive([-a for a in l]))
    return sorted(found)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


# ---------------------------------------------------------------------------------------
# 5. Leibniz sums, minors, exact rational SO(4)
# ---------------------------------------------------------------------------------------


def leibniz_weights(O):
    out = []
    for p, s in zip(PERMS, SGN):
        out.append(s * O[0][p[0]] * O[1][p[1]] * O[2][p[2]] * O[3][p[3]])
    return out


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


def mat_inv(M):
    n = len(M)
    aug = [list(M[i]) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    R, piv = rref(aug)
    assert piv[:n] == list(range(n))
    return [row[n:] for row in R]


def random_rational_SO4(rng, size=4):
    """(I - S)(I + S)^{-1} for a random rational skew S: exactly special orthogonal."""
    while True:
        S = [[Fr(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i + 1, 4):
                val = Fr(rng.randint(-size, size), rng.randint(1, size))
                S[i][j], S[j][i] = val, -val
        A = [[Fr(int(i == j)) - S[i][j] for j in range(4)] for i in range(4)]
        B = [[Fr(int(i == j)) + S[i][j] for j in range(4)] for i in range(4)]
        O = matmul(A, mat_inv(B))
        if all(x != 0 for row in O for x in row):
            return O


def is_SO4(O):
    T = [[O[j][i] for j in range(4)] for i in range(4)]
    return matmul(O, T) == [[Fr(int(i == j)) for j in range(4)] for i in range(4)] and det4(O) == 1


def transfer_to_SO4(g, O):
    """For O in SO(4) return O'' in SO(4) with Phi_{g.c}(O'') = Phi_c(O) for all c:
    g.O, with the first column negated when det(g.O) = -1."""
    M = g_matrix(g, O)
    if det4(M) == -1:
        M = [[-row[0]] + list(row[1:]) for row in M]
    return M


# ---------------------------------------------------------------------------------------
# 6. quaternion model of SO(4)
# ---------------------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------------------
# 7. cross-checks against the supplement's numbering (data) and the audit library
# ---------------------------------------------------------------------------------------


def supplement_facet_values():
    """24-value vectors of the supplement's facet list (read as DATA), or None if
    MOC4_SUPP is not set / the file is missing."""
    root = os.environ.get("MOC4_SUPP")
    if not root:
        return None
    path = os.path.join(root, "data", "root_permutation_polytope_complete_20261022.json")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        data = json.load(fh)
    perms = [tuple(x - 1 for x in p) for p in data["permutations_1_based"]]
    return perms, [tuple(f["values"]) for f in data["facets"]]


def crosscheck_numbering(ck):
    """Own Appendix-D transcription versus (a) the supplement's data numbering used by the
    documents and (b) the audit library."""
    ck.section("Facet numbering: own transcription of Appendix D vs supplement data vs audit lib")
    prof = sorted(tuple(sorted(v)) for v in FVAL)
    ck.check("82 distinct value vectors; 34 minor (16 with 18 zeros, 18 with 20 zeros), 48 h with profile 13x0, 10x1, 1x2",
             len(set(FVAL)) == 82 and sum(IS_MINOR) == 34
             and sorted(v.count(0) for j, v in enumerate(FVAL) if IS_MINOR[j]) == [18] * 16 + [20] * 18
             and all((v.count(0), v.count(1), v.count(2)) == (13, 10, 1)
                     for j, v in enumerate(FVAL) if not IS_MINOR[j]))
    okG = len(set(G_MAPS)) == 1152
    for gm in G_MAPS:
        img = [FIDX.get(g_values(gm, FVAL[j])) for j in range(82)]
        okG = okG and None not in img and sorted(img) == list(range(82)) \
            and all(IS_MINOR[j] == IS_MINOR[img[j]] for j in range(82))
    ck.check("the 1152 elements of G act as distinct permutations of the 24 indices, permute the 82 facet vectors and "
             "map minors to minors", okG)
    sup = supplement_facet_values()
    if sup is None:
        ck.info("MOC4_SUPP not set or data file missing: supplement numbering NOT compared")
    else:
        perms, vals = sup
        ck.check("supplement data: permutations listed lexicographically", perms == PERMS)
        ck.check("supplement data: facet j has the values of the paper's f_j for all 82 j "
                 "(documents' numbering = Appendix D numbering)", vals == FVAL)
    lib = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "..", "zheng-n4-audit", "lib"))
    if os.path.isfile(os.path.join(lib, "moc4lib.py")):
        sys.path.insert(0, lib)
        try:
            import moc4lib as L
            ck.check("audit lib: PERMS identical", [tuple(p) for p in L.PERMS] == PERMS)
            ck.check("audit lib: F_VALUES identical to own values for all 82 facets",
                     [tuple(v) for v in L.F_VALUES] == FVAL)
            ck.check("audit lib: its 1152 group index maps are the same set as mine",
                     set(L._G_MAPS) == set(G_MAPS) and len(set(G_MAPS)) == 1152)
            Pcols = [[L.P[s_][a_] for s_ in range(24)] for a_ in range(14)]
            mins = [FVAL[j] for j in range(82) if IS_MINOR[j]]
            ck.check("audit lib: the column space of its matrix P (rank 14) is the span of my 34 minor value vectors",
                     rank(Pcols) == 14 and rank(mins) == 14 and rank(Pcols + mins) == 14)
            for nu, Tset in ((3, {0, 1, 2, 3, 7, 8}), (4, {0, 1, 4, 5, 7, 8}),
                             (5, {1, 2, 3, 4, 6, 8}), (6, {0, 1, 2, 4, 6, 8})):
                ck.check("audit lib: CONTACT_SETS[%d] = %s and generators agree" % (nu, sorted(Tset)),
                         set(L.CONTACT_SETS[nu]) == Tset and L.generators(Tset) == generators(Tset))
        finally:
            sys.path.pop(0)
    else:
        ck.info("audit lib not found: not compared")


def orbit_check(ck, family, rep, target, sym_positions=(0, 1), expected_orbit=None):
    """Check that `target` (tuple of facet indices with roles) is in the G-orbit of `rep`,
    where the roles in `sym_positions` may be exchanged (the inequality is symmetric in
    them).  Returns a group element g (as tuple) with g.rep = target (possibly after the
    exchange) or None."""
    reps = [tuple(rep)]
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
    ck.check("%s: roles %s lie in the G-orbit of the proved representative %s"
             % (family, target, tuple(rep)), hit is not None,
             "" if hit is None else "g = (rho=%s, rho'=%s, transpose=%d), %d group elements, source roles %s"
             % (hit[0][0], hit[0][1], hit[0][2], hit[2], hit[1]))
    if expected_orbit is not None:
        orb = orbit_of_tuple(rep)
        unordered = {tuple(sorted((t[i], t[j]))) + tuple(t[m] for m in range(len(t)) if m not in (i, j))
                     for t in orb}
        ck.check("%s: orbit size (pair in positions %s unordered) = %d" % (family, sym_positions, expected_orbit),
                 len(unordered) == expected_orbit, len(unordered))
    return None if hit is None else hit[0]


def covariance_check(ck, label, g, pairs, rng, trials=6):
    """Exact: for random rational O in SO(4), O'' = transfer_to_SO4(g, O) is in SO(4) and
    Phi_{f_target}(O'') = Phi_{f_rep}(O) for each (rep, target) in pairs."""
    ok = True
    for _ in range(trials):
        O = random_rational_SO4(rng)
        O2 = transfer_to_SO4(g, O)
        ok = ok and is_SO4(O) and is_SO4(O2)
        for a, b in pairs:
            ok = ok and phi(FVAL[a], O) == phi(FVAL[b], O2)
    ck.check("%s: Phi_{g.f}(O'') = Phi_f(O) with O'' in SO(4), exact on %d random rational O"
             % (label, trials), ok)


# ---------------------------------------------------------------------------------------
# 8. replays of the two universal inequalities used by T_3 .. T_6
# ---------------------------------------------------------------------------------------


def replay_SOC_representative(ck):
    """SOC representative (f1, f9; f22):  a f1 + b f9 + d f22 has Phi >= 0 on SO(4) when
    a, b, d >= 0 and ab >= d^2.  Statement and proof: audit, foundations/so4-atoms.md 3.1.
    Here: every polynomial identity that proof uses, exactly, in the model x -> p x conj(q)."""
    from sympy.polys.rings import ring
    from sympy.polys.domains import ZZ
    ck.section("SOC representative (f1, f9; f22): identities of the audit's proof (so4-atoms.md 3.1), own code")
    R, p0, p1, p2, p3, q0, q1, q2, q3 = ring("p0,p1,p2,p3,q0,q1,q2,q3", ZZ)
    P, Q = (p0, p1, p2, p3), (q0, q1, q2, q3)
    pp, qq = sum(x * x for x in P), sum(x * x for x in Q)
    M = quat_matrix(P, Q, R(1), R(0), conj_q=True)
    MT = [[M[j][i] for j in range(4)] for i in range(4)]
    MM = matmul(MT, M)
    ck.check("model x -> p x conj(q): M^t M = |p|^2|q|^2 I and det M = (|p|^2|q|^2)^2 (so M in SO(4) for unit p, q)",
             all(MM[i][j] == (pp * qq if i == j else 0) for i in range(4) for j in range(4))
             and det4(M) == (pp * qq) ** 2)
    a = [x * x for x in P]
    A = [x * x for x in Q]
    vh = [a[i] * qq - A[i] * pp for i in range(4)]          # v_i = a_i - A_i, homogenised
    x1 = minor_det(M, *minor_spec(1))
    x9 = minor_det(M, *minor_spec(9))
    ck.check("Phi_{f1} = x1^2, Phi_{f9} = x9^2 with x1 = det M[12,34], x9 = det M[14,23] (polynomial identities)",
             phi(FVAL[1], M) == x1 ** 2 and phi(FVAL[9], M) == x9 ** 2)
    ck.check("x1 = v0 + v1 and x9 = v0 + v3 (homogenised with |p|^2, |q|^2)",
             x1 == vh[0] + vh[1] and x9 == vh[0] + vh[3])
    rhs = (2 * (a[0] * A[2] + a[2] * A[0] + a[1] * A[3] + a[3] * A[1]) * pp * qq
           - 4 * p0 * p1 * p2 * p3 * qq ** 2 - 4 * q0 * q1 * q2 * q3 * pp ** 2)
    ck.check("Phi_{f22} = 2(a0A2 + a2A0 + a1A3 + a3A1) - 4 p0p1p2p3 - 4 q0q1q2q3 on S^3 x S^3 (homogenised)",
             phi(FVAL[22], M) == rhs)
    # the supplement's scalar-compensation note, Theorem 3, is written in the model x -> p x q:
    M2 = quat_matrix(P, Q, R(1), R(0), conj_q=False)
    S3 = (p2 ** 2 + p3 ** 2) * qq - (q2 ** 2 + q3 ** 2) * pp
    T3 = (p1 ** 2 + p2 ** 2) * qq - (q1 ** 2 + q2 ** 2) * pp
    H3 = (2 * sum(P[i] ** 2 * Q[(i + 2) % 4] ** 2 for i in range(4)) * pp * qq
          - 4 * p0 * p1 * p2 * p3 * qq ** 2 + 4 * q0 * q1 * q2 * q3 * pp ** 2)
    y1, y9 = minor_det(M2, *minor_spec(1)), minor_det(M2, *minor_spec(9))
    ck.check("bridge for the supplement's Theorem 3 (scalar-compensation note; not in any PDF), model x -> p x q: "
             "Phi_{f1} = S^2, Phi_{f9} = T^2, Phi_{f22} = H (homogenised)",
             y1 ** 2 == S3 ** 2 and y9 ** 2 == T3 ** 2 and phi(FVAL[1], M2) == S3 ** 2
             and phi(FVAL[9], M2) == T3 ** 2 and phi(FVAL[22], M2) == H3)
    # scalar theorem: with a_i = al_i^2, A_i = be_i^2 (al, be >= 0)
    S, a0, a1, a2, a3, b0, b1, b2, b3 = ring("al0,al1,al2,al3,be0,be1,be2,be3", ZZ)
    al, be = (a0, a1, a2, a3), (b0, b1, b2, b3)
    aa, AA = [x * x for x in al], [x * x for x in be]
    v = [aa[i] - AA[i] for i in range(4)]
    Rr = (2 * (aa[0] * AA[2] + aa[1] * AA[3] + aa[2] * AA[0] + aa[3] * AA[1])
          - 4 * a0 * a1 * a2 * a3 - 4 * b0 * b1 * b2 * b3)
    de = a0 * a2 - a1 * a3
    De = b0 * b2 - b1 * b3
    ck.check("scalar step 1: R = 2 delta^2 + 2 Delta^2 - 2(v0v2 + v1v3) identically",
             Rr == 2 * de ** 2 + 2 * De ** 2 - 2 * (v[0] * v[2] + v[1] * v[3]))
    s_ = v[2] + v[3]
    t_ = v[1] + v[2]
    summ = v[0] + v[1] + v[2] + v[3]
    ck.check("scalar step 1: st = (v2+v3)(v1+v2) = v1v3 - v0v2 + v2 sum(v), and (v0+v1)(v0+v3) = st + (v0-v2) sum(v)",
             s_ * t_ == v[1] * v[3] - v[0] * v[2] + v[2] * summ
             and (v[0] + v[1]) * (v[0] + v[3]) == s_ * t_ + (v[0] - v[2]) * summ)
    U, r_, s2, u_, w_ = ring("r,s,u,w", ZZ)
    ck.check("increment lemma: (r^2+u^2)(s^2+w^2) - (rs + uw)^2 = (rw - us)^2",
             (r_ ** 2 + u_ ** 2) * (s2 ** 2 + w_ ** 2) - (r_ * s2 + u_ * w_) ** 2 == (r_ * w_ - u_ * s2) ** 2)
    ck.info("PROOF-READ (scalar-compensation note, Theorem 2 = so4-atoms.md 3.1 steps 1-5): R + 2|st| = 2delta^2 + 2Delta^2 - 4 min(v0v2, v1v3);")
    ck.info("  if the min is > 0, after a <-> A, Delta - delta >= sqrt(alpha gamma) + sqrt(MN), so 2delta^2 + 2Delta^2 >= (Delta-delta)^2")
    ck.info("  >= 4 sqrt(alpha gamma M N) >= 4 min.  With -4p0p1p2p3 >= -4 sqrt(a0a1a2a3): Phi_22 + 2|x1 x9| >= 0 on SO(4);")
    ck.info("  then a x1^2 + b x9^2 + d Phi_22 >= d(2|x1x9| + Phi_22) >= 0 for ab >= d^2, d >= 0.")


def replay_four_coefficient_note(ck):
    """The four-coefficient inequality for its representative (f0, f6; f10; f54), i.e.
    Theorem 1 of the supplement's final_referee_T5_four_budget_capacity note, replayed:
    the bridge to Leibniz sums in the model x -> p x q, the certificate identity (4) modulo
    N = |p|^2 - |q|^2, and the determinant identities behind 'Z is PSD iff capacity'."""
    from sympy.polys.rings import ring
    from sympy.polys.domains import ZZ
    ck.section("Four-coefficient representative (f0, f6; f10; f54): replay of the four-budget note, own code")
    R, p0, p1, p2, p3, q0, q1, q2, q3, a, b, c, d = ring("p0,p1,p2,p3,q0,q1,q2,q3,a,b,c,d", ZZ)
    P, Q = (p0, p1, p2, p3), (q0, q1, q2, q3)
    pp, qq = sum(x * x for x in P), sum(x * x for x in Q)
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
    M = quat_matrix(P, Q, R(1), R(0), conj_q=False)
    MT = [[M[j][i] for j in range(4)] for i in range(4)]
    MM = matmul(MT, M)
    ck.check("model x -> p x q: M^t M = |p|^2|q|^2 I and det M = (|p|^2|q|^2)^2",
             all(MM[i][j] == (pp * qq if i == j else 0) for i in range(4) for j in range(4))
             and det4(M) == (pp * qq) ** 2)
    quart = {0: (A - B) ** 2, 6: (A + B) ** 2, 10: D ** 2, 54: H}
    for j in (0, 6, 10, 54):
        ck.check("bridge: 4 Phi_{f%d}(M(p,q)) = n^2 Q mod N, Q = %s"
                 % (j, {0: "(A-B)^2", 6: "(A+B)^2", 10: "D^2", 54: "H"}[j]),
                 red(4 * phi(FVAL[j], M)) == red(n * n * quart[j]))
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

    S = [y,
         (-p1, p0, -p3, p2, q1, -q0, -q3, q2),
         (-p2, p3, p0, -p1, q2, q3, -q0, -q1),
         (-p3, -p2, p1, p0, q3, -q2, q1, -q0)]
    rhs = R(0)
    for Sy in S:
        l, m = lm(Sy)
        rhs += qform(l) + qform(m)
    lhs = n * (a * (A - B) ** 2 + b * (A + B) ** 2 + c * D ** 2 + d * H)
    ck.check("certificate (4): n[a(A-B)^2 + b(A+B)^2 + cD^2 + dH] = sum_j l(S_j y)^t Z l(S_j y) + m(S_j y)^t Z m(S_j y) mod N",
             red(lhs) == red(rhs))
    ck.check("each S_j preserves |p|^2 and |q|^2 (signed coordinate permutations within p and within q)",
             all(sum(t * t for t in Sy[:4]) == pp and sum(t * t for t in Sy[4:]) == qq for Sy in S))
    # Z >= 0  <=>  capacity (for a, b, c, d >= 0): the Schur-complement identities
    cap = a * b * (c - d) - d * d * (a + b + c + d)
    S3 = [[a + b + 2 * d, a - b, R(0)], [a - b, a + b - 2 * d, -2 * d], [R(0), -2 * d, c - d]]
    ck.check("Schur complement of the entry Z44 = d: d * (Z[:3,:3] - col col^t / d) = d * S3 with S3 as in the note",
             all(d * Z[i][j] - Z[i][3] * Z[3][j] == d * S3[i][j] for i in range(3) for j in range(3)))
    det2 = (a + b + 2 * d) * ((a + b - 2 * d) * (c - d) - 4 * d * d) - (a - b) ** 2 * (c - d)
    ck.check("second Schur complement (of c - d): (c-d) * det of the 2x2 block = 4[ab(c-d) - d^2(a+b+c+d)]",
             det2 == 4 * cap)

    def det_sym(Mx):
        nn = len(Mx)
        if nn == 1:
            return Mx[0][0]
        tot = R(0)
        for k in range(nn):
            if Mx[0][k] == 0:
                continue
            sub = [[Mx[i][j] for j in range(nn) if j != k] for i in range(1, nn)]
            tot += (-1) ** k * Mx[0][k] * det_sym(sub)
        return tot

    ck.check("det Z = 4 d [ab(c-d) - d^2(a+b+c+d)]", det_sym(Z) == 4 * d * cap)
    ck.info("PROOF-READ: d = 0: Z = [[a+b, a-b],[a-b, a+b]] + [c] + [0] is PSD.  d > 0 and capacity: ab(c-d) >= d^2(..) > 0")
    ck.info("  forces c > d; Z44 = d > 0, S3 has (3,3) entry c-d > 0, the 2x2 complement has positive (1,1) entry and det >= 0.")
    ck.info("  So Z >= 0, each quadratic form in (4) is >= 0, and n = 2 on S^3 x S^3: a Phi_0 + b Phi_6 + c Phi_10 + d Phi_54 >= 0.")
    ck.info("  CITED (classical): every O in SO(4) is the matrix of x -> p x q for unit quaternions p, q.")


def infinity_convention_check(ck, kind, idx):
    """Exact: sending the point A_idx (kind 'A') or C_idx (kind 'C') to infinity by
    s -> 1/(s - P) and dropping its factor changes all 24 products by ONE common nonzero
    rational factor.  Symbols are generic finite points."""
    from sympy.polys.fields import field
    from sympy.polys.domains import QQ
    K, a0, a1, a2, a3, c0, c1, c2, c3 = field("a0,a1,a2,a3,c0,c1,c2,c3", QQ)
    A, Cc = [a0, a1, a2, a3], [c0, c1, c2, c3]
    Pinf = A[idx] if kind == "A" else Cc[idx]
    mA = [None if (kind == "A" and i == idx) else 1 / (A[i] - Pinf) for i in range(4)]
    mC = [None if (kind == "C" and j == idx) else 1 / (Cc[j] - Pinf) for j in range(4)]
    fin, cha = [], []
    for s in PERMS:
        f, g = K(1), K(1)
        for i in range(4):
            f = f * (A[i] - Cc[s[i]])
            if mA[i] is None or mC[s[i]] is None:
                continue
            g = g * (mA[i] - mC[s[i]])
        fin.append(f)
        cha.append(g)
    fac = fin[0] / cha[0]
    ck.check("point %s_%d at infinity with its factor dropped: all 24 products change by one common factor"
             % (kind, idx), all(cha[k] * fac == fin[k] for k in range(24)))


# ---------------------------------------------------------------------------------------
# 9. exact instances: gap vectors from FINITE Gaussian-rational spectra (paper (3.1)-(3.3))
# ---------------------------------------------------------------------------------------


class GQ:
    """Exact Gaussian rational re + i im."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re, self.im = Fr(re), Fr(im)

    @staticmethod
    def of(x):
        return x if isinstance(x, GQ) else GQ(x, 0)

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

    def __truediv__(self, o):
        o = GQ.of(o)
        n = o.re * o.re + o.im * o.im
        return self * GQ(o.re / n, -o.im / n)

    def __rtruediv__(self, o):
        return GQ.of(o) / self

    def __eq__(self, o):
        o = GQ.of(o)
        return self.re == o.re and self.im == o.im

    def __hash__(self):
        return hash((self.re, self.im))

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = GQ(1)
        for _ in range(n):
            out = out * self
        return out

    def conj(self):
        return GQ(self.re, -self.im)

    def norm2(self):
        return self.re * self.re + self.im * self.im

    def __repr__(self):
        return "(%s%+si)" % (self.re, self.im)


def exact_paper_gaps(Apts, Cpts, ref, other):
    """Apts, Cpts: four points of P^1(C) each, as GQ or None (= infinity), all distinct.
    A fixed Moebius map s -> 1/(s - s0) makes all eight finite; then a_i = m(A_i),
    b_j = -m(C_j), z_sigma = prod_i (a_i + b_sigma(i)) as in the paper, and the gap vector is
    Re(alpha (z_sigma - z_ref)) for the direction alpha = +- i conj(z_other - z_ref), the
    unique real line of directions making `other` a contact (`other` may be a tuple of
    candidates; the first with z_other != z_ref is used); the sign is chosen so that the
    gaps sum to a positive number.  Returns (a, b, alpha, gaps) with exact Fractions."""
    s0 = GQ(Fr(7, 3), Fr(5, 7))
    pts = [p for p in list(Apts) + list(Cpts) if p is not None]
    assert all(not (p == s0) for p in pts)
    mob = lambda p: GQ(0) if p is None else 1 / (GQ.of(p) - s0)
    a = [mob(p) for p in Apts]
    b = [-mob(p) for p in Cpts]
    assert len({(x.re, x.im) for x in a}) == 4 and len({(x.re, x.im) for x in b}) == 4
    zs = []
    for s in PERMS:
        p = GQ(1)
        for i in range(4):
            p = p * (a[i] + b[s[i]])
        zs.append(p)
    if isinstance(other, int):
        other = (other,)
    other = next(k for k in other if not (zs[k] == zs[ref]))
    alpha = GQ(0, 1) * (zs[other] - zs[ref]).conj()
    gaps = [(alpha * (zz - zs[ref])).re for zz in zs]
    if sum(gaps) < 0:
        alpha = -alpha
        gaps = [-x for x in gaps]
    return a, b, alpha, gaps


def exact_instance_report(ck, label, gaps, T):
    """Common checks on an exact gap vector: nonnegative, contact set, normalisation."""
    contact = [k for k in range(24) if gaps[k] == 0]
    ck.check("%s: all 24 gaps >= 0 (exact) and the contact set contains T; contact set = %s"
             % (label, contact), all(x >= 0 for x in gaps) and set(T) <= set(contact))
    return contact


def fl(q):
    """Short float rendering of an exact rational, for log lines only."""
    return "%.6g" % float(q)


def verify_contact_cone(ck, T):
    """Exact, independent of the facet classification: the cone
        {g in span(minor value vectors) = im P : g = 0 on T, g >= 0 on all 24 permutations}
    has as extreme rays exactly the value vectors of the f_j vanishing on T.  Exhaustive
    over all (d-1)-subsets of the 18 inequality constraints (python-flint integer nullspaces)."""
    import flint
    T = sorted(T)
    mins = [FVAL[j] for j in range(82) if IS_MINOR[j]]
    basis = rref(mins)[0][:14]
    cons = [[b[s] for b in basis] for s in T]
    W = [[sum(co[m] * basis[m][s] for m in range(14)) for s in range(24)] for co in nullspace(cons)]
    d = len(W)
    others = [s for s in range(24) if s not in T]
    rows = [primitive([w[s] for w in W]) for s in others]
    found = set()
    for S in itertools.combinations(range(len(rows)), d - 1):
        X, nul = flint.fmpz_mat([list(rows[i]) for i in S]).nullspace()
        if nul != 1:
            continue
        c = [int(X[i, 0]) for i in range(d)]
        vals = [sum(a * b for a, b in zip(c, r)) for r in rows]
        if all(x <= 0 for x in vals):
            c = [-x for x in c]
        elif not all(x >= 0 for x in vals):
            continue
        found.add(primitive([sum(Fr(c[m]) * W[m][s] for m in range(d)) for s in range(24)]))
    gens = generators(T)
    ck.check("N_T directly: {g in im P : g|T = 0} has dimension %d and the cone {g >= 0} in it has exactly %d extreme rays, "
             "the value vectors of f_j, j in %s (own exhaustive enumeration; Lemma 4.2 for this T without Prop. 4.1)"
             % (d, len(found), gens), found == {tuple(FVAL[j]) for j in gens})
    return d
