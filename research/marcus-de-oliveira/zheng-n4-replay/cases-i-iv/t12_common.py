"""Shared exact objects for the independent replay of the contact cases T_1 and T_2.

Everything here is rebuilt from the DEFINITIONS in Zheng's paper (Sections 3.1, 4.1, 4.2 and the numbering
of Appendix D), by code written for this replay.  Nothing from the supplement is executed or imported.
The audit library (zheng-n4-audit/lib/moc4lib.py) is used ONLY in `crosscheck_with_audit_lib`, to compare
objects; no result below depends on it.

All objects live in the space of "profiles": integer vectors indexed by the 24 permutations of {0,1,2,3}
in lexicographic order of one-line notation (index 0 = identity).  This is the image of the audit's map
c -> Pc, which is injective, so nothing is lost:

  * W        = P(V_R) = linear span of the minor profiles  sigma -> 1{sigma(I) = J}   (dimension 14);
  * f_j      = the j-th of the 82 facet profiles in the numbering of the paper's Appendix D;
  * N_T      = {g in W : g >= 0, g = 0 on T}, the contact cone of a set T of permutations;
  * Phi_g(O) = sum_sigma sgn(sigma) g_sigma prod_i O[i][sigma(i)], the Leibniz sum.

Why W is the span of the minor profiles: prod_i (a_i + b_sigma(i)) = sum_S prod_{i in S} a_i prod_{i not in S}
b_sigma(i), and the coefficient of the monomial a^S b^J (|J| = 4 - |S|) is 1{sigma(S^c) = J}.  The vectors
(prod_i (a_i + b_sigma(i)))_sigma span P(V_C) because the matching basis polynomials are linearly independent.

Python packages: standard library only (fractions, itertools), plus sympy for the polynomial identities.
"""
from fractions import Fraction
from itertools import permutations, combinations
from math import gcd
import os
import sys

# ----------------------------------------------------------------------------------------------------
# permutations (0-based, sigma[i] = image of i), lexicographic order
# ----------------------------------------------------------------------------------------------------
PERMS = [tuple(p) for p in permutations(range(4))]
PIDX = {p: k for k, p in enumerate(PERMS)}
IDENT = (0, 1, 2, 3)


def sgn(p):
    return -1 if sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) % 2 else 1


SGN = [sgn(p) for p in PERMS]


def compose(p, q):
    """(p o q)(i) = p(q(i))."""
    return tuple(p[q[i]] for i in range(4))


def inverse(p):
    r = [0] * 4
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def oneline(p):
    return "".join(str(v) for v in p)


# ----------------------------------------------------------------------------------------------------
# profiles of the 82 facet vectors (paper Section 4.1), 1-based index sets as in the paper
# ----------------------------------------------------------------------------------------------------
def minor_profile(I, J):
    """sigma -> 1{sigma(I) = J};  I, J are 1-based index collections of equal size."""
    I0 = [i - 1 for i in I]
    J0 = set(j - 1 for j in J)
    return tuple(1 if set(s[i] for i in I0) == J0 else 0 for s in PERMS)


def vadd(u, v, c=1):
    return tuple(a + c * b for a, b in zip(u, v))


def t_profile(pi):
    """t_pi = sum_i n_{i, pi(i)};  pi is a 1-based one-line tuple (pi(1),..,pi(4))."""
    out = (0,) * 24
    for i in range(1, 5):
        out = vadd(out, minor_profile((i,), (pi[i - 1],)))
    return out


def d_profile(pi):
    """d_pi = n_{12,{pi1,pi3}} + n_{12,{pi2,pi4}} - n_{13,{pi1,pi2}} - n_{13,{pi3,pi4}}."""
    p1, p2, p3, p4 = pi
    out = minor_profile((1, 2), (p1, p3))
    out = vadd(out, minor_profile((1, 2), (p2, p4)))
    out = vadd(out, minor_profile((1, 3), (p1, p2)), -1)
    out = vadd(out, minor_profile((1, 3), (p3, p4)), -1)
    return out


def h_profile(pi, eps):
    """h_{pi,eps} = (t_pi + eps d_pi)/2; the result is integral (asserted)."""
    s = vadd(t_profile(pi), d_profile(pi), eps)
    assert all(v % 2 == 0 for v in s)
    return tuple(v // 2 for v in s)


# Numbering of Appendix D of the paper (p. 59), typed in for this replay from the page image.
# "nI,J" is a minor vector (digits = elements of I and J), "hPQRS+/-" is h_{pi,eps} with pi = PQRS.
_APPENDIX_D = """
n1,4 n12,34 n13,34 n14,34 n1,3 n12,24 n4,1 n3,1 n12,23 n14,23
n13,24 n13,23 n14,24 n2,1 n1,2 n12,14 n2,4 h3421+ h3412- h2413-
h2431+ n4,2 h3412+ n14,12 h4312+ h3142- h4132- n3,2 h3421- n13,12
h4321- h4123+ h3124+ n12,13 h4321+ h4312- n2,3 h2341+ h2314- n14,13
h4213+ n4,3 h4123- h2413+ h2143- n13,14 h3241- n3,4 h3142+ h2143+
h2341- h1423- h1342+ n13,13 h4231- n3,3 h4132+ h2431- h2134+ n14,14
h3214+ n4,4 h3124- h2134- h2314+ h1432+ h1324- n2,2 h4213- h3241+
h4231+ h3214- n12,12 h1243+ h1432- h1423+ h1234- h1342- h1234+ h1243-
h1324+ n1,1
""".split()
assert len(_APPENDIX_D) == 82


def _parse_spec(s):
    if s[0] == "n":
        a, b = s[1:].split(",")
        return ("n", tuple(int(ch) for ch in a), tuple(int(ch) for ch in b))
    return ("h", tuple(int(ch) for ch in s[1:5]), 1 if s[5] == "+" else -1)


F_SPEC = [_parse_spec(s) for s in _APPENDIX_D]
F_NAME = list(_APPENDIX_D)
F_IS_MINOR = [sp[0] == "n" for sp in F_SPEC]
F_PROF = [minor_profile(sp[1], sp[2]) if sp[0] == "n" else h_profile(sp[1], sp[2]) for sp in F_SPEC]
F_ZERO = [frozenset(k for k in range(24) if f[k] == 0) for f in F_PROF]
PROF_INDEX = {f: j for j, f in enumerate(F_PROF)}

# the 34 minor profiles and the 48 h profiles straight from the definitions (no numbering involved)
ALL_MINOR_PROFILES = set()
for _k in (1, 2):
    for _I in combinations(range(1, 5), _k):
        for _J in combinations(range(1, 5), _k):
            ALL_MINOR_PROFILES.add(minor_profile(_I, _J))
ALL_H_PROFILES = set(h_profile(tuple(v + 1 for v in p), e) for p in PERMS for e in (1, -1))


def generators(T):
    """Indices j (Appendix D numbering) with f_j = 0 on every permutation index of T."""
    T = set(T)
    return [j for j in range(82) if T <= F_ZERO[j]]


# ----------------------------------------------------------------------------------------------------
# exact linear algebra over Q
# ----------------------------------------------------------------------------------------------------
def rref(rows):
    A = [[Fraction(x) for x in r] for r in rows]
    piv, r = [], 0
    if not A:
        return A, piv
    nc = len(A[0])
    for c in range(nc):
        k = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    return A, piv


def rank(rows):
    return len(rref(rows)[1])


def nullspace(rows, ncols):
    """Basis of {x in Q^ncols : rows . x = 0}."""
    if not rows:
        return [[Fraction(int(i == j)) for j in range(ncols)] for i in range(ncols)]
    R, piv = rref(rows)
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for f in free:
        x = [Fraction(0)] * ncols
        x[f] = Fraction(1)
        for r, p in enumerate(piv):
            x[p] = -R[r][f]
        out.append(x)
    return out


def primitive(v):
    """Scale a rational vector to a primitive integer vector (same direction)."""
    den = 1
    for x in v:
        den = den * Fraction(x).denominator // gcd(den, Fraction(x).denominator)
    w = [int(Fraction(x) * den) for x in v]
    g = 0
    for x in w:
        g = gcd(g, abs(x))
    return tuple(x // g for x in w) if g else tuple(w)


# W = span of all minor profiles; a basis as rows
_W_ROWS, _W_PIV = rref([list(p) for p in sorted(ALL_MINOR_PROFILES)] + [[1] * 24])
W_BASIS = [r for r in _W_ROWS if any(x != 0 for x in r)]
W_DIM = len(W_BASIS)


def in_W(g):
    return rank(W_BASIS + [list(g)]) == W_DIM


# ----------------------------------------------------------------------------------------------------
# contact cone N_T by direct enumeration of extreme rays (no use of the 82-facet theorem)
# ----------------------------------------------------------------------------------------------------
def contact_cone_extreme_rays(T, extra_eqs=()):
    """Extreme rays of N_T = {g in W : g >= 0, g|_T = 0}, as primitive integer profiles.  With `extra_eqs`
    (a list of 24-vectors e) the cone is cut further by the equations e . g = 0.

    Method: N_T is a pointed cone in the d-dimensional space W_T = {g in W : g|_T = 0}.  A nonzero g in N_T
    spans an extreme ray iff the inequalities g_sigma >= 0 that are active at g have rank d - 1 on W_T.
    So it suffices to run through every linearly independent (d-1)-subset of the inequalities, take the
    line they cut out, and keep it if it meets the cone.  The search is a depth-first enumeration that
    carries the kernel of the rows chosen so far; rows that do not lower the dimension are skipped, which
    loses nothing because only independent subsets are needed.  Returns (d, sorted list of rays)."""
    T = sorted(set(T))
    # coordinates y in Q^14 on W: g = sum_k y_k W_BASIS[k]
    eqs = [[W_BASIS[k][s] for k in range(W_DIM)] for s in T]
    eqs += [[sum(e[s] * W_BASIS[k][s] for s in range(24)) for k in range(W_DIM)] for e in extra_eqs]
    K = nullspace(eqs, W_DIM)                       # basis of W_T in y-coordinates
    d = len(K)
    G = [[sum(K[a][k] * W_BASIS[k][s] for k in range(W_DIM)) for s in range(24)] for a in range(d)]
    assert all(G[a][s] == 0 for a in range(d) for s in T)
    others = [s for s in range(24) if s not in T]
    rows = [primitive([G[a][s] for a in range(d)]) for s in others]     # inequality rows on Q^d
    rows = [r for r in rows if any(r)]             # a zero row is no constraint
    rows = sorted(set(rows))
    m = len(rows)
    assert rank([list(r) for r in rows]) == d, "cone not pointed"
    rays = set()

    def rec(start, N):
        if len(N) == 1:
            y = N[0]
            vals = [sum(a * b for a, b in zip(r, y)) for r in rows]
            if all(v >= 0 for v in vals):
                pass
            elif all(v <= 0 for v in vals):
                y = [-a for a in y]
            else:
                return
            g = [sum(y[a] * G[a][s] for a in range(d)) for s in range(24)]
            rays.add(primitive(g))
            return
        need = len(N) - 1
        for i in range(start, m - need + 1):
            dots = [sum(a * b for a, b in zip(rows[i], v)) for v in N]
            p = next((j for j, x in enumerate(dots) if x != 0), None)
            if p is None:
                continue
            newN = []
            for j in range(len(N)):
                if j == p:
                    continue
                v = [dots[p] * a - dots[j] * b for a, b in zip(N[j], N[p])]
                gg = 0
                for x in v:
                    gg = gcd(gg, abs(x))
                newN.append([x // gg for x in v] if gg > 1 else v)
            rec(i + 1, newN)

    ident = [[int(i == j) for j in range(d)] for i in range(d)]
    rec(0, ident)
    return d, sorted(rays)


def exact_cone_decomposition(g, gen_profiles):
    """Exact nonnegative coefficients x with sum_k x_k gen_profiles[k] = g, or None.  Tries every linearly
    independent subset of the generators of full rank (a basic feasible solution exists if any solution does)."""
    n = len(gen_profiles)
    r = rank([list(f) for f in gen_profiles])
    for S in combinations(range(n), r):
        cols = [gen_profiles[k] for k in S]
        if rank([list(f) for f in cols]) < r:
            continue
        aug = [[Fraction(cols[a][s]) for a in range(r)] + [Fraction(g[s])] for s in range(24)]
        R, piv = rref(aug)
        if r in piv or piv != list(range(r)):
            continue                                   # inconsistent
        x = [R[a][r] for a in range(r)]
        if all(v >= 0 for v in x):
            out = [Fraction(0)] * n
            for a, k in enumerate(S):
                out[k] = x[a]
            assert all(sum(out[k] * gen_profiles[k][s] for k in range(n)) == g[s] for s in range(24))
            return out
    return None


# ----------------------------------------------------------------------------------------------------
# the symmetry group of order 1152: sigma -> rho o sigma o rho'  or  rho o sigma^{-1} o rho'
# ----------------------------------------------------------------------------------------------------
G_ELEMENTS = [(r1, r2, e) for e in (0, 1) for r1 in PERMS for r2 in PERMS]


def g_perm(g, s):
    r1, r2, e = g
    s = inverse(s) if e else s
    return compose(r1, compose(s, r2))


G_MAPS = [tuple(PIDX[g_perm(g, s)] for s in PERMS) for g in G_ELEMENTS]


def g_image(k, T):
    return frozenset(G_MAPS[k][j] for j in T)


def orbit_of(T):
    return set(g_image(k, T) for k in range(len(G_ELEMENTS)))


def orbit_min(T):
    return min(tuple(sorted(G_MAPS[k][j] for j in T)) for k in range(len(G_ELEMENTS)))


def g_profile(k, g):
    """Push a profile forward: (g.f)_{g.sigma} = f_sigma."""
    out = [0] * 24
    for j in range(24):
        out[G_MAPS[k][j]] = g[j]
    return tuple(out)


def entry_face(i, j):
    """Indices of the six permutations with sigma(i) = j (0-based)."""
    return frozenset(k for k, s in enumerate(PERMS) if s[i] == j)


# ----------------------------------------------------------------------------------------------------
# spectral side, exact over any commutative ring:  F_sigma = lam * prod_i (A_i - C_sigma(i))
# ----------------------------------------------------------------------------------------------------
def F_products(A, C, lam=1):
    out = []
    for s in PERMS:
        t = lam
        for i in range(4):
            t = t * (A[i] - C[s[i]])
        out.append(t)
    return out


def classify_factors(expr, A, C, lam):
    """(sympy)  Factor a polynomial into lam, cross factors A_i - C_j, and within-spectrum differences.
    Returns (cross, withinA, withinC, has_lam, leftover) with cross a set of pairs (i, j)."""
    import sympy as sp
    coeff, facs = sp.factor_list(sp.expand(expr))
    cross, wA, wC, has_lam, left = set(), set(), set(), False, []
    for f, mult in facs:
        hit = False
        if sp.expand(f - lam) == 0 or sp.expand(f + lam) == 0:
            has_lam, hit = True, True
        for i in range(4):
            for j in range(4):
                if sp.expand(f - (A[i] - C[j])) == 0 or sp.expand(f + (A[i] - C[j])) == 0:
                    cross.add((i, j))
                    hit = True
                if i < j and (sp.expand(f - (A[i] - A[j])) == 0 or sp.expand(f + (A[i] - A[j])) == 0):
                    wA.add((i, j))
                    hit = True
                if i < j and (sp.expand(f - (C[i] - C[j])) == 0 or sp.expand(f + (C[i] - C[j])) == 0):
                    wC.add((i, j))
                    hit = True
        if not hit:
            left.append(f)
    return cross, wA, wC, has_lam, left, coeff


# ----------------------------------------------------------------------------------------------------
# check harness
# ----------------------------------------------------------------------------------------------------
class Checker:
    def __init__(self, title):
        self.n = 0
        self.failed = []
        print("=" * 100)
        print(title)
        print("=" * 100)

    def section(self, name):
        print()
        print("--- " + name)

    def check(self, name, ok, detail=""):
        self.n += 1
        tag = "PASS" if ok else "FAIL"
        if not ok:
            self.failed.append(name)
        print("[%s] %s%s" % (tag, name, (" :: " + detail) if detail else ""))
        return ok

    def info(self, text):
        print("       " + text)

    def finish(self):
        print()
        if self.failed:
            print("FAILED checks:")
            for f in self.failed:
                print("   " + f)
            print("FINAL: FAIL (%d of %d checks failed)" % (len(self.failed), self.n))
            return 1
        print("FINAL: PASS (%d checks)" % self.n)
        return 0


# ----------------------------------------------------------------------------------------------------
# checks shared by both cases
# ----------------------------------------------------------------------------------------------------
def check_basic_objects(ck):
    ck.section("Own objects: the 82 facet profiles, the space W, the group")
    ck.check("W = span of minor profiles has dimension 14", W_DIM == 14, "dim = %d" % W_DIM)
    ck.check("34 distinct minor profiles (16 entries + 18 two-by-two)", len(ALL_MINOR_PROFILES) == 34)
    ck.check("48 distinct h profiles, none of them a minor profile",
             len(ALL_H_PROFILES) == 48 and not (ALL_H_PROFILES & ALL_MINOR_PROFILES))
    ck.check("Appendix D table: 82 distinct profiles = the 34 minor profiles and the 48 h profiles",
             len(set(F_PROF)) == 82 and set(F_PROF) == ALL_MINOR_PROFILES | ALL_H_PROFILES)
    ck.check("every f_j is nonnegative and lies in W",
             all(min(f) >= 0 for f in F_PROF) and all(in_W(f) for f in F_PROF))
    ok = True
    for j in range(82):
        cnt = (F_PROF[j].count(0), F_PROF[j].count(1), F_PROF[j].count(2))
        if F_IS_MINOR[j]:
            ok = ok and cnt in ((18, 6, 0), (20, 4, 0))
        else:
            ok = ok and cnt == (13, 10, 1)
    ck.check("value counts: minors 18 or 20 zeros, each h has 13 zeros, ten 1s and one 2", ok)
    ck.check("group: 1152 distinct maps of the 24 permutation indices", len(set(G_MAPS)) == 1152)
    okg = all(g_profile(k, f) in PROF_INDEX for k in range(1152) for f in F_PROF)
    okm = all((g_profile(k, f) in ALL_MINOR_PROFILES) == (f in ALL_MINOR_PROFILES)
              for k in range(1152) for f in F_PROF)
    ck.check("every group element permutes the 82 profiles and maps minor profiles to minor profiles "
             "(all 1152 elements, all 82 profiles)", okg and okm)


def check_lemma_6_2(ck):
    """Phi_{n_{I,J}}(O) = det(O[I,J])^2 on SO(4): three polynomial identities in 16 free entries, then a
    two-line deduction.  (L) Laplace: sum_{sigma(I)=J} sgn prod M[i,sigma(i)] = eps det M[I,J] det M[I^c,J^c];
    (A) M cof(M)^t = det(M) Id;  (J) Jacobi: det cof(M)[I,J] = eps det(M) det M[I^c,J^c] for |I| = 2.
    On SO(4): det = 1 and M^{-1} = M^t, so (A) gives cof(O) = O, (J) gives det O[I,J] = eps det O[I^c,J^c],
    and (L) gives Phi = det O[I,J]^2 (for |I| = 1: Phi = O_ij cof_ij(O) = O_ij^2)."""
    import sympy as sp
    ck.section("Lemma 6.2 (minor vectors have Leibniz sum det(O[I,J])^2 on SO(4)): exact polynomial identities")
    M = sp.Matrix(4, 4, lambda i, j: sp.Symbol("m%d%d" % (i, j)))

    def minor(rows, cols):
        return M.extract(list(rows), list(cols)).det(method="berkowitz")

    detM = sp.expand(M.det(method="berkowitz"))
    cof = sp.Matrix(4, 4, lambda i, j: (-1) ** (i + j) * minor([r for r in range(4) if r != i],
                                                                [c for c in range(4) if c != j]))
    okA = all(sp.expand((M * cof.T)[i, j] - (detM if i == j else 0)) == 0 for i in range(4) for j in range(4))
    ck.check("(A) M cof(M)^t = det(M) Id as a polynomial identity", okA)
    okL, okJ, nL, nJ = True, True, 0, 0
    for k in (1, 2):
        for I in combinations(range(4), k):
            for J in combinations(range(4), k):
                Ic = [r for r in range(4) if r not in I]
                Jc = [c for c in range(4) if c not in J]
                eps = (-1) ** (sum(I) + sum(J))
                lhs = 0
                for s, sg in zip(PERMS, SGN):
                    if set(s[i] for i in I) == set(J):
                        lhs += sg * M[0, s[0]] * M[1, s[1]] * M[2, s[2]] * M[3, s[3]]
                okL = okL and sp.expand(lhs - eps * minor(I, J) * minor(Ic, Jc)) == 0
                nL += 1
                if k == 2:
                    cm = cof.extract(list(I), list(J)).det(method="berkowitz")
                    okJ = okJ and sp.expand(cm - eps * detM * minor(Ic, Jc)) == 0
                    nJ += 1
    ck.check("(L) Laplace form of the Leibniz sum of n_{I,J}, all %d pairs (I,J) with |I| = 1, 2" % nL, okL)
    ck.check("(J) Jacobi complementary-minor identity for the cofactor matrix, all %d pairs with |I| = 2" % nJ,
             okJ)
    # exact spot check on rational SO(4) matrices (Cayley transforms), all 34 minor vectors
    import random
    rng = random.Random(20261125)
    okS = True
    for _ in range(6):
        S = sp.zeros(4, 4)
        for i in range(4):
            for j in range(i + 1, 4):
                v = sp.Rational(rng.randint(-5, 5), rng.randint(1, 4))
                S[i, j], S[j, i] = v, -v
        O = (sp.eye(4) - S) * (sp.eye(4) + S).inv()
        okS = okS and (O * O.T == sp.eye(4)) and O.det() == 1
        w = [sg * O[0, s[0]] * O[1, s[1]] * O[2, s[2]] * O[3, s[3]] for s, sg in zip(PERMS, SGN)]
        for k in (1, 2):
            for I in combinations(range(4), k):
                for J in combinations(range(4), k):
                    prof = minor_profile([i + 1 for i in I], [j + 1 for j in J])
                    phi = sum(a * b for a, b in zip(prof, w))
                    okS = okS and phi == O.extract(list(I), list(J)).det() ** 2
    ck.check("spot check: Phi_{n_{I,J}}(O) = det(O[I,J])^2 on 6 exact rational SO(4) matrices, all 52 (I,J)", okS)


def check_entry_faces(ck):
    """The contact cone of each entry face {sigma : sigma(i) = j} is generated by 15 minor profiles.
    Direct extreme-ray enumeration is done for the face sigma(0) = 0; the other fifteen faces are images of
    it under the group, which acts on W by coordinate permutations (checked), so their cones are images."""
    ck.section("Entry faces {sigma(i) = j}: contact cone generated by minor vectors only")
    ok_gen, ok_cnt = True, True
    for i in range(4):
        for j in range(4):
            gens = generators(entry_face(i, j))
            ok_gen = ok_gen and len(gens) == 15 and all(F_IS_MINOR[k] for k in gens)
            ok_cnt = ok_cnt and sorted(len(F_SPEC[k][1]) for k in gens) == [1] * 6 + [2] * 9
    ck.check("all 16 entry faces: exactly 15 of the f_j vanish on the face, all of them minor vectors", ok_gen)
    ck.check("all 16 entry faces: the 15 are 6 entry vectors and 9 two-by-two vectors", ok_cnt)
    E0 = entry_face(0, 0)
    d, rays0 = contact_cone_extreme_rays(E0)
    ck.check("face sigma(0)=0 = indices {0,..,5}: direct extreme-ray enumeration of N_E (dimension %d) returns "
             "exactly the 15 minor profiles f_0..f_14 (no use of the 82-facet theorem)" % d,
             sorted(E0) == [0, 1, 2, 3, 4, 5] and generators(E0) == list(range(15))
             and set(rays0) == set(F_PROF[k] for k in range(15)))
    ok_tr = True
    for i in range(4):
        for j in range(4):
            E = entry_face(i, j)
            k = next(k for k in range(len(G_ELEMENTS)) if g_image(k, E0) == E)
            keepsW = all(g_profile(k, f) in ALL_MINOR_PROFILES for f in ALL_MINOR_PROFILES)
            img = set(g_profile(k, r) for r in rays0)
            ok_tr = ok_tr and keepsW and img == set(F_PROF[m] for m in generators(E))
    ck.check("the other entry faces: a group element carries {sigma(0)=0} to the face, preserves W, and carries "
             "the 15 extreme rays onto the 15 minor generators of that face", ok_tr)
    ck.info("generators of the face sigma(0)=0: %s"
            % ", ".join("f%d=%s" % (k, F_NAME[k]) for k in generators(E0)))


def crosscheck_with_audit_lib(ck, cases=(1, 2)):
    """Compare own objects with the audit library.  The library is not used anywhere else."""
    ck.section("Cross-check of own objects against the audit library (zheng-n4-audit/lib/moc4lib.py)")
    here = os.path.dirname(os.path.abspath(__file__))
    libdir = os.path.normpath(os.path.join(here, "..", "..", "zheng-n4-audit", "lib"))
    if not os.path.isfile(os.path.join(libdir, "moc4lib.py")):
        ck.check("audit library present", False, "not found at " + libdir)
        return None
    sys.path.insert(0, libdir)
    try:
        import moc4lib as L
    finally:
        sys.path.remove(libdir)
    ck.check("permutation order agrees", [tuple(p) for p in L.PERMS] == PERMS)
    ck.check("all 82 profiles agree with the library's F_VALUES (same numbering)",
             [tuple(v) for v in L.F_VALUES] == F_PROF)
    ck.check("minor flags agree", list(L.F_IS_MINOR) == F_IS_MINOR)
    # W equals the column space of the library's P
    cols = [[L.P[s][a] for s in range(24)] for a in range(14)]
    ck.check("W equals the column space of the library's matrix P",
             rank(cols) == 14 and rank(cols + W_BASIS) == 14)
    for nu in cases:
        ck.check("library CONTACT_SETS[%d] = own T_%d" % (nu, nu), set(L.CONTACT_SETS[nu]) == set(PAPER_T[nu]))
        ck.check("library generators(T_%d) = own generators" % nu,
                 L.generators(L.CONTACT_SETS[nu]) == generators(PAPER_T[nu]))
    ck.check("group action on permutations agrees with the library (all 1152 elements)",
             [tuple(m) for m in L._G_MAPS] == G_MAPS)
    return L


# Appendix E of the paper, entries T_1 and T_2 only (typed in from p. 60).
PAPER_T = {1: frozenset({0, 1, 2, 3, 6, 7}), 2: frozenset({0, 1, 2, 3, 6, 8})}


def supplement_dir():
    d = os.environ.get("MOC4_SUPP", "")
    return d if d and os.path.isdir(d) else None


def compare_with_supplement_data(ck, T):
    """Optional comparison with the supplement's JSON tables, read as DATA (never executed).  This is not
    evidence for anything; it only records whether the document's numbering and generator list agree with
    the recomputation.  Skipped when MOC4_SUPP is not set."""
    import json
    ck.section("Comparison with the supplement's data tables (untrusted data, comparison only)")
    base = supplement_dir()
    if base is None:
        ck.info("SKIPPED: environment variable MOC4_SUPP not set or not a directory")
        return
    p1 = os.path.join(base, "proofs", "principal_minor_result_20261022",
                      "root_permutation_polytope_complete_20261022.json")
    p2 = os.path.join(base, "proofs", "manin_transfer_20261021", "root_moc4_six_contact_orbits_20261022.json")
    if not (os.path.isfile(p1) and os.path.isfile(p2)):
        ck.info("SKIPPED: data files not found under MOC4_SUPP")
        return
    with open(p1) as fh:
        poly = json.load(fh)
    with open(p2) as fh:
        inv = json.load(fh)
    ck.check("supplement facet table: 82 value rows equal own f_0..f_81 in the same order",
             [tuple(f["values"]) for f in poly["facets"]] == F_PROF)
    row = [r for r in inv["rows"] if r["representative"] == sorted(T)]
    ck.check("supplement orbit table has a row for %s" % sorted(T), len(row) == 1)
    if len(row) == 1:
        row = row[0]
        gens = generators(T)
        ck.check("its generator list equals own generators(T)", row["source_facets"] == gens,
                 "supplement %s" % row["source_facets"])
        ck.check("its non-minor list equals own non-minor generators",
                 row["bad_source_facets"] == [j for j in gens if not F_IS_MINOR[j]],
                 "supplement %s" % row["bad_source_facets"])
        ck.check("its orbit size equals own orbit size", row["orbit_size"] == len(orbit_of(T)),
                 "own %d" % len(orbit_of(T)))
