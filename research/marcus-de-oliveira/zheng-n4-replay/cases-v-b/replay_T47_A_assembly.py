"""Independent replay, T_47 part A, document 1: the top-level assembly of the paper's case
T_47 = {1,3,5,6,7,8} (group V) of D. Zheng's proof of the Marcus-de Oliveira conjecture, n = 4.

Proof document (read as PDF, never executed): root_moc4_T46_whole_assembly_20261023.pdf, 3 pages
(the document's "T46" is the paper's T_47).

Run:   MOC4_SUPP=/path/to/supplement  python replay_T47_A_assembly.py [--no-numeric]

Exit code 0 iff every check passes.  EXACT checks use integers, rationals and exact polynomial
arithmetic (sympy, python-flint); lines marked NUMERIC are floating-point falsification attempts.
MOC4_SUPP is needed for: the sha256 comparison of the files that were read, the comparison of
stored objects with the ones rebuilt here, and the sixteen cubic features of the sum-of-squares
certificate used for the direct proof of capacity (5) (they exist only in JSON data of the
supplement).  Without MOC4_SUPP these parts are skipped and the script FAILS.
"""
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T47_A_common as K          # noqa: E402

NUMERIC = "--no-numeric" not in sys.argv
ck = K.Checker("Replay of T_47 part A, document 1: root_moc4_T46_whole_assembly_20261023.pdf  (T_47 = {1,3,5,6,7,8})")
K.print_versions(ck)

T = K.PAPER_T[47]
RAYS = [0, 1, 2, 3, 4, 5, 6, 27, 28, 59]
CANON = [27, 1, 59, 5, 2, 28, 6, 4]                   # coefficients (q, p, j, c, e, d, b, a)
CNAMES = ["q", "p", "j", "c", "e", "d", "b", "a"]
EXCL = [18, 22, 0, 20, 19, 10, 9, 14]                 # permutation where the coefficient is read off
cv = K.combo_values

# ----------------------------------------------------------------------------------------
ck.section("Provenance")
have_supp = K.check_file_identity(ck, ["root_moc4_T46_whole_assembly_20261023.pdf",
                                       "root_moc4_T46_full_radial_20261022.json",
                                       "final_referee_T46_closed_budget_faces_20261023.json",
                                       "independent_moc4_T46_preclosure_allowlist_20261023.json",
                                       "final_referee_T46_budget_PSD8_20261022.json",
                                       "final_referee_T46_cubic_reduction_20261022.json"])
ck.check("EXACT", "MOC4_SUPP is set and the six files are the ones that were read", have_supp)

K.check_foundations(ck)
K.check_group_covariance(ck)
K.crosscheck_with_supplement_data(ck)
K.crosscheck_with_audit_library(ck, [1, 7, 29, 36, 39, 42, 47, 48, 65, 67, 77])

# ----------------------------------------------------------------------------------------
K.check_contact_set(ck, 47, (1, 3, 5, 6, 7, 8), "T46", 65)
ck.check("EXACT", "the six cited contact cases and T_65 lie in orbits different from the orbit of T_47 and from each other",
         len({K.orbit_rep(K.PAPER_T[k]) for k in (1, 7, 29, 36, 39, 42, 47, 65)}) == 8)
gens = K.check_generators(ck, T, RAYS, [28])
K.check_minor_leibniz(ck, [j for j in RAYS if K.F_IS_MINOR[j]])

# ----------------------------------------------------------------------------------------
ck.section("The two linear identities (2) and the absorption into the normal form (3)  [p. 1]")
ck.check("EXACT", "f28 + f0 = f1 + f2 + f5 + f27 and f28 + f3 = f4 + f5 + f27 (all 24 values; P is injective, so these "
         "are identities in V_R)",
         cv([(28, 1), (0, 1)]) == cv([(1, 1), (2, 1), (5, 1), (27, 1)])
         and cv([(28, 1), (3, 1)]) == cv([(4, 1), (5, 1), (27, 1)]))
d8, _ = K.extreme_rays_of_contact_cone(T)
ck.check("EXACT", "the eight vectors f27, f1, f59, f5, f2, f28, f6, f4 are linearly independent and N_T spans a space of "
         "dimension 8: every vector of the span L of N_T has unique coordinates (q,p,j,c,e,d,b,a)",
         K.rank([list(K.F_VALUES[j]) for j in CANON]) == 8 and d8 == 8)
vm = [[K.F_VALUES[j][s] for j in CANON] for s in EXCL]
want = [[int(i == k) for k in range(8)] for i in range(8)]
want[5][3] = 1                                             # at permutation 10 the value is c + d
ck.check("EXACT", "values of the eight vectors at the permutations (18,22,0,20,19,10,9,14): identity matrix except that "
         "the value at 10 is c + d; so q, p, j, c, e, b, a are the gaps at 18, 22, 0, 20, 19, 9, 14 (hence >= 0 on N_T) "
         "and d = gap(10) - gap(20)", vm == want)
import sympy as sp
lam = dict(zip(RAYS, sp.symbols("l0 l1 l2 l3 l4 l5 l6 l27 l28 l59")))
su, sv = sp.symbols("u v")
old = [sum(lam[j] * K.F_VALUES[j][s] for j in RAYS) for s in range(24)]
newc = {0: lam[0] - su, 3: lam[3] - sv, 1: lam[1] + su, 2: lam[2] + su, 4: lam[4] + sv, 5: lam[5] + su + sv,
        6: lam[6], 27: lam[27] + su + sv, 28: lam[28] - su - sv, 59: lam[59]}
new = [sum(newc[j] * K.F_VALUES[j][s] for j in RAYS) for s in range(24)]
ck.check("EXACT", "absorption identity, for all real l_j, u, v: sum l_j f_j = the same sum with f0, f3, f28 reduced by "
         "u, v, u+v and q, c raised by u+v, p, e by u, a by v (24 linear identities; the whole cost is rewritten, nothing "
         "is discarded)", all(sp.expand(x - y) == 0 for x, y in zip(old, new)))
rng = random.Random(4701)
ok = True
n_pos = n_min = 0
for _ in range(4000):
    lv = {j: rng.randint(0, 6) for j in RAYS}
    g = cv(list(lv.items()))
    co = [g[18], g[22], g[0], g[20], g[19], g[10] - g[20], g[9], g[14]]
    ok &= cv(list(zip(CANON, co))) == g
    dd = lv[28] - lv[0] - lv[3]
    ok &= co[5] == dd
    if dd > 0:
        n_pos += 1
        ok &= co == [lv[27] + lv[0] + lv[3], lv[1] + lv[0], lv[59], lv[5] + lv[0] + lv[3], lv[2] + lv[0], dd,
                     lv[6], lv[4] + lv[3]] and min(co) >= 0
    else:
        n_min += 1
        u_ = min(lv[0], lv[28])
        v_ = lv[28] - u_
        mc = {0: lv[0] - u_, 3: lv[3] - v_, 1: lv[1] + u_, 2: lv[2] + u_, 4: lv[4] + v_, 5: lv[5] + u_ + v_,
              6: lv[6], 27: lv[27] + u_ + v_, 59: lv[59]}
        ok &= min(mc.values()) >= 0 and cv(list(mc.items())) == g
ck.check("EXACT", "dichotomy on 4000 random nonnegative integer combinations of the ten rays (seed 4701): the coordinate d "
         "equals l28 - l0 - l3; if it is > 0 the eight coordinates are the nonnegative numbers of the absorption rule "
         "(%d samples); if it is <= 0 the cost equals an explicit nonnegative combination of the nine minor rays "
         "(%d samples)" % (n_pos, n_min), ok and n_pos > 500 and n_min > 500)
ck.info("PROOF-READ: N_T is a pointed cone, hence the hull of its extreme rays, which are the ten vectors (direct "
        "enumeration above).  So every c in N_T is sum l_j f_j with l_j >= 0; with u <= l0, v <= l3, u + v = min(l28, l0 + l3) "
        "the identity gives either a combination of minor vectors (if l0 + l3 >= l28) or the form (3) with d = l28 - l0 - l3 > 0 "
        "and all coefficients >= 0.  Equivalently: for c in N_T, d := gap(10) - gap(20); if d > 0, c has the form (3).")

# ----------------------------------------------------------------------------------------
ck.section("Closed faces where a coefficient vanishes (4)  [pp. 1-2]")
FACES = [("q", 29, "T28"), ("p", 36, "T35"), ("j", 1, "T0"), ("c", 39, "T38"), ("e", 42, "T41"), ("b", 7, "T6")]
stored = {}
registry = {}
if have_supp:
    with open(os.path.join(K.supp_dir(), "final_referee_T46_closed_budget_faces_20261023.json")) as fh:
        stored = {c["budget_zero"]: c for c in json.load(fh)["checks"]}
    with open(os.path.join(K.supp_dir(), "independent_moc4_T46_preclosure_allowlist_20261023.json")) as fh:
        registry = json.load(fh)["registry"]["families"]
G_SET = set(K.G_MAPS)
for nm, nu, label in FACES:
    jz = CANON[CNAMES.index(nm)]
    rest = [j for j in CANON if j != jz]
    Z = [s for s in range(24) if all(K.F_VALUES[j][s] == 0 for j in rest)]
    extra = [s for s in Z if s not in T]
    Tk = K.PAPER_T[nu]
    hits = [gm for gm in K.G_MAPS if set(gm[s] for s in Tk) <= set(Z)]
    ok = bool(hits) and extra == [EXCL[CNAMES.index(nm)]]
    detail = ""
    if hits:
        gm = hits[0]
        g = K.G_ELEMENTS[K.G_MAPS.index(gm)]
        img = K.g_set(gm, Tk)
        ginv = [0] * 24
        for s in range(24):
            ginv[gm[s]] = s
        back = [K.g_facet(tuple(ginv), j) for j in rest]            # g^{-1}.f_j
        gk = K.generators(Tk)
        ok &= all(x in gk for x in back) and K.orbit_rep(img) == K.orbit_rep(Tk)
        detail = ("zero set of the face = T + {%d}; g = (rho=%s, rho'=%s, transpose=%d) maps T_%d onto %s; %d group "
                  "elements; g^-1 maps the seven remaining generators to %s, all generators of T_%d"
                  % (extra[0], g[0], g[1], g[2], nu, list(img), len(hits), back, nu))
    ck.check("EXACT", "face %s = 0 (all other coefficients arbitrary >= 0, a included): the cost vanishes on a seven-element "
             "set that contains a G-image of the paper's T_%d (document label %s); so c(x) lies in g.N_{T_%d} and Case T_%d "
             "applies to g^-1.x" % (nm, nu, label, nu, nu), ok, detail)
    if have_supp:
        c = stored[nm]
        a_map = tuple(c["all24_action"])
        conv = K.g_index_map((tuple(c["row"]), tuple(c["column"]), 1 if c["transpose"] else 0))
        ck.check("EXACT", "   comparison with the supplement's data for this face: its stored index map is the map "
                 "sigma -> row o sigma^(+-1) o column, it is one of the group elements found here, its family label %s names "
                 "the set %s = the paper's T_%d, and its zero set equals mine"
                 % (label, list(Tk), nu),
                 a_map == conv and a_map in G_SET and a_map in hits and c["family"] == label
                 and tuple(registry.get(label, ())) == Tk and c["actual_zero_set"] == Z)
ck.info("the six cited cases: T_29, T_36, T_1, T_39, T_42, T_7 (documents' labels T28, T35, T0, T38, T41, T6: the label is "
        "the paper's number minus one).  All six are in stage 1 of the paper's Table 6; T_47 is in stage 3.  "
        "T_1, T_7, T_29, T_36 were replayed by this owner earlier; T_39 and T_42 belong to the present batch.")
Z14 = [s for s in range(24) if all(K.F_VALUES[j][s] == 0 for j in CANON if j != 4)]
d7, rays7 = K.extreme_rays_of_contact_cone(tuple(T) + (14,))
r7 = sorted(K.F_INDEX[r] for r in rays7)
ck.check("EXACT", "face a = 0: the cost vanishes exactly on T + {14}; N_{T+{14}} has dimension 7 and its extreme rays "
         "(direct enumeration) are f0, f1, f2, f5, f6, f27, f28, f59; with the first identity (2) every element of it is a "
         "combination of minor vectors or has the form (3) with a = 0",
         Z14 == sorted(tuple(T) + (14,)) and d7 == 7 and r7 == [0, 1, 2, 5, 6, 27, 28, 59])
known = {K.orbit_rep(v): k for k, v in K.PAPER_T.items()}
subs14 = [(sub, known.get(K.orbit_rep(sub))) for sub in itertools.combinations(sorted(tuple(T) + (14,)), 6)]
lab14 = sorted(k for _, k in subs14 if k is not None)
ck.check("EXACT", "the seven six-element subsets of T + {14} lie in the orbits of T_47 (two of them), T_48 (two), T_67 (two) "
         "and T_77 (one); none lies in the orbit of one of the six cited cases.  T_48, T_67, T_77 are in stages 8, 8, 9 of "
         "the paper's Table 6, after T_47 (stage 3): the face a = 0 cannot be reduced to an earlier six-contact case in "
         "this way; it is the seven-contact statement of part C",
         lab14 == [47, 47, 48, 48, 67, 67, 77], "%s" % [(list(sb), k) for sb, k in subs14])
ck.info("CITED-INTERNAL (part C, final_referee_T46_a0_uniform): if x in X and c(x) = q f27 + p f1 + j f59 + c f5 + e f2 + "
        "d f28 + b f6 with d > 0 and q, p, j, c, e, b > 0 (a = 0), then Phi_{c(x)} >= 0 on SO(4).  (The document cites it for "
        "the whole cone N_{T+{14}}; the assembly needs only the positive part, the other zero faces being covered above.)")
ck.info("d = 0: the cost is a nonnegative combination of the minor vectors f27, f1, f59, f5, f2, f6, f4 (Lemma 6.2).")

# ----------------------------------------------------------------------------------------
K.check_radial_necessity(ck)
ck.section("The cubic on the lift, C0 and Q  [p. 2]")
t0 = time.time()
rad = K.Radial(CANON, CNAMES)
R = rad.R
q, p, j, c, e, d, b, a, h = R.gens
byh = {}
for ex, co in R.terms(rad.F):
    byh.setdefault(ex[8], []).append((ex[:8] + (0,), co))
C0 = R.from_terms(byh.get(0, []))
Q = R.from_terms(byh.get(1, []))
Dd = a + d
M = (c + d) * (a + d + e)
L = M + b * Dd
ck.check("EXACT", "C(g + h kappa) = C0 + h Q (degree 1 in h) with C0 = -(q+d) L - j M, M = (c+d)(a+d+e), L = M + b(a+d), "
         "Q = j p - j(c+e+d) - (c+d)(e+d) = (m0 m6 - m1 m5)(g); at d = 1 these are the document's formulas",
         set(byh) == {0, 1} and C0 == -(q + d) * L - j * M and Q == j * p - j * (c + e + d) - (c + d) * (e + d)
         and Q == rad.u[0] * rad.u[6] - rad.u[1] * rad.u[5])
c0t = R.terms(C0)
ck.check("EXACT", "every coefficient of C0 is negative and C0 contains -d^3 (%d terms): C0 <= -d^3 < 0 on the closed orthant "
         "with d > 0; so C0 + h Q = 0 forces Q != 0 and h = -C0/Q" % len(c0t),
         all(co < 0 for _, co in c0t) and dict(c0t).get((0, 0, 0, 0, 0, 3, 0, 0, 0)) == -1)

ck.section("The radial identity and the degree-nine polynomial G9  [p. 2]")
det = rad.det()
import flint
rng = random.Random(4702)
okd = True
for _ in range(6):
    x = [rng.randint(-9, 9) for _ in range(9)]
    Mi = flint.fmpz_mat(14, 14, [int(R.evaluate(rad.H[i][k], x)) for i in range(14) for k in range(14)])
    okd &= int(Mi.det()) == R.evaluate(det, x)
ck.check("EXACT", "det Hess C(g + h kappa) as a polynomial in nine variables: %d terms (Laplace expansion in flint "
         "polynomials, %d minors); equal to the integer determinant (fmpz_mat) at 6 random integer points (seed 4702); "
         "degree %d in h" % (len(det), rad.n_minors, R.degrees(det)[8]), okd and R.degrees(det)[8] <= 4,
         "%.1f s" % (time.time() - t0))
okr, G, S = rad.radial_factor(-C0, Q, "j")
N = rad.cleared(-C0, Q)
ck.check("EXACT", "N := Q^4 det Hess C(g - (C0/Q) kappa) has %d terms and N = -4 Q^2 G^2 as a polynomial identity "
         "(exact division, G = gcd(S, dS/dj) with G^2 = S verified, and the product multiplied out again)" % len(N),
         okr and N == Q * Q * G * G * (-4))
gt = R.terms(G)
from math import gcd
cont = 0
for _, co in gt:
    cont = gcd(cont, abs(co))
ck.check("EXACT", "G has 1556 terms (the document's count), is homogeneous of degree 9, has degree 2 in j, and is "
         "primitive (integer content 1); at d = 1 it is the document's G9 (dehomogenisation keeps the number of terms)",
         len(G) == 1556 and K.total_degrees(R, G) == [9] and R.degrees(G)[2] == 2 and cont == 1
         and len({ex[:5] + ex[6:] for ex, _ in gt}) == 1556)
if have_supp:
    with open(os.path.join(K.supp_dir(), "root_moc4_T46_full_radial_20261022.json")) as fh:
        data = json.load(fh)

    def load(terms):
        return R.from_terms([(tuple(ex) + (0,) * (9 - len(ex)), int(co)) for ex, co in terms])
    facs = {f["degree"]: load(f["terms"]) for f in data["radial_factors"]}
    ck.check("EXACT", "supplement data root_moc4_T46_full_radial (JSON, comparison only): same coefficient names, its C and Q "
             "equal C0 and Q, its degree-nine factor equals G with the same sign",
             data["variables"][:8] == CNAMES and load(data["C_terms"]) == C0 and load(data["Q_terms"]) == Q
             and facs.get(9) == G)
ck.info("PROOF-READ: if x in X and c(x) = g with d > 0, then u = Re x = g + h kappa (h = Re F_pi*(x)), C(u) = 0 and "
        "det Hess C(u) = 0 (7.3); so Q(g) != 0, h = -C0/Q, N(g) = 0 and therefore G(g) = 0.  No converse is used.")

# ----------------------------------------------------------------------------------------
ck.section("The second capacity (5): a b (e - 1) >= a + b + e + 1  [p. 2]")
g_four = K.orbit_check(ck, "four-coefficient: (a: f4, b: f6; c = e: f2; d: f28)", (0, 6, 10, 54), (4, 6, 2, 28),
                       sym_positions=(0, 1))
sa, sb, se, sd, ss = sp.symbols("a b e d s", positive=True)
ck.check("EXACT", "the coefficient condition a b (c - d) >= d^2 (a + b + c + d) of PROOF.md 7.6 with c = e, d = 1 is (5); "
         "(5) with a, b, e >= 0 forces e > 1, a, b > 0 and ab > 1: (5) <=> (e-1)(ab-1) >= a + b + 2",
         sp.expand((sa * sb * (se - 1) - (sa + sb + se + 1)) - ((se - 1) * (sa * sb - 1) - (sa + sb + 2))) == 0)
ck.info("CITED (PROOF.md 7.6, four-coefficient inequality): a f0 + b f6 + c f10 + d f54 has Phi >= 0 on SO(4) if a,b,c,d >= 0 "
        "and ab(c-d) >= d^2(a+b+c+d), on its whole G-orbit.  The orbit member (f4, f6; f2; f28) gives (5).")
gf = None
if K.supp_dir() is not None:
    gf = K.GramFace(K.supp_dir(), "final_referee_T46")
    ck.check("EXACT", "certificate data (JSON): source order (27,1,59,5,2,28,6,4); the four signed coordinate maps permute "
             "p and q separately; quaternion model M(p,q): M^t M = |p|^2|q|^2 I and det M = (|p|^2|q|^2)^2",
             gf.order == CANON and gf.signed_perms_ok and gf.model_ok)
    key = {k: i for i, k in enumerate(gf.keys)}

    def unit(i, k):
        v = [Fr(0)] * 36
        v[key[(i, k)]] = Fr(1)
        return v
    wr = gf.weight_rows
    w3 = [x - y for x, y in zip(wr[3], unit(3, 3))]
    w5 = [x - y + z for x, y, z in zip(wr[5], unit(5, 5), unit(3, 3))]
    ck.check("EXACT", "own linear algebra on the identity (*) of GramFace (%d coefficient equations, rank %d): the solutions "
             "(Gram matrix, weights) are: 10 independent linear relations on the Gram matrix, and weights "
             "(G00, G11, G22, G33, G44, G55 - G33, G66, G77) modulo the relations, i.e. Gram diagonal (q,p,j,c,e,c+d,b,a); "
             "the ten stored relations span the same space" % (gf.n_equations, gf.rank),
             len(gf.relations) == 10 and all(wr[m] == unit(m, m) for m in (0, 1, 2, 4, 6, 7))
             and K.rank(gf.relations + [w3]) == 10 and K.rank(gf.relations + [w5]) == 10
             and K.rank(gf.stored_R) == 10 and K.rank(gf.stored_R + gf.relations) == 10)
    emin = 1 + (sa + sb + 2) / (sa * sb - 1)
    K.sym_gram_check(ck, gf, "direct certificate of (5) at d = 1 (no PDF contains it; the supplement has only a program): "
                     "vectors (E, D, B, A) with |D|^2 = 1, |B|^2 = b, |A|^2 = a, B.A = 1, E.A = -1, E.B = 1, E.D = -1, "
                     "|E|^2 = 1 + (a+b+2)/(ab-1) + s",
                     {(4, 5): -1, (4, 6): 1, (4, 7): -1, (6, 7): 1},
                     [0, 0, 0, 0, emin + ss, 1, sb, sa], [0, 0, 0, 0, emin + ss, 1, sb, sa],
                     [5, 6, 7, 4], [1, sb, (sa * sb - 1) / sb, ss],
                     {sa: Fr(5, 2), sb: Fr(2), ss: Fr(1, 7)})
    ck.check("EXACT", "e >= 1 + (a+b+2)/(ab-1) with ab > 1 is exactly (5)",
             sp.cancel((se - emin) * (sa * sb - 1) - (sa * sb * (se - 1) - (sa + sb + se + 1))) == 0)
    ck.info("PROOF-READ: for ab > 1, s >= 0 all pivots are >= 0, so the Gram matrix is PSD and (*) gives "
            "Phi_{a f4 + b f6 + e f2 + f28} >= 0 on SO(4) (surjectivity of (p,q) -> M(p,q) onto SO(4): CITED, classical); "
            "the other minor vectors are added with Phi >= 0.  This re-proves (5) on this face without the citation.")
else:
    ck.check("EXACT", "direct certificate of (5) replayed", False, "MOC4_SUPP not available")

# ----------------------------------------------------------------------------------------
ck.section("Exhaustiveness of the three sectors (6)  [p. 2]")


def sectors(qv, bv):
    out = []
    if (qv <= 1 and bv <= 1) or qv == bv:
        out.append(1)
    if qv * bv <= 1 and max(qv, bv) > 1:
        out.append(2)
    if qv * bv > 1 and qv != bv:
        out.append(3)
    return out


grid = [Fr(n, m) for n in range(1, 13) for m in range(1, 7)]
cnt = {}
okc = True
for qv in grid:
    for bv in grid:
        sc = sectors(qv, bv)
        okc &= len(sc) == 1
        cnt[tuple(sc)] = cnt.get(tuple(sc), 0) + 1
ck.check("EXACT", "on the grid of %d rational pairs (q, b) (boundaries qb = 1, q = b, q = 1, b = 1 included) every pair "
         "lies in exactly one of: [q,b <= 1 or q = b], [qb <= 1, max > 1], [qb > 1, q != b]" % (len(grid) ** 2), okc,
         "counts per sector %s" % sorted(cnt.items()))
ck.info("PROOF-READ: for q, b > 0: if qb <= 1 then not both exceed 1, so either both are <= 1 (sector 1) or exactly one "
        "exceeds 1 (sector 2); if qb > 1 then q = b (sector 1) or q != b (sector 3).  The three sets are disjoint.")
ck.info("CITED-INTERNAL (part B), needed for each sector S: for all positive q, p, j, c, e, b, a with d = 1 and (q, b) in S: "
        "if x in X and c(x) = g, then Phi_g >= 0 on SO(4)  [in the notes: a sufficient capacity holds for g, or G9(g) != 0].")

# ----------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence only)")
    Gt, C0tt, Qt = R.terms(G), R.terms(C0), R.terms(Q)
    A8 = np.array([K.F_VALUES[jx] for jx in CANON], float).T
    # (i) the necessary equations on generic points of X satisfying the five contact equations
    pts = K.numeric_equality_points(T, 47005, 60)
    worst = [0.0, 0.0, 0.0, 1e9]
    for pnt in pts:
        co, *_ = np.linalg.lstsq(A8, pnt["gaps"], rcond=None)
        worst[0] = max(worst[0], float(np.max(np.abs(A8 @ co - pnt["gaps"]))))
        c0v, qv, gv = K.float_eval(C0tt, co), K.float_eval(Qt, co), K.float_eval(Gt, co)
        worst[1] = max(worst[1], abs(c0v + pnt["H"] * qv) / (K.float_eval(C0tt, co, True) + abs(pnt["H"]) * K.float_eval(Qt, co, True)))
        worst[2] = max(worst[2], abs(gv) / K.float_eval(Gt, co, True))
        worst[3] = min(worst[3], abs(qv) / K.float_eval(Qt, co, True))
    ck.check("NUMERIC", "%d generic points x = alpha Psi(a,b) whose real parts satisfy the five contact equations (seed "
             "47005; the other gaps have either sign): the gap vector is in the span of the eight generators (fit %.1e), "
             "|C0 + h Q| <= %.1e (relative), |G| <= %.1e (relative to the sum of |monomials|), smallest relative |Q| = %.1e"
             % (len(pts), worst[0], worst[1], worst[2], worst[3]),
             len(pts) >= 40 and worst[0] < 1e-9 and worst[1] < 1e-9 and worst[2] < 1e-9)
    # (ii) stratum points: unsteered, then steered towards d > 0
    t1 = time.time()
    free = K.numeric_stratum_search(T, 47100, 250)
    dfree = [pt["gaps"][10] - pt["gaps"][20] for pt in free]
    ck.check("NUMERIC", "unsteered search for points with all gaps >= 0, zero on T, sum of gaps 1 (seed 47100, 250 starts): "
             "%d points found; the coordinate d ranges over [%.3g, %.3g]: none has d > 1e-6 (they are combinations of "
             "minor vectors)" % (len(free), min(dfree), max(dfree)), len(free) > 50, "%.0f s" % (time.time() - t1))
    found = []
    for dmin, seed in ((0.002, 47200), (0.002, 47201), (0.01, 47200), (0.01, 47201)):
        found += [(dmin, pt) for pt in K.numeric_stratum_search(
            T, seed, 300, extra=lambda g, dmin=dmin: [10 * min(0.0, (g[10] - g[20]) - dmin)])]
    rows = []
    for dmin, pt in found:
        g = pt["gaps"]
        co = np.array([g[18], g[22], g[0], g[20], g[19], g[10] - g[20], g[9], g[14]])
        x = co / co[5]
        relG = abs(K.float_eval(Gt, x)) / K.float_eval(Gt, x, True)
        Hn = pt["H"] / co[5]
        relC = abs(K.float_eval(C0tt, x) + Hn * K.float_eval(Qt, x)) / (K.float_eval(C0tt, x, True) + abs(Hn) * K.float_eval(Qt, x, True))
        qq, pp, jj, cc, ee, _, bb, aa = x
        rows.append(dict(d=co[5], ncont=int(np.sum(np.abs(g) < 1e-9)), sep=pt["sep"], relG=relG, relC=relC,
                         four=aa * bb * (ee - 1) - (aa + bb + ee + 1), minc=float(np.min(x)), x=x,
                         sec=sectors(qq, bb)))
    good = [r for r in rows if r["minc"] > 1e-3 and r["sep"] > 0.05]
    for r in sorted(rows, key=lambda r: -r["sep"])[:12]:
        ck.info("   d = %.3g, contacts %d, label separation %.3f, |G| rel %.1e, |C0+hQ| rel %.1e, (5) margin %.3g, sector %s, "
                "(q,p,j,c,e,b,a)/d = %s" % (r["d"], r["ncont"], r["sep"], r["relG"], r["relC"], r["four"], r["sec"],
                                             " ".join("%.3g" % v for i, v in enumerate(r["x"]) if i != 5)))
    ck.check("NUMERIC", "search steered to d >= 0.002 and d >= 0.01 (seeds 47200, 47201; 300 starts each): %d stratum points "
             "with d > 0, of which %d have exactly the six contacts, all seven other coefficients > 1e-3 d and separated "
             "labels; on all %d points |G| <= %.1e and |C0 + hQ| <= %.1e (relative), and the capacity (5) holds at every "
             "one of them (so the case is not vacuous, and on the points found the conclusion comes through (5))"
             % (len(rows), len(good), len(rows), max(r["relG"] for r in rows), max(r["relC"] for r in rows)),
             len(rows) >= 5 and len(good) >= 2 and max(r["relG"] for r in rows) < 1e-7
             and all(r["four"] > 0 for r in rows), "%.0f s" % (time.time() - t1))
    # (iii) capacity (5) at its threshold
    rng2 = np.random.default_rng(47300)
    worst5 = float("inf")
    ctrl = float("inf")
    for _ in range(10):
        av, bv = np.exp(rng2.normal(size=2))
        if av * bv <= 1.05:
            av = 1.5 / bv + av
        em = 1 + (av + bv + 2) / (av * bv - 1)
        vals = np.array(cv([(4, 1)]), float) * av + np.array(cv([(6, 1)]), float) * bv \
            + np.array(cv([(2, 1)]), float) * em + np.array(cv([(28, 1)]), float)
        worst5 = min(worst5, K.numeric_min_phi(vals, rng2, restarts=8))
        vals2 = vals - 0.3 * (em - 1) * np.array(cv([(2, 1)]), float)
        ctrl = min(ctrl, K.numeric_min_phi(vals2, rng2, restarts=8))
    ck.check("NUMERIC", "capacity (5) at equality for 10 random (a, b) with ab > 1 (seed 47300): smallest Phi found over "
             "SO(4) by multistart BFGS is %.3e (>= -1e-9); control with e lowered by 30%% of (e - 1): %.3g (< 0: the search "
             "detects a violation below the threshold)" % (worst5, ctrl), worst5 > -1e-9 and ctrl < -1e-7)

ck.finish()
