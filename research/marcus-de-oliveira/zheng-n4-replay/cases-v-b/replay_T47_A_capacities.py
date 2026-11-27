"""Independent replay, T_47 part A, documents 2 and 3: the two Gram capacities on the cone of the
paper's T_47 = {1,3,5,6,7,8}.

Proof documents (read as PDF, never executed):
  final_referee_T46_six_ray_coherent_capacity_20261023.pdf   (4 pages; Theorem 1, Corollaries 2-4)
  final_referee_T46_core_coherent_capacity_20261023.pdf      (5 pages; Theorems 1-2, Corollary 3)

Run:   MOC4_SUPP=/path/to/supplement  python replay_T47_A_capacities.py [--no-numeric]

Exit code 0 iff every check passes.  The sum-of-squares certificate behind both notes (sixteen
cubic polynomials, four signed coordinate permutations, an 8 x 8 change of basis) exists only in
JSON data of the supplement; it is read as data and the identity it is claimed to satisfy is
verified here against Leibniz sums rebuilt from the paper's definitions (class GramFace of the
helper).  Without MOC4_SUPP the script FAILS, because nothing can then be replayed.
"""
import json
import os
import sys
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T47_A_common as K          # noqa: E402

NUMERIC = "--no-numeric" not in sys.argv
ck = K.Checker("Replay of T_47 part A, documents 2 and 3: the six-ray and the five-ray Gram capacities")
K.print_versions(ck)
import sympy as sp                        # noqa: E402

CANON = [27, 1, 59, 5, 2, 28, 6, 4]                   # (q, p, j, c, e, d, b, a); Gram vectors (Q,P,J,C,E,D,B,A)
T47, T65 = K.PAPER_T[47], K.PAPER_T[65]

ck.section("Provenance")
have_supp = K.check_file_identity(ck, ["final_referee_T46_six_ray_coherent_capacity_20261023.pdf",
                                       "final_referee_T46_core_coherent_capacity_20261023.pdf",
                                       "final_referee_T46_budget_PSD8_20261022.json",
                                       "final_referee_T46_cubic_reduction_20261022.json",
                                       "final_referee_T46_core_T47_transports_20261023.json",
                                       "independent_moc4_T47_budget_PSD8_20261022.json"])
ck.check("EXACT", "MOC4_SUPP is set and the six files are the ones that were read", have_supp)
if not have_supp:
    ck.finish()
K.check_foundations(ck)
K.check_group_covariance(ck)

# ----------------------------------------------------------------------------------------
ck.section("The sum-of-squares representation on SO(4) (GramFace), rebuilt")
gf = K.GramFace(K.supp_dir(), "final_referee_T46")
ck.check("EXACT", "certificate data (JSON): source order (27,1,59,5,2,28,6,4), the same eight generators and names as the "
         "top-level note; the four signed coordinate maps permute p and q separately (they preserve |p|^2, |q|^2); "
         "quaternion model M(p,q) = matrix of x -> p x q: M^t M = |p|^2|q|^2 I, det M = (|p|^2|q|^2)^2",
         gf.order == CANON and gf.signed_perms_ok and gf.model_ok)
key = {k: i for i, k in enumerate(gf.keys)}


def unit(i, k):
    v = [Fr(0)] * 36
    v[key[(i, k)]] = Fr(1)
    return v


wr = gf.weight_rows
w3 = [x - y for x, y in zip(wr[3], unit(3, 3))]
w5 = [x - y + z for x, y, z in zip(wr[5], unit(5, 5), unit(3, 3))]
ck.check("EXACT", "own linear algebra on (*): (|p|^2+|q|^2) Sigma(Gm) = 4 Phi_w(M(p,q)) mod (|p|^2-|q|^2)  (%d coefficient "
         "equations, rank %d) holds exactly when Gm satisfies 10 independent linear relations and "
         "w = (G00, G11, G22, G33, G44, G55 - G33, G66, G77) modulo them: the Gram diagonal is (q, p, j, c, e, c + d, b, a); "
         "the ten relations stored in the supplement span the same space" % (gf.n_equations, gf.rank),
         len(gf.relations) == 10 and all(wr[m] == unit(m, m) for m in (0, 1, 2, 4, 6, 7))
         and K.rank(gf.relations + [w3]) == 10 and K.rank(gf.relations + [w5]) == 10
         and K.rank(gf.stored_R) == 10 and K.rank(gf.stored_R + gf.relations) == 10)
ck.info("relations derived here (Gik = Gram entries of the vectors (Q,P,J,C,E,D,B,A)):")
for r in gf.relations:
    ck.info("   " + " ".join("%+d*G%d%d" % (int(cf), gf.keys[i][0], gf.keys[i][1]) for i, cf in enumerate(r) if cf) + " = 0")
ck.info("consequence (PROOF-READ + CITED surjectivity of unit (p,q) -> M(p,q) onto SO(4)): on unit quaternions (*) reads "
        "Sigma(Gm) = 2 Phi_w; if a PSD matrix Gm satisfies the ten relations, Phi of the cost with weights w is >= 0 on SO(4)")


def comb(coeffs, ent):
    return sum(sp.Rational(c.numerator, c.denominator) * e for c, e in zip(coeffs, ent) if c != 0)


def gram_ok(Gm, weights):
    """relations vanish and the weights are the claimed ones (symbolic, exact)."""
    ent = [Gm[i, k] for (i, k) in gf.keys]
    rel = all(sp.cancel(comb(r, ent)) == 0 for r in gf.relations)
    w = [comb(gf.weight_rows[m], ent) for m in range(8)]
    return rel and all(sp.cancel(x - y) == 0 for x, y in zip(w, weights))


def direct(Gm, weights, inst):
    """identity (*) multiplied out directly at a rational parameter point."""
    sub = {s: sp.Rational(v.numerator, v.denominator) for s, v in inst.items()}
    Gn = sp.Matrix(Gm).subs(sub).applyfunc(sp.cancel)
    wn = [sp.cancel(sp.sympify(w).subs(sub)) for w in weights]
    den = 1
    for x in list(Gn) + wn:
        x = sp.cancel(x)
        assert x.is_Rational, x
        den = sp.ilcm(den, x.q)
    Gi = [[int(Gn[i, k] * den) for k in range(8)] for i in range(8)]
    return gf.direct_identity(Gi, [int(w * den) for w in wn])


def minors_positive(Mm, inst):
    sub = {s: sp.Rational(v.numerator, v.denominator) for s, v in inst.items()}
    Mn = sp.Matrix(Mm).subs(sub)
    return all(Mn[:k, :k].det() > 0 for k in range(1, Mn.shape[0] + 1))


# ========================================================================================
ck.section("Document 2 (six-ray note), Theorem 1  [pp. 1-2]")
q, b, j, beta, t, v = sp.symbols("q b j beta t v", real=True)
s_ = q + b + 2 * beta
h_ = 1 + beta
ell = q * b - beta ** 2
R0 = b * (j - t * t) - v * v
Delta = s_ * (j - t * t) - (t + v) ** 2
lam = (-h_ * (j - t * t) + (t + v) * (1 + v)) / Delta
mu = (h_ * (t + v) - s_ * (1 + v)) / Delta
metric = sp.Matrix([[q, beta, t, 0], [beta, b, v, 0], [t, v, j - t * t, 0], [0, 0, 0, 1]])     # basis (Q, B, J0, Z)
cfac = (t - v) / R0
VV = sp.Matrix([[1, 0, 0, 0], [lam, lam, mu, 0], [0, 0, 1, t], [0, -v * cfac, b * cfac, 0], [0, 0, 0, 0],
                [0, -v * cfac, b * cfac, -1], [0, 1, 0, 0], [0, 0, 0, 0]])                      # rows Q,P,J,C,E,D,B,A
Gm6 = (VV * metric * VV.T).applyfunc(sp.cancel)
pstar = ((j - t * t) * h_ ** 2 - 2 * h_ * (t + v) * (1 + v) + s_ * (1 + v) ** 2) / Delta
cstar = (t - v) ** 2 / (j - t * t - v * v / b)
W6 = [q, pstar, j, cstar, 0, 1, b, 0]
ck.check("EXACT", "leaf Schur complement: det of the 3 x 3 metric of (Q, B, J0) = ell (j - t^2) - (b t^2 - 2 beta t v + q v^2), "
         "ell = qb - beta^2; so the hypotheses q, b > 0, ell > 0, j - t^2 > (b t^2 - 2 beta t v + q v^2)/ell say exactly "
         "that the 4 x 4 metric of (Q, B, J0, Z) is positive definite; Delta = Gram determinant of (Q + B, J0) and "
         "b(j - t^2) - v^2 = Gram determinant of (B, J0) are then > 0",
         sp.expand(metric[:3, :3].det() - (ell * (j - t * t) - (b * t * t - 2 * beta * t * v + q * v * v))) == 0
         and sp.expand(Delta - sp.Matrix([[s_, t + v], [t + v, j - t * t]]).det()) == 0
         and sp.expand(R0 - sp.Matrix([[b, v], [v, j - t * t]]).det()) == 0
         and sp.expand((sp.Matrix([[1, 1, 0, 0]]) * metric * sp.Matrix([1, 1, 0, 0]))[0] - s_) == 0)
ck.check("EXACT", "the eight vectors (Q, P, J, C, 0, D, B, 0) with P = lambda (Q+B) + mu J0, J = J0 + t Z, "
         "C = (t-v)/(j-t^2-v^2/b) (J0 - (v/b) B), D = C - Z have squared lengths (q, p*, j, c*, 0, c* + 1, b, 0) with the "
         "document's p* and c*; P.(Q+B) = -h, P.J0 = -1-v, C.B = C.Z = 0, C.J = t - v",
         all(sp.cancel(Gm6[i, i] - z) == 0 for i, z in enumerate([q, pstar, j, cstar, 0, cstar + 1, b, 0]))
         and sp.cancel(Gm6[1, 0] + Gm6[1, 6] + h_) == 0
         and sp.cancel((VV[1, :] * metric * sp.Matrix([0, 0, 1, 0]))[0] + 1 + v) == 0
         and sp.cancel(Gm6[3, 6]) == 0 and sp.cancel(Gm6[3, 3] - Gm6[3, 5]) == 0 and sp.cancel(Gm6[3, 2] - (t - v)) == 0)
ck.check("EXACT", "p* = (-h, -1-v) Gram(Q+B, J0)^-1 (-h, -1-v)^t (least squared length with the two prescribed products) "
         "and c* = (t-v)^2 / |J0 - (v/b)B|^2",
         sp.cancel((sp.Matrix([[-h_, -1 - v]]) * sp.Matrix([[s_, t + v], [t + v, j - t * t]]).inv()
                    * sp.Matrix([-h_, -1 - v]))[0] - pstar) == 0
         and sp.cancel(cstar - (t - v) ** 2 * b / R0) == 0)
ck.check("EXACT", "THEOREM 1 (sufficiency): this Gram matrix satisfies the ten relations identically in (q, b, j, beta, t, v) "
         "and its weights are (q, p*, j, c*, 0, 1, b, 0): (*) holds for the cost q f27 + p* f1 + j f59 + c* f5 + f28 + b f6",
         gram_ok(Gm6, W6))
inst6 = {q: Fr(2), b: Fr(3), beta: Fr(1, 2), t: Fr(1, 3), v: Fr(-1, 2), j: Fr(5)}
inst6b = {q: Fr(1, 3), b: Fr(5, 2), beta: Fr(-3, 4), t: Fr(-2, 5), v: Fr(1, 4), j: Fr(3, 2)}
ck.check("EXACT", "identity (*) multiplied out directly (no linear algebra) at two rational parameter points with positive "
         "definite metric: (q,b,beta,t,v,j) = (2,3,1/2,1/3,-1/2,5) and (1/3,5/2,-3/4,-2/5,1/4,3/2)",
         minors_positive(metric, inst6) and minors_positive(metric, inst6b)
         and direct(Gm6, W6, inst6) and direct(Gm6, W6, inst6b))
s1, s2 = sp.symbols("s1 s2", nonnegative=True)
Gsur = sp.Matrix(Gm6)
Gsur[1, 1] += s1
Gsur[3, 3] += s2
Gsur[3, 5] += s2
Gsur[5, 3] += s2
Gsur[5, 5] += s2
ck.check("EXACT", "surplus: adding s1 to |P|^2 and one vector of squared length s2, orthogonal to everything, to both C and "
         "D keeps the ten relations and changes the weights to (q, p* + s1, j, c* + s2, 0, 1, b, 0)",
         gram_ok(Gsur, [q, pstar + s1, j, cstar + s2, 0, 1, b, 0]))
ck.info("PROOF-READ: Gm = V metric V^t is PSD when the metric is; so Phi >= 0 on SO(4) for p >= p*, c >= c*.  Further "
        "nonnegative multiples of the minor vectors f2, f4 (and of any minor vector) may be added.  The set of costs with "
        "Phi >= 0 on SO(4) is closed, so every finite limit of admissible (q, p, j, c, b) is admissible.")
ck.info("REMARK (EXACT): the identity is rational in the parameters, so it holds wherever Delta and b(j-t^2)-v^2 are nonzero; "
        "the Gram matrix is PSD whenever the 3 x 3 metric is PSD.  Corollary 2 uses this with ell = 0.")

ck.section("Document 2, the classification (converse) part of Theorem 1  [p. 2]")
core = [0, 1, 2, 3, 5, 6]                                   # Q, P, J, C, D, B
zz = sp.symbols("z0:21")
ce = [(i, k) for i in range(6) for k in range(i, 6)]
cm = sp.zeros(6, 6)
for (i, k), x in zip(ce, zz):
    cm[i, k] = cm[k, i] = x
emb = sp.zeros(8, 8)
for i, oi in enumerate(core):
    for k, ok_ in enumerate(core):
        emb[oi, ok_] = cm[i, k]
ent = [emb[i, k] for (i, k) in gf.keys]
restricted = sp.Matrix([comb(r, ent) for r in gf.relations]).jacobian(zz)
bad = cm[4, 4] - cm[3, 3]                                   # |D|^2 - |C|^2  (= d)
nine = sp.Matrix([cm[0, 1] + cm[0, 5] + cm[1, 5] + bad,     # QP + QB + PB + d = 0
                  cm[0, 3] - cm[0, 4],                      # Q.Z = 0   (Z = C - D)
                  cm[1, 3] - cm[1, 4],                      # P.Z = 0
                  cm[3, 5] - cm[4, 5],                      # B.Z = 0
                  cm[3, 3] - cm[3, 4],                      # C.Z = 0
                  cm[3, 5],                                 # C.B = 0
                  cm[0, 2] - cm[2, 3] + cm[2, 4],           # J.Q = J.Z
                  cm[1, 2] + cm[2, 5] + bad,                # JP + JB + d = 0
                  cm[2, 5] + cm[2, 4]]).jacobian(zz)        # JC + JB = JZ  (with C.. : J.B + J.D = ... see report)
ck.check("EXACT", "when e = a = 0 (rows and columns E, A of a PSD Gram matrix vanish) the ten relations restricted to the six "
         "vectors (Q, P, J, C, D, B) have rank 9 and the same row space as the nine equations displayed in the document "
         "(written with Z = C - D and d = |D|^2 - |C|^2)",
         restricted.rank() == 9 and nine.rank() == 9 and restricted.col_join(nine).rank() == 9)
tv = sp.Matrix([[1, 0, 0, 0, 0, 0], [sp.Rational(-1, 2), sp.Rational(-1, 2), -1, 0, 1, 0], [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1], [0, 0, 0, -1, 0, 1], [0, 1, 0, 0, 0, 0]])           # rows Q, P, J, C, D, B
tmpl = tv * tv.T
ck.check("EXACT", "the document's template (six explicit vectors) has a positive definite Gram matrix (leading minors %s), "
         "satisfies the nine equations and has |D|^2 - |C|^2 = 1" % [str(tmpl[:k, :k].det()) for k in range(1, 7)],
         all(tmpl[:k, :k].det() > 0 for k in range(1, 7)) and tmpl[4, 4] - tmpl[3, 3] == 1
         and nine * sp.Matrix([tmpl[i, k] for i, k in ce]) == sp.zeros(9, 1))
ck.info("PROOF-READ (not needed for any support inequality): for a positive definite completion, Z = C - D is a unit vector "
        "orthogonal to Q, P, B, C; with beta = Q.B, t = J.Z, v = J.B, J0 = J - tZ the equations give P.(Q+B) = -1-beta, "
        "P.J0 = -1-v, C.B = 0, C.J0 = t-v, so |P|^2 >= p* and |C|^2 >= c* by projection; mixing with the template gives "
        "density of positive definite completions.  That every degree-six certificate has this one-block form is a claim of "
        "the supplement's data ('proof' field of the cubic-reduction file) and was NOT checked here.")

ck.section("Document 2, the transport to the cone of the paper's T_65 (document label T47)  [p. 2]")
img = K.doc_action_facets((1, 2, 0, 3), (2, 1, 3, 0), True, CANON)
g65 = K.generators(T65)
ck.check("EXACT", "the action with row (1,2,0,3), column (2,1,3,0) and transposition, in the supplement's convention "
         "(pullback by sigma -> row o sigma^-1 o column), sends the roles [27,1,59,5,2,28,6,4] to [0,45,5,9,29,35,7,36]",
         img == [0, 45, 5, 9, 29, 35, 7, 36], "computed %s" % img)
ck.check("EXACT", "the six images f0, f45, f5, f9, f35, f7 of the core (q, p, j, c, d, b) are generators of T_65 = "
         "{0,4,5,6,9,12}, whose generators are f0, f1, f2, f5, f7, f9, f35, f45 (only f35 non-minor); the images f29, f36 of "
         "e, a are not: only the six-ray core is transported, f1 and f2 are added as minor vectors",
         g65 == [0, 1, 2, 5, 7, 9, 35, 45] and all(x in g65 for x in (0, 45, 5, 9, 35, 7))
         and 29 not in g65 and 36 not in g65 and [x for x in g65 if not K.F_IS_MINOR[x]] == [35])
ck.info("roles: (q, b, j, p, c) of T46 = coefficients of (f27, f6, f59, f1, f5) correspond to the coefficients of "
        "(f0, f7, f5, f45, f9) on the T_65 cone, which the document calls (b, a, e, c, p)_47; d_47 is the coefficient of f35.")
with open(os.path.join(K.supp_dir(), "independent_moc4_T47_budget_PSD8_20261022.json")) as fh:
    so65 = json.load(fh)["source_order"]
ck.check("EXACT", "comparison with the supplement's data for its 'T47' cone: the coefficient names (q,p,j,c,e,d,b,a)_47 "
         "belong to (f2, f9, f1, f45, f5, f35, f0, f7); so (b, a, e, c, p)_47 are the coefficients of (f0, f7, f5, f45, f9), "
         "as the document says, and q_47, j_47 are those of f2, f1",
         so65 == [2, 9, 1, 45, 5, 35, 0, 7] and sorted(so65) == g65)
with open(os.path.join(K.supp_dir(), "final_referee_T46_core_T47_transports_20261023.json")) as fh:
    tr = json.load(fh)
okw = tr["T46_core"] == [1, 59, 6, 27, 28] and tr["T47_core"] == [45, 5, 7, 0, 35]
for w in tr["witnesses"]:
    okw &= K.doc_action_facets(w["row"], w["column"], w["transpose"], CANON) == w["full_T46_canonical_image"]
    okw &= K.doc_action_facets(w["row"], w["column"], w["transpose"], [1, 59, 6, 27, 28]) == [45, 5, 7, 0, 35]
ck.check("EXACT", "comparison with the supplement's transport data (%d stored witnesses): each stored action, in the same "
         "convention, maps the eight roles as stored and the five-ray core [1,59,6,27,28] to [45,5,7,0,35]"
         % len(tr["witnesses"]), okw)

ck.section("Document 2, Corollary 2 (singular leaf)  [p. 3]")
r, Bq, E, Pp = sp.symbols("r B E P", positive=True)          # T47 roles: B = b_47, r = sqrt(a b), E = e_47
mp = {q: Bq, b: r ** 2 / Bq, j: E, beta: -r, v: -r * t / Bq}
kap = 1 + 1 / Bq
rho = 1 + r / Bq
eta = (1 - r ** 2 / Bq) / (Bq - r)
kk = Bq * (1 - r) ** 2 / (Bq - r) ** 2
rem = E - kap * t * t
msub = metric.subs(mp, simultaneous=True)
schur = (msub[1:3, 1:3] - msub[1:3, 0] * msub[0, 1:3] / msub[0, 0]).applyfunc(sp.cancel)
ck.check("EXACT", "with q = b_47 = B, b = a_47 = r^2/B, beta = -r, v = -r t/B, j = e_47 = E the leaf block has determinant 0 "
         "and the Schur complement of the 3 x 3 metric with respect to the entry q = B is diag(0, E - kappa t^2), "
         "kappa = 1 + 1/B: the metric is PSD iff E >= kappa t^2; Delta = (B-r)^2 (E - kappa t^2)/B and "
         "b(j-t^2)-v^2 = r^2 (E - kappa t^2)/B, so the construction is defined for B != r, r > 0, E > kappa t^2",
         sp.cancel(ell.subs(mp, simultaneous=True)) == 0
         and schur == sp.Matrix([[0, 0], [0, sp.cancel(rem)]])
         and sp.cancel(Delta.subs(mp, simultaneous=True) - (Bq - r) ** 2 * rem / Bq) == 0
         and sp.cancel(R0.subs(mp, simultaneous=True) - r ** 2 * rem / Bq) == 0)
p2 = kk + (1 - eta * t) ** 2 / rem
c2 = rho ** 2 * t * t / rem
ck.check("EXACT", "p* becomes k + (1 - eta t)^2/(E - kappa t^2) and c* becomes rho^2 t^2/(E - kappa t^2) with "
         "rho = 1 + r/B, eta = (1 - a)/(B - r), k = B(1-r)^2/(B-r)^2 (a = r^2/B)",
         sp.cancel(pstar.subs(mp, simultaneous=True) - p2) == 0 and sp.cancel(cstar.subs(mp, simultaneous=True) - c2) == 0)
G2 = Gm6.subs(mp, simultaneous=True).applyfunc(sp.cancel)
W2 = [Bq, p2, E, c2, 0, 1, r ** 2 / Bq, 0]
inst2 = {Bq: Fr(3), r: Fr(1, 2), E: Fr(2), t: Fr(1, 2)}
ck.check("EXACT", "COROLLARY 2, first statement: the substituted Gram matrix still satisfies the ten relations with weights "
         "(B, k + (1-eta t)^2/(E-kappa t^2), E, rho^2 t^2/(E-kappa t^2), 0, 1, r^2/B, 0); identity (*) multiplied out "
         "directly at (B, r, E, t) = (3, 1/2, 2, 1/2), i.e. b_47 = 3, a_47 = 1/12",
         gram_ok(G2, W2) and direct(G2, W2, inst2))
ts = eta * E / kap
ck.check("EXACT", "stationary choice t = eta E/kappa: E - kappa t^2 = E (kappa - eta^2 E)/kappa; thresholds become "
         "c_47 >= 1/E + k - eta^2/kappa and p_47 >= rho^2 eta^2 E/(kappa(kappa - eta^2 E)); and "
         "P kappa (kappa - eta^2 E) - rho^2 eta^2 E = P kappa^2 - eta^2 E (P kappa + rho^2), so for eta^2 E < kappa the "
         "second is E <= P kappa^2/(eta^2 (P kappa + rho^2)); this bound is itself < kappa/eta^2",
         sp.cancel(rem.subs(t, ts) - E * (kap - eta ** 2 * E) / kap) == 0
         and sp.cancel(p2.subs(t, ts) - (kk + 1 / E - eta ** 2 / kap)) == 0
         and sp.cancel(c2.subs(t, ts) - rho ** 2 * eta ** 2 * E / (kap * (kap - eta ** 2 * E))) == 0
         and sp.expand(Pp * kap * (kap - eta ** 2 * E) - rho ** 2 * eta ** 2 * E
                       - (Pp * kap ** 2 - eta ** 2 * E * (Pp * kap + rho ** 2))) == 0
         and sp.cancel(kap / eta ** 2 - Pp * kap ** 2 / (eta ** 2 * (Pp * kap + rho ** 2))
                       - kap * rho ** 2 / (eta ** 2 * (Pp * kap + rho ** 2))) == 0)

ck.section("Document 2, Corollaries 3 and 4 (rational sections)  [pp. 3-4]")
KAP, RHO, ETA, EP, PC = sp.symbols("kappa rho eta E P", positive=True)
D0 = KAP + RHO ** 2 / PC
M0 = D0 + ETA ** 2 * EP
NN = M0 ** 2 - KAP * ETA ** 2 * EP
t3 = ETA * EP / M0
ck.check("EXACT", "Corollary 3, t = eta e/M0: e - kappa t^2 = e N/M0^2; consumed c = k + D0^2/(e N); consumed p = "
         "rho^2 eta^2 e/N; N = D0^2 + eta^2 e (2 D0 - kappa + eta^2 e); N - rho^2 eta^2 e/p = D0^2 + D0 eta^2 e + eta^4 e^2; "
         "2 D0 - kappa = kappa + 2 rho^2/p > 0, so N > D0^2 > 0 (eta != 0) and the consumed p is < p",
         sp.cancel(EP - KAP * t3 ** 2 - EP * NN / M0 ** 2) == 0
         and sp.cancel((1 - ETA * t3) ** 2 / (EP - KAP * t3 ** 2) - D0 ** 2 / (EP * NN)) == 0
         and sp.cancel(RHO ** 2 * t3 ** 2 / (EP - KAP * t3 ** 2) - RHO ** 2 * ETA ** 2 * EP / NN) == 0
         and sp.expand(NN - D0 ** 2 - ETA ** 2 * EP * (2 * D0 - KAP + ETA ** 2 * EP)) == 0
         and sp.expand(NN - RHO ** 2 * ETA ** 2 * EP / PC - (D0 ** 2 + D0 * ETA ** 2 * EP + ETA ** 4 * EP ** 2)) == 0
         and sp.expand(2 * D0 - KAP - (KAP + 2 * RHO ** 2 / PC)) == 0)
AP = ETA ** 2 * EP
t4 = ETA * EP / D0
N0 = D0 ** 2 - KAP * AP
ck.check("EXACT", "Corollary 4, A = eta^2 e.  Branch A <= D0, t = eta e/D0: e - kappa t^2 = e N0/D0^2 with "
         "N0 = D0^2 - kappa A = D0 (D0 - kappa) + kappa (D0 - A) > 0; consumed c = k + (D0-A)^2/(e N0); unused "
         "p = p D0 (D0-A)/N0 >= 0.  Branch A >= D0, t = 1/eta: e - kappa t^2 = (A - kappa)/eta^2 > 0, consumed "
         "p = rho^2/(A - kappa) with p - rho^2/(A-kappa) = p (A - D0)/(A - kappa) >= 0, consumed c = k.  The branches "
         "agree at A = D0",
         sp.cancel(EP - KAP * t4 ** 2 - EP * N0 / D0 ** 2) == 0
         and sp.expand(N0 - (D0 * (D0 - KAP) + KAP * (D0 - AP))) == 0
         and sp.cancel((1 - ETA * t4) ** 2 / (EP - KAP * t4 ** 2) - (D0 - AP) ** 2 / (EP * N0)) == 0
         and sp.cancel(PC - RHO ** 2 * t4 ** 2 / (EP - KAP * t4 ** 2) - PC * D0 * (D0 - AP) / N0) == 0
         and sp.cancel(EP - KAP / ETA ** 2 - (AP - KAP) / ETA ** 2) == 0
         and sp.cancel(RHO ** 2 * (1 / ETA) ** 2 / (EP - KAP / ETA ** 2) - RHO ** 2 / (AP - KAP)) == 0
         and sp.cancel(PC - RHO ** 2 / (AP - KAP) - PC * (AP - D0) / (AP - KAP)) == 0
         and (1 - ETA * (1 / ETA)) == 0)
Aa, Dz, Kz = sp.symbols("A D0 kappa", positive=True)
brace = 2 * Dz * (Dz - Kz) + Kz * Aa - Aa ** 2
ck.check("EXACT", "comparison of the two sections for 0 < A <= D0: D0^2/(e N) - (D0-A)^2/(e N0) = "
         "A^2 {2 D0 (D0-kappa) + kappa A - A^2}/(e N0 N), N = (D0+A)^2 - kappa A; the brace equals 2 D0 (D0 - kappa) at "
         "A = 0 and D0 (D0 - kappa) at A = D0 and is concave in A, hence positive on [0, D0] (D0 > kappa)",
         sp.cancel(Dz ** 2 / ((Dz + Aa) ** 2 - Kz * Aa) - (Dz - Aa) ** 2 / (Dz ** 2 - Kz * Aa)
                   - Aa ** 2 * brace / ((Dz ** 2 - Kz * Aa) * ((Dz + Aa) ** 2 - Kz * Aa))) == 0
         and sp.expand(brace.subs(Aa, 0) - 2 * Dz * (Dz - Kz)) == 0
         and sp.expand(brace.subs(Aa, Dz) - Dz * (Dz - Kz)) == 0 and sp.diff(brace, Aa, 2) == -2)
ck.info("PROOF-READ: in Corollaries 2-4 the hypotheses b_47 > 1, 0 < a_47 < 1/b_47 give 0 < r < 1 < B, so B != r, "
        "eta = (1-a)/(B-r) > 0, kappa > 1, D0 > kappa; every t used satisfies E > kappa t^2 strictly, so the Gram matrix "
        "of Corollary 2 is defined and PSD; the unused p and c are added as multiples of minor vectors.")

# ========================================================================================
ck.section("Document 3 (five-ray core note), Theorem 1  [pp. 1-2]")
D3 = s_ * j - (s_ + 4) * t * t
lam3 = (-h_ * (j - t * t) + 2 * t * (t - 1)) / D3
mu3 = (-2 * t * h_ + s_ * (t - 1)) / D3
metric3 = sp.Matrix([[q, beta, -t, 0], [beta, b, -t, 0], [-t, -t, j - t * t, 0], [0, 0, 0, 1]])   # basis (Q, B, J0, D)
V3 = sp.Matrix([[1, 0, 0, 0], [lam3, lam3, mu3, 0], [0, 0, 1, t], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1],
                [0, 1, 0, 0], [0, 0, 0, 0]])
Gm5 = (V3 * metric3 * V3.T).applyfunc(sp.cancel)
p0 = ((j - t * t) * h_ ** 2 - 4 * t * h_ * (t - 1) + s_ * (t - 1) ** 2) / D3
W5 = [q, p0, j, 0, 0, 1, b, 0]
ck.check("EXACT", "det of the 3 x 3 metric of (Q, B, J0) = ell j - (ell + q + b - 2 beta) t^2, so the hypotheses q, b > 0, "
         "ell > 0, j > (1 + (q+b-2beta)/ell) t^2 say that the metric is positive definite; Delta = s j - (s+4) t^2 is the "
         "Gram determinant of (Q+B, J0); squared lengths (q, p0, j, 0, 0, 1, b, 0) with p0 the right side of (1)",
         sp.expand(metric3[:3, :3].det() - (ell * j - (ell + q + b - 2 * beta) * t * t)) == 0
         and sp.expand(D3 - sp.Matrix([[s_, -2 * t], [-2 * t, j - t * t]]).det()) == 0
         and all(sp.cancel(Gm5[i, i] - z) == 0 for i, z in enumerate(W5)))
ck.check("EXACT", "THEOREM 1 (sufficiency): the Gram matrix of (Q, P, J, 0, 0, D, B, 0) satisfies the ten relations "
         "identically in (q, b, j, beta, t) with weights (q, p0, j, 0, 0, 1, b, 0): (*) holds for "
         "p0 f1 + j f59 + q f27 + b f6 + f28", gram_ok(Gm5, W5))
inst5 = {q: Fr(1, 8), b: Fr(4), j: Fr(4), beta: Fr(-1, 2), t: Fr(1, 4)}
ck.check("EXACT", "identity (*) multiplied out directly at the document's example (q,b,j,beta,t) = (1/8,4,4,-1/2,1/4), where "
         "the metric is positive definite and p0 = 399/1543",
         minors_positive(metric3, inst5) and direct(Gm5, W5, inst5)
         and sp.cancel(p0.subs({k_: sp.Rational(x.numerator, x.denominator) for k_, x in inst5.items()}))
         == sp.Rational(399, 1543))
ck.check("EXACT", "it is the slice v = t of the six-ray family after t -> -t: p*(t -> -t, v -> -t) = p0 and c* = 0 there",
         sp.cancel(pstar.subs({t: -t, v: -t}, simultaneous=True) - p0) == 0
         and sp.cancel(cstar.subs({t: -t, v: -t}, simultaneous=True)) == 0)
ck.check("EXACT", "remarks on p. 2: p0 at t = 0 is 1/j + (1+beta)^2/(q+b+2beta); d p0/dt at t = 0 is "
         "2(2-q-b)/((q+b+2beta) j); and 399/1543 < 13/50 < 1/4 + (1 - 1/sqrt2)^2/(33/8 - sqrt2) "
         "[the last number is (25839 - 16800 sqrt2)/96100 + 13/50 and 25839^2 > 2*16800^2]",
         sp.cancel(p0.subs(t, 0) - (1 / j + h_ ** 2 / s_)) == 0
         and sp.cancel(sp.diff(p0, t).subs(t, 0) - 2 * (2 - q - b) / (s_ * j)) == 0
         and Fr(399, 1543) < Fr(13, 50) and 25839 ** 2 > 2 * 16800 ** 2
         and sp.simplify(sp.Rational(1, 4) + (1 - 1 / sp.sqrt(2)) ** 2 / (sp.Rational(33, 8) - sp.sqrt(2))
                         - sp.Rational(13, 50) - (25839 - 16800 * sp.sqrt(2)) / 96100) == 0)

ck.section("Document 3, Theorem 2 (classification) and the optimisation in t  [pp. 2-4]")
core5 = [0, 1, 2, 5, 6]                                     # Q, P, J, D, B
y5 = sp.symbols("y0:15")
c5 = [(i, k) for i in range(5) for k in range(i, 5)]
m5 = sp.zeros(5, 5)
for (i, k), x in zip(c5, y5):
    m5[i, k] = m5[k, i] = x
emb5 = sp.zeros(8, 8)
for i, oi in enumerate(core5):
    for k, ok_ in enumerate(core5):
        emb5[oi, ok_] = m5[i, k]
ent5 = [emb5[i, k] for (i, k) in gf.keys]
restr5 = sp.Matrix([comb(r_, ent5) for r_ in gf.relations]).jacobian(y5)
seven = sp.Matrix([m5[0, 1] + m5[0, 4] + m5[1, 4] + m5[3, 3], m5[0, 2] + m5[2, 3], m5[0, 3],
                   m5[1, 2] + m5[2, 4] + m5[3, 3], m5[1, 3], m5[2, 3] + m5[2, 4], m5[3, 4]]).jacobian(y5)
tm5 = sp.Matrix([[1, sp.Rational(-1, 2), 0, 0, 0], [sp.Rational(-1, 2), 2, -1, 0, sp.Rational(-1, 2)],
                 [0, -1, 1, 0, 0], [0, 0, 0, 1, 0], [0, sp.Rational(-1, 2), 0, 0, 1]])
ck.check("EXACT", "when c = e = a = 0 the ten relations restricted to the five vectors (Q, P, J, D, B) have rank 7 and the "
         "same row space as the seven equations displayed in the document; its 5 x 5 template is positive definite "
         "(leading minors %s), has |D|^2 = 1 and satisfies the seven equations"
         % [str(tm5[:k, :k].det()) for k in range(1, 6)],
         restr5.rank() == 7 and seven.rank() == 7 and restr5.col_join(seven).rank() == 7
         and all(tm5[:k, :k].det() > 0 for k in range(1, 6)) and tm5[3, 3] == 1
         and seven * sp.Matrix([tm5[i, k] for i, k in c5]) == sp.zeros(7, 1))
A_, B_, C_, Dn, En = s_ - h_ ** 2 - 4 * h_, 4 * h_ - 2 * s_, j * h_ ** 2 + s_, s_ * j, s_ + 4
kq = 2 - q - b
form2 = h_ ** 2 / s_ + (s_ + kq * t) ** 2 / (s_ * D3)
uu = (q + b - 2 * beta) / ell
tst = -kq * j / (s_ + 4)
ck.check("EXACT", "p0 = (A t^2 + B t + C)/(D - E t^2) with the document's A..E; the numerator of its t-derivative is "
         "B E t^2 + 2 (A D + E C) t + B D; p0 = h^2/s + (s + k t)^2/(s (s j - (s+4) t^2)) with k = 2 - q - b; "
         "d p0/dt = 2 (s + k t)(k j + (s+4) t)/(s j - (s+4) t^2)^2; u - 4/s = (q-b)^2/(s ell) for u = (q+b-2beta)/ell; "
         "at t = -k j/(s+4) the value is 1/j + (s+2k)/4 + k^2/(4(s+4))",
         sp.cancel(p0 - (A_ * t * t + B_ * t + C_) / (Dn - En * t * t)) == 0
         and sp.expand((2 * A_ * t + B_) * (Dn - En * t * t) + 2 * En * t * (A_ * t * t + B_ * t + C_)
                       - (B_ * En * t * t + 2 * (A_ * Dn + En * C_) * t + B_ * Dn)) == 0
         and sp.cancel(p0 - form2) == 0
         and sp.cancel(sp.diff(p0, t) - 2 * (s_ + kq * t) * (kq * j + (s_ + 4) * t) / D3 ** 2) == 0
         and sp.cancel(uu - 4 / s_ - (q - b) ** 2 / (s_ * ell)) == 0
         and sp.cancel(p0.subs(t, tst) - (1 / j + (s_ + 2 * kq) / 4 + kq ** 2 / (4 * (s_ + 4)))) == 0)
kS = sp.Symbol("k", positive=True)
ck.check("EXACT", "with beta = 1 - q - b (so s = k): 1/j + 3k/4 + k^2/(4(k+4)) = 1/j + k(k+3)/(k+4) < 1/j + k, and "
         "|t*| < T = sqrt(j/(1+u)) is j < (k+4)^2/(k^2 (1+u))",
         sp.cancel(sp.Rational(3, 4) * kS + kS ** 2 / (4 * (kS + 4)) - kS * (kS + 3) / (kS + 4)) == 0
         and sp.cancel(kS - kS * (kS + 3) / (kS + 4) - kS / (kS + 4)) == 0
         and sp.cancel(s_.subs(beta, 1 - q - b) - kq) == 0)
ck.info("PROOF-READ (not needed for any support inequality): the projection argument of Theorem 2 and the choice of t_opt "
        "from the sign of the derivative are correct as written; they concern only the size of the certified region.")

ck.section("Document 3, the transports  [p. 4]")
im1 = K.doc_action_facets((0, 3, 1, 2), (3, 1, 0, 2), True, [59, 53, 55, 67, 54])
im2 = K.doc_action_facets((0, 2, 3, 1), (3, 1, 0, 2), False, [1, 59, 6, 27, 28])
g25 = K.generators(K.PAPER_T[25])
ck.check("EXACT", "row (0,3,1,2), column (3,1,0,2), transposed: [59,53,55,67,54] -> [1,59,6,27,28]; row (0,2,3,1), column "
         "(3,1,0,2), not transposed: [1,59,6,27,28] -> [45,5,7,0,35] (supplement's convention); f59,f53,f55,f67,f54 are "
         "generators of the paper's T_25 (document label T24), f45,f5,f7,f0,f35 of the paper's T_65",
         im1 == [1, 59, 6, 27, 28] and im2 == [45, 5, 7, 0, 35] and all(x in g25 for x in [59, 53, 55, 67, 54])
         and all(x in g65 for x in im2), "computed %s, %s" % (im1, im2))

ck.section("Document 3, Corollary 3 (interior leaf-product capacity)  [pp. 4-5]")
sub3 = {t: 0, beta: 1 - q - b}
G3 = Gm5.subs(sub3, simultaneous=True).applyfunc(sp.cancel)
W3 = [q, 1 / j + 2 - q - b, j, 0, 0, 1, b, 0]
m3 = metric3.subs(sub3, simultaneous=True)
ck.check("EXACT", "COROLLARY 3 at d = 1, t = 0, beta = 1 - q - b: the metric is diag-block [[q, 1-S],[1-S, b]] + (j) + (1), "
         "PSD iff q, b >= 0, j >= 0 and (S-1)^2 <= qb (S = q + b); s = h = 2 - S, Delta = (2-S) j; p0 = 1/j + 2 - S; the "
         "Gram matrix satisfies the ten relations with weights (q, 1/j + 2 - q - b, j, 0, 0, 1, b, 0) wherever S != 2, j != 0: "
         "the certificate is direct on the closed leaf condition, the point S = 2 excepted",
         m3[0, 2] == 0 and m3[1, 2] == 0 and sp.expand(m3[:2, :2].det() - (q * b - (q + b - 1) ** 2)) == 0
         and sp.cancel(p0.subs(sub3, simultaneous=True) - (1 / j + 2 - q - b)) == 0 and gram_ok(G3, W3)
         and direct(G3, W3, {q: Fr(1, 2), b: Fr(1), j: Fr(3)})              # interior: (S-1)^2 = 1/4 < 1/2
         and direct(G3, W3, {q: Fr(1, 3), b: Fr(1, 3), j: Fr(2)})           # boundary: (S-1)^2 = 1/9 = qb
         and direct(G3, W3, {q: Fr(0), b: Fr(1), j: Fr(5, 2)}))             # one leaf zero, the other one
jj = sp.Symbol("j", positive=True)
Gb = sp.zeros(8, 8)                                         # q = b = 1: B = -Q, P = -J/j, J orthogonal to Q, D
Gb[0, 0] = Gb[6, 6] = 1
Gb[0, 6] = Gb[6, 0] = -1
Gb[2, 2] = jj
Gb[1, 1] = 1 / jj
Gb[1, 2] = Gb[2, 1] = -1
Gb[5, 5] = 1
piv, resid = K.ldl_certificate(Gb, [0, 2, 5])
ck.check("EXACT", "the excepted point q = b = 1 (S = 2): the Gram matrix with B = -Q, |Q| = 1, J orthogonal to Q and D, "
         "|J|^2 = j, P = -J/j, |D| = 1 satisfies the ten relations with weights (1, 1/j, j, 0, 0, 1, 1, 0) and is PSD "
         "(pivots 1, j, 1, zero residual): p >= 1/j = 1/j + 2 - S there, without a limit",
         gram_ok(Gb, [1, 1 / jj, jj, 0, 0, 1, 1, 0]) and piv == [1, jj, 1] and resid == sp.zeros(8, 8)
         and direct(Gb, [1, 1 / jj, jj, 0, 0, 1, 1, 0], {jj: Fr(3, 2)}))
dS = sp.Symbol("d", positive=True)
ck.check("EXACT", "homogeneity: with (q, b, j, p) -> (q, b, j, p)/d the condition and the threshold become "
         "(q + b - d)^2 <= qb and p >= d^2/j + 2d - q - b; and (S-1)^2 <= qb <= S^2/4 forces 2/3 <= S <= 2",
         sp.cancel(dS * (1 / (j / dS) + 2 - q / dS - b / dS) - (dS ** 2 / j + 2 * dS - q - b)) == 0
         and sp.expand(sp.Symbol("S") ** 2 / 4 - (sp.Symbol("S") - 1) ** 2
                       - sp.Rational(3, 4) * (sp.Symbol("S") - sp.Rational(2, 3)) * (2 - sp.Symbol("S"))) == 0)
Ss = sp.Symbol("S")
ck.check("EXACT", "closing remarks: d/dbeta (1+beta)^2/(S+2beta) = 2(1+beta)(S+beta-1)/(S+2beta)^2; the quadratic "
         "(q+b-1)^2 - qb in q has discriminant b(4-3b)",
         sp.cancel(sp.diff((1 + beta) ** 2 / (Ss + 2 * beta), beta)
                   - 2 * (1 + beta) * (Ss + beta - 1) / (Ss + 2 * beta) ** 2) == 0
         and sp.expand(sp.discriminant(sp.expand((q + b - 1) ** 2 - q * b), q) - b * (4 - 3 * b)) == 0)
ck.info("STATEMENT (Corollary 3): d > 0, j > 0, q, b >= 0, (q + b - d)^2 <= qb, p >= d^2/j + 2d - q - b  =>  "
        "Phi >= 0 on SO(4) for p f1 + j f59 + q f27 + b f6 + d f28 (plus any minor vectors).  On the T_65 cone: "
        "(a + b - d)^2 <= ab, e > 0, c >= d^2/e + 2d - a - b  =>  Phi >= 0 for c f45 + e f5 + a f7 + b f0 + d f35.")

# ----------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence only)")
    cvv = {jx: np.array(K.F_VALUES[jx], float) for jx in CANON}
    rng = np.random.default_rng(47400)
    f_p = sp.lambdify((q, b, j, beta, t, v), pstar, "numpy")
    f_c = sp.lambdify((q, b, j, beta, t, v), cstar, "numpy")
    worst = float("inf")
    n = 0
    while n < 12:
        qv, bv = np.exp(rng.normal(size=2))
        bt = rng.uniform(-0.9, 0.9) * np.sqrt(qv * bv)
        tv_, vv_ = rng.normal(size=2) * 0.6
        lf = qv * bv - bt ** 2
        jv = tv_ ** 2 + (bv * tv_ ** 2 - 2 * bt * tv_ * vv_ + qv * vv_ ** 2) / lf + np.exp(rng.normal())
        ps, cs = float(f_p(qv, bv, jv, bt, tv_, vv_)), float(f_c(qv, bv, jv, bt, tv_, vv_))
        vals = qv * cvv[27] + ps * cvv[1] + jv * cvv[59] + cs * cvv[5] + cvv[28] + bv * cvv[6]
        worst = min(worst, K.numeric_min_phi(vals, rng, restarts=8))
        n += 1
    ctrl = K.numeric_min_phi(cvv[28], rng, restarts=8)
    ck.check("NUMERIC", "six-ray Theorem 1 at its thresholds (p = p*, c = c*) for 12 random admissible (q, b, beta, t, v, j) "
             "(seed 47400): smallest Phi found over SO(4) by multistart BFGS is %.3e (>= -1e-9); control: the same search "
             "finds %.4f for f28 alone" % (worst, ctrl), worst > -1e-9 and ctrl < -0.01)
    worst3 = float("inf")
    for _ in range(10):
        Sv = rng.uniform(0.7, 1.95)
        # q + b = S and (S-1)^2 <= qb: q in the allowed interval
        disc = Sv ** 2 - 4 * (Sv - 1) ** 2
        lo, hi = (Sv - np.sqrt(disc)) / 2, (Sv + np.sqrt(disc)) / 2
        qv = rng.uniform(lo, hi)
        bv = Sv - qv
        jv = np.exp(rng.normal())
        vals = qv * cvv[27] + (1 / jv + 2 - Sv) * cvv[1] + jv * cvv[59] + cvv[28] + bv * cvv[6]
        worst3 = min(worst3, K.numeric_min_phi(vals, rng, restarts=8))
    ck.check("NUMERIC", "Corollary 3 of the core note at p = 1/j + 2 - q - b for 10 random (q, b, j) with (q+b-1)^2 <= qb "
             "(boundary included): smallest Phi found %.3e (>= -1e-9)" % worst3, worst3 > -1e-9)

ck.finish()
