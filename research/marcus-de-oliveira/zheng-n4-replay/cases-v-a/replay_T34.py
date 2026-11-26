#!/usr/bin/env python3
"""Independent replay of case T_34 = {2,3,5,9,10,12} (group V) of D. Zheng, "The Marcus-de
Oliveira conjecture for n = 4" (October 2026).

Proof document replayed (author's supplement, read as a PDF; never executed):
    proofs/manin_transfer_20261021/blind_moc4_bad70_full_safe_disjunction_20261022.pdf  (3 pages)

Run:   MOC4_SUPP=/path/to/supplement  python replay_T34.py
Exit code 0 iff every check passes.  Lines tagged [EXACT] come from exact arithmetic (Python
ints, Fractions, python-flint integer polynomials, sympy); lines tagged [NUMERIC] are floating
point falsification attempts (evidence, not proof).

What is read from MOC4_SUPP, and how
  * the PDF and four JSON files are hashed (provenance);
  * DATA, parsed with json.load and a strict term parser (no eval, nothing executed):
      - the sixteen cubic "features", four signed coordinate maps and the 8 x 8 matrix T of
        the degree-six Gram certificate behind the document's conditions (PB), (CQ)
        (the certificate exists only as data: the PDF, p. 1, describes the Gram matrices
        and says the identities are in its verifier).  Every polynomial the certificate is
        an identity ABOUT (the Leibniz sums on the quaternion model of SO(4)) is rebuilt
        here from the definitions; the features are the witness of the identity;
      - the stored polynomial R9 and the stored term counts / digests of the twelve
        Bernstein numerators, for COMPARISON with the polynomials rebuilt here.
The supplement's programs are not imported or executed.
"""
import json
import os
import random
import re
import sys
import time
from fractions import Fraction
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay_T34_common as K   # noqa: E402

NU = 34
T34 = (2, 3, 5, 9, 10, 12)                     # Appendix E of the paper (p. 60), entry T_34
DOC_T = (2, 3, 5, 9, 10, 12)                   # the document, p. 1
DOC_SRC = [0, 1, 2, 12, 39, 67, 70, 72]        # the document's generators, p. 1
NAMES = ["p", "r", "a", "b", "c", "q", "d", "s"]   # the document's coefficients, eq. (1)
POS = {n: i for i, n in enumerate(NAMES)}
GSTAR = ((2, 3, 0, 1), (1, 0, 3, 2), 0)        # the non-trivial element of G fixing T (found below)

MANIN = os.path.join("proofs", "manin_transfer_20261021")
FILES = {
    "pdf": (os.path.join(MANIN, "blind_moc4_bad70_full_safe_disjunction_20261022.pdf"),
            "479f5594c498be1fbaf2cc637b35fbffd281c185d3266ae24da0bdb13e4db037"),
    "bern": (os.path.join(MANIN, "blind_moc4_bad70_full_safe_disjunction_20261022.json"),
             "d59cc38bc1a3a9c0835312aefc5285972b345097724bc2bd7ab112da852d0c9a"),
    "r9": (os.path.join(MANIN, "blind_moc4_next_singlebad_R9_restriction_20261022.json"),
           "c7e4e885addcbaedcd22d5cfcf0a20e6d6a9dd8a8af7948d9504a67c392c3910"),
    "gramT": (os.path.join(MANIN, "blind_moc4_next_singlebad_budget_PSD8_20261022.json"),
              "3d5ecfb9b6cb8f2e54dd410ffa0a46e6373a71caf7cee022958694c3a481ac0c"),
    "gramF": (os.path.join(MANIN, "blind_moc4_next_singlebad_cubic_reduction_20261022.json"),
              "c8974ea5f4bb3c5bae8b8de3387855a8945860e501dbbfd9495701006cf37f96"),
}


# ----------------------------------------------------------------------------------------
# 0. provenance
# ----------------------------------------------------------------------------------------


def load_supplement(ck):
    ck.section("Provenance: files of the supplement that were read (hashed; data only)")
    import flint
    import sympy
    import numpy
    import scipy
    ck.info("python %s, python-flint %s, sympy %s (exact part); numpy %s, scipy %s (NUMERIC part only)"
            % (sys.version.split()[0], flint.__version__, sympy.__version__, numpy.__version__,
               scipy.__version__))
    root = os.environ.get("MOC4_SUPP")
    if not ck.check("EXACT", "MOC4_SUPP is set (the Gram certificate of (PB), (CQ) exists only as "
                    "data of the supplement and is needed)", bool(root)):
        return None
    data = {}
    for key, (rel, sha) in FILES.items():
        path = os.path.join(root, rel)
        if not os.path.isfile(path):
            ck.check("EXACT", "file present: " + rel, False)
            return None
        h = K.sha256_file(path)
        ck.check("EXACT", "sha256 of %s is the one of the file read for this replay"
                 % os.path.basename(rel), h == sha, h)
        if key != "pdf":
            with open(path) as fh:
                data[key] = json.load(fh)
    return data


# ----------------------------------------------------------------------------------------
# 1. contact set, generators, scope
# ----------------------------------------------------------------------------------------


def check_contact_set_and_generators(ck):
    ck.section("Contact set, generators of the contact cone, scope (document p. 1)")
    ck.check("EXACT", "the document's contact set %s equals the paper's T_34 (Appendix E) as a set of "
             "lexicographic indices: group element = identity (the file name 'bad70' and the labels "
             "'T33' are historical; the PDF is titled with T_34)" % list(DOC_T),
             tuple(sorted(DOC_T)) == T34)
    rep = K.orbit_rep(T34)
    ck.check("EXACT", "lexicographically smallest set in the G-orbit of T_34 is {0,1,2,6,10,23}",
             rep == (0, 1, 2, 6, 10, 23), rep)
    gens = K.generators(T34)
    ck.check("EXACT", "F_T = {j : f_j = 0 on T} = %s equals the document's list" % gens,
             gens == DOC_SRC)
    nm = [j for j in gens if not K.F_IS_MINOR[j]]
    ck.check("EXACT", "the only non-minor generator is f_70 = h_{4231,+}",
             nm == [70] and K.f_name(70) == "h_{4231,+}",
             ", ".join("f_%d = %s" % (j, K.f_name(j)) for j in gens))
    t0 = time.time()
    d, rays = K.extreme_rays_of_contact_cone(T34)
    want = sorted(K.F_VALUES[j] for j in gens)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} from P alone "
             "(no use of the 82-facet classification): the subspace {(Pc)_T = 0} has dimension %d and "
             "the cone has exactly the 8 rays f_j, j in F_T" % d, d == 8 and rays == want,
             "%d rays, %.1f s" % (len(rays), time.time() - t0))
    # coefficients are gaps
    excl = []
    for j in gens:
        e = [s for s in range(24) if K.F_VALUES[j][s] == 1
             and all(K.F_VALUES[k][s] == 0 for k in gens if k != j)]
        excl.append(e)
    ck.check("EXACT", "the cone is simplicial and every coefficient is a gap: for each generator there is a "
             "permutation at which it is 1 and the other seven vanish, so "
             "(p, r, a, b, c, q, d, s) = (gap_18, gap_16, gap_13, gap_8, gap_4, gap_14, gap_11, gap_7)",
             [e[0] if e else None for e in excl] == [18, 16, 13, 8, 4, 14, 11, 7], excl)
    ck.info("so g in N_T  <=>  g = p f0 + r f1 + a f2 + b f12 + c f39 + q f67 + d f70 + s f72 with all eight "
            "coefficients >= 0 (unique); a zero coefficient is exactly an extra contact, and non-strict "
            "contacts are included in the document's statement")
    stab = [K.G_ELEMENTS[k] for k, gm in enumerate(K.G_MAPS) if sorted(gm[s] for s in T34) == list(T34)]
    gm = K.G_MAPS[K.G_ELEMENTS.index(GSTAR)]
    img = [K.g_facet(gm, j) for j in gens]
    ck.check("EXACT", "stabiliser of T_34 in G = {identity, g*}, g* = (rho, rho', eps) = %s; g* permutes the "
             "generators as f0<->f67, f1<->f72, f12<->f39 and fixes f2, f70 (coefficients: p<->q, r<->s, b<->c)"
             % (GSTAR,), sorted(stab) == sorted([((0, 1, 2, 3), (0, 1, 2, 3), 0), GSTAR])
             and img == [67, 72, 2, 39, 12, 0, 70, 1], img)
    return gens


def check_covariance_gstar(ck):
    """Phi_{g*.c}(g*.O) = Phi_c(O) as an identity in a generic 4 x 4 matrix, and g*.O is a row
    and column permutation of O of sign +1."""
    R = K.PR(["o%d%d" % (i, j) for i in range(4) for j in range(4)])
    O = [[R.gens[4 * i + j] for j in range(4)] for i in range(4)]
    gO = K.g_matrix(GSTAR, O)
    gm = K.G_MAPS[K.G_ELEMENTS.index(GSTAR)]
    w, wg = K.leibniz_weights(O), K.leibniz_weights(gO)
    ok = all(wg[gm[s]] == w[s] for s in range(24))
    ck.check("EXACT", "covariance for g*: the Leibniz weights of a generic 4 x 4 matrix satisfy "
             "w_{g*.sigma}(g*.O) = w_sigma(O) for all 24 sigma (identity in 16 variables); g*.O is O with "
             "rows and columns permuted by even permutations, so g* maps SO(4) onto SO(4) and "
             "Phi_{g*.c}(g*.O) = Phi_c(O) for every c", ok and K.sgn(GSTAR[0]) * K.sgn(GSTAR[1]) == 1)


# ----------------------------------------------------------------------------------------
# 2. the three sufficient conditions (document Section 1)
# ----------------------------------------------------------------------------------------


def check_old_H(ck):
    ck.section("Condition (H) (document p. 1): the two 'endpoint inequalities' are instances of the old-H "
               "inequality, which the document neither states nor proves")
    k1 = K.find_group_elements((12, 2, 70), (39, 2, 70))
    ck.check("EXACT", "(f12; f2; f70) is the representative of old-H in PROOF.md 7.6; (f39; f2; f70) is its "
             "image under %d elements of G, among them the stabiliser element g*" % len(k1),
             len(k1) > 0 and K.G_ELEMENTS.index(GSTAR) in k1, [K.G_ELEMENTS[k] for k in k1])
    orb = K.orbit_of_tuple((12, 2, 70))
    ck.check("EXACT", "orbit of the ordered triple (aux f12, main f2; f70) under G has 288 elements "
             "(as stated in PROOF.md 7.6)", len(orb) == 288, len(orb))
    # the inequality proved in the supplement's "scalar compensation" note is SOC, a different family
    soc = K.orbit_of_tuple((1, 9, 22))
    rem = True
    partners = {}
    for aux in (12, 39):
        ks = sorted({t[2] for t in soc if set(t[:2]) == {aux, 2}})
        partners[aux] = ks
        for kk in ks:
            diff = [x + y - z for x, y, z in zip(K.F_VALUES[70], K.F_VALUES[2], K.F_VALUES[kk])]
            rem &= min(diff) < 0
    ck.check("EXACT", "remark: the SOC family (orbit of (f1, f9; f22), the inequality of the supplement's "
             "scalar-compensation note) contains triples with the pair {f12, f2} only for the non-minor "
             "vectors f26, f56, and with {f39, f2} only for f46, f69; for each such k the vector "
             "f70 + f2 - f_k has a negative value, so it is not a nonnegative combination of minor vectors: "
             "the endpoint inequalities of (H) are not obtained from SOC by adding minors",
             partners == {12: [26, 56], 39: [46, 69]} and rem, partners)
    # --- replay of the audit's proof of the representative (foundations/so4-atoms.md 3.2)
    R = K.PR(["P0", "P1", "P2", "P3", "Q0", "Q1", "Q2", "Q3"])
    y = R.gens
    Pv, Qv = y[:4], y[4:]
    W2 = [[1, 1, 1, 1], [1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1]]     # 2 W, W the Hadamard matrix
    ck.check("EXACT", "W = (1/2)[[1,1,1,1],[1,-1,-1,1],[1,-1,1,-1],[1,1,-1,-1]] is symmetric with W^2 = I "
             "(so (P, Q) = (Wp, Wq) ranges over S^3 x S^3 with (p, q))",
             all(sum(W2[i][k] * W2[k][j] for k in range(4)) == (4 if i == j else 0)
                 for i in range(4) for j in range(4))
             and all(W2[i][j] == W2[j][i] for i in range(4) for j in range(4)))

    def lin(row, vec):
        s = R.zero
        for c_, v in zip(row, vec):
            s = s + v * c_
        return s
    p = tuple(lin(W2[i], Pv) for i in range(4))            # p = 2 W P, |p|^2 = 4 |P|^2
    q = tuple(lin(W2[i], Qv) for i in range(4))
    nP = Pv[0] * Pv[0] + Pv[1] * Pv[1] + Pv[2] * Pv[2] + Pv[3] * Pv[3]
    nQ = Qv[0] * Qv[0] + Qv[1] * Qv[1] + Qv[2] * Qv[2] + Qv[3] * Qv[3]
    M = K.quat_matrix(p, q, R.one, R.zero, conj_q=True)    # x -> p x conj(q)
    MtM = K.ring_matmul([[M[j][i] for j in range(4)] for i in range(4)], M)
    w = K.leibniz_weights(M)
    detM = K.det4(M)
    n2 = nP * nQ * 16
    ck.check("EXACT", "model x -> p x conj(q) with p = 2WP, q = 2WQ: M^t M = 16|P|^2|Q|^2 I and "
             "det M = (16|P|^2|Q|^2)^2, so O = M/4 is in SO(4) for unit P, Q",
             all((MtM[i][j] == n2) if i == j else (len(MtM[i][j]) == 0) for i in range(4) for j in range(4))
             and detM == n2 * n2)
    vh = [Pv[i] * Pv[i] * nQ - Qv[i] * Qv[i] * nP for i in range(4)]   # v_i = P_i^2 - Q_i^2, homogenised
    x2 = K.minor_det(M, *K.minor_spec(2))
    x12 = K.minor_det(M, *K.minor_spec(12))
    Phi2 = K.phi(K.F_VALUES[2], M, w)
    Phi12 = K.phi(K.F_VALUES[12], M, w)
    Phi70 = K.phi(K.F_VALUES[70], M, w)
    ck.check("EXACT", "Phi_2 = x_2^2 and Phi_12 = x_12^2 with x_2 = det M[13,34], x_12 = det M[14,24] "
             "(polynomial identities in the model)", Phi2 == x2 * x2 and Phi12 == x12 * x12)
    ck.check("EXACT", "bridge: x_2 = v_0 + v_1 and x_12 = v_1 + v_2 for O = M/4 (homogenised identities "
             "x(M) = 16 (v_i + v_j)); sum v_i = 0",
             x2 == (vh[0] + vh[1]) * 16 and x12 == (vh[1] + vh[2]) * 16
             and len(vh[0] + vh[1] + vh[2] + vh[3]) == 0)
    delta = Pv[0] * Qv[3] + Pv[1] * Qv[2]
    Delta = Pv[2] * Qv[1] + Pv[3] * Qv[0]
    rhs = ((delta * delta + Delta * Delta) * nP * nQ - (vh[0] * vh[2] + vh[1] * vh[3])) * 512
    ck.check("EXACT", "bridge: Phi_2 + Phi_70 = 2[delta^2 + Delta^2 - (v_0 v_2 + v_1 v_3)] on S^3 x S^3, "
             "delta = P0Q3 + P1Q2, Delta = P2Q1 + P3Q0 (homogenised identity of bidegree (4,4))",
             Phi2 + Phi70 == rhs)
    ck.check("EXACT", "(v_0 + v_1)(v_1 + v_2) = v_0 v_2 - v_1 v_3 + v_1 (v_0+v_1+v_2+v_3): on sum v = 0, "
             "x_2 x_12 = v_0 v_2 - v_1 v_3",
             (vh[0] + vh[1]) * (vh[1] + vh[2]) == vh[0] * vh[2] - vh[1] * vh[3]
             + vh[1] * (vh[0] + vh[1] + vh[2] + vh[3]))
    R4 = K.PR(["r", "s", "u", "w"])
    r_, s_, u_, w_ = R4.gens
    ck.check("EXACT", "increment lemma: (r^2+u^2)(s^2+w^2) - (rs + uw)^2 = (rw - us)^2, i.e. "
             "sqrt((R+U)(S+W)) - sqrt(RS) >= sqrt(UW) for R, S, U, W >= 0",
             (r_ * r_ + u_ * u_) * (s_ * s_ + w_ * w_) - (r_ * s_ + u_ * w_) ** 2 == (r_ * w_ - u_ * s_) ** 2)
    ck.info("PROOF-READ (so4-atoms.md 3.2): Phi_2 + Phi_70 + 2|x_12 x_2| = 2[delta^2 + Delta^2 - 2 min(v0v2, v1v3)]. "
            "If the min is <= 0, done.  Otherwise (sum v = 0) v0, v2 < 0 < v1, v3 or the reverse; in the first "
            "case with a_i = P_i^2, A_i = Q_i^2: |delta| >= |X|, |Delta| >= |Y|, X = sqrt(a1A2) - sqrt(a0A3), "
            "Y = sqrt(a3A0) - sqrt(a2A1), and X + Y = [sqrt((A1+v1)(a2-v2)) - sqrt(A1a2)] + "
            "[sqrt((A3+v3)(a0-v0)) - sqrt(A3a0)] >= sqrt(-v1v2) + sqrt(-v3v0) > 0, hence delta^2 + Delta^2 >= "
            "(X+Y)^2/2 >= 2 sqrt(v0v1v2v3) >= 2 min(v0v2, v1v3).  The reverse case exchanges a and A.")
    # --- identification with inequality (star) of ../../star-orthogonal (Theorem O)
    import itertools
    PM = [[[1 if s[i] == j else 0 for j in range(4)] for i in range(4)] for s in K.PERMS]

    def relabel(Mx, tr, rp, cp):
        if tr:
            Mx = [[Mx[j][i] for j in range(4)] for i in range(4)]
        return [[Mx[rp[i]][cp[j]] for j in range(4)] for i in range(4)]

    def D(Mx, S, Tt):
        return (Mx[S[0]][Tt[0]] * Mx[S[1]][Tt[1]] - Mx[S[0]][Tt[1]] * Mx[S[1]][Tt[0]]) ** 2

    def g0(Mx):
        return (D(Mx, (0, 2), (0, 2)) + D(Mx, (0, 3), (0, 3)) + D(Mx, (0, 2), (0, 3))
                + D(Mx, (0, 2), (1, 2)) - D(Mx, (0, 1), (2, 3)))
    hits = {12: [], 39: []}
    for tr in (0, 1):
        for rp in itertools.permutations(range(4)):
            for cp in itertools.permutations(range(4)):
                rel = [relabel(Mx, tr, rp, cp) for Mx in PM]
                if tuple(g0(Mx) for Mx in rel) != K.F_VALUES[70]:
                    continue
                if tuple(D(Mx, (0, 1), (0, 1)) for Mx in rel) != K.F_VALUES[2]:
                    continue
                side = tuple(D(Mx, (0, 2), (0, 3)) for Mx in rel)
                for aux in (12, 39):
                    if side == K.F_VALUES[aux]:
                        hits[aux].append((tr, rp, cp))
    ck.check("EXACT", "identification with inequality (star) of ../../star-orthogonal (Theorem O): there are "
             "relabellings O -> O~ (row and column permutations, with or without transposition) with "
             "Phi_70 = g0(O~), Phi_2 = D[12,12](O~) and Phi_12 = D[13,14](O~) (resp. Phi_39 = D[13,14](O~)) "
             "at all 24 permutation matrices, hence on SO(4) (both sides are Leibniz sums of vectors of V_R "
             "with the same 24 values); (star) for O~ in O(4) is then the old-H inequality for (f12; f2; f70) "
             "(resp. (f39; f2; f70))", len(hits[12]) > 0 and len(hits[39]) > 0,
             "%d and %d relabellings, e.g. %s and %s" % (len(hits[12]), len(hits[39]),
                                                         hits[12][:1], hits[39][:1]))
    # --- allocation algebra of (H)
    import sympy
    b, c, d = sympy.symbols("b c d", positive=True)
    db, dc = d * b / (b + c), d * c / (b + c)
    ck.check("EXACT", "allocation of (H): with d_b = d b/(b+c), d_c = d c/(b+c): d_b + d_c = d and "
             "d_b + d_b^2/b + d_c + d_c^2/c = d + d^2/(b+c)",
             sympy.cancel(db + dc - d) == 0
             and sympy.cancel(db + db ** 2 / b + dc + dc ** 2 / c - d - d ** 2 / (b + c)) == 0)
    ck.info("PROOF-READ: if b, c > 0 and a >= d + d^2/(b+c), put a_b = d_b + d_b^2/b, a_c = d_c + d_c^2/c "
            "(a_b + a_c <= a); then g = [a_b f2 + b f12 + d_b f70] + [a_c f2 + c f39 + d_c f70] + "
            "(a - a_b - a_c) f2 + p f0 + r f1 + q f67 + s f72, the two brackets satisfy the old-H "
            "coefficient condition t(c' - d') >= d'^2 with equality, and the rest are minors.  If c = 0 < b "
            "(or b = 0 < c), (H) is the old-H condition for (f12; f2; f70) (or (f39; f2; f70)) itself: "
            "no continuity argument is needed.  If b = c = 0, (H) is void.")


def parse_cubic(s):
    """Strict parser for a sum of signed monomials  [int*]v[**k]*v[**k]...[/int]  in the variables
    p0..p3, q0..q3.  Returns {exponent 8-tuple: Fraction}.  Raises on anything else."""
    toks = s.replace(" ", "")
    if not re.fullmatch(r"[pq0-9+\-*/]+", toks):
        raise ValueError("unexpected character in feature string")
    terms = re.findall(r"[+-]?[^+-]+", toks)
    if "".join(terms) != toks:
        raise ValueError("cannot split feature string")
    term_re = re.compile(r"(?:(\d+)\*)?((?:[pq][0-3](?:\*\*\d)?)(?:\*[pq][0-3](?:\*\*\d)?)*)(?:/(\d+))?")
    out = {}
    for t in terms:
        sign = 1
        if t[0] == "+":
            t = t[1:]
        elif t[0] == "-":
            sign, t = -1, t[1:]
        m = term_re.fullmatch(t)
        if not m:
            raise ValueError("bad term %r" % t)
        co = Fraction(int(m.group(1) or 1), int(m.group(3) or 1)) * sign
        ex = [0] * 8
        for var, idx, pw in re.findall(r"([pq])([0-3])(?:\*\*(\d))?", m.group(2)):
            ex[(0 if var == "p" else 4) + int(idx)] += int(pw or 1)
        ex = tuple(ex)
        out[ex] = out.get(ex, 0) + co
    return {e: c for e, c in out.items() if c}


def check_chain_gram(ck, data):
    ck.section("Conditions (PB), (CQ) (document p. 1): the degree-six Gram identities, verified exactly "
               "from the supplement's certificate data")
    cub, bud = data["gramF"], data["gramT"]
    ck.check("EXACT", "certificate data refer to this face: contacts %s, source order %s"
             % (cub.get("contacts"), cub.get("source_facets")),
             cub.get("contacts") == list(T34) and cub.get("source_facets") == DOC_SRC
             and bud.get("contacts") == list(T34) and bud.get("source_order") == DOC_SRC)
    R = K.PR(["p0", "p1", "p2", "p3", "q0", "q1", "q2", "q3"])
    y = R.gens
    feats = [[parse_cubic(s) for s in cub["features"][e]] for e in range(2)]
    ck.check("EXACT", "data: 2 x 8 features, each a homogeneous cubic in (p0..p3, q0..q3) with coefficients "
             "in (1/2)Z (parsed by a strict term parser)",
             len(feats) == 2 and all(len(f) == 8 for f in feats)
             and all(sum(ex) == 3 and (2 * co).denominator == 1
                     for e in range(2) for f in feats[e] for ex, co in f.items()))
    F2 = [[R.from_terms([(ex, int(co * 2)) for ex, co in f.items()]) for f in feats[e]] for e in range(2)]
    Tm = [[Fraction(x) for x in row] for row in bud["T"]]
    Ti = K.mat_inv(Tm)
    ck.check("EXACT", "data: T is an invertible 8 x 8 rational matrix and T^{-1} is an integer matrix",
             len(Tm) == 8 and all(x.denominator == 1 for row in Ti for x in row))
    Ti = [[int(x) for x in row] for row in Ti]
    maps = cub["left_four"]
    okm = len(maps) == 4 and all(sorted(m[0] for m in mp) == list(range(8)) and all(m[1] in (1, -1) for m in mp)
                                 for mp in maps)
    ck.check("EXACT", "data: four signed coordinate permutations S_0..S_3 of (p0..p3, q0..q3)", okm)
    # the S_j are (p, q) -> (e p, q conj(e)), e = 1, i, j, k
    units = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    oku = True
    for mp, e in zip(maps, units):
        ep = K.quat_mul(tuple(R.const(x) for x in e), tuple(y[:4]))
        qe = K.quat_mul(tuple(y[4:]), tuple(R.const(x) for x in (e[0], -e[1], -e[2], -e[3])))
        img = list(ep) + list(qe)
        oku &= all(img[i] == y[mp[i][0]] * mp[i][1] for i in range(8))
    ck.check("EXACT", "(information) S_j is (p, q) -> (e p, q conj(e)) for e = 1, i, j, k", oku)

    def apply_map(poly, mp):
        out = R.zero
        for ex, co in R.terms(poly):
            t = R.const(co)
            for i, e in enumerate(ex):
                if e:
                    t = t * (y[mp[i][0]] * mp[i][1]) ** e
            out = out + t
        return out
    # U[(e, j)][k] = 2 u_k^{(e,j)}(y),   u^{(e,j)}(y) = T^{-1} F^{(e)}(S_j y)
    U = {}
    for e in range(2):
        for j in range(4):
            Fs = [apply_map(F2[e][m], maps[j]) for m in range(8)]
            vec = []
            for k in range(8):
                s = R.zero
                for m in range(8):
                    if Ti[k][m]:
                        s = s + Fs[m] * Ti[k][m]
                vec.append(s)
            U[(e, j)] = vec

    def Gam(k, l):                      # 4 * Gamma_kl,  Gamma_kl = sum_{e,j} u_k u_l
        s = R.zero
        for key in U:
            s = s + U[key][k] * U[key][l]
        return s
    np_ = y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3]
    nq_ = y[4] * y[4] + y[5] * y[5] + y[6] * y[6] + y[7] * y[7]

    def reduce_modN(f):
        """Normal form modulo N = |p|^2 - |q|^2 (monic in p0): p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2.
        f is in the ideal (N) iff the normal form is 0."""
        cs = R.coeff_in(f, 0)
        Rr = nq_ - y[1] * y[1] - y[2] * y[2] - y[3] * y[3]
        out, pw = R.zero, R.one
        for m in range((len(cs) + 1) // 2):
            c0 = cs[2 * m]
            c1 = cs[2 * m + 1] if 2 * m + 1 < len(cs) else R.zero
            out = out + (c0 + c1 * y[0]) * pw
            pw = pw * Rr
        return out
    chk = reduce_modN(np_ - nq_)
    ck.check("EXACT", "reduction modulo N = |p|^2 - |q|^2 is a normal form (N -> 0, p0^3 -> p0(|q|^2 - ...))",
             len(chk) == 0 and reduce_modN(y[0] ** 3) == y[0] * (nq_ - y[1] * y[1] - y[2] * y[2] - y[3] * y[3]))
    M = K.quat_matrix(tuple(y[:4]), tuple(y[4:]), R.one, R.zero, conj_q=False)   # x -> p x q
    MtM = K.ring_matmul([[M[j][i] for j in range(4)] for i in range(4)], M)
    w = K.leibniz_weights(M)
    n2 = np_ * nq_
    ck.check("EXACT", "model x -> p x q (own quaternion multiplication): M^t M = |p|^2|q|^2 I and "
             "det M = (|p|^2|q|^2)^2; so O = M is in SO(4) for unit p, q (and every O in SO(4) arises: "
             "double cover S^3 x S^3 -> SO(4), CITED)",
             all((MtM[i][j] == n2) if i == j else (len(MtM[i][j]) == 0) for i in range(4) for j in range(4))
             and K.det4(M) == n2 * n2)
    Phi = [K.phi(K.F_VALUES[j], M, w) for j in DOC_SRC]          # rebuilt from the definitions
    okmin = True
    for pos in (0, 1, 2, 3, 4, 5, 7):
        okmin &= len(reduce_modN(Phi[pos] * 8 - np_ * Gam(pos, pos))) == 0
    ck.check("EXACT", "Gram identity, minors: 2 Phi_k(M(y)) = |p|^2 Gamma_kk(y) modulo N for the seven minor "
             "generators k (Gamma_kl = sum over 2 parities and 4 maps of u_k u_l, u = T^{-1} F(S_j y))", okmin)
    pb = Gam(6, 6) - Gam(0, 2) * 2 - Gam(2, 3) * 2
    cq = Gam(6, 6) + Gam(2, 4) * 2 - Gam(2, 5) * 2
    ck.check("EXACT", "Gram identity (PB): 2 Phi_70(M(y)) = |p|^2 (Gamma_66 - 2 Gamma_02 - 2 Gamma_23) modulo N "
             "(positions: 0 = f0 (p), 2 = f2 (a), 3 = f12 (b), 6 = f70 (d)): the document's Gram matrix with "
             "G_02 = G_23 = -d", len(reduce_modN(Phi[6] * 8 - np_ * pb)) == 0)
    ck.check("EXACT", "Gram identity (CQ): 2 Phi_70(M(y)) = |p|^2 (Gamma_66 + 2 Gamma_24 - 2 Gamma_25) modulo N "
             "(positions: 4 = f39 (c), 5 = f67 (q)): the document's Gram matrix with G_24 = +d, G_25 = -d",
             len(reduce_modN(Phi[6] * 8 - np_ * cq)) == 0)
    # control: the identity is sensitive (a wrong sign must fail)
    bad = Gam(6, 6) - Gam(0, 2) * 2 + Gam(2, 3) * 2
    ck.check("EXACT", "control: with the sign of G_23 flipped the (PB) identity fails (the check is not vacuous)",
             len(reduce_modN(Phi[6] * 8 - np_ * bad)) != 0)
    # the quadratic forms
    import sympy
    xs = sympy.symbols("x0:8")
    p, r, a, b, c, q, d, s = sympy.symbols("p r a b c q d s", positive=True)
    diag = p * xs[0] ** 2 + r * xs[1] ** 2 + a * xs[2] ** 2 + b * xs[3] ** 2 + c * xs[4] ** 2 \
        + q * xs[5] ** 2 + d * xs[6] ** 2 + s * xs[7] ** 2
    form_pb = diag - 2 * d * xs[0] * xs[2] - 2 * d * xs[2] * xs[3]
    form_cq = diag + 2 * d * xs[2] * xs[4] - 2 * d * xs[2] * xs[5]
    sos_pb = (p * (xs[0] - d * xs[2] / p) ** 2 + b * (xs[3] - d * xs[2] / b) ** 2
              + (a - d ** 2 / p - d ** 2 / b) * xs[2] ** 2
              + r * xs[1] ** 2 + c * xs[4] ** 2 + q * xs[5] ** 2 + d * xs[6] ** 2 + s * xs[7] ** 2)
    sos_cq = (c * (xs[4] + d * xs[2] / c) ** 2 + q * (xs[5] - d * xs[2] / q) ** 2
              + (a - d ** 2 / c - d ** 2 / q) * xs[2] ** 2
              + p * xs[0] ** 2 + r * xs[1] ** 2 + b * xs[3] ** 2 + d * xs[6] ** 2 + s * xs[7] ** 2)
    ck.check("EXACT", "x^t G_PB x = p(x0 - d x2/p)^2 + b(x3 - d x2/b)^2 + (a - d^2/p - d^2/b) x2^2 + (squares "
             "with coefficients r, c, q, d, s): G_PB is positive semidefinite when p, b > 0 and (PB) holds",
             sympy.expand(form_pb - sos_pb) == 0)
    ck.check("EXACT", "x^t G_CQ x = c(x4 + d x2/c)^2 + q(x5 - d x2/q)^2 + (a - d^2/c - d^2/q) x2^2 + (squares "
             "with coefficients p, r, b, d, s): G_CQ is positive semidefinite when c, q > 0 and (CQ) holds",
             sympy.expand(form_cq - sos_cq) == 0)
    ck.info("PROOF-READ: for unit p, q (N = 0, |p| = 1) and g as in (1): 2 Phi_g(O) = sum over the 8 pairs "
            "(parity, map) of u^t G u with G = G_PB or G_CQ (linear in the coefficients; follows from the "
            "identities above), which is >= 0 when G is positive semidefinite.  Pointwise form: "
            "2 Phi_70 >= -2<u_0,u_2> - 2<u_2,u_3> >= -4|x_2|(|x_0| + |x_12|) since |u_k|^2 = 2 x_k^2.")
    k2 = K.find_group_elements((0, 12, 2, 70), (67, 39, 2, 70))
    orb = K.orbit_of_tuple((2, 0, 12, 70))
    unordered = {(t[0], tuple(sorted(t[1:3])), t[3]) for t in orb}
    ck.check("EXACT", "consistency: (leaves f67, f39; centre f2; f70) of (CQ) is the image under g* of "
             "(leaves f0, f12; centre f2; f70) of (PB); the orbit of the chain configuration has 576 elements "
             "(leaves unordered), as stated in PROOF.md 7.6",
             K.G_ELEMENTS.index(GSTAR) in k2 and len(unordered) == 576, "%d" % len(unordered))


# ----------------------------------------------------------------------------------------
# 3. the radial equation (document Section 2)
# ----------------------------------------------------------------------------------------


def check_radial(ck, data):
    ck.section("The necessary radial equation (document p. 2): restricted cubic, exclusion of Q = 0, "
               "Hessian determinant, the polynomial R9")
    R = K.PR(NAMES + ["h"])
    g = R.gens
    p, r, a, b, c, q, d, s, h = g
    coord = [K.F_COORD[j] for j in DOC_SRC]
    ck.check("EXACT", "the lift v of g = p f0 + ... + s f72 to V_R (unique preimage under P) has integer "
             "coordinates in (p, ..., s)", all(x.denominator == 1 for cc in coord for x in cc))
    v = []
    for m in range(14):
        sm = R.zero
        for k in range(8):
            if coord[k][m]:
                sm = sm + g[k] * int(coord[k][m])
        v.append(sm)
    vh = list(v)
    vh[13] = vh[13] + h                                   # v + h kappa
    Cv = K.cubic_C(vh)
    U = d + r + s
    Mm = b + c
    Qd = (a - U) * Mm - d * U
    Br = b * (d + r) * (d + q + s) + c * (d + s) * (d + p + r) + d * (d + p + r) * (d + q + s)
    Cd = b * c * (a - U) - Br
    ck.check("EXACT", "C(v + h kappa) = C + h Q identically, with the document's Q = (a-U)M - dU and "
             "C = bc(a-U) - b(d+r)(d+q+s) - c(d+s)(d+p+r) - d(d+p+r)(d+q+s) (U = d+r+s, M = b+c), for the "
             "paper's cubic (A.2) with its sign", Cv == Cd + h * Qd)
    ck.check("EXACT", "Q = v_0 v_6 - v_1 v_5 = (1/2) v^t Hess C(kappa) v", Qd == v[0] * v[6] - v[1] * v[5])
    # Q = 0 exclusion
    lhs = Mm * Cd - b * c * Qd
    ck.check("EXACT", "M C - bc Q = bc dU - M [b(d+r)(d+q+s) + c(d+s)(d+p+r) + d(d+p+r)(d+q+s)] identically",
             lhs == b * c * d * U - Mm * Br)
    pos = lhs * (-1)
    td = R.to_dict(pos)
    ck.check("EXACT", "-(M C - bc Q) is a polynomial with nonnegative coefficients (%d terms) containing "
             "b d^3 and c d^3: so Q = 0, M > 0, d > 0 force M C < 0, i.e. C != 0" % len(td),
             all(x > 0 for x in td.values())
             and td.get((0, 0, 0, 1, 0, 0, 3, 0, 0)) == 1 and td.get((0, 0, 0, 0, 1, 0, 3, 0, 0)) == 1)
    Q00 = R.from_terms([(e, co) for e, co in R.terms(Qd) if e[3] == 0 and e[4] == 0])
    ck.check("EXACT", "at b = c = 0: Q = -d(d + r + s) < 0 for d > 0", Q00 == d * U * (-1))
    ck.info("PROOF-READ: for x in X with c(x) = g and d > 0, u = Re x = v + H kappa (H = Re F_pi*(x)) "
            "satisfies C(u) = 0, so C + HQ = 0; by the two checks above Q != 0, hence H = -C/Q")
    # Hessian determinant
    Hm = K.hess_at(vh, R.zero)
    t0 = time.time()
    det = K.det_laplace(Hm, R.zero, R.one)
    Dk = R.coeff_in(det, 8)
    ck.check("EXACT", "det Hess C(v + h kappa): computed by division-free Laplace expansion (14 x 14, entries "
             "linear in 9 variables): %d terms, homogeneous of degree 14, degree %d in h"
             % (len(det), len(Dk) - 1), len(det) > 0 and R.total_degree_set(det) == [14] and len(Dk) - 1 <= 4,
             "%.1f s" % (time.time() - t0))
    # control: evaluate the determinant at a random integer point by Fraction elimination
    rng = random.Random(3434)
    okdet = True
    for _ in range(3):
        pt = [rng.randint(1, 7) for _ in range(9)]
        val = 0
        for e, co in R.terms(det):
            t = co
            for x, k in zip(pt, e):
                t *= x ** k
            val += t
        vnum = [sum(int(coord[k][m]) * pt[k] for k in range(8)) for m in range(14)]
        vnum[13] += pt[8]
        Hn = K.hess_at([Fraction(x) for x in vnum], Fraction(0))
        A = [row[:] for row in Hn]
        dv = Fraction(1)
        for i in range(14):
            piv = next((k for k in range(i, 14) if A[k][i] != 0), None)
            if piv is None:
                dv = Fraction(0)
                break
            if piv != i:
                A[i], A[piv] = A[piv], A[i]
                dv = -dv
            dv *= A[i][i]
            for k in range(i + 1, 14):
                if A[k][i] != 0:
                    f = A[k][i] / A[i][i]
                    A[k] = [x - f * z for x, z in zip(A[k], A[i])]
        okdet &= dv == val
    ck.check("EXACT", "control: the determinant polynomial agrees with Gaussian elimination over Q at 3 random "
             "integer points (independent of the Laplace code)", okdet)
    Sx = R.zero
    for k, Dc in enumerate(Dk):
        Sx = Sx + Dc * (Cd * (-1)) ** k * Qd ** (4 - k)      # = Q^4 det Hess C(v - (C/Q) kappa)
    neg = Sx * (-1)
    ok4 = all(x % 4 == 0 for x in R.coeffs(neg))
    quarter = R.from_terms([(e, co // 4) for e, co in R.terms(neg)]) if ok4 else R.zero
    t0 = time.time()
    try:
        root = R.exact_sqrt(quarter)
        R9 = R.exact_div(root, Qd)
    except ValueError as ex:
        ck.check("EXACT", "-(1/4) Q^4 det Hess C(v - (C/Q) kappa) is the square of Q times a polynomial", False,
                 str(ex))
        return None
    if R.to_dict(R9).get((0, 0, 0, 0, 0, 0, 9, 0, 0), 0) < 0:
        R9 = R9 * (-1)
    ck.check("EXACT", "Hessian identity: Q^4 det Hess C(v - (C/Q) kappa) = -4 Q^2 R9^2 (the document's "
             "Q^3 det = -4 Q R9^2), with R9 computed here as an exact square root (%d terms) and exact "
             "quotient by Q" % len(root), ok4 and Sx == (Qd * Qd * R9 * R9) * (-4),
             "%.1f s" % (time.time() - t0))
    ck.check("EXACT", "R9 (own): %d terms, homogeneous of degree 9, free of h, degree 2 in a, coefficient of "
             "d^9 equal to +1" % len(R9),
             len(R9) == 1381 and R.total_degree_set(R9) == [9] and R.degrees(R9)[8] == 0
             and R.degrees(R9)[2] == 2 and R.to_dict(R9)[(0, 0, 0, 0, 0, 0, 9, 0, 0)] == 1)
    rad = data["r9"]
    R9s = R.from_terms([(tuple(e), int(co)) for e, co in rad["factors"][0]["terms"]])
    Cs = R.from_terms([(tuple(e), int(co)) for e, co in rad["C_terms"]])
    Qs = R.from_terms([(tuple(e), int(co)) for e, co in rad["Q_terms"]])
    ck.check("EXACT", "comparison with the supplement's data (variables %s): its stored R9 equals the R9 "
             "computed here (same sign normalisation), its stored C, Q equal the document's, its stored "
             "lift equals the lift computed here" % rad["variables"],
             rad["variables"] == NAMES + ["h"] and R9s == R9 and Cs == Cd and Qs == Qd
             and [[int(x) for x in row] for row in rad["costlift"]]
             == [[int(coord[k][m]) for k in range(8)] for m in range(14)])
    ck.info("PROOF-READ: det Hess C(u) = 0 at u = v - (C/Q) kappa (radial necessity; u != 0 since g != 0), "
            "so R9(p, r, a, b, c, q, d, s) = 0 for every x in X with c(x) = g in N_T and d > 0")
    return R, R9, Cd, Qd


# ----------------------------------------------------------------------------------------
# 4. the polynomial implication (document Section 3, Lemma 2)
# ----------------------------------------------------------------------------------------


def to_sympy(Rg, poly, syms):
    import sympy
    tot = sympy.Integer(0)
    for e, co in Rg.terms(poly):
        t = sympy.Integer(co)
        for sy, k in zip(syms, e):
            t *= sy ** k
        tot += t
    return tot


def chart_data(Rt, case):
    """Numerators and denominators (polynomials of Rt = Z[M, T, u, W, s, r]) of z, M - z, P and the
    other primary coefficient X in the three charts of the document (theta = T/(1+T),
    X = P W/(1+W))."""
    M, T, u, W, ss, rr = Rt.gens
    one = Rt.one
    if case == "III":
        zn, zd = M * T, (T + 1) * (M + 1)
        on, od = M * ((M + 1) * (T + 1) - T), zd
        Pn, Pd = u, one
    else:
        zn, zd = M * (T + 1 + M * T), (T + 1) * (M + 1)
        on, od = M * M, zd
        P0n, P0d = T + 1 + M * T, (M + 1) * T
        if case == "I":
            Pn, Pd = P0n + u * P0d, P0d
        else:
            Pn, Pd = P0n * u, P0d * (u + 1)
    Xn, Xd = Pn * W, Pd * (W + 1)
    return zn, zd, on, od, Pn, Pd, Xn, Xd


def check_lemma2(ck, data, Rrad, R9):
    ck.section("The polynomial implication (document pp. 2-3, Lemma 2): R9 = 0, d > 0  =>  (H) or (PB) or (CQ)")
    # d = 1; ring order (p, b, c, q, s, r)
    Rk = K.PR(["p", "b", "c", "q", "s", "r"])
    p, b, c, q, s, r = Rk.gens
    buckets = [[], [], []]
    for e, co in Rrad.terms(R9):
        buckets[e[2]].append(((e[0], e[3], e[4], e[5], e[7], e[1]), co))
    k = [Rk.from_terms(t) for t in buckets]
    ck.check("EXACT", "d = 1 (homogeneity): R9 = k2 a^2 + k1 a + k0 with k0, k1, k2 in Z[p, r, b, c, q, s] "
             "(%d, %d, %d terms)" % (len(k[0]), len(k[1]), len(k[2])), len(k[0]) + len(k[1]) + len(k[2]) == 1381)
    ck.check("EXACT", "every coefficient of k0 is positive and its constant term is 1 (so k0 >= 1)",
             all(x > 0 for x in Rk.coeffs(k[0])) and Rk.to_dict(k[0]).get((0,) * 6) == 1)
    k1_00 = [(e, co) for e, co in Rk.terms(k[1]) if e[1] == 0 and e[2] == 0]
    ck.check("EXACT", "at b = c = 0: k2 = 0 and every coefficient of k1 is positive (%d terms): R9(a) > 0 "
             "for a >= 0, so no realizable g has b = c = 0 < d" % len(k1_00),
             all(e[1] or e[2] for e, co in Rk.terms(k[2])) and len(k1_00) > 0
             and all(co > 0 for e, co in k1_00))

    def swap(f):
        return Rk.from_terms([((e[3], e[2], e[1], e[0], e[5], e[4]), co) for e, co in Rk.terms(f)])
    ck.check("EXACT", "k0, k1, k2 are invariant under (p, b, r) <-> (q, c, s) (the action of g*): the two "
             "orders p >= q and q >= p are mirror images (both are expanded below anyway)",
             all(swap(x) == x for x in k))
    # Bernstein form of a quadratic
    import sympy
    K0, K1, K2, A, t = sympy.symbols("k0 k1 k2 A t")
    b0, b1, b2 = K0, K0 + K1 * A / 2, K0 + K1 * A + K2 * A ** 2
    ck.check("EXACT", "Bernstein form: k0 + k1 (tA) + k2 (tA)^2 = b0 (1-t)^2 + 2 b1 t(1-t) + b2 t^2 with "
             "b0 = k0, b1 = k0 + k1 A/2, b2 = k0 + k1 A + k2 A^2; for A = n/e: b1 = (2 k0 e + k1 n)/(2e), "
             "b2 = (k0 e^2 + k1 n e + k2 n^2)/e^2",
             sympy.expand(K0 + K1 * t * A + K2 * (t * A) ** 2
                          - (b0 * (1 - t) ** 2 + 2 * b1 * t * (1 - t) + b2 * t ** 2)) == 0)
    ck.info("PROOF-READ: if b1, b2 >= 0 and b0 = k0 > 0 then R9(a) > 0 for 0 <= a < A (the term b0 (1-t)^2 "
            "is positive for t < 1)")
    # chart identities
    Rt = K.PR(["M", "T", "u", "W", "s", "r"])
    Ms, Ts, us, Ws, ss_, rs_ = sympy.symbols("M T u W s r", positive=True)
    syms = (Ms, Ts, us, Ws, ss_, rs_)
    th = Ts / (1 + Ts)
    okc = True
    # charts I, II
    zn, zd, on, od, Pn, Pd, Xn, Xd = chart_data(Rt, "I")
    z1 = to_sympy(Rt, zn, syms) / to_sympy(Rt, zd, syms)
    o1 = to_sympy(Rt, on, syms) / to_sympy(Rt, od, syms)
    P1 = to_sympy(Rt, Pn, syms) / to_sympy(Rt, Pd, syms)
    P0 = (1 + Ms * th) / ((1 + Ms) * th)
    okc &= sympy.cancel(z1 - Ms * (1 + Ms * th) / (1 + Ms)) == 0
    okc &= sympy.cancel(o1 - (Ms - z1)) == 0 and sympy.cancel(o1 - Ms ** 2 * (1 - th) / (1 + Ms)) == 0
    okc &= sympy.cancel(z1 * (Ms + 1) - Ms - Ms ** 2 * th) == 0                 # V0 = M^2 theta > 0
    okc &= sympy.cancel(1 / P0 + 1 / z1 - (1 + 1 / Ms)) == 0                    # A1 = A0 at P = P0
    okc &= sympy.cancel(P1 - (P0 + us)) == 0
    zn, zd, on, od, Pn, Pd, Xn, Xd = chart_data(Rt, "II")
    P2 = to_sympy(Rt, Pn, syms) / to_sympy(Rt, Pd, syms)
    okc &= sympy.cancel(P2 - P0 * us / (1 + us)) == 0
    X2 = to_sympy(Rt, Xn, syms) / to_sympy(Rt, Xd, syms)
    okc &= sympy.cancel(X2 - P2 * Ws / (1 + Ws)) == 0
    zn, zd, on, od, Pn, Pd, Xn, Xd = chart_data(Rt, "III")
    z3 = to_sympy(Rt, zn, syms) / to_sympy(Rt, zd, syms)
    o3 = to_sympy(Rt, on, syms) / to_sympy(Rt, od, syms)
    okc &= sympy.cancel(z3 - Ms * th / (1 + Ms)) == 0 and sympy.cancel(o3 - (Ms - z3)) == 0
    okc &= sympy.cancel(z3 * (Ms + 1) - Ms - Ms * (th - 1)) == 0                # V0 = M(theta - 1) <= 0
    okc &= sympy.cancel(o3 - Ms * (1 + Ms - th) / (1 + Ms)) == 0
    ck.check("EXACT", "chart identities (theta = T/(1+T)): I, II: z = M(1+M theta)/(1+M), M - z = "
             "M^2(1-theta)/(1+M), V0 = z(M+1) - M = M^2 theta, P0 = (1+M theta)/((1+M) theta) solves "
             "1/P0 + 1/z = 1 + 1/M, P = P0 + u (I), P = P0 u/(1+u) (II); III: z = M theta/(1+M), "
             "M - z = M(1+M-theta)/(1+M), V0 = M(theta - 1); other primary = P W/(1+W)", okc)
    ck.info("PROOF-READ (coverage, order p >= q; the other order is the mirror image): put P = p, z = b, "
            "X = q in [0, P], c = M - z, A0 = 1 + 1/M, A1 = 1/P + 1/z (+infinity if Pz = 0), V0 = z(M+1) - M. "
            "V0 > 0 <=> 1/z < A0, and then A1 <= A0 <=> P >= P0.  Region I = {V0 > 0, P >= P0} uses A = A1; "
            "region II = {V0 > 0, 0 <= P <= P0} and region III = {V0 <= 0} use A = A0; in each region "
            "A = min(A0, A1) and A is finite and continuous.  The three charts with M > 0, T > 0 (T >= 0 in "
            "III), u, W >= 0 cover the regions except z = M (I, II), z = M/(M+1) (III), X = P > 0 and, for "
            "II, P = P0; these lie in the closure of the chart images inside the same region, and b1, b2 are "
            "continuous there, so b1, b2 >= 0 extend to the closed regions.")
    # the twelve Bernstein numerators
    rec = {(x["order"], x["case"], x["coefficient"]): x for x in data["bern"]["Bernstein_tables"]}
    pdf_counts = {("I", "middle"): 24805, ("I", "endpoint"): 37747, ("II", "middle"): 38066,
                  ("II", "endpoint"): 52638, ("III", "middle"): 10902, ("III", "endpoint"): 16934}
    names = {"I": "chain_min", "II": "Hmin_positive_V", "III": "Hmin_nonpositive_V"}
    ss, rr = Rt.gens[4], Rt.gens[5]
    one = Rt.one
    spot_rng = random.Random(3402)
    spot_ok = True
    for order in ("p>=q", "q>=p"):
        primary, endpoint = (p, b) if order == "p>=q" else (q, c)
        for case in ("I", "II", "III"):
            zn, zd, on, od, Pn, Pd, Xn, Xd = chart_data(Rt, case)
            if order == "p>=q":           # ring order p, b, c, q, s, r :  p = P, b = z, c = M - z, q = X
                nums, dens = [Pn, zn, on, Xn, ss, rr], [Pd, zd, od, Xd, one, one]
            else:                         # q = P, c = z, b = M - z, p = X
                nums, dens = [Xn, on, zn, Pn, ss, rr], [Xd, od, zd, Pd, one, one]
            if case == "I":
                n_, e_ = primary + endpoint, primary * endpoint       # A1 = 1/P + 1/z = n/e
            else:
                n_, e_ = b + c + 1, b + c                             # A0 = 1 + 1/M = n/e
            mid = k[0] * e_ * 2 + k[1] * n_
            end = k[0] * e_ * e_ + k[1] * n_ * e_ + k[2] * n_ * n_
            for label, poly in (("middle", mid), ("endpoint", end)):
                t0 = time.time()
                out = Rk.subs_clear(poly, nums, dens)
                co = Rt.coeffs(out)
                x = rec[(order, names[case], label)]
                # independent control of the substitution: exact evaluation at a random rational point
                pt = [Fraction(spot_rng.randint(1, 7), spot_rng.randint(1, 3)) for _ in range(6)]
                nv = [eval_terms(Rt.terms(n1), pt) for n1 in nums]
                dv = [eval_terms(Rt.terms(d1), pt) for d1 in dens]
                degs = Rk.degrees(poly)
                mult = Fraction(1)
                for dv1, dg in zip(dv, degs):
                    mult *= dv1 ** dg
                spot_ok &= (eval_terms(Rt.terms(out), pt)
                            == mult * eval_terms(Rk.terms(poly), [a1 / b1 for a1, b1 in zip(nv, dv)]))
                spot_ok &= all(v > 0 for v in dv) and mult > 0
                ck.check("EXACT", "Bernstein certificate, order %s, chart %s, %s numerator: after the "
                         "substitution and multiplication by prod_i den_i^{deg_i} (a positive polynomial) "
                         "all %d coefficients are positive integers (min %d); the count equals the document's "
                         "table (%d) and the canonical digest equals the one stored in the supplement's data"
                         % (order, case, label, len(co), min(co) if co else 0, pdf_counts[(case, label)]),
                         len(co) > 0 and all(v > 0 for v in co) and len(co) == pdf_counts[(case, label)]
                         and len(co) == x["terms"] and Rt.digest(out) == x["sha256"],
                         "%.1f s" % (time.time() - t0))
    ck.check("EXACT", "control of the twelve substitutions: at one random positive rational point each "
             "(seed 3402) the expanded numerator equals prod_i den_i^{deg_i} times the Bernstein numerator "
             "evaluated at the chart point, and the multiplier is positive (exact rational evaluation, "
             "independent of the expansion code and of the supplement's digests)", spot_ok)
    return Rk, k


def eval_terms(terms, pt):
    """Exact value of sum coef * prod pt_i^e_i for a term list [(exponents, coef)]."""
    tot = Fraction(0)
    for e, co in terms:
        t = Fraction(co)
        for x, kk in zip(pt, e):
            if kk:
                t *= x ** kk
        tot += t
    return tot


def exact_boundary_spot_checks(ck, Rk, k, seed=3401, samples=40):
    """Exact evaluation of the Bernstein coefficients b1, b2 at rational points ON the boundaries of the
    three closed regions (which the charts reach only in the limit) and at interior points."""
    rng = random.Random(seed)
    terms = [Rk.terms(x) for x in k]

    def ev(tl, pt):
        tot = Fraction(0)
        for e, co in tl:
            t = Fraction(co)
            for x, kk in zip(pt, e):
                if kk:
                    t *= x ** kk
            tot += t
        return tot

    def rq(lo=1, hi=9):
        return Fraction(rng.randint(lo, hi), rng.randint(1, 6))
    worst = None
    n = 0
    ok = True
    kinds = ["I:z=M", "I:X=P", "I:P=P0", "I:z=M,X=P", "II:P=0", "II:z=M", "II:X=P", "III:z=0", "III:V0=0",
             "III:P=0", "III:X=P", "III:z=0,P=0", "I:int", "II:int", "III:int"]
    for order in ("p>=q", "q>=p"):
        for kind in kinds:
            for _ in range(samples):
                M = rq()
                reg, what = kind.split(":")
                th = Fraction(rng.randint(1, 9), 10)
                if reg in ("I", "II"):
                    if "z=M" in what:
                        th = Fraction(1)
                    z = M * (1 + M * th) / (1 + M)
                    P0 = (1 + M * th) / ((1 + M) * th)
                    if reg == "I":
                        P = P0 if "P=P0" in what else P0 + rq()
                    else:
                        P = Fraction(0) if "P=0" in what else P0 * Fraction(rng.randint(1, 9), 10)
                else:
                    if "V0=0" in what:
                        th = Fraction(1)
                    if "z=0" in what:
                        th = Fraction(0)
                    z = M * th / (1 + M)
                    P = Fraction(0) if "P=0" in what else rq()
                X = P if "X=P" in what else P * Fraction(rng.randint(0, 9), 10)
                rr = rq(0, 5) if rng.random() < 0.7 else Fraction(0)
                sv = rq(0, 5) if rng.random() < 0.7 else Fraction(0)
                A0 = 1 + 1 / M
                A1 = (1 / P + 1 / z) if (P > 0 and z > 0) else None
                A = A0 if (A1 is None or A0 <= A1) else A1
                # consistency of the region with the threshold selection
                V0 = z * (M + 1) - M
                if reg == "I":
                    ok &= V0 > 0 and A1 is not None and A1 <= A0 and A == A1
                else:
                    ok &= (A == A0) and ((V0 > 0) == (reg == "II"))
                if order == "p>=q":
                    pt = (P, z, M - z, X, sv, rr)          # p, b, c, q, s, r
                else:
                    pt = (X, M - z, z, P, sv, rr)
                v0, v1, v2 = ev(terms[0], pt), ev(terms[1], pt), ev(terms[2], pt)
                b1 = v0 + v1 * A / 2
                b2 = v0 + v1 * A + v2 * A * A
                tt = Fraction(rng.randint(0, 99), 100)
                val = v0 + v1 * (tt * A) + v2 * (tt * A) ** 2
                ok &= v0 >= 1 and b1 >= 0 and b2 >= 0 and val > 0
                m_ = min(b1, b2)
                worst = m_ if worst is None or m_ < worst else worst
                n += 1
    ck.check("EXACT", "closed regions: at %d rational points (seed %d; both orders; boundaries z = M, X = P, "
             "P = P0, P = 0, z = 0, V0 = 0 and combinations, r or s = 0 allowed, plus interior points) the "
             "selected threshold is min(A0, A1), b1 >= 0, b2 >= 0, k0 >= 1 and R9(tA) > 0 for a random "
             "t in [0, 1) (exact rational evaluation of own k0, k1, k2)" % (n, seed), ok,
             "smallest Bernstein coefficient met: %s" % (float(worst) if worst is not None else None))


# ----------------------------------------------------------------------------------------
# 5. NUMERIC falsification attempts
# ----------------------------------------------------------------------------------------


def numeric_part(ck, Rrad, R9, Cd, Qd):
    import numpy as np
    from scipy.optimize import minimize, least_squares
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    S = np.array(K.PERMS)
    sg = np.array(K.SGN, float)
    FV = {j: np.array(K.F_VALUES[j], float) for j in DOC_SRC}

    def so4_from_quats(pq):
        """rows of pq: 8 numbers; returns (n, 4, 4) rotation matrices x -> p x q (normalised)."""
        p = pq[:, :4] / np.linalg.norm(pq[:, :4], axis=1, keepdims=True)
        q = pq[:, 4:] / np.linalg.norm(pq[:, 4:], axis=1, keepdims=True)
        cols = []
        for kk in range(4):
            e = np.zeros((len(pq), 4))
            e[:, kk] = 1.0
            pe = np.stack(K.quat_mul(tuple(p.T), tuple(e.T)), axis=1)
            cols.append(np.stack(K.quat_mul(tuple(pe.T), tuple(q.T)), axis=1))
        return np.stack(cols, axis=2)

    def weights(O):
        return sg[None, :] * O[:, 0, S[:, 0]] * O[:, 1, S[:, 1]] * O[:, 2, S[:, 2]] * O[:, 3, S[:, 3]]

    def minor(O, j):
        I, J = K.minor_spec(j)
        if len(I) == 1:
            return O[:, I[0], J[0]]
        return O[:, I[0], J[0]] * O[:, I[1], J[1]] - O[:, I[0], J[1]] * O[:, I[1], J[0]]

    def pointwise(O):
        w = weights(O)
        ph70 = w @ FV[70]
        x = {j: np.abs(minor(O, j)) for j in (0, 2, 12, 39, 67)}
        return {
            "old-H (f12; f2; f70)": ph70 + x[2] ** 2 + 2 * x[12] * x[2],
            "old-H (f39; f2; f70)": ph70 + x[2] ** 2 + 2 * x[39] * x[2],
            "chain (f0, f12; f2; f70) = (PB)": ph70 + 2 * x[2] * (x[0] + x[12]),
            "chain (f67, f39; f2; f70) = (CQ)": ph70 + 2 * x[2] * (x[67] + x[39]),
        }, w
    rng = np.random.default_rng(20261134)
    nsamp = 300000
    pq = rng.normal(size=(nsamp, 8))
    O = so4_from_quats(pq)
    orth = float(np.max(np.abs(np.einsum("nij,nkj->nik", O, O) - np.eye(4)[None])))
    vals, w = pointwise(O)
    okm = float(np.max(np.abs((w @ FV[2]) - minor(O, 2) ** 2))) < 1e-12
    ck.check("NUMERIC", "sampler: %d matrices from unit quaternion pairs are orthogonal (max defect %.1e) "
             "and Phi_2 = x_2^2 (tolerance 1e-12)" % (nsamp, orth), orth < 1e-12 and okm)
    for name, v in vals.items():
        def fobj(x, name=name):
            return float(pointwise(so4_from_quats(x[None, :]))[0][name][0])
        best = float(v.min())
        idx = np.argsort(v)[:40]
        for i0 in idx:
            res = minimize(fobj, pq[i0], method="Nelder-Mead",
                           options={"xatol": 1e-12, "fatol": 1e-15, "maxiter": 4000})
            best = min(best, float(res.fun))
        ck.check("NUMERIC", "pointwise %s >= 0 on SO(4): %d random samples (seed 20261134), then Nelder-Mead "
                 "from the 40 smallest; minimum found %.3e (tolerance -1e-10)" % (name, nsamp, best),
                 best > -1e-10)
    # coefficient forms at the thresholds, zero coefficients included
    worst = {"H": np.inf, "PB": np.inf, "CQ": np.inf}
    polished = {"H": np.inf, "PB": np.inf, "CQ": np.inf}
    rr = np.random.default_rng(20261135)
    nsub = 40000
    wsub = w[:nsub]
    PH = {j: wsub @ FV[j] for j in DOC_SRC}
    ntrial, npolish = 300, 12
    for trial in range(ntrial):
        co = 10.0 ** rr.uniform(-2, 2, size=8)
        co[rr.random(8) < 0.25] = 0.0
        p, r, a, b, c, q, d, s = co
        d = d if d > 0 else 1.0
        for cond in ("H", "PB", "CQ"):
            cc = dict(p=p, r=r, b=b, c=c, q=q, d=d, s=s)
            if cond == "H":
                if cc["b"] + cc["c"] == 0:
                    cc["b"] = 1.0
                cc["a"] = d + d * d / (cc["b"] + cc["c"])
            elif cond == "PB":
                if cc["p"] == 0:
                    cc["p"] = 0.5
                if cc["b"] == 0:
                    cc["b"] = 2.0
                cc["a"] = d * d / cc["p"] + d * d / cc["b"]
            else:
                if cc["c"] == 0:
                    cc["c"] = 0.5
                if cc["q"] == 0:
                    cc["q"] = 2.0
                cc["a"] = d * d / cc["c"] + d * d / cc["q"]
            scale = sum(abs(cc[nm]) for nm in NAMES)
            tot = sum(cc[nm] * PH[j] for nm, j in zip(NAMES, DOC_SRC)) / scale
            worst[cond] = min(worst[cond], float(tot.min()))
            if trial < npolish:
                def fobj2(x, cc=cc, scale=scale):
                    ww = weights(so4_from_quats(x[None, :]))[0]
                    return float(sum(cc[nm] * (ww @ FV[j]) for nm, j in zip(NAMES, DOC_SRC)) / scale)
                for i0 in np.argsort(tot)[:2]:
                    res = minimize(fobj2, pq[i0], method="Nelder-Mead",
                                   options={"xatol": 1e-12, "fatol": 1e-15, "maxiter": 3000})
                    polished[cond] = min(polished[cond], float(res.fun))
    ck.check("NUMERIC", "coefficient forms with equality in (H), resp. (PB), resp. (CQ): %d random coefficient "
             "vectors each (seed 20261135; log-uniform in [1e-2, 1e2], each coefficient set to 0 with "
             "probability 1/4 where the condition allows) on %d rotations: smallest Phi_g / (sum of "
             "coefficients) = %.2e, %.2e, %.2e; after Nelder-Mead from the 2 smallest samples of the first %d "
             "vectors: %.2e, %.2e, %.2e (tolerance -1e-10)"
             % (ntrial, nsub, worst["H"], worst["PB"], worst["CQ"], npolish,
                polished["H"], polished["PB"], polished["CQ"]),
             min(min(worst.values()), min(polished.values())) > -1e-10)
    # Lemma 2 by sampling
    terms = Rrad.terms(R9)
    E = np.array([e[:8] for e, co in terms], dtype=np.int64)
    Cf = np.array([float(co) for e, co in terms])

    def R9f(wv):
        wv = np.asarray(wv, dtype=np.longdouble)
        out = np.zeros(wv.shape[0], dtype=np.longdouble)
        mag = np.zeros(wv.shape[0], dtype=np.longdouble)
        # powers table
        pw = [np.stack([wv[:, i] ** kk for kk in range(10)], axis=1) for i in range(8)]
        for e, co in zip(E, Cf):
            t = np.ones(wv.shape[0], dtype=np.longdouble)
            for i in range(8):
                if e[i]:
                    t = t * pw[i][:, e[i]]
            out += co * t
            mag += abs(co) * np.abs(t)
        return out, mag
    r3 = np.random.default_rng(20261136)
    n3 = 60000
    wv = 10.0 ** r3.uniform(-2.5, 2.5, size=(n3, 8))
    wv[r3.random((n3, 8)) < 0.12] = 0.0
    wv[:, 6] = 1.0                                           # d = 1
    p, r, a, b, c, q, d, s = [wv[:, i] for i in range(8)]
    M = b + c
    keep = M > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        A0 = 1 + 1 / M
        Pm = np.where(p >= q, p, q)
        zm = np.where(p >= q, b, c)
        A1 = np.where((Pm > 0) & (zm > 0), 1 / Pm + 1 / zm, np.inf)
        Amin = np.minimum(A0, A1)
    tt = r3.random(n3) * 0.999999
    wv[:, 2] = np.where(keep, tt * Amin, r3.random(n3) * 5)   # a < min(A0, A1), or a >= 0 if b = c = 0
    val, mag = R9f(wv)
    rel = val / mag
    ck.check("NUMERIC", "Lemma 2: %d random points (seed 20261136; d = 1, other coefficients log-uniform in "
             "[10^-2.5, 10^2.5] or 0 with probability 0.12; a = t min(A0, A1), t in [0, 1), or a in [0, 5) when "
             "b = c = 0): R9 > 0 at all of them; smallest R9 / (sum of |terms|) = %.2e" % (n3, float(rel.min())),
             bool(np.all(val > 0)))
    # the roots of the quadratic R9(a) directly
    kq = [np.zeros(n3, dtype=np.longdouble) for _ in range(3)]
    wl = np.asarray(wv, dtype=np.longdouble)
    pw = [np.stack([wl[:, i] ** kk for kk in range(10)], axis=1) for i in range(8)]
    for e, co in zip(E, Cf):
        tm = np.ones(n3, dtype=np.longdouble)
        for i in range(8):
            if i != 2 and e[i]:
                tm = tm * pw[i][:, e[i]]
        kq[e[2]] += co * tm
    k0f, k1f, k2f = [np.asarray(x, dtype=float) for x in kq]
    disc = k1f * k1f - 4 * k2f * k0f
    with np.errstate(divide="ignore", invalid="ignore"):
        sq = np.sqrt(np.where(disc >= 0, disc, np.nan))
        # numerically stable roots of k2 a^2 + k1 a + k0
        qq = -0.5 * (k1f + np.where(k1f >= 0, 1.0, -1.0) * sq)
        ra = np.where(k2f != 0, qq / k2f, np.nan)
        rb = np.where(qq != 0, k0f / qq, np.nan)
        lin = np.where((k2f == 0) & (k1f != 0), -k0f / k1f, np.nan)
        cand = np.stack([ra, rb, lin], axis=1)
        cand = np.where(cand > 0, cand, np.inf)
        small = np.nanmin(np.where(np.isnan(cand), np.inf, cand), axis=1)
        ratio = small / Amin
    has = np.isfinite(small) & keep
    ck.check("NUMERIC", "Lemma 2 through the roots: on the same coefficient vectors, %d have a positive root a "
             "of R9 (d = 1); the smallest positive root divided by min(A0, A1) is >= %.6f at all of them "
             "(the lemma asserts >= 1; tolerance 1e-9); no positive root when b = c = 0 (%d vectors)"
             % (int(has.sum()), float(ratio[has].min()) if has.any() else float("nan"), int((~keep).sum())),
             has.sum() > 1000 and bool(np.all(ratio[has] >= 1 - 1e-9))
             and not bool(np.any(np.isfinite(small[~keep]))))
    # Q = 0 exclusion numerically
    U = d + r + s
    with np.errstate(divide="ignore", invalid="ignore"):
        a0 = U + d * U / M
        Cq = b * c * (a0 - U) - b * (d + r) * (d + q + s) - c * (d + s) * (d + p + r) \
            - d * (d + p + r) * (d + q + s)
    ck.check("NUMERIC", "Q = 0 exclusion on the same vectors with M > 0 (a chosen so that Q = 0): C < 0 at all "
             "%d of them" % int(keep.sum()), bool(np.all(Cq[keep] < 0)))
    # points of X with the five equalities only: C + HQ = 0 and R9 = 0
    T = list(T34)
    others = [x for x in range(24) if x not in T]
    excl = [18, 16, 13, 8, 4, 14, 11, 7]

    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])

    def gaps(x):
        a_, b_, al = unpack(x)
        F = (al * np.prod(a_[None, :] + b_[S], axis=1)).real
        return F, F - F[T[0]]
    r4 = np.random.default_rng(20261137)
    pts = []
    for _ in range(120):
        x0 = r4.normal(size=17)
        try:
            sol = least_squares(lambda x: np.r_[gaps(x)[1][T[1:]], np.sum(gaps(x)[1] ** 2) - 1.0,
                                                np.zeros(11)], x0,
                                method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception:
            continue
        F, g = gaps(sol.x)
        if np.max(np.abs(g[T])) < 1e-12:
            pts.append((F[T[0]], g))
    chq = Rrad.terms(Cd + Rrad.gens[8] * Qd)            # C + h Q in the variables (p, ..., s, h)
    worst_r = worst_c = 0.0
    nd = nneg = ntriv = 0
    for H, g in pts:
        wv1 = np.array([[g[x] for x in excl]])
        val, mag = R9f(wv1)
        pt9 = list(wv1[0]) + [H]
        tot = mg = 0.0
        for e, co in chq:
            tm = float(co)
            for xx, kk in zip(pt9, e):
                if kk:
                    tm *= xx ** kk
            tot += tm
            mg += abs(tm)
        if mg > 0:
            worst_c = max(worst_c, abs(tot) / mg)
        if float(mag[0]) > 0:
            worst_r = max(worst_r, float(abs(val[0]) / mag[0]))
        else:
            ntriv += 1                    # every monomial of R9 vanishes (several coefficients are 0)
        nd += int(abs(wv1[0][6]) > 1e-3)
        nneg += int(np.min(wv1[0]) < -1e-6)
    ck.check("NUMERIC", "Section 2 on actual points of X: %d points x = alpha Psi(a, b) (Levenberg-Marquardt "
             "from 120 random starts, seed 20261137) whose gaps vanish on T, WITHOUT sign conditions on the "
             "other gaps (the identities of Section 2 are polynomial, so C + HQ = 0 and, where Q != 0, R9 = 0 "
             "must hold there for coefficients of either sign; %d of the points have a negative coefficient, "
             "%d have |d| > 1e-3 in the normalisation sum of squared gaps = 1, %d are degenerate points at "
             "which every monomial of R9 vanishes): max |C + HQ| / (sum of |terms|) = %.1e, "
             "max |R9| / (sum of |terms|) = %.1e over the other points (tolerance 1e-9)"
             % (len(pts), nneg, nd, ntriv, worst_c, worst_r),
             len(pts) >= 30 and nd >= 10 and worst_c < 1e-9 and worst_r < 1e-9)
    # genuine points of the stratum: all gaps outside T at least `margin` times their mean
    margin = 0.01
    r5 = np.random.default_rng(20261138)
    feas = []

    def res(x):
        F, g = gaps(x)
        return np.r_[g[T[1:]], np.minimum(0.0, g[others] - margin), np.mean(g[others]) - 1.0]
    nstart = 300
    for _ in range(nstart):
        x0 = r5.normal(size=17)
        try:
            sol = least_squares(res, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        if np.all(np.isfinite(sol.fun)) and np.max(np.abs(sol.fun)) < 1e-11:
            feas.append(gaps(sol.x))
    okf = True
    worst_r = worst_c = 0.0
    slack = []
    dvals = []
    nall3 = 0
    for F, g in feas:
        wv1 = np.array([[g[x] for x in excl]])
        p, r, a, b, c, q, d, s = wv1[0]
        okf &= bool(np.min(g[others]) >= margin - 1e-9) and d > 0
        val, mag = R9f(wv1)
        worst_r = max(worst_r, float(abs(val[0]) / mag[0]))
        pt9 = list(wv1[0]) + [F[T[0]]]
        tot = mg = 0.0
        for e, co in chq:
            tm = float(co)
            for xx, kk in zip(pt9, e):
                if kk:
                    tm *= xx ** kk
            tot += tm
            mg += abs(tm)
        worst_c = max(worst_c, abs(tot) / mg)
        th = [d + d * d / (b + c), d * d / p + d * d / b, d * d / c + d * d / q]
        okf &= a >= min(th) * (1 - 1e-9)
        nall3 += int(a >= max(th))
        slack.append(a / min(th))
        dvals.append(d)
    if not feas:
        ck.info("NUMERIC: the search for points of the stratum itself (%d starts, seed 20261138, every gap "
                "outside T >= %.2f times the mean gap) did not converge on this platform: NOT tested"
                % (nstart, margin))
        return
    ck.check("NUMERIC", "points of the stratum itself: Levenberg-Marquardt search (%d starts, seed 20261138) for "
             "x = alpha Psi(a, b) with gaps 0 on T and every other gap >= %.2f times the mean gap: %d points "
             "found (contact set exactly T; d between %.3f and %.3f of the mean gap).  At all of them "
             "|C + HQ| / (sum of |terms|) <= %.1e, |R9| / (sum of |terms|) <= %.1e (tolerance 1e-9), and "
             "a >= the smallest of the three thresholds (ratio between %.0f and %.0f; all three conditions "
             "hold at %d of the points)"
             % (nstart, margin, len(feas), min(dvals), max(dvals), worst_c, worst_r,
                min(slack), max(slack), nall3),
             okf and worst_c < 1e-9 and worst_r < 1e-9)
    ck.info("the stratum is therefore not empty for d > 0; at the points found d is at or near the imposed "
            "lower bound and the conditions hold with a large margin, so these points do not probe the "
            "thresholds (the thresholds are probed by the algebraic tests above)")


# ----------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------


def main():
    ck = K.Checker("Replay of case T_34 = {2,3,5,9,10,12} (group V): document "
                   "blind_moc4_bad70_full_safe_disjunction_20261022.pdf")
    data = load_supplement(ck)
    if data is None:
        ck.finish()
    K.check_foundations(ck)
    K.crosscheck_with_audit_library(ck, NU, T34)
    gens = check_contact_set_and_generators(ck)
    check_covariance_gstar(ck)
    K.check_minor_leibniz(ck, [j for j in gens if K.F_IS_MINOR[j]])
    K.check_radial_necessity(ck)
    check_old_H(ck)
    check_chain_gram(ck, data)
    out = check_radial(ck, data)
    if out is not None:
        Rrad, R9, Cd, Qd = out
        Rk, k = check_lemma2(ck, data, Rrad, R9)
        exact_boundary_spot_checks(ck, Rk, k)
        numeric_part(ck, Rrad, R9, Cd, Qd)
    ck.finish()


if __name__ == "__main__":
    main()
