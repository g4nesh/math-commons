"""Independent replay of Case T_20 = {0,1,5,6,7,8} (Zheng, MOC n = 4; group II).

Document replayed: root_moc4_T19_uniform_20261022.pdf (3 pages; historical label "T19").
Nothing from the supplement is imported or executed.  If MOC4_SUPP is set, the PDF is hashed
(sha256) to pin the document; no other use.  All objects are re-derived in replay_T7_common.py.

Run:  MOC4_SUPP=/path/to/supplement python replay_T20.py     (exit code 0 iff every check passes)
Packages: sympy (exact part), numpy + scipy (NUMERIC part only).

Sections 1-7 are EXACT.  Section 8 is NUMERIC (falsification attempts; seeds printed).

Method for the kernel identities.  In the chart kappa F_0 = i mu and R_j = F_j/F_0 is a rational
function of the slope parameter q with real coefficients (x, h, k, l real), so
        g_j / (mu Im q) = -(R_j(Q) - R_j(Qb)) / (Q - Qb)     at Q = q, Qb = conj(q),
computed in Q(x,h,k,l,Q,Qb) with Q, Qb independent.  |p(q)|^2 = p(Q) p(Qb) for real polynomials p.
"""
import math
import os
import random
import sys
from fractions import Fraction

sys.dont_write_bytecode = True           # keep the output directory free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T7_common as H            # noqa: E402
from replay_T7_common import GQ         # noqa: E402

ck = H.Checks("replay_T20: Case T_20 = {0,1,5,6,7,8}, document root_moc4_T19_uniform_20261022.pdf")
T = H.PAPER_T[20]
DOC_T = (0, 1, 5, 6, 7, 8)                       # the document's "T19"
DOC_GENS = [0, 1, 2, 3, 4, 5, 6, 27, 28, 33, 34]
DOC_COORD = [3, 4, 9, 10, 14, 16, 19, 22]
A_COEF = {3: 1, 2: -1, 18: 1}                    # A = g3 - g2 + g18
B_COEF = {13: 1, 3: -1}                          # B = g13 - g3

# ---------------------------------------------------------------------------------------------
ck.section("1. Objects (EXACT)")
H.document_hash_check(ck, "root_moc4_T19_uniform_20261022.pdf")
H.facet_table_sanity(ck)
H.crosscheck_with_audit_lib(ck)
ck.check("document's contact set 'T19' = paper's T_20 (identical set, same indexing; group element = identity)",
         tuple(sorted(DOC_T)) == tuple(sorted(T)))

# ---------------------------------------------------------------------------------------------
ck.section("2. Contact cone and complete minor dual, doc p.1 (EXACT)")
gens, minor, nonminor, facets, cols = H.minor_cone_report(
    ck, T, DOC_GENS, [28, 34], DOC_COORD, {"A": A_COEF, "B": B_COEF}, 10)
ck.check("the minor cone has nine generators", len(minor) == 9)

# ---------------------------------------------------------------------------------------------
ck.section("3. Rectangle argument, the curve (2) and its parametrization (3), doc p.1-2 (EXACT)")
K, A, C, F = H.rectangle_checks(ck)
ck.check("F5/F0 = (A1-C3)(A3-C1)/((A1-C1)(A3-C3))",
         F[5] / F[0] == (A[1] - C[3]) * (A[3] - C[1]) / ((A[1] - C[1]) * (A[3] - C[3])))
ck.check("F5 - F0 = (A0-C0)(A2-C2)(A1-A3)(C3-C1)   (nonzero for distinct labels, so h != 1)",
         F[5] - F[0] == (A[0] - C[0]) * (A[2] - C[2]) * (A[1] - A[3]) * (C[3] - C[1]))

Kn, gn = H.make_field("x,A2,A3,v,w")
x, A2, A3, v, w = [gn[n] for n in ["x", "A2", "A3", "v", "w"]]
Fn = H.products([x, Kn(1), A2, A3], [None, Kn(0), v, w], inf_C=0)
ck.check("normalized labels A = (x,1,A2,A3), C = (inf,0,v,w): F6/F0 = x, F8/F6 = (1-v)/(A2-v)",
         Fn[6] / Fn[0] == x and Fn[8] / Fn[6] == (1 - v) / (A2 - v))
ck.check("F5/F0 = (1 - w) A3/(A3 - w)", Fn[5] / Fn[0] == (1 - w) * A3 / (A3 - w))
ck.check("F1/F0 = (A2 - w)(A3 - v)/((A2 - v)(A3 - w))",
         Fn[1] / Fn[0] == (A2 - w) * (A3 - v) / ((A2 - v) * (A3 - w)))

Ke, ge = H.make_field("h,k,l,v,w,q")
h, k, l, v, w, q = [ge[n] for n in ["h", "k", "l", "v", "w", "q"]]
a2 = 1 + k * (v - 1)
a3 = h * w / (h - 1 + w)
ck.check("A3 = hw/(h-1+w) is the solution of F5/F0 = h", (1 - w) * a3 / (a3 - w) == h)
lhs2 = l * (a2 - v) * w * (1 - w)
rhs2 = (a2 - w) * ((h - v) * w - v * (h - 1))
ck.check("curve (2): F1/F0 - l = -[l(a2-v)w(1-w) - (a2-w){(h-v)w - v(h-1)}] / ((a2-v) w (1-w))",
         (a2 - w) * (a3 - v) / ((a2 - v) * (a3 - w)) - l == -(lhs2 - rhs2) / ((a2 - v) * w * (1 - w)))
c = 1 + l * (k - 1)
Tq = c * q - k
wq = 1 + q * (v - 1)
a2q = a2
lhs2q = l * (a2q - v) * wq * (1 - wq)
rhs2q = (a2q - wq) * ((h - v) * wq - v * (h - 1))
ck.check("with w = 1 + q(v-1): lhs(2) - rhs(2) = -(v-1)^2 [ v q T(q) - (q-1){h(q-k) + l(k-1)q} ], linear in v",
         lhs2q - rhs2q == -(v - 1) ** 2 * (v * q * Tq - (q - 1) * (h * (q - k) + l * (k - 1) * q)))
Kq, gq_ = H.make_field("h,k,l")
h_, k_, l_ = [gq_[n] for n in ["h", "k", "l"]]
c_ = 1 + l_ * (k_ - 1)
q0 = k_ / c_
ck.check("at T(q) = 0 (q = k/c): (q-1){h(q-k)+l(k-1)q} = k l (h-1)(k-1)^2 (l-1)/c^2  (nonzero, so T(q) != 0)",
         (q0 - 1) * (h_ * (q0 - k_) + l_ * (k_ - 1) * q0) == k_ * l_ * (h_ - 1) * (k_ - 1) ** 2 * (l_ - 1) / c_ ** 2)

# ---------------------------------------------------------------------------------------------
ck.section("4. Ratios, kernels (4), Pick identities (5), (6), obstruction functions (8), doc p.2-3 (EXACT)")
K3, g3 = H.make_field("x,h,k,l,Q,Qb")
NAMES = ["x", "h", "k", "l", "Q", "Qb"]
x, h, k, l, Q, Qb = [g3[n] for n in NAMES]
c = 1 + l * (k - 1)


def Tpol(q):
    return c * q - k


def Jpol(q):
    return (h - 1) * q ** 2 + (-h * k - h - k * l + k + l) * q + h * k


def Spol(q):
    return (h * k + c - k) * q - h * k


def chart(q):
    v = (q - 1) * (h * (q - k) + l * (k - 1) * q) / (q * Tpol(q))
    w = 1 + q * (v - 1)
    a2 = 1 + k * (v - 1)
    a3 = h * w / (h - 1 + w)
    Fv = H.products([x, K3(1), a2, a3], [None, K3(0), v, w], inf_C=0)
    return [f / Fv[0] for f in Fv], v, w


R1, v1, w1 = chart(Q)
R2, v2, w2 = chart(Qb)


def dd(a, b):
    """divided difference (a - b)/(Q - Qb) = Im a / Im q for a = f(q), b = f(conj q), f real."""
    return (a - b) / (Q - Qb)


G = [-dd(a, b) for a, b in zip(R1, R2)]                  # G_j = g_j / (mu Im q)
ck.check("(R0,R1,R5,R6,R7,R8) = (1, l, h, x, lx, -x/(k-1))",
         [R1[j] for j in (0, 1, 5, 6, 7, 8)] == [K3(1), l, h, x, l * x, -x / (k - 1)])
ck.check("the chart satisfies all six contacts and no other gap vanishes identically",
         all(G[j] == 0 for j in T) and all(G[j] != 0 for j in range(24) if j not in T))
ck.check("R14 = (x - v)/((k-1)(v-1)) and R16 = h R14",
         R1[14] == (x - v1) / ((k - 1) * (v1 - 1)) and R1[16] == h * R1[14])
ck.check("hence g14 = mu (x-1) Im v / ((k-1)|v-1|^2) and g16 = h g14",
         G[14] == (x - 1) * dd(v1, v2) / ((k - 1) * (v1 - 1) * (v2 - 1)) and G[16] == h * G[14])
T1, T2 = Tpol(Q), Tpol(Qb)
ck.check("(4) g3 = mu h k Im q / ((k-1)|q|^2)", G[3] == h * k / ((k - 1) * Q * Qb))
ck.check("(4) g4/g3 = -l (h-1)(k-1)^2 (l-1)|q|^2 / (h |T(q)|^2)",
         G[4] / G[3] == -l * (h - 1) * (k - 1) ** 2 * (l - 1) * Q * Qb / (h * T1 * T2))
ck.check("(4) g9/g3 = -x (l-1)(k-1)|q|^2 / (h k |q-1|^2)",
         G[9] / G[3] == -x * (l - 1) * (k - 1) * Q * Qb / (h * k * (Q - 1) * (Qb - 1)))
ck.check("(4) g10/g3 = l x (k-1)|q|^2 / (h |q-k|^2)",
         G[10] / G[3] == l * x * (k - 1) * Q * Qb / (h * (Q - k) * (Qb - k)))
ck.check("(5) Pick identity: Im v / Im q = h/|q|^2 - l (h-1)(k-1)^2 (l-1)/|T(q)|^2",
         dd(v1, v2) == h / (Q * Qb) - l * (h - 1) * (k - 1) ** 2 * (l - 1) / (T1 * T2))
ck.check("g14/g3 = (x-1)|q|^2 (Im v/Im q) / (h k |v-1|^2)",
         G[14] / G[3] == (x - 1) * Q * Qb * dd(v1, v2) / (h * k * (v1 - 1) * (v2 - 1)))
ck.check("(6) Im w / Im q = ((h-1)/c) {1 - k l (k-1)^2 (l-1)/|T(q)|^2}",
         dd(w1, w2) == (h - 1) / c * (1 - k * l * (k - 1) ** 2 * (l - 1) / (T1 * T2)))
ck.check("R19 = h - h(x-1)/(w-1), hence g19 = -mu h (x-1) Im w / |w-1|^2",
         R1[19] == h - h * (x - 1) / (w1 - 1)
         and G[19] == -h * (x - 1) * dd(w1, w2) / ((w1 - 1) * (w2 - 1)))
J1, J2 = Jpol(Q), Jpol(Qb)
ck.check("v - 1 = J(q)/(q T(q))", v1 - 1 == J1 / (Q * T1))
AR = (1 - x) * l * Spol(Q) / J1 + h + l - 1
BR = -h * (1 - x) * (Q - k) * T1 / ((k - 1) * J1)
ck.check("(8) R3 - R2 + R18 = A_R exactly (A = g3 - g2 + g18 = -mu Im A_R)", R1[3] - R1[2] + R1[18] == AR)
ck.check("(8) R13 - R3 = B_R exactly (B = g13 - g3 = -mu Im B_R)", R1[13] - R1[3] == BR)
p = k / c
ck.check("J(0) = hk,  J(p) = k l (h-1)(k-1)^2 (l-1)/c^2,  J(1) = -(k-1)(l-1)   (p = k/c)",
         Jpol(K3(0)) == h * k and Jpol(p) == k * l * (h - 1) * (k - 1) ** 2 * (l - 1) / c ** 2
         and Jpol(K3(1)) == -(k - 1) * (l - 1))
ck.check("slope of S = hk + (k-1)(l-1),  S(0) = -hk,  S(p) = k (h-1)(k-c)/c",
         h * k + c - k == h * k + (k - 1) * (l - 1) and Spol(K3(0)) == -h * k
         and Spol(p) == k * (h - 1) * (k - c) / c)
z0 = h * k / (h * k + c - k)
ck.check("zero z0 of S:  J(z0) = -h k (h-1)(k-1)(l-1)/(hk + c - k)^2",
         Spol(z0) == 0 and Jpol(z0) == -h * k * (h - 1) * (k - 1) * (l - 1) / (h * k + c - k) ** 2)
ck.check("T(q) = c (q - p)", T1 == c * (Q - p))

# ---------------------------------------------------------------------------------------------
ck.section("5. Sign analysis, doc p.2-3 (EXACT: complete enumeration of sign cells; certificates)")
# After removing the squared moduli each premise is a rational function of (x,h,k,l) alone; its
# sign is constant on each of the 81 cells of (k,x,l,h) relative to 0 and 1 and is certified.
Kp, gp = H.make_field("a,b,c,d")
pa, pb, pc, pd = [gp[n] for n in ["a", "b", "c", "d"]]
premises = [h,                                                   # h > 0 (g16 = h g14, g14 != 0)
            G[4] / G[3] * T1 * T2 / (Q * Qb),                    # g4/g3  * |T|^2/|q|^2
            G[9] / G[3] * (Q - 1) * (Qb - 1) / (Q * Qb),         # g9/g3  * |q-1|^2/|q|^2
            G[10] / G[3] * (Q - k) * (Qb - k) / (Q * Qb),        # g10/g3 * |q-k|^2/|q|^2
            (x - 1) / (h * k)]                                   # g14/g3 divided by positive factors
table = H.cell_signs(premises, NAMES, ["k", "x", "l", "h"], Kp, [pb, pd, pc, pa])
ck.check("the premises, stripped of squared moduli, are free of q and have a certified sign on each of the "
         "81 cells", len(table) == 81 and all(0 not in sg for _, sg in table))
rows = [(lab, sg[4]) for lab, sg in table if min(sg[:4]) > 0]
doc_rows = [  # (k, sign of x, l, h) of the document's table
    (">1", 1, "(0,1)", ">1"), ("(0,1)", -1, "(0,1)", ">1"), ("<0", 1, "<0", "(0,1)"), ("<0", -1, ">1", "(0,1)")]
got = set((lab["k"], -1 if lab["x"] == "<0" else 1, lab["l"], lab["h"]) for lab, _ in rows)
ck.check("the sign cells compatible with h > 0 and the three positive ratios of (4) are exactly the four "
         "rows of the document's table (each x > 0 row splits into 0<x<1 and x>1)",
         got == set(doc_rows) and len(rows) == 6, "%d cells" % len(rows))
surv = [lab for lab, s5 in rows if s5 > 0]
exp_surv = [{"k": ">1", "x": ">1", "l": "(0,1)", "h": ">1"},
            {"k": "<0", "x": "(0,1)", "l": "<0", "h": "(0,1)"},
            {"k": "<0", "x": "<0", "l": ">1", "h": "(0,1)"}]
ck.check("adding (x-1)/k > 0 (from g14/g3 > 0 and the Pick identity) leaves exactly: "
         "(k>1, x>1, 0<l<1, h>1), (k<0, 0<x<1, l<0, 0<h<1), (k<0, x<0, l>1, 0<h<1)",
         len(surv) == 3 and all(e in surv for e in exp_surv))


def sub(elem, val):
    """Substitute a region parametrization into an element of Q(x,h,k,l,Q,Qb) free of Q, Qb."""
    return H.substitute(elem, NAMES, val, Kp)


def sgn(elem, val):
    return H.certified_sign(sub(K3(elem), val))


row1 = dict(k=1 + pb, x=1 + pd, l=1 / (1 + pc), h=1 + pa)          # k>1, x>1, 0<l<1, h>1
row3 = dict(k=-pb, x=1 / (1 + pd), l=-pc, h=1 / (1 + pa))          # k<0, 0<x<1, l<0, 0<h<1
row4 = dict(k=-pb, x=-pd, l=1 + pc, h=1 / (1 + pa))                # k<0, x<0, l>1, 0<h<1  = region (7)
for name, val in (("first", row1), ("third", row3)):
    ck.check("%s case: c > 1 and -k l (l-1) > 0 (so the bracket of (6) exceeds 1); h k/(k-1) > 0 (mu Im q > 0)"
             % name, sgn(c - 1, val) == 1 and sgn(-k * l * (l - 1), val) == 1 and sgn(h * k / (k - 1), val) == 1)
    ck.check("%s case: g19/g3 = -(x-1)(k-1)|q|^2 (Im w/Im q)/(k |w-1|^2) has the sign of -(x-1)(k-1)(h-1)/(k c) < 0"
             % name,
             G[19] / G[3] == -(x - 1) * (k - 1) * Q * Qb * dd(w1, w2) / (k * (w1 - 1) * (w2 - 1))
             and sgn(-(x - 1) * (k - 1) * (h - 1) / (k * c), val) == -1)
ck.info("=> only region (7): 0 < h < 1, k < 0, l > 1, x < 0, mu Im q > 0 can be supported")
ck.check("region (7): c < 0, 0 < p < 1 (p = k/c), h k/(k-1) > 0",
         sgn(c, row4) == -1 and sgn(p, row4) == 1 and sgn(1 - p, row4) == 1 and sgn(h * k / (k - 1), row4) == 1)
ck.check("region (7): J(0) < 0, J(p) > 0, J(1) > 0, leading coefficient h - 1 < 0",
         sgn(h * k, row4) == -1 and sgn(k * l * (h - 1) * (k - 1) ** 2 * (l - 1) / c ** 2, row4) == 1
         and sgn(-(k - 1) * (l - 1), row4) == 1 and sgn(h - 1, row4) == -1)
ck.check("region (7): slope of S < 0, S(0) > 0, S(p) < 0, J(z0) > 0",
         sgn(h * k + c - k, row4) == -1 and sgn(-h * k, row4) == 1
         and sgn(k * (h - 1) * (k - c) / c, row4) == -1
         and sgn(-h * k * (h - 1) * (k - 1) * (l - 1) / (h * k + c - k) ** 2, row4) == 1)
ck.check("region (7): prefactors (1-x) l > 0 and -h(1-x)/(k-1) > 0",
         sgn((1 - x) * l, row4) == 1 and sgn(-h * (1 - x) / (k - 1), row4) == 1)

# ---------------------------------------------------------------------------------------------
ck.section("6. Second, root-free proof of A > 0 and B > 0 on region (7) (EXACT; own argument)")
# A/(mu Im q) |J|^2 and B/(mu Im q) |J|^2 are real quadratic forms a|q|^2 + b(q + conj q) + c0.
# If a > 0 and a c0 - b^2 > 0 the form is positive for every complex q.


def quad_coeffs(elem):
    """Coefficients (a, b, c0) of elem = a Q Qb + b (Q + Qb) + c0; raises unless elem has this form."""
    assert all(mon[4] == 0 and mon[5] == 0 for mon, _ in elem.denom.terms())   # denominator free of Q, Qb
    den = K3(elem.denom)
    ring5 = elem.numer.ring
    co = {}
    for mon, cf in elem.numer.terms():
        key = (mon[4], mon[5])
        co[key] = co.get(key, ring5.zero) + ring5({(mon[0], mon[1], mon[2], mon[3], 0, 0): cf})
    assert set(co) <= {(1, 1), (1, 0), (0, 1), (0, 0)} and co.get((1, 0), 0) == co.get((0, 1), 0)
    return K3(co.get((1, 1), 0)) / den, K3(co.get((1, 0), 0)) / den, K3(co.get((0, 0), 0)) / den


GA = (G[3] - G[2] + G[18]) * J1 * J2
GB = (G[13] - G[3]) * J1 * J2
for name, form in (("A", GA), ("B", GB)):
    qa, qb, qc = quad_coeffs(form)
    ck.check("region (7): (%s/(mu Im q)) |J(q)|^2 = a|q|^2 + b(q+conj q) + c0 with a > 0 and a c0 - b^2 > 0 "
             "(certificates), so %s > 0 for every nonreal q with mu Im q > 0" % (name, name),
             sgn(qa, row4) == 1 and sgn(qa * qc - qb ** 2, row4) == 1)

# ---------------------------------------------------------------------------------------------
ck.section("6b. Observation, not in the document: v real forces q real, i.e. all eight labels real (EXACT)")
ck.check("Im v/Im q = h [ 1/|q|^2 + (g4/g3)/|q|^2 ]:  it is nonzero whenever g3, g4 > 0",
         dd(v1, v2) == h * (1 / (Q * Qb) + (G[4] / G[3]) / (Q * Qb)))
ck.info("=> for distinct labels and all gaps >= 0: Im q != 0 gives g3 != 0 and g4 != 0 by (4), hence Im v != 0.")
ck.info("   So v real implies q real; then w = 1 + q(v-1), A2 and A3 are real: all eight labels are real,")
ck.info("   every ratio R_j is real and g_j = -mu Im R_j = 0 for all j.  No circle theorem is needed there.")

# ---------------------------------------------------------------------------------------------
ck.section("7. Exact evaluation at Gaussian-rational points, directly from the 24 products (EXACT)")


def labels_q(xq, hq, kq, lq, qq):
    cq = 1 + lq * (kq - 1)
    vq = (qq - 1) * (hq * (qq - kq) + lq * (kq - 1) * qq) / (qq * (cq * qq - kq))
    wq = 1 + qq * (vq - 1)
    return [GQ(xq), GQ(1), 1 + kq * (vq - 1), hq * wq / (hq - 1 + wq)], [None, GQ(0), vq, wq]


exact_points = [
    # (x, h, k, l, Re q, Im q, mu): three points in region (7), three outside
    (Fraction(-3, 2), Fraction(1, 2), Fraction(-1, 3), Fraction(5, 2), Fraction(1, 4), Fraction(3, 5), 1),
    (Fraction(-1, 5), Fraction(4, 5), Fraction(-3, 1), Fraction(8, 7), Fraction(-2, 1), Fraction(1, 3), 1),
    (Fraction(-7, 1), Fraction(1, 9), Fraction(-1, 10), Fraction(11, 1), Fraction(3, 2), Fraction(-2, 7), -1),
    (Fraction(3, 2), Fraction(5, 2), Fraction(7, 3), Fraction(2, 5), Fraction(1, 4), Fraction(3, 5), 1),
    (Fraction(1, 2), Fraction(1, 3), Fraction(-2, 1), Fraction(-5, 2), Fraction(1, 3), Fraction(4, 7), 1),
    (Fraction(-4, 3), Fraction(5, 2), Fraction(2, 7), Fraction(-3, 5), Fraction(1, 3), Fraction(2, 9), 1),
]
ok_all = True
supp_pts = []
for (xq, hq, kq, lq, pr, pi_, mu) in exact_points:
    qq = GQ(pr, pi_)
    Aq, Cq = labels_q(xq, hq, kq, lq, qq)
    gq = H.exact_ratio_gaps(Aq, Cq, mu)
    pt = {"x": GQ(xq), "h": GQ(hq), "k": GQ(kq), "l": GQ(lq), "Q": qq, "Qb": H.gq_conj(qq)}
    for j in range(24):
        val = H.eval_frac(G[j], NAMES, pt)
        if not (val.im == 0 and val.re * mu * pi_ == gq[j]):
            ok_all = False
    if not all(gq[j] == 0 for j in T):
        ok_all = False
    labs = Aq + Cq[1:]
    if len(set((u.re, u.im) for u in labs)) != 7:
        ok_all = False
    if min(gq) >= 0:
        Aval = gq[3] - gq[2] + gq[18]
        Bval = gq[13] - gq[3]
        supp_pts.append(Aval > 0 and Bval > 0 and 0 < hq < 1 and kq < 0 and lq > 1 and xq < 0)
ck.check("at 6 Gaussian-rational chart points: 24 gaps from the products = mu Im q * (polarized formula), "
         "contacts vanish exactly, labels distinct", ok_all)
ok_real = True
for (xq, hq, kq, lq, qr, mu) in [
        (Fraction(-3, 2), Fraction(1, 2), Fraction(-1, 3), Fraction(5, 2), Fraction(-4, 7), 1),
        (Fraction(3, 2), Fraction(5, 2), Fraction(7, 3), Fraction(2, 5), Fraction(9, 4), -1)]:
    Aq, Cq = labels_q(xq, hq, kq, lq, GQ(qr, 0))
    gq = H.exact_ratio_gaps(Aq, Cq, mu)
    labs = Aq + Cq[1:]
    ok_real = ok_real and all(u == 0 for u in gq) and all(u.im == 0 for u in labs) \
        and len(set((u.re, u.im) for u in labs)) == 7
ck.check("at 2 chart points with q real (all labels real, distinct): all 24 gaps are exactly 0", ok_real)
ck.check("exactly the 3 points of region (7) are supported; there A > 0 and B > 0",
         len(supp_pts) == 3 and all(supp_pts))

# ---------------------------------------------------------------------------------------------
ck.section("8. NUMERIC falsification attempts (floating point prefilter, exact confirmation)")
import numpy as np   # noqa: E402

np.seterr(all="ignore")     # degenerate random samples are filtered explicitly below
SEED = 20261125
rng = random.Random(SEED)


def labels_f(xf, hf, kf, lf, qf):
    cf = 1 + lf * (kf - 1)
    vf = (qf - 1) * (hf * (qf - kf) + lf * (kf - 1) * qf) / (qf * (cf * qf - kf))
    wf = 1 + qf * (vf - 1)
    return [xf, 1.0, 1 + kf * (vf - 1), hf * wf / (hf - 1 + wf)], [None, 0.0, vf, wf]


def rnd_real():
    u = rng.random()
    if u < 0.4:
        return rng.choice([-1, 1]) * math.exp(rng.uniform(-3, 3))
    if u < 0.5:
        return rng.choice([0.0, 1.0]) + rng.choice([-1, 1]) * 10 ** rng.uniform(-6, -1)
    return rng.uniform(-3, 3)


def rnd_unit():
    u = rng.random()
    if u < 0.2:
        return 10 ** rng.uniform(-6, -1)
    if u < 0.4:
        return 1 - 10 ** rng.uniform(-6, -1)
    return rng.uniform(0.001, 0.999)


def rnd_above_one():
    u = rng.random()
    if u < 0.3:
        return 1 + 10 ** rng.uniform(-6, -1)
    return 1 + math.exp(rng.uniform(-3, 4))


def rnd_neg():
    return -math.exp(rng.uniform(-7, 4))


def to_q(val):
    return Fraction(val).limit_denominator(10 ** 7)


stats = dict(supported=0, region=0, region_supported=0, violations=0, contact_fail=0, minA=None, minB=None)


def exact_tests(xq, hq, kq, lq, pr, pi_, mu):
    try:
        Aq, Cq = labels_q(xq, hq, kq, lq, GQ(pr, pi_))
        gq = H.exact_ratio_gaps(Aq, Cq, mu)
    except ZeroDivisionError:
        return
    if not all(gq[j] == 0 for j in T):
        stats["contact_fail"] += 1
        return
    in7 = 0 < hq < 1 and kq < 0 and lq > 1 and xq < 0
    supported = min(gq) >= 0
    Aval = gq[3] - gq[2] + gq[18]
    Bval = gq[13] - gq[3]
    ok = True
    if supported:
        stats["supported"] += 1
        ok = in7 and mu * pi_ > 0 and Aval > 0 and Bval > 0
        ok = ok and all(gq[j] > 0 for j in (3, 4, 9, 10, 14, 16))
        for key, nrm in facets.items():
            if sum(a * gq[cc] for a, cc in zip(nrm, cols)) < 0:
                ok = False
    if in7:
        # the document's intermediate claims on region (7) (they do not need support)
        stats["region"] += 1
        g3_ = gq[3]
        ok = ok and Aval / g3_ > 0 and Bval / g3_ > 0 and all(gq[j] / g3_ > 0 for j in (4, 9, 10))
        ok = ok and (g3_ > 0) == (mu * pi_ > 0)
        if supported:
            stats["region_supported"] += 1
        ra, rb = Aval / g3_, Bval / g3_
        stats["minA"] = ra if stats["minA"] is None else min(stats["minA"], ra)
        stats["minB"] = rb if stats["minB"] is None else min(stats["minB"], rb)
    if not ok:
        stats["violations"] += 1


def run_batch(sampler, n, prefilter):
    cand = 0
    for _ in range(n):
        xf, hf, kf, lf = sampler()
        qf = complex(rng.uniform(-3, 3), rng.choice([-1, 1]) * 10 ** rng.uniform(-4, 0.5))
        mu = rng.choice([-1, 1])
        if prefilter:
            try:
                Af, Cf = labels_f(xf, hf, kf, lf, qf)
                gf = H.float_ratio_gaps(Af, Cf, mu)
            except ZeroDivisionError:
                continue
            sc = np.max(np.abs(gf))
            if not np.isfinite(sc) or sc == 0 or gf.min() < -1e-9 * sc:
                continue
        cand += 1
        exact_tests(to_q(xf), to_q(hf), to_q(kf), to_q(lf), to_q(qf.real), to_q(qf.imag), mu)
    return cand


N1 = 100000
c1 = run_batch(lambda: (rnd_real(), rnd_real(), rnd_real(), rnd_real()), N1, True)
ck.info("N1a: %d unrestricted chart samples (seed %d): %d float-supported candidates re-evaluated exactly"
        % (N1, SEED, c1))
N1b = 6000
c2 = run_batch(lambda: (rnd_neg(), rnd_unit(), rnd_neg(), rnd_above_one()), N1b, False)
ck.info("N1b: %d exact samples in region (7) (with mass near the boundary values)" % c2)
ck.check("N1 (EXACT at sampled points): contacts vanish; every supported sample lies in region (7) with "
         "mu Im q > 0, g3,g4,g9,g10,g14,g16 > 0, A > 0, B > 0 and all 10 facet inequalities; on region (7) "
         "A/g3 > 0, B/g3 > 0 and the kernel ratios are positive",
         stats["violations"] == 0 and stats["contact_fail"] == 0 and stats["supported"] > 500
         and stats["region"] > 3000,
         "supported %d; region (7) samples %d (supported %d); min A/g3 = %.3e, min B/g3 = %.3e"
         % (stats["supported"], stats["region"], stats["region_supported"],
            float(stats["minA"]), float(stats["minB"])))

# N2: chart-free solutions of the six contact equalities
nprng = np.random.default_rng(SEED)
N2 = 300
found = used = 0
worst_H = worst_real = worst_rebuild = 0.0
for it in range(N2):
    res = H.find_contact_point(T, nprng, tries=3)
    if res is None:
        continue
    found += 1
    a, cc, kap, gf = res
    if not H.labels_well_separated(a, cc):
        continue
    Fv = H.float_products(a, cc)
    mob = H.moebius_to(cc[0], cc[1], a[1])           # C0 -> infinity, C1 -> 0, A1 -> 1
    xn, _, a2n, a3n = [mob(u) for u in a]
    vn, wn = mob(cc[2]), mob(cc[3])
    if max(abs(u) for u in (xn, a2n, a3n, vn, wn)) > 1e4 or abs(vn - 1) < 1e-3:
        continue
    used += 1
    worst_H = max(worst_H, abs((kap * Fv[0]).real) / abs(Fv[0]))
    kn = (a2n - 1) / (vn - 1)
    hn = Fv[5] / Fv[0]
    ln = Fv[1] / Fv[0]
    qn = (wn - 1) / (vn - 1)
    worst_real = max(worst_real, max(abs(u.imag) / (1 + abs(u)) for u in (xn, kn, hn, ln)))
    mu_n = (kap * Fv[0]).imag
    Af, Cf = labels_f(xn.real, hn.real, kn.real, ln.real, qn)
    g_chart = H.float_ratio_gaps(Af, Cf, mu_n)
    worst_rebuild = max(worst_rebuild, np.max(np.abs(gf - g_chart)) / np.max(np.abs(gf)))
ck.check("N2 (NUMERIC): at chart-free solutions of the six contact equalities with well separated labels: "
         "Re(kappa F0) = 0, the parameters x, h, k, l are real, and the chart (with q = (w-1)/(v-1)) "
         "reproduces all 24 gaps",
         used >= 30 and worst_H < 1e-7 and worst_real < 1e-6 and worst_rebuild < 1e-6,
         "%d searches, %d converged, %d usable; max |Re(kappa F0)|/|F0| = %.1e; max |Im(x,h,k,l)|/(1+|.|) = %.1e; "
         "max relative gap mismatch = %.1e" % (N2, found, used, worst_H, worst_real, worst_rebuild))

# N3: the dual criterion against nonnegative least squares
N3 = 3000
mism = nin = 0
Mall = np.array([H.F_VALUES[j] for j in gens], float)
for it in range(N3):
    lam = nprng.exponential(size=len(gens)) * (nprng.random(len(gens)) < 0.6)
    for j in nonminor:
        lam[gens.index(j)] *= nprng.choice([0.0, 0.3, 1.0, 3.0])
    gv = lam @ Mall
    if gv.max() == 0:
        continue
    gv = gv / gv.max()
    Aval = sum(cf * gv[j] for j, cf in A_COEF.items())
    Bval = sum(cf * gv[j] for j, cf in B_COEF.items())
    if abs(Aval) < 1e-7 or abs(Bval) < 1e-7:
        continue
    inside_dual = (Aval >= 0 and Bval >= 0)
    nin += inside_dual
    if inside_dual != (H.minor_cone_nnls_residual(gv, minor) < 1e-9):
        mism += 1
ck.check("N3 (NUMERIC): on %d random points of the contact cone, 'A >= 0 and B >= 0' agrees with "
         "a nonnegative-least-squares test of membership in the minor cone" % N3, mism == 0,
         "%d inside, mismatches %d" % (nin, mism))

ck.finish()
