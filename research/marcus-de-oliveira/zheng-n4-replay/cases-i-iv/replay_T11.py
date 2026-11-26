"""Independent exact replay of the case T_11 = {0,1,3,5,6,8} of Zheng's n = 4 proof.

Proof documents replayed (read as PDFs; historical number "T10", contact set written as {0,1,2,4,7,9}):
    final_referee_moc4_T10_uniform_20261022.pdf   (2 pages, assembly; relies on BOTH branch documents)  -> sections 1-3
    root_moc4_T10_B_branch_20261022.pdf           (4 pages, branch B < 0, four-coefficient inequality)  -> sections 4, 5
    final_referee_T10_A_branch_20261022.pdf       (4 pages, branch A < 0, SOC inequality)               -> sections 4, 6
No program or data file of the supplement is read, imported or executed (MOC4_SUPP is not used).
Run:  python replay_T11.py     (exit code 0 iff every check passes; about one minute).

Chart of the documents: labels 0..3, A_sp = (infinity, 0, 1, z), C_sp = -B = (w, v, r, c),
    F_sigma = prod_i [A_i, C_sigma(i)]  (bracket of the point at infinity with a finite point = +1),
    R_j = (F_j - F_0)/((z-c)(v-r)),  kappa = i mu/((z-c)(v-r)),  g_j = -mu Im R_j = -rho ghat_j,  rho = mu Im v,
    ghat_j = (R_j(v) - R_j(vb))/(v - vb)  (v, vb independent symbols; h, t, s, l real).
"""
import sys

import sympy as sp
from sympy.polys.fields import field

import replay_T10_common as C
from replay_T10_common import feq, fzero

ck = C.Checker("Replay of T_11 = {0,1,3,5,6,8}  (documents: final_referee_moc4_T10_uniform, root_moc4_T10_B_branch, "
               "final_referee_T10_A_branch)")
F = C.F
T11 = [0, 1, 3, 5, 6, 8]
TD = [0, 1, 2, 4, 7, 9]            # the contact set as written in the three documents
T7 = [0, 1, 2, 4, 7, 8]

# ---------------------------------------------------------------------------------------------------------------
ck.section("0. objects (own construction) and cross-check with the audit library")
C.basic_selfchecks(ck)
C.crosscheck_library(ck, {11: set(T11), 7: set(T7)})

# ---------------------------------------------------------------------------------------------------------------
ck.section("1. correspondence with the paper's T_11; generators; minor case; canonical decompositions (uniform p.1-2)")
ks = C.find_set_maps(TD, T11)
ck.check("exactly one element of G maps the document's set {0,1,2,4,7,9} to the paper's T_11 = {0,1,3,5,6,8}: "
         "sigma -> sigma o (3 4) (a relabelling of the rows, i.e. of the A-spectrum)",
         len(ks) == 1 and C.G_ELEMENTS[ks[0]] == ((0, 1, 2, 3), (0, 1, 3, 2), 0), [C.G_ELEMENTS[k_] for k_ in ks])
kg = ks[0]
gens = C.generators(TD)
ck.check("generators vanishing on {0,1,2,4,7,9} are exactly {0,1,2,3,4,5,7,11,12,21,22,26}",
         gens == [0, 1, 2, 3, 4, 5, 7, 11, 12, 21, 22, 26], gens)
ck.check("this group element maps them onto the generators of the paper's T_11",
         sorted(C.g_facet(C.G_MAPS[kg], j) for j in gens) == C.generators(T11))
ck.check("Leibniz covariance for this group element on all 82 vectors: Phi_{g.f}(O') = sgn(g) Phi_f(O) identically in O",
         C.covariance_identity(kg, range(82)))
ck.check("f_22 and f_26 are the only non-minor generators", [j for j in gens if not C.IS_MINOR[j]] == [22, 26])
Afun = C.functional({16: 1, 3: -1})
Bfun = C.functional({3: -1, 6: 1, 12: -1, 13: 1, 18: 1})
va = {j: Afun(F[j]) for j in gens}
vbv = {j: Bfun(F[j]) for j in gens}
ck.check("A = g16 - g3 is 1 on f_1, f_4, f_7, f_22, -1 on f_26, 0 on the other generators",
         va == {0: 0, 1: 1, 2: 0, 3: 0, 4: 1, 5: 0, 7: 1, 11: 0, 12: 0, 21: 0, 22: 1, 26: -1}, va)
ck.check("B = -g3 + g6 - g12 + g13 + g18 is 1 on f_0, f_2, f_12, f_26, -1 on f_22, 0 on the other generators",
         vbv == {0: 1, 1: 0, 2: 1, 3: 0, 4: 0, 5: 0, 7: 0, 11: 0, 12: 1, 21: 0, 22: -1, 26: 1}, vbv)


def eqp(lhs, rhs):
    return C.lin_comb(lhs) == C.lin_comb(rhs)


ids = [({22: 1, 26: 1}, {3: 1, 5: 1, 11: 1, 21: 2}), ({26: 1, 1: 1}, {0: 1, 11: 1, 21: 1}),
       ({26: 1, 4: 1}, {2: 1, 3: 1, 11: 1, 21: 1}), ({26: 1, 7: 1}, {3: 1, 5: 1, 12: 1, 21: 1}),
       ({22: 1, 0: 1}, {1: 1, 3: 1, 5: 1, 21: 1}), ({22: 1, 2: 1}, {4: 1, 5: 1, 21: 1}),
       ({22: 1, 12: 1}, {7: 1, 11: 1, 21: 1})]
ck.check("the seven minor replacement identities of p.1 hold, with minor vectors only on the right",
         all(eqp(a_, b_) for a_, b_ in ids) and all(C.IS_MINOR[j] for _, b_ in ids for j in b_))
eA = {16: -1, 3: 1}
coefA = {0: {23: 1}, 2: {15: 1}, 3: {14: 1}, 5: {10: 1}, 11: {17: 1}, 12: {8: 1}, 21: {16: 1}, 26: eA}
dB = {3: 1, 6: -1, 12: 1, 13: -1, 18: -1}
coefB = {1: {23: 1}, 3: {18: 1}, 4: {15: 1}, 5: {21: 1}, 7: {8: 1}, 11: {6: 1}, 21: {19: 1}, 22: dB}
ck.check("A-branch identity: g = g23 f0 + g15 f2 + g14 f3 + g10 f5 + g17 f11 + g8 f12 + g16 f21 + e f26 (e = -A) on the "
         "span of the generators", all(C.lin_comb({k_: C.functional(e_)(F[j]) for k_, e_ in coefA.items()}) == F[j] for j in gens))
ck.check("B-branch identity: g = g23 f1 + g18 f3 + g15 f4 + g21 f5 + g8 f7 + g6 f11 + g19 f21 + d f22 (d = -B) on the "
         "span of the generators", all(C.lin_comb({k_: C.functional(e_)(F[j]) for k_, e_ in coefB.items()}) == F[j] for j in gens))
ck.check("the span of the generators has dimension 8", sp.Matrix([F[j] for j in gens]).rank() == 8)


def lin0(*terms):
    acc = {}
    for dct, cf in terms:
        for k_, v_ in dct.items():
            acc[k_] = acc.get(k_, 0) + cf * v_
    return all(C.functional(acc)(F[j]) == 0 for j in gens)


ck.check("linear identities on the span: g19 = g3 - d, g21 = g10 - g8 - d, g16 = g3 - e, g17 = g6 - g8 - e, "
         "d = -(g15 + g8 + g23 + e)",
         lin0(({19: 1}, 1), ({3: 1}, -1), (dB, 1)) and lin0(({21: 1}, 1), ({10: 1}, -1), ({8: 1}, 1), (dB, 1))
         and lin0(({16: 1}, 1), ({3: 1}, -1), (eA, 1)) and lin0(({17: 1}, 1), ({6: 1}, -1), ({8: 1}, 1), (eA, 1))
         and lin0((dB, 1), ({15: 1}, 1), ({8: 1}, 1), ({23: 1}, 1), (eA, 1)))

# ---------------------------------------------------------------------------------------------------------------
ck.section("2. the universal inequalities invoked (uniform p.2; A-branch p.1; B-branch p.1)")
ck.check("orbit sizes: SOC 144 (pair unordered), four-coefficient 288 (a, b unordered)",
         C.orbit_size(C.SOC_REP, 2) == 144 and C.orbit_size(C.FOUR_REP, 2) == 288)
C.check_atom_use(ck, "SOC (f_2, f_12; f_26)", C.SOC_REP, (2, 12, 26), 2)
C.check_atom_use(ck, "four-coefficient (f_4, f_7; f_3; f_22)", C.FOUR_REP, (4, 7, 3, 22), 2)
ck.check("roles: A-branch a = g15 on f_2, b = g8 on f_12, e on f_26; B-branch a = g15 on f_4, b = g8 on f_7, "
         "c = g18 on f_3, d on f_22",
         coefA[2] == {15: 1} and coefA[12] == {8: 1} and coefB[4] == {15: 1} and coefB[7] == {8: 1} and coefB[3] == {18: 1})
Os = C.rational_so4_samples(6, 1101)
okp = True
for O in Os:
    x2, x12, p26 = C.minor_value(2, O), C.minor_value(12, O), C.leibniz(F[26], O)
    x4, x7, x3, p22 = C.minor_value(4, O), C.minor_value(7, O), C.minor_value(3, O), C.leibniz(F[22], O)
    okp &= (p26 + 2 * abs(x2 * x12) >= 0)
    okp &= (p22 + x3 ** 2 + 2 * (abs(x4 * x7) + abs(x4 * x3) + abs(x7 * x3)) >= 0)
ck.check("sanity (6 exact rational SO(4) points, not a proof): Phi_26 + 2|x_2 x_12| >= 0 and "
         "Phi_22 + Phi_3 + 2(|x_4 x_7| + |x_4 x_3| + |x_7 x_3|) >= 0", okp)

# ---------------------------------------------------------------------------------------------------------------
ck.section("3. the boundary b = g_8 = 0 (both branch documents): the cited family {0,1,2,4,7,8}")
ck.check("{0,1,2,4,7,9} + {8} contains {0,1,2,4,7,8}, which is literally the paper's T_7 (Appendix E)",
         set(T7) <= set(TD) | {8} and sorted(T7) == [0, 1, 2, 4, 7, 8])

# ---------------------------------------------------------------------------------------------------------------
ck.section("4. the chart (B-branch p.1-2; A-branch p.1)")
Kg, zg, wg, vg, rg, cg = field("z w v r c", sp.ZZ)
Fg = C.field_products(Kg.one, [None, Kg.zero, Kg.one, zg], [wg, vg, rg, cg])
Rg = [(Fg[j] - Fg[0]) / ((zg - cg) * (vg - rg)) for j in range(24)]
h_ = zg * (vg - 1) / (vg * (zg - 1))
t_ = (zg - cg) * (vg - rg) / ((rg - cg) * (vg - zg))
s_ = (wg - vg) * (1 - cg) * (zg - rg) / ((zg - cg) * (vg - rg))
l_ = (wg - rg) * zg / ((wg - vg) * (zg - rg))
H_, k_ = h_ - 1, 1 + t_ * (h_ - 1)
ck.check("for EVERY configuration, with h, t, s, l defined by the cross-ratio formulas of B-branch p.1: R_1 = -1/(t(h-1)), "
         "R_2 = 1, R_4 = 1 + h R_1, R_7 = R_1 - s, R_9 = R_1 + s(l-1); hence real contacts force h, t, s, l real",
         feq(Rg[1], -1 / (t_ * H_)) and feq(Rg[2], 1) and feq(Rg[4], 1 + h_ * Rg[1]) and feq(Rg[7], Rg[1] - s_)
         and feq(Rg[9], Rg[1] + s_ * (l_ - 1)))
ck.check("so the four parameters are quotients of contact differences: h = (F_4 - F_2)/(F_1 - F_0), "
         "t (h-1) = -(F_2 - F_0)/(F_1 - F_0), s = (F_1 - F_7)/(F_2 - F_0), l - 1 = (F_9 - F_1)/(F_1 - F_7)",
         feq(h_, (Fg[4] - Fg[2]) / (Fg[1] - Fg[0])) and feq(t_ * H_, -(Fg[2] - Fg[0]) / (Fg[1] - Fg[0]))
         and feq(s_, (Fg[1] - Fg[7]) / (Fg[2] - Fg[0])) and feq(l_ - 1, (Fg[9] - Fg[1]) / (Fg[1] - Fg[7])))
ck.check("for EVERY configuration the chart formulas hold identically: z = h v/(1+(h-1)v), c = (z(v-r) - t r(v-z))/"
         "((v-r) - t(v-z)), w = v + s(z-c)(v-r)/((1-c)(z-r)), r = h v((k-1)s(l-1)-k)/(((h-1) l s (k-1) - h(k-1)) v + l s (k-1) - h)",
         feq(zg, h_ * vg / (1 + H_ * vg)) and feq(cg, (zg * (vg - rg) - t_ * rg * (vg - zg)) / ((vg - rg) - t_ * (vg - zg)))
         and feq(wg, vg + s_ * (zg - cg) * (vg - rg) / ((1 - cg) * (zg - rg)))
         and feq(rg, h_ * vg * ((k_ - 1) * s_ * (l_ - 1) - k_) / ((H_ * l_ * s_ * (k_ - 1) - h_ * (k_ - 1)) * vg + l_ * s_ * (k_ - 1) - h_)))
zs_, ws_, vs_, rs_, cs_ = sp.symbols("z w v r c")
ck.check("degenerations are collisions: h = 1 iff z = v; t = 1 iff (v-c)(z-r) = 0; k = 1 iff (z-c)(v-r) = 0 "
         "(numerators of h-1, t-1, k-1 factor accordingly)",
         C.factor_strings([(h_ - 1).numer.as_expr()]) <= C.factor_strings([vs_ - zs_])
         and C.factor_strings([(t_ - 1).numer.as_expr()]) <= C.factor_strings([vs_ - cs_, zs_ - rs_])
         and C.factor_strings([(k_ - 1).numer.as_expr()]) <= C.factor_strings([zs_ - cs_, vs_ - rs_]))

K, v, vb, h, t, s, l = field("v vb h t s l", sp.ZZ)
H, k = h - 1, 1 + t * (h - 1)
D = h - k


def chart(vv, ll):
    zz = h * vv / (1 + H * vv)
    rr = h * vv * ((k - 1) * s * (ll - 1) - k) / ((H * ll * s * (k - 1) - h * (k - 1)) * vv + ll * s * (k - 1) - h)
    cc = (zz * (vv - rr) - t * rr * (vv - zz)) / ((vv - rr) - t * (vv - zz))
    ww = vv + s * (zz - cc) * (vv - rr) / ((1 - cc) * (zz - rr))
    Fp = C.field_products(K.one, [None, K.zero, K.one, zz], [ww, vv, rr, cc])
    return [(Fp[j] - Fp[0]) / ((zz - cc) * (vv - rr)) for j in range(24)]


R = chart(v, l)
Rb = chart(vb, l)
gh = [(R[j] - Rb[j]) / (v - vb) for j in range(24)]
ck.check("on the chart (R_0,R_1,R_2,R_4,R_7,R_9) = (0, -1/(tH), 1, 1 - h/(tH), -s - 1/(tH), s(l-1) - 1/(tH)): all real",
         all(feq(a1, b1) for a1, b1 in zip([R[j] for j in TD],
                                            [0, -1 / (t * H), 1, 1 - h / (t * H), -s - 1 / (t * H), s * (l - 1) - 1 / (t * H)])))
Cc = h * (k - 1) + l * D
J0 = l * D - h
N0 = k - s * (k - 1) * (l - 1)
M1 = lambda vv: (D + s * H * (k - 1)) * vv + s * (k - 1) - D
M2 = lambda vv: Cc * vv - J0
hs, ts, ss, ls, vs = sp.symbols("h t s l v")
M1s, M2s = C.fexpr(M1(v)), C.fexpr(M2(v))
found = C.denominator_factors(R)
allowed = C.factor_strings([hs, hs - 1, ts, vs - 1, M1s, M2s])
ck.check("every R_j is in Q(h,t,s,l,v) with denominator a product of h, h-1, t, v-1, M_1(v), M_2(v): for real v with "
         "M_1(v) M_2(v) != 0 every R_j is real, i.e. the whole gap vector vanishes", found <= allowed, sorted(found - allowed))
kS_, lS_, sS_ = sp.symbols("k l s")
ck.check("C - J0 = h k, so M_2 = C v - J0 vanishes identically in v only if k = 0 and l D = h, i.e. (D = h at k = 0) "
         "l = 1, where N0 = k - s(k-1)(l-1) = 0; and M_1 = 0 identically would need D(1 + H) = D h = 0",
         feq(Cc - J0, h * k) and (kS_ - sS_ * (kS_ - 1) * (lS_ - 1)).subs({kS_: 0, lS_: 1}) == 0
         and feq((D + s * H * (k - 1)) + H * (s * (k - 1) - D) * (-1) - (D + H * D), 0))
ah, bh, Gh, Vh, qh = gh[15], gh[8], gh[3], gh[10], gh[23]
dh = gh[3] - gh[6] + gh[12] - gh[13] - gh[18]            # d = -B
eh = gh[3] - gh[16]                                      # e = -A
ck.check("B (1): g18 = (t-1)(a + d), g19 = G - d, g21 = V - b - d  (a = g15, b = g8, G = g3, V = g10, d = -B)",
         feq(gh[18], (t - 1) * (ah + dh)) and feq(gh[19], Gh - dh) and feq(gh[21], Vh - bh - dh))
ck.check("B (2) = A (3): q = (h-k) a - k d - h b and J0 g6 = -h b  (q = g23)",
         feq(qh, (h - k) * ah - k * dh - h * bh) and feq(J0 * gh[6], -h * bh))
ck.check("B (3) = A (4): a = -rho l s (h-1)/h, G = rho h s D N0/|M1|^2, b = -rho h J0 N0/|M2|^2  (g = -rho ghat)",
         feq(ah, l * s * (h - 1) / h) and feq(Gh, -h * s * D * N0 / (M1(v) * M1(vb))) and feq(bh, h * J0 * N0 / (M2(v) * M2(vb))))
ck.check("B (4) = A (5): k V = -rho h s/|v-1|^2 - D b (cross-multiplied, no division by k)",
         feq(k * Vh, h * s / ((v - 1) * (vb - 1)) - D * bh))
ck.check("A (2): g16 = G - e, g17 = g6 - b - e, g14 = b - t(b + q + e), g18 = (1-t)(b + q + e); A (3): d = -(a + b + q + e)",
         feq(gh[16], Gh - eh) and feq(gh[17], gh[6] - bh - eh) and feq(gh[14], bh - t * (bh + qh + eh))
         and feq(gh[18], (1 - t) * (bh + qh + eh)) and feq(dh, -(ah + bh + qh + eh)))

# ---------------------------------------------------------------------------------------------------------------
ck.section("5. branch B < 0 (root_moc4_T10_B_branch): four-coefficient capacity")
Yb_ = h * (h - 1) / Cc
NP = h + s * (k - 1) * (h - l)
Pb = D * NP / (s * Cc * (k - 1))
ck.check("B p.2: D = (t-1)(1-h); q/a = (t-1)(1+delta) - h(t-1+eta+t delta) for a = 1 (delta = d/a, eta = b/a)",
         feq(D, (t - 1) * (1 - h)) and feq(qh / ah, (t - 1) * (1 + dh / ah) - h * (t - 1 + bh / ah + t * dh / ah)))
ck.check("B p.2: C = h k + J0;  d = P a - Y b with P = D(h + s(k-1)(h-l))/(s C (k-1)), Y = h(h-1)/C;  "
         "(l-1) N_P + (h-l) N0 = C", feq(Cc, h * k + J0) and feq(dh, Pb * ah - Yb_ * bh) and feq((l - 1) * NP + (h - l) * N0, Cc))
ck.check("B p.2: Y < 1/t iff l D < 0 given C < 0:  t h (h-1) - C = -l D", feq(t * h * (h - 1) - Cc, -l * D))
eta, dl = bh / ah, dh / ah
p_, T0 = 1 - h, t - 1
Zb = 1 - t * Yb_
Eb = p_ + k * Yb_
Jb = 1 + dl + Yb_ * eta
Dn = Jb - Yb_
Wb = T0 * p_ - k * dl - k * Yb_ * eta
B0 = h * eta + k * dl
xi, xib = p_ * v - Eb, p_ * vb - Eb
ck.check("B p.3: E = p + k Y = p Z + Y;  W = T0 p - k delta - k Y eta = -N0 Dn;  h Z = -l T0 Y",
         feq(Eb, p_ * Zb + Yb_) and feq(Wb, -N0 * Dn) and feq(h * Zb, -l * T0 * Yb_))
ck.check("B p.3: R = |xi|^2 = t Y E W/(Z eta) with xi = p v - E;  G/a = h Y W Dn/(Z |J xi - Y W|^2)",
         feq(xi * xib, t * Yb_ * Eb * Wb / (Zb * eta)) and feq(Gh / ah, h * Yb_ * Wb * Dn / (Zb * (Jb * xi - Yb_ * Wb) * (Jb * xib - Yb_ * Wb))))
ck.check("B p.3: k g21/a = h^2 p/((-l)|xi + k Y|^2) - B0, B0 = h eta + k delta; and xi + k Y = p (v - 1)",
         feq(k * gh[21] / ah, h * h * p_ / ((-l) * (xi + k * Yb_) * (xib + k * Yb_)) - B0) and feq(xi + k * Yb_, p_ * (v - 1)))
# abstract interval identities in the independent variables (h, t, Y, eta, delta, Re xi)
Ka, ha, ta, Ya, ea, da, xr = field("h t Y eta delta xr", sp.QQ)
pa, T0a, ka = 1 - ha, ta - 1, 1 - ta * (1 - ha)
Za, Ea, Ja = 1 - ta * Ya, (1 - ha) + ka * Ya, 1 + da + Ya * ea
Dna = Ja - Ya
Wa = T0a * pa - ka * da - ka * Ya * ea
B0a = ha * ea + ka * da
Ra_ = ta * Ya * Ea * Wa / (Za * ea)
SG, SV = ha * Ya * Wa * Dna / (Za * da), ha * pa * T0a * Ya / (Za * B0a)
LG = (Ja * Ja * Ra_ + Ya * Ya * Wa * Wa - SG) / (2 * Ja * Ya * Wa)
LV = (SV - Ra_ - ka * ka * Ya * Ya) / (2 * ka * Ya)
hS, tS, YS, eS, dS = sp.symbols("h t Y eta delta")
B0S = hS * eS + (1 - tS * (1 - hS)) * dS
ck.check("B p.3: the reduced denominator of L_V is 2 eta (tY - 1) B0 up to sign: no factor k",
         sp.expand(LV.denom.as_expr() - 2 * eS * (tS * YS - 1) * B0S) == 0
         or sp.expand(LV.denom.as_expr() + 2 * eS * (tS * YS - 1) * B0S) == 0, sp.factor(LV.denom.as_expr()))
nJ = Ja * Ja * Ra_ - 2 * Ja * Ya * Wa * xr + Ya * Ya * Wa * Wa
nK = Ra_ + 2 * ka * Ya * xr + ka * ka * Ya * Ya
ck.check("B p.3: g19 = G - d = 2 delta J Y W (Re xi - L_G)/|J xi - Y W|^2 and g21 = 2 Y B0 (L_V - Re xi)/|xi + k Y|^2 "
         "(a = 1; |J xi - YW|^2 = J^2 R - 2 J Y W Re xi + Y^2 W^2, |xi + kY|^2 = R + 2 k Y Re xi + k^2 Y^2)",
         feq(ha * Ya * Wa * Dna / (Za * nJ) - da, 2 * da * Ja * Ya * Wa / nJ * (xr - LG))
         and feq((ha * pa * T0a * Ya / (Za * nK) - B0a) / ka, 2 * Ya * B0a / nK * (LV - xr)))
Nf = -(LV - LG) * 2 * ea * da * Za * Ja * B0a / ha
ck.check("B p.3: L_V - L_G = -h N(h,Y)/(2 eta delta Z J B0) with N a polynomial", Nf.denom == 1)
Ne = sp.expand(Nf.numer.as_expr())
Kcap = eS * (tS - 1 - dS) - tS * dS ** 2
hstar = (tS - 1) * (1 + dS) / (tS - 1 + eS + tS * dS)
z0 = lambda e_: sp.simplify(sp.together(e_)) == 0
Nh = sp.diff(Ne, hS)
ck.check("B p.3: N is affine in h and affine in Y", sp.degree(Ne, hS) == 1 and sp.degree(Ne, YS) == 1)
ck.check("B p.3: N_h(h,0) = -(eta + t delta)(eta + delta)(delta + 1);  N_h(h,1/t) = -eta (eta + t delta)(eta + t(delta+1) - 1)/t",
         z0(Nh.subs(YS, 0) + (eS + tS * dS) * (eS + dS) * (dS + 1))
         and z0(Nh.subs(YS, 1 / tS) + eS * (eS + tS * dS) * (eS + tS * (dS + 1) - 1) / tS))
ck.check("B p.3: N(h*,0) = -(eta + delta)(delta + 1) K/(eta + t delta + t - 1);  N(h*,1/t) = (eta + delta)(t delta - K)/t, "
         "K = eta (t - 1 - delta) - t delta^2, h* = (t-1)(1+delta)/(t-1+eta+t delta)",
         z0(Ne.subs({hS: hstar, YS: 0}) + (eS + dS) * (dS + 1) * Kcap / (eS + tS * dS + tS - 1))
         and z0(Ne.subs({hS: hstar, YS: 1 / tS}) - (eS + dS) * (tS * dS - Kcap) / tS))
a_, b_, d_ = sp.symbols("a b d")
c_ = (tS - 1) * (a_ + d_)
ck.check("B p.3: a b (c - d) - d^2 (a + b + c + d) = (a + d)(b((t-1)a - d) - t d^2) for c = (t-1)(a + d)",
         sp.expand(a_ * b_ * (c_ - d_) - d_ ** 2 * (a_ + b_ + c_ + d_) - (a_ + d_) * (b_ * ((tS - 1) * a_ - d_) - tS * d_ ** 2)) == 0)

# ---------------------------------------------------------------------------------------------------------------
ck.section("6. branch A < 0 (final_referee_T10_A_branch): SOC bound g15 g8 >= e^2")
b, G, V, q, e, f, g6 = bh / ah, Gh / ah, Vh / ah, qh / ah, eh / ah, gh[14] / ah, gh[6] / ah     # normalisation a = 1
Mm = b + q + e
th = l * D / h
be = 1 - th
ck.check("A p.2: t = (b - f)/M (M = b + q + e, f = g14); g6 = b/beta (theta = l D/h, beta = 1 - theta); "
         "h(1-f) + e + f = 0; H(1-f) = -(1+e); D = H(1-t); J0 = -h beta",
         feq(t * Mm, b - f) and feq(g6 * be, b) and fzero(h * (1 - f) + e + f) and feq(H * (1 - f), -(1 + e))
         and feq(D, H * (1 - t)) and feq(J0, -h * be))
Dl = 1 + e - th * Mm
al = t * H + th
ga = Dl + H
Wa_ = -k * Dl - D + th * h
et = h / H
xiA, xiAb = al * v + be, al * vb + be
ck.check("A p.3: alpha = tH + theta = C/h;  (8): (1-t)/s = t Delta and W = -N0 Delta  (Delta = 1 + e - theta M, "
         "W = -k Delta - D + theta h)", feq(al, Cc / h) and feq((1 - t) / s, t * Dl) and feq(Wa_, -N0 * Dl))
ck.check("A p.3: R = |xi|^2 = t beta W/(theta b) with xi = alpha v + beta", feq(xiA * xiAb, t * be * Wa_ / (th * b)))
ck.check("A p.3 (9): M_1 = (D/(Delta alpha))(gamma xi + W);  G = eta Delta alpha^2 W/(theta |gamma xi + W|^2);  "
         "k V = D (eta alpha^2/(theta |xi - k|^2) - b)   (gamma = Delta + H, eta = h/H)",
         feq(M1(v), D / (Dl * al) * (ga * xiA + Wa_))
         and feq(G, et * Dl * al * al * Wa_ / (th * (ga * xiA + Wa_) * (ga * xiAb + Wa_)))
         and feq(k * V, D * (et * al * al / (th * (xiA - k) * (xiAb - k)) - b)))
ck.check("A p.3-4: gamma = -(f(1+e) + theta M (1-f))/(1-f);  W = -Delta - H + theta h + t H (1 - Delta)",
         feq(ga, -(f * (1 + e) + th * Mm * (1 - f)) / (1 - f)) and feq(Wa_, -Dl - H + th * h + t * H * (1 - Dl)))
# abstract identities in the independent variables (b, e, f, theta, M, Re xi)
Kq, bq, eq_, fq, thq, Mq, xq = field("b e f theta M xr", sp.QQ)
hq, Hq = -(eq_ + fq) / (1 - fq), -(1 + eq_) / (1 - fq)
tq = (bq - fq) / Mq
kq, Dq, beq = 1 + tq * Hq, Hq * (1 - tq), 1 - thq
Dlq, alq = 1 + eq_ - thq * Mq, tq * Hq + thq
gaq = Dlq + Hq
Wq = -kq * Dlq - Dq + thq * hq
etq = hq / Hq
Rq = tq * beq * Wq / (thq * bq)
SGq, SVq = etq * Dlq * alq * alq * Wq / (thq * eq_), etq * alq * alq / (thq * bq)
LGq = (SGq - gaq * gaq * Rq - Wq * Wq) / (2 * gaq * Wq)
LVq = (Rq + kq * kq - SVq) / (2 * kq)
bS, eS2, fS, thS, MS, dS2 = sp.symbols("b e f theta M delta")
ck.check("A p.3: the reduced denominator of L_V is 2 b theta (1+e)(f-1) M up to sign: no factor k",
         sp.expand(LVq.denom.as_expr() - 2 * bS * thS * (1 + eS2) * (fS - 1) * MS) == 0
         or sp.expand(LVq.denom.as_expr() + 2 * bS * thS * (1 + eS2) * (fS - 1) * MS) == 0, sp.factor(LVq.denom.as_expr()))
nG = gaq * gaq * Rq + 2 * gaq * Wq * xq + Wq * Wq
nKq = Rq - 2 * kq * xq + kq * kq
ck.check("A (10): G - e = 2 e gamma W (L_G - Re xi)/|gamma xi + W|^2 and V = 2 D b (Re xi - L_V)/|xi - k|^2",
         feq(etq * Dlq * alq * alq * Wq / (thq * nG) - eq_, 2 * eq_ * gaq * Wq / nG * (LGq - xq))
         and feq(Dq * (etq * alq * alq / (thq * nKq) - bq) / kq, 2 * Dq * bq / nKq * (xq - LVq)))
Lq = thq * (bq + eq_) * Mq - bq * (1 + eq_)
ck.check("A (11): L_G - L_V = -eta L alpha^2/(2 b e theta gamma), L = theta (b+e) M - b (1+e);  "
         "L - (e^2 - b) = (theta (b+e) - e) M + e q, q = M - b - e",
         feq(LGq - LVq, -etq * Lq * alq * alq / (2 * bq * eq_ * thq * gaq))
         and feq(Lq - (eq_ * eq_ - bq), (thq * (bq + eq_) - eq_) * Mq + eq_ * (Mq - bq - eq_)))
# negative t: M = E/theta, E = 1 + e + delta, Delta = -delta
K2, b2, e2, f2, th2, dl2 = field("b e f theta delta", sp.QQ)
E2 = 1 + e2 + dl2
h2, H2, M2_ = -(e2 + f2) / (1 - f2), -(1 + e2) / (1 - f2), E2 / th2
t2 = (b2 - f2) / M2_
k2, D2_, be2 = 1 + t2 * H2, H2 * (1 - t2), 1 - th2
Dl2, al2 = 1 + e2 - th2 * M2_, t2 * H2 + th2
W2 = -k2 * Dl2 - D2_ + th2 * h2
R2_ = t2 * be2 * W2 / (th2 * b2)
SV2 = (h2 / H2) * al2 * al2 / (th2 * b2)
LV2 = (R2_ + k2 * k2 - SV2) / (2 * k2)
ck.check("A p.4: with Delta = -delta: t = -theta (f-b)/E, W = delta - H + theta h - theta (f-b) H (1+delta)/E, E = 1 + e + delta",
         feq(Dl2, -dl2) and feq(t2, -th2 * (f2 - b2) / E2)
         and feq(W2, dl2 - H2 + th2 * h2 - th2 * (f2 - b2) * H2 * (1 + dl2) / E2))
We = C.fexpr(W2)
HS, hS2 = (1 + eS2) / (fS - 1), (eS2 + fS) / (fS - 1)
d0 = (1 - bS) * HS
ck.check("A p.4: E = f H at delta = H; W(delta = H) = theta h b/f; dW/ddelta = 1 - theta (f-b) H e/E^2",
         z0((1 + eS2 + HS) - fS * HS) and z0(We.subs(dS2, HS) - thS * hS2 * bS / fS)
         and z0(sp.diff(We, dS2) - (1 - thS * (fS - bS) * HS * eS2 / (1 + eS2 + dS2) ** 2)))
ck.check("A p.4: at delta0 = (1-b) H: E = (f-b) H and W = -beta b H",
         z0(1 + eS2 + d0 - (fS - bS) * HS) and z0(We.subs(dS2, d0) + (1 - thS) * bS * HS))
WE = sp.expand(sp.cancel(We * (1 + eS2 + dS2)))
ck.check("A p.4: W E = (1+delta) E - theta (f-b) H (1+delta) - beta h E is a quadratic in delta with leading coefficient 1",
         sp.Poly(WE, dS2).degree() == 2 and sp.Poly(WE, dS2).LC() == 1
         and z0(WE - ((1 + dS2) * (1 + eS2 + dS2) - thS * (fS - bS) * HS * (1 + dS2) - (1 - thS) * hS2 * (1 + eS2 + dS2))))
Pcal = ((1 + e2 - e2 * th2) * (f2 - b2) - (b2 * (1 + e2) - th2 * (b2 + e2))) ** 2 + 4 * b2 * (1 + e2) ** 2 * (1 - th2) * (f2 - b2)
ck.check("A p.4 (SOS): L_V^2 - R = alpha^2 Pcal/(4 b^2 theta^2 (1+e)^2), Pcal = [(1+e-e theta)(f-b) - (b(1+e) - theta(b+e))]^2 "
         "+ 4 b (1+e)^2 (1-theta)(f-b)", feq(LV2 * LV2 - R2_, al2 * al2 * Pcal / (4 * b2 * b2 * th2 * th2 * (1 + e2) ** 2)))
ck.check("A p.4: along the auxiliary segment the reduced denominator of L_V is 2 b (1+e)(f-1) E up to sign",
         sp.expand(LV2.denom.as_expr() - 2 * bS * (1 + eS2) * (fS - 1) * (1 + eS2 + dS2)) == 0
         or sp.expand(LV2.denom.as_expr() + 2 * bS * (1 + eS2) * (fS - 1) * (1 + eS2 + dS2)) == 0)
ck.check("A p.4: at delta0: alpha = 0, R = beta^2, L_V = beta",
         z0(C.fexpr(al2).subs(dS2, d0)) and z0(C.fexpr(R2_).subs(dS2, d0) - (1 - thS) ** 2)
         and z0(C.fexpr(LV2).subs(dS2, d0) - (1 - thS)))
ck.check("A p.4: R = -beta (f-b) W/(b E) on the segment (positive when W < 0 < E, f > b)", feq(R2_, -be2 * (f2 - b2) * W2 / (b2 * E2)))

# merged stratum C = 0: l = h(1-k)/D
lm = h * (1 - k) / D
ck.check("merged stratum: C = h(k-1) + l D = 0 for l = h(1-k)/D, and then theta = l D/h = 1 - k", fzero(h * (k - 1) + lm * D) and feq(lm * D / h, 1 - k))
Rm = chart(v, lm)
Rmb = chart(vb, lm)
gm = [(Rm[j] - Rmb[j]) / (v - vb) for j in range(24)]
ck.check("merged stratum: the six contacts still hold", all(fzero(gm[j]) for j in TD))
am = gm[15]
bm, Gm, Vm, qm, em = gm[8] / am, gm[3] / am, gm[10] / am, gm[23] / am, (gm[3] - gm[16]) / am
thm = 1 - k
bem = 1 - thm
ck.check("merged (7): D = s H theta (1-b);  k e = b H theta (1 - 2 Re v) - h;  k (e + (H theta/D) V |v-1|^2) = -b H theta |v|^2",
         feq(D, s * H * thm * (1 - bm)) and feq(k * em, bm * H * thm * (1 - (v + vb)) - h)
         and feq(k * (em + (H * thm / D) * Vm * (v - 1) * (vb - 1)), -bm * H * thm * v * vb))
pm = -H
wm, wmb = pm * bm * v + pm * (1 - bm) - 1, pm * bm * vb + pm * (1 - bm) - 1
ck.check("merged: G = p (p-1) beta b (1-b)/(theta N_w);  2 theta Re w = beta e + (2 theta - 1)(p-1) - p b theta;  "
         "theta q = -(p-1)(1-b) + beta (b+e)   (p = -H, w = p b v + p(1-b) - 1)",
         feq(Gm, pm * (pm - 1) * bem * bm * (1 - bm) / (thm * wm * wmb))
         and feq(thm * (wm + wmb), bem * em + (2 * thm - 1) * (pm - 1) - pm * bm * thm)
         and feq(thm * qm, -(pm - 1) * (1 - bm) + bem * (bm + em)))
ck.check("merged: w - (p-1) = p b (v-1);  D = -(p - theta) for h = 1 - p, k = beta",
         feq(wm - (pm - 1), pm * bm * (v - 1)) and feq(D, -(pm - thm)))
pS, bS3, eS3, thS3, NwS, RwS = sp.symbols("p b e theta Nw Rew")
beS = 1 - thS3
# |w-(p-1)|^2 >= p b (p-1)/theta  <=>  N_w >= ((p-1) beta/theta)(e - p(1-b) + 1), using 2 theta Re w = ...
Rew = (beS * eS3 + (2 * thS3 - 1) * (pS - 1) - pS * bS3 * thS3) / (2 * thS3)
ck.check("merged: N_w - 2(p-1) Re w + (p-1)^2 - p b (p-1)/theta = N_w - ((p-1) beta/theta)(e - p(1-b) + 1)",
         z0((NwS - 2 * (pS - 1) * Rew + (pS - 1) ** 2 - pS * bS3 * (pS - 1) / thS3)
            - (NwS - (pS - 1) * beS / thS3 * (eS3 - pS * (1 - bS3) + 1))))
ck.check("merged: e (e - p(1-b) + 1) <= p b (1-b) is e^2 + e <= p (1-b)(e + b)",
         sp.expand((pS * bS3 * (1 - bS3) - eS3 * (eS3 - pS * (1 - bS3) + 1)) - (pS * (1 - bS3) * (eS3 + bS3) - eS3 ** 2 - eS3)) == 0)

# ---------------------------------------------------------------------------------------------------------------
ck.section("7. NUMERIC falsification attempts (floating point; evidence only)")
import numpy as np

npr = np.random.default_rng(20261127)
TOL = 1e-9
CH = 500000


def gapsN(hv, tv, sv, lv, vv, mu):
    Hn, kn = hv - 1, 1 + tv * (hv - 1)
    zv = hv * vv / (1 + Hn * vv)
    rv = hv * vv * ((kn - 1) * sv * (lv - 1) - kn) / ((Hn * lv * sv * (kn - 1) - hv * (kn - 1)) * vv + lv * sv * (kn - 1) - hv)
    cv = (zv * (vv - rv) - tv * rv * (vv - zv)) / ((vv - rv) - tv * (vv - zv))
    wv = vv + sv * (zv - cv) * (vv - rv) / ((1 - cv) * (zv - rv))
    kap = 1j * mu / ((zv - cv) * (vv - rv))
    gN = C.numeric_gaps([None, 0.0 * vv, 1.0 + 0.0 * vv, zv], [wv, vv, rv, cv], kap)
    sc = np.abs(gN).max(axis=1)
    return gN / sc[:, None]


def mag(n, sd=1.1):
    return np.exp(npr.normal(0, sd, n))


def stats(gN, hv, tv, sv, lv, vv, mu, acc):
    ok = np.isfinite(gN).all(axis=1)
    acc["contact"] = max(acc["contact"], np.abs(gN[ok][:, TD]).max())
    supp = ok & (gN.min(axis=1) >= -TOL)
    An = gN[:, 16] - gN[:, 3]
    Bn = -gN[:, 3] + gN[:, 6] - gN[:, 12] + gN[:, 13] + gN[:, 18]
    acc["nsupp"] += int(supp.sum())
    acc["minAB"] = min(acc["minAB"], (An + Bn)[supp].min()) if supp.any() else acc["minAB"]
    mA = supp & (An < -1e-7) & (gN[:, 8] > 1e-9)
    mB = supp & (Bn < -1e-7) & (gN[:, 8] > 1e-9)
    acc["nA"] += int(mA.sum())
    acc["nB"] += int(mB.sum())
    rho = mu * vv.imag
    Dn_ = hv - (1 + tv * (hv - 1))
    theta = lv * Dn_ / hv
    if mA.any():
        acc["socA"] = min(acc["socA"], (gN[:, 15] * gN[:, 8] - An ** 2)[mA].min())
        acc["badA"] += int((mA & ~((tv > 0) & (tv < 1) & (theta > 0) & (theta < 1) & (sv > 0) & (lv > 0) & (rho < 0))).sum())
    if mB.any():
        dN = -Bn
        cap = gN[:, 15] * gN[:, 8] * (gN[:, 18] - dN) - dN ** 2 * (gN[:, 15] + gN[:, 8] + gN[:, 18] + dN)
        acc["capB"] = min(acc["capB"], cap[mB].min())
        acc["badB"] += int((mB & ~((hv > 0) & (hv < 1) & (tv > 1) & (lv < 0) & (sv > 0) & (rho < 0))).sum())


def new_acc():
    return dict(contact=0.0, nsupp=0, nA=0, nB=0, socA=np.inf, capB=np.inf, badA=0, badB=0, minAB=np.inf)


acc = new_acc()
for _ in range(8):
    sg = lambda: npr.choice([-1.0, 1.0], CH)
    hv, tv, sv, lv = sg() * mag(CH), sg() * mag(CH), sg() * mag(CH), sg() * mag(CH)
    vv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
    mu = sg()
    stats(gapsN(hv, tv, sv, lv, vv, mu), hv, tv, sv, lv, vv, mu, acc)
ck.check("NUMERIC global: 4,000,000 chart points (seed 20261127; h, t, s, l = random sign * lognormal(0, 1.1); v complex "
         "normal(0, 1.5); mu = +-1; gaps normalised to max 1): max |g| on the contact set = %.1e; supporting: %d; min(A+B) "
         "there = %.2e; supporting with A < -1e-7, b > 0: %d (outside the A-branch sign region 0<t<1, 0<theta<1, s,l>0, "
         "rho<0: %d; min(g15 g8 - e^2) = %.3e); with B < -1e-7, b > 0: %d (outside 0<h<1, t>1, l<0, s>0, rho<0: %d; "
         "min capacity = %.3e)" % (acc["contact"], acc["nsupp"], acc["minAB"], acc["nA"], acc["badA"], acc["socA"],
                                   acc["nB"], acc["badB"], acc["capB"]),
         acc["contact"] < 1e-6 and acc["badA"] == 0 and acc["badB"] == 0 and acc["socA"] > -1e-9 and acc["capB"] > -1e-12
         and acc["minAB"] > -1e-9 and acc["nsupp"] > 0)

accB = new_acc()
for _ in range(8):                                  # B-branch sign region
    hv, tv, sv, lv = npr.uniform(0, 1, CH), 1 + mag(CH), mag(CH), -mag(CH)
    vv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
    mu = -np.sign(vv.imag)
    stats(gapsN(hv, tv, sv, lv, vv, mu), hv, tv, sv, lv, vv, mu, accB)
ck.check("NUMERIC B-branch region: 4,000,000 points with 0<h<1, t>1, l<0, s>0, rho<0: supporting with B < -1e-7: %d; "
         "min of a b (c-d) - d^2 (a+b+c+d) = %.3e (claim: > 0)" % (accB["nB"], accB["capB"]),
         accB["nB"] > 50 and accB["capB"] > -1e-12)
accA = new_acc()
for _ in range(8):                                  # A-branch sign region: t in (0,1), theta in (0,1), h < 0 or h > 1
    hv = np.where(npr.random(CH) < 0.5, -mag(CH), 1 + mag(CH))
    tv = npr.uniform(0, 1, CH)
    theta = npr.uniform(0, 1, CH)
    Dn_ = hv - (1 + tv * (hv - 1))
    lv = theta * hv / Dn_
    sv = mag(CH)
    vv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
    mu = -np.sign(vv.imag)
    stats(gapsN(hv, tv, sv, lv, vv, mu), hv, tv, sv, lv, vv, mu, accA)
ck.check("NUMERIC A-branch region: 4,000,000 points with 0<t<1, 0<theta<1, s>0, rho<0, h<0 or h>1: supporting with "
         "A < -1e-7: %d; min(g15 g8 - e^2) = %.3e (tolerance -1e-9)" % (accA["nA"], accA["socA"]),
         accA["nA"] > 50 and accA["socA"] > -1e-9)
accN = new_acc()
for _ in range(8):                                  # the excluded sign t < 0 of the A-branch
    hv = np.where(npr.random(CH) < 0.5, -mag(CH), 1 + mag(CH))
    tv = -mag(CH)
    theta = npr.uniform(0, 1, CH)
    Dn_ = hv - (1 + tv * (hv - 1))
    lv = theta * hv / Dn_
    sv = mag(CH)
    vv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
    mu = -np.sign(vv.imag)
    stats(gapsN(hv, tv, sv, lv, vv, mu), hv, tv, sv, lv, vv, mu, accN)
ck.check("NUMERIC A-branch, t < 0: 4,000,000 points with t<0, 0<theta<1, s>0, rho<0: supporting with A < -1e-7 and "
         "b > 0: %d (claim: none)" % accN["nA"], accN["nA"] == 0)
# merged stratum C = 0
accM = new_acc()
for _ in range(4):
    hv = np.where(npr.random(CH) < 0.5, -mag(CH), 1 + mag(CH))
    tv = np.where(hv < 0, 1.0, -1.0) * npr.uniform(0, 1, CH) / np.abs(hv - 1)      # theta = -t H in (0, 1)
    kn = 1 + tv * (hv - 1)
    lv = hv * (1 - kn) / (hv - kn)
    sv = mag(CH)
    vv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
    mu = -np.sign(vv.imag)
    stats(gapsN(hv, tv, sv, lv, vv, mu), hv, tv, sv, lv, vv, mu, accM)
ck.check("NUMERIC merged stratum C = 0: 2,000,000 points with theta = 1 - k in (0,1), s > 0, rho < 0: supporting with "
         "A < -1e-7: %d, all with 0 < t < 1: %s; min(g15 g8 - e^2) = %.3e" % (accM["nA"], accM["badA"] == 0, accM["socA"]),
         accM["badA"] == 0 and (accM["nA"] == 0 or accM["socA"] > -1e-9))
# real v
hv, tv, sv, lv = [npr.choice([-1.0, 1.0], 200000) * mag(200000) for _ in range(4)]
vr = npr.normal(0, 1.5, 200000) + 0j
Hn, kn = hv - 1, 1 + tv * (hv - 1)
zv = hv * vr / (1 + Hn * vr)
rv = hv * vr * ((kn - 1) * sv * (lv - 1) - kn) / ((Hn * lv * sv * (kn - 1) - hv * (kn - 1)) * vr + lv * sv * (kn - 1) - hv)
cv = (zv * (vr - rv) - tv * rv * (vr - zv)) / ((vr - rv) - tv * (vr - zv))
wv = vr + sv * (zv - cv) * (vr - rv) / ((1 - cv) * (zv - rv))
gR = C.numeric_gaps([None, 0 * vr, 1 + 0 * vr, zv], [wv, vr, rv, cv], 1j / ((zv - cv) * (vr - rv)))
Rsc = np.abs(C.numeric_gaps([None, 0 * vr, 1 + 0 * vr, zv], [wv, vr, rv, cv], 1 / ((zv - cv) * (vr - rv)))).max(axis=1)
ratio = np.nanmax(np.abs(gR).max(axis=1) / Rsc)
ck.check("NUMERIC real v: 200,000 chart points with v real: max |gap| / max |R_j| = %.1e (the gap vector is zero)" % ratio,
         ratio < 1e-9)
b1 = C.numeric_atom_attack("SOC", (2, 12, 26), 20000, 6, 1102)
b2 = C.numeric_atom_attack("FOUR", (4, 7, 3, 22), 20000, 6, 1103)
ck.check("NUMERIC atoms on SO(4) (20000 random points + 6 Nelder-Mead runs each, seeds 1102, 1103): min of "
         "Phi_26 + 2|x_2 x_12| = %.3e; min of Phi_22 + Phi_3 + 2(|x_4 x_7| + |x_4 x_3| + |x_7 x_3|) = %.3e (tolerance -1e-9)"
         % (b1, b2), b1 > -1e-9 and b2 > -1e-9)

ck.finish()
