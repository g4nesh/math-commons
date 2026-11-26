"""Independent exact replay of the case T_10 = {0,1,3,5,7,8} of Zheng's n = 4 proof.

Proof documents replayed (read as PDFs; historical number "T9"):
    final_referee_moc4_T9_uniform_20261022.pdf        (4 pages, assembly)                 -> sections 1-4, 7
    final_referee_T9_first_branch_20261022.pdf        (4 pages, region I,  SOC bound)     -> section 5
    final_referee_T9_second_branch_20261022.pdf       (3 pages, region II, four-coeff.)   -> section 6
No program or data file of the supplement is read, imported or executed (MOC4_SUPP is not used).
Run:  python replay_T10.py     (exit code 0 iff every check passes; about one minute).

Chart of the documents: labels 0..3, A = (z, 0, 1, x), C = -B = (w, infinity, v, q), q = 1 + t(v-1),
    F_sigma = prod_i [A_i, C_sigma(i)]  (bracket with the point at infinity = -1),
    R_j = (F_j - F_0)/((z-w)(v-1)),  kappa = i mu/((z-w)(v-1)),  g_j = Re(kappa(F_j - F_0)) = -mu Im R_j = -rho ghat_j,
    rho = mu Im v,  ghat_j = Im R_j / Im v = (R_j(v) - R_j(vb))/(v - vb)   (v, vb independent symbols; x, t, r, u real).
"""
import sys

import sympy as sp
from sympy.polys.fields import field

import replay_T10_common as C
from replay_T10_common import feq, fzero

ck = C.Checker("Replay of T_10 = {0,1,3,5,7,8}  (documents: final_referee_moc4_T9_uniform, T9_first_branch, "
               "T9_second_branch)")
F = C.F
T10 = [0, 1, 3, 5, 7, 8]

# ---------------------------------------------------------------------------------------------------------------
ck.section("0. objects (own construction) and cross-check with the audit library")
C.basic_selfchecks(ck)
C.crosscheck_library(ck, {10: set(T10)})

# ---------------------------------------------------------------------------------------------------------------
ck.section("1. contact set, generators, minor case, canonical decomposition (uniform p.1)")
gens = C.generators(T10)
ck.check("the document's contact set {0,1,3,5,7,8} is literally the paper's T_10 (Appendix E); identity group element",
         sorted(T10) == [0, 1, 3, 5, 7, 8])
ck.check("generators vanishing on T_10 are exactly {0,1,2,3,4,5,6,11,27,28}", gens == [0, 1, 2, 3, 4, 5, 6, 11, 27, 28], gens)
ck.check("f_28 is the only non-minor generator", [j for j in gens if not C.IS_MINOR[j]] == [28])
dfun = C.functional({2: 1, 18: -1})
vals = {j: dfun(F[j]) for j in gens}
ck.check("d = g_2 - g_18 equals 1 on f_28, -1 on f_0 and f_3, 0 on the other generators",
         vals == {0: -1, 1: 0, 2: 0, 3: -1, 4: 0, 5: 0, 6: 0, 11: 0, 27: 0, 28: 1}, vals)


def eqp(lhs, rhs):
    return C.lin_comb(lhs) == C.lin_comb(rhs)


ck.check("f_28 + f_0 = f_1 + f_2 + f_5 + f_27", eqp({28: 1, 0: 1}, {1: 1, 2: 1, 5: 1, 27: 1}))
ck.check("f_28 + f_3 = f_4 + f_5 + f_27", eqp({28: 1, 3: 1}, {4: 1, 5: 1, 27: 1}))
ck.check("the right-hand sides use minor vectors only", all(C.IS_MINOR[j] for j in (1, 2, 5, 27, 4)))
coef3 = {1: {22: 1}, 2: {19: 1}, 4: {14: 1}, 5: {20: 1}, 6: {9: 1}, 11: {6: 1}, 27: {18: 1}, 28: {2: 1, 18: -1}}
ck.check("(3): g = g22 f1 + g19 f2 + g14 f4 + g20 f5 + g9 f6 + g6 f11 + g18 f27 + d f28 on the span of the generators",
         all(C.lin_comb({k: C.functional(e)(F[j]) for k, e in coef3.items()}) == F[j] for j in gens))
ck.check("f_28(pi_2) = 1 > 0 and f_28(pi_18) = 0", F[28][2] == 1 and F[28][18] == 0)
ck.check("(7c): g_10 = g_20 + d on the span of the generators",
         all(C.functional({10: 1, 20: -1, 2: -1, 18: 1})(F[j]) == 0 for j in gens))
ck.check("the span of the generators has dimension 8", sp.Matrix([F[j] for j in gens]).rank() == 8)

# ---------------------------------------------------------------------------------------------------------------
ck.section("2. the universal inequalities invoked (uniform p.1 eqs.(1),(2); p.3)")
ck.check("orbit sizes: SOC 144 (pair unordered), four-coefficient 288 (a, b unordered)",
         C.orbit_size(C.SOC_REP, 2) == 144 and C.orbit_size(C.FOUR_REP, 2) == 288)
C.check_atom_use(ck, "SOC (f_1, f_11; f_28)", C.SOC_REP, (1, 11, 28), 2)
C.check_atom_use(ck, "four-coefficient (f_4, f_6; f_2; f_28)", C.FOUR_REP, (4, 6, 2, 28), 2)
ck.check("roles in (3): region I uses s = g22 (on f_1), b = g6 (on f_11); region II uses a = g14 (on f_4), "
         "B = g9 (on f_6), c = g19 (on f_2)",
         coef3[1] == {22: 1} and coef3[11] == {6: 1} and coef3[4] == {14: 1} and coef3[6] == {9: 1} and coef3[2] == {19: 1})
a_, b_, c_, d_ = sp.symbols("a b c d")
M3 = sp.Matrix([[a_, -d_, -d_], [-d_, b_, -d_], [-d_, -d_, c_ - d_]])
ck.check("pointwise four-coefficient form => coefficient form: det [[a,-d,-d],[-d,b,-d],[-d,-d,c-d]] = "
         "a b (c-d) - d^2 (a+b+c+d)", sp.expand(M3.det() - (a_ * b_ * (c_ - d_) - d_ ** 2 * (a_ + b_ + c_ + d_))) == 0)
Os = C.rational_so4_samples(6, 1001)
okp = True
for O in Os:
    x1, x11, p28 = C.minor_value(1, O), C.minor_value(11, O), C.leibniz(F[28], O)
    x4, x6, x2 = C.minor_value(4, O), C.minor_value(6, O), C.minor_value(2, O)
    okp &= (p28 + 2 * abs(x1 * x11) >= 0)
    okp &= (p28 + x2 ** 2 + 2 * (abs(x4 * x6) + abs(x4 * x2) + abs(x6 * x2)) >= 0)
ck.check("sanity (6 exact rational SO(4) points, not a proof): Phi_28 + 2|x_1 x_11| >= 0 and "
         "Phi_28 + Phi_2 + 2(|x_4 x_6| + |x_4 x_2| + |x_6 x_2|) >= 0", okp)

# ---------------------------------------------------------------------------------------------------------------
ck.section("3. the chart (uniform p.2)")
Kg, z, w, v, x, t, q = field("z w v x t q", sp.ZZ)
Fg = C.field_products(Kg.one, [z, Kg.zero, Kg.one, x], [w, None, v, q])
ck.check("F_1 - F_0 = (z-w)(q-v)(x-1), F_3 - F_5 = -(z-w)(q-v), F_5 - F_0 = -x(z-w)(v-1)",
         feq(Fg[1] - Fg[0], (z - w) * (q - v) * (x - 1)) and feq(Fg[3] - Fg[5], -(z - w) * (q - v))
         and feq(Fg[5] - Fg[0], -x * (z - w) * (v - 1)))
qq = 1 + t * (v - 1)
Yg = x - qq
Fg2 = C.field_products(Kg.one, [z, Kg.zero, Kg.one, x], [w, None, v, qq])
Rg = [(Fg2[j] - Fg2[0]) / ((z - w) * (v - 1)) for j in range(24)]
r_def = Rg[7] - Rg[1]
u_def = Rg[8]
ck.check("r := R_7 - R_1 = t z (v-x)/(z-w) and u := R_8 = -Y (v z - v + w - z)/((v-1)(z-w)), Y = x - q",
         feq(r_def, t * z * (v - x) / (z - w)) and feq(u_def, -Yg * (v * z - v + w - z) / ((v - 1) * (z - w))))
Qg = u_def * t * (v - 1) * (v - x) + Yg * (v * (r_def - t) + t * x)
ck.check("(4) holds for EVERY configuration (z, w free): z Q = r v Y and w = z [1 - t (v-x)/r] with r, u as defined",
         feq(z * Qg, r_def * v * Yg) and feq(w, z * (1 - t * (v - x) / r_def)))
ck.check("general contact responses: R_1 = (t-1)(x-1), R_3 = 1-t-x, R_5 = -x (z, w free)",
         feq(Rg[1], (t - 1) * (x - 1)) and feq(Rg[3], 1 - t - x) and feq(Rg[5], -x))

K, v, vb, x, t, r, u = field("v vb x t r u", sp.ZZ)


def Yf(vv):
    return x - 1 - t * (vv - 1)


def Qf(vv):
    return u * t * (vv - 1) * (vv - x) + Yf(vv) * (vv * (r - t) + t * x)


def chart(vv):
    qv = 1 + t * (vv - 1)
    zz = r * vv * Yf(vv) / Qf(vv)
    ww = zz * (1 - t * (vv - x) / r)
    Fp = C.field_products(K.one, [zz, K.zero, K.one, x], [ww, None, vv, qv])
    return [(Fp[j] - Fp[0]) / ((zz - ww) * (vv - 1)) for j in range(24)]


R = chart(v)
Rb = chart(vb)
gh = [(R[j] - Rb[j]) / (v - vb) for j in range(24)]
ck.check("on the chart (R_0,R_1,R_3,R_5,R_7,R_8) = (0, (t-1)(x-1), 1-t-x, -x, r+(t-1)(x-1), u): all real, so the six "
         "contacts hold", all(feq(a1, b1) for a1, b1 in zip([R[j] for j in T10],
                                                            [0, (t - 1) * (x - 1), 1 - t - x, -x, r + (t - 1) * (x - 1), u])))
xs, ts, rs, us, vs = sp.symbols("x t r u v")
Ys = xs - 1 - ts * (vs - 1)
Qs = us * ts * (vs - 1) * (vs - xs) + Ys * (vs * (rs - ts) + ts * xs)
allowed = C.factor_strings([ts, vs, vs - 1, vs - xs, Ys, Qs])
found = C.denominator_factors(R)
ck.check("every R_j is in Q(x,t,r,u,v) with denominator a product of t, v, v-1, v-x, Y, Q (all nonzero for eight "
         "distinct points): for real v every R_j is real, i.e. the whole gap vector vanishes", found <= allowed,
         sorted(found - allowed))
y, h = x - 1, 1 - t
Y, Yb, Q, Qb = Yf(v), Yf(vb), Qf(v), Qf(vb)
dh = gh[2] - gh[18]
ck.check("(5): ghat_2 = -y/|v-1|^2, ghat_6 = r y h/(t |v-x|^2), ghat_9 = x y (r/|v-x|^2 - t u/|Y|^2)",
         feq(gh[2], -y / ((v - 1) * (vb - 1))) and feq(gh[6], r * y * h / (t * (v - x) * (vb - x)))
         and feq(gh[9], x * y * (r / ((v - x) * (vb - x)) - t * u / (Y * Yb))))
ck.check("(6): (x+t-1) ghat_10 = -x h ((x+t+u-1)/|v|^2 + t u y/|Y|^2)",
         feq((x + t - 1) * gh[10], -x * h * ((x + t + u - 1) / (v * vb) + t * u * y / (Y * Yb))))
ck.check("(7a),(7b): g14 + d + t g6 + h g19 = 0 and g9 + d + x g19 - y g22 = 0 on the chart",
         fzero(gh[14] + dh + t * gh[6] + h * gh[19]) and fzero(gh[9] + dh + x * gh[19] - y * gh[22]))
ck.check("the chart gap vector satisfies the linear identity g10 = g20 + d (consistency of numbering)",
         feq(gh[10], gh[20] + dh))
ck.check("ghat_6/ghat_2 = -r h |v-1|^2/(t |v-x|^2)", feq(gh[6] / gh[2], -r * h * (v - 1) * (vb - 1) / (t * (v - x) * (vb - x))))

# ---------------------------------------------------------------------------------------------------------------
ck.section("4. exclusion of x < 0, t < 0 (uniform p.3) and of x < 0, t > 1, u >= 0 (second branch, Lemma 1)")
Cv = lambda vv: -r * (Yf(vv) + u) / Qf(vv)
m = (Cv(v) - Cv(vb)) / (v - vb)
ck.check("with rho = -1: g14 = -(1+X)(1+T) m + T d and g22 = -d + X m - X T u/|Y|^2  (x = -X, t = -T, m = Im C(v)/Im v, "
         "C(v) = -r(Y+u)/Q)",
         feq(gh[14], -(1 - x) * (1 - t) * m - t * dh) and feq(gh[22], -dh - x * m - x * t * u / (Y * Yb)))
Xp, Tp = sp.symbols("X T", positive=True)
ck.check("d/X <= T d/((1+X)(1+T)) is impossible for d > 0: (1+X)(1+T) - T X = 1 + X + T > 0",
         sp.expand((1 + Xp) * (1 + Tp) - Tp * Xp - (1 + Xp + Tp)) == 0)
Qpoly = sp.Poly(Qs, vs)
ck.check("Q is a quadratic in v with leading coefficient t(u - r + t) and constant term t x (u + x - 1 + t)",
         Qpoly.degree() == 2 and sp.expand(Qpoly.all_coeffs()[0] - ts * (us - rs + ts)) == 0
         and sp.expand(Qpoly.all_coeffs()[2] - ts * xs * (us + xs - 1 + ts)) == 0)
a0 = 1 + (xs - 1) / ts
ck.check("Y(a0) = 0 for a0 = 1 + (x-1)/t; Q(x) = r x Y(x); Q(a0) = u t (a0-1)(a0-x)",
         sp.simplify(Ys.subs(vs, a0)) == 0 and sp.expand(Qs.subs(vs, xs) - rs * xs * Ys.subs(vs, xs)) == 0
         and sp.simplify(Qs.subs(vs, a0) - us * ts * (a0 - 1) * (a0 - xs)) == 0)
Pf = lambda vv: (1 + t * (vv - 1)) * u * (vv - x) + Yf(vv) * ((r + 1 - t) * vv - x * (1 - t))
e1 = R[2] - R[18] + r * Pf(v) / Q
ck.check("R_2 - R_18 = -r P(v)/Q(v) + (a real constant, namely x(1-t))", feq(e1, x * (1 - t)))
ck.check("R_14 = u - r Y (Y+u)/Q", feq(R[14], u - r * Y * (Y + u) / Q))
ck.check("P(v) - (Y(v)+u)(v-x) = Q(v) identically, so P(theta) = (Y(theta)+u)(theta-x) at the roots of Q",
         feq(Pf(v) - (Y + u) * (v - x), Q))
fth = Ys / (vs - xs)
ck.check("f(theta) = Y(theta)/(theta - x) has derivative (1-x)(1-t)/(theta-x)^2 (> 0 for x, t < 0), and f = "
         "t (a0 - theta)/(theta - x)",
         sp.simplify(sp.diff(fth, vs) - (1 - xs) * (1 - ts) / (vs - xs) ** 2) == 0
         and sp.simplify(fth - ts * (a0 - vs) / (vs - xs)) == 0)
ck.check("Y(theta) + u vanishes at theta = a0 + u/t (= a0 + U/T)", sp.simplify((Ys + us).subs(vs, a0 + us / ts)) == 0)
# second branch Lemma 1
k = x + t - 1
L = k + u
N0, Nx, NY = v * vb, (v - x) * (vb - x), (k - t * v) * (k - t * vb)
Mq = gh[10] * N0 * NY
D2 = -x * L * k + t * x * L * (v + vb) - t * (t * x + u) * N0
ck.check("second branch L.1: N_Y = |k - t v|^2 = |Y|^2; M := ghat_10 N_0 N_Y = -x h (k L - t L (v+vb) + t (t+u) N_0)",
         feq(NY, Y * Yb) and feq(Mq, -x * h * (k * L - t * L * (v + vb) + t * (t + u) * N0)))
ck.check("second branch L.1: (ghat_2 - ghat_18) |Q|^2 = r^2 D_2 + r t u L N_x", feq(dh * Q * Qb, r * r * D2 + r * t * u * L * Nx))
ck.check("second branch L.1: D_2 = M/h - t u (1-x) N_0 and D_2 + L N_Y = M/x",
         feq(D2, Mq / h - t * u * (1 - x) * N0) and feq(D2 + L * NY, Mq / x))
ck.check("second branch L.1: d|Q|^2/r = (r - r0) D_2 + r0 (D_2 + L N_Y) with r0 = t u N_x/N_Y",
         feq(dh * Q * Qb / r, (r - t * u * Nx / NY) * D2 + (t * u * Nx / NY) * (D2 + L * NY)))

# ---------------------------------------------------------------------------------------------------------------
ck.section("5. region I: x > 1, t < 0, r > 0, u < 0 (first branch document)")
g = [gh[j] / gh[2] for j in range(24)]          # normalisation g_2 = 1
ze, zb = (v - x) / (v - 1), (vb - x) / (vb - 1)
X1, T1, U1 = x - 1, -t, -u                       # x = 1 + X, t = -T, u = -U
A0, h1, N = 1 + X1, 1 + T1, ze * zb
Qz = lambda zz: r * (zz - h1) * (zz - A0) - T1 * X1 * zz * (zz - h1 - U1)
Kz = lambda zz: r * ((X1 - U1) * zz - (h1 * X1 - U1)) / Qz(zz)
G1, H1 = r / N, T1 * U1 / ((ze - h1) * (zb - h1))
n1 = -(Kz(ze) - Kz(zb)) / (ze - zb)
d1 = -(ze * Kz(ze) - zb * Kz(zb)) / (ze - zb)
ck.check("I: Q(zeta) = (x-1) Q_v/(v-1)^2 (so Q(zeta) != 0 on the chart)", feq(Qz(ze), (x - 1) * Q / (v - 1) ** 2))
ck.check("I: d := -Im(zeta K)/Im zeta equals g_2 - g_18 (normalised g_2 = 1)", feq(d1, g[2] - g[18]))
ck.check("I (1): g14 = h n - d, g19 = G - n, g9 = A0 (H - G), g18 = 1 - d, g6 = h G/T, g22 = (d + A0 (H - n))/X",
         feq(g[14], h1 * n1 - d1) and feq(g[19], G1 - n1) and feq(g[9], A0 * (H1 - G1)) and feq(g[18], 1 - d1)
         and feq(g[6], h1 * G1 / T1) and feq(g[22], (d1 + A0 * (H1 - n1)) / X1))
ck.check("I (2): g20 = A0 h (T X (1+H) - H |zeta-1|^2)/(T |zeta-A0|^2) - d",
         feq(g[20], A0 * h1 * (T1 * X1 * (1 + H1) - H1 * (ze - 1) * (zb - 1)) / (T1 * (ze - A0) * (zb - A0)) - d1))
# the two 'forms' behind (3), in a small field of the zeta-coordinates
Kz_, zt, ztb, Xz, Tz, rz, Uz = field("zeta zetab X T r U", sp.ZZ)
A0z, hz, Nz, xi2z = 1 + Xz, 1 + Tz, zt * ztb, zt + ztb
Gz, Hz = rz / Nz, Tz * Uz / ((zt - hz) * (ztb - hz))
Qzz = rz * (zt - hz) * (zt - A0z) - Tz * Xz * zt * (zt - hz - Uz)
Sz = Xz * Hz - A0z * Gz
ck.check("I: Q/(zeta(zeta-h)) = G N - X (T + h H) + S conj(zeta), S = X H - A0 G",
         feq(Qzz / (zt * (zt - hz)), Gz * Nz - Xz * (Tz + hz * Hz) + Sz * ztb))
ck.check("I: r((X-U) zeta - (hX-U))/(zeta(zeta-h)) = G (X - H (N - 2 xi + h)/T) conj(zeta) + G N H",
         feq(rz * ((Xz - Uz) * zt - (hz * Xz - Uz)) / (zt * (zt - hz)),
             Gz * (Xz - Hz * (Nz - xi2z + hz) / Tz) * ztb + Gz * Nz * Hz))
ck.info("K = -d + n conj(zeta) by the definitions of n, d (K - conj K = n (conj zeta - zeta), zeta K - conj(zeta K) = "
        "-d (zeta - conj zeta)); comparing the coefficients of 1 and conj(zeta) (zeta nonreal) gives two real equations.")
Ra, Ga, Ha, na, da, Xa, Ta, Na, xia = sp.ring("G H n d X T N xi", sp.QQ)


def defs(G_, H_, n_, d_, X_, T_):
    A0_, h_ = 1 + X_, 1 + T_
    S_ = X_ * H_ - A0_ * G_
    L_ = G_ * (d_ + H_ - A0_ * n_) + X_ * n_ * H_
    M_ = d_ * X_ * (T_ + h_ * H_)
    B0_ = T_ * n_ * S_ - G_ * H_
    Ec_ = T_ * G_ * X_ - G_ * H_ * h_ + T_ * d_ * S_ + T_ * n_ * X_ * (T_ + h_ * H_)
    EL_ = Ec_ * L_ - G_ * M_ * (T_ * n_ + H_)
    return A0_, h_, S_, L_, M_, B0_, Ec_, EL_


A0a, ha, Sa, La, Ma, B0a, Eca, ELa = defs(Ga, Ha, na, da, Xa, Ta)
C1a = Ga * Na - Xa * (Ta + ha * Ha)
eq1 = da * C1a + na * Sa * Na + Ga * Na * Ha                       # constant coefficient
TE1 = Ta * Ga * Xa - Ga * Ha * (Na - 2 * xia + ha)
Teq2 = Ta * (na * C1a - da * Sa + 2 * xia * na * Sa) - TE1       # T times the coefficient of conj(zeta)
ck.check("I (3): N L - M and 2 xi B0 - (Ec - G N (T n + H)) are exactly the two coefficient equations of "
         "(-d + n conj zeta)(G N - X(T+hH) + S conj zeta) = E1 conj zeta + G N H (polynomial identities in G,H,n,d,X,T,N,xi)",
         Na * La - Ma == eq1 and 2 * xia * B0a - (Eca - Ga * Na * (Ta * na + Ha)) == Teq2)
Fa = ha * ELa - B0a * (Ma + ha ** 2 * La)
ck.check("I (4): F + L B0 |zeta-h|^2 = h (E_L - 2 xi L B0) - B0 (M - L N), so (3) gives F = -L B0 |zeta-h|^2",
         Fa + La * B0a * (Na - 2 * ha * xia + ha ** 2) == ha * (ELa - 2 * xia * La * B0a) - B0a * (Ma - La * Na))
# exact spot check of (3) on the chart itself at rational points (v, vb independent rationals)
import random

rng = random.Random(77)
spot = True
for _ in range(3):
    pt = {s_: sp.Rational(rng.randint(-9, 9), rng.randint(1, 5)) + sp.Rational(1, 7) for s_ in (sp.Symbol("v"), sp.Symbol("vb"))}
    pt.update({sp.Symbol("x"): sp.Rational(rng.randint(2, 9), 3) + 1, sp.Symbol("t"): -sp.Rational(rng.randint(1, 9), 4),
               sp.Symbol("r"): sp.Rational(rng.randint(1, 9), 5), sp.Symbol("u"): -sp.Rational(rng.randint(1, 9), 2)})
    ev = lambda e_: sp.Rational(C.fexpr(e_).subs(pt)) if hasattr(e_, "numer") else sp.Rational(e_)
    Gv, Hv, nv, dv, Xv, Tv, Nv, xiv = ev(G1), ev(H1), ev(n1), ev(d1), ev(X1), ev(T1), ev(N), ev(ze + zb) / 2
    A0v, hv, Sv, Lv, Mv, B0v, Ecv, ELv = defs(Gv, Hv, nv, dv, Xv, Tv)
    spot &= (Nv * Lv == Mv) and (2 * xiv * B0v == Ecv - Gv * Nv * (Tv * nv + Hv)) and (2 * xiv * Lv * B0v == ELv)
ck.check("I (3): N L = M, 2 xi B0 = Ec - G N (T n + H), 2 xi L B0 = E_L at 3 exact rational chart points", spot)

Kf, a, c, d, e, X, T = field("a c d e X T", sp.QQ)
hh, AA, kk = 1 + T, 1 + X, T * X - 1
n = (a + d) / hh
G = c + n
H = G + e
A0f, hf, Sf, Lf, Mf, B0f, Ecf, ELf = defs(G, H, n, d, X, T)
ck.check("I L.1: L = G d + G c + e (G + X n) and B0 = e (k n - c) - G (c + h n)  (variables a, c, d, e, X, T; "
         "n = (a+d)/h, G = c + n, H = G + e, k = T X - 1)",
         feq(Lf, G * d + G * c + e * (G + X * n)) and feq(B0f, e * (kk * n - c) - G * (c + hh * n)))
Ff = hf * ELf - B0f * (Mf + hf ** 2 * Lf)
Ncal = (AA ** 2 * a ** 3 + hh * AA * (X + 2) * c * a ** 2
        + a ** 2 * ((1 - d) * X * hh * AA + d * (T * X ** 2 + 3 * X ** 2 + 4 * X + 2)) + hh ** 2 * AA * c ** 2 * a
        + c * a * ((1 - d) * X * hh * (T + X + 2) + 2 * d * hh * AA ** 2)
        + a * (d * (1 - d) * X * (T + 2) * AA + d ** 2 * (T * X ** 2 + 3 * X ** 2 + 3 * X + 1))
        + X * hh ** 2 * c ** 2 + X * d * hh * (X + 2) * c + X * d ** 2 * AA)
F1 = Ff.diff(e)
F0 = Ff - e * F1
qf = kk * n - c
estar = G * (c + hh * n) / qf
ck.check("I L.1: F is affine in e with dF/de = (T/h) Ncal (Ncal as displayed in the document)",
         fzero(F1.diff(e)) and fzero(F0.diff(e)) and feq(F1, (T / hh) * Ncal))
ck.check("I L.1: F(e*) = T h G^2 (A0 a + X d)(A0 G a + X (a + d + c))/q, e* = G (c + h n)/q, q = k n - c",
         feq(F0 + estar * F1, T * hh * G ** 2 * (AA * a + X * d) * (AA * G * a + X * (a + d + c)) / qf))
# sign-explicit form of Ncal: nonnegative coefficients in (a, c, X, T, d, d1) with d1 = 1 - d
sa, sc, sd, sd1, sX, sT = sp.symbols("a c d d1 X T")
shh, sAA = 1 + sT, 1 + sX
Nexp = (sAA ** 2 * sa ** 3 + shh * sAA * (sX + 2) * sc * sa ** 2
        + sa ** 2 * (sd1 * sX * shh * sAA + sd * (sT * sX ** 2 + 3 * sX ** 2 + 4 * sX + 2)) + shh ** 2 * sAA * sc ** 2 * sa
        + sc * sa * (sd1 * sX * shh * (sT + sX + 2) + 2 * sd * shh * sAA ** 2)
        + sa * (sd * sd1 * sX * (sT + 2) * sAA + sd ** 2 * (sT * sX ** 2 + 3 * sX ** 2 + 3 * sX + 1))
        + sX * shh ** 2 * sc ** 2 + sX * sd * shh * (sX + 2) * sc + sX * sd ** 2 * sAA)
ck.check("I L.1: Ncal has nonnegative coefficients as a polynomial in a, c, X, T, d, d1 = 1 - d, and contains X d^2",
         C.poly_coeffs_nonneg(Nexp, [sa, sc, sd, sd1, sX, sT])
         and sp.Poly(sp.expand(Nexp), sa, sc, sd, sd1, sX, sT).coeff_monomial(sX * sd ** 2) > 0)

Kn, n_, c_, e_, d_, X_, T_ = field("n c e d X T", sp.QQ)
Gn = n_ + c_
Hn = Gn + e_
A0n, hn, Sn, Ln, Mn, B0n, Ecn, ELn = defs(Gn, Hn, n_, d_, X_, T_)
kn = T_ * X_ - 1


def Ppoly(cc, ee, nn, dd, XX, TT):
    A0_, h_, k_ = 1 + XX, 1 + TT, TT * XX - 1
    p2 = 2 * A0_ * h_ * (ee + nn) + (k_ + 2 * h_) * dd
    p1 = (A0_ * h_ * (ee * ee + (2 - k_) * nn * ee + nn * nn) + (k_ + 2 * h_) * dd * ee + TT * A0_ * dd * dd
          - (k_ + 1) * dd + (2 * h_ - TT * k_) * dd * nn)
    p0 = -nn * (A0_ * ee + dd) * (h_ * k_ * (ee + nn) + TT * (k_ + A0_ * (1 - dd)))
    return A0_ * h_ * cc ** 3 + p2 * cc ** 2 + p1 * cc + p0


Rr = Gn * (A0n * c_ + d_) + e_ * (A0n * A0n * n_ + A0n * c_ - X_ * d_)
lhs6 = (A0n * hn * B0n * (T_ * X_ * (1 + Hn) * Ln - Hn * (Mn + Ln)) + A0n * hn * Hn * ELn
        - d_ * T_ * B0n * (Mn + A0n ** 2 * Ln) + A0n * d_ * T_ * ELn)
ck.check("I (6): L B0 F20 = -T R P(c), where L B0 F20 is expanded with N L = M and 2 xi L B0 = E_L "
         "(F20 = A0 h (T X (1+H) - H (N - 2 xi + 1)) - d T (N - 2 A0 xi + A0^2)); polynomial identity in n, c, e, d, X, T",
         feq(lhs6, -T_ * Rr * Ppoly(c_, e_, n_, d_, X_, T_)))
Kg2, a2, c2, d2, X2, T2 = field("a c d X T", sp.QQ)
h2, A2, k2 = 1 + T2, 1 + X2, T2 * X2 - 1
n2 = (a2 + d2) / h2
G2 = c2 + n2
q2 = k2 * n2 - c2
es2 = G2 * (c2 + h2 * n2) / q2
ck.check("I Thm 1: R(e*) = A0 G^2 (A0 a + X d)/q",
         feq(G2 * (A2 * c2 + d2) + es2 * (A2 * A2 * n2 + A2 * c2 - X2 * d2), A2 * G2 ** 2 * (A2 * a2 + X2 * d2) / q2))
Kh, a3, d3, e3, X3, k3 = field("a d e X k", sp.QQ)
T3 = (k3 + 1) / X3
h3, A3, W3 = 1 + T3, 1 + X3, X3 + k3 + 1
n3 = (a3 + d3) / h3
cs3 = k3 * d3 / (A3 * h3)
R0 = X3 ** 2 + 2 * X3 * k3 + 2 * X3 + k3 + 1
Q0 = (1 - d3) * A3 ** 2 * (k3 + 1) * W3 + d3 * X3 * k3 * R0
Q1 = (1 - d3) * A3 * (k3 + 1) * W3 + 2 * d3 * k3 * (X3 ** 2 + X3 * k3 + X3 + k3 + 1)
Q2 = (1 - d3) * A3 * W3 ** 2 + d3 * k3 * R0
rhsP = (-k3 * (A3 * a3 + X3 * d3) * e3 ** 2 - (d3 * Q0 + A3 * a3 * Q1 + X3 * k3 * A3 ** 2 * a3 ** 2) / (A3 * W3) * e3
        - (d3 ** 2 * (1 - d3) * (k3 + 1) * W3 * R0 + d3 * (k3 + 1) * a3 * Q2 + X3 * d3 * k3 * A3 * (k3 + 1) * a3 ** 2)
        / (A3 * W3 ** 2))
ck.check("I Thm 1: W = X h, and P(c*) equals the displayed quadratic in e with coefficients built from Q0, Q1, Q2, R0 "
         "(c* = k d/(A0 h), T = (k+1)/X, a = h n - d)",
         feq(W3, X3 * h3) and feq(Ppoly(cs3, e3, n3, d3, X3, T3), rhsP))
Ki, c4, d4, X4, T4 = field("c d X T", sp.QQ)
ck.check("I Thm 2: (d + h c)(d + A0 c) = T X d^2 + (c + d)(A0 h c - (T X - 1) d)",
         feq((d4 + (1 + T4) * c4) * (d4 + (1 + X4) * c4),
             T4 * X4 * d4 ** 2 + (c4 + d4) * ((1 + X4) * (1 + T4) * c4 - (T4 * X4 - 1) * d4)))
ck.check("uniform p.3: x (1-t) g19 - (-t(x-1) - 1) d is A0 h g19 - (T X - 1) d",
         feq(x * (1 - t), A0 * h1) and feq(-t * (x - 1) - 1, T1 * X1 - 1))

# ---------------------------------------------------------------------------------------------------------------
ck.section("6. region II: x < 0, t > 1, r > 0, u < 0 (second branch document)")
XII, TII, UII, AII = -x, t - 1, -u, 1 - x                          # x = -X, t = T + 1, u = -U, A0 = X + 1
QzII = lambda zz: r * (zz + TII) * (zz + XII) - t * AII * zz * (zz + TII - UII)
KzII = lambda zz: -r * ((AII + UII) * zz + TII * AII - UII) / QzII(zz)
GII, HII = r / N, t * UII / ((ze + TII) * (zb + TII))
mII = (KzII(ze) - KzII(zb)) / (ze - zb)
dII = -(ze * KzII(ze) - zb * KzII(zb)) / (ze - zb)
ck.check("II: d := -Im(zeta K)/Im zeta equals g_2 - g_18 (normalised g_2 = 1)", feq(dII, g[2] - g[18]))
ck.check("II (2): g14 = T m - d, g19 = m + G, g9 = X (G + H)",
         feq(g[14], TII * mII - dII) and feq(g[19], mII + GII) and feq(g[9], XII * (GII + HII)))
ck.check("II (3): g18 = 1 - d, g22 = (X m - X H - d)/A0",
         feq(g[18], 1 - dII) and feq(g[22], (XII * mII - XII * HII - dII) / AII))
tz = Tz + 1
HzII = tz * Uz / ((zt + Tz) * (ztb + Tz))
QzzII = rz * (zt + Tz) * (zt + Xz) - tz * A0z * zt * (zt + Tz - Uz)
VzII = Tz * HzII - tz
C0z, S0z = Gz * Nz + A0z * VzII, Xz * Gz + A0z * HzII
E0z = Gz * (HzII * (xi2z + Tz - Nz) / tz - A0z)
ck.check("II L.2: Q/(zeta(zeta+T)) = r - t A0 + X r/zeta + A0 t U/(zeta+T) = C0 + S0 conj(zeta)",
         feq(QzzII / (zt * (zt + Tz)), rz - tz * A0z + Xz * rz / zt + A0z * tz * Uz / (zt + Tz))
         and feq(QzzII / (zt * (zt + Tz)), C0z + S0z * ztb))
ck.check("II L.2: -r((A0+U) zeta + T A0 - U)/(zeta(zeta+T)) = E0 conj(zeta) - G N H, E0 = G (H (2 Re zeta + T - N)/t - A0)",
         feq(-rz * ((A0z + Uz) * zt + Tz * A0z - Uz) / (zt * (zt + Tz)), E0z * ztb - Gz * Nz * HzII))
Kc, dc, mc, C0c, S0c, E0c, Gc, Nc, Hc, xic, zbc = field("d m C0 S0 E0 G N H xi zb", sp.QQ)
expr = (-dc - mc * zbc) * (C0c + S0c * zbc) - (E0c * zbc - Gc * Nc * Hc)
red = expr + mc * S0c * (zbc ** 2 - (2 * xic * zbc - Nc))              # reduce conj(zeta)^2 = 2 xi conj(zeta) - N
c1 = red.diff(zbc)
ck.check("II L.2 (4): the constant coefficient of (-d - m conj zeta)(C0 + S0 conj zeta) - (E0 conj zeta - G N H) is "
         "-d C0 + m S0 N + G N H, i.e. d C0 = N (m S0 + G H)",
         fzero(c1.diff(zbc)) and feq(red - c1 * zbc, -dc * C0c + mc * S0c * Nc + Gc * Nc * Hc))
spot = True
for _ in range(3):
    pt = {s_: sp.Rational(rng.randint(-9, 9), rng.randint(1, 5)) + sp.Rational(1, 7) for s_ in (sp.Symbol("v"), sp.Symbol("vb"))}
    pt.update({sp.Symbol("x"): -sp.Rational(rng.randint(1, 9), 3), sp.Symbol("t"): 1 + sp.Rational(rng.randint(1, 9), 4),
               sp.Symbol("r"): sp.Rational(rng.randint(1, 9), 5), sp.Symbol("u"): -sp.Rational(rng.randint(1, 9), 2)})
    ev = lambda e_: sp.Rational(C.fexpr(e_).subs(pt)) if hasattr(e_, "numer") else sp.Rational(e_)
    Gv, Hv, mv, dv, Xv, Tv, Nv, tv, Av = ev(GII), ev(HII), ev(mII), ev(dII), ev(XII), ev(TII), ev(N), ev(t), ev(AII)
    spot &= (dv * (Gv * Nv + Av * (Tv * Hv - tv)) == Nv * (mv * (Xv * Gv + Av * Hv) + Gv * Hv))
ck.check("II L.2 (4): d C0 = N (m S0 + G H) at 3 exact rational chart points", spot)
Kb, Gb, Hb, mb, db, Nb, Xb, Tb = field("G H m d N X T", sp.QQ)
tb, Ab = Tb + 1, Xb + 1
Vb = Tb * Hb - tb
C0b, S0b = Gb * Nb + Ab * Vb, Xb * Gb + Ab * Hb
ck.check("II Thm 1 (5): (X m - X H) C0 - N (m S0 + G H) = A0 (m (X V - N H) - H (G N + X V))",
         feq((Xb * mb - Xb * Hb) * C0b - Nb * (mb * S0b + Gb * Hb),
             Ab * (mb * (Xb * Vb - Nb * Hb) - Hb * (Gb * Nb + Xb * Vb))))
aB, BB, cB = Tb * mb - db, Xb * (Gb + Hb), mb + Gb
L0 = tb * mb + Gb
D0 = Tb * mb * (mb + Gb) - L0 * db
ck.check("II Thm 1: T m (G + H) - L0 = G (T m - 1) + m (T H - t)",
         feq(Tb * mb * (Gb + Hb) - L0, Gb * (Tb * mb - 1) + mb * (Tb * Hb - tb)))
ck.check("II Thm 1: a B (c - d) - d^2 (a + B + c + d) = X (G + H) D0 - L0 d^2",
         feq(aB * BB * (cB - db) - db ** 2 * (aB + BB + cB + db), Xb * (Gb + Hb) * D0 - L0 * db ** 2))
ck.check("II Thm 1: (G + H) D0 - L0 d (m - H) = (m + G)(T m (G + H) - L0 d)",
         feq((Gb + Hb) * D0 - L0 * db * (mb - Hb), (mb + Gb) * (Tb * mb * (Gb + Hb) - L0 * db)))
ck.check("II Thm 1: D0 - (T m (G+H) - L0) = T m (m - H) + L0 (1 - d)  (so m > H, d <= 1 give D0 > 0)",
         feq(D0 - (Tb * mb * (Gb + Hb) - L0), Tb * mb * (mb - Hb) + L0 * (1 - db)))

# ---------------------------------------------------------------------------------------------------------------
ck.section("7. NUMERIC falsification attempts (floating point; evidence only)")
import numpy as np

npr = np.random.default_rng(20261126)
TOL = 1e-9
CH = 500000


def sample(n, sx=None, st=None, sr=None, su=None):
    def sgn_(s_):
        return npr.choice([-1.0, 1.0], n) if s_ is None else s_ * np.ones(n)

    def mag():
        return np.exp(npr.normal(0, 1.1, n))

    return (sgn_(sx) * mag(), sgn_(st) * mag(), sgn_(sr) * mag(), sgn_(su) * mag(),
            npr.normal(0, 1.5, n) + 1j * npr.normal(0, 1.5, n))


def gapsN(xv, tv, rv, uv, vv, mu):
    qv = 1 + tv * (vv - 1)
    Yv = xv - qv
    Qv = uv * tv * (vv - 1) * (vv - xv) + Yv * (vv * (rv - tv) + tv * xv)
    zv = rv * vv * Yv / Qv
    wv = zv * (1 - tv * (vv - xv) / rv)
    kap = 1j * mu / ((zv - wv) * (vv - 1))
    gN = C.numeric_gaps([zv, 0.0 * vv, 1.0 + 0.0 * vv, xv + 0.0 * vv], [wv, None, vv, qv], kap)
    sc = np.abs(gN).max(axis=1)
    return gN / sc[:, None]


# ---- global sampling: every supporting point with d > 0 must lie in region I or II and satisfy its capacity
st_ = dict(contact=0.0, nsupp=0, nhard=0, nout=0, nI=0, nII=0, socI=np.inf, linI=np.inf, capII=np.inf)
for _ in range(6):
    xv, tv, rv, uv, vv = sample(CH)
    mu = npr.choice([-1.0, 1.0], CH)
    gN = gapsN(xv, tv, rv, uv, vv, mu)
    ok = np.isfinite(gN).all(axis=1)
    st_["contact"] = max(st_["contact"], np.abs(gN[ok][:, T10]).max())
    supp = ok & (gN.min(axis=1) >= -TOL)
    dN = gN[:, 2] - gN[:, 18]
    hard = supp & (dN > 1e-7)
    regI = (xv > 1) & (tv < 0) & (rv > 0) & (uv < 0)
    regII = (xv < 0) & (tv > 1) & (rv > 0) & (uv < 0)
    st_["nsupp"] += int(supp.sum())
    st_["nhard"] += int(hard.sum())
    st_["nout"] += int((hard & ~(regI | regII)).sum())
    mI, mII_ = hard & regI, hard & regII
    st_["nI"] += int(mI.sum())
    st_["nII"] += int(mII_.sum())
    if mI.any():
        st_["socI"] = min(st_["socI"], (gN[:, 6] * gN[:, 22] - dN ** 2)[mI].min())
        st_["linI"] = min(st_["linI"], ((xv * (1 - tv)) * gN[:, 19] - (-tv * (xv - 1) - 1) * dN)[mI].min())
    if mII_.any():
        st_["capII"] = min(st_["capII"], (gN[:, 14] * gN[:, 9] * (gN[:, 19] - dN)
                                          - dN ** 2 * (gN[:, 14] + gN[:, 9] + gN[:, 19] + dN))[mII_].min())
ck.check("NUMERIC global: 3,000,000 chart points (seed 20261126; x, t, r, u = random sign * lognormal(0, 1.1); v complex "
         "normal(0, 1.5); mu = +-1; gaps normalised to max 1): max |g| on T_10 = %.1e; supporting (all gaps >= -1e-9): %d; "
         "of these with d > 1e-7: %d (region I: %d, region II: %d, elsewhere: %d); min(g6 g22 - d^2) in I = %.3e, "
         "min(A0 h g19 - (TX-1) d) in I = %.3e, min capacity in II = %.3e"
         % (st_["contact"], st_["nsupp"], st_["nhard"], st_["nI"], st_["nII"], st_["nout"], st_["socI"], st_["linI"], st_["capII"]),
         st_["contact"] < 1e-6 and st_["nout"] == 0 and st_["nsupp"] > 0
         and st_["socI"] > -1e-9 and st_["linI"] > -1e-9 and st_["capII"] > -1e-12)


def region_points(region):
    if region == "I":
        xv, tv, rv, uv, vv = sample(CH, 1, -1, 1, -1)
        xv = 1 + xv
    else:
        xv, tv, rv, uv, vv = sample(CH, -1, 1, 1, -1)
        tv = 1 + tv
    mu = np.sign(xv - 1) * np.sign(vv.imag)                 # rho = mu Im v has the sign of y = x - 1 (g_2 > 0)
    return xv, tv, rv, uv, vv, gapsN(xv, tv, rv, uv, vv, mu)


nsI, socm, linm, nhyp, b0max = 0, np.inf, np.inf, 0, -np.inf
for _ in range(8):
    xv, tv, rv, uv, vv, gN = region_points("I")
    dN = gN[:, 2] - gN[:, 18]
    supp = (gN.min(axis=1) >= -TOL) & (dN > 1e-7)
    nsI += int(supp.sum())
    if supp.any():
        socm = min(socm, (gN[:, 6] * gN[:, 22] - dN ** 2)[supp].min())
        linm = min(linm, ((xv * (1 - tv)) * gN[:, 19] - (-tv * (xv - 1) - 1) * dN)[supp].min())
    # Lemma 1 of the first branch under its own hypotheses only: g14, g19, g9, g18 >= 0, d > 0  =>  B0 < 0
    g2 = gN[:, 2]
    hyp = (g2 > 1e-6) & (gN[:, 14] >= 0) & (gN[:, 19] >= 0) & (gN[:, 9] >= 0) & (gN[:, 18] >= 0) & (dN > 1e-7 * g2)
    Xn, Tn = xv - 1, -tv
    gn = gN / g2[:, None]
    Gn_ = Tn * gn[:, 6] / (1 + Tn)
    nn_ = Gn_ - gn[:, 19]
    Hn_ = Gn_ + gn[:, 9] / (1 + Xn)
    B0N = Tn * nn_ * (Xn * Hn_ - (1 + Xn) * Gn_) - Gn_ * Hn_
    scaleB = Tn * np.abs(nn_) * (Xn * Hn_ + (1 + Xn) * Gn_) + Gn_ * Hn_
    nhyp += int(hyp.sum())
    if hyp.any():
        b0max = max(b0max, (B0N / scaleB)[hyp].max())
ck.check("NUMERIC region I: 4,000,000 points with x > 1, t < 0, r > 0, u < 0: %d supporting with d > 1e-7; "
         "min(g6 g22 - d^2) = %.3e, min(A0 h g19 - (TX-1) d) = %.3e (tolerance -1e-9)" % (nsI, socm, linm),
         nsI > 50 and socm > -1e-9 and linm > -1e-9)
ck.check("NUMERIC first branch Lemma 1 under its own hypotheses only (g14, g19, g9, g18 >= 0, d > 0; %d of these "
         "points): max of B0/(sum of |terms|) = %.3e (claim: B0 < 0)" % (nhyp, b0max), nhyp > 1000 and b0max < 1e-9)

nsII, capm, nhyp, relm, posm = 0, np.inf, 0, np.inf, np.inf
for _ in range(16):
    xv, tv, rv, uv, vv, gN = region_points("II")
    dN = gN[:, 2] - gN[:, 18]
    supp = (gN.min(axis=1) >= -TOL) & (dN > 1e-7)
    cap = gN[:, 14] * gN[:, 9] * (gN[:, 19] - dN) - dN ** 2 * (gN[:, 14] + gN[:, 9] + gN[:, 19] + dN)
    nsII += int(supp.sum())
    if supp.any():
        capm = min(capm, cap[supp].min())
    hyp = (gN[:, 2] > 1e-6) & (dN > 1e-7 * gN[:, 2]) & (gN[:, 18] >= 0) & (gN[:, 22] >= 0)
    relcap = cap / (np.abs(gN[:, 14] * gN[:, 9] * gN[:, 19]) + np.abs(gN[:, 14] * gN[:, 9]) * dN
                    + dN ** 2 * (np.abs(gN[:, 14]) + np.abs(gN[:, 9]) + np.abs(gN[:, 19]) + dN))
    nhyp += int(hyp.sum())
    if hyp.any():
        relm = min(relm, relcap[hyp].min())
        posm = min(posm, gN[hyp][:, [14, 9, 19]].min())
ck.check("NUMERIC region II: 8,000,000 points with x < 0, t > 1, r > 0, u < 0: %d supporting with d > 1e-7; "
         "min of a B (c-d) - d^2 (a+B+c+d) = %.3e (claim: > 0)" % (nsII, capm), nsII > 50 and capm > -1e-12)
ck.check("NUMERIC second branch Theorem 1 under its own hypotheses only (d > 0, g18 >= 0, g22 >= 0; %d points): "
         "min relative capacity = %.3e and min(g14, g9, g19) = %.3e (claims: all > 0)" % (nhyp, relm, posm),
         nhyp > 50 and relm > -1e-9 and posm > -1e-9)


def excluded(nchunks, sx, st, su, cond):
    bad = 0
    for _ in range(nchunks):
        xv, tv, rv, uv, vv = sample(CH, sx, st, 1, su)
        if st > 0:
            tv = 1 + tv
        if sx > 0:
            xv = 1 + xv
        mu = np.sign(xv - 1) * np.sign(vv.imag)
        gN = gapsN(xv, tv, rv, uv, vv, mu)
        dN = gN[:, 2] - gN[:, 18]
        bad += int(cond(gN, dN).sum())
    return bad


nb1 = excluded(4, -1, -1, 1, lambda gN, dN: (dN > 1e-7) & (gN[:, 14] >= -TOL) & (gN[:, 22] >= -TOL))
nb2 = excluded(4, -1, -1, -1, lambda gN, dN: (dN > 1e-7) & (gN[:, 14] >= -TOL))
nb3 = excluded(4, -1, 1, 1, lambda gN, dN: (dN > 1e-7) & (gN[:, 9] >= -TOL) & (gN[:, 10] > 1e-7))
nb4 = excluded(4, 1, 1, None, lambda gN, dN: (gN[:, 9] >= -TOL) & (gN[:, 10] >= -TOL) & (dN > 1e-7))
ck.check("NUMERIC exclusions, 2,000,000 points each, r > 0, sign of rho fixed by g_2 > 0: "
         "(x<0, t<0, u>0): d > 0, g14 >= 0, g22 >= 0 together: %d; (x<0, t<0, u<0): d > 0 and g14 >= 0: %d; "
         "(x<0, t>1, u>0): d > 0, g9 >= 0, g10 > 0: %d; (x>1, t>1): g9 >= 0, g10 >= 0, d > 0: %d  (claims: all 0)"
         % (nb1, nb2, nb3, nb4), nb1 == nb2 == nb3 == nb4 == 0)
# real v: the gap vector vanishes
xv, tv, rv, uv, vv = sample(200000)
vr = vv.real + 0j
qv = 1 + tv * (vr - 1)
Yv = xv - qv
Qv = uv * tv * (vr - 1) * (vr - xv) + Yv * (vr * (rv - tv) + tv * xv)
zv = rv * vr * Yv / Qv
wv = zv * (1 - tv * (vr - xv) / rv)
pts_ = ([zv, 0 * vr, 1 + 0 * vr, xv + 0 * vr], [wv, None, vr, qv])
gR = C.numeric_gaps(pts_[0], pts_[1], 1j / ((zv - wv) * (vr - 1)))
Rscale = np.abs(C.numeric_gaps(pts_[0], pts_[1], 1 / ((zv - wv) * (vr - 1)))).max(axis=1)
ratio = np.nanmax(np.abs(gR).max(axis=1) / Rscale)
ck.check("NUMERIC real v: 200,000 chart points with v real: max |gap| / max |R_j| = %.1e (the gap vector is zero)" % ratio,
         ratio < 1e-9)
b1 = C.numeric_atom_attack("SOC", (1, 11, 28), 20000, 6, 1002)
b2 = C.numeric_atom_attack("FOUR", (4, 6, 2, 28), 20000, 6, 1003)
ck.check("NUMERIC atoms on SO(4) (20000 random points + 6 Nelder-Mead runs each, seeds 1002, 1003): min of "
         "Phi_28 + 2|x_1 x_11| = %.3e; min of Phi_28 + Phi_2 + 2(|x_4 x_6| + |x_4 x_2| + |x_6 x_2|) = %.3e (tolerance -1e-9)"
         % (b1, b2), b1 > -1e-9 and b2 > -1e-9)

ck.finish()
