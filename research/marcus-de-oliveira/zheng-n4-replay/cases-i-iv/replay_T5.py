"""Independent replay of the proof of Case T_5 = {1,2,3,4,6,8}.

Documents replayed: supplement, proofs/manin_transfer_20261021/
    final_referee_moc4_T5_uniform_20261022.pdf           (5 pages; read, nothing executed)
    final_referee_T5_four_budget_capacity_20261022.pdf   (3 pages; read)

Run:  MOC4_SUPP=/path/to/supplement python replay_T5.py
(MOC4_SUPP is used only to read one JSON data file for the facet-numbering comparison.)
Exit status 0 iff every check passes.

Exact parts: sympy rational-function fields / polynomial rings, Fractions.
Section NUMERIC: numpy floating point, falsification attempt only.
"""
import random
import sys

sys.dont_write_bytecode = True

import replay_T3_common as C
from sympy.polys.domains import QQ
from sympy.polys.fields import field

T5 = (1, 2, 3, 4, 6, 8)
ck = C.Checker("replay_T5: Case T_5 = {1,2,3,4,6,8} (final_referee_moc4_T5_uniform_20261022.pdf)")
V = C.FVAL

# ---------------------------------------------------------------------------------------
# A. objects
# ---------------------------------------------------------------------------------------
C.crosscheck_numbering(ck)

ck.section("A. contact set, generators, A0, the three replacements (document pp. 1-2)")
ck.check("document's T_5 = {1,2,3,4,6,8} is literally the paper's T_5 (identity of G); reference permutation pi_1 is in T_5",
         set(T5) == {1, 2, 3, 4, 6, 8})
gens = C.generators(T5)
ck.check("generators vanishing on T_5: {0,1,2,3,4,5,6,9,10,53,54}", gens == [0, 1, 2, 3, 4, 5, 6, 9, 10, 53, 54], gens)
ck.check("only f54 is not a minor vector", [j for j in gens if not C.IS_MINOR[j]] == [54])
C.verify_contact_cone(ck, T5)
ck.check("span of the generators has dimension 8", C.rank([V[j] for j in gens]) == 8)
A0_of = lambda g: g[16] - g[0]
ck.check("A0 = g16 - g0 is -1 on f54, +1 on f1, f4, f9 and 0 on the other generators",
         {j: A0_of(V[j]) for j in gens} == {0: 0, 1: 1, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 9: 1, 10: 0, 53: 0, 54: -1})
vec = lambda *js: tuple(sum(V[j][s] for j in js) for s in range(24))
ck.check("f54 + f1 = f0 + f6 + f53", vec(54, 1) == vec(0, 6, 53))
ck.check("f54 + f4 = f2 + f3 + f6 + f53", vec(54, 4) == vec(2, 3, 6, 53))
ck.check("f54 + f9 = f2 + f5 + f10 + f53", vec(54, 9) == vec(2, 5, 10, 53))
ck.info("PROOF-READ: g = sum lam_j f_j, lam >= 0 (Lemma 4.2); A0 = lam1 + lam4 + lam9 - lam54 >= 0 lets lam54 be split")
ck.info("  against lam1, lam4, lam9 and replaced by minors: A0 >= 0 => minor decomposition.")
ck.check("f54 is positive at the permutation indices 0, 9, 11 (so A0 < 0 forces g0, g9, g11 > 0)",
         all(V[54][s] > 0 for s in (0, 9, 11)), [V[54][s] for s in (0, 9, 11)])
d_of = lambda g: g[0] - g[16]
ck.check("(3): g = g22 f0 + g13 f2 + g12 f3 + g10 f5 + g17 f6 + g7 f10 + g16 f53 + d f54, d = g0 - g16, on the span",
         all(tuple(g[22] * V[0][s] + g[13] * V[2][s] + g[12] * V[3][s] + g[10] * V[5][s] + g[17] * V[6][s]
                   + g[7] * V[10][s] + g[16] * V[53][s] + d_of(g) * V[54][s] for s in range(24)) == g
             for g in (V[k] for k in gens)))
ck.info("(4) is the SOC condition for (f5: g10, f2: g13; f54: d); (5) is the four-coefficient condition for")
ck.info("(f0: g22, f6: g17; f10: g7; f54: d).  The remaining terms of (3) are minors with gap coefficients.")

# ---------------------------------------------------------------------------------------
# B. the universal inequalities used
# ---------------------------------------------------------------------------------------
ck.section("B. SOC member a f5 + b f2 + d f54 (document p. 1, eq. (1))")
rng = random.Random(20261125)
g = C.orbit_check(ck, "SOC", (1, 9, 22), (5, 2, 54), expected_orbit=144)
C.covariance_check(ck, "SOC (1,9;22) -> (5,2;54)", g, [(1, 5), (9, 2), (22, 54)], rng)
C.replay_SOC_representative(ck)
ck.section("B. four-coefficient member a f0 + b f6 + c f10 + d f54 (document p. 1, eq. (2)): the representative itself")
C.orbit_check(ck, "four-coefficient", (0, 6, 10, 54), (0, 6, 10, 54), expected_orbit=288)
C.replay_four_coefficient_note(ck)

# ---------------------------------------------------------------------------------------
# C. the chart
# ---------------------------------------------------------------------------------------
ck.section("C. normalisation and contact ratios (document p. 2)")
C.infinity_convention_check(ck, "C", 2)


def products(A, Cc):
    """F_sigma = - prod_{sigma(i) != 2} (A_i - C_sigma(i)).  The sign -1 is the bracket of a
    finite point with the point at infinity; with it the identities (7) hold as printed."""
    one = A[1] * 0 + 1
    F = []
    for s in C.PERMS:
        p = -one
        for i in range(4):
            if Cc[s[i]] is None:
                continue
            p = p * (A[i] - Cc[s[i]])
        F.append(p)
    return F


K0, z, w, v, x, t = field("z,w,v,x,t", QQ)
F = products([z, K0(1), K0(0), x], [w, v, None, t * v])
r_def = z * (v - w) * (x - t * v) / (v * (z - w))
u_def = (z - v) * (x - t * v) / (v * (z - w))
ck.check("contact ratios R_sigma = (F_sigma - F_1)/(v(z-w)) on T_5 = (1,2,3,4,6,8): 0, x-t, t(x-1), 1-t, x-t-r-u, x-t-r "
         "(identities in complex z, w, v, x, t; r, u as defined on p. 2)",
         [(F[s] - F[1]) / (v * (z - w)) for s in T5]
         == [K0(0), x - t, t * (x - 1), 1 - t, x - t - r_def - u_def, x - t - r_def])
ck.check("the four transposition differences of p. 2, up to one common sign: F_3-F_1, F_4-F_1, F_3-F_2, F_4-F_2 = "
         "+-(-t v (z-w)(x-1), v (t-1)(z-w), -v x (t-1)(z-w), v (x-1)(z-w))",
         [F[3] - F[1], F[4] - F[1], F[3] - F[2], F[4] - F[2]]
         == [-(-t * v * (z - w) * (x - 1)), -(v * (t - 1) * (z - w)), -(-v * x * (t - 1) * (z - w)), -(v * (x - 1) * (z - w))])
ck.info("PROOF-READ: kappa = i m/(v(z-w)), m in C*: contacts 2,3,4 <=> m(x-t), m t(x-1), m(1-t) real; the first minus")
ck.info("  the third gives m(x-1) real, so t real (x != 1), m real (t != 1), x real; contacts 8, 6 <=> r, u real.")
ck.check("exhaustiveness: z (x - (t+u) v) = r v and u w = x - r - t v identically for r, u defined from z, w",
         z * (x - (t + u_def) * v) == r_def * v and u_def * w == x - r_def - t * v)

K, x, t, u, r, v, vb = field("x,t,u,r,v,vb", QQ)


def chart(x, t, u, r, v):
    return r * v / (x - (t + u) * v), (x - r - t * v) / u


def gdoc(x, t, u, r, v, vb):
    """gamma_sigma with g_sigma = rho * gamma_sigma (rho = mu Im v, kappa = i mu/(v(z-w)))."""
    out = []
    for vv in (v, vb):
        z, w = chart(x, t, u, r, vv)
        F = products([z, vv * 0 + 1, vv * 0, x], [w, vv, None, t * vv])
        out.append([(f - F[1]) / (vv * (z - w)) for f in F])
    return [-(p - q) / (v - vb) for p, q in zip(*out)]


zc, wc = chart(x, t, u, r, v)
ck.check("the chart (6) inverts the definitions of r and u",
         zc * (v - wc) * (x - t * v) / (v * (zc - wc)) == r and (zc - v) * (x - t * v) / (v * (zc - wc)) == u)
gam = gdoc(x, t, u, r, v, vb)
ck.check("on the chart the gaps vanish identically exactly on T_5", [s for s in range(24) if gam[s] == 0] == list(T5))
ck.check("the chart's gap vector satisfies (3) in all 24 entries (so it lies in the span)",
         [gam[22] * V[0][s] + gam[13] * V[2][s] + gam[12] * V[3][s] + gam[10] * V[5][s] + gam[17] * V[6][s]
          + gam[7] * V[10][s] + gam[16] * V[53][s] + d_of(gam) * V[54][s] for s in range(24)] == gam)
z3, w3 = chart(x, t, -t, x, v)
ck.check("F = 0 with nonreal v needs t + u = 0, r = x, and then z = v: collision", z3 == v)

ck.section("C. identities (7) (document p. 3); g_sigma = rho * gamma_sigma, sigma = rho/E")
N = v * vb
En = (x - t * v) * (x - t * vb)
Fn = (x - r - (t + u) * v) * (x - r - (t + u) * vb)
sig = 1 / En                                  # sigma / rho
L = x - r - t - u
D0 = r * t + u * x
J = r + (t + u - 1) * x
H = -r * t + t * (x - 1) - u
W = t * (x - 1) - r - u
A0 = A0_of(gam)
ck.check("g0 = -rho x/N", gam[0] == -x / N)
ck.check("g7 = -sigma t x (r+u)", gam[7] == -sig * t * x * (r + u))
ck.check("g9 = -sigma t x (r + u x)", gam[9] == -sig * t * x * (r + u * x))
ck.check("g10 = -sigma r t (x-1)", gam[10] == -sig * r * t * (x - 1))
ck.check("x g11 = (x-1) g9", x * gam[11] == (x - 1) * gam[9])
ck.check("g12 = -rho r u L/F", gam[12] == -r * u * L / Fn)
ck.check("g17 = A0 - sigma t u x (x-1)", gam[17] == A0 - sig * t * u * x * (x - 1))
ck.check("H g13 = t L A0 + sigma t u x (t-1) W", H * gam[13] == t * L * A0 + sig * t * u * x * (t - 1) * W)
ck.check("H g22 = (t-1)(x-r) A0 + sigma r t x (t-1)(x-r-u-1)",
         H * gam[22] == (t - 1) * (x - r) * A0 + sig * r * t * x * (t - 1) * (x - r - u - 1))
ck.check("D0 A0 = sigma r u [x t (t-1)(x-1) - J H E/F]",
         D0 * A0 == sig * r * u * (x * t * (t - 1) * (x - 1) - J * H * En / Fn))
ck.check("H = t L + (t-1)(t+u) and W = H + r(t-1)", H == t * L + (t - 1) * (t + u) and W == H + r * (t - 1))

ck.section("C. negative branch x = -X, t = -T, r = -R, u = U, sigma = 1 (document pp. 3-4)")
KN, X, T, R, U, d = field("X,T,R,U,d", QQ)
Ln = R + T - X - U
Hn = T * (X + 1 - R) - U
Wn = R - U + T * (X + 1)
Bn = T * U * X * (X + 1)
x_, t_, r_, u_ = -X, -T, -R, U
ck.check("specialisation: L, H, W of (7) become R+T-X-U, T(X+1-R)-U, R-U+T(X+1)",
         x_ - r_ - t_ - u_ == Ln and -r_ * t_ + t_ * (x_ - 1) - u_ == Hn and t_ * (x_ - 1) - r_ - u_ == Wn)
# (7) at sigma = 1, A0 = -d:
g10n = -r_ * t_ * (x_ - 1)
g17n = -d - t_ * u_ * x_ * (x_ - 1)
Hg13n = t_ * Ln * (-d) + t_ * u_ * x_ * (t_ - 1) * Wn
Hg22n = (t_ - 1) * (x_ - r_) * (-d) + r_ * t_ * x_ * (t_ - 1) * (x_ - r_ - u_ - 1)
ck.check("(8): g10 = R T (X+1), g17 = B - d, H g13 = T[L d - U X (T+1) W], H g22 = (T+1)[(R-X) d + R T X (R-X-U-1)]",
         g10n == R * T * (X + 1) and g17n == Bn - d and Hg13n == T * (Ln * d - U * X * (T + 1) * Wn)
         and Hg22n == (T + 1) * ((R - X) * d + R * T * X * (R - X - U - 1)))
ck.check("g7 = T X (R - U) and g12 = (rho/F) R U L at sigma = 1 (signs of the support conditions R >= U, L >= 0)",
         -t_ * x_ * (r_ + u_) == T * X * (R - U) and -r_ * u_ * Ln == R * U * Ln)
da = R * T * X * (X + U + 1 - R) / (R - X)
ck.check("d_a - B = T X (R + U X)(X + 1 - R)/(R - X)", da - Bn == T * X * (R + U * X) * (X + 1 - R) / (R - X))
ck.check("H - (T - U) = T (X - R)   (so R <= X would give H >= T - U >= 0)", Hn - (T - U) == T * (X - R))
Psi = lambda dd: R * T * (X + 1) * (T / Hn) * (Ln * dd - U * X * (T + 1) * Wn) - dd * dd
Qn = -(X + 1) * Hn - R * (T + 1)
Pn = R * (R - U) + (X + 1) * (R * T - (T + 1) * X)
ck.check("Psi(d) = g10 g13 - d^2 with g13 from (8)", Psi(d) == g10n * (Hg13n / Hn) - d * d)
ck.check("Psi(0) = -R T^2 (X+1) U X (T+1) W/H  (> 0 for H < 0, W > 0)", Psi(KN(0)) == -R * T ** 2 * (X + 1) * U * X * (T + 1) * Wn / Hn)
ck.check("Psi(B) = T^2 U X (R + U X)(X+1) Q/H, Q = -(X+1)H - R(T+1)", Psi(Bn) == T ** 2 * U * X * (R + U * X) * (X + 1) * Qn / Hn)
ck.check("Psi(d_a) = R T^2 X (R + U X) P/(R-X)^2, P = R(R-U) + (X+1)[R T - (T+1) X]",
         Psi(da) == R * T ** 2 * X * (R + U * X) * Pn / (R - X) ** 2)
ck.check("Q = -X(T+1)(X+1-R) - (X+1) L   (so L >= 0 gives Q <= -X(T+1)(X+1-R))", Qn == -X * (T + 1) * (X + 1 - R) - (X + 1) * Ln)
ck.check("P = X(T+1)(R-X-1) + R L        (so L >= 0 gives P >= X(T+1)(R-X-1))", Pn == X * (T + 1) * (R - X - 1) + R * Ln)

ck.section("C. positive branch x > 1, sigma = -1 (document pp. 4-5)")
KP, x, t, r, u, d, Wv = field("x,t,r,u,d,W", QQ)
# (a) all three negative is impossible: t = -T, r = -R, u = -U  (use symbols t, r, u for T, R, U here)
Tt, Rr, Uu = t, r, u
tm, rm, um = -Tt, -Rr, -Uu
Hm = -rm * tm + tm * (x - 1) - um
D0m = rm * tm + um * x
Jm = rm + (tm + um - 1) * x
Hg22m = (tm - 1) * (x - rm) * (-d) + (-1) * rm * tm * x * (tm - 1) * (x - rm - um - 1)
ck.check("t = -T, r = -R, u = -U, sigma = -1: H g22 = (T+1)[(x+R) d + R T x (x+R+U-1)], H = U - T(R+x-1)",
         Hg22m == (Tt + 1) * ((x + Rr) * d + Rr * Tt * x * (x + Rr + Uu - 1)) and Hm == Uu - Tt * (Rr + x - 1))
ck.check("D0 = R T - U x = -T(x-1)(R+x) - x H and J = -R - (T+U+1) x",
         D0m == Rr * Tt - Uu * x and D0m == -Tt * (x - 1) * (Rr + x) - x * Hm and Jm == -Rr - (Tt + Uu + 1) * x)
ck.info("(7) at sigma = -1, A0 = -d reads D0 d = r u [x t (t-1)(x-1) - J H E/F]; with t = -T etc.: d = (R U/D0)[x T (T+1)(x-1) - J H E/F].")
ck.check("sign bookkeeping: r u x t (t-1)(x-1) at t=-T, r=-R, u=-U equals R U x T (T+1)(x-1)",
         rm * um * x * tm * (tm - 1) * (x - 1) == Rr * Uu * x * Tt * (Tt + 1) * (x - 1))
# (b) t, r, u > 0
y, h = x - 1, 1 - t
Lp = x - r - t - u
Hp = -r * t + t * (x - 1) - u
Wp = t * (x - 1) - r - u
Bp = t * u * x * y
Cp = t * x * (r + u)
g10p = r * t * y
Hg13p = t * Lp * (-d) + (-1) * t * u * x * (t - 1) * Wp
Hap = (t - 1) * (x - r) * (-d) + (-1) * r * t * x * (t - 1) * (x - r - u - 1)
ck.check("(9): H g13 = t[u x h W - L d], H a = h[(x-r) d + r t x (y-r-u)], g17 = B - d, B = t u x y, g7 = C = t x (r+u), g10 = r t y",
         Hg13p == t * (u * x * h * Wp - Lp * d) and Hap == h * ((x - r) * d + r * t * x * (y - r - u))
         and (-d) - (-1) * t * u * x * (x - 1) == Bp - d and -(-1) * t * x * (r + u) == Cp and -(-1) * r * t * (x - 1) == g10p)
ck.check("H = t L - h (t+u), W = H - r h, W = t L + (t-1)(r+t+u)",
         Hp == t * Lp - h * (t + u) and Wp == Hp - r * h and Wp == t * Lp + (t - 1) * (r + t + u))
ck.info("H = 0 (sigma = -1): (7) gives D0 d = r u x t (t-1)(x-1) < 0 with D0 = r t + u x > 0: impossible.")
dap = r * t * x * (r + u - y) / (x - r)
PsiP = lambda dd: r * t * y * (t / Hp) * (u * x * h * Wp - Lp * dd) - dd * dd
Pp = r * (r + u) - y * (r * t + h * x)
Qp = y * Hp + r * h
ck.check("H < 0: Psi(d) = g10 g13 - d^2; Psi(0) = r t^2 y u x h W/H", PsiP(d) == g10p * (Hg13p / Hp) - d * d
         and PsiP(KP(0)) == r * t ** 2 * y * u * x * h * Wp / Hp)
ck.check("Psi(d_a) = -r t^2 x (r + u x) P/(x-r)^2, P = r(r+u) - y(r t + h x)", PsiP(dap) == -r * t ** 2 * x * (r + u * x) * Pp / (x - r) ** 2)
ck.check("Psi(B) = -t^2 u x (r + u x) y Q/H, Q = y H + r h", PsiP(Bp) == -t ** 2 * u * x * (r + u * x) * y * Qp / Hp)
ck.check("d_a - B = t x (r - y)(r + u x)/(x - r)", dap - Bp == t * x * (r - y) * (r + u * x) / (x - r))
ck.check("P = -h x (y - r) - r L   (L >= 0, r < y => P < 0)", Pp == -h * x * (y - r) - r * Lp)
ck.check("Q = (r - y) h (y+1) + y L (L >= 0, r >= y => Q >= 0)", Qp == (r - y) * h * (y + 1) + y * Lp)
# H > 0: new capacity
Kc = t * u * x * y * (r + u) / (r + u * x)
a_of = lambda dd: (h / Hp) * ((x - r) * dd + r * t * x * (y - r - u))
chi = lambda dd: a_of(dd) * (Kc - dd) - dd * dd
ck.check("K = B C/(B + C) = t u x y (r+u)/(r + u x)", Kc == Bp * Cp / (Bp + Cp))
ck.check("the difference in (5): a (B-d)(C-d) - d^2 (a + (B-d) + C + d) = (B + C) chi(d), chi = a (K - d) - d^2",
         a_of(d) * (Bp - d) * (Cp - d) - d * d * (a_of(d) + (Bp - d) + Cp + d) == (Bp + Cp) * chi(d))
ck.check("h (x - r) + H = L, hence the d^2-coefficient of chi is -L/H", h * (x - r) + Hp == Lp
         and chi(d) - chi(KP(0)) - d * (chi(KP(1)) - chi(KP(0)) + Lp / Hp) == -(Lp / Hp) * d * d)
ck.check("chi(0) = (h/H) r t x (y-r-u) K", chi(KP(0)) == (h / Hp) * r * t * x * (y - r - u) * Kc)
# substitute x = 1 + (r+u+W)/t
xs = 1 + (r + u + Wv) / t
# recompute the branch quantities with x replaced by xs:
def branch_quantities(xv):
    yv = xv - 1
    Lv = xv - r - t - u
    Hv = -r * t + t * (xv - 1) - u
    Wq = t * (xv - 1) - r - u
    Kv = t * u * xv * yv * (r + u) / (r + u * xv)
    av = lambda dd: (h / Hv) * ((xv - r) * dd + r * t * xv * (yv - r - u))
    chiv = lambda dd: av(dd) * (Kv - dd) - dd * dd
    dstar = u * xv * h * Wq / Lv
    return yv, Lv, Hv, Wq, chiv, dstar


yv, Lv, Hv, Wq, chiv, dstar = branch_quantities(xs)
tPstar = r * Wv ** 2 + (r + u) * (r * (2 - t) + u * (1 - t)) * Wv + (r + u) ** 2 * (1 - t) * (r + t + u)
ck.check("with x = 1 + (r+u+W)/t: W(x) = W and y - r - u = ((r+u)(1-t) + W)/t", Wq == Wv and yv - r - u == ((r + u) * (1 - t) + Wv) / t)
ck.check("chi(d*) = u x^2 (r+u) h P*/((r + u x) L) with t P* = r W^2 + (r+u)[r(2-t) + u(1-t)] W + (r+u)^2 (1-t)(r+t+u), d* = u x h W/L",
         chiv(dstar) == u * xs ** 2 * (r + u) * h * (tPstar / t) / ((r + u * xs) * Lv))
ck.info("PROOF-READ: all coefficients of t P* are > 0 for r, u, W > 0, 0 < t < 1; chi concave, chi(0) > 0, chi(d*) > 0.")

# ---------------------------------------------------------------------------------------
# C'. EXACT instances
# ---------------------------------------------------------------------------------------
ck.section("C'. EXACT instances built from the paper's definitions (finite Gaussian-rational spectra)")
from fractions import Fraction as Fr

GQ = C.GQ


def exact_T5(label, xq, tq, uq, rq, vq, branch):
    zq = rq * vq / (xq - (tq + uq) * vq)
    wq = (xq - rq - tq * vq) / uq
    a_s, b_s, alpha, gp = C.exact_paper_gaps([zq, GQ(1), GQ(0), GQ(xq)], [wq, vq, None, tq * vq], 1, (2, 3, 4))
    C.exact_instance_report(ck, label, gp, T5)
    dq = gp[0] - gp[16]
    cf = {0: gp[22], 2: gp[13], 3: gp[12], 5: gp[10], 6: gp[17], 10: gp[7], 53: gp[16], 54: dq}
    soc = gp[10] * gp[13] - dq * dq
    four = gp[22] * gp[17] * (gp[7] - dq) - dq * dq * (gp[22] + gp[17] + gp[7] + dq)
    Hq = -rq * tq + tq * (xq - 1) - uq
    ck.check("%s: d = -A0 = %s > 0 (not in the minor cone); (3) reproduces all 24 gaps with coefficients >= 0"
             % (label, C.fl(dq)), dq > 0 and all(val >= 0 for val in cf.values())
             and [sum(cf[j] * V[j][s] for j in cf) for s in range(24)] == gp)
    want = {"x<0": xq < 0 and Hq < 0 and soc >= 0, "x>1,H<0": xq > 1 and Hq < 0 and soc >= 0,
            "x>1,H>0": xq > 1 and Hq > 0 and four >= 0}[branch]
    ck.check("%s: branch %s, H = %s; SOC margin g10 g13 - d^2 = %s, four-coefficient margin = %s; the document's "
             "inequality for this branch holds" % (label, branch, Hq, C.fl(soc), C.fl(four)), want)
    Eq = (xq - tq * vq).norm2()
    rho = -gp[0] * vq.norm2() / xq
    sg = rho / Eq
    ck.check("%s: with rho := -g0 N/x, sigma = rho/E: formulas (7) for g7, g9, g10, g17 hold at this point" % label,
             gp[7] == -sg * tq * xq * (rq + uq) and gp[9] == -sg * tq * xq * (rq + uq * xq)
             and gp[10] == -sg * rq * tq * (xq - 1) and gp[17] == -dq - sg * tq * uq * xq * (xq - 1))
    return soc, four


exact_T5("instance (x,t,u,r,v) = (-3,-1,2,-4,-3/2+i/2)", Fr(-3), Fr(-1), Fr(2), Fr(-4), GQ(Fr(-3, 2), Fr(1, 2)), "x<0")
exact_T5("instance (x,t,u,r,v) = (2,1/4,3/4,1,-1-5i/4)", Fr(2), Fr(1, 4), Fr(3, 4), Fr(1), GQ(-1, Fr(-5, 4)), "x>1,H<0")
socH, fourH = exact_T5("instance (x,t,u,r,v) = (6,1/2,1/4,1/4,6+i/4)", Fr(6), Fr(1, 2), Fr(1, 4), Fr(1, 4),
                       GQ(6, Fr(1, 4)), "x>1,H>0")
ck.info("in the H > 0 instance the SOC margin is %s: %s" % (C.fl(socH), "the SOC condition (4) does not hold there, so the "
        "four-coefficient inequality is genuinely needed for this decomposition" if socH < 0 else "SOC (4) also holds there"))

# ---------------------------------------------------------------------------------------
# D. NUMERIC falsification attempts
# ---------------------------------------------------------------------------------------
ck.section("D. NUMERIC falsification attempts on the intermediate claims (numpy, seeds printed)")
import numpy as np


def sample(xs_, ts_, us_, rs_, vs_, rngn):
    """Gaps from FINITE spectra: random Moebius image of (A, C) = ((z,1,0,x), (w,v,inf,tv)),
    a_i = m(A_i), b_j = -m(C_j), F_sigma = prod (a_i + b_sigma(i)); reference pi_1; direction with
    Re kappa(F_2 - F_1) = 0."""
    zs = rs_ * vs_ / (xs_ - (ts_ + us_) * vs_)
    ws = (xs_ - rs_ - ts_ * vs_) / us_
    n = len(xs_)
    al, be, ga, de = [rngn.normal(size=n) + 1j * rngn.normal(size=n) for _ in range(4)]
    mob = lambda s: (al * s + be) / (ga * s + de)
    a = np.stack([mob(zs), mob(0 * zs + 1), mob(0 * zs), mob(xs_ + 0j)])
    b = -np.stack([mob(ws), mob(vs_), al / ga, mob(ts_ * vs_)])
    Fs = np.ones((24, n), dtype=complex)
    for k, s in enumerate(C.PERMS):
        for i in range(4):
            Fs[k] *= a[i] + b[s[i]]
    d1 = Fs[2] - Fs[1]
    kap = 1j * np.conj(d1) / np.abs(d1)
    g = (kap * (Fs - Fs[1])).real
    return g / np.max(np.abs(g), axis=0)


def run_family(name, seed, n, gen):
    rngn = np.random.default_rng(seed)
    xs_, ts_, us_, rs_, vs_ = gen(rngn, n)
    g = sample(xs_, ts_, us_, rs_, vs_, rngn)
    cerr = np.max(np.abs(g[list(T5)]), axis=0)
    good = np.isfinite(cerr) & (cerr < 1e-6)
    st = dict(n=int(good.sum()), supp=0, neg=0, posHneg=0, posHpos=0)
    worst = dict(soc_neg=np.inf, soc_pos=np.inf, four=np.inf)
    viol = []
    for sign in (1.0, -1.0):
        gg = sign * g
        idx = np.nonzero(good & (np.min(gg, axis=0) >= -1e-9))[0]
        st["supp"] += len(idx)
        for i in idx:
            hv = gg[:, i]
            dv = hv[0] - hv[16]
            if dv <= 1e-9:
                continue
            xv, tv, uv, rv = xs_[i], ts_[i], us_[i], rs_[i]
            Hv_ = -rv * tv + tv * (xv - 1) - uv
            Lv_ = xv - rv - tv - uv
            soc = hv[10] * hv[13] - dv * dv
            four = hv[22] * hv[17] * (hv[7] - dv) - dv * dv * (hv[22] + hv[17] + hv[7] + dv)
            if xv < 0:
                st["neg"] += 1
                if not (tv < 0 and rv < 0 and uv > 0 and -rv >= uv - 1e-9 and Lv_ >= -1e-9 and Hv_ < 0 and -rv > -xv):
                    viol.append(("negative branch sign/order claims fail", i, xv, tv, rv, uv))
                worst["soc_neg"] = min(worst["soc_neg"], soc)
                if soc < -1e-10:
                    viol.append(("negative branch: d^2 > g10 g13", i))
            elif xv > 1:
                if not (tv > 0 and rv > 0 and uv > 0 and tv < 1 and Lv_ >= -1e-9 and abs(Hv_) > 0):
                    viol.append(("positive branch sign claims fail", i, xv, tv, rv, uv))
                if Hv_ < 0:
                    st["posHneg"] += 1
                    worst["soc_pos"] = min(worst["soc_pos"], soc)
                    if soc < -1e-10:
                        viol.append(("positive branch H<0: d^2 > g10 g13", i))
                else:
                    st["posHpos"] += 1
                    Wv_ = tv * (xv - 1) - rv - uv
                    worst["four"] = min(worst["four"], four)
                    if not (Wv_ > 0 and Lv_ > 0):
                        viol.append(("positive branch H>0: W>0, L>0 fail", i))
                    if four < -1e-12:
                        viol.append(("positive branch H>0: four-coefficient condition fails", i))
            else:
                viol.append(("supporting source with A0<0 and 0<=x<=1", i, xv))
    ck.check("NUMERIC %s: seed %d, %d samples on the contact locus, %d supporting; with A0 < 0: %d (x<0), %d (x>1, H<0), "
             "%d (x>1, H>0); no violation" % (name, seed, st["n"], st["supp"], st["neg"], st["posHneg"], st["posHpos"]),
             not viol, viol[:3])
    ck.info("   min g10 g13 - d^2: %.3e (x<0), %.3e (x>1, H<0) | min four-coefficient margin (x>1, H>0): %.3e"
            % (worst["soc_neg"], worst["soc_pos"], worst["four"]))
    return st


def gen_broad(rngn, n):
    sg = lambda: rngn.choice([-1.0, 1.0], size=n)
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    return (sg() * mag(), sg() * mag(), sg() * mag(), sg() * mag(), rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_neg(rngn, n):
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    X_, T_, U_ = mag(), mag(), mag()
    R_ = U_ * (1 + mag() * rngn.integers(0, 2, size=n))            # R = U exactly half the time
    k = n // 4
    T_[:k] = (X_ + U_ - R_)[:k]                                     # L = 0 exactly on a quarter (kept only if T > 0)
    T_ = np.abs(T_) + 1e-12
    return (-X_, -T_, U_, -R_, 2 * rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_pos(rngn, n):
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    xs_ = 1 + mag()
    rs_ = mag()
    k = n // 5
    rs_[:k] = (xs_ - 1)[:k]                                         # r = y exactly on a fifth
    return (xs_, rngn.uniform(0, 1, size=n), mag(), rs_, 2 * rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_posH(rngn, n):
    """Region of the 'new capacity': H > 0, i.e. x = 1 + (r+u+W)/t with W > 0; Re v of the order of x .. x/t."""
    ts_ = rngn.uniform(0, 1, size=n)
    rs_ = np.exp(rngn.normal(loc=-1.0, scale=1.2, size=n))
    us_ = np.exp(rngn.normal(loc=-1.0, scale=1.2, size=n))
    Ws_ = np.exp(rngn.normal(scale=1.5, size=n))
    xs_ = 1 + (rs_ + us_ + Ws_) / ts_
    vre = xs_ * np.exp(rngn.uniform(-0.6, 1, size=n) * (np.log(1 / ts_) + 0.3))
    vim = xs_ * np.exp(rngn.normal(loc=-2.0, scale=1.5, size=n)) * rngn.choice([-1.0, 1.0], size=n)
    return (xs_, ts_, us_, rs_, vre + 1j * vim)


s1 = run_family("broad chart sample", 51001, 300000, gen_broad)
s2 = run_family("negative branch region (incl. R = U, L = 0)", 51002, 300000, gen_neg)
s3 = run_family("positive branch region (incl. r = y)", 51003, 300000, gen_pos)
s4 = run_family("positive branch, H > 0 region", 51004, 400000, gen_posH)
ck.check("NUMERIC: tests not vacuous (supporting samples with A0 < 0 in all three branches)",
         s2["neg"] > 200 and s3["posHneg"] > 200 and s4["posHpos"] > 200, (s2["neg"], s3["posHneg"], s4["posHpos"]))

ck.finish()
