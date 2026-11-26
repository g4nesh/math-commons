"""Independent exact replay of the case T_15 = {1,2,3,4,6,7} of Zheng's n = 4 proof.

Proof document replayed: final_referee_moc4_T14_uniform_20261022.pdf (3 pages; historical number "T14").
The document was read as a PDF; none of its programs was read for this script, imported or executed, and no data of the
supplement is needed (MOC4_SUPP is not used).  Run:  python replay_T15.py   (exit code 0 iff every check passes).

Notation of the document: labels 0..3, A = (z, 1, 0, x), C = -B = (w, v, infinity, t v), reference permutation pi_1,
    P_pi = prod_i [A_i, C_pi(i)]  (bracket with the point at infinity = -1),   R_pi = (P_pi - P_1) / (v (z - w)),
    kappa = i mu / (v (z - w)),   g_pi = Re(kappa (P_pi - P_1)) = -mu Im R_pi.
Complex parameters are represented by pairs of independent symbols (p, pb), (q, qb), ...; conjugation swaps them.
"""
import sys

import sympy as sp

import replay_T10_common as C

ck = C.Checker("Replay of T_15 = {1,2,3,4,6,7}  (document: final_referee_moc4_T14_uniform)")
F = C.F
T15 = [1, 2, 3, 4, 6, 7]

# ---------------------------------------------------------------------------------------------------------------
ck.section("0. objects (own construction) and cross-check with the audit library")
C.basic_selfchecks(ck)
C.crosscheck_library(ck, {15: set(T15)})

# ---------------------------------------------------------------------------------------------------------------
ck.section("1. contact set, generators, minor case (document p.1)")
gens = C.generators(T15)
ck.check("the document's contact set {1,2,3,4,6,7} is literally the paper's T_15 (Appendix E); identity group element",
         sorted(T15) == [1, 2, 3, 4, 6, 7])
ck.check("generators vanishing on T_15 are exactly {0,...,8,53,54}", gens == [0, 1, 2, 3, 4, 5, 6, 7, 8, 53, 54], gens)
ck.check("f_54 is the only non-minor generator", [j for j in gens if not C.IS_MINOR[j]] == [54])
mfun = C.functional({16: 1, 0: -1})          # g_16 - g_0 = -d
vals = {j: mfun(F[j]) for j in gens}
ck.check("g_16 - g_0 equals 1 on f_1, f_4, f_7, -1 on f_54 and 0 on the other generators",
         vals == {0: 0, 1: 1, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 7: 1, 8: 0, 53: 0, 54: -1}, vals)


def eqp(lhs, rhs):
    return C.lin_comb(lhs) == C.lin_comb(rhs)


ck.check("f_54 + f_1 = f_0 + f_6 + f_53", eqp({54: 1, 1: 1}, {0: 1, 6: 1, 53: 1}))
ck.check("f_54 + f_4 = f_2 + f_3 + f_6 + f_53", eqp({54: 1, 4: 1}, {2: 1, 3: 1, 6: 1, 53: 1}))
ck.check("f_54 + f_7 = f_0 + f_5 + f_8 + f_53", eqp({54: 1, 7: 1}, {0: 1, 5: 1, 8: 1, 53: 1}))
ck.check("the right-hand sides of the three identities use minor vectors only",
         all(C.IS_MINOR[j] for j in (0, 6, 53, 2, 3, 5, 8)))
# canonical decomposition (2):  g = g22 f0 + g13 f2 + g12 f3 + g10 f5 + g17 f6 + g8 f8 + g16 f53 + d f54, d = g0 - g16
coef2 = {0: {22: 1}, 2: {13: 1}, 3: {12: 1}, 5: {10: 1}, 6: {17: 1}, 8: {8: 1}, 53: {16: 1}, 54: {0: 1, 16: -1}}
ck.check("identity (2) holds on the linear span of the 11 generators (checked on each generator)",
         all(C.lin_comb({k: C.functional(e)(F[j]) for k, e in coef2.items()}) == F[j] for j in gens))
ck.check("(3): g_11 = g_10 + g_17 + d and g_19 = g_13 + g_22 + d on the span of the generators",
         all(C.functional({11: 1, 10: -1, 17: -1, 0: -1, 16: 1})(F[j]) == 0 and
             C.functional({19: 1, 13: -1, 22: -1, 0: -1, 16: 1})(F[j]) == 0 for j in gens))
ck.check("f_54 is positive at pi_0 and zero at pi_16 (so d = g_0 - g_16 has coefficient +1 on f_54)",
         F[54][0] == 1 and F[54][16] == 0)

# ---------------------------------------------------------------------------------------------------------------
ck.section("2. the SOC inequality invoked: (f_2, f_5; f_54) with a b >= d^2 (document p.1 eq.(1), p.3)")
k = C.check_atom_use(ck, "SOC", C.SOC_REP, (2, 5, 54), 2)
ck.check("the SOC orbit has 144 elements (pair unordered)", C.orbit_size(C.SOC_REP, 2) == 144)
ck.check("in (2) the coefficients of f_2, f_5, f_54 are a = g_13, b = g_10, d: the pair (a, b) of the document",
         coef2[2] == {13: 1} and coef2[5] == {10: 1})
Os = C.rational_so4_samples(6, 1501)
lam = [sp.Rational(1, 3), sp.Rational(1), sp.Rational(5, 2)]
soc_exact = True
for O in Os:
    x2, x5, p54 = C.minor_value(2, O), C.minor_value(5, O), C.leibniz(F[54], O)
    # pointwise form squared-free test: Phi_54 + 2|x_2 x_5| >= 0
    soc_exact &= (p54 + 2 * abs(x2 * x5) >= 0)
    for tau in lam:                 # coefficient form with a b = d^2: a = tau, b = 1/tau, d = 1
        soc_exact &= (tau * x2 ** 2 + x5 ** 2 / tau + p54 >= 0)
ck.check("sanity (exact rational SO(4) points, not a proof): Phi_54 + 2|x_2 x_5| >= 0 and its coefficient form at 6 "
         "matrices", soc_exact)

# ---------------------------------------------------------------------------------------------------------------
ck.section("3. the chart (document p.2)")
z, w, v, x, t = sp.symbols("z w v x t")
A = [z, sp.Integer(1), sp.Integer(0), x]
Cc = [w, v, C.INF, t * v]
P = C.products(A, Cc)
ck.check("with C_2 at infinity, P_pi = - prod_{i: pi(i) != 2} (A_i - C_pi(i))",
         all(C.same(P[j], -sp.Mul(*[(A[i] - Cc[s[i]]) for i in range(4) if s[i] != 2])) for j, s in enumerate(C.PERMS)))
R = [sp.cancel((P[i] - P[1]) / (v * (z - w))) for i in range(24)]
ck.check("(4): R_2 = x - t, R_3 = t(x - 1), R_4 = 1 - t (z, w, v, x, t free)",
         C.same(R[2], x - t) and C.same(R[3], t * (x - 1)) and C.same(R[4], 1 - t))
ck.check("(P_3 - P_1)/(P_2 - P_4) = t", C.same(R[3] / (R[2] - R[4]), t))
J = x - t - R[6]
ck.check("R_7 + t = t v/(x - t v) * J with J = x - t - R_6", C.same(R[7] + t, t * v / (x - t * v) * J))
ck.check("J = (x - t v)[z(v + 1 - w) - v] / (v (z - w))", C.same(J, (x - t * v) * (z * (v + 1 - w) - v) / (v * (z - w))))
vb = sp.Symbol("vb")
ck.check("Im(t v/(x - t v)) = t x Im v / |x - t v|^2 (nonzero for real t, x != 0 and nonreal v)",
         C.same(C.im_ratio(t * v / (x - t * v), v, vb), t * x / ((x - t * v) * (x - t * vb))))

# chart coordinates p = v/(w-1), q = 1/v  (both complex);  v = 1/q, w = 1 + 1/(p q), z = p/(p-1)
p, q, pb, qb = sp.symbols("p q pb qb")
CONJ = {p: pb, pb: p, q: qb, qb: q}


def cj(e):
    return e.xreplace(CONJ)


def im2i(e):
    """2 i Im(e) for a rational function of (p, q) with real coefficients."""
    return e - cj(e)


vv = 1 / q
ww = 1 + 1 / (p * q)
zz = sp.cancel(vv / (vv + 1 - ww))
ck.check("(5) z = v/(v + 1 - w) equals p/(p - 1); z - w = (p q - p + 1)/(q p (p - 1))",
         C.same(zz, p / (p - 1)) and C.same(zz - ww, (p * q - p + 1) / (q * p * (p - 1))))
Rc = [sp.cancel(r.subs({z: zz, w: ww, v: vv})) for r in R]
ck.check("on the chart all six contact responses are real: (R_1,R_2,R_3,R_4,R_6,R_7) = (0, x-t, t(x-1), 1-t, x-t, -t)",
         [Rc[j] for j in (1, 2, 3, 4, 6, 7)] == [0, x - t, sp.expand(t * (x - 1)), 1 - t, x - t, -t]
         or all(C.same(a_, b_) for a_, b_ in zip([Rc[j] for j in (1, 2, 3, 4, 6, 7)],
                                                 [0, x - t, t * (x - 1), 1 - t, x - t, -t])))
ck.check("(6): Im R_11 = (1 - x) Im R_10 and Im R_19 = ((t - 1)/t) Im R_13 on the chart (hence g_11 = (1-x) g_10, "
         "g_19 = ((t-1)/t) g_13)",
         C.is_zero(im2i(Rc[11]) - (1 - x) * im2i(Rc[10])) and C.is_zero(im2i(Rc[19]) - (t - 1) / t * im2i(Rc[13])))
K = p * q - p * t - t + 1
Lq = ((x - 1) * p - x) / (p * q - p + 1)
ck.check("K := R_10 = p q - p t - t + 1 and R_13 = -t (L + 1), L = ((x-1)p - x)/(p q - p + 1)",
         C.same(Rc[10], K) and C.same(Rc[13], -t * (Lq + 1)))

# ---------------------------------------------------------------------------------------------------------------
ck.section("4. the reciprocal obstruction (document p.3)")
X, Tt = sp.symbols("X T")
neg = {x: -X, t: -Tt}
Hh, Cx, kk = 1 + Tt, 1 + X, Tt * X - 1
U = Hh - K.subs(neg)
V = X - Lq.subs(neg)
ck.check("U = H - K = -p (T + q) and V = X - L = p (X q + 1)/(p q - p + 1)  (x = -X, t = -T)",
         C.same(U, -p * (Tt + q)) and C.same(V, p * (X * q + 1) / (p * q - p + 1)))
D1 = Cx * U + kk
D2 = kk - Hh * V
Z = (Tt - X) * U * V - Hh * X * U - Tt * Cx * V
ck.check("(8): R_12 = (T - X) - Z/D_1 and R_8 = (T - X) + Z/D_2 (identities in p, q, X, T)",
         C.same(Rc[12].subs(neg), (Tt - X) - Z / D1) and C.same(Rc[8].subs(neg), (Tt - X) + Z / D2))
Us, Vs, Ub, Vb = sp.symbols("U V Ub Vb")
D1s, D2s = Cx * Us + kk, kk - Hh * Vs
Zs = (Tt - X) * Us * Vs - Hh * X * Us - Tt * Cx * Vs
ck.check("(9): Z = (X-T)/(C H) D_1 D_2 - (T C/H) D_1 + (X H/C) D_2 (identity in U, V, X, T)",
         C.same(Zs, (X - Tt) / (Cx * Hh) * D1s * D2s - Tt * Cx / Hh * D1s + X * Hh / Cx * D2s))
CONJ2 = {Us: Ub, Ub: Us, Vs: Vb, Vb: Vs}
cj2 = lambda e: e.xreplace(CONJ2)
im2 = lambda e: e - cj2(e)
ck.check("Im D_1 = C Im U and Im D_2 = -H Im V", C.is_zero(im2(D1s) - Cx * im2(Us)) and C.is_zero(im2(D2s) + Hh * im2(Vs)))
# positive-weight identity, multiplied by 2i:  m = Im(Z/D1), n = -Im(Z/D2), A = Im V, B = Im U
n1, n2 = D1s * cj2(D1s), D2s * cj2(D2s)
lhs = Tt * Cx ** 2 * n1 * im2(Zs / D1s) - X * Hh ** 2 * n2 * im2(Zs / D2s)
rhs = (Tt - X) * (Tt * Cx * im2(Vs) * n1 + X * Hh * im2(Us) * n2)
ck.check("T C^2 |D_1|^2 m + X H^2 |D_2|^2 n = (T - X)(T C A |D_1|^2 + X H B |D_2|^2), with m = Im(Z/D_1), "
         "n = -Im(Z/D_2), A = Im V, B = Im U (polynomial identity in U, Ub, V, Vb, X, T)", C.is_zero(lhs - rhs))
ck.check("gap dictionary: g_10 = -mu Im K = mu Im U, g_13 = -mu Im R_13 = -mu T Im L = mu T Im V; g_12 = mu Im(Z/D_1), "
         "g_8 = -mu Im(Z/D_2)",
         C.is_zero(im2i(K.subs(neg)) + im2i(U)) and C.is_zero(im2i(Rc[13].subs(neg)) - Tt * im2i(Lq.subs(neg)))
         and C.is_zero(im2i(Lq.subs(neg)) + im2i(V))
         and C.is_zero(im2i(Rc[12].subs(neg)) + im2i(Z / D1)) and C.is_zero(im2i(Rc[8].subs(neg)) - im2i(Z / D2)))
a_, b_, d_, r_, s_ = sp.symbols("a b d r s")
# the arithmetic of 'Spectral signs': (1-x) b = b + r + d  and  ((t-1)/t) a = a + s + d
ck.check("(3) and (6) give -x b = r + d and -a/t = s + d",
         C.same(sp.solve(sp.Eq((1 - x) * b_, b_ + r_ + d_), r_)[0] + d_, -x * b_)
         and C.same(sp.solve(sp.Eq((t - 1) / t * a_, a_ + s_ + d_), s_)[0] + d_, -a_ / t))
ck.info("PROOF-READ (see T15.md): d > 0 and r, s >= 0 force a, b > 0, x = -X < 0, t = -T < 0, d <= X b, d <= a/T; "
        "m, n >= 0 and the identity force T >= X; hence a b >= (T/X) d^2 >= d^2.")

# ---------------------------------------------------------------------------------------------------------------
ck.section("5. NUMERIC falsification attempts (floating point; evidence only)")
import numpy as np

npr = np.random.default_rng(20261125)
TOL = 1e-9
CH = 500000


def gapsN(xv, tv, pv, qv, mu):
    v_ = 1 / qv
    w_ = 1 + 1 / (pv * qv)
    z_ = pv / (pv - 1)
    kap = 1j * mu / (v_ * (z_ - w_))
    gN = C.numeric_gaps([z_, 1.0 + 0 * v_, 0.0 * v_, xv + 0 * v_], [w_, v_, None, tv * v_], kap, ref=1)
    sc = np.abs(gN).max(axis=1)
    return gN / sc[:, None]


def mag(n):
    return np.exp(npr.normal(0, 1.2, n))


def run(nchunks, sx, st):
    acc = dict(contact=0.0, nsupp=0, nhard=0, nout=0, soc=np.inf, TX=np.inf)
    for _ in range(nchunks):
        sg = lambda s_: npr.choice([-1.0, 1.0], CH) if s_ is None else s_ * np.ones(CH)
        xv, tv = sg(sx) * mag(CH), sg(st) * mag(CH)
        pv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
        qv = npr.normal(0, 1.5, CH) + 1j * npr.normal(0, 1.5, CH)
        mu = npr.choice([-1.0, 1.0], CH)
        gN = gapsN(xv, tv, pv, qv, mu)
        ok = np.isfinite(gN).all(axis=1)
        acc["contact"] = max(acc["contact"], np.abs(gN[ok][:, T15]).max())
        supp = ok & (gN.min(axis=1) >= -TOL)
        dN = gN[:, 0] - gN[:, 16]
        hard = supp & (dN > 1e-7)
        acc["nsupp"] += int(supp.sum())
        acc["nhard"] += int(hard.sum())
        acc["nout"] += int((hard & ~((xv < 0) & (tv < 0))).sum())
        if hard.any():
            acc["soc"] = min(acc["soc"], (gN[:, 13] * gN[:, 10] - dN ** 2)[hard].min())
            acc["TX"] = min(acc["TX"], (xv - tv)[hard & (xv < 0) & (tv < 0)].min())       # T - X = x - t
    return acc


acc = run(6, None, None)
ck.check("NUMERIC global: 3,000,000 chart points (seed 20261125; x, t = random sign * lognormal(0, 1.2); p, q complex "
         "normal(0, 1.5); mu = +-1; gaps normalised to max 1): max |g| on T_15 = %.1e; supporting (all gaps >= -1e-9): %d; "
         "of these with d > 1e-7: %d, outside {x < 0, t < 0}: %d; min(g13 g10 - d^2) = %.3e, min(T - X) = %.3e"
         % (acc["contact"], acc["nsupp"], acc["nhard"], acc["nout"], acc["soc"], acc["TX"]),
         acc["contact"] < 1e-6 and acc["nout"] == 0 and acc["nhard"] > 50 and acc["soc"] > -1e-9 and acc["TX"] > -1e-9)
acc = run(6, -1, -1)
ck.check("NUMERIC region x < 0, t < 0: 3,000,000 points: supporting with d > 1e-7: %d; min(g13 g10 - d^2) = %.3e, "
         "min(T - X) = %.3e (tolerance -1e-9)" % (acc["nhard"], acc["soc"], acc["TX"]),
         acc["nhard"] > 50 and acc["soc"] > -1e-9 and acc["TX"] > -1e-9)
# the intermediate claim 'm, n >= 0, Im U, Im V > 0 imply T >= X' on its own region (U, V free in the upper half-plane)
n_ = 2000000
Xv, Tv = np.exp(npr.normal(0, 1.0, n_)), np.exp(npr.normal(0, 1.0, n_))
Uv = npr.normal(0, 2, n_) + 1j * (np.abs(npr.normal(0, 2, n_)) + 1e-6)
Vv = npr.normal(0, 2, n_) + 1j * (np.abs(npr.normal(0, 2, n_)) + 1e-6)
Hn, Cn, kn = 1 + Tv, 1 + Xv, Tv * Xv - 1
D1n, D2n = Cn * Uv + kn, kn - Hn * Vv
Zn = (Tv - Xv) * Uv * Vv - Hn * Xv * Uv - Tv * Cn * Vv
m_, nn_ = (Zn / D1n).imag, -(Zn / D2n).imag
sel = Xv > Tv
both = int((sel & (m_ >= 0) & (nn_ >= 0)).sum())
ck.check("NUMERIC: %d samples with X > T and U, V in the upper half-plane (seed as above): m >= 0 and n >= 0 together "
         "occurs %d times (claim: never); largest min(m, n) = %.3e" % (sel.sum(), both, np.minimum(m_, nn_)[sel].max()),
         both == 0 and sel.sum() > 1000)
best = C.numeric_atom_attack("SOC", (2, 5, 54), 20000, 6, 1502)
ck.check("NUMERIC: min over SO(4) of Phi_54 + 2|x_2 x_5| found = %.3e (20000 random points + 6 Nelder-Mead runs, "
         "seed 1502); tolerance -1e-9" % best, best > -1e-9)
# real v (the stratum delegated to the three-plus-three theorem): the gap vector is NOT forced to vanish here
xv, tv = -mag(4), -mag(4)
vr = npr.normal(0, 1.5, 4) + 0j
zc = npr.normal(0, 1.5, 4) + 1j * npr.normal(0, 1.5, 4)
wc = npr.normal(0, 1.5, 4) + 1j * npr.normal(0, 1.5, 4)
gR = C.numeric_gaps([zc, 1.0 + 0 * vr, 0.0 * vr, xv + 0 * vr], [wc, vr, None, tv * vr], 1j / (vr * (zc - wc)), ref=1)
ck.check("NUMERIC remark: for real v and generic complex z, w the responses of pi_2, pi_3, pi_4 are real but other "
         "gaps are not zero (so the real-v stratum of T_15 is not trivial): max |g| on {1,2,3,4} = %.1e, max |g| = %.2f"
         % (np.abs(gR[:, [1, 2, 3, 4]]).max(), np.abs(gR).max()), np.abs(gR[:, [1, 2, 3, 4]]).max() < 1e-9 and np.abs(gR).max() > 1e-3)

ck.finish()
