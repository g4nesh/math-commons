# Auditor replay for independent_moc4_T39_high_j_uniform_20261022 (own code; nothing from the supplement is executed or imported).
# Run from this directory: python -I -B independent_moc4_T39_high_j_uniform_20261022__check.py   (needs ../lib/moc4lib.py, sympy, python-flint, numpy, scipy, mpmath)
# Env: NLINES / NSEEDS / NPTS control sample sizes.
# ---- shared auditor helper (own code; rebuilt from paper definitions via moc4lib) ----
import os, sys, random
from fractions import Fraction as Fr
_here = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(_here, "..", "lib"), os.environ.get("MOC4LIB", "")):
    if _p and os.path.isdir(_p): sys.path.insert(0, _p)
from moc4lib import *
import sympy as sp
import flint
FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(info)) if info != "" else ""))
    if not ok: FAIL.append(name)
def info(name, msg): print("INFO " + name + " : " + str(msg))
def fvec(j): return [Fr(x) for x in F_VECS[j]]
def comb(pairs):
    v = [Fr(0)] * 14
    for cf, j in pairs:
        for r in range(14): v[r] += Fr(cf) * F_VECS[j][r]
    return v
def vals_of(v): return [Fr(x) for x in values(v)]
# third derivative tensor of the paper cubic C (A.2)
_ms, _Cexpr = cubic_sympy()
TENS = {}
for _i in range(14):
    for _j in range(_i, 14):
        _d2 = sp.diff(_Cexpr, _ms[_i], _ms[_j])
        for _k in range(14):
            _v = sp.diff(_d2, _ms[_k])
            if _v != 0: TENS[(_i, _j, _k)] = int(_v)
def hessian(w):
    H = [[Fr(0)] * 14 for _ in range(14)]
    for (i, j, k), t in TENS.items():
        H[i][j] += t * w[k]
        if i != j: H[j][i] += t * w[k]
    return H
def fdet(M):
    return Fr(str(flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator) for x in row] for row in M]).det()))
def CQ(g):
    """C(g + h kappa) = C0 + h*Q exactly (checks h^2, h^3 terms vanish)."""
    w0 = list(g); w1 = list(g); w1[13] += 1; w2 = list(g); w2[13] += 2; w3 = list(g); w3[13] += 3
    v = [cubic_C(x) for x in (w0, w1, w2, w3)]
    C0 = v[0]; Q = v[1] - v[0]
    assert v[2] - v[0] == 2 * Q and v[3] - v[0] == 3 * Q, "cubic not affine in h"
    return C0, Q
def G2_at(g):
    """-(Q^2/4) det Hess C(g - (C/Q) kappa); claimed to equal G^2."""
    C0, Q = CQ(g)
    w = list(g); w[13] -= C0 / Q
    return -Q * Q * fdet(hessian(w)) / 4, C0, Q
def is_rat_square(x):
    if x < 0: return False
    n, d = x.numerator, x.denominator
    from math import isqrt
    return isqrt(n) ** 2 == n and isqrt(d) ** 2 == d
def interp_poly(ts, ys):
    t = sp.Symbol('t')
    return sp.Poly(sp.interpolate(list(zip([sp.Rational(x.numerator, x.denominator) for x in ts],
                                          [sp.Rational(y.numerator, y.denominator) for y in ys])), t), t)
def G2_on_line(gfun, deg, lo, hi, extra=4):
    """Exact G^2(t) along t -> gfun(t) (coefficients affine in t), interpolated with
    deg+1+extra points; the extra points certify the degree bound."""
    n = deg + 1 + extra
    ts = [Fr(lo) + (Fr(hi) - Fr(lo)) * Fr(k, n - 1) for k in range(n)]
    ys = [G2_at(gfun(t))[0] for t in ts]
    t = sp.Symbol('t')
    Pfit = interp_poly(ts[:deg + 1], ys[:deg + 1])
    ok = all(Pfit.eval(sp.Rational(x.numerator, x.denominator)) == sp.Rational(y.numerator, y.denominator) for x, y in zip(ts[deg + 1:], ys[deg + 1:]))
    return Pfit, ok
def roots_in(Pfit, lo, hi):
    sq = sp.Poly(sp.sqf_part(Pfit.as_expr()), Pfit.gen)
    if sq.degree() <= 0: return 0
    return sq.count_roots(sp.Rational(lo.numerator, lo.denominator), sp.Rational(hi.numerator, hi.denominator))
# ===== independent_moc4_T39_high_j_uniform_20261022 : high-j chamber of paper T_40 =====
import mpmath as mp
mp.mp.dps = 50
q, p, j, c, e, d, b, a, h = sp.symbols('q p j c e d b a h')
SY = [q, p, j, c, e, d, b, a]; S35 = [5, 1, 12, 9, 2, 35, 33, 29]
v = [0] * 14
for s_, jj in zip(SY, S35):
    for k in range(14): v[k] += s_ * sp.Rational(F_VECS[jj][k].numerator, F_VECS[jj][k].denominator)
v[13] += h
Cp = sp.Poly(sp.expand(cubic_C(v)), h); Qs, Cs = Cp.all_coeffs()
A_ = 1 + a + b; B_ = 1 + q + e; N = (q + 1) * (e + 1) * A_ + B_ * a * b; Nj = N - j * e * a; S_ = a + e; J_ = j - 1
check("C = -(ABc + Nj) at d=1", sp.expand(Cs.subs(d, 1) + (A_ * B_ * c + Nj)) == 0)
check("-Q formula at d=1", sp.expand(-Qs.subs(d, 1) - ((q + b + 1 + S_ - p) * c + 1 - J_ * S_ + (q + b) * (1 + S_) + p * j)) == 0)
R1 = sp.Poly(sp.expand(N - (1 + S_ + e * a)), q, e, a, b)
check("N - (1 + S + ea) has nonnegative coefficients", all(x >= 0 for x in R1.coeffs()))
# chart (5) inverse: P = pj/(1-pj), V = JS/(1-JS), L:M = a:e reproduces (p, S, a, e)
P_, V_, L_, M_, JJ = sp.symbols('P V L M J', positive=True)
jj_ = 1 + JJ
pp = P_ / (jj_ * (1 + P_)); SS = V_ / (JJ * (1 + V_))
check("chart (5): p*j = P/(1+P) and J*S = V/(1+V), i.e. P = pj/(1-pj), V = JS/(1-JS)", sp.simplify(pp * jj_ - P_ / (1 + P_)) == 0 and sp.simplify(JJ * SS - V_ / (1 + V_)) == 0)
# pointwise tests at random rational points of chamber (1)
def src(x, cc):
    q_, p_, j_, e_, b_, a_ = x
    return comb([(q_, 5), (p_, 1), (j_, 12), (cc, 9), (e_, 2), (1, 35), (b_, 33), (a_, 29)])
rng = random.Random(40)
R = lambda lo, hi: Fr(rng.randint(int(lo * 1000), int(hi * 1000)), 1000)
def rnd_chamber():
    while True:
        j_ = 1 + R(0, 3); p_ = R(0, 0.999) / j_; a_ = R(0, 3); e_ = R(0, 3); q_ = R(0, 3) * rng.choice([0, 1]); b_ = R(0, 3) * rng.choice([0, 1])
        if (j_ - 1) * (a_ + e_) < 1: return [q_, p_, j_, e_, b_, a_]
H_np = None
def hess_mp(w):
    H = mp.matrix(14, 14)
    for (i, jx, k), t in TENS.items():
        H[i, jx] += t * w[k]
        if i != jx: H[jx, i] += t * w[k]
    return H
npt = 0; okall = True; stats = dict(n=0, ineq2=0, gamma_neg=0, alpha_pos=0)
# chamber points with alpha > 0 found by an auditor search over 23668 random chamber points (558 had alpha > 0)
APOS = [['6780', '93800/100839', '100839/100000', '18705/839', '7100', '10295/839'], ['7580', '4940000/10000681', '10000681/10000000', '309680/681', '533/100', '5220320/681'], ['857/10', '2425000/2500061', '2500061/2500000', '1250', '777/100', '8750'], ['520', '1862/2123', '2123/2000', '8103/2050', '765', '48691/6150'], ['27/2', '867000/1000111', '1000111/1000000', '32818/111', '711/10', '9394/37'], ['591/100', '11740/20199', '20199/20000', '483723/9950', '119/20', '273277/9950'], ['7160', '50200/100213', '100213/100000', '90117/710', '143/2', '62883/710']]
NP = int(os.environ.get("NPTS", "6"))
for k in range(NP + len(APOS)):
    if k < NP:
        x = rnd_chamber()
        if k == 0: x[2] = Fr(1)   # include the wall j = 1
    else:
        x = [Fr(z) for z in APOS[k - NP]]
        assert x[1] * x[2] < 1 and x[2] >= 1 and (x[2] - 1) * (x[3] + x[5]) < 1
    Pf, okd = G2_on_line(lambda t: src(x, t), 4, Fr(0), Fr(5))
    t = Pf.gen
    Gq = sp.Poly(sp.sqf_part(Pf.as_expr()), t)   # G up to a constant factor
    isq = sp.expand(Gq.as_expr() ** 2 * Pf.LC() / Gq.LC() ** 2 - Pf.as_expr()) == 0
    al, be, ga = [Gq.coeff_monomial(t ** i) for i in (2, 1, 0)]
    if ga > 0: al, be, ga = -al, -be, -ga    # sign of G fixed by the claim gamma < 0 (G is only defined up to sign by G^2)
    q_, p_, j_, e_, b_, a_ = x
    AB = (1 + a_ + b_) * (1 + q_ + e_); Njv = (q_ + 1) * (e_ + 1) * (1 + a_ + b_) + (1 + q_ + e_) * a_ * b_ - j_ * e_ * a_
    ineq2 = AB * be - Njv * al < 0
    stats["n"] += 1; stats["ineq2"] += bool(ineq2); stats["gamma_neg"] += bool(ga < 0); stats["alpha_pos"] += bool(al > 0)
    Cn, Qn = CQ(src(x, Fr(0)))
    rts = sorted([r_ for r_ in sp.Poly(al * t**2 + be * t + ga, t).nroots(n=40) if r_.is_real], key=lambda z: float(z))
    pos = [r_ for r_ in rts if r_ > 0]
    good = okd and isq and Gq.degree() == 2 and ga < 0 and ineq2 and (len(pos) == 1 if al > 0 else len(pos) == 0)
    # at c+: Hessian rank, P*, bordered minors, binary form
    cp = mp.mpf(str(pos[0])) if len(pos) == 1 else None
    if cp is not None:
        g0 = [mp.mpf(z.numerator) / z.denominator for z in src(x, Fr(0))]
        f9 = [mp.mpf(z.numerator) / z.denominator for z in F_VECS[9]]
        C0, Q0 = CQ(src(x, Fr(0))); C1, Q1 = CQ(src(x, Fr(1)))
        Ca, Cb = mp.mpf(C0.numerator) / C0.denominator, mp.mpf((C1 - C0).numerator) / (C1 - C0).denominator
        Qa, Qb = mp.mpf(Q0.numerator) / Q0.denominator, mp.mpf((Q1 - Q0).numerator) / (Q1 - Q0).denominator
        Cc, Qc = Ca + Cb * cp, Qa + Qb * cp
        hh = -Cc / Qc; hd = -(Cb * Qc - Cc * Qb) / Qc ** 2
        w = [g0[i] + cp * f9[i] for i in range(14)]; w[13] += hh
        vdir = list(f9); vdir[13] += hd
        H = hess_mp(w)
        U, Sv, Vt = mp.svd_r(H)
        sv = sorted([abs(Sv[i]) for i in range(14)])
        Pstar = mp.det(H[2:14, 2:14])
        bord0 = mp.det(mp.matrix([[H[i, jx] for jx in [0] + list(range(2, 14))] for i in [0] + list(range(2, 14))]))
        bord1 = mp.det(mp.matrix([[H[i, jx] for jx in [1] + list(range(2, 14))] for i in [1] + list(range(2, 14))]))
        Hv = hess_mp(vdir)
        idx = sorted(range(14), key=lambda i: abs(Sv[i]))[:2]
        K = [[Vt[i, r] for r in range(14)] for i in idx]
        Bf = [[sum(K[u][r] * Hv[r, s] * K[w_][s] for r in range(14) for s in range(14)) for w_ in range(2)] for u in range(2)]
        detB = Bf[0][0] * Bf[1][1] - Bf[0][1] * Bf[1][0]
        scale = abs(Pstar) * max(abs(H[i, jx]) for i in range(14) for jx in range(14))  # relative tolerance for 13-minors
        good = good and sv[1] < sv[2] * mp.mpf(10) ** -30 and Pstar < 0 and abs(bord0) < scale * mp.mpf(10) ** -25 and abs(bord1) < scale * mp.mpf(10) ** -25 and detB > 0 and Qc < 0
        info(f"point {k}", f"x(q,p,j,e,b,a)={[float(z) for z in x]} G deg {Gq.degree()}, gamma<0 {ga < 0}, c+={float(cp):.5g}, sing.vals {mp.nstr(sv[0],3)},{mp.nstr(sv[1],3)},{mp.nstr(sv[2],3)}, P*={mp.nstr(Pstar,5)}, bordered {mp.nstr(bord0,3)},{mp.nstr(bord1,3)}, det binary form {mp.nstr(detB,5)}, Q={mp.nstr(Qc,5)}")
    else:
        info(f"point {k}", f"x={[float(z) for z in x]} alpha<=0, no nonnegative root (as claimed: G<0 on c>=0)")
    npt += 1; okall &= bool(good)
info("sign statistics", stats)
check(f"NUMERIC: at {npt} rational chamber points (one on j=1): G is quadratic in c with gamma<0, unique simple positive root c+; at c+ the Hessian has rank 12, P*<0, both bordered 13-minors vanish, K^T H(w_c) K is definite (NUMERIC, 50 digits; G from exact interpolation)", okall)
ex = [f for f in FAIL if not f.startswith("NUMERIC")]
print(f"SUMMARY: {len(FAIL)} failed checks ({len(ex)} exact)")
sys.exit(1 if ex else 0)
