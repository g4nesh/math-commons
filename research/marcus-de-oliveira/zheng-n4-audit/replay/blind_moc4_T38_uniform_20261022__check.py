# Auditor replay for blind_moc4_T38_uniform_20261022 (own code; nothing from the supplement is executed or imported).
# Run from this directory: python -I -B blind_moc4_T38_uniform_20261022__check.py   (needs ../lib/moc4lib.py, sympy, python-flint, numpy, scipy, mpmath)
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
# ===== blind_moc4_T38_uniform_20261022 : paper case T_39 =====
import numpy as np
T = CONTACT_SETS[39]
check("contact set {2,3,5,7,9,12} == paper T_39", T == frozenset({2, 3, 5, 7, 9, 12}))
check("generators(T_39) == {0,1,2,5,7,12,39,40,67,70}", generators(T) == [0, 1, 2, 5, 7, 12, 39, 40, 67, 70])
check("non-minor generators are exactly 40, 70", [j for j in generators(T) if not F_IS_MINOR[j]] == [40, 70])
def vv(lst): return [sum(F_VALUES[j][s] * k for k, j in lst) for s in range(24)]
check("identity f40+f70 = f0+f5+f12+f39+f67", vv([(1, 40), (1, 70)]) == vv([(1, 0), (1, 5), (1, 12), (1, 39), (1, 67)]))
check("identity f1+f40 = f0+f7+f39", vv([(1, 1), (1, 40)]) == vv([(1, 0), (1, 7), (1, 39)]))
check("identity f7+f70 = f1+f5+f12+f67", vv([(1, 7), (1, 70)]) == vv([(1, 1), (1, 5), (1, 12), (1, 67)]))
B40 = [39, 12, 5, 7, 2, 40, 67, 0]; B70 = [0, 1, 2, 12, 39, 67, 70, 5]
check("source bases (1) and (2) have rank 8", rank([fvec(j) for j in B40]) == 8 and rank([fvec(j) for j in B70]) == 8)
# T24 core action: some g in G with g.f59,f53,f55,f67,f54 = f5,f2,f67,f0,f40 (as value vectors)
def gmap_f(g, j):
    vals = g_act_values(g, list(F_VALUES[j])); return [k for k in range(82) if list(F_VALUES[k]) == [int(x) for x in vals]]
found = None
for g in G_ELEMENTS:
    if all(gmap_f(g, s) == [t] for s, t in zip([59, 53, 55, 67, 54], [5, 2, 67, 0, 40])): found = g; break
check("some g in G maps (f59,f53,f55,f67,f54) to (f5,f2,f67,f0,f40) in order", found is not None, found)
info("generators of paper T_25 (supplement T24)", generators(CONTACT_SETS[25]))
# radial formulas
q, p, j, c, e, d, b, a, r, s, h = sp.symbols('q p j c e d b a r s h')
def CQsym(pairs):
    v = [0] * 14
    for s_, jj in pairs:
        for k in range(14): v[k] += s_ * sp.Rational(F_VECS[jj][k].numerator, F_VECS[jj][k].denominator)
    v[13] += h
    Cp = sp.Poly(sp.expand(cubic_C(v)), h); return Cp.degree(), Cp.all_coeffs()
dg, (Q40, C40) = CQsym([(q, 39), (p, 12), (j, 5), (c, 7), (e, 2), (d, 40), (b, 67), (a, 0)])
check("g40: Q = d^2+de+(d+e)(p+q)-ej", dg == 1 and sp.expand(Q40 - (d**2 + d*e + (d + e)*(p + q) - e*j)) == 0)
CP = sp.Poly(C40, q, p, j, c, e, d, b, a)
check("g40: C has 21 positive coefficients incl. d^3", len(CP.terms()) == 21 and all(x > 0 for x in CP.coeffs()) and CP.coeff_monomial(d**3) == 1)
dg, (Q70, C70) = CQsym([(p, 0), (r, 1), (a, 2), (b, 12), (c, 39), (q, 67), (d, 70), (s, 5)])
M = b + c; U = d + r; T0 = b*U*(d + q) + c*d*(d + p + r) + d*(d + p + r)*(d + q)
check("g70: C = bc(a-U) - T0", sp.expand(C70 - (b*c*(a - U) - T0)) == 0)
check("g70: Q = (a-U)M - dU - as", sp.expand(Q70 - ((a - U)*M - d*U - a*s)) == 0)
ex = sp.expand((M*T0 - b*c).subs(d, 1))
check("g70: -MC = bc(M+1-aM) + rMbc + (MT0 - bc) at d=1", sp.expand((-M*C70 - (b*c*(M + 1 - a*M) + r*M*b*c + (M*T0 - b*c))).subs(d, 1)) == 0)
check("g70: MT0 - bc (d=1) has nonnegative coefficients and contains b and c", all(x >= 0 for x in sp.Poly(ex, b, c, p, q, r).coeffs()) and sp.Poly(ex, b, c, p, q, r).coeff_monomial(b) > 0 and sp.Poly(ex, b, c, p, q, r).coeff_monomial(c) > 0)
check("g70: M=0 gives C = -(1+p+r)(1+q) at d=1", sp.expand(C70.subs({b: 0, c: 0, d: 1}) + (1 + p + r)*(1 + q)) == 0)
def s40(x):  # x = q,p,j,c,e,b,a ; d=1
    q_, p_, j_, c_, e_, b_, a_ = x
    return comb([(q_, 39), (p_, 12), (j_, 5), (c_, 7), (e_, 2), (1, 40), (b_, 67), (a_, 0)])
def s70(x):  # x = p,r,a,b,c,q,s ; d=1
    p_, r_, a_, b_, c_, q_, s_ = x
    return comb([(p_, 0), (r_, 1), (a_, 2), (b_, 12), (c_, 39), (q_, 67), (1, 70), (s_, 5)])
rng = random.Random(39)
R = lambda lo, hi: Fr(rng.randint(int(lo * 1000), int(hi * 1000)), 1000)
nsq = 0
for _ in range(4):
    x = [R(0.01, 3) for _ in range(7)]
    G2a = G2_at(s40(x))[0]; G2b = G2_at(s70(x))[0]
    nsq += is_rat_square(G2a) + is_rat_square(G2b)
check("-(Q^2/4) det Hess at h=-C/Q is a rational square (4 points per source)", nsq == 8)
# g40: exact root count of G40 in j on [0, J] in domains 1 and 7
bad = 0; n1 = 0; n7 = 0
for k in range(int(os.environ.get("NLINES", "10"))):
    while True:   # domain 1: a,b<=1, ab<1 ; ac>1 -> p < Pac
        a_, b_, c_, e_, q_ = R(0.01, 1), R(0.01, 1), R(0, 4), R(0.01, 4), R(0, 3)
        if a_ * b_ >= 1: continue
        p_ = R(0, 3)
        if a_ * c_ > 1 and p_ >= (a_ + 1) * (c_ + 1) / (a_ * c_ - 1): continue
        break
    J = 1 + 1 / e_
    Pf, okd = G2_on_line(lambda t: s40([q_, p_, t, c_, e_, b_, a_]), 4, Fr(0), J)
    nr = roots_in(Pf, Fr(0), J); bad += (nr != 0 or not okd); n1 += 1
    while True:   # domain 7: ab>1, ac<1, e < Eab
        a_ = R(0.05, 4); b_ = (1 + R(0.001, 3)) / a_; c_ = R(0, 0.999) / a_; Eab = (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)
        e_ = Eab * R(0.001, 0.999); q_ = R(0, 3); p_ = R(0, 3); break
    J = 1 / e_ - 1 / Eab
    Pf, okd = G2_on_line(lambda t: s40([q_, p_, t, c_, e_, b_, a_]), 4, Fr(0), J)
    nr2 = roots_in(Pf, Fr(0), J); bad += (nr2 != 0 or not okd); n7 += 1
    info(f"g40 lines {k}", f"domain1 roots in [0,J]: {nr}; domain7 roots in [0,J]: {nr2}")
check(f"NUMERIC(random lines, exact Sturm): g40: G40 has no root in j on [0,J] ({n1} domain-1 and {n7} domain-7 random lines; exact Sturm)", bad == 0)
# g70: exact root count in a on [0, min(J, 1/s))
def Jcap(b_, c_, p_, q_):
    best = None
    for k in range(0, 20001):
        b1 = k / 20000; b2 = 1 - b1
        d1 = min(p_ / 2, b1); d2 = min(q_ / 2, b2)
        val = 1 - d1 - d2 + (b1 ** 2 / b_ if b_ > 0 else (0 if b1 == 0 else np.inf)) + (b2 ** 2 / c_ if c_ > 0 else (0 if b2 == 0 else np.inf)) + (d1 ** 2 / p_ if p_ > 0 else 0) + (d2 ** 2 / q_ if q_ > 0 else 0)
        best = val if best is None else min(best, val)
    return best
bad = 0; n = 0
for k in range(int(os.environ.get("NLINES", "10"))):
    p_, r_, b_, c_, q_, s_ = R(0, 3), R(0, 3), R(0.01, 3), R(0.01, 3), R(0, 3), R(0, 2)
    J = Jcap(float(b_), float(c_), float(p_), float(q_))
    top = min(J, 1 / float(s_) if s_ > 0 else np.inf) * (1 - 1e-6)
    top = Fr(top).limit_denominator(10 ** 6)
    Pf, okd = G2_on_line(lambda t: s70([p_, r_, t, b_, c_, q_, s_]), 4, Fr(0), top)
    nr = roots_in(Pf, Fr(0), top); bad += (nr != 0 or not okd); n += 1
    info(f"g70 line {k}", f"J~{J:.4f}, as<1 cap {float(1/s_) if s_>0 else 'inf'}, roots of G70 in a on [0,{float(top):.4f}] = {nr}")
check(f"NUMERIC(random lines, J numeric, exact Sturm): g70: G70 has no root in a on [0, min(J,1/s)) on {n} random lines (J by numeric 1-D minimization; exact Sturm)", bad == 0)
ex = [f for f in FAIL if not f.startswith("NUMERIC")]
print(f"SUMMARY: {len(FAIL)} failed checks ({len(ex)} exact)")
sys.exit(1 if ex else 0)
