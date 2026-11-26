# Auditor replay for independent_moc4_T39_uniform_20261022 (own code; nothing from the supplement is executed or imported).
# Run from this directory: python -I -B independent_moc4_T39_uniform_20261022__check.py   (needs ../lib/moc4lib.py, sympy, python-flint, numpy, scipy, mpmath)
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
# ===== independent_moc4_T39_uniform_20261022 : paper case T_40 =====
import numpy as np
T = CONTACT_SETS[40]
check("contact set {0,1,5,9,14,18} == paper T_40", T == frozenset({0, 1, 5, 9, 14, 18}))
check("generators(T_40) == {1,2,5,9,11,12,28,29,33,35}", generators(T) == [1, 2, 5, 9, 11, 12, 28, 29, 33, 35])
check("non-minor generators are exactly 28, 35", [j for j in generators(T) if not F_IS_MINOR[j]] == [28, 35])
rows = [[P[s][r] for r in range(14)] for s in T]
check("contact nullspace has dimension 8", 14 - rank(rows) == 8)
def vv(lst): return [sum(F_VALUES[j][s] * k for k, j in lst) for s in range(24)]
check("(1a) f28+f12 = f1+f2+f11+f29", vv([(1, 28), (1, 12)]) == vv([(1, 1), (1, 2), (1, 11), (1, 29)]))
check("(1b) f35+f11 = f1+f5+f12+f33", vv([(1, 35), (1, 11)]) == vv([(1, 1), (1, 5), (1, 12), (1, 33)]))
check("(1c) f28+f35 = 2f1+f2+f5+f29+f33", vv([(1, 28), (1, 35)]) == vv([(2, 1), (1, 2), (1, 5), (1, 29), (1, 33)]))
S35 = [5, 1, 12, 9, 2, 35, 33, 29]; S28 = [2, 1, 11, 9, 5, 28, 29, 33]   # q p j c e d b a
check("both canonical sources have rank 8", rank([fvec(j) for j in S35]) == 8 and rank([fvec(j) for j in S28]) == 8)
# Appendix A: (q,p,j,c,e,d,b,a) = (g11,g17,g6,g7,g15,g3-g12,g12,g4) on the whole nullspace
def coeffs_A(vals): return [vals[11], vals[17], vals[6], vals[7], vals[15], vals[3] - vals[12], vals[12], vals[4]]
ok = True
for k in generators(T):
    cf = coeffs_A(list(F_VALUES[k]))
    ok &= vv(list(zip(cf, S35))) == list(F_VALUES[k])
check("App. A: c = g11 f5+g17 f1+g6 f12+g7 f9+g15 f2+(g3-g12) f35+g12 f33+g4 f29 on the contact nullspace (checked on a spanning set)", ok)
q, p, j, c, e, d, b, a, m = sp.symbols('q p j c e d b a m')
SY = [q, p, j, c, e, d, b, a]
def symvals(src, cf):
    return [sp.expand(sum(k * F_VALUES[jj][s] for k, jj in zip(cf, src))) for s in range(24)]
V = symvals(S35, SY)
check("App. A gap identities g3=b+d, g20=q+d, g23=p+d, g8=a+j+d, g19=e+j+d, g22=p+j+2d",
      [V[3], V[20], V[23], V[8], V[19], V[22]] == [b + d, q + d, p + d, a + j + d, e + j + d, p + j + 2 * d])
Vm = symvals(S35, [q, p, j, c, e, -m, b, a])
W1 = [sp.expand(x) for x in symvals([5, 1, 12, 9, 2, 11, 33, 29], [q - m, p - m, j - m, c, e, m, b - m, a])]
check("App. A case j>=m: rewriting into minor sum is an identity", all(sp.expand(x - y) == 0 for x, y in zip(Vm, W1)))
dl = m - j
W2 = symvals(S28, [e - dl, p - m - dl, j, c, q - m, dl, a - dl, b - m])
check("App. A case j<m: rewriting into second source is an identity", all(sp.expand(x - y) == 0 for x, y in zip(Vm, W2)))
# transpose: inverse-permutation action
g_inv = ((0, 1, 2, 3), (0, 1, 2, 3), 1)
check("inverse action preserves T_40", g_act_set(g_inv, T) == frozenset(T))
def gmap_f(g, jj):
    vals = g_act_values(g, list(F_VALUES[jj])); return [k for k in range(82) if list(F_VALUES[k]) == [int(x) for x in vals]]
check("inverse action maps second source (2,1,11,9,5,28,29,33) to first (5,1,12,9,2,35,33,29) term by term",
      [gmap_f(g_inv, jj) for jj in S28] == [[k] for k in S35])
# capacity (4) algebra
da = d * a / (a + e); de = d * e / (a + e)
check("split d=da+de gives j_a+j_e = d + d^2/(a+e)", sp.simplify(da + da**2 / a + de + de**2 / e - (d + d**2 / (a + e))) == 0)
# radial C, Q
h = sp.Symbol('h')
v = [0] * 14
for s_, jj in zip(SY, S35):
    for k in range(14): v[k] += s_ * sp.Rational(F_VECS[jj][k].numerator, F_VECS[jj][k].denominator)
v[13] += h
Cp = sp.Poly(sp.expand(cubic_C(v)), h); Qs, Cs = Cp.all_coeffs()
A_ = 1 + a + b; B_ = 1 + q + e; N = (q + 1) * (e + 1) * A_ + B_ * a * b; Nj = N - j * e * a
check("C = -(ABc + Nj) at d=1", sp.expand(Cs.subs(d, 1) + (A_ * B_ * c + Nj)) == 0)
S_ = a + e; J_ = j - 1
check("-Q = (q+b+1+S-p)c + 1 - JS + (q+b)(1+S) + pj at d=1", sp.expand(-Qs.subs(d, 1) - ((q + b + 1 + S_ - p) * c + 1 - J_ * S_ + (q + b) * (1 + S_) + p * j)) == 0)
# NUMERIC: stratum sampling for T_40 and the uniform disjunction
T = CONTACT_SETS[40]
nf, nv, nminor = 0, 0, 0
for seed in range(int(os.environ.get("NSEEDS", "25"))):
    r = random_point_on_contact_stratum(T, seed=seed, tries=20, margin=1e-2)
    if r is None: continue
    g = r["gaps"]; x = [g[11], g[17], g[6], g[7], g[15], g[3] - g[12], g[12], g[4]]
    q_, p_, j_, c_, e_, d_, b_, a_ = x
    if d_ < 0:   # second source; use transposed coefficients (e-dl, p-m-dl, j, c, q-m, dl, a-dl, b-m)
        m_ = -d_
        if j_ >= m_: nminor += 1; nf += 1; info(f"seed {seed}", "minor sum"); continue
        dl = m_ - j_
        q_, p_, j_, c_, e_, d_, b_, a_ = e_ - dl, p_ - m_ - dl, j_, c_, q_ - m_, dl, a_ - dl, b_ - m_   # transposed: g28 -> g35 coefficientwise
    nf += 1
    ok = p_ * j_ >= d_ * d_ or (a_ + e_ > 0 and j_ >= d_ + d_ * d_ / (a_ + e_))
    nv += (not ok)
    info(f"seed {seed}", f"normalized q,p,j,c,e,b,a = {[round(t / d_, 4) for t in (q_, p_, j_, c_, e_, b_, a_)]}; pj>=1: {p_*j_>=d_*d_}; Hsum: {a_+e_>0 and j_>=d_+d_*d_/(a_+e_)}")
check(f"NUMERIC: all {nf} stratum points (contact set exactly T_40) satisfy the uniform disjunction ({nminor} minor sums)", nv == 0 and nf > 0)
def phimin(pairs): return min_phi_over_SO4([float(t) for t in comb(pairs)], restarts=30, seed=1)[0]
m1 = phimin([(1, 1), (1, 12), (1, 35)]); m2 = phimin([(Fr(9, 10), 1), (Fr(9, 10), 12), (1, 35)])
check("NUMERIC SOC (3) at p=j=1: min Phi >= -1e-9", m1 >= -1e-9, f"{m1:.3e}; 10% below: {m2:.3e}")
m1 = phimin([(Fr(1, 2), 29), (Fr(1, 2), 2), (2, 12), (1, 35)]); m2 = phimin([(Fr(1, 2), 29), (Fr(1, 2), 2), (Fr(18, 10), 12), (1, 35)])
check("NUMERIC old-H sum (4) at a=e=1/2, j=2: min Phi >= -1e-9", m1 >= -1e-9, f"{m1:.3e}; 10% below: {m2:.3e}")
ex = [f for f in FAIL if not f.startswith("NUMERIC")]
print(f"SUMMARY: {len(FAIL)} failed checks ({len(ex)} exact)")
sys.exit(1 if ex else 0)
