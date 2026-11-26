# Auditor replay for root_moc4_T37_uniform_20261022 (own code; nothing from the supplement is executed or imported).
# Run from this directory: python -I -B root_moc4_T37_uniform_20261022__check.py   (needs ../lib/moc4lib.py, sympy, python-flint, numpy, scipy, mpmath)
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
# ===== root_moc4_T37_uniform_20261022 : paper case T_38 =====
import numpy as np
T = CONTACT_SETS[38]
check("contact set {1,2,3,8,11,12} == paper T_38", T == frozenset({1, 2, 3, 8, 11, 12}))
GEN = [1, 2, 15, 9, 10, 56, 0, 53]  # q p j c e d b a
check("generators(T_38) == {0,1,2,9,10,15,53,56}", sorted(generators(T)) == sorted(GEN), generators(T))
check("only f56 non-minor", [j for j in GEN if not F_IS_MINOR[j]] == [56])
PER = [17, 15, 4, 10, 9, 6, 20, 14]
Mev = [[F_VALUES[j][s] for j in GEN] for s in PER]
check("evaluation of (f1,f2,f15,f9,f10,f56,f0,f53) at perms (17,15,4,10,9,6,20,14) is the identity", Mev == [[int(i == k) for k in range(8)] for i in range(8)])
# contact nullspace has dim 8 and is spanned by these generators
rows = [[P[s][r] for r in range(14)] for s in T]
check("contact equations have rank 6", rank(rows) == 6)
check("the 8 generators lie in the nullspace and have rank 8", all(F_VALUES[j][s] == 0 for j in GEN for s in T) and rank([fvec(j) for j in GEN]) == 8)
# zero-face table: (zero coeff, its gap perm, supplement label, listed image)
ZT = [("q", 17, "T45", {1, 2, 8, 11, 12, 17}), ("j", 4, "T4", {1, 2, 3, 4, 8, 11}),
      ("e", 9, "T7", {1, 2, 3, 8, 9, 11}), ("a", 14, "T10", {1, 2, 3, 8, 12, 14})]
for nm, s, lab, img in ZT:
    inside = img <= (set(T) | {s})
    nus = [nu for nu in range(1, 91) if find_g(CONTACT_SETS[nu], img) is not None]
    check(f"zero {nm}: listed image {sorted(img)} lies in T_38 u {{{s}}}", inside)
    LAB = {"T45": 46, "T4": 4, "T7": 8, "T10": 11}  # supplement historical label -> paper index (cases/index.tsv)
    check(f"zero {nm}: image is in the orbit of paper T_{LAB[lab]} (supplement historical label {lab})", nus == [LAB[lab]], nus)
# radial formulas
q, p, j, c, e, d, b, a, h = sp.symbols('q p j c e d b a h')
syms = [q, p, j, c, e, d, b, a]
v = [0] * 14
for s_, jj in zip(syms, GEN):
    for r in range(14): v[r] += s_ * sp.Rational(F_VECS[jj][r].numerator, F_VECS[jj][r].denominator)
v[13] += h
Cp = sp.Poly(sp.expand(cubic_C(v)), h)
Q1, C1 = Cp.all_coeffs()
S = a + e
Pp = sp.Poly(sp.expand(S * C1 - a * e * (d * S + d ** 2 - j * S)), *syms)
check("cubic is affine in h along kappa", Cp.degree() == 1)
check("S*C - ae(dS + d^2 - jS) has nonnegative coefficients", all(x >= 0 for x in Pp.coeffs()), f"{len(Pp.terms())} terms")
check("... and contains a d^3 and e d^3", Pp.coeff_monomial(a * d ** 3) > 0 and Pp.coeff_monomial(e * d ** 3) > 0)
info("Q (h-coefficient)", Q1)
# Hessian square structure at random rational points of the hard region
def hard(x):
    q_, p_, j_, c_, e_, b_, a_ = x
    if p_ * j_ >= 1: return False
    if c_ > 0 and p_ >= 1 + 1 / c_: return False
    if j_ >= 1 + 1 / (a_ + e_): return False
    if b_ > 0 and j_ >= 1 / e_ + 1 / b_: return False
    return min(q_, j_, e_, a_) > 0 and min(p_, c_, b_) >= 0
def src(x):
    q_, p_, j_, c_, e_, b_, a_ = x
    return comb([(q_, 1), (p_, 2), (j_, 15), (c_, 9), (e_, 10), (1, 56), (b_, 0), (a_, 53)])
rng = random.Random(38)
def rnd_hard():
    while True:
        x = [Fr(rng.randint(1, 400), rng.randint(1, 200)) for _ in range(7)]
        if hard(x): return x
nsq = 0
for _ in range(6):
    G2, _, _ = G2_at(src(rnd_hard()))
    nsq += is_rat_square(G2) and G2 != 0
check("-(Q^2/4) det Hess at h=-C/Q is a nonzero rational square at 6 hard points", nsq == 6)
# exact root count of G on random lines inside the hard region
nl, bad = 0, 0
for k in range(int(os.environ.get("NLINES", "12"))):
    x0 = rnd_hard(); dv = [Fr(rng.randint(-50, 50), 50) for _ in range(7)]
    # find interval in hard region (sampled), then shrink
    ts = np.linspace(-20, 20, 4001)
    ok = [hard([x0[i] + Fr(float(t)).limit_denominator(1000) * dv[i] for i in range(7)]) for t in ts]
    i0 = 2000; lo = i0; hi = i0
    while lo > 0 and ok[lo - 1]: lo -= 1
    while hi < 4000 and ok[hi + 1]: hi += 1
    tlo = Fr(float(ts[lo])).limit_denominator(1000); thi = Fr(float(ts[hi])).limit_denominator(1000)
    if not (hard([x0[i] + tlo * dv[i] for i in range(7)]) and hard([x0[i] + thi * dv[i] for i in range(7)])): continue
    Pf, okdeg = G2_on_line(lambda t: src([x0[i] + t * dv[i] for i in range(7)]), 18, tlo, thi)
    nr = roots_in(Pf, tlo, thi)
    nl += 1; bad += (nr != 0 or not okdeg)
    info(f"line {k}", f"t in [{float(tlo):.3f},{float(thi):.3f}], deg G^2 = {Pf.degree()}, degree certified {okdeg}, roots of G = {nr}")
check(f"NUMERIC(segments chosen numerically, root count exact): G has no real root on {nl} random segments inside the hard region (exact Sturm count; segments chosen numerically)", bad == 0 and nl > 0)
# NUMERIC part for root_moc4_T37_uniform_20261022
import numpy as np
T = CONTACT_SETS[38]
GEN = [1, 2, 15, 9, 10, 56, 0, 53]; PER = [17, 15, 4, 10, 9, 6, 20, 14]
def capacities(x):
    q_, p_, j_, c_, e_, d_, b_, a_ = x
    out = []
    if p_ * j_ >= d_ * d_: out.append("SOC pj>=1")
    if c_ > 0 and p_ >= d_ + d_ * d_ / c_: out.append("H p>=1+1/c")
    if a_ + e_ > 0 and j_ >= d_ + d_ * d_ / (a_ + e_): out.append("Hsum j>=1+1/(a+e)")
    if e_ > 0 and b_ > 0 and j_ >= d_ * d_ / e_ + d_ * d_ / b_: out.append("chain j>=1/e+1/b")
    return out
nfound, nviol = 0, 0
for seed in range(int(os.environ.get("NSEEDS", "25"))):
    r = random_point_on_contact_stratum(T, seed=seed, tries=20, margin=1e-2)
    if r is None: continue
    g = r["gaps"]; x = [g[s] for s in PER]
    if x[5] < 1e-9: continue
    nfound += 1
    caps = capacities(x)
    if not caps: nviol += 1
    info(f"stratum seed {seed}", f"d={x[5]:.3g} q,p,j,c,e,b,a/d={[round(v / x[5], 4) for v in x[:5] + x[6:]]} capacities={caps}")
check(f"NUMERIC: every one of {nfound} stratum points (contact set exactly T_38) satisfies one of (2)-(5)", nviol == 0 and nfound > 0)
def phimin(pairs, seed=0):
    c = [float(v) for v in comb(pairs)]
    return min_phi_over_SO4(c, restarts=30, seed=seed)[0]
tests = [("(2) at p=j=1", [(1, 2), (1, 15), (1, 56)], [(0.9, 2), (0.9, 15), (1, 56)]),
         ("(3) c=1,p=2", [(1, 9), (2, 2), (1, 56)], [(1, 9), (1.8, 2), (1, 56)]),
         ("(4) a=e=1/2,j=2", [(0.5, 53), (0.5, 10), (2, 15), (1, 56)], [(0.5, 53), (0.5, 10), (1.8, 15), (1, 56)]),
         ("(5) e=b=1,j=2", [(1, 10), (1, 0), (2, 15), (1, 56)], [(1, 10), (1, 0), (1.8, 15), (1, 56)])]
for nm, at, below in tests:
    m1 = phimin(at); m2 = phimin(below)
    check(f"NUMERIC capacity {nm}: min Phi over SO(4) at threshold >= -1e-9", m1 >= -1e-9, f"{m1:.3e}; 10% below threshold: {m2:.3e}")
ex = [f for f in FAIL if not f.startswith("NUMERIC")]
print(f"SUMMARY: {len(FAIL)} failed checks ({len(ex)} exact)")
sys.exit(1 if ex else 0)
