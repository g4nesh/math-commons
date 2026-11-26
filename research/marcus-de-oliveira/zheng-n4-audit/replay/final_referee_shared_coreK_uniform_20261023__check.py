"""Exact replay of final_referee_shared_coreK_uniform_20261023 (seven-contact set K = AUX7[16]).

Own code only (../lib/moc4lib.py + python-flint + sympy). Nothing from the supplement is read
or executed; in particular the degree-nine polynomial G is re-derived here from the 14x14
Hessian determinant, not loaded from their certificate.
Prints one line per check; exits 1 if an EXACT check fails.
Run: python -I -B final_referee_shared_coreK_uniform_20261023__check.py [n_numeric]
"""
import os, sys, itertools, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc4lib import *            # noqa
from fractions import Fraction as Fr
import flint, sympy

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + str(info)) if info != "" else ""))
    if not ok:
        FAIL.append(name)

def hessian_forms():
    ms, Cx = cubic_sympy()
    H = []
    for i in range(14):
        row = []
        for j in range(14):
            h = sympy.expand(sympy.diff(Cx, ms[i], ms[j]))
            lf = {k: int(h.coeff(ms[k])) for k in range(14) if h.coeff(ms[k]) != 0}
            assert sympy.expand(h - sum(cc * ms[k] for k, cc in lf.items())) == 0
            row.append(lf)
        H.append(row)
    return H

def bareiss_det(M, zero, one):
    M = [list(r) for r in M]; n = len(M); sign = 1; prev = one; ndiv = 0
    for k in range(n - 1):
        if M[k][k] == zero:
            sw = next((r for r in range(k + 1, n) if M[r][k] != zero), None)
            if sw is None:
                return zero, ndiv
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev
                ndiv += 1
        prev = M[k][k]
    return (M[n - 1][n - 1] if sign == 1 else -M[n - 1][n - 1]), ndiv

# ================= source cone =================
K = [0, 1, 2, 9, 10, 14, 17]
check("K equals paper Table 5 seven-set #16", [k for k, v in AUX7.items() if v == frozenset(K)] == [16])
order = [0, 2, 12, 13, 15, 21, 26]
check("generators(K) = [0,2,12,13,15,21,26] (simplicial)", generators(K) == order, generators(K))
check("f26 is the only non-minor generator", [j for j in order if not F_IS_MINOR[j]] == [26])
check("rank of P rows at K = 7", rank([P[s] for s in K]) == 7)
rows = [23, 15, 8, 7, 4, 16, 11]
check("profile values at rows [23,15,8,7,4,16,11] form I7",
      [[F_VALUES[j][s] for j in order] for s in rows] == [[int(i == k) for k in range(7)] for i in range(7)])
check("f1 = n12,34, f9 = n14,23, f22 = h3412,+", (f_name(1), f_name(9), f_name(22)) == ("n12,34", "n14,23", "h3412,+"))

def fimg(G_, j):
    vv = g_act_values(G_, list(F_VALUES[j])); r = [k for k in range(82) if list(F_VALUES[k]) == vv]
    return r[0] if r else None
socg = {}
for pair in [(2, 12), (2, 15), (12, 15)]:
    found = [G_ for G_ in G_ELEMENTS if fimg(G_, 22) == 26 and {fimg(G_, 1), fimg(G_, 9)} == set(pair)]
    socg[pair] = found[0] if found else None
check("(2): group elements carry (f1,f9;f22) to (fi,fj;f26) for all three pairs in {2,12,15}",
      all(socg[p] is not None for p in socg), {p: socg[p] for p in socg})
stab = [G_ for G_ in G_ELEMENTS if g_act_set(G_, K) == frozenset(K)]
acts = set(tuple(fimg(G_, j) for j in order) for G_ in stab)
ok = len(stab) == 6 and all(p[6] == 26 for p in acts) and len(set((p[1], p[2], p[4]) for p in acts)) == 6
# pairing a<->q (f2<->f0), b<->v (f12<->f13), c<->w (f15<->f21)
pair = {2: 0, 12: 13, 15: 21}
ok = ok and all(p[order.index(pair[m])] == pair[p[order.index(m)]] for p in acts for m in (2, 12, 15))
check("stabilizer of K: 6 elements, all of S3 on (a,b,c), same permutation on (q,v,w), d fixed", ok)

# ================= cubic, lift and G =================
names = ('q', 'a', 'b', 'v', 'c', 'w', 'd', 'z')
ctx = flint.fmpz_mpoly_ctx.get(names, 'degrevlex')
q, a, b, v, c, w, d, z = ctx.gens(); zero = ctx.from_dict({}); one = zero + 1
coef = {0: q, 2: a, 12: b, 13: v, 15: c, 21: w, 26: d}
x = [sum((int(F_VECS[j][k]) * coef[j] for j in coef), zero) for k in range(14)]
x[13] += z
Cp = (q + d) * (v + d) * (w + d) + a * b * c; Qp = a * b + a * c + b * c + d * d
check("(5): cubic(g + z kappa) = C + z Q", cubic_C(x) - (Cp + z * Qp) == zero)
Hf = hessian_forms()
H = [[sum((cc * x[k] for k, cc in Hf[i][j].items()), zero) for j in range(14)] for i in range(14)]
t0 = time.time()
D, nd = bareiss_det(H, zero, one)
check("Bareiss: 819 exact divisions", nd == 819, "time %.0fs" % (time.time() - t0))
zc = {}
for mon, cf in zip(D.monoms(), D.coeffs()):
    m2 = list(mon); k = m2[7]; m2[7] = 0; zc.setdefault(k, {})[tuple(m2)] = cf
zc = {k: ctx.from_dict(vv) for k, vv in zc.items()}
check("det H(g + z kappa) has z-degree <= 4", max(zc) <= 4, max(zc))
N = sum((zc[k] * (-Cp) ** k * Qp ** (4 - k) for k in zc), zero)
try:
    G2 = N / (-4 * Qp ** 2); G = G2.sqrt()
    mono_d9 = (0, 0, 0, 0, 0, 0, 9, 0)
    if G.coefficient(G.monoms().index(mono_d9)) < 0:
        G = -G
    ok6 = (-4 * Qp ** 2 * G ** 2 == N)
except Exception as ex:
    ok6 = False; G = None; print("  ", ex)
check("(6): Q^4 det H(g - C kappa/Q) = -4 Q^2 G^2 for a polynomial G", ok6)
check("G: 746 terms, multidegrees (3,2,2,3,2,3,9), d^9 coefficient +1, homogeneous degree 9",
      len(G) == 746 and tuple(G.degrees())[:7] == (3, 2, 2, 3, 2, 3, 9)
      and G.coefficient(G.monoms().index(mono_d9)) == 1 and all(sum(m) == 9 for m in G.monoms()),
      (len(G), G.degrees()))
# symmetry of G under the induced S3
def sub(p, perm):   # perm: dict var index -> new var
    gens_ = list(ctx.gens()); return p.compose(*[gens_[perm.get(i, i)] for i in range(8)])
sw_ab = {1: 2, 2: 1, 0: 3, 3: 0}; sw_bc = {2: 4, 4: 2, 3: 5, 5: 3}
check("G invariant under (a b)(q v) and (b c)(v w)", sub(G, sw_ab) == G and sub(G, sw_bc) == G)
check("(10): G(q,a,0,v,0,w,d) = {a(v+w+d)+(v+d)(w+d)}(v+d)^2(w+d)^2(q+d)^3",
      G.compose(q, a, zero, v, zero, w, d, z) ==
      (a * (v + w + d) + (v + d) * (w + d)) * (v + d) ** 2 * (w + d) ** 2 * (q + d) ** 3)

# ================= the two charts =================
cx = flint.fmpz_mpoly_ctx.get(('Q0', 'V', 'W', 'X', 'Y', 'Z'), 'degrevlex')
Q0, V, W, Xc, Yc, Zc = cx.gens(); cz = cx.from_dict({}); c1 = cz + 1
Gd = {}   # group G (with d = 1) by (deg_a, deg_b, deg_c): polynomial in Q0, V, W
for mon, cf in zip(G.monoms(), G.coeffs()):
    iq, ia, ib, iv, ic, iw, idd, _ = mon
    key = (ia, ib, ic)
    Gd[key] = Gd.get(key, cz) + int(cf) * Q0 ** iq * V ** iv * W ** iw
kabc = max(sum(k) for k in Gd)
# balanced chart (8)
Dd = 1 + Xc + Yc + Zc
Pb = sum((Gd[k] * (Yc + Zc) ** k[0] * (Xc + Zc) ** k[1] * (Xc + Yc) ** k[2] * Dd ** (kabc - sum(k)) for k in Gd), cz)
cfs = [int(t) for t in Pb.coeffs()]
check("(8) balanced chart: D^%d * G has all coefficients > 0, constant coefficient 1" % kabc,
      all(t > 0 for t in cfs) and Pb.coefficient(Pb.monoms().index((0,) * 6)) == 1,
      "nonzero coefficients: %d (document: 5373)" % len(cfs))
# dominant chart (9)
L = 1 + Xc
ka = max(k[0] for k in Gd); kbc = max(k[1] + k[2] for k in Gd)
Pd = sum((Gd[k] * (1 + Yc * L ** 2) ** k[0] * (1 + Yc) ** (ka - k[0]) * Zc ** k[1] * (1 + Zc) ** (kbc - k[1] - k[2])
          * L ** (kabc - sum(k)) for k in Gd), cz)
cfs = [int(t) for t in Pd.coeffs()]
mX = (0, 0, 0, 1, 0, 0)
coefX = Pd.coefficient(Pd.monoms().index(mX)) if mX in Pd.monoms() else 0
const = Pd.coefficient(Pd.monoms().index((0,) * 6)) if (0,) * 6 in Pd.monoms() else 0
check("(9) dominant chart: L^%d (1+Y)^%d (1+Z)^%d * G has all coefficients > 0 and a pure-X monomial"
      % (kabc, ka, kbc), all(t > 0 for t in cfs) and coefX > 0,
      "nonzero coefficients: %d (document: 7297); coeff of X = %s; constant = %s" % (len(cfs), coefX, const))
# chart inverses (sympy, exact)
xs, ys, zs_, aS = sympy.symbols('x y z0 a', positive=True)
Xs, Ys, Zs = [t / (1 - xs - ys - zs_) for t in (xs, ys, zs_)]
Ds = 1 + Xs + Ys + Zs
check("(8) inverse: (Y+Z)/D = y+z0 etc. (so a,b,c recovered)",
      all(sympy.simplify(e1 - e2) == 0 for e1, e2 in [((Ys + Zs) / Ds, ys + zs_), ((Xs + Zs) / Ds, xs + zs_), ((Xs + Ys) / Ds, xs + ys)]))
S_, bq, cq = sympy.symbols('S b c', positive=True)
Xd = 1 / S_ - 1; Yd = (aS * S_ - S_ ** 2) / (1 - aS * S_); Zd = bq / cq; Ld = 1 + Xd
ok9 = (sympy.simplify((1 + Yd * Ld ** 2) / (Ld * (1 + Yd)) - aS) == 0 and
       sympy.simplify(Zd / (Ld * (1 + Zd)) - bq * S_ / (bq + cq)) == 0 and
       sympy.simplify(1 / (Ld * (1 + Zd)) - cq * S_ / (bq + cq)) == 0)
check("(9) inverse: chart returns a, and b,c with b+c = S, b/c = Z", ok9)

# ================= (3) resource allocations =================
A_, B_, C_ = sympy.symbols('a b c', nonnegative=True)
dab, dac, dbc = (A_ + B_ - C_) / 2, (A_ + C_ - B_) / 2, (B_ + C_ - A_) / 2
ok3 = (sympy.expand(dab + dac - A_) == 0 and sympy.expand(dab + dbc - B_) == 0 and sympy.expand(dac + dbc - C_) == 0
       and sympy.expand(dab + dac + dbc - (A_ + B_ + C_) / 2) == 0)
s_ = B_ + C_
xa, xb = A_ * B_ / s_, A_ * C_ / s_
ok3 = ok3 and sympy.simplify(xa + xb - A_) == 0 and sympy.simplify(sympy.sqrt(xa * B_) + sympy.sqrt(xb * C_) - sympy.sqrt(A_ * s_)) == 0
check("(3) balanced allocations sum to a,b,c with total s/2; dominant split gives sqrt(a(b+c))", ok3)

# ================= NUMERIC =================
import numpy as np
from scipy.optimize import minimize
nn = int(sys.argv[1]) if len(sys.argv) > 1 else 40
f1, f9, f22 = [np.array([float(t) for t in F_VECS[j]]) for j in (1, 9, 22)]
def soc_obj(th):
    O = so4_from_params(th)
    p1, p9, p22 = Phi_np(f1, O), Phi_np(f9, O), Phi_np(f22, O)
    return p22 + 2 * math.sqrt(max(p1, 0) * max(p9, 0))
rng = np.random.default_rng(7)
best = 1e9
for r in range(nn):
    res = minimize(soc_obj, rng.normal(size=6) * 2, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-12, fatol=1e-14))
    best = min(best, res.fun)
print("NUMERIC SOC atom (f1,f9;f22): min over SO(4) of Phi_f22 + 2 sqrt(Phi_f1 Phi_f9) = %.3e (%d Nelder-Mead restarts, seed 7)" % (best, nn))
# the conic claim: a f2 + b f12 + c f15 + H3(a,b,c) f26 has Phi >= 0
def H3(aa, bb, cc):
    s = aa + bb + cc; m = max(aa, bb, cc)
    return s / 2 if 2 * m <= s else math.sqrt(m * (s - m))
FV = {j: np.array([float(t) for t in F_VECS[j]]) for j in order}
worst = 1e9
for r in range(12):
    abc = rng.exponential(size=3)
    if r % 3 == 0:
        abc[0] = abc[1] + abc[2] + rng.exponential()
    cvec = abc[0] * FV[2] + abc[1] * FV[12] + abc[2] * FV[15] + H3(*abc) * FV[26]
    m = min_phi_over_SO4(cvec, restarts=10, seed=r)
    m = m[0] if isinstance(m, tuple) else m
    worst = min(worst, m / (abc.sum()))
print("NUMERIC conic claim (3): min Phi of a f2+b f12+c f15+H3 f26 over SO(4), 12 random (a,b,c), 10 restarts each: %.3e (normalized by a+b+c)" % worst)
print("SUMMARY: %d EXACT failures" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
