"""Exact replay of blind_shared_coreJ_uniform_20261023 (seven-contact set J = AUX7[22]).

Own code only. Uses the shared library ../lib/moc4lib.py (paper objects rebuilt from the
paper's definitions) and python-flint for exact multivariate polynomials. Nothing from the
supplement is read or executed. Prints one line per check; exits 1 if an EXACT check fails.
Run: python -I -B blind_shared_coreJ_uniform_20261023__check.py [numeric_samples]
"""
import os, sys, itertools, time
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

# ---------- Hessian of the paper's cubic (A.2) as integer linear forms ----------
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

def hessian_at(Hf, w, zero):
    return [[sum((cc * w[k] for k, cc in Hf[i][j].items()), zero) for j in range(14)] for i in range(14)]

def bareiss_det(M, zero, one):
    """Fraction-free determinant; every division is exact or raises."""
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

def zcoeffs(D, zi, nvars):
    """coefficients of D in the variable with index zi"""
    out = {}
    for mon, cf in zip(D.monoms(), D.coeffs()):
        k = mon[zi]; m2 = list(mon); m2[zi] = 0
        out.setdefault(k, {})[tuple(m2)] = cf
    return {k: D.context().from_dict(v) for k, v in out.items()}

# ================= A/B: source cone =================
J = [0, 2, 3, 7, 11, 12, 23]
check("J equals paper Table 5 seven-set #22", [k for k, v in AUX7.items() if v == frozenset(J)] == [22])
gens = generators(J)
check("generators(J) = [2,7,8,12,15,18,39,40]", gens == [2, 7, 8, 12, 15, 18, 39, 40], gens)
check("minor flags: all minor except f18, f40",
      [j for j in gens if not F_IS_MINOR[j]] == [18, 40])
check("rank of P rows at J = 7", rank([P[s] for s in J]) == 7)
seven = [2, 7, 8, 12, 15, 39, 18]
check("the seven profiles (2) are independent (rank 7)", rank([list(F_VECS[j]) for j in seven]) == 7)
check("all seven profiles vanish on J", all(F_VALUES[j][s] == 0 for j in seven for s in J))

names = ('a', 'b', 'c', 'e', 'f', 'g', 'd', 'z')
ctx = flint.fmpz_mpoly_ctx.get(names, 'degrevlex')
a, b, c, e, f, g, d, z = ctx.gens(); zero = ctx.from_dict({}); one = zero + 1
coef = {2: a, 7: b, 8: c, 12: e, 15: f, 39: g, 18: d}
vals = [sum((F_VALUES[j][s] * coef[j] for j in coef), zero) for s in range(24)]
check("measurements (3), first row", [vals[s] for s in (21, 20, 9, 6, 18, 1)] == [a, b, c, e, f, g])
check("measurements (3), second row", [vals[s] for s in (13, 16, 5, 17)] == [a + d, b + 2 * d, f + d, g + d])
circ = all(F_VALUES[18][s] + F_VALUES[40][s] == F_VALUES[2][s] + 2 * F_VALUES[7][s] + F_VALUES[15][s] + F_VALUES[39][s]
           for s in range(24))
check("circuit (4) f18+f40 = f2+2f7+f15+f39 on all 24 values", circ)
gel = ((1, 2, 0, 3), (1, 2, 0, 3), 1)   # pi -> column[pi^{-1}(row[i])]
ok_conv = all(perm_index(tuple((1, 2, 0, 3)[inverse(p)[(1, 2, 0, 3)[i]]] for i in range(4))) == g_index_map(gel)[s]
              for s, p in enumerate(PERMS))
check("stated stabilizer convention = lib element (rho=col, rho'=row, eps=1)", ok_conv)
def fimg(G_, j):
    v = g_act_values(G_, list(F_VALUES[j])); r = [k for k in range(82) if list(F_VALUES[k]) == v]
    return r[0] if r else None
check("stabilizer fixes J, fixes six minors, swaps f18<->f40",
      g_act_set(gel, J) == frozenset(J) and all(fimg(gel, j) == j for j in (2, 7, 8, 12, 15, 39))
      and fimg(gel, 18) == 40 and fimg(gel, 40) == 18)

# ================= cubic and lift =================
for j in coef:
    assert all(x.denominator == 1 for x in F_VECS[j])
w = [sum((int(F_VECS[j][k]) * coef[j] for j in coef), zero) for k in range(14)]
check("w0 has no kappa component (m13 = 0)", w[13] == zero)
w[13] = w[13] + z
X = c + e; Y = f + g; Pp = c * e; R = f * g; T = (f + d) * (g + d)
Cp = (b + X + d) * (a * (Y + d) + T) + Pp * (Y + d)
Qp = a * (c + e + f + g + d) + c * e + c * g + c * d + e * f + e * d + f * g + f * d + g * d + d * d
check("(7)-(8): cubic(w0 + z kappa) = C + z Q", cubic_C(w) - (Cp + z * Qp) == zero)
check("C, Q have positive coefficients; C contains d^3, Q contains d^2",
      all(x > 0 for x in Cp.coeffs()) and all(x > 0 for x in Qp.coeffs())
      and Cp.coefficient(Cp.monoms().index((0,) * 6 + (3, 0))) == 1 and (0,) * 6 + (2, 0) in Qp.monoms())

A = (Y + d) * (b * b * Pp + d * d * (b + d) * (b + d + X)) - 4 * Pp * R * (b + d)
B = b * b * (Pp + d * d) * T + b * d * d * ((X + 2 * d) * T + Pp * (Y + d)) + d ** 3 * ((X + d) * T + Pp * (Y + d))
G = a * A + B
check("B has 38 terms, all positive, d^6 coefficient 1",
      len(B) == 38 and all(x > 0 for x in B.coeffs()) and B.coefficient(B.monoms().index((0,) * 6 + (6, 0))) == 1, len(B))
negG = zero + sum((ctx.from_dict({m: cf}) for m, cf in zip(G.monoms(), G.coeffs()) if cf < 0), zero)
check("G has 64 terms; negative part = -4acefg(b+d)", len(G) == 64 and negG == -4 * a * c * e * f * g * (b + d), len(G))

# ================= (10) exact determinant identity =================
Hf = hessian_forms()
H = hessian_at(Hf, w, zero)
t0 = time.time()
D, nd = bareiss_det(H, zero, one)
check("Bareiss: 819 exact divisions for the 14x14 determinant", nd == 819, "time %.1fs" % (time.time() - t0))
Dz = zcoeffs(D, 7, 8)
check("det H(w0 + z kappa) has degree <= 4 in z", max(Dz) <= 4, max(Dz))
lhs = sum((Dz[k] * (-Cp) ** k * Qp ** (4 - k) for k in Dz), zero)
check("(10): Q^4 det H(w0 - C kappa/Q) = -4 Q^2 C^2 G^2", lhs == -4 * Qp ** 2 * Cp ** 2 * G ** 2)

# ================= (11) Schur borders =================
I = list(range(12))
borders = [(I + [12], I + [12]), (I + [13], I + [13]), (I + [12], I + [13])]
allok = True; ndiv_tot = 0
for rs, cs in borders:
    Mb = [[H[r][s] for s in cs] for r in rs]
    Db, ndb = bareiss_det(Mb, zero, one); ndiv_tot += ndb
    Dbz = zcoeffs(Db, 7, 8) if Db != zero else {0: zero}
    if max(Dbz) > 4:
        allok = False; print("  border z-degree", max(Dbz)); continue
    num = sum((Dbz[k] * (-Cp) ** k * Qp ** (4 - k) for k in Dbz), zero)
    try:
        _ = num / G
    except Exception:
        allok = False
check("(11): Q^4 * each of the three 13x13 border determinants is divisible by G (exact)", allok,
      "border divisions %d" % ndiv_tot)
gAB = A.gcd(B)
check("gcd(A, B) = 1", gAB == one, str(gAB))
fac = G.factor()
check("G irreducible over Q (flint factor)", len(fac[1]) == 1 and fac[1][0][1] == 1, str(fac[0]))

# ================= (12), (19) sample =================
smp = dict(a=Fr(5, 3), b=Fr(0), c=Fr(2), e=Fr(2), f=Fr(2), g=Fr(2), d=Fr(1))
def ev(p):
    return sum(Fr(int(cf)) * eval_monomial(m) for m, cf in zip(p.monoms(), p.coeffs()))
order = ['a', 'b', 'c', 'e', 'f', 'g', 'd', 'z']
def eval_monomial(m, zval=None):
    out = Fr(1)
    for k, ex in enumerate(m):
        if ex:
            out *= (smp[order[k]] if order[k] != 'z' else zval_g[0]) ** int(ex)
    return out
zval_g = [Fr(0)]
Cs, Qs = ev(Cp), ev(Qp); zs = -Cs / Qs; zval_g[0] = zs
check("(12): C=320/3, Q=40, z=-8/3, G=0 at the sample", (Cs, Qs, zs, ev(G)) == (Fr(320, 3), Fr(40), Fr(-8, 3), 0))
wnum = [ev(x) for x in w]
stated = [Fr(-4), Fr(20, 3), Fr(3), Fr(-14, 3), Fr(2), Fr(-5), Fr(-5, 3), Fr(2), Fr(2), Fr(5, 3), Fr(2), Fr(-14, 3), Fr(3), Fr(-8, 3)]
check("(19): full vector at the sample", wnum == stated)
Hs = [[sum(Fr(cc) * wnum[k] for k, cc in Hf[i][j].items()) for j in range(14)] for i in range(14)]
H3 = [[3 * x for x in r] for r in Hs]
check("3H(w) integral", all(x.denominator == 1 for r in H3 for x in r))
dI = det([r[:12] for r in H3[:12]])
check("(19): det(3H(w)[I,I]) = -310542336000", dI == -310542336000, dI)
check("rank H(w) = 12 at the sample", rank(Hs) == 12, rank(Hs))
mins12 = [det([[Hs[r][s] for s in S_] for r in S_]) for S_ in itertools.combinations(range(14), 12)]
check("all 91 principal 12-minors at the sample are <= 0, some < 0",
      all(x <= 0 for x in mins12) and any(x < 0 for x in mins12), "nonzero: %d" % sum(1 for x in mins12 if x != 0))

# ================= (17)-(18) branch parametrization (d = 1) =================
bs, cs_, es, xi, ze = sympy.symbols('b c e xi zeta', positive=True)
k0 = bs ** 2 / (4 * (bs + 1)) + (bs + 1 + cs_ + es) / (4 * cs_ * es)
fs = k0 + xi; gs = k0 * (k0 + xi + 1) / xi + ze
Xs, Ys, Ps, Rs = cs_ + es, fs + gs, cs_ * es, fs * gs
As = (Ys + 1) * (bs ** 2 * Ps + (bs + 1) * (bs + 1 + Xs)) - 4 * Ps * Rs * (bs + 1)
check("(18): A = -4ce(b+1) xi zeta on the parametrization", sympy.simplify(As + 4 * cs_ * es * (bs + 1) * xi * ze) == 0)
ff, gg = sympy.symbols('f g', positive=True)
A_fg = (ff + gg + 1) * (bs ** 2 * Ps + (bs + 1) * (bs + 1 + Xs)) - 4 * Ps * ff * gg * (bs + 1)
check("(17): A = -4ce(b+1)[fg - kappa0(f+g+1)] (so A<0 iff fg > kappa0(f+g+1))",
      sympy.simplify(A_fg + 4 * cs_ * es * (bs + 1) * (ff * gg - k0 * (ff + gg + 1))) == 0)

# ================= NUMERIC: cited Hessian inertia theorem (1) =================
import numpy as np
Hk = np.zeros((14, 14, 14))
for i in range(14):
    for j in range(14):
        for k, cc in Hf[i][j].items():
            Hk[k, i, j] = cc
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rng = np.random.default_rng(20261125)
subsets = list(itertools.combinations(range(14), 12))
worst = 1.0; ranks = {}
for t in range(N):
    av = rng.normal(size=4) + 1j * rng.normal(size=4); bv = rng.normal(size=4) + 1j * rng.normal(size=4)
    if t % 3 == 0:
        av = av.real + 0j
    al = np.exp(1j * rng.uniform(0, 2 * np.pi))
    wv = 2 * np.real(matching_vector_np(av, bv, al))
    Hn = np.tensordot(wv, Hk, axes=1)
    sv = np.linalg.svd(Hn, compute_uv=False)
    r = int(np.sum(sv > 1e-9 * sv[0])); ranks[r] = ranks.get(r, 0) + 1
    sc = np.prod(sv[:12])
    worst = min(worst, min(np.linalg.det(Hn[np.ix_(S_, S_)]) / sc for S_ in subsets))
print("NUMERIC Hessian-inertia (1) on w = 2Re x, x = matching vector of random complex spectra: samples %d, seed 20261125, "
      "min normalized 12-principal minor %.3e, numerical ranks %s" % (N, worst, ranks))

# ================= NUMERIC: realizable points on the J stratum have rho = 0 =================
inv7 = None
found = 0; maxrho = 0.0
for sd in range(8):
    pt = random_point_on_contact_stratum(J, seed=sd, tries=20)
    if pt is None:
        continue
    found += 1
    vv = np.array(pt["gaps"], dtype=float)
    # rho = g13 - g21 (= (a+rho) - a)
    maxrho = max(maxrho, abs(vv[13] - vv[21]) / max(1e-300, np.max(np.abs(vv))))
print("NUMERIC stratum search on J: points found %d of 8 seeds, max |rho|/max gap = %.2e" % (found, maxrho))

print("SUMMARY: %d EXACT failures" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
