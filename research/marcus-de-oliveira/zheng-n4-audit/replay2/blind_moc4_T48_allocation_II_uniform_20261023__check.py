"""Second-pass exact replay for blind_moc4_T48_allocation_II_uniform_20261023 (paper case T_49 cell).

Everything is rebuilt from the paper's definitions through the audit library moc4lib
(lib/moc4lib.py): generators f_j (Appendix D), the cubic C (A.2), kappa. Nothing from the
supplement is executed, imported or read. Run from this directory:
    python -I -B blind_moc4_T48_allocation_II_uniform_20261023__check.py
Prints one line per check, exits nonzero if an EXACT check fails. About 90 s.
"""
import os, sys, time, itertools
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
import flint

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" : " + str(info) if info != "" else ""), flush=True)
    if not ok:
        FAIL.append(name)

t0 = time.time()
T = frozenset({0, 3, 5, 6, 9, 12})
roles = dict(q=7, p=5, j=29, c=39, e=2, d=40, b=41, a=0)
order = "qpjcedba"

# ---- A. contact set and roles -------------------------------------------------------
check("contact set equals paper T_49 (identity group element)", L.CONTACT_SETS[49] == T)
gens = L.generators(T)
check("generators(T_49) = [0,1,2,5,7,9,29,39,40,41]", gens == [0, 1, 2, 5, 7, 9, 29, 39, 40, 41], gens)
check("the eight roles all vanish on T", all(roles[k] in gens for k in order))
check("only f40 (d) is non-minor among the 10 generators",
      [j for j in gens if not L.F_IS_MINOR[j]] == [40])
M8 = [L.F_VECS[roles[k]] for k in order]
check("the eight role vectors are linearly independent (rank 8)", L.rank(M8) == 8)
def comb(co):
    return [sum(co.get(k, 0) * L.F_VECS[roles[k]][i] for k in order) for i in range(14)]
check("f1 = q + c + a - d (in role vectors)", comb(dict(q=1, c=1, a=1, d=-1)) == list(L.F_VECS[1]))
check("f9 = q + e + b - d (in role vectors)", comb(dict(q=1, e=1, b=1, d=-1)) == list(L.F_VECS[9]))

# ---- C. pointwise allocation algebra (sympy, exact) ----------------------------------
a, e, j, q, b, c, p = sp.symbols("a e j q b c p", positive=True)
u = (1 + 2 / j) / (2 * (1 / a + 1 / e + 1 / j)); v = 1 - u; Q = q + b
h = v - Q / 2
jh, jq, jb = j * h / v, j * (q / 2) / v, j * (b / 2) / v
check("j resources sum to j", sp.simplify(jh + jq + jb - j) == 0)
need = (u**2 / a + u**2 / e) + (h + h**2 / jh) + ((q / 2)**2 / q + (q / 2)**2 / jq) + ((b / 2)**2 / b + (b / 2)**2 / jb)
pstar_doc = u**2 * (1 / a + 1 / e) + v - Q / 4 + v**2 / j
check("sum of atom requirements = p* (first display)", sp.simplify(need - pstar_doc) == 0)
check("non-minor currents sum to d = 1", sp.simplify(u + h + q / 2 + b / 2 - 1) == 0)
check("u(1/a+1/e) - v/j = 1/2", sp.simplify(u * (1 / a + 1 / e) - v / j - sp.Rational(1, 2)) == 0)
check("p* = 1 - u/2 + v/j - Q/4", sp.simplify(pstar_doc - (1 - u / 2 + v / j - Q / 4)) == 0)
check("p* - (1/2 + v/j) = (2v - Q)/4", sp.simplify(pstar_doc - (sp.Rational(1, 2) + v / j) - (2 * v - Q) / 4) == 0)

# chart
Ls, As, Bs, C1s, Ss, Ps = sp.symbols("L A B C1 S P", positive=True)
uc = Ls / (1 + Ls); vc = 1 / (1 + Ls)
ac = 2 * uc * (1 + As); jc = Ss
ec = 2 * uc * Ss * (1 + As) / (As * Ss + 2 * vc * (1 + As))
Qc = 2 * vc * Bs / (1 + Bs); qc = Qc * C1s / (1 + C1s); bc = Qc / (1 + C1s)
psc = 1 - uc / 2 + vc / Ss - Qc / 4
pc = Ps / (1 + Ps) * psc
check("chart: u recomputed from (a,e,j) equals L/(1+L)",
      sp.simplify(u.subs({a: ac, e: ec, j: jc}) - uc) == 0)
check("chart: e denominator times (1+L) = E_N",
      sp.simplify((As * Ss + 2 * vc * (1 + As)) * (1 + Ls) - (As * Ss * (1 + Ls) + 2 * (1 + As))) == 0)
check("chart inverse: L=u/v, A=a/(2u)-1, B=Q/(2v-Q), C1=q/b, P=p/(p*-p)",
      all(sp.simplify(x) == 0 for x in [uc / vc - Ls, ac / (2 * uc) - 1 - As, Qc / (2 * vc - Qc) - Bs,
                                         qc / bc - C1s, pc / (psc - pc) - Ps]))

# ---- C. radial necessary condition: C(c + r kappa) = C + r Q0, G --------------------
names = ['q', 'p', 'j', 'c', 'e', 'b', 'a', 'r']
ctx = flint.fmpz_mpoly_ctx.get(names, 'lex')
V = dict(zip(names, ctx.gens())); one = ctx.from_dict({(0,) * 8: 1}); zero = ctx.from_dict({})
coef = {k: (V[k] if k != 'd' else one) for k in order}
m = [zero] * 14
for k in order:
    for i in range(14):
        f = L.F_VECS[roles[k]][i]
        assert f.denominator == 1
        if f:
            m[i] = m[i] + int(f) * coef[k]
mm = m[:13] + [m[13] + V['r']]
ms, Cx = L.cubic_sympy()
Hs = sp.hessian(Cx, ms)
def tofl(expr, vals):
    expr = sp.expand(expr); out = zero
    if expr == 0:
        return out
    for mon, cf in sp.Poly(expr, *ms).terms():
        t = one * int(cf)
        for ex, vv in zip(mon, vals):
            for _ in range(ex):
                t = t * vv
        out = out + t
    return out
Cc = tofl(Cx, m); Q0 = tofl(sp.diff(Cx, ms[13]), m)
qq, pp, jj, cc, ee = V['q'], V['p'], V['j'], V['c'], V['e']
check("Q0 = m0 m6 - m1 m5 at the source = (c+1)(1+e+j) - p(e+j)", Q0 == (cc + 1) * (1 + ee + jj) - pp * (ee + jj))
check("C(c + r kappa) = C(c) + r Q0", tofl(Cx, mm) == Cc + V['r'] * Q0)
W = (ee + jj) * Cc - ee * jj * Q0
check("(e+j)C - ejQ0 has 42 terms, all positive, containing e and j",
      len(W) == 42 and all(x > 0 for x in W.coeffs()) and W.to_dict().get((0, 0, 1, 0, 0, 0, 0, 0)) and W.to_dict().get((0, 0, 0, 0, 1, 0, 0, 0)), len(W))
H = [[tofl(Hs[i, k], mm) for k in range(14)] for i in range(14)]
def bareiss(Mx):
    Mx = [row[:] for row in Mx]; n = len(Mx); sign = 1; prev = one
    for k in range(n - 1):
        if Mx[k][k] == 0:
            for i in range(k + 1, n):
                if Mx[i][k] != 0:
                    Mx[k], Mx[i] = Mx[i], Mx[k]; sign = -sign; break
            else:
                return zero
        for i in range(k + 1, n):
            for l in range(k + 1, n):
                num = Mx[i][l] * Mx[k][k] - Mx[i][k] * Mx[k][l]
                Mx[i][l] = num / prev if prev != one else num   # exact division (Bareiss)
            Mx[i][k] = zero
        prev = Mx[k][k]
    return sign * Mx[n - 1][n - 1]
Dt = bareiss(H)
check("det H(c + r kappa) has degree 3 in r", Dt.degrees()[7] == 3, Dt.degrees())
Dk = [zero] * 4
for mon, cf in Dt.to_dict().items():
    mon2 = list(mon); k = mon2[7]; mon2[7] = 0
    Dk[k] = Dk[k] + ctx.from_dict({tuple(mon2): cf})
N = Q0**3 * Dk[0] - Q0**2 * Cc * Dk[1] + Q0 * Cc**2 * Dk[2] - Cc**3 * Dk[3]   # Q0^3 det H(c - (C/Q0) kappa)
N1, r1 = divmod(N, Q0)
N2, r2 = divmod(N1, -4 * one)
G = N2.sqrt() if (r1 == 0 and r2 == 0) else None
ok = G is not None and G * G == N2
check("Q0^2 det H(c - (C/Q0) kappa) = -4 G^2 for a polynomial G", ok)
check("G has 1772 terms and content 1", ok and len(G) == 1772 and G.content() == 1, len(G) if ok else None)

# ---- C. chart substitution and the coefficient certificate --------------------------
cn = ['L', 'A', 'B', 'C1', 'S', 'P', 'c']
cx = flint.fmpz_mpoly_ctx.get(cn, 'lex')
Lx, Ax, Bx, Cx1, Sx, Px, cxv = cx.gens(); o = cx.from_dict({(0,) * 7: 1})
EN = Ax * Sx * (1 + Lx) + 2 * (1 + Ax)
psn = 2 * Sx * (1 + Lx) * (1 + Bx) - Lx * Sx * (1 + Bx) + 2 * (1 + Bx) - Sx * Bx
sub = {'q': (2 * Bx * Cx1, (1 + Lx) * (1 + Bx) * (1 + Cx1)), 'b': (2 * Bx, (1 + Lx) * (1 + Bx) * (1 + Cx1)),
       'a': (2 * Lx * (1 + Ax), 1 + Lx), 'e': (2 * Lx * Sx * (1 + Ax), EN), 'j': (Sx, o), 'c': (cxv, o),
       'p': (Px * psn, (1 + Px) * 2 * Sx * (1 + Lx) * (1 + Bx))}
degs = dict(zip(names, G.degrees()))
Num = cx.from_dict({}); cache = {}
def pw(k, w, ex):
    if (k, w, ex) not in cache:
        cache[(k, w, ex)] = sub[k][w]**ex
    return cache[(k, w, ex)]
for mon, cf in G.to_dict().items():
    t = o * int(cf)
    for n_, ex in zip(names, mon):
        if n_ != 'r':
            t = t * pw(n_, 0, ex) * pw(n_, 1, degs[n_] - ex)
    Num = Num + t
Dfull = o
for n_ in names:
    if n_ != 'r':
        Dfull = Dfull * sub[n_][1]**degs[n_]
Dclear = 4 * Sx**2 * (1 + Lx)**9 * (1 + Bx)**6 * (1 + Cx1)**6 * (1 + Px)**2 * EN**2
DG, rem = divmod(Dclear * Num, Dfull)
check("D*G is a polynomial in (L,A,B,C1,S,P,c)", rem == 0)
check("D*G has 149933 monomials", len(DG) == 149933, len(DG))
mDG = -DG
groups = defaultdict(dict)
for mon, cf in mDG.to_dict().items():
    groups[tuple(mon[:3] + mon[4:])][mon[3]] = int(cf)
neg = [k for k, gr in groups.items() if any(x < 0 for x in gr.values())]
check("exactly 18 C1-groups of -DG contain a negative coefficient", len(neg) == 18, len(neg))
sq = flint.fmpz_poly([1, -2, 1]); allsq = True
for k in neg:
    gr = groups[k]; pol = flint.fmpz_poly([gr.get(i, 0) for i in range(max(gr) + 1)])
    qt, rt = divmod(pol, sq)
    allsq = allsq and rt == 0 and all(x >= 0 for x in qt.coeffs())
check("each such group = (C1-1)^2 R(C1), R with nonnegative integer coefficients", allsq)
check("-DG coefficient of B^6 S^2 is 432", mDG.to_dict().get((0, 0, 6, 0, 2, 0, 0)) == 432)
check("-DG has degree <= 6 in B (homogenization with W_B of degree 6 is valid)", mDG.degrees()[2] <= 6)
negtop = [k for k in neg if k[2] == 6]
check("B^6 slice: its negative coefficients lie in complete square groups (18 of them)", len(negtop) == 18)

# ---- B. zero faces -> paper cases ---------------------------------------------------
want = dict(q=39, p=1, j=3, e=24, b=43, a=44)   # doc's T38,T0,T3,T23,T42,T43 = historical labels
for k, nu in want.items():
    rest = [roles[x] for x in order if x != k]
    Z = frozenset(s for s in range(24) if all(L.F_VALUES[f][s] == 0 for f in rest))
    hit = any(L.find_g(L.CONTACT_SETS[nu], set(S6)) is not None for S6 in itertools.combinations(sorted(Z), 6))
    check("zero face %s=0: contact set %s contains a G-image of paper T_%d" % (k, sorted(Z), nu), hit)

print("elapsed %.1f s" % (time.time() - t0))
print("ALL EXACT CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)
sys.exit(1 if FAIL else 0)
