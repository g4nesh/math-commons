"""Auditor replay for blind_inventory84_rank_drop_uniform_20261023 (paper case T_85).

Own code only; nothing from the supplement is imported or executed, and no supplement data is
read. Run from any directory:  python -I -B <this file>   (optionally --no-numeric).
Prints one line per check; exits 1 if any EXACT check fails.
"""
import itertools, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + "EXACT " + name + ((" | " + str(info)) if info else ""))
    if not ok: FAIL.append(name)

# ---------- paper proof order (Appendix G, Table 6) ----------
ST6 = {1:[1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
       2:[10,26,46,90],3:[9,11,12,47],4:[8,13,17,49],5:[14,16,24,33,38,70],6:[52,57,58,65,68,71,76,78],
       7:[45,51,66,69,75,81],8:[48,50,53,59,62,67,73,79,86],9:[30,41,54,55,60,61,63,74,77,80,83,84],
       10:[56,64,82,85,87],11:[32],12:[89]}
STAGE6 = {n: s for s, l in ST6.items() for n in l}
ST7 = {1:[1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24],2:[3,4,23,26],4:[6,18],5:[20,27],7:[19,25]}
STAGE7 = {n: s for s, l in ST7.items() for n in l}

def covering(Tj, before_stage):
    """Earlier-stage paper statements (6-sets T_n or 7-sets T^(7)_n) some g in G maps INTO Tj."""
    out = []
    for name, S, st in ([("T%d" % n, S, STAGE6[n]) for n, S in L.CONTACT_SETS.items()] +
                        [("T7_%d" % n, S, STAGE7[n]) for n, S in L.AUX7.items()]):
        if st >= before_stage: continue
        for g, tm in zip(L.G_ELEMENTS, L._G_MAPS):
            if frozenset(tm[j] for j in S) <= Tj:
                out.append((name, st, g)); break
    return out

T = frozenset({0, 5, 7, 9, 14, 18})
MYSTAGE = 10
from fractions import Fraction as Fr
check("contact set equals paper T_85", T == L.CONTACT_SETS[85])
check("six contact rows have rank 6; kernel dimension 8",
      L.rank([L.P[j] for j in sorted(T)]) == 6 and len(L.nullspace([L.P[j] for j in sorted(T)])) == 8)
basis = [1, 2, 5, 12, 29, 33, 45, 35]; names = list("abcefgkd")
check("eq.(2): basis vectors vanish on T and are independent",
      all(L.F_VALUES[j][i] == 0 for j in basis for i in T) and L.rank([L.F_VECS[j] for j in basis]) == 8)
check("first seven basis vectors minor, f_35 non-minor",
      all(L.F_IS_MINOR[j] for j in basis[:7]) and not L.F_IS_MINOR[35])
excl = [17, 15, 11, 6, 4, 12, 1]
check("exclusive rows (17,15,11,6,4,12,1): coefficient matrix I7 with zero d-column",
      all([L.F_VALUES[basis[p]][sg] for p in range(8)] == [1 if p == q else 0 for p in range(8)] for q, sg in enumerate(excl)))
gens = L.generators(T)
check("full generator list has 11 elements (1,2,5,11,12,28,29,33,35,39,45)", gens == [1, 2, 5, 11, 12, 28, 29, 33, 35, 39, 45], gens)
def vadd(*vs):
    out = [Fr(0)]*14
    for c, vv in vs: out = [a + c*b for a, b in zip(out, vv)]
    return out
F = L.F_VECS
check("eq.(3a): f35 + f11 = f1 + f5 + f12 + f33", vadd((1,F[35]),(1,F[11]),(-1,F[1]),(-1,F[5]),(-1,F[12]),(-1,F[33])) == [0]*14)
check("eq.(3b): f35 + f39 = f1 + f2 + f29 + f45", vadd((1,F[35]),(1,F[39]),(-1,F[1]),(-1,F[2]),(-1,F[29]),(-1,F[45])) == [0]*14)
check("f11, f39 are minor vectors", L.F_IS_MINOR[11] and L.F_IS_MINOR[39])
a_, b_, c_, e_, f_, g_, k_, dl = sp.symbols('a b c e f g k delta')
coef = [a_, b_, c_, e_, f_, g_, k_, -dl]
rowsd = {sp.expand(sum(coef[q]*L.F_VALUES[basis[q]][sg] for q in range(8))) for sg in range(24) if sg not in T}
need = {a_ - dl, c_ + b_ - dl, c_ + f_ - dl, c_ + k_ - dl, e_ + b_ - dl, e_ + f_ - dl, g_ + b_ - dl, g_ + f_ - dl, g_ + k_ - dl}
check("negative branch: gap rows contain a >= delta and the eight cross inequalities", need <= rowsd)
check("negative branch: every other gap row is implied (single coefficients or a+e+k-2delta)",
      rowsd - need <= {a_, b_, c_, e_, f_, g_, k_, a_ + e_ + k_ - 2*dl}, rowsd - need)
th = sp.Symbol('theta')
lhs = [a_, b_, c_, e_, f_, g_, k_, -dl]
vec_old = vadd(*[(lhs[q], F[basis[q]]) for q in range(8)])
newc = [a_ - dl, b_ - dl + th, c_ - th, e_ - th, f_ - dl + th, g_ - th, k_ - dl + th]
vec_new = vadd(*([(newc[q], F[basis[q]]) for q in range(7)] + [(th, F[11]), (dl - th, F[39])]))
check("minor rewriting: old vector = new seven coefficients + theta f11 + (delta-theta) f39",
      all(sp.expand(p - q) == 0 for p, q in zip(vec_old, vec_new)))
gT = ((0, 1, 2, 3), (0, 1, 2, 3), 1)
check("transpose sigma -> sigma^-1 preserves T", L.g_act_set(gT, T) == T)
syms = sp.symbols('a b c e f g k d')
Bm = [F[j] for j in basis]
cvec = [sum(syms[q]*Bm[q][i] for q in range(8)) for i in range(14)]
vals = [sum(L.P[sg][i]*cvec[i] for i in range(14)) for sg in range(24)]
newvals = L.g_act_values(gT, vals)
cp = sp.symbols('A0:8')
sol = sp.solve([sp.expand(newvals[sg] - sum(cp[q]*L.F_VALUES[basis[q]][sg] for q in range(8))) for sg in range(24)], cp, dict=True)
check("transpose: new d' = -d - e - k (all 24 rows consistent)", len(sol) == 1 and sp.expand(sol[0][cp[7]] + syms[7] + syms[3] + syms[6]) == 0)
check("transpose: new minor coefficients are the transformed gaps at the exclusive rows",
      all(sp.expand(sol[0][cp[q]] - newvals[excl[q]]) == 0 for q in range(7)))
fz = [j for j in range(24) if j not in T and L.F_VALUES[35][j] == 0]
fp = [j for j in range(24) if j not in T and L.F_VALUES[35][j] > 0]
check("boundary table: f35 vanishes off T exactly at 1,4,6,11,12,15,17", fz == [1, 4, 6, 11, 12, 15, 17], fz)
check("boundary table: f35 > 0 at 2,3,8,10,13,16,19,20,21,22,23 (d=0 rows)", fp == [2, 3, 8, 10, 13, 16, 19, 20, 21, 22, 23], fp)
for js, n in [((1, 6), 40), ((4, 11, 12, 15), 55), ((17,), 81)]:
    for j in js:
        hit = [gg for gg, tm in zip(L.G_ELEMENTS, L._G_MAPS) if frozenset(tm[q] for q in L.CONTACT_SETS[n]) <= T | {j}]
        check("boundary +%d: some g maps paper T_%d (stage %d < 10) into T u {%d}" % (j, n, STAGE6[n], j),
              len(hit) > 0 and STAGE6[n] < MYSTAGE, hit[0] if hit else None)
for row in L.TABLE3:
    j, n, rho, rho2, eps, img = row
    gg = (L.perm_from_oneline(rho), L.perm_from_oneline(rho2), eps)
    im = L.g_act_set(gg, L.CONTACT_SETS[n])
    check("paper Table 3 row j=%d: image of T_%d equals listed set and lies in T u {j}" % (j, n),
          im == frozenset(img) and im <= T | {j})
p, r, s, u, h, zz = sp.symbols('p r s u h z')
Lm = [[Bm[q][i] for q in range(8)] for i in range(14)]
v0 = [-r, 2*r - s, -h + p + 2*r - s + u, -r + s - u, 0, 0, 0, -p - r + s, 0, 0, 0, 0, 0, 0]
m = [sum(Lm[i][q]*syms[q] for q in range(8)) + zz*L.KAPPA[i] + v0[i] for i in range(14)]
Pm = [sp.expand(sum(L.P[sg][i]*m[i] for i in range(14))) for sg in range(24)]
check("eq.(4): (P5-P0)m=p, (P7-P0)m=r, (P9-P0)m=s, (P14-P0)m=u, (P18-P0)m=h",
      all(sp.expand(Pm[j] - Pm[0] - dv) == 0 for j, dv in {5: p, 7: r, 9: s, 14: u, 18: h}.items()))
allv = list(syms) + [zz, p, r, s, u, h]
Jm = sp.Matrix([[sp.diff(m[i], vv) for vv in allv] for i in range(14)])
check("eq.(4): coordinate change is invertible (det = -1)", Jm.det() == -1)
check("eq.(5): d(m) = m12 - m4 - m8", sp.expand(m[12] - m[4] - m[8] - syms[7]) == 0)
check("gaps: (P_sigma - P_0)m = sum coeff f_j(sigma) + real-linear(p,r,s,u,h) for all sigma",
      all(sp.expand(Pm[sg] - Pm[0] - sum(syms[q]*L.F_VALUES[basis[q]][sg] for q in range(8))).free_symbols <= {p, r, s, u, h} for sg in range(24)))
msy = sp.symbols('m0:14')
Cub = sp.expand(L.cubic_C(msy))
Gs = [sp.expand(sp.diff(Cub, mi).subs(dict(zip(msy, m)), simultaneous=True)) for mi in msy]
a, b, c, e, f, g, k, d = syms
check("eq.(6): -G2+G3-G4+G7-G8-G9+G10+G11-G12+G13 = -r(a+d) + h(s-r) - pu",
      sp.expand(-Gs[2]+Gs[3]-Gs[4]+Gs[7]-Gs[8]-Gs[9]+Gs[10]+Gs[11]-Gs[12]+Gs[13] - (-r*(a+d) + h*(s-r) - p*u)) == 0)
check("eq.(8): G13 - G2 = -(d+f)(d+g) + r(d+f) + (s-r)g", sp.expand(Gs[13] - Gs[2] - (-(d+f)*(d+g) + r*(d+f) + (s-r)*g)) == 0)
G1 = sp.Poly(Gs[1], zz)
check("eq.(9): G1 affine in z with z-coefficient -(c+d+g+k)", G1.degree() == 1 and sp.expand(G1.coeff_monomial(zz) + c + d + g + k) == 0)
check("rank [L8 | kappa] = 9", L.rank([[Fr(Lm[i][q]) for q in range(8)] + [L.KAPPA[i]] for i in range(14)]) == 9)
check("source differences (rows minus row 0) have rank 8",
      L.rank([[Fr(L.F_VALUES[basis[q]][sg] - L.F_VALUES[basis[q]][0]) for q in range(8)] for sg in range(24)]) == 8)
# Lemma 2 charts: on {r=0} the relation hs = pu recovers the fifth parameter
for fixed, rec, val in [(p, u, h*s/p), (u, p, h*s/u), (s, h, p*u/s), (h, s, p*u/h)]:
    check("Lemma 2 chart (nonzero %s): r=0 and hs=pu give %s = %s" % (fixed, rec, val),
          sp.simplify((h*s - p*u).subs(rec, val)) == 0)
# ---------- D. numeric falsification ----------
if "--no-numeric" not in sys.argv:
    import numpy as np
    Pn = L.P_numpy(); found = 0; worst = np.inf; t0 = time.time()
    for sd in range(500, 700):
        if found >= 6 or time.time() - t0 > 150: break
        pt = L.random_point_on_contact_stratum(T, seed=sd, tries=8, margin=1e-3)
        if pt is None: continue
        found += 1
        m = L.matching_vector_np(pt['a'], pt['b'], pt['alpha'])
        c = m.real.copy(); c[13] -= (Pn @ m)[pt['pi_star']].real
        vals = Pn @ c; val, _ = L.min_phi_over_SO4(c, restarts=10, seed=sd)
        worst = min(worst, val/np.max(np.abs(vals)))
        mi = 1j*m
        print('INFO NUMERIC seed %d: Im w = %.3e, |gamma7| = %.2e (expected 0 when Im w > 0)' % (sd, (mi[12]-mi[4]-mi[8]).imag, abs((Pn[7]-Pn[0]) @ mi)))
    print("INFO NUMERIC T_85 stratum: %d points (seeds 500.., margin 1e-3), 10 BFGS restarts each; "
          "worst min Phi / max gap = %.3e" % (found, worst))
print("SUMMARY: %d EXACT failures" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
