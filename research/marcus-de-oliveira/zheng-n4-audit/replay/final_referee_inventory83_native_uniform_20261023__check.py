"""Auditor replay for final_referee_inventory83_native_uniform_20261023 (paper case T_84).

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

T = frozenset({2, 5, 9, 11, 12, 13})
MYSTAGE = 9
T0 = frozenset({3, 4, 8, 10, 12, 13})
check("working contact set equals paper T_84", T == L.CONTACT_SETS[84])
g = L.find_g(T0, T)
check("stated set {3,4,8,10,12,13} maps onto T by a group element", g is not None, g)
gens = L.generators(T)
check("source S = (0,1,7,12,39,40,41,45,67,68,72)", gens == [0, 1, 7, 12, 39, 40, 41, 45, 67, 68, 72], gens)
check("non-minor generators are exactly f_40, f_68", [j for j in gens if not L.F_IS_MINOR[j]] == [40, 68])
check("rank of generators 8", L.rank([L.F_VECS[j] for j in gens]) == 8)
check("eq.(2): f_40 and f_68 are >= 1 at sigma = 10, 14, 18 (so Theta <= g_10, g_14, g_18)",
      all(L.F_VALUES[j][s] >= 1 for j in (40, 68) for s in (10, 14, 18)))
extra = [j for j in range(24) if j not in T and (L.F_VALUES[40][j] == 0 or L.F_VALUES[68][j] == 0)]
print("INFO extra ties where f_40 or f_68 vanishes (list NOT given in the PDF):", extra)
for j in extra:
    cv = covering(T | {j}, MYSTAGE)
    check("boundary +%d covered by an earlier-stage statement (auditor-supplied)" % j, len(cv) > 0,
          [(c[0], "stage %d" % c[1]) for c in cv][:4])
other = [j for j in range(24) if j not in T and j not in extra]
check("at every other extra tie both f_40 and f_68 are positive (so Theta = 0, minor sum)",
      all(L.F_VALUES[40][j] > 0 and L.F_VALUES[68][j] > 0 for j in other), other)

PERMS = list(itertools.permutations(range(4)))
x, z, w, v, t, r, s, u, h, tau, sig = sp.symbols('x z w v t r s u h tau sigma')
A = [0, 1, z, w]; C = [None, v, x, t]
def Fp(i):
    p = PERMS[i]; e = sp.Integer(1)
    for a in range(4): e *= (-1 if p[a] == 0 else (A[a] - C[p[a]]))   # [A, inf] = -1
    return sp.expand(e)
Fs = [Fp(i) for i in range(24)]
E = sp.expand(Fs[12] - Fs[2])
check("E = F12 - F2 = (z-v)(w-t)", sp.expand(E - (z - v)*(w - t)) == 0)
check("height F12/E = x", sp.cancel(Fs[12]/E - x) == 0)
H = {i: sp.cancel((Fs[i] - Fs[12])/E) for i in range(24)}
check("H2 = -1 and H12 = 0 identically", H[2] == -1 and H[12] == 0)
check("F5/F13 and F9/F13 elementary forms",
      sp.cancel(Fs[5]/Fs[13] + (1-t)*(z-x)/(x*(z-t))) == 0 and sp.cancel(Fs[9]/Fs[13] - v*(1-x)/(x*(w-v))) == 0)
V = r + x*(1 + u - r); Td = z*(x + r) + x*(s - r); Tn = z*(x**2 + x*(s + 1) + r) - x*(x + r)
vs = x*(x + u)*w/V; ts = Tn/Td; sub = {v: vs, t: ts}
check("eq.(4): v, t solve F5(x+r) = F13(x+s) and F9(x+r) = F13(x+u)",
      sp.cancel((Fs[5]*(x+r) - Fs[13]*(x+s)).subs(sub)) == 0 and sp.cancel((Fs[9]*(x+r) - Fs[13]*(x+u)).subs(sub)) == 0)
check("Td - Tn = x(x+s)(1-z)", sp.expand(Td - Tn - x*(x + s)*(1 - z)) == 0)
n13 = sp.factor_list(sp.numer(sp.cancel((H[13] - r).subs(sub))))
bigs = [f for f, m in n13[1] if len(sp.Poly(f, z, w, x, r, s, u).terms()) >= 10]
lins = [f for f, m in n13[1] if f not in bigs]
check("H13 = r becomes (x+r) C = 0 up to a nonzero constant/monomial factor",
      len(bigs) == 1 and (x + r) in lins, lins)
Cc = sp.expand(bigs[0])
lc = sp.Poly(Cc, z, w).coeff_monomial(z*w**2); Cc = sp.expand(Cc*sp.cancel(-x*(x+u)*(x+r)/lc))
check("C has 53 terms, bidegree (2,2) in (z,w)", len(sp.Poly(Cc, z, w, x, r, s, u).terms()) == 53 and
      sp.Poly(Cc, z, w).degree(z) == 2 and sp.Poly(Cc, z, w).degree(w) == 2)
check("normalization of C is a constant (independent of x,r,s,u)", sp.cancel(-x*(x+u)*(x+r)/lc).is_number)
print("INFO eq.(7) [z w^2]C = -x(x+u)(x+r) is used as the normalization of C; it is a genuine identity because the normalizing factor checked above is a pure number")
check("eq.(13): [w^2] C = -x(x+u) Td", sp.expand(sp.Poly(Cc, w).coeff_monomial(w**2) + x*(x+u)*Td) == 0)
check("eq.(11): C(x,w) = x^2(x+s)(x-w)((x+u)w - V)", sp.expand(Cc.subs(z, x) - x**2*(x+s)*(x-w)*((x+u)*w - V)) == 0)
def modC(e):
    n = sp.expand(sp.numer(sp.cancel(sp.together(e))))
    return n == 0 or sp.prem(sp.Poly(n, w), sp.Poly(Cc, w)).is_zero
K = x*(x + s)*(x + u)/((x - 1)*(x + r)) - x
check("eq.(6): H11 = K(x) on C = 0", modC(H[11].subs(sub) - K))
J = (x + K)*(w - x)/(z - x) - x
check("eq.(8): H10 = J on C = 0", modC(H[10].subs(sub) - J))
check("x + K(x) = x(x+s)(x+u)/((x-1)(x+r))", sp.cancel(x + K - x*(x+s)*(x+u)/((x-1)*(x+r))) == 0)
check("eq.(9): H14 = x(1-z)/(z-v) = x(1-z)V/(zV - x(x+u)w) (no curve needed)",
      sp.cancel(H[14] - x*(1-z)/(z-v)) == 0 and sp.cancel(H[14].subs(sub) - x*(1-z)*V/(z*V - x*(x+u)*w)) == 0)
check("eq.(10): H18 = t(w-x)/(w-t) - x = w(Tn - x Td)/(Td w - Tn) (no curve needed)",
      sp.cancel(H[18] - (t*(w-x)/(w-t) - x)) == 0 and sp.cancel(H[18].subs(sub) - w*(Tn - x*Td)/(Td*w - Tn)) == 0)
w0 = V/(x + u)
check("eq.(12): w0 - x = (1-x)(x+r)/(x+u)", sp.cancel(w0 - x - (1-x)*(x+r)/(x+u)) == 0)
check("H14 denominator at (x,x) is x(1-x)(x+r); value V/(x+r)",
      sp.expand((z*V - x*(x+u)*w).subs({z: x, w: x}) - x*(1-x)*(x+r)) == 0 and
      sp.cancel((x*(1-z)*V/(z*V - x*(x+u)*w)).subs({z: x, w: x}) - V/(x+r)) == 0)
check("at (x,w0): v = t = x, H18 = 0, (x+u)(Td w0 - Tn) = x(x+s)(1-x)(x+r)",
      sp.cancel(vs.subs({z: x, w: w0}) - x) == 0 and sp.cancel(ts.subs({z: x, w: w0}) - x) == 0 and
      sp.cancel((w*(Tn - x*Td)/(Td*w - Tn)).subs({z: x, w: w0})) == 0 and
      sp.cancel(((x+u)*(Td*w0 - Tn)).subs(z, x) - x*(x+s)*(1-x)*(x+r)) == 0)
check("w = infinity only at Td = 0, z = x(r-s)/(x+r) != x", sp.solve(Td, z) == [x*(r - s)/(r + x)])
Ct = sp.expand(tau**2*sig**2*Cc.subs({z: 1/tau, w: 1/sig}))
check("corner: Ctilde(0,0) = 0, d_tau Ctilde(0,0) = -x(x+u)(x+r)",
      Ct.subs({tau: 0, sig: 0}) == 0 and sp.expand(sp.diff(Ct, tau).subs({tau: 0, sig: 0}) + x*(x+u)*(x+r)) == 0)
check("x F0 = (x - z) F14 (raw products)", sp.expand(x*Fs[0] - (x - z)*Fs[14]) == 0)

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
    print("INFO NUMERIC T_84 stratum: %d points (seeds 500.., margin 1e-3), 10 BFGS restarts each; "
          "worst min Phi / max gap = %.3e" % (found, worst))
print("SUMMARY: %d EXACT failures" % len(FAIL), FAIL)
sys.exit(1 if FAIL else 0)
