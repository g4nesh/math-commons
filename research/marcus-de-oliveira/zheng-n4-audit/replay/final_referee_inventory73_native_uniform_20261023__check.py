"""Own-code replay for the Zheng n=4 audit, document final_referee_inventory73_native_uniform_20261023
(primary proof of paper case T_74 = {3,4,5,7,8,12}).
Sigma-points from the definitions with A = (y,z,w,x), C = (t,0,1,inf), c = -b, minus brackets
A_i - C_j (the common factor of the point at infinity is dropped; it only changes the overall sign).
Nothing from the supplement is read, imported or executed.  Exit code 1 if an EXACT check fails.
The radial polynomial of the K face (eq. (17)) is checked in the T_72 replay
(final_referee_inventory71_uniform_20261023__check.py), which rebuilds G9 from the cubic."""
import sys, os, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp
from fractions import Fraction

FAIL = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name, flush=True)
    if not cond:
        FAIL.append(name)
ok = lambda e: sp.expand(sp.fraction(sp.cancel(sp.together(e)))[0]) == 0
t0 = time.time()
FV = L.F_VALUES
IDX = {tuple(v): j for j, v in enumerate(FV)}
def gf(gi, j):
    return IDX[tuple(L.g_act_values(L.G_ELEMENTS[gi], list(FV[j])))]

# ---------- A/B ----------
T = set(L.CONTACT_SETS[74]); gen = L.generators(T)
p3 = [0, 1, 2, 5, 6, 59, 67, 70]; r3 = [18, 16, 13, 10, 9, 2, 14, 6]
check("T_74 = {3,4,5,7,8,12}; generators = [0,1,2,5,6,59,67,70] (doc (1)), rank 8",
      sorted(T) == [3, 4, 5, 7, 8, 12] and gen == p3 and L.rank([list(L.F_VECS[j]) for j in gen]) == 8)
check("all minor except f70 = h4231,+", [L.F_IS_MINOR[j] for j in gen] == [True] * 7 + [False])
check("identity rows [18,16,13,10,9,2,14,6] (coefficients = those gaps; d = g6)",
      all(FV[p3[i]][r3[k]] == (1 if i == k else 0) for k in range(8) for i in range(8)))
check("f70 = 1 at rows 0 and 1, all generators >= 0: d <= g0, d <= g1",
      FV[70][0] >= 1 and FV[70][1] >= 1 and all(FV[j][s] >= 0 for j in gen for s in range(24)))
ext = [j for j in range(24) if j not in T and FV[70][j] == 0]
check("extra ties with f70(j) = 0 are exactly 2, 9, 10, 13, 14, 16, 18", ext == [2, 9, 10, 13, 14, 16, 18])
for j in ext:
    S = T | {j}
    cv = [k for k, Tk in L.CONTACT_SETS.items() if k != 74 and any(L.g_act_set(gi, Tk) <= S for gi in range(1152))]
    a7 = [k for k, A in L.AUX7.items() if L.find_g(A, S)]
    print(f"info tie {j}: paper cases T_k with g.T_k inside T u {{{j}}}: {cv}; Table-5 sets equal to an image: {a7}")
    if j != 16:
        check(f"tie {j} covered by some paper case T_k with k < 74", any(k < 74 for k in cv))
g = (L.perm_from_oneline("1234"), L.perm_from_oneline("1432"), 0); gi = L.G_ELEMENTS.index(g)
K = L.AUX7[16]
check("K = {0,1,2,9,10,14,17} = Table 5 set #16 = T'(of T_72 doc) u {9}; rows identity, columns (0,3,2,1) map K onto T u {16}",
      K == frozenset({0, 1, 2, 9, 10, 14, 17}) and L.g_act_set(gi, K) == T | {16})
check("(3) profile action [0,2,12,13,15,21,26] -> [0,2,5,6,59,67,70] in order",
      [gf(gi, j) for j in [0, 2, 12, 13, 15, 21, 26]] == [0, 2, 5, 6, 59, 67, 70] and L.generators(K) == [0, 2, 12, 13, 15, 21, 26])
check("tie 16 is also covered by paper case T_72 (the K argument can equally cite T_72)",
      any(L.g_act_set(x, L.CONTACT_SETS[72]) <= T | {16} for x in range(1152)))

# ---------- C: coordinates ----------
y, z, w, x, t = sp.symbols('y z w x t')
u, Z, r, s, p, h, X, Tt = sp.symbols('u Z r s p h x t_')
A = [y, z, w, x]; Cc = [t, 0, 1, None]
def Fk(k):
    sig = L.PERMS[k]; pr = 1
    for i in range(4):
        if Cc[sig[i]] is not None:
            pr *= (A[i] - Cc[sig[i]])
    return sp.expand(pr)
F = [Fk(k) for k in range(24)]
E = F[5] - F[3]
check("E = F5 - F3 = x(y-t)(z-w) up to the global sign of the dropped infinite factor", ok(E + x * (y - t) * (z - w)))
check("F4 - F5 = (w-x)(y-t) up to the same global sign", ok(F[4] - F[5] + (w - x) * (y - t)))
H = [sp.cancel((F[i] - F[3]) / E) for i in range(24)]
rr, ss, pp, hh, uu = H[7], H[8], H[4], H[12], H[10]
Dd = rr - uu; kk = (ss - uu) / Dd
check("(5) D = r - u = -y(x-1)/(x(y-t)), k = (s-u)/D = (t-w)(x-z)/((x-1)(z-w))",
      ok(Dd + y * (x - 1) / (x * (y - t))) and ok(kk - (t - w) * (x - z) / ((x - 1) * (z - w))))
check("(6) y = Dxt/(x(D+1)-1), t = w + k(x-1)(z-w)/(x-z), w = z + (x-z)(z-1)/(s(1-x)+u(z-1))",
      ok(y - Dd * x * t / (x * (Dd + 1) - 1)) and ok(t - (w + kk * (x - 1) * (z - w) / (x - z)))
      and ok(w - (z + (x - z) * (z - 1) / (ss * (1 - x) + uu * (z - 1)))))
YY = 1 / (x - 1); ZZ = (z - 1) / (x - 1); LL = uu * ZZ - ss; JJ = (uu - pp) * ZZ - ss
check("H4 = (YL+Z)/(Z(Y+1)) and (7) Y = Z(p-1)/J", ok(pp - (YY * LL + ZZ) / (ZZ * (YY + 1))) and ok(YY - ZZ * (pp - 1) / JJ))
check("(9) H6 = -L + (r+1)Z + rZ(1-Z)/L, H1 = (p-1)L/(L-Z) (in source coordinates)",
      ok(H[6] - (-LL + (rr + 1) * ZZ + rr * ZZ * (1 - ZZ) / LL)) and ok(H[1] - (pp - 1) * LL / (LL - ZZ)))
Lq = u * Z - s; Dq = r - u; Jq = (u - p) * Z - s
VN = Z * ((u - 1) * Z - s + 1); TN = Dq * VN - (s - u) * Z; YN = Z * (p - 1)
Asp = [((YN + Jq) * TN - YN * Lq, Lq * (Dq * YN + (Dq + 1) * Jq)), (Z, 1), (VN, Lq), (1, 1)]
Csp = [(TN, Dq * Lq), (-YN, Jq), (0, 1), (1, 0)]
br = lambda P, Q: P[0] * Q[1] - P[1] * Q[0]
G = [sp.expand(sp.Mul(*[br(Asp[i], Csp[L.PERMS[k][i]]) for i in range(4)])) for k in range(24)]
Eg = G[5] - G[3]
Hg = [sp.cancel((G[i] - G[3]) / Eg) for i in range(24)]
check("spinor inverse: H3 = 0, H5 = 1, H7 = r, H8 = s, H4 = p, H10 = u identically",
      Hg[3] == 0 and Hg[5] == 1 and ok(Hg[7] - r) and ok(Hg[8] - s) and ok(Hg[4] - p) and ok(Hg[10] - u))
check("(9) in spinor coordinates: H6 and H1 formulas", ok(Hg[6] - (-Lq + (r + 1) * Z + r * Z * (1 - Z) / Lq)) and ok(Hg[1] - (p - 1) * Lq / (Lq - Z)))
num, den = sp.fraction(sp.cancel(Hg[12] - h))
C = sp.expand(-num)    # global sign chosen so that [u^5]C = +Z^4 as in the document
Pc = sp.Poly(C, u, Z)
check("(8) C(u,Z): degrees (5,4), 144 terms, [u^5]C = Z^4",
      Pc.degree(u) == 5 and Pc.degree(Z) == 4 and len(sp.Poly(C, u, Z, r, s, p, h).terms()) == 144 and ok(C.coeff(u, 5) - Z**4))
check("(11) [Z^4]C = -(p-u)(u-1)^2(r-u)^2", ok(C.coeff(Z, 4) + (p - u) * (u - 1)**2 * (r - u)**2))
# left-inverse property at random rational sources
random.seed(7); good = True
for _ in range(4):
    vals = {y: Fraction(random.randint(-90, 90), random.randint(1, 9)) for y in []}
    sv = {y: sp.Rational(random.randint(-50, 50), random.randint(1, 7)), z: sp.Rational(random.randint(-50, 50), random.randint(1, 7)),
          w: sp.Rational(random.randint(-50, 50), random.randint(1, 7)), x: sp.Rational(random.randint(-50, 50), random.randint(1, 7)),
          t: sp.Rational(random.randint(-50, 50), random.randint(1, 7))}
    Hn = [Hi.subs(sv) for Hi in H]
    pv = {u: Hn[10], Z: ZZ.subs(sv), r: Hn[7], s: Hn[8], p: Hn[4], h: Hn[12]}
    Hb = [Hi.subs(pv) for Hi in Hg]
    good &= all(sp.simplify(Hb[i] - Hn[i]) == 0 for i in range(24)) and C.subs(pv) == 0
check("spinor inverse is a left inverse: at 4 random rational sources all 24 H_i are reproduced and C(u,Z) = 0", good)
# (12) u = inf, Z = 0
Ct = sp.Poly(sp.expand(Tt**5 * C.subs({u: 1 / Tt, Z: X})), X, Tt)
md = min(sum(m_) for m_ in Ct.monoms())
low = sum(c_ * X**m_[0] * Tt**m_[1] for m_, c_ in Ct.terms() if sum(m_) == md)
check("(12) at (u,Z) = (inf,0): lowest form has degree 4 and equals (st-ht-x)(st-pt-x+t)(st-x)^2, value x^4 at t=0",
      md == 4 and ok(low - (s * Tt - h * Tt - X) * (s * Tt - p * Tt - X + Tt) * (s * Tt - X)**2))
check("no point of the closure at (u,Z) = (inf,inf): [u^5 Z^4]C = 1 != 0", Pc.coeff_monomial(u**5 * Z**4) == 1)
# (13)-(15) u = 1, Z = inf
Cr = sp.Poly(sp.expand(Tt**4 * C.subs({u: 1 + X, Z: 1 / Tt})), X, Tt)
wd = min(m_[0] + 2 * m_[1] for m_ in Cr.monoms())
lowr = sum(c_ * X**m_[0] * Tt**m_[1] for m_, c_ in Cr.terms() if m_[0] + 2 * m_[1] == wd)
check("(13) weighted (1,2) lowest part = (r-1)(p-1)^2 t - (p-1)(r-1)^2 x^2; dt-derivative at 0 = (r-1)(p-1)^2",
      ok(lowr - ((r - 1) * (p - 1)**2 * Tt - (p - 1) * (r - 1)**2 * X**2))
      and ok(sp.diff(Cr.as_expr(), Tt).subs({X: 0, Tt: 0}) - (r - 1) * (p - 1)**2) and Cr.as_expr().subs({X: 0, Tt: 0}) == 0)
c2, c3, c4 = sp.symbols('c2 c3 c4'); ex = sp.expand(Cr.as_expr().subs(Tt, c2 * X**2 + c3 * X**3 + c4 * X**4)); sol = {}
for k_, cc in [(2, c2), (3, c3), (4, c4)]:
    sol[cc] = sp.solve(sp.expand(ex.coeff(X, k_)).subs(sol), cc)[0]
check("(14) branch t = (r-1)/(p-1) x^2 + O(x^3)", ok(sol[c2] - (r - 1) / (p - 1)))
brn = {u: 1 + X, Z: 1 / (sol[c2] * X**2 + sol[c3] * X**3 + sol[c4] * X**4)}
def lim(e):
    nu, de = sp.fraction(sp.cancel(sp.together(e.subs(brn))))
    nu = sp.Poly(sp.expand(nu), X); de = sp.Poly(sp.expand(de), X)
    on = min(m_[0] for m_ in nu.monoms()); od = min(m_[0] for m_ in de.monoms())
    return on, od, sp.cancel(nu.as_expr().coeff(X, on) / de.as_expr().coeff(X, od))
res = {i: lim(Hg[i]) for i in (13, 16, 10, 6)}
print("info orders (num, den, ratio of leading coefficients) along the ramified branch:", {i: (a_, b_) for i, (a_, b_, _) in res.items()})
check("(15) H13 -> 1, H16 -> 1, H10 -> 1 on the ramified branch", all(res[i][0] == res[i][1] and ok(res[i][2] - 1) for i in (13, 16, 10)))
check("H6 has (at most) a simple pole there", res[6][1] - res[6][0] == 1)
check("Z -> inf growth of H6 at constant u: coefficient -u + (r+1) - r/u = (r-u)(u-1)/u",
      ok(sp.limit((-Lq + (r + 1) * Z + r * Z * (1 - Z) / Lq) / Z, Z, sp.oo) - (r - u) * (u - 1) / u))
print("info exact part time", round(time.time() - t0, 1))

# ---------- NUMERIC (brief) ----------
import numpy as np
found = 0; minphi = np.inf
for seed in range(6):
    pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=30, margin=1e-2)
    if pt is None:
        continue
    found += 1
    c_vec = np.linalg.lstsq(np.array(L.P_numpy(), float), np.array(pt['gaps'], float), rcond=None)[0]
    mv, _ = L.min_phi_over_SO4(c_vec, restarts=8, seed=seed)
    minphi = min(minphi, mv)
print(f"NUMERIC T_74 stratum points found {found}/6 (seeds 0-5, margin 1e-2): min Phi over SO(4) (8 restarts each) = {minphi:.4g}")
print("time", round(time.time() - t0, 1))
print("SUMMARY", "all EXACT checks passed" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
