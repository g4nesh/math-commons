"""Own-code replay for the Zheng n=4 audit, document independent_inventory74_native_harmonic_uniform_20261023\n(primary proof of paper case T_75). Exact checks by default; floating-point falsification with --numeric.\nNothing from the supplement is read, imported or executed."""
import sys, os, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
NUMERIC = "--numeric" in sys.argv
FAIL = []
def chk(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(extra)) if extra else ""), flush=True)
    if not ok: FAIL.append(name)
# Table 6 of the paper (order of proof), transcribed from the paper.
_ST6 = {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
        2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
        6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81], 8: [48,50,53,59,62,67,73,79,86],
        9: [30,41,54,55,60,61,63,74,77,80,83,84], 10: [56,64,82,85,87], 11: [32], 12: [89]}
_ST7 = {1: [1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24], 2: [3,4,23,26], 4: [6,18], 5: [20,27], 7: [19,25]}
STAGE6 = {k: s for s, ks in _ST6.items() for k in ks}
STAGE7 = {k: s for s, ks in _ST7.items() for k in ks}
REP6 = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
REP7 = {L.orbit_rep(T): k for k, T in L.AUX7.items()}
def closers(S, max_stage):
    """Paper cases (stage <= max_stage) or trivial (all-minor) orbits among the 6-subsets of S,
    and Table 5 sets equal to S up to G."""
    out = set()
    for sub in itertools.combinations(sorted(S), 6):
        if all(L.F_IS_MINOR[j] for j in L.generators(sub)):
            out.add(("trivial", 0)); continue
        k = REP6[L.orbit_rep(sub)]
        if STAGE6[k] <= max_stage: out.add(("T_%d" % k, STAGE6[k]))
    if len(S) == 7 and L.orbit_rep(S) in REP7:
        k = REP7[L.orbit_rep(S)]
        if STAGE7[k] <= max_stage: out.add(("T7_%d" % k, STAGE7[k]))
    return sorted(out)
def exclusive(gens, j):
    return [s for s in range(24) if L.F_VALUES[j][s] != 0 and all(L.F_VALUES[i][s] == 0 for i in gens if i != j)]

import sympy as sp
PERMS = L.PERMS
print("== Part A: contact set, complete source, seventh ties (paper definitions, library)")
T = frozenset({3, 5, 7, 8, 10, 12})
chk("document tuple equals paper T_75 literally", T == L.CONTACT_SETS[75])
g = L.generators(T)
chk("complete source = [0,1,2,6,39,59,67,70]", g == [0, 1, 2, 6, 39, 59, 67, 70], g)
chk("only f70 non-minor", [j for j in g if not L.F_IS_MINOR[j]] == [70])
ex = [exclusive(g, j) for j in g]
chk("exclusive gaps (value 1) g18,g16,g13,g9,g4,g2,g14,g6", ex == [[18], [16], [13], [9], [4], [2], [14], [6]] and all(L.F_VALUES[j][e[0]] == 1 for j, e in zip(g, ex)), ex)
chk("(1) f70(pi_0) >= 1 and f70(pi_1) >= 1, so g0 >= d and g1 >= d", L.F_VALUES[70][0] >= 1 and L.F_VALUES[70][1] >= 1)
extra = sorted(L.F_ZERO[70] - T)
chk("(2) zeros of f70 outside T = 2,4,9,13,14,16,18", extra == [2, 4, 9, 13, 14, 16, 18], extra)
for jx in extra:
    cl = [c for c in closers(T | {jx}, 6)]
    chk("seventh tie %d: T+{%d} contains a closed case of stage <= 6 (T_75 is stage 7)" % (jx, jx), len(cl) > 0, cl)
print("== Part B: algebra of the charts, curve, atlas (sympy, exact)")
def br(P, Q):
    return P[0]*Q[1] - P[1]*Q[0]
def Fs(A, B):
    return [sp.prod([br(A[i], B[s[i]]) for i in range(4)]) for s in PERMS]
def fin(v): return (v, sp.Integer(1))
INF = (sp.Integer(1), sp.Integer(0))
def z0(e): return sp.cancel(sp.together(e)) == 0

y, z, w, x, t = sp.symbols('y z w x t')
D, k, u, h, Y, Z = sp.symbols('D k u h Y Z')
s_ = D*k + u; r_ = D + u
# ---------- original chart
A = [fin(y), fin(z), fin(w), fin(x)]; B = [fin(t), fin(sp.Integer(0)), fin(sp.Integer(1)), INF]
F = Fs(A, B)
E = sp.expand(F[5] - F[3])
chk("orig E = x(y-t)(z-w)", sp.expand(E - x*(y-t)*(z-w)) == 0)
H = [sp.cancel((F[i] - F[3]) / E) for i in range(24)]
Dex = sp.cancel(H[7] - H[10]); kex = sp.cancel((H[8] - H[10]) / Dex)
chk("D = -y(x-1)/(x(y-t))", z0(Dex + y*(x-1)/(x*(y-t))))
chk("k = (t-w)(x-z)/((x-1)(z-w))", z0(kex - (t-w)*(x-z)/((x-1)*(z-w))))
ysol = D*x*t/(x*(D+1)-1)
chk("y-formula solves D-equation", z0(Dex.subs(y, ysol) - D))
tsol = w + k*(x-1)*(z-w)/(x-z)
kk = sp.cancel(kex.subs(y, ysol))  # k-expression does not involve y
chk("k-expression independent of y", sp.cancel(kex - kk) == 0 or not kex.has(y), kex.has(y))
chk("t-formula solves k-equation", z0(kex.subs(t, tsol) - k))
Lold = s_*(1-x) + u*(z-1)
raw = sp.cancel((H[10] - u).subs(y, ysol).subs(t, tsol))
rawE = sp.cancel(((H[10] - u) * E).subs(y, ysol).subs(t, tsol))
claim = (x*(x-1)*t/((x-z)*(x*(D+1)-1)) * ((w-z)*Lold - (x-z)*(z-1))).subs(t, tsol)
chk("raw contact equation E(H10-u) = x(x-1)t/((x-z)(x(D+1)-1)) [(w-z)L-(x-z)(z-1)] after y,t substitution", z0(rawE - claim))
wsol = z + (x-z)*(z-1)/Lold
chk("w-formula solves H10=u", z0(raw.subs(w, wsol)))

# ---------- new chart
Lc = u*(Z-1) - D*k
VN = Z*((u-1)*(Z-1) - D*k)
TN = VN - k*Z
K = Y*Lc + TN
YC = D*Y + D + 1
A2 = [(D*(Y+1)*TN - Y*Lc, Lc*YC), (Z, sp.Integer(1)), (VN, Lc), (sp.Integer(1), sp.Integer(1))]
B2 = [(TN, Lc), (-Y, sp.Integer(1)), (sp.Integer(0), sp.Integer(1)), INF]
F2 = [sp.expand(f) for f in Fs(A2, B2)]
E2 = sp.expand(F2[5] - F2[3])
chk("(5) E = -Z(Z-1)L(Y+1)K", sp.expand(E2 + Z*(Z-1)*Lc*(Y+1)*K) == 0, sp.factor(E2))
H2 = [sp.cancel((F2[i] - F2[3]) / E2) for i in range(24)]
chk("(5) H3=0,H5=1,H7=D+u,H8=Dk+u,H10=u", [H2[3], H2[5]] == [0, 1] and z0(H2[7]-(D+u)) and z0(H2[8]-(D*k+u)) and z0(H2[10]-u))
chk("height F3/E = -L/(Z-1)", z0(F2[3]/E2 + Lc/(Z-1)))
AZ = (u-1)*Z - s_
N6 = D*(u-1)*Z**2 + (-D**2*k + D*k*(u-1) - D*(u-1) + u**2)*Z - (D**2*k**2 + 2*D*k*u + u**2)
claims7 = {0: (Y-AZ)/(Y+1), 1: Y*Lc/(Z*(Y+1)), 2: -AZ/(Y+1), 4: (Y*Lc+Z)/(Z*(Y+1)), 9: -D*(u-1)*Z + (D+1)*s_, 6: N6/Lc}
for i, c in claims7.items():
    chk("(7) H%d" % i, z0(H2[i] - c), sp.factor(H2[i]))
chk("H6 leading term D(u-1)Z/u", sp.limit(H2[6]/Z, Z, sp.oo) == sp.cancel(D*(u-1)/u) or z0(sp.limit(H2[6]/Z, Z, sp.oo) - D*(u-1)/u))
den_all = Z*Lc*(Y+1)*K
bad = []
for i in range(24):
    n_, d_ = sp.fraction(sp.cancel(H2[i]))
    q, rem = sp.div(sp.Poly(den_all, Y, Z, D, k, u), sp.Poly(d_, Y, Z, D, k, u))
    # divisibility: den_all / d_ must be polynomial
    if not sp.cancel(den_all / d_).is_polynomial(Y, Z, D, k, u):
        bad.append((i, sp.factor(d_)))
chk("all 24 H denominators divide Z L (Y+1) K", not bad, bad)
# matching coordinates v = P^{-1}(H) : denominators
# coverage of the new chart: original config -> (D,k,u,Y,Z) -> inverse == Mobius image
subsd = {D: Dex, k: kex, u: H[10], Y: 1/(x-1), Z: (z-1)/(x-1)}
mob = lambda xi: (xi - 1)/(x - 1)
origA = [y, z, w, x]; origB = [t, 0, 1, None]
okc = True
for P_, o in list(zip(A2, origA)) + list(zip(B2, origB)):
    P0 = sp.cancel(P_[0].subs(subsd)); P1 = sp.cancel(P_[1].subs(subsd))
    if o is None:
        okc &= z0(P1)
    else:
        okc &= z0(P0 - mob(o)*P1) and not z0(P1)
chk("new chart inverse reproduces Mobius image (x->Y(xi-1)) of every original 8-point config", okc)

# ---------- the curve
n12, d12 = sp.fraction(sp.cancel(H2[12]))
C = sp.expand(n12 - h*d12)
Cp = sp.Poly(C, Y, Z)
print("INFO H12 =", sp.factor(H2[12]))
chk("curve bidegree (2,4)", (sp.degree(C, Y), sp.degree(C, Z)) == (2, 4), (sp.degree(C, Y), sp.degree(C, Z)))
nterms = len(sp.Poly(C, Y, Z, D, k, u, h).terms())
print("INFO curve terms in (Y,Z,D,k,u,h):", nterms, " (doc: 160)")
cZ4 = sp.factor(sp.Poly(C, Z).coeff_monomial(Z**4))
print("INFO [Z^4]C =", cZ4)
cnorm = sp.cancel(cZ4 / (-D*(u-1)**2*(Y+1)))
print("INFO normalisation constant of C relative to doc ([Z^4] = -D(u-1)^2(Y+1)):", cnorm)
Cn = sp.expand(sp.cancel(C / cnorm))
chk("[Z^4]C = -D(u-1)^2(Y+1) (after constant normalisation)", cnorm.is_number)
CP = sp.Poly(Cn, Y, Z)
def co(i, j): return sp.factor(CP.coeff_monomial(Y**i * Z**j))
print("INFO [Y^2Z^2] =", co(2, 2), " [Y^0Z^2] =", co(0, 2), " [Y^1Z^2] =", co(1, 2))
chk("Lemma2: u=1 => [Y^2Z^2] forces h=-Dk(D+1)", z0(co(2, 2).subs(u, 1).subs(h, -D*k*(D+1))) and sp.degree(sp.expand(co(2,2).subs(u,1)), h) == 1)
c02 = sp.factor(co(0, 2).subs({u: 1, h: -D*k*(D+1)}))
chk("Lemma2: then [Y^0Z^2] = k(D+1)(Dk+1) (up to sign)", z0(c02 - k*(D+1)*(D*k+1)) or z0(c02 + k*(D+1)*(D*k+1)), c02)
c12 = sp.factor(co(1, 2).subs({u: 1, D: -1, h: 0}))
chk("Lemma2: D=-1,h=0 => [Y^1Z^2] = k-1 (up to sign)", z0(c12 - (k-1)) or z0(c12 + (k-1)), c12)
chk("VN - L = (Z-1)((u-1)Z - s)", sp.expand(VN - Lc - (Z-1)*AZ) == 0)
# u=1,s=0 => is C identically zero?
C_deg = sp.expand(Cn.subs({u: 1, k: -1/D}))
print("INFO C at u=1, s=0 (k=-1/D):", sp.factor(C_deg))
# Lemma 2 exhaustiveness: solve C == 0 identically
coeffs = sp.Poly(Cn, Y, Z).coeffs()
sol = sp.solve(coeffs, [h, u, k], dict=True)
print("INFO parameter solutions making C identically 0:", sol)

# ---------- resultants (8)
def res(a, b, var): return sp.factor(sp.resultant(sp.Poly(a, var), sp.Poly(b, var)))
def same_upto_const(a, b):
    q = sp.cancel(a / b); return q.is_number and q != 0, q
r8 = [("ResZ(C,Z)", res(Cn, Z, Z), -Y*(Y*(D*k+k+u-h-1) + s_ - h)*s_**2),
      ("ResY(C,Y+1)", res(Cn, Y+1, Y), -(Z+k-1)*Lc*AZ),
      ("ResZ(C,L)", res(Cn, Lc, Z), -D*(Y+1)*k**2*r_**2*s_**2),
      ("ResY(C,K)", res(Cn, K, Y), k*(Z+k-1)*Lc*sp.cancel(TN/Z)*(-claims7[9])*Z**2)]
for nm, got, exp in r8:
    ok, q = same_upto_const(got, exp)
    chk("(8) %s matches up to a constant" % nm, ok, "ratio %s ; got %s" % (q, got))
# point L=0, Y=-1, Z=s/u
pt = {Y: -1, Z: s_/u}
CY = sp.factor(sp.diff(Cn, Y).subs(pt)); CZ = sp.factor(sp.diff(Cn, Z).subs(pt))
chk("point (Y,Z)=(-1,s/u) on curve", z0(Cn.subs(pt).subs(h, h)), sp.factor(Cn.subs(pt)))
print("INFO C_Y, C_Z at that point:", CY, CZ)
chk("C_Y = -Dk^2r^2s^2/u^4", z0(CY + D*k**2*r_**2*s_**2/u**4))
chk("C_Z = krs/u", z0(CZ - k*r_*s_/u))
lim_ratio = sp.cancel(-CZ / (CY * sp.diff(Lc, Z)))
chk("(Y+1)/L -> u^2/(Dkrs)", z0(lim_ratio - u**2/(D*k*r_*s_)))
H1lim = sp.cancel((Y*Lc/(Z*(Y+1))).subs(Y, -1).subs(Z, s_/u) if False else (-1)*(1/lim_ratio)/(s_/u))
chk("H1 -> -Dkr/u", z0(H1lim + D*k*r_/u))
chk("H6(0) finite real (L(0) != 0)", z0(Lc.subs(Z, 0) + s_), sp.cancel(H2[6].subs(Z, 0)))
chk("L(1) = -Dk", z0(Lc.subs(Z, 1) + D*k))
chk("Y=oo: H0 -> 1", sp.limit(claims7[0], Y, sp.oo) == 1)
# Z = oo
chk("[Z^3]C at Y=-1 = -u(u-1)", z0(sp.Poly(Cn, Z).coeff_monomial(Z**3).subs(Y, -1) + u*(u-1)))
V, tau = sp.symbols('V tau')
eq10 = sp.factor(sp.expand(tau**5 * Cn.subs({Y: V/tau, Z: 1/tau})).subs(tau, 0))
chk("(10) tau^5 C(V/tau,1/tau)|_0 = -D(u-1)V(uV+u-1)", z0(eq10 + D*(u-1)*V*(u*V+u-1)), eq10)
chk("second branch: H0 -> 1+u", z0(sp.limit(claims7[0].subs(Y, V*Z).subs(V, (1-u)/u), Z, sp.oo) - (1+u)))
ell = -u/(D*(u-1))
# implicit branch Y = -1 + ell tau + c2 tau^2 + c3 tau^3
c2, c3 = sp.symbols('c2 c3')
Ys = -1 + ell*tau + c2*tau**2 + c3*tau**3
ser = sp.expand(sp.cancel(tau**4 * Cn.subs({Y: Ys, Z: 1/tau})))
ser = sp.Poly(ser, tau)
co_ = [sp.factor(ser.coeff_monomial(tau**i)) for i in range(0, 4)]
chk("branch (9): tau^0 coefficient vanishes", z0(co_[0]))
chk("branch (9): tau^1 coefficient vanishes for ell=-u/(D(u-1))", z0(co_[1]), co_[1])
c2v = sp.solve(co_[2], c2); print("INFO c2 =", c2v)
c3v = sp.solve(co_[3].subs(c2, c2v[0]), c3)
Ysr = Ys.subs({c2: c2v[0], c3: c3v[0]})
Hcomb = (claims7[4] + (1-u)*claims7[6]).subs({Y: Ysr, Z: 1/tau})
lser = sp.series(sp.cancel(Hcomb), tau, 0, 1).removeO()
lser = sp.expand(lser)
print("INFO H4+(1-u)H6 on branch (9), series to O(tau):", sp.simplify(lser))
chk("H4+(1-u)H6 has no pole on branch (9)", sp.limit(sp.cancel(Hcomb*tau), tau, 0) == 0)
H4p = sp.limit(sp.cancel(claims7[4].subs({Y: Ysr, Z: 1/tau})*tau), tau, 0)
H6p = sp.limit(sp.cancel(claims7[6].subs({Z: 1/tau})*tau), tau, 0)
chk("1/tau coeffs: H4 -> (1-u)/ell, H6 -> -1/ell", z0(H4p - (1-u)/ell) and z0(H6p + 1/ell), (sp.factor(H4p), sp.factor(H6p)))
lim69 = sp.limit(claims7[6] + claims7[9]/u, Z, sp.oo)
chk("H6 + H9/u finite at Z=oo", lim69.is_finite is not False and not lim69.has(sp.oo), sp.factor(lim69))
# generic specialisation (2,3,4,5)
spec = {D: 2, k: 3, u: 4, h: 5}
vals = [e.subs(spec) for e in [D, k, u, u-1, r_, s_, D*k+k+u-h-1]]
c_Y2 = sp.Poly(Cn, Y).coeff_monomial(Y**2).subs(spec); c_Y0 = sp.Poly(Cn, Y).coeff_monomial(Y**0).subs(spec)
chk("generic exclusions at (2,3,4,5)", all(v != 0 for v in vals) and sp.degree(sp.gcd(c_Y2, c_Y0), Z) == 0, (vals, sp.gcd(c_Y2, c_Y0)))
print("INFO [Y^2]C =", sp.factor(sp.Poly(Cn, Y).coeff_monomial(Y**2)))
print("INFO [Y^0]C =", sp.factor(sp.Poly(Cn, Y).coeff_monomial(Y**0)))
print("INFO [Y^1]C =", sp.factor(sp.Poly(Cn, Y).coeff_monomial(Y**1)))

if NUMERIC:
    import numpy as np
    rng = np.random.default_rng(20261125)
    Pn = L.P_numpy()
    def F_hom(A, B):
        out = []
        for s_ in PERMS:
            v = 1.0 + 0j
            for i in range(4):
                P_, Q_ = A[i], B[s_[i]]
                v *= P_[0]*Q_[1] - P_[1]*Q_[0]
            out.append(v)
        return np.array(out)
    cf = [sp.lambdify((Z, D, k, u, h), sp.Poly(Cn, Y).coeff_monomial(Y**i), "numpy") for i in (2, 1, 0)]
    Af = sp.lambdify((Y, Z, D, k, u), [list(p) for p in A2], "numpy"); Bf = sp.lambdify((Y, Z, D, k, u), [list(p) for p in B2], "numpy")
    NS = 120000; nsup = 0; nphi = 0; worst = np.inf; res = 0.0
    for it in range(NS):
        Dv, kv, uv, hv = rng.normal(size=4) * 2
        Zv = complex(rng.normal(), rng.normal()) * 1.5
        for Yv in np.roots([complex(f(Zv, Dv, kv, uv, hv)) for f in cf]):
            A = [tuple(complex(q) for q in p) for p in Af(Yv, Zv, Dv, kv, uv)]
            B = [tuple(complex(q) for q in p) for p in Bf(Yv, Zv, Dv, kv, uv)]
            F_ = F_hom(A, B); E_ = F_[5] - F_[3]
            if abs(E_) < 1e-8: continue
            H_ = (F_ - F_[3]) / E_
            res = max(res, max(abs(H_[s_].imag) for s_ in T))
            for eta in (1.0, -1.0):
                gg = -eta * H_.imag
                if gg[6] < 0 or min(gg[s_] for s_ in range(24) if s_ not in T) < 1e-9: continue
                nsup += 1
                if nphi < 40:
                    c_ = np.linalg.lstsq(Pn, gg, rcond=None)[0]
                    val, O_ = L.min_phi_over_SO4(c_, restarts=6, seed=it); nphi += 1
                    worst = min(worst, val / np.max(np.abs(gg)))
    print("NUMERIC T_75 curve sampling: seed 20261125, %d parameter draws, %d strictly supported points, contact residual %.1e, min over SO(4) (6 BFGS restarts) on %d of them, worst min Phi / max gap = %.2e" % (NS, nsup, res, nphi, worst))

print("SUMMARY failures:", FAIL)
sys.exit(1 if FAIL else 0)
