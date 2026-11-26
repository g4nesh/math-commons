"""Own-code replay for the Zheng n=4 audit, document root_inventory72_native_uniform_20261023
(primary proof of paper case T_73 = {0,4,5,7,9,12}).
Sigma-points are rebuilt from the definitions: A = (a,z,x,w), B = (t,0,1,inf) with plus brackets
[i,4+j] = A_i + B_j (the factor of the point at infinity is a common constant and is dropped).
Nothing from the supplement is read, imported or executed.  Exit code 1 if an EXACT check fails."""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name, flush=True)
    if not cond:
        FAIL.append(name)
ok = lambda e: sp.expand(sp.fraction(sp.cancel(sp.together(e)))[0]) == 0
t0 = time.time()
FV = L.F_VALUES

# ---------- A/B: generators, coefficient bounds, extra ties ----------
T = set(L.CONTACT_SETS[73])
gen = L.generators(T)
check("T_73 = {0,4,5,7,9,12}; generators = [0,1,2,5,7,12,35,45] (doc (1))", sorted(T) == [0, 4, 5, 7, 9, 12] and gen == [0, 1, 2, 5, 7, 12, 35, 45])
check("rank 8, all minor except f35 = h4312,-", L.rank([list(L.F_VECS[j]) for j in gen]) == 8
      and [L.F_IS_MINOR[j] for j in gen] == [True] * 6 + [False, True])
check("(2) f35 = 1 at rows 2, 8, 10 and every generator is >= 0 on all rows, so 0 <= theta <= g2, g8, g10",
      [FV[35][s] for s in (2, 8, 10)] == [1, 1, 1] and all(FV[j][s] >= 0 for j in gen for s in range(24)))
ext = [j for j in range(24) if j not in T and FV[35][j] == 0]
check("extra ties not forcing theta = 0 are exactly 1, 6, 11, 14, 15, 17, 18", ext == [1, 6, 11, 14, 15, 17, 18])
def cover(S):
    out = []
    for k, Tk in L.CONTACT_SETS.items():
        if k == 73:
            continue
        if any(L.g_act_set(gi, Tk) <= S for gi in range(1152)):
            out.append(k)
    return out
for j in ext:
    S = T | {j}
    cv = cover(S); a7 = [k for k, A in L.AUX7.items() if L.find_g(A, S)]
    print(f"info tie {j}: paper cases T_k with g.T_k inside T u {{{j}}}: {cv}; Table-5 seven-sets equal to an image: {a7}")
    if j != 15:
        check(f"tie {j} is covered by some paper case T_k with k < 73", any(k < 73 for k in cv))
check("tie 15: T u {15} is an image of E = {0,1,2,9,10,16,17} = Table 5 set #17",
      L.AUX7[17] == frozenset({0, 1, 2, 9, 10, 16, 17}) and L.find_g(L.AUX7[17], T | {15}) is not None)

# ---------- C: coordinates ----------
a, z, x, w, t = sp.symbols('a z x w t')
J, r, s, u, h = sp.symbols('J r s u h')
A = [a, z, x, w]; B = [t, 0, 1, None]
def Fk(k):
    sig = L.PERMS[k]; p = 1
    for i in range(4):
        if B[sig[i]] is not None:
            p *= (A[i] + B[sig[i]])
    return sp.expand(p)
F = [Fk(k) for k in range(24)]
E = F[5] - F[0]
check("(1) E = F5 - F0 = (w-z)(x+1)(a+t)", ok(E - (w - z) * (x + 1) * (a + t)))
check("(2) F4 - F5 = (x-w)(a+t)", ok(F[4] - F[5] - (x - w) * (a + t)))
check("(4) F9 - F7 = a(t-1)(z-w)", ok(F[9] - F[7] - a * (t - 1) * (z - w)))
check("F7 = a(z+t)(w+1), F0 = (a+t)z(x+1)", ok(F[7] - a * (z + t) * (w + 1)) and ok(F[0] - (a + t) * z * (x + 1)))
H = [sp.cancel((F[i] - F[0]) / E) for i in range(24)]
rr = H[4] - 1; uu = H[7]; ss = H[9] - H[7]
check("(6) J = H2 = (x-z)/((w-z)(x+1)) and J - r = 1/(x+1)", ok(H[2] - (x - z) / ((w - z) * (x + 1))) and ok(H[2] - rr - 1 / (x + 1)))
qq = 1 / (x + 1); D = qq * (w + 1) * (z + 1); N = z * (1 - uu) + uu * w + ss * (w + 1)
check("raw identity a(t-1)/(a+t) = -s/q and N/D = a/(a+t) (7)", ok(a * (t - 1) / (a + t) + ss / qq) and ok(N / D - a / (a + t)))
q_ = J - r; x_ = 1 / q_ - 1; z_ = (J * (w + 1) - r - 1) / r
D_ = q_ * (w + 1) * (z_ + 1); N_ = z_ * (1 - u) + u * w + s * (w + 1)
sub = {a: (q_ * N_ - s * D_) / (q_ * (D_ - N_)), z: z_, x: x_, t: (q_ * N_ - s * D_) / (q_ * N_)}
Hs = [sp.cancel(sp.together(Hi.subs(sub, simultaneous=True))) for Hi in H]
check("inverse (6),(8): H0 = 0, H5 = 1, H4 = 1+r, H7 = u, H9 = u+s identically",
      Hs[0] == 0 and Hs[5] == 1 and ok(Hs[4] - 1 - r) and ok(Hs[7] - u) and ok(Hs[9] - u - s))
check("(3) H2 = J, H3 = 1 - wJ", ok(Hs[2] - J) and ok(Hs[3] - (1 - w * J)))
H8c = (-J * (s + 1) * (w + 1) + r - u + 1) / ((r - J) * (w + 1))
H10c = (-J**2 * (w + 1) - J * r * s * (w + 1) + J * (r + 1) - r * u) / ((r - J) * (J * (w + 1) - 1))
check("(9) H8 formula", ok(Hs[8] - H8c))
check("(10) H10 formula", ok(Hs[10] - H10c))
num, den = sp.fraction(sp.cancel(Hs[12] - h))
Pn = sp.Poly(sp.expand(num), J, w)
c3 = sp.expand(Pn.as_expr().coeff(J, 3))
sgn = 1 if ok(c3 - (u - 1) * (w + 1)) else (-1 if ok(c3 + (u - 1) * (w + 1)) else 0)
Cjw = sp.expand(sgn * num)
check("last equation C(J,w): bidegree (3,2), 60 terms, [J^3]C = (u-1)(w+1)",
      sgn != 0 and Pn.degree(J) == 3 and Pn.degree(w) == 2 and len(sp.Poly(Cjw, J, w, r, s, u, h).terms()) == 60)
check("(11) [w^2]C = J s ((1-h)J + h r)", ok(Cjw.coeff(w, 2) - J * s * ((1 - h) * J + h * r)))
check("[J^2]C at w = -1 equals r - u + 1", ok(Cjw.coeff(J, 2).subs(w, -1) - (r - u + 1)))
V0 = (u - r - 1) / (u - 1)
check("(13) V0 - 1 = -r/(u-1)", ok(V0 - 1 + r / (u - 1)))
check("(14) (V0 - r - 1)/(V0 - 1) = u", ok((V0 - r - 1) / (V0 - 1) - u))
# (12) and (14) as limits along the w = -1 branch, from the series of w(J)
tau, v1, v2 = sp.symbols('tau v1 v2')
Ct = sp.expand(tau**3 * Cjw.subs({J: 1 / tau}))
check("C~(tau,w) = tau^3 C(1/tau,w) vanishes at (0,-1) with dC~/dw(0,-1) = u-1",
      ok(Ct.subs({tau: 0, w: -1})) and ok(sp.diff(Ct, w).subs({tau: 0, w: -1}) - (u - 1)))
ser = sp.expand(Ct.subs(w, -1 + v1 * tau + v2 * tau**2))
s1 = sp.solve(ser.coeff(tau, 1), v1)[0]
check("branch at (J,w) = (inf,-1): w + 1 = V0/J + O(J^-2)", ok(s1 - V0))
s2 = sp.solve(sp.expand(ser.coeff(tau, 2).subs(v1, s1)), v2)[0]
H10b = sp.cancel(sp.together(H10c.subs({J: 1 / tau, w: -1 + s1 * tau + s2 * tau**2})))
nb, db = sp.fraction(H10b)
check("H10 -> u along that branch (leading terms of numerator/denominator)",
      ok(sp.expand(nb).subs(tau, 0) - u * sp.expand(db).subs(tau, 0)) and not ok(sp.expand(db).subs(tau, 0)))
lam = sp.symbols('lam')   # (12): J, w -> inf with arbitrary relation: H8 = (-(s+1) + (r-u+1)/(J(w+1)))/(r/J - 1)
check("(12) rewritten H8 = (-J(s+1) + (r-u+1)/(w+1))/(r-J)", ok(H8c - (-J * (s + 1) + (r - u + 1) / (w + 1)) / (r - J)))
print("info exact part time", round(time.time() - t0, 1))

# ---------- NUMERIC (brief) ----------
import numpy as np
from fractions import Fraction
found = 0; minphi = np.inf; worst = 0.0
for seed in range(6):
    pt = L.random_point_on_contact_stratum(sorted(T), seed=seed, tries=30, margin=1e-2)
    if pt is None:
        continue
    found += 1
    c_vec = np.linalg.lstsq(np.array(L.P_numpy(), float), np.array(pt['gaps'], float), rcond=None)[0]
    mv, _ = L.min_phi_over_SO4(c_vec, restarts=8, seed=seed)
    minphi = min(minphi, mv)
    # parametrization consistency: Moebius-normalize b2,b3,b4 -> 0,1,inf and evaluate C(J,w) at the point
    aa = np.array(pt['a']); bb = np.array(pt['b'])
    pts = list(aa) + list(-bb)          # projective points z_i; brackets z_i - z_{4+j}
    zb2, zb3, zb4 = -bb[1], -bb[2], -bb[3]
    mob = lambda zz: -((zz - zb2) * (zb3 - zb4)) / ((zz - zb4) * (zb3 - zb2))   # sends -b2,-b3,-b4 to 0,-1,inf (B = (t,0,1,inf))
    An = [mob(zz) for zz in aa]; tn = -mob(-bb[0])
    vals = {a: An[0], z: An[1], x: An[2], w: An[3], t: tn}
    Hn = [complex(sp.N(H[i].subs(vals))) for i in range(24)]
    Jn, wn = Hn[2], An[3]
    rn, un, sn, hn = (Hn[4] - 1).real, Hn[7].real, (Hn[9] - Hn[7]).real, Hn[12].real
    resid = complex(sp.N(Cjw.subs({J: Jn, w: wn, r: rn, u: un, s: sn, h: hn})))
    scale = sum(abs(complex(sp.N(term.subs({J: Jn, w: wn, r: rn, u: un, s: sn, h: hn})))) for term in Cjw.args)
    worst = max(worst, abs(resid) / scale)
    # contact values real?
    worst = max(worst, max(abs(Hn[i].imag) for i in (4, 7, 9, 12)) / (1 + max(abs(Hn[i]) for i in (4, 7, 9, 12))))
print(f"NUMERIC T_73 stratum points found {found}/6 (seeds 0-5, margin 1e-2): min Phi over SO(4) (8 restarts each) = {minphi:.4g}; "
      f"max relative residual of C(J,w)=0 and of Im(contact H) at these points = {worst:.2e}")
print("time", round(time.time() - t0, 1))
print("SUMMARY", "all EXACT checks passed" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
