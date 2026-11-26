"""Audit replay (batch C01-T1-T3).  Exact replay for independent_moc4_six_contact_reduction_20261022 (paper T_1, T_2)."""
import os, sys, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("MOC4LIB", os.path.join(HERE, "..", "lib")))
# Own code only.  Nothing from the supplement is read or executed.
# Run: python -I -B THIS.py [--numeric]   (--numeric adds the floating-point tests)
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)

T1 = frozenset({0, 1, 2, 3, 6, 7}); T2 = frozenset({0, 1, 2, 3, 6, 8})
check("doc_T_equals_paper_T1", L.CONTACT_SETS[1] == T1)
check("doc_T1_equals_paper_T2", L.CONTACT_SETS[2] == T2)

# ---------- (M2), (M3), Lemma 4 ingredients ----------
m = sp.symbols("m0:14"); t = sp.Symbol("t")
syms, C = L.cubic_sympy()
syms = list(syms)
grad = [sp.diff(C, s) for s in syms]
kap = [0] * 13 + [1]
gk = [g.subs(dict(zip(syms, kap))) for g in grad]
check("M2_gradF_kappa_zero", all(sp.simplify(v) == 0 for v in gk))
H = sp.Matrix([[sp.diff(C, si, sj).subs(dict(zip(syms, kap))) for sj in syms] for si in syms])
xs = sp.symbols("x0:14")
shift = {syms[i]: xs[i] + t * kap[i] for i in range(14)}
ok = True
for i in range(14):
    lhs = sp.expand(grad[i].subs(shift, simultaneous=True))
    rhs = sp.expand(grad[i].subs(dict(zip(syms, xs))) + t * sum(H[i, j] * xs[j] for j in range(14)))
    ok &= sp.expand(lhs - rhs) == 0
check("M2_eq2_gradF(x+t kappa)=gradF(x)+t H x", ok)
check("M3_rank_H_kappa_4", H.rank() == 4)
E = [L.n_vec((i,), (j,)) for i in range(1, 5) for j in range(1, 5)]
check("M3_singleton_vectors_rank_10", L.rank(E) == 10)
check("M3_H_e_ij_zero", all(all(v == 0 for v in H * sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in e])) for e in E))
check("M3_P_e_ij_profile", all(L.values(L.n_vec((i,), (j,)))[s] == (1 if L.PERMS[s][i - 1] == j - 1 else 0)
                                for i in range(1, 5) for j in range(1, 5) for s in range(24)))
check("M1_rank_P_14_and_P_kappa_1", L.rank(L.P) == 14 and L.values(L.KAPPA) == [1] * 24)
diffs = [[L.P[s][r] - L.P[0][r] for r in range(14)] for s in range(24)]
check("Lemma4_gap_forms_common_kernel_is_C_kappa(rank 13)", L.rank(diffs) == 13)

# ---------- Proposition 8 inventory ----------
orbs = L.enumerate_six_set_orbits()
check("Prop8_214_orbits", len(orbs) == 214 and sum(s for _, s in orbs) == 134596)
def face_vertices(gens):
    Z = set(range(24))
    for j in gens:
        Z &= L.F_ZERO[j]
    return Z
n90 = 0; nface = 0; dist = {}; reps90_face = True; nfour = 0
four = {2, 7, 19, 20}
four_imgs = set(L.g_act_set(g, four) for g in range(1152))
for T, size in orbs:
    gens = L.generators(T)
    nm = [j for j in gens if not L.F_IS_MINOR[j]]
    exact_face = bool(gens) and face_vertices(gens) == set(T)
    nface += exact_face
    if nm:
        n90 += 1
        dist[len(nm)] = dist.get(len(nm), 0) + 1
        if not exact_face:
            reps90_face = False
        if any(F <= set(T) for F in four_imgs):
            nfour += 1
check("Prop8_90_orbits_with_nonminor_generator", n90 == 90, str(n90))
check("Prop8_113_orbits_exact_six_vertex_faces", nface == 113, str(nface))
check("Prop8_nonminor_count_distribution_51_35_3_1", dist == {1: 51, 2: 35, 3: 3, 4: 1}, str(dist))
check("Prop8_all_90_reps_exact_six_vertex_faces", reps90_face)
check("Prop8_9_retained_orbits_contain_image_of_{2,7,19,20}", nfour == 9, str(nfour))

# ---------- entry faces ----------
ok = True
for i in range(4):
    for j in range(4):
        Fc = [s for s in range(24) if L.PERMS[s][i] == j]
        g = L.generators(Fc)
        n1 = sum(1 for k in g if len(L.F_TABLE[k][1]) == 1) if False else None
        ok &= len(g) == 15 and all(L.F_IS_MINOR[k] for k in g)
check("entry_faces_have_15_generators_all_minor", ok)
Fc = [s for s in range(24) if L.PERMS[s][0] == 0]
g = L.generators(Fc)
names = [L.f_name(k) for k in g]
print("INFO entry face (0,0) generators:", names)

# ---------- Proposition 10 / 13 factor identities ----------
A = sp.symbols("A0:4"); Cc = sp.symbols("C0:4"); lam = sp.Symbol("lam")
def Fs(s):
    p = L.PERMS[s]
    return lam * sp.Mul(*[A[i] - Cc[p[i]] for i in range(4)])
def D(r, s):
    return sp.expand(Fs(r) - Fs(s))
z = lambda e: sp.expand(e) == 0
check("P10_D06_factor", z(D(0, 6) - lam * (A[0] - A[1]) * (Cc[0] - Cc[1]) * (A[2] - Cc[2]) * (A[3] - Cc[3])))
check("P10_D17_factor", z(D(1, 7) - lam * (A[0] - A[1]) * (Cc[0] - Cc[1]) * (A[2] - Cc[3]) * (A[3] - Cc[2])))
check("P10_D02_factor", z(D(0, 2) - lam * (A[0] - Cc[0]) * (A[3] - Cc[3]) * (A[1] - A[2]) * (Cc[1] - Cc[2])))
check("P10_D13_factor", z(D(1, 3) - lam * (A[0] - Cc[0]) * (A[2] - Cc[3]) * (A[1] - A[3]) * (Cc[1] - Cc[2])))
check("P10_D01_factor", z(D(0, 1) - lam * (A[0] - Cc[0]) * (A[1] - Cc[1]) * (A[2] - A[3]) * (Cc[2] - Cc[3])))
check("P10_D23_factor", z(D(2, 3) - lam * (A[0] - Cc[0]) * (A[1] - Cc[2]) * (A[2] - A[3]) * (Cc[1] - Cc[3])))
check("P10_ratio1", z(D(0, 6) * (A[2] - Cc[3]) * (A[3] - Cc[2]) - D(1, 7) * (A[2] - Cc[2]) * (A[3] - Cc[3])))
check("P10_ratio2", z(D(0, 2) * (A[2] - Cc[3]) * (A[1] - A[3]) - D(1, 3) * (A[3] - Cc[3]) * (A[1] - A[2])))
check("P10_ratio3", z(D(0, 1) * (A[1] - Cc[2]) * (Cc[1] - Cc[3]) - D(2, 3) * (A[1] - Cc[1]) * (Cc[2] - Cc[3])))
# imaginary numerators
x, y, tt, u, v, p, q = sp.symbols("x y t u v p q", real=True)
A1 = u + sp.I * v
r2 = (y * (A1 - x)) * sp.conjugate(x * (A1 - y))     # ratio2 numerator * conj(denominator), C3=0
check("P10_imag_numerator_ratio2_eq_xy(x-y)v", z(sp.im(sp.expand(r2)) - x * y * (x - y) * v))
beta = p + sp.I * q
r3 = ((u - beta) * tt) * sp.conjugate((u - tt) * beta)  # ratio3 with A1=u, C1=beta, C2=t, C3=0
check("P10_imag_numerator_ratio3_eq_-tu(u-t)q", z(sp.im(sp.expand(r3)) + tt * u * (u - tt) * q))
# exceptional collision u = 0 (A1 = C3 = 0)
Lsym = sp.Symbol("L")
sub = {A[1]: 0, A[2]: x, A[3]: y, Cc[1]: sp.Symbol("beta"), Cc[2]: tt, Cc[3]: 0}
b_ = sp.Symbol("beta")
d01 = (D(0, 1).subs(sub)); d02 = (D(0, 2).subs(sub))
Lval = lam * (A[0] - Cc[0])
check("P10_exceptional_D01", z(d01 - (-Lval * b_ * tt * (x - y))))
check("P10_exceptional_D02", z(d02 - (-Lval * x * y * (b_ - tt))))
# Prop 13 identity D79 D02 = D68 D13
check("P13_D79_D02_eq_D68_D13", z(D(7, 9) * D(0, 2) - D(6, 8) * D(1, 3)))
# Prop 13 facet identities (paper numbering = supplement numbering, checked separately)
J = [0, 1, 2, 3, 4, 5, 6, 9, 10, 15, 16, 17]
check("P13_generators_T2_equal_J", L.generators(T2) == J, str(L.generators(T2)))
check("P13_only_f17_nonminor", [j for j in J if not L.F_IS_MINOR[j]] == [17])
V = L.F_VALUES
check("P13_f17+f9=f4+f10+f16", all(V[17][s] + V[9][s] == V[4][s] + V[10][s] + V[16][s] for s in range(24)))
ell = {j: V[j][7] - V[j][9] for j in J}
check("P13_ell_values", ell == {j: (-1 if j in (6, 17) else 1 if j == 9 else 0) for j in J}, str(ell))
# coefficient transfer with formal nonnegative coefficients
a = {j: sp.Symbol(f"a{j}", nonnegative=True) for j in J}
sol_a9 = a[6] + a[17]
beta_c = dict(a); beta_c.pop(17); beta_c[9] = a[6]; beta_c[4] = a[4] + a[17]; beta_c[10] = a[10] + a[17]; beta_c[16] = a[16] + a[17]
ok = True
for s in range(24):
    lhs = sum(a[j] * V[j][s] for j in J).subs(a[9], sol_a9)
    rhs = sum(beta_c[j] * V[j][s] for j in beta_c)
    ok &= sp.expand(lhs - rhs) == 0
check("P13_coefficient_transfer_exact", ok)
# T1 generators
g1 = L.generators(T1)
print("INFO T1 generators", g1, "nonminor", [j for j in g1 if not L.F_IS_MINOR[j]])

if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    for which in ("T12",):

        PN = L.P_numpy()
        FV = np.array(L.F_VALUES, float)


        def phi_vals(vals, O):
            return float(np.asarray(vals) @ L.leibniz_weights_np(O))


        def min_phi_vals(vals, restarts, rng):
            f = lambda th: phi_vals(vals, L.so4_from_params(th))
            best = np.inf
            for _ in range(restarts):
                r = minimize(f, rng.uniform(-np.pi, np.pi, 6), method="BFGS")
                best = min(best, r.fun)
            return best


        if which in ("soc", "all"):
            # SOC input of the T_3 note: Phi(a f1 + b f11 + c f17) >= 0 when ab >= c^2; worst case ab = c^2, c = 1
            for trip, seed in (((1, 11, 17), 11), ((9, 1, 22), 12)):
                X, Y, H = trip
                rng = np.random.default_rng(seed)
                def f(z):
                    O = L.so4_from_params(z[:6]); tau = np.exp(z[6])
                    return phi_vals(tau * FV[X] + FV[Y] / tau + FV[H], O)
                best = np.inf
                for _ in range(80):
                    z0 = np.concatenate([rng.uniform(-np.pi, np.pi, 6), rng.normal(0, 1.5, 1)])
                    r = minimize(f, z0, method="BFGS")
                    best = min(best, r.fun)
                # direct sampling of Phi(fH) + 2 sqrt(Phi fX Phi fY)
                worst = np.inf
                for _ in range(5000):
                    O = L.so4_from_params(rng.uniform(-np.pi, np.pi, 6))
                    px, py, ph = (phi_vals(FV[k], O) for k in trip)
                    worst = min(worst, ph + 2 * np.sqrt(max(px, 0) * max(py, 0)))
                print(f"NUMERIC SOC triple {trip}: min over SO(4) x tau of Phi(tau f{X} + f{Y}/tau + f{H}) = {best:.3e} "
                      f"(80 BFGS restarts, seed {seed}); min of Phi(f{H})+2sqrt(Phi f{X} Phi f{Y}) on 5000 random O = {worst:.3e}")

        if which in ("scalar", "all"):
            rng = np.random.default_rng(5)
            worst = np.inf; worst3 = np.inf
            for k in range(50000):
                x = rng.exponential(1, 4) * (rng.random(4) > 0.25)
                y = rng.exponential(1, 4) * (rng.random(4) > 0.25)
                if x.sum() == 0 or y.sum() == 0:
                    continue
                y = y * x.sum() / y.sum()
                a, b, c, d = x; A, B, C, D = y
                s = c + d - C - D; t = b + c - B - C
                R = 2 * (a * C + b * D + c * A + d * B) - 4 * np.sqrt(a * b * c * d) - 4 * np.sqrt(A * B * C * D)
                worst = min(worst, (R + 2 * abs(s * t)) / (x.sum() ** 2))
                p = rng.normal(size=4); q = rng.normal(size=4); q *= np.linalg.norm(p) / np.linalg.norm(q)
                S = p[2]**2 + p[3]**2 - q[2]**2 - q[3]**2; T = p[1]**2 + p[2]**2 - q[1]**2 - q[2]**2
                H = 2 * sum(p[i]**2 * q[(i + 2) % 4]**2 for i in range(4)) - 4 * np.prod(p) + 4 * np.prod(q)
                worst3 = min(worst3, (H + 2 * abs(S * T)) / np.dot(p, p) ** 2)
            print(f"NUMERIC Thm2 min (R+2|st|)/mass^2 over 50000 samples (seed 5, ~25% zero coords) = {worst:.3e}; "
                  f"Thm3 min (H+2|ST|)/|p|^4 = {worst3:.3e}")

        if which in ("T3", "all"):
            rng = np.random.default_rng(7)
            feas = 0; Aneg = 0; worst_soc = np.inf; worst_phi = np.inf; tried = 0; nphi = 0
            while tried < 400000 and feas < 400:
                tried += 1
                xr, t, u, r = rng.normal(0, 2, 4)
                v = complex(*rng.normal(0, 1.5, 2)); mu = rng.normal()
                z = u * v / ((u + t) * v - xr)
                w = (v * (z - 1) - r * z) / (z - 1 - r)
                kap = 1j * mu / (v * (z - w))
                Av = [z, 0, 1, xr]; Cv = [w, v, t * v, None]
                F = np.array([-kap * np.prod([Av[i] - Cv[p[i]] for i in range(4) if p[i] != 3]) for p in L.PERMS])
                g = F.real - F.real[0]
                sc = np.abs(g).max()
                if sc == 0 or not np.isfinite(sc):
                    continue
                g = g / sc
                if g.min() < -1e-12:
                    continue
                feas += 1
                A = -g[4] + g[6] + g[10] - g[12] + g[13] + g[18]
                if A < 0:
                    Aneg += 1
                    worst_soc = min(worst_soc, g[6] * g[22] - A * A)
                if nphi < 25:
                    nphi += 1
                    worst_phi = min(worst_phi, min_phi_vals(g, 8, rng))
            print(f"NUMERIC T3 chart: {tried} random chart points (seed 7), {feas} with all gaps >= 0, {Aneg} with A<0; "
                  f"min (g6 g22 - A^2) over A<0 = {worst_soc:.3e}; min Phi over SO(4) at {nphi} feasible points "
                  f"(8 BFGS restarts each) = {worst_phi:.3e}")

        if which in ("T12", "all"):
            for nu, seeds in ((2, range(0, 6)), (3, range(0, 4))):
                T = sorted(L.CONTACT_SETS[nu])
                for seed in seeds:
                    rp = L.random_point_on_contact_stratum(T, seed=seed, margin=1e-2)
                    if rp is None:
                        print(f"NUMERIC T_{nu} seed {seed}: no stratum point found"); continue
                    a, b, al = np.array(rp["a"]), np.array(rp["b"]), rp["alpha"]
                    F = L.F_values_np(a, b, al); g = F.real - F.real[T[0]]
                    mind_a = min(abs(a[i] - a[j]) for i in range(4) for j in range(i))
                    mind_b = min(abs(b[i] - b[j]) for i in range(4) for j in range(i))
                    mph = min_phi_vals(g, 20, np.random.default_rng(seed))
                    extra = f" g7-g9={g[7]-g[9]:.2e}" if nu == 2 else ""
                    print(f"NUMERIC T_{nu} seed {seed}: contact {rp['contact']} min|a_i-a_j|={mind_a:.2e} "
                          f"min|b_i-b_j|={mind_b:.2e}{extra} minPhi={mph:.3e}")

print("RESULT", "FAIL" if FAIL else "ALL EXACT CHECKS PASS", FAIL)
sys.exit(1 if FAIL else 0)
