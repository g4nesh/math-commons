"""Audit replay (batch C01-T1-T3).  Exact replay for final_referee_T2_scalar_compensation_20261022 (scalar SOC inequality)."""
import os, sys, random
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

# square roots made polynomial: a = a1^2 etc., all >= 0
a1, b1, c1, d1, A1, B1, C1, D1 = sp.symbols("a1 b1 c1 d1 A1 B1 C1 D1", nonnegative=True)
a, b, c, d, A, B, C, D = [s**2 for s in (a1, b1, c1, d1, A1, B1, C1, D1)]
s = c + d - C - D; t = b + c - B - C
R = 2 * (a * C + b * D + c * A + d * B) - 4 * a1 * b1 * c1 * d1 - 4 * A1 * B1 * C1 * D1
va, vb, vc, vd = a - A, b - B, c - C, d - D
delta = a1 * c1 - b1 * d1; Delta = A1 * C1 - B1 * D1
check("Thm2_identity_R", sp.expand(R - (2 * delta**2 + 2 * Delta**2 - 2 * (va * vc + vb * vd))) == 0)
# st identity needs sum v = 0: substitute A = a+b+c+d-B-C-D (A1 enters only through A and Delta)
Asub = a + b + c + d - B - C - D
check("Thm2_identity_st_given_equal_mass",
      sp.expand((s * t - (vb * vd - va * vc)).subs(A1**2, Asub)) == 0)
# checkerboard identity F = R + s^2 + t^2 = 2 delta^2 + 2 Delta^2 + (va-vc)^2 - 4 vb vd under equal mass
Fexpr = R + s**2 + t**2 - (2 * delta**2 + 2 * Delta**2 + (va - vc)**2 - 4 * vb * vd)
check("Remark_checkerboard_identity_given_equal_mass", sp.expand(Fexpr.subs(A1**2, Asub)) == 0)
# Lemma 1: ((r+u)(s+v)) - (sqrt(rs)+sqrt(uv))^2 = (sqrt(rv)-sqrt(su))^2
r1, s1, u1, v1 = sp.symbols("r1 s1 u1 v1", nonnegative=True)
check("Lemma1_square_difference", sp.expand((r1**2 + u1**2) * (s1**2 + v1**2) - (r1 * s1 + u1 * v1)**2 - (r1 * v1 - s1 * u1)**2) == 0)
# tau completion
tau = sp.Symbol("tau", positive=True); S_, T_ = sp.symbols("S T", real=True)
check("Thm2_(2)_completion", sp.expand(tau * S_**2 + tau**-1 * T_**2 - 2 * S_ * T_ - (sp.sqrt(tau) * S_ - T_ / sp.sqrt(tau))**2) == 0)
# equality family
rr, mm = sp.symbols("r m", positive=True)
sub = {a1: sp.sqrt(rr), b1: sp.sqrt(rr + mm), c1: sp.sqrt(rr), d1: sp.sqrt(rr + mm),
       A1: sp.sqrt(rr + mm), B1: sp.sqrt(rr), C1: sp.sqrt(rr + mm), D1: sp.sqrt(rr)}
check("Remark_equality_family_s=t=R=0", all(sp.simplify(e.subs(sub)) == 0 for e in (s, t, R)))
# Theorem 3: H - R >= 0 with a=p0^2.. ; check H formula with the (i, i+2 mod 4) pairing equals R-part
p = sp.symbols("p0:4", real=True); q = sp.symbols("q0:4", real=True)
Hq = 2 * sum(p[i]**2 * q[(i + 2) % 4]**2 for i in range(4)) - 4 * p[0] * p[1] * p[2] * p[3] + 4 * q[0] * q[1] * q[2] * q[3]
Rq = 2 * (p[0]**2 * q[2]**2 + p[1]**2 * q[3]**2 + p[2]**2 * q[0]**2 + p[3]**2 * q[1]**2)
check("Thm3_H_minus_signed_products_equals_R_quadratic_part",
      sp.expand(Hq + 4 * p[0] * p[1] * p[2] * p[3] - 4 * q[0] * q[1] * q[2] * q[3] - Rq) == 0)

# ---------- orbit claim: [11,39,17] -> [9,1,22] ----------
V = L.F_VALUES
idx = {tuple(v): j for j, v in enumerate(V)}
def img(g, j):
    return idx[tuple(L.g_act_values(g, list(V[j])))]
src, dst = (11, 39, 17), (9, 1, 22)
hits_ordered = [g for g in L.G_ELEMENTS if tuple(img(g, j) for j in src) == dst]
hits_set = [g for g in L.G_ELEMENTS if {img(g, j) for j in src} == set(dst) and img(g, 17) == 22]
print("INFO group elements mapping (11,39,17)->(9,1,22) ordered:", len(hits_ordered), "as sets with h->h:", len(hits_set))
rowp = (0, 2, 3, 1); colp = (3, 1, 2, 0)
inv = L.inverse
cands = {}
for nr, R_ in (("row", rowp), ("row^-1", inv(rowp))):
    for nc, C_ in (("col", colp), ("col^-1", inv(colp))):
        g = (C_, R_, 0)
        cands[(nr, nc)] = tuple(img(g, j) for j in src)
print("INFO images under the stated row/column permutations in the 4 conventions:", cands)
check("orbit_claim_some_convention_maps_(11,39,17)_to_(9,1,22)", any(v == dst for v in cands.values()))
# the T_3 note needs the triple (1,11,17); is it in the same orbit, with h -> h?
src2 = (1, 11, 17)
hit2 = [g for g in L.G_ELEMENTS if {img(g, j) for j in src2} == set(dst) and img(g, 17) == 22]
check("T3_triple_(1,11,17)_in_same_orbit_as_(9,1,22)_with_h_to_h", len(hit2) > 0, f"{len(hit2)} elements")

if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    for which in ("scalar", "soc"):

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
            for trip, seed in (((9, 1, 22), 12),):
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
