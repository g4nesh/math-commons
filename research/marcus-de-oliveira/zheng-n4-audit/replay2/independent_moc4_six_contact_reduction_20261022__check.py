"""Second-pass independent replay for the supplement document
independent_moc4_six_contact_reduction_20261022 (8 pages).

Usage:  MOC4_SUPP=<supplement root> python -I -B this_file.py [--no-numeric]
Only the supplement's facet JSON (root_permutation_polytope_complete_20261022.json) is
loaded, as data, to translate the document's facet numbering f_j into profiles.
Nothing from the supplement is executed or imported.

EXACT checks print 'EXACT ... PASS/FAIL'; NUMERIC checks print 'NUMERIC ...'.
Exit code is nonzero iff some EXACT check fails.
"""
import itertools
import json
import os
import random
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"EXACT {name}: {'PASS' if ok else 'FAIL'} {info}".rstrip())
    if not ok:
        FAIL.append(name)


PERMS = L.PERMS
IDX = {p: i for i, p in enumerate(PERMS)}

# ---------------------------------------------------------------------------------
# (M1)-(M4)
# ---------------------------------------------------------------------------------
check("M1_P_rank14", L.rank([list(map(Fraction, r)) for r in L.P]) == 14)
check("M2_P_kappa_is_ones", L.values(L.KAPPA) == [1] * 24)

# P m = lambda prod(a_i + b_sigma(i)) symbolically in 8 dehomogenized variables
a_s = sp.symbols("a0:4")
b_s = sp.symbols("b0:4")
m_sym = L.matching_vector(list(a_s), list(b_s))
ok = True
for k, s in enumerate(PERMS):
    lhs = sp.expand(sum(L.P[k][r] * m_sym[r] for r in range(14)))
    rhs = sp.expand(sp.Mul(*[a_s[i] + b_s[s[i]] for i in range(4)]))
    ok &= sp.expand(lhs - rhs) == 0
check("M1_P_rows_are_sigma_products_symbolic", ok)

ms, Cexpr = L.cubic_sympy()
grad = [sp.diff(Cexpr, v) for v in ms]
kap = {v: (1 if i == 13 else 0) for i, v in enumerate(ms)}
check("M2_gradC_kappa_zero", all(sp.expand(g.subs(kap)) == 0 for g in grad))
H = sp.Matrix([[sp.diff(Cexpr, u, v) for v in ms] for u in ms]).subs(kap)
check("M2_Hkappa_matches_paper_3.10", H == sp.Matrix(L.hess_kappa()))
t = sp.Symbol("t")
shift = {v: v + t * (1 if i == 13 else 0) for i, v in enumerate(ms)}
xvec = sp.Matrix(ms)
ok = all(sp.expand(grad[i].subs(shift, simultaneous=True) - grad[i] - t * (H * xvec)[i]) == 0
         for i in range(14))
check("M2_eq2_polarization_identity", ok)

# X subset Sing(C): grad C vanishes on the symbolic matching vector (dehomogenized,
# 8 variables; X is the closure of the cone over these, and the identity is homogeneous)
gens = list(a_s) + list(b_s)
mpoly = [sp.Poly(x, *gens) for x in m_sym]
sub = dict(zip(ms, mpoly))
ok = True
for g in grad:
    pg = sp.Poly(g, *ms)
    tot = sp.Poly(0, *gens)
    for mon, coeff in pg.terms():
        term = sp.Poly(coeff, *gens)
        for v, e in zip(ms, mon):
            for _ in range(e):
                term = term * sub[v]
        tot = tot + term
    ok &= tot.is_zero
check("M2_X_in_Sing_C_symbolic_8vars", ok)

# (M3): singleton response vectors e_ij, kernel of H_kappa
E = [L.solve_P([1 if s[i] == j else 0 for s in PERMS]) for i in range(4) for j in range(4)]
check("M3_e_ij_span_rank10", L.rank(E) == 10)
Hk = [[Fraction(x) for x in r] for r in L.hess_kappa()]
check("M3_Hkappa_rank4", L.rank(Hk) == 4)
check("M3_Hkappa_kills_e_ij", all(all(x == 0 for x in L.matvec(Hk, e)) for e in E))
# kernel dimension 10 and span(e_ij) of dim 10 inside it => equal
check("M3_ker_Hkappa_equals_span_e_ij", len(L.nullspace(Hk)) == 10)

# (M4): Phi_O(kappa) = 1, Phi_O(e_ij) = O_ij^2 on SO(4); the det -1 case via O*diag(1,1,1,-1)
rng = random.Random(20261125)
R = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]
okk = oke = okd = True
for trial in range(8):
    O = L.random_rational_SO4(rng)
    okk &= L.Phi(L.KAPPA, O) == 1
    oke &= all(L.Phi(E[4 * i + j], O) == O[i][j] ** 2 for i in range(4) for j in range(4))
    Om = L.matmul(O, R)  # det -1
    a = [L.QI(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
    b = [L.QI(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(4)]
    lam = L.QI(rng.randint(-3, 3), rng.randint(1, 3))
    m = L.matching_vector(a, b, lam)
    F = [sum((L.P[k][r] * m[r] for r in range(14)), L.QI(0)) for k in range(24)]
    for OO in (O, Om):
        Ophi = O  # functional Phi_O := Leibniz sum at O (det 1) or at Om*R = O (det -1)
        lhs = lam * L.det_Da_plus_ODbOt(a, b, OO)
        rhs = sum((F[k] * w for k, w in enumerate(L.leibniz_weights(Ophi))), L.QI(0))
        okd &= lhs == rhs
check("M4_Phi_kappa_1", okk)
check("M4_Phi_e_ij_is_Oij_squared", oke)
check("M4_Phi_gives_lambda_det_for_det_pm1", okd, "(8 Gaussian-rational spectra, both det signs)")

# Lemma 4 / Theorem 5: common kernel of gap forms = C kappa; h_{O,tau} in their real span
for tau in (0, 7, 23):
    rows = [[Fraction(L.P[s][r] - L.P[tau][r]) for r in range(14)] for s in range(24)]
    ns = L.nullspace(rows)
    ok = len(ns) == 1 and all(ns[0][r] == 0 for r in range(13))
    check(f"Thm5_gap_forms_common_kernel_is_kappa_tau{tau}", ok)
    O = L.random_rational_SO4(rng)
    w = L.leibniz_weights(O)
    hrow = [sum(Fraction(w[s]) * L.P[s][r] for s in range(24)) - L.P[tau][r] for r in range(14)]
    check(f"Thm5_h_in_real_span_of_gap_forms_tau{tau}", L.rank(rows + [hrow]) == L.rank(rows))

# ---------------------------------------------------------------------------------
# Section 5: facet profiles and Proposition 8 inventory
# ---------------------------------------------------------------------------------
FV = [tuple(v) for v in L.F_VALUES]
nonminor = [j for j in range(82) if not L.F_IS_MINOR[j]]
check("Sec5_82_facets_34_minor_48_nonminor", len(FV) == 82 and len(nonminor) == 48)
check("Sec5_nonminor_profile_13zero_10one_1two",
      all(sorted(FV[j]).count(0) == 13 and sorted(FV[j]).count(1) == 10 and sorted(FV[j]).count(2) == 1
          for j in nonminor))

orbits = L.enumerate_six_set_orbits()
check("Prop8_214_orbits", len(orbits) == 214 and sum(s for _, s in orbits) == 134596)


def face_vertices(T):
    gs = L.generators(T)
    if not gs:
        return frozenset(range(24))
    out = set(range(24))
    for j in gs:
        out &= L.F_ZERO[j]
    return frozenset(out)


nm_counts = {}
n90 = 0
n_six_face = 0
all90_six = True
empty_gen = 0
for T, _ in orbits:
    gs = L.generators(T)
    k = sum(1 for j in gs if not L.F_IS_MINOR[j])
    if not gs:
        empty_gen += 1
    sixface = face_vertices(T) == frozenset(T)
    n_six_face += sixface
    if k > 0:
        n90 += 1
        nm_counts[k] = nm_counts.get(k, 0) + 1
        all90_six &= sixface
check("Prop8_90_orbits_with_nonminor_generator", n90 == 90)
check("Prop8_nonminor_count_distribution_51_35_3_1",
      nm_counts == {1: 51, 2: 35, 3: 3, 4: 1}, str(dict(sorted(nm_counts.items()))))
check("Prop8_113_orbits_are_exact_six_vertex_faces", n_six_face == 113, f"(found {n_six_face})")
check("Prop8_all_90_reps_exact_six_vertex_faces", all90_six)
print(f"INFO Prop8 orbits with empty generator set: {empty_gen}")
# "retained orbits containing the specified four-contact type" (Sect. 6 names {2,7,19,20})
F4 = frozenset({2, 7, 19, 20})
G_MAPS = [L.g_index_map(g) for g in L.G_ELEMENTS]
orb4 = {frozenset(tm[j] for j in F4) for tm in G_MAPS}
cnt4 = 0
for T, _ in orbits:
    gs = L.generators(T)
    if any(not L.F_IS_MINOR[j] for j in gs):
        if any(S <= frozenset(T) for S in orb4) or any(
                frozenset(tm[j] for j in T) >= F4 for tm in G_MAPS):
            cnt4 += 1
check("Prop8_9_retained_orbits_contain_type_{2,7,19,20}", cnt4 == 9, f"(found {cnt4})")
check("Prop8_paper_T1..T90_in_distinct_retained_orbits",
      len({L.orbit_rep(L.CONTACT_SETS[k]) for k in range(1, 91)}) == 90)

# ---------------------------------------------------------------------------------
# Proposition 10 (paper T_1) and Proposition 13 (paper T_2)
# ---------------------------------------------------------------------------------
T1doc = frozenset({0, 1, 2, 3, 6, 7})
T2doc = frozenset({0, 1, 2, 3, 6, 8})
check("Prop10_set_is_paper_T1", T1doc == L.CONTACT_SETS[1])
check("Prop13_set_is_paper_T2", T2doc == L.CONTACT_SETS[2])

# entry faces: all 16 have 15 generators, all minor, 6 singleton + 9 two-by-two
ok = True
for i in range(4):
    for j in range(4):
        Ef = frozenset(k for k, s in enumerate(PERMS) if s[i] == j)
        gs = L.generators(Ef)
        names = [L.F_TABLE[g] for g in gs]
        sing = sum(1 for g in gs if L.F_IS_MINOR[g] and len(L.F_TABLE[g][1]) == 1)
        ok &= len(gs) == 15 and all(L.F_IS_MINOR[g] for g in gs) and sing == 6
check("Prop10_entry_faces_15_minor_generators_6_singleton_9_twobytwo", ok)

A = sp.symbols("A0:4")
Cc = sp.symbols("C0:4")
lam = sp.Symbol("lam")


def Fs(k):
    s = PERMS[k]
    return lam * sp.Mul(*[A[i] - Cc[s[i]] for i in range(4)])


def D(r, s):
    return Fs(r) - Fs(s)


fac = {
    (0, 6): lam * (A[0] - A[1]) * (Cc[0] - Cc[1]) * (A[2] - Cc[2]) * (A[3] - Cc[3]),
    (1, 7): lam * (A[0] - A[1]) * (Cc[0] - Cc[1]) * (A[2] - Cc[3]) * (A[3] - Cc[2]),
    (0, 2): lam * (A[0] - Cc[0]) * (A[3] - Cc[3]) * (A[1] - A[2]) * (Cc[1] - Cc[2]),
    (1, 3): lam * (A[0] - Cc[0]) * (A[2] - Cc[3]) * (A[1] - A[3]) * (Cc[1] - Cc[2]),
    (0, 1): lam * (A[0] - Cc[0]) * (A[1] - Cc[1]) * (A[2] - A[3]) * (Cc[2] - Cc[3]),
    (2, 3): lam * (A[0] - Cc[0]) * (A[1] - Cc[2]) * (A[2] - A[3]) * (Cc[1] - Cc[3]),
    (6, 8): lam * (A[0] - Cc[1]) * (A[3] - Cc[3]) * (A[1] - A[2]) * (Cc[0] - Cc[2]),
    (7, 9): None,
}
for (r, s), f in fac.items():
    if f is not None:
        check(f"Prop10_13_factorization_D{r}{s}", sp.expand(D(r, s) - f) == 0)
# the three ratios as displayed
rat = [((0, 6), (1, 7), (A[2] - Cc[2]) * (A[3] - Cc[3]) / ((A[2] - Cc[3]) * (A[3] - Cc[2]))),
       ((0, 2), (1, 3), (A[3] - Cc[3]) * (A[1] - A[2]) / ((A[2] - Cc[3]) * (A[1] - A[3]))),
       ((0, 1), (2, 3), (A[1] - Cc[1]) * (Cc[2] - Cc[3]) / ((A[1] - Cc[2]) * (Cc[1] - Cc[3])))]
for (p, q, r_) in rat:
    check(f"Prop10_ratio_D{p[0]}{p[1]}_over_D{q[0]}{q[1]}", sp.simplify(fac[p] / fac[q] - r_) == 0)
# imaginary cross-multiplied numerators
x, y, tt, u, v, pp, qq = sp.symbols("x y t u v p q", real=True)
A1 = u + sp.I * v
num2 = sp.expand(sp.im(sp.expand(y * (A1 - x) * sp.conjugate(x * (A1 - y)))))
check("Prop10_second_ratio_imag_numerator_xy(x-y)v", sp.expand(num2 - x * y * (x - y) * v) == 0)
beta = pp + sp.I * qq
num3 = sp.expand(sp.im(sp.expand((u - beta) * tt * sp.conjugate((u - tt) * beta))))
check("Prop10_third_ratio_imag_numerator_-tu(u-t)q", sp.expand(num3 + tt * u * (u - tt) * qq) == 0)
# exceptional collision u = 0 (A1 = C3 = 0)
Lsym = sp.Symbol("L")
subs0 = {A[1]: 0, A[2]: x, A[3]: y, Cc[1]: beta, Cc[2]: tt, Cc[3]: 0}
d01 = (fac[(0, 1)] / (lam * (A[0] - Cc[0]))).subs(subs0) * Lsym
d02 = (fac[(0, 2)] / (lam * (A[0] - Cc[0]))).subs(subs0) * Lsym
check("Prop10_collision_D01=-L*beta*t*(x-y)", sp.expand(d01 + Lsym * beta * tt * (x - y)) == 0)
check("Prop10_collision_D02=-L*x*y*(beta-t)", sp.expand(d02 + Lsym * x * y * (beta - tt)) == 0)
# with A1..A3, C1..C3 real and L = i*l, every product with sigma(0)=0 is L * real
check("Prop10_sigma0_fixed_products_are_L_times_(A0-C0)_free",
      all(sp.expand(Fs(k) / (lam * (A[0] - Cc[0]))).has(A[0]) is False for k in range(6)))

# identity (6) D79*D02 = D68*D13
check("Prop13_identity_6_D79D02=D68D13", sp.expand(D(7, 9) * D(0, 2) - D(6, 8) * D(1, 3)) == 0)

# facet numbering of the document = supplement JSON (data only)
supp = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else None)
if supp:
    jp = os.path.join(supp, "proofs", "manin_transfer_20261021", "root_permutation_polytope_complete_20261022.json")
    dat = json.load(open(jp))
    # verify JSON permutation order equals ours
    check("JSON_perm_order_lexicographic", [tuple(x - 1 for x in p) for p in dat["permutations_1_based"]] == PERMS)
    JV = [tuple(f["values"]) for f in dat["facets"]]
    check("JSON_facet_profiles_equal_lib_set", len(JV) == 82 and set(JV) == set(FV))
    lib_of = {FV[j]: j for j in range(82)}
    Jdoc = [j for j in range(82) if all(JV[j][s] == 0 for s in T2doc)]
    check("Prop13_J_equals_doc_list", Jdoc == [0, 1, 2, 3, 4, 5, 6, 9, 10, 15, 16, 17], str(Jdoc))
    check("Prop13_J_matches_lib_generators_T2",
          sorted(lib_of[JV[j]] for j in Jdoc) == sorted(L.generators(L.CONTACT_SETS[2])))
    minor = {j: L.F_IS_MINOR[lib_of[JV[j]]] for j in Jdoc}
    check("Prop13_only_f17_nonminor", [j for j in Jdoc if not minor[j]] == [17])
    f = {j: JV[j] for j in Jdoc}
    check("Prop13_f17+f9=f4+f10+f16",
          all(f[17][s] + f[9][s] == f[4][s] + f[10][s] + f[16][s] for s in range(24)))
    ell = {j: f[j][7] - f[j][9] for j in Jdoc}
    check("Prop13_ell_values", ell[6] == -1 and ell[9] == 1 and ell[17] == -1 and
          all(ell[j] == 0 for j in Jdoc if j not in (6, 9, 17)), str(ell))
    # formal coefficient transfer: g = sum a_j f_j with a9 = a6 + a17 equals sum beta_j f_j
    asy = {j: sp.Symbol(f"a{j}", nonnegative=True) for j in Jdoc}
    asy9 = asy[6] + asy[17]
    g = [sum((asy9 if j == 9 else asy[j]) * f[j][s] for j in Jdoc) for s in range(24)]
    beta_ = {j: asy[j] for j in Jdoc if j != 17}
    beta_[9] = asy[6]
    beta_[4] = asy[4] + asy[17]
    beta_[10] = asy[10] + asy[17]
    beta_[16] = asy[16] + asy[17]
    g2 = [sum(beta_[j] * f[j][s] for j in beta_) for s in range(24)]
    check("Prop13_formal_coefficient_transfer", all(sp.expand(g[s] - g2[s]) == 0 for s in range(24)))
    # the constraint a9 = a6 + a17 is exactly ell(g)=0
    gen_g = [sum(asy[j] * f[j][s] for j in Jdoc) for s in range(24)]
    check("Prop13_ell(g)=0_iff_a9=a6+a17", sp.expand(gen_g[7] - gen_g[9] - (asy[9] - asy[6] - asy[17])) == 0)
else:
    print("SKIP JSON-based checks: set MOC4_SUPP")

# ---------------------------------------------------------------------------------
# NUMERIC falsification
# ---------------------------------------------------------------------------------
if "--no-numeric" not in sys.argv:
    import numpy as np
    for name, T in (("T1", T1doc), ("T2", T2doc)):
        worst = np.inf
        found = 0
        extra = []
        degen = []
        g79 = []
        for seed in range(12):
            pt = L.random_point_on_contact_stratum(T, seed=seed, tries=15, tol=1e-9)
            if pt is None:
                continue
            found += 1
            a, b, al = pt["a"], pt["b"], pt["alpha"]
            mind = min(min(abs(a[i] - a[j]), abs(b[i] - b[j])) for i in range(4) for j in range(i + 1, 4))
            minprod = min(abs(a[i] + b[j]) for i in range(4) for j in range(4))
            degen.append((mind, minprod))
            gaps = pt["gaps"]
            extra.append(sorted(set(pt["contact"]) - set(T)))
            g79.append(abs(gaps[7] - gaps[9]))
            m = L.matching_vector_np(a, b, al)
            c = np.real(m).astype(float)
            Fps = np.real(L.F_values_np(a, b, al))[pt["pi_star"]]
            c[13] -= Fps
            val, _ = L.min_phi_over_SO4(c, restarts=12, seed=seed)
            worst = min(worst, val)
        print(f"NUMERIC {name}: {found}/12 stratum points (seeds 0..11, tries 15, tol 1e-9); "
              f"worst min Phi over SO(4) (12 BFGS restarts) = {worst:.3e}")
        print(f"NUMERIC {name}: extra contacts per point = {extra}")
        print(f"NUMERIC {name}: (min within-spectrum distance, min |a_i+b_j|) per point = "
              f"{[(round(float(d), 4), round(float(p), 4)) for d, p in degen]}")
        if name == "T2":
            print(f"NUMERIC T2: |g7-g9| per point = {[round(float(z), 6) for z in g79]}")
    # Prop 13 prediction g7 = g9 on NONDEGENERATE T2 points (distinct within each spectrum,
    # all six contact products nonzero); degenerate samples are filtered out.
    nd = []
    for seed in range(12, 60):
        pt = L.random_point_on_contact_stratum(T2doc, seed=seed, tries=4, tol=1e-9)
        if pt is None:
            continue
        a, b, al = pt["a"], pt["b"], pt["alpha"]
        mind = min(min(abs(a[i] - a[j]), abs(b[i] - b[j])) for i in range(4) for j in range(i + 1, 4))
        Fv = L.F_values_np(a, b, al)
        minF = min(abs(Fv[k]) for k in T2doc)
        if mind > 1e-4 and minF > 1e-4:
            nd.append((seed, round(float(mind), 4), round(float(minF), 4),
                       float(abs(pt["gaps"][7] - pt["gaps"][9]))))
    print(f"NUMERIC T2 nondegenerate points (seeds 12..59, tries 4): (seed, min dist, min |F_T|, |g7-g9|) = {nd}")
    # Prop 10 structural prediction: with margin > 0 (contact set exactly T1) a point must be
    # degenerate (repeated eigenvalue within a spectrum).
    res = []
    for seed in range(6):
        pt = L.random_point_on_contact_stratum(T1doc, seed=seed, tries=15, tol=1e-9, margin=1e-2)
        if pt is None:
            res.append(None)
            continue
        a, b = pt["a"], pt["b"]
        mind = min(min(abs(a[i] - a[j]), abs(b[i] - b[j])) for i in range(4) for j in range(i + 1, 4))
        res.append(round(float(mind), 6))
    print(f"NUMERIC T1 exact-contact (margin 1e-2) points: min within-spectrum distance per seed = {res}")

print(f"SUMMARY: {len(FAIL)} EXACT failures {FAIL}")
sys.exit(1 if FAIL else 0)
