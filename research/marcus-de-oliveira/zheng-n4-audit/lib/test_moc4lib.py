"""Exact tests for moc4lib.  Run from this directory:  python -I -B test_moc4lib.py
Prints one line per check; exit code 0 iff every check passes.
Every check is exact (integers, Fractions, Gaussian rationals, sympy) unless its name
starts with numeric_.  Lines starting with INFO are informational, not checks."""
import itertools
import os
import random
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moc4lib as L  # noqa: E402

RESULTS = []
T0 = time.time()


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(f"{'PASS' if cond else 'FAIL'} {name}" + (f" : {detail}" if detail else ""), flush=True)


def info(msg):
    print("INFO " + msg, flush=True)


# ---------------------------------------------------------------- permutations
def test_perms():
    check("perms_24_lexicographic", len(L.PERMS) == 24 and L.PERMS == sorted(L.PERMS)
          and L.PERMS[0] == (0, 1, 2, 3) and L.PERMS[23] == (3, 2, 1, 0))
    ok = all(L.compose(p, L.inverse(p)) == L.IDENTITY for p in L.PERMS)
    ok &= all(L.sgn(L.compose(p, q)) == L.sgn(p) * L.sgn(q) for p in L.PERMS for q in L.PERMS)
    ok &= L.sgn((1, 0, 2, 3)) == -1 and sum(L.SGN) == 0
    ok &= all(L.perm_index(L.perm_from_oneline(L.perm_to_oneline(p))) == i
              for i, p in enumerate(L.PERMS))
    check("perms_sgn_compose_inverse", ok)


# ---------------------------------------------------------------- matching basis and P
def test_basis_and_P():
    # all perfect matchings of 8 points, noncrossing ones
    def matchings(pts):
        if not pts:
            yield []
            return
        a = pts[0]
        for k in range(1, len(pts)):
            rest = pts[1:k] + pts[k + 1:]
            for m in matchings(rest):
                yield [(a, pts[k])] + m
    allm = list(matchings(list(range(1, 9))))
    nc = {frozenset(m) for m in allm if L.is_noncrossing(m)}
    listed = {frozenset(m) for m in L.MATCHINGS}
    check("A1_basis_is_the_14_noncrossing_matchings", len(allm) == 105 and len(nc) == 14
          and nc == listed and len(L.MATCHINGS) == 14)
    masks = sorted({m for p in L.B_POLYS for m in p})
    check("A1_B_linearly_independent", L.rank([[p.get(m, 0) for m in masks] for p in L.B_POLYS]) == 14)
    # every [ij][kl] product of any perfect matching lies in span (dim R_1 = 14)
    ok = True
    for m in allm:
        try:
            L.express_in_basis(L.matching_poly(m))
        except ValueError:
            ok = False
    check("R1_all_105_matching_products_in_span_of_B", ok)
    # P identity as a polynomial identity in eight symbolic 2-vectors (16 variables)
    ok = True
    for s, row in zip(L.PERMS, L.P):
        lhs = L.p_sigma_poly(s)
        rhs = {}
        for a in range(14):
            rhs = L.poly_add(rhs, L.B_POLYS[a], row[a])
        ok &= lhs == rhs
    check("P_identity_prod_bracket_eq_sum_P_B_16vars", ok, "all 24 rows, coefficientwise")
    # independent check of (4.5) with sympy, w_i = (z_i, 1)
    import sympy
    z = sympy.symbols("z1:9")
    br = lambda i, j: z[i - 1] - z[j - 1]          # det((z_i,1),(z_j,1))
    Bz = [sympy.expand(sympy.Mul(*[br(i, j) for i, j in m])) for m in L.MATCHINGS]
    ok = True
    for s, row in zip(L.PERMS, L.P):
        lhs = sympy.expand(sympy.Mul(*[(z[4 + s[i]] - z[i]) for i in range(4)]))
        ok &= sympy.expand(lhs - sum(row[a] * Bz[a] for a in range(14))) == 0
    check("eq4.5_sympy_prod_(z_{4+s(i)}-z_i)_eq_sum_P_B(w)", ok)
    check("P_integer_entries", all(isinstance(x, int) for r in L.P for x in r),
          f"entry range [{min(min(r) for r in L.P)}, {max(max(r) for r in L.P)}]")
    check("P_rank_14", L.rank(L.P) == 14)
    d = L.det([L.P[r] for r in L.P_PIVOT_ROWS])
    check("P_minor_rows_0,1,2,3,4,6,7,8,9,10,12,13,16,18_det_is_1", d == 1, f"det = {d}")
    check("P_kappa_column13_all_ones", all(r[13] == 1 for r in L.P)
          and L.values(L.KAPPA) == [1] * 24)


# ---------------------------------------------------------------- the cubic
def hpoly_from_mask(poly):
    """R_1 element -> dict over exponent tuples of x1..x8,y1..y8."""
    out = {}
    for mask, c in poly.items():
        e = tuple([1 if mask >> l & 1 else 0 for l in range(8)] +
                  [0 if mask >> l & 1 else 1 for l in range(8)])
        out[e] = out.get(e, 0) + c
    return out


def hmul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            e = tuple(x + y for x, y in zip(a, b))
            out[e] = out.get(e, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def test_cubic():
    import sympy
    ms, C = L.cubic_sympy()
    grads = [sympy.expand(sympy.diff(C, m)) for m in ms]
    # (a) partials vanish identically on matching vectors of 8 symbolic 2-vectors
    HB = [hpoly_from_mask(p) for p in L.B_POLYS]
    ok = True
    for k, gk in enumerate(grads):
        tot = {}
        for mon, coeff in sympy.Poly(gk, *ms).terms():
            idx = [i for i, e in enumerate(mon) for _ in range(e)]
            assert len(idx) == 2
            prod = hmul(HB[idx[0]], HB[idx[1]])
            for e, v in prod.items():
                tot[e] = tot.get(e, 0) + int(coeff) * v
        ok &= all(v == 0 for v in tot.values())
    check("A2a_grad_C_vanishes_on_X_16vars", ok, "X subset Sing(C), all 14 partials")
    # (b) adjacent transpositions
    Ts = [L.transposition_matrix(i) for i in range(1, 8)]
    I14 = [[int(i == j) for j in range(14)] for i in range(14)]
    ok_sq = all(L.matmul(T, T) == I14 for T in Ts)
    ok_sign = True
    for T in Ts:
        Tm = [sum(T[r][s] * ms[s] for s in range(14)) for r in range(14)]
        ok_sign &= sympy.expand(L.cubic_C(Tm) + C) == 0
    check("A2b_T_i_squared_is_identity_i=1..7", ok_sq)
    check("A2b_C(T_i m)_eq_-C(m)_i=1..7", ok_sign)
    # (c) independence of partials via the listed monomials
    M = []
    for gk in grads:
        pk = sympy.Poly(gk, *ms)
        row = []
        for (i, j) in L.C_PARTIAL_MONOMIALS:
            e = [0] * 14; e[i] += 1; e[j] += 1
            row.append(int(pk.coeff_monomial(tuple(e))))
        M.append(row)
    d = L.det(M)
    check("A2c_partials_coefficient_minor_det_is_1", d == 1,
          f"det = {d} (rows dC/dm0..dC/dm13, columns in the paper's monomial order)")
    # (d) Hessian at kappa
    kap = {m: (1 if m == ms[13] else 0) for m in ms}
    H = [[int(sympy.diff(C, ms[i], ms[j]).subs(kap)) for j in range(14)] for i in range(14)]
    check("eq3.10_Hess_C(kappa)_eq_e0e6+e6e0-e1e5-e5e1", H == L.hess_kappa())
    check("eq3.10_rank_H_kappa_is_4", L.rank(H) == 4)
    ker = L.nullspace(H)
    span_ok = L.rank(ker) == 10 and all(v[i] == 0 for v in ker for i in (0, 1, 5, 6))
    check("eq3.11_ker_H_kappa_eq_{m0=m1=m5=m6=0}", span_ok)
    # E = ker H : coordinates of E(Q) span exactly this 10-dim space
    Es = []
    for i in range(4):
        for j in range(4):
            Q = [[int(a == i and b == j) for b in range(4)] for a in range(4)]
            Es.append(L.E_coords(Q))
    check("eq3.11_E_eq_ker_H_kappa", L.rank(Es) == 10 and L.rank(Es + ker) == 10)
    # (e) grad C(kappa) = 0 and polarization (3.12)
    check("eq3.12_grad_C(kappa)_is_0", all(g.subs(kap) == 0 for g in grads))
    t = sympy.Symbol("t")
    shifted = {m: (m + t if m == ms[13] else m) for m in ms}
    Hx = [sum(L.hess_kappa()[i][j] * ms[j] for j in range(14)) for i in range(14)]
    ok = all(sympy.expand(grads[i].xreplace(shifted) - grads[i] - t * Hx[i]) == 0 for i in range(14))
    check("eq3.12_grad_C(x+t kappa)_eq_grad_C(x)+t H x", ok)
    # (f) dim X = 6: Jacobian of (lambda, z1..z8) -> lambda * B(z)
    rng = random.Random(11)
    ranks = []
    for _ in range(3):
        z = [Fraction(rng.randint(-30, 30), rng.randint(1, 9)) for _ in range(8)]
        lam = Fraction(rng.randint(1, 9), rng.randint(1, 9))
        rows = [[L.eval_poly(L.B_POLYS[r], z) for r in range(14)]]
        for k in range(8):
            der = []
            for r in range(14):
                dp = {m & ~(1 << k): c for m, c in L.B_POLYS[r].items() if m >> k & 1}
                der.append(lam * L.eval_poly(dp, z))
            rows.append(der)
        ranks.append(L.rank(rows))
    check("dimX_jacobian_rank_6_at_rational_points", ranks == [6, 6, 6], f"ranks {ranks}")
    # kappa in X: first four columns (0,1), last four (1,1)
    mk = L.matching_vector_from_points([(0, 1)] * 4 + [(1, 1)] * 4)
    check("kappa_in_X_matching_vector_of_(0,1)^4(1,1)^4", mk == L.KAPPA)


# ---------------------------------------------------------------- entry subspace
def test_entry():
    import sympy
    q = sympy.symbols("q11:15 q21:25 q31:35 q41:45")
    Q = [[q[4 * i + j] for j in range(4)] for i in range(4)]
    E = L.E_coords(Q)
    vals = [sympy.expand(sum(L.P[s][a] * E[a] for a in range(14))) for s in range(24)]
    ok = all(sympy.expand(vals[k] - sum(Q[i][s[i]] for i in range(4))) == 0
             for k, s in enumerate(L.PERMS))
    check("eq3.7_A2_(P E(Q))_sigma_eq_sum_q_i,sigma(i)_symbolic", ok)
    J = [[sympy.diff(e, v) for v in q] for e in E]
    J = [[int(x) for x in r] for r in J]
    check("entry_subspace_dim_10", L.rank(J) == 10)


# ---------------------------------------------------------------- minor vectors and facets
def test_minors_and_facets():
    singles = {}
    for i in range(1, 5):
        for j in range(1, 5):
            singles[(i, j)] = tuple(L.n_vec((i,), (j,)))     # solve_P checks all 24 rows
    pairs = {}
    comp_ok = True
    for I in itertools.combinations(range(1, 5), 2):
        for J in itertools.combinations(range(1, 5), 2):
            v = tuple(L.n_vec(I, J))
            pairs[(I, J)] = v
            Ic = tuple(sorted(set(range(1, 5)) - set(I)))
            Jc = tuple(sorted(set(range(1, 5)) - set(J)))
            comp_ok &= v == tuple(L.n_vec(Ic, Jc))
    check("eq4.1_minor_vectors_exist_all_24_rows", True, "52 index pairs solved exactly")
    check("n_IJ_eq_n_IcJc_for_2sets", comp_ok)
    allv = set(singles.values()) | set(pairs.values())
    check("minor_vectors_16+18_distinct", len(set(singles.values())) == 16
          and len(set(pairs.values())) == 18 and len(allv) == 34)
    z1 = {sum(1 for x in L.values(list(v)) if x == 0) for v in singles.values()}
    z2 = {sum(1 for x in L.values(list(v)) if x == 0) for v in pairs.values()}
    check("minor_vectors_zero_counts_18_and_20", z1 == {18} and z2 == {20})
    # t, d, h
    ok_t = ok_d = ok_h = True
    hs = set()
    for pi in L.PERMS:
        tv, dv = L.values(L.t_vec(pi)), L.values(L.d_vec(pi))
        for k, s in enumerate(L.PERMS):
            agree = sum(1 for i in range(4) if s[i] == pi[i])
            ok_t &= tv[k] == agree
            ok_d &= (dv[k] != 0) == (agree == 1)
        ok_d &= sorted(dv).count(1) == 4 and sorted(dv).count(-1) == 4
        for eps in (1, -1):
            hv = L.values(L.h_vec(pi, eps))
            ok_h &= sorted(hv) == [0] * 13 + [1] * 10 + [2] and hv[L.perm_index(pi)] == 2
            hs.add(tuple(hv))
    check("t_pi_values_count_agreements", ok_t)
    check("d_pi_nonzero_exactly_at_8_relative_3cycles_4plus_4minus", ok_d)
    check("h_profile_13x0_10x1_1x2_and_max_at_pi", ok_h)
    check("h_48_distinct", len(hs) == 48)
    # Appendix D list
    mset = {tuple(L.values(list(v))) for v in allv}
    fv = [tuple(v) for v in L.F_VALUES]
    check("D_f_list_is_bijection_onto_34_minor+48_h",
          len(set(fv)) == 82 and set(fv[j] for j in range(82) if L.F_IS_MINOR[j]) == mset
          and set(fv[j] for j in range(82) if not L.F_IS_MINOR[j]) == hs)
    prof_ok = True
    for j in range(82):
        cnt = sorted(fv[j])
        k, a, b = L.F_TABLE[j]
        if k == "n":
            prof_ok &= cnt == [0] * (18 if len(a) == 1 else 20) + [1] * (6 if len(a) == 1 else 4)
        else:
            prof_ok &= cnt == [0] * 13 + [1] * 10 + [2]
    check("D_f_value_profiles", prof_ok)
    check("D_f_nonnegative_at_all_vertices", all(x >= 0 for v in fv for x in v))
    check("D_zero_sets_pairwise_distinct", len(set(L.F_ZERO)) == 82)
    rk = [L.rank([L.P[s] for s in sorted(L.F_ZERO[j])]) for j in range(82)]
    check("D_each_zero_set_affine_rank_12_(facet)", all(r == 13 for r in rk),
          "linear rank 13 of zero-set vertices in R^14 = affine dim 12")
    # pairwise distinct as hyperplanes is implied by distinct zero sets + facet property


# ---------------------------------------------------------------- Prop 4.1 completeness (exact certificate)
def test_prop41_completeness():
    """qhull (floating point) proposes a triangulation of the boundary; every claim used
    is then checked exactly: (i) each 13-vertex simplex lies in a listed facet zero set;
    (ii) the chain of simplices is a mod-2 cycle (each 12-vertex face occurs an even
    number of times); (iii) a rational point p in the relative interior of one facet F is
    interior to an odd number of the simplices in F and on the boundary of none.  Then
    [z] != 0 in H_12(dQ; Z/2) = Z/2, so supp z = dQ (else [z] would come from a disc),
    so dQ is covered by listed facets and there is no unlisted facet."""
    try:
        import numpy as np
        from scipy.spatial import ConvexHull
    except ImportError:
        check("prop4.1_completeness_exact_certificate", False, "scipy missing")
        return
    X = np.array(L.P, float)[:, :13]
    h = ConvexHull(X - X.mean(0))
    simp = [tuple(sorted(int(v) for v in s)) for s in h.simplices]
    in_listed = all(any(set(s) <= L.F_ZERO[j] for j in range(82)) for s in simp)
    ridges = {}
    for s in simp:
        for k in range(13):
            r = s[:k] + s[k + 1:]
            ridges[r] = ridges.get(r, 0) + 1
    cycle = all(c % 2 == 0 for c in ridges.values())
    # point count in one facet
    j0 = next(j for j in range(82) if not L.F_IS_MINOR[j])
    inF = [s for s in simp if set(s) <= L.F_ZERO[j0]]
    rng = random.Random(5)
    s0 = inF[0]
    w = [Fraction(rng.randint(1, 1000)) for _ in range(13)]
    tw = sum(w)
    p = [sum(w[k] * L.P[s0[k]][a] for k in range(13)) / tw for a in range(14)]
    odd, clean = 0, True
    for s in inF:
        A = [[L.P[s[k]][a] for k in range(13)] for a in range(14)]   # 14 x 13
        if L.rank(A) != 13:
            clean = False
            continue
        R, piv = L.rref([row + [p[a]] for a, row in enumerate(A)])
        if 13 in piv:
            clean = False   # p not in affine span (cannot happen inside one facet)
            continue
        lam = [R[i][13] for i in range(13)]
        if all(x > 0 for x in lam):
            odd += 1
        elif all(x >= 0 for x in lam):
            clean = False   # p on a simplex boundary: choose another point
    check("prop4.1_completeness_exact_certificate",
          in_listed and cycle and clean and odd % 2 == 1,
          f"{len(simp)} simplices, all in listed facets: {in_listed}; mod-2 cycle: {cycle}; "
          f"point in facet f_{j0} covered {odd} time(s)")


# ---------------------------------------------------------------- Leibniz sums
def test_leibniz():
    rng = random.Random(2026)
    Os = [L.random_rational_SO4(rng) for _ in range(24)]
    check("random_rational_O_in_SO4_nonzero_entries", all(L.is_special_orthogonal(O) for O in Os)
          and all(x != 0 for O in Os for r in O for x in r), f"{len(Os)} matrices")
    keys = [((i,), (j,)) for i in range(1, 5) for j in range(1, 5)] + \
        [(I, J) for I in itertools.combinations(range(1, 5), 2)
         for J in itertools.combinations(range(1, 5), 2)]
    nv = {k: L.n_vec(*k) for k in keys}
    ok = all(L.Phi(nv[k], O) == L.minor_det(O, *k) ** 2 for O in Os for k in keys)
    check("eq4.2_Phi_n_IJ_eq_det_O[I,J]^2", ok, f"{len(keys)} (I,J) x {len(Os)} matrices")
    check("Phi_kappa_eq_1", all(L.Phi(L.KAPPA, O) == 1 for O in Os))
    ok = True
    for O in Os:
        Q = [[Fraction(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(4)] for _ in range(4)]
        ok &= L.Phi(L.E_coords(Q), O) == sum(Q[i][j] * O[i][j] ** 2 for i in range(4) for j in range(4))
    check("eq3.8_Phi_E(Q)_eq_sum_q_ij_O_ij^2", ok)
    ok = True
    for O in Os[:6]:
        for pi in L.PERMS:
            for eps in (1, -1):
                lhs = L.Phi(L.h_vec(pi, eps), O)
                rhs = (sum(O[i][pi[i]] ** 2 for i in range(4)) + eps * L.Phi(L.d_vec(pi), O)) / 2
                ok &= lhs == rhs
    check("eq4.4_Phi_h_formula", ok)
    # signed permutation matrices with det 1
    ok = True
    for k, s in enumerate(L.PERMS):
        for signs in itertools.product((1, -1), repeat=4):
            O = [[signs[i] if j == s[i] else 0 for j in range(4)] for i in range(4)]
            if L.det(O) != 1:
                continue
            c = [Fraction(rng.randint(-5, 5)) for _ in range(14)]
            ok &= L.Phi(c, O) == L.values(c)[k]
    check("Phi_at_signed_permutation_eq_value", ok)
    Os_ = L.O_STAR
    check("O_star_in_SO4", L.is_special_orthogonal(Os_))
    idp = L.IDENTITY
    a, b, c = L.Phi(L.t_vec(idp), Os_), L.Phi(L.d_vec(idp), Os_), L.Phi(L.h_vec(idp, -1), Os_)
    check("eq4.3_O_star_values", (a, b, c) == (Fraction(4, 49), Fraction(512, 2401), Fraction(-158, 2401)),
          f"Phi_t={a}, Phi_d={b}, Phi_h(id,-1)={c}")
    # Example 1.3: psi(P_sigma) multiplicities 13, 10, 1 = value profile of h_{id,-1}
    hv = L.values(L.h_vec(idp, -1))
    check("ex1.3_psi_values_on_permutations_13_10_1", sorted(hv) == [0] * 13 + [1] * 10 + [2])


# ---------------------------------------------------------------- spectral side
def rq(rng, n=9):
    return Fraction(rng.randint(-n, n), rng.randint(1, n))


def rqi(rng):
    return L.QI(rq(rng), rq(rng))


def test_spectral():
    rng = random.Random(77)
    ok_pm = ok_35 = ok_gap = ok_hom = True
    for trial in range(12):
        a = [rqi(rng) for _ in range(4)]
        b = [rqi(rng) for _ in range(4)]
        al = rqi(rng)
        m = L.matching_vector(a, b, al)
        F = L.values(m)
        ok_pm &= F == [al * z for z in L.sigma_points(a, b)]
        pts = [(x, 1) for x in a] + [(-x, 1) for x in b]
        ok_hom &= L.matching_vector_from_points(pts, al) == m
        O = L.random_rational_SO4(rng)
        ps = rng.randrange(24)
        c = L.gap_vector(m, ps)
        v = L.values(c)
        ok_gap &= v == [L.re_part(F[k] - F[ps]) for k in range(24)] and v[ps] == 0
        lhs = L.Phi(c, O)
        rhs = L.re_part(al * L.det_Da_plus_ODbOt(a, b, O) - F[ps])
        ok_35 &= lhs == rhs
    check("eq3.1_P_matching_vector_eq_alpha_prod(a_i+b_sigma(i))_gaussian_rational", ok_pm)
    check("matching_vector_homogeneous_eq_dehomogenized", ok_hom)
    check("eq3.3_gap_vector_values_eq_Re_g_sigma", ok_gap)
    check("eq3.5_Phi_c(x)(O)_eq_Re(alpha det(Da+ODbOt)-F_pi*)", ok_35, "12 Gaussian-rational spectra")
    # matching vectors satisfy grad C = 0 (consistency of library pieces)
    import sympy
    ms, C = L.cubic_sympy()
    grads = [sympy.lambdify(ms, sympy.diff(C, m_)) for m_ in ms]
    a = [rqi(rng) for _ in range(4)]; b = [rqi(rng) for _ in range(4)]
    m = L.matching_vector(a, b, rqi(rng))
    check("grad_C_zero_at_gaussian_rational_matching_vector", all(g(*m) == 0 for g in grads))
    # contact set on a real example
    a = [Fraction(x) for x in (0, 1, 3, 7)]; b = [Fraction(x) for x in (0, 2, 5, 6)]
    m = L.matching_vector(a, b)
    F = L.values(m)
    k = min(range(24), key=lambda s: F[s])
    T = L.contact_set(L.gap_vector(m, k))
    check("contact_set_real_example", k in T and all(F[s] == F[k] for s in T), f"T = {sorted(T)}")


# ---------------------------------------------------------------- symmetry group
def test_group():
    maps = set(L._G_MAPS)
    check("G_order_1152_faithful", len(L.G_ELEMENTS) == 1152 and len(maps) == 1152)
    comp_closed = all(tuple(t2[t1[j]] for j in range(24)) in maps
                      for t1 in L._G_MAPS[::37] for t2 in L._G_MAPS)
    check("G_closed_under_composition_(sampled_32x1152)", comp_closed)
    T85 = L.CONTACT_SETS[85]
    ok_rows = []
    for (j, nu, r1, r2, e, img) in L.TABLE3:
        g = (L.perm_from_oneline(r1), L.perm_from_oneline(r2), e)
        im = L.g_act_set(g, L.CONTACT_SETS[nu])
        ok_rows.append(im == frozenset(img) and im <= T85 | {j})
    check("table3_all_7_rows_with_convention_rho.sigma.rho2", all(ok_rows), f"{ok_rows}")
    # other conventions, for the record
    alt = {}
    for name, fn in [("rho2.sigma.rho", lambda g, s: L.compose(g[1], L.compose(L.inverse(s) if g[2] else s, g[0]))),
                     ("rho^-1.sigma.rho2^-1", lambda g, s: L.compose(L.inverse(g[0]), L.compose(L.inverse(s) if g[2] else s, L.inverse(g[1]))))]:
        cnt = 0
        for (j, nu, r1, r2, e, img) in L.TABLE3:
            g = (L.perm_from_oneline(r1), L.perm_from_oneline(r2), e)
            im = frozenset(L.perm_index(fn(g, L.PERMS[x])) for x in L.CONTACT_SETS[nu])
            cnt += im == frozenset(img)
        alt[name] = cnt
    info(f"Table 3 rows matched under alternative conventions: {alt}")
    # covariance (4.7)
    rng = random.Random(4242)
    O = L.random_rational_SO4(rng)
    c = [rq(rng) for _ in range(14)]
    base = L.Phi(c, O)
    w = L.leibniz_weights(O)
    ok_w = ok_phi = ok_flip = True
    for g in L.G_ELEMENTS:
        gO = L.g_act_O(g, O)
        dg = L.det(gO)
        wg = L.leibniz_weights(gO)
        t = L.g_index_map(g)
        ok_w &= all(wg[t[k]] / dg == w[k] for k in range(24))
        gc = L.g_act_c(g, c)
        ok_phi &= L.Phi_O4(gc, gO) == base
        if dg == -1:
            gO2 = [[-x if j == 0 else x for j, x in enumerate(r)] for r in gO]
            ok_flip &= L.Phi(gc, gO2) == base
    check("eq4.7_w_g.sigma(g.O)_eq_w_sigma(O)_all_g", ok_w)
    check("eq4.7_Phi_g.c(g.O)_eq_Phi_c(O)_all_1152_g", ok_phi)
    check("eq4.7_det-1_column_flip_returns_to_SO4", ok_flip)
    ok = True
    for _ in range(20):
        O = L.random_rational_SO4(rng); c = [rq(rng) for _ in range(14)]
        g = rng.choice(L.G_ELEMENTS)
        ok &= L.Phi_O4(L.g_act_c(g, c), L.g_act_O(g, O)) == L.Phi(c, O)
    check("eq4.7_random_20", ok)
    # G maps matching vectors to matching vectors (relabelled / swapped spectra)
    ok = True
    for _ in range(30):
        a = [rqi(rng) for _ in range(4)]; b = [rqi(rng) for _ in range(4)]; al = rqi(rng)
        g = rng.choice(L.G_ELEMENTS)
        a2, b2 = L.g_act_spectra(g, a, b)
        ok &= L.solve_P(L.g_act_values(g, L.values(L.matching_vector(a, b, al)))) == \
            L.matching_vector(a2, b2, al)
    check("G_preserves_X_via_relabelled_spectra", ok, "30 random g, Gaussian-rational spectra")
    fset = set(L.F_VALUES)
    ok = all(tuple(L.g_act_values(g, list(v))) in fset for g in L.G_ELEMENTS for v in L.F_VALUES)
    check("G_permutes_f0..f81", ok)
    mset = {v for j, v in enumerate(L.F_VALUES) if L.F_IS_MINOR[j]}
    ok = all(tuple(L.g_act_values(g, list(v))) in mset for g in L.G_ELEMENTS for v in mset)
    check("G_permutes_minor_vectors", ok)


# ---------------------------------------------------------------- contact sets
def test_contact_sets():
    check("E_90_six_sets", len(L.CONTACT_SETS) == 90 and all(len(T) == 6 for T in L.CONTACT_SETS.values()))
    check("table5_27_seven_sets_distinct", len(L.AUX7) == 27 and all(len(T) == 7 for T in L.AUX7.values())
          and len(set(L.AUX7.values())) == 27)
    t = time.time()
    orbits = L.enumerate_six_set_orbits()
    check("six_sets_214_orbits", len(orbits) == 214 and sum(n for _, n in orbits) == 134596,
          f"{len(orbits)} orbits, {time.time() - t:.1f}s")
    nonminor = {rep for rep, _ in orbits if not all(L.F_IS_MINOR[j] for j in L.generators(rep))}
    check("orbits_124_all_minor_90_nonminor", len(orbits) - len(nonminor) == 124 and len(nonminor) == 90)
    reps = {nu: L.orbit_rep(T) for nu, T in L.CONTACT_SETS.items()}
    check("E_T1..T90_are_the_90_nonminor_orbits_one_each", set(reps.values()) == nonminor
          and len(set(reps.values())) == 90)
    ok = all(L.find_g(T, reps[nu]) is not None for nu, T in list(L.CONTACT_SETS.items())[:10])
    check("find_g_maps_T_to_orbit_rep_(T1..T10)", ok)
    info(f"F_T85 = {[L.f_name(j) for j in L.generators(L.CONTACT_SETS[85])]}")
    info("orbits of the 27 seven-sets: " + ", ".join(
        f"T7_{k}:{'nonminor' if not all(L.F_IS_MINOR[j] for j in L.generators(T)) else 'minor-only'}"
        for k, T in L.AUX7.items()))


# ---------------------------------------------------------------- numeric helpers
def test_numeric():
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        check("numeric_skipped", False, "numpy missing")
        return
    O = L.so4_from_params([0.3, -1.1, 0.7, 0.2, 2.0, -0.4])
    check("numeric_so4_param_orthogonal", float(abs(O @ O.T - np.eye(4)).max()) < 1e-12
          and abs(np.linalg.det(O) - 1) < 1e-12)
    Oq = L.random_rational_SO4(random.Random(9))
    c = [rq(random.Random(10)) for _ in range(14)]
    check("numeric_Phi_np_matches_exact", abs(L.Phi_np(c, np.array(Oq, float)) - float(L.Phi(c, Oq))) < 1e-12)
    v, _ = L.min_phi_over_SO4(L.n_vec((1, 2), (3, 4)), restarts=10, seed=1)
    check("numeric_min_phi_minor_vector_nonnegative", v > -1e-10, f"min = {v:.3e} (10 restarts, seed 1)")
    v, Om = L.min_phi_over_SO4(L.h_vec(L.IDENTITY, -1), restarts=20, seed=1)
    check("numeric_min_phi_h_id_minus_below_O_star_value", v <= -158 / 2401 + 1e-12,
          f"min = {v:.6f} vs Phi(O*) = {-158 / 2401:.6f} (20 restarts, seed 1)")
    r = L.random_point_on_contact_stratum(L.CONTACT_SETS[85], seed=3, tries=40)
    if r is None:
        info("random_point_on_contact_stratum(T85, seed=3) found no point (best effort)")
    else:
        check("numeric_contact_stratum_point_T85", set(L.CONTACT_SETS[85]) <= set(r["contact"])
              and min(r["gaps"]) > -1e-9, f"contact found = {r['contact']}")
    r = L.random_point_on_contact_stratum(L.CONTACT_SETS[1], seed=3, tries=40)
    if r is None:
        info("random_point_on_contact_stratum(T1, seed=3) found no point (best effort)")
    else:
        check("numeric_contact_stratum_point_T1", set(L.CONTACT_SETS[1]) <= set(r["contact"])
              and min(r["gaps"]) > -1e-9, f"contact found = {r['contact']}")
    r = L.random_point_on_contact_stratum(L.CONTACT_SETS[1], seed=3, tries=40, margin=1e-2)
    if r is None:
        info("random_point_on_contact_stratum(T1, margin=1e-2) found no point (best effort)")
    else:
        check("numeric_contact_stratum_point_T1_exact_contact", r["contact"] == sorted(L.CONTACT_SETS[1])
              and min(r["gaps"]) > -1e-9, f"contact found = {r['contact']}")


if __name__ == "__main__":
    for fn in [test_perms, test_basis_and_P, test_cubic, test_entry, test_minors_and_facets,
               test_prop41_completeness, test_leibniz, test_spectral, test_group,
               test_contact_sets, test_numeric]:
        t = time.time()
        try:
            fn()
        except Exception as ex:  # a crash is a failure, never a pass
            import traceback
            traceback.print_exc()
            check(fn.__name__ + "_crashed", False, repr(ex))
        info(f"{fn.__name__} took {time.time() - t:.1f}s")
    n, k = len(RESULTS), sum(RESULTS)
    print(f"SUMMARY {k}/{n} checks passed in {time.time() - T0:.1f}s")
    sys.exit(0 if k == n else 1)
