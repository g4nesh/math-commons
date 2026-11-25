"""Generate contact-configs.json: numerical seeds for the 90 contact sets T_1..T_90 (n = 4).

Evidence level: NUMERIC (floating point data, not certificates).

Usage (from any directory; paths are relative to this script):
    python -I -B make_contact_configs.py part K1 K2 SCRATCHDIR   # cases K1..K2 -> SCRATCHDIR/part_K1_K2.json
    python -I -B make_contact_configs.py merge SCRATCHDIR        # merge parts -> contact-configs.json
"""
import time, sys, os, json, itertools, glob, math
import numpy as np
from scipy.optimize import nnls

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import moc4lib as L

PER_CASE = 3
MARGINS = [0.05, 0.02, 0.01]       # tried in this order for each seed
SEEDS_PER_CASE = 12                # seeds base = 1000*k + s
TRIES = 25
PERMS = list(itertools.permutations(range(4)))   # lexicographic, independent of the lib


def F_direct(a, b, alpha):
    return np.array([alpha * np.prod([a[i] + b[p[i]] for i in range(4)]) for p in PERMS])


def polish(a, b, al, T, ps):
    """Gauss-Newton with minimum-norm steps on the 5 contact equations (18 real unknowns:
    a, b, log-alpha; alpha kept as a complex number, scaled later)."""
    T = sorted(T)
    x = np.concatenate([a.real, a.imag, b.real, b.imag, [al.real, al.imag]])

    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], x[16] + 1j * x[17]

    def res(x):
        A, B, AL = unpack(x)
        F = F_direct(A, B, AL).real
        return np.array([F[s] - F[ps] for s in T if s != ps])

    for _ in range(30):
        r = res(x)
        if np.max(np.abs(r)) < 1e-16:
            break
        J = np.zeros((len(r), 18))
        h = 1e-7
        for k in range(18):
            e = np.zeros(18); e[k] = h
            J[:, k] = (res(x + e) - res(x - e)) / (2 * h)
        x = x - np.linalg.pinv(J, rcond=1e-12) @ r
    return unpack(x)


def build_record(k, T, sol, margin, seed):
    T = sorted(T)
    ps = T[0]
    a, b, al = polish(np.array(sol["a"], complex), np.array(sol["b"], complex),
                      complex(sol["alpha"]), T, ps)
    F = F_direct(a, b, al).real
    gaps = F - F[ps]
    if gaps.min() < -1e-9 * max(1.0, abs(gaps).max()):
        return None
    tot = gaps.sum()
    al = al / tot                       # positive real rescaling of alpha: gaps scale linearly
    F = F_direct(a, b, al).real
    gaps = F - F[ps]
    gaps = np.where(np.isin(np.arange(24), T), gaps, gaps)
    mx = gaps.max()
    resid_abs = max(abs(gaps[s]) for s in T)
    non = [s for s in range(24) if s not in T]
    smallest = float(min(gaps[s] for s in non))
    if smallest <= 1e-4:
        return None
    # decomposition
    gens = L.generators(T)
    A = np.array([[L.F_VALUES[j][s] for j in gens] for s in range(24)], float)
    coef, rn = nnls(A, gaps, maxiter=10000)
    gens_info = [dict(index=j, name=L.f_name(j), is_minor=bool(L.F_IS_MINOR[j]),
                      coefficient=float(coef[i])) for i, j in enumerate(gens)]
    # O_min
    vmin, Omin = L.min_phi_over_SO4(_c_from_gaps(gaps), restarts=40, seed=7)
    rec = dict(case=k, contact_set=T,
               contact_set_oneline=["".join(str(i + 1) for i in PERMS[s]) for s in T],
               a=[[float(z.real), float(z.imag)] for z in a],
               b=[[float(z.real), float(z.imag)] for z in b],
               alpha=[float(al.real), float(al.imag)],
               reference_permutation_index=ps,
               gaps=[float(g) for g in gaps],
               contact_residual=float(resid_abs), contact_residual_relative=float(resid_abs / mx),
               smallest_noncontact_gap=smallest,
               search=dict(seed=seed, margin=margin),
               generators=gens_info,
               decomposition=dict(coefficients_in_generator_order=[float(c) for c in coef],
                                  residual_2norm=float(rn)),
               O_min=dict(matrix=Omin.tolist(), min_phi=float(vmin)))
    # O_hard
    nm = [(coef[i], j) for i, j in enumerate(gens) if not L.F_IS_MINOR[j] and coef[i] > 1e-12]
    if nm:
        cbest, jh = max(nm)
        vh, Oh = L.min_phi_over_SO4(L.F_VECS[jh], restarts=40, seed=11)
        w = L.leibniz_weights_np(Oh)
        phic = float(np.dot(gaps, w))
        D = L.det_Da_plus_ODbOt(list(a), list(b), Oh.tolist())
        z = [complex(np.prod([a[i] + b[p[i]] for i in range(4)])) for p in PERMS]
        rec["O_hard"] = dict(generator_index=jh, generator_name=L.f_name(jh),
                             coefficient=float(cbest), matrix=Oh.tolist(),
                             phi_generator_at_O_hard=float(vh), phi_c_at_O_hard=phic,
                             det_Da_plus_O_Db_Ot=[float(complex(D).real), float(complex(D).imag)],
                             sigma_points=[[float(q.real), float(q.imag)] for q in z])
    else:
        rec["O_hard"] = None
    return rec


def _c_from_gaps(gaps):
    """A vector c in V_R with P c = gaps (gaps lie in im P by construction)."""
    # least squares in floating point (exact solve_P would reject rounding noise)
    Pm = np.array(L.P, float)
    c, *_ = np.linalg.lstsq(Pm, np.asarray(gaps, float), rcond=None)
    return c


def fast_search(T, seed, tries, margin, tol=1e-8, max_nfev=120):
    """Same formulation as lib random_point_on_contact_stratum (mean off-T gap = 1,
    gaps off T >= margin) but with a small evaluation cap per start, so hopeless starts fail fast."""
    from scipy.optimize import least_squares
    T = sorted(T); ps = T[0]
    others = [s for s in range(24) if s not in T]
    rng = np.random.default_rng(seed)

    def unpack(x):
        return x[0:4] + 1j * x[4:8], x[8:12] + 1j * x[12:16], np.exp(1j * x[16])

    def res(x):
        a, b, al = unpack(x)
        F = L.F_values_np(a, b, al).real
        g = F - F[ps]
        r = [g[s] for s in T if s != ps]
        r += [min(0.0, g[s] - margin) for s in others]
        r.append(np.mean(g[others]) - 1.0)
        return np.array(r)

    for _ in range(tries):
        x0 = rng.normal(size=17)
        sol = least_squares(res, x0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=max_nfev)
        if np.max(np.abs(res(sol.x))) < tol:
            a, b, al = unpack(sol.x)
            return dict(a=a, b=b, alpha=al)
    return None


def run_case(k, T):
    out = []
    taken = []
    for s in range(SEEDS_PER_CASE):
        if len(out) >= PER_CASE:
            break
        seed = 1000 * k + s
        margin = MARGINS[s % len(MARGINS)]
        sol = fast_search(T, seed, 40, margin)
        if sol is None:
            continue
        rec = build_record(k, T, sol, margin, seed)
        if rec is None:
            continue
        g = np.array(rec["gaps"])
        if any(np.linalg.norm(g - t) < 1e-3 for t in taken):
            continue
        taken.append(g)
        out.append(rec)
    return out


def part(k1, k2, scratch):
    res = {}
    t0 = time.time()
    for k in range(k1, k2 + 1):
        recs = run_case(k, sorted(L.CONTACT_SETS[k]))
        res[k] = recs
        print(k, len(recs), round(time.time() - t0, 1), flush=True)
        with open(os.path.join(scratch, f"part_{k1}_{k2}.json"), "w") as f:
            json.dump(res, f)


def independent_contact_check(rec):
    a = np.array([complex(*z) for z in rec["a"]]); b = np.array([complex(*z) for z in rec["b"]])
    al = complex(*rec["alpha"])
    F = F_direct(a, b, al).real
    g = F - F.min()
    tol = 1e-10 * g.max()
    return sorted(int(s) for s in range(24) if g[s] <= tol)


def merge(scratch):
    parts = {}
    for fn in sorted(glob.glob(os.path.join(scratch, "part_*.json"))):
        parts.update({int(k): v for k, v in json.load(open(fn)).items()})
    configs = []
    counts = {}
    worst_res = 0.0; min_phi = 1e9; bad = []
    for k in range(1, 91):
        recs = parts.get(k, [])
        counts[k] = len(recs)
        for i, r in enumerate(recs):
            r["config_id"] = f"T{k}-{i}"
            r["sanity"] = dict(independent_contact_set_equals_T=(independent_contact_check(r) == r["contact_set"]))
            if not r["sanity"]["independent_contact_set_equals_T"]:
                bad.append(r["config_id"])
            worst_res = max(worst_res, r["contact_residual_relative"])
            min_phi = min(min_phi, r["O_min"]["min_phi"])
            configs.append(r)
    doc = dict(
        evidence_level="NUMERIC",
        conventions=dict(
            permutations="Index 0..23 = lexicographic order of one-line notation of S_4 (index 0 = identity 1234). "
                         "contact_set_oneline gives the images of 1,2,3,4. P_sigma has a 1 in row i, column sigma(i).",
            F_sigma="F_sigma = alpha * prod_{i=1..4} (a_i + b_{sigma(i)}); complex numbers are stored as [re, im].",
            gaps="gaps[sigma] = Re F_sigma - Re F_ref with ref = reference_permutation_index = the minimizer, so all gaps >= 0 "
                 "(this is the vector P c(x) of the paper). Rescaled (alpha multiplied by a positive real) so that the 24 gaps sum to 1.",
            contact="contact set = {sigma : gap = 0}; contact_residual is max |gap| over the contact set (absolute, equals relative "
                    "to the largest gap times max gap); smallest_noncontact_gap is the smallest gap outside T.",
            generators="generators = f_j of Appendix D vanishing on T (lib generators(T)), name n<I>,<J> (minor) or h<pi>,<+/-> (non-minor); "
                       "decomposition: gaps ~ sum_j coefficient_j * (P f_j), coefficients >= 0 by NNLS.",
            O_min="multistart (40) BFGS minimizer over SO(4) of Phi_c(O) = sum_sigma sgn(sigma) gaps_sigma prod_i O[i][sigma(i)].",
            O_hard="minimizer over SO(4) of Phi_h for the non-minor generator h with the largest positive coefficient; "
                   "phi_c_at_O_hard = Phi_c(O_hard) for the actual gap vector; det = det(D_a + O D_b O^t) (complex); "
                   "sigma_points = prod_i (a_i + b_sigma(i)) in the 24-permutation order. Phi_c(O) = Re(alpha*det - F_ref) before rescaling, "
                   "after rescaling use the stored alpha.",
        ),
        provenance=dict(
            script="make_contact_configs.py",
            library="research/marcus-de-oliveira/zheng-n4-audit/lib/moc4lib.py",
            seeds=f"configuration search seeds 1000*k + s, s = 0..{SEEDS_PER_CASE - 1}; margins cycle {MARGINS}; "
                  f"40 random starts per seed; min_phi seeds 7 (O_min) and 11 (O_hard)",
            tolerances="stratum solver (own copy of the lib formulation, max_nfev 120, 40 starts) tol 1e-8, then Gauss-Newton polish to ~1e-16; contact set recheck tol 1e-10*max gap; "
                       "configurations kept only if smallest non-contact gap > 1e-4 (sum of gaps = 1)",
            note="All values are floating point numerical data, NOT certificates. min_phi values are upper bounds of true minima.",
        ),
        counts_per_case={str(k): v for k, v in counts.items()},
        configs=configs)
    with open(os.path.join(HERE, "contact-configs.json"), "w") as f:
        json.dump(doc, f, indent=1)
    from collections import Counter
    print("count histogram", sorted(Counter(counts.values()).items()))
    print("zero cases", [k for k, v in counts.items() if v == 0])
    print("short cases", {k: v for k, v in counts.items() if v < 3})
    print("worst relative residual", worst_res, "smallest min phi", min_phi, "bad independent", bad)
    print("n configs", len(configs))


if __name__ == "__main__":
    if sys.argv[1] == "part":
        part(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    else:
        merge(sys.argv[2])
