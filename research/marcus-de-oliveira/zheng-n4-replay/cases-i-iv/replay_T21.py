"""Independent exact replay of the proof of Case T_21 = {0,3,5,9,11,12} (group IV) of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4".

Document replayed (read as a PDF, nothing executed):
    MOC4_SUPP/proofs/manin_transfer_20261021/root_moc4_T20_uniform_20261022.pdf   (2 pages)
The document calls the set "T20" (its numbering is one less than the paper's).

Run:   MOC4_SUPP=/path/to/supplement  python replay_T21.py  [--no-numeric]
MOC4_SUPP is used only (i) to compare the sha256 of the PDF with the one that was read and
(ii) to compare two JSON data files (P, the cubic, the 82 facet vectors) with the objects rebuilt
here.  Nothing from the supplement is executed; without MOC4_SUPP these comparisons are skipped.
Exit status 0 iff every check passes.  Packages: sympy (exact), numpy + scipy (NUMERIC part).
Report: T21.md in this directory.
"""
import os
import random
import sys
from fractions import Fraction

sys.dont_write_bytecode = True          # nothing is written to disk by this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_iv_common as R   # noqa: E402

NU = 21
DOC = "proofs/manin_transfer_20261021/root_moc4_T20_uniform_20261022.pdf"
T = (0, 3, 5, 9, 11, 12)                      # the document's "T20"

# ---- statements of the document, restated as data -------------------------------------------
DOC_SOURCES = [0, 1, 2, 7, 9, 12, 29, 39, 40, 41]            # p. 2, "complete source facets"
DOC_NONMINOR = [40]
RECTANGLE = (3, 5, 9, 11)                                     # p. 1: F_3 F_11 = F_5 F_9
# p. 1: F_3 - F_5 = (A0-C0)(A1-A2)(A3-C1)(C2-C3),  F_3 - F_9 = (A0-A3)(A1-C2)(A2-C3)(C0-C1),
# with A_i = a_i = label i+1 and C_j = -b_j = label 5+j; a difference of two labels is a bracket.
SWAPS = [(3, 5, [(1, 5), (2, 3), (4, 6), (7, 8)]),
         (3, 9, [(1, 4), (2, 7), (3, 8), (5, 6)])]
DELTA = {8: 1, 2: -1, 6: -1, 16: -1}                          # p. 2: delta = g8 - g2 - g6 - g16
DELTA_ON_SOURCES = {40: 1, 1: -1, 9: -1}                      # = coeff f40 - coeff f1 - coeff f9
ABSORB = {40: [(1, [0, 7, 39]), (9, [2, 7, 41])]}            # p. 2: f40+f1 = f0+f7+f39, f40+f9 = f2+f7+f41
NORMAL_ORDER = (0, 2, 7, 12, 29, 39, 40, 41)                  # (1): coefficients q,a,b,c,e,f,d,h


def claimed_cubic(q, a, b, c, e, f, d, h):                    # (2) of the document
    M = (b + d) * (q + a + d) + c * (a + d)
    return (f + d + h) * M + e * ((q + a + c + d) * (f + d) + a * c)


NAMED = sorted(set(tuple(sorted(pr)) for _, _, pairs in SWAPS for pr in pairs))   # the eight brackets


def main():
    numeric = "--no-numeric" not in sys.argv
    ck = R.Checker("Replay of Case T_21 = %s  (document root_moc4_T20_uniform_20261022)" % (list(T),))
    R.print_versions(ck, numeric)
    R.check_document_identity(ck, DOC)
    R.check_foundations(ck)

    # (1) contact set and scope
    R.check_contact_set(ck, NU, T, "T20")
    ck.info("scope: the argument uses only c(x) in N_T (gaps >= 0, gaps = 0 on T), so further "
            "contacts are allowed; coincident labels are delegated (see T21.md, citations)")

    # (2) generators
    gens = R.check_generators(ck, T, DOC_SOURCES, DOC_NONMINOR)

    # (3a) multiplicative relation and the two differences
    ck.section("Multiplicative relation among the contact products (document p. 1)")
    ck.check("EXACT", "the four permutations of the relation are contacts", set(RECTANGLE) <= set(T))
    R.check_rectangle(ck, RECTANGLE, SWAPS)
    ck.info("F_sigma(x) = alpha p_sigma(v); with alpha != 0 and eight distinct labels each bracket is "
            "nonzero, hence F_3 != F_5 and F_3 != F_9")

    R.check_collision_delegation(ck, T, SWAPS)

    # (3b) common real part vanishes
    R.check_real_part_lemma(ck)
    ck.info("applied with (a,b,c,d) = (3,5,9,11): H != 0 would give F_3 = F_5 or F_3 = F_9; so H = 0, "
            "and Re x = c(x) + H kappa = c(x)")

    # (3c) radial necessity
    R.check_radial_necessity(ck)

    # (3d) the functional delta and the two replacements
    ck.section("Removing the non-minor generator (document p. 2)")
    ok = True
    for j in gens:
        val = sum(cf * R.F_VALUES[j][s] for s, cf in DELTA.items())
        ok &= val == DELTA_ON_SOURCES.get(j, 0)
    ck.check("EXACT", "delta = g_8 - g_2 - g_6 - g_16 takes the value +1 on f_40, -1 on f_1 and f_9, and 0 "
             "on the other seven generators; so delta(c) = lambda_40 - lambda_1 - lambda_9 for EVERY "
             "nonnegative representation of c", ok)
    R.check_linear_identities(ck, [([40, r], right) for r, right in ABSORB[40]], "replacements")
    R.check_elimination(ck, gens, DOC_NONMINOR, {}, ABSORB, {40: list(NORMAL_ORDER)},
                        seed=2101, samples=4000)
    rng = random.Random(2102)
    ok = True
    npos = 0
    for _ in range(4000):
        mode = rng.random()
        c = {j: Fraction(rng.randint(0, 40), rng.randint(1, 6)) for j in gens if rng.random() < mode}
        vals = R.combo_values(c)
        delta = sum(cf * vals[s] for s, cf in DELTA.items())
        out = R.eliminate(c, DOC_NONMINOR, {}, ABSORB)
        if delta > 0:
            npos += 1
            ok &= out.get(40, 0) == delta and 1 not in out and 9 not in out
        else:
            ok &= 40 not in out
    ck.check("EXACT", "dichotomy on 4000 random rational coefficient vectors (seed 2102): delta <= 0 gives "
             "a minor decomposition; delta > 0 (%d times) gives the form (1) with d = delta and no f_1, f_9"
             % npos, ok and npos > 0)

    # (3e) cubic identity and its sign
    ck.section("The cubic on the form (1) (document p. 2, identity (2))")
    R.check_cubic_identity(ck, "T_21", NORMAL_ORDER, +1, claimed_cubic, d_index=6)
    R.exact_sampling_of_normal_form(ck, "T_21", NORMAL_ORDER, +1, 6, seed=2103, samples=20000)
    ck.info("so C(c(x)) >= d^3 > 0 on the form (1), against C(Re x) = C(c(x)) = 0: d > 0 is impossible "
            "for eight distinct labels, and c(x) is a nonnegative combination of minor vectors")

    # final step: minor vectors have nonnegative Leibniz sums on SO(4)
    ck.section("Minor decomposition implies Phi >= 0 on SO(4)")
    R.check_minor_leibniz(ck, [j for j in gens if R.F_IS_MINOR[j]])

    # observation (mine): no further contact is forced by the eighteen two-transposition relations
    ck.section("Observation, not part of the document: relations with three corners in T")
    ck.check("EXACT", "there are 18 relations p_a p_d = p_b p_c from two disjoint transpositions; "
             "T_21 is closed under completing a relation with three corners in the set",
             len(R.all_rectangles()) == 18 and R.rectangle_closure(T) == sorted(T))

    R.crosscheck_with_audit_library(ck)
    R.crosscheck_with_supplement_data(ck)

    if numeric:
        def extra(pts, cov):
            def delta_of(p):
                return sum(cf * p["gaps"][s] for s, cf in DELTA.items())
            if cov:
                ck.check("NUMERIC", "T_21 (c): no covered point has delta > 1e-7 (claim: delta <= 0)",
                         all(delta_of(p) <= 1e-7 for p in cov),
                         "%d points tested; delta in [%.3e, %.3e]"
                         % (len(cov), min(delta_of(p) for p in cov), max(delta_of(p) for p in cov)))
            else:
                ck.info("T_21 (c): delta <= 0 is NOT tested (no covered point)")
            if pts:
                ck.info("T_21: max delta over all converged points (delegated ones included) = %.3e"
                        % max(delta_of(p) for p in pts))
        R.numeric_falsification(ck, NU, T, gens, seed=21001, named_pairs=NAMED, starts_free=500,
                                eps_list=(0.02, 0.1, 0.5), extra=extra)
        # is the layer of eight distinct labels of the T_21 stratum empty?  (positive control: T_22)
        for nu2, T2 in ((21, T), (22, R.PAPER_T[22])):
            for s0 in (0.05, 0.01):
                n, best, hmax = R.numeric_hinge_search(T2, seed=21200 + nu2, starts=80, s0=s0)
                ck.info("search with all 28 label distances kept >= %g (normalised), T_%d, 80 starts "
                        "(seed %d): %d starts reached residual < 1e-10; smallest residual %.3e%s"
                        % (s0, nu2, 21200 + nu2, n, best,
                           "" if hmax is None else "; max |H| at the zeros %.3e" % hmax))
        ck.info("so no point of the T_21 stratum with eight distinct labels was found; whether that layer "
                "is empty is not decided here and is not needed by the proof")
        R.numeric_common_weights(ck, seed=21900, samples=60)
    else:
        ck.info("NUMERIC part skipped (--no-numeric)")
    ck.finish()


if __name__ == "__main__":
    main()
