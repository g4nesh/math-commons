"""Independent exact replay of the proof of Case T_22 = {2,3,4,8,10,12} (group IV) of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4".

Document replayed (read as a PDF, nothing executed):
    MOC4_SUPP/proofs/manin_transfer_20261021/root_moc4_T21_T22_uniform_20261022.pdf   (3 pages)
The document calls the set "T21" (its numbering is one less than the paper's); the same
document also treats the paper's T_23 (replay_T23.py).

Run:   MOC4_SUPP=/path/to/supplement  python replay_T22.py  [--no-numeric]
MOC4_SUPP is used only (i) to compare the sha256 of the PDF with the one that was read and
(ii) to compare two JSON data files (P, the cubic, the 82 facet vectors) with the objects rebuilt
here.  Nothing from the supplement is executed; without MOC4_SUPP these comparisons are skipped.
Exit status 0 iff every check passes.  Packages: sympy (exact), numpy + scipy (NUMERIC part).
Report: T22.md in this directory.
"""
import os
import sys

sys.dont_write_bytecode = True          # nothing is written to disk by this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_iv_common as R   # noqa: E402

NU = 22
DOC = "proofs/manin_transfer_20261021/root_moc4_T21_T22_uniform_20261022.pdf"
T = (2, 3, 4, 8, 10, 12)                      # the document's "T21"

# ---- statements of the document, restated as data -------------------------------------------
DOC_SOURCES = [0, 1, 2, 6, 10, 53, 54, 55, 56, 67, 69, 70, 72]      # p. 1, list "S21"
DOC_NONMINOR = [54, 56, 69, 70]                                      # p. 2
RECTANGLE = (2, 4, 8, 10)                                            # p. 1: F_2 F_10 = F_4 F_8
# p. 1 only says that F_2 - F_4 and F_2 - F_8 factor into two differences inside one spectrum
# and two cross differences.  The explicit factors below are mine and are verified.
SWAPS = [(2, 4, [(1, 5), (3, 6), (2, 4), (7, 8)]),
         (2, 8, [(2, 7), (4, 8), (1, 3), (5, 6)])]
PAIR_REPL = {(54, 56): [0, 2, 10, 53, 55], (54, 69): [2, 6, 10, 53, 67], (54, 70): [0, 6, 55, 67],
             (56, 69): [1, 2, 2, 10, 53, 72], (56, 70): [0, 1, 2, 55, 72], (69, 70): [1, 2, 6, 67, 72]}
ABSORB = {54: [(1, [0, 6, 53]), (72, [10, 55, 67])],
          56: [(6, [1, 2, 10, 55]), (67, [0, 2, 53, 72])],
          69: [(0, [1, 2, 10, 67]), (55, [2, 6, 53, 72])],
          70: [(10, [0, 6, 72]), (53, [1, 55, 67])]}
# p. 3 table: surviving non-minor generator -> (sign sigma, ordered sources lambda_0..lambda_7)
NORMAL = {54: (-1, (0, 2, 6, 10, 53, 55, 67, 54)),
          56: (+1, (1, 2, 10, 0, 55, 53, 72, 56)),
          69: (+1, (1, 2, 10, 6, 67, 53, 72, 69)),
          70: (-1, (0, 2, 6, 1, 72, 55, 67, 70))}
D_INDEX = 7                                                          # d = r_7


def P24(r0, r1, r2, r3, r4, r5, r6, r7):                             # p. 2 of the document
    d = r7
    return (d ** 3 + d ** 2 * (r0 + r2 + r3 + r4 + r5 + r6)
            + d * (r5 * r2 + r5 * r3 + r5 * r6 + r5 * r0 + r2 * r4 + r2 * r6 + r2 * r0 + r3 * r4
                   + r3 * r6 + r4 * r0 + r6 * r0)
            + r5 * r2 * r6 + r5 * r2 * r0 + r5 * r3 * r6 + r5 * r6 * r0 + r2 * r4 * r0 + r2 * r6 * r0)


NAMED = sorted(set(tuple(sorted(pr)) for _, _, pairs in SWAPS for pr in pairs))   # the eight brackets


def main():
    numeric = "--no-numeric" not in sys.argv
    ck = R.Checker("Replay of Case T_22 = %s  (document root_moc4_T21_T22_uniform_20261022, family '21')"
                   % (list(T),))
    R.print_versions(ck, numeric)
    R.check_document_identity(ck, DOC)
    R.check_foundations(ck)

    # (1) contact set and scope
    R.check_contact_set(ck, NU, T, "T21")
    ck.info("scope: the argument uses only c(x) in N_T (gaps >= 0, gaps = 0 on T), so further "
            "contacts are allowed; coincident labels are delegated (see T22.md, citations)")

    # (2) generators
    gens = R.check_generators(ck, T, DOC_SOURCES, DOC_NONMINOR)

    # (3a) multiplicative relation and the two differences
    ck.section("Multiplicative relation among the contact products (document p. 1)")
    ck.check("EXACT", "the four permutations of the relation are contacts", set(RECTANGLE) <= set(T))
    R.check_rectangle(ck, RECTANGLE, SWAPS)
    ck.info("F_sigma(x) = alpha p_sigma(v); with alpha != 0 and eight distinct labels each bracket is "
            "nonzero, hence F_2 != F_4 and F_2 != F_8")

    R.check_collision_delegation(ck, T, SWAPS)

    # (3b) common real part vanishes
    R.check_real_part_lemma(ck)
    ck.info("applied with (a,b,c,d) = (2,4,8,10): H != 0 would give F_2 = F_4 or F_2 = F_8; so H = 0, "
            "and Re x = c(x) + H kappa = c(x)")

    # (3c) radial necessity
    R.check_radial_necessity(ck)

    # (3d) finite elimination
    ck.section("Finite elimination of the non-minor generators (document p. 2)")
    R.check_linear_identities(ck, [(list(k), v) for k, v in PAIR_REPL.items()], "pair replacements")
    R.check_linear_identities(ck, [([k, r], right) for k in ABSORB for r, right in ABSORB[k]],
                              "absorption replacements")
    R.check_elimination(ck, gens, DOC_NONMINOR, PAIR_REPL, ABSORB,
                        {k: list(v[1]) for k, v in NORMAL.items()}, seed=2201, samples=6000)

    # (3e) cubic identities and their signs
    ck.section("The cubic on the four normal forms (document pp. 2-3)")
    for k in sorted(NORMAL):
        sigma, order = NORMAL[k]
        ck.check("EXACT", "normal form %d: the designated coefficient d = r_%d multiplies the non-minor "
                 "generator f_%d" % (k, D_INDEX, k), order[D_INDEX] == k)
        R.check_cubic_identity(ck, "T_22, surviving f_%d, P24" % k, order, sigma, P24, D_INDEX)
        R.exact_sampling_of_normal_form(ck, "T_22, surviving f_%d" % k, order, sigma, D_INDEX,
                                        seed=2200 + k, samples=10000)
        r, Ck = R.cubic_on_normal_form(order)
        ck.check("EXACT", "normal form %d: C does not depend on r_1, the coefficient of f_%d (the document's "
                 "remark that P24 is independent of r_1)" % (k, order[1]), not Ck.has(r[1]))
    ck.info("so |C(c(x))| >= d^3 > 0 on every normal form, against C(Re x) = C(c(x)) = 0: no non-minor "
            "generator survives for eight distinct labels, and c(x) is a nonnegative combination of "
            "minor vectors")

    # final step
    ck.section("Minor decomposition implies Phi >= 0 on SO(4)")
    R.check_minor_leibniz(ck, [j for j in gens if R.F_IS_MINOR[j]])

    # observation (mine): on the layer of eight distinct labels a seventh contact is forced
    ck.section("Observation, not part of the document: the relation p_2 p_13 = p_3 p_12")
    ck.check("EXACT", "p_2 p_13 = p_3 p_12 identically; 2, 3, 12 are contacts and 13 is not in T; the "
             "closure of T under the 18 two-transposition relations is T + {13}",
             R.pmul(R.P_SIGMA[2], R.P_SIGMA[13]) == R.pmul(R.P_SIGMA[3], R.P_SIGMA[12])
             and {2, 3, 12} <= set(T) and 13 not in T
             and R.rectangle_closure(T) == sorted(set(T) | {13}))
    ck.info("so for eight distinct labels (H = 0, all products nonzero) alpha z_13 = (alpha z_3)(alpha z_12)"
            "/(alpha z_2) is purely imaginary and 13 is a contact as well: the layer of eight distinct "
            "labels has no point with exactly the six contacts T")

    R.crosscheck_with_audit_library(ck)
    R.crosscheck_with_supplement_data(ck)

    if numeric:
        def extra(pts, cov):
            eight = [p for p in pts if p["sep"] > 1e-3]
            if eight:
                ck.check("NUMERIC", "T_%d: every converged point with eight distinct labels has the "
                         "seventh contact 13 (the observation above)" % NU,
                         all(13 in p["contacts"] for p in eight), "%d points tested" % len(eight))
        R.numeric_falsification(ck, NU, T, gens, seed=22001, named_pairs=NAMED, starts_free=400,
                                extra=extra)
        R.numeric_common_weights(ck, seed=22900, samples=60)
    else:
        ck.info("NUMERIC part skipped (--no-numeric)")
    ck.finish()


if __name__ == "__main__":
    main()
