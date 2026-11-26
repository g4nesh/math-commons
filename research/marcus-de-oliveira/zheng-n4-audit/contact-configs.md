Evidence level: NUMERIC

# contact-configs.json: numerical seeds for the ninety contact sets (n = 4)

## Deviations and limits

- All values are 64-bit floating point data, not certificates. `min_phi` values are upper bounds of the true minimum over SO(4) (multistart BFGS, 40 restarts).
- Eight of the ninety contact sets got no configuration: T_19, T_21, T_22, T_23, T_33, T_79, T_89, T_90. For these the solver (12 seeds x 40 random starts, margin 0.05/0.02/0.01, 120 evaluations per start) found nothing. This says nothing about whether those strata are empty. The CPU budget was spent before a second, different solver could be tried.
- The stratum search is a copy of the library's formulation (`random_point_on_contact_stratum`: mean off-T gap = 1, gaps off T at least `margin`) with a cap on evaluations per start. The library function itself was too slow on some cases. Cases are searched independently, so counts depend on seeds.
- The smallest non-contact gap in the file is 5.6e-4 (gaps sum to 1), so some configurations are close to the boundary of their stratum. The per-configuration field `smallest_noncontact_gap` is the margin.
- Configurations of one case are required to differ (gap vectors at distance at least 1e-3) but may lie on the same connected piece of the stratum.

## Counts

90 cases: 82 have 3 configurations, 8 have 0 (above). Total 246 configurations. `counts_per_case` in the file lists every case.

## Format

Top level: `evidence_level`, `conventions` (indexing, formulas), `provenance` (script, library path, seeds, tolerances), `counts_per_case`, `configs` (list).

Each config: `config_id` (T<k>-<i>), `case` k, `contact_set` (indices 0..23, lexicographic order of S_4), `contact_set_oneline` (images of 1,2,3,4), `a`, `b`, `alpha` ([re, im] pairs), `reference_permutation_index` (a contact permutation; its F is the minimum), `gaps` (24 values Re F_sigma - Re F_ref, scaled to sum 1), `contact_residual` (max |gap| on T), `contact_residual_relative`, `smallest_noncontact_gap`, `search` (seed, margin), `generators` (index j, name, `is_minor`, coefficient), `decomposition` (NNLS coefficients in generator order, 2-norm residual), `O_min` (matrix, min_phi), `O_hard` (null if no non-minor generator has positive coefficient; else generator, coefficient, matrix, Phi_h at that matrix, Phi_c at that matrix for the actual gaps, det(D_a + O D_b O^t) as [re, im], the 24 sigma-points prod_i (a_i + b_sigma(i))), `sanity`.

F_sigma = alpha * prod_i (a_i + b_sigma(i)); P_sigma has a 1 in row i, column sigma(i). Because alpha was rescaled to make the gaps sum to 1, Phi_c(O) = Re(alpha det(D_a + O D_b O^t) - F_ref) holds with the stored alpha.

## Sanity checks run

- Contact set recomputed from a, b, alpha with an independent numpy computation (itertools permutations, no library): equals T for all 246 configurations (tolerance 1e-10 times the largest gap).
- Worst contact residual 5.3e-15 relative to the largest gap (requirement 1e-12). Largest NNLS decomposition residual 1.5e-15 (all gaps lie in the cone, as expected).
- Smallest `O_min.min_phi` over all configurations: 6.7e-16. None is below -1e-9, so no counterexample candidate and no mpmath re-evaluation was needed.
- 239 configurations have an `O_hard`; the smallest Phi_h at the stored matrices is -0.067442 (the known negative value of a non-minor generator), and the smallest Phi_c(O_hard) for an actual gap vector is 0.0052 (positive). 7 configurations have only minor generators with positive coefficient.

## Regenerate

`python -I -B make_contact_configs.py part K1 K2 SCRATCHDIR` for case ranges, then `... merge SCRATCHDIR`. The script imports `lib/moc4lib.py` next to it. Full runtime is about 15 minutes of CPU (most of it spent on the eight failing cases and on cases 19 to 33).
