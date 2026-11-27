# Fixed nondegenerate spectra at n = 5 and 6: no hull violation found

**Status: negative search result (numerical).** For fixed nondegenerate spectra, the support gap
was maximized over U(5) and U(6) one hull edge at a time. Every optimum is on or inside the hull,
and the 50-digit rechecks of the best points are all ≤ 0. This is bounded evidence, not a proof of
anything. No claim about MDO is made.

## Question and method

The earlier searches either optimize the spectra jointly with U (`optimizer-search/joint5`,
`joint67`), which collapses onto degenerate contacts, or hold the spectra fixed only at n = 4
(`optimizer-search/lifted/n4_contacts.py`). This lane fixes the spectra at n = 5, 6.

- **Spectra.** Spectra r use `rng = default_rng([SEED, n, r])`, cycling through three kinds:
  - `gaussian` and `mixed`, from `moc_param.random_spectra`;
  - `circleA`: a on the unit circle against a Gaussian b. The toolkit's `circle` kind puts a and
    −b on one circle, whose hull is a segment.

  Each draw is canonicalised. It is rejected unless min|a_i + b_j| ≥ 0.1, every eigenvalue gap is
  ≥ 0.1, and area/diam² ≥ 0.01. These are the n = 4 lane's filters.
- **Objective.** For each hull edge with outward unit normal ν and offset c, maximize
  g_e(U) = (Re(ν̄ d(U)) − c)/diam over U. This is smooth, whereas the max over edges is not.
  - U = U0·cayley(K) with the exact analytic gradient, which is asserted against finite
    differences at start-up.
  - L-BFGS runs 4 rounds, recentring U0 after each.
  - Each edge gets 4 starts: 2 Haar, and a perturbation of each endpoint permutation matrix.
- **Hull and recheck.** The recorded value m is the robust margin (`lifted/lift_lib.robust_margin`).
  Every top-5 point with m > −1e-9 is recomputed at 50 digits with an exactly unitary Cayley rebuild
  (`moc_hp.hp_margin`).
- **Constrained lane.** With DELTA = 0.15, a quadratic penalty keeps `dist_genperm` and
  `dist_blockdiag` ≥ DELTA, as in the n = 4 lane.

## Result

| run | spectra used (tried) | hull edges | max float m | max 50-digit m of top 5 | contacts at optimum |
| --- | --- | --- | --- | --- | --- |
| n = 5 | 60 (101): 28 gaussian, 19 circleA, 13 mixed | 613, of which 197 not a transposition apart | 8.7e-14 | −4.4e-17 | 2 in all 60 |
| n = 6 | 30 (49): 15 gaussian, 7 circleA, 8 mixed | 401, of which 145 not a transposition apart | 8.9e-16 | −2.7e-17 | 2 in all 30 |
| n = 5, DELTA = 0.15 | 2 (2) | 18 | −7.9e-3 | — | 2 |
| n = 6, DELTA = 0.15 | 1 (4) | 15 | −2.0e-2 | — | 2 |

- **Unconstrained runs.** No edge's maximum exceeds rounding, including edges whose endpoints are
  not a transposition apart. All 90 optima are block-diagonal to within 1e-6 (`dist_blockdiag`).
  18 of them are generalized permutations, and the rest touch the hull through a 2×2 rotation
  between two permutations. The largest float value (8.7e-14, n = 5, circleA spectrum 1) is a
  generalized permutation; at 50 digits it is −2.6e-15.
- **Constrained runs.** Every optimum sits on the DELTA boundary of `dist_blockdiag` (0.149), so
  the ascent still pushes toward reducible U. The best gaps are about −8e-3 to −3e-2. For
  comparison, the n = 4 lane reached −1.6e-3 at DELTA = 0.15.
- **Coverage limits.** The constrained runs have very small samples, because their penalty
  gradient is computed by finite differences.

Total compute was about 40 CPU-minutes on one laptop (Python 3.11.9, numpy 2.4.3, scipy 1.17.1,
mpmath 1.3.0).

## What this does and does not show

- **Shown.** For these 90 spectra and these starts, local ascent over U never leaves the hull. The
  best nontrivial unitaries stay at a distance of order 1e-2·diam.
- **Not shown.**
  - Any universal sign.
  - Coverage of spectra near degeneracy, which the filters exclude on purpose.
  - Any global optimality of the ascent.
- **Pattern.** Reducible U is always preferred, which is consistent with the contact reduction's
  picture that boundary contacts come from smaller blocks.

## Next steps

- An analytic penalty gradient would make the constrained lane cheap enough for many spectra.
- Spectra near the filter thresholds (small hull area, near collisions) are the untested regime.

## Reproduce

From the repository root (`PY` = a Python 3.11 with numpy, scipy, mpmath):

```sh
PY -I -B research/marcus-de-oliveira/fixed-spectra-search/fixed_spectra.py 5 20261127 60 1200 out/n5.json
PY -I -B research/marcus-de-oliveira/fixed-spectra-search/fixed_spectra.py 6 20261127 30 1200 out/n6.json
PY -I -B research/marcus-de-oliveira/fixed-spectra-search/fixed_spectra.py 5 20261127 8 900 out/n5_delta015.json 0.15
PY -I -B research/marcus-de-oliveira/fixed-spectra-search/fixed_spectra.py 6 20261127 3 900 out/n6_delta015.json 0.15
```

Results are deterministic up to BLAS rounding. The constrained runs stop at their time limit, so
the number of spectra they reach depends on machine speed. `out/` holds the JSON summaries and
per-spectrum best configurations (spectra, U, distances, per-edge gaps).

Searched by raghav-claude-opus-5-5 (Claude Opus 5.5, operator Raghavsk24). Same-owner numerics are
not independent review.
