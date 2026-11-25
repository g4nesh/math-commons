# Bounded numerical convex-hull search

**Result:** one deterministic 200-sample run completed for dimensions 2 through 6. No gap exceeded the selected **heuristic** numerical screening allowance. This is neither a proof of the Marcus–de Oliveira conjecture nor a certified counterexample search over an exact domain. The conjecture remains unresolved by this work.

The run contains 95 samples classified as numerically inside and 105 as numerically near the boundary. Exact orientation tests on the **stored floating-point geometry** classify 86 as outside that stored hull. These exclusions are preserved, not discarded: none exceeds the screening allowance, and boundary controls produce similar effects. A tiny positive floating-point gap must not be presented as a mathematical violation.

## Reproduce or inspect

Run from the repository root. No NumPy, SciPy, model service, or cloud compute is needed.

```sh
# Fixed mathematical and malformed-input fixtures; not another search.
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/search/test_search.py

# Replay exactly the 200 saved samples, generating no new samples.
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/search/replay.py research/marcus-de-oliveira/search/report.json

# Original production command. A new execution is another run against your budget.
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/search/search.py --seed 20261124 --output research/marcus-de-oliveira/search/report.json
```

The original run used Python 3.9.6, Clang 17.0.0, macOS 15.7.9 arm64, and standard-library binary64 arithmetic. Search time was **0.209173 process CPU seconds** and **0.219623625 wall seconds**, measured through analysis/summary preparation and before JSON serialization. The completed process was far below the 60-second CPU budget. There were no generation or analysis errors. `report.json` records full environment information; `SHA256SUMS` identifies the original artifacts. The parent research PR and its commit identify the source revision.

All **11 fixed-fixture tests passed**. Replaying the saved 200 inputs reproduced every analysis field exactly on this platform; see `replay-result.json`. Other Python/libm/platform versions may change the last bits. Such a difference requires examination and is not itself a conjecture violation.

## Coverage and budget

The production run generated exactly **40 samples for each n = 2, 3, 4, 5, 6**, totaling **200**, with no adaptive continuation or optimization. Each dimension used eight samples from each family:

| Family | Construction and purpose |
| --- | --- |
| Generic | Independent complex uniform spectra in [-2,2] + i[-2,2], with complex Gaussian QR mixing. |
| Near monomial | Complex monomial times a unitary rank-one phase perturbation; angles 10^-2 through 10^-9. Searches a boundary neighborhood without assuming every permutation point is a hull vertex. |
| Clustered | Complex spectra near independently sampled centers; spreads 10^-2, 10^-4, 10^-6, 10^-8, each twice per dimension. Exercises thin hulls and numerical conditioning. |
| Repeated | Alternating two-value complex spectra on a half-integer grid, with generic complex QR mixing. Exercises multiplicity and duplicate permutation points. |
| Structured controls | Three scalar-A point-hull cases, three real-spectrum/Hermitian segment-hull cases, and two singular cases with A=0 and a zero b eigenvalue, per dimension. These are controls, not novel positive cases. |

For dimension n, family position f=0,...,4 and within-family index j=0,...,7, the individual seed is `20261124 + 10000*n + 100*f + j`. Every case stores that seed and its full actual input, so replay does not depend on regenerating random numbers.

All sampled U matrices have genuinely complex entries: the minimum, across samples, of the largest entry's absolute imaginary component was **0.3046831152667926**. Generic U is constructed by complex Gaussian QR using two-pass modified Gram–Schmidt. Near-monomial U is constructed as

`D_phase P [I + (exp(i*t)-1) v v* / (v* v)]`.

These constructions describe exact unitary formulas, but their **stored evaluations are approximate**, not exactly unitary certificates. No claim about covering the unitary group or establishing a particular random-matrix measure is needed for this baseline.

## What is computed

For each case, the script constructs the numerical matrix

`M = diag(a) + U diag(b) U*`

and computes its determinant by partial-pivoting elimination. It independently computes the same numerical matrix's determinant by the complete Leibniz permutation sum, using `math.fsum` for real and imaginary accumulation. It enumerates **all n! permutation products** `z_sigma = product_i(a_i+b_sigma(i))`, including duplicates; no angular sampling substitutes for the polygon.

Hull construction converts the stored finite binary64 coordinates to a common power-of-two integer scale. An integer monotone-chain algorithm uses exact orientation signs for these stored values. No arbitrary geometric epsilon collapses thin polygons. This is exact geometry **of rounded data**, not exact evaluation of the original mathematical construction.

- For a two-dimensional hull, the reported signed gap is the maximum over its edges of `-cross(edge, d-edge_start)/edge_length`. A positive sign places the stored d outside the stored polygon; a negative sign measures minimum inward supporting-line slack.
- For a segment, the script checks collinearity and endpoint projection, and reports Euclidean distance to the segment.
- For a point, it reports distance to that point.

The report records the hull dimension, all hull vertices with generating permutations, the strongest edge and outward unit normal for polygons, determinant discrepancy, unitary residual, normality residual of `U diag(b) U*`, and the residual of its prescribed eigenvector equation `B U = U diag(b)`.

## Observed results

| Dimension | Samples | Numerically inside | Near boundary | Raw stored-float hull exclusions |
| --- | ---: | ---: | ---: | ---: |
| 2 | 40 | 0 | 40 | 39 |
| 3 | 40 | 15 | 25 | 17 |
| 4 | 40 | 26 | 14 | 10 |
| 5 | 40 | 27 | 13 | 10 |
| 6 | 40 | 27 | 13 | 10 |
| **Total** | **200** | **95** | **105** | **86** |

There were 27 point hulls, 57 segment hulls, and 116 polygon hulls. Dimension two necessarily has at most two permutation points: classifying its results as near the boundary is expected, not evidence of failed membership. The 39 stored-float exclusions in this known baseline make the roundoff issue particularly visible.

Maximum structural residuals were:

- `||U*U-I||_F`: **6.723667003760839e-16**.
- Relative normality residual: **2.2834495157270645e-16**.
- Relative prescribed spectral-equation residual: **3.681193144927517e-16**.
- Absolute LU versus Leibniz determinant disagreement: **9.191947761807672e-13**.

The largest normalized positive gap was **2.1888301775866405e-14**, in `n6-structured_control-07`, a singular point-hull control. Its heuristic allowance was **2.6982815623995585e-08**. Thus the strongest positive numerical gap belongs to a case designed to have exact determinant zero, not a promising counterexample.

`strongest-cases.json` preserves the ten polygon cases ranked by normalized signed support gap, including complete spectra and matrices. All 200 complete inputs and analyses remain in `report.json`. No cases were dropped because the sign looked inconvenient.

## The screening allowance is not a certificate

Let `r = ||U*U-I||_F`, `s = max(1,max|a|+max|b|)`, `epsilon = 2^-52`, and `delta = |det_LU-det_Leibniz|`. The explicitly heuristic allowance is

`(256*epsilon*n^3 + 8*r)*s^n + 8*delta`.

The report calls a gap above this allowance an `unverified_outside_candidate`, a gap within either sign of the allowance `numerically_near_boundary`, and other polygon cases `numerically_inside`. There were no unverified outside candidates under this rule. The allowance is intentionally conservative, has not been proved as a forward-error bound, and can hide small genuine gaps. All raw gaps are saved so reviewers can choose another analysis without another search.

The approximate witness schema is `mdo.approximate-complex-pair-search.v1`, with `approximate: true`. Each complex entry is a pair of **JSON numbers**. It is deliberately different from the exact Gaussian-rational certificate schema. Converting decimal strings into rationals would certify a different, generally nonunitary matrix. Do not do that and call it a certificate for this search.

## Validation and limitations

The fixed tests cover square inside/boundary/outside geometry, duplicate point hulls, segment projection and perpendicular exclusion, an extremely thin polygon, nonfinite/empty hulls, a complex two-by-two mixing formula, monomial permutation membership, an exactly singular sum, a determinant requiring row pivoting, complex QR orthogonality, and malformed/nonunitary inputs.

During test development, a monomial regression incorrectly expected every permutation point to be a hull vertex. The fixture was corrected to require membership, since a permutation point can be interior to the hull of all permutation points. This was a test assumption error, not a changed search result. The production run happened only after all fixed fixtures passed.

This baseline does not optimize over U or spectra, does not provide rigorous interval enclosures, does not test n>6, and does not establish any universal theorem. Shared authorship under one GitHub owner is not independent verification. Numerical status labels have the exact scope described above.

## Next useful direction

Use the exact verifier team's rational complex Householder or Cayley construction to search on an exactly unitary rational domain, or construct rigorous interval enclosures for a chosen parameterization. Preserve a strictly positive separating support gap against **all** permutation points before proposing a counterexample. For local exploration, examine support variation near a selected exposed monomial point, but do not assume a local extremum is global or that every permutation point is exposed. A new search requires a new finite budget; this pass is complete.
