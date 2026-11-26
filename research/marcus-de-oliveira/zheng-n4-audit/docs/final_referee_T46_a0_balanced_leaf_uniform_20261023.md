Evidence level: EXACT for the a = 0 source, zero faces, height formulas, the radial identity with G8, the lambda and threshold algebra and the three sign tables (own code, 18 checks, 0 failures). PROOF-READ for the closure argument at j = lambda and the case split. CITED for the beta = -1 coherent Gram construction (audit foundations so4-atoms-part2b G3), radial necessity and the paper cases T_29, T_36, T_1, T_39, T_7. NUMERIC as support only.

# final_referee_T46_a0_balanced_leaf_uniform_20261023 ("A support inequality when two leaf coefficients are equal", 2 pp.)

This is a branch document, used for paper case T_47 through T^(7)_26 = {1, 3, 5, 6, 7, 8, 14}.

Replay: `replay/final_referee_T46_a0_balanced_leaf_uniform_20261023__check.py`, output in `__check_output.txt`.

## Deviations and limits

1. G8 was rebuilt by my own code. Term counts match the document: 415, 60, 664. The strict monomials I found differ from the printed ones because of the clearing powers: -32 V^5 rather than the constant -1, and -V^9 J^3 rather than -64 V J^3. Both kinds serve the same purpose.
2. The beta = -1 Gram construction is CITED (foundations G3: "balanced q=b=x>1: lambda = K = (x+1)/(x-1) and t = j/lambda gives 1/j - (x-1)/(x+1)", EXACT there).

## A. Statements

- **Lemma 1 (p. 1).** On the 8-ray T46 face, take d = 1, q = b = x > 1 and j > 0. If p >= max{1/j - (x-1)/(x+1), 0}, then Phi >= 0 on SO(4), with any a, c, e >= 0 added. This is a universal SO(4) inequality.
- **Theorem 2 (p. 1).** For x in X with c(x) on the a = 0 face and q = b: Phi >= 0 on SO(4). All zero coefficients are included, and spectra are unrestricted.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| five-ray coherent Gram family, beta = -1 | common lemma | CITED foundations G3 (PROVED there) |
| q (= b), p, j, c = 0 lead to paper T_29, T_7, T_36, T_1, T_39 | paper cases | CITED (Stage 1); EXACT containment |
| radial necessity | common lemma | CITED |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | lambda = (s+4)/s = (x+1)/(x-1) = 1 + (2x+2)/(x^2-1) | identity | EXACT | sympy |
| 2 | t = j/lambda gives threshold 1/j - 1/lambda; t = 1 gives 0 when j > lambda | identity | EXACT | - |
| 3 | wall j = lambda by closedness of the safe cone | analytic | PROOF-READ | the safe set {c : Phi_c >= 0 on SO(4)} is closed, and the threshold tends to 0 as t tends to 1^- |
| 4 | C0 < 0, Q != 0, G8 = 0 is necessary | identity | EXACT | - |
| 5 | 0 < q = b < 1: 415 negative coefficients | certificate | EXACT | - |
| 6 | q = b = 1: 60 negative coefficients, constant -32 | certificate | EXACT | - |
| 7 | q = b > 1 strict region chart: 664 negative coefficients; the chart covers exactly j < lambda, p < 1/j - 1/lambda | certificate + parametrization | EXACT | 1/j - 1/lambda = V/((V+2)J) checked |
| 8 | case split 0 < x < 1, x = 1, x > 1 with j >= lambda or j < lambda, plus zero cases | case split | PROOF-READ | exhaustive |

## D. Falsification

- NUMERIC: 12 samples of Lemma 1 at its threshold (a = c = e = 0, seed 12, x in [1.02, 4], j in [0.05, 3]). Worst min Phi / max coef is 1.1e-15, so there is no violation.

## E. Verdict

**confirmed-modulo-citations** (the beta = -1 Gram family from the foundations, radial necessity, and Stage-1 cases). No gaps.
