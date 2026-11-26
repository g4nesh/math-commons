Evidence level: EXACT for the seven-contact source, identity (1), zero faces, height formulas, the radial identity with G8, the Lemma 2 algebra and chart (6), and the four sign tables (own code, 23 checks, 0 failures). PROOF-READ for the j-range coverage, the closure limits and the assembly. CITED for the beta = -1 Gram construction (audit foundations G3), the two a = 0 sector notes (audited in this batch), radial necessity and the paper cases. NUMERIC as support only.

# final_referee_T46_a0_uniform_20261023 ("The seven-contact face obtained by adding contact 14 to T46", 3 pp.)

This document proves Table 5's auxiliary proposition T^(7)_26 = {1, 3, 5, 6, 7, 8, 14}, which paper case T_47 uses for its a = 0 boundary.

Replay: `replay/final_referee_T46_a0_uniform_20261023__check.py`, output in `__check_output.txt`.

## Deviations and limits

1. G8 was rebuilt by my own code. The four table term counts equal the document's exactly: 53537, 215772, 50282, 207473, all strictly negative. My strict monomials differ from the printed ones by clearing powers. Mine are -A^27 V^4 J^4 C1^2 (low j) and -A^31 V^12 C1^2 (middle j), and the middle one contains no J, S or E. They serve the same purpose.
2. The beta = -1 Gram construction is CITED to the foundations (G3: threshold (t-1)^2/(j - lambda t^2) and leaf Schur condition j > K t^2).

## A. Statements

- **Theorem 1 (whole, p. 1).** For x in X with P c(x) >= 0 vanishing on T* = {1,3,5,6,7,8,14}: Phi_{c(x)} >= 0 on SO(4). All spectra are included, and the contact set is exactly T^(7)_26 of Table 5.
- **Lemma 2 (p. 2).** On the canonical a = 0 source with qb > 1 and q != b, every realizable cost satisfies the support inequality.

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| generators(T*) = {f0,f1,f2,f5,f6,f27,f28,f59}; (1) f28 + f0 = f27 + f1 + f5 + f2 | core | EXACT | - |
| q, p, j, c, e, b = 0 lead to paper T_29, T_36, T_1, T_39, T_42, T_7 | paper cases | CITED (Stage 1) | EXACT containment on the a = 0 face |
| final_referee_T46_a0_all_small_product_uniform_20261023 (qb <= 1) | branch | CITED; audited in this batch, confirmed-modulo-citations | statement matches |
| final_referee_T46_a0_balanced_leaf_uniform_20261023 (q = b) | branch | CITED; audited in this batch, confirmed-modulo-citations | statement matches |
| beta = -1 coherent Gram family | common lemma | CITED foundations G3 | lambda, K, K - lambda re-checked EXACT |
| radial necessity | common lemma | CITED | - |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | absorption via (1) to form (2), or a minor sum | identity | EXACT | - |
| 2 | zero faces | citation | EXACT + CITED | - |
| 3 | C0 < 0, Q as printed, (q+d) G8 = 0 necessary, G8 with 380 terms | identity | EXACT | - |
| 4 | K - lambda = (q-b)^2/((qb-1)(q+b-2)) > 0 | identity | EXACT | - |
| 5 | j >= K: t = 1, threshold 0; j = K by closedness | analytic | PROOF-READ | - |
| 6 | j < lambda^2/K: t = j/lambda satisfies j > K t^2, threshold 1/j - 1/lambda | identity | EXACT | - |
| 7 | lambda^2/K <= j < K: boundary t = sqrt(j/K); for t' < t the thresholds tend to (1-t)^2/(j - lambda t^2), whose denominator is (K - lambda) t^2 > 0; closedness | analytic | PROOF-READ + EXACT formula | - |
| 8 | chart (6): t^2 = j/K, P* = (K-lambda)/(lambda+KJ)^2, j increasing from lambda^2/K (J = 0) towards K (J -> oo) | parametrization | EXACT | sympy, including dj/dJ |
| 9 | low-j chart j = lambda^2 J/(K(1+J)) covers (0, lambda^2/K) | parametrization | PROOF-READ | - |
| 10 | four sign tables (q > b and b > q, low and middle j) all strictly negative, with strict monomials free of S and E (and of J in the middle table) | certificate | EXACT | denominators certified positive |
| 11 | assembly: qb <= 1 (cited), q = b (cited), qb > 1 with q != b (Lemma 2) | case split | PROOF-READ | exhaustive |

## D. Falsification

- NUMERIC: 12 samples of the Lemma 2 capacity at p = P* (c = e = a = 0, seed 13, 20 BFGS restarts). Worst min Phi / max coef is 9.6e-17, so there is no violation.

## E. Verdict

**confirmed-modulo-citations.** The citations are the foundations beta = -1 Gram family, radial necessity, Stage-1 paper cases, and the two a = 0 sector notes. Both notes are confirmed-modulo-citations in this batch. No gaps.
