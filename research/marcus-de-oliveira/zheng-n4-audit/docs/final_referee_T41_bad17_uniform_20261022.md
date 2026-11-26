# Audit: final_referee_T41_bad17_uniform_20261022

Evidence level: EXACT for the cubic restriction, the radial reduction, the Heron/radial identities, the p=0 factorization, the small-product parametrization and its four positivity certificates (term counts 1360/14373/1330/10831 reproduced), the symmetry and transport group elements, and the Lemma 1 sub-identities; EXACT (randomized identity test, 8 random rational points) for det H_F = -4(c+1)^2 G^2/Q^2; PROOF-READ for the sector logic and Lemma 1; CITED for the T24 coefficient bounds (3), (4), the large-product certificate, and the negative-Heron proposition (audited separately in this batch); NUMERIC for falsification of Lemma 1 and of the large-product inequality.

PDF: proofs/manin_transfer_20261021/final_referee_T41_bad17_uniform_20261022.pdf (6 pages). Auditor label C17-T42. Replay: `replay/final_referee_T41_bad17_uniform_20261022__check.py` (31 checks, all pass, about 60 s).

## Deviations and limits

- The determinant identity (step 3) was tested by exact evaluation at random rational points, not by symbolic expansion. The 14x14 determinant is in 8 variables.
- G is not displayed in full in the note. I built it from (1): G = V D^2 + R0 j D + A0 j^2, with R0 = XY r(e/Y, q/X) and A0 = (R0^2 - Fe Fq)/(4V). This polynomial (1191 terms, degree 8) satisfies the determinant identity for the paper's cubic (A.2). That check confirms (1) together with "R0/(XY) = r" and "R0^2 - 4 V A0 = Fe Fq".
- Lemma 1 (pp.3-5): I re-derived the identities it uses (see C). The case analysis was read, not fully re-derived. The H = 0 slope formulas and the a = 0 sub-case were read only. Lemma 1 was also tested numerically: 3.26 million points in D, edges included, and no violation.
- (7) and (8): the stationary-value closed forms were compared with the exact discriminant roots at one rational parameter point each, in floating point. The sign argument after (8) is EXACT.
- The large-product inequality (A >= 1) is proved in a separate note that I did not audit. It is CITED here and was tested numerically: 4 million points, no violation.
- I did not test the main statement numerically. Spectrally realizable sources are not easy to sample, and the main statement is a special case of the conjecture.

## A. Statement

**Theorem 1 (p.1).** Let x be in X, and suppose c(x), up to a positive scale, equals g = q f3 + p f1 + j f53 + c f6 + e f15 + d f17 + b f16 + a f4 with all coefficients >= 0. Then Phi_{c(x)}(O) >= 0 for every O in SO(4).

All eight generators are generators of paper T_42 = {1,2,3,6,7,8} (EXACT). They are linearly independent, so this is a simplicial 8-dimensional subcone of N_{T_42}. Contacts are non-strict, coefficient zeros are allowed, and every spectrum is covered, because the argument uses only grad C(x) = 0. This note is a branch: covering the full cone N_{T_42}, which has 12 generators, needs the f54 branch and the "elementary generator cancellations" of the primary T41 note. That is not checked here.

**Lemma 1 (p.3).** For a, b >= 0, p > 0 and a + b < 2: sqrt(f(x)k(y)) <= r(x,y) on D = {(x,y) in [0,1]^2 : f(x) >= 0, k(y) >= 0}.

## B. Dependencies

- The T24 "coherent coefficient bounds", transported by the profile map [59,53,55,67,54] -> [1,53,4,16,17]: (3), 1/j <= p + (A-1)/R for A >= 1 (at A = 1 this is the "double-SOC Gram"); and (4), 1/j <= p - kappa0 for A = z^2 < 1, S >= 1 + z, a >= b. **CITED.** A group element mapping those generators as stated exists: (rho, rho', eps) = (0321, 2103, 0) in the library's convention sigma -> rho o sigma o rho' (EXACT). The transported scalar statements themselves were not re-read.
- final_referee_T41_bad17_large_product_uniform_20261022 (h + r >= 0 and (h + r)^2 >= fk for A >= 1). **CITED**, and numerically consistent.
- final_referee_T41_negative_Heron_region_20261022 (Fe < 0, Fq < 0 is not realizable). Audited in this batch: confirmed.
- The paper: X ⊂ Sing C, cubic (A.2) (lib, EXACT).

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | d = 0 gives a minor sum; otherwise normalise d = 1 (p.1) | case split | PROOF-READ | the 7 other generators are minors (lib) |
| 2 | C + hQ, with C = -(c+1)(...) < 0 and Q = -ej-eq-e+jp-jq-j-q-1 (p.1) | identity | EXACT | own restriction of (A.2) |
| 3 | realizable implies F(w) = 0 and det H_F(w) = 0 (auditor: x = w + i delta, H(w) delta = 0, delta != 0 since Q != 0); det H_F = -4(c+1)^2 G^2/Q^2 (p.1) | identity | EXACT (randomized, 8 pts) + PROOF-READ | |
| 4 | (1): G = V D^2 + R0 j D + A0 j^2; R0 and A0 depend only on a, b, p, q, e; G(0,0) > 0 (p.1) | identity | EXACT | |
| 5 | J = jXY/D: 0 < J <= j and G(0,J) = 0; adding (j - J) f53 + c f6 (minors) recovers the source (p.1) | identity | EXACT | G(c,j) = (D/XY)^2 G(0, jXY/D) |
| 6 | R0/(XY) = r; R0^2 - 4 V A0 = Fe Fq with Fe = Y^2 f(x), Fq = X^2 k(y); (2) is the quadratic formula (p.2) | identity | EXACT + PROOF-READ | |
| 7 | (3) and (4) via the T24 transport; the ordering a >= b uses a stabilizer that swaps (a,b) and (q,e) (p.2) | citation | CITED (bounds); EXACT (stabilizer (3210, 3210, 1) exists and fixes T_42) | |
| 8 | A >= 1: h + r >= 0 and (h + r)^2 >= fk, with h = 2V[p + (A-1)/R] (p.2) | citation | CITED; NUMERIC consistent | |
| 9 | small-product chart: 1 - z, S - 1 - z, kappa0 forms and inverse Z; the chart covers 0 < z < 1, S >= 1 + z, a > b > 0 (p.3) | identity | EXACT | a = b is impossible here (auditor) |
| 10 | four positivity certificates H = R(h+r) and E = H^2 - R fN kN in (U, Z, P, xi, eta), both layers (p.3) | certificate | EXACT | numerator and denominator have positive coefficients; term counts match |
| 11 | positive root gives p > kappa0 and 1/j <= 1/J <= p - kappa0, hence (4) (p.3) | inequality | PROOF-READ | |
| 12 | ab = 0, S >= 1: closure of the region; limiting threshold 1/S (p.3) | analytic | PROOF-READ (auditor supplied the approach sequence at a = 1, b = 0: b = delta^2/4, a = 1 + delta) | the threshold itself is in (4), CITED |
| 13 | S < 1 + sqrt(A), p = 0: G(0,j) = (XY + jT)(RXY + jB0) with B0 > 0, so no root (p.5) | identity | EXACT (B0 > 0 read: b^2 + b + 1 - ab > 0 when ab < 1) | |
| 14 | Lemma 1, origin: r00^2 - f(0)k(0) = -4VW/R with W < 0 (p.3) | identity | EXACT | |
| 15 | Lemma 1: f is decreasing (vertex formula); H = f(1) = k(1) = r(1,1); both discriminants > 0 (p.3) | identity | EXACT | disc_x f = 16 a b^2 p V/(b+1)^2 |
| 16 | Lemma 1: zero-fiber signs, f(x_r) = -4(Ap-ap-a-p-1)HV/((b+1)D1^2) (p.4) | identity | EXACT | |
| 17 | Lemma 1: no interior maximum; (6) J22 = -4[a-(b+1)p]^2 V^2/(b+1)^2; the x-matrix determinant equals k(y) y^T J y (p.4) | identity + analytic | EXACT + PROOF-READ | see note N1 |
| 18 | Lemma 1 edges: (7), (8), edge derivatives, the identity after (8), H = 0 slopes, a = 0 (pp.4-5) | inequality | EXACT for the derivatives and the (8) identity; NUMERIC (one point) for the closed forms in (7) and (8); PROOF-READ for the remainder | |
| 19 | zero/negative Heron boundary: r > 0 (p.5) | case split | PROOF-READ | when H = 0, f > 0 on [0,1) and x = e/Y < 1, so the case cannot occur |
| 20 | opposite strict signs: negative discriminant; both strictly negative: the negative-Heron proposition (p.5) | case split | PROOF-READ / CITED (batch doc, confirmed) | |
| 21 | the sectors are exhaustive (p.5) | case split | PROOF-READ | A >= 1; A < 1 with S >= 1 + sqrt A; A < 1 with S < 1 + sqrt A, then p = 0 or p > 0 with the four sign classes |

## D. Falsification

- Lemma 1 (5): 3,260,739 points of D (seed 2026), with a + b < 2, p in [1e-4, 1e3], 5% each of a = 0 and b = 0, and edges x, y in {0, 1}. The maximum of (sqrt(fk) - r)/scale was -1.9e-8, attained near the corner (1,1), where B = 0. No violation.
- Large-product (cited): 4,000,000 points with A >= 1. min (h + r)/scale = 0.018 and min ((h + r)^2 - fk)/scale^2 = 2.5e-8 >= 0. No violation.
- An earlier scratch run of the small-product inequality on 4.3 million points (A < 1, S >= 1 + z) found no violation. It was superseded by the exact certificate (step 10).

## E. Verdict

**confirmed-modulo-citations.** Citations: the T24 coherent coefficient bounds (3) and (4), including the A = 1 double-SOC Gram and the z -> 0 threshold; and the large-product note. The negative-Heron note was confirmed in this batch.

No GAP was found.

Presentation note N1 (p.4, the line "Strict concavity of f gives E_xx > 0"; step 17 is marked PROOF-READ, argument supplied by auditor): f is strictly convex in x, not concave. The needed fact is that sqrt(f) is strictly concave where f > 0, which holds because disc f > 0. Equivalently, at a stationary point E_xx = (k/2f)(f_x^2 - 4 b^2 f) = (k/2f) disc f > 0. The conclusion is right.
