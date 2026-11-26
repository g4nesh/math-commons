Evidence level: EXACT for every identity and finite assertion in the document (own code, `replay/blind_inventory88_matching_curve_uniform_20261023__check.py`, all EXACT checks pass); PROOF-READ for the curve, pole and harmonic arguments; CITED for T_32, T_79, the seven-element set T^(7)_24, the paper's Lemma 3.5, HMSV, the GIT lifting and Appendix B; established for radial necessity (det Hess C(Re x) = 0); NUMERIC for the falsification runs.

# Audit: blind_inventory88_matching_curve_uniform_20261023 (primary proof of paper case T_89)

PDF: `proofs/manin_transfer_20261021/blind_inventory88_matching_curve_uniform_20261023.pdf`, 6 pages. Auditor label C44-T87-T89.

## Deviations and limits

- None from DOC_AUDIT. I did not open the certificate `blind_inventory88_matching_curve_exact_20261023.py/json` or `root_scout_moc_spinor_20261021.json`. The 211-term C211, the 47-term R47 and every solve were recomputed from the paper's definitions.
- The coordinates are the paper's. The displayed offset v0 of (7), read in the paper's coordinates m_0..m_13, gives exactly the five fixed differences (C2). The document's F is the paper's cubic C, with sign +1 in (13) and C1. The f_j numbering is the paper's Appendix D numbering (B1).
- Statement label `wholeSing` (p.1) is the identification (1) {grad F = 0} = M^ together with the necessity (2). Both are inputs, not new theorems of this document. They are CITED (HMSV) and established (radial necessity) respectively.
- The cited earlier cases, Lemma 3.5, HMSV and Appendix B were not audited here.

## A. Statements

**Theorem 1 (p.1).** For every x ∈ X with P c(x) ≥ 0 and (P c(x))_σ = 0 for σ ∈ T = {2, 5, 9, 10, 13, 14} (non-strict contacts; every spectrum, including repeated eigenvalues and cross-collisions), Φ_{c(x)}(O) ≥ 0 for every O ∈ O(4). T is literally the paper's T_89. The statement is (*) for T_89.

**(1)/(2), label wholeSing (p.1).** (1) is {m : grad C(m) = 0} = X̂, quoted from HMSV with the identification checked (CITED). (2) is that C(w) = 0 and det Hess C(w) = 0 for every real w = Re(ν m), m ∈ X. This is the established radial necessity (item (5) of the audit brief), applied to ν m ∈ X.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| extra contact 3, 8, 16 -> {0,2,3,7,11,12} | paper case | = T_32; an image lies in T ∪ {j} for each j (EXACT) | CITED (stage 11 < 12) |
| extra contact 4, 7, 12 -> {1,3,4,8,11,12} | paper case | = T_79 (EXACT) | CITED (stage 8) |
| extra contact 18 -> N = {0,3,7,11,12,16,20} | seven-contact | = Table 5, j = 24 (EXACT) | CITED (stage 1) |
| (2) det Hess = 0 for real matching vectors ([2] global Hessian inertia) | common lemma | only det = 0 is used, which is the established radial necessity; no inertia statement is used | established |
| [4] bounded-negative-part harmonic lemma | common lemma | the paper's Lemma 3.5 (p.19); the hypotheses as used match | CITED (PROOF-READ in reduction audit) |
| HMSV ideal generation; GIT lift; same-color collisions via Lemma B.2 | external/core | as in the paper | CITED |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | kernel dim 8, basis (3), coefficients = gaps (18,16,12,8,4,3,7) | identity | EXACT | B3, B4 |
| 2 | generators(T_89) = {0,1,11,12,39,42,45,70,72}, non-minor f42 and f70 | finite | EXACT | B1, B2 |
| 3 | circuit (4) f42 + f70 = 2f0 + f11 + f39 + f72; gaps 21, 11, 15, 0 = a+2d, c+d, f+d, k+d | identity | EXACT | B5, B6 |
| 4 | the action row=(0,2,1,3), col=(0,3,2,1) (σ ↦ row∘σ∘col) fixes T, swaps f42 and f70, permutes the minors | finite | EXACT | B7 (the printed data, literally) |
| 5 | extra-contact table complete; three reductions; f70 > 0 on the other 11 | finite + citation | EXACT + CITED | B9 to B12 |
| 6 | perturb O so the minors are > 0 (dense open, real-analytic); (5): a_i/d ≤ -B/f_i(O), -L ≤ -Bd | inequality | PROOF-READ | correct |
| 7 | stabilizer 12; six fix f70 and map (2,5) onto a 6-cycle of edges | finite | EXACT | B8 |
| 8 | full-complex-equal: three gradient identities give e = g = b and then b = 0 (gap 16 = 0, so T_32 applies) | identity + logic | EXACT + PROOF-READ | C1: all three identities hold for the paper C. Re ν(2c+b+d) > 0 because these real coefficients are > 0 at a strict point |
| 9 | v0 of (7) gives the differences 0, 1, r, s, u, h; v0 real, so the real coefficients are -η Im x | identity | EXACT | C2, C2b |
| 10 | solve order: each equation is linear in its variable with raw coefficient (8); after substitution it is α·var + β with α, β ∈ Q[r,s,u,h,b,d] | identity | EXACT | C3, C4 (all seven coefficients match (8) up to the constant factor from the normalization of F) |
| 11 | the coefficients (8) are nonzero at a strict point (imaginary parts); reduced denominators are nonzero there | logic | PROOF-READ | the reduced denominator divides the cleared α, which is a product of nonzero factors |
| 12 | (9): the nonzero substituted numerators (indices 0,5,6,9,10,12,13) have the common factor C211 with 211 terms, deg_b 3, deg_d 2, [b^3 d^2] = 16 | identity | EXACT | C6, C7, C7b (the gcd is exactly C211) |
| 13 | one numerator = const·F6·C211 (index 9); one = C211·R47 with R47 47 terms (index 13) | identity | EXACT | C8, C9 |
| 14 | (10): (2d-u+h)^2 R47(d(r-s-1)/(2d-u+h), d) = const·d^2·(four linear factors); all nonzero when Im d ≠ 0, so F6 and R47 have no common zero, and C211 = 0 at every strict realizable point | identity + logic | EXACT + PROOF-READ | C10; deg_b R47 = 2 |
| 15 | conversely, C211 = 0 maps into X where defined; the GIT lift makes these points spectral | identity + citation | EXACT + CITED | the gcd in step 12 covers all 14 partials |
| 16 | irreducible component of C211 = 0 through the point; b and d are not both constant; [L8 \| κ] has rank 9, so some H_i is nonconstant | logic | PROOF-READ | |
| 17 | (11) closed forms for e and g | identity | EXACT | C5 |
| 18 | poles of d have finite b: the b^3 d^2 term dominates at (∞,∞) | logic | PROOF-READ | Newton-polygon argument: deg_b ≤ 3, deg_d ≤ 2, and b^3 d^2 has a nonzero constant coefficient |
| 19 | at a pole of d: e → b - r + u + 1, g → b + h - r + s, so b, e, g are bounded | identity | PROOF-READ | from (11) |
| 20 | unbounded real d: the normalized limit has b = e = g = 0, a, c, f, k finite (by (5)), and normalized height → -(a+d) since Q = -K2 ≤ -1 | identity + analytic | EXACT + PROOF-READ | C11 gives (13), with sign +1 |
| 21 | (14) det Hess C(L8(a,0,c,0,f,0,k,d) - (a+d)κ) = -4 d^4 K2^2 (a+d)^6 < 0, contradicting det = 0 along the sequence (continuity, homogeneity) | identity + analytic | EXACT + PROOF-READ | C12: exact Bareiss, 819 divisions |
| 22 | completion: Lemma 3.5 hypotheses (-L bounded near all punctures by (5) plus steps 18 to 21; boundary handled by step 5; d cannot cross 0 inside Ω) | citation | CITED + PROOF-READ | matches p.19 |

## D. Falsification

- NUMERIC D1. I took 3000 random complex spectra (seed 3) and computed (r,s,u,h,b,d) from E = F5 - F2 and the leaf decomposition. C211 vanishes at every one: the maximum relative residual is 4.3e-15. A control with b perturbed by 0.1 gives a minimum relative residual of 1.8e-8, more than 10^6 times larger. This is an independent check that every spectral point lies on C211 = 0.
- NUMERIC D9. I took 6 stratum points (seeds 0-5, margin 0) and minimized Φ over SO(4) with 8-start BFGS. The minimum is -4.7e-15, which is zero within rounding. This is a brief test only.

## E. Verdict

**confirmed-modulo-citations.** Every displayed identity was replayed exactly with the paper's C and coordinates. They are (3), (4), the action, the table, the three full-complex-equal identities, v0, the solve chain (8), C211 with its degree data, F6 and R47, (10), (11), (13) and (14). The pole argument and the harmonic completion are correct. Citations: T_32, T_79, T^(7)_24, Lemma 3.5, HMSV, the GIT lift and Lemma B.2; radial necessity is established.

Presentation remarks (not gaps):

- p.1: the document justifies (2) by "HMSV's secant/Hessian theorem" and a "global-Hessian note". The elementary radial-necessity argument (paper p.27) gives (2) directly.
- p.3: "every reduced denominator is nonzero at such a point". This is true because the reduced denominator divides the cleared coefficient. The text's induction says the same thing tersely.
