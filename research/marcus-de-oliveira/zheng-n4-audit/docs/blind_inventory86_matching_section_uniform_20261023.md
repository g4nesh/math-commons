Evidence level: EXACT for every identity and finite assertion in the document (own code, `replay/blind_inventory86_matching_section_uniform_20261023__check.py`, all checks pass); PROOF-READ for the curve and harmonic argument; CITED for the earlier cases T_68, T_51, T_61, the paper's Lemma 3.5, HMSV (X = {grad C = 0}), the GIT lifting and Appendix B; NUMERIC for the falsification runs.

# Audit: blind_inventory86_matching_section_uniform_20261023 (primary proof of paper case T_87)

PDF: `proofs/manin_transfer_20261021/blind_inventory86_matching_section_uniform_20261023.pdf`, 4 pages. Auditor label C44-T87-T89.

## Deviations and limits

- None from DOC_AUDIT. I did not open the certificate `blind_inventory86_matching_section_exact_20261023.py/json`. Every displayed identity was re-derived from the paper's definitions with the shared library. No document numbers or arrays were copied in.
- The document writes F for "the real cubic" and uses matching coordinates from the supplement. I checked that they are the paper's coordinates m_0..m_13 and the paper's cubic C (A.2), with sign +1. Evidence: with my own particular solution v0 of the leaf equations, dC/dm13 - dC/dm2 equals the document's (8) term by term, and the z-coefficient of dC/dm1 equals the document's (7) coefficient -d-e-g+h-u term by term. The generator numbering f_j is the paper's Appendix D numbering: `generators(T_87)` is exactly the document's nine profiles.
- I did not audit the cited earlier cases, Lemma 3.5, HMSV or Appendix B. They are CITED.
- The main-statement search is brief, as instructed (section D).

## A. Statements

**Theorem 1 (p.1).** Let x ∈ X be any matching vector, i.e. alpha times the sigma-points of spectra a, b ∈ C^4, with all spectra allowed, including repeated eigenvalues and a_i + b_j = 0. Suppose the real gap vector c(x) satisfies (P c)_σ ≥ 0 for all σ and (P c)_σ = 0 for σ ∈ T = {4, 5, 7, 9, 12, 14}, with contacts non-strict. Then Φ_{c(x)}(O) ≥ 0 for every O ∈ O(4), hence for every O ∈ SO(4).

- T is literally the paper's T_87, 0-based (lib `CONTACT_SETS[87]`). This is exactly (*) for T_87.
- Degenerate configurations are claimed in the statement. They are handled in the proof: extra contacts by the table, d = 0 as a sum of minors, d < 0 by the swap, full-complex-equal by (5), same-color collisions by Appendix B through the GIT-lift paragraph. Nothing is restricted to a generic stratum.

## B. Dependencies

| dependency | kind | resolution | status |
|---|---|---|---|
| extra contact 2, 10, 13 -> {1,4,5,7,9,12} | paper case | = T_68 (literal); some g ∈ G maps T_68 into T ∪ {j} for each j (EXACT) | CITED (stage 6 < 10) |
| extra contact 3, 8, 16 -> {1,3,5,10,11,12} | paper case | = T_51 (literal); g found for each j (EXACT) | CITED (stage 7) |
| extra contact 18 -> {1,2,3,9,11,12} | paper case | = T_61 (literal); g found (EXACT) | CITED (stage 9) |
| [3] bounded-coefficient harmonic lemma | common lemma | the paper's Lemma 3.5. The hypotheses used are a finite puncture set, meromorphic H_j that are not all constant, the extra-contact hypothesis at boundary points, and -L bounded near the punctures. All match the statement on p.19 | CITED (PROOF-READ in `reduction/contact-reduction.md`) |
| HMSV: {grad C = 0} = X | external | paper Sect. 3 / App. A (CITED [15] there) | CITED |
| GIT lift of nonzero points of X; same-color repetitions and strictly semistable points | core | paper's spectral representation of X. Repeated eigenvalues give MOC by Lemma B.2 (App. B, established) | CITED |
| [2] six-contact reduction "Lifting and scope" | core doc | only used for the lift above | CITED |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | kernel of P[T,:] has dim 8; the eight profiles (2) are a basis (p.1) | identity | EXACT | B3, B4 |
| 2 | the first seven coefficients equal the gaps at 18, 16, 13, 10, 8, 3, 2 (p.2) | identity | EXACT | B5: unit-matrix pattern, f70 = 0 there |
| 3 | generators(T) = {0,1,2,5,12,35,45,59,70}, non-minor f35 and f70 | finite | EXACT | B1, B2 |
| 4 | circuit (3) f35 + f70 = f1+f2+f5+f12+f45+f59 | identity | EXACT | B6 |
| 5 | d = -δ < 0: the gaps at 17, 15, 11, 6, 1, 0 are b+d, ..., k+d, so b..k ≥ δ; absorbing gives a source on f35 | identity | EXACT | B7 |
| 6 | a stabilizer element exchanges f35 and f70 and permutes the seven minors | finite | EXACT | B9 (6 such elements) |
| 7 | extra-contact table covers all 18 noncontacts; three rows reduce to T_68, T_51, T_61; for the 11 others f70 > 0, so gap = 0 forces d = 0 | finite + citation | EXACT + CITED | B13 to B16 |
| 8 | perturb O so that the minors are > 0; B = f70(O) < 0 and 0 < -L ≤ -Bd (4) | inequality | PROOF-READ | the minors are ≥ 0 and d > 0; the perturbation is not even needed for (4) |
| 9 | the stabilizer has 12 elements; the 6 that fix f70 permute the minors and send (4,9) to the 6 listed edges, which form a cycle on T | finite | EXACT | B8, B10 to B12 |
| 10 | if E = 0 on every edge, all six F_i are equal | logic | PROOF-READ | connectivity |
| 11 | (5): (dC/dm13 - dC/dm2)(L8 x + zκ) = -d^2, so complex d = 0, so the real d = 0 | identity | EXACT | C5, C5b (paper C, coefficient for coefficient) |
| 12 | E ≠ 0 normalization: H_i real; the gaps are -η Im H_i; the leaf (6) has [L8 \| κ] of rank 9 and 5 independent constraint rows; v0 is real; the coefficients are -η Im x | identity + logic | EXACT + PROOF-READ | C6, C6b; v0 real gives Re(iη v0) = 0 |
| 13 | X ∩ leaf has local dim ≥ 1 at the point; choose an irreducible curve | analytic | PROOF-READ | affine dimension theorem on the 6-dim irreducible cone |
| 14 | (7): z-coefficient of G1 is -d-e-g+h-u, nonzero at a strict point; if the x-coordinates are constant then z is fixed, so some H_i is nonconstant | identity + logic | EXACT + PROOF-READ | C7a: grad C(m + zκ) = grad C(m) + z H_κ m. INFO line: the coefficient of partial 1 is exactly -d-e-g+h-u. The rank-9 faithfulness argument is correct |
| 15 | (8): G13 - G2 = -d^2 + (h+r+s-u-1)d + r(u-h) on the whole leaf | identity | EXACT | C8, C8b, INFO: my v0 reproduces (8) verbatim. It is also the only constant combination of partials that involves d alone (1-dim solution space) |
| 16 | monic quadratic with constant coefficients on an irreducible curve, so d is constant; a real root is incompatible with Im d ≠ 0 | logic | PROOF-READ | finitely many values on a connected curve |
| 17 | -L ≤ -B·(constant real d) near every puncture | inequality | PROOF-READ | from (4) |
| 18 | boundary points: L = 0 or an extra gap vanishes, and then L ≥ 0 by step 7 (spectral lift, collisions via App. B) | citation | CITED | matches the hypothesis of Lemma 3.5 |
| 19 | Lemma 3.5: the subharmonic extension is constant; a positive constant means no poles, so the H_i are constant, contradicting step 14 | citation + logic | CITED + PROOF-READ | the statement on p.19 matches |
| 20 | exhaustion: non-strict, d ≤ 0, full-complex-equal | case split | PROOF-READ | the cases are exhaustive under the hypotheses |

## D. Falsification

- NUMERIC D1. I took 4000 random complex spectra (seed 7). For each I computed m, E = F9 - F4, r, s, u, h, and the decomposition m/E - v0 = L8 x + zκ. I then evaluated the d-only quadratic of (8). Maximum normalized residual: 9.7e-14. (The identity is polynomial, so it also holds for complex r, s, u, h.)
- NUMERIC D9. Brief main-statement test. I took 6 stratum points from lib `random_point_on_contact_stratum` (seeds 0-5, margin 1e-3, so the contact set is exactly T) and minimized Φ over SO(4) with 8-start BFGS. The minimum is 1.1e-13 ≥ 0. This test proves little.
- I did not try to break the constancy of d on non-spectral points. (8) holds identically on the leaf ∩ {grad C = 0}, which contains all of X ∩ leaf. It is therefore stronger than needed, and it was replayed exactly.

## E. Verdict

**confirmed-modulo-citations.** Every displayed identity was replayed exactly in the paper's coordinates. They are (2), (3), the negative-side rows, the stabilizer and edge claims, the table, (5), the leaf ranks, (7) and (8). The argument is correct as written. Citations: T_68, T_51, T_61; the paper's Lemma 3.5; HMSV generation; the GIT lift and Lemma B.2 for repeated eigenvalues.

Presentation remarks (not gaps):

- p.2: "Perturb only O ... until every good minor f_i(O) > 0" is not needed for (4). It is harmless.
- p.1: the claim that a same-color repetition is safe "by the established repeated-spectrum/one-simple result" rests on Lemma B.2 (repeated eigenvalue in one spectrum), which is established by this audit. Strictly semistable points always contain a same-color repetition, by pigeonhole, so they are covered by the same lemma.
