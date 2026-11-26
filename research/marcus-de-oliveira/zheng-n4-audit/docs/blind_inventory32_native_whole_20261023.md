# Audit: blind_inventory32_native_whole_20261023 (primary proof of paper case T_33)

Evidence level: EXACT for every algebraic identity, cone computation and group action the proof uses (own code, `replay/blind_inventory32_native_whole_20261023__check.py`, all checks pass); PROOF-READ for the analytic steps; CITED for the boundary cases, Appendix B/C theorems, the harmonic lemma and Theorem 1.2. Numeric falsification was attempted and produced no admissible sample (see D).

PDF: `proofs/manin_transfer_20261021/blind_inventory32_native_whole_20261023.pdf`, 4 pages. Auditor label C13-T32-T34.

## Deviations and limits

- The document names its boundary targets only through the certificate ("the 69 preceding contact-case theorems"). I resolved each added contact myself with the symmetry group.
- Numeric falsification failed to produce samples. The library sampler found no point on the T_33 stratum in 6 seeds with margin 0.01, and my own sampler in the document's chart found none in 30 s. So the main inequality has **no** numeric evidence from this audit. Only the exact replays and the reading support it.
- Boundary cases, Thm C.1, Thm B.1, the harmonic lemma and Thm 1.2 are CITED, not audited.

## A. Statements

**Theorem 1 (p. 1).** Let T = {0, 3, 6, 7, 11, 12} (one-line 1234, 1342, 2134, 2143, 2431, 3124). For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 on T, Phi_{c(x)}(O) >= 0 for all O in SO(4). Extra contacts, zero coefficients, collisions and nongeneric parameters are all included. T equals the paper's T_33 literally (identity group element, EXACT), so the statement is (*) for T_33 with all spectra.

## B. Dependencies

| added contact j | resolution (own code; 6-subset of T ∪ {j} mapped onto) | paper stage |
|---|---|---|
| 1 | T_19 or T_22 | 1 |
| 2 | T_29, T_7, T_36 | 1 |
| 5 | T_7, T_29, T_36 | 1 |
| 9 | T_7, T_29, T_36 | 1 |
| 10 | T_4, T_31 | 1 |
| 13 | T_29, T_7, T_36 | 1 |
| 14 | T_31, T_4 | 1 |
| 16 | T ∪ {16} is in the orbit of T^(7)_22 (the document says so); every 6-subset other than T is in the orbit of T_32, T_89 or T_33 itself | 1 (7-set) |
| 17 | T_34, T_21, T_22 | 1 |
| 18 | T_31, T_4 | 1 |
| 20 | T_40 | 1 |
| 21 | T_4, T_31 | 1 |
| 23 | T^(7)_22 (orbit) | 1 |
| 4, 8, 15, 19, 22 | f18, f30, f40 all > 0 there: pure minor sum (EXACT) | none |

T_33 is at stage 5 in Table 6, and every input above is at stage 1, so the order is non-circular. Other inputs: Thm C.1, used when Z is real (A = (∞,0,1,Z) concyclic) or a spectrum repeats; Thm B.1 for cross-spectrum common points; the harmonic lemma; Thm 1.2. All CITED.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = T_33 | finite | EXACT | find_g |
| 2 | Contact rank 6; kernel basis f0,f1,f2,f7,f8,f15,f29,f39; gap formulas (1); identities (2) (three) | identity | EXACT | profiles |
| 3 | Nonnegative representation on the eleven-list (3) (= `generators(T)`) | case split | PROOF-READ | b < 0: g14, g10, g20, g16 give c, ℓ, x0, x1 >= δ, then the third identity. b >= 0: the first or second identity as in T_32. Exhaustive. |
| 4 | f18, f30, f40 equal 1 at indices 8 and 19, so -L_O <= M min(g8, g19) | inequality | EXACT + PROOF-READ | |
| 5 | Added contacts | finite + citation | EXACT / CITED | table B |
| 6 | Chart A = (∞,0,1,Z), B = (b,c,d,e); E, rE, sE, uE | identity | EXACT | |
| 7 | b(V+1)T0 = N0 + V N1; the inverse (7) recovers b, d, c, e from (Z,V,r,s,u), and the chart recovers V, r, s, u | identity | EXACT | both directions |
| 8 | T0 ≠ 0 for nonreal Z; De ≠ 0; all of B finite | argument | PROOF-READ | |
| 9 | H6,H0,H7,H11,H12 = 0,1,r,s,u; H3 - h = C41/D with 41 terms and bidegree (2,2); H8 = -u/Y; H19 as in (8) | identity | EXACT | |
| 10 | P75(Z,Y) has 75 terms and bidegree (4,2); the numerator of P75(Z, Y(Z,V)) is divisible by C41 | identity | EXACT | exact polynomial division, remainder 0 |
| 11 | Disc_V C41 at (2,3,5,7) = 25(2704Z^4 - … + 225), squarefree; (9) ≠ 0 there (= 21840) | identity | EXACT | |
| 12 | P75(Z,0) = Z r u^2 (Z-1)((s-1)Z-s)((r-su)Z+u(s+1)) (same normalisation); its Z^4 coefficient r u^2 (s-1)(r-su); the derivatives (11), (12) | identity | EXACT | |
| 13 | Atlas: Y ≠ 0 or Y = ∞: H8 bounded. Y = 0: Z is one of four distinct roots (distinctness ⇔ (9) factors, checked by hand). H19 → 0 at the two roots other than 0, 1. At 0 and 1, P75 is smooth with both partials nonzero, so Y/Z (resp. Y/(Z-1)) is bounded on every branch. The corner (Z,Y) = (∞,0) is absent. | analytic | PROOF-READ | |
| 14 | H8 nonconstant; C41 absolutely irreducible under the generic conditions | argument | PROOF-READ | |
| 15 | Harmonic lemma on the normalisation of C41 = 0 | citation | CITED | hypotheses as in T_32; see note (ii) |
| 16 | C41 ≡ 0 iff r = 1, s = u = -1, h = 1 (given rsu ≠ 0); there H22 ≡ 0, so g22 = 0 and the cost is a pure minor sum | identity | EXACT | Groebner basis [λ-1, r-1, s+1, u+1, h-1] |
| 17 | Rouché restoration to generic parameters, keeping Z nonreal | analytic | PROOF-READ (details supplied by auditor) | as in T_32: isolated zero on a complex line, Hurwitz, density of the generic set, all conditions open |

## D. Falsification

Attempted, no admissible sample (see Deviations). This is a limit of the samplers, not evidence either way. The cost and chart identities were exercised exactly in place of the numeric test.

## E. Verdict

**confirmed-modulo-citations.** Citations: T_4, T_7, T_19, T_21, T_22, T_29, T_31, T_34, T_36, T_40 (stage 1); T^(7)_22; Thm C.1; Thm B.1; the harmonic lemma; Thm 1.2.

Presentation notes: (i) the boundary targets are named only in the certificate (resolved above). (ii) The document never states that the denominator D of H3 - h = C41/D is nonzero at spectrally realizable points. D is u·(VZ^2(r-su) + 2VZsu - Vsu + Z^2(r-u) + Zu), and it cannot vanish identically on the irreducible curve because its V-degree is 1. So the curve points with D = 0 are finitely many and can be added to the punctures. They are covered by the same Y-atlas, and H3 is holomorphic at genuine spectral points in any case. (iii) The Rouché paragraph omits the density argument (supplied above).
