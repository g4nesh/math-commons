Evidence level: EXACT for the relabelling, source, decomposition, tie covering, chart, homogeneous inverse, coordinate functions, curve, resultants and specialization (own sympy code `replay/independent_inventory70_native_harmonic_uniform_20261023__check.py`, all checks pass); PROOF-READ for the analytic curve and harmonic argument; CITED for earlier cases, the seven-contact theorem on Table 5 set #27, Appendix B, Theorem 1.2 and the Birkhoff theorem; NUMERIC for a brief falsification.

# Audit of independent_inventory70_native_harmonic_uniform_20261023 (primary proof of paper case T_71)

PDF: `proofs/manin_transfer_20261021/independent_inventory70_native_harmonic_uniform_20261023.pdf`, 4 pages. Audit batch C37-T69-T71.

## Deviations and limits

- I worked from the text extract. Every formula I rely on was replayed from the definitions and agrees with the extract, apart from the sign convention noted below. I did not open the page images.
- Convention. This document uses difference brackets: F_σ = ∏ (A_i − C_σ(i)), where C = (t,0,1,∞) is c = −b. In the affine form, with factor 1 at ∞, I get E = F7 − F1 = −t(y−z)(x−1). The printed value t(y−z)(x−1) is the homogeneous one, where the factor at ∞ is −A_den. All 24 products change sign together, so no ratio H_σ changes. Harmless.
- The seven-contact theorem on {1,4,5,7,9,12,14} (`root_inventory67_g14_uniform_20261023`, Table 5 #27) is CITED, not audited here.
- The falsification is brief: 156 supported points.

## A. Statements

**Theorem 1 (p.1).** Let x ∈ X with P c(x) ≥ 0 and zero gaps on T_71 = {1,3,4,10,11,12}, with non-strict contacts. Then Φ_{c(x)}(O) ≥ 0 for every O ∈ SO(4), and so for unitaries. The literal tuple equals Appendix E's T_71.

The proof runs on the working tuple T = {1,4,5,7,9,14}. It is obtained by relabelling a_2 ↔ a_3 (0-based rows 1, 2), which sends σ ↦ σ∘τ. This is the group element (id, τ, 0) of G, and Φ, X and the cone are equivariant under it (paper (4.7)). EXACT: the image of T_71 is T, and the literal source maps to the working source.

Coverage: coincidences by Thm B.1 / Lemma B.2; seventh ties as in section B; the remainder by the curve argument. So this is (*) for T_71 on all of X.

Internal claims:
- Source and decomposition. Working source {0,1,2,5,11,12,33,35,59}, sole non-minor f35, and f35 + f11 = f1 + f5 + f12 + f33 (2).
- Readouts. (q, p, a, b, j, e, d, c) = (g18, g17, g15, g11, g6, g12, g3 − g12, g0), hence 0 ≤ d ≤ g3 (4).
- Extra ties. Those with f35 = 0 are 0, 6, 11, 12, 15, 17, 18. Tie 12 is the seven-contact theorem.
- Chart. Ratios r, k, v; inverse (5); the homogeneous inverse in (Z,Q); E = −QZKN, H3 = r/Q, F1/E = −r(Z+1)/K (6).
- Curve. C = N14 − hD14 has bidegree (2,2) and 33 terms, with Z² coefficient k(Qv−r)(Qk+r−1) (8).
- Edge identity. H14 − H12 = (r/Q − v)(kQ − 1), so d = −η Im(−vkQ) = ηvk Im Q on C (9).
- Resultants: Res_Z(C,Z) and Res_Q(C,N).
- Z = ∞ projections, and the height points K = 0.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| G-equivariance of (*) under (id, τ, 0) | core (paper (4.7)) | EXACT in lib |
| Lemma 4.2 / Prop 4.1 | core | CITED (generators EXACT) |
| Theorem 1.2; Thm B.1 / Lemma B.2; (3.11); Birkhoff | core / common / external | CITED / EXACT (lib) / classical |
| Seven-contact set {1,4,5,7,9,12,14} = T ∪ {12} literally = Table 5 #27, stage 5 < 6 | seven-contact | CITED |
| tie 0: T_5, T_9, T_12 or T_37 (all earlier) | paper cases | covering EXACT, CITED |
| tie 6: T_29 or T_31 | paper cases | covering EXACT, CITED |
| tie 11: T_42 | paper case | covering EXACT, CITED |
| tie 15: T_5, T_9, T_12 or T_37 | paper cases | covering EXACT, CITED |
| tie 17: T_29 or T_31 | paper cases | covering EXACT, CITED |
| tie 18: T_9, T_14 or T_40 | paper cases | covering EXACT, CITED |

For tie 12, no six-subset of T ∪ {12} lies in an earlier six-contact orbit. Its six-subsets are G-equivalent to T_87, T_76, T_67 and T_68, all of stage ≥ 6. So the seven-contact theorem really is needed, as the document says. All ties with f35 > 0 give all-minor cones (EXACT).

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | relabelling image and source pullback | identity | EXACT |
| 2 | source, sole non-minor generator, (2) in 24 and 14 coordinates; kept 8 independent with span {c : (Pc)_T = 0}; readouts; d ≤ g3 | identity | EXACT |
| 3 | extra ties and their covering, including T ∪ {12} = Table 5 #27 | finite | EXACT (covering), CITED (cases) |
| 4 | E (up to the convention sign), r, k, v formulas; H7 = 1 | identity | EXACT |
| 5 | inverse (5) recovers (x, y, t) from (z, w, r, k, v) | identity | EXACT on 20 random rational spectra, plus symbolic agreement of the homogeneous inverse with (5) |
| 6 | the denominators in (5) are nonzero for eight distinct labels | case analysis | PROOF-READ |
| 7 | E_hom = −QZKN; H3 = r/Q; F1/E = −r(Z+1)/K; the five contacts | identity | EXACT |
| 8 | every denominator factor of the 24 H_σ lies in {Q, Z, N} | identity | EXACT |
| 9 | C: bidegree (2,2), 33 terms, Z² coefficient (8), independent of h, so H14 is nonconstant at every start point | identity | EXACT |
| 10 | (9): H14 − H12 identity; H3 − H12 = rk + v − h − vkQ on C | identity | EXACT |
| 11 | Q = 0: d → 0 by (9). Q = ∞: g3 = −η Im(r/Q) → 0, so d → 0 by (4) | analytic | PROOF-READ |
| 12 | Z = 0: C(0,Q) = −(Qv+Q−r)(Qkv−rk−v+h), real Q | identity | EXACT |
| 13 | N = 0: Res_Q(C,N) as printed, real Z-roots; N's coefficients have a common zero only if 1+v−r(k+1) = 0; so Q is real | identity | EXACT |
| 14 | Z = ∞ with finite Q: Q ∈ {r/v, (1−r)/k}, real; [Z²Q²] = k²v | identity | EXACT |
| 15 | K = 0: no H_σ has K in its denominator; Z = −1 and K = 0 give Q = 1/(k+1) and N = 0; height lemma | identity + analytic | EXACT + PROOF-READ |
| 16 | genericity: (2,3,4,5) satisfies the nine listed nonvanishing conditions, and the Z² and Z⁰ coefficients are coprime (no constant-Q component) | identity | EXACT |
| 17 | Rouché approximation; Q nonconstant so it reaches ∞; minimum principle; relabel back | analytic | PROOF-READ |

Completeness of punctures. By step 8, a pole of some H_σ lies over Q = 0, Z = 0, N = 0, Z = ∞ or Q = ∞. Each case has d → 0, by steps 11-14. Points where only f diverges (K = 0) are handled by step 15. No puncture is missing.

## D. Falsification (NUMERIC)

Seed 20261125, 40000 trials, 156 supported points. On all of them, d = ηvk Im Q to 1e−7 and d ≤ g3. The worst normalized Leibniz minimum over SO(4) at 25 points (8 restarts) was +4.2e−14. No violation.

## E. Verdict: confirmed-modulo-citations

Presentation gaps only:
- G1 (p.2). The ties 0, 6, 11, 15, 17, 18 are said to "contain earlier closed six-tuples" without naming them. They are supplied above (EXACT).
- G2 (p.2). The printed sign of E belongs to the homogeneous convention. It is harmless.
- G3 (p.4). The harmonic completion is compressed, as in the sibling documents. The completeness argument above fills it in.

Citations: the seven-contact theorem on Table 5 #27; T_5, T_9, T_12, T_14, T_29, T_31, T_37, T_40, T_42 (one per tie suffices); Lemma 4.2; Theorem 1.2; Thm B.1 / Lemma B.2; Birkhoff.
