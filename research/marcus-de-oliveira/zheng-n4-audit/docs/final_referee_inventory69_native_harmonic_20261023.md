Evidence level: EXACT for every displayed identity, factorization, resultant, node blow-up, Newton-polygon and division identity, for the generator decomposition and the tie table (own sympy code `replay/final_referee_inventory69_native_harmonic_20261023__check.py`, all checks pass); PROOF-READ for the analytic curve and harmonic argument; CITED for earlier cases, Appendix B, Theorem 1.2 and the Birkhoff theorem; NUMERIC for a brief falsification.

# Audit of final_referee_inventory69_native_harmonic_20261023 (primary proof of paper case T_70)

PDF: `proofs/manin_transfer_20261021/final_referee_inventory69_native_harmonic_20261023.pdf`, 5 pages. Audit batch C37-T69-T71.

## Deviations and limits

- I worked from the text extract. Every formula I rely on was replayed from the definitions and agrees with the extract, so I did not open the page images.
- The supplement's companion program was not used.
- Historical labels in the tie table are read as paper index = label + 1, except T4 = T_4. These are the conventions of `deps/dependency-map.md`, and I checked the reading exactly.
- The falsification is brief. The supported set on this chart is thin, and random sampling found only 23 points.

## A. Statements

**Theorem 1 (p.1).** Let x ∈ X with P c(x) ≥ 0 and (P c(x))_π = 0 on T = {1,3,4,8,9,12}, with non-strict contacts. Then Φ_{c(x)}(O) ≥ 0 for all O ∈ SO(4), and so for every unitary by Theorem 1.2. T equals Appendix E's T_70 literally. Coverage of degenerate spectra:
- Coincidences: pointwise by Thm B.1 / Lemma B.2.
- A seventh tie: section B.
- Labels at ∞: Möbius change of chart.

So the statement is (*) for T_70 on all of X.

Internal claims:
- Source and decomposition. Source {0,1,2,5,9,53,55,56,59}, sole non-minor f56, and f56 + f5 = f0 + f9 + f55.
- Readouts. (C0, C1, C2, C9, C53, C55, C59) = (g20, g17, g15, g10, g14, g11, g2) and d = g7 − g10.
- Tie table. 2, 10, 11, 14, 15, 17, 20 go to T_13, T_42, T_29, T_44, T_4, T_44, T_29.
- Chart. E = F3 − F1 = (a+t)(w−z), the ratios r, s, u, and the inverse (1). Also L ≠ 0, uV/L = a/(a+t) and P = tuV.
- Coordinate functions (4). D = H7 − H10 = ruV/(r−W), H10 = V(s−u) + W, f = F1/E = (W−V)/(V−1). Also H6 − H11 = D − u and H18 − H20 = D + r.
- Curve. C has bidegree (2,4), 80 terms, and A = (r−1)K².
- The r = 1 derivative, E_hom = VWL(V−1)P (5), and Res_V(C,P) (6).
- W = 0 slice. The W = r node and its blow-up, with the limiting spinors, A*, and the 24 limits.
- W = ∞. The leading factorization.
- V = ∞ over W* = rs/(s−u). The barrier G* = H15 + H7 and the division identity (8).
- Height slice V = 1, and the specialization (2,3,5,7).

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| Lemma 4.2 / Prop 4.1 | core | CITED (generators EXACT via lib) |
| Theorem 1.2 | core | CITED |
| Thm B.1 / Lemma B.2 (collisions; also the shared-zero node limit A1 = B3 = ∞) | common lemma | CITED, pointwise use matches |
| (3.10)-(3.11) ker H_κ = E | core | EXACT (lib) |
| Birkhoff | external | hypotheses checked |
| T_13, T_42, T_29, T_44, T_4 (ties 2, 10, 11, 14, 15, 17, 20), all of stage < 5 = stage(T_70) | paper cases | covering EXACT, cases CITED |

Every tie σ ∉ T with f56(σ) > 0 has an all-minor cone (EXACT). The ties with f56(σ) = 0 are exactly the seven printed ones (EXACT).

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | source, sole non-minor generator, relation (1) in 24 and 14 coordinates | identity | EXACT |
| 2 | kept 8 independent, span = {c : (Pc)_T = 0}, readouts | identity | EXACT |
| 3 | tie table as printed and stage order | finite | EXACT (cases CITED) |
| 4 | ratio numerators; inverse (1); forward identities L = (w+1)(w−x)/(w−z), uV/L = a/(a+t), P = tuV; round trip | identity | EXACT |
| 5 | chart covers every eight-distinct start (V, W, L, L−uV, P ≠ 0) | case analysis | PROOF-READ (uses step 4) |
| 6 | (4) and the two other non-minor measurements | identity | EXACT |
| 7 | all 24 H_σ have denominator factors only in {V, W, W−r, P} (and r) | identity | EXACT |
| 8 | C: bidegree (2,4), 80 terms, A(W) = (r−1)K(W)², deg_W(A,B,C) = (2,3,4) | identity | EXACT |
| 9 | nonconstancy for r ≠ 1; the printed ∂_V H12 at r = 1, nonzero unless s = −1, u = 0 | identity | EXACT |
| 10 | E_hom = VWL(V−1)P | identity | EXACT |
| 11 | Res_V(C,P) as printed (6); its W-roots are real; V real from P = p1(W)V + p0(W), p1 = (s−u)W + u − rs, generically nonzero; hence D has a real limit | identity + analytic | EXACT + PROOF-READ |
| 12 | V = 0, W ≠ r: D → 0 | analytic | PROOF-READ |
| 13 | C(V,0) = r²s(Vr−V−r)(Vs−s−1), real V-roots; A(0) ≠ 0 | identity | EXACT |
| 14 | node: C(V,r) = r²u²(r−1)V²; (0,r) is a double point; the lowest-order τ² coefficient is a quadratic in Y with nonzero end coefficients; all 24 H_σ have finite limits along V = τY, W = r+τ, equal to the values of the printed limit spinors identically in Y; reference difference −rA*Y; bracket [A1, B3] = 0 (a shared zero) | identity | EXACT |
| 15 | W = ∞: the Newton polygon has the single slope 1 (max W-degree 4, 3, 2 for V⁰, V¹, V²); leading polynomial {Y(s−u)+1}{(r−1)Y(s−u)−h+r+s−u}; D → −ruY ∈ ℝ | identity + analytic | EXACT + PROOF-READ |
| 16 | V = ∞ with finite W: only W* = rs/(s−u) (double root of A); −η Im G* = C2 + C9 + d ≥ d; G* = a*V + b* + c*/(d1V + d0) with a* as printed; d1 = W(r−W)(K−u); (r−1)K²V → −B(W*); so G* → a real limit on every branch | identity + analytic | EXACT + PROOF-READ |
| 17 | height V = 1: C(1,W) = (W−1)·(cubic); H holomorphic there except where W ∈ {0, r} or P(1,W) = −(W−1)(W−r−u) = 0, generically disjoint from the cubic; height lemma as in T_69 | identity + analytic | EXACT + PROOF-READ |
| 18 | (2,3,5,7): node equation 10(10Y²+16Y+3), W* = −3, C(1,W) = −(W−1)(7W³−50W²+16W+12), coprimality, disjointness, d1(W*) ≠ 0 | identity | EXACT |
| 19 | Rouché approximation; harmonic completion | analytic | PROOF-READ |

Note on steps 7 and 11-17: these settle completeness of the punctures. A pole of some H_σ on the compact curve lies over V = 0, W = 0, W = r, P = 0, V = ∞ or W = ∞. Each of these is listed, and at each d → 0, or G* → ℝ, or, at the node, all H_σ stay finite. At the node, L_O therefore extends harmonically, and a boundary point of S there lies in X because X is closed. The printed shared-zero identification is correct (EXACT), but it is not even needed. Points where only f diverges (V = 1) are handled by the height lemma. I found no unlisted puncture.

## D. Falsification (NUMERIC)

Seed 20261125, 40000 trials on the chart, 23 supported points (η = ±1). On all of them, d = −η Im(ruV/(r−W)) to 1e−7, and d ≤ g15 + g7. The worst normalized Leibniz minimum over SO(4) (8 restarts) was +3.8e−14. No violation.

## E. Verdict: confirmed-modulo-citations

Presentation gaps only:
- G1 (p.4-5). The harmonic completion is compressed, as in T_69. The separation of poles of the H_σ from points where only the inverse degenerates is supplied in the note above.
- G2 (p.2). The genericity conditions ("r ≠ 1, s ≠ u", etc.) are not listed in full. They are proper algebraic conditions, and one exact specialization satisfies all of those used (EXACT).

Citations: T_4, T_13, T_29, T_42, T_44; Lemma 4.2; Theorem 1.2; Thm B.1 / Lemma B.2; Birkhoff.
