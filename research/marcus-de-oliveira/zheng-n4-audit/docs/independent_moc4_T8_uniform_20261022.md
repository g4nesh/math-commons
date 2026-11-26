Evidence level: EXACT for all displayed identities of Sections 2-3 except Lemma 4's (8) and the -N identity (UNCHECKED); EXACT-randomized (exact rational arithmetic at random points, Schwartz-Zippel) for the Hessian-determinant identity that defines G; PROOF-READ for the sign arguments of Lemmas 2-3 and the complexified-cubic lemma; Section 1 as in the polar-reduction report; CITED for Appendix B, SOC (f1,f9,f18), T_4..T_7, T^(7)_1.

# independent_moc4_T8_uniform_20261022 (8 pages). Primary proof of paper T_9 = {0,1,2,3,11,12}

## Deviations and limits

- Section 1 has the same text as `independent_moc4_T8_polar_reduction_20261022`. See that report and its replay (57 EXACT checks).
- Sections 2-3 are replayed by `replay/independent_moc4_T8_uniform_20261022__check.py` (33 checks, all pass), which imports the chart helpers of the double-boundary replay.
- **G is not printed in the document.** I reconstructed it from the stated identity Q2^4 det Hess C(c0 - (aQ1/Q2)κ) = -4a^4 Q1^2 Q2^2 G^2. Method: interpolation mod p = 2^61-1 along 820 lines through one anchor, where the ratio was a perfect square of a quintic on every line. The identity was then checked in exact rational arithmetic at 25 random rational points with h ≠ 0. This is a randomized identity test, not a symbolic expansion. Every claim about G (independence of h, affine in B, coefficient signs, the J = 0 identity, Lemma 2's J-quadratic data) was then checked EXACTLY on the reconstructed G. If the reconstruction were wrong, the randomized test would have failed with overwhelming probability.
- Lemma 4's Vieta–norm identities (8) and the "-N" positive identity were **not replayed** (time). I tested Lemma 4's conclusion numerically instead. On 990 chart points of T_9∪{18} in the branch 0 < x, t < 1, r = -t-R, 0 < M < R(1-x) (seed 4), none was supported. The best normalized min gap was -0.51.
- The SOC lemma for (f1,f9,f18) is CITED. It was tested numerically in the double-boundary replay (min ≈ 7e-6 >= 0, 20000 samples).

## A. Statements

- **Theorem 1** (Section 1) = the polar reduction: T_9 reduces to T_9∪{18} or T_9∪{20}.
- **Theorem 2.** For every x in X with P c(x) >= 0 and contacts ⊇ T_9∪{20} (= G-image of Table 5 T^(7)_4): Phi >= 0 on SO(4). The proof does not depend on the spectral chart. It uses only X ⊂ Sing C, so it covers degenerate spectra as well.
- **Lemma 1.** For x in X ∩ N_{T_9} with A = g4 - g18 > 0: G(A,B,C,Q,J,K,L) = 0, where G = B G_B + G_0, and G_0 < 0 on {A, C, K > 0, CK < A^2, Q, J, L >= 0}. This holds for all of N_{T_9}, not only on J = 0.
- **Theorem 3.** The same conclusion for contacts ⊇ T_9∪{18} (= T^(7)_3). Eight distinct points via the chart. Degenerate spectra via Appendix B.
- **Lemmas 2, 3, 4.** In the negative region (A, C, K > 0, CK < A^2) of T_9∪{18}, the branches x < 0, then (0 < x < 1, t < 0), then (0 < x, t < 1) are impossible.
- Combined (p.5-8): Phi >= 0 for every x in X with contacts ⊇ T_9 and every O in SO(4). This is (*) for T_9 (the same set, g = identity), with non-strict contacts and all spectra (degenerate ones through Appendix B).

## B. Dependencies

| ref | kind | status |
|---|---|---|
| X ⊂ Sing C (paper A.1; lib-verified) | core | EXACT in lib |
| complexified-cubic lemma (∇C(c) = ∇C(d), H(c)d = 0 ⇒ C(c) = 0, det H(c) = 0) | — | PROOF-READ (checked line by line; Euler identity EXACT) |
| SOC (f1,f9,f18) | common-lemma | CITED, NUMERIC support |
| T_4 (C = 0 or K = 0 in Thms 2/3) | paper-case | containment EXACT for T_9∪{20,7}, {20,23}, {18,7}, {18,23}; case CITED (Stage 1) |
| T^(7)_1, T_4..T_7, Appendix B | via Section 1 | CITED |
| Appendix B (repeated points and cross-collisions in Thm 3) | common-lemma | CITED |
| two-orthogonal theorem (Thm 1.2) | core | CITED (only for the unitary remark) |

## C. Steps (Sections 2-3; Section 1 as in the polar report)

| # | claim | status |
|---|---|---|
| 2.1 | profile identity with f1,f2,f7,f8,f9,f12,f15,f18; A <= 0 ⇒ minor | EXACT (polar replay, auditor's f0/f10 decomposition) |
| 2.2 | c0 = P^{-1} g (zero-height lift); C(c0 + zκ) = aQ1 + zQ2 (affine in z); Q1 = be - cd + df = BQ + (A+L)(B+Q+A+J) | EXACT |
| 2.3 | complex matching vector ⇒ C(Re x) = 0 and det H(Re x) = 0; Re x = c0 + zκ; Q2 ≠ 0 since aQ1 > 0 | PROOF-READ |
| 2.4 | Q2^4 det H(c0 - (aQ1/Q2)κ) = -4a^4 Q1^2 Q2^2 G^2 | EXACT-randomized (see limits) |
| 2.5 | G independent of h, affine in B | EXACT (on the reconstructed G) |
| 2.6 | (C+D)^2 G_0(K = A^2/(C+D)): all 81 coefficients negative; value -A^4 D(A+C)(A+C+D) at Q=J=L=0 | EXACT |
| 2.7 | -G_B at J = 0 equals the displayed sum of nonnegative terms | EXACT |
| 2.8 | conclusion CK >= A^2, then SOC | PROOF-READ + CITED |
| 3.1 | chart: inverse formulas; (3), (4), (5) | EXACT |
| 3.2 | x = (Q-C)/(A+Q+J), t = (Q-K)/(A+Q) | EXACT (they hold as x(A+Q+J) - (Q-C) = x·g11 and t(A+Q) - (Q-K) = t·g11, so the contact 11 is needed, and it is present) |
| 3.3 | signs sgn ρ = sgn B = sgn x, r opposite, -PH/B > 0 (identity F/C - 1 = -PH/B, checked by hand), tM > 0 | PROOF-READ |
| 3.4 | Lemma 2: rR2 + PFK = 0 | EXACT |
| 3.5 | G_B quadratic in J, leading (K-Q)(A+L); (A+L)^2 G_B(J*) = -A(A+L)(A+Q+L)F(A+K)^2; linear coefficient = -H1 identically | EXACT |
| 3.6 | t > 0 case: "C > Q > K and A > K". A > K is not justified in the text. It follows from CK < A^2 and C > K (if A <= K < C then CK > A^2) | PROOF-READ (argument supplied by auditor) |
| 3.7 | H9 = 1 + H11 + H11/(tv-1) ⇒ Q = ρtu/\|tv-1\|^2 | EXACT |
| 3.8 | Π real cubic; (6), (7) | EXACT |
| 3.9 | Lemma 3 sign chain; z - c0 formula; H11(w) formula | EXACT identities + PROOF-READ signs |
| 3.10 | Lemma 4: (8) and the -N positive identity; I = (x-tw)D/[w(w-1)Zn] | UNCHECKED (numeric support only) |

## D. Falsification

- Theorem 3 CK >= A^2: random chart points of T_9∪{18}, seed 2026. Only 6 supported points with A > 0 were found. Min (CK - A^2)/A^2 = 6.7. This is weak, because the supported set is thin under this sampler.
- Lemma 4 branch: 990 points, 0 supported (seed 4).
- SOC: see the double-boundary report.

## E. Verdict

**confirmed-modulo-citations**, with two caveats:
1. The defining determinant identity for G (Lemma 1) was checked by an exact-arithmetic randomized test, not by symbolic expansion.
2. Lemma 4's identities (8) and -N were not replayed; only numeric support exists for them.

Citations: Appendix B, SOC (f1,f9,f18), T_4..T_7, T^(7)_1, and Theorem 1.2 for the unitary remark only.

Presentation gaps:
- p.7: "A > K" is not justified in the text (auditor's one-line argument above).
- p.5: G is never displayed. The proof rests on a companion script, and I reconstructed G myself.
- p.3: the labels T2/T4/T5 for the endpoints are inconsistent with the paper's indices (see the polar report).
- p.4: Theorem 2's remark that degenerate spectra need the global boundary theorems is unnecessary, because Lemma 1 holds on all of X.

The statement (*) for T_9 is proved at least as strongly as needed (non-strict contacts, all O in SO(4), all spectra modulo Appendix B), provided Lemma 4 holds as stated.
