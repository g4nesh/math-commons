Evidence level: EXACT for all coordinate identities, the spinor inverse, the 144-term last equation, the fibre and branch computations and the K-face radial polynomial (own code, 32 checks here plus the G9 face check in the T_72 replay, 0 failures); PROOF-READ for the puncture analysis, the closure argument, the subharmonic argument and the wall step; CITED for the prior cases, Appendix B, the compact pasting lemma and Theorem 1.2. NUMERIC only for a brief falsification.

# final_referee_inventory73_native_uniform_20261023: primary proof of paper case T_74

PDF: `proofs/manin_transfer_20261021/final_referee_inventory73_native_uniform_20261023.pdf`, 4 pages. Title: "A ramified coefficient pole and a coercive boundary identity". Replay: `replay/final_referee_inventory73_native_uniform_20261023__check.py` (about 75 s). Eq. (17) is replayed in `replay/final_referee_inventory71_uniform_20261023__check.py`, where G9 is rebuilt.

## Deviations and limits

- None in scope. Nothing from the supplement was executed or read.
- Convention: with minus brackets and the infinite factor dropped, E = F5 - F3 equals -x(y-t)(z-w). The document prints +x(y-t)(z-w). This is the sign of the dropped constant, it does not affect any H_i, and it is a presentation point only.
- The overall sign of C was normalized to [u^5]C = +Z^4, as in the document.
- The left-inverse property of the spinor inverse was checked exactly at 4 random rational sources rather than symbolically. The symbolic identities (5)-(7) already show that the source is determined by (u, Z, r, s, p).

## A. Statement

**Theorem 1 (p. 1).** Let x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T = {3,4,5,7,8,12} (non-strict). Then Phi_{c(x)}(O) >= 0 for all O in SO(4).

- T is literally the paper's T_74.
- All spectra are covered. Non-eight-distinct spectra go to Appendix B (B.1, B.2), extra ties to earlier cases, and the generic restrictions (10) are restored by perturbation.

## B. Dependencies

| dependency | kind | used for | check |
|---|---|---|---|
| Lemma 4.2 | core | generators [0,1,2,5,6,59,67,70]; identity rows [18,16,13,10,9,2,14,6]; f70 = h4231,+ is the only non-minor | EXACT |
| ties 2, 9, 10, 13, 14, 18 | paper cases | extra ties with f70(j) = 0 | covering cases computed. Tie 2: T_10, T_13, T_26, T_27. Ties 9 and 18: T_50, T_52, T_56, T_62, T_72, T_73. Ties 10 and 13: T_24, T_25, T_26, T_48, T_66, T_75. Tie 14: T_73, T_80 (or Table 5 #17). Each has a case with index < 74 (EXACT); stage order CITED |
| tie 16 | seven-set K = Table 5 #16 = {0,1,2,9,10,14,17} | T u {16} | the stated action (rows identity, columns (0,3,2,1), no transpose) maps K onto T u {16} and the profiles [0,2,12,13,15,21,26] to [0,2,5,6,59,67,70] in order (EXACT). K = T' u {9}, where T' is the working set of the T_72 document; T_72 itself also covers tie 16 (EXACT) |
| T_72 document's polynomial G9 | other-doc (audited in this batch, confirmed) | (16)-(17) | G9 at a = b = t = 0, d = 1 equals ((q+1)(v+1) + c(q+v+1))(q+1)^2(v+1)^2(w+1)^3 (EXACT); C(g + z kappa) = C + zQ gives Q_K = 1, C_K = (q+1)(v+1)(w+1) |
| radial necessity (paper p. 27) | established | closure argument | applies to g^{-1} x_n in X directly, so no semi-invariance of C is needed; kappa is G-fixed |
| Theorem B.1, Lemma B.2 | Appendix B | collisions | pointwise |
| compact pasting / bounded-negative-part lemma | common lemma | curve argument | CITED |
| Theorem 1.2 | core | unitary | CITED |

## C. Steps

| # | claim (location) | kind | status | check |
|---|---|---|---|---|
| 1 | (1), (2): canonical coefficients are the gaps; d = g6 <= g0, g1; -L_O <= M d | identity | EXACT | |
| 2 | perturb O so that all seven minor values are > 0; (4): 0 <= x_i/d <= -B/f_i(O) | analytic | PROOF-READ | Phi < 0 is open in O; each minor squared is a nonzero analytic function on SO(4) |
| 3 | (5), (6) inverse with the raw residual; (7) Y = Z(p-1)/J, H4 formula, p != 1 | identity | EXACT | |
| 4 | spinor inverse A0..C3: H3 = 0, H5 = 1, H7 = r, H8 = s, H4 = p, H10 = u identically | identity | EXACT | also a left inverse at 4 random rational sources (all 24 H_i) |
| 5 | (8): C(u,Z) has degrees (5,4), 144 terms, [u^5]C = Z^4 | identity | EXACT | |
| 6 | (9) formulas for H6 and H1 | identity | EXACT | |
| 7 | finite punctures: L != 0 gives H6 regular; L = 0 with s != 0 gives Z != 0 and H1 -> 0, so d <= g1 is bounded | analytic | PROOF-READ | |
| 8 | (11): the Z = inf fibre is u in {p, 1, r}; at u = r, p, H1 -> (p-1)u/(u-1) | algebraic + limit | EXACT + PROOF-READ | the limit holds for any approach to (u0, inf) with u0 != 1 |
| 9 | u = inf: only Z = 0, and (inf, inf) is not on the closure; (12) the lowest form has degree 4, with value x^4 at t = 0; so uZ and L are bounded; the |L/Z| disjunction bounds H6 or H1 | algebraic + analytic | EXACT + PROOF-READ | tangent lines x = lambda t with lambda finite |
| 10 | (13), (14): one smooth branch at (u,Z) = (1, inf), t = (r-1)/(p-1) x^2 + O(x^3) | analytic | EXACT | the t-derivative (r-1)(p-1)^2 is nonzero under (10) |
| 11 | (15): H13, H16, H10 -> 1 on that branch; H6 has a simple pole | limit | EXACT | series to x^4 |
| 12 | (16), (17) closure argument: if d is unbounded, the normalized limit is a K-source with a = b = 0, d = 1, and z_n -> -C_K; F = 0 and det HF = 0 pass to the limit, but G_K > 0. So d is bounded | analytic | PROOF-READ + EXACT identities | Q_n -> 1 != 0; the limit has det Hess = -4 G_K^2/Q^2 != 0, contradiction |
| 13 | nonconstancy: u varies, or u is constant nonreal and H6 has linear growth coefficient (r-u)(u-1)/u != 0; then the pasting lemma | analytic | EXACT + PROOF-READ + CITED lemma | |
| 14 | walls (10): [u^5]C = Z^4, so C is never identically zero; Rouché gives nearby roots for perturbed real (r, s, p, h) satisfying (10); all strict conditions persist | analytic | PROOF-READ | the same argument as in T_73 |

## D. Falsification (NUMERIC)

- Six floating-point points with contact set exactly T_74 (seeds 0-5, margin 1e-2).
- Minimum of Phi over SO(4) (8 restarts each): 2.1e-14. No violation.

## E. Verdict

**confirmed-modulo-citations.** The citations are the prior cases and K for the extra ties (reductions checked exactly; stage order not checked), Appendix B, the compact pasting lemma and Theorem 1.2. The K-face polynomial was replayed exactly from G9, whose identity (9) is proved exactly in the T_72 replay.

Presentation remarks:
- The sign of E (see Deviations).
- "Finite intersection is dense" in step 2 means the complement of finitely many proper analytic subsets. That is fine.
