Evidence level: EXACT for all coordinate identities, the inverse, the last-equation coefficients and the fibre and limit computations (own code, 30 checks, 0 failures); PROOF-READ for the puncture analysis, the subharmonic argument and the wall (Rouché) step; CITED for the prior cases, Appendix B, the bounded-negative-part (compact pasting) lemma and Theorem 1.2. NUMERIC only for a brief falsification.

# root_inventory72_native_uniform_20261023: primary proof of paper case T_73

PDF: `proofs/manin_transfer_20261021/root_inventory72_native_uniform_20261023.pdf`, 4 pages. Title: "A six-contact support inequality from two coefficient bounds". Replay: `replay/root_inventory72_native_uniform_20261023__check.py` (about 45 s).

## Deviations and limits

- None in scope. Nothing from the supplement was executed or read. The sigma-points are rebuilt from A = (a,z,x,w), B = (t,0,1,inf) with plus brackets.
- The Appendix G stage order was not checked; it is not in the library. For each extra tie I list every paper case that covers it and check that one with a smaller index exists.

## A. Statement

**Theorem 1 (p. 1).** Let x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T = {0,4,5,7,9,12} (non-strict). Then Phi_{c(x)}(O) >= 0 for all O in SO(4).

- T is literally the paper's T_73.
- All spectra are covered. Configurations that are not eight-distinct (a_i = a_j, b_i = b_j, a_i + b_j = 0) satisfy MOC by Theorem B.1 and Lemma B.2. Those hold pointwise for every U, so the reduction to eight-distinct sources is legitimate. Extra ties go to earlier cases (section B). The wall u = 1 is restored by a perturbation argument (step 14).

## B. Dependencies

| dependency | kind | used for | check |
|---|---|---|---|
| Lemma 4.2 | core | generators [0,1,2,5,7,12,35,45], all minor except f35 = h4312,- | EXACT (lib) |
| ties 1, 6, 11, 14, 17, 18 | paper cases | an extra tie that does not force theta = 0 | covering cases computed with the group. Tie 1: T_5, T_10, T_17, T_26, T_57, T_68. Tie 6: T_24, T_26, T_37, T_44, T_65, T_67. Tie 11: T_50, T_52, T_56, T_62, T_72, T_74. Tie 14: T_51, T_52, T_63, T_87. Tie 17: T_65, T_66, T_67, T_69, T_75, T_86. Tie 18: T_50, T_51, T_55, T_61, T_62, T_63. Each has a case with index < 73 (EXACT); the stage order is CITED |
| tie 15 | seven-set | T u {15} | equals an image of E = {0,1,2,9,10,16,17} = Table 5 #17 (EXACT); the only paper cases inside are T_74 and T_80 |
| other extra ties | - | f35(j) > 0 forces theta = 0, so the source is a minor sum | EXACT: the ties with f35(j) = 0 are exactly {1,6,11,14,15,17,18} |
| Theorem B.1, Lemma B.2 | Appendix B (established) | non-eight-distinct spectra | pointwise in U; no single-collision compensation is used |
| bounded-negative-part / compact pasting lemma (Lemma 3.5 type) | common lemma | curve argument | CITED |
| Theorem 1.2 | core | unitary assertion | CITED |
| "one-circle result ... degenerate normalizations" | common lemma | unclear | not needed: the normalization B = (t,0,1,inf) only needs b_2, b_3, b_4 distinct |

## C. Steps

| # | claim (location) | kind | status | check |
|---|---|---|---|---|
| 1 | (1)-(3) p. 1: theta <= g2, g8, g10, so -L_O <= M min(g2, g8, g10) | identity + inequality | EXACT | f35 = 1 at rows 2, 8, 10; every generator >= 0 on all rows; the representation is unique (rank 8) |
| 2 | extra ties (p. 1) | citation | EXACT reductions + CITED | section B |
| 3 | (1), (2), (4) raw products; F7, F0 | identity | EXACT | sympy |
| 4 | (5), (6): J = H2, q = J - r = 1/(x+1), x and z inverse | identity | EXACT | |
| 5 | (7): N/D = a/(a+t); the inverse (8) | identity | EXACT | |
| 6 | the inverse reproduces H0 = 0, H5 = 1, H4 = 1+r, H7 = u, H9 = u+s identically; (3) H2 = J, H3 = 1 - wJ; (9) H8; (10) H10 | identity | EXACT | |
| 7 | nonvanishing of r, s, q, J, D, N, D - N at eight-distinct sources (p. 2) | inequality | PROOF-READ | each is a product of cross differences (e.g. s = a(t-1)(z-w)/E with a = a_1 + b_2) |
| 8 | the last equation C(J,w): bidegree (3,2), 60 terms, [J^3]C = (u-1)(w+1), [w^2]C = Js((1-h)J + hr) (11) | identity | EXACT | |
| 9 | punctures with J finite: g2 = -eta Im J is bounded | analytic | PROOF-READ | |
| 10 | J = inf with u != 1: only w = -1 or w = inf lie above it | algebraic | EXACT + PROOF-READ | C/J^3 -> (u-1)(w+1) on any branch with J -> inf and w finite |
| 11 | (12) at (inf, inf): H8 -> s+1, for any relative pole orders | limit | EXACT + PROOF-READ | rewritten form verified |
| 12 | (13), (14) at (inf, -1): dC~/dw = u-1 != 0 gives one smooth branch; w+1 = V0/J + O(J^-2); H10 -> u | limit | EXACT | series to second order; V0 - 1 = -r/(u-1) != 0 |
| 13 | nonconstancy (J varies, or J constant nonreal and H3 = 1 - wJ varies; real constant J gives theta = 0); pasting, a subharmonic function, constancy, contradiction | analytic | PROOF-READ + CITED lemma | boundary of the supported set: extra tie, so Phi >= 0 by step 2; a positive constant forces Im J of one sign, so J is constant, so contradiction |
| 14 | wall u = 1 (p. 4): C is never identically zero; take a complex line through the root, the zero is isolated of finite order; for real u' near 1 (r, s, h fixed) a nearby root exists (Rouché/Hurwitz); the inverse, distinctness, strict gaps and Phi < 0 persist; this contradicts the u != 1 case | analytic | PROOF-READ (argument re-derived by auditor; the document's text is adequate) | the perturbed point keeps H4, H7, H9, H12 real, so the contacts are kept |

## D. Falsification (NUMERIC)

- Six floating-point points with contact set exactly T_73 (seeds 0-5, margin 1e-2).
- After a Möbius renormalization to B = (t,0,1,inf), C(J,w) = 0 and the reality of H4, H7, H9, H12 hold to relative residual 6.5e-13. So the document's coordinates describe genuine sources.
- Minimum of Phi over SO(4) (8 restarts each): 1.0e-13. No violation.

## E. Verdict

**confirmed-modulo-citations.** The citations are the six prior cases and Table 5 #17 for the extra ties (reductions checked exactly; stage order not checked), Theorem B.1 and Lemma B.2, the compact bounded-negative-part lemma, and Theorem 1.2.

Presentation remarks:
- The appeal to the "one-circle result" is unexplained and unnecessary.
- The Rouché step does not say that the complex line must be chosen so that C restricted to it is not identically zero. This is standard.
