Evidence level: EXACT for the source decomposition, zero faces, restricted cubic (7), identity (8), radial identity (9) and every positivity certificate of Lemma 1 and Section 3.1 (own code, 48 checks, 0 failures). PROOF-READ for the region bookkeeping (10)-(14) and the logic. CITED for the capacities (2)-(6) (audit foundations, so4-atoms part 1 and 2b), radial necessity (audit item 5), Theorem 1.2, and the paper cases T_7, T_15, T_34. NUMERIC only as supporting evidence.

# independent_moc4_T45_uniform_20261022 ("The support inequality for T45 by product-domain polynomial identities", 6 pp.)

Primary proof of paper case **T_46 = {3, 4, 5, 9, 10, 12}**. The supplement's label is historical T45.

Replay: `replay/independent_moc4_T45_uniform_20261022__check.py`, output in `__check_output.txt`. The script runs in stages: `setup`, `direct LO HI`, `square`, `numeric`. Set MOC4_CACHE to reuse the G9 that the script itself computed.

## Deviations and limits

1. I did not use the supplement's verifier, its JSON, or its 1158-term G9. I rebuilt G9 myself as the 1158-term squared factor of the full 14x14 Hessian determinant of the paper's cubic (A.2) on the face. It agrees with the document's description: 1158 terms, total degree 9, and identity (9) holds exactly.
2. In Lemma 1 I cleared denominators my own way. For the 17 domains without substitution (14), my coefficient counts equal the document's table exactly, for example 7229, 2719, 13644, 20860, 36870. For the nine domains with (14) the counts differ, because the document first strips an M-monomial (for example 1443 vs 1258, 1299921 vs 499377). The claim being certified is the same: nonnegative coefficients plus a positive monomial in the strictly positive variables.
3. Section 3.1: I replayed the Bernstein certificate in an equivalent form. On [0, L] I substitute Z = L u/(1+u), so the coefficients are the Bernstein coefficients up to positive binomial factors. On the PR side my four coefficient lists equal the document's exactly (4140, 7956, 11076, 13139, 13756 and 5400, ..., 12912). On the PB side I removed an M^1 factor, which swaps parities, so one list matches (792, ..., 13756) and the other is an equivalent clearing.
4. The capacities (2)-(6) are not re-proved here. The audit foundations prove them: SOC(2,59;70), chain(2;0,12;70), chain(59;0,1;70), chain(59;67,72;70), and the old-H split with H(1,59;70) and H(72,59;70), all in the proved orbits with the roles as printed (so4-atoms-part2b.md, G2). The literal (row, col) old-H table of p. 2 was not reproduced (foundations note G2/limit 5). This is a notation point, not a gap.
5. Process: some first exploratory runs were in my scratch area; every result reported here comes from the replay script.

## A. Statements

- **Theorem 1 (thm:whole, p. 1).** Let x in X, pi* any contact, c = c(x) with P c >= 0 and (P c)_pi = 0 for pi in {3, 4, 5, 9, 10, 12}, with extra contacts allowed. Then Phi_c(O) >= 0 for every O in O(4), and through Theorem 1.2 the support inequality holds for every U in U(4). The proof works only with the generator coefficients and the radial necessity, which holds for every x in X. So all spectra are covered, including repeated eigenvalues and a_i + b_j = 0. The contact set is literally the paper's T_46 (identity group element), and the statement is at least as strong as (*).
- **Lemma 1 (lem:26, p. 4).** Normalize d = 1 and assume a, c, r > 0. In the 26 sign domains of (x-1, y-1, z-1) other than (1,1,-1), with the charts (11)-(14), G9 times a positive denominator is a polynomial with nonnegative coefficients and a positive monomial in the strictly positive variables. Hence G9 > 0 on those domains, including p, b, q, s = 0 and x, y, z = 1.
- **Section 3.1 (unlabelled).** G9 > 0 on the remaining domain x > 1, y > 1, z < 1 of the failure region (10).

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| Lemma 4.2 / generators(T) = {f0,f1,f2,f12,f59,f67,f70,f72} | core | EXACT (lib) | matches |
| SOC, chain x3, old-H split (2)-(6) | common lemma | CITED, foundations so4-atoms-part2b G2 | orbit and roles checked there |
| zero face a=0: T ∪ {13} contains g.T_34 (historical T33) | paper case | CITED T_34 (Stage 1) | EXACT: g = ((1,3,2,0),(1,3,2,0),1); printed image {3,5,9,10,12,13} |
| zero face c=0: T ∪ {2} contains g.T_7 (historical T6) | paper case | CITED T_7 (Stage 1) | EXACT, image {2,3,4,5,9,10} |
| zero face r=0: T ∪ {16} contains g.T_15 (historical T14) | paper case | CITED T_15 (Stage 1) | EXACT, image {3,4,5,9,10,16} |
| radial necessity: C(u) = 0, det Hess C(u) = 0 for u = c + Re(F_pi*) kappa | common lemma | CITED (audit item 5); the document's own polarization argument p. 2-3 is also correct | - |
| Theorem 1.2 (O(4) to U(4)) | core | CITED (audit foundations) | - |

T_46 is in Stage 2 of Appendix G and every cited case is in Stage 1, so the order is acyclic.

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | nullspace has dimension 8; the eight profiles evaluated at (18,16,13,8,2,14,11,7) form the identity, so coefficients = gaps >= 0 | identity | EXACT | lib values |
| 2 | d = 0 gives a sum of minors | trivial | PROOF-READ | only f70 is non-minor (EXACT) |
| 3 | capacities (2)-(6) imply the support inequality | inequality | CITED | foundations G2 |
| 4 | a = 0, c = 0, r = 0 add contacts 13, 2, 16, and each face contains an image of T_34, T_7, T_15 | citation | EXACT + CITED | own G search |
| 5 | (7) C(g + h kappa) = C + hQ with the printed C and Q | identity | EXACT | flint, paper cubic (A.2) |
| 6 | (8) RC = (cR - dR - d^2) rs - P, and P > 0 for d, r > 0 | identity | EXACT | P has positive coefficients and contains d^3 r |
| 7 | if (6) fails then C < 0, so Q != 0 and h = -C/Q | logic | PROOF-READ | uses R = r + s > 0 |
| 8 | (9) Q^4 det Hess C(g - (C/Q) kappa) = -4 Q^2 G9^2, G9 with 1158 terms | identity | EXACT | Bareiss over Z[p,r,a,b,c,q,d,s,h], then factorization |
| 9 | failure of all of (2)-(6) gives exactly region (10) | case split | PROOF-READ | including the "or p = b = 0" clauses |
| 10 | charts (11)-(14) cover region (10) bijectively, with positive denominators | parametrization | PROOF-READ + EXACT (denominators certified positive in code) | h = T/(1+T) gives ac < 1; bounds (13) are the failed chains; (14) is the failed old-H when y, z >= 1; elsewhere M is free (superset) |
| 11 | Lemma 1, 26 domains | certificate | EXACT | all 26 pass: nonnegative coefficients and a strict monomial |
| 12 | (15) b/U <= r/V iff M^2 <= L | identity | EXACT | b/U - r/V = (M^2 - L)(1+T)(1+U)/(TUM) |
| 13 | Section 3.1 Bernstein certificate on both sides | certificate | EXACT | all lists nonnegative; the last list has a positive monomial in T, U, V |
| 14 | conclusion: G9 > 0 contradicts G9 = 0 | logic | PROOF-READ | - |

Observation (not a gap): G9 vanishes identically on {r = d = 0}. The proof never uses G9 there, because r > 0 and d > 0 have been imposed by then.

## D. Falsification

- NUMERIC necessity check, seed 7. At 10 matching vectors whose gaps vanish on T, |G9|/sum|terms| <= 1.9e-14; at random coefficients the ratio is >= 0.26. Five further solver points had r or d near 0. These were skipped because G9 vanishes there identically; in the one I inspected, a_2 + b_4 = 0.
- NUMERIC main statement. At 4 lib stratum points (seeds 0-3, 15 BFGS restarts), the worst min Phi / max gap is 9.0e-17. No violation. As the brief says, this test is weak.
- The capacities were attacked in the foundations (G2: minimum 1.7e-18 for (6)).

## E. Verdict

**confirmed-modulo-citations.** Every certificate claimed in Sections 2-3 holds exactly: (7), (8), (9), the 26 domains, and the square interval. The case logic is complete. Citations: the capacities (2)-(6) (foundations), radial necessity, Theorem 1.2, and paper cases T_7, T_15, T_34.

Gaps: none of substance. Presentation: the literal (row, col) old-H action table on p. 2 could not be reproduced under the conventions tried (foundations G2). Orbit membership with the printed roles holds, so the step stands.
