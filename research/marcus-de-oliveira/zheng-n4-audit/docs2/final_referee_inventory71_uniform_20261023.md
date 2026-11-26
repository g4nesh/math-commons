# Second-pass audit: final_referee_inventory71_uniform_20261023 (paper case T_72)

Evidence level: EXACT for the transport, the generator representation, the signed side, the lifted cubic (8), the Hessian identity (9) (with G9 re-derived independently as an exact square root), the symmetry (10), all 17 positivity tables (17)-(19) with their strict markers, and the chart inverses. PROOF-READ for the allocation arithmetic, the failure conditions (6)-(7), the polarization necessity, the chart coverage and the conclusion. CITED, with NUMERIC support, for the three "reference inequalities" (sharp SOC, old-H, ordinary chain), which the document uses but neither states in full nor proves.

Auditor label: C38-T72-T74-second (independent second pass; the first auditor's report was not opened).
Document: proofs/manin_transfer_20261021/final_referee_inventory71_uniform_20261023.pdf, 5 pages, read in full (text extract, checked against the displayed formulas).
Replay: replay2/final_referee_inventory71_uniform_20261023__check.py (modes `exact` and `numeric`); output in replay2/final_referee_inventory71_uniform_20261023__check_output.txt.

## Deviations and limits

- None of the supplement's scripts or JSON files was executed, imported or read. Every polynomial (C, Q, the 14x14 Hessian determinant, G9 and the 17 cleared tables) was rebuilt from the paper's definitions (Appendix D generators, cubic (A.2)) through lib/moc4lib.py. G9 was not copied from anywhere. It was obtained as the exact square root of Q^4 det H / (-4Q^2), so (9) is checked as an identity and not just evaluated.
- The three reference atoms were **not proved**. Section D reduces each one to a sharp scalar inequality on SO(4) and tests it numerically. They are CITED dependencies, and the document does not name the documents that prove them.
- My cleared tables use my own numerator/denominator representation of each chart. The term counts nevertheless agree exactly with the paper's table (18) and with the counts for (19), 496798 coefficients in all. This suggests, but does not prove, that I used the same representation.
- The batch prompt file for this run contains a line saying second passes are paused. The task prompt assigned this audit, so I carried it out and report the line here for the orchestrator.

## A. Statements

**Theorem 1 (p. 1).** Let T = {1,3,4,7,8,12} (zero-based lexicographic indices), which is paper T_72 exactly (EXACT: equals CONTACT_SETS[72]). The claim: for every x in X (all spectra a, b in C^4 and alpha, including repeated eigenvalues, a_i + b_j = 0 and extra contacts) with P c(x) >= 0 and (P c(x))_pi = 0 for all pi in T, Phi_{c(x)}(O) >= 0 for every O in O(4) (in particular every O in SO(4)). This is at least as strong as (*) for T_72. Contacts are non-strict ("contain T"), and no genericity is assumed anywhere in the proof (see C). The added sentence about the unitary assertion is a remark that rests on the paper's Theorem 1.2. It is not needed for (*).

Internal claims (all for the working set T' = {0,1,2,10,14,17}):
- (1) Every supported cost has the unique form g = q f0 + a f2 + b f12 + v f13 + c f15 + w f21 + t f23 + d f26. Here q, a, b, v, c, w, t are the gaps at 23, 15, 8, 7, 4, 16, 9, and d = g11 - g9.
- (2) f10 + f26 = f0 + f13 + f23. Also g20 = q+d, g12 = v+d, g11 = t+d.
- (3)-(5) Sufficient conditions from the atoms: d <= H3(a,b,c); t(c-d) >= d^2; and, for w, t > 0, c >= d^2/w + d^2/t.
- (8) C(g + z kappa) = C + zQ with C = (q+d)(v+d)(w+t+d) + abc and Q = ab+ac+bc+d^2 - t(c-d).
- (9) Q^4 det Hess C(g - C kappa/Q) = -4 Q^2 G9^2. G9 is homogeneous of degree 9 with 1615 terms and [d^9] = 1.
- (10) G9 is invariant under (q,a,b,v) -> (v,b,a,q).
- (17)-(19) G9 > 0 on the strict-failure domain (6)-(7), by 15 + 2 positive cleared tables.

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| G-action (rho, chi) = (0213, 0231), no transpose, maps T' to T and is applied to spectra and O together; Phi is invariant and matching vectors go to matching vectors | core (paper Sect. 4.2) | EXACT for the set map (lib); the equivariance itself is a lib EXACT test of paper (4.7) | `g_act_set` gives T exactly; the action permutes all 82 f_j |
| X is contained in Sing C (grad C = 0 on X, X closed under conjugation) | core (paper App. A.1) | EXACT in lib (16-variable identity) | used in the polarization step |
| Phi_{n_{I,J}}(O) = det O[I,J]^2 >= 0 | core (paper (4.2)) | EXACT/PROOF-READ in lib | minor remainders are nonnegative |
| Sharp SOC atom (f1, f9; f22) = (n12,34; n14,23; h3412,+): x f1 + y f9 + z f22 is safe when z^2 <= xy | other-doc (unnamed "reference inequality") | CITED + NUMERIC | images onto (f2,f12;f26), (f2,f15;f26), (f12,f15;f26) EXACT (4 group elements each). Image inequality Phi_f26 >= -2 abs(m_i m_j) holds numerically (min ~ 0, sharp) |
| Old-H atom (f3, f15; f26): x f3 + y f15 + z f26 is safe when x(y - z) >= z^2 | other-doc (unnamed) | CITED + NUMERIC | image (f23, f15; f26) with f3 -> f23 and f15 -> f15 EXACT. Phi_f26 + m15^2 + 2 abs(m23 m15) >= 0 numerically (min ~ 0) |
| Ordinary chain atom (center f2; leaves f0, f12; f70) | other-doc (unnamed) | CITED + NUMERIC | image (center f15; leaves f21, f23; f26) role by role EXACT (2 group elements). Phi_f26 + 2 abs(m15)(abs(m21) + abs(m23)) >= 0 numerically (min ~ 0) |
| Theorem 1.2 (unitary to orthogonal) | core (paper) | CITED | only for the remark about the unitary assertion |

No other contact case, branch document, Appendix B/C result or external theorem is used. Lemma 4.2 is not needed either, because representation (1) follows from linear algebra on the 8-dimensional contact kernel.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | The literal action sends T' to T; proving the claim on T' suffices | citation + finite | EXACT (set map); equivariance CITED to paper 4.2 / lib test | `g_act_set`, permutation of the 82 f_j |
| 2 | generators(T') = {0,2,10,12,13,15,21,23,26}; only f26 = h4132,- is non-minor | finite | EXACT | lib `generators`, `F_IS_MINOR` |
| 3 | The contact rows have rank 6; the 8 profiles of (1) are independent and therefore span the kernel; the coefficients are the gaps at 23,15,8,7,4,16,9; d = g11 - g9 | identity | EXACT | rank computations, value table |
| 4 | (2): f10 + f26 = f0 + f13 + f23; g20 = q+d, g12 = v+d, g11 = t+d | identity | EXACT | 24 values |
| 5 | If d <= 0, then g = (q+d)f0 + a f2 + b f12 + (v+d)f13 + c f15 + w f21 + (t+d)f23 - d f10, with all coefficients gaps or -d >= 0 and all profiles minor, so Phi >= 0 | identity + inequality | EXACT + PROOF-READ | from step 4 and (4.2) |
| 6 | Image claims for the three atoms | finite | EXACT | search over G for role-preserving images |
| 7 | The atom inequalities themselves | inequality on SO(4) | CITED (NUMERIC support) | Section D |
| 8 | Allocation: balanced case, pair coefficients (a+b-c)/2 etc. give total (a+b+c)/2; dominant case, a proportional split gives sqrt(m(s-m)). So d <= H3 is sufficient (smaller nonminor coefficients keep z^2 <= xy) | arithmetic | PROOF-READ | sum over pairs of sqrt(x y) computed by hand |
| 9 | (4) and (5) follow from the atom images with the unused coefficients as minor remainders | arithmetic | PROOF-READ | |
| 10 | Failure of (3)-(5) with d = 1 gives (6) H3 < 1 and t(c-1) < 1, and (7) 0 <= w < t/(ct-1) when ct > 1 (w = 0 included; t > 0 automatic) | arithmetic | PROOF-READ | |
| 11 | (8): C(g + z kappa) = C + zQ, linear in z | identity | EXACT | fmpz_mpoly, paper cubic (A.2) |
| 12 | C >= d^3 > 0, and under (6) Q > ab+ac+bc >= 0, so z = -C/Q | inequality | PROOF-READ | |
| 13 | For realizable g, Re x = g + z kappa with z = Re F_pi*; Re x = (x + conj x)/2 with x, conj x in X. Polarization gives C(Re x) = 0 and Hess C(Re x)(Im x) = 0; if Im x = 0, Euler gives Hess C(x)x = 2 grad C(x) = 0; x = 0 is trivial. Hence det Hess C(Re x) = 0 at every realizable point, degenerate ones included | analytic (algebraic) | PROOF-READ | uses only grad C = 0 on X (lib EXACT) |
| 14 | (9): Q^4 det H(g - C kappa/Q) = -4 Q^2 G9^2, G9 with 1615 terms, degree 9, [d^9] = 1, degree profile (3,2,2,3,2,3,2,9) | identity | EXACT | 14x14 Bareiss determinant of Hess C(g + z kappa) (27680 terms, z-degree 3). Then z = -C/Q, exact division by -4Q^2 and exact square root |
| 15 | So at a negative-response realizable source (d = 1, conditions (6)-(7)) G9 = 0 | logic | PROOF-READ | Q != 0 by step 12 |
| 16 | (10) symmetry of G9, C, Q and of the failure conditions | identity | EXACT (G9, C, Q) + PROOF-READ (conditions symmetric in a,b, with q,v free) | compose |
| 17 | The charts (11)-(14), the solo chart, the interfaces (15)-(16) and the c = 0 face (19) cover the whole strict-failure domain | case split | PROOF-READ + EXACT (the inverse maps recover the parameters) | balanced/dominant split. c > 1 forces dominance; c = 1 balanced is impossible; c = S = 1 is safe; ct = 1 has its own row; for (14) the old-H bound gives tau = ct - 1 in (0, 1/X) |
| 18 | Each of the 17 cleared polynomials has only positive coefficients; term counts 28693, 10012, 93570, 27186, 10665, 70122, 36773, 13175, 70650, 38081, 12700, 79550, 703, 271, 1674, 1293, 1680 (total 496798) | finite | EXACT | own substitution with multiplicities (3,2,2,3,2,3,2) and d = 1; all counts equal the document's |
| 19 | Strict positivity: a constant term, or marker monomials X^2, Y^2 / T^3X^3, T^3Y^3 / X / T^3 / X, Y / T^5 present, each positive on its chart's domain; the denominators (products of 1+..., X > 0 in (14), T > 0 in high rows, X+Y > 0 in (11)) are positive | finite + logic | EXACT + PROOF-READ | monomial membership |
| 20 | Conclusion: G9 > 0 on the failure domain contradicts step 15, so no realizable source has Phi < 0 | logic | PROOF-READ | |

## D. Falsification

- Sharp forms of the atoms. Write m_j = det O[I_j, J_j] for the minor generator f_j. Optimizing out the minor coefficients turns each atom into an equivalent pointwise inequality: SOC images Phi_f26 + 2|m_i m_j| >= 0 for {i,j} in {2,12,15}; old-H Phi_f26 + m15^2 + 2|m23 m15| >= 0; chain Phi_f26 + 2|m15|(|m21| + |m23|) >= 0. With 100 BFGS restarts plus Nelder-Mead polishing each (seed 1), the minima were -2.6e-28, -3.9e-28, -1.9e-28, -1.7e-28 and -1.7e-31. The atoms look sharp (minimum 0) and were not violated. Phi_f26 alone reaches -0.0674, so the atoms are not vacuous.
- Realizable points. I sampled points on the T' stratum with lib `random_point_on_contact_stratum` (seeds 0-29, margins 0 and 1e-2). At each point I checked that the representation (1) reproduces the 24 gaps, that C(Re x) = 0 and C + zQ = 0 to rounding, and that G9 vanishes at the exact-contact points (relative size about 1e-14 or less). The sampler found 60 points. G9 vanished at the exact-contact points to a relative size of at most 1.2e-13. All 8 points with d > 0 satisfied one of (3)-(5). The worst minimum of Phi over SO(4) was -7.2e-10. The values of order -1e-10 come from degenerate margin-0 points with 14 to 20 contacts, where the gaps are accurate only to about 1e-10, so they are floating-point noise and not violations. No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Every step is EXACT or PROOF-READ, except the three reference-atom inequalities. Those are CITED. The document neither states them precisely nor names where they are proved, but their role-preserving group images are EXACT and the sharp forms of the images hold numerically. The algebraic core is replayed exactly and independently: the lifted cubic, the Hessian-determinant square identity with G9 re-derived, the 17 positivity tables with the stated term counts, and the chart coverage. It needs no genericity, so it covers repeated eigenvalues, a_i + b_j = 0 and extra contacts.

Gap (presentation): p. 1-2 cite the "sharp SOC", "old-H" and "ordinary chain" reference inequalities without a document reference and without stating their exact hypotheses (z^2 <= xy; x(y - z) >= z^2; chain: center coefficient >= z^2/x1 + z^2/x2). A reader can only infer the forms from (3)-(5). Repair: name the source documents and state each atom. I believe the atoms are true.
