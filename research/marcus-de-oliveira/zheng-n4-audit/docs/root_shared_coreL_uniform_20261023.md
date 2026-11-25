# Audit: root_shared_coreL_uniform_20261023 (seven-contact set L)

Evidence level: EXACT for every identity, generator decomposition, group relabelling, the Hessian identity (10) and all 12 coefficient tables; PROOF-READ for the capacity allocations, the real-part bridge, the C > 0 argument and the chart coverage; CITED for the two families of coefficient atoms (sharp pair-SOC, old-H) and for the paper's Theorem 1.2. NUMERIC only for falsification.

PDF: `proofs/manin_transfer_20261021/root_shared_coreL_uniform_20261023.pdf`, 4 pages. Replay: `replay/root_shared_coreL_uniform_20261023__check.py` (output in `__check_output.txt`, 34 EXACT checks, all pass).

## Deviations and limits

- None from the audit plan. Nothing from the supplement was executed or loaded, not even its JSON: the degree-8 polynomial G8 that the document defines through its "companion coefficient array" was **re-derived** from the paper's cubic (A.2) by my own fraction-free determinant, and every one of the 12 tables was rebuilt from that.
- The atoms "sharp pair-SOC (f1, f9; f22)" and "old-H (f3, f15; f26)" are pointwise inequalities proved elsewhere in the supplement. I checked only that the group relabellings the document uses exist (EXACT) and ran a brief numerical non-violation test. They remain CITED.
- Time spent: about 30 minutes, shared with the M document.

## A. Statements

Contact set L = {0, 1, 2, 7, 10, 16, 23} (0-based lexicographic indices). This is the paper's Table 5 set #14 (EXACT: `AUX7[14]`). It contains G-copies of the six-sets T_55, T_63, T_64 and T_85 (remove 16 or 23; 0 or 7; 2 or 10; 1). It extends T_85 by one contact for j in {3, 4, 8, 11, 12, 15, 19, 20}, and T_55 or T_63 for j = 20.

**Theorem 1 (p. 1).** Let x be in X and pi* in L, let c = c(x), and suppose P c >= 0 and (P c)_sigma = 0 for every sigma in L. Extra contacts are allowed. Then Phi_c(O) >= 0 for every O in SO(4), and hence on every real orthogonal response in the paper's O(4) normalization. The unitary form follows from the paper's Theorem 1.2. No distinctness of eigenvalues is assumed, so all spectra are covered, including repeated eigenvalues and a_i + b_j = 0. The proof uses only x in X, a subset of Sing C (covers_all_spectra: yes).

Internal claims:
- (1)-(2): the contact kernel {c : (Pc)_s = 0 for s in L} is 7-dimensional, with basis f2, f3, f11, f12, f15, f23, f26. A cost is c = a f2 + b f3 + q f11 + c f12 + e f15 + v f23 + d f26, and (a, b, q, c, e, v) = (g15, g14, g17, g8, g4, g9).
- Relation f8 + f26 = f2 + f3 + f12 + f23, with g13 = a + d, g20 = b + d, g22 = c + d and g3 = v + d.
- (3)-(5): if d > 0 (normalized to d = 1) and Phi is negative somewhere, then H3(a, c, e) < 1 and xi(e - 1) < 1, where xi = b + v.
- (6)-(7): C(w0 + z kappa) = C + z Q for the paper's cubic C.
- (8)-(9): C > 0 on (5).
- (10): Q^4 det H_C(w0 - C kappa / Q) = -4 (a + d)^2 Q^2 G8^2, where G8 is homogeneous of degree 8 with 327 terms and [d^8] G8 = 1.
- (11)-(16): G8 > 0 on all of region (5), by 12 cleared tables with nonnegative coefficients.

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| sharp pair-SOC atom (f1, f9; f22): Phi(x f1 + y f9 + delta f22) >= 0 whenever xy >= delta^2 | common lemma (supplement, "entire-coreD package") | CITED | EXACT: G elements map (f1, f9; f22) to (f2, f12; f26), (f2, f15; f26) and (f12, f15; f26), 8 elements each. Phi is G-invariant (paper (4.7), replayed in lib). NUMERIC: no violation at xy = delta^2 |
| old-H atom (f3, f15; f26): Phi(b f3 + e f15 + delta f26) >= 0 whenever b(e - delta) >= delta^2 | common lemma | CITED | EXACT: G maps it to (f23, f15; f26). NUMERIC: no violation |
| Phi(n_{I,J}) = det O[I,J]^2 >= 0 (paper (4.2)) | paper | EXACT/PROOF-READ in lib | used for the d <= 0 side and for the minor remainders |
| X is contained in Sing C (paper App. A.1) | paper | EXACT in lib | the bridge F(Re x) = 0, det H(Re x) = 0 |
| Theorem 1.2 (two real orthogonal responses) | paper | CITED | used only for the unitary sentence |

No other contact case, branch or Appendix B/C result is used. Corollary B.3 and "single collision compensation" are not invoked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | generators(L) = {2, 3, 8, 11, 12, 15, 23, 26}; contact rows have rank 7; the 7 listed profiles are a basis; coefficients equal the gaps g15, g14, g17, g8, g4, g9 | identity | EXACT | lib `generators`, `rank`; symbolic gap vector |
| 2 | f8 + f26 = f2 + f3 + f12 + f23; g13, g20, g22, g3 = a+d, b+d, c+d, v+d; for d <= 0 the cost is a nonnegative minor combination (f8 = n_{12,23} is minor) | identity | EXACT | symbolic |
| 3 | SOC atoms on the three pairs give ac, ae, ce >= delta^2 as sufficient | citation + relabel | CITED (relabel EXACT) | group search |
| 4 | the allocations give the sufficient capacity H3 in (3) | inequality | PROOF-READ | balanced case: pair amounts (a+c-e)/2, etc. use a, c, e exactly and give total delta = S/2. Dominant case: proportional split gives sqrt(M(S - M)). If H3 >= 1, scale every delta by 1/H3 |
| 5 | old-H split: xi(e - 1) >= 1 is sufficient | inequality | PROOF-READ (atom CITED, relabel EXACT) | b(eb/xi - b/xi) >= (b/xi)^2 is equivalent to xi(e - 1) >= 1; zero leaves are omitted |
| 6 | a negative response at d = 1 forces (5) | logic | PROOF-READ | Phi is linear in c; leftovers are minor |
| 7 | F(w0 + z kappa) = C + z Q with (6), (7) | identity | EXACT | the paper's cubic (A.2) on the lifted cost. The match is exact, with no scalar factor |
| 8 | C = ace + (q+1)(xi+1) + P(q+1-e); identity (9); hence C > 0 on (5) | identity + inequality | EXACT + PROOF-READ | symbolic. RHS(9) >= 0 because xi(e-1) < 1, so C >= 1 + 3xi/4 + ... > 0 when xi > 0, and C = ace + q + 1 when xi = 0 |
| 9 | for realizable c, w = Re x = c + Re F_pi*(x) kappa satisfies C(w) = 0 and det H(w) = 0 | analytic | PROOF-READ (argument checked by auditor) | x̄ is in X; x, x̄ are in Sing C; trilinear polarization gives C(x + x̄) = 0 and H(x + x̄)(x - x̄) = 0; if x is real, Euler gives H(x)x = 0 with x ≠ 0, since d = 1 > 0 |
| 10 | C > 0 forces Q ≠ 0 and z = -C/Q | logic | PROOF-READ | linear in z |
| 11 | (10): Q^4 det H(w0 - C kappa/Q) = -4(a+d)^2 Q^2 G8^2; G8 has 327 terms, max degrees (1,2,2,2,2,2,8), [d^8] = 1 | identity | EXACT | own Bareiss on the 14x14 Hessian along w0 + z kappa (d = 1, 819 exact divisions); z-degree 3; exact division and an exact square root |
| 12 | hence G8 = 0 at any unsafe realizable source, since a + d > 0 | logic | PROOF-READ | |
| 13 | region (5) is covered by: balanced (12); dominant a or c (13) in both split orders, plus the S = 0 tables; dominant e <= 1 (14) with E + T > 0, plus the S = 0 table; dominant e > 1 (15)-(16) using b<->v | case split | PROOF-READ (auditor rederived every inverse and range) | e.g. (13): S = 1/(1+X) in (0,1) and m in [S, 1/S); (14): the excluded point e = S = 1 has H3 = 1, which is safe; (15): t-ranges and (16) cover xi > 0 with v > 0; xi = 0 is B = 0 |
| 14 | G8 is invariant under b <-> v, and (5) is symmetric | identity | EXACT + PROOF-READ | |
| 15 | the 12 cleared tables (11) have nonnegative integer coefficients with the stated markers (constant 1; 4X; 8X/4X; 27-, 72- and 80-term tables with constants 1, 1, 2; 8E + 4T and 4E + 2T; 2E^4) and 29584 coefficients in total | finite certificate | EXACT | every table and every term count reproduced exactly from my own G8 |
| 16 | contradiction, so Phi >= 0 | logic | PROOF-READ | |
| 17 | unitary form | citation | CITED | paper Theorem 1.2 |

## D. Falsification (NUMERIC)

- The three SOC atoms at xy = delta^2 and the two old-H atoms at b(e - delta) = delta^2, with s in {0.05, 0.2, 0.5, 1, 2, 5, 20}, minimized over SO(4) (12 BFGS restarts each, seed 1): minimum >= -3e-16. No violation.
- Realizable points on the stratum with contact set exactly L (lib `random_point_on_contact_stratum`, margin 0.01, seeds 0-2; all three have d > 0): (a + d) G8 / scale^9 is at most 2e-21 in absolute value, consistent with the necessary equation G8 = 0. The capacity test succeeds at all three points. min Phi over SO(4) >= -1e-14 (Phi = 0 is attained at the contact permutation matrices).
- Earlier seeds without a margin gave degenerate points (many extra contacts, d <= 0): all had min Phi >= -3e-10, which is numerical zero.

## E. Verdict

**confirmed-modulo-citations.** Every algebraic and finite step was replayed exactly from the paper's definitions, with term counts identical to the document's. G8 was re-derived independently, not taken from the supplement. The analytic bridge (step 9) is correct and needs only X ⊂ Sing C. The case split covers region (5) together with its boundary walls. The citations are the sharp pair-SOC atom (f1, f9; f22), the old-H atom (f3, f15; f26) and paper Theorem 1.2.

Gaps: none. Presentation remarks:
- The text says the seven profiles "give the entire contact kernel". The generator set of L also contains f8, which lies in their span by the displayed relation, so the claim stands.
- (8) is stated with the reason "4P <= xi^2". The rigorous form is (9), which is exact.
