# Audit: root_moc4_T17_uniform_20261022 (paper case T_18)

Evidence level: EXACT for every identity, the generator list, the minor-cone dual and the sign-cell analysis; PROOF-READ for the chart derivation, the H = 0 argument and the region inequalities; CITED for Appendix B (collision layer), the paper's two-orthogonal reduction and (4.2); the real-v layer rests on an argument supplied by the auditor (PROOF-READ plus an EXACT generator check). NUMERIC only for the falsification run.

Document: `proofs/manin_transfer_20261021/root_moc4_T17_uniform_20261022.pdf`, 3 pages, dated 22 October 2026. Supplement label T17, which is the paper's T_18 = {0, 1, 3, 6, 7, 8} (identical as a set, not just up to the group G).

## Deviations and limits

- The replay script was written from the PDF and from the definitions. It uses the shared library `lib/moc4lib.py` for P, the 82 facet profiles and the generator lists. Nothing from the supplement was executed or read as data, including the companion checker the document names.
- Following the document, gaps are computed as g_j = Re{kappa (F_j - F_0)}, with F_sigma = prod_i (A_i - C_sigma(i)) and C = -b. A point at infinity contributes the same constant factor to every product, and that factor is absorbed into kappa. This matches the paper's gap vector for pi* = pi_0, with alpha = kappa.
- Imaginary parts are computed symbolically: v and its conjugate vb are independent symbols, and Im f = (f(v) - f(vb))/(2i) for f with real coefficients. Every displayed identity is checked as an identity of rational functions in (x, h, k, l, v, vb).
- The "three-A/three-C common-circle theorem" (or "six-circle theorem") that the document cites for real v is NOT the paper's Theorem C.1, which needs a whole spectrum on one generalized circle. The supplement does not state this theorem as a standalone result. The auditor supplies the needed argument in full (step S9 below). It is elementary.

## A. Statements

**Theorem 1 (p. 1).** Let x be in X with P c(x) >= 0 and (P c(x))_pi = 0 for every pi in T_18 = {0, 1, 3, 6, 7, 8}. Extra contacts are allowed, and the spectra may be arbitrary finite complex ones, with repeated eigenvalues and a_i + b_j = 0 allowed. Then Phi_{c(x)}(O) >= 0 for every O in SO(4), equivalently every O in O(4) and, by the two-orthogonal reduction, every unitary U. The stronger claim is that, in the layer of eight distinct projective labels with v nonreal, c(x) lies in the cone of the ten minor generators. This covers (*) for T_18 in full, with non-strict contacts and all spectra.

Contact set as permutation indices: [0, 1, 3, 6, 7, 8]. `find_g` matches it with paper T_18 (EXACT).

Intermediate claims, each strong enough to fail:
- (1) The ten minor generators {0,1,2,3,4,5,6,15,16,27} span the same 8-dimensional space as all twelve generators. Their cone has exactly 12 facets: g_j >= 0 for j in {2,5,9,10,13,14,16,19,20,22}, A = g18 - g2 >= 0, and B = g2 - g4 + g10 - g12 + g13 + g18 >= 0.
- Chart (2): A = (x, 1, 1+k(v-1), hv/(v+h-1)), C = (inf, 0, v, w), with the stated w. Contact ratios (R0,R1,R3,R6,R7,R8) = (1, l, hl, x, lx, -x/(k-1)).
- (3) Six sign-kernel identities.
- (4) The supported nonreal layer lies in region I (0<h,k,l<1, x>1) or region II (h,l>1, k<0, 0<x<1).
- (5) Two "Poisson" identities for A/g2 and B/g2, with every summand positive in I and II.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Lemma 4.2 / Prop. 4.1 (82 facets, F_T) | core (paper) | generator list | The library's exact replay is used. Only that g lies in im P with g_T = 0 is needed, and the check that this space has dimension 8 and equals the span of the minor generators is EXACT |
| (4.2): Phi of a minor vector is a squared minor | core (paper) | minor cone implies Phi >= 0 | EXACT in the library and PROOF-READ (library README) |
| Appendix B: Theorem B.1 and Lemma B.2 | common lemma (paper, audited by the project) | collision layer: a_i = c_j, or a repeat within one spectrum | CITED. Common convex weights give Phi = sum mu_sigma g_sigma >= 0 |
| "three-A/three-C common-circle theorem" | common lemma (supplement, not stated) | real-v layer | Replaced by an auditor-supplied argument (S9). Its combinatorial input, that the entry face {sigma(4)=4} together with T has only minor generators [0,1,2,5,6,15,16], is EXACT |
| two-orthogonal reduction (Thm 1.2 / Cor. 2.2) | core (paper) | from O(4) to U(4) | CITED |
| projective (Mobius) covariance | elementary | normalization A0=x, A1=1, C0=inf, C1=0 | PROOF-READ. Each label occurs exactly once in each F_sigma, so a Mobius map multiplies all 24 products by a common factor, which kappa absorbs |

No other contact case and no seven-element set is used.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | F_T = {0,1,2,3,4,5,6,15,16,17,27,28}; nonminor {17,28} | identity | EXACT | `generators(T)` |
| S2 | the minor cone is 8-dimensional; the span of all generators equals the span of the minor ones, which equals im P intersected with {g_T = 0} | identity | EXACT | ranks 8, 8; 14 - rank P[T] = 8 |
| S3 | the dual of the minor cone is exactly the 10 coordinates plus A and B (1) | finite assertion | EXACT | own facet enumeration over 7-subsets in basis coordinates, with a ray-by-ray match up to positive scaling. INFO: A = 1, B = -1 on f17, and A = -1, B = 1 on f28 |
| S4 | kappa = 0 is trivial; collisions go to Appendix B | citation | CITED | Theorem B.1, Lemma B.2 |
| S5 | F0 F7 = F1 F6, hence H = 0, kappa F0 = i mu with mu real and nonzero, and all six contact ratios are real | identity + algebra | PROOF-READ | both sides expand to the same 8 factors. If H != 0, then q0+q7 = q1+q6 and q0 q7 = q1 q6, so -(q0-q1)(q0-q6) = 0, which contradicts F0 != F1, F6 |
| S6 | chart (2): F6/F0 = x, F8/F6 = 1/(1-k), F3/F1 = h, l-equation linear in C3 with the stated unique solution | identity | EXACT | script lines "F6/F0", "F8/F6", "F3/F1", "F1/F0 = l is linear" |
| S7 | the chart is exhaustive on the eight-distinct layer, and the excluded parameter values are collisions | analytic | PROOF-READ | each normalization step is an invertible Mobius or cross-ratio solve. A zero denominator would put A2, A3 or C3 at C0 = inf |
| S8 | contact ratios (1, l, hl, x, lx, -x/(k-1)) | identity | EXACT | |
| S9 | real v: A0,A1,A2 and C0,C1,C2 are real, so the face {sigma(4)=4} = {0,2,6,8,12,14} joins the contact set and the cone becomes minor-only | analytic + finite | PROOF-READ (argument supplied by auditor) + EXACT | Each F_sigma with sigma(4)=4 is (A3-C3) r_sigma with r_sigma real. Contacts 0 and 6 lie in the face with r_0 != r_6, so Re(kappa(A3-C3)) = 0, and every gap on the face is 0. The generators of the face united with T are [0,1,2,5,6,15,16], all minor |
| S10 | six sign kernels (3) | identity | EXACT | 6 script lines |
| S11 | the kernels are nonzero in the nonreal layer (D has only the real zero) | algebra | EXACT / PROOF-READ | |
| S12 | the sign constraints leave exactly regions I and II (4) | case split | EXACT | Every constraint is a product of x, x-1, h, h-1, k, k-1, l, l-1, so feasibility depends only on the 81 cells. The enumeration returns exactly I and II |
| S13 | in I: e, L, P > 0; in II: e, L, P < 0 | inequality | PROOF-READ | I: e = (1-h)+hk, L = (1-hl)+hlk, P = (x-1)+(1-h)+hk. II: h(k-1) < -1 and hl(k-1) < -1 |
| S14 | the two identities (5) | identity | EXACT | |
| S15 | each summand of (5) is positive in I and II | inequality | PROOF-READ | sign bookkeeping in the report text |
| S16 | A, B > 0 puts g in the minor cone, so Phi >= 0 | conclusion | EXACT (S2, S3) + (4.2) | |

## D. Falsification

The script draws 4000 random chart points (seed 20261125). Half are unbiased, with x, h, k, l ~ N(0, 4) and v complex Gaussian. The other half are biased into regions I and II. Each point is tried with both signs of mu. A point counts as supported when every gap is >= -1e-9 max|g|. 2040 points were supported. None lay outside regions I and II, and min(A, B)/max|g| = 3.7e-5 > 0. No violation was found. The main statement is a special case of MOC, and minor-cone membership already implies it, so no SO(4) minimization was run.

## E. Verdict

**confirmed-modulo-citations.** Every identity and finite assertion in the document was replayed exactly and passes. The case analysis and the positivity argument are complete. Citations: Appendix B (Theorem B.1, Lemma B.2) for collisions; the paper's two-orthogonal reduction; (4.2).

Presentation gap (1): p. 2, "the separate six-circle theorem applies". The cited "three-A/three-C common-circle theorem" is neither stated in the document nor a result of the paper. Theorem C.1 needs a whole spectrum on one circle. The claim is true and the auditor supplies the argument (S9). Severity: presentation.

Corollary B.3 / single-collision compensation is not invoked.
