# Audit: root_moc4_T16_uniform_20261022 (paper case T_17)

Evidence level: EXACT for the generator expansion, the extra-contact reductions (group elements found), the chart identities, the curve, the z = w divisor, the fibre at v = 1, the branch series, the residues and the compensators. PROOF-READ for the analytic steps: the curve, the boundary behaviour, the harmonic minimum principle and the perturbation at a1 = 0. CITED for the six earlier cases used at extra contacts, Appendix B, (4.2) and the two-orthogonal reduction. The real-v layer rests on an argument supplied by the auditor. NUMERIC for the falsification run.

Document: `proofs/manin_transfer_20261021/root_moc4_T16_uniform_20261022.pdf` (4 pages). Supplement label T16 = paper T_17 = {0, 1, 4, 5, 9, 12}, the same set literally.

## Deviations and limits

- Own code only; the companion checker named on p. 4 was not run or read.
- The document's normalization puts C2 at infinity. A homogeneous bracket with the point at infinity is -1. That sign is common to every product and is kept in the replay, so R5 = 1 - t and F0/(v(z-w)(1-t)) = ((h-1)v - h)/h come out exactly as printed.
- The fibre over v = 1 was replayed with exact truncated Laurent series in s = v - 1, with coefficients in the rational function field QQ(x, h, r, u). The branches come from exact Newton iteration on curve (3), accurate to O(s^6). An earlier run with sympy's `series` reached the same w = 0 residues but was too slow, so the final script uses the field-coefficient series.
- The "three-A/three-C common-circle theorem" cited for real v is not a paper theorem (Theorem C.1 needs a whole spectrum on one circle). The auditor supplies the argument (S9).
- Corollary B.3 / single-collision compensation is not invoked. The v = 1 compensation is internal and done by explicit residues, so the three-point checklist for B.3 does not apply.

## A. Statements

**Theorem 1 (p. 1).** Let x be in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T_17 = [0, 1, 4, 5, 9, 12]. Extra contacts are allowed, and so are all finite spectra, including repeats and a_i + b_j = 0. Then Phi_{c(x)}(O) >= 0 for every O in SO(4), hence for every O in O(4) and every unitary U. This is (*) for T_17, in full.

Intermediate claims:
- The source is simplicial: the generators are {0,1,2,5,7,9,12,35}, f35 = h_{4312,-} is the only nonminor one, and g has the displayed expansion with coefficients that are linear forms in g. Hence (1): L_O(g) >= -C_O d with d = g2.
- Every extra contact with f35 > 0 removes f35. Each of the seven nontied zeros {6,7,11,14,15,17,18} of f35 reduces to an earlier closed case.
- Chart A = (z, 1, x, 0), C = (w, v, inf, tv), with x, h, r, u real and tv = (1-h)v + h. The contacts are H1 = x/h, H4 = x/h - x + 1, H5 = 1, H9 = r, H12 = u.
- (2): H2 and d = mu(x-1) Im v/(h|v-1|^2).
- E9, E12, the factorization of a, the formula for D, curve (3), and no zero of z - w over nonreal v.
- Fibre v = 1 with a1 = 1 - hu != 0: the curve there is (x-1)^2 w(w-1). On the w = 0 branch, res H2 = (x-1)/h = -res H7. On the w = 1 branch, w = 1 + a1/(x-1) (v-1)^2 + ..., z = 1 + (v-1) + ..., and res H2 = -res H14.
- Harmonic descent, and the perturbation argument for a1 = 0.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Prop 4.1 / Lemma 4.2, (4.2) | core | generators; minor terms are >= 0 | EXACT (library) |
| T_4 = {0,1,4,5,7,8} | paper case | extra contacts j = 6, 17. A G-image of T_4 equals {0,1,4,5,6,9}, resp. {0,1,4,5,9,17}, inside T u {j} | CITED; embedding EXACT (find_g) |
| T_5 = {1,2,3,4,6,8} | paper case | j = 7, 15. Images {0,1,4,5,7,9} and {0,1,4,5,9,15} | CITED; embedding EXACT |
| T_9 = {0,1,2,3,11,12} | paper case | j = 14, 18. Images {0,1,4,5,9,14} and {0,1,4,5,9,18} | CITED; embedding EXACT |
| T_15 = {1,2,3,4,6,7} | paper case | j = 11. Image {0,1,4,5,9,11} | CITED; embedding EXACT |
| alternatives T_10, T_11, T_16 | paper cases | each j also embeds one of these | not needed |
| Appendix B (Theorem B.1, Lemma B.2) | common lemma | collisions, including boundary points of the curve that have collisions | CITED |
| common-circle theorem (supplement, unstated) | common lemma | real v | replaced by the auditor's argument. The face {sigma(1)=1} united with T has generators [0,1,2,5,7,9,12], all minor (EXACT) |
| two-orthogonal reduction | core | U(4) | CITED |
| minimum principle for harmonic functions on an open subset of a compact Riemann surface | external, standard | descent | PROOF-READ |

The paper's proof order (Table 6) puts T_17 at Stage 4. The cases used (T_4, T_5, T_15 at Stage 1, T_9 at Stage 3) all come earlier. The document's input list also names T_16 (historical T15, paper Stage 5) and T_8 (Stage 4). Neither is needed: every j has a reduction that uses only earlier stages. So the dependency graph is consistent with Table 6.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | F_T = {0,1,2,5,7,9,12,35}, f35 the only nonminor generator, rank 8 = dim(im P intersected with {g_T = 0}) | identity | EXACT | |
| S2 | the displayed expansion: its coefficient forms are the dual basis of the 8 generators | identity | EXACT | each form evaluated on each generator profile gives the identity matrix |
| S3 | the coefficients are >= 0 on N_T, which gives (1) | inequality | PROOF-READ | simplicial cone; (4.2) |
| S4 | an extra contact j with f35(j) > 0 forces the f35 coefficient to 0, so g is minor | algebra | PROOF-READ | f35 > 0 at {2,3,8,10,13,16,19,20,21,22,23} (EXACT) |
| S5 | each j in {6,7,11,14,15,17,18} reduces to an earlier case | citation | EXACT embedding + CITED | table above |
| S6 | a hypothetical counterexample is strict outside T | logic | PROOF-READ | |
| S7 | chart: R1 = x(v-1)/v, R4 = R1 + (t-1)(x-1), R5 = 1-t, so x is real; h = v(1-t)/(v-1) is real and nonzero; h != 1 | identity + logic | EXACT + PROOF-READ | |
| S8 | kappa = i mu/[v(z-w)(1-t)]; H1, H4, H5 as stated; g_j = -mu Im H_j | identity | EXACT | |
| S9 | real v: (1,x,0) and (v,inf,tv) are real, so the face {sigma(1)=1} = {0,...,5} becomes all contacts and the cone is minor-only | analytic + finite | PROOF-READ (argument supplied by auditor) + EXACT | F_sigma = (A0-C0) r_sigma with r_sigma real. Contacts 0 and 1 have r0 != r1, so Re(kappa(A0-C0)) = 0 |
| S10 | (2): H2 and d | identity | EXACT | |
| S11 | H9 = r iff E9 = 0, H12 = u iff E12 = 0 (the numerators are equal); a and D as stated; a E9(z = (Dw + t(v-x))/a) = curve (3) | identity | EXACT | |
| S12 | at nonreal v, a != 0, D(x - tv) != 0, w and z are finite, and v is nonconstant on any component | analytic | PROOF-READ | the zeros of a, D and x - tv are real. If u = 1 - x then D = x - 1 (EXACT) |
| S13 | no zero of z - w at nonreal v | identity + logic | EXACT + PROOF-READ | |
| S14 | F0/(v(z-w)(1-t)) = ((h-1)v - h)/h | identity | EXACT | |
| S15 | approaching any real v != 1, or infinity, inside the supported region gives liminf L_O >= -C_O lim d = 0 | analytic | PROOF-READ | d depends only on v |
| S16 | the fibre at v = 1 is (x-1)^2 w(w-1): two simple roots, so exactly two smooth branches with v as local parameter, and no point over v = 1 at infinity | identity + analytic | EXACT + PROOF-READ | the leading coefficient at v = 1 is (x-1)^2 != 0 |
| S17 | residues: on w = 0, H2 and H7 have simple poles with residues +-(x-1)/h; on w = 1, H2 and H14 likewise; w and z expansions; the compensator sums are regular with real coefficients (values (-xh + x + hr + hu - 1)/h and (hr + hu - 1)/h) | identity | EXACT | Laurent series over QQ(x,h,r,u) |
| S18 | on the supported region 0 < d <= g2 + g7, resp. g2 + g14, -> 0 at the two points, so liminf L_O >= 0 | inequality | PROOF-READ | |
| S19 | harmonic descent: L_O = Im of a holomorphic function on U. Every boundary point of U in the compact normalization is over nonreal v (an extra-contact configuration, handled by S4/S5 or Appendix B) or over real v / infinity (S15, S18). U is not the whole curve, so L_O >= 0 on U | analytic | PROOF-READ | |
| S20 | the stratum a1 = 0 is reduced to a1 != 0 by perturbing (x,h,r,u) with v fixed | analytic | PROOF-READ (details supplied by auditor) | The roots of (3) move continuously, since the leading coefficient is nonzero. z is continuous because a != 0. The six contacts are built into the chart. Strictness, distinctness and L_O < 0 are open conditions |
| S21 | conclusion for O(4), then U(4) | citation | CITED | |

## D. Falsification

3000 chart samples were drawn (seed 20261125): x, h, r, u ~ N(0, 4), v complex Gaussian, both roots w of (3), and both signs of mu. A sample counts as supported if the six contacts hold to 1e-8 relative and all other gaps are >= -1e-9 relative. On those samples the eight expansion coefficients were checked to be >= 0, and at the first 12 supported points Phi was minimized over SO(4) (6 BFGS restarts each). 33 samples were supported with exactly the six contacts. The smallest expansion coefficient over max|g| was 3.8e-4 > 0. The worst minimum of Phi over SO(4) at 12 supported points was 9.5e-14, which is 0 up to rounding (the minimum sits on a contact). No violation was found. The main statement is a special case of MOC, so this test proves little. The exact series checks carry the weight.

## E. Verdict

**confirmed-modulo-citations.** Citations: paper cases T_4, T_5, T_9 and T_15, each through an exactly verified G-embedding into T_17 u {j}; Appendix B (Theorem B.1, Lemma B.2); (4.2); the two-orthogonal reduction.

Gaps, all of presentation severity:
1. p. 1, "the exact checker supplies a previously closed family image inside T16 u {j} for each". The document names no case for any j. The auditor resolved them (table in B), and they are consistent with the paper's Table 6.
2. p. 2, real v: "satisfy the stated support inequality" rests on an unstated three-A/three-C common-circle theorem. The auditor supplies the argument (S9).
3. p. 3, a1 = 0: the perturbation argument is brief but correct. The auditor wrote out the details (S20).
