Evidence level: EXACT for the generator bookkeeping, the two substitution identities, the zero-face and seventh-tie reductions, the height equation at c = 0 and the positivity of (e+j)C0 on Q0 = 0; EXACT pointwise (20 random rational points, own Hessian, supplement polynomial G read as JSON data) for det H = -4G^2/Q^2; PROOF-READ for Lemma 2 and the assembly logic; CITED for the four sub-theorems that carry the analytic and inequality content.

# Audit: root_true_T48_whole_assembly_20261023 (primary proof of paper T_49)

PDF: proofs/manin_transfer_20261021/root_true_T48_whole_assembly_20261023.pdf, 4 pages. Replay: `replay/root_true_T48_whole_assembly_20261023__check.py` (output in `__check_output.txt`).

## Deviations and limits

- This document is an assembly. Its analytic content (harmonic curve argument, Cramer/orientation theorem, cubic interlacing, all-j chain capacity with an 84 925-coefficient table) lives in four other documents that are not in this batch. I re-checked only what this document itself asserts and whether the cited statements, as restated here, fit together. They are CITED.
- The determinant identity (3) was checked exactly at 20 random rational points, with my own 14x14 Hessian of the paper's cubic (A.2) and the supplement's 1772-term G loaded as data. This proves the identity at those points; as a polynomial identity it holds with overwhelming probability (Schwartz-Zippel) but I did not expand it symbolically.
- The Cramer determinant, the relations mN(Rv+1) = tH, L/y^3, B/y^2, F6/y^6 and the interlacing sign pattern (p.3) are restated from the cited documents in coordinates (J, Z, r, s, v, t, u, h) that this document does not define completely. I did not replay them (UNCHECKED here; CITED).
- No numerical search for violations of the main statement (it is a special case of the conjecture; per instructions).

## A. Statements

**Theorem 1 (p.1).** For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in T = {0,3,5,6,9,12}, Phi_{c(x)}(O) >= 0 for all O in SO(4); the unitary case by Theorem 1.2/Cor. 2.2. Contacts non-strict, all spectra. T is literally paper T_49 (EXACT). Degenerate spectra: repeated eigenvalues and cross collisions a_i + b_j = 0 are delegated to "global collision/repetition safety", i.e. paper Appendix B (Theorem B.1, Lemma B.2), which this audit has established; this matches (MOC holds outright at such spectra). The statement is (*) for T_49.

**Lemma 2 (p.2).** F real homogeneous cubic, x, y in Sing F: F(x+y) = 0 and, if x+y != 0, H(x+y) is singular. This is the established "radial necessity" (item 5 of the audit's established list) in a slightly different form; the proof given (polarization, H(x+y)(x-y) = 0, and H(2x)(2x) = 0 when x = y) is correct (PROOF-READ). Applied with x = alpha m/2, y = conj(alpha m)/2, both in X, which is conjugation-invariant.

## B. Dependencies

| ref | kind | used for | status |
|---|---|---|---|
| independent_true_T48_native_harmonic_conditional_20261023 | branch | reduces the whole six-contact family to the frontier c = 0 (contact 17), full pole atlas, Rouche restoration | CITED |
| independent_true_T48_c0_native_Heron_20261023 | branch | F6 >= 0, sign L = -sign r at strict c = 0 sources (Cramer) | CITED |
| blind_moc4_T48_c0_native_interlacing_20261023 | branch | excludes r > 0 | CITED |
| blind_moc4_T48_c0_all_j_chain_20261023 | branch | G = 0, L >= 0, F6 >= 0 imply pae >= d^2(a+e) for every j > 0 | CITED |
| chain atom (centre f5 = p; leaves f0 = a, f2 = e; non-minor f40) | common lemma (so4-atoms) | pae >= d^2(a+e) is safe | tuple verified in the proved 576-orbit with roles as used (EXACT, in the allocation-III replay) |
| Appendix B (Thm B.1, Lemma B.2) | core | repeated / cross-colliding spectra | established by this audit; use matches |
| Lemma 3.4 + Birkhoff ("singleton kernel") | core | common-height punctures | established (paper Sect. 3) |
| bounded-negative-part / maximum principle lemma | common lemma | harmonic step | CITED |
| T_49 + {17} = Table 5 set 23 | seven-contact | the frontier | proved here (via the cited sub-theorems), paper stage 2 |
| six other seventh ties {2,7,11,13,16,23} and zero faces | earlier cases | boundaries | EXACT: each T_49+{j} contains a six-set of paper stage < 4 (stage of T_49) |

## C. Steps

| # | claim (page) | kind | status | how |
|---|---|---|---|---|
| 1 | generator set {0,1,2,5,7,9,29,39,40,41}, only f40 non-minor (p.1) | identity | EXACT | lib `generators` |
| 2 | f40+f1 = f7+f39+f0, f40+f9 = f7+f2+f41 (p.1) | identity | EXACT | 24 values (P injective) |
| 3 | the substitutions either exhaust d (minor decomposition) or give canonical form (1) with q,p,j,c,e,b,a >= 0, d > 0 | case split | PROOF-READ | replacing min(d, coeff f1) f1 then min(d', coeff f9) f9 by the identities keeps all coefficients >= 0; exhaustive |
| 4 | canonical eight profiles independent | identity | EXACT | rank 8 |
| 5 | Lemma 2, w = g + rho kappa realizable lies in Sing-F secant | analytic | PROOF-READ | matches established radial necessity |
| 6 | F(g + rho kappa) = C0 + rho Q0, Q0 = d(d+e+j) - p(e+j) at c = 0 (p.2) | identity | EXACT | own cubic; F affine in rho since kappa in X |
| 7 | Q0 = 0 forces (e+j)C0 > 0 (27 positive terms, d^3 j present), contradiction | identity | EXACT | own expansion |
| 8 | det H(g - C0 kappa/Q0) = -4G^2/Q0^2, G|c=0 has 759 terms (p.2) | identity | EXACT at 20 points (pointwise), see limits | own Hessian; supplement G as data |
| 9 | c = 0 frontier: d = 0 minor; zero coefficients of q,p,j,e,b,a reduce to earlier six-sets (p.3) | reduction | EXACT | each zero face adds one contact (16, 11, 2, 13, 7, 23) and T_49 plus it contains a six-set of stage <= 3 |
| 10 | Cramer analysis, m != 0, sign L = -sign r, F6 >= 0 (p.3) | citation | CITED (UNCHECKED here) | independent_true_T48_c0_native_Heron |
| 11 | interlacing excludes r > 0; r = 0 excluded by eight-distinctness (p.3) | citation | CITED | blind_moc4_T48_c0_native_interlacing |
| 12 | all-j chain: below pae = d^2(a+e) no realizable point; at or above, chain atom (p.3) | citation + atom | CITED; atom tuple EXACT | the chain tuple (5;0,2;40) is in the proved orbit |
| 13 | harmonic return to the six-contact family, all punctures, Rouche (p.4) | citation | CITED | independent_true_T48_native_harmonic_conditional |
| 14 | stage consistency: everything used precedes T_49 | bookkeeping | EXACT for the six-sets; the four branch documents are T_49-internal | Appendix G stages |

## D. Falsification

Lemma 2 sanity check on 50 random realizable w = Re(alpha m) (NUMERIC, seed 2): max scaled residual 7.9e-16. No search on the main statement. The determinant identity (step 8) is the only strong intermediate claim made in this document and it passed exactly at every tested point.

## E. Verdict

**confirmed-modulo-citations.** Every step that this document itself carries is EXACT or PROOF-READ. Its conclusion rests on the four cited branch theorems, on the bounded-negative-part lemma and on the chain atom (whose tuple is in the proved orbit). No gap found in this document. Presentation note: the document says the frontier is "c = 0, or contact 17"; the seven-contact cone of T_49+{17} also contains f9, which is eliminated by the second substitution identity, so the two descriptions agree (EXACT).
