# Audit: blind_moc4_bad70_full_safe_disjunction_20261022 (primary proof of paper case T_34)

> **Second pass:** one step (the endpoint inequalities behind sufficient condition (H)) is not proved in the document; a proof is supplied by the audit; the statement holds. See [the second-pass note](../docs2/blind_moc4_bad70_full_safe_disjunction_20261022__gapcheck.md).

Evidence level: EXACT for the cone, the restricted cubic, the Q ≠ 0 argument, the (PB) and (CQ) Gram identities and Schur complements, and all twelve Bernstein coefficient-sign expansions of Lemma 2 (own code, `replay/blind_moc4_bad70_full_safe_disjunction_20261022__check.py`, all pass). The Hessian identity defining R9 was checked in exact arithmetic at 12 random rational points only (a randomised identity test, not a symbolic proof). GAP: the sufficiency of condition (H) is asserted without proof. NUMERIC: (H) on 25 samples and the main inequality on 5 stratum points.

PDF: `proofs/manin_transfer_20261021/blind_moc4_bad70_full_safe_disjunction_20261022.pdf`, 3 pages. Auditor label C13-T32-T34.

## Deviations and limits

- Three supplement JSON files were loaded as **data** (json.load only): R9 (`data/blind_moc4_next_singlebad_R9_restriction_20261022.json`), the cubic Gram feature polynomials with their four signed variable permutations (`data/blind_moc4_next_singlebad_cubic_reduction_20261022.json`), and the 8×8 change of basis T (`data/blind_moc4_next_singlebad_budget_PSD8_20261022.json`). Every identity involving them was checked by my own code against Leibniz sums I built from the library's generator profiles. The document itself gives neither the Gram vectors nor R9.
- The Hessian identity Q^3 det Hess C(v - (C/Q)κ) = -4Q R9^2 was not expanded symbolically. Each side was evaluated exactly at 12 random rational points.
- Term counts: my Bernstein numerators agree with the document's table for region I (24805 and 37747). They differ for II and III (mine: 25370/25342 and 6166/6166) because I cleared denominators differently. Only coefficient signs matter, and all are nonnegative.

## A. Statements

**Theorem 1 (p. 1).** T = {2, 3, 5, 9, 10, 12} = paper T_34 literally (EXACT). For every x in X with P c(x) in the cone N_T, i.e. g = p f0 + r f1 + a f2 + b f12 + c f39 + q f67 + d f70 + s f72 with all coefficients >= 0 (generator list = `generators(T)`, EXACT; only f70 = h_{4231,+} is non-minor), Phi_{c(x)}(O) >= 0 for all O in SO(4). The proof uses only the fact that Re x lies on {C = 0, det Hess C = 0}, which holds for every x in X. So the statement covers all spectra (repeated eigenvalues, a_i + b_j = 0) and all extra contacts. It is at least as strong as (*) for T_34.

**Lemma 2 (p. 2).** For nonnegative p, r, a, b, c, q, s and d > 0: R9 = 0 implies (H) a >= d + d^2/(b+c), or (PB) a >= d^2/p + d^2/b, or (CQ) a >= d^2/c + d^2/q.

## B. Dependencies

No other contact case and no boundary input (consistent with stage 1 of Table 6). Inputs:

- X ⊂ Sing C: paper (HMSV / Appendix A). The library checks exactly that ∇C vanishes identically on matching vectors.
- The paper's argument on pp. 26–27 that Re x satisfies C = 0 and det Hess C = 0: re-derived and PROOF-READ. I also checked exactly that C(Re x) = 0 and det Hess C(Re x) = 0 at a Gaussian-rational matching point.
- The quaternion parametrisation of SO(4) (standard fact that x ↦ p x q is onto SO(4)): PROOF-READ. A^T A = |p|^2|q|^2 I and det A = (|p|^2|q|^2)^2 are EXACT.
- The two endpoint inequalities behind (H): not stated or proved in the document. This is the GAP below.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = T_34; generators (f0,f1,f2,f12,f39,f67,f70,f72); only f70 non-minor | finite | EXACT | |
| 2 | (PB): Phi_O(p f0 + a f2 + b f12 + d f70) is a Gram form with support {0,2,3,6}, diagonal (p,a,b,d), off-diagonal -d at (0,2) and (2,3) | identity | EXACT | 2·Phi(f_j) = Σ_blocks y_j^2 for j = f0, f2, f12, and 2·Phi(f70) = Σ(y6^2 - 2y0y2 - 2y2y3), modulo the unit-sphere ideal (normal form), with Leibniz sums from moc4lib |
| 3 | (CQ): support {2,4,5,6}, entries +d at (2,4), -d at (2,5) | identity | EXACT | 2·Phi(f70) = Σ(y6^2 + 2y2y4 - 2y2y5); diagonals for f39, f67 |
| 4 | The PSD thresholds are the stated Schur complements | identity | EXACT | determinants of the 3×3 path matrices |
| 5 | **(H) is sufficient**: split d in ratio b:c between two "endpoint inequalities", then Cauchy | inequality | **GAP** | The splitting identity d_b + d_b^2/b + d_c + d_c^2/c = d + d^2/(b+c) is correct (PROOF-READ). The endpoint inequalities themselves are never stated or proved: Phi_O(a' f2 + b f12 + d' f70) >= 0 for a' >= d' + d'^2/b, and the same with f39, c. They are not consequences of the (PB)/(CQ) Gram forms, since those need p > 0 or q > 0. NUMERIC: 25 random (b, c), a = 1 + 1/(b+c), d = 1, other coefficients 0: worst min Phi over SO(4) = -4.6e-13, i.e. zero to rounding. |
| 6 | (2), (3): C(v + Hκ) = C + HQ with the displayed Q, C | identity | EXACT | v = Σ coef · f_j, library cubic (A.2) |
| 7 | Re x = c(x) + Re F_{π*} κ satisfies C = 0 and det Hess = 0 | analytic | PROOF-READ | ∇C(u+iv) = 0 gives ∇C(u) = ∇C(v) and H(u)v = 0; then Euler and polarisation |
| 8 | Q^3 det Hess C(v - (C/Q)κ) = -4Q R9^2 | identity | NUMERIC-exact (12 random rational points; randomised test) | own Hessian of (A.2), exact Fraction determinant |
| 9 | d > 0 ⇒ Q ≠ 0, so R9 = 0 | inequality | PROOF-READ + EXACT | the identity MC - bcQ = bcdU - M[…] is EXACT; the bound uses M^2 - bc = b^2 + bc + c^2 > 0 and the strictly positive last bracket |
| 10 | R9 is homogeneous of degree 9 and quadratic in a; k0 > 0 coefficientwise with constant term 1; for b = c = 0, k2 = 0 and k1 > 0 | identity | EXACT | |
| 11 | Threshold logic: R9 > 0 on [0, min(A0, A1)) gives the disjunction | argument | PROOF-READ | |
| 12 | Regions I, II, III cover all (P, z, M) with M > 0, 0 <= z <= M, P >= 0 | case split | PROOF-READ | I/II: z ∈ (M/(M+1), M) via θ ∈ (0,1), with A1 <= A0 ⇔ P >= P0. III: z ∈ [0, M/(M+1)], A0 <= A1. θ = 1, p = q, Pz = 0 are limits. |
| 13 | 12 Bernstein numerators (3 regions × 2 orders × middle/endpoint) have only nonnegative coefficients after clearing positive denominators | identity | EXACT | flint substitution, 0 negative coefficients in all 12 |
| 14 | Boundary limits: the Bernstein coefficients are continuous polynomials in the original coefficients, and k0 >= 1 | argument | PROOF-READ | |
| 15 | Conclusion of Theorem 1 | assembly | PROOF-READ, except that it inherits the GAP of step 5 | |

## D. Falsification

- NUMERIC: `random_point_on_contact_stratum(T_34, margin 0.01)`, seeds 300–307. Five points found; worst min Phi / sum(gaps) = 1.9e-14 (no violation).
- NUMERIC: (H) test as in step 5 (seed 7, 25 samples, 12 BFGS restarts each). No violation.

## E. Verdict

**gap-found** (one substantive gap; the claim is believed true).

- Gap (p. 1, Section 1, the paragraph on (H)): the two endpoint inequalities that (H) is split into are neither stated nor proved, and they do not follow from the (PB)/(CQ) Gram identities in the document. Belief: true (numerically tight at zero on 25 samples). Repair: state the endpoint inequalities and prove them, for example by an explicit SOS or Gram identity for Phi(a' f2 + b f12 + d' f70) with threshold a' = d' + d'^2/b. (The supplement's check program refers to a "scalar compensation" argument in spin coordinates. I neither ran nor audited it.)
- Everything else is EXACT or PROOF-READ. The one exception is step 8, which is only a randomised exact-arithmetic test; a symbolic expansion would close it.
