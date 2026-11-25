Evidence level: EXACT (eq. (1), eq. (2), rank H_κ = 4, the single-entry profile and Φ formula, injectivity of P, from the paper's cubic (A.2) through moc4lib); PROOF-READ (Lemma 1, Lemma 2, Theorem 3, Theorem 4, line by line); CITED (the Birkhoff–von Neumann theorem; the external theorems invoked only in the "Lifting and scope" commentary); NUMERIC (one sanity test).

# Audit: final_referee_moc4_six_contact_reduction_20261022

PDF: `proofs/manin_transfer_20261021/final_referee_moc4_six_contact_reduction_20261022.pdf` (4 pages). Batch F4-common-lemmas-b.
Replay: `replay/final_referee_moc4_six_contact_reduction_20261022__check.py`, with output in `..._check_output.txt`: 14/14 EXACT checks pass and 1 NUMERIC check passes.

## Deviations and limits

- The analytic and algebraic-geometric steps were read, not machine-checked. Standard facts used without re-proof: affine dimension theorem; through a point of an irreducible variety on which a regular function is nonconstant there passes an irreducible curve on which it is nonconstant; normalization of an irreducible curve is a connected compact Riemann surface minus finitely many points; removable singularities for bounded subharmonic functions.
- The "Lifting and scope" paragraph (p. 3) is commentary that the theorems do not use. It was read and assessed, but the external theorems it names (a "one-simple three-value theorem", "projection DPPs") were not audited. See E.
- This document is essentially the paper's Theorem 3.6 (with Lemma 3.4 and Lemma 3.5). This project audited that theorem independently in an earlier commit, which I did not consult. My verdict here rests only on this document and on the replay.

## A. Statements

Setting: the paper's V_R, P, X (dim 6, irreducible), κ, reference τ, g_σ = F_σ − F_τ, c(x) = Re(x − F_τ(x)κ), Φ_c(O) = Σ_σ sgn σ (Pc)_σ Π_i O_{i,σ(i)}. K ⊂ V_R is any closed convex cone that contains every single-entry cost c = Σ a_ij n_ij with Pc ≥ 0, for example K = {c : Φ_c ≥ 0 on SO(4)}.

**Eq. (1)–(2) (p. 1).** ∇C(x + tκ) = ∇C(x) + t H_κ x, and ker H_κ = span_C{n_{i,j}} (dim 10, rank H_κ = 4).

**Lemma 1 (p. 2).** If x + Cκ ⊂ X and Pc(x) ≥ 0, then c(x) is a single-entry cost and c(x) ∈ K; Φ_{c(x)}(O) = Σ a_ij O_ij^2 ≥ 0.

**Lemma 2 (p. 2).** C is an irreducible affine curve and g_1..g_N are regular on C, not all constant, with Re Σ_j g_j ≡ 1. c is a real-linear function of (Re g_j) and K is closed and convex. Suppose c ∈ K wherever all Re g_j ≥ 0 and some Re g_j = 0. Then c(p) ∈ K at every p where all Re g_j > 0.

**Theorem 3 (p. 2).** If c(x) ∈ K for every x ∈ X with Pc(x) ≥ 0 and at least six zero entries, then c(x) ∈ K for every x ∈ X with Pc(x) ≥ 0. (These are closed faces: "at least six" contacts, any x ∈ X, which includes all spectra since the statement is on X.)

**Theorem 4 (p. 3).** Fix O_0 ∈ SO(4). If some x ∈ X has Pc(x) ≥ 0 and Φ_{c(x)}(O_0) < 0, then some x' ∈ X has Pc(x') ≥ 0, at least six zero entries, the old contacts kept, and Φ_{c(x')}(O_0) < 0.

There are no contact-specific statements, so no T_k is involved.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| X ⊂ Sing C, ∇C(κ) = 0, P injective, P κ = 1, dim X = 6, X irreducible | paper / moc4lib | EXACT (moc4lib + this script); irreducibility PROOF-READ (closure of the image of an irreducible variety) |
| Birkhoff–von Neumann theorem | external classical | CITED |
| standard algebraic geometry and potential theory (see Deviations) | standard | PROOF-READ |
| Lifting paragraph: semistable representatives (HMSV ordered-points), a "one-simple three-value theorem", "projection DPPs" | external | not needed by Theorems 3–4; UNCHECKED |
| Appendix B, Corollary B.3, single-collision compensation | — | not invoked |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | eq. (1) | identity | EXACT | Symbolic in all 14 coordinates and t. Behind it: ∇C(x+tκ) = ∇C(x) + t H(x)κ + t^2 ∇C(κ), with H(x)κ = H(κ)x (checked) and ∇C(κ) = 0 (checked). |
| 2 | eq. (2) | identity | EXACT | Hess C(κ) computed from (A.2) equals (3.10). Its rank is 4, H_κ κ = 0, and the 16 n_{i,j} span a 10-dimensional subspace inside the kernel, which therefore equals the kernel. |
| 3 | Lemma 1 | identity + inequality | EXACT + PROOF-READ + CITED (Birkhoff) | A line in Sing C forces H_κ x = 0. H_κ is real and kills κ, so c(x) = Re x − Re F_τ κ lies in the real kernel, the real span of the n_ij. (Pc)_σ = Σ_i a_{i,σ(i)} and Φ = Σ a_ij O_ij^2 were checked EXACT on 20 random cases each. (O_ij^2) is doubly stochastic (EXACT on 20 exact O). |
| 4 | Lemma 2: on D the real gaps lie in the simplex, so u = β − L(c) is bounded, positive and harmonic | analytic | PROOF-READ | Re Σ g_j = 1 with all Re g_j > 0 |
| 5 | Lemma 2: boundary values | analytic | PROOF-READ | A boundary point with all gaps > 0 and L < β would lie in D. Otherwise a gap vanishes, so c ∈ K, L ≥ β, and with continuity L = β. |
| 6 | Lemma 2: u extended by 0 is subharmonic, bounded, removable at the points at infinity, hence constant | analytic | PROOF-READ | Same mechanism as the bounded-bad-harmonic lemma; boundedness comes free from the normalization |
| 7 | Lemma 2: then all g_j are bounded in real part near every puncture, hence holomorphic and constant | analytic | PROOF-READ | A pole has a real part unbounded in both directions |
| 8 | Theorem 3: x with r ≤ 5 zeros; freeze the r − 1 other contact gaps and the complex sum; each component of the fibre has dim ≥ 6 − r ≥ 1 | dimension | PROOF-READ | r affine-linear equations on the irreducible 6-dimensional X |
| 9 | Theorem 3: if all gaps are constant on the component, it is the line x + Cκ ⊂ X | algebraic | PROOF-READ + EXACT (P injective, Pκ = 1) | Lemma 1 then contradicts c(x) ∉ K |
| 10 | Theorem 3: the curve keeps the contacts and the normalized sum; boundary points have ≥ r + 1 zeros and fall under the induction hypothesis; Lemma 2 applies to the noncontact gaps | analytic | PROOF-READ | The frozen contacts have fixed complex values with real part 0, so the sum of the noncontact gaps is a constant with real part 1 |
| 11 | Theorem 4: the same, with K' = {Φ_c(O_0) ≥ 0}, a closed half-space | analytic | PROOF-READ | The constant-gap case gives Φ ≥ 0, a contradiction. At most 5 iterations. |
| 12 | Lifting paragraph | commentary | not needed | see E |

## D. Falsification

The document contains no inequality that could be tested independently of the 90 cases. The only checkable algebra is EXACT (C.1–C.3, C.9). NUMERIC sanity: Σ a_ij O_ij^2 ≥ 0 whenever all permutation sums are ≥ 0, on 200 random cases (seed 7); the minimum was 0.56.

## E. Verdict

**confirmed-modulo-citations.** The one citation is the Birkhoff–von Neumann theorem in Lemma 1; the rest is standard algebraic geometry and potential theory, PROOF-READ. Theorems 3 and 4 are correct as stated, for closed contact faces ("at least six" contacts), on all of X.

Presentation remarks (not gaps):
1. The "Lifting and scope" paragraph is unnecessary and partly muddled. Every nonzero x ∈ X has the form α·m(v_1..v_8) for a semistable configuration: the degree-one invariants generate the invariant ring (Kempe; paper Sect. 3) and the quotient map is surjective. A Möbius change then puts all eight points in the affine line whatever their coincidences, and rescaling the v_l changes all matchings by a common factor. Repeated eigenvalues and a_i + b_j = 0 are therefore representable directly. Detouring through a "one-simple three-value theorem" and "projection DPPs" is not needed. The sentence that strictly semistable points "satisfy the support inequality for the same reason" is not justified as written, but nothing depends on it.
2. As the document says, Theorem 3 needs the six-contact cases on CLOSED faces (contact set ⊇ T). Each case proof must state its contacts as non-strict.
