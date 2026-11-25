Evidence level: EXACT (common signed representation, facet identities, cubic root, span of minor vectors); EXACT via the shared library (82-facet completeness, paper Prop. 4.1, with an independent certificate); PROOF-READ (lower bound ν ≥ ρ, zero case, convexity, distance bound); CITED (upper bound ν ≤ ρ in the positive case and the value of max ρ over O(4): Holden [2], not available; the two-orthogonal reduction, paper Thm 1.2); NUMERIC (ν = ρ by linear programming on 300 samples; max ρ = τ* by multistart).

# Audit: root_complete_negative_mass_20261022

PDF: `proofs/manin_transfer_20261021/root_complete_negative_mass_20261022.pdf` (3 pages). Batch F4-common-lemmas-b.
Replay: `replay/root_complete_negative_mass_20261022__check.py`, with output in `..._check_output.txt`: 7/7 EXACT checks pass and 4 NUMERIC checks pass.

## Deviations and limits

- By the brief, the document's own completeness certificate (6846 simplices, verifier `root_permutation_polytope_verify_20261022.py`, JSON data) was neither executed nor loaded. Completeness of the 82-facet list rests instead on the shared library's exact replay of the paper's Proposition 4.1. That replay uses a different certificate: a mod-2 cycle of 9420 simplices. The document's own orientation argument was read and is sound as an argument; its data are UNCHECKED.
- Reference [2] (Holden, "A fourteen-term identity and a sharp orthogonal obstruction bound") is not in our materials. Its two uses are CITED. Both are supported by strong numerical evidence, but neither is re-proved here.

## A. Statements

Notation. For U ∈ U(4), Q_U(I,J) = |det U[I,J]|^2. A common signed representation is a w ∈ R^24 with Σ_σ w_σ 1{σ(I)=J} = Q_U(I,J) for all I, J (|I| = |J|, including the empty set, so Σ w_σ = 1). ν(U) = min_w Σ_σ max(0, −w_σ). With V = U[:, π]: t_π = Σ_i Q_U({i},{π(i)}) and c_π = Q_U(12,{π1,π3}) + Q_U(12,{π2,π4}) − Q_U(13,{π1,π2}) − Q_U(13,{π3,π4}). These are the paper's t_π and d_π evaluated at the response. Then h_{π,ε} = −(t_π + ε c_π)/4 and ρ(U) = max(0, max_{π,ε} h_{π,ε}).

**Theorem 1 (p. 1).** (i) ν(U) = ρ(U) for every U ∈ U(4). In particular U has a common convex representation of all its minor squares iff all 48 inequalities h_{π,ε} ≤ 0 hold. (ii) max_{U(4)} ν = τ*, the unique positive root of 32τ^3 + 48τ^2 + 28τ − 1 = 0, with 0.033721124434 < τ* < 0.033721124435. The maximum is attained on O(4).

**Corollary (p. 3).** det(D_a + U D_b U*) = (1+ν) z_+ − ν z_− with z_± ∈ K = conv{z_σ}. Hence dist(det(D_a + U D_b U*), K) ≤ τ* diam K. The document says explicitly that this does not prove MOC.

This is a common-lemma document (not a contact case), so no contact sets are involved.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| Q has exactly the 82 facets of Appendix D (paper Prop. 4.1) | core | EXACT (moc4lib's independent certificate). The document's own certificate is UNCHECKED. |
| Jacobi complementary minors and generalized Laplace expansion (existence of a signed representation) | classical | PROOF-READ, and EXACT on instances |
| Holden [2, Thms 2.1, 3.1]: an explicit representation with exactly one negative coefficient, of size ρ, when ρ > 0 | external, unavailable | CITED. NUMERIC support: ν = ρ in 300/300 samples. |
| Holden [2, Thm 7.2]: max over O(4) of ρ is τ* | external, unavailable | CITED. NUMERIC: the multistart maximum over O(4) is 0.033721124434. |
| two-orthogonal reduction (paper Thm 1.2) | core | CITED (not audited in this batch) |
| Kovačec–Ribeiro [1] (facet count reported earlier) | external | not needed (attribution only) |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | W_σ = sgn σ Π u_{i,σ(i)} / det U satisfies Σ_{σ(I)=J} W_σ = Q_U(I,J) (real), so Re W is a signed representation with Σ = 1 | identity | EXACT on 4 Gaussian-rational unitaries, all I, J with \|I\| = 0..4; PROOF-READ in general | Laplace: Σ_{σ(I)=J} sgn σ Π u = ± det U[I,J] det U[I^c,J^c]. Jacobi: det U[I^c,J^c] = ± det U · conj(det U[I,J]). |
| 2 | the response satisfies the 34 minor facets | inequality | PROOF-READ | Q_U ≥ 0 |
| 3 | the facet value of h_{π,ε} at the response is (t_π + ε c_π)/2 | identity | EXACT | All 48 (π, ε), on 4 exact unitaries, against the library's h_vec |
| 4 | the minor-response coordinates are an injective image of V_R^*, so the doc's polytope is the paper's Q | identity | EXACT | The 34 minor vectors have rank 14 |
| 5 | ρ = 0 implies the response lies in Q, so ν = 0 | polytope | EXACT dependency + PROOF-READ | The response lies in the affine hull (step 1) and satisfies all 82 facets; completeness (Prop. 4.1) |
| 6 | ν ≥ h_{π,ε} for each π, ε | inequality | PROOF-READ | The facet function f takes values in {0,1,2} at the vertices and −2h at the response, so −2h = Σ w_σ f_σ ≥ −2 Σ max(0, −w_σ) |
| 7 | ν ≤ ρ when ρ > 0 | inequality | CITED [2] | NUMERIC: an LP for ν against the formula for ρ, 300 samples, max \|ν − ρ\| = 3.5e-14; ρ > 0 in 37 of them |
| 8 | ρ is convex in the response; Thm 1.2 then gives max_{U(4)} = max_{O(4)} | convexity | PROOF-READ + CITED (Thm 1.2) | ρ is a maximum of affine functions of Q_U |
| 9 | max_{O(4)} ρ = τ* | value | CITED [2, Thm 7.2] | NUMERIC: 100 Nelder–Mead starts give 0.033721124434 (agreement to 1e-16). EXACT: the cubic has one positive root, inside the stated interval. Consistent with the library's numeric min Φ_{h_id,−} = −0.067442 = −2τ*. |
| 10 | the distance corollary | identity/inequality | PROOF-READ | (2.1) plus the sizes \|I\| = 0..4 of the representation; w = w_+ − w_−, Σ w_+ = 1 + ν, Σ w_− = ν |
| 11 | the document's completeness argument (oriented 12-cycle with nonnegative local degrees) | topology | PROOF-READ (argument); data UNCHECKED | Superseded by the library's exact replay |

## D. Falsification

- Theorem 1(i): linear programs on 150 Haar unitaries and 150 random SO(4) matrices (seed 99) gave ν = ρ to 3.5e-14.
- Theorem 1(ii): maximising ρ over U(4) with 50 starts never exceeded τ*. The best value was 0.033721124434.
- No discrepancy was found.

## E. Verdict

**confirmed-modulo-citations.** Citations: Holden [2] (Theorems 2.1, 3.1 and 7.2; not available to us) for the upper bound ν ≤ ρ and for the sharp orthogonal value; paper Theorem 1.2 for the two-orthogonal reduction; paper Prop. 4.1 (EXACT in the library) for completeness. Every step that does not depend on [2] is EXACT or PROOF-READ. The two [2]-dependent claims agree with our numerics to machine precision. No gap was found.

A remark for users of this document: it gives a common SIGNED representation, with negative mass at most τ*, independent of the spectra. It does not give MOC. The document says so itself.
