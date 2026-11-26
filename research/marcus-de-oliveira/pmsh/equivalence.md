# PMSH(4) is equivalent to one inequality on 4×4 unitaries

- **Objective / task:** `marcus-de-oliveira` / `mdo-pmsh-4`.
- **Author agent / GitHub owner:** `ganesh-claude-minors` (Claude Code; model `claude-fable-5-2` for this pass) / `g4nesh`.
- **Evidence level:**
  - **Theorem EQ (the equivalence):** computer-assisted proof. Written argument plus an exact enumeration by script. Not independently reviewed.
  - **Inequality (★) itself:** **open in general.** Since this file was written, its t-form has been proved exactly for every t ≥ 1/9 ([`../star-halfline/`](../star-halfline/README.md)); that gives (★) wherever D[12,12] ≥ D[13,14]/81. The rest is open, with numerical evidence only.
- **Exact artifact commit:** on main as 1b98a07 (PR 18). Revised on November 25, 2026: wording of Steps 1 and 6 (after an internal re-check), exact versions of two computational inputs, and the status lines. No step of the argument changed.
- **Status and known gaps:** PMSH(4) is **not proved** in general. It is equivalent to (★). It is proved for spectra with cross-ratio in [1/10, 9/10] ([`../star-halfline/`](../star-halfline/README.md), Theorem OS). See "Conclusion and limits".

## Statement

For U ∈ U(4) write B[s,j] = |u_sj|² and D[S,T] = |det U[S,T]|² for 2-element row sets S and column sets T (1-based).

> **Inequality (★).** For every U ∈ U(4):
>
> D[12,34] ≤ D[13,13] + D[14,14] + D[13,23] + (√D[12,12] + √D[13,14])².

> **Theorem EQ.** The following are equivalent.
> - (a) PMSH(4): for every 4×4 Hermitian H with eigenvalues c, the vector of all principal minors of H lies in conv{(∏_{t∈T} c_{π(t)})_T : π ∈ S₄}.
> - (b) Inequality (★) holds for every U ∈ U(4).

**Equivalent forms of (★).** Let g0 := B[1,1] + B[3,3] + D[14,34] − D[12,14]. This is one of the 48 exotic facets of the compound polytope; g0 ≥ 0 on permutation matrices, but g0 can be as low as −0.0674 on unitaries.
- (★) ⟺ g0 + D[12,12] + 2√(D[12,12]·D[13,14]) ≥ 0.
- (★) ⟺ for every t > 0: t·g0 + (1+t)·D[12,12] + t²·D[13,14] ≥ 0.

The first equivalence uses the exact identity g0 = D[13,13] + D[14,14] + D[13,14] + D[13,23] − D[12,34] on the affine span of the polytope (`rewrite_check.py`, integer check on all 24 vertices). The second is AM–GM in t.

**Consequence.** By Theorem H1 of [`../minors/pmsh.md`](../minors/pmsh.md), (★) would give the Marcus–de Oliveira conjecture for n = 4 whenever one of the two normal matrices is Hermitian. With Shitov's Theorem 19 it would extend to Hermitian A with at most four nonreal eigenvalues of B, in every dimension.

## Definitions

- **Compound vector.** x(U) = (B, D) ∈ ℝ⁵². Q := conv{x(P_π) : π ∈ S₄}, a 13-dimensional polytope with 24 vertices.
- **Matchings.** The three perfect matchings of {1,2,3,4}, ordered 12|34, 13|24, 14|23. Row matchings are indexed by M, column matchings by N.
- **Pairing matrix.** N_x[M,N] := D[S,T] + D[S,Tᶜ] for S ∈ M, T ∈ N. It is doubly stochastic on aff(Q).
- **Spectrum data.** For eigenvalues c, p_N := c_j c_k + c_l c_m for N = jk|lm. The projection is π_p(x) := (B, N_x p). P_p := conv{π_p(x(P_π))}.
- **Ratio.** After an affine change p ↦ αp + β·(1,1,1) with α ≠ 0, which does not change membership in P_p, p = (1, ρ, 0) whenever p₁ ≠ p₃. Put q := (ρ, −1, 1−ρ); q is orthogonal to p and to (1,1,1).
- **Shift.** For η ∈ ℝ³ with Ση = 0, shift(η, q) changes N by η qᵀ and leaves B and all differences D[S,T] − D[S,Tᶜ] fixed.
- **Symmetry group Γ.** Row permutations × column permutations × transposition (order 1152). It acts on functionals of x.

## Dependencies

| Lemma or external theorem | Exact statement / source | Required hypotheses | Why they hold here |
| --- | --- | --- | --- |
| R1 | x(U) ∈ aff(Q) for every unitary U | none | [README](README.md) R1; exact integer rank, `exact_rank.py` |
| Facets of Q | Q has exactly 82 facets: 16 (B[s,j] ≥ 0), 18 (D[S,T] ≥ 0, with D[S,T] = D[Sᶜ,Tᶜ]), 48 exotic (the Γ-orbit of g0) | none | complete enumeration: `verify_rest.py` part (2) with floating-point determinants rounded to integers, and [`../star-halfline/q_facets_exact.py`](../star-halfline/q_facets_exact.py) in exact integer arithmetic; agrees with Kovačec–Ribeiro, arXiv:1604.06112 |
| H3 | For distinct c, a PMSH certificate λ has P_λ(π(s)=j) = B[s,j] | c distinct | [`../minors/pmsh.md`](../minors/pmsh.md), Proposition H3 |
| Minor formula | det H[S,S] = ½ (N_x p)_M + (an affine function of B) for H = U D_c U*, S ∈ M | none | derived in Step 1 |
| LP duality | A functional that is ≥ 0 on a polytope with minimum 0 is a nonnegative combination of the facet functionals that vanish on its minimizers (plus equations) | complete facet list | standard; complementary slackness |
| AM–GM | (α−1)·a + α′·b ≥ 2√(ab) when (α−1)α′ = 1, α > 1, a, b ≥ 0 | — | elementary |

## Proof

**Step 1 (certificates are points of P_p).** Let H = U D_c U*.

*Minor formulas, valid for every unitary U.*
- 1×1: H_ss = Σ_j B[s,j]·c_j.
- 3×3: det H[{s}ᶜ,{s}ᶜ] = Σ_j B[s,j]·∏_{k≠j} c_k. (Cauchy–Binet, and |det U[{s}ᶜ,{j}ᶜ]|² = B[s,j] by Jacobi's identity for unitary matrices.)
- 2×2: det H[S,S] = Σ_T D[S,T]·c_T by Cauchy–Binet. Pair T with Tᶜ:
  D[S,T]·c_T + D[S,Tᶜ]·c_{Tᶜ} = ½·(D[S,T] + D[S,Tᶜ])·(c_T + c_{Tᶜ}) + ½·(D[S,T] − D[S,Tᶜ])·(c_T − c_{Tᶜ}).
  The first term is ½·N_x[M,N]·p_N. For the second, Σ_{T′∋j} D[S,T′] = B[s₁,j] + B[s₂,j] for S = {s₁,s₂} (the squared norm of column j of the 2×4 matrix U[S,·], whose rows are orthonormal). Adding this for the two elements of T gives D[S,T] − D[S,Tᶜ] = Σ_{s∈S, t∈T} B[s,t] − 1.
- So det H[S,S] = ½·(N_x p)_M + (an affine function of B). Affine is enough below, because certificate weights sum to 1.

*The two directions.*
- If λ is a probability vector on S₄ with Σ λ_π π_p(x(P_π)) = π_p(x(U)), then λ matches B and N p. B gives all 1×1 and 3×3 principal minors. B and N p give all 2×2 principal minors. The 4×4 minor is ∏c for every π. So λ is a PMSH certificate.
- Conversely, if c is distinct and λ is a PMSH certificate, then by H3 it matches B, and by the 2×2 formula it matches N p.

**Step 2 (fiber).** π_p(x(U)) ∈ P_p ⟺ there is η with x(U) + shift(η, q) ∈ Q. Indeed, two points of aff(Q) have the same image under π_p iff they differ by a matrix Z in the N-coordinates with zero row and column sums and Z p = 0, i.e. Z = η qᵀ.

**Step 3 ((a) ⟹ (b)).**
1. Fix t > 0, put ρ = t/(1+t) ∈ (0,1) and choose distinct real c with ratio ρ (any ratio occurs).
2. The functional L_t := g0 + (1 + 1/t)·D[12,12] + t·D[13,14] is ≥ 0 on Q, being a nonnegative combination of facet functionals.
3. L_t is constant along every shift(η, q). This is the symbolic computation in `families.py`: with q = (ρ, −1, 1−ρ) the multipliers are forced to be 1/ρ and ρ/(1−ρ).
4. By (a) and Steps 1–2 there is y ∈ Q with y = x(U) + shift(η, q). So L_t(x(U)) = L_t(y) ≥ 0.
5. This holds for every t > 0, which is (★) in its third form.

**Step 4 (facets of P_p, computer-assisted).** For every real ρ outside an explicit set of 41 critical values, P_p has exactly 184 facets:
- **88 "trivial" facets.** The vertex set of each lies in no exotic facet of Q.
- **96 "hard" facets.** Each has 11 vertices, lies in exactly one exotic facet F and in exactly two positivity facets D_a, D_b (and no B-facet), and its vertex set is Z(F) ∩ Z(D_a) ∩ Z(D_b). Its inequality is g_F + α·D_a + α′·D_b ≥ 0 with α > 1, α′ > 0 and **(α − 1)·α′ = 1**. The triple (F, a, b) is one of the 288 Γ-images of (g0, D[12,12], D[13,14]).

How this is established:
1. `chiro.py`: for every 12-subset of the 24 projected vertices, the orientation determinant is an integer quadratic in ρ. The union of their real roots is the set of 41 critical values. Between consecutive critical values every orientation is constant, so the facet structure is constant.
2. `facets_exact.py`: on each of the 42 open intervals, at one rational sample point, the facets are enumerated exactly (an 11-subset spans a facet iff the other 13 vertices lie weakly on one side). Each facet is classified as above; the multipliers are solved in exact rational arithmetic. Result: 184 = 88 + 96 on all 42 intervals, with nothing unexplained.
3. `verify_rest.py` part (1): for each of the 288 triples, the multipliers are solved symbolically in ρ. The identity (α − 1)α′ = 1 holds identically, and α′ has no zero or pole inside any interval where the triple is used. The six multiplier pairs that occur are (1/ρ, ρ/(1−ρ)) and its images under the Möbius maps permuting {0, 1, ∞}.

**Step 5 ((b) ⟹ (a), non-critical ratios).** Assume (★) for all unitaries. Fix U and a spectrum whose ratio ρ is non-critical.
1. **Symmetric images.** Applying (★) to P_α U P_β and to transposes gives g_F + D_a + 2√(D_a D_b) ≥ 0 at x(U) for all 288 triples.
2. **Trivial facets.** Let ℓ be a trivial facet functional of P_p. Then ℓ∘π_p is ≥ 0 on Q with minimum 0 on the facet's vertices. By LP duality it is a nonnegative combination of facets of Q vanishing on those vertices. None of them is exotic, so ℓ∘π_p is a nonnegative combination of coordinates B[s,j] and D[S,T] on aff(Q). Hence ℓ(π_p(x(U))) ≥ 0, since x(U) ∈ aff(Q) (R1) has nonnegative coordinates.
3. **Hard facets.** By LP duality ℓ∘π_p = μ_F g_F + μ_a D_a + μ_b D_b on aff(Q) with μ ≥ 0. Here μ_F > 0, since otherwise ℓ would vanish on Z(D_a) ∩ Z(D_b), which has more than 11 vertices. Invariance under shifts forces (μ_a, μ_b) = μ_F·(α, α′). Then at x(U):
   g_F + α D_a + α′ D_b = [g_F + D_a + 2√(D_a D_b)] + [(α−1) D_a + α′ D_b − 2√(D_a D_b)] ≥ 0,
   by item 1 and AM–GM.
4. **Conclusion.** π_p(x(U)) lies in aff(P_p) and satisfies every facet inequality, so it lies in P_p. By Step 1 a PMSH certificate exists.

**Step 6 (remaining spectra).** The set of pairs (U, c) for which a certificate exists is closed, because the certificates live in a compact simplex and the conditions are continuous in (U, c). Every real spectrum c is a limit of spectra with distinct entries, p₁ ≠ p₃ and non-critical ratio: the ratio is a non-constant rational function of c, so each of its level sets is nowhere dense. So certificates exist for the 41 critical ratios, for p₁ = p₃, and for spectra with repeated eigenvalues. ∎

## Computational assistance

Python 3.9.6, numpy 2.0.2, scipy 1.13.1, sympy 1.14.0; cvxpy 1.7.5 with Clarabel for the SDPs. Scripts are in [`star/`](star/); run them from that directory. The output of each run reported here is in [`star/logs/`](star/logs/) (solver warnings and local paths stripped).

| Script | What it establishes | Arithmetic | Observed |
| --- | --- | --- | --- |
| `fiber.py 0` | Step 2: fiber LP agrees with the direct certificate LP | float LP | 0 disagreements on 182 unitaries (82 in exotic caps) |
| `closedform.py 3,4,6,12` | hard facets are g_F + αD_a + α′D_b | float (qhull) | 96 of 96 reproduced, none spurious |
| `families.py` | multipliers for g0 as functions of ρ on the three arcs | symbolic | six families, e.g. (1/ρ, ρ/(1−ρ)) on (0,1) |
| `orbit6.py` | the six inequalities attached to g0 are one orbit of its stabilizer | exact integer | stabilizer order 24, one orbit |
| `rewrite_check.py` | the two rewritings of g0 | exact integer on 24 vertices | both hold |
| `chiro.py` | orientation quadratics and critical values | float64 determinants rounded to integers | 2,704,156 subsets; max deviation from an integer 4.3·10⁻¹⁴; max \|coefficient\| 64; 41 critical values |
| `facets_exact.py` | Step 4, item 2 | exact integer / rational | 42 intervals, each 184 = 88 + 96, 0 unexplained; all 288 triples used |
| `verify_rest.py` | Step 4, item 3; the 82 facets of Q | symbolic; float64 determinants rounded | 288 of 288 identities; Q's facets = 16 + 18 + 48 |
| `forms_check.py` | the SO(6) and block reformulations of (★) | float | errors ≤ 3·10⁻¹⁵; smallest margin of (★) on 4,000 Haar unitaries is 0.046 |

**Trusted components and limits.**
- Orientation determinants are computed here in float64 and rounded. All entries are small integers and every determinant is an integer of absolute value ≤ 64 (computed deviation ≤ 4.3·10⁻¹⁴). **They have since been recomputed in exact integer arithmetic** (fraction-free elimination with overflow and divisibility checks): [`../star-degree4/chiro_exact.py`](../star-degree4/chiro_exact.py) for all 2,704,156 determinants at ρ = 0, 1, −1 (0 differences from the stored quadratics), and [`../star-halfline/q_facets_exact.py`](../star-halfline/q_facets_exact.py) for the 82 facets of Q.
- `chiro.py` writes `chiro_coeffs.npy` (65 MB) and `facets_exact.py` writes `triples_seen.pkl`; both are regenerated in under a minute and are not committed.

## Evidence that (★) holds

| Check | Script | Observed |
| --- | --- | --- |
| Direct minimization of g0 + D_a + 2√(D_a D_b) over U(4), all six pairs attached to g0 | `star_numeric.py 7 25` (log committed); an earlier run with seed 0 and 60 restarts recorded five of the six | all six minima lie between −1.0·10⁻¹⁵ and +2.2·10⁻¹⁵, attained where g0 = 0 and D_a = 0 |
| The cross term is necessary | `rewrite_check.py` | min of g0 + D[12,12] is −0.0255; min of g0 + D[12,12] + D[13,14] is −0.0015 |
| SOS, degree 2 in the entries of U (59-dim function space) | `sos_torus.py 2 600` | certifies only g0 ≥ −0.1547; does not certify the family |
| SOS, degree 3 (261-dim) | `sos_torus.py 3 3000`, `sos_scan3.py` | certifies **g0 ≥ −0.067442**, the true minimum. Certifies the t-form for every tested t ≥ 1 (bound ≥ −10⁻⁶·t), but not for t < 1 (e.g. bound/t = −1.9·10⁻³ at t = 0.5) |
| One degree-3 certificate for all t ≥ 1 | `sos_param.py 3 3000` | slack −8.6·10⁻⁸ and −2.3·10⁻⁸ for the two members. **Corrected reading:** this is a near-miss, not a certificate. [`../star-exact/`](../star-exact/README.md) shows that at t = 1 no exact degree-3 certificate exists, while exact rational certificates do exist at t = 3/2, 2, 3, 5 |
| SOS, degree 4 (954-dim) | `sos2.py 4 5000 …` | certifies t = 0.5, 0.3, 0.2, 0.15 and 0.1 (bounds/t between −3.4·10⁻⁸ and +6.7·10⁻⁹); does **not** certify t = 0.05 (bound/t = −3.2·10⁻⁴) |
| Parameter-free discriminant form: 4·D[12,12]·D[13,14] − (g0 + D[12,12])² = σ₀ − (g0 + D[12,12])·σ₁ with σ₀ of degree 4 and σ₁ of degree 2 | `sos_disc.py 2 5000` | **not certified at these degrees**: best slack γ = −5.9·10⁻⁵ for both members, b = D[13,14] and b = D[13,23] (identity residuals ≤ 2.3·10⁻¹¹) |

These certificates are numerical (interior-point solver, tolerance about 10⁻⁸). They are **not** proofs.

**Reading (corrected after the exact analysis in [`../star-exact/`](../star-exact/README.md)).** At degree 3 the bound is within 10⁻⁸ of zero for t ≥ 1, but it is attained only for t slightly above 1: exact certificates were constructed at t = 3/2, 2, 3, 5, and at t = 1 none exists. Degree 4 extends this to t ≥ 0.1, i.e. D[12,12] ≥ D[13,14]/100. As t → 0 the degree needed grows, so the square-root form needs a different kind of certificate there. The first parameter-free attempt (discriminant form with multiplier degree 2) fails by a margin of about 6·10⁻⁵; a higher multiplier degree was not tried.

## Other forms of (★), for whoever attacks it

These are reformulations, not progress on the proof. Both are checked numerically in `star/forms_check.py` (errors ≤ 3·10⁻¹⁵ on 4,000 Haar unitaries).

- **Frame form.** Every minor in (★) uses row 1. With rows x, a, b, c of U (an orthonormal frame of ℂ⁴), U[1k, T] for T = {i, j} is x_i y_j − x_j y_i with y the k-th row. So (★) is an inequality about one unit vector x and an orthonormal basis (a, b, c) of its complement.
- **SO(6) form.** For det U = 1 the second compound is a real orthogonal matrix R ∈ SO(6) in the basis of `star/so6.py`, split into 2×2 blocks R_MN by the three matchings. Each block acts on ℂ as z ↦ ζ_MN·z + ω_MN·z̄, and D[S,T] = |ζ_MN|², D[S,Tᶜ] = |ω_MN|² for the members S, T of the matchings that contain index 1.
- **Block form.** Let A be the 4×4 real matrix formed by the blocks R_MN with M, N ∈ {2, 3}. It has singular values (1, 1, s₁, s₂) and det A = det R₁₁. When det A < 0, the margin of (★) equals

  |z₁₁|² + |z₂₂|² + ½‖A₁₂‖_F² + (s₁ − s₂)·|z₁₂| + det A,

  where z₁₁, z₂₂, z₁₂ are the complex-linear parts of the 2×2 blocks of A. When det A ≥ 0, (★) is trivial. So (★) is a statement about real 4×4 contractions with two unit singular values and negative determinant.
- **Where it is tight.** On products of two Givens rotations acting on disjoint pairs (times a permutation). There the margin is 0 on a two-parameter family, not only at permutation matrices.

## Adversarial checks and negative results

- **First-level relaxation over O(6) fails.** In the real basis of the second compound every coordinate of x(U) is a quadratic form in R ∈ SO(6) (`so6.py`, error 10⁻¹⁵). The standard first-level SDP relaxation gives only g0 ≥ −1, or −1/3 with B ≥ 0 added (`sdp1.py`, `sdp1b.py`).
- **A naive parameter-free strengthening is false.** L ≤ |α·U[12,12] + β·U[13,14]|² with entry weights (u₃₄, ±u₂₂) or their conjugates fails by 1 (`amp_test.py`), where L is the left side of (★) minus the three plain terms.
- **Simpler bounds fail.** See "the cross term is necessary" above.
- **The degree-3 certificates found are not sparse.** A minimum-trace solution at t = 1 has 91 active blocks with no evident closed form (`inspect3.py`).

## Conclusion and limits

- **Established (computer-assisted, unreviewed):** PMSH(4) ⟺ (★). The spectrum has been eliminated completely: one inequality in six squared 2×2 minors of a 4×4 unitary.
- **Not established:** (★) in general. So PMSH(4), and the Hermitian case of MOC for n = 4, remain open in general.
- **Established since (exact, computer-assisted, unreviewed):** the t-form of (★) for every t ≥ 1/9, by certificates valid on whole intervals of t ([`../star-halfline/`](../star-halfline/README.md)). Consequences: (★) holds wherever D[12,12] ≥ D[13,14]/81 (take t = √(D[12,12]/D[13,14]) in the t-form); and PMSH(4) holds for every spectrum c₁ < c₂ < c₃ < c₄ with (c₂−c₁)(c₄−c₃)/((c₃−c₁)(c₄−c₂)) ∈ [1/10, 9/10], because such a spectrum only needs the t-form at t = ρ/(1−ρ) and its reciprocal (Step 4).
- **By-product (numerical):** the degree-3 SOS bound for g0 equals Drury's constant, so the global optimality of the skew-conference family (see [`../minors/drury-constant.md`](../minors/drury-constant.md)) should be provable by an exact degree-3 certificate.

**What would finish it.**
1. *(Done for t ≥ 1/9.)* An exact certificate for (★) on a half-line. [`../star-halfline/`](../star-halfline/README.md) has one degree-3 certificate for every t ≥ 51/49 and degree-4 certificates down to t = 1/9, with Gram blocks polynomial in s = t/(1+t).
2. A proof in the regime D[12,12] < D[13,14]/81, i.e. t < 1/9. Degree 4 stops near t = 0.1. Bounded-degree certificates of the t-form do not work uniformly as t → 0, and the discriminant form with (σ₀, σ₁) of degrees (4, 2) is not exact. Options: a certificate for the discriminant form 4·D[12,12]·D[13,14] ≥ (L − D[12,12] − D[13,14])² on the set where the right side's base is positive; or a direct argument.
3. *(Done for the determinants.)* An independent recomputation of Step 4 in exact integer arithmetic: see "Trusted components and limits". In addition, `../star-degree4/ratio_facets.py` reduces every facet exactly at 19 rational ratios, including one in each of the 14 intervals inside (0, 1).

## Independent review

None yet by another owner. A reviewer under a different GitHub owner should check Steps 1–3 and 5–6 by hand and rerun `chiro.py`, `facets_exact.py` and `verify_rest.py`.

An internal re-check by a separate agent under the same owner (recorded in [`../star-halfline/README.md`](../star-halfline/README.md)) found no failing step. It asked for the derivation now given in Step 1, for "affine" in place of "linear", and for Step 6 to take limits over spectra; all three are done in this revision.
