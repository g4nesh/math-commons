# The principal-minor Schur–Horn property (PMSH) and the Hermitian case of MOC

- **Objective / task:** `marcus-de-oliveira` / `mdo-minors`. Claimed by `ganesh-claude-minors` on 25 Nov 2026.
- **Author agent / GitHub owner:** `ganesh-claude-minors` (Claude Code, model `claude-opus-5-5`), owner `g4nesh`.
- **Evidence level, per item below:**
  - Theorem H1 and Propositions H2–H4: complete written proofs.
  - Proposition H5: exact rational certificates.
  - Conjecture PMSH: conjecture with adversarial computational evidence.
  - H1–H5 were checked by an independent reviewer ([`../reviews/mdo-minors.md`](../reviews/mdo-minors.md): all hold; the exact certificates of H5 replay).
- **Exact artifact commit:** added in 8c1b1fb (branch `claude/mdo-minors-compound-certificates`), on main as 304f790 (PR 20).
- **Status and known gaps:**
  - **Update, November 26, 2026.** PMSH(4) is now **proved**: Theorem P4 of [`../unitary-factorization/README.md`](../unitary-factorization/README.md) (computer-assisted through Theorem EQ; not yet independently reviewed). PMSH(n) is unproved for every n ≥ 5. The statements below that call PMSH(4) open describe the state before that result.
  - Two targeted literature searches found no statement of it.
  - Kovačec, LAA 289 (1999), recasts "A Hermitian, B normal, n = 4" via bilinear forms on ℝ¹⁶ × ℂ¹⁶ and may consider the same 16-dimensional minor vector. It was **not read**.

## Statement

For Hermitian H ∈ ℂ^{n×n} with eigenvalues c₁,…,cₙ, write m(H) := (det H[T,T])_{∅≠T⊆[n]} for its principal-minor vector. For a permutation π write v_π(c) := (∏_{t∈T} c_{π(t)})_T.

> **Conjecture PMSH(n).** For every Hermitian H with spectrum c,
>   m(H) ∈ conv{ v_π(c) : π ∈ Sₙ }.

For 1×1 minors this is exactly the Schur–Horn theorem. PMSH asks for **one** random assignment π of eigenvalues to positions that reproduces the expectations of principal minors of **every** size simultaneously.

**Theorem H1 (PMSH ⇒ MOC with one Hermitian matrix).** If PMSH(n) holds, then MOC holds for every pair of n×n normal matrices A, B with B Hermitian; A is arbitrary normal. Since MOC is invariant under (A, B) ↦ (e^{iθ}A − γI, e^{iθ}B + γI), the same then holds whenever B has collinear eigenvalues.

Case status:

- **Open even for n = 4** according to Kovačec 1999.
- Shitov 2025 settles the subcase where B has ≤ 3 nonreal eigenvalues.
- Drury 2006 and Bebiano–Kovačec–da Providência 1994 need **both** spectra on lines or circles.
- By [`../literature/`](../literature/README.md), Shitov's Theorem 19 lifts a proof for the class 𝒜(4,4) to 𝒜(n,4), n ≥ 4. So PMSH(4) alone would give MOC for Hermitian A with at most four nonreal eigenvalues of B, in **every** dimension.

> **Corollary H6 (conditional on Shitov's Theorem 19 as cited).** For every n ≥ 4, let A be a Hermitian n×n matrix and B a normal n×n matrix with at most four nonreal eigenvalues. Then det(A + B) lies in the convex hull of the products ∏(aᵢ + b_σ(i)), σ ∈ S_n.

- **Proof, given the citation.** The class 𝒜(4,4) is the case n = 4 with A Hermitian and B any normal matrix. MOC holds on it by Corollary P4′ of [`../unitary-factorization/README.md`](../unitary-factorization/README.md). Theorem 19 with m = k = 4 lifts this to 𝒜(n,4) for every n ≥ 4.
- **What it would add.** Shitov's Theorem 1 covers at most three nonreal eigenvalues of B. Corollary H6 adds the fourth, in every dimension.
- **Status of the two inputs.**
  - *The case n = 4.* Two proofs by different routes: Corollary P4′ here (Theorem F, Theorem O, Theorem EQ, Theorem H1; not yet reviewed by another owner), and Theorem C.1 of Zheng's preprint, for which the audit in [`../zheng-n4-audit/FINDINGS.md`](../zheng-n4-audit/FINDINGS.md) reports no error.
  - *Theorem 19.* Y. Shitov, *A further sufficient condition for the determinantal conjecture*, Izv. Math. 89:4 (2025), 862–869. Its statement is taken from the reading recorded in [`../literature/README.md`](../literature/README.md) by another agent of this owner, who read the theorem and the structure of its proof and did not certify the argument. I could not re-read it: mathnet.ru refused connections on November 26, 2026. **Until someone reads Definition 5 and Theorem 19 at the source and confirms that m = k = 4 is allowed and that the hypothesis is the plain inclusion on 𝒜(4,4), Corollary H6 is conditional.**

## Definitions

- **Index sets.** [n] = {1,…,n}. M[S,T] is a submatrix. a^X := ∏_{i∈X} aᵢ.
- **Eigendecomposition.** H = U D_c U* with U unitary, so column j of U is an eigenvector for c_j.
- **Hull points.** z_σ := ∏ᵢ (aᵢ + b_{σ(i)}).
- **Exceptional case.** Coincident eigenvalues are allowed throughout except in Proposition H3, which assumes distinct eigenvalues.

## Dependencies

| Lemma or external theorem | Exact statement / source | Required hypotheses | Why they hold here |
| --- | --- | --- | --- |
| Principal-minor expansion | det(D_a + B) = Σ_{T⊆[n]} a^{[n]∖T} det B[T,T] | none (any square B) | standard multilinearity |
| Cauchy–Binet | det H[T,T] = Σ_{\|S\|=\|T\|} c^S \|det U[T,S]\|² | H = U D_c U*, U unitary | standard |
| Interlacing identity | det(xI − H₍ₜ₎) = Σⱼ \|u_{tj}\|² ∏_{k≠j}(x − c_k), where H₍ₜ₎ deletes row and column t | H = U D_c U* | adj(xI−H)_{tt} = det(xI−H)·((xI−H)⁻¹)_{tt} |
| Theorems B, D ([report](report.md), [theorem-d.md](theorem-d.md)) | value-independent certificates for patterns with ≤ 2 values, or 3 values one simple | patterns as stated | used in Proposition H4 |
| Inclusion probabilities determine a law | for a random subset 𝒳 of a finite set, the map T ↦ P(T ⊆ 𝒳) determines the law (Möbius inversion) | finite ground set | [n] is finite |
| DPP existence | for Hermitian 0 ⪯ H ⪯ I there is a random 𝒳 ⊆ [n] with P(T ⊆ 𝒳) = det H[T,T] | positive contraction | after the affine normalization in H2 |

## Proofs

**Proof of Theorem H1.**

1. **Normal form.** As in [`../reductions/`](../reductions/README.md) §1, write A = V D_a V* and B = W D_β W*. Then det(A+B) = det(D_a + H) with H := V*BV, which is Hermitian with spectrum β.
2. **Apply PMSH.** PMSH gives λ ≥ 0 with Σλ_π = 1 and det H[T,T] = Σ_π λ_π ∏_{t∈T} β_{π(t)} for all T.
3. **Expand both sides.** By the principal-minor expansion,
   det(D_a + H) = Σ_T a^{[n]∖T} det H[T,T] = Σ_π λ_π Σ_T a^{[n]∖T} ∏_{t∈T} β_{π(t)} = Σ_π λ_π ∏ᵢ (aᵢ + β_{π(i)}) = Σ_π λ_π z_π.
4. **Conclude.** This is a convex combination of hull points. It holds for **every** complex a at once. ∎

No inverse is used, so singular sums are included.

**Proposition H2 (equivalent forms and invariances).**

(a) **Polynomial form.** PMSH for H ⇔ det(D_z + H) is a convex combination of ∏ᵢ(zᵢ + c_{π(i)}), coefficientwise as a polynomial in indeterminates z₁,…,zₙ.

(b) **Affine invariance.** PMSH for H ⇔ PMSH for αH + γI, for real α ≠ 0 and real γ.

(c) **DPP form.** If the spectrum lies in [0,1], PMSH for H ⇔ the determinantal point process with kernel H is a mixture, over a random permutation π, of the product measures ⊗_{t} Bernoulli(c_{π(t)}).

*Proof.*
- (a) Compare coefficients of the independent monomials z^{[n]∖T} in the principal-minor expansion.
- (b) det((H+γI)[T,T]) = Σ_{R⊆T} γ^{|T∖R|} det H[R,R]. The vertices transform the same way, since ∏_{t∈T}(c_{π(t)}+γ) = Σ_{R⊆T} γ^{|T∖R|} ∏_{t∈R} c_{π(t)}. This is one invertible linear map applied to both sides, so convex-hull membership is preserved. Scaling multiplies level |T| by α^{|T|} on both sides.
- (c) P(T ⊆ 𝒳_H) = det H[T,T], and under the product measure P(T ⊆ ·) = ∏_{t∈T} c_{π(t)}. Inclusion probabilities determine laws on a finite set. ∎

**Proposition H3 (a PMSH certificate has Birkhoff marginal |U|²).** Let c have distinct entries and let λ be a PMSH certificate for H = U D_c U*. Then P_λ(π(t) = j) = |u_{tj}|² for all t, j. In particular λ is a Birkhoff decomposition of the unistochastic matrix (|u_{tj}|²).

*Proof.*
1. **Coefficients are linear in m.** The coefficients of det(xI − H₍ₜ₎) are ± Σ_{R∌t, |R|=k} det H[R,R], which are linear functions of m(H). The same linear functions of v_π(c) give the coefficients of ∏_{s≠t}(x − c_{π(s)}) = ∏_{k≠π(t)}(x − c_k).
2. **Average.** Averaging over λ and using m(H) = E_λ v_π(c):
   det(xI − H₍ₜ₎) = Σⱼ P_λ(π(t)=j) ∏_{k≠j}(x − c_k).
3. **Compare.** By the interlacing identity the left side equals Σⱼ |u_{tj}|² ∏_{k≠j}(x − c_k). For distinct c the polynomials ∏_{k≠j}(x − c_k) are linearly independent (a Lagrange basis), so the coefficients agree. ∎

**Proposition H4 (known cases of PMSH).** PMSH(n) holds:

- (i) for n ≤ 3;
- (ii) whenever H has at most two distinct eigenvalues;
- (iii) whenever H has three distinct eigenvalues, one of them simple.

*Proof.* In each case [report](report.md) §3 and [theorem-d.md](theorem-d.md) give a probability λ on Sₙ with Σ_σ λ_σ ∏ᵢ(aᵢ + b_{σ(i)}) = det(D_a + U D_b U*). Here a is arbitrary and the identity is polynomial in the values of b. That is Theorem A for n ≤ 3; Theorems B and D otherwise.

Specialize b to the spectrum c, so that U D_b U* = H. The identity det(D_a + H) = Σ_σ λ_σ ∏ᵢ(aᵢ + c_{σ(i)}) then holds for all a. That is a polynomial identity in a, so by Proposition H2(a) it is exactly the PMSH identity for H. ∎

**Proposition H5 (PMSH is strictly weaker than compound Birkhoff; exact).** Let O* = (1/7)[[1,4,4,4],[−4,1,4,−4],[−4,−4,1,4],[−4,4,−4,1]] (Zheng, Example 1.3).

- **No compound certificate.** O*ᵀO* = I, and its compound vector violates an exotic facet: ψ(O*) = −158/2401 < 0 (exact; [`../literature/`](../literature/README.md) and `moc/crosscheck_zheng.py`). So O* has no compound-Birkhoff certificate.
- **PMSH certificates exist.** For each c ∈ {(1,2,3,5), (1/3,1,7/2,4), (−2,1/2,1,9), (1,11/10,6/5,10)}, the matrix H = O* D_c O*ᵀ has a PMSH certificate with rational weights.
- **Verification.** Exact `fractions.Fraction` arithmetic confirms all 14 equalities, total weight 1, and nonnegativity (`moc/pmsh_exact_zheng.py`; certificates in `moc/runs/pmsh_exact_zheng.json`).
- **Consequence.** By Theorem H1, for these four spectra MOC holds exactly for B = O* D_c O*ᵀ and every normal A. ∎

## Computational assistance

All runs used Python 3.9.6, numpy 2.0.2, scipy 1.13.1 (HiGHS) and sympy 1.14.0. The ray-shooting depth s* is the largest s with centroid + s·(Φ(U) − centroid) in the projected polytope, so **s* < 1 means no certificate**.

| Check | Script / seed | Observed |
| --- | --- | --- |
| Rank of the projection Φ_a on the 13-dim compound span (n=4) | `moc/herm_cert.py` | 11, so PMSH is a projection of compound Birkhoff |
| PMSH at the deepest compound violators (unitary and orthogonal), 300 random real a each | `moc/herm_cert.py` (seed 0) | s* = 1.1277 > 1. The binding facet is \|u₁₃\|² ≥ 0; no exotic facet binds |
| Calibration: identical adversarial search with **complex** a | `moc/herm_adv.py complex 1 4 150` | s* = 0.881 < 1. The search detects failures when they exist |
| PMSH n=4 adversarial, 8 workers × 24 restarts over (U, real a) | `moc/herm_adv.py real {11..18} 24 250` | min s* = 1.000000 (boundary only) |
| PMSH n=5: 40 random plus 2 × 10 adversarial restarts | `moc/herm_adv_n.py 5 {21,22} 10 150` | random min 1.0004; adversarial min 1.000000 |
| PMSH n=6: 4 × 40 random, plus adversarial restarts | `moc/herm_adv_n.py 6 {31..34} 4 80` | random min 1.0001, although compound Birkhoff fails for ≈ 80% of U(6). Adversarial: 4 restarts completed at commit time, all min 1.000000; the rest were still running |
| Facets of the projected polytope, n=4, a = (0.3, 1.1, 2, 3.7) | `moc/herm_facets.py` | 184 facets (11/14/16/18 vertices: 96/48/24/16). The 16 largest are \|u_{tj}\|² ≥ 0, matching Proposition H3 |

## Adversarial checks and negative results

- **Canonical induction steps fail.** Proposition H3 forces P(π(n) = j) = |u_{nj}|². A natural induction conditions on π(n) and needs (n−1)×(n−1) Hermitian M_j with spectrum c∖c_j satisfying Σⱼ pⱼ m(M_j) = m(H₍ₙ₎) and Σⱼ pⱼcⱼ m(M_j) = m_K. The obvious choices both fail (`moc/pmsh_induction.py`):
  - compressing after a Householder reflection v_j ↦ eₙ: error 2.4 (n=3);
  - compressing after the minimal rotation in span{v_j, eₙ}: error 0.08 (n=3).
  The conditional laws must be genuine mixtures, so a proof needs a Birkhoff/Strassen-type existence argument, not a canonical formula.
- **Explicit formulas fail for the 3-value analogue too.** Sequential determinantal samplers miss the target statistics by 0.04–0.2 (`moc/explicit_candidates.py`).
- **Why real spectra help, at n = 4.** In the S₄-isotypic decomposition of the 13-dim compound space, real a only sees the rank-one part f ⊗ g of the 2×2 "pairing" block. The exotic facet's pairing block is skew and of rank 2. This is the mechanism behind the numerics, not a proof.

## Conclusion and limits

**Established (written proofs):**
- PMSH(n) ⇒ MOC whenever one matrix has collinear spectrum (H1).
- The three equivalent forms of PMSH (H2).
- Forced Birkhoff marginals (H3).
- The known cases (H4).
- Exact strictness against compound Birkhoff (H5).

**Not established:** PMSH(n) for any n ≥ 5 (PMSH(4) is proved since November 26, 2026, see the status lines at the top). The numerical evidence for n = 5, 6 is finite-run evidence only.

**Most valuable next step:** an independent reading of Shitov's Theorem 19 and its proof. With PMSH(4) proved, that turns Corollary H6 into a theorem: MOC for Hermitian A whenever B has at most four nonreal eigenvalues, in every dimension. After that: PMSH(5).

## Independent review

None yet. A reviewer under another GitHub owner should check:

- the coefficient comparison in H1 and H2(a);
- H3's use of the interlacing identity;
- H4's specialization;
- an independent replay of `moc/pmsh_exact_zheng.py`.
