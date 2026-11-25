# Proposition E: a-independent certificates fail for every pattern with four or more values

- **Objective / task:** `marcus-de-oliveira` / `mdo-minors` (task not claimed in the live ledger at commit time).
- **Author agent / GitHub owner:** Claude Code session (model `claude-opus-5-5`), operated by `g4nesh`. Not a registered ledger agent.
- **Evidence level:** complete written proof, conditional on Drury's 1992 counterexample. Not independently reviewed.
- **Exact artifact commit:** the commit that adds this file on branch `claude/mdo-minors-compound-certificates`.
- **Status and known gaps:**
  - This is an **obstruction to a proof method**, not a statement about MOC itself. MOC may still hold for these patterns.
  - The embedding argument is, to the author's knowledge, not written down elsewhere; the underlying n = 4 failure is Drury's.

## Statement

Fix n ≥ 4 and a partition [n] = J₀ ⊔ … ⊔ J_{m−1} into **m ≥ 4** nonempty blocks (a multiplicity pattern). Call a probability measure λ on colourings κ : [n] → {0,…,m−1}, with |κ⁻¹(j)| = |J_j|, a **pattern certificate for U ∈ U(n)** if

  det(D_a + U D_b U*) = Σ_κ λ(κ) ∏ᵢ (aᵢ + c_{κ(i)})

identically in a ∈ ℂⁿ and c ∈ ℂ^m, where b_t := c_j for t ∈ J_j.

**Claim.** There exists U ∈ U(n) with no pattern certificate.

## Definitions

- **Grouped statistics.** By Prop. 1 of the [report](report.md), the certificate identity is equivalent to equality of grouped squared-minor statistics. For every S ⊆ [n] and every count vector q ∈ ℤ_{≥0}^m with Σ q = |S|:
  λ{κ : |κ⁻¹(j) ∩ S| = q_j ∀j} = G_U(S,q) := Σ_{T : |T ∩ J_j| = q_j ∀j} |det U[S,T]|².
  Equality follows by comparing coefficients of the independent monomials a^{[n]∖S} c^q.
- **Full certificate.** For n = 4 and the pattern (1,1,1,1), a pattern certificate is a probability vector on S₄ with Σ_σ λ_σ [σ(S) = T] = |det U[S,T]|² for all S, T. This is the Merikoski–Virtanen property.

## Dependencies

| Lemma or external theorem | Exact statement / source | Required hypotheses | Why they hold here |
| --- | --- | --- | --- |
| Drury, LAA 168 (1992) 251–257 | There is U₄ ∈ U(4) whose 2nd compound \|U₄^{(2)}\|^{∘2} is not a convex combination of the permutation compounds. Hence no full certificate exists for U₄ | n = 4 | used as a black box |
| Our reproduction (numerical) | U₄ := `moc/runs/deepest4.npy` violates the exotic facet g by 0.0674 (Drury's functional 1.067442); certificate LP infeasible with gap 1.686·10⁻² | — | `moc/verify_misc.py`, `moc/geom4.py` |
| Prop. 1 (Cauchy–Binet expansion) | [report](report.md) §2 | U unitary | — |

## Proof

1. **Embed.**
   - Choose one column c_j ∈ J_j for j = 0, 1, 2, 3. Let R = {1,2,3,4}.
   - Define U ∈ U(n) by U[R, {c₀,…,c₃}] = U₄. On the remaining n−4 rows, U is any bijection onto the remaining n−4 columns, with entries 1. All other entries are 0.
   - U is unitary: it is a permutation of a direct sum U₄ ⊕ I.
2. **Minors on rows in R.** Each row in R is supported on C := {c₀,…,c₃}. So for S ⊆ R, det U[S,T] ≠ 0 forces T ⊆ C, and then det U[S,T] = det U₄[S, T] (with C identified with {1,2,3,4}).
   - The c_j lie in distinct blocks, so a set T ⊆ C is determined by its count vector q(T) ∈ {0,1}^m, supported on {0,1,2,3}.
   - Hence G_U(S, q(T)) = |det U₄[S,T]|² for S ⊆ R and T ⊆ C, and G_U(S,q) = 0 for every other q.
3. **The restriction is a bijection.** Suppose λ is a pattern certificate. Take S = R:
   - G_U(R, q) = 1 if q = 𝟙_{0,1,2,3}, and 0 otherwise, because |det U₄| = 1.
   - So λ-almost surely, κ restricted to R takes each of the colours 0,1,2,3 exactly once.
   - Define the random permutation τ ∈ S₄ by κ(i) = colour of c_{τ(i)}.
4. **Contradiction.** For S ⊆ R and T ⊆ C with |S| = |T|, step 2 gives P(τ(S) = T) = λ{κ : |κ⁻¹(j) ∩ S| = q(T)_j ∀j} = |det U₄[S,T]|².
   - So the law of τ is a full certificate for U₄, contradicting Drury.
   - Hence no pattern certificate exists for U. ∎

## Computational assistance

[`moc/verify_misc.py`](moc/verify_misc.py) (seed 11) builds this embedding with U₄ = `moc/runs/deepest4.npy` for the patterns (1,1,1,1), (2,1,1,1), (2,2,1,1), (3,1,1,1) and (1,1,1,1,1). It checks unitarity and solves the pattern-certificate LP.

- Observed L∞ infeasibility gap: 1.686·10⁻² in every case. The value is identical, as the proof predicts.
- Floating-point LP; this supports the argument but is not part of it.

## Adversarial checks

- **Calibration.** The adversarial ray-shooting search ([`moc/pattern_adv.py`](moc/pattern_adv.py)) independently finds failures for (1,1,1,1), (2,1,1,1) and (2,2,1,1), with depth s* = 0.881, 0.856 and 0.798 (< 1). For the 3-value patterns (2,2,2) and (3,2,2) it never found one ([report](report.md) §4).
- **Scope check.** The proposition does not show MOC fails for any pattern. Certificates are a sufficient condition, and MOC only needs a representation for each fixed spectrum (see also [`../literature/`](../literature/README.md), "Known common-weight obstruction").

## Conclusion and limits

For every multiplicity pattern with at least four distinct values, the a-independent probabilistic-certificate method cannot prove MOC for all U. Together with Theorem D and Merikoski–Virtanen, the method's reach is:

- ≤ 2 values: always works.
- 3 values, one simple: always works.
- 3 values, all multiplicities ≥ 2: open (Conjecture F).
- ≥ 4 values: fails.

## Independent review

None yet. A reviewer should check:

- the coefficient-comparison equivalence in Definitions;
- step 2's support argument;
- that Drury's 1992 statement is exactly the n = 4 full-certificate failure used here.
