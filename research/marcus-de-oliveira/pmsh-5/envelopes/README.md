# (E20) and (E21): one 3×3 block, and (★) one dimension up

Task `mdo-pmsh-5-envelopes` · agent `ganesh-claude-pmsh` · owner g4nesh · 27 November 2026. Continues [`../facets/README.md`](../facets/README.md), R5.

## Status

- **(E20) and (E21) are not proved on U(5).** They remain numerically true.
- **Proved (hand proof plus exact check; not independently reviewed):**
  - **Proposition B.** Each inequality depends only on one 3×3 block W of U. Its two "outside" terms are defect terms of W. So each is an inequality on 3×3 contractions.
  - **Proposition S.** On 3×3 blocks of 4×4 unitaries, both inequalities are images of (★) under row and column relabellings.
  - **Corollary S.** Both hold whenever W has at least two singular values equal to 1. In particular they hold whenever U ∈ U(5) splits off a 1×1 block outside the block W, after unitary changes of the two outside rows and of the two outside columns. This uses (★) on U(4), which is Corollary F2 of [`../../unitary-factorization/`](../../unitary-factorization/README.md) (written proof, not reviewed).
  - **Proposition M.** A closure-disc identity, which proves either inequality at W whenever a certain quadratic form is ≤ 0 in some admissible direction. Numerically this holds at 189 of 200 random (U, t) for (E20); it fails elsewhere.
- **So (E20) and (E21) are the n = 5 lifts of (★).** The open case is a 3×3 block with exactly one singular value equal to 1. That is precisely the gap between U(4) and U(5).
- **Numerical:**
  - Both inequalities hold on all 3×3 contractions, not only on blocks of U(5).
  - They have zeros at blocks of U(5) that are neither blocks of U(4) nor singular, for example with singular values (1, 0.892, 0.854).
  - A hull property (H) that would prove both is supported by a support-function test.

## Setting

- **Notation** as in [`../facets/README.md`](../facets/README.md): x_{K→J} = |det U[K,J]|². Rows are positions and columns are eigenvalue labels, numbered 0–4.
- **The two inequalities** (R5 there), with g the E-facet of each decomposition in its minimal form ([`../facets/facet_forms.py`](../facets/facet_forms.py)):
  - (E20) g₂₀ + Y + 2√(x01>12·Y) ≥ 0, with g₂₀ = x1>3 + x01>12 + x01>13 + x04>13 − x23>04 and Y = x04>01 + x04>14;
  - (E21) g₂₁ + A + 2√(A·Y) ≥ 0, with g₂₁ = x1>2 + x01>13 + x01>23 + x03>23 − x24>04, A = x01>23 and Y = x02>12 + x04>12.
  - Each is equivalent to the whole spectral family of its orbit of facets of P₅(c).
- **The block.** Rows R = (a, p, q) and columns C = (a′, γ, δ):
  - for (E20), R = (1, 0, 4) and C = (3, 1, 2);
  - for (E21), R = (1, 0, 3) and C = (2, 3, 1).
  - In both cases C = {1, 2, 3}: the three middle eigenvalues, which collide in the limit that defines the inequalities.
  - The other two rows and columns are "outside". Write W = U[R, C] and M(rows; cols) for squared minors of W.
- **The block forms.**
  - E(W) = |w_aa′|² + M(ap; a′γ) + M(ap; γδ) + M(pq; a′γ) − |det W|².
  - A₂₀ = M(ap; γδ) and Y₂₀ = |w_pγ|² + |w_qγ|² − M(pq; a′γ) − M(pq; γδ).
  - A₂₁ = M(ap; a′γ) and Y₂₁ = |w_pa′|² + |w_pδ|² − M(ap; a′δ) − M(pq; a′δ).
- **The defect forms.** Y₂₀ = f*(I − WW*)f with f = w̄_pγ·e_q − w̄_qγ·e_p. Y₂₁ = h*(I − W*W)h with h = w_pa′·e_δ − w_pδ·e_a′.
  - Both are ≥ 0 on every contraction.
  - They vanish when W is isometric along f, resp. h.

## Results

### Proposition B: each inequality is a function of W alone

> For every U ∈ U(5),
>
>   (E20) = E(W) + Y₂₀ + 2√(A₂₀·Y₂₀),  (E21) = E(W) + A₂₁ + 2√(A₂₁·Y₂₁),
>
> with W = U[R, C] for the respective R and C.

*Proof.*
1. **The E-facet.** By Jacobi's identity, |det W|² = |det U[Rᶜ, Cᶜ]|², the coordinate x23>04, resp. x24>04. With that, g₂₀ − E(W) and g₂₁ − E(W) are affine functions on aff(Q₅) that vanish on all 120 vertices. D(U) ∈ aff(Q₅) by R1 of the facets README, so both vanish at D(U).
2. **The outside terms.** By the identity Σ_{T∋j} |det U[S,T]|² = Σ_{s∈S} |u_sj|² for |S| = 2 (proof in the facets README, R4), with S = {p, q} and j = γ, Y = Σ_{o∉C} |det U[pq, γo]|² equals Y₂₀. For (E21), apply the transposed identity with the column pair {a′, δ} and row p.
3. **The square-root partners.** x01>12 is literally M(ap; γδ), and x01>23 is literally M(ap; a′γ).
4. **The defect forms.** Σ_{o∉C} |(f̃*U)_o|² = ‖U*f̃‖² − ‖W*f‖² = f*(I − WW*)f for f̃ = conj(u_pγ)·e_q − conj(u_qγ)·e_p. The same computation applied to U·h gives Y₂₁. ∎

- `reduce_check.py` checks steps 1–3 on all 120 vertices in integer arithmetic.
- It also compares the original and block forms on 2,000 Haar unitaries, with largest difference 2·10⁻¹⁵.

**Consequences.**
- Every 3×3 contraction is a block of some unitary of size at most 6. So the block forms make sense, and are numerically ≥ 0, on all contractions (`contraction_min.py`).
- The 3×3 block W of U ∈ U(5) has singular values (1, s₂, s₃), where (s₂, s₃) are the singular values of the outside block U[Rᶜ, Cᶜ] (CS decomposition).

### Proposition S: on blocks of 4×4 unitaries they are (★)

Write (★) as the statement that the margin

  D[13,13] + D[14,14] + D[13,23] + (√D[12,12] + √D[13,14])² − D[12,34]

is ≥ 0 for every V ∈ U(4), where D[S,T] = |det V[S,T]|². This is the form proved in Corollary F2 of [`../../unitary-factorization/README.md`](../../unitary-factorization/README.md).

> Let V ∈ U(4) with rows (a, p, q, o) and columns (a′, γ, δ, o′), and let W = V[apq, a′γδ]. Then
> - the block form of (E20) at W equals the (★) margin of V with rows 1, 2, 3, 4 ↦ p, q, a, o and columns 1, 2, 3, 4 ↦ γ, o′, a′, δ;
> - the block form of (E21) at W equals the (★) margin of V with rows ↦ q, o, a, p and columns ↦ o′, δ, a′, γ.

*Proof.*
1. **(E20), square-root partners.** In U(4), Σ_{T∋γ} D[pq, T] = |v_pγ|² + |v_qγ|² leaves one term: Y₂₀ = D[pq; γo′]. Under the relabelling this is D[12,12]. Likewise A₂₀ = D[pa; γδ] is D[13,14].
2. **(E21), square-root partners.** The transposed identity gives Y₂₁ = D[po; a′δ]. Jacobi's identity in U(4) gives D[qa; o′γ] = D[po; a′δ] and D[qo; o′δ] = D[ap; a′γ] = A₂₁. These are the relabelled D[13,14] and D[12,12].
3. **Linear parts.** In both cases the two linear parts agree on all 24 vertices of Q₄. So they agree on aff(Q₄), which contains every x(V) (R1 of [`../../pmsh/equivalence.md`](../../pmsh/equivalence.md)). ∎

- `reduce_check.py` checks steps 1–3 on the 24 vertices in integer arithmetic.
- On 2,000 Haar V ∈ U(4) the two sides agree to 5·10⁻¹⁴.

> **Corollary S.** (E20) and (E21) hold at every U ∈ U(5) whose block W has two singular values equal to 1, that is, whenever ‖U[Rᶜ, Cᶜ]‖ = 1. More generally they hold at every 3×3 contraction with at least two unit singular values, since such a matrix is a block of a 4×4 unitary.

*Proof.* Proposition B reduces both to W. Proposition S identifies them with (★) at a 4×4 unitary dilation of W. (★) holds on U(4) by Corollary F2 of `unitary-factorization`. ∎

**Corollary S′.** Both also hold when W is singular. Then −|det W|² = 0, and every other term is ≥ 0 on contractions.

### Proposition M: closure discs

- **Closing one outside channel.** Let U ∈ U(5). Change the two outside rows and the two outside columns by unitaries; this leaves W, hence both inequalities, unchanged. Then close the last outside row into the last outside column with a phase ζ:
  - V(ζ) = U[A, A] + U[A, o′]·ζ(1 − uζ)⁻¹·U[o, A], where A is the other four indices and u = U[o, o′], with |u| < 1.
  - For |ζ| = 1, V(ζ) is a 4×4 unitary.
  - Its block is F(ζ) = W + z·x yᵀ, with x = U[R, o′], y = U[o, C] and z = ζ/(1 − uζ).
- **Admissible directions.** As the outside changes vary, x runs over the range of I − WW* and yᵀ over that of I − Wᵀ·conj(W). These are the directions that keep the unit singular pair of W.
- **β.** Each minor is affine in z: det F[S,T] = det W[S,T] + z·β_{S,T}.

> **Proposition M.** For every linear function L of the 19 squared minors of a 3×3 matrix,
>
>   L(v(W)) = mean over |ζ| = 1 of L(v(F(ζ))) − L(|β|²)/(1 − |u|²).
>
> So if L ≥ 0 on blocks of U(4) and L(|β|²) ≤ 0 for one admissible (x, y), then L(v(W)) ≥ 0.

*Proof.* z = Σ_{k≥0} u^k ζ^{k+1}, so over the circle the mean of z is 0 and the mean of |z|² is 1/(1 − |u|²). Hence the mean of |α + zβ|² is |α|² + |β|²/(1 − |u|²). ∎

**Application.**
- (E20) = min_{t>0} L_t with L_t = E + (1 + t)Y₂₀ + A₂₀/t. (E21) = min_{t>0} L_t with L_t = E + (1 + 1/t)A₂₁ + t·Y₂₁. (AM–GM.)
- Each L_t is ≥ 0 on blocks of U(4) by Corollary S.
- So Proposition M proves L_t(W) ≥ 0 wherever some admissible direction gives L_t(|β|²) ≤ 0.

**Coverage** (`closure.py`, 40 Haar unitaries × 5 values of t):
- Such a direction exists at 189 of 200 (U, t) for (E20).
- At the other 11, L_t(W) ≥ 0.31: they are far from failure, but not covered.
- So closure discs alone do not prove either inequality.

### What is left, and a hull property that would finish it

- **The open case.** By Corollaries S and S′, what remains is W with singular values (1, s₂, s₃), 0 < s₃ ≤ s₂ < 1. These are blocks of U(5) that are neither blocks of U(4) nor singular.
- **Hull property (H).** For every such W, the vector v(W) of its 19 squared minors lies in the convex hull of v(F₂), where F₂ is the set of blocks of U(4) together with the singular contractions.
  - (H) implies (E20) and (E21). Both are of the form ℓ₀ + 2√(ℓ₁ℓ₂), with ℓᵢ linear in v and ℓ₁, ℓ₂ ≥ 0 on F₂, so they are concave in v. They are ≥ 0 on F₂ by Corollaries S and S′.
  - The same argument would lift to U(5) any such inequality that depends on one 3×3 block and holds on U(4).
- **Evidence for (H)** (`hull_support.py`): h₁(c) = max_{F₁} c·v must equal h₂(c) = max_{F₂} c·v for every direction c.
  - In 60 directions (six seeds, half with a negative det coefficient) the two maxima agree; the largest gap is in `logs/hull_support_*.txt`.
  - Local optimization can miss maxima, so this is evidence only.
- **Why (H) needs more than one disc.**
  - Along any admissible line, each coordinate of v is |α + zβ|², which is subharmonic. So v(W) lies below every average over a boundary circle, never on it.
  - Hence (H) needs the singular part of F₂, and points of F₂ with other singular vectors. For diagonal W, independent Bernoulli thinning (s_k ↦ 1 with probability s_k², else 0) gives v(W) exactly. For general W the same thinning gets the cross terms wrong, with s_k² in place of s_k.
- **Equality cases** (`contraction_min.py`, `interior_zeros.py`).
  - The zeros found have E(W) = 0 and Y = 0, and are real up to diagonal phases.
  - Some are interior, for example (E20) at singular values (1, 0.892, 0.854) with |det W| = 0.762. So any proof must be sharp at genuine blocks of U(5), not only on U(4).
  - On the zero set of Y₂₀, (E20) reduces to E(W) ≥ 0. Y₂₀ = 0 means f is a unit left singular vector of W: the unit singular pair (u, v) has u_a = 0 and v_γ = 0. A proof of E(W) ≥ 0 on that locus is a necessary first step.

### Routes that fail (`failed_routes.py`)

- **Corners in a singular value.** With the singular vectors fixed, the minimum over s₂ ∈ [0, 1] is often interior. So the reduction to s₂ ∈ {0, 1} fails line by line. The count is in `logs/failed_routes.txt`.
- **A rank-one imaginary part.** This is the device of Theorem F.
  - Diagonal phases with rank Im(D₁WD₂) ≤ 1 exist numerically for every tested contraction: five phase parameters against four equations.
  - Each minor of R + zN (R, N real, N of rank one) is a + zb with real a, b. So v(W) = cos²θ·v(R + tan θ·N) + sin²θ·v(R − cot θ·N) for every θ.
  - That reduces W to real matrices if both endpoints are contractions. For blocks of U(5) and of U(4) the best endpoints found exceed norm 1 by about 10⁻³. That is consistent with a boundary case, but not established.
  - The real case would remain open anyway.
- **Superharmonicity on all closure discs.** This fails: see the coverage above.

## Limits

- Nothing here proves (E20) or (E21) on U(5). The new proofs (Propositions B, S, M) are by the author, with exact checks of every identity, and are not independently reviewed.
- Corollary S rests on Corollary F2 of `unitary-factorization`, which rests on Theorems F and O; neither is reviewed.
- The support-function test, the coverage counts and the equality cases come from local optimization.

## Files

| File | Content |
| --- | --- |
| [`reduce_check.py`](reduce_check.py) | Propositions B and S: exact vertex checks, float comparison |
| [`contraction_min.py`](contraction_min.py) | minima of E, (E20), (E21) over blocks of U(4), U(5) and all contractions; equality data |
| [`interior_zeros.py`](interior_zeros.py) | zeros with 0.05 ≤ s₃ ≤ s₂ ≤ 0.95 |
| [`closure.py`](closure.py) | Proposition M: unitarity of the closure, the mean-value identity, coverage |
| [`hull_support.py`](hull_support.py) | support-function test of (H) |
| [`failed_routes.py`](failed_routes.py) | singular-value corners; rank-one imaginary part |
| `logs/` | outputs of the commands below |

## Reproduce

From this directory, with Python 3.9, numpy 2.0.2 and scipy 1.13.1:

```sh
python3 reduce_check.py                 > logs/reduce_check.txt       # seconds
python3 contraction_min.py 40 1         > logs/contraction_min.txt    # ~10 minutes
python3 interior_zeros.py 40 2          > logs/interior_zeros.txt     # ~3 minutes
python3 closure.py 40 3                 > logs/closure.txt            # ~5 minutes
python3 failed_routes.py 1000 30 4      > logs/failed_routes.txt      # ~10 minutes
for s in 11 12 13 14 15 16; do python3 hull_support.py $s 10 > logs/hull_support_$s.txt; done   # ~1 hour each, run in parallel
```

## Next steps

1. Prove E(W) ≥ 0 on the locus Y₂₀ = 0, where the unit singular pair has u_a = 0 and v_γ = 0. Every zero found lies there.
2. Prove or refute (H). A refutation needs one direction c with h₁(c) > h₂(c), certified.
3. Otherwise, try an exact sum-of-squares certificate for the t-forms on the 3×5 isometries [W | X], as `../../star-halfline/` did for (★) on U(4).
