# Rotation baselines and an explicit forest certificate

Task `mdo-baselines` · Agent `ganesh-reductions` · Owner `g4nesh` · November 25, 2026.

**Result:** for arbitrary complex diagonal spectra, the Marcus–de Oliveira inclusion holds when the relative unitary is a product of arbitrary two-coordinate unitary gates, on distinct edges of a forest, in any order. Explicit nonnegative permutation weights match **every squared minor simultaneously**. Disjoint two-coordinate blocks are the elementary special case. Left and right monomial unitary factors are allowed.

**Evidence level:** a complete written derivation for review, supported by fixed exact rational regression checks. This is not formal verification, has not been independently reviewed under another GitHub owner, and makes no novelty claim. It neither proves the general conjecture nor recovers all unitaries in dimension three. The [existing minors report](../minors/report.md), §3, already derives the full known `n ≤ 3` case by Birkhoff decomposition and complementary minors; that work is not repeated here. The single-Givens identity already appears in the [adversarial notes](../adversarial/README.md), §4. The present contribution gives its full product weights, a larger precisely delimited family, and an exact failure of extending these weights across a cycle.

The exact target and normal-pair reduction are in the [objective](../README.md) and [reductions](../reductions/README.md). No claims from the unreviewed Theorem D or Proposition E are needed below.

## 1. Conventions and elementary baseline

Write `D_a = diag(a₁,…,aₙ)`, `D_b = diag(b₁,…,bₙ)`, and

\[
d_{a,b}(U)=\det(D_a+UD_bU^*),\qquad
z_\sigma=\prod_{i=1}^n(a_i+b_{\sigma(i)}).
\]

Every spectral entry is an arbitrary complex number. No inverse of a spectral value or of a matrix sum will be used. For a permutation `p`, use the **column-to-row** convention `P_p e_j = e_{p(j)}`. Thus

\[
d_{a,b}(P_p)=z_{p^{-1}}=\prod_j(a_{p(j)}+b_j).
\tag{1}
\]

This inverse matters: the repository's `z_σ` convention sends a row to its assigned `b` index.

Let `G` equal the identity off coordinates `{i,j}`, and let its block `g` be any complex `2 × 2` unitary. Row and column norms imply

\[
|g_{11}|^2=|g_{22}|^2=t,\qquad
|g_{12}|^2=|g_{21}|^2=1-t,\qquad 0\le t\le1.
\]

The determinant of the corresponding block of `D_a+GD_bG*` is

\[
t(a_i+b_i)(a_j+b_j)
 +(1-t)(a_i+b_j)(a_j+b_i).
\tag{2}
\]

To see this directly, the diagonal entries of `g diag(b_i,b_j) g*` are `t b_i+(1−t)b_j` and `(1−t)b_i+t b_j`, while its determinant is `b_i b_j`. Expanding the two-by-two sum gives (2). These statements use complex unitarity, not real orthogonality; arbitrary entry phases are included. The full determinant is (2) times `∏_{k≠i,j}(a_k+b_k)`.

Now take pairwise disjoint pairs `e₁,…,eₘ`. The sum is block diagonal, so multiplication of (2) gives

\[
d_{a,b}(G_m\cdots G_1)
 =\sum_{\varepsilon\in\{0,1\}^m}
 \left(\prod_{r:\varepsilon_r=0}t_r
       \prod_{r:\varepsilon_r=1}(1-t_r)\right)
 z_{p_\varepsilon^{-1}},
\quad
p_\varepsilon=\tau_m^{\varepsilon_m}\cdots\tau_1^{\varepsilon_1},
\tag{3}
\]

where `τᵣ` exchanges the endpoints of `eᵣ`. The weights are nonnegative and sum to `∏ᵣ(tᵣ+1−tᵣ)=1`. In this disjoint case the transpositions commute and `p_ε` is its own inverse. Isolated coordinates simply contribute their one-dimensional factor. This proves the inclusion for this family, including `n=1` and every `n=2` relative unitary.

## 2. Forest-product theorem

**Theorem.** Let `e₁,…,eₘ` be distinct edges of an undirected forest on `{1,…,n}`. Let `Gᵣ` be the identity off `eᵣ`, with an arbitrary complex unitary two-by-two block `gᵣ` on that pair. Set

\[
V=G_m\cdots G_1,\quad q_r=|(g_r)_{12}|^2,\quad
w_\varepsilon=\prod_{r:\varepsilon_r=0}(1-q_r)
                     \prod_{r:\varepsilon_r=1}q_r,
\quad p_\varepsilon=\tau_m^{\varepsilon_m}\cdots\tau_1^{\varepsilon_1}.
\tag{4}
\]

For any equal-sized row and column sets `S,T`,

\[
|\det V[S,T]|^2
 =\sum_\varepsilon w_\varepsilon\,
   \mathbf1\{p_\varepsilon(T)=S\}.
\tag{5}
\]

Consequently, simultaneously for every `a,b ∈ ℂⁿ`,

\[
d_{a,b}(V)=\sum_\varepsilon w_\varepsilon z_{p_\varepsilon^{-1}}.
\tag{6}
\]

Order is arbitrary but fixed, and **each forest edge is used at most once**. The two-by-two gates may overlap in coordinates. A path or star therefore gives instances beyond disjoint blocks. Different bit strings that yield the same permutation may be combined by adding weights; it is not necessary to assert distinctness. Empty products have weight one.

### Proof of the squared-minor identity

Fix a subset size `k`. The `k`th compound matrix of `G` has entries `det G[S,T]` with subsets in increasing order. Cauchy–Binet says compounds multiply. Thus `det V[S,T]` is the sum over intermediate subsets

\[
X_0=T, X_1,\ldots,X_m=S
\]

of products `∏ᵣ det Gᵣ[Xᵣ,Xᵣ₋₁]`. For one gate on `{i,j}`, a nonzero entry can only preserve the subset, or exchange `i` and `j` when exactly one is occupied. The squared magnitudes of the possibilities are:

| Occupied endpoints before the gate | Subset transition | Squared magnitude |
| --- | --- | --- |
| Neither | Stay | 1 |
| Both | Stay | `|det gᵣ|² = 1` |
| Exactly one | Stay | `1−qᵣ` |
| Exactly one | Exchange endpoints | `qᵣ` |

Signs caused by ordering exterior basis vectors do not change these magnitudes. Zero entries simply remove a possible path.

**Uniqueness of the subset path.** Delete an edge `eᵣ` from its forest component and choose one of the two resulting vertex sets `C`. No gate except `Gᵣ` crosses from `C` to its complement: other edges cannot cross this cut without making a cycle, and `eᵣ` is not repeated. Hence a path from `T` to `S` must change the number of occupied vertices in `C` at step `r` by exactly

\[
\delta_r=|S\cap C|-|T\cap C|.
\]

If `δᵣ=1` or `−1`, that step must exchange the endpoints in the uniquely specified direction. If `δᵣ=0`, the subset must stay; a one-particle exchange would change the count. If `|δᵣ|>1`, no path exists. Starting from `T`, these prescribed moves determine every intermediate subset. They might be incompatible with the current occupation or with a zero gate entry, but there can never be two nonzero paths. The same argument applies separately to every component and includes isolated vertices.

There is consequently no interference between different products in the compound expansion: the squared modulus of the sum equals the sum of the products of the squared moduli, because at most one summand is nonzero.

Now independently choose the transposition `τᵣ` with probability `qᵣ` and the identity with probability `1−qᵣ`, in the same gate order. The induced subset process has exactly the transition table above. When both or neither endpoints are occupied, both choices give the same subset, and their probabilities sum to one. Its probability of ending at `S`, starting at `T`, is therefore both the left side of (5) and the right side of (5). This proves (5), including `k=0,n`.

### From all squared minors to the determinant

For completeness, expand the determinant by choosing diagonal entries of `D_a` on the complement of a set `S`, then apply Cauchy–Binet to the principal submatrix `(VD_bV*)[S,S]`:

\[
d_{a,b}(V)
 =\sum_{k=0}^n\sum_{|S|=|T|=k}
   |\det V[S,T]|^2
   \left(\prod_{i\notin S}a_i\right)
   \left(\prod_{j\in T}b_j\right).
\tag{7}
\]

The conjugate minor from `V*` has the same row/column ordering, so it multiplies the first minor to its squared modulus. This derivation needs no nonsingularity or distinctness assumption. Substitute (5) into (7) and interchange finite sums. For each `p_ε`, the inner sum is exactly (7) for the permutation matrix `P_{p_ε}`, hence equals (1). This proves (6). Nonnegativity and total mass one follow directly from (4). ∎

## 3. Phases, permutations and degeneracies

The theorem already allows every complex `U(2)` block, including determinant phases. It does not assume that a general complex unitary can be made real by diagonal phases.

It also applies to

\[
U=M_\ell V M_r
\]

for arbitrary monomial unitaries, with `M_ℓ e_j = phase_j e_{ℓ(j)}` and `M_r e_j = phase'_j e_{r(j)}`. Multiplying minors on either side permutes their index sets and multiplies by modulus-one scalars. Thus push the same weights forward to

\[
p'_\varepsilon=\ell\circ p_\varepsilon\circ r,
\qquad
d_{a,b}(U)=\sum_\varepsilon w_\varepsilon z_{(p'_\varepsilon)^{-1}}.
\tag{8}
\]

This also follows by permuting the two spectra and cancelling their commuting diagonal phases. The executable checks include nontrivial permutations and phases on both sides; (8) uses the actual ordered compositions, not a silent relabeling of the output.

All formulas are polynomial identities in the spectral entries. They remain valid with zeros, repeated eigenvalues, singular matrix sums, repeated permutation points, point or segment hulls, and gate probabilities zero or one. They do not divide by a determinant or remove boundary points. If a weight is zero, that term may simply be omitted.

## 4. A cycle breaks these weights, not the conjecture

Let `Rᵢⱼ` have real block `[[c,s],[-s,c]]` on coordinates `{i,j}`, with `c=3/5`, `s=4/5`, and put

\[
U=R_{13}R_{23}R_{12}.
\]

There are now two paths from coordinate `1` back to `1`: staying at the first and last gates, or transferring around the triangle. Direct multiplication gives

\[
U_{11}=c^2+s^3=109/125,\qquad
|U_{11}|^2=11881/15625.
\]

The independent-transposition construction instead gives return probability

\[
c^4+s^6=6121/15625.
\]

The difference is the nonzero interference term `2c²s³=1152/3125`. With `a=(0,1,1)` and `b=(1,0,0)`, (7) isolates this coefficient, so the actual determinant is `11881/15625` while the proposed weights give `6121/15625`. The MOC permutation points are just `0` and `1`, and the actual value lies in their segment. This is therefore **only a counterexample to extending this particular product-weight formula**. Other convex weights do exist for this case.

Repeating an edge also destroys the proof's unique-crossing argument, even if the set of distinct edges is a tree. Some repeated or cyclic products may nevertheless have a valid certificate, possibly after simplification or with different weights. Neither a cycle nor a repeated edge is a necessary obstruction to MOC or to all certificates.

## 5. Reproduction and bounded evidence

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/baselines/forest.py --output research/marcus-de-oliveira/baselines/results.json
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/marcus-de-oliveira/baselines -p 'test_*.py' -v
```

The suite uses Python's `fractions.Fraction` and the adjacent [exact verifier](../certificates/verify.py)'s Gaussian-rational arithmetic. It uses no NumPy, floating-point tolerance, optimization, model invocation or random sampling. Its fixed fixtures define every spectrum, gate and permutation in source; the complete resulting unitaries and all nonzero permutation weights are stored in [results.json](results.json). Gaussian-rational entries are `[real rational string, imaginary rational string]`; permutations use zero-based column-to-row indices. These are exact constructions, never rounded numerical unitaries. The recorded Python version and dependency hashes pin the computation.

The first run used Python 3.9.6 and **0.861063 process CPU seconds**, checking 11 forest fixtures in dimensions 1–6, with at most five gates, and all **3,332 compound entries** across every subset size. Both unitarity products, nonnegative weights of total one, every compound identity, and the determinant identity passed. Determinants were also checked by a separate Leibniz expansion sharing only scalar arithmetic. One additional triangle fixture exactly reproduced the obstruction in §4. The command has a 60-second POSIX CPU limit. All **11 unittest checks passed in 0.807 seconds**; these replay the same fixed cases and additionally reject repeated edges, cycles, invalid endpoints, nonunitary gates and inexact scalar input. No new sampled inputs were generated by the replay.

The saved run is an author-side check of finitely many instances. The general statement rests on the written proof, which still needs independent mathematical review. A same-owner agent read the cut/path lemma and found no gap; that is not an independent-owner ledger review. No runtime token counter is exposed to this child agent. Token usage remains unknown here and must not be estimated or counted again if included in a parent receipt.

## 6. What remains and suggested review

The forest restriction is substantive: this construction uses at most `n−1` distinct coordinate-pair gates, whereas arbitrary unitary decompositions generally involve cycles or repeated edges. It supplies a concrete certificate family useful for testing an exact checker and separating phase/index errors from the harder global problem. It does not close the original objective, give a general support bound, validate Zheng's dimension-four claim, or supply a new counterexample.

An independent reviewer should check the no-interference cut argument for arbitrary gate order, the compound transition table in the double-occupation sector, the permutation inverse in (1)/(6)/(8), and the stated degenerate cases. A next research question is whether controlled interference for a single cycle can be absorbed into **different** nonnegative permutation weights. The exact triangle computation is a regression target for that question, not evidence that a general cycle extension is true.
