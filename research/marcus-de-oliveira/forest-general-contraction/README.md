# The actual rational forest determinant law

**The matrix-side bridge is now proved at every subset order.**
`ForestGeneralContraction.forestCommonDeterminantLaw` proves that one explicit
nonnegative normalized law on genuine permutations represents the **actual**
determinant of $D_a+U D_b U^*$ for every pair of Gaussian-rational spectra,
when $U$ is the existing finite forest unitary gate product.

Let $U=G_{m-1}\cdots G_0$. Each gate acts on its specified pair of coordinates
with block $(\alpha_k,\beta_k;\gamma_k,\delta_k)$. Put $q_k=|\beta_k|^2$.
The finite law independently chooses the identity or the corresponding
transposition with factors $1-q_k$ and $q_k$, respectively. It composes the
newest transposition on the left. Different choices producing the same
permutation can remain separate entries. For each entry $(p,w_p)$,

$$
w_p\ge0,\qquad \sum_p w_p=1,\qquad
\det(D_a+U D_b U^*)=\sum_p w_p\prod_j(a_{p(j)}+b_j).
$$

The law depends only on the gate data and prefix length, **not on the spectra**.
The final declaration places `∀ a b` after nonnegativity and normalization in
its conclusion. The forward permutation maps columns to rows. The equivalent
row-indexed product is $\prod_i(a_i+b_{p^{-1}(i)})$, as already proved in
`ForestSpectralExpansion.pairingProductAsRows`.

This is a kernel-checked result over `Gaussian = Rat × Rat`, with the
repository's explicit finite matrix product and signed Laplace determinant.
It is not an unrestricted MDO proof or a Lean theorem over arbitrary complex
matrices. The previously written [complex forest theorem](../baselines/README.md#2-forest-product-theorem)
already supplies the broader mathematical statement; this artifact closes the
formal rational matrix bridge rather than claiming a new complex theorem.

## Exact declarations and assumptions

All **16** public declarations are in namespace
`MathCommons.ForestGeneralContraction`. Their exact Lean types are pinned in
[statements.txt](statements.txt), with the inventory in
[declarations.json](declarations.json).

| Declaration | Proved scope |
| --- | --- |
| `firstColumnExpansion` | Actual first-column Laplace expansion, derived from the row-recursive determinant. Arbitrary rectangular lists and repeated indices are allowed. |
| `transposeMinor` | Transposition exchanges the supplied row and column lists in every actual minor; no square or uniqueness premise. |
| `adjointMinor` | Every actual adjoint minor equals the conjugate of the reversed-index minor. |
| `matrixMulAssoc` | Associativity of the actual finite matrix sum, proved by finite-sum interchange. |
| `gateProductContraction` | Arbitrary-order Cauchy–Binet for an actual **left product of two-coordinate gates** against an unrestricted right matrix. Neither forest nor unitarity is needed. |
| `leftDiagonalMinor` | Actual left multiplication by a diagonal matrix scales a selected minor by the product of its row factors. |
| `gateConjugatedMinor` | Every principal minor of $U D_b U^*$ has the squared-$U$-minor spectral expansion for an actual gate product $U$, with no forest or unitarity premise. |
| `gateSpectralDeterminant` | The full actual determinant of $D_a+U D_b U^*$ equals the canonical double subset sum in the exact start/finish convention. |
| `forestDeterminantLaw` | The actual determinant equals the nonnegative normalized weighted pairing-product sum for supplied spectra. |
| `forestCommonDeterminantLaw` | The same explicit normalized law works simultaneously for **all** spectra. |

The last theorem assumes exactly:

- a finite coordinate list provided by `edges : Nat → Edge vertices`; each
  `Edge` certifies the same complete duplicate-free enumeration and two
  distinct endpoints;
- a finite used prefix of length `m`;
- `EveryEdgeIsBridge edges m`, the existing explicit indexed-edge reachability
  condition; and
- `UnitaryBlock` for each used two-by-two Gaussian-rational block.

These are the previously formalized forest and local-unitarity hypotheses.
There is no matrix-to-history equality, transpose identity, Cauchy–Binet
identity, determinant spectral expansion or convex-weight representation
supplied as an extra premise. Spectra may contain zeros or repetitions, and
the actual matrix sum may be singular.

The inherited `Nat → Edge vertices` interface requires at least two ambient
coordinates even for a zero-length schedule. Empty and full **selected
subsets** and unequal selected cardinalities are included. The transpose,
adjoint and scalar matrix identities do not have that `Edge` restriction.
A zero-length gate product is the actual identity matrix.

## Derivation of the missing matrix side

### 1. Transpose and adjoint

Expand a minor with first column $c$ using the existing first-row recursion.
Its first term uses the top-left entry. In the remaining terms, apply the
induction hypothesis to the first column of the smaller determinant.
Interchange the finite row-choice and column-choice sums. `signedCross`
checks the two sign orders, including both parity values, and multiplication
commutativity identifies the summands. This derives `firstColumnExpansion`.

Now induct on the row list of the transpose. Its defining row expansion is
the just-proved column expansion of the original matrix, with the smaller
transpose identities supplied by induction. Combine `transposeMinor` with the
previous arbitrary-order entrywise `conjugateMinor` to obtain `adjointMinor`.
Rectangular and duplicate-index lists require no division or exceptional case.

### 2. Contraction using sparse gates

The proof does **not** assume generic Cauchy–Binet for arbitrary matrix factors.
For a zero-length left product, the actual identity-minor delta and the
complete duplicate-free Boolean subset enumeration give the contraction.
For one more gate, actual matrix associativity and the already proved
`gateMinorPredecessors` reduce the output minor to its one or two predecessor
subsets. Apply the induction hypothesis, interchange the two finite sums,
and recombine the same sparse recurrence on the left product minor.

This proves, for the actual gate product $U$, every retained subset $K$, every
right matrix $N$, and every supplied column list $C$,

$$
\det(UN)[K,C]=\sum_{J\subseteq V}\det U[K,J]\det N[J,C].
$$

The statement covers all minor orders and uses no forest or unitary property.
It is explicitly restricted to the left gate-product factor; it does not
claim every arbitrary matrix has been expressed as such a product.

### 3. Conjugated diagonal and determinant

Reassociate $(U D_b)U^*=U(D_bU^*)$, apply the derived contraction, and use
`leftDiagonalMinor` and `adjointMinor`. The two matching minor factors become
the squared norm. Thus

$$
\det(U D_b U^*)[K,K]
=\sum_{J\subseteq V}|\det U[K,J]|^2\prod_{j\in J}b_j.
$$

Substitute into the previously proved actual principal-minor expansion of
$D_a+M$, then interchange the finite subset sums. `gateSpectralDeterminant`
now matches `ForestSpectralExpansion.forestWeightedSubsetExpansion` exactly:
**minor rows = finish $K$, columns = start $J$**. Their equality proves the
actual forest determinant law and its common-for-all-spectra version.

## Evidence and remaining boundary

Reproduce with Python's standard library and Lean **4.24.0**:

```sh
python3 -I -B research/marcus-de-oliveira/forest-general-contraction/check.py /path/to/lean
```

The checker freshly builds **16 local modules** into a temporary directory,
treats warnings as errors, audits **360 transitive declaration reports**,
compares all **16 new exact theorem types**, and checks every dependency source
SHA-256 against [audit.json](audit.json). No cached `.olean`, Mathlib, network,
`sorry`, `native_decide`, or new axiom is used. Only `propext`,
`Classical.choice`, and `Quot.sound` are allowed.

The remaining generic arbitrary-matrix Cauchy–Binet theorem is **not needed for
this proved forest conclusion**. A standard-complex matrix-library interface
and a Lean-checked continuity/density transfer are still outside this artifact.
The mathematical complex forest result is already written in
[the baseline proof](../baselines/README.md), and the rational-unitary density
argument is already written in [the reductions](../reductions/README.md#4-exact-rational-data-suffice-for-finding-a-counterexample).
Those written results are separate evidence; they are not silently included in
the Lean audit. No new density theorem or unrestricted MDO claim is made here.

The source builds on Raghav's path-probability module and the existing support,
magnitude, path normalization, endpoint uniqueness, actual product-minor,
permutation-law and spectral modules. Source merging and same-owner semantic
feedback do not replace independent-owner mathematical review.
