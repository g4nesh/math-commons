# Actual conjugated-minor contraction through order two

This artifact proves the **order-zero, order-one and order-two matrix-side
contraction**, in arbitrary ambient dimension, for arbitrary Gaussian-rational
matrices $U$ and spectra $b$. It also proves entrywise conjugation and diagonal
scaling identities for actual minors of **every order**. It does not prove
higher-order Cauchy–Binet, the arbitrary-order adjoint transpose identity, or a
full formal Marcus–de Oliveira theorem.

Write $C=U D_b U^*$ using the literal finite matrix operations of
[ForestMatrixSumExpansion](../forest-matrix-sum-expansion/README.md). For any
ordered pair of retained indices $r,s$, the checked principal identity is

$$
\det C[\{r,s\},\{r,s\}]
=\sum_{j<k}\left|\det U[(r,s),(j,k)]\right|^2 b_jb_k.
$$

The notation $j<k$ means **earlier in the supplied coordinate list**, not an
order on vertex labels. Each pair of distinct intermediate coordinates occurs
once. The literal Lean statement retains the ordered lists `[r,s]` and
`[p.1,p.2]`; no sorting sign is silently removed. Repeated queried indices are
allowed in the generic order-two identities and give the corresponding zero
minor. The finite multiplication list must be duplicate-free.

The entries of $b$ may be complex Gaussian rationals, zero, or repeated.
Neither unitarity, a forest schedule, nonsingularity, nor a determinant
expansion premise is needed for this generic contraction.

## Main checked declarations

All declarations are in namespace `MathCommons.ForestAdjointContraction`.
The exact inventory of **20** declarations and all quantified types are in
[declarations.json](declarations.json) and [statements.txt](statements.txt).

| Declaration | Scope |
| --- | --- |
| `conjugateMinor` | Entrywise conjugation commutes with the actual minor in every order, including rectangular lists. This is not a transpose identity. |
| `columnScaleMinor` | Multiplying each column by its scalar multiplies any actual minor by the product of the selected column scalars; rectangular lists and zero factors are allowed. |
| `minorCongrEntries` | Equality only on the selected rectangle suffices for equality of the actual minors. |
| `rightDiagonalMinor` | Actual multiplication by $D_b$ has the corresponding minor-scaling effect, when selected columns lie in the duplicate-free multiplication universe. |
| `adjointMinorZero`, `adjointMinorOne`, `adjointMinorTwo` | Actual adjoint minors are conjugates of the reversed-index minors, in these three orders, with both list orders explicit. |
| `conjugatedMinorZero`, `conjugatedMinorOne`, `conjugatedMinorTwo` | Actual principal-minor contraction in orders zero, one and two. |
| `conjugatedMinorTwoMixed` | The off-principal identity $\det C[(r,s),(c,d)]=\sum_{j<k}b_jb_k\det U[(r,s),(j,k)]\overline{\det U[(c,d),(j,k)]}$. |
| `conjugatedTwoSubset` | The order-two principal identity on the actual filtered list of any retained Boolean subset of size two. |
| `twoSpectralPencil` | A closed expansion of the actual two-coordinate principal determinant of $D_a+U D_b U^*$ into squared actual $U$ minors. |
| `forestTwoContraction` | The previously constructed **same full-permutation law**, nonnegative and normalized, supplies the coefficients of every actual retained two-subset minor for a forest unitary gate product. |

The complete two-coordinate spectral formula is

$$
\begin{aligned}
\det(D_a+C)[(r,s),(r,s)]
={}&a_ra_s
+a_r\sum_j|U_{sj}|^2b_j
+a_s\sum_j|U_{rj}|^2b_j\\
&+\sum_{j<k}|\det U[(r,s),(j,k)]|^2b_jb_k.
\end{aligned}
$$

`twoSpectralPencil` assumes $r\ne s$ for the diagonal-matrix convention. This
is a two-coordinate principal determinant in an arbitrary ambient system. It
is **not** the full ambient determinant when the ambient dimension exceeds
two. The empty determinant is checked directly without requiring a gate or
nonempty coordinate type.

## Why the identities are derived

All determinants are the existing signed recursive
`minorDet gaussianArithmetic M rows cols`. `conjugatedDiagonal` is the
ordinary finite matrix product `(U * D_b) * adjoint U`, not a formal polynomial
chosen to have the desired coefficients.

1. Induct on the actual Laplace recursion to prove entrywise conjugation and
   column scaling. The `choicePerm` certificate identifies the product of all
   column factors with the chosen factor times those of the remaining columns.
   No division by a column factor occurs.
2. Prove that the actual matrix product $U D_b$ equals column scaling on the
   selected rectangle. `minorCongrEntries` transports this entry equality to
   an arbitrary-order minor.
3. Derive the actual adjoint identities directly in orders zero, one and two.
   The order-two proof uses the existing explicit two-by-two determinant and
   Gaussian-rational arithmetic. No arbitrary-order transpose theorem is used.
4. Apply the previously **proved** `ForestProductMinor.cauchyBinetTwo` to
   $(U D_b)U^*$. Substitute the derived scaling and adjoint identities. For a
   principal minor, $z\overline z$ becomes its rational squared norm. The mixed
   identity keeps both distinct minor factors.
5. Combine these coefficients with the literal two-by-two diagonal-addition
   identity to obtain `twoSpectralPencil`.
6. For the forest specialization, apply
   `ForestPermutationLaw.commonSquaredMinorLaw`. The proved `pairFilter`
   identifies each canonical intermediate pair with its Boolean subset, so the
   same law can be substituted into the actual order-two contraction.

`forestTwoContraction` puts nonnegativity and normalization before the universal
quantifier over retained size-two subsets. Its law depends only on the gate
schedule, its off-diagonal entries and the prefix length, never on $b$ or the
chosen subset. The actual minor convention remains **rows = finish, columns =
start**. Its forest assumptions are precisely unitary used Gaussian-rational
blocks and `EveryEdgeIsBridge`. The inherited `Nat → Edge vertices` API still
requires two ambient endpoints even when the used prefix is empty; the generic
contraction results have no such restriction.

## Remaining obligation

For retained sets of size greater than two, the general contraction

$$
\det(U D_b U^*)[K,K]
=\sum_{J\subseteq V}|\det U[K,J]|^2\prod_{j\in J}b_j
$$

is still **unproved in this formal library**. Arbitrary-order entrywise
conjugation alone does not justify replacing an adjoint minor by a conjugate
minor with its index lists reversed. The general transpose identity and
ordered-tuple grouping into distinct intermediate subsets still need derived
proofs, or a clearly delimited forest-specific replacement.

Once that gap is closed, substitute it into the already proved
`ForestMatrixSumExpansion.spectralPrincipalExpansion` and compare with
`ForestSpectralExpansion.forestWeightedSubsetExpansion`. This artifact closes
the contraction only through order two. It does not claim that the resulting
low-order formula alone proves the full forest determinant identity in every
ambient dimension. Transfer from Gaussian rationals to arbitrary complex
matrices remains separate.

These are formalizations of standard algebraic identities, not novelty claims.
No `sorry`, assumed target identity, `native_decide`, or new axiom is used.

## Reproduce

Use Python's standard library and Lean **4.24.0**:

```sh
python3 -I -B research/marcus-de-oliveira/forest-adjoint-contraction/check.py /path/to/lean
```

The checker freshly builds **15 local dependency modules** into a temporary
directory, treats warnings as errors, checks **344 transitive axiom reports**,
compares the exact types of all **20 new declarations**, and verifies every
source SHA-256 hash against [audit.json](audit.json). No cached `.olean` files,
network or Mathlib are used. Only `propext`, `Classical.choice`, and `Quot.sound`
are allowed. After the separate registry integration, `python3 formal/check.py`
includes the module and its exact audit alongside all previous checks.

The kernel proves the displayed statements over the specified scalar domain.
Source merging, CI, same-owner semantic feedback and independent-owner
mathematical review remain distinct evidence categories.
