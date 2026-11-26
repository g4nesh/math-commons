# Actual diagonal-update determinant expansion

This is the **proved principal-minor stage** of the matrix side of the forest
spectral bridge. It does not complete the squared-minor expansion of
$U D_b U^*$, and it does not prove the general Marcus–de Oliveira conjecture.

`formal/ForestMatrixSumExpansion.lean` derives, from the existing recursive
signed determinant, the identity

$$
\det(D_a+M)=\sum_{K\subseteq V}\left(\prod_{i\notin K}a_i\right)\det M[K,K].
$$

The matrix $M$ and spectra $a$ are arbitrary **Gaussian rationals**
(`Rat × Rat`). The finite vertex list is duplicate-free. Neither invertibility,
normality, positivity, nonzero entries nor distinct spectral values are assumed.
Empty and full subsets are included. The empty determinant and empty product
are one. No premise is the expansion being proved.

## Exact statements and concrete definitions

The public declaration inventory contains **21 proved declarations**. Its exact
Lean types are in [statements.txt](statements.txt); source hashes and axiom
reports are in [audit.json](audit.json).

| Declaration | Proved content |
| --- | --- |
| `diagonalRowMinor` | Replacing one row by a diagonal delta gives its positive principal cofactor, at any list position. |
| `singleDiagonalUpdate` | Actual determinant is affine in an arbitrary diagonal entry, with the smaller principal minor as coefficient. |
| `partialDiagonalExpansion` | Every distinct finite collection of diagonal updates expands into principal minors; update coordinates must lie in the duplicate-free row list. |
| `diagonalSumExpansion` | Full diagonal-plus-arbitrary-matrix expansion, indexed by the coordinates chosen from the diagonal. |
| `retainedDiagonalExpansion` | The same expansion indexed by the coordinates retained in $M$, exactly as displayed above. |
| `spectralPrincipalExpansion` | Applies the proved formula to the **actual matrix** $D_a+U D_b U^*$, leaving the actual principal minors of $U D_b U^*$ as coefficients. |
| `conjugatedDiagonalEntry` | Derives $(U D_b U^*)_{ij}=\sum_k U_{ik}b_k\overline{U_{jk}}$ from actual finite matrix multiplication. |
| `diagonalDeterminant` | The actual determinant of a diagonal matrix equals its product of entries. |
| `emptySpectralDeterminant`, `singletonSpectralDeterminant` | Direct order-zero and order-one checks, including the coefficient $|U_{xx}|^2b_x$ in order one. |

All names are in namespace `MathCommons.ForestMatrixSumExpansion`.

`matrixAdd`, `diagonal`, and `adjoint` are entrywise addition, the diagonal
Kronecker-delta matrix, and conjugate transpose. `conjugatedDiagonal vertices U b`
is literally

```lean
matrixMul vertices (matrixMul vertices U (diagonal b)) (adjoint U)
```

Here `matrixMul` is the prior explicit finite sum, and every determinant is
`ForestGateSupport.minorDet gaussianArithmetic M rows cols`. The target
determinant is **not defined as a subset sum**. In a finite-dimensional
interpretation, `vertices` enumerates the coordinate space. The principal-minor
identities themselves do not require that this list exhaust the ambient Lean
type: only the listed coordinates and finite matrix-product sums are used.

The retained-coordinate theorem has this concrete shape:

```lean
minorDet gaussianArithmetic (matrixAdd (diagonal a) M) rows rows =
  scalarSum gaussianArithmetic ((subsets rows).map fun s =>
    gaussianArithmetic.mul
      (product rows (fun x => if s x then gaussianArithmetic.one else a x))
      (minorDet gaussianArithmetic M (rows.filter s) (rows.filter s)))
```

Its sole mathematical hypothesis is `rows.Nodup`. `subsets` is the existing
proved complete finite enumeration of supported Boolean subsets, with no
repeated functions when the coordinate list is duplicate-free. A local
complement map is false outside the coordinate list; `complementsPerm` proves
that complement permutes the enumeration without introducing multiplicities.

## Derivation

1. Move a distinguished delta row to the first row using the already proved
   actual row-movement parity. Expand in its sole supported column using the
   actual Laplace recursion. Both signs are the parity of the same prefix and
   cancel. The two remaining lists keep their original common order.
2. Apply the previously proved row multilinearity to add one diagonal entry.
   This derives the single-update formula, including zero updates.
3. Induct on an arbitrary duplicate-free update list. The branch using the
   original entry keeps the full row list; the branch using the diagonal entry
   removes that coordinate. The remaining update coordinates lie in both row
   lists. The induction proves the complete principal-minor sum, with every
   coefficient a product of the chosen diagonal entries.
4. Prove that supported complementation permutes the canonical subset list.
   Reindexing gives the retained-coordinate convention used by the existing
   forest squared-minor law.
5. Substitute the concrete conjugated diagonal matrix. No property of $U$ is
   needed for this principal-minor stage.

The proof uses no divisions, spectral genericity, numerical search, theorem
placeholders, `sorry`, `native_decide`, or newly declared axioms. It is an
implementation of the standard diagonal-update identity, not a claim of a new
mathematical theorem.

## Exact remaining bridge

For each retained subset $K$, the next obligation is

$$
\det (U D_b U^*)[K,K]
=\sum_{J\subseteq V}\left|\det U[K,J]\right|^2\prod_{j\in J}b_j.
$$

This requires a derived arbitrary-order contraction/Cauchy–Binet identity and
the appropriate adjoint-minor identity, or a genuinely derived restriction to
forest gate products. The actual recursive minor is already zero for unequal
selected cardinalities, but that fact alone is not the missing contraction.
This file neither assumes nor claims the displayed remaining identity.

Once that step is proved, substitute it into `spectralPrincipalExpansion`,
interchange the finite sums, and compare with
`ForestSpectralExpansion.forestWeightedSubsetExpansion`. That prior scalar
result already supplies a single nonnegative normalized permutation law for
actual forest gate products and all subset pairs. The convention is
**minor rows = finish $K$, columns = start $J$**; the corresponding pairing
product uses the column-to-row permutation (or its inverse in row indexing).

No arbitrary-complex transfer is included. Unlike the concrete forest API,
these principal-minor results need no `Edge`, so zero- and one-dimensional
coordinate spaces are directly available.

## Reproduce

Use Lean **4.24.0** and Python's standard library:

```sh
python3 -I -B research/marcus-de-oliveira/forest-matrix-sum-expansion/check.py /path/to/lean
```

The checker freshly builds all **14 transitive local modules** in a temporary
directory, treats warnings as errors, checks **324 transitive axiom reports**,
compares the exact types of the **21 new declarations**, and verifies all
source SHA-256 hashes. It does not use cached `.olean` files, Mathlib or the
network. Only `propext`, `Classical.choice`, and `Quot.sound` are allowed in the
audit. Re-running `python3 formal/check.py` also exercises the separately
registered module and this exact-type audit after registry integration.

Kernel checking verifies the stated rational identities. Independent review
of their intended mathematical interpretation remains separate from source
merging and CI. The artifact makes no independent-owner review claim.
