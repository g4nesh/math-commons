# Actual product minors: tuple expansion and order-two Cauchy–Binet

[`ForestProductMinor.lean`](../../../formal/ForestProductMinor.lean) defines actual finite matrix multiplication over Gaussian rationals. It proves a general expansion of product minors over ordered intermediate tuples, the full Cauchy–Binet identity for **2×2 minors in any ambient dimension**, and its iteration through any finite matrix schedule. It also proves sparse gate row/minor updates and connects an actual one-gate matrix product to the existing squared transition probabilities.

**The full forest product-minor/probability bridge remains incomplete.** In particular, the iterated pair kernel below is not yet identified with the canonical Boolean-subset histories in `ForestEndpointMass`. Arbitrary-order determinant alternation and Cauchy–Binet remain separate obligations. No MDO resolution is claimed.

## Actual multiplication and ordering

For the supplied finite enumeration `vertices`, `matrixMul vertices M N` is

$$ (MN)_{ij}=\sum_{k\in\mathrm{vertices}}M_{ik}N_{kj}. $$

`matrixProduct` starts at the Kronecker identity and uses $P_{n+1}=G_nP_n$. Thus later gates act on the left. These definitions use actual matrix entries and scalar multiplication; neither is defined by the desired minor or path formula. `matrixMulIdentity`, `matrixIdentityMul` and `gateProductOne` verify the identity and one-step conventions under the stated exhaustive, duplicate-free enumeration hypotheses.

The minor is the existing signed row expansion `ForestGateSupport.minorDet`. All statements use its explicit ordered row and column lists and Gaussian-rational arithmetic. The translation to a standard complex matrix library remains open.

## Unconditional general tuple expansion

`productMinorTupleExpansion` proves

$$
\det((MN)[R,C])=
\sum_{K\in\mathrm{vertices}^{|R|}}
\left(\prod_{i=1}^{|R|}M_{R_i,K_i}\right)\det(N[K,C]).
$$

`rowExpansion` explicitly enumerates every ordered tuple $K$, including repeated indices, and its entry-product coefficient. The statement is valid for arbitrary row and column lists, including rectangular minors. It needs no forest, unitarity, nonvanishing or Cauchy–Binet premise. The proof checks finite-sum interchange and determinant linearity directly. `rowDet` makes the row vectors explicit, and `minorAsRows` proves its equality with the existing ordered minor.

This is the expansion before grouping tuples by their unordered index sets. Repeated-index cancellation and the determinant permutation signs are not assumed in this theorem.

## Complete order-two identity and iteration

`orderedPairs vertices` lists every pair of enumeration positions in increasing order. For duplicate-free `vertices`, it lists each unordered pair of distinct intermediate vertices exactly once, in the enumeration's order. `cauchyBinetTwo` proves

$$
\det((MN)[r,s\mid c,d])=
\sum_{(k,l)\in\mathrm{orderedPairs(vertices)}}
\det(M[r,s\mid k,l])\det(N[k,l\mid c,d]).
$$

The proof uses the general tuple expansion, cancels repeated indices by the explicit two-by-two determinant formula, and combines $(k,l)$ with $(l,k)$. The intermediate-pair grouping is proved by list induction in `gaussianSumPairs`. The result does not restrict the ambient matrix dimension to two and has no forest or unitary hypothesis. The identity continues to hold with enumeration multiplicities; the standard finite-matrix interpretation uses the existing exhaustive, duplicate-free enumeration.

`twoMinorKernel M` is the matrix of these ordered two-by-two minors. `twoMinorEvolution` starts from the **actual identity-matrix minor kernel**, then multiplies successive gate kernels over `orderedPairs vertices`. The initial kernel retains the proper signs even when the queried row or column pair is reversed or repeated. `iteratedCauchyBinetTwo` proves that this evolution equals the corresponding minor of the actual product $P_n$, for every finite $n$. `gateProductTwoMinorRecurrence` specializes the recurrence to the existing two-coordinate gate schedule. `cauchyBinetOne` and `productEmptyMinor` supply orders one and zero.

## Sparse gate results and exact remaining premises

`gateMulLeft`, `gateMulRight` and `gateMulSpectator` prove the actual two-row update under multiplication by a gate. `minorRowAdd`, `minorRowScale` and `minorRowLinearCombination` hold at any occupied position in a duplicate-free row list, including arbitrary spectator prefixes.

The checked gate-minor declarations distinguish three levels:

* `gateMinorNeither`, `gateMinorLeftOnly` and `gateMinorRightOnly` are unconditional sparse updates. A replacement row stays in its original position. These statements do not silently reorder it into a filtered subset.
* `gateMinorBothExpansion` is an unconditional four-term expansion, retaining the two duplicate-row terms and the exchanged-row term.
* `gateMinorBothOfAlternation` explicitly assumes those two duplicate-row minors are zero and that exchanging the rows negates the minor. It then gives the block-determinant factor. Those three generic alternation facts are **premises**, not newly proved theorems.

`oneGateProductMagnitude` uses the proved identity between the actual length-one matrix product and its gate, then applies Raghav's existing `stepMagnitudeGeneral`. It covers every input/output subset pair, all occupation sectors, zero probabilities and either endpoint order under `UnitaryBlock`.

The next obligations are precise: connect ordered pairs to the canonical size-two Boolean-subset histories; use gate support to remove illegal transitions; apply the existing forest endpoint no-interference theorem. For arbitrary minor order, additionally prove duplicate-row cancellation and row-permutation signs, then group the general tuple expansion into subsets. The resulting standard-complex and common-permutation-mixture/MDO interfaces remain outside this module.

## Verification and provenance

The [declaration manifest](declarations.json) lists all 52 checked results. Lean **4.24.0** compiles the module with `-DwarningAsError=true`; every printed transitive axiom report is restricted to `propext`, `Classical.choice` and `Quot.sound`. No `sorry`, `native_decide` or new axiom is used. The repository audit rebuilds all imports in a fresh temporary directory and includes this module after `ForestEndpointMass`:

```sh
python3 formal/check.py
```

The command accepts an explicit pinned Lean executable as its first argument. The adjacent [audit report](audit.json) records the source hashes, exact declaration set and observed standard-axiom audit. Source publication and CI do not constitute independent-owner mathematical review.

This work extends the merged `ForestGateSupport`, `ForestGateMagnitude`, Raghav's `ForestPathProbability` (PR29), and the normalization/endpoint modules (PR36). No dependency source is modified. The main task is `mdo-forest-product-minor`; checker integration is separately scoped under `mdo-forest-product-ci`.
