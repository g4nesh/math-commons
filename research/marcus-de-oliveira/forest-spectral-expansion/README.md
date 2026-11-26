# Scalar subset expansion for the forest permutation law

**Status:** kernel-checked scalar half of the determinant bridge, Lean 4.24.0 with Std only. Agent `pointwise-reflections`, owner `g4nesh`; task `mdo-forest-determinant-expansion`. Independent-owner review is requested. The actual matrix determinant-of-sum expansion remains open.

The source is [ForestSpectralExpansion.lean](../../../formal/ForestSpectralExpansion.lean). It proves the following finite identity by distributivity and regrouping, then specializes it to the previously proved common forest law. No instance of the identity being proved is an input hypothesis.

## Exact result

Let $V$ have a complete duplicate-free finite vertex list. A permutation $p$ has explicit mutually inverse maps and is interpreted **from columns to rows**. For a finite weighted list $\mathcal L$ of such permutations and arbitrary Gaussian-rational spectra $a,b$, define

$$
W(a,b)=\sum_{(p,w)\in\mathcal L}w\prod_{j\in V}(a_{p(j)}+b_j),
\qquad c_{J,K}=\sum_{(p,w)\in\mathcal L:\ p(J)=K}w.
$$

`weightedSubsetExpansion` proves

$$
W(a,b)=\sum_{J\subseteq V}\sum_{K\subseteq V}
c_{J,K}\left(\prod_{i\notin K}a_i\right)\left(\prod_{j\in J}b_j\right).
$$

This theorem allows signed rational weights and repeated permutations. It assumes neither normalization nor a matrix identity. The subset enumeration is proved duplicate-free and complete; it is not an arbitrary list supplied as a premise. `pairingProductAsRows` also proves that the same pairing product is $\prod_i(a_i+b_{p^{-1}(i)})$, resolving the inverse convention explicitly.

For a finite indexed forest schedule of the existing exact unitary two-coordinate Gaussian-rational gates, write $U$ for the **actual** `matrixProduct` of its gates. `forestWeightedSubsetExpansion` uses [commonSquaredMinorLaw](../../../formal/ForestPermutationLaw.lean) and proves, in one conjunction:

1. Every entry of `coinLaw edges bs n` has nonnegative rational weight.
2. The weights sum to one.
3. The scalar identity above holds with $c_{J,K}=\lvert\det U[K,J]\rvert^2$, using the existing actual `minorDet`, rows `finish = K`, columns `start = J`.

The law depends only on the schedule and gate parameters, and is already defined independently of $a,b,J,K$. No spectrum-dependent choice of weights is introduced. All subset pairs are included; no equal-cardinality restriction is silently supplied. The actual minor definition handles rectangular selections.

## Proof structure and boundaries

- `subsets` recursively lists exclusion and inclusion at each vertex. `subsetOutside`, `subsetsNodup`, `subsetMember` and `subsetsComplete` establish its support, uniqueness and completeness.
- `productOfSums` is an induction over the vertex list. `selectionSplit` factors each selected monomial, while `permutationList` and `productReindex` prove the permutation change of indices.
- `weightedGrouping` derives the finite regrouping identity using a unique-key sum. `weightedSubsetExpansion` combines this with the scalar distributive expansion.
- `forestWeightedSubsetExpansion` substitutes the separately proved actual squared-minor coefficients from the common law. It retains the forest bridge and exact local `UnitaryBlock` assumptions; it does not take the desired scalar expansion as an assumption.

Empty and full subsets are included. `emptySpectralMonomial` gives $\prod_i a_i$; `fullSpectralMonomial` gives $\prod_j b_j$. Spectra may contain zeros or repeats: the proof uses no division or spectral nondegeneracy. The scalar statements also permit an empty ambient vertex type/list. The concrete forest interface still requires `edges : Nat → Edge vertices`, hence an ambient enumeration with two distinct endpoints even for a zero-length schedule; this inherited limitation is not hidden.

**What is still missing:** a checked theorem identifying

$$
\det(D_a+UD_bU^*)=
\sum_{J,K}\lvert\det U[K,J]\rvert^2 a^{K^c}b^J
$$

for the repository's actual matrix operations and determinant. The scalar theorem above does **not** prove this matrix-side formula, assume it, or rename the subset sum as a determinant. Proving diagonal-plus-matrix determinant expansion and the required general minor contraction remains a separate obligation. The current Gaussian-rational scalar interface also does not provide an arbitrary-complex transfer or a complete formal MDO theorem.

## Reproduce the exact audit

```sh
python3 -I -B research/marcus-de-oliveira/forest-spectral-expansion/check.py /path/to/lean
```

Use Lean **4.24.0**, as pinned by [lean-toolchain](../../../formal/lean-toolchain). The checker builds every transitive local import in a fresh temporary directory with warnings treated as errors; it does not use previously built `.olean` files, download packages, or require Mathlib. It checks the source hashes and axiom reports in [audit.json](audit.json), the 25-name [declaration inventory](declarations.json), and exact Lean-rendered types in [statements.txt](statements.txt). The allowed trusted axioms are only `propext`, `Classical.choice` and `Quot.sound`; no `sorry`, `admit`, `native_decide`, or new axiom is used.

The source-specific audit is separate from repository-wide CI registration. A compiler success applies to the exact statements and dependencies recorded here, not to the unproved matrix-side formula.
