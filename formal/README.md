# Lean verification

Lean 4.24.0 checks an exact three-dimensional instance of the Marcus–de Oliveira conjecture. The proof starts from the spectra and relative eigenbasis, computes the determinant and the six permutation products, and constructs a convex combination for the determinant. The tree on the website draws these data. The [Lean source](MathCommons.lean) uses Lean's standard library.

## Reproduce it

From the repository root, with [elan](https://lean-lang.org/install/) installed, run:

```sh
python3 formal/check.py
```

The script uses the toolchain pinned in this directory and treats warnings as errors. It checks the transitive axiom report for each named theorem. The allowlist contains Lean's standard propositional extensionality, classical choice and quotient soundness. A missing report or an extra axiom fails the check. The same command runs in the required GitHub check. The executable path can be supplied as the first argument when using a standalone installation of the pinned version.

## What is checked

Write the spectra as $a=(0,1,i)$ and $b=(-1-i,1,-1+i)$. The relative eigenbasis is the real symmetric matrix with diagonal entries −1/3 and off-diagonal entries 2/3. The formalization represents complex rationals as pairs of rationals with explicit complex multiplication.

| Lean declaration | Result |
| --- | --- |
| `allPairings` | Exhausts all triples of distinct indices from a three-element set. Each triple defines one pairing of the two spectra. |
| `orthogonal` | Checks all nine entries of the row Gram matrix against the identity. |
| `determinantValue` | Evaluates the three-by-three Leibniz determinant of the diagonal matrix plus the conjugated second diagonal matrix. Its coordinates are (2/3, 2/9). |
| `witnessInHull` | Constructs two binary convex combinations of three computed endpoints. The endpoints (2,1), (0,−1), (0,2) have weights 9/27, 13/27 and 5/27. |
| `inEveryConvexSet` | Proves by induction that every constructed certificate belongs to every set containing the vertices and closed under rational convex combinations. |

The rational coordinates embed in the real plane. Addition, multiplication and nonnegative rational interpolation have their usual meaning under this embedding. The code formalizes the rational operations and interpolation predicate directly. The bridge to Mathlib's complex matrices, spectral theorem, determinant and real convex hull is still a separate formalization task. The [independent arithmetic implementation](../research/marcus-de-oliveira/literature/check_source_examples.py) reproduces the example.

## Contributing a theorem

1. Claim the specific lemma or translation in the live ledger. State its quantifiers, all hypotheses, the Lean declaration you intend to prove and the source file you own. Separate mathematical hypotheses from any unproved bridge from the paper's statement.
2. Pin dependencies. For Mathlib, commit its revision and Lake manifest with a matching Lean toolchain. Agree on matrix indices, normality, conjugate transpose, the scalar field for convexity and the direction of each permutation before building dependent results.
3. Cover repeated eigenvalues, singular sums, point or segment hulls and the intended dimension range. If a proof divides by a quantity, prove it nonzero or give the zero case. A strict-separation argument needs its exact sign and quantifiers.
4. Add a reproducible checker to the required CI job. For this project add a named `#print axioms` report, list the declaration in `check.py`, and keep warnings fatal. Inspect the statement and imports as well as its axiom report. Drafts with placeholders belong in clearly marked notes until they pass. Avoid `native_decide` or introducing an axiom to discharge the claim. The current finite computations use kernel reduction.
5. Checkpoint the declaration names, exact commit, local checker result, and unresolved dependencies. Submit the evidence to a different GitHub owner, who should rebuild it and compare the formal statement against the mathematical claim. Merging the source makes it available to collaborators. A ledger review records the independent verdict.

A natural next lemma is the relation between the rational interpolation constructor and convexity of subsets of the complex plane over the reals. A separate, larger direction is to formalize the coefficient identity in the [forest-rotation argument](../research/marcus-de-oliveira/baselines/README.md). Neither the finite witness nor that restricted family establishes the unrestricted conjecture.
