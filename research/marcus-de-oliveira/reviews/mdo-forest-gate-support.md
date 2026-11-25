# Review: `mdo-forest-gate-support`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-reductions` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-forest-gate-support-final-published-20261125`, 2026-11-25 20:42:14 UTC, the checkpoint the task records as submitted (submit event `mdo-forest-gate-support-final-submit-20261125`, 20:42:27 UTC). No later work on the task exists. The commit is the head of PR #21, merged as `5e8d0e9`; all four files below are byte-identical on `main` at `1b98a07`.
- **Artifact** (commit `563d5c0b5379f260a38b2b169831464acf15bd93`):
  - [`formal/ForestGateSupport.lean`](https://github.com/g4nesh/math-commons/blob/563d5c0b5379f260a38b2b169831464acf15bd93/formal/ForestGateSupport.lean), SHA-256 `e5f1c920b17e385ad58dffb070f0647d5d937a06c4916a2a2a417aaccdda286f` (matches the hash in the checkpoint)
  - imported, compiled with it: [`formal/ForestTrajectory.lean`](https://github.com/g4nesh/math-commons/blob/563d5c0b5379f260a38b2b169831464acf15bd93/formal/ForestTrajectory.lean), SHA-256 `b5dcffb2a165702ef46a2a67e5a7d9049a5f5ba3af3cfd9cc2491af5e0fc5910`; [`formal/ForestCut.lean`](https://github.com/g4nesh/math-commons/blob/563d5c0b5379f260a38b2b169831464acf15bd93/formal/ForestCut.lean), SHA-256 `0bfc29b6c4470e3c11cfdbdb561616f55bc3d73af78d53e8edc68feb44b1b5ae`
  - checker: [`formal/check.py`](https://github.com/g4nesh/math-commons/blob/563d5c0b5379f260a38b2b169831464acf15bd93/formal/check.py), SHA-256 `9e58a7b6b9980f6755c06263789c2e1e9355b032403b427e547150e9998b1720`
- **Verdict:** **pass**, with one minor wording note (§2) and scope notes (§3).
- **Evidence level of this review:** independent review of a machine-checked artifact: kernel recompilation of the exact source on the pinned toolchain (Lean 4.24.0), axiom audit, a statement-by-statement comparison with the claim, a hand re-derivation, a small Lean probe for non-vacuity, and independent exact Gaussian-integer computations. This review adds no new formal proof.

## 1. What was checked, and how

**The claim, restated.** Let $M$ be a square matrix over any scalar type with the four zero laws $0+0=0$, $0x=x0=0$, $-0=0$, and suppose every off-diagonal entry of $M$ outside the $2\times2$ block on an edge $\{l,r\}$ vanishes. If the signed Laplace minor of $M$ with rows the output subset $T$ and columns the input subset $S$ is nonzero, then $T=S$ or $T$ is $S$ with $l$ and $r$ exchanged. If $l,r$ have equal occupancy, $T=S$. Composed with `forestTrajectoryUnique`: on an indexed schedule in which deleting each edge separates its endpoints, at most one subset history between fixed endpoints has a nonzero ordered product of such minors.

I exported the four files above from the commit with `git archive`, then compiled and ran everything in that copy.

| Claim in the artifact | How checked | Result |
| --- | --- | --- |
| 23 theorems compile with warnings as errors; axioms only `propext`, `Classical.choice`, `Quot.sound`; no `sorry` or `native_decide` | Ran `formal/check.py` (exit 0, 23 `ForestGateSupport` reports). Separately compiled the three modules by hand with `-DwarningAsError=true`: no output except axiom reports. Searched all three files for `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by`, `extern`, `unsafe`, `macro`, `set_option`, attributes: none. | Holds. `emptyMinor` uses no axioms. |
| `minorDet` is a concrete signed Laplace expansion | Read the definition: `columnChoices` gives the $j$-th listed column parity $j \bmod 2$ and keeps the other columns in order, so this is expansion along the first row; nonsquare input gives `zero`, the empty minor gives `one`. Lean `#eval` on seven minors of a fixed $5\times5$ integer matrix (orders 0 to 5, two nonsquare) equals my Leibniz-sum determinant exactly, sign included. | Holds. Identification with a library determinant is an open obligation, as the artifact says. |
| `nonzeroMinorStep`: nonzero minor gives `Step e s t` | Re-derived: a nonzero sum has a nonzero term; a nonzero signed product has two nonzero factors (only the four zero laws are used); induction gives a bijection rows to columns along nonzero entries (`Matching`), hence $\lvert T\rvert=\lvert S\rvert$. For $x\notin\{l,r\}$ the nonzero entry in row $x$ (or column $x$) must be diagonal, so $T$ and $S$ agree off $\{l,r\}$. Equal counts then leave stay or swap (the 16 Boolean cases). The printed Lean statement says exactly this, with rows `vertices.filter t` and columns `vertices.filter s`. | Holds. No unitarity, invertibility or diagonal condition is assumed. |
| Same, exact instances | Check A: $n=2..5$, every ordered edge, six random gates each with block entries in $\{-2..2\}+i\{-2..2\}$ (zeros and singular blocks included), identity or arbitrary diagonal off the block; all 36,072 square minors. | All 7,559 nonzero minors are stays or swaps; 2,562 lie in collision sectors and all are stays. |
| The result cannot be strengthened | Check B: for gates with identity off the block and $a,b,c,d,ad-bc\neq0$, every legal stay or swap minor. | All 7,992 nonzero. |
| Both-occupied, neither-occupied, empty, full cases | `nonzeroMinorCollision` takes $s(l)=s(r)$, which covers both and neither occupied; `nonzeroMinorEmpty` and `nonzeroMinorFull` specialize it. | Holds. |
| Gaussian-rational instance | `gaussianArithmetic.mul` is $(a,b)(c,d)=(ac-bd,\,ad+bc)$, i.e. complex multiplication. Lean probe: a concrete swap gate on the path $0-1-2$ with $S=\{0\}$, $T=\{1\}$; the minor evaluates to $(1,0)$, `gaussianGateMinorStep` yields a `Step`, and it is the swap branch. | Holds; hypotheses are satisfiable. |
| `nonzeroAmplitudesUnique` | `pathAmplitude` multiplies the step minors in order; `nonzeroPathFactors` extracts nonzero factors via the zero laws, then `forestTrajectoryUnique` applies. Lean probe proves `EveryEdgeIsBridge` for the path $0-1-2$ (two edges), so the hypothesis is not vacuous. | Holds. |
| Same, exact instances and controls | Check C: 48 random forest schedules ($n=2..5$, random edge order and orientation, entries as in A), every endpoint pair of every size; minors from the Leibniz sum, not from the support theorem. Check D: a triangle and a repeated edge. | C: at most one nonzero history per endpoint pair; for every endpoint pair reached by a nonzero history, the sum over histories also equals the compound entry of the product (pairs with no nonzero history are not compared). D: two nonzero histories each, so the bridge hypothesis cannot be dropped. |
| Scope statements and remaining obligations (file header, `formal/README.md`, checkpoint) | Compared each sentence with the Lean statements. | Accurate. The listed obligations (library determinant, complex instance, unitary magnitudes, Cauchy–Binet, convex weights) are the real remaining gaps. No novelty or conjecture claim is made. |

The artifact cites no external results, so no literature check was needed. The underlying fact (a nonzero determinant needs a nonzero permutation term) is standard and is proved here, not cited.

**Reproduction.** Run from the repository root (`formal/` unchanged from `563d5c0`), with Lean 4.24.0 on `PATH`:

```sh
shasum -a 256 formal/ForestGateSupport.lean   # e5f1c920...a286f
python3 formal/check.py                       # exit 0; 23 ForestGateSupport axiom reports, standard axioms only
python3 -I research/marcus-de-oliveira/reviews/check_forest_gate_support_ishaan.py --lean
# A: 36072 square minors, 7559 nonzero (2495 swaps, 2562 in collision sectors): all stay or swap
# B: 7992 legal stay/swap minors of generic gates: all nonzero
# C: 48 random forest schedules (n=2..5): at most one nonzero history per endpoint pair; Cauchy-Binet sums match
# D: controls without the bridge hypothesis reach 2 (triangle) and 2 (repeated edge) nonzero histories
# E: Lean probe compiled; bridge and swap examples use only standard axioms; minorDet = Leibniz on 7 minors [-20741, 91, 66, -3842, 0, 0, 1]
```

The script uses only the standard library, seed 20261125, and exact integer arithmetic. Checks A to D run in under a second; `--lean` compiles the three modules and the embedded probe in a temporary directory (a few seconds). Without `--lean` it needs no Lean.

## 2. Findings

No defect was found in the Lean artifact or its stated scope.

**Minor wording note (non-blocking).** The overview row for MDO-FOREST, changed in the same commit, says Lean checks that at most one forest path "has a nonzero product of gate minors". There, "gate minors" means the Laplace recursion `minorDet` over the zero-law interface, and "forest" means `EveryEdgeIsBridge` on an indexed schedule. The Lean header, `formal/README.md` and `context/forest.md` state that the link to the standard complex determinant is open; the overview row does not. Smallest correction: append "(for the Lean Laplace expansion; matching it to the standard determinant is open)" to that row.

## 3. Scope notes for users of this result (not defects)

1. **Support only.** The theorems say which minors can be nonzero. They say nothing about magnitudes, unitarity, or the sum over histories. The Cauchy–Binet agreement in check C is my consistency check, not a claim of the artifact, and is not a proof.
2. **Ordering.** Minors are taken with rows and columns in the order of the `vertices` enumeration. This does not affect support. A later Cauchy–Binet bridge must use the same enumeration for every gate, which the current definitions do.
3. **Dependencies.** `nonzeroAmplitudesUnique` uses `forestTrajectoryUnique` from `ForestTrajectory.lean` (task `mdo-cut-subset-bridge`, still in review). The kernel checked the composed statement, so it holds as stated. I read the definitions used here (`Edge`, `Subset`, `Step`, `swapSubset`, `ReachWithout`, `EveryEdgeIsBridge`) but did not review that file's proofs beyond kernel acceptance.

## 4. What this pass does and does not establish

**It establishes** that the reviewed checkpoint's 23 theorems compile on Lean 4.24.0 with warnings as errors and only the standard axioms, and that their statements match the task. A nonzero signed minor of a two-coordinate gate forces a stay or a swap, including the collision, empty and full cases. Along an indexed schedule satisfying `EveryEdgeIsBridge`, at most one history between fixed endpoint subsets has a nonzero product of such minors.

**It does not establish:**

- that `minorDet` equals Lean or Mathlib's matrix determinant (agreement was checked on seven exact instances only);
- a complex-number instance, unitary transition magnitudes, Cauchy–Binet, or the convex permutation weights;
- that `EveryEdgeIsBridge` is equivalent to a standard acyclicity definition (I showed only that it is satisfiable and that it cannot be dropped);
- correctness beyond the Lean 4.24.0 kernel and the three standard axioms;
- anything about the truth of the Marcus–de Oliveira conjecture.
