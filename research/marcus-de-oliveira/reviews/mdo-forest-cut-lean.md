# Review: `mdo-forest-cut-lean`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-reductions` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-forest-cut-lean-published-20261125`, recorded 2026-11-25 18:43:47 UTC and submitted by `mdo-forest-cut-lean-submit-20261125` at 18:43:58 UTC. The ledger has no later event for this task.
- **Artifact:** [`formal/ForestCut.lean` at `9123d54`](https://github.com/g4nesh/math-commons/blob/9123d545f43045ab192703edaf468262b4a7b688/formal/ForestCut.lean), SHA-256 `0bfc29b6…b1b5ae`. This matches the hash in the checkpoint. The file is byte-identical on `main` at the PR #15 squash merge `8c01009` and at `1b98a07`. I read it together with the files beside it at the same commit: [`formal/check.py`](https://github.com/g4nesh/math-commons/blob/9123d545f43045ab192703edaf468262b4a7b688/formal/check.py) (`6ad71168…5fea79`), [`formal/README.md`](https://github.com/g4nesh/math-commons/blob/9123d545f43045ab192703edaf468262b4a7b688/formal/README.md) (`d7c705b4…fabfadf`), `lakefile.toml` and `lean-toolchain` (`v4.24.0`).
- **Verdict:** **pass**. Four minor notes follow in §4. None of them affects a claim in the checkpoint.
- **Evidence level of this review:**
  - a kernel recheck of the exact file with Lean 4.24.0;
  - Lean statement pins and concrete instantiations that I wrote;
  - a re-derivation by hand;
  - an exhaustive finite test of the unformalized cut bridge on 304 small forests.

  This is not a formalization of the forest-product theorem.

## 1. Claim restated

Take any $f:\mathbb N\to\mathbb Z$ and any $j<n$. Suppose $f(k+1)=f(k)$ for every $k<n$ with $k\ne j$ (`SingleOpportunity f n j`). Then:

- `singleCutChange`: $f(j+1)-f(j)=f(n)-f(0)$.
- `cutHistoryUnique`: if $g$ satisfies the same hypothesis with the same $n,j$, and $f(0)=g(0)$ and $f(n)=g(n)$, then $f(k)=g(k)$ for all $k\le n$.
- `crossingDecisionUnique`: if decisions $x,y\in\{\text{stay},\text{enter},\text{leave}\}$, encoded as $0,+1,-1$, equal the actual jumps of $f$ and $g$, then $x=y$.

The supporting lemmas are `telescoping`, `constantOnInterval`, `totalChangeAtOpportunity`, `sameEndpointsSameCrossing` and `crossingDeltaInjective`.

**Hand check.** By induction, $\sum_{k<n}(f(k+1)-f(k))=f(n)-f(0)$. Under the hypothesis every term except $k=j$ vanishes, so the sum is the single jump. For $k\le j$, $f$ is constant on $[0,k]$, so $f(k)=f(0)=g(0)=g(k)$. For $k>j$, $f$ is constant on $[k,n]$, so $f(k)=f(n)=g(n)=g(k)$. The map $0,1,-1$ is injective on the three decisions. The Lean proofs follow exactly this structure.

## 2. What was checked, and how

| Check | Method | Result |
| --- | --- | --- |
| Exact file | SHA-256 of `git archive 9123d54 formal` | Matches the checkpoint. |
| Kernel acceptance | `lean -DwarningAsError=true ForestCut.lean` with Lean 4.24.0 (commit `797c613`) on the exported copy | Exit 0 in 0.6 s, with no warnings. |
| Axioms | The file's eight `#print axioms` reports | Seven theorems use `[propext, Quot.sound]`. `crossingDeltaInjective` uses `[propext]`. There is no `sorryAx` and no `Lean.ofReduceBool`, so no `sorry` and no `native_decide`. This matches the checkpoint. |
| Source hygiene | I read all 153 lines. I grepped for `sorry`, `axiom`, `admit`, `native_decide`, `set_option`, `macro`, `notation`, `instance`, `#eval`, `unsafe` and `extern`. I scanned the non-ASCII characters. | Nothing was found apart from the `#print axioms` lines and one docstring sentence that says the hypotheses are not axioms. The only import is `Std`. The non-ASCII characters are standard Lean notation (`→`, `≤`, `∀`, `∧`, `≠`, `·`), one en dash in the docstring, and the bound type variable `α`. All declared names are ASCII, so there are no homoglyphs. |
| Statements say what is claimed | I wrote eight `theorem pinX : ⟨my restatement⟩ := X` declarations. Separate `rfl`/`Iff.rfl` checks fix `totalChange`, `SingleOpportunity` and `crossingDelta`. In a negative control, a pin with the sides of the conclusion swapped is rejected. | All pins are accepted. |
| Hypotheses satisfiable (non-vacuity) | I used the concrete history 2,2,2,3,3,3 with $n=5$, $j=2$, plus a second, definitionally different history with the same values. Every hypothesis of `cutHistoryUnique` and `crossingDecisionUnique` is discharged. The decision is forced to `enter`. The opportunity at the last step ($j=n-1$) also works. | Accepted, with only `propext` and `Quot.sound`. |
| Hypothesis not redundant | The histories 0,1,0 and 0,0,0 have equal endpoints and different first jumps. The history 0,1,0 satisfies neither `SingleOpportunity … 2 0` nor `… 2 1`. | Accepted. |
| Boundary | With $n=0$, `SingleOpportunity f 0 j` is false. | This is harmless. With zero gates there is no crossing to determine. |
| Unformalized bridge, finite test | This test is outside Lean, in my own Python. It covers 4 fixed forests (empty, a path in two orders, a star) and 300 seeded random forests with 2 to 7 vertices, random gate order and random orientation. For every start subset I enumerated every stay/exchange trajectory under the transition table of `baselines/README.md` §2. For every edge and both sides of its cut, I tested the count history. | All 608,844 cut histories satisfy `SingleOpportunity` with a jump in $\{-1,0,1\}$. No two trajectories share both endpoints. A triangle and a repeated edge give 8 and 4 clashing endpoint pairs, so there the hypothesis correctly fails. |
| Repo checker at the reviewed commit | `python3 -I formal/check.py` in the export | Exit 0, but only the 5 `MathCommons.lean` reports appear. `ForestCut.lean` is not compiled (§4, note 1). |
| Lake | `lake build ForestCut` in a scratch copy | `unknown target` (§4, note 2). |
| README | At `9123d54` the README does not mention `ForestCut`. I read the section "An unbounded cut-history lemma", added in `8c01009`, against the Lean statements. | At `8c01009` the section is accurate. It names the same three conclusions and states that the exterior-power trajectory bridge remains a written argument. On `main` (`1b98a07`) the heading and first paragraph are unchanged, but that second paragraph has been replaced by a section on a later file, `ForestTrajectory.lean`, which I did not review. |

## 3. How much of the argument this covers

The target is the paragraph "Uniqueness of the subset path" in [`baselines/README.md` §2](https://github.com/g4nesh/math-commons/blob/9123d545f43045ab192703edaf468262b4a7b688/research/marcus-de-oliveira/baselines/README.md). Its steps are:

1. Deleting forest edge $e_r$ gives a cut that no other gate crosses. **Not formal.** It enters only as the hypothesis `SingleOpportunity`.
2. The occupied-count history on side $C$ is an integer history. **Not formal.** The file has no graphs and no subsets.
3. The count therefore changes only at step $r$, by $\delta_r=\lvert S\cap C\rvert-\lvert T\cap C\rvert$. **Formal**: `singleCutChange` and `totalChangeAtOpportunity`.
4. $\delta_r$ fixes the stay/enter/leave move at step $r$. **Formal**, but only given that a decision exists and represents the jump (`crossingDecisionUnique`). The case $\lvert\delta_r\rvert>1$ ("no path exists") is not stated.
5. Together, the decisions for all edges determine every intermediate subset. **Not formal.** `cutHistoryUnique` gives agreement of *counts* on one cut, not of subsets. The file says this.

The amplitude and probability table, compound minors, Cauchy–Binet, identities (5) to (8), and the convex-combination conclusion are also not formal. The file header, the checkpoint detail and the README at `8c01009` all state these exclusions correctly.

"Not formal" here means not formalized in `ForestCut.lean` at the reviewed commit. The README on `main` describes a later file, `ForestTrajectory.lean`, that builds on `ForestCut` for steps 1, 2 and 5. That file belongs to other work and is not covered by this review.

The formal content is elementary telescoping. It is correct and stated at full generality: every finite length, every opportunity index, and arbitrary integer values. Its mathematical weight is one inference in the proof of a restricted family, namely forests of two-coordinate gates. It says nothing about the conjecture in general.

My finite test supports the claim that the hypotheses are the right ones for the intended forest application. It is evidence on small instances, not a proof of steps 1, 2 or 5.

## 4. Notes (non-blocking)

1. **The checker did not cover the file at the submitted commit.** At `9123d54`, `formal/check.py`, and so the CI step "Check the formal statements and their axiom dependencies", compiles only `MathCommons.lean`. The checkpoint's "required command" is the direct `lean` call, which I reproduced. The squash merge `8c01009` adds all eight declarations to `check.py` and the README section. `ForestCut.lean` is byte-identical there. I read that `check.py` but did not run it. The task description assigns CI and docs to the coordinator, so this checkpoint needs no correction.
2. **No lake target.** `formal/lakefile.toml` declares only `lean_lib MathCommons`, so `lake build ForestCut` fails with `unknown target`. The lakefile is unchanged on `main`. Verification does not depend on lake. Smallest fix: add a `[[lean_lib]]` entry for `ForestCut`, or say in the README that only `check.py` covers it.
3. **The checker does not pin statements.** Both versions of `check.py` parse only the names and axiom lists from `#print axioms` output. A weakened theorem under the same name would still pass. Smallest fix: add type-ascribed pins like those in my script (`theorem pin : ⟨statement⟩ := name`) to a checked file.
4. **Wording.** The task title, "Formalize the unique forest-cut crossing lemma in Lean", reads as if the forest cut were formalized. Only the integer cut-count step is. The forest cut is a hypothesis. The file header and checkpoint detail are accurate. Downstream citations should say "cut-count arithmetic machine-checked", not "forest-cut lemma formalized".

## 5. Reproduction

From the repository root (`ForestCut.lean` is byte-identical at `9123d54` and on `main`):

```sh
lean --version                      # Lean (version 4.24.0, ...)
(cd formal && lean -DwarningAsError=true ForestCut.lean)
# eight axiom reports: seven [propext, Quot.sound], crossingDeltaInjective [propext]
python3 -I research/marcus-de-oliveira/reviews/check_forest_cut_lean_ishaan.py
# Lean file: accepted, 8 axiom reports, all within ['Quot.sound', 'propext']
# Lean file+review: accepted, 11 axiom reports, all within ['Quot.sound', 'propext']
# Bridge: 304 forests, 608844 cut histories satisfy SingleOpportunity, no endpoint clash
# Control triangle: 8 endpoint pairs with two trajectories (hypothesis fails, as expected)
# Control repeated edge: 4 endpoint pairs with two trajectories (hypothesis fails, as expected)
# All checks passed.
```

The script uses only the standard library and seed 20261125. It takes about 3 s, and repeated runs give identical output. It refuses any `ForestCut.lean` whose SHA-256 differs from the reviewed one. It writes the Lean test file only to a temporary directory. To see note 1, check out `9123d54` and run `python3 formal/check.py`: it prints only the five `MathCommons.*` reports.

## 6. Limits of this review

- The trusted base is the Lean 4.24.0 binary and its `Std`. I did not replay the proofs in an independent kernel checker.
- The bridge test is finite: up to 7 vertices, 304 forests, the combinatorial stay/exchange model only. It does not touch matrix amplitudes. It does not establish steps 1, 2 or 5 of §3.
- I did not review other files in `formal/`, or the merged versions of `check.py` and the README, beyond the comparison in §2 and note 1.
- No external sources are cited by the artifact, and telescoping is standard, so there was no attribution to check. The artifact claims no novelty.
- This pass covers the exact file above. It says nothing about the forest-product theorem as a whole or about the Marcus–de Oliveira conjecture.
