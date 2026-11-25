# Review: `mdo-cut-subset-bridge`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-reductions` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-cut-subset-bridge-final-published-20261125`, recorded 2026-11-25 19:52:31 UTC. The submit `mdo-cut-subset-bridge-final-submit-20261125` followed at 19:52:42 UTC. The ledger shows no later checkpoint for this task.
- **Artifact:** [`formal/ForestTrajectory.lean` at `bbbf34f`](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/formal/ForestTrajectory.lean), SHA-256 `b5dcffb2a165702ef46a2a67e5a7d9049a5f5ba3af3cfd9cc2491af5e0fc5910`. This matches the hash in the checkpoint. Its only import is [`formal/ForestCut.lean` at the same commit](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/formal/ForestCut.lean), SHA-256 `0bfc29b6c4470e3c11cfdbdb561616f55bc3d73af78d53e8edc68feb44b1b5ae`. Both files are byte-identical on `main` at `1b98a07`. The toolchain pin is `leanprover/lean4:v4.24.0`.
- **Verdict:** **pass**, with three minor notes (§2). None of them affects the formal content.
- **Evidence level of this review:** an independent check of each statement against the claim it stands for, plus a kernel re-check of the exact source. It also includes my own Lean restatement, a non-vacuous instance and a negative control, an exhaustive finite model check, and a seeded floating-point check of the gate-minor support that the model assumes. The forest-product theorem as a whole is not formally verified here.

## 1. What was checked, and how

**The claim, restated.** Fix a finite vertex type `V`, listed without repeats, and indexed edges `e_0,…,e_{n-1}` with distinct endpoints. Suppose that for each `j < n`, the endpoints of `e_j` are not connected by the other `n-1` edges. Take two histories `a, b : ℕ → 2^V` in which each step `k < n` either keeps the subset or exchanges the memberships of the endpoints of `e_k`. If `a 0 = b 0` and `a n = b n`, then `a k = b k` for every `k ≤ n`. This is the "Uniqueness of the subset path" paragraph of `research/marcus-de-oliveira/baselines/README.md` §2, at the same commit.

**My derivation.** I argue by induction on `k`, assuming `a k = b k`.
- Let `C` be the set of vertices reachable from the left end of `e_k` without using `e_k`. By hypothesis, the right end of `e_k` is not in `C`.
- Every other edge has both ends in `C` or both outside it, because one end in `C` puts the other in `C`. Steps `t ≠ k` therefore preserve `|x ∩ C|`.
- Hence `|a(k+1) ∩ C| = |a n ∩ C| = |b n ∩ C| = |b(k+1) ∩ C|`.
- If `a(k+1) ≠ b(k+1)`, then one of them is `a k` and the other is its swap, and they differ. That happens only when exactly one endpoint is occupied, and then the swap changes `|· ∩ C|` by ±1, a contradiction.

The Lean proof has the same structure: `cutsOfBridges`, `subsetCountsSingleOpportunity` with `ForestCut.cutHistoryUnique`, and `nextSubsetUnique`.

| Item in the artifact or checkpoint | How I checked it | Result |
| --- | --- | --- |
| Identity: the SHA-256 above, 17 theorems, 17 `#print axioms` lines | `shasum -a 256`; read the file; counted declarations and report lines | Holds |
| Compiles with `-DwarningAsError=true` and uses only `propext`, `Classical.choice` and `Quot.sound` | Compiled `ForestCut` and then `ForestTrajectory` into a fresh `LEAN_PATH` in an exported copy, with `lean` 4.24.0 | Holds. Exit 0. All 17 reports and the 8 from `ForestCut` lie within the allowlist, and none mentions `sorryAx` |
| Nothing could fake the result | Searched both files for `sorry`, `admit`, `axiom` declarations, `native_decide`, `unsafe`, `implemented_by`, `extern`, `opaque`, `set_option`, `macro`, `syntax`, `notation`, `instance`, `#eval` and `IO.` | No matches |
| `Edge`: a duplicate-free, exhaustive enumeration split as `left :: right :: rest` | Read `Edge`, `countSplit` and `restCountSwap` | Holds. `cutCount` counts each vertex once, so "cardinality" in `trajectoryPreservesCardinality` means `countP` over the full enumeration |
| `Step` (stay, or swap the two endpoints) models the support of a two-coordinate gate | By hand: outside `{p,q}` the rows and columns of `G` are unit vectors, so `det G[R',R] ≠ 0` forces `R'∖{p,q} = R∖{p,q}`, and equal sizes then force `R' ∈ {R, swap R}`. Numerically: 12,024 minors of random gates for `d ≤ 5` and every pair (seed 20261125) | Holds. The squared magnitudes are 1, `1−q` and `q`, as in the baselines table. This step is not formalized in the file, which says so |
| `EveryEdgeIsBridge` is the forest hypothesis and is not vacuous | In Lean, I proved that `ReachWithout` implies an ordinary walk along edges `k < n`, `k ≠ j` (`toWalk`). By hand: every edge is a bridge exactly when the multigraph has no cycle, a parallel pair counting as a cycle. Exhaustive check of all 9,701 edge sequences with `d ≤ 4` and `m ≤ 5`: bridge equals acyclic in every case. In Lean, the path 0–1–2 satisfies it (`pathBridges`) | Holds |
| Cut construction: `edgeSide`, `edgeSideSeparates`, `survivingEdgeSameSide`, `cutsOfBridges` | Matched them against the derivation above | Holds. The cut is constructed from the hypothesis, not assumed |
| Collisions and the empty and full subsets (`swapCollision`, `swapEmpty`, `swapFull`, `equalCountSwapCollision`) | Read and compiled them | Holds. The general theorem covers these cases without special treatment |
| `forestTrajectoryUnique` states the claim | My theorem `restated` spells out the swap rule as an `if` expression, takes the walk-based bridge hypothesis, and is proved directly by the author's theorem. `pathForcesMiddle` uses it to show that any legal path from `{0}` to `{2}` passes through `{1}` | Holds |
| The hypothesis is needed | `triangleControl` (Lean): on the triangle, the bridge condition fails, and two legal histories from `{0}` to `{2}` differ at time 1. Exhaustively, all 145 forests in range have unique histories, and all 9,556 non-forests have two histories with the same endpoints | Holds |
| Scope statements in the file docstring and checkpoint | Read them | Accurate. The file contains no matrices, exterior powers, Cauchy–Binet, amplitudes or convex weights. The equivalence with a library acyclicity predicate and the automatic construction of the enumeration are stated as not formalized |
| Attribution | Read both files and baselines §2 | The file uses no external results and claims no novelty. It formalizes the repository's own written argument |

**Reproduction.** Run from the repository root. The two Lean files are identical at `bbbf34f` and on `main`. The check script is part of this review's commit.

```sh
shasum -a 256 formal/ForestTrajectory.lean formal/ForestCut.lean
B=$(mktemp -d)
(cd formal && LEAN_PATH=$B lean -DwarningAsError=true -o $B/ForestCut.olean ForestCut.lean \
  && LEAN_PATH=$B lean -DwarningAsError=true -o $B/ForestTrajectory.olean ForestTrajectory.lean)
# 17 lines "'MathCommons.ForestTrajectory.…' depends on axioms: [...]", each within
# propext, Classical.choice, Quot.sound

python3 -I research/marcus-de-oliveira/reviews/check_cut_subset_bridge_ishaan.py --lean
# [A] gate support and transition table: 12024 minors checked, all as stated
# [B] 9701 edge sequences (d <= 4, m <= 5): bridge == forest in all; 145 forests, all with
#     unique histories; 9556 of 9556 non-forests have two histories with equal endpoints
# [L] sha256 pinned; ForestCut and ForestTrajectory compile with warnings as errors; 17
#     ForestTrajectory axiom reports within [...]; reviewer audit (...) compiles
```

The script needs Python 3 with NumPy, plus `lean` 4.24.0 for `--lean`, and runs in a few seconds. In `--lean` mode it checks both SHA-256 hashes and compiles the two modules in a temporary directory. It then compiles my own Lean check file, which is embedded in the script: `Walk`, `toWalk`, `restated`, `pathBridges`, `pathForcesMiddle` and `triangleControl`. Its four axiom reports are within the same allowlist.

## 2. Findings (all minor, none blocking)

1. **At the reviewed commit, neither the repository checker nor CI covers this file.** At `bbbf34f`, `formal/check.py` lists only the `MathCommons` and `ForestCut` declarations. I ran it in an exported copy: exit 0, 13 axiom reports, none from `ForestTrajectory`. The Lean step of CI runs only that script. `formal/README.md` at that commit still describes the step from the count-history hypotheses to actual trajectories on a forest as a written argument.

   The checkpoint is accurate here: it says the compile was run by hand. Commit `97adf40` (19:53 UTC, after the submit, merged to `main` as `abad852`) adds the 17 names to `check.py` and a README section. I read that change but did not run it.

   *Correction:* none needed on `main`. Anyone citing `bbbf34f` itself should use the commands above.
2. **The checker does not pin statements.** `check.py` matches declaration names and axiom lists only (the regex at lines 32–38). A weakened theorem under the same name would still pass.

   *Smallest fix:* commit the expected `#check @MathCommons.ForestTrajectory.forestTrajectoryUnique` output, or a separate restatement like `restated`, and have CI compare against it.
3. **The remaining matrix-minor bridge is named but not stated precisely in the reviewed artifact.** The task asks for it to be specified. The checkpoint lists the remaining obligations: the gate-minor, exterior-power, Cauchy–Binet and convex-weight bridge, the library acyclicity translation, and the enumeration construction. The fact itself appears at the reviewed commit only as a step of the written proof in baselines §2 (a nonzero gate entry either keeps the subset or exchanges the endpoints when exactly one is occupied), and neither the Lean file nor the checkpoint points to it. It is stated as the next lemma ("a nonzero minor of a two-coordinate gate forces a legal `Step`") only after the checkpoint, in request `mdo-followups-ask-gate` (19:53:50 UTC). `proof-obligations.md` at `97adf40` names it more loosely, as connecting the support of the gate minors to the legal transitions.

   I checked that the lemma is true, by hand and on 12,024 minors (row [A]), and that it matches `Step`.

   *Optional fix:* one sentence in the module docstring stating that lemma and the transition table.

## 3. Limits of this review

- **Trusted base.** The trusted base is the Lean 4.24.0 kernel and elaborator from the pinned toolchain. I used no independent kernel replay and no second Lean version.
- **The `Edge` encoding.** `restated` still relies on the author's `Edge` structure, a list enumeration of a finite type. No translation to a library graph, `Fintype` or `Finset.card` notion is formalized, as the artifact itself says. That "every edge is a bridge" equals "forest" is checked by hand and on a finite domain only.
- **Floating-point checks.** Row [A] is seeded floating point with tolerance `1e-12` and four random gates per pair. It supports the specification of the next step. It is not part of what this pass certifies.
- **Not reviewed.**
  - `ForestCut.lean` as its own task (`mdo-forest-cut-lean`). I only compiled it and checked the statements used here: `SingleOpportunity`, `constantOnInterval` and `cutHistoryUnique`.
  - Later files on `main`, such as `ForestGateSupport.lean`.
  - The single-cycle part of PR #17.
  - The written Cauchy–Binet and probability steps of baselines §2.
- **What this pass establishes and does not.**
  - It establishes that the subset-path uniqueness step of the forest argument is machine-checked as stated, under a hypothesis equivalent to "the indexed edges form a forest".
  - It does not establish the forest-product theorem, any case of the Marcus–de Oliveira conjecture, or the conjecture itself.
