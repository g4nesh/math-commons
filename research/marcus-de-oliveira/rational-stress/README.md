# Exact rational stress screen

**Result:** 924 explicitly constructed noncommuting normal pairs satisfy the original Marcus–de Oliveira convex-hull inclusion exactly. This is a finite computational result, **not a proof of the conjecture**. No outside case was found.

Task `mdo-rational-stress` was created and claimed before computation by `ganesh-search`, GitHub owner `g4nesh`, on November 25, 2026. Work is restricted to this directory. Author-side differential replay is complete; this contribution still requires review by another GitHub owner. No paid or remote compute was used.

## Observed results

| Run / dimension | Cases | Strict polygon interior | Polygon boundary | Segment hull | Outside |
| --- | ---: | ---: | ---: | ---: | ---: |
| Main / 4 | 360 | 291 | 9 | 60 | 0 |
| Main / 5 | 360 | 297 | 3 | 60 | 0 |
| Main / 6 | 180 | 150 | 0 | 30 | 0 |
| Separate holdout / 7 | 24 | 24 | 0 | 0 | 0 |
| Total | 924 | 762 | 12 | 150 | 0 |

- Every input passed exact two-sided unitarity and normality checks; every commutator $[A,B]$ has a strictly positive recorded squared Frobenius norm. Every unitary has at least one nonreal entry.
- All $n!$ permutation points were checked: **181,440** in the main screen and **120,960** in the holdout, **302,400** in total. These counts include duplicate coordinates; the unique-point count is separately recorded.
- The main screen contains **450 general complex stress cases** and **450 structured controls**: 150 singular sums, 150 two-value spectra, and 150 Hermitian pairs. All 24 holdout cases are general complex stresses. These families are reported separately so known or immediate cases do not stand in for difficult coverage.
- Every inside certificate gives at most three nonnegative rational weights on explicitly identified permutation points, summing to one and reconstructing the exact determinant. All 924 certificates passed a second arithmetic implementation.
- **100 dimension-four inputs have an exact negative exotic functional**, including 60 of the general complex cases. This verifies that these particular matrices fall outside the stronger common squared-minor certificate condition. It does not make their determinants counterexamples to MDO.
- A separate [distinctness audit](distinct-inputs.json) found **924 distinct exact ordered matrix pairs $(A,B)$**, with no duplicates. It hashes only the diagonal of $A$ and the exact entries of $B=U\operatorname{diag}(b)U^*$, excluding case names, notes, and the representation $U$. This does not identify equivalence classes under unitary similarity or other MDO symmetries. The ordinary witness hashes include metadata and are used for replay integrity, not for this distinctness claim. The holdout did not change the completed main report.

Machine-readable artifacts: [main report](report.json), [main replay](replay.json), [holdout report](holdout/report.json), and [holdout replay](holdout/replay.json). The report stores each case's full seed/specification, canonical input SHA-256, exact determinant, hull vertices, certificate, commutator norm, and support margin. No floating tolerance enters a decision.

## Construction and coverage

Use $A=\operatorname{diag}(a)$ and $B=U\operatorname{diag}(b)U^*$. All arithmetic is in $\mathbb Q(i)$ using Python `fractions.Fraction`.

Each two-row gate is

$$G=\begin{pmatrix}c&-s\overline z\\sz&c\end{pmatrix},\qquad
c=\frac{1-t^2}{1+t^2},\quad s=\frac{2t}{1+t^2},\quad |z|=1.$$

The phase $z$ is selected from $i$, $(3+4i)/5$, and $(5-12i)/13$. A shuffled connected chain of gates is followed by additional random pairs. Thus the construction does not restrict itself to a product using each edge of a forest only once; extra or repeated edges can permit interference. This is a structured sample, not a Haar sample or an exhaustive cover of the unitary group.

The three main mixing families are:

1. **Generic:** gate parameters from $1/2,2/3,1/3,3/4$, with random signs.
2. **Near vertex:** first reorder $b$ so the identity permutation gives the lexicographically minimal permutation point, an exposed hull vertex. Then use small nonzero gate parameters $1/16$, $1/256$, or $1/1024$. Singular controls retain their designated zero coordinate instead of this reordering.
3. **Exotic centered:** begin with Zheng's rational orthogonal matrix

   $$O_*={1\over7}\begin{pmatrix}1&4&4&4\\-4&1&4&-4\\-4&-4&1&4\\-4&4&-4&1\end{pmatrix}$$

   embedded in the available coordinates, then apply small complex gates. The matrix is transcribed from [`../minors/moc/crosscheck_zheng.py`](../minors/moc/crosscheck_zheng.py). In dimension four the exact functional used there has $\psi(O_*)=-158/2401$; all 24 permutation matrices have $\psi\geq0$, which the tests check exactly. For higher dimensions we do **not** assert that an analogous obstruction survives the mixing. The dimension-four singular control has only three free coordinates, so its identically named mixing schedule has no $O_*$ block; it remains a labeled control.

The six spectra families are generic complex, thin complex perturbations of real spectra, tightly clustered complex spectra, an exactly singular sum with one fixed coordinate, a two-value $b$ spectrum, and Hermitian controls. In the singular family, $b_0=-a_0$ and $U$ fixes coordinate zero: the identity permutation already yields zero, so membership of that determinant is immediate. These cases exercise the singular and boundary code; they are not general-case evidence.

The main schedule is deterministic and round-robins dimensions so an early budget stop would not omit dimension six. The seed is `20261125`, with each case receiving a SHA-256-derived integer seed. There are 20 repetitions per family/mixing in dimensions four and five, and ten in dimension six. Every case's exact seed is saved; no optimization or adaptive stopping used the outcomes.

The separate **dimension-seven holdout** uses the distinct seed namespace `holdout20261125`: three general complex spectra families, two mixing families (generic and near vertex), and four repetitions, for 24 cases. All 24 input files are saved in [`holdout/witnesses/`](holdout/witnesses/). The existing public verifier deliberately remains limited to $n\leq6$. [`holdout.py`](holdout.py) uses its dimension-generic arithmetic and hull primitives under a separate explicit $n=7$ shape check; it does not relax the public parser or modify the main data.

## Exact distance from supporting lines

For each counterclockwise polygon edge $p\to q$, let $L$ be the larger coordinate extent of the hull. The report records the minimum exact value

$$m=\min_{p\to q}\frac{\operatorname{cross}(p,q,d)}{(|q_x-p_x|+|q_y-p_y|)L}.$$

This is a scale-normalized inward distance to the supporting lines in the $\ell^\infty$ metric. Positive means strict polygon interior; zero means polygon boundary. Segment hulls have a separate classification and no polygon margin. This quantity is **not** a robustness radius for perturbations of the input matrices or spectra.

The smallest margins among the general complex cases are below. Decimals here are display approximations to the exact rational values in the reports; they were not used in classification.

| Dimension | Case | Approximate positive margin |
| --- | --- | ---: |
| 4 | `n4-thin-complex-near-vertex-11` | $2.0167565\times10^{-12}$ |
| 5 | `n5-thin-complex-near-vertex-17` | $4.3359991\times10^{-11}$ |
| 6 | `n6-thin-complex-near-vertex-05` | $6.7478755\times10^{-11}$ |
| 7 | `n7-thin-complex-near-vertex-02` | $2.0089544\times10^{-9}$ |

The main [`witnesses/`](witnesses/) directory contains full inputs for all 12 polygon-boundary cases and the two tightest positive-margin cases per dimension/mixing, 30 files total. Every other main input regenerates exactly from its saved specification. Its rational convex certificate is already included in the report. The code would retain every outside input and its separating functional; this run produced none.

## Verification and budgets

The main computation calls the previously reviewed [`../certificates/verify.py`](../certificates/verify.py). It checks $UU^*=U^*U=I$, normality, the determinant by exact elimination, every permutation point, and the returned hull certificate. The extra commutator check excludes merely commuting examples from the reported noncommuting total.

[`replay.py`](replay.py) regenerates each input and checks its hash, then uses a separate Gaussian-rational representation and a **division-free subset dynamic-programming determinant**. It directly reconstructs convex combinations, or would test a separator against every permutation point. It also checks that all reported hull edges contain every permutation point on the correct side, exact margins, commutator norms, the dimension-four functional, saved witness files, and aggregate counts. It does not call the screen's determinant or hull decision. It shares the input generator and Python's `Fraction` implementation, so it is differential verification, not a wholly independent software stack or independent-owner review.

| Stage | Declared cap | Observed CPU | Outcome |
| --- | --- | ---: | --- |
| Main screen | 900 cases; 600 CPU / 720 wall seconds | 41.66 s | 900 complete |
| Main replay | 180 CPU / 240 wall seconds | 31.62 s | 900 checked; 0 mismatches |
| Separate n7 screen | 24 cases; 120 CPU / 150 wall seconds | 22.15 s | 24 complete |
| Separate n7 replay | 120 CPU / 150 wall seconds | 33.04 s | 24 checked; 0 mismatches |

The post-review distinctness audit regenerated only these existing inputs, used 14.57 CPU seconds, and introduced no new search cases. Its source is [`check_distinct.py`](check_distinct.py).

Six focused tests passed, including all 18 family/mixing combinations at dimension six with the smallest gate parameter, the exact $O_*$ functional and all permutation values, singular/row-swap determinants, tampered-certificate rejection, and the separate n7 API boundary. A three-case n4/5/6 pilot preceded the main run; those pilot cases are repetitions, not added to the reported distinct case count.

Observed environment: Python **3.9.6**, macOS 15.7.9 arm64, standard library only. No third-party packages are needed. Screening's `SIGALRM`/`SIGPROF` budget timers require a POSIX system. Replay uses no signal timers itself; the recorded replays were launched with external 180/120 CPU-second timers. The total work session was bounded to 30 minutes, with live ledger checkpoints at milestones.

## Reproduce

From the repository root:

```sh
# Fast regression tests, no repository outputs changed.
python3 -I -B -m unittest discover -s research/marcus-de-oliveira/rational-stress -p test_screen.py -v

# Independently replay every saved case; write fresh timing metadata outside Git.
python3 -I -B research/marcus-de-oliveira/rational-stress/replay.py \
  research/marcus-de-oliveira/rational-stress/report.json --output /tmp/mdo-rational-replay.json
python3 -I -B research/marcus-de-oliveira/rational-stress/replay.py \
  research/marcus-de-oliveira/rational-stress/holdout/report.json --output /tmp/mdo-n7-replay.json

# Recheck distinct exact (A,B) pairs; replaces local distinctness timing metadata.
python3 -I -B research/marcus-de-oliveira/rational-stress/check_distinct.py

# Regenerate the bounded main screen in a separate output folder.
python3 -I -B research/marcus-de-oliveira/rational-stress/screen.py \
  --output /tmp/mdo-rational-screen

# Optional full n7 regeneration. This intentionally replaces local holdout outputs.
python3 -I -B research/marcus-de-oliveira/rational-stress/holdout.py

# Verify committed artifact bytes in this directory (LF checkout).
cd research/marcus-de-oliveira/rational-stress
shasum -a 256 -c SHA256SUMS
```

Timing/platform fields vary on rerun. The exact case records, deterministic input hashes, classifications, and certificates should agree. For a main selected input, the existing `certificates/verify.py` CLI also verifies its JSON directly. For n7 use the explicitly scoped holdout/replay harness; the existing CLI correctly rejects that dimension.

## What remains open

No point outside the original hull was found, so there is no proposed MDO counterexample. These samples do not cover arbitrary spectra, unitaries, dimensions, or neighborhoods. The near-vertex schedules are structured, not an adversarial optimizer over all unitary directions. The negative exotic functional obstructs a stronger common-certificate approach, not the original conjecture. No novelty claim, universal theorem, formal verification, or resolution of the recent dimension-four literature claims follows from this run.

The immediate handoff is independent-owner replay of the immutable artifacts, especially the n7 extension and the thin polygons. A useful later search could target the exact coefficient geometry at the exotic matrices or introduce an explicitly budgeted optimization over rational parameters, while keeping candidate generation separate from exact certification.
