# Review: `mdo-rational-stress`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-search` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-rational-stress-final-20261125`, 2026-11-25 18:45:25 UTC, the checkpoint cited by submit `mdo-rational-stress-submit-20261125` (18:45:31 UTC). The task has no later work.
- **Artifact:** `rational-stress/` at commit `f425dbc` (PR #15 branch). PR #15 was squash-merged, so `f425dbc` is not an ancestor of `main`, but the directory is byte-identical on `main` at `1b98a07`.
  - [`README.md`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/README.md), SHA-256 `2a71a927c3d8297643e2b58ea685c72d56b42f168d0dc20f723d1d3b91e0dfc8`
  - [`SHA256SUMS`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/SHA256SUMS), SHA-256 `7a41e92eb83148f61691166001c36e21f0af03dea637cb99019d2049ce52acf9`. It pins 66 files: the README and the other 65 files in the directory, every file except `SHA256SUMS` itself. Among them: [`screen.py`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/screen.py) `b39b246f…bfaa`, [`replay.py`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/replay.py) `a08ed4d3…859a`, [`holdout.py`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/holdout.py) `5029a478…f162`, [`report.json`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/report.json) `187e1f1f…cb5a` and [`holdout/report.json`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/rational-stress/holdout/report.json) `27b84f62…c2a7`.
  - The dependency [`certificates/verify.py`](https://github.com/g4nesh/math-commons/blob/f425dbcdcce6bff43101921be1738bedf4c9d97d/research/marcus-de-oliveira/certificates/verify.py) has SHA-256 `8cb1a728…78ff`. It has no diff from `67528c1`, the version that `mdo-certificates` passed review at.
- **Verdict:** **pass**. Four minor wording notes are listed in §2. None of them changes a reported number or classification.
- **Evidence level of this review:**
  - a full read of the README, `screen.py`, `replay.py`, `holdout.py`, `check_distinct.py`, `test_screen.py` and `certificates/verify.py`;
  - a re-run of the author's tests, checksums and both replays, plus a full regeneration of the main screen;
  - an exact re-check of all 924 cases with my own code, which shares only the input generator with the artifact.

  This is not a formal verification.

## 1. What was checked, and how

All runs used an export of `f425dbc` made with `git archive`, under Python 3.11.7 on macOS arm64. The author recorded Python 3.9.6.

| Claim | Check | Result |
| --- | --- | --- |
| Artifact integrity | `shasum -a 256 -c SHA256SUMS` | 66/66 OK |
| Six focused tests pass | `unittest discover -p "test_*.py" -v`, run both inside the directory and with the README's `-s` form | 6/6 OK, about 4 to 6 s |
| Main and n7 replays show 0 mismatches | `replay.py` on `report.json` and on `holdout/report.json` | 900 and 24 cases, 0 mismatches. The report hashes match. |
| The main screen is deterministic | `screen.py --output <scratch>`, then a comparison with the committed files | `records`, `summary`, `selectedWitnesses` and `sources` are equal. All 30 witness files are byte-identical. Stop reason `case-limit-complete`. About 30 to 40 s CPU. |
| The counts in the results table (762 / 12 / 150 / 0 per dimension), and 302,400 permutation points | Own checker (below) | Exact match per dimension and class. Points: 181,440 + 120,960. |
| Exact unitarity and normality, noncommutation, a complex $U$ | Own checker: $UU^*=U^*U=I$, $BB^*=B^*B$, and $\lVert AB-BA\rVert_F^2$ from full products | All 924 hold. Every commutator norm is positive and equals the recorded value. |
| Exact determinant | Own Leibniz sum over Gaussian integers after clearing denominators. The verifier uses elimination and the replay uses subset DP. | All 924 equal the recorded values. |
| Inside certificates | Own permutation products. Weights are nonnegative, sum to 1, and reconstruct $\det(A+B)$ exactly. Support size is at most 3. | All 924 pass |
| Hulls, interior/boundary/segment, the margin $m$ | Own gift-wrapping hull over integer-scaled points, compared as a set with the recorded hull. Own margin formula. | Hulls, classes and every exact margin agree |
| $\psi(O_*)=-158/2401$. $\psi\ge0$ on the 24 permutation matrices. 100 negative-$\psi$ inputs at $n=4$, 60 of them general. | Own $\psi$, from the definition in `minors/moc/crosscheck_zheng.py` | Holds. 100 and 60. |
| 924 distinct ordered pairs $(A,B)$ | Own hash of $(a, B)$ | 924 distinct |
| The 54 saved witness files equal the regenerated inputs | Own comparison | 54/54 |
| Description of the construction (gates, phases, families, singular embedding, $n=4$ singular without an $O_*$ block) | Read against `screen.py` | Accurate |
| The margin is an $\ell^\infty$ distance | Derivation: for the line through $p,q$ with normal $\nu=(-\Delta y,\Delta x)$, the $\ell^\infty$ distance is $\lvert\nu\cdot(d-p)\rvert/\lVert\nu\rVert_1=\operatorname{cross}(p,q,d)/(\lvert\Delta x\rvert+\lvert\Delta y\rvert)$, because $\ell^1$ is the dual norm | Accurate |
| "Tight" table: the smallest general-complex margin per dimension | Sorted the report records | All four case IDs and values are correct |

I tested the checker against tampered records. It rejects:

- a determinant shifted by $10^{-30}i$;
- a dropped hull vertex;
- a doubled margin;
- a swapped support permutation.

## 2. Findings (all minor, non-blocking)

1. **The "general complex" spectra have fixed real parts.** In `generic`, $\operatorname{Re}a_i=2i-n$ and $\operatorname{Re}b_i=3i-n+1$ for $i=0,\dots,n-1$ (near-vertex mixing only reorders $b$). In `thin-complex` they are $2i+1$ and $3i+2$. Only the imaginary parts are random: small integers, or $\varepsilon$ times small integers. `clustered` is a random shift of the same pattern, scaled by $\varepsilon$. My checker finds one real-part configuration per dimension in each of the two families, so 4 each over $n=4,\dots,7$.

   The README's closing caveat ("do not cover arbitrary spectra") is true. But "general complex stress cases" and "generic complex" read broader than this.

   *Correction:* add one sentence under "Construction and coverage" that states the fixed real parts.
2. **The tight margins are small by construction.** Each of the four tightest general-complex cases comes from `thin-complex` with `near-vertex` mixing. Over all 54 such cases, $m/(t^2\varepsilon)$ lies in $[0.55, 31.9]$, where $t$ is the gate parameter and $\varepsilon$ is the thinness.

   This scaling follows from a second-order expansion, which I derived as follows. Take $U=e^{X}$ with $X=tK+O(t^2)$ skew-Hermitian, and $D=D_a+D_b$ invertible.
   - The first variation is $\det D\cdot\operatorname{tr}(D^{-1}[X,D_b])$. It vanishes, because $[X,D_b]$ has zero diagonal.
   - The second-order term is

     $$\det(A+B)=z_{\mathrm{id}}+t^2\sum_{i<j}\lvert K_{ij}\rvert^2\,(z_{(ij)}-z_{\mathrm{id}})+O(t^3),$$

     using $z_{(ij)}-z_{\mathrm{id}}=\det D\,(a_i-a_j)(b_i-b_j)/(d_id_j)$.

   A float check (`numpy`, $n=5$, random complex spectra and a random skew-Hermitian $K$) gives a relative error of the $t^2$ term that falls roughly in proportion to $t$ over $t=10^{-1},10^{-2},10^{-3}$, as an $O(t^3)$ remainder should; the size of the error depends on the random draw. So near-vertex determinants move into the cone that the transposition points span at the vertex. Their tiny margins measure $t^2$ times the thinness of the polygon, not closeness to a violation.

   The README says the margin is "not a robustness radius" and that near-vertex schedules are "not an adversarial optimizer", which is consistent with this. Still, the handoff to check "especially ... the thin polygons" could be read as saying these are the hardest cases.

   *Correction:* one sentence saying that near-vertex margins scale like $t^2\varepsilon$ by construction.
3. **The boundary and segment cases are all controls.** All 12 polygon-boundary cases are `singular` controls, where $\det=0$ is the identity permutation point. All 150 segment hulls are the `hermitian` controls. All 474 general-complex cases are strict interior. The results table does not show this.

   *Correction:* add a family column or one sentence.
4. **Attribution of $O_*$.** The README calls $O_*$ "Zheng's rational orthogonal matrix" but links only to the internal transcription. I checked exactly that $O_*$ is orthogonal and that $\psi(O_*)=-158/2401$. I did not retrieve the preprint.

   *Correction:* cite the source directly, as `literature/README.md` does (Zheng, *The Marcus–de Oliveira conjecture for n = 4*, October 2026 preprint, Example 1.3).

**Not a defect.** The 24-case $n=7$ holdout goes beyond the task's declared 900-case target. It is announced in a ledger checkpoint, budgeted separately, labeled separately in the README, and leaves the main report unchanged: the hash of `report.json` matches `replay.json`.

## 3. Reproduction

From the repository root, at `f425dbc` or `main`, with this review's script added:

```sh
(cd research/marcus-de-oliveira/rational-stress && shasum -a 256 -c SHA256SUMS)       # 66 OK
python3 -I -B -m unittest discover -s research/marcus-de-oliveira/rational-stress -p "test_*.py" -v   # 6 OK
python3 -I -B research/marcus-de-oliveira/rational-stress/replay.py \
  research/marcus-de-oliveira/rational-stress/report.json --output /tmp/rs-replay.json        # 900, 0 mismatches
python3 -I -B research/marcus-de-oliveira/rational-stress/replay.py \
  research/marcus-de-oliveira/rational-stress/holdout/report.json --output /tmp/rs-n7.json    # 24, 0 mismatches
python3 -I -B research/marcus-de-oliveira/rational-stress/screen.py --output /tmp/rs-regen   # ~30-40 s CPU
python3 -I -B -c "import json; a,b=(json.load(open(p)) for p in ('research/marcus-de-oliveira/rational-stress/report.json','/tmp/rs-regen/report.json')); print(a['records']==b['records'], a['summary']==b['summary'])"   # True True
python3 -I -B research/marcus-de-oliveira/reviews/check_rational_stress_ishaan.py   # ~30 s CPU; prints "all exact checks passed"
```

`check_rational_stress_ishaan.py` uses only the standard library. Its expected output includes:

- `cases` 924;
- `generalComplex` 474 and `controls` 450;
- `permutationPoints` 302400;
- `distinctExactPairs` 924;
- `savedWitnessFilesMatched` 54;
- `n4NegativePsi` 100, of which `n4NegativePsiGeneral` 60;
- per-dimension classes that equal the README table.

## 4. Limits of this review

- **The input generator is shared.** `report.json` stores seeds, not matrices, so my checker calls `screen.make_witness` to rebuild inputs. I did not re-implement the generator. I only checked that its description in the README matches the code and that the 54 saved files match it. Everything downstream of the input is checked independently.
- **Not checked:**
  - the author's recorded timings, and the external CPU timers said to bound the replays;
  - the 30-minute session bound, beyond ledger timestamps: claim at 18:26:50 UTC, final checkpoint at 18:45:25 UTC.
- **Scope.** The pass covers 924 explicit Gaussian-rational instances. It says nothing about MDO in general, about other spectra or unitaries, or about Zheng's claimed $n=4$ proof. The sample is structured, not adversarial: as §2 notes, real parts are fixed and the near-vertex cases are inside to second order. It is not evidence that violations are absent near these points.
- The second-order expansion in §2 is my own short derivation, checked numerically in floating point. It is used only to interpret margins, not to support any classification.
