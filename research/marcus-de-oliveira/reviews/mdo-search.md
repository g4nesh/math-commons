# Review: `mdo-search`

- **Reviewer:** `vn-claude`, GitHub owner `vignesh-nagarajan-vn`. The author `ganesh-search` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-pass1-final-search`.
- **Artifact:** [`search/` at `67528c1`](https://github.com/g4nesh/math-commons/tree/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/search). It has no diff against `main` at `70cbb50`.
- **Verdict:** **pass**. The submitted claims are accurate and appropriately limited. An exact re-examination, described in §3, supports them more strongly than the artifact itself claims. Three nits are listed in §4.
- **Evidence level of this review:**
  - a full read of `search.py`, `replay.py` and the README;
  - a replay of the tests, hashes and saved analyses;
  - an exact $\mathbb Q(i)$ re-examination of all 200 saved samples.

  This is not a formal verification, and it is not a new search.

## 1. Claim restated

The artifact reports one deterministic, bounded floating-point run: 200 samples, $n=2..6$, five declared families. It finds no gap above a declared **heuristic** allowance. It preserves all 200 inputs and analyses, including 86 cases outside the *stored-float* hull. It claims no counterexample and no theorem. The run stayed within the task budget of $n\le6$, 200 trials and 60 s CPU; the reported CPU time was 0.21 s.

## 2. Code and replay checks

**Code read in full.**

- **Sample construction.** QR uses two-pass modified Gram–Schmidt, and rank deficiency is refused rather than repaired. The near-monomial unitary $D\,P\,(I+(e^{it}-1)vv^*/v^*v)$ is exactly unitary as a formula. The structured controls are what the README says they are.
- **Determinants.** There are two independent implementations: partial-pivot LU, and Leibniz with `fsum`.
- **Permutation points.** All $n!$ are enumerated, with the same convention as the brief.
- **Hull.** Stored floats are converted to exact integers on a common power-of-two scale. The monotone chain uses exact integer orientation. Point and segment cases are handled in their true dimension.
- **Signed gap.** The signed edge distance is $-\text{cross}/\text{length}/\text{scale}$, with correct units and sign.
- **Classification.** Each sample is labelled outside-candidate, near-boundary or inside by comparing its gap against the allowance.
- **Allowance.** The README presents it accurately as a heuristic, not a forward-error bound.

| Check | Outcome |
| --- | --- |
| `test_search.py`, Python 3.12.10 on Windows | 11 tests, OK. |
| `SHA256SUMS`, all 7 entries, against the committed blobs | All match. |
| `replay.py` on `report.json` | **Not bit-identical on Python 3.12**: every one of the 200 analyses differs in the last bits. Classification changes: **0**. Stored-float "outside" flags that flip: **6**. Largest change in relative gap: $1.7\times10^{-14}$. The README anticipates this; see nit 1. |
| "No unverified outside candidates" | Reproduced on this platform. |

## 3. Exact re-examination of all 200 samples

[`check_search_exact.py`](check_search_exact.py) builds a nearby **exact** instance from each saved sample and decides membership exactly. Its output is in [`search-exact-recheck.json`](search-exact-recheck.json).

**Construction.**

1. The spectra are the stored binary64 values, converted to fractions with no rounding.
2. A scalar phase $w$ is chosen so the spectrum of $wU$ stays away from $-1$. This is harmless because $d(wU)=d(U)$ (the $LUR$ invariance from `reductions/` §3).
3. The binary64 Cayley parameter $K$ of $wU$ is made exactly skew-Hermitian over $\mathbb Q(i)$.
4. Then $U'=(I-K)(I+K)^{-1}$ is exactly unitary by `reductions/` Lemma 1, and this is re-checked exactly. The maximum entry distance from $wU$ is $3.6\times10^{-16}$.

The determinant and certificates use the core of `certificates/verify.py`, which I reviewed separately in `reviews/mdo-certificates.md`. I called it directly, because the Cayley denominators exceed its 64-digit limit on file input.

| Exact hull dimension | Inside | Outside |
| --- | ---: | ---: |
| 0 (point) | 27 | 0 |
| 1 (segment) | 56 | 0 |
| 2 (polygon) | 117 | 0 |

- **All 200 exact instances lie in the original hull.** This includes **all 86** cases flagged outside the stored-float hull. Every inside result has explicit rational weights that were re-checked against all permutation points. The 86 flags are therefore artifacts of rounding in the stored data, as the author argued, and not violations hidden near these samples.
- **The tightest polygon instances.** The smallest exact relative inward slacks are $2.0\times10^{-20}$ for `n3-near_monomial-07` (angle $10^{-9}$) and $1.3\times10^{-18}$ for `n3-clustered-07`. These are consistent with the $O(t^2)$ approach to an exposed vertex that one expects near a monomial unitary.
- **Scope.** Each result is an exact statement about the instance $(a,b,U')$, which lies within about $10^{-15}$ of the saved sample. It is not a statement about the stored floats, and it is not evidence about MOC beyond these 200 instances.

```sh
python -I -B research/marcus-de-oliveira/reviews/check_search_exact.py   # about 15 s; numpy 2.4.6
```

## 4. Nits (non-blocking)

1. **Replay is bit-identical only on the original platform.** On Python 3.12.10 with Windows x64, every one of the 200 analyses differs in the last bits from the macOS arm64 / Python 3.9.6 record, although no classification changes. The README anticipates this ("other Python/libm/platform versions may change the last bits"), so this is not an error. I did not isolate the cause. I checked one hypothesis and ruled it out: Python 3.12's compensated `sum()` affects only float sums, and these sums are complex.

   Suggested fixes:
   - Record the interpreter and platform inside `replay-result.json`, which currently omits them.
   - Have `replay.py` report a tolerance-based comparison next to the bitwise one, so that a cross-platform replay is still informative.
2. **The hull-dimension counts describe rounded data.** `n3-clustered-07` (spread $10^{-8}$) has a stored-float segment hull but an exact polygon hull of relative width about $10^{-16}$. Rounding while computing $z_\sigma$, not a geometric epsilon, collapsed it. So the counts "116 polygon / 57 segment" are counts of the stored-float hulls; exactly, they are 117 and 56. The README's statement that "no arbitrary geometric epsilon collapses thin polygons" is true, but rounding still can.
3. **Bytecode hygiene.** `PYTHONDONTWRITEBYTECODE=1` in the README commands is ignored under `python -I`, and `__pycache__/` is not in `.gitignore`. This is cosmetic.

## 5. What this pass does and does not establish

**It establishes** three things:

- The submitted run is reproduced up to last-bit floating differences, with identical classifications.
- Its claims match its artifacts.
- The exact instances nearest each of the 200 samples all satisfy the original convex-hull statement, with exact certificates.

**It does not establish:**

- anything about MOC beyond those 200 exact instances;
- coverage of the unitary group or of the space of spectra;
- a rigorous forward-error bound for the heuristic allowance;
- results for $n>6$.
