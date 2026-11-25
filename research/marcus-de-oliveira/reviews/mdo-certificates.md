# Review: `mdo-certificates`

- **Reviewer:** `vn-claude`, GitHub owner `vignesh-nagarajan-vn`. The author `ganesh-certificates` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `mdo-pass1-final-certificates`, submitted 2026-11-25.
- **Artifact:** [`certificates/` at `67528c1`](https://github.com/g4nesh/math-commons/tree/67528c1f93ba64bf474c8a5954f8c029e58ae268/research/marcus-de-oliveira/certificates). It has no diff against `main` at `70cbb50`.
- **Verdict:** **pass**. No correctness defect was found. Two robustness nits are listed in §3 and one scope note in §4.
- **Evidence level of this review:** a full read of `verify.py` and the README, a replay of the author's tests and results, and differential testing against oracle code written independently. This is not a formal verification of the program.

## 1. Claim restated

For one input $(a,b,U)$ with Gaussian-rational entries and $n\le6$, the program does five things:

1. It checks $UU^*=U^*U=I$ exactly.
2. It forms $A=\operatorname{diag}(a)$ and $B=U\operatorname{diag}(b)U^*$, then computes $d=\det(A+B)$ exactly.
3. It enumerates all $n!$ points $z_\sigma=\prod_i(a_i+b_{\sigma(i)})$.
4. It returns either
   - **inside**: at most three nonnegative rational weights on permutation points that sum to 1 and reconstruct $d$ exactly, or
   - **outside**: a rational functional whose value at $d$ exceeds its value at every permutation point by a strictly positive rational gap.
5. A direct check re-validates that certificate against the full point set.

This matches the convention in the objective brief.

## 2. What was checked

**Code read in full: `verify.py`, 401 lines.** I traced every branch by hand.

- **Parsing.** The regex admits only integer or `p/q` strings with a positive denominator and at most 64 digits per integer. `bool` is rejected for `schemaVersion`. Duplicate keys, `NaN` and unknown fields are rejected. Float, bool and non-`Fraction` scalars are rejected in both the generator and the direct checker.
- **Unitarity.** Both products are compared with $I$ exactly. Nothing is rounded or repaired.
- **Determinant.** Gaussian elimination over $\mathbb Q(i)$ uses the first nonzero pivot and tracks the sign through row swaps. A missing pivot returns $0$, so singular sums are handled without division by zero.
- **Hull geometry.**
  - The monotone chain removes collinear points (`<= 0`) and returns counterclockwise vertices.
  - **Point hull:** equality test, or separator $t-h$.
  - **Segment hull:** a nonzero orientation gives a separator by the left or right normal. A collinear target gives a parameter, and the parameter lies in $[0,1]$ exactly when the target is inside. Otherwise the separator is $\pm$ the displacement, and its sign is correct for both overshoot directions.
  - **Polygon:** an edge with negative orientation gives the outward normal $(\Delta y,-\Delta x)$. Otherwise fan-triangle barycentric weights $\operatorname{area}(t,v_k,v_{k+1})/\operatorname{area}(v_0,v_k,v_{k+1})$ are used. Every polygon area is nonzero, because collinear points are removed.
- **Direct check.** It does not reuse the hull. It requires support points in the full point set, nonnegative weights summing to one and exact reconstruction. For outside certificates it requires a nonzero normal, a bound that holds at every point, and a strictly positive gap that equals its claimed value.
- **Hypotheses.** Normality and the prescribed spectra follow from exact unitarity, as the README says. The extra normality checks are redundant but harmless.

**Replays.**

| Check | Outcome |
| --- | --- |
| `python -I -m unittest discover …`, Python 3.12.10 on Windows | 16 tests, OK. |
| `SHA256SUMS`, all 12 entries, against the committed blobs | All match. |
| Rerun of `verify.py` on all three fixtures | The output equals the committed `results/*.json` byte for byte once line endings are normalized (see nit 2). |
| The two printed weight identities in the README | Confirmed by hand. They sum to 1 and give $2/3+2i/9$ and $(2864+1840i)/49$. |
| Determinants of the Goyal and Zhang–Ma inputs, and their hull membership | Agree with my separate exact code: `adversarial/polygon_enclosure.py` on PR #11. |

**Differential and adversarial tests.** These are in [`check_certificates_adversarial.py`](check_certificates_adversarial.py), with seed 20261125. The oracle decides membership by brute force over Carathéodory supports (points, segments, closed triangles) and builds no hull. It computes determinants by Leibniz expansion.

- **A. Geometry: 3000 cases**, biased toward degeneracies: single and repeated points, collinear sets, and targets on vertices, edges and chords. The results were 1554 inside and 1446 outside, with **0 mismatches** against the oracle. Every returned certificate also passed my own re-check.
- **B. End to end: 150 exact witnesses.** These used rational Cayley unitaries, $n=2..4$, and Gaussian-integer spectra, including repeated eigenvalues and 3 singular sums. There were **0 determinant mismatches** against Leibniz and **0 membership mismatches** against the oracle. All 150 were inside. This bounded run is **not** evidence about MOC at scale.
- **C. Malformed files: 14 of 14 rejected** with exit 2 and nothing on stdout. The files covered:
  - float, decimal, exponent and padded strings;
  - zero and negative denominators;
  - a 65-digit integer;
  - a bool `schemaVersion`;
  - an unknown field and a duplicate key;
  - `NaN`;
  - a unitary off by $10^{-30}$;
  - an expected-determinant mismatch.

```sh
python -I research/marcus-de-oliveira/reviews/check_certificates_adversarial.py   # about 5 s
```

## 3. Nits (non-blocking)

1. **Uncaught recursion on hostile JSON.** A file of 100,000 nested brackets, well under the 1 MiB limit, makes `json.loads` raise `RecursionError`. `main()` does not catch it, so the program exits with status 1 and a traceback, not the documented status 2. It emits no certificate, so this is not a soundness issue. Fix: catch `RecursionError`, or `ValueError` generally, in `load_witness`.
2. **Platform-dependent replay bytes.** `--output` uses `Path.write_text` without `newline="\n"`, so Windows writes CRLF. Separately, with `core.autocrlf=true` a Windows checkout fails `sha256sum -c SHA256SUMS` on every file, even though the committed blobs are correct. Fix: add a `.gitattributes` rule such as `research/**/certificates/** text eol=lf`, or note that the hash check must run on committed blobs.

## 4. Scope note

The task description asked for "exact **or** rigorous interval certificates" and asked that roundoff be covered. The implementation is exact-only over $\mathbb Q(i)$. It handles roundoff by **refusing** approximate inputs, not by enclosing them, and the README says so plainly. This satisfies the task's either/or wording. It also means this tool cannot certify a binary64 candidate from `search/`. Such a candidate needs an exact rational reconstruction first, for example a rational Cayley unitary as in `reductions/` §4, or a separate interval argument.

## 5. What this pass does and does not establish

**It establishes** that for any accepted input, an `inside` or `outside` output is backed by an exact certificate. That certificate is checked against all $n!$ permutation points and against the stated construction. The three committed fixtures are inside the original hull.

**It does not establish:**

- anything about MOC beyond those single instances;
- correctness for inputs the parser rejects;
- performance bounds;
- the soundness of Python's `fractions` module or interpreter, which form the trusted base.
