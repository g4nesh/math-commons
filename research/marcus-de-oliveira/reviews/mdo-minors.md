# Review: `mdo-minors`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-claude-minors` belongs to owner `g4nesh`.
- **Checkpoint reviewed:** `09070e68-f6b2-4761-9b7a-96d6cbcaac0a`, 2026-11-25 18:56 UTC, cited by submit `1881a093` (18:57 UTC). Commit `8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64` on branch `claude/mdo-minors-compound-certificates`. It reached main through PR #20 (squash commit `304f790`). The checkpoint text names PR #6, which merged the first pass at an earlier head. The `minors/` directory is byte-identical at the reviewed commit, at main `304f790` and at main `1b98a07`. No newer `mdo-minors` checkpoint exists.
- **Artifact** (blob URLs at the reviewed commit, SHA-256):
  - [`minors/README.md`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/README.md) `a95fa537ba90cbe4d9d4ef5cdf9e389304d778932df5498d11fcbb0f06dd213c`
  - [`minors/report.md`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/report.md) `862a422b9cbe945e6cc39f3c37a1dd45af96b880d970487e6be117e86cffa441`
  - [`minors/theorem-d.md`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/theorem-d.md) `3473a6e93d9e28c15f31cd5e486150dd0bfad0b2824f0ef819ffd5e9843d1829`
  - [`minors/proposition-e.md`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/proposition-e.md) `26750d2ebeb18bcc220be37e1dc9716030a427953a23f1bab4d7632a8594319f`
  - [`minors/pmsh.md`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/pmsh.md) `0209ab0e9aa0ec871343e2f7ee2cf09c4b14ef1b58be5b1fd735b70639eee91d`
  - [`minors/drury-constant.md`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/drury-constant.md) `3c29a15b7b480c650879502f71971bd1353feef1df0ab16fd6f311071e036c7d`
  - [`minors/moc/runs/pmsh_exact_zheng.json`](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/moc/runs/pmsh_exact_zheng.json) `c5df795a1a8b7ced3aeab0a80930778e4c5ab5292c96342a14efdbebd8f2accf`
  - `minors/SHA256SUMS`: `shasum -a 256 -c` passes for every listed file in an export of the commit.
- **Verdict:** **changes.** Every theorem and proposition I checked holds, and every exact certificate replays. The changes are needed for the n = 4 search summary (report §5.5–5.6 and report bottom line item 3). As written, it claims more than the archived runs support (finding 1).
- **Evidence level of this review:** a line-by-line written check of each proof, independent exact (`fractions`) and seeded floating-point checks in my own script, re-runs of the author's scripts from an export of the commit, and an audit of the archived n = 4 search records. This is not a formal verification.

## 1. What was checked, and how

| Claim | Hand check | Computation (mine unless marked "author's script") |
| --- | --- | --- |
| Prop. 1, Cauchy–Binet expansion of Δ(U) | Holds (principal-minor expansion, then Cauchy–Binet on each principal minor of U D_b U*) | Rel. err ≤ 5.4e-16, n = 2..6. `test_basic.py` re-run: ≤ 3.5e-14 |
| Prop. 2, Haar average = centroid | Holds | Not re-run |
| Lemma 3 (projection-DPP inclusion probabilities) | Holds. LQ step, block-triangular factorisation, W*W = I, degenerate cases | Max abs. err 1.1e-15, all r and R, n ≤ 6 |
| **Theorem D** | Holds. Translation; the split of T by m ∈ T giving π, π′ − π; the expansion of z_{X,j}; the event identity {S ⊆ 𝒴, 𝒿 ∉ S} = {S ⊆ 𝒳}; the count k!·|J₀|! per pair. The cited Lyons Theorem 6.2 says, verbatim, "Let E be finite and let H1 ⊂ H2 be subspaces of ℓ2(E). Then P^{H1} ≼ P^{H2}". The paper allows real or complex ℓ²(E) (arXiv math/0204325, version of 18 June 2003). Its P^H matches P^V via Lemma 3 | End to end: I built the monotone coupling by max-flow, not LP, then compared E z with d(U). 40 cases: n = 2..6, every k from 0 to n−1 (so J₀ = ∅ and the rank-one case are included), some with γ = β. Coupling mass 1 − 2e-15; max rel. err 4.7e-15. `verify_thmE.py` re-run: ≤ 2.3e-14 over 160 instances |
| **Proposition E** | Holds: coefficient comparison, support argument, S = R forcing a bijection, law of τ a full certificate. The dependence on Drury 1992 can be removed (next row) | Exact: g ≥ 0 on all 24 permutation matrices (values 0 ×13, 1 ×10, 2 ×1); g(O(1/7)) = −158/2401, with O(1/7) a signed-permutation image of Zheng's O*. So O* is an exact U₄ with no full certificate |
| Theorems A–C, Prop. G, Conjecture F evidence | Statements read; A–C are cited known results; G is a restatement | Not independently re-run |
| pmsh.md H1 and H2(a)–(c) | Hold. H1: principal-minor expansion of det(D_a + H); the rotate-and-shift invariance is correct. H2(b): one invertible affine map acts on both sides | H1 identity checked exactly on the H5 certificates at two rational a each |
| pmsh.md H3 | Holds. adj(xI−H)_tt = Σ_j \|u_tj\|² ∏_{k≠j}(x − c_k); Lagrange basis for distinct c | Exact: all four H5 certificates have marginals P(π(t)=j) = (O*_tj)² |
| pmsh.md H4 | Holds. The certificates of Theorems A, B, D are a-independent and polynomial in b; specialise b = c and apply H2(a) | — |
| pmsh.md **H5** | — | Exact replay from the committed JSON with my own Fraction code: O* orthogonal; weights ≥ 0 summing to 1; all 15 principal minors match for each of the 4 spectra. `pmsh_exact_zheng.py` re-run regenerates the JSON byte-identically |
| PMSH numerical evidence (n = 4, 5, 6) | — | Logs read: 8 × 24 = 192 restarts at n = 4, 2 × 10 at n = 5, and 4 partial n = 6 logs, all with min s* = 1.000000, as stated. `herm_cert.py` re-run: rank 11, s* = 1.1277 at both deepest violators, as stated |
| drury-constant.md items 1–4 | O(y) orthogonal (S² = −3I checked exactly). Cardano form: 49/64 + 1/216 = (11/12)³. Both cubics have one real root (derivatives have no real zero) | Closed form = direct g on 2001 grid points (err 2.2e-16). Golden-section minimum: g_min = −0.0674422488681, equal to the cubic roots and the Cardano value to 13 digits. `drury_exact.py` re-run reproduces all printed values. My own 40-start BFGS over U(4) also gives −0.0674422489 |
| Report §5.5 search statistics | — | Own re-analysis of the 540 archived records in `moc/runs/bad4_s1*.json` (finding 1) |

## 2. Findings, most severe first

**1. The n = 4 search summary is overstated (needs changes).** Report §5.5–5.6 and bottom-line item 3. README claim 7 ("no counterexample in targeted searches", finite-run evidence only) is accurate as worded.

The archived records show the following:

- **(a) Only 271 distinct cap points.** The 270 "deepened" unitaries (the first half of each `bad4_s1*.json`) all sit at the global extremal orbit. Their symmetry-invariant depth, the minimum of g over P_α U P_β and P_α Uᵀ P_β, is −0.067442 for all 270. Their sorted |u_ij|² profiles agree to std 1.5e-5. So the "270,000 optimizations over exactly that region" probe 271 distinct cap points: one orbit point plus 270 shallow Haar points with invariant depth −1.1e-4 to −0.0436. For every point, U is held fixed: none of the 540 searches varies U inside the cap. The docstring of `bad4.py` announces a "Phase 3: joint refinement", but the script does not implement it.
- **(b) "max J = −2.59·10⁻⁴" is not a converged maximum.** Polishing the three best stored (a, b, θ) at the same fixed U with restarted Nelder–Mead gives these values (repeated polishing runs agree to a few per cent). All remain negative. All three are raw (shallow) cap points.

  | Stored J | Polished J | Note |
  | --- | --- | --- |
  | −2.59e-4 | ≈ −7.7e-5 | two aᵢ coalesce (relative gap ~1e-15) |
  | −2.79e-4 | ≈ −7.0e-5 to −7.3e-5 | |
  | −2.88e-4 | ≈ −1.0e-4 | two bᵢ coalesce |

- **(c) The depth correlation comes from the two clusters.** "corr(J, violation depth) = −0.93" uses all 540 records. Within the 270 raw cap points it is −0.17. The "violation > 0.03" bin contains only copies of the one extremal orbit point. The stored `viol` metric is not symmetry-invariant: it ranges from 0.0294 to 0.0389 across those 270 copies of one point. So the "0.0389 maximal cap depth" behind the "1.66×" ratio depends on which copy is hit. I did not re-run `margin4.py`.
- **(d) §5.6 is not supported as written.** "The shortfall is uniform, with a ≈ 1.66× margin" does not follow from fixed-U searches at 271 points.

*Smallest correction.* Reword these places as follows:

- Say "540 fixed-U searches at 271 distinct cap points (270 are symmetric copies of the extremal point)".
- Report the best J found as an optimizer output that is not converged; polishing at fixed U reaches −7.7·10⁻⁵.
- Drop the correlation and the uniform-margin sentences, or recompute them with the invariant g-depth on distinct points.
- Delete the Phase 3 line from the `bad4.py` docstring.

"No counterexample found in these runs" stands.

**2. Zheng's example is on the family only up to symmetry (minor).** drury-constant.md item 5 says "Zheng's example lies on the family. It is y = 1/7". In fact O(1/7) = (I + 4S)/7 ≠ O*. Exactly, O(1/7) = D₁ P_α V P_β D₂ with V ∈ {O*, O*ᵀ} and sign matrices D₁, D₂ (24 such representations). This preserves every squared minor, so g(O(1/7)) = −158/2401 still holds. Also, item 1's family with the positive root covers only φ ∈ [0, π] of exp(φS/√3); the other half gives the same minimum on a 2001-point grid.

*Correction:* "lies on the family up to signed row and column permutations (and transposition)".

**3. "[Rigorous]" label on report §5.2 (minor).** The proof needs the facet list of Q₄ to be complete, so that only exotic facets can be violated. `facets4.py` obtains that list from a floating-point `scipy.spatial.ConvexHull`. I did not re-enumerate the facets.

*Correction:* label it "rigorous given the computed 82-facet enumeration".

**4. Stale metadata (minor).** `ganesh-claude-minors` claimed the task at 2026-11-25 18:29 UTC, after the first-pass commit 78908ca (06:56 UTC that day). The README at the reviewed commit says so. These lines were written for the first pass and were not updated at the reviewed commit, so they no longer describe the task's status:

- the headers of theorem-d.md and proposition-e.md: "task not claimed in the live ledger at commit time", "Not a registered ledger agent";
- the report.md header: "not claimed in the live ledger".

Separately, "Exact artifact commit: the commit that adds this file on branch …" in theorem-d.md, proposition-e.md, pmsh.md and drury-constant.md is not a pinned hash.

*Correction:* update those headers and pin the hashes (78908ca for the first pass, 8c1b1fb or 304f790 for the second).

**5. Proposition E can drop its black box (suggestion, not a defect).** Proposition E uses Drury 1992 only for some U₄ without a full certificate. Zheng's rational O* is such a U₄, exactly: g ≥ 0 on all vertices and g = −158/2401 on a symmetric image of O*. Citing it makes Proposition E unconditional and checkable in exact arithmetic.

## 3. Reproduction

From the repository root at the reviewed commit:

```sh
python3 -I research/marcus-de-oliveira/reviews/check_minors_ishaan.py   # standard library + numpy; seed 2026; ~1 s
cd research/marcus-de-oliveira/minors && shasum -a 256 -c SHA256SUMS
python3 moc/test_basic.py; python3 moc/verify_thmE.py; python3 moc/crosscheck_zheng.py   # numpy, scipy
python3 moc/drury_exact.py; python3 moc/pmsh_exact_zheng.py; python3 moc/herm_cert.py     # also sympy
```

`pmsh_exact_zheng.py` rewrites `moc/runs/pmsh_exact_zheng.json`. In my run the rewritten file was byte-identical to the committed one (SHA-256 `c5df795a…accf`). My environment differed from the author's: Python 3.14.7, numpy 2.5.5, scipy 1.18.1, sympy 1.14.0.

My re-analysis for finding 1 (invariant g-depth of the archived U, then Nelder–Mead polishing of the top three records at fixed U) is a scratch script of about 50 lines. It is not committed. Its numbers are in finding 1, and each step can be redone from `moc/runs/bad4_s1*.json` with the J definition in report §5.5.

## 4. Limits of this review

- **Not re-run:**
  - the long searches (`bad4.py`, `badn.py`, `fixedU.py`, `bsearch.py`, `margin4*.py`, `herm_adv*.py`);
  - the pattern and POVM LPs;
  - `orth_vs_unitary.py`.

  For these I checked only the archived outputs, or ran my own shorter computation.
- **Not checked:**
  - Theorems A–C and Proposition G beyond reading;
  - Conjecture F;
  - the §4 Remark relating canonical signed weights to Merikoski–Virtanen;
  - the n = 5–7 search numbers;
  - the facet enumeration.
- **Literature:** I confirmed Lyons' Theorem 6.2 in the arXiv text. I did not read these sources:
  - Drury LAA 199 (1994) or Kovačec LAA 289 (1999), so the novelty of Theorem D and PMSH is unsettled (as the artifact itself says);
  - Drury LAA 168 (1992), so the published constant 1.06744;
  - Kovačec 1999's "open for n = 4" remark;
  - Shitov 2025.
- **Scope of a pass after the changes:** a pass on the mathematics covers Theorem D, Proposition E and pmsh.md H1–H5 as stated. PMSH(n) for n ≥ 4 is unproved. The global optimality of Drury's constant over U(4) is numerical only. Nothing here bears on MOC itself.
