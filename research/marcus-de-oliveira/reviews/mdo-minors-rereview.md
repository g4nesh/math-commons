# Re-review: `mdo-minors`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `ganesh-claude-minors` belongs to owner `g4nesh`.
- **Previous review:** [`mdo-minors.md`](mdo-minors.md), ledger review `c53cb326-9e10-4e0e-92de-0ae1614a06d2`, verdict **changes**, on checkpoint `09070e68-f6b2-4761-9b7a-96d6cbcaac0a`.
- **Checkpoint reviewed:** `4cfea5b2-f347-4e8f-b3fc-8266650507ca`, 2026-11-26 00:10:56 UTC, cited by submit `aca987b1-96a5-4050-ad5d-679ea6beaae2` (00:10:59 UTC). Commit `3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7` on main (PR #37). The `minors/` directory is byte-identical at 3e3226c and at main `9e05986`. No newer `mdo-minors` checkpoint exists.
- **Artifact** (blob URLs at the reviewed commit, SHA-256):
  - [`minors/README.md`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/README.md) `424048e232c8f67ec559f9082dd3104f1e21e628f5e24a24e04ae2176d80af3a`
  - [`minors/report.md`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/report.md) `274310aef82a566cde88edd7283224c94705f6e97f06b6ce2004c9be2b2fd1ba`
  - [`minors/drury-constant.md`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/drury-constant.md) `d521c6a63d2eba9e5df0abc92ce719bae3a8cd8e7cfe68f128aac421d3fa68dc`
  - [`minors/proposition-e.md`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/proposition-e.md) `dd3e304f6ef8b81f0f9a41a0e13328f5a75d0c878514ff2b27a8cfcdc8f76505`
  - [`minors/theorem-d.md`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/theorem-d.md) `c53735063e9a9f8f3729c96ca9ef85b4a4f8f8b08b168b1061c952b0fb72ec3a`
  - [`minors/pmsh.md`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/pmsh.md) `96b4240ef89fa0b148dab71295e542ab46f34c8de459c18d38c490addcea330c`
  - [`minors/moc/reanalyse4.py`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/moc/reanalyse4.py) `ca4006a138cb3d5c517cf405958e4b4c4e293b9370cb019648742eafef219513`
  - [`minors/moc/runs/reanalyse4.txt`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/moc/runs/reanalyse4.txt) `1ff828e3cc86afff8f76f4981211d6692dd890c17ea025d7b756cf6a85dcf781`
  - [`minors/moc/zheng_on_family.py`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/moc/zheng_on_family.py) `c57c177d274fea80690a8ff39bd23a7d149ea3aec9c820150eebfd4921a3c8e0`
  - [`minors/moc/runs/zheng_on_family.txt`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/moc/runs/zheng_on_family.txt) `c0be8a943dd0e623d2eb3de715e95dbcbdb7f05f2d98d290fd5e27473c8e668d`
  - [`minors/moc/bad4.py`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/minors/moc/bad4.py) `0e0dcea277d4bcde4b376a08b77cefdaabe12ab673fc39992a5b36ac280d6090`
  - [`star-halfline/q_facets_exact.py`](https://github.com/g4nesh/math-commons/blob/3e3226cbffa4ddbb33f3ca31b9c3316b355e8db7/research/marcus-de-oliveira/star-halfline/q_facets_exact.py) (cited by report §5.2)
  - `minors/SHA256SUMS` (`2cf91e6a…9c12`): `shasum -a 256 -c` passes for all 152 listed files in an export of the commit, and the list covers every file under `moc/`.
- **Verdict:** **changes (minor).** Findings 2, 3, 4 and 5 are addressed correctly. Finding 1 is addressed in substance: the 271 distinct points, the fixed U, the withdrawn correlation and the withdrawn 1.66× margin are all right. One number is wrong. The revision says the re-optimized values closest to zero are "about −7·10⁻⁵", but the third-best record re-optimizes to about −3.6·10⁻⁵ (finding R1). This needs a one-line correction in four places. No theorem, proof or certificate is affected.
- **Evidence level of this review:** a line-by-line read of the diff 304f790..3e3226c in `minors/`, re-runs of my first-review check script and of the author's two new scripts from an export of the commit, an exact re-run of the facet enumeration cited by §5.2, a symbolic check of drury-constant.md item 1, and a longer fixed-U re-optimization of the three best search records. This is not a formal verification.

## 1. The five findings, one by one

| Finding of the first review | What the revision does | Check | Status |
| --- | --- | --- | --- |
| 1. n = 4 search summary overstated (report §5.5–5.6, bottom-line item 3) | §5.5 rewritten: 540 fixed-U searches at 271 distinct cap points, 270 copies of the extremal point; stored best J flagged as not converged; re-optimized table; correlation and the 1.66× margin withdrawn; `margin4.py` described as one line experiment with its limits. §5.6 and item 3 reworded. README claim 7 updated. Phase 3 line removed from `bad4.py` | Re-ran `moc/reanalyse4.py`: lines (a) and (b) reproduce exactly (invariant depth −0.067442 for all 270 deepened unitaries; corr −0.93 overall, −0.17 raw, −0.06 deepened; stored violation 0.0294–0.0389). Line (c) does not reproduce: see R1 | **Addressed except R1** |
| 2. Zheng's example on the family only up to symmetry (drury-constant.md item 5) | Item 5 now says "up to symmetry", with `moc/zheng_on_family.py`. Item 1 now says O(y) is the half φ ∈ [0, π] | `zheng_on_family.py` re-run: output byte-identical to the committed `runs/zheng_on_family.txt`. My own exact check agrees: O(1/7) ≠ O*, 24 representations D₁P_αVP_βD₂, g(O(1/7)) = −158/2401, minimum over the 1,152 images = −158/2401. Item 1 parametrization: O(y) = cos φ·I + sin φ·S/√3 with y = cos φ and the nonnegative root, so φ ∈ [0, π]; correct | **Addressed** |
| 3. "[Rigorous]" on report §5.2 | Label now "[Rigorous, given the facet list of Q₄]", plus a dependence note citing an exact enumeration in `../star-halfline/q_facets_exact.py`. README claim 6 updated | Re-ran `q_facets_exact.py` from the export: all 1,961,256 orientation determinants of 14-subsets in int64 Bareiss elimination with overflow and divisibility asserts (largest intermediate 64), result "facets of Q: 82; equal to the union of 16 (B = 0) + 18 (D = 0) + 48 (exotic): True", 13 s. I read the script and its two helpers (`../pmsh/star/base.py`, `bareiss` in `../star-degree4/chiro_exact.py`). A nonzero determinant shows the 13 coordinates are an affine chart of the 13-dimensional hull, and the facet test (other vertices weakly on one side, not all on the hyperplane) is the standard one | **Addressed.** With the exact list, §5.2 is rigorous modulo that script |
| 4. Stale metadata | Headers of report.md, theorem-d.md, proposition-e.md updated (claim time, registered agent, review status). Artifact commits pinned in theorem-d.md, proposition-e.md, pmsh.md, drury-constant.md | 78908ca, 8c1b1fb and 304f790 are commits of the repository. theorem-d.md and proposition-e.md are identical at 78908ca and 304f790; pmsh.md and drury-constant.md are identical at 8c1b1fb and 304f790. `git grep` at 3e3226c finds no remaining "not claimed", "not a registered ledger agent", "not independently reviewed" or "the commit that adds this file" in `minors/` | **Addressed** (nit R3) |
| 5. Proposition E can drop Drury 1992 | Gap note in proposition-e.md: the n = 4 input can be Zheng's rational O*, since g ≥ 0 on all 24 permutation matrices and g = −158/2401 at a symmetric image of O* | Exact (`fractions`): g takes values 0 (×13), 1 (×10), 2 (×1) on the permutation matrices, and g(O(1/7)) = −158/2401 with O(1/7) rational orthogonal. Since g is linear in the squared minors, this separates r(O(1/7)) from the convex hull of the permutation vectors, which is all Proposition E needs from Drury 1992 | **Addressed** |

## 2. Remaining findings

**R1. The re-optimized J closest to zero is about −3.6·10⁻⁵, not about −7·10⁻⁵ (needs a correction).**

The revision prints the re-optimized values as −7.7·10⁻⁵, −6.6·10⁻⁵ and −1.0·10⁻⁴ (report §5.5 table, README "Response" row 1). It then summarizes them as "about −7·10⁻⁵" (bottom-line item 3) and "down to about −7·10⁻⁵" (§5.6).

- Re-running `moc/reanalyse4.py` unchanged (seed 0) in my environment gives −7.731·10⁻⁵, **−4.063·10⁻⁵** and −9.988·10⁻⁵. So line (c) of `runs/reanalyse4.txt` is environment-dependent; lines (a) and (b) reproduce exactly.
- A longer re-optimization at the same fixed U has 40 restarts, perturbation scales 10⁻⁶ to 10⁻², adaptive Nelder–Mead and 40,000 evaluations per restart. It converges to the same values for two seeds:

  | record | stored J | revision's value | my longer re-optimization (seeds 1 and 2) | at the optimum |
  | --- | --- | --- | --- | --- |
  | s19, raw | −2.59·10⁻⁴ | −7.7·10⁻⁵ | −7.695·10⁻⁵ | two aᵢ coincide (relative gap ~10⁻¹⁴) |
  | s17, raw | −2.79·10⁻⁴ | −6.6·10⁻⁵ | **−3.648·10⁻⁵** | two aᵢ coincide (relative gap ~10⁻¹⁴) |
  | s11, raw | −2.88·10⁻⁴ | −1.0·10⁻⁴ | −9.973·10⁻⁵ | two bᵢ coincide (relative gap ~10⁻¹⁵) |

- My own first review had the same problem: it gave −7.0·10⁻⁵ to −7.3·10⁻⁵ for the s17 record. That was also not converged, and I correct it here.
- The qualitative statement still holds. All values are negative, and the optima sit at a repeated eigenvalue. There the spectrum {α, α, β, γ} has three distinct values with a simple one, so Theorem D already gives J ≤ 0.

*Smallest correction.*

- Replace "about −7·10⁻⁵" (bottom-line item 3, §5.6) by "about −4·10⁻⁵".
- Replace −6.6·10⁻⁵ in the §5.5 table and in the README response row by "−3.6·10⁻⁵ (longer re-optimization by the re-reviewer; −4.1·10⁻⁵ to −6.6·10⁻⁵ from `reanalyse4.py` depending on environment)", or by "≈ −4·10⁻⁵".
- Add once in §5.5 that the re-optimized values are optimizer outputs, not certified maxima. The checkpoint text says so, but the report does not.

**R2. Stale file-table line (nit).** The `moc/` table near the end of report.md still describes `bad4.py`, `analyse4.py`, `margin4*.py` as "270k targeted search, depth correlation, safety margin". Suggested: "270k fixed-U search; line experiment (§5.5)", and add a row for `reanalyse4.py` and `zheng_on_family.py`.

**R3. "Unchanged on main at 304f790 apart from this header" (nit).** In theorem-d.md and proposition-e.md the file at 304f790 is byte-identical to 78908ca. The header was changed later, in 3e3226c. Suggested: "78908ca (first pass); unchanged through 304f790; header updated in 3e3226c".

**R4. drury-constant.md item 1 can be exact (suggestion).** The note "the other half gives the same minimum (numerically, on a 2001-point grid by the reviewer)" can cite an exact one-line argument instead. S is skew, so the other half is {O(y)ᵀ} = {−O(−y)}. Every term of g is even in U, so g(O(y)ᵀ) = g(O(−y)). I confirmed the two closed forms symbolically: 2y² ∓ (4√3/9)·y·(1−y²)^{3/2}.

## 3. Reproduction

From the repository root at the reviewed commit (an export with `git archive` works):

```sh
python3 -I research/marcus-de-oliveira/reviews/check_minors_ishaan.py      # first-review checks; all exact parts unchanged
cd research/marcus-de-oliveira/minors && shasum -a 256 -c SHA256SUMS       # 152 files OK
python3 moc/zheng_on_family.py        # byte-identical to moc/runs/zheng_on_family.txt
python3 moc/reanalyse4.py             # (a), (b) identical; (c) environment-dependent (R1); 9 s
cd ../star-halfline && python3 q_facets_exact.py   # 82 facets, exact; 13 s, about 4.7 GB peak memory
```

The longer re-optimization behind R1 is a scratch script of about 30 lines. It is not committed. It loads the three best records from `moc/runs/bad4_s1*.json` and evaluates J with `Problem.evaluate` from `moc/bsearch.py`, exactly as `reanalyse4.py` does. It then runs 40 restarts of `scipy.optimize.minimize(method='Nelder-Mead', options={'adaptive': True, 'maxiter': 40000, 'maxfev': 40000, 'xatol': 1e-13, 'fatol': 1e-16})` from the best point so far, perturbed by `10**U(-6,-2)` times a standard normal vector, for seeds 1 and 2. Run time: 200 s.

Environment: Python 3.14.7, numpy 2.5.5, scipy 1.18.1, sympy 1.14.0.

## 4. Limits of this re-review

- I checked only what changed between 304f790 and 3e3226c in `minors/`, plus the facet script cited by §5.2. The scope and limits of the first review (its section 4) still apply to everything else.
- I did not check the ledger registration time of `ganesh-claude-minors` stated in the new headers.
- The re-optimized J values, mine included, are optimizer outputs. They are not certified maxima.
- The `margin4.py` line experiment was not re-run.
- **Scope of a pass after R1:** report §5.5–5.6 and bottom-line item 3 would then state only what the archived runs support. The mathematical content keeps the scope given in the first review: Theorem D, Proposition E (now with an exact n = 4 input) and pmsh.md H1–H5 as stated, and report §5.2 given the exact facet list. PMSH(n) for n ≥ 4 is unproved. The global optimality of Drury's constant is numerical. Nothing here bears on MOC itself.
