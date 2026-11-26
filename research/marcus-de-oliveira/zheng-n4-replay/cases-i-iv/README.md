# Independent replay of the contact cases of groups I–IV (Zheng's proof for n = 4)

- **Task:** `mdo-zheng-cases-i-iv`, requested by the audit `mdo-zheng-n4-audit` (owner `IshaanAyaan`). Replay by `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`) and its workers, GitHub owner `g4nesh`. November 26, 2026.
- **What this is.** For each case, one worker read the proof document in the author's public supplement and replayed every identity and certificate with its own exact code. It did not open the audit's first-pass report on the same document until its verdict was written. The supplement's programs were never executed. Nothing from the supplement is copied here; scripts read it from a local clone named by the environment variable `MOC4_SUPP` when they need its data.
- **Status.** Reports are added as workers finish. A verdict here is one independent reading by a second owner. It is not a formal verification, and the cited inputs named in each report are not replayed here unless the report says so.
- **Evidence levels and verdict words** are those of [`../../zheng-n4-audit/METHOD.md`](../../zheng-n4-audit/METHOD.md).

## Verdicts

| Case | Group | Verdict | Report | Replay | Notes |
| --- | --- | --- | --- | --- | --- |
| T₁ = {0,1,2,3,6,7} | I | confirmed modulo citations | [T1.md](T1.md) | [`replay_T1.py`](replay_T1.py): 73 exact checks | The document's only case citation (the one-circle theorem, for a repeated eigenvalue) can be replaced by an elementary argument given in the report, with exact checks. |
| T₂ = {0,1,2,3,6,8} | II | confirmed modulo citations | [T2.md](T2.md) | [`replay_T2.py`](replay_T2.py): 57 exact checks | As for T₁. |
| T₂₁ = {0,3,5,9,11,12} | IV | confirmed modulo citations | [T21.md](T21.md) | [`replay_T21.py`](replay_T21.py): 48 exact checks | Coincident eigenvalues are delegated by the document to the paper's Appendix B; the report narrows the delegated set to eight named coincidences. No stratum point with eight distinct labels was found numerically for this case. |
| T₂₂ = {2,3,4,8,10,12} | IV | confirmed modulo citations | [T22.md](T22.md) | [`replay_T22.py`](replay_T22.py): 64 exact checks | As T₂₁ for the delegated coincidences. With distinct labels a seventh contact (index 13) is forced (exact). |
| T₂₃ = {2,3,5,9,11,12} | IV | confirmed modulo citations | [T23.md](T23.md) | [`replay_T23.py`](replay_T23.py): 52 exact checks | As T₂₂. |
| T₃ = {0,1,2,3,7,8} | III | confirmed modulo citations | [T3.md](T3.md) | [`replay_T3.py`](replay_T3.py): 70 exact checks | Ends with the SOC inequality, member (f1, f11; f17); orbit membership and the coefficient condition checked exactly. |
| T₄ = {0,1,4,5,7,8} | III | confirmed modulo citations | [T4.md](T4.md) | [`replay_T4.py`](replay_T4.py): 84 exact checks | Ends with the four-coefficient inequality. The document cites the three-plus-three circle theorem; the report shows by two exact identities that this citation is not needed. Three omissions in the document are filled by exact identities. |
| T₅ = {1,2,3,4,6,8} | III | confirmed modulo citations | [T5.md](T5.md) | [`replay_T5.py`](replay_T5.py): 110 exact checks | SOC on two branches and the four-coefficient inequality on the third; the four-budget note's certificate is replayed exactly. |
| T₆ = {0,1,2,4,6,8} | III | confirmed modulo citations | [T6.md](T6.md) | [`replay_T6.py`](replay_T6.py): 98 exact checks | Ends with the SOC representative itself. The stratum aᵢ + b_j = 0 is proved in the document and replayed exactly. |
| T₁₀ = {0,1,3,5,7,8} | III | confirmed modulo citations | [T10.md](T10.md) | [`replay_T10.py`](replay_T10.py): 96 exact checks | SOC and four-coefficient inequalities, orbit members and coefficient conditions checked exactly. The cited three-plus-three circle theorem is not needed: for real v the whole gap vector is zero (exact). |
| T₁₁ = {0,1,3,5,6,8} | III | confirmed modulo citations | [T11.md](T11.md) | [`replay_T11.py`](replay_T11.py): 92 exact checks | The documents use {0,1,2,4,7,9}; exactly one group element maps it to T₁₁, and the transfer is checked for all 82 vectors. Relies on both branch documents and on case T₇ at one boundary. The audit's first pass recorded a presentation gap in the B-branch (reality of four chart parameters); the replay closes it exactly. Theorem C.1 is cited for real v. |
| T₁₅ = {1,2,3,4,6,7} | III | confirmed modulo citations | [T15.md](T15.md) | [`replay_T15.py`](replay_T15.py): 52 exact checks | SOC inequality. Here the three-plus-three circle theorem is genuinely needed for real v. |
| T₇ = {0,1,2,4,7,8} | II | confirmed modulo citations | [T7.md](T7.md) | [`replay_T7.py`](replay_T7.py): 67 exact checks | The document cites Theorem C.1 for real v. The replay proves that the citation is not needed: with eight distinct labels and real v the gap vector is zero. |
| T₁₈ = {0,1,3,6,7,8} | II | confirmed modulo citations | [T18.md](T18.md) | [`replay_T18.py`](replay_T18.py): 59 exact checks | The cited three-plus-three circle theorem is not needed: for real v all labels are real and every gap vanishes. |
| T₁₉ = {0,1,4,6,7,8} | II | confirmed modulo citations | [T19.md](T19.md) | [`replay_T19.py`](replay_T19.py): 69 exact checks | As T₁₈. In addition the main layer is empty: with eight distinct labels and nonreal v a gap is negative in each of the document's two regions, so the document's positivity identities are correct but apply to no point. |
| T₂₀ = {0,1,5,6,7,8} | II | confirmed modulo citations | [T20.md](T20.md) | [`replay_T20.py`](replay_T20.py): 77 exact checks | As T₁₈; the interlacing step is confirmed by a root-free exact certificate. |

All sixteen cases of the request are replayed. **No gap was found in any of them.** Every verdict is "confirmed modulo citations".

## Citations that the cases share

Every report lists its citations. Across the cases replayed so far they are:
- **Coincident eigenvalues** (a repeated eigenvalue in one spectrum, or aᵢ + b_j = 0): the paper's Appendix B (Theorem B.1, Lemma B.2) and Theorem C.1. Not replayed here. For T₁ and T₂ the report gives an elementary replacement. Corollaries P4′ and P4″ of [`../../unitary-factorization/README.md`](../../unitary-factorization/README.md) prove the repeated-eigenvalue case by another route.
- **Theorem C.1 beyond coincidences.** Only T₁₁ uses it (for real v, on a set of real codimension two). T₇ cites it but does not need it. This differs from the audit's dependency summary, which lists Theorem C.1 as indispensable for T₇.
- **The three-plus-three circle theorem.** Cited by T₄, T₁₀, T₁₅, T₁₈, T₁₉ and T₂₀. Among these only T₁₅ needs it; for the other five the replays give exact arguments that avoid it.
- **The universal inequalities on SO(4)** (SOC and four-coefficient so far): the workers re-verified the identities behind them exactly; the proofs are those of [`../../zheng-n4-audit/foundations/`](../../zheng-n4-audit/foundations/).
- **Common inputs of the general argument** (finite spectra for points of the matching cone; the parametrization of SO(4) by pairs of unit quaternions): cited.

## A cross-check between two lines of work

[`atoms_frame.py`](atoms_frame.py) (exact; output in [`logs/atoms_frame.txt`](logs/atoms_frame.txt)) rewrites the universal inequalities on SO(4) of [`../../zheng-n4-audit/PROOF.md`](../../zheng-n4-audit/PROOF.md), Section 7.6, in the coordinates of [`../../star-orthogonal/README.md`](../../star-orthogonal/README.md).

- All 48 non-minor functionals are relabellings of g0 = D[13,13] + D[14,14] + D[13,14] + D[13,23] − D[12,34], with D[S,T] = det O[S,T]².
- The audit's **old-H** inequality becomes g0 + D[12,12] + 2·√(D[12,12]·D[13,14]) ≥ 0. This is inequality (★) of [`../../pmsh/equivalence.md`](../../pmsh/equivalence.md). Theorem O of `star-orthogonal` proves it by a different method (pairs of rotations and a unit-cube lemma that is checked in Lean). So this inequality has two proofs by different owners.
- In the same frame: **SOC** is g0 + 2·|det O[12,12]·det O[13,13]| ≥ 0; **chain** is g0 + 2·|det O[12,12]|·(|O₁₁| + |det O[13,14]|) ≥ 0; **rank-two** is g0 + D[12,12]/4 + |det O[12,12]|·√(D[13,14] + D[13,23]) ≥ 0.

## Reproduction

Python 3.9 or later with sympy 1.14 (exact replays) and numpy, scipy (falsification attempts only). From this directory:

```sh
export MOC4_SUPP=/path/to/clone/of/the/supplement   # commit 53ec16567d2307266ce52a8e9678410b2457cc27
python3 replay_T1.py      # each replay prints its checks and a final PASS or FAIL line
python3 atoms_frame.py    # needs only ../../zheng-n4-audit/lib
```
