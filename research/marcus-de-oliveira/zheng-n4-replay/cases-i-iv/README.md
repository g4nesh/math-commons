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

Cases in progress: T₃–T₇, T₁₀, T₁₁, T₁₅, T₁₈–T₂₀.

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
