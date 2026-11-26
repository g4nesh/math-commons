# Independent replay of the contact cases of group V, first half (Zheng's proof for n = 4)

- **Task:** `mdo-zheng-cases-v-a`, requested by the audit `mdo-zheng-n4-audit` (owner `IshaanAyaan`). Replay by `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`) and its workers, GitHub owner `g4nesh`. November 26, 2026.
- **What this is.** The same procedure as [`../cases-i-iv/`](../cases-i-iv/README.md): for each case one worker read the proof document in the author's public supplement and replayed every identity and certificate with its own exact code, without opening the audit's first-pass report on the same document until its verdict was written. The supplement's programs were never executed. Nothing from the supplement is copied here; scripts read it from a local clone named by the environment variable `MOC4_SUPP`.
- **Group V.** The real part of a point of the matching cone satisfies C = 0 and a Hessian equation. Written in the generators of the contact cone, these exclude, by polynomial sign certificates, the region where every universal inequality on SO(4) fails.
- **Status.** Reports are added as workers finish. A verdict here is one independent reading by a second owner, not a formal verification. The coordinator re-ran every replay script before committing it.
- **Evidence levels and verdict words** are those of [`../../zheng-n4-audit/METHOD.md`](../../zheng-n4-audit/METHOD.md).

## Verdicts

| Case | Verdict | Report | Replay | Notes |
| --- | --- | --- | --- | --- |
| T₃₄ = {2,3,5,9,10,12} | confirmed modulo citations | [T34.md](T34.md) | [`replay_T34.py`](replay_T34.py): 90 exact and 11 numeric checks | The paper's representative of the group. The document does not state or prove the two "endpoint inequalities" behind its condition (H); they are the old-H inequality on SO(4), which is inequality (★) of [`../../star-orthogonal/README.md`](../../star-orthogonal/README.md) (Theorem O) and is also proved in the audit's notes. The chain inequality is verified exactly from the supplement's certificate data; no PDF contains it. The degree-9 polynomial of the exclusion is rebuilt here from the 14×14 Hessian determinant by a full symbolic expansion (1,381 terms), and the twelve Bernstein numerators have positive integer coefficients. |

Cases in progress: T₂₅–T₂₉, T₃₁, T₃₅–T₃₇.

## Reproduction

Python 3.9 or later with sympy 1.14 and python-flint 0.6 (exact replays), numpy and scipy (falsification attempts only). From this directory:

```sh
export MOC4_SUPP=/path/to/clone/of/the/supplement   # commit 53ec16567d2307266ce52a8e9678410b2457cc27
python3 replay_T34.py     # prints its checks and a final PASS or FAIL line (about 3 minutes)
```
