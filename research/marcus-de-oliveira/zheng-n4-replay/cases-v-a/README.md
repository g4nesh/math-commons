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
| T₂₇ = {2,3,4,7,8,12} | confirmed modulo citations | [T27.md](T27.md) | [`replay_T27.py`](replay_T27.py): 54 exact, 6 numeric, 5 data comparisons | A direct proof; the two boundary documents of the supplement are not needed by the main one. The factorization of the 14×14 Hessian determinant is verified with the document's degree-6 polynomial (401 terms). Ends with the SOC and four-coefficient inequalities (cited; orbit members and coefficient conditions exact). |
| T₂₈ = {2,3,5,6,9,12} | confirmed modulo citations | [T28.md](T28.md) | [`replay_T28.py`](replay_T28.py): 84 exact, 5 numeric, 5 data comparisons | The degree-9 polynomial (2,124 terms), all 16 constant-term tables, all 24 Bernstein tables (up to 100,982 terms) and 302 Sturm groups are recomputed; every term count matches the document. It cites Lemma 1 of the document for the paper's T₂₅ (coherent, parallel and double-SOC bounds); that lemma is tested numerically only in this report, and is verified exactly in [T25.md](T25.md). No point of the matching cone was found in the hard chamber, so the endpoint inequalities were tested on the coefficient region. |
| T₂₅ = {3,4,7,8,10,12} | confirmed modulo citations | [T25.md](T25.md) | [`replay_T25.py`](replay_T25.py): 71 exact and 9 numeric checks | Uses the old-H inequality on SO(4) in one region; an exact rank check shows that the document's own Gram representation cannot supply it there. The document's Lemma 1 (coherent, parallel and double-SOC bounds) is verified exactly: its sum-of-squares identity is checked against Leibniz sums rebuilt from the paper, with the 16 cubic features taken from the supplement's data, since no PDF contains them. The degree-9 polynomial (2,095 terms) is computed here from the 14×14 Hessian determinant. |
| T₂₆ = {2,3,5,7,8,12} | confirmed modulo citations | [T26.md](T26.md) | [`replay_T26.py`](replay_T26.py): 96 exact and 12 numeric checks | Old-H in one region, as for T₂₅. All tables recomputed (the largest has 1,066,354 terms); every count stated in the document matches. The data file that the document says records eight binomial squares is not in the public supplement; the squares were found here. **T₂₆ does not need T₂₅**: its bounds are verified directly on its own face. |
| T₃₅ = {2,3,4,6,9,12} | confirmed modulo citations | [T35.md](T35.md) | [`replay_T35.py`](replay_T35.py): 82 exact and 11 numeric checks | Self-contained: nothing is delegated to Appendix B, Theorem C.1 or another case. 19 tables regenerated, the largest with 289,848 terms. Cites old-H, chain and four-coefficient inequalities; orbit members exact. |
| T₃₆ = {2,3,4,6,11,12} | confirmed modulo citations | [T36.md](T36.md) | [`replay_T36.py`](replay_T36.py): 92 exact and 14 numeric checks | Self-contained. 10 tables regenerated. Cites four-coefficient, SOC, old-H and chain inequalities. One inequality written as strict in the document holds only non-strictly on a boundary; the conclusion is unaffected. |
| T₃₇ = {0,5,6,7,9,12} | confirmed modulo citations | [T37.md](T37.md) | [`replay_T37.py`](replay_T37.py): 115 exact and 19 numeric checks | Self-contained. 85 tables regenerated. Cites SOC, four-coefficient, old-H and chain inequalities and Lemma 1 of the T₂₅ document (verified exactly in [T25.md](T25.md)). A statement written for O(4) holds on SO(4) only. |

Cases in progress: T₂₉, T₃₁.

## Reproduction

Python 3.9 or later with sympy 1.14 and python-flint 0.6 (exact replays), numpy and scipy (falsification attempts only). From this directory:

```sh
export MOC4_SUPP=/path/to/clone/of/the/supplement   # commit 53ec16567d2307266ce52a8e9678410b2457cc27
python3 replay_T34.py     # prints its checks and a final PASS or FAIL line (about 3 minutes)
```
