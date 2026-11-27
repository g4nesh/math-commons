# Independent replay of the contact cases of group V, second half (Zheng's proof for n = 4)

This page tracks eleven cases. The verdict table lists completed replays and links to their reports and scripts as they become available.

- **Task:** `mdo-zheng-cases-v-b`, requested by the audit `mdo-zheng-n4-audit` (owner `IshaanAyaan`). Replay by `ganesh-claude-minors` (Claude Code, model `claude-fable-5-2`) and its workers, GitHub owner `g4nesh`. November 26, 2026.
- **Procedure.** As in [the first half of group V](../cases-v-a/README.md), each worker reads the proof documents in the author’s public supplement, replays the identities and certificates with its own exact code, and writes its verdict before opening the audit’s first-pass report. Scripts read the supplement as data from a local clone named by `MOC4_SUPP`. The supplement’s programs remain outside the replay.
- **Group V.** The real part of a point of the matching cone satisfies C = 0 and a Hessian equation. Written in the generators of the contact cone, these exclude, by polynomial sign certificates, the region where every universal inequality on SO(4) fails.
- **Status.** Reports are added as workers finish. A verdict here is one independent reading by a second owner, not a formal verification. The coordinator re-runs every replay script before committing it.
- **Evidence levels and verdict words** are those of [`../../zheng-n4-audit/METHOD.md`](../../zheng-n4-audit/METHOD.md).

## Verdicts

| Case | Verdict | Report | Replay | Notes |
| --- | --- | --- | --- | --- |
| T₇₂ = {1,3,4,7,8,12} | confirmed modulo citations | [T72.md](T72.md) | [`replay_T72.py`](replay_T72.py): 92 exact and 8 numeric checks | The document works on the image {0,1,2,10,14,17} of the set under a group element; the transport is verified, including that it preserves the matching cone and transforms the Leibniz sums correctly. Nine generators (the cone is not simplicial); both signs of the non-minor coefficient are treated. The degree-9 polynomial (1,615 terms) is computed here and the radial identity holds for every value of that coefficient; all 17 tables (496,798 coefficients, all positive) are rebuilt with the document's term counts. Cites the SOC, old-H and chain inequalities and nothing else. All 19 finite assertions of the supplement's portable check `checks/fresh_T72_audit.py` are re-established with own code. |
| T₈₈ = {1,5,6,10,13,14} | confirmed modulo citations | [T88.md](T88.md) | [`replay_T88.py`](replay_T88.py): 79 exact and 6 numeric checks | A negative non-minor coefficient is reduced to a positive one by the inversion σ ↦ σ⁻¹ (checked exactly). The degree-9 polynomial (2,285 terms) is computed here; all 10 tables (185,823 coefficients, all positive) are rebuilt with the document's term counts. Cites only the SOC inequality, at three pairs. All 17 finite assertions of the portable check `checks/fresh_T88_audit.py` are re-established. The audit's stratum test had found only points where the non-minor coefficient is about 0; 30 points where it is not are found here, and they satisfy the radial equation and the conclusion. Shares the helper [`replay_T72_common.py`](replay_T72_common.py). |

Cases in progress: T₃₈, T₃₉, T₄₀, T₄₂, T₄₃, T₄₄, T₄₆ (one proof document each), and T₄₇ and T₆₅, whose proofs are assemblies of 11 and 22 documents. Those two are replayed in parts (`T47_A`–`T47_C`, `T65_A`–`T65_D`): every document belongs to exactly one part, each part states the exact region and conclusion of every note it replays, and the coordinator checks that the parts fit together.

## Reproduction

The planned replay environment uses Python 3.9 or later with sympy 1.14 and python-flint 0.6 for exact work, and numpy and scipy for falsification attempts. Set `MOC4_SUPP` to a local clone of the supplement at commit `53ec16567d2307266ce52a8e9678410b2457cc27`. Each published replay will include its invocation and expected result.
