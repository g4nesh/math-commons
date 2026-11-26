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

Cases in progress: T₃₈, T₃₉, T₄₀, T₄₂, T₄₃, T₄₄, T₄₆, T₇₂, T₈₈ (one proof document each), and T₄₇ and T₆₅, whose proofs are assemblies of 11 and 22 documents. Those two are replayed in parts (`T47_A`–`T47_C`, `T65_A`–`T65_D`): every document belongs to exactly one part, each part states the exact region and conclusion of every note it replays, and the coordinator checks that the parts fit together.

## Reproduction

The planned replay environment uses Python 3.9 or later with sympy 1.14 and python-flint 0.6 for exact work, and numpy and scipy for falsification attempts. Set `MOC4_SUPP` to a local clone of the supplement at commit `53ec16567d2307266ce52a8e9678410b2457cc27`. Each published replay will include its invocation and expected result.
