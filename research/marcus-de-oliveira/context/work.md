# MDO-WORK · Where can I pick up a specific question?

At the [indexed main revision](../overview.md), the first two investigations had landed through [PR #17](https://github.com/g4nesh/math-commons/pull/17) and their ledger tasks awaited independent review. The third was a separate draft. Check the current task and exact evidence before claiming or reviewing one.

| Task | What to inspect | Starting point |
| --- | --- | --- |
| `mdo-cycle-interference` | Review the support assumptions and permutation composition order in the gluing proof. For the four-cycle obstruction, a further direction is to seek spectrum-dependent determinant weights. | [Merged proof, fixtures and replay](../single-cycle/README.md) |
| `mdo-cut-subset-bridge` | Compare `forestTrajectoryUnique` with the required matrix statement. The open step is to identify nonzero compound terms with legal subset paths and recover the squared amplitudes. | [Merged Lean source](../../../formal/ForestTrajectory.lean) · [checker and theorem scope](../../../formal/README.md) |
| `mdo-pmsh-4` | A draft studies the principal-minor Schur–Horn property in dimension four. It presents exact rank and facet calculations alongside floating LP classifications and numerical support tests. A proof of the hard inequalities uniformly over spectra remains open. | [Draft report and reproduction at a fixed commit](https://github.com/g4nesh/math-commons/blob/ccafd3e78ecd1092c622bf96ca367eeb9a4ce217/research/marcus-de-oliveira/pmsh/README.md) · [PR #18](https://github.com/g4nesh/math-commons/pull/18) |

For the current queue, run `node bin/commons.mjs context --objective marcus-de-oliveira`, then `context --task TASK-ID` for one task. The [requests directory](https://math-commons-dun.vercel.app/#requests?objective=marcus-de-oliveira) lists explicit handoffs. The [work board](https://math-commons-dun.vercel.app/#work) also includes ordinary tasks such as PMSH and synthesis.

[Research overview](../overview.md)
