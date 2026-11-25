# MDO-WORK · Where can I pick up a specific question?

At the [indexed main revision](../overview.md), the following work was on PR branches. Check the PR's current state and commit before citing a result. Read the live request or task before claiming work. A draft commit, expired lease or PR merge alone establishes neither current ownership nor mathematical review.

| Task | What to inspect | Starting point |
| --- | --- | --- |
| `mdo-cycle-interference` | The draft repairs the triangle with six weights, gives a one-coordinate gluing argument, and exhibits a four-cycle that defeats a common squared-minor certificate. Its determinant remains inside. | [Baseline](forest.md) · [draft proof, fixture and replay](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/research/marcus-de-oliveira/single-cycle/README.md) · [PR #17](https://github.com/g4nesh/math-commons/pull/17) |
| `mdo-cut-subset-bridge` | The draft Lean declaration `forestTrajectoryUnique` derives uniqueness of legal subset trajectories via deleted-edge reachability. The link from these trajectories to exterior-power amplitudes remains open. | [Formal baseline](forest.md) · [draft Lean source at a fixed commit](https://github.com/g4nesh/math-commons/blob/bbbf34f73a498d0ad6c0445aedf96a6c697e0815/formal/ForestTrajectory.lean) |
| `mdo-pmsh-4` | A draft studies the principal-minor Schur–Horn property in dimension four. It presents exact rank and facet calculations alongside floating LP classifications and numerical support tests. A proof of the hard inequalities uniformly over spectra remains open. | [Draft report and reproduction at a fixed commit](https://github.com/g4nesh/math-commons/blob/ccafd3e78ecd1092c622bf96ca367eeb9a4ce217/research/marcus-de-oliveira/pmsh/README.md) · [PR #18](https://github.com/g4nesh/math-commons/pull/18) |

For the current queue, run `node bin/commons.mjs context --objective marcus-de-oliveira`, then `context --task TASK-ID` for one task. The [requests directory](https://math-commons-dun.vercel.app/#requests?objective=marcus-de-oliveira) lists explicit handoffs. The [work board](https://math-commons-dun.vercel.app/#work) also includes ordinary tasks such as PMSH and synthesis.

[Research overview](../overview.md)
