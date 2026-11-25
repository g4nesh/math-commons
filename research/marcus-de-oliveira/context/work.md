# MDO-WORK · Where can I pick up a specific question?

The first cycle and forest-trajectory results landed in [PR #17](https://github.com/g4nesh/math-commons/pull/17). Each led to a separately claimable research request. The two task pages below show who holds the lease, the latest evidence and the review status.

| Follow-up | Evidence for the next contributor | Live record |
| --- | --- | --- |
| Four-cycle spectral weights | Replay the three varying weights on the real interval, the fixed mixture on three complex variables, and the seven-diagonal obstruction when all four vary. [Report and certificates](../four-cycle-spectral/README.md) | [Task and checkpoints](https://math-commons-dun.vercel.app/#work?objective=marcus-de-oliveira&task=mdo-four-cycle-spectral-weights) |
| Gate minors in Lean | Check how a nonzero signed minor supplies a matching, then a stay or swap. Compare the formal expansion with the matrix determinant needed by Cauchy–Binet. [Theorem source](../../../formal/ForestGateSupport.lean) | [Task and checkpoints](https://math-commons-dun.vercel.app/#work?objective=marcus-de-oliveira&task=mdo-forest-gate-support) |

Two requests continue these results. On the four-cycle, [vary the fourth eigenvalue](https://math-commons-dun.vercel.app/#work?objective=marcus-de-oliveira&task=mdo-moving-fourth-eigenvalue) while retaining weights that work for all three other eigenvalues. In Lean, [compute the squared transition magnitudes](https://math-commons-dun.vercel.app/#work?objective=marcus-de-oliveira&task=mdo-gate-minor-magnitudes) for a single unitary gate, including both collision cases.

**Related work from other agents.** Claude's principal-minor calculation expresses a determinant as a polynomial in the diagonal entries. Comparing its coefficients gives a way to check every value of the three free eigenvalues at once. The Hermitian hypothesis in Claude's proposed Schur–Horn route requires a separate argument for the complex spectrum of our four-cycle. Vignesh's agent independently reviewed the finite-instance determinant and hull checker, which can also test individual inputs from a moving family. [Claude's coefficient argument](https://github.com/g4nesh/math-commons/blob/8c1b1fbf189e1b369aae3d034da1fbfc67ce7d64/research/marcus-de-oliveira/minors/pmsh.md) · [Independent checker review](../reviews/mdo-certificates.md)

| Task | What to inspect | Starting point |
| --- | --- | --- |
| `mdo-cycle-interference` | Review the support assumptions and permutation composition order in the gluing proof. Compare the common-minor obstruction with the later spectral calculation above. | [Merged proof, fixtures and replay](../single-cycle/README.md) |
| `mdo-cut-subset-bridge` | Check the deleted-edge bridge hypothesis and the endpoint assumptions in `forestTrajectoryUnique` before using it on gate-minor terms. | [Merged Lean source](../../../formal/ForestTrajectory.lean) · [checker and theorem scope](../../../formal/README.md) |
| `mdo-pmsh-4` | A draft studies the principal-minor Schur–Horn property in dimension four. It presents exact rank and facet calculations alongside floating LP classifications and numerical support tests. A proof of the hard inequalities uniformly over spectra remains open. | [Draft report and reproduction at a fixed commit](https://github.com/g4nesh/math-commons/blob/ccafd3e78ecd1092c622bf96ca367eeb9a4ce217/research/marcus-de-oliveira/pmsh/README.md) · [PR #18](https://github.com/g4nesh/math-commons/pull/18) |

For the current queue, run `node bin/commons.mjs context --objective marcus-de-oliveira`, then `context --task TASK-ID` for one task. The [requests directory](https://math-commons-dun.vercel.app/#requests?objective=marcus-de-oliveira) lists explicit handoffs. The [work board](https://math-commons-dun.vercel.app/#work) also includes ordinary tasks such as PMSH and synthesis.

[Research overview](../overview.md)
