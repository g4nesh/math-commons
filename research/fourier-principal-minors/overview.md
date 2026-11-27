# Fourier principal minors · Start here

**Current objective.** Determine whether every principal minor of the complex Fourier matrix of order 154 is nonzero, independently establish what the literature already covers, then audit the prime-multiplier family of orders 14p. The kickoff contains contributor evidence awaiting independent review.

For a primitive complex Nth root of unity ζ, the entry in row i and column j of the Fourier matrix is ζ^(ij), with indices 0 through N−1. A principal minor uses the *same* set of indices for rows and columns. The larger conjecture concerns all square-free N. [Conjecture and small cases](https://arxiv.org/abs/2409.09793v4)

| When | Task | Start with this page |
|---|---|---|
| Now, parallel | `fourier-154-replay` | [Independently reconstruct the order-154 proof](order-154/replay.md) |
| Now, parallel | `fourier-154-prior-art` | [Check priority and earlier theorems](order-154/prior-art.md) |
| After both reports receive independent review | `fourier-154-synthesis` | [Resolve discrepancies and record the order-154 decision](order-154/synthesis.md) |
| After reviewed synthesis clears the handoff | `fourier-14p-audit` | [Audit the 14p family](family-14p/README.md) |

Run `node bin/commons.mjs context --objective fourier-principal-minors` from the repository root for ownership and checkpoints. Pick one available task and read its linked page. Cold reviewers start from the replay brief and record their own results before consulting comparison artifacts. [GitHub tracking issue](https://github.com/g4nesh/math-commons/issues/130) · [Full project index](README.md) · [Independent review rules](../../docs/verification.md#independent-review-procedure)
