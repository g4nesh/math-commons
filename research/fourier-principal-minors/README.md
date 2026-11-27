# Square-free Fourier principal minors

**Project brief.** For a square-free positive integer N, define P(N) to mean that for every S ⊆ {0,…,N−1}, the determinant of the submatrix (ζ_N^(ij)) indexed by i,j ∈ S is nonzero. Fix a primitive Nth root in ℂ. The empty determinant is one. Changing the primitive root is a Galois conjugation, and the customary nonzero Fourier normalization preserves vanishing. The first milestone is P(154). The second is an audit of P(14p) for prime p with 14p square-free.

New workers should begin with the [short overview](overview.md) and the live `fourier-principal-minors` ledger. The order-154 replay and literature assignments are deliberately parallel so that a numerical outcome cannot steer the priority search.

| Evidence or decision | Location |
|---|---|
| Cold reconstruction brief | [Order-154 independent replay](order-154/replay.md) |
| Literature scope, dates and entailment | [Prior-art assignment](order-154/prior-art.md) |
| Owner-supplied proof sketch and reproducible baseline, for comparison after a cold run | [Order-154 comparison packet](order-154/seed.md) |
| Gate for the next stage, including known results and disagreements | [Order-154 synthesis](order-154/synthesis.md) |
| Preliminary family argument and exact code | [14p audit packet](family-14p/README.md) |

The seeded computations and proof sketches were supplied by the kickoff contributor. Report review and priority separately. Attribute an earlier result when it covers a case, even when the project's reconstruction succeeds independently. An independent review in this workspace requires a GitHub owner distinct from every contributor to the submitted checkpoint. Follow the [review procedure](../../docs/verification.md#independent-review-procedure) when recording a verdict.

The earlier Marcus–de Oliveira Fourier work is available for historical method comparison. It studies distributions of squared minors and a different conjecture. [Historical Fourier workspace](../marcus-de-oliveira/fourier-composite-laws/README.md)
