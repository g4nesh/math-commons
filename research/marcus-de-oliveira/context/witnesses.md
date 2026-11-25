# MDO-TEST · What has been tested, and how is a witness certified?

| Evidence | Result and limits | Full evidence |
| --- | --- | --- |
| Numerical screen, independently reviewed | 200 floating-point samples in dimensions 2–6 produced no gap above the heuristic allowance. The reviewer constructed 200 **nearby exactly unitary inputs** and certified each inside. Exact membership applies to those reconstructed inputs. | [Screen and replay](../search/README.md) · [review and reconstruction](../reviews/mdo-search.md) |
| Exact stress screen, author-side differential replay | All 924 distinct ordered pairs in dimensions 4–7 are inside, including 24 dimension-seven holdouts. Distinctness is by exact matrix entries. Search coverage is finite. | [Report, separate holdout and replay](../rational-stress/README.md) · [distinctness audit](../rational-stress/distinct-inputs.json) |
| Exact witness verifier, independently reviewed | For an accepted input, it checks unitarity and normality, then returns explicit convex weights or a strict rational separator against every permutation point. It handles points, segments and polygons. | [Schema, limits and CLI](../certificates/README.md) · [reviewed checkpoint and tests](../reviews/mdo-certificates.md) |

**If you have a candidate**, preserve its exact spectra and unitary construction. Rounding approximate entries to fractions can destroy unitarity. The general file verifier accepts dimensions 1–6. The dimension-seven holdout uses its own harness. Reproduce a strict outside gap independently before claiming a counterexample.

[Research overview](../overview.md)
