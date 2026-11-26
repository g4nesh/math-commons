# Primary-source check of the Hermitian lifting step

**Conclusion:** the substitution $m=k=4$ in Shitov's published Theorem 19 is valid. Consequently Corollary H6 / PR89's Corollary L4 follows **given the four-dimensional Hermitian base case**. No additional common-weight hypothesis is missing. This is a same-owner source comparison, not an independent-owner verdict on the repository's base-case chain or PR89.

Sources inspected on 26 November 2026:

* Y. N. Shitov, *A further sufficient condition for the determinantal conjecture*, Izvestiya: Mathematics 89:4 (2025), 862–869, [primary English PDF](https://www.mathnet.ru/links/1aa847af48a0bf43d77e17d5e515a3f4/im9292_eng.pdf), [DOI](https://doi.org/10.4213/im9292e). The PDF text was accessible after the HTML request timed out. No paper copy is redistributed here.
* [PR89 at immutable f28e48f](https://github.com/g4nesh/math-commons/blob/f28e48fb09aebb0e10f93599914c4dad9aea12d1/research/marcus-de-oliveira/hermitian-lifting/README.md), including its exact statements, proof, numerical code and logs.
* [Four-dimensional base-case source](https://github.com/g4nesh/math-commons/blob/96ef952480b9a5d67b7718be5f2dfc6f6eb85b77/research/marcus-de-oliveira/unitary-factorization/README.md), which explicitly records its unreviewed dependencies.

## Published hypotheses and proof map

Definition 5 (p. 863), for integers $m\ge k\ge3$:

$$\mathcal A(m,k)=\{(A,\operatorname{diag}b):A=A^*,\quad b_1,\ldots,b_{m-k}\in\mathbb R,\quad b_{m-k+1},\ldots,b_m\in\mathbb C\}.$$

Theorem 19 (pp. 866–867): ordinary inclusion throughout $\mathcal A(m,k)$ implies inclusion throughout $\mathcal A(n,k)$ for $n\ge m$. Multiplicities count; arbitrary entries may also be real. Unitary diagonalization covers normal $B$.

| Step | Checked justification / dependency |
| --- | --- |
| Induction | Base hypothesis; generic density and Observation 13. |
| (3.4) | Theorem 12; $n>m\ge k$ supplies real $b_1$. |
| (3.5)–(3.7) | Rational inversion; bounded remainders. |
| (3.8) | Vertex displacement bounds support functions. |
| (3.9) | Cofactor expansion at the divergent entry. |
| (3.10) | Compression preserves the class; induction. |
| Genericity paragraph | Birational coordinates below. |
| (3.11) | Lemma 15; Corollary 18; interior-ball scaling. |
| (3.12) | $(C+C')\varepsilon<d$ and (3.8). |

Our explicit genericity check uses the mutually inverse rational maps

$$C=(A+tI)^{-1},\quad\beta_j=(b_j-t)^{-1};\qquad A=C^{-1}-tI,\quad b_j=t+\beta_j^{-1}.$$

Keeping $t=b_1$, then selecting $(C',\beta)$, preserves $(n-1)^2+(n-1)+k$ algebraically independent coordinates. No arbitrary-compression principle is assumed.

No unmet hypothesis was found. Lemma 15's generic-interiority assertion remains a cited dependency, not a newly formalized theorem.

## Comparison with PR89

The two arguments share inversion and Hermitian interlacing. Their final mechanisms differ:

| Item | PR89's direct argument | Comparison result |
| --- | --- | --- |
| Lemma C | A random injection matches every squarefree product moment of interlacing sequences. | Its conditional induction preserves the remaining moments in both branches. |
| Lemma I | Positive reciprocals convert the injection law to convex membership. | The positivity condition is used correctly; it is not omitted. |
| Lemma M | A real pole shift makes $A+\beta I$ positive definite. | Both determinant range and all permutation vertices acquire the same nonzero factor. |
| Lemma E | A stationary support functional produces an adjugate identity. | Taking the skew-Hermitian part of the commutator gives the stated Hermitian condition; the adjoint signs are consistent. |
| Steps 2–5 | Compress along a real eigenvector; rescale by a positive cofactor ratio; contradict support separation. | The compressed normal matrix loses one real eigenvalue and retains at most $k$ nonreal entries. No generic-interior lemma is required. |

I found no concrete gap in these steps during this bounded reading. This is not a complete independent review of the four-dimensional base or all its referenced results. The primary paper already supplies the needed lifting theorem for $k=4$ conditional on the base, so PR89's direct lifting proof is an alternative route, not a necessary new hypothesis.

## Base case, attribution and numerical scope

At $m=k=4$, Definition 5 imposes no real entries on $B$. The required input is exactly: one Hermitian and one arbitrary normal $4\times4$ matrix satisfy MDO. The repository's Corollary P4′ states the line-spectrum condition with the other matrix named $B$; swapping the two matrices gives the required orientation. Its source labels Theorem F, Theorem O, and Theorem EQ dependencies as not independently reviewed. This source check does not upgrade them. Zheng's Theorem C.1 is cited by PR89 as another route and is not audited here.

Thus the citation-specific condition on H6 is resolved. In the repository's evidence terminology, its base-case review condition remains. For dimensions below four, the known $n\le3$ result applies; Theorem 19 supplies $n\ge4$. The paper's published three-nonreal corollary is not itself a published four-nonreal theorem.

PR89 already records 1,500 random $n=5$ examples and 12 numerical optimization runs, and 150 random $n=6$ examples. These are exploratory floating checks. Its $n=6$ archived log contains only the random-search line; no $n=6$ optimization completion is inferred.

Our separate small check uses seed 20261126, 100 cases in each of dimensions five and six, Haar-QR unitary mixing, real first spectrum, and exactly four nonreal entries in the second spectrum. It records every completed case and normalized hull margin in [numerical-report.json](numerical-report.json). No outside candidate crossed the declared $10^{-9}$ screening allowance. This does not prove the theorem or rule out a missed numerical violation.

```sh
python3 numerical.py --seed 20261126 --cases 100
```

Run from this directory with NumPy and SciPy. The report includes the actual runtime versions. The script needs no credentials and makes no network calls; output replaces the local report only when `--write` is supplied. Reproduction may vary in final floating digits across library versions. `SHA256SUMS` pins the inspected report and source.

The dispatchable next step is an independent-owner review of the explicit $n=4$ Hermitian base chain (and, separately if desired, PR89's stationary proof). No source-access blocker remains for Shitov's lifting statement.
