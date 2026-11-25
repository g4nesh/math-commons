# Literature audit and exact finite checks

Task `mdo-literature` · agent `ganesh-literature` · GitHub owner `g4nesh` · November 24, 2026 (America/Phoenix), November 25 UTC.

**Status: research artifact awaiting review by another owner.** This pass establishes neither the general Marcus–de Oliveira conjecture nor the claimed four-dimensional proof. It identifies a reproducibility gap in that proof's accessible materials and checks two finite examples with exact arithmetic. The [objective brief](../README.md) states the original problem.

## Evidence and access

| Primary source | Precise location and relevant scope | Evidence obtained here |
| --- | --- | --- |
| [Zheng, *The Marcus–de Oliveira conjecture for n = 4*](https://doi.org/10.13140/RG.2.2.36790.25923), October 2026 preprint | Theorem 1.1 claims the full normal $4\times4$ case; Theorem 1.2 represents all squared minors using two orthogonal endpoints and one common coefficient. | All 66 pages readable through the web PDF reader. Direct byte download returned HTTP 403. No PDF hash or independent proof certification. |
| [Kovačec, DMUC p2611](https://www.mat.uc.pt/preprints/ps/p2611.pdf), **v2, November 4, 2026** | Current title: *A large polygonal region contained in a simply connected Marcus de Oliveira region*. The unnumbered theorem on p. 5 assumes simple connectedness and a Jordan boundary; its conclusion is an inner-region inclusion. The former stronger outer-enclosure conjecture has counterexamples. | Full 10-page PDF downloaded and hashed; page 4 rendered to check notation. This is not the April version indexed by some searches. |
| [Shitov, *A further sufficient condition for the determinantal conjecture*](https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=9292&wshow=paper), 2025 | Theorem 1: $A$ Hermitian, $B$ normal with at most three nonreal eigenvalues. Theorem 19 lifts a fixed $\mathcal A(m,k)$ case to $\mathcal A(n,k)$, $n\ge m$. | Full English HTML including proof read. English publication date August 25, 2025; receipt date November 21, 2021. PDF download timed out; no byte hash. |
| [Drury, 2006](https://journals.uwyo.edu/index.php/ela/article/view/403) | Hermitian $A$, with the other matrix a nonreal scalar direction times a Hermitian matrix. Translate the paper's subtraction convention using $B\mapsto-B$. | Primary article and [full PDF](https://ftp.gwdg.de/pub/EMIS/journals/ELA/ela-articles/articles/vol15_pp285-296.pdf) screened; no full proof audit. |
| [Rodtes, arXiv:2006.14846v1](https://arxiv.org/abs/2006.14846v1), 2020 | Particular normal block-matrix families and direct-sum closure. | [Full text](https://arxiv.org/html/2006.14846v1) screened. This does not establish the arbitrary normal-pair case. |
| [Bebiano and da Providência, 2025](https://www.mdpi.com/2227-7390/13/5/711) | Survey and established special-case context. | Read as a map to primary results, not evidence about developments after February 22, 2025. |
| [Mulero-Martínez, 2026](https://doi.org/10.1016/j.laa.2026.03.019) | Variational analysis of the conjecture; structured positive cases and limitations of first-order arguments. | Publisher abstract only; full text was inaccessible in this pass. No theorem-level use in our derivations. |

The machine-readable [source manifest](sources.json) distinguishes access from verification. Kovačec v2 acknowledges Zheng's claim tentatively; neither that acknowledgment nor the claim itself substitutes for checking the proof. No source screened here establishes a general resolution.

## Four-dimensional proof: the concrete audit blocker

[Zheng's full preprint](https://www.researchgate.net/profile/Dongzhe-Zheng/publication/416773999_The_Marcus-de_Oliveira_conjecture_for_n_4/links/6add31918aa8df110de33813/The-Marcus-de-Oliveira-conjecture-for-n-4.pdf), §4.3 and Appendix F, sends 90 contact-case proofs to supplementary material: S0 notation, Table S1 proof index, S2 programs/data. Section 5 develops T85, using T40/T55/T81. Appendix G gives a dependency order. The reported checks include 82 generators, 6846 simplices and 214 contact orbits, over integers/rationals. Reported environment: Python 3.14.2, SymPy 1.14.0, python-flint 0.9.0 / FLINT 3.6.0. We did not locate these supplements in the linked landing page, author site or targeted searches; this is a bounded search result, not proof of absence.

The following are **proposed audit tasks**, not findings that the argument is wrong:

1. Obtain a stable complete source package and pin its hashes. Map every claimed certificate to a proposition, inputs and expected output. Do not accept screenshots of successful runs as a certificate inventory.
2. Separate the two mathematical obligations: unitary-to-orthogonal reduction and the orthogonal spectral inequality. The first cannot be replaced by a common-permutation-weight assertion; the exact obstruction below forbids that shortcut.
3. Check genericity, limiting cases and the endpoints in the geometric reduction, then inspect every pole and boundary used in contact accumulation. Preserve spectral realizability; enlarging the parameter domain could change the theorem.
4. Replay the finite classification independently, including coverage and dependency order. An exact identity checker only establishes the identity it is given; it does not establish that all cases or parameter regions were supplied.

Until those obligations are met, use **“claimed $n=4$ proof; audit pending”**, not “verified $n=4$” or a blanket “no $n=4$ proof has been claimed.” We have not established peer-reviewed acceptance or a complete formal verification.

## Known common-weight obstruction, reproduced

The input matrix and functional are [Zheng, Example 1.3, p. 3](https://www.researchgate.net/profile/Dongzhe-Zheng/publication/416773999_The_Marcus-de_Oliveira_conjecture_for_n_4/links/6add31918aa8df110de33813/The-Marcus-de-Oliveira-conjecture-for-n-4.pdf). This is a reproduction of a known obstruction, not a new counterexample.

$$
O_*={1\over7}\begin{pmatrix}
1&4&4&4\\-4&1&4&-4\\-4&-4&1&4\\-4&4&-4&1
\end{pmatrix}.
$$

For a matrix $V$, with $V[I,J]$ denoting its indicated minor, set

$$
2\psi(V)=\sum_i|V_{ii}|^2
-|\det V[12,13]|^2-|\det V[12,24]|^2
+|\det V[13,12]|^2+|\det V[13,34]|^2.
$$

Our independent standard-library calculation verifies $O_*^TO_*=I$, $\det O_*=1$, and $\psi(O_*)=-158/2401$. Exhausting all 24 permutation matrices gives $\psi(P)\in\{0,1,2\}$, with multiplicities 13, 10, 1.

Here is the separating-functional argument. Let $r(V)$ be the vector of squared absolute values of all minors. The displayed functional is linear in $r(V)$. If nonnegative weights $w_\sigma$, summing to one, represented **all coordinates** of $r(O_*)$ as $\sum_\sigma w_\sigma r(P_\sigma)$, linearity would give

$$
-158/2401=\psi(O_*)=\sum_\sigma w_\sigma\psi(P_\sigma)\ge0,
$$

a contradiction. This rules out those simultaneous weights. For a fixed pair of spectra, MOC asks only for a convex representation of one complex determinant; its weights may depend on those spectra. The contradiction above therefore does not refute MOC. It also does not refute a representation whose endpoints are arbitrary orthogonal matrices.

## The stronger-polygon witness stays inside the original hull

[Kovačec v2, §4, pp. 8–9](https://www.mat.uc.pt/preprints/ps/p2611.pdf) records Goyal's example

$$
a=(0,1,i),\quad b=(-1-i,1,-1+i),\quad
U={1\over3}\begin{pmatrix}-1&2&2\\2&-1&2\\2&2&-1\end{pmatrix}.
$$

Our exact computation gives $U^*U=I$ and

$$
\det(D_a+UD_bU^*)={2\over3}+{2\over9}i.
$$

In lexicographic permutation order $123,132,213,231,312,321$, the original comparison points are

$$
(6-2i,\;2,\;2+i,\;-i,\;2i,\;2-2i).
$$

The independently supplied weights

$$
(0,\;0,\;1/3,\;13/27,\;5/27,\;0)
$$

are nonnegative, sum to one, and recover the determinant exactly. This certifies membership in the original MOC hull for this witness. We have **not** reproduced the smaller polygon's boundary/topology; its claimed failure is recorded as the source's result. Confusing the two regions would turn a valid stronger-conjecture counterexample into an invalid claim about MOC.

There is also a local displayed-formula discrepancy in v2 Corollary 3, p. 4: its outside product retains $b_\nu$ after permuting the diagonal, where direct expansion requires $b_{\pi(\nu)}$. With $a=(1,2,3),b=(4,5,6),\pi=(3,2,1),i=1,j=2,G=I$, the actual determinant is 343; the displayed unpermuted outside factor gives 441, and the corrected factor gives 343. The adjacent segment statement survives this correction. This is a notation/algebra correction, not a refutation of that result or MOC.

## A conditional route worth separating from the full conjecture

[Shitov, Definition 5 and Theorem 19](https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=9292&wshow=paper) concerns $\mathcal A(m,k)$: Hermitian $A$, diagonal $B$, with at most $k$ nonreal entries. A proof for every pair in $\mathcal A(m,k)$ lifts to every $\mathcal A(n,k)$ with $n\ge m$. We read the theorem and proof structure but did not independently certify the full argument.

**Our conditional inference:** setting $m=k=4$, a valid four-dimensional proof would extend the case where $A$ is Hermitian and $B$ has at most four nonreal eigenvalues to every $n\ge4$. In fact, that four-dimensional Hermitian subcase alone would suffice. This suggests a smaller research target even if the arbitrary $4\times4$ proof remains under audit. It does not extend the full arbitrary-normal-pair conjecture to higher dimensions.

## Reproduction and handoff

From the repository root:

```sh
python3 research/marcus-de-oliveira/literature/check_source_examples.py
```

The run used Python 3.9.6 with only `fractions`, `itertools`, `dataclasses` and `json`; all arithmetic and comparisons are exact. [Recorded output](source-example-results.json) includes the three checks. No downloaded program was executed. These checks certify finite inputs and the stated elementary contradiction; they do not certify the source papers' universal theorems.

Priority handoff: first recover and pin the complete $n=4$ supplementary package; separately have another owner rerun these finite checks and inspect the distinction between common weights and a single spectral evaluation. A geometry specialist can audit the unitary reduction while a different worker checks classification coverage. Preserve finite budgets and checkpoint incomplete audits explicitly. No independent-owner review has happened in this artifact.
