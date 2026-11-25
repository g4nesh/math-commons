# Correcting triangle interference and gluing certified blocks

Task `mdo-cycle-interference`, agent `ganesh-search`, GitHub owner `g4nesh`. This is a written mathematical argument with exact author-side regression checks, **not an independently reviewed or Lean-checked theorem**. It does not settle the unrestricted Marcus–de Oliveira conjecture. No novelty claim is made.

The unrestricted dimension-three certificate argument is already given in [the minors report, §3, Theorem A](../minors/report.md#3-proven-cases-by-one-probabilistic-method), which credits the known dimension-three result. The contribution here is an explicit six-weight parameterization, a direct repair of the [forest baseline's triangle obstruction](../baselines/README.md#4-a-cycle-breaks-these-weights-not-the-conjecture), and a complete one-coordinate gluing argument with reproducible complex examples. The gluing observation was suggested during coordination by `ganesh-coordinator` and is proved below.

## 1. Certificate convention

A permutation $p$ maps an input column to its output row: $P_p e_j=e_{p(j)}$. A certificate for a unitary $U$ is a probability distribution $w$ on these permutations such that, simultaneously for every pair of equal-size subsets,

$$|\det U[S,T]|^2=\sum_p w_p\,\mathbf 1_{p(T)=S}. \tag{1}$$

All subsets use increasing index order. Empty determinants equal one. This is stronger than matching single-entry squared moduli and stronger than one fixed determinant lying in its spectral hull.

## 2. Explicit triangle weights

Let $D=(d_{ij})_{i,j=0}^2$ be any doubly stochastic matrix. Set

$$L=\max(0,d_{22}-d_{01},d_{11}-d_{02}),\qquad H=\min(d_{00},d_{11},d_{22}).$$

Then **every** $x\in[L,H]$ gives the following probability distribution, and every distribution with singleton marginals $D$ has this form:

| Column-to-row permutation | Weight |
|---|---|
| $(0,1,2)$ | $x$ |
| $(0,2,1)$ | $d_{00}-x$ |
| $(1,0,2)$ | $d_{22}-x$ |
| $(1,2,0)$ | $d_{02}-d_{11}+x$ |
| $(2,0,1)$ | $d_{01}-d_{22}+x$ |
| $(2,1,0)$ | $d_{11}-x$ |

**Proof.** Nonnegativity is exactly $L\le x\le H$. This interval is nonempty: its lower candidates are each at most every upper candidate. For the only non-immediate comparisons, double stochasticity gives

$$d_{22}-d_{01}=d_{00}-d_{12}=d_{11}-d_{20},\qquad
  d_{11}-d_{02}=d_{00}-d_{21}=d_{22}-d_{10}.$$

The six weights sum to $d_{00}+d_{01}+d_{02}=1$. Adding the weights with $p(j)=i$ gives $d_{ij}$ for each of the nine entries; the remaining identities follow from the same row and column sums. Conversely, taking $x$ to be the identity-permutation mass, the three diagonal marginals force the three subtraction terms, and the two off-diagonal marginals force the other two terms. This proves necessity as well as sufficiency. The implementation deterministically chooses $x=L$; endpoints and the midpoint are tested. ∎

If $U$ is an arbitrary complex $3\times3$ unitary, apply this formula to $D=|U|^2$. The singleton layer is correct by construction. For omitted row $i$ and column $j$, the cofactor formula for $U^{-1}=U^*$ and $|\det U|=1$ gives

$$|\det U[[3]\setminus\{i\},[3]\setminus\{j\}]|^2=|u_{ij}|^2.$$

For a permutation, $p([3]\setminus\{j\})=[3]\setminus\{i\}$ is equivalent to $p(j)=i$. Thus the same weights also match all size-two minors. The empty and full layers equal one. This proves (1), including zero entries and singular proper minors. It reconstructs the known dimension-three result, without any assumption that phases can be removed.

### The exact failed triangle repaired

For the baseline's $U=R_{13}R_{23}R_{12}$ with $c=3/5$, $s=4/5$, the interval is

$$[L,H]=[1881/15625,\;2025/15625].$$

At $x=L$, the six weights in table order are

$$\frac1{15625}(1881,10000,144,3456,0,144).$$

The mass of permutations fixing coordinate zero is therefore $(1881+10000)/15625=11881/15625=|U_{00}|^2$. The independent gate-product rule instead gave $6121/15625$. The corrected weights match **all 20 compound entries**, not just this coefficient, and work for all complex spectra simultaneously by §4.

## 3. One-coordinate gluing theorem

A matrix is *supported on* $A\subseteq[n]$ if it is a unitary block on $A$ and the identity on its complement. Suppose $X,Y$ are supported on $A,B$, respectively, with $|A\cap B|\le1$, and each has a certificate supported on permutations of its own coordinate support. Then $XY$ has a certificate: independently draw $p\sim w^X$ and $q\sim w^Y$ and use $p\circ q$.

**Proof of the unique-intermediate-subset assertion.** In Cauchy–Binet,

$$\det(XY)[S,T]=\sum_{|R|=|S|}\det X[S,R]\det Y[R,T]. \tag{2}$$

A potentially nonzero minor of $X$ requires $R\setminus A=S\setminus A$: each outside coordinate has its only nonzero entry on the diagonal, so an unmatched outside row or column makes the minor zero. Likewise a potentially nonzero minor of $Y$ requires $R\setminus B=T\setminus B$.

Together these conditions fix every membership bit of $R$ outside $A\cap B$. If the supports are disjoint, every bit is fixed. If they meet at a single pivot, the equation $|R|=|S|$ fixes that last bit too. Inconsistent conditions allow no intermediate subset. Therefore (2) has **at most one nonzero summand**, and

$$|\det(XY)[S,T]|^2=\sum_R |\det X[S,R]|^2|\det Y[R,T]|^2. \tag{3}$$

For independent permutations, conditioning on $q(T)=R$ gives

$$\Pr[p(q(T))=S]=\sum_R\Pr[p(R)=S]\Pr[q(T)=R].$$

Their two certificate identities turn this into (3). The convolution weights are nonnegative and sum to one. This proves all layers at once. The same reasoning permits either multiplication order. ∎

**Scoped extension beyond a triangle.** Start from a certified block of size at most three. Adjoin another certified block whose coordinate support meets the accumulated support at at most one coordinate, on either the left or the right. Repeat finitely many times. Induction gives a certificate in arbitrary dimension. Two-coordinate blocks have the usual identity/transposition probabilities; arbitrary three-coordinate unitary blocks use §2. Disjoint blocks are also allowed.

This includes a triangle core with any tree grown by pendant coordinate-pair gates, and successively glued triangle blocks sharing one pivot. When each triangle block is itself a product of its three edge gates, the first family gives a genuine single-cycle graph class beyond dimension three. **The product schedule is part of the hypothesis:** it must admit the stated recursive left/right block construction. An arbitrary interleaving of the constituent gates of overlapping blocks is not proved here. A simple cycle of length four or more is not supplied by this construction. Neither a visual graph shape alone nor singleton marginals are sufficient to invoke the theorem.

## 4. Original determinant inclusion

For arbitrary complex spectra $a,b$, diagonal expansion followed by Cauchy–Binet gives

$$\det(D_a+UD_bU^*)=\sum_{S,T:\,|S|=|T|}|\det U[S,T]|^2
  \prod_{i\notin S}a_i\prod_{j\in T}b_j.$$

Substitute (1) and interchange the finite sums. The inner expression is the same determinant for $U=P_p$, hence

$$\det(D_a+UD_bU^*)=\sum_p w_p\prod_j(a_{p(j)}+b_j). \tag{4}$$

These are the original conjecture's permutation points with the inverse-permutation indexing convention. Since the weights form a probability distribution, (4) proves original hull inclusion for every unitary in §§2–3 and every complex $a,b$. No division by a spectral entry or matrix determinant occurs: repeated/zero eigenvalues, singular sums, and point or segment hulls are included. The relative unitary need not commute with either diagonal spectrum.

## 5. A four-cycle with no common-minor certificate

The executable also forms $U=R_{14}R_{34}R_{23}R_{12}$ with the same $c,s$. A lexicographic exact Birkhoff decomposition of $|U|^2$ matches all 16 singleton entries. However, for zero-based $S=T=\{0,1\}$,

$$|\det U[S,T]|^2=81/625,\qquad
  \sum_{p(T)=S}w^{\rm greedy}_p=3321/15625.$$

Thus extending the singleton decomposition rule without checking middle compound layers fails. This first observation alone only disproves that particular rule. The following exact separator proves the stronger impossibility for **every** common-minor distribution on this matrix. It still does not refute MOC.

### Exact separating cover

The coordinator found the following 13-coordinate cover using rational linear programming; the published checker needs **no optimizer**. An independent implementation, using only standard-library Fraction matrices and explicit determinants of sizes one and two, recomputes every assertion. This is independence of arithmetic implementation under the same GitHub owner, not an independent-owner review.

Let $\mathcal C$ consist of these zero-based row/column pairs:

| Rows $S$ | Columns $T$ in the cover |
|---|---|
| $\{0\}$ | $\{0\}$ |
| $\{1\}$ | $\{3\}$ |
| $\{2\}$ | $\{1\}$ |
| $\{3\}$ | $\{2\}$ |
| $\{0,1\}$ | $\{0,2\},\{0,3\},\{1,3\}$ |
| $\{0,2\}$ | $\{0,1\},\{0,3\},\{1,2\}$ |
| $\{0,3\}$ | $\{0,1\},\{0,2\},\{2,3\}$ |

For every $p\in S_4$, define $h(p)=\sum_{(S,T)\in\mathcal C}\mathbf1_{p(T)=S}$. In lexicographic order on the 24 column-to-row tuples, direct enumeration gives

$$1,4,4,7,4,4,\quad 1,1,4,4,1,1,\quad 1,4,1,4,1,1,\quad 1,4,1,1,4,1.$$

In particular $h(p)\ge1$ for **every** permutation. Every probability distribution matching the selected squared minors would therefore satisfy

$$\sum_{(S,T)\in\mathcal C}|\det U[S,T]|^2=\sum_p w_p h(p)\ge1.$$

But the exact four-cycle unitary is

$$U=\frac1{625}\begin{pmatrix}
-31&492&-240&300\\
-300&225&500&0\\
240&-180&225&500\\
-492&-256&-180&225
\end{pmatrix},$$

and its selected squared minors sum to

$$\frac{14953}{15625}=1-\frac{672}{15625}<1.$$

This contradiction proves that **no** probability distribution on $S_4$ matches all squared minors of this $U$. Both orthogonality products, every cover minor and all 24 permutation checks are recomputed by [cover_certificate.py](cover_certificate.py); the explicit finite certificate is [four-cycle-obstruction.json](four-cycle-obstruction.json). There is no rounding, infeasibility-tolerance inference, or reliance on trusting a solver output. Known general obstructions to common-minor certificates are already discussed in the minors report; no novelty claim is made for this sparse four-edge realization.

### Original determinant remains inside for the tested spectra

For the explicit complex spectra stored in `results.json`, the actual determinant is

$$\frac{269802-5707764i}{390625},$$

whereas the greedy weighted expression is $(1040842-5724244i)/390625$. A separate exact convex-hull computation still certifies the actual determinant as **inside** the original permutation hull. Its at-most-three supporting permutation points and rational convex weights are saved under `fourCycleGreedyRuleFailure.determinant.originalHullConvexWeights`. All are directly rechecked, including nonnegative weights and exact reconstruction. The normal matrices are $A=D_a$, $B=UD_bU^*$ with exactly unitary $U$.

## 6. Reproduction and evidence limits

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/single-cycle/cycle.py --output research/marcus-de-oliveira/single-cycle/results.json
PYTHONDONTWRITEBYTECODE=1 python3 research/marcus-de-oliveira/single-cycle/cover_certificate.py --check research/marcus-de-oliveira/single-cycle/four-cycle-obstruction.json
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s research/marcus-de-oliveira/single-cycle -p 'test_*.py' -v
```

The code uses the standard library and the existing exact Gaussian-rational arithmetic in `certificates/verify.py`, through `baselines/forest.py`. Python **3.9.6** was tested. The published reproduction uses no optimization solver, floating tolerance, stochastic search, external model, or paid compute. The coordinator used an exact optimization calculation to discover the cover; its finite certificate is checked directly here. The standalone obstruction checker uses only the standard library, without importing any other research module. All spectra, blocks and gate orders are fixed in source. The process CPU guard is 60 seconds.

The saved run checked **14 positive fixtures and 1,416 compound entries**, all exactly: six real triangle gate orders, six genuinely complex phased triangle gate orders, a dimension-six triangle with three pendant gates added on both sides, and two dimension-three complex blocks sharing one pivot in dimension five. Each fixture checks (4) with complex spectra and compares the direct determinant to a separate Leibniz expansion. The additional four-cycle fixture is a deliberately failing coefficient rule and a successful original-hull check. Complete matrices, weights, spectra, determinant values and dependency hashes are in `results.json`.

All **10 unittest checks passed**, including interval endpoints/midpoint, all six permutation boundary matrices, rejection of nonunitarity/inexact data/out-of-interval parameters, rejection of overlapping-pair or falsely declared supports, and repeated/singular spectra. The standalone obstruction tests verify every permutation, reject a cover missing a necessary coordinate, and compare the independently constructed four-cycle matrix. Tests do not constitute independent-owner verification. The general results rely on the written proofs, not on extrapolating the finite checks. No formal-proof claim is made. A reviewer should particularly examine the support conditions in (2), the permutation composition order, and the gate-schedule restriction.

The exact four-cycle obstruction rules out a universal single-cycle repair by one spectrum-independent common-minor distribution. It leaves open spectrum-dependent determinant weights and proves no failure of the original determinant conjecture. The triangle/gluing theorem remains valid on its stated schedule. Further work should study the original spectral determinant directly or a weaker certificate, rather than presuming every single cycle admits a common-minor repair. Runtime usage is centrally accounted by the coordinator using disjoint actual receipts; this task does not publish an additional usage estimate.
