# Rank-one unitary phases: explicit common-minor certificates

**Written algebraic proof and exact regression tests; independent-owner review pending.** These results concern the relative unitary $U=I+(z-1)vv^*$, where $v^*v=1$ and $|z|=1$. The two spectra are arbitrary complex numbers. A common-minor certificate works simultaneously for every pair of spectra and is stronger than the original Marcus–de Oliveira inclusion. Its failure is **not** a counterexample to that inclusion.

This note gives explicit weights for a dimension-independent phase range, an exact triangle-packing criterion for every vector and phase, and sharp flat-vector thresholds. The pointwise Householder case is attributed to Drury's 1999 paper; no novelty claim is made for the present constructions. [Reference status](#reference-status-and-overlap)

## 1. Exact phase dependence

Put $p_i=|v_i|^2$, $c=|z-1|^2\in[0,4]$, $e_2=\sum_{i<j}p_ip_j$, and

$$P_\sigma=\prod_i(a_i+b_{\sigma(i)}).$$

In every dimension, including repeated or zero eigenvalues,

$$\boxed{\det(D_a+UD_bU^*)=(1-ce_2)P_{\mathrm{id}}
       +c\sum_{i<j}p_ip_jP_{(ij)}.}\tag{1}$$

For a rank-two derivation set $\alpha=z-1$, $\mu=\sum_i p_i b_i$, and $D=D_a+D_b$. Then

$$UD_bU^*-D_b=[v,D_bv]
 \begin{pmatrix}c\mu&\alpha\\\overline\alpha&0\end{pmatrix}
 \begin{bmatrix}v^*\\v^*D_b\end{bmatrix}.$$

When the entries of $D$ are nonzero, define $s_j=\sum_i p_i b_i^j/(a_i+b_i)$. The determinant lemma, using $\alpha+\overline\alpha=-c$, gives

$$\frac{\det(D_a+UD_bU^*)}{\det D}
 =1+c(\mu s_0-s_1-s_0s_2+s_1^2).$$

The pairwise identity

$$\frac{P_{(ij)}-P_{\mathrm{id}}}{P_{\mathrm{id}}}
 =\frac{(a_i-a_j)(b_i-b_j)}{(a_i+b_i)(a_j+b_j)}$$

turns the preceding expression into (1). Multiplying denominators gives a polynomial identity, so it remains true when any $a_i+b_i=0$. Neither $\mu$ nor $s_j$ assumes that $D_b$ is Hermitian.

Consequently the entire phase curve is a line segment. With $H=I-2vv^*$,

$$\det(D_a+UD_bU^*)=(1-c/4)P_{\mathrm{id}}
                 +(c/4)\det(D_a+HD_bH^*).\tag{2}$$

A pointwise Householder MOC theorem therefore implies the same inclusion for every phase. Equation (2) does not assert a common-minor certificate. The particular identity-and-swaps weights in (1) are nonnegative exactly when $ce_2\le1$.

## 2. Every squared minor

Let $I,J$ have the same cardinality and write $s_I=\sum_{i\in I}p_i$. Then

$$|\det U[I,J]|^2=
\begin{cases}
1-cs_I(1-s_I),&I=J,\\
cp_ip_j,&I\setminus J=\{i\},\ J\setminus I=\{j\},\\
0,&|I\setminus J|\ge2.
\end{cases}\tag{3}$$

To see this, each minor is a submatrix of the identity plus one rank-one term. By multilinearity, its determinant has degree at most one in $\alpha$: terms using the perturbation in two columns vanish. A principal minor is $1+\alpha s_I$. A one-exchange minor is, up to the ordering sign, $\alpha v_i\overline v_j$. With two or more exchanges, the unperturbed submatrix has rank at most $|I|-2$, so both its determinant and every first-order cofactor vanish. Taking squared moduli proves (3), including empty and full minors.

A common-minor certificate is a probability distribution on permutations satisfying

$$|\det U[I,J]|^2=\mathbb P\{\sigma(J)=I\}\quad\text{for every }I,J.$$

The usual determinant expansion then proves MOC simultaneously for all spectra. Our distributions are invariant under permutation inversion, so the column-to-row convention in this definition and the row-to-column convention in $P_\sigma$ give the same weights.

## 3. Explicit triangle weights

For every unordered triple assign $t_{ijk}\ge0$ as the weight of **each** of its two oriented 3-cycles. Set

$$s_{ij}=cp_ip_j-\sum_{k\ne i,j}t_{ijk},\qquad
w_0=1-ce_2+\sum_{i<j<k}t_{ijk}.\tag{4}$$

Use $s_{ij}$ for transposition $(ij)$ and $w_0$ for the identity. These weights sum to one. If they are nonnegative, they reproduce (3): for a transition replacing $j$ by $i$, a paired triangle containing $i,j$ contributes exactly one of its orientations, regardless of whether the third vertex is occupied. Only one-exchange transitions occur; diagonal probabilities follow by subtracting outgoing mass. At the determinant level the same substitution follows from

$$P_{(ij)}+P_{(ik)}+P_{(jk)}
 =P_{\mathrm{id}}+P_{(ijk)}+P_{(ikj)},$$

the vanishing $3\times3$ determinant of the rank-at-most-two matrix $(a_i+b_j)$, multiplied by the unchanged factors.

Choose $t_{ijk}=cp_ip_jp_k$. Then

$$s_{ij}=cp_ip_j(p_i+p_j),\qquad
w_0=1-\frac c3\left(1-\sum_i p_i^3\right),\tag{5}$$

where $e_2-e_3=(1-\sum_i p_i^3)/3$ was used. Thus:

> A common-minor certificate exists whenever $c(1-\sum_i p_i^3)\le3$. In particular, it exists in every dimension whenever $c\le3$, equivalently $|\arg z|\le2\pi/3$ for the principal argument.

For a positive support of size $m\ge3$, the explicit choice $t_{ijk}=\lambda cp_ip_jp_k$ improves this sufficient condition, where

$$\lambda=\frac1{1-p_{(1)}-p_{(2)}}$$

and $p_{(1)},p_{(2)}$ are the two smallest **positive** coordinates. All swap weights remain nonnegative, because $1-p_i-p_j\le1-p_{(1)}-p_{(2)}$. The remaining condition is

$$c(e_2-\lambda e_3)\le1.\tag{6}$$

Coordinates outside the support stay fixed. Support size one is trivial; support size two permits every phase by (1). Condition (6) is sufficient, not asserted necessary for a non-flat vector.

## 4. Exact finite criterion and its dual

For $c>0$, a common-minor certificate exists **if and only if** there are triangle weights with

$$t_{ijk}\ge0,\quad \sum_{k\ne i,j}t_{ijk}\le cp_ip_j,
\quad \sum_{i<j<k}t_{ijk}\ge ce_2-1.\tag{7}$$

Sufficiency is (4). For necessity, any positively weighted permutation must satisfy $|I\setminus\sigma(I)|\le1$ for every subset, since all two-exchange minors vanish. A permutation moving at least four elements has two moving edges with disjoint endpoints: choose one from each of two nontrivial cycles, or two separated edges of a cycle of length at least four. The set of their sources has two images outside it, contradicting (3). Thus only the identity, one transposition, or one 3-cycle can have positive weight.

Equation (3) is symmetric in $I,J$. Inverting every permutation therefore gives another certificate. Average the original distribution with its inverse; the two orientations of every triangle now have equal weights. Its singleton off-diagonal entries and total mass give exactly (4), hence (7).

Let $\tau(p)$ be the maximum total triangle weight with edge capacities $p_ip_j$. Then (7) is equivalent to

$$c\bigl(e_2-\tau(p)\bigr)\le1.\tag{8}$$

The triangle-packing LP is feasible and bounded. Its dual minimizes $\sum_{i<j}p_ip_jy_{ij}$ with $y_{ij}\ge0$ and $y_{ij}+y_{ik}+y_{jk}\ge1$. Truncating any $y_{ij}>1$ to one preserves feasibility and never increases the objective. Substituting $x_{ij}=1-y_{ij}$ yields the exact equivalent weighted fractional triangle-free problem

$$\begin{aligned}
\beta(p)&=\max\sum_{i<j}p_ip_jx_{ij},\\
0\le x_{ij}&\le1,\qquad x_{ij}+x_{ik}+x_{jk}\le2,
\end{aligned}\tag{9}$$

with $\beta(p)=e_2-\tau(p)$. The necessary-and-sufficient condition is $c\beta(p)\le1$. This is a finite LP characterization, not a closed-form formula for arbitrary $p$.

### A tempting shortcut fails

Constant edge weights $x_{ij}=2/3$ and cut indicators give

$$\beta(p)\ge\max\left\{\frac{1-\sum_i p_i^2}{3},\ 
\max_S p(S)(1-p(S))\right\}.$$

Equality need not hold. For

$$p=(4,2,6,1,19,8,3)/43,$$

the exact primal and dual witnesses in [packing-results.json](packing-results.json) give

$$e_2=679/1849,\quad \tau=649/5547,\quad
\beta=1388/5547,\quad c_{\max}=5547/1388<4.$$

The cut bound is $462/1849$ and the constant-edge bound is $1358/5547$; both are strictly below $\beta$. At the half-turn, the exact obstruction is $c\beta=5552/5547>1$, although every cut bound is below one and the expected number of moved points is below three. This is a failed simplification of the **common-minor** criterion. It does not challenge the pointwise Householder theorem.

## 5. Sharp flat threshold

If $p_i=1/m$ on its support, take $t_{ijk}=c/[m^2(m-2)]$ for $m\ge3$. All swap weights vanish and

$$w_0=1-\frac{c(m-1)}{3m}.$$

The exact threshold in the physical phase range is therefore

$$\boxed{c\le\min\{4,3m/(m-1)\}.}\tag{10}$$

For necessity, every permitted permutation moves at most three points. But the singleton minors prescribe expected moved-point count

$$\sum_i(1-|U_{ii}|^2)=c\left(1-\sum_i p_i^2\right)=c(m-1)/m,$$

which cannot exceed three. The formal threshold $9/2$ for $m=3$ lies beyond the actual range $c\le4$. At $m=4$ all phases work. At $m=5,c=4$, the expected count is $16/5>3$, so no common certificate exists. As $m\to\infty$, (10) approaches three, proving that the uniform constant in (5) cannot be increased over all dimensions for common-minor certificates.

## 6. Reproduction and evidence limits

From the repository root:

```sh
python3 research/marcus-de-oliveira/reflection-angles/check.py
python3 research/marcus-de-oliveira/reflection-angles/check_packing.py
```

Both use only Python's standard library and exact rational arithmetic. `check.py` independently implements Gaussian-rational determinants and verifies (1)–(5), phase interpolation, all compound entries for its selected fixtures, zero-support coordinates, singular $D_a+D_b$, the flat threshold formulas through dimension seven, and the permutation-support obstruction through dimension seven. It checks 33 unitary fixtures, 18,910 compound entries, and 5,710 permutations moving at least four coordinates. [Recorded results](results.json)

`check_packing.py` solves six bounded exact rational LP fixtures and **separately verifies** nonnegative primal weights, edge capacities, every dual triangle inequality, and equality of primal/dual objectives. The witnesses certify the stated finite optima without trusting floating-point optimization. Its report preserves the non-flat failed-shortcut example. [Exact witnesses](packing-results.json)

The finite checks support arithmetic and indexing; the written arguments establish the general statements. These are same-owner calculations, not independent-owner review, and no Lean theorem is claimed. No original-MOC counterexample or resolution follows.

## Reference status and overlap

S. W. Drury, [*OMC for Householder reflections*](https://www.sciencedirect.com/science/article/pii/S0024379599001731), *Linear Algebra and its Applications* **298** (1999), 159–169; [DOI](https://doi.org/10.1016/S0024-3795(99)00173-1). Its accessible primary abstract states the Householder MOC theorem without a dimension restriction and separately reports a dimension-five failure of the stronger Merikoski–Virtanen condition. **The full primary proof was not accessible and was not audited in this pass.** Equation (2) is the direct phase-reduction argument supplied here; the pointwise endpoint is attributed to Drury, not presented as a new result.

The [2025 review](https://www.mdpi.com/2227-7390/13/5/711) mentions the dimension-four Householder case. That narrower statement does not establish that higher-dimensional pointwise Householder cases are open. Precise overlap of the triangle construction with Drury's full paper remains to be checked. The dimension-five flat obstruction is consistent with his reported stronger-conjecture failure and is presented as a reproduction, not a new counterexample.

Within this repository, [the search](../search/README.md) already samples near-monomial rank-one phases numerically, [the forest proof](../baselines/README.md) constructs common weights for forest products of two-coordinate gates, and [the minors report](../minors/report.md) discusses general common-minor certificates. Focused text inspection found no existing explicit reflection-angle weights or triangle-packing criterion. This is not an exhaustive novelty search. The next bounded questions are an audit of Drury's full proof and a useful closed-form evaluation of (9) beyond the flat case; the trace-plus-cut shortcut is already ruled out by the exact witness above.
