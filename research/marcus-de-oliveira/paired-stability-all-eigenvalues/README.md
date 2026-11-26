# Breaking all the eigenvalue pairs

In dimension six, choose a base pair of spectra with three double eigenvalues on each side. Suppose their 21 indexed permutation products form a polygon with positive area. **For an open dense set of relative unitaries, the determinant inclusion persists when all twelve eigenvalues and the relative unitary move independently by sufficiently small amounts.** The neighborhood may depend on the center. It contains pairs with six distinct eigenvalues in each matrix. A computable lower bound for the radius follows below.

The prerequisite is the [universal paired-spectrum law](../double-paired-five-term/GRASSMANN.md), a written proof with an exact polynomial certificate submitted for independent review. This extension does not establish the conjecture for arbitrary simple spectra.

## A strictly positive common distribution

Let $q(U)$ be the 21 canonical signed coefficients from the [dual classification](../double-paired-contingencies/README.md). For each of the 189 extreme dual rays $y$, the preceding theorem gives $f_y(U)=y\cdot q(U)\ge0$. This is a real analytic function on connected $U(6)$. It is positive at some permutation matrix: choose a table with $y_j>0$, then a permutation realizing that table. At the corresponding matrix the canonical coefficient vector puts all its mass there. Each zero set therefore has empty interior. The finite intersection of the strictly positive sets is open and dense and is independent of the spectral values.

At any such $U$ take
$$
\alpha=\min\left\{\frac1{42},\frac12\min_y\frac{y\cdot q}{\sum_jy_j}\right\}>0.
$$
The normalized vector $(q-\alpha\mathbf 1)/(1-21\alpha)$ satisfies every dual inequality. The kernel repair from the classification gives it nonnegative weights. Scale those weights by $1-21\alpha$ and add $\alpha$ at each of the 21 tables. The result represents the original determinant for every paired spectrum, with every weight at least $\alpha$.

## A quantitative margin and radius

Write $v_j$ for the 21 product points at the center and $z$ for its determinant. Choose any triangle of these points with doubled area $A>0$ and longest side at most $D$. In every unit direction the triangle has width at least $A/D$. One elementary proof takes the vertices of least and greatest projection onto the perpendicular axis. Twice the area is their projection separation, which is at most $D$, times the distance in the support coordinate between the third vertex and the interpolated opposite side, which is at most the width. Thus the whole product hull has width at least $A/D$. In each direction, the gap between the support of the hull and the projection of $z$ is a sum of nonnegative weighted deficits. The opposite extreme contributes at least $\alpha$ times the width. A uniform support margin is therefore $m=\alpha A/D$.

Let $R\ge \max_i|a_i|+\max_j|b_j|$ and $B\ge\max_j|b_j|$. Perturb each of the twelve individual complex eigenvalues by at most $\delta\le1$ in modulus. Allow the relative unitary to vary by at most $\delta$ in operator norm. Each genuine six-factor product changes by at most $12(R+2)^5\delta$ by telescoping. Every indexed product at the center still has at least one actual six-coordinate permutation carrying it after the split.

For the matrix sum, the norm change is at most $(2+2B)\delta$, by separating the two diagonal perturbations and the two unitary factors. Its old and new operator norms are at most $R+2$. Replacing columns one at a time and using Hadamard's inequality bounds the determinant change by $6(R+2)^5(2+2B)\delta$. The largest possible loss of support clearance is consequently at most
$$
12(2+B)(R+2)^5\delta.
$$
Any $\delta\le1$ for which this is smaller than $m$ works for **every** simultaneous perturbation described above. The 21 weights are used at the center. The split matrices may have different convex weights.

## Exact example with two fully simple spectra

Take base triples $(6-i,1+6i,-6-8i)$ and $(-6-8i,-1+2i,3-i)$. These are also the spectra for the [failure of the elementary 45-generator approach](../paired-pointwise-cone/README.md). Use the relative unitary $U=H_3H_2H_1$ given by the rational Householder reflections $H_j=I-2v_jv_j^*/(v_j^*v_j)$ with

$$
v_1=(1,1,i,1,1+i,-i)^T,\quad
v_2=(1+i,-1,2,i,1,1)^T,\quad
v_3=(1,2i,-1,2+i,1,-2)^T.
$$

Exact arithmetic gives $\alpha=277/94080$, support margin at least $49072489/16611$, and $R=20$, $B=10$. The entire ball of radius **$1/262144$** satisfies the strict margin inequality. Independently splitting all repeated values along distinct offsets and rotating the relative unitary by a small complex gate gives an explicit pair in this ball. The checker directly computes its determinant and verifies strict barycentric weights on three actual permutation products.

Run `python3 -I -B research/marcus-de-oliveira/paired-stability-all-eigenvalues/check.py`. The [checker](check.py) rebuilds the dual rays, strict common law, triangle width estimate, perturbation bound and the exact split example. These are local neighborhoods in the unrestricted spectral problem. The pointwise statement about simple spectra is separate from the question of one common law that works for every possible splitting.
