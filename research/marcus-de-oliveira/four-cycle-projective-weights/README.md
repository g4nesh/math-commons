# Four-cycle weights on the exterior of the second-eigenvalue plane

Fix the real orthogonal matrix and the opposing spectrum from the [four-cycle fixture](../full-four-cycle/README.md),

$$
U=\frac1{625}\begin{pmatrix}-31&492&-240&300\\-300&225&500&0\\240&-180&225&500\\-492&-256&-180&225\end{pmatrix},\qquad
b=(3-2i,2,1+2i,-i).
$$

Let $B=U\operatorname{diag}(b)U^*$ and $a=(x,\eta,y,z)$. For every complex $\eta$ with $|\eta|\ge384$ and every complex $z$, there are 24 nonnegative weights of sum one such that
$$
\det(\operatorname{diag}(a)+B)=\sum_p \lambda_p(\eta,z)\prod_{j=1}^4(a_{p(j)}+b_j)
$$
**simultaneously for all complex $x,y$.** Eighteen weights are strictly positive. The other six, indexed by $p(4)=2$ in one-based notation, vanish. This extends the archived fixture into an unbounded region, alongside the [finite rectangle](../complex-second-region/README.md). The intervening region and the general conjecture remain open.

## A chart at infinity

Write $z=s+it$ and $D=1+s^2+t^2$. Give each of the 18 allowed permutations a real quadratic numerator in $s,t$. Stack its six coefficients into a vector, obtaining 108 coordinates in all. For fixed $\eta$, the desired polynomial identity before normalization is the vanishing of the weighted sum of the differences between each pairing and the determinant. Equating the coefficients of $1,x,y,xy$, then the ten monomials of total degree at most three in $s,t$, gives four complex blocks of ten equations. Stack their real parts above their imaginary parts. Call the resulting $80\times108$ matrix $C(\eta)$. Its entries are reconstructed from the original principal minors.

Set $A=C(1)-C(0)$ and $B_0=C(0)$. Let $J$ send the real/imaginary row pair $(r,i)$ to $(-i,r)$. For $w=1/\eta=u+iv$, multiplying every complex equation by $w$ gives the real matrix
$$M(w)=A+uB_0+vJB_0.$$
This follows from complex affine dependence of each unstacked residual on $\eta$. For $w\ne0$ this operation preserves rank and kernel. The [26 polynomial row relations](../complex-second-region/README.md#global-rank-control-proved-by-polynomial-identities) give $\operatorname{rank} C(\eta)\le54$ for every complex $\eta$. Every size-55 minor of $M$ also vanishes at zero by continuity. The checker repeats the exact polynomial relation test.

## A positive center and its correction

For each allowed permutation start with one copy of every square $(a+bs+ct)^2$ with a nonzero triple in $\{-1,0,1\}^3$, taking one representative of each pair of opposite signs. The 13 squares sum to $9D$. The 40 additional rational masses in [center.json](center.json) are positive. Expand all squares to form a vector $c$. Exact multiplication gives $Ac=0$, and every allowed numerator is at least $9D$.

Let $R,K$ be the 54 chosen rows and columns in that file. The checker inverts $P=A_{R,K}$ over the rationals and computes
$$
X=P^{-1}(B_0)_{R,K},\quad Y=P^{-1}(JB_0)_{R,K},\quad
f=P^{-1}(B_0c)_R,\quad h=P^{-1}((JB_0)c)_R .
$$
The induced maximum row-sum norms satisfy $\|X\|_\infty+\|Y\|_\infty<145$. For the two vectors, $\|f\|_\infty+\|h\|_\infty<814$.

If $|w|\le1/384$, then $I+uX+vY$ is invertible by the geometric series. Add
$$d=-(I+uX+vY)^{-1}(uf+vh)$$
to the entries of $c$ indexed by $K$. The selected rows now vanish. They span all rows since their pivot stays invertible and the full matrix has rank at most 54. The bound
$$\|d\|_\infty < 814/(384-145)<7/2$$
controls the correction to every polynomial coefficient. The sum of the absolute values of the six monomials is at most $2D$, by $2|s|\le1+s^2$, $2|t|\le1+t^2$ and $2|st|\le s^2+t^2$. Every allowed numerator therefore remains at least $2D$. Divide by the sum of the numerators to obtain the weights. This argument includes the circle $|\eta|=384$ and has no bound on $z$.

## Reproduction and handoff

From the repository root run

```sh
python3 -I -B -O research/marcus-de-oliveira/four-cycle-projective-weights/check.py
```

The checker uses Python's standard library. It reconstructs the center, verifies the exact rank identities and both perturbation bounds, and compares the returned weights against four original determinants. The inequalities prove the continuum statement. To evaluate rational inputs directly, import `weights_at(second_real, second_imag, fourth_real, fourth_imag)` from [check.py](check.py). Pass integers or `fractions.Fraction` values. Results follow lexicographic column-to-row permutation order.

A follow-up is to construct positive centers that join this exterior to the [finite rectangle](../complex-second-region/README.md). Weights may depend on both moving entries. A failed choice of numerator degree only rules out that choice if infeasibility is itself certified.
