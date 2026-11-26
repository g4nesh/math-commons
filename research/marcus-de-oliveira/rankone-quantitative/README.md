# An explicit complex square with pointwise inclusion and no common law

For the archived four-cycle unitary $U$ and spectrum $b=(3-2i,2,1+2i,-i)$, put

$$A(x,y,z)=\operatorname{diag}(x,y,2-i/4,z),\qquad B=U\operatorname{diag}(b)U^*.$$

**For every $x,y\in\mathbb C$ and every $z$ in the closed square**

$$|\Re z|\le 1/32768,\qquad |\Im z|\le 1/32768,$$

**the original Marcus–de Oliveira determinant inclusion holds. At each fixed $z$ in this square, no single probability distribution on permutations works simultaneously for all $x,y\in\mathbb C$.** These are two different assertions: the pointwise convex weights can depend on $x$ and $y$.

This is a quantitative extension of the coordinator's [exact two-free-eigenvalue certificate](https://github.com/g4nesh/math-commons/blob/cb253d9c3c7746fd2278db9d41a5eb1668ed52b4/research/marcus-de-oliveira/two-free-rank-one/README.md). That result proves the center case $z=0$ and qualitative stability under several simultaneous perturbations. Here only the fourth eigenvalue moves, but the neighborhood radius and both strict margins are explicit. The exact fixed data are

$$U=\frac1{625}\begin{pmatrix}-31&492&-240&300\\-300&225&500&0\\240&-180&225&500\\-492&-256&-180&225\end{pmatrix}.$$

Indices in this report are one-based: $a_1=x,a_2=y,a_3=2-i/4,a_4=z$. No genericity assumption on $x,y,z$, nonzero determinant assumption, or nondegenerate pairing hull assumption is needed. The square is a sufficient region, not a sharp boundary.

[Exact bound](bound.json) · [Replay](radius.py) · [Replay output](verification.json) · [Tests](test_radius.py)

## The identity that proves pointwise inclusion

For a column-to-row permutation $p$, let

$$D_p(x,y;z)=\prod_{j=1}^4(a_{p(j)}+b_j)-\det(A(x,y,z)+B).$$

At fixed $z$ this has four complex coefficients, in the order $(1,x,y,xy)$. Write

$$\ell_p(L;z)=\Re\sum_{k=0}^3L_k(D_p)_k(z),\qquad \Delta(L)=L_0L_3-L_1L_2.$$

The base certificate provides 35 permutation pairs and positive coefficients $w^0_{pq}$ summing to one, together with real $h_0,k_0$, such that

$$\sum_{p,q}w^0_{pq}\ell_p(L;0)\ell_q(L;0)=h_0\Re\Delta(L)+k_0\Im\Delta(L).$$

Our replay uses precisely those 35 pairs. It proves that for every $z$ in the square the same support admits positive coefficients $w_{pq}(z)$ summing to one and real $h(z),k(z)$ satisfying the identity with $z$ in place of zero. These are coefficients of a quadratic identity, not convex masses directly attached to the 24 pairing vertices.

If the determinant lay outside its pairing hull at some $x,y$, strict separation from that finite planar convex hull would give $\lambda\in\mathbb C$ with all $\Re(\lambda D_p(x,y;z))>0$. Taking $L=\lambda(1,x,y,xy)$ makes all factors on the left positive, but $\Delta(L)=0$. The identity is a contradiction. This proves inclusion for all complex $x,y$, including coincident eigenvalues and degenerate hulls.

## Exact matrix perturbation bound

Write $z=u+iv$. Each coefficient row of $D_p$ is affine in $(u,v)$, because the determinant and pairing products are affine in the fourth diagonal entry. Therefore the 36 coefficients of each quadratic product are polynomials of degree at most two in $(u,v)$.

Let $C(u,v)$ be the $37\times37$ real matrix with the 35 quadratic-product columns, the columns for $-\Re\Delta,-\Im\Delta$, and the normalization row. Then

$$C(u,v)=C_0+uC_u+vC_v+u^2C_{uu}+uvC_{uv}+v^2C_{vv}.$$

The replay constructs all six matrices exactly from determinants and products. It verifies the base solution $b_0=(w^0,h_0,k_0)$, and computes the exact rational inverse of $C_0$, checking both residuals. Let $S$ be diagonal with the first 35 entries $w^0_{pq}$ and final entries $\max(1,|h_0|),\max(1,|k_0|)$. Set

$$E(u,v)=S^{-1}C_0^{-1}(C(u,v)-C_0)S,\qquad y_0=S^{-1}b_0.$$

The first 35 entries of $y_0$ equal one. For each of the five nonconstant matrices the replay records both the exact infinity norm of $S^{-1}C_0^{-1}C_jS$ and the exact infinity norm of its action on $y_0$. Triangle inequalities on the closed square give bounds $\rho,\nu$ satisfying

$$\|E(u,v)\|_\infty\le\rho<1/5,\qquad \|E(u,v)y_0\|_\infty\le\nu<3/100.$$

All rational values, including the chosen radius, are independently recomputed in [the bound file](bound.json). In particular $\rho<1$, so the Neumann series makes $I+E$ invertible throughout the whole square. The unique solution is $b(u,v)=S(I+E)^{-1}y_0$, and

$$\|S^{-1}b(u,v)-y_0\|_\infty\le\frac{\nu}{1-\rho}<3/80<1/2.$$

Consequently every product coefficient is at least half its positive base value. The uniform lower bound recorded in the certificate is greater than $1/4000$. The normalization row enforces sum one, and the other 36 rows enforce the quadratic identity. These analytic inequalities prove the result at irrational coordinates as well as rational ones. The square's corners are regression checks, not a sampling argument for its interior.

## A strict obstruction to common weights throughout the square

The base certificate also gives a fixed eight-dimensional real coefficient separator $v_*$. At $z=0$, its minimum value on the 24 gap rows is

$$m=267131743/156385125000>0.$$

Write a gap row as $r_p(u,v)=r_p(0,0)+u r_{p,u}+v r_{p,v}$. The replay computes, separately for each permutation,

$$m_p=v_*\cdot r_p(0,0),\qquad d_p=|v_*\cdot r_{p,u}|+|v_*\cdot r_{p,v}|.$$

At radius $r=1/32768$, exact arithmetic gives

$$\min_p(m_p-rd_p)=\frac{27177805772189}{25622138880000000}>\frac m2>0.$$

Thus every one of the 24 gap coefficient vectors remains strictly on the same side of this functional. Their convex hull cannot contain zero, which rules out any common distribution valid for every $x,y$. This is not a pointwise determinant counterexample; the preceding identity proves pointwise inclusion everywhere in the family.

## Reproduction and limits

From the repository root, using Python 3.9 or later and only its standard library:

```sh
python3 -B research/marcus-de-oliveira/rankone-quantitative/radius.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/rankone-quantitative -p 'test_*.py'
```

To regenerate the exact recorded bound and report:

```sh
python3 -B research/marcus-de-oliveira/rankone-quantitative/radius.py --derive --report research/marcus-de-oliveira/rankone-quantitative/verification.json
```

`coefficients_at(u, v)` returns the 35 positive identity coefficients followed by $h,k$, solving over rational numbers and checking the residual. Its parameters must be integers or `Fraction` values in the certified square. The file-level `verify` recomputes every bound rather than trusting the JSON. Imports use the repository's archived arithmetic and base certificate; no optimization package or floating-point solver is used during replay. The numerical displays above are rational inequalities, not tolerance tests.

The tests replay the complete bound, reconstruct an off-center quadratic identity, check its rank-one specialization against direct physical determinants, verify the separator and closed endpoints, and reject malformed or out-of-square parameters and modified bound data. This is an exact, family-specific certificate and same-owner reproducibility work. Independent-owner review is pending. It does not establish the general MDO conjecture, nor give an explicit formula for the pointwise vertex masses.
