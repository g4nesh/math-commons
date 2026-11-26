# Three arbitrary complex eigenvalues on the archived four-cycle

Fix the archived four-cycle $U$, $b=(3-2i,2,1+2i,-i)$, and $B=U\operatorname{diag}(b)U^*$. For **every** $x,z,y\in\mathbb C$, put $A=\operatorname{diag}(x,z,y,0)$. We prove that $\det(A+B)$ lies in its permutation pairing hull. For every $z$ the certificate constructs one probability vector that works for both other complex variables simultaneously. An earlier exact separator proves that no constant vector works for the full three-variable family. This result concerns the indicated relative eigenbasis, fixed spectrum $b$ and fixed fourth eigenvalue. [Original four-cycle and obstruction](../moving-fourth-spectrum/README.md) · [Real-line precursor](../two-free-real-atlas/README.md)

Write the second diagonal entry in one of two charts:

$$z=u(1+ik)\quad\text{or}\quad z=u(k+i),\qquad u\in\mathbb R,\quad -1\le k\le1.$$

For $z=a+ib\ne0$, use the first chart with $u=a,k=b/a$ if $|a|\ge|b|$, and the second with $u=b,k=a/b$ otherwise. For $z=0$ take $u=k=0$. Both charts have four closed slope intervals with endpoints $-1,-1/2,0,1/2,1$. There are no gaps between them.

## How to read a patch

Let $D_p(x,z,y)$ be the pairing vertex minus the actual determinant. Its coefficients come from the exact polynomial expansion of the unitary fixture. On a patch $[l,h]$ write $k=(1-\lambda)l+\lambda h$, with $0\le\lambda\le1$. Define $\omega(k)=1+ik$ for the horizontal chart and $\omega(k)=k+i$ for the vertical chart. Let $B_{j,d}(\lambda)=\binom dj(1-\lambda)^{d-j}\lambda^j$ and $H_R(u)=\prod_{r\in R}(u-r)^2$, including the empty product. We construct
$$N_p(\lambda,u)=\sum_{j,R} c_{j,p,R}\,B_{j,d}(\lambda)\,H_R(u).$$
Every stored $c$ is a positive rational. For each control index $j$, all its radial products have no common real root: the intersection of their sets of centers is empty. At least one Bernstein factor is positive for every $\lambda$, so $H=\sum_p N_p$ is positive for every finite $u$, including both endpoints in $\lambda$. The weights are $N_p/H$. The executable function <code>weights_at(real, imag)</code> returns them as exact fractions for any Gaussian-rational $z$. Permutation words map columns to one-based diagonal rows.

The gap $D_p(x,\omega(k)u,y)$ is affine in $\lambda$. Thus $\sum_p N_pD_p$ has angular degree at most $d+1$. Denote the gap operators at the two ends by $M_l,M_h$, and the coefficient blocks of the numerators by $C_0,\ldots,C_d$. Multiplication in the Bernstein basis gives the cleared coefficients
$$ (d+1-j)M_l C_j+jM_h C_{j-1}=0,\qquad 0\le j\le d+1,$$
omitting out-of-range blocks. Each equality is a polynomial identity in $x,y,u$. The checker expands every coefficient over Gaussian rationals, including the mixed terms. At degree one this includes the middle equation $M_lC_1+M_hC_0=0$.

| Slopes in each chart | Angular degree | Maximum radial degree | Positive coefficients per patch |
|---|---:|---:|---:|
| $[-1,-1/2],[-1/2,0],[0,1/2]$ | 1 | 2 | 64 |
| $[1/2,1]$ | 4 | 4 | 247 |

The last two patches use products with centers in $\{-2,-1,0,1,2\}$. The construction grew out of a quartic identity along $z=(1+i)u$ supplied by the quantitative-weight agent. Numerical linear programming selected supports. We solved each resulting system over the rationals and retained only supports with positive exact solutions. The public replay checks all **878** positive coefficients directly using Python's standard library. It also checks normalization, absence of a common radial root, all polynomial identities and full interval coverage. [Positive coefficients](patches.json) · [Exact replay](check.py) · [Independent determinant examples](test_cover.py)

**Second implementation.** A separate checker reconstructs $B$ from the archived entries and uses Gaussian elimination on the original dense matrices, with its own arithmetic and weight decoder. A patch of angular degree $d$ and radial degree $R$ has weighted gap of degree at most $(d+1,R+1,1,1)$ in $(\lambda,u,x,y)$. The grid contains $d+2$, $R+2$, two and two distinct values, respectively. Repeated univariate interpolation therefore certifies the entire identity, not just its sampled points. The replay evaluates **576 dense determinants and 13,824 pairing vertices**. The second implementation was produced under the same GitHub owner. Independent-owner mathematical review remains open. [Grid replay](grid_replay.py)

## Changing spectral coordinates

The same weights apply in any common Möbius coordinate on the signed spectra. This standard symmetry gives versions where the fourth eigenvalue is nonzero and the radial lines become generalized circles. Write $C=-B$ and $f(t)=(pt+q)/(rt+s)$ with $ps-qr\ne0$. Assume the spectra of $A,C$ avoid the pole. Without requiring $AC=CA$,
$$
(rA+sI)[f(A)-f(C)](rC+sI)
=(pA+qI)(rC+sI)-(rA+sI)(pC+qI)
=(ps-qr)(A-C).
$$
After taking determinants, $T$ and **every** pairing vertex acquire the same nonzero complex multiplier $(ps-qr)^4/[\det(rA+sI)\det(rC+sI)]$. Every pointwise weight identity therefore transports to the pair with spectra $f(a_i),-f(-b_j)$ and the same relative eigenbasis. For fixed $z$, the weights remain independent of the first and third complex eigenvalues. This transports the entire configuration together. Applying the map to the four diagonal entries alone would change the problem.

~~~sh
python3 -B research/marcus-de-oliveira/complex-direction-cover/check.py
python3 -B research/marcus-de-oliveira/complex-direction-cover/grid_replay.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/complex-direction-cover -p 'test_*.py' -v
~~~
