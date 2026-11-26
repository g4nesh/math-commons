# Common radial weights over complex directions

Fix the archived four-cycle $U$, $b=(3-2i,2,1+2i,-i)$, and $A=\operatorname{diag}(x,z,y,0)$. The first and third eigenvalues $x,y$ can be arbitrary complex numbers. For the second eigenvalue $z$, this report proves the inclusion for every radius along the following real lines:

$$z=u(1+ik)\quad\text{or}\quad z=u(k+i),
\qquad u\in\mathbb R,\quad -1\le k\le\tfrac12.$$

The covered directions include the real and imaginary axes and both orientations of every indicated line. In ordinary slope coordinates these six patches cover every line outside the open positive-slope range $(1/2,2)$. The statement holds also at $z=0$. Thus three eigenvalues here can vary over a region with nonempty interior in $\mathbb C^3$. On the full three-variable family, the earlier exact separator at fourth eigenvalue zero rules out one common law. The weights below can change with $z$ and work for every $x,y$. [Three-variable obstruction](../moving-fourth-spectrum/README.md) · [Real-line precursor](../two-free-real-atlas/README.md)

## How to read a patch

Let $D_p(x,z,y)$ be the pairing vertex minus the actual determinant. Its coefficients come from the exact polynomial expansion of the unitary fixture. On a patch $[l,h]$ write $k=(1-\lambda)l+\lambda h$, with $0\le\lambda\le1$. Define $\omega(k)=1+ik$ for the horizontal chart and $\omega(k)=k+i$ for the vertical chart. We construct nonnegative polynomial numerators
$$N_p(\lambda,u)=(1-\lambda)\sum_r c_{0,p,r}(u-r)^2+
                      \lambda\sum_r c_{1,p,r}(u-r)^2.$$
Every stored $c$ is a positive rational. Within each endpoint list the centers do not have a common root, so $H=\sum_p N_p$ is positive for every finite $u$, including either endpoint in $\lambda$. The weights are $N_p/H$.

The gap $D_p(x,\omega(k)u,y)$ is itself affine in $\lambda$. Consequently the identity $\sum_p N_p D_p=0$ can be checked by expanding a degree-two polynomial in the Bernstein basis in $\lambda$. For the two endpoint coefficient lists $C_0,C_1$ and their gap operators $M_l,M_h$, its three coefficients vanish exactly when
$$M_l C_0=0,\qquad M_l C_1+M_h C_0=0,\qquad M_h C_1=0.$$
Each equation here is a polynomial identity in $x,y,u$. The checker expands every one over Gaussian rationals, including the mixed angular coefficient. Agreement at just the two ends would miss the middle equation.

The six half-unit patches each use 64 primitive coefficients indexed by endpoint, pairing and real center. A numerical linear program found small supports. On each support the exact linear equations were solved over the rationals. The archived replay requires only Python's standard library. It verifies the positivity, normalization, common-root condition, coefficient identities and interval coverage. [Positive coefficients](patches.json) · [Exact replay](check.py)

~~~sh
python3 -B research/marcus-de-oliveira/complex-direction-cover/check.py
~~~
