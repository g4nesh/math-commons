# Three arbitrary complex eigenvalues while the fourth moves

Let $U_0$ be the archived four-cycle relative eigenbasis and let $b^0=(3-2i,2,1+2i,-i)$. There are a number $r>0$, a neighborhood of $U_0$ in $U(4)$, and a neighborhood of $b^0$ in $\mathbb C^4$ with the following property. For every $U,b$ in these neighborhoods, every $|\delta|<r$, and every $z\in\mathbb C$, there exist **24 strictly positive real probabilities**, independent of $x,y$, satisfying
$$
\det\!\left(\operatorname{diag}(x,z,y,\delta)+U D_b U^*\right)
=\sum_p m_p\prod_{t=1}^4(a_{p(t)}+b_t)
\quad\text{for all }x,y\in\mathbb C,
\qquad a=(x,z,y,\delta).
$$
The probabilities can depend on $U,b,z,\delta$. The proof gives existence of a common positive radius for every middle eigenvalue. The explicit [four-collision bounds](../README.md) give radius $1/1024$ for both spectral motions there. The full theorem leaves the numerical value of its global radius open.

## One fixed set of nine real equations

Translate the first and third coordinates by writing $h=x-\delta$, $k=y-\delta$, $w=z-\delta$, $b'=b+\delta\mathbf1$. The translated diagonal is $(h,w,k,0)$. When a pairing assigns columns $(i,j,l,t)$ to successive rows, its product is
$$(h+b'_i)(w+b'_j)(k+b'_l)b'_t.$$
The target has the same constant value as each pairing at $h=w=k=0$. Subtract that scalar constant. In the value at $h=k=0$ cancel the factor $w$ and retain the resulting equation at $w=0$ as well. With mass one enforced, the four complex entries of the resulting raw pairing column are
$$
 b'_i b'_l b'_t,\quad
 (w+b'_j)b'_l b'_t,\quad
 (w+b'_j)b'_i b'_t,\quad
 (w+b'_j)b'_t.
$$
The first entry gives the coefficient of $w$ in the constant term. The other three are the coefficients of $h,k,hk$. Form a real system from the mass equation and the eight real components, using the corresponding exact determinant coefficients as the right side. Any positive solution certifies every $h,k$ simultaneously.

On the reciprocal chart $v=1/w$, keep the first complex equation and multiply the last three by $v$, including the same factors on the right. Replace each $(w+b'_j)$ there by $(1+v b'_j)$. The equations extend continuously to $v=0$. The determinant coefficients depend continuously on $U,b,\delta$ in both charts. The point at infinity is used only to control arbitrarily large finite $w$.

## Strict probabilities at the base parameters

The [eight angular identities](../../complex-direction-cover/README.md) cover all $w$. Their original masses can vanish for individual permutations. In each patch write $u$ for its real radial variable and $R$ for the largest radial degree, two or four. At every Bernstein control index add an empty radial product and $u^R$ to the numerator for every permutation, both with coefficient $10^{-7}$. Subtract $10^{-7}$ times the corresponding stored [correction](corrections.json) from each original positive coefficient. Exact replay verifies that each adjusted coefficient stays positive, all normalization equations still hold and the complete mixed polynomial gap still vanishes.

The Bernstein functions are nonnegative and sum to one. Every numerator is therefore at least $10^{-7}(1+u^R)$. At finite $u$ this gives strict positivity. After division by $u^R$ each numerator has a positive finite leading coefficient and the normalized law has a strictly positive limit at radial infinity. The same limit satisfies the reciprocal system by continuity. At $w=0$ the equation obtained by canceling $w$ holds by taking a radial limit, or by the polynomial identity. This constructs a strictly positive base law at every point of the compactified complex plane.

## The four exceptional points

Set $w=X+iY$ at the base parameters. Each nine-column determinant of the real system is a polynomial in $X,Y$ of total degree at most six, since its mass row and two middle-singleton rows are constant. The six minors recorded in [rank.json](rank.json) admit explicit quadratic polynomial multipliers producing the two identities
$$
 6X+5Y^4-20Y^2+3Y+12=\sum_m P_m(X,Y)D_m(X,Y),\qquad
 Y^2(Y-2)(Y-1)(Y+2)=\sum_m Q_m(X,Y)D_m(X,Y).
$$
A rank drop forces both left sides to vanish. Their only joint zeros are $-b^0_j$, and the [collision replay](../check.py) gives rank eight at all four. At infinity a nine-column determinant of the reciprocal matrix equals $100000$. All other parameter values have full row rank. The rank data use row-to-column words, the inverse convention to the archived angular weights. The rank proof uses all 24 columns, so this reordering preserves its conclusion.

The checker evaluates the two identities by integer determinants and exact rational multipliers on $X,Y\ge0$, $X+Y\le8$, with integer coordinates. These 45 values determine a polynomial of total degree at most eight: for example the Newton products $\binom Xi\binom Yj$ for $i+j\le8$ form a triangular basis. Polynomial identity gives the claimed rank conclusion for the whole complex plane.

## Obtaining one radius for every middle eigenvalue

At a regular point of the compactified plane choose any nonsingular nine-column block. Fix the other 15 masses at their positive base values. The nine adjustable masses are the inverse of that block applied to the right side minus the contribution of the fixed masses. This formula is continuous in the chart coordinate, $\delta,U,b$, and gives the original strictly positive law at the center. Thus it works on a product of open neighborhoods in these variables. This includes the reciprocal-chart point at infinity. At a collision use the [ten-equation fourth-row marginal construction](../README.md), whose positive center and invertible block give the same local conclusion. Translation between $w$ and $z=w+\delta$ is continuous there.

These open parameter patches cover the compact Riemann sphere at the base values of $\delta,U,b$. Choose finitely many. Intersect their neighborhoods of those three parameters; the intersection contains a disk $|\delta|<r$ and neighborhoods of the archived $U,b$. Every $w$ lies in at least one retained patch. This proves the theorem with the same $r$ for the entire sphere, hence for every actual $z$. The proof requires neither a common law across different $z$ nor a single global choice of a nine-column basis.

## Calculate a shared weight vector

For Gaussian-rational inputs at the archived $U,b$, the calculator returns all 24 exact masses. Pass the real and imaginary parts of the actual second and fourth eigenvalues, respectively.

~~~sh
python3 -B research/marcus-de-oliveira/projective-fiber-stability/uniform/weights.py 3/5 7/10 1/2048 -1/4096
~~~

The output keys map columns to rows using one-based indices. It checks the original coefficient equations, so the returned law works for every first and third eigenvalue. The exact simplex calculation first enforces feasibility and then maximizes the smallest of the 24 probabilities. With no arguments it checks four covers against dense matrix determinants, including a collision and a large middle eigenvalue.

## Reproduction and scope

Run from the repository root:

~~~sh
python3 -B research/marcus-de-oliveira/projective-fiber-stability/check.py
python3 -B research/marcus-de-oliveira/projective-fiber-stability/uniform/verify.py
~~~

Both commands use Python's standard library. The first checks the exceptional centers and quantitative margins. The second checks all angular corrections, both polynomial identities and the infinity determinant. The compactness deduction is written above. The claims concern neighborhoods of one specific relative unitary and opposing spectrum in dimension four. Independent-owner mathematical review is open.
