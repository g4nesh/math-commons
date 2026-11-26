# A determinant family with no common permutation weights

Fix the real orthogonal four-cycle $U$ and the complex spectrum $b=(3-2i,2,1+2i,-i)$ from the [archived example](../fourth-boundary/README.md). Set
$$A(x,y)=\operatorname{diag}(x,y,2-i/4,0),\qquad B=U\operatorname{diag}(b)U^*.$$

**For every $(x,y)\in\mathbb C^2$, $\det(A(x,y)+B)$ belongs to its MDO pairing hull. No probability distribution on the 24 permutations works for all $(x,y)$. Both conclusions persist in an open neighborhood of the fixed data.** The argument below uses exact rational arithmetic. It gives a sufficient pointwise proof method when a direct coefficient-matching approach has a strict obstruction. [Rational evidence](certificate.json) · [Executable verifier](check.py)

## Why two different conclusions can hold

Write $T(x,y)$ for the determinant, $V_p(x,y)=\prod_j(a_{p(j)}+b_j)$ for its 24 pairing vertices, and $D_p=V_p-T$. Each $D_p$ has four complex coefficients in the order $(1,x,y,xy)$. The checker reconstructs them by evaluating determinants and vertices at $(0,0),(1,0),(0,1),(1,1)$. Multi-affinity in the two diagonal entries makes these four evaluations exact.

For $L=(L_0,L_1,L_2,L_3)\in\mathbb C^4$, define the real linear forms
$$\ell_p(L)=\Re\sum_{k=0}^3 L_k(D_p)_k.$$
A common law would put the zero coefficient vector in the real convex hull of the 24 coefficient vectors. The explicit $L$ stored as eight real coordinates in the certificate satisfies
$$\min_p\ell_p(L)=267131743/156385125000>0.$$
Taking a convex combination immediately rules out that law.

A violation of the **original** inclusion at any particular $(x,y)$ would give a complex number $\lambda$ with $\Re(\lambda D_p(x,y))>0$ for all $p$, by strict separation of a point from a compact planar convex hull. Such a separator has the special form
$$L=(\lambda,\lambda x,\lambda y,\lambda xy),\qquad L_0L_3-L_1L_2=0.$$
We can therefore prove pointwise inclusion by excluding only strict separators satisfying this last equation. This is the rank-one condition on a complex two-by-two matrix. Degenerate pairing hulls, zero determinants and coincident eigenvalues need no separate nondegeneracy assumption in this separation argument.

## The exact quadratic identity

The certificate lists 35 positive rational numbers $w_{pq}$ summing to one and two real rational numbers $h,k$. With $\Delta=L_0L_3-L_1L_2$, they satisfy
$$\sum_{p\le q} w_{pq}\,\ell_p(L)\ell_q(L)=h\,\Re\Delta+k\,\Im\Delta.$$
Only the listed pairs occur in the sum. The checker expands both sides as quadratic forms in the eight real coordinates of **arbitrary** $L$. All 36 coefficients agree exactly. At a hypothetical pointwise violation, every factor on the left is strictly positive and $\Delta=0$ on the right. The contradiction proves MDO for the entire family. The 35 numbers certify an identity between separators. The convex weights on permutation vertices at a given $(x,y)$ follow by separation and may depend on both eigenvalues.

To find the certificate, we considered the 300 unordered pairs of permutations including repetitions. Each pair supplies one known quadratic form. We asked a linear program for nonnegative coefficients whose weighted sum lies in the real span of $\Re\Delta,\Im\Delta$, with total weight one. We retained its support, solved the resulting 37 equations over $\mathbb Q$, and checked positivity and both sides of the identity independently. No floating-point equality or sample in $(x,y)$ enters the proof.

## Stability and limits

For these 35 columns, append the columns for $-\Re\Delta$ and $-\Im\Delta$, and a final normalization row equal to one on the product columns and zero on the two determinant columns. The resulting $37\times37$ rational matrix has determinant $396862542\bmod 1000000007$. Thus it is nonsingular. Its entries are continuous in the two fixed diagonal entries, all four eigenvalues of $B$, and the relative unitary $U$. The unique solution for the same support remains strictly positive on its first 35 coordinates in a neighborhood. The original strict separating margin is also continuous and stays positive there with the same separator. Both assertions therefore persist under simultaneous small perturbations of all these data, with $x$ and $y$ still completely free. This argument proves existence of a neighborhood; it does not specify its radius.

For more free complex variables the evaluation vector consists of products of the selected variables. Its rank-one tensor still satisfies quadratic relations. Searching for identities in their span is a separate question. Failure to find a quadratic certificate would leave the original membership question open.

Run the replay using the Python standard library:

~~~sh
python3 -B research/marcus-de-oliveira/two-free-rank-one/check.py
~~~
