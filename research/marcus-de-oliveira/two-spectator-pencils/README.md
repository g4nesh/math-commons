# One physical two-spectator pencil, with four fixed weights

**Written polynomial proof and exact coefficient check; independent-owner review pending.** This family allows two arbitrary complex outside eigenvalues simultaneously. It is an explicit certificate for a particular pair of overlapping normal blocks, not a theorem that arbitrary two-variable pointwise inclusions admit common weights.

Let

$$P=\begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix},\qquad
G=(I+P)^{-1}=\frac{I-P+P^2}{2}
=\frac12\begin{pmatrix}1&-1&1\\1&1&-1\\-1&1&1\end{pmatrix}.$$

As a polynomial in the unitary $P$, $G$ is normal. Its eigenvalues are

$$a=\frac12,\qquad b=\frac{1+i\sqrt3}{2},\qquad c=\frac{1-i\sqrt3}{2};
\qquad b+c=1,\quad bc=1.$$

On four coordinates put $A(h)=G\oplus[h]$ and $B(k)=[k]\oplus G$. Both are normal for **every** $(h,k)\in\mathbb C^2$. Their sum is

$$\begin{pmatrix}
k+1/2&-1/2&1/2&0\\
1/2&1&-1&1/2\\
-1/2&1&1&-1/2\\
0&-1/2&1/2&h+1/2
\end{pmatrix},$$

whose determinant is

$$D(h,k)=1+\frac32(h+k)+2hk.\tag{1}$$

Order the spectra as $\alpha=(h,a,b,c)$ and $\beta=(k,a,b,c)$, and define $P_\sigma=\prod_{i=0}^3(\alpha_i+\beta_{\sigma(i)})$. The following fixed weights certify the original MDO inclusion:

| Zero-based permutation $\sigma$ | Weight | Its vertex polynomial |
|---|---:|---|
| $(0,2,3,1)$ | $2/7$ | $\frac74(h+k)$ |
| $(1,0,2,3)$ | $1/3$ | $4(h+1/2)(k+1/2)$ |
| $(2,3,1,0)$ | $4/21$ | $\frac74(h+b)(k+c)$ |
| $(3,2,0,1)$ | $4/21$ | $\frac74(h+c)(k+b)$ |

The masses are positive and sum to one. Their weighted sum is

$$\frac12(h+k)+\frac43(h+1/2)(k+1/2)
 +\frac13\bigl[(h+b)(k+c)+(h+c)(k+b)\bigr].$$

Using $b+c=bc=1$ gives exactly (1). This proves the entire complex two-parameter family, including repeated eigenvalues and singular matrix sums. It does not use division by $h$, $k$, a determinant, or a hull vertex.

## A sufficient interface proxy can fail inside this family

Set $k=0$, retain the middle two coordinates as the interface, and eliminate the two outside coordinates. For $q=h+1/2\ne0$ the left and right effective matrices are

$$M=\begin{pmatrix}1&-1\\0&1\end{pmatrix},\quad
K=\begin{pmatrix}d&-d\\e&d\end{pmatrix},\qquad
d=\frac12+\frac1{4q},\quad e=\frac12-\frac1{4q}.$$

Their exterior determinant factor is $q/2$, and direct calculation again gives $\det(A(h)+B(0))=1+3h/2$.

The centered numerical range of $N=M-I$ is a disk of radius $1/2$. For $T=K-dI$, its numerical radius is $(|d|+|e|)/2$, while $\operatorname{tr}(NT)=-e$. The centered-ellipse-product proxy from [the interface analysis](../pointwise-gluing/larger-interface.md) therefore passes exactly when $|e|\le|d|$. Here

$$|d|^2-|e|^2=\frac{\operatorname{Re}q}{2|q|^2},$$

so that proxy covers $\operatorname{Re}h\ge-1/2$ when $q\ne0$. The singular endpoint follows by continuity for the original determinant inclusion. For $\operatorname{Re}h>-1/2$ and $h\ne0$, both effective matrices are nonnormal and have no common eigenvector: $M$ has unique eigenline $\mathbb C e_1$, whereas $e\ne0$ prevents its invariance under $K$.

At $h=-1,k=0$, the proxy fails: $d=0$, $e=1$, so the desired centered value has modulus $1/2$ but the product-range disk has radius $1/4$. Nevertheless (1) equals $-1/2$ and the same four weights above certify it. Thus the proxy is sufficient, not necessary for the original inclusion. This is an explicit physical example, distinct from an abstract counterexample to multivariable coefficient separation.

## Reproduction and limits

```sh
python3 research/marcus-de-oliveira/two-spectator-pencils/check.py
```

The standalone standard-library checker uses exact arithmetic in $\mathbb Q(i\sqrt3)$ and bivariate polynomials. It verifies normality and the three eigenvalues of $G$, builds the original sum from its two physical blocks, expands its determinant by all 24 Leibniz terms, independently constructs the four spectral vertices, and matches all polynomial coefficients and probability masses. It also checks the strict proxy failure at $h=-1$. [Recorded output](results.json)

No finite sampling is being promoted to a universal claim: equality of the full polynomial coefficients supplies the universal quantifier over $h,k$. The checker and derivation share one GitHub owner; independent review remains pending. These calculations do not prove a general two-coordinate gluing theorem or decide whether arbitrary local normal blocks have a common two-spectator certificate. That is a separate investigation. No novelty claim is made for this particular family.
