# A spectral triangle that works for every relative unitary

Take $a_0,a_1,a_2$ and $b_0,b_1,b_2$ complex, each eigenvalue repeated twice in its respective six-dimensional normal matrix. For $\pi\in S_3$ set $p_\pi=\prod_{r=0}^2(a_r+b_{\pi(r)})$. These six complex numbers are simple to compute from the spectra alone.

**Theorem.** If $\operatorname{conv}\{p_\pi\}$ has at most three vertices, the Marcus–de Oliveira inclusion holds for **every** relative unitary $U\in U(6)$. More precisely, let $v_1,\ldots,v_h$ be the actual distinct endpoints or vertices, $1\le h\le3$. The determinant is a convex combination of the at most $h(h+1)/2$ points $v_i v_j$ with $i\le j$. These are genuine six-coordinate pairing products. An open family of noncollinear pairs of spectra satisfying the condition is given below. The weights may depend on both spectra and $U$.

## Two elementary positive cones

Index the permutations by $012,120,201,021,102,210$, with sign vector $s=(1,1,1,-1,-1,-1)$. Since $(a_r+b_j)_{rj}$ has rank at most two, $s^\mathsf T p=0$. The constant vector $\mathbf1$ also satisfies this balance.

A real vector $f\ge0$ with $s^\mathsf T f=0$ is a nonnegative combination of the nine indicators $g_{eo}$ of an even index $e$ and an odd index $o$. Indeed the even and odd coordinates have the same total. Transport the former to the latter with these row and column sums to obtain $f=\sum c_{eo}g_{eo}$ with $c_{eo}\ge0$.

Write the determinant as $\sum_{i\le j}q_{ij}p_i p_j$, using the real signed weights obtained by grouping the Leibniz expansion. Define a symmetric matrix $Q$ by $Q_{ii}=q_{ii}$ and $Q_{ij}=Q_{ji}=q_{ij}/2$ for $i<j$. The following finite positivity statement holds for every physical $U$:
$$g_{eo}^{\mathsf T}Q g_{e'o'}\ \ge\ 0.$$
To check it, the 21 coefficients of twice this expression as a linear form on $q$ are $g_{eo}(i)g_{e'o'}(j)+g_{eo}(j)g_{e'o'}(i)$ for $i\le j$, with the same convention on the diagonal. Divide by their positive integer content. Each of the 45 unordered indicator pairs then equals a primitive coefficient row of the 21 polynomial products in the variables $a,b$. This is an exhaustive finite integer identity, rechecked from polynomial expansion in the verifier. The corresponding coefficient in the physical determinant is a sum of squared unitary minors by Cauchy–Binet. Consequently
$$f^\mathsf T Qg\ge0\quad\text{for all } f,g\ge0\text{ with }s^\mathsf T f=s^\mathsf T g=0.$$
The bilinear form is used on real vectors; its complex-bilinear extension gives the determinant as $p^\mathsf TQp$. [Two-sided grouped coefficient proof](../double-paired-contingencies/README.md)

## Positive barycentric weights

In a triangle, its three barycentric coordinates are real affine functions of the point. On a segment use its two endpoint coordinates, and at a point use the constant one. Let $f_k(i)$ be the coordinate of $p_i$ at $v_k$. Then $f_k\ge0$, $\sum_k f_k=\mathbf1$ and $p=\sum_k v_k f_k$. Being affine functions of $\Re p,\Im p$, these vectors are balanced: the even and odd triples have identical vector sums and equal cardinalities. Define $\gamma_{kl}=f_k^\mathsf T Q f_l$. The preceding positivity gives $\gamma_{kl}\ge0$ and
$$\sum_{k,l}\gamma_{kl}=\mathbf1^\mathsf T Q\mathbf1=\sum_{i\le j}q_{ij}=1,\qquad
\det(D_a+UD_bU^*)=\sum_{k,l}\gamma_{kl}v_k v_l.$$
The displayed construction gives weights at the original product vertices and works for singular matrices, endpoint coincidences and degenerate hulls.

This common barycentric construction has a precise geometric limit. Suppose one tries to use balanced nonnegative coordinates on the actual extreme points of an arbitrary hull, with a probability vector in each row and the correct barycenter $p_i$. Every extreme point forces its own row to be concentrated at itself. The equal even and odd column sums must therefore each assign at least unit mass to every distinct extreme point. Each side has total mass three, so at most three extreme points are possible. This limits that construction only. Spectra with four or more extreme points need a different certificate.

## An open set of entirely complex examples

At the following center all three values on each side are distinct and noncollinear:
$$a=(-2-2i,-1+3i,3+3i),\qquad b=(-1-3i,-3-i,-3-2i).$$
The six products in parity order are
$$(-14+22i,\ 46+14i,\ -16+20i,\ -34+34i,\ -6+10i,\ 56+12i).$$
The last three form a triangle. The first three have barycentric matrix
$$\frac1{193}\begin{pmatrix}95&80&18\\18&5&170\\80&108&5\end{pmatrix}.$$
Each entry is strictly positive. For any complex perturbations of the six spectral values by at most $1/16384$ in modulus, every relative unitary still satisfies MDO.

Here are conservative bounds for that radius. The original factors $|a_r+b_j|$ are less than ten. After perturbations of size $\delta\le1$ the factors have modulus at most twelve. Telescoping a three-factor product bounds its movement by $6\cdot12^2\delta\le 1/16$ at the stated radius. At the center every pair of products has distance less than 100. In a signed doubled triangle area $\det(y-x,z-x)$ each difference moves by at most $2/16$. Expanding the determinant bounds the area change by $2(100)(2/16)+(2/16)^2<40$. The oriented doubled area of the outer triangle is 1544; all nine inner barycentric numerators are at least 40. Their signs survive uniformly. The proof therefore applies throughout a spectral neighborhood independent of $U$.

Run `python3 -I -B research/marcus-de-oliveira/paired-triangular-spectra/check.py`. The replay checks the 45 universal coefficient identities, all barycentric numbers and perturbation bounds, and directly reconstructs a nonreal determinant on a dense complex unitary from the six claimed hull products. [Exact checker](check.py)
