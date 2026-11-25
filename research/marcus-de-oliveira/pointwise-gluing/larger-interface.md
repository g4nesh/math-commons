# A normal effective matrix allows a larger interface

Use coordinates $L,J,R$, with $J$ of arbitrary positive size. Let $X$ act on $L,J$ and $Y$ on $J,R$. In the conjugated global sum set $B'=Y D_{b_{J,R}}Y^*$, $T=D_{a_R}+B'_{RR}$. Assume $T$ is invertible and that its Schur complement at the interface,
$$K=B'_{JJ}-B'_{JR}T^{-1}B'_{RJ},$$
is normal. Write $K=W D_\kappa W^*$. Suppose the ordinary MDO inclusion holds for the required local spectra at the rotated relative unitaries $X(I_L\oplus W)$ and $(W^*\oplus I_R)Y$. In particular it suffices that both local dimensions are known MDO dimensions. Then the global MDO inclusion holds for these data.

**Proof.** Eliminate $R$. The remaining determinant is that of $X^* D_{a_{L,J}}X+\operatorname{diag}(b_L,K)$, multiplied by $\det T$. The first local inclusion therefore expresses the full determinant as a convex mixture of terms
$$C\det(T)\prod_{j\in J}(a_{i(j)}+\kappa_j),$$
where $C$ contains the factors assigned to $b_L$ and the $a_{i(j)}$ are precisely the unused eigenvalues. For each term form the right local sum of $B'$ and the normal matrix $\operatorname{diag}(W D_{(a_{i(j)})}W^*,D_{a_R})$. Eliminating $R$ shows that its determinant is exactly $\det(T)\prod_j(a_{i(j)}+\kappa_j)$. Apply the second local inclusion to that right sum. Each product multiplied by $C$ pairs every original eigenvalue once. The second distribution may depend on the first pairing. Their resulting nonnegative masses sum to one. This argument works also for repeated eigenvalues of $K$ and zero products.

**A checkable spectral hypothesis.** Suppose $b_{J,R}=\eta+\omega t$ and $a_R=-\eta+\omega h$ coordinatewise for real vectors $t,h$, with $|\omega|=1$. Then $B'=\eta I+\omega H$ for a Hermitian $H$. Whenever $T$ is invertible, $K$ is a scalar shift and rotation of the Hermitian Schur complement of $H_{RR}+D_h$ and is normal. If the local dimensions are unconditionally covered, the singular case follows by adding a common real increment to $h$ and taking a limit. Finite hulls vary continuously with their vertices. The spectra $a_{L,J}$ and $b_L$ can still be arbitrary complex numbers.

For example, a four-dimensional relative unitary with a prescribed zero can be factored across an interface of size two using two three-dimensional blocks. The chosen zero row and column are the two outside coordinates after independent relabeling. Both local dimensions are known. The stated line hypothesis on three right $b$-values and the negative of the outside $a$-value suffices. For general $r$ rows and $c$ columns forming a zero rectangle, the same orthonormal completion gives interface size $n-r-c$, when this number is positive.

**A second sufficient condition.** Assume the left exterior block also is invertible. Eliminate it to obtain the left effective interface matrix $M$. Suppose $M$ and $K$ have a common invariant complete flag, so there is a unitary $W$ in whose basis both are upper triangular. The same two local inclusions suffice. Indeed replacing $K$ by $W D_\kappa W^*$ preserves the full determinant, since the two interface sums have the same triangular diagonal. On the right, adding a matrix diagonal in this basis to the original upper triangular $K$ gives a triangular sum, so its determinant is again the product of diagonal sums. In size two a common eigenvector is enough. Normality of $K$ alone in the first theorem needs no hypothesis on $M$.

## What fails for a general effective matrix

Even the Schur complement of a genuine normal block can carry an off-diagonal entry that changes the global determinant. Let $P$ be the $3\times3$ cyclic shift and $F=(I+P)^{-1}=(I-P+P^2)/2$. Both $F$ and $F^T$ are real normal. Put one on $L,J$ and the other on $J,R$, oriented so their effective triangular entries are on opposite sides of the diagonal. Set the left outside $b$ and the right outside $a$ to zero. Both outside blocks are $[1/2]$. The two effective matrices are
$$\begin{pmatrix}1&-1\\0&1\end{pmatrix},\qquad
\begin{pmatrix}1&0\\-1&1\end{pmatrix}.$$
Each is obtained by inverting the corresponding $2\times2$ principal block of $I+P$ or its transpose. Their sum has determinant three. Replacing either with the diagonal matrix of its eigenvalues would give four. After including the exterior factors the actual determinant is $3/4$, and the replacement gives one.

The original MDO inclusion holds in this example. Both spectra are $0,1/2,(1+i\sqrt3)/2,(1-i\sqrt3)/2$. Pair $0$ and $1/2$ across the two spectra in both orders. Pairing the two remaining numbers with their conjugates gives vertex $1/4$, and pairing them with themselves gives one. The convex weights $1/3$ and $2/3$ give $3/4$. Resolving generic interfaces requires retaining or controlling the off-diagonal contribution. [Exact rational replay](interface.py)

## A convex test for two nonnormal effective matrices

Suppose $|J|=2$, both exterior blocks are invertible, and denote their determinant product by $C$. Let $M,K$ be the left and right effective matrices. Here we assume each local block, in each interface basis used below, has weights matching its determinant polynomial in the **two independent diagonal interface variables**. This holds, for example, when both block dimensions are at most three, by the common squared-minor certificate. It is a stronger input than the one-variable lemma. The bases on the left and right may be chosen independently.

Write $M=m_0 I+l\cdot\sigma$ and $K=k_0 I+r\cdot\sigma$ in the usual Pauli matrices. Here $l,r\in\mathbb C^3$, the dot is complex bilinear, and $m_0,k_0$ are half-traces. Define the centered numerical ranges
$$E_l=\{l\cdot n:n\in\mathbb R^3,\ \|n\|=1\},\qquad E_r=\{r\cdot m:\|m\|=1,\ m\in\mathbb R^3\}.$$
The following sufficient condition includes pairs having no common invariant line:
$$l\cdot r\in\operatorname{conv}\{uv:u\in E_l,\ v\in E_r\}. \tag{2}$$

**Proof of the certificate.** In a chosen left basis, for each local permutation let $c$ be its exterior product and $x_1,x_2$ its two unused $a$ eigenvalues. Write $d_L$ for the left exterior determinant. Equating the four coefficients of the local diagonal pencil gives
$$\mathbb E[c]=d_L,\quad \mathbb E[cx_i]=d_L M_{ii},\quad
\mathbb E[cx_1x_2]=d_L\det M,$$
with $M_{ii}$ taken in that basis. For instance the coefficient of the second interface variable gives the first diagonal entry. The analogous right equalities have $e,y_i,d_R,K$. Pair the unused eigenvalues at corresponding positions across independent local permutations. The expected product of the global factors is
$$\mathbb E[ce(x_1+y_1)(x_2+y_2)]
=C\big(\det M+\det K+M_{11}K_{22}+M_{22}K_{11}\big).$$
In the independent bases the two diagonals are $m_0\pm l\cdot n$ and $k_0\pm r\cdot m$. Every real unit axis occurs as the basis of the Hermitian matrix $n\cdot\sigma$. The constructed points therefore have the form $D_0-2C(l\cdot n)(r\cdot m)$, where $D_0=C(\det M+\det K+2m_0k_0)$. The exact determinant $C\det(M+K)$ equals $D_0-2C(l\cdot r)$ by the Pauli trace identities. Condition (2) makes this a convex combination of the constructed points. All local weights used here are probabilities, even though $c,e,d_L,d_R$ may be complex. ∎

An exact support-function test for (2) is
$$\big|\operatorname{tr} H_\theta\big|\le \|H_\theta\|_{\mathrm{op}}
\quad(0\le\theta<2\pi),\qquad H_\theta=\Re(e^{-i\theta}l r^T).$$
The matrix is real. The largest directional support of the product set is $\max_{n,m}n^T H_\theta m=\|H_\theta\|_{\mathrm{op}}$, and its central symmetry permits the absolute value. Failure of (2) excludes only this construction.

**Defective interfaces.** Suppose the traceless parts are nonzero nilpotents $M_0=\alpha u v^*$ and $K_0=\beta w x^*$, with $u,v$ and $w,x$ orthonormal pairs. Each centered numerical range is a disk of radius $|\alpha|/2$ or $|\beta|/2$: for a unit test vector the expectation is the product of its two basis coordinates with freely chosen relative phase. Moreover $l\cdot r=\operatorname{tr}(M_0K_0)/2$, and $|v^*w|=|x^*u|$ in dimension two. Thus (2) is exactly
$$|v^*w|^2\le 1/2.$$
For the rational example above, rotate the right interface by the real rotation with cosine $3/5$ and sine $4/5$. The squared overlap is $9/25$. Both effective matrices remain defective, their eigenlines differ, and (2) holds strictly. The direct determinant becomes $91/100$. More generally every relative interface rotation satisfying the displayed inequality is covered under the stated local certificate hypothesis.


**Evaluating the convex set without sampling axes.** The image of a real unit sphere in three dimensions under a real linear map to the plane is a filled ellipse (allowing a segment or a point). Each centered range $E_l$ therefore has boundary parameter $e^{i\phi}(a e^{it}+b e^{-it})$ with $a\ge b\ge0$: if its semiaxes are $s_1\ge s_2$, take $a=(s_1+s_2)/2$, $b=(s_1-s_2)/2$. Use $c,d,\psi$ for the second ellipse. The convex hull of products is the rotation through $\phi+\psi$ of the Minkowski sum of the two filled ellipses with boundary curves
$$ac e^{iu}+bd e^{-iu},\qquad ad e^{iv}+bc e^{-iv}.$$
Indeed the product of boundary parameters gives precisely their sum with $u=t+s$, $v=t-s$. Every pair of phases $u,v$ is attainable, and taking convex hulls commutes with Minkowski addition. Products of interior points already lie in this hull by separate real linearity. This also covers a circular or degenerate range. The expression reduces (2) to planar convex geometry with four semiaxis lengths and one rotation.
