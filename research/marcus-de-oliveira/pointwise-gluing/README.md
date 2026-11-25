# One-coordinate gluing for pointwise determinant inclusions

For a fixed unitary $U$, let $P(U)$ mean that $\det(D_a+UD_bU^*)$ lies in the convex hull of the permutation products for **every** pair of complex spectra $a,b$. The weights can depend on both spectra. We prove that $P$ is closed under multiplication of two blocks meeting in one coordinate. This extends the repo's [gluing of simultaneous squared-minor distributions](../single-cycle/README.md#3-one-coordinate-gluing-theorem). The proof below uses the ordinary, pointwise MDO inclusion for each block.

## The affine-pencil lemma

Let $p(z)=p_0+zp_1$ and $v_j(z)=v_{j0}+zv_{j1}$ be finitely many complex affine pencils. Then
$$
p(z)\in\operatorname{conv}_\mathbb R\{v_j(z)\}\quad\text{for all }z\in\mathbb C
\quad\Longleftrightarrow\quad
(p_0,p_1)\in\operatorname{conv}_\mathbb R\{(v_{j0},v_{j1})\}.
$$
A common convex combination of the coefficients suffices by evaluation. Conversely every real linear functional on $\mathbb C^2$ has the form $\Re(\alpha x+\beta y)$. When $\alpha\ne0$, pointwise inclusion at $z=\beta/\alpha$, followed by the real functional $w\mapsto\Re(\alpha w)$, bounds the functional by its maximum on the proposed vertices. For $\alpha=0$ take a sequence of nonzero $\alpha$ tending to zero. A maximum of finitely many continuous functions is continuous. Real convex separation now proves the assertion, including degenerate hulls. Carathéodory's theorem gives a certificate with at most five vertices. A strict separating functional with $\alpha=0$ remains strict after a sufficiently small perturbation of $\alpha$, so failure of coefficient membership always yields a **finite** complex scalar witness.

In particular, all choices of **one** free eigenvalue in MDO are equivalent to a single distribution uniform over that eigenvalue, with every other eigenvalue and the unitary fixed. This implication uses the whole complex plane.

## The gluing theorem

Partition the coordinates as $L,\{v\},R$. Extend a unitary $X$ on $L\cup\{v\}$ and a unitary $Y$ on $\{v\}\cup R$ by the identity. If $P(X)$ and $P(Y)$ hold on their respective blocks, then $P(XY)$ holds. The conclusion also applies to $YX$ by exchanging the names of the sides. Independent row and column permutations and diagonal unitary rephasings preserve $P$. Disjoint supports give the usual product of independent local certificates.

**Proof.** Fix arbitrary spectra $a,b$. Conjugating by $X^*$ gives the matrix $X^*D_aX+YD_bY^*$, whose $L,R$ and $R,L$ blocks are zero. Introduce two local pencils
$$
l(t)=\det\!\left(X^*D_{a_{L,v}}X+\operatorname{diag}(b_L,t)\right)=l_0+t l_1,
\qquad
r(t)=\det\!\left(\operatorname{diag}(t,a_R)+YD_{b_{v,R}}Y^*\right)=r_0+t r_1 .
$$
Their slopes are respectively the determinants of the two outside diagonal blocks of the global matrix. Its determinant is
$$
D=l_0r_1+l_1r_0. \tag{1}
$$
In a Leibniz term, the cycle through $v$ can visit only one side, since distinct sides have no edge between them. If $v$ is a fixed point, its diagonal entry is the sum of the two local entries. Collecting the two contributions gives (1). This argument also covers singular exterior blocks, a zero slope and an empty side.

The left pencil has pointwise MDO for every $t$ by conjugating it back by $X$ and applying $P(X)$. The right pencil has it by $P(Y)$. For each left permutation $\pi$, denote by $i$ the position of the $a$-eigenvalue paired with $t$, and by $C_\pi$ the product of the other factors. Its pencil is $C_\pi(a_i+t)$. Similarly a right permutation $\rho$ leaves one $b_j$ and gives $E_\rho(t+b_j)$. The lemma supplies probability distributions $\lambda,\mu$ reproducing both coefficients of their respective pencils. Equation (1) gives
$$
D=\sum_{\pi,\rho}\lambda_\pi\mu_\rho C_\pi E_\rho(a_i+b_j).
$$
Each product pairs every original $a$ and $b$ exactly once. The left assignments exhaust $b_L$, the right assignments exhaust $a_R$, and the leftover pair supplies the last factor. This is a convex combination of global permutation vertices. No division by a determinant or by a vertex is involved. Each step uses at most 25 local pairings, allowing repeated global vertices. ∎

**Zero-block recognition.** After independent row and column relabelings, $U=XY$ as above if and only if $U[R,L]=0$. Indeed the columns indexed by $L$ are orthonormal in the $(|L|+1)$-dimensional space on $L,v$. Complete those columns to an orthonormal basis of that space, and extend the resulting $X$ by the identity on $R$. Then $X^*U$ fixes every coordinate of $L$ and restricts to a unitary on $v,R$. Thus an $r$-by-$c$ zero rectangle with $r+c=n-1$ splits the fixed-unitary problem into two pointwise problems of sizes $c+1$ and $r+1$. The two smaller unitaries still need proofs of $P$. For example, an order-five unitary with a two-by-two zero rectangle splits into two order-three factors, covered by the known order-three result.

**Recursive schedules and Householder factors.** The argument can be iterated whenever each new block meets the previously assembled coordinates only at one vertex. Drury proved the pointwise determinant inclusion for Householder reflections. For $U_z=I+(z-1)vv^*$ with $v^*v=1$ and $|z|=1$, the determinant interpolates between the identity pairing and the Householder value with parameter $|1-z|^2/4$. Hence products assembled in this schedule can include arbitrary-angle rank-one unitary rotations. A five-coordinate flat Householder already fails the stronger simultaneous squared-minor property. This gives a concrete application for which requiring that stronger property of each factor would exclude a valid factor. No novelty claim is made. [Drury, *OMC for Householder reflections*](https://doi.org/10.1016/S0024-3795(99)00173-1)

## Why the proof stops at one coordinate

A simple example isolates the limitation of the separation step. The polynomial $1+zw$ is in the pointwise convex hull of the twelve polynomials $\eta(1-zw),\eta z,\eta w$ with $\eta\in\{4,-4,4i,-4i\}$. Indeed with $M=\max(|1-zw|,|z|,|w|)$ the square formed by the four multiples of a term of modulus $M$ contains the disk of radius $2\sqrt2 M$. Also $(1+zw)^2=(1-zw)^2+4zw$, so $|1+zw|\le\sqrt5 M$. But all twelve coefficient vectors satisfy “constant coefficient plus $zw$ coefficient equals zero,” while the target sum is two. Thus two free complex coordinates can have pointwise inclusion without any common coefficient distribution. This example concerns abstract pencils. It does not decide the spectrum-constrained interface of two MDO blocks.

The theorem addresses relative unitaries with the indicated factorization. A dense irreducible unitary need not admit it. A two-coordinate interface needs additional information about cross terms or a different proof.

## Exact order-nine example

Two flat order-five Householder matrices meeting at one coordinate give an order-nine unitary. For the nonreal spectra in the checker, each side has a five-permutation certificate for both pencil coefficients. Their 25 products give the global value $(64384+7470336i)/625$. The verifier reconstructs the two pencils, checks every positive rational mass, assembles each full permutation and evaluates the order-nine determinant directly. It also checks (1) with zero left, right and simultaneous slopes.

Run `python3 -B research/marcus-de-oliveira/pointwise-gluing/check.py` from the repository root. [Stored rational certificates](pencil-weights.json) · [Exact checker](check.py)

**A neighborhood of unrestricted relative unitaries and spectra.** For this example the 25 selected points have a six-vertex hull. Orient its edges counterclockwise. The cross product of an edge vector with the vector to the target, divided by the sum of the absolute edge coordinates, is a lower bound on the distance to that edge's supporting line. The least such lower bound is $m=18001547776/550625>0$.

Let $a',b',U'$ vary, with $U'$ any unitary. Bound both spectral changes coordinatewise in complex modulus and the change in $U$ in operator norm by
$$\epsilon=35159273/23229492187500.$$
This includes dense unitaries without the zero rectangle. The original spectra have maximum modulus at most four, so both original and perturbed sums have norm at most ten. The difference between the two matrix sums has norm at most $2\epsilon+2\|D_b\|\epsilon\le10\epsilon$ by adding and subtracting the mixed conjugations. For matrices of norm at most ten, multilinearity in columns and Hadamard's inequality bound the determinant change by $9\cdot 10^8$ times the norm change. Applying the same bound to each diagonal pairing gives vertex movement at most $9\cdot10^8\cdot2\epsilon$. Consequently the slack in every unit support direction drops by at most $9\cdot10^8\cdot12\epsilon=m/2$. The perturbed determinant stays strictly in the hull of the same 25 indexed pairings. This supplies a rigorous joint neighborhood; its certified radius is conservative.
