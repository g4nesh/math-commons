# Which pole deletions preserve normality?

In the eigenbasis of a normal matrix $B$, choose a simple eigenvalue $b$ with unit eigenvector $e$. Let $A$ be normal and $A+bI$ invertible. Deletion at $b$ expresses the determinant as a scalar times a determinant one dimension smaller. The following criterion decides whether the smaller matrix is normal, using only the eigenvalues of $A$ visible from $e$. The order-four consequence gives one set of weights that works as the other three eigenvalues of $B$ move independently in the complex plane. [Pole identity and injective hull](../pole-compression/README.md)

## The deleted direction

Write $G=(A+bI)^{-1}$, $Q=I-ee^*$ and $C=QG|_{e^\perp}$. For each distinct eigenvalue $\gamma$ of $G$ let $E_\gamma$ be its orthogonal spectral projection and $p_\gamma=\|E_\gamma e\|^2$. Define

$$
\mu=\sum_\gamma p_\gamma\gamma,\quad
V=\sum_\gamma p_\gamma|\gamma-\mu|^2,\qquad
P=\sum_\gamma p_\gamma(\gamma-\mu)^2 .
$$

**Criterion.** $C$ is normal exactly when $\{\gamma:p_\gamma>0\}$ lies on a real affine line in $\mathbb C$. For every deletion, including the nonnormal case,
$$\|CC^*-C^*C\|_{\rm op}^2=V^2-|P|^2.$$

Indeed normality of $G$ gives $CC^*-C^*C=vv^*-uu^*$, where $u=QGe$ and $v=QG^*e$. Both have squared norm $V$ and $v^*u=P$. The difference has rank at most two and trace zero. Its squared Frobenius norm is $2(V^2-|P|^2)$, giving the formula even if both vectors vanish. For $V>0$ equality to zero means $u=\zeta v$ for $|\zeta|=1$. Projecting this identity onto each eigenspace gives $\gamma-\mu=\zeta(\overline\gamma-\overline\mu)$ on the support, precisely the equation of a line through $\mu$. Conversely any such line gives the identity. A one-point support has $V=0$.

Möbius inversion therefore gives the criterion directly on the original spectra. All eigenvalues $a$ of $A$ whose eigenspaces meet $e$ in a nonzero projection must lie on one generalized circle through $-b$. A generalized circle may be a line. Eigenvalues orthogonal to $e$ are unrestricted. When normality fails, the two nonzero eigenvalues of the self-commutator are $\pm\sqrt{V^2-|P|^2}$. This number measures a failure of the reduction hypothesis, with no implication that the original hull inclusion fails.

The supporting-line geometry is classical. Fan and Pall's codimension-one normal compression theorem describes interlacing on a line after common eigenvalues are removed. The spectral measure identifies the appropriate reducing summand for the specified vector, including shared eigenvalues. Their higher-codimension observation on p. 299 matters for iteration. A normal final compression alone does not guarantee a chain of normal deletions. Our application combines the codimension-one geometry with the determinant and injection formulas. [Fan and Pall, Theorem 2, p. 298](https://doi.org/10.4153/CJM-1957-036-1) · [Holbrook, Mudalige and Pereira, Proposition 8](https://files.ele-math.com/articles/oam-07-47.pdf)

## Transferring the weights

The cyclic space $K=\operatorname{span}\{E_\gamma e:p_\gamma>0\}$ reduces $G$. On $K$, write $G=\eta I+\omega H$ with $H$ Hermitian, choosing $|\omega|=1$ along the support line. On $K^\perp$ the matrix passes unchanged into $C$, including the remaining copies of any multiple eigenvalue. On $K\cap e^\perp$ the reduced eigenvalues are affine images of an interlacing Hermitian list. This includes the one-point case, where that part disappears.

For the active list in increasing real line coordinates, assign the largest reduced value to either of the last two original indices with probabilities giving its mean. Delete that index. Both remaining lists interlace, so induction gives an injective assignment with the prescribed mixed products for the smaller list in each conditional branch. Multiplying by the last assigned factor and taking expectations completes the induction. Pair the inactive indices deterministically. This constructs a law on injections $J$ such that for the eigenvalues $\alpha_i$ of $C$,
$$\mathbb E\prod_i(\gamma_{J(i)}+t_i)=\prod_i(\alpha_i+t_i)\qquad(t_i\in\mathbb C).$$
The affine change of coordinate only rescales the identity and translates its free parameters. For any lower-dimensional pairing, permute the labels in the identity, then average over the lower-dimensional certificate. The unused original eigenvalue pairs with $b$. The determinant factorization proves MDO upstairs whenever the reduced pair has a certificate. This uses conditional mixed products, which are stronger than interpolating each marginal. [Interlacing and lifting in Shitov, Lemma 17](https://www.mathnet.ru/eng/im9292)

**Order four.** Fix $A,e,b$ and fix the three other eigenvectors of $B$. Under the criterion, the reduced $3\times3$ pair has one permutation law valid at all three remaining spectral parameters. To see this, decompose the doubly stochastic matrix of squared entries of a unitary eigenbasis of $C$ by Birkhoff's theorem. Squared order-two minors equal the complementary squared entries by unitarity. Expanding the determinant in its three diagonal variables therefore matches every coefficient with this same law. The interpolation above lifts it without using any of those three variables. When an additional $b_j$ equals $b$, the identity extends by polynomial continuity. The eigenbasis on $e^\perp$ can also be chosen arbitrarily before constructing the law.

**Moving the fourth eigenvalue on a circle.** If a row of the relative order-four eigenbasis contains a zero, the corresponding direction $e$ sees at most three distinct eigenvalues of $A$. Choose a generalized circle $\Gamma$ through those eigenvalues. For each finite $b\in-\Gamma$, the preceding construction gives a permutation vector valid for every choice of the other three complex eigenvalues of $B$. When $e$ sees at most two eigenvalues, every $b\in\mathbb C$ belongs to such a reflected circle. At a pole $a_i+b=0$, approach along the circle and take a convergent subsequence of the weights in the compact probability simplex. They match each of the finitely many polynomial coefficients in the limit, so the limiting weights still work for all three free variables. For three distinct active points the circle is unique. We assert sufficiency along this circle; points off it can admit other weights. In the fixture below all three active eigenvalues of $A$ are real, hence the moving parameter $b$ runs over the entire real axis. Each $b$ can have its own weights.

**Higher dimensions.** Both the criterion and the transfer are dimension free. Each successive deletion requires its own normality test and an invertible shifted matrix. Claude's [merged sum-of-squares certificate](https://github.com/g4nesh/math-commons/pull/34) enlarges the last four-dimensional step in the [earlier lift](../pole-compression/README.md) to a sorted line-spectrum cross-ratio in $[1/10,9/10]$. This consequence inherits the mathematical review status of the cross-ratio argument. The criterion here also handles intermediate normal matrices with noncollinear full spectra. An arbitrary normal four-dimensional remainder with noncollinear spectrum still requires its own proof.

## The archived four-cycle

Use the ordering in the [four-cycle report](../fourth-boundary/README.md), with $A$ diagonal and $B=U\operatorname{diag}(3-2i,2,1+2i,-i)U^*$. Exchange the roles of $A$ and $B$ when applying the theorem. Row two of $U$ sees precisely the first three eigenvalues, which lie on the line $2+\mathbb R(1-2i)$. Therefore, for every real $t$, fixing $a_2=-2+t(1-2i)$ (one-based index) gives weights for **all** $(a_1,a_3,a_4)\in\mathbb C^3$. This allows the fourth eigenvalue to move anywhere in the plane. It is complementary to fixing the fourth coordinate and leaving the other three unrestricted. The scopes of the two claims are different.

Five exact laws are included at $t=-1,0,1,-1/2,2$. The first three are singular pole values and have, respectively, nine, nine and seven positive masses. The other two have eleven and ten. Each archived permutation word maps the four columns of $U$ to the one-based diagonal rows of $A$. We found the masses by rational phase-I linear programming and replay them independently against the original determinant on a unisolvent eight-point grid. [Rational masses](archive-slices.json) · [Independent replay](archive_slices.py)

## An exact nonconcyclic example and a failure control

Start with the resolvent eigenvalues $1,3,8,4+i$. Apply the three-coordinate Householder $I-\tfrac23\mathbf{1}\mathbf{1}^T$ on coordinates $(1,2,3)$ and then on $(2,3,4)$. The deleted first coordinate sees squared weights $1/9,4/9,4/9,0$. The compression roots are $4/3,17/3,4+i$. Use $A=G^{-1}$ and $B=\operatorname{diag}(0,1+i,2-i,-1+3i)$. Their spectral cross-ratios in the displayed ordering are $(-3+5i)/7$ and $1/2+i$, so neither spectrum lies on a generalized circle.

The checker constructs sixteen positive rational masses from spectral projectors, a Birkhoff decomposition and conditional interlacing. It evaluates both sides at all eight vertices of $\{0,1\}^3$ to certify the multi-affine identity for every triple of complex eigenvalues with the same weights. At the indicated $B$ the determinant is $-19955/16524+(116525/33048)i$. Changing the conjugating matrix to the four-coordinate flat Householder makes the deleted coordinate see all four points. The squared operator norm of the resulting self-commutator is $39/8$. [Exact checker](check.py) · [All sixteen masses](certificate.json)

From the repository root, run

~~~sh
python3 -B research/marcus-de-oliveira/normal-pole-step/check.py | diff - research/marcus-de-oliveira/normal-pole-step/certificate.json
python3 -B research/marcus-de-oliveira/normal-pole-step/archive_slices.py
~~~
