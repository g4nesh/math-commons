# A Pfaffian certificate for three quaternionic coordinates

Let both six-dimensional spectra consist of three arbitrary complex values repeated twice. Write $D_a=\operatorname{diag}(a_0 I_2,a_1 I_2,a_2 I_2)$ and likewise $D_b$. Put $J=\operatorname{diag}(J_2,J_2,J_2)$ with $J_2=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)$. If the relative unitary $U$ is symplectic, $UJU^\mathsf T=J$, **one nonnegative distribution on the 21 distinct pairing products works for all six spectral values simultaneously**. The result also holds after independent unitary changes of basis inside any of the row or column pairs. At every symplectic center with all nine two-by-two blocks nonzero, an explicit neighborhood of unrestricted complex unitaries has the same existence property.

## Pfaffian expansion in any even dimension

Take $m$ pairs and $U\in U(2m)$ with $UJU^\mathsf T=J$. From $U^*J=JU^\mathsf T$ we get
$$K=(D_a+UD_bU^*)J = D_a J+U D_b J U^\mathsf T.$$
This is skew-symmetric over $\mathbb C$ for arbitrary complex $a,b$. Its determinant is the determinant of the original sum since $\det J=1$. The identity $\det K=\operatorname{Pf}(K)^2$ supplies a polynomial square with no choice of branch.

For $I\subset\{0,\ldots,m-1\}$ let $I^\uparrow$ contain both complex coordinates of each pair in $I$. With columns $u_k=Ue_k$, the two-form of $K$ is
$$\sum_i a_i\,e_{2i}\wedge e_{2i+1} + \sum_j b_j\,u_{2j}\wedge u_{2j+1}.$$
The coefficient of the standard volume form in its $m$th exterior power divided by $m!$ is the Pfaffian. Selecting the standard basis pairs outside $I$ and the transformed pairs on $T$ gives
$$\operatorname{Pf}(K)=\sum_{|I|=|T|} a^{I^c}b^T \det U[I^\uparrow,T^\uparrow].$$
Moving entire pairs across each other introduces even signs. The empty minor is one. This is the paired-block instance of Pfaffian minor summation. [Ishikawa, *Minor summation formula of Pfaffians and Schur function identities*](https://fpsac-archive.github.io/FPSAC95/ARTICLES/38.ISHIKAWA.pdf)

## A positive square root on three pairs

The unitary symplectic condition also reads $UJ=J\overline U$. Every $2\times2$ block $U_{rs}$ therefore has the form $\left(\begin{smallmatrix}z&w\\-\bar w&\bar z\end{smallmatrix}\right)$. Set $P_{rs}=\det U_{rs}=|z|^2+|w|^2$. The sum over a block row is one by the norm of its first ordinary row, and the column sums are one by the first ordinary column. Hence $P$ is doubly stochastic.

Pfaffian congruence gives $\mathrm{Pf}(UJU^\mathsf T)=(\det U)\mathrm{Pf}(J)$ and $\mathrm{Pf}(J)=1$, so $\det U=1$. By the complementary-minor identity, a two-pair minor for $m=3$ equals the complex conjugate of the one-block minor on the complementary pair of indices. Its sign is positive because the sums of the row and column index lists have the same parity. Since $P$ is real, the two-pair coefficients are precisely the complementary entries of $P$.

Choose nonnegative Birkhoff weights $\lambda_\pi$ with $\sum\lambda_\pi=1$ and $P_{rs}=\sum_{\pi:\pi(r)=s}\lambda_\pi$. One permutation then matches every coefficient of the three-pair Pfaffian in expectation: sizes one and two follow from $P$ and its complement, and sizes zero and three have mass one. Thus, with $p_\pi=\prod_r(a_r+b_{\pi(r)})$,
$$\operatorname{Pf}(K)=\sum_{\pi\in S_3}\lambda_\pi p_\pi,\qquad
\det(D_a+UD_bU^*)=\sum_{\pi,\rho\in S_3}\lambda_\pi\lambda_\rho\,p_\pi p_\rho.$$
Pair one of the two rows in each block with a column of block $\pi(r)$ and the other with a column of block $\rho(r)$. This uses each column block twice, and its two slots can be assigned bijectively. Hence every summand is a genuine six-coordinate pairing. Collecting unordered pairs gives diagonal weights $\lambda_\pi^2$ and off-diagonal weights $2\lambda_\pi\lambda_\rho$. In particular all four [remaining universal signs](https://github.com/g4nesh/math-commons/pull/81) hold on this class. Row and column block changes of basis commute with $D_a$ and $D_b$ respectively, which proves the stated invariance.

## A neighborhood beyond symplectic matrices

Let $h=\min_{r,s}P_{rs}>0$. Reserve mass $h/2$ on each of the six permutations. Each matrix entry is used by two permutations, leaving the nonnegative matrix $P-h\mathbf1\mathbf1^\mathsf T$ with row and column sums $1-3h$. Decompose its normalized residual by Birkhoff if that total is positive. This gives $\lambda_\pi\ge h/2$ and a base law $w^U$ with every unordered-pair weight at least $h^2/4$.

For an arbitrary six-dimensional unitary $V$, aggregate its canonical signed three-color masses into the 21 row/column contingency tables, giving $q(V)$. Set $w^V=w^U+q(V)-q(U)$. This works coefficientwise for both paired spectra. For an entrywise maximum change $\eta$, each canonical coloring mass varies by at most $30\eta$: its formula contains one six-minor and three two-minors of modulus at most one, and a $k$-minor varies by at most $k\sqrt{k}\eta$ by column telescoping. At most eight colorings aggregate to any table, since the two rows in each of the three pairs have at most two orderings. Thus each table mass varies by at most $240\eta$. Every unitary with
$$\max_{i,j}|V_{ij}-U_{ij}|\ \le\ h^2/1920$$
admits the displayed law with each weight at least $h^2/8$. A nonsymplectic matrix may belong to this neighborhood. [Canonical signed coloring formula](../paired-matching-boundary/README.md#4-keeping-the-labels-when-symmetry-breaks)

## The first limit of the square-root argument

For four pairs take $U=R\otimes I_2$ with the exact orthogonal four-cycle $R$ in the repository. This $8\times8$ matrix is symplectic. Its paired minors are $\det(R[I,T])^2$. A common four-coordinate Pfaffian law for arbitrary $a,b$ would match all these squared minors by coefficient comparison. The archived cover is at least one on each of the 24 permutations and equals $14953/15625$ on these minors. Thus the Pfaffian's common-law conclusion already fails at four pairs. The squared determinant might still admit a different law. [Exact four-cycle obstruction](../single-cycle/README.md#5-a-four-cycle-with-no-common-minor-certificate)

## Replay

Run `python3 -I -B research/marcus-de-oliveira/quaternionic-paired-spectra/check.py`. The standard-library checker verifies a product of gates with noncommuting quaternionic phases, the Pfaffian against the ordinary determinant, all paired coefficients against Birkhoff weights, and the full symplectic and perturbed laws against 924 squared minors in each fixture. It also replays the four-pair obstruction directly. [Exact checker](check.py)

The center has $h=334864/6076225$ and certified radius $876046082/553807653759375$. The perturbation has a two-by-two block with distinct singular values, which rules out even the blockwise change-of-basis form of a symplectic matrix. Both matrices pass all 141 distinct full-determinant coefficient checks.
