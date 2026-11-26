# A common determinant law when both spectra have three double eigenvalues

**Theorem.** For every complex unitary $U$ of size six there are 21 nonnegative real weights summing to one, depending on $U$ alone, such that for all $a,b\in\mathbb C^3$
$$
\det(D_a+U D_b U^*)=\sum_{i\le j}w_{ij}p_i p_j,\qquad
D_a=\mathrm{diag}(a_0,a_0,a_1,a_1,a_2,a_2),\qquad
p_i=\prod_{r=0}^2(a_r+b_{\pi_i(r)}).
$$
Here $D_b$ repeats the entries of $b$ in the same way and $\pi_i$ ranges over $012,120,201,021,102,210$. The entries within either triple may coincide. Each $p_i p_j$ is an actual six-dimensional pairing product. This proves the Marcus–de Oliveira inclusion for these two multiplicity patterns. The unrestricted conjecture and the case where only one spectrum has this pattern remain open.

The mathematical evidence is a written reduction and a finite exact polynomial proof. A different GitHub owner should independently review the unitary chart, the Hodge orientation, and the four-facet reduction. The repository's Lean forest results concern a different family.

## The four signs that suffice

Use the 21 real canonical signed coefficients $q_{ij}$ obtained by grouping the real parts of the signed Leibniz terms $\overline{\det U}\,\mathrm{sgn}(\sigma)\prod_k U_{k,\sigma(k)}$. They sum to one. Their pairing with the 21 formal product polynomials gives the determinant. The coefficient of each spectral monomial is real by the squared-minor expansion, which justifies taking real parts before choosing complex spectra.

The earlier [contingency-table theorem and shortest-path construction](../double-paired-contingencies/README.md) classifies the nonnegative dual rays of this polynomial-weight problem. The 189 rays form seven orbits under independent row- and column-pair relabeling and transposition of the table. These operations preserve the set of physical problems, with transposition exchanging $U$ and $U^*$. Three orbits are nonnegative by Cauchy–Binet. It suffices to prove one representative of each of the four other orbits for **every** $U$. With $ij$ denoting $q_{ij}$, these representatives are

| Orbit | Functional |
| --- | --- |
| 4 | $12+15+24+44+55$ |
| 5 | $2(12)+2(15)+23+24+34+2(55)$ |
| 6 | $02+05+11+14+23+34+55$ |
| 7 | $02+05+12+14+2(23)+2(33)+45$ |

For any such functional, let $h(I,J,K)$ be the coefficient it assigns to the table of row counts of the three ordered column pairs. Here $I,J,K$ partition the six rows into sets of size two, and each set is listed in increasing order. Let $\epsilon(I,J,K)$ be the sign of their concatenation. Splitting the determinant into the three pairs gives the functional as the real part of a sum of $h\epsilon$ times the three two-column minors and $\overline{\det U}$. This identity also shows that the grouped functional is invariant under independent unitary changes of basis within each of the three column pairs. On the row side, sum first over terms with prescribed column sets for each row pair. Each row change of basis then contributes one block determinant, canceled by $\overline{\det U}$. Equivalently, the grouped coefficients are those of the six-by-six determinant with nine independent variables scaling its two-by-two blocks.

## Put the first column plane in real coordinates

First assume its projection onto the first row pair is invertible. Choose a basis for that plane to write it as the column span of $[I;D;M]$, initially with two complex two-by-two lower blocks. Left actions from the three row-pair unitaries act on both lower blocks on the right by the same unitary and on the left independently. An SVD diagonalizes $D$ with real nonnegative entries. A left QR change makes $M$ upper triangular with real nonnegative diagonal, including by a unitary completion when singular. Simultaneous conjugation by a diagonal unitary fixes $D$ and the diagonal of $M$ and makes its off-diagonal entry real. These basis choices leave every grouped functional unchanged.

Write the resulting real frames as
$$
D=\begin{pmatrix}a&0\\0&b\end{pmatrix},\quad
M=\begin{pmatrix}c&d\\0&e\end{pmatrix},\quad
X=\begin{bmatrix}I\\D\\M\end{bmatrix},\quad
Y=\begin{bmatrix}-D^T&-M^T\\ I&0\\0&I\end{bmatrix}.
$$
The five real letters in this chart are independent of the spectral triples in the theorem. Put $S=X^T X$, $G=Y^T Y$. The columns of $Y$ span the orthogonal complement of $X$. Writing $W=[D;M]$ gives $S=I+W^TW$, $G=I+WW^T$, and $\det[X,Y]=\det S=\det G>0$ by a block determinant. Normalize the first pair to $XS^{-1/2}$ and write the other pairs $Y V_1,Y V_2$. Then $C=G^{1/2}[V_1,V_2]$ is unitary of size four. Changing the phase of a column in the last pair allows $\det C=1$. In this gauge $\det U=1$ as well.

## Compress to a quadratic form on six coordinates

Use the increasing pairs $01,02,03,12,13,23$ of the four complement coordinates. Write $A^{[2]}$ for the matrix of all ordered two-by-two minors of $A$. The real Hodge matrix, extended complex-linearly, has $*_{P,Q}=\mathrm{sgn}(Q,P)$ if $P=Q^c$ and is zero otherwise. Let $z$ be the vector of two-by-two minors of $V_1$. Jacobi's complementary-minor identity on $C\in SU(4)$ identifies the last pair with the Hodge image of the conjugate of the first pair. Pulling it back through $G^{1/2}$ gives
$$
\bigwedge^2 V_2=(\det G)^{-1/2}\,*\,G^{[2]}\overline z .
$$
Indeed for any invertible real $H$ the identity $(H^{-1})^{[2]}*=(\det H)^{-1}*(H^T)^{[2]}$ follows from the complementary-minor formula. Take $H=G^{1/2}$ and use the multiplicativity of two-by-two compounds. This step retains the physical first column plane, unlike the earlier relaxation of the entire third compound.

For a row pair $I$, set $x_I=\det X[I,:]$, $v_I=Y^{[2]}[I,:]$ and $t_I=(Y^{[2]}*G^{[2]})[I,:]$. Define
$$
L=\sum_{I,J,K} h(I,J,K)\epsilon(I,J,K)x_I\,v_J^T t_K,\qquad T=L+L^T.
$$
The original grouped functional equals $z^*Tz/(2\det S)$. To see the factors, the first minor contributes $x_I/\sqrt{\det S}$, and the other two contribute $v_Jz$ and $t_K\bar z/\sqrt{\det G}$. Their real part symmetrizes $L$. The denominator is positive. It will suffice to show that $T$ is positive semidefinite on the whole six-dimensional space.

## Exact polynomial certificate

For each representative, form $T$ as above in the polynomial ring $\mathbb Z[a,b,c,d,e]$ and expand its six leading principal minors. Every monomial of each of the 24 polynomials has positive integer coefficient and even exponents in all five variables. Here are the counts of distinct monomials after collecting all terms.

| Orbit | Order 1 | Order 2 | Order 3 | Order 4 | Order 5 | Order 6 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 7 | 47 | 166 | 475 | 1,111 | 1,545 |
| 5 | 3 | 33 | 147 | 558 | 1,386 | 2,834 |
| 6 | 14 | 99 | 282 | 755 | 1,344 | 2,557 |
| 7 | 15 | 120 | 381 | 1,329 | 2,964 | 5,568 |

The [standard-library verifier](grassmann.py) constructs the matrix from signed set partitions and calculates each minor by Laplace expansion. It inspects every coefficient and exponent and checks the [fingerprints of the full polynomials](grassmann-results.json). The maximum total degree of an entry is ten, so a determinant intermediate has total degree at most sixty. The verifier's five base-64 exponent digits cannot carry into adjoining coordinates.

When all five variables are nonzero, all six minors for each matrix are **strictly** positive. Sylvester's criterion makes each real symmetric matrix positive definite there. Polynomial continuity gives positive semidefiniteness everywhere on the real chart. Finally, unitaries with invertible upper two-by-two block are dense in $U(6)$, since the squared modulus of its determinant is a real analytic function on a connected real analytic manifold and is nonzero at the identity. The original functionals are continuous on all unitaries. This extends all four signs across the missing chart boundary.

Applying the physical symmetries and the earlier dual classification now proves feasibility of nonnegative weights for every relative unitary. The [12-vertex solver](../double-paired-contingencies/check.py) constructs them from the canonical coefficients and uses no spectral parameters. A representative permutation realizing each table transfers its mass to an actual six-coordinate pairing if a distribution on permutations is desired.

## Reproduce and challenge

From the repository root run:

```sh
python3 -I -B research/marcus-de-oliveira/double-paired-contingencies/check.py
python3 -I -B research/marcus-de-oliveira/double-paired-five-term/grassmann.py
python3 -I -B research/marcus-de-oliveira/double-paired-five-term/grassmann_fixtures.py
```

The [complex frame checks](grassmann_fixtures.py) compare the expression with all 720 signed Leibniz terms on four Gaussian-rational graph frames, including a singular chart parameter and nonreal determinant phases. Normalizing a frame $Z$ on the second plane divides the quadratic form by $\det(Z^*GZ)$. The direct calculation groups by table, divides by the complex frame determinant, and then takes real parts. Block orthogonality ensures that positive Gram normalizations cancel between the numerator and determinant. Each example has negative canonical entries and reconstructs a nonnegative law with all 141 polynomial coefficients checked. These fixtures test the conventions. The universal sign proof uses the 24 exact polynomials. Both calculations also replay under Python's optimized mode.
