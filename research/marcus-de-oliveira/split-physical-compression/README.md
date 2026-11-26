# A split eigenvalue pair with an arbitrary complementary flag

Let the first spectrum have three pairs of equal eigenvalues and the second have two pairs and two singletons, in dimension six. All seven spectral values are arbitrary complex numbers. In the basis of the first matrix, write its three two-dimensional coordinate blocks as $B_0,B_1,B_2$. Designate by $E$ the plane of the two singleton eigenvectors of the second matrix. Call this plane balanced if, for the orthogonal projections,
$$
 P_E P_{B_j}|_E=I_E/3\qquad (j=0,1,2).
$$

**Theorem.** For every balanced $E$, every orthonormal ordered basis of $E$, and every orthogonal decomposition of $E^\perp$ into two two-dimensional eigenspaces, the determinant has a convex law on the 33 pairing tables. The law depends only on the relative eigenbasis and works simultaneously for all seven spectral values. The same conclusion holds for every unitary within entrywise distance $1/57480192$ of this entire family. In particular the fourth eigenvalue can move arbitrarily across the complex plane.

This includes the entire complex Grassmannian of two-planes in a four-dimensional complement, for every orientation of the two singleton eigenvectors. The proof allows the weights to change with the orientation. It leaves the general case of an arbitrary singleton eigenspace open.

## Reduce the flag using the balanced plane

The map $v\mapsto\sqrt3 P_{B_j}v$ from $E$ to each $B_j$ is unitary. A separate choice of basis in each $B_j$ therefore sends any prescribed ordered orthonormal split basis to the columns of $X/\sqrt3$, with $X=[I_2;I_2;I_2]$. Such changes preserve the canonical table coefficients: the determinant terms grouped by physical blocks are coefficients of a block-scaled determinant, and block determinants cancel in the ratio. Order the columns as singleton, singleton, first pair, second pair. This column permutation has positive sign; any shared sign would cancel in the normalized determinant in any event.

Take $Y=[-W^T;I_4]$ with $W=[I_2;I_2]$. It spans $E^\perp$. Put $G=Y^TY=I_4+WW^T$; $\det G=9$ and the eigenvalues of $G$ are $3,3,1,1$. Choose any full-rank four-by-two matrix $Z$ to represent the first complementary eigenspace by $YZ$. Let $z$ be its six two-minors in increasing coordinate order $01,02,03,12,13,23$. Write $g_z=\det(Z^*GZ)=z^*G^{[2]}z$, so $0<g_z\le9\|z\|^2$. Superscript $[2]$ takes all two-by-two minors.

The real symmetric Hodge matrix $S$ has $S_{I,J}=\operatorname{sgn}(J,I)$ when $I=J^c$ and zero otherwise. Every such $z$ satisfies the complex Plücker equation $z^TSz=0$. Set $v_K=(Y^{[2]})_{K,:}$ and $u_L=(Y^{[2]}S G^{[2]})_{L,:}$ for row pairs $K,L$ of the six-dimensional matrix.

For a real table score $h$, write $h(K,L,i,j)$ for its value when the two paired columns take row sets $K,L$ and the singleton columns take distinct rows $i,j$. Sum over all ordered partitions of the six rows of this form, with increasing $K,L$. Form the integer matrix
$$
 L_h=\sum \operatorname{sgn}(i,j,K,L)\,h(K,L,i,j) X_{i0}X_{j1}\,v_K^T u_L,
 \qquad T_h=3(L_h+L_h^T).
$$
The physical value on the canonical signed table coefficients $q$ is exactly
$$
 h\cdot q=\frac{z^*T_h z}{54g_z}.
$$
Indeed the two normalized singleton columns supply $X_{i0}X_{j1}/3$. Normalize the first pair in the $G$ metric, which replaces $z$ by $z/\sqrt{g_z}$ after a column phase choice. Choose the phase of the last pair so that the induced four-by-four unitary has determinant one. Jacobi's complementary-minor identity then gives the last pair's minors as $S G^{[2]}\bar z/(3\sqrt{g_z})$ before applying $Y^{[2]}$. The determinant of the whole orthonormal frame is one because $\det[X,Y]=9$. Multiplying the factors and taking real parts gives the displayed expression. The derivation of the metric Hodge identity is in the [paired-spectrum Grassmann proof](../double-paired-five-term/GRASSMANN.md). The verifier also compares each of the 33 basis scores with the original 720-term determinant on three nonreal flags.

## A phase choice on the Plücker quadric

For a real symmetric matrix $T$, choose a real $s$ and a coordinate $k$. If
$$
 T+sS\succeq cI,\qquad
 (T-sS)|_{e_k^\perp}\succeq cI,
$$
then $z^*Tz\ge c\|z\|^2$ for every complex $z$ with $z^TSz=0$. Multiply $z$ by a scalar phase so that $z_k$ is real, and write $z=x+iy$. Then $y_k=0$ and $x^TSx=y^TSy$. Adding $s$ times their difference leaves the value unchanged and gives the two lower bounds. This proof also covers $z_k=0$.

The [complete split-pattern dual classification](../split-dual-rays/README.md) contains 935 rays in 62 physical orbits. Ten representatives are positive multiples of grouped squared-minor coefficients and hence are always nonnegative. For representative 3, merge physical row colors 1 and 2 in the original table convention. The verifier checks that the resulting rank-12 coefficient row space contains this nonnegative score. Thus it is constant on each coarse table fiber and its evaluation is determined by the coarse spectral polynomial. The merged spectrum has multiplicities $(2,3,1)$; [Claude's independently reviewed Theorem D](../minors/theorem-d.md) provides a nonnegative common law for that polynomial. Evaluating the coarse score against that law proves the desired sign for every unitary.

For each of the 51 remaining representatives the [exact checker](check.py) proves the phase criterion with $c=1/4$ by checking positive integer leading principal minors after subtracting the margin. The five nonzero shifts, in the enumeration's zero-based IDs and coordinates, are below. The other representatives use $s=0$.

| Orbit | Shift $s$ | Real coordinate $k$ |
| ---: | ---: | ---: |
| 8 | -4 | 0 |
| 11 | 8 | 0 |
| 14 | -16 | 5 |
| 46 | 21 | 0 |
| 47 | -29 | 5 |

Ordinary positive semidefiniteness fails for exactly 14, 46 and 47 on this plane. The phase argument proves actual positivity there. The physical lower bound for all 51 is $1/(4\cdot54\cdot9)=1/1944$. All physical relabelings preserve the family: the three row pairs can permute, the two singleton columns can swap, and the two complementary pair colors can swap, exchanging the chosen two-plane with its orthogonal complement. Thus all 935 signs hold. Finite-dimensional cone duality in the classification gives nonnegative masses $m$ with $Am=Aq$, where $A$ is the spectral coefficient matrix. Normalization is among the consequences. Assigning each table's mass to any realizing permutation produces the claimed spectral law.

## An open neighborhood and an intrinsic sufficient condition

Take a relative unitary $V$ in this family and another relative unitary $U$ whose entries all differ by at most $\epsilon$ in modulus of the difference. Each canonical table coefficient changes by at most $672\epsilon$ by the [split-pair Lipschitz estimate](../split-one-eigenvalue-pair/README.md). The largest primitive ray has coordinate sum 22, giving loss at most $22\cdot672\epsilon$. At $\epsilon=1/57480192$ every strictly controlled ray retains at least $1/3888$. The other eleven orbits have global nonnegativity from the preceding argument. This proves the entrywise tube, including across changes of determinant phase. Transposition or adjoint of the relative unitary merely exchanges the two spectral sides and preserves the entrywise metric.

It also suffices that the orthogonal projection onto the actual split plane is within operator distance $1/114960384$ of some balanced split-plane projection. To see this, rotate principal vector pairs to carry the actual plane to the balanced one. For the largest principal angle $\theta$, the rotating unitary differs from the identity in operator norm by $2\sin(\theta/2)\le2\sin\theta$, and $\sin\theta$ is the projection distance. Applying this rotation to the entire eigenframe produces a member of the balanced family within the required entrywise distance. This leaves both complementary eigenspaces arbitrary.

## Exact replay

Run from the repository root. The verifier uses the standard library.

```sh
python3 -O -I -B research/marcus-de-oliveira/split-physical-compression/check.py
```

Independent mathematical review is requested, particularly for the sign conventions and the use of coarse weights for orbit 3. The input ray enumeration and the separate paired-spectrum proof have their own review records. The exact computations here establish a split-plane fiber and a neighborhood; they do not establish positivity on every split plane.
