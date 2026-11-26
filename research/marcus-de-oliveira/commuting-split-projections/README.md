# Different projection proportions for the two singleton eigenvectors

In dimension six let one normal matrix have a spectral pattern $(2,2,1,1)$ and the other $(2,2,2)$. Write $u,v$ for the two singleton eigenvectors and $P_j$ for the projections onto the three opposite spectral blocks. All seven eigenvalues may range independently over the complex plane. The following conditions concern their relative eigenspaces.

1. **All aligned pairs.** If $\langle P_j u,P_j v\rangle=0$ for every $j$, a single law on the 33 pairing tables gives the determinant for all spectral values. The two probability triples $(\|P_j u\|^2)_j$ and $(\|P_j v\|^2)_j$ are arbitrary and independent. The orthogonal complement of the singleton plane can be split into its two eigenspaces in every possible way.
2. **All orientations at a nonscalar plane.** Suppose the singleton plane has orthonormal axes with opposite projection triples $(1,1,9)/11$ and $(1,4,16)/21$, and the axes project orthogonally in each block. Then every choice of singleton eigenvectors in that plane and every split of its complement admits such weights. The weights may depend on the orientation. This also holds on an open unitary neighborhood of the whole family.

The aligned condition says the singleton eigenvectors diagonalize the three compressed projections simultaneously. It includes the [isoclinic theorem](../isoclinic-split-planes/README.md). The second result includes singleton directions that do not diagonalize those projections. The combination leaves the general case of arbitrary commuting compressed projections and arbitrary singleton orientation open.

## Four independent parameters

Choose a basis within each opposite block so the two projected axes are nonnegative multiples of its first and second coordinate vectors. When both have nonzero projection in block zero, use the two-column frame
$$
 X=[I;\operatorname{diag}(a,b);\operatorname{diag}(c,d)],\qquad
 r_0=1+a^2+c^2,\quad r_1=1+b^2+d^2.
$$
The two axes are the columns of $X$ divided by $\sqrt{r_0}$ and $\sqrt{r_1}$ respectively. For $W=X[2{:},:]$, take $Y=[-W^T;I_4]$, $G=Y^TY$. The block determinant identities give $\det[X,Y]=\det G=r_0r_1$. As in the [balanced-plane reduction](../split-physical-compression/README.md), write $z$ for the six minors of an arbitrary full-rank four-by-two complement frame $Z$, and $g_z=\det(Z^*GZ)>0$. The physical ray functional is
$$
 h\cdot q=\frac{z^*M_hz}{2r_0r_1g_z}.
$$
Here $M_h=L_h+L_h^T$ uses exactly the ordered partitions and Hodge matrix $S$ from the isoclinic reduction, now with $X_{i0}X_{j1}$ from these four independent variables. Normalizing the singleton columns supplies $1/\sqrt{r_0r_1}$; the metric complementary-minor identity supplies the identical factor a second time. The [checker](check.py) builds every matrix from these partitions and also checks all 33 coefficients directly against the original 720-term determinant on two aligned complex flags and one rotated flag.

The [dual classification](../split-dual-rays/README.md) has 62 physical orbits covering 935 rays. Ten are always nonnegative by the squared-minor expansion. For orbit 3, the checker verifies nonnegativity of all 63 principal minors. For each of the other 51 non-coefficient orbits, it proves $M_h+sS$ positive definite and $M_h-sS$ positive definite after deleting one row and the same column. The coefficient of $abcd$ in $s$, and the coordinate deleted, are listed in the checker's `shifts` table. The remaining orbits use the unshifted full form. Both matrices split into a four-dimensional block on indices $(0,2,3,5)$ and a two-dimensional block on $(1,4)$, in the order $01,02,03,12,13,23$.

These two positive forms suffice on the complex Grassmannian. Choose the global phase of a decomposable $z=x+iy$ so that its deleted coordinate is real. Then $z^TSz=0$ implies $x^TSx=y^TSy$, and the deleted coordinate of $y$ vanishes. Applying the two forms gives $z^*M_hz\ge0$. The strict certificates give strict positivity for every nonzero $z$ when $a,b,c,d>0$.

### The polynomial certificates

There are 411 principal-minor checks across the 52 representatives. All but three polynomials have nonnegative monomial coefficients already. The three remaining polynomials are the order-four mixed and order-two pure determinants for orbit 47, and the unshifted order-four mixed determinant for orbit 51. The 15 terms in [squares.json](squares.json) subtract positive multiples of
$$
 (x^{u/2}-x^{v/2})^2=x^u-2x^{(u+v)/2}+x^v,\qquad x=(a,b,c,d).
$$
For example the exceptional pure determinant uses $4(x^{(0,2,1,1)}-x^{(2,0,1,1)})^2$. Each stored endpoint exponent is even. The verifier reconstructs the determinants over the integers, subtracts the exact squares, and checks every residual coefficient. Every strict test retains a positive monomial. No numerical optimization is needed to replay the proof. To bound the exponent packing, each matrix entry has total degree at most eight and no product in a principal minor has more than six factors. Each exponent is therefore below 64.

Permuting the three opposite blocks, exchanging the singleton labels, or exchanging the complementary eigenspaces preserves the family as the two projection triples range freely. The complete dual cone then supplies $m\ge0$ with the same spectral polynomial as $q$; the law normalizes to one. Degenerate triples follow by continuity even if their supports are disjoint. Approximate the two unit amplitude vectors in their separate three-dimensional coordinate spaces by strictly positive unit vectors. Near-identity rotations in those two spaces transport the singleton axes and the whole complement flag to the approximants. The finite family of dual inequalities is closed.

## Rotation of the singleton basis

For a plane with commuting projections in the above chart, common diagonal phases in all three coordinate blocks rotate the relative phase between the two axes without changing a physical score. Independent phases of the singleton columns also cancel. Thus every orientation, including the limit at infinity, can be represented by two orthogonal columns proportional to $X(1,t)^T$ and $X(-r_1t,r_0)^T$ for $0\le t<\infty$. Replace $X_{i0}X_{j1}$ in the partition formula by their ordered entries; use denominator $2r_0r_1(r_0+r_1t^2)g_z$. This follows from the product of their norms, $(r_0+r_1t^2)\sqrt{r_0r_1}$. The resulting matrix need not split into the previous four- and two-dimensional blocks.

For the specific two triples in result 2, [orientations.py](orientations.py) evaluates all six row-block permutations. The amplitudes are $\alpha=(1,1,3)$ and $\beta=(1,2,4)$ before normalization. For pivot $i$ and the remaining blocks $j,k$, the four chart values are $a=\alpha_j/\alpha_i,b=\beta_j/\beta_i,c=\alpha_k/\alpha_i,d=\beta_k/\beta_i$. Exchanging the singleton and complementary labels preserves the family of all orientations and flags. All 51 strict representatives use the full six-dimensional matrix and a shift $s=C abcd(r_0+r_1t^2)$ with $C$ from the aligned table. Orbit 3 is globally nonnegative by the [coarsening via Claude's reviewed Theorem D](../split-physical-compression/README.md).

Put $x=t^2$ and $R=r_0/r_1$. For a principal minor of order $k$ the executable checks directly that only even powers of $t$ occur and that $P(x)$ has degree $k$. On $0\le x\le R$ it checks every coefficient of $(1+y)^kP(Ry/(1+y))$ is strictly positive. On $R\le x$ it checks every coefficient of $P(R+y)$ is strictly positive. The leading coefficient covers the normalized limit at infinity. Each side chooses a permissible real coordinate for the phase lemma. This simultaneously certifies all angles and every complex complement flag; the calculation does not sample orientations.

Compactness gives the neighborhood assertion. On this compact orbit of planes and flags the 51 strict scores (and all their physical relabelings) have a strictly positive minimum. The remaining ten coefficient signs and orbit 3 are globally nonnegative. By continuity all signs persist in an open unitary neighborhood. The claim does not fix a shared parent across different orientations.

## Replay and review

```sh
python3 -O -I -B research/marcus-de-oliveira/commuting-split-projections/check.py
python3 -O -I -B research/marcus-de-oliveira/commuting-split-projections/orientations.py
```

Both verifiers use standard-library exact arithmetic. The polynomial searches that suggested the squares and phase choices are outside the proof. Independent mathematical review is requested, in particular for the changes of basis, paired complementary-minor signs and completeness transfer from the 62 orbits.
