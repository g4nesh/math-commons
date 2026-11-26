# A split plane with arbitrary projection proportions

Consider the six-dimensional Marcus–de Oliveira problem with four spectral values of multiplicities $(2,2,1,1)$ against three of multiplicities $(2,2,2)$. All seven eigenvalues can vary independently in $\mathbb C$. Write $E$ for the plane spanned by the two simple eigenvectors and $B_0,B_1,B_2$ for the opposite two-dimensional spectral blocks. Suppose there exist $p_j\ge0$ summing to one such that
$$
P_EP_{B_j}|_E=p_j I_E\quad(j=0,1,2).
$$
This says that every unit vector in $E$ has the same squared projection $p_j$ onto $B_j$. In geometric terms the plane is isoclinic to each block.

**Theorem.** For every such triple, every ordered orthonormal basis of $E$, and every orthogonal splitting of $E^\perp$ into the two remaining spectral pairs, there are 33 nonnegative table weights summing to one which give the determinant simultaneously for every choice of the spectral values. The relative eigenbasis determines the weights. This extends the [balanced-plane theorem](../split-physical-compression/README.md), where all three proportions are $1/3$. No positive lower bound on a proportion is required.

## Physical form on an entire flag

First take all $p_j>0$. In an arbitrary chosen ordered singleton basis, the maps from $E$ to the three blocks are scalar multiples of isometries. Independent basis changes inside the opposite blocks put the singleton columns in the form $X/\sqrt r$ with
$$
 X=[I_2;aI_2;bI_2],\quad a=\sqrt{p_1/p_0},\quad b=\sqrt{p_2/p_0},\quad r=1+a^2+b^2.
$$
This gauge makes no restriction on the singleton orientation. With $W=[aI_2;bI_2]$, set $Y=[-W^T;I_4]$ and $G=Y^TY$. Then $X^TX=rI$, $\det G=r^2$, and $\det[X,Y]=r^2$. The complementary planes are parametrized by the arbitrary full-rank four-by-two frame $Z$ of the first one in these $Y$ coordinates. Write $z=\wedge^2 Z$ for its six minors and $g_z=z^*G^{[2]}z>0$. The ordered physical row table in the original convention assigns paired groups $K,L$ and singleton groups $i,j$. Use the Hodge matrix $S$, order of coordinates, and row-pair notation of the [balanced proof](../split-physical-compression/README.md). Put
$$
 L_h=\sum \operatorname{sgn}(i,j,K,L)h(K,L,i,j)X_{i0}X_{j1}(Y^{[2]}_{K,:})^T(Y^{[2]}S G^{[2]})_{L,:},\qquad M_h=L_h+L_h^T.
$$
The same complementary-minor calculation yields the physical ray value $h\cdot q=z^*M_h z/(2 r^2 g_z)$. The two singleton columns contribute $1/r$, and the final normalized pair contributes $1/(r\sqrt{g_z})$ in addition to the first pair's $1/\sqrt{g_z}$. In this note $M_h$ omits the balanced proof's extra factor of three. All entries are integer polynomials of total degree at most eight in $a,b$. The [exact checker](check.py) also evaluates the original 720-term determinant on two unequal positive triples and one boundary case with complex frames, matching all 33 coefficients in each case.

## Polynomial proof of all remaining signs

The [exact classification](../split-dual-rays/README.md) has 935 rays in 62 physical orbits. Ten orbits are grouped squared-minor coefficients and are nonnegative for any unitary. The calculation here deals directly with every one of the other 52 representatives. In particular, it does not require the separate 17-orbit coarsening theorem. [Claude's Theorem D](../minors/theorem-d.md) and the [coarsening proof](https://github.com/g4nesh/math-commons/pull/108) give an alternative sign argument for some of them.

Use the order $01,02,03,12,13,23$. In the following table, $c$ is a scalar and $k$ is an index in this zero-based ordering. For each row take $s=c\,a^2b^2$ and choose a global phase of $z$ making $z_k$ real.

| Orbit | $c$ | $k$ |
| ---: | ---: | ---: |
| 8 | -2 | 0 |
| 10 | 2 | 0 |
| 11 | 2 | 0 |
| 14 | -2 | 5 |
| 16 | -2 | 2 |
| 17 | -4 | 2 |
| 29 | -4 | 2 |
| 36 | -2 | 2 |
| 46 | 4 | 0 |
| 47 | -4 | 5 |
| 49 | 4 | 0 |

The matrices have zero interaction between indices $(0,2,3,5)$ and $(1,4)$, as do their Hodge shifts. For the eleven listed orbits the checker proves $M_h+sS$ is positive definite and $(M_h-sS)|_{e_k^\perp}$ is positive definite. It checks the leading principal minors in each block. The two-by-two block on $(1,4)$ has the same leading minor and determinant for both signs because its unshifted off-diagonal entry vanishes; the implementation checks this structural zero. For 40 other orbits it proves positive definiteness of $M_h$ by the same block test without a shift. For orbit 3 it proves positive semidefiniteness by checking all 63 nonempty principal minors, including identically zero polynomials.

Each nonzero required minor has only positive integer coefficients and even powers of both parameters. The strict tests also require the polynomial to contain at least one monomial. The 402 polynomial checks are recorded per orbit with SHA256 fingerprints of every coefficient and exponent. Our checker reconstructs the matrix entries as exact integer polynomials from the table scores, before inspecting or hashing any minors. Six base-64 exponent bits per variable suffice because total degree never exceeds 48. A separate FLINT calculation gave the same sign decisions.

Write $z=x+iy$ in the chosen phase. The Plücker equation implies $x^TSx=y^TSy$ and $y_k=0$. Thus $z^*M_hz=x^T(M_h+sS)x+y^T(M_h-sS)y\ge0$. The positive denominator now gives the sign of the original determinant functional. Physical relabelings permute the proportions and the two complementary planes; because **all** proportions and both complementary planes were allowed, each representative controls its entire orbit. The classification's finite-dimensional duality then gives the 33 nonnegative weights.

For boundary proportions, choose an opposite block of positive weight as block 0 and keep the gauge above. Approximate vanishing $a$ or $b$ by positive values with the chosen chart $Z$ fixed. The matrices, metric and positive denominator are continuous, and $G$ stays positive definite. At a block with $p_j=0$ either basis is permissible since the projected singleton vectors vanish. This proves the closed probability simplex including its vertices.

For a fixed triple of strictly positive proportions, the 51 strict scores have a positive common minimum over all unit singleton frames and complementary flags, by compactness and the finite number of labels. Orbit 3 is globally nonnegative by the coarsening through Claude's reviewed Theorem D in the [balanced proof](../split-physical-compression/README.md). The other ten global inequalities come from squared minors. Consequently a unitary neighborhood of the whole fiber for each strictly positive triple also satisfies the spectral law. The explicit radius proved at the balanced triple remains available there.

## Reproduce

```sh
python3 -O -I -B research/marcus-de-oliveira/isoclinic-split-planes/check.py
```

The standard-library replay takes under a second on the development machine. Its [polynomial manifest](polynomials.json) can be regenerated to stdout with `--print-certificate`. Independent mathematical review of the physical reduction and finite sign proof is pending. The open distinction is a split plane for which the squared projection to an opposite block varies with the vector in the plane.
