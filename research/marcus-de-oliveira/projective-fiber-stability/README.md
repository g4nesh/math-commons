# Positive pairing weights near four spectral collisions

Take the same real orthogonal four-cycle matrix $U$ and $b=(3-2i,2,1+2i,-i)$ as the [three-free-eigenvalue theorem](../complex-direction-cover/README.md). Set $B=U D_b U^*$ and choose any index $j$. Write

$$A=\operatorname{diag}(x,-b_j+\eta,y,\delta).$$

For arbitrary $x,y\in\mathbb C$ and
$$\max(|\Re\eta|,|\Im\eta|,|\Re\delta|,|\Im\delta|)\le 1/1024,$$
one probability vector on all 24 pairings reconstructs $\det(A+B)$ simultaneously for **every** $x,y$. Every probability exceeds $1/4000$. The vector may depend on $j,\eta,\delta$. In each neighborhood 14 probabilities stay at $1/2000$ and the other ten are specified by a rational linear solve. The four centers in the second diagonal entry are $-3+2i,-2,-1-2i,i$. The fourth diagonal entry varies independently of it, while $U$ and $b$ stay fixed.

## Removing a rank failure at a collision

Use $h=x-\delta$, $k=y-\delta$, $z=-b_j+\eta$, $w=z-\delta$. For a column-to-row permutation $p$, write the pairing product minus the determinant as
$$D_p=w C_p+h L_p+k M_p+hk N_p.$$
The constant of the gap vanishes when all four entries of $A$ equal $\delta$. Holding $\delta$ fixed, $C_p$ is independent of $w$, and $L_p,M_p,N_p$ are affine in $w$. The usual equations for a probability law require normalization and both real components of the four weighted coefficients to vanish. Their nine-row real matrix has rank eight at each of these four centers with $\delta=0$. Directly inverting nine rows there would fail.

Require instead that the probability of assigning column $t$ to row four is $|U_{4t}|^2$. Three independent marginals and normalization suffice. When $x=y=z$, the actual determinant is
$$\prod_i(z+b_i)+(\delta-z)\sum_t |U_{4t}|^2\prod_{i\ne t}(z+b_i).$$
This follows by adding a rank-one diagonal perturbation to $zI+B$ and expanding the cofactor in its eigenbasis. Each pairing product has the same expression with its own row-four assignment indicator. The formula is polynomial, so it holds at a collision too. Consequently for any law with these marginals, evaluating its averaged gap at $h=k=w$ gives the polynomial identity
$$w\sum_p m_p(C_p+L_p+M_p+wN_p)=0.$$
For a fixed law and $\delta$ the identity holds for every $w$. Cancel the polynomial factor $w$ before evaluating. Requiring the weighted $C,L,N$ to vanish therefore forces the weighted $M$ to vanish, even at $w=0$. These give six real equations, ten in total with the probability and marginal equations. Each of the four certified bases is invertible at its center.

## Why the same bases work throughout the squares

Let $A_0$ be the ten selected center columns, $m^0$ the full positive center solution and $S$ the diagonal matrix of its ten adjustable probabilities. Hold the other probabilities fixed. Write $\Delta A$ for the change in all 24 columns and $\Delta A_B$ for the change in the selected ten. Then the relative changes $e$ of the adjustable probabilities solve
$$(I+S^{-1}A_0^{-1}\Delta A_B S)e=-S^{-1}A_0^{-1}\Delta A\,m^0.$$
The entries are real and imaginary parts of polynomials in $\eta,\delta$. For every nonconstant monomial $c\delta^d\eta^k$, its absolute value on the stated box is at most $(|\Re c|+|\Im c|)(2/1024)^{d+k}$. The verifier computes these bounds and propagates them through the absolute entries of the exact inverse. For each center it proves that the matrix infinity norm is bounded by $\rho<1$ and the right-hand side by $\nu$ with $\nu/(1-\rho)<1/2$. The Neumann series gives invertibility everywhere in the box and bounds every relative change by less than one half. Every center mass is at least $1/2000$. This proves strict positivity and the advertised common lower bound, including at irrational points and on the boundary.

## A neighborhood in the relative eigenbasis and opposing spectrum

The conclusion also has qualitative stability under small simultaneous changes of $U$ within $U(4)$ and of $b$ in $\mathbb C^4$, using the new centers $-b_j$. For an arbitrary complex unitary use $|U_{4t}|^2$ in the marginal equations. The ten equations vary continuously with all entries and their conjugates. For each center, the determinant of the chosen ten-column block has nonzero absolute minimum on the compact four-real-dimensional parameter box. Its inverse and the associated weights are continuous on an open neighborhood of that box. At the original fixture all 24 weights have a strict uniform margin above $1/4000$. Compactness and the finite number of centers preserve invertibility and that strict lower bound for all sufficiently small changes of $U,b$. The polynomial marginal argument above is valid for every such unitary. The checker gives the stated explicit radius for spectral motion at the archived fixture. This continuity argument supplies existence of a neighborhood in the additional parameters without assigning it a numerical radius.

## Reproduce the certificate

From the repository root run `python3 -B research/marcus-de-oliveira/projective-fiber-stability/check.py`. The checker constructs all coefficient polynomials from determinant minors, computes the four rank drops, reconstructs and bounds the four bases over the rationals, and checks off-center weights against the original dense determinant. The 10-column indices in [bases.json](bases.json) refer to the lexicographic order of column-to-row permutations of $(0,1,2,3)$. A linear program selected these indices. The replay uses Python's standard library and checks the indices directly.

This is a family-specific proof with public exact arithmetic. Its mathematical statement awaits independent-owner review. The result leaves open simultaneous motion of the fourth eigenvalue for every value of the middle eigenvalue.
