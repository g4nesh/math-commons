# Three double eigenvalues on both sides

Assume both spectra have multiplicities $(2,2,2)$ in dimension six. This note reduces the stronger, spectrum-independent determinant certificate for a fixed unitary to **four universal inequality types**. Each type has 36 row/column relabelings. Three additional types hold for every unitary by Cauchy–Binet. We give an exact algorithm that produces the 21 nonnegative weights whenever the inequalities pass. All four signs are now proved for every complex $U(6)$ by the [Grassmann compression and exact polynomial certificates](../double-paired-five-term/GRASSMANN.md). That written proof awaits independent-owner mathematical review. The reduction and solver below are unchanged.

## The 21 possible products and their complete kernel

Index $S_3$ by $0,\ldots,5$ in this order: $012,120,201,021,102,210$. Write $s_i$ for the parity of the $i$th permutation and $p_i=\prod_{r=0}^2(a_r+b_{\pi_i(r)})$. Each six-row permutation gives a $3\times3$ nonnegative integer table with row and column sums two. Such a table decomposes into two permutation matrices: the corresponding bipartite multigraph is a union of even cycles and doubled edges. Alternate the edges along each component. There can be at most one cycle of length four or more on six vertices, so the **unordered** pair of permutations is unique. Conversely every such pair occurs. The 21 vertices to average are $p_i p_j$ for $i\le j$.

The matrix $(a_r+b_s)_{r,s}$ has rank at most two, hence $L=\sum_i s_i p_i=0$. Multiplying by each $p_i$ gives six dependencies between the 21 products. Their quadratic coefficient vectors are independent since each contains its own diagonal term. They span the entire kernel: expand the 21 products in the monomials with six individual exponents at most two and total degree six. The resulting integer matrix has rank 15 by the exact elimination in the replay. This also checks that 141 monomials occur. These identities concern the whole polynomial in the six eigenvalues. Weights for a single numerical spectrum may have more freedom.

Obtain signed starting coefficients $q_{ij}$ by grouping any canonical signed coloring law according to its three row-block counts. Equivalently, group the real parts of the Leibniz coefficients $\mathrm{sgn}(\sigma)\overline{\det U}\prod_k U_{k,\sigma(k)}$ by the same tables. Their sum is one. [Signed determinant expansion](../n5-contact-program/contact-reduction/README.md)

All alternative coefficient vectors arise by adding $\sum_i s_i u_i p_iL$ to the identity, for reals $u_i$. Thus a common law exists exactly when
$$w_{ii}=q_{ii}+u_i\ge0,\qquad
w_{ij}=q_{ij}+s_i s_j(u_i+u_j)\ge0\quad(i<j).$$
Their sum stays one since $\sum s_i=0$. This provides a small constructive test. Introduce two vertices per index, with potentials $v_i^+,v_i^-$. Turn $u_i+u_j\le c$ into the two difference constraints $v_i^+-v_j^-\le c,\ v_j^+-v_i^-\le c$; reverse both signs for a bound on $-u_i-u_j$. The diagonal constraint is $v_i^--v_i^+\le 2q_{ii}$. A feasible system gives $u_i=(v_i^+-v_i^-)/2$, by adding the two reflected constraints. A feasible $u$ gives $v_i^\pm=\pm u_i$. Shortest-path relaxation on these 12 vertices therefore gives exact weights or detects a negative cycle.

## All dual obstructions

By finite-dimensional linear-programming duality, the above inequalities are feasible exactly when $\sum_{i\le j}y_{ij}q_{ij}\ge0$ for every nonnegative symmetric array $y$ with $\sum_j s_j y_{ij}=0$ at each $i$. The pairing uses each unordered pair once. Multiplying the balance equation by $s_i$ gives a six-vertex graph: within-parity edges and loops have one color, opposite-parity edges have the other. At every vertex the total of each color must be the same, with a loop counted once.

Lift the graph to the 12 states $(i,c)$ recording the next color. An edge of color $c$ from $i$ to $j$ gives an arc $(i,c)\to(j,1-c)$, in both directions when $i\ne j$. A simple directed cycle projects to a feasible vector by counting each edge traversal once and each loop traversal twice. Conversely a feasible $y$ lifts to a circulation by putting $y_{ij}$ on each of the two arcs of a nonloop edge and $y_{ii}$ on its single loop arc. Flow balance is exactly the original equation. The projection gives $2y$. Decomposing the circulation into simple directed cycles proves that their projections generate every obstruction.

For a nonnegative homogeneous null cone a feasible ray is extreme exactly when its support is inclusion-minimal. A proper sub-support permits subtraction of a maximal multiple of its feasible vector. For a minimal support any two positive vectors on it must be proportional by the same subtraction argument. Consequently enumerating simple cycles, dividing projected vectors by their greatest common divisors, and retaining inclusion-minimal supports is exhaustive. The exact replay gives 573 distinct primitive projections and 189 extreme rays.

Independent row and column relabelings of the three values, together with transposing the table (exchange the two normal matrices), induce precisely the 72 permutations of the six indices preserving the unordered parity partition. The replay constructs the group directly from the maps $\pi\mapsto\rho\pi\gamma^{-1}$ and $\pi\mapsto\pi^{-1}$. Transposition in the determinant problem exchanges $U$ with $U^*$. We use these symmetries to count universal inequalities. For a particular fixed $U$, all 189 rays must be tested unless its extra symmetry has been proved.

## The seven representatives

In the table $ij$ means $q_{ij}$, and each displayed sum must be nonnegative. Evaluate every row/column relabeling of its representative.

| Number of rays | Representative | Universal status |
| ---: | --- | --- |
| 9 | $22+25+55$ | A grouped squared-minor coefficient |
| 18 | $2(22)+24+25+45$ | A grouped squared-minor coefficient |
| 18 | $12+15+24+45$ | A grouped squared-minor coefficient |
| 36 | $12+15+24+44+55$ | Grassmann polynomial certificate |
| 36 | $2(12)+2(15)+23+24+34+2(55)$ | Grassmann polynomial certificate |
| 36 | $02+05+11+14+23+34+55$ | Grassmann polynomial certificate |
| 36 | $02+05+12+14+2(23)+2(33)+45$ | Grassmann polynomial certificate |

To verify the first three rows, compute the vector of coefficients contributed by the 21 pairing polynomials to each monomial. Its dot product with $q$ is that monomial's coefficient in the physical determinant, a sum of squared minors of $U$ by Cauchy–Binet. After division by their positive integer content, the 141 rows give 71 distinct vectors. Exactly 45 of them are extreme rays, consisting of the first three entire orbits.

The four non-Cauchy–Binet types are independent obstacles in the abstract data. For each of their representatives the checker constructs a normalized signed coefficient vector that satisfies every other extreme inequality and every monomial-positivity test, but violates the chosen inequality. It does this by placing initial mass on all coordinates missing from the chosen support and making one coordinate of the support slightly negative. These artificial vectors do not establish that any unitary violates an inequality. The Grassmann proof supplies the additional properties of unitary minors needed along this route.

## Exact checks and handoff

Run `python3 -I -B research/marcus-de-oliveira/double-paired-contingencies/check.py` using Python 3 and its standard library. The verifier rebuilds the polynomial matrix, kernel, graph, symmetry group, extreme supports and seven orbits. It tests proof independence on abstract signed data and reconstructs nonnegative pairing laws on complex unitary fixtures, checking the polynomial coefficients directly against all 924 squared minors. [Exact classifier and constructive solver](check.py)

The five-term representative $q_{12}+q_{15}+q_{24}+q_{44}+q_{55}$ and the other three signs are covered by the [universal two-sided proof](../double-paired-five-term/GRASSMANN.md). Its separate `grassmann.py` and `grassmann_fixtures.py` replays verify the 24 leading-minor polynomial certificates and the complex-frame conventions. Review of that proof is pending. The case of only one $(2,2,2)$ spectrum opposite an unrestricted spectrum remains open; this reduction concerns both spectra paired.
