# Five triple correlations beyond the colored pairs

Let $U\in U(6)$ have its columns grouped into three pairs $J_0,J_1,J_2$. A successful spectrum-independent certificate chooses a random partition of the six rows into three labeled pairs. For every row subset $S$, its color counts must have the determinantal histogram $h_U(S)$ obtained by choosing columns $T$ with mass $|\det U[S,T]|^2$. This supplies the convex combination in Marcus–de Oliveira for every complex opposing spectrum and every choice of the three block eigenvalues. [Polynomial identity and matching boundary](../paired-matching-boundary/README.md)

**Reduction.** It is enough to reproduce the 45 individual colored-edge probabilities and five real triple correlations. The edge statistics have linear rank 38. All five additional correlations are needed to distinguish arbitrary signed or positive laws with those same edge data. For cyclically invariant data, we also obtain a necessary-and-sufficient test involving only 15 signed matching variables and elementary inequalities.

## 1. Recovering all the histograms

For an unordered edge $e$ and color $j$ let $D_j(e)=|\det U[e,J_j]|^2$. For a triple $S$ set
$$A_S=\mathbb E_{h_U(S)}\,\alpha(N_0,N_1,N_2),\qquad
\alpha(n)=\tfrac12(n_0-n_1)(n_1-n_2)(n_0-n_2).$$
Here $\alpha$ is zero for a rainbow, and is $+1$ or $-1$ on the six permutations of $(2,1,0)$. An equivalent way to compute the target uses $\zeta=e^{2\pi i/3}$ and the unitary $H=U\operatorname{diag}(I_2,\zeta I_2,\zeta^2 I_2)U^*$. Then
$$A_S=\frac{2}{\sqrt3}\Im\det H[S,S].$$
Indeed Cauchy–Binet gives the expectation of $\zeta^{N_1+2N_2}$, whose imaginary part has the stated sign on each histogram.

**Recovering the counts from edge data.** For color $j$ the random two rows carrying that color must have law $D_j$. In the target data the chance of both hits in $S$ is $\sum_{e\subset S}D_j(e)$ by the determinantal inclusion identity, and the chance of neither is $\sum_{e\subset S^c}D_j(e)$ by complementary minors. The mean count follows from the endpoint marginals of $D_j$. These quantities have the same interpretation in any candidate coloring law.

On one row the three means recover the histogram. On two rows the three double probabilities are known. Subtract twice those probabilities from the means. What remains gives the three pairwise sums of the three mixed-color probabilities, which determine them individually. On three rows each non-rainbow outcome has a unique doubled color $i$ and missing color $j\ne i$. Put these six probabilities in an off-diagonal $3\times3$ table. Its row sums are the three double probabilities and its column sums are the three missing probabilities. The remaining probability is the rainbow. The off-diagonal positions form one six-cycle in the bipartite row-column graph. Tables with zero row and column sums on that support are scalar multiples of the alternating signs around the cycle. The scalar $A_S$ distinguishes exactly this remaining degree of freedom. Larger row sets follow by complementation, valid for both candidate colorings and unitary minors.

**Why five triples suffice.** For any fixed coloring, its 20 triple signs $a_S=\alpha(N(S))$ satisfy
$$\sum_{S\supset e,\ |S|=3}a_S=0\qquad (|e|=2).$$
If the two rows of $e$ have the same color, the remaining four have two copies of each other color, giving opposite signs. If they differ, adding another copy of the first or second color gives opposite signs and the other two choices are rainbows. The same relations hold after averaging. They also hold for the target because every unitary admits the canonical signed coloring law by the determinant expansion. [Canonical signed law](../paired-matching-boundary/README.md#4-keeping-the-labels-when-symmetry-breaks)

The 15 incidence equations are independent. If edge coefficients $b_{ij}$ sum to zero on each triple, subtract the equations on $(i,j,k)$ and $(i,j,l)$. The four differences $b_{ik}-b_{il}$, for $i$ outside $\{k,l\}$, have all pairwise sums zero, so all are zero. Hence the edges at each vertex have a common coefficient, and symmetry and a triple equation force that coefficient to vanish. The space $\mathcal H$ of admissible sign vectors has dimension $20-15=5$. With zero-based rows the triples $012,013,014,023,024$ give coordinates on it. Appending those five coordinate rows to the incidence matrix gives rank 20 in the [integer replay](check.py).

For clarity, the matrix of 45 colored-edge indicators on the 90 possible colorings has rank 38. Its column span consists of the tables for which each color has the same total mass $t$ and the combined edge degrees at every vertex equal $t$. Counting gives dimension at most 38: the nine constraints in the 46 coordinates $(D,t)$ have one redundancy, and $t$ is determined by $D$. Exact elimination gives the matching lower bound. Together with the five triples the rank is 43, the full rank of the 358 grouped histogram rows. [Earlier exchange-kernel proof](../three-color-exchanges/README.md)

## 2. How to change only the five correlations

Write $M$ for one of the 15 uncolored perfect matchings. Give its three lexicographically ordered edges the labels $0,1,2$ and define $B_{S,M}$ to be the resulting triple sign. Relabeling by $\sigma$ multiplies it by $\mathrm{sgn}(\sigma)$. Changing the mass of that labeling by $\mathrm{sgn}(\sigma)t_M$ preserves normalization and every individual colored-edge marginal, even conditional on the matching. It changes the triple sign vector by $6B_{:,M}t_M$.

The columns form a tight frame on $\mathcal H$:
$$BB^\mathsf T d=36d\qquad(d\in\mathcal H).$$
Here is a direct count. $G=BB^\mathsf T$ depends only on $|S\cap T|$ because renaming sites can change the orientation of a matching only by a common sign in the two factors. $G_{S,S}=9$: choose the internal edge on $S$ and the internal edge on $S^c$ in three ways each, then join the leftover vertices. Complementation reverses a sign. Applying the pair-incidence relation to a fixed pair inside $S$ shows that $G_{S,T}=-3$ for $|S\cap T|=2$. The four values in increasing order of intersection size are therefore $-9,3,-3,9$.

For $d\in\mathcal H$, all single-site sums and the total sum of $d$ are zero by adding the pair equations. Inclusion–exclusion gives $d_{S^c}=-d_S$. Summing the constraints on the three pairs inside $S$ gives $\sum_{|S\cap T|=2}d_T=-3d_S$. The sum for intersection size one is $3d_S$ by complementation. Substituting these four sums into $Gd$ proves the identity.

Consequently any discrepancy $\Delta$ in the five triple correlations can be corrected algebraically by taking $t=B^\mathsf T\Delta/216$. If a base law $u$ already has the correct colored-edge marginals, the explicit adjustment is a valid probability law whenever
$$-\min_{\sigma\ \mathrm{even}}u_{M,\sigma}\ \le\ t_M\ \le\
\min_{\sigma\ \mathrm{odd}}u_{M,\sigma}\qquad\text{for every }M.$$
One may add any vector in $\ker B$ to $t$. For fixed per-matching colored-edge tables these intervals and $6Bt=\Delta$ are **necessary and sufficient**: each such table records the three-by-three edge-to-label assignments, a sum of six permutation matrices. Their unique linear relation is the sum of the even permutations minus the sum of the odd permutations. This follows directly since fixing an edge and its color pairs one even labeling with one odd labeling, and the nine such pairs form a connected bipartite graph. Allowing the tables to vary can give additional solutions. Around the strictly positive uniform base on 90 colorings one may perturb in any of these five directions while preserving all 45 edge probabilities. Thus edge marginals alone cannot determine any of the five remaining degrees of freedom.

## 3. Exact criterion under cyclic color symmetry

Suppose all target histograms are invariant under $0\to1\to2\to0$. Write $x_e=\sum_j D_j(e)$, so $D_j(e)=x_e/3$. There is a full coloring certificate **if and only if** there exist 15 reals $\delta_M$ satisfying
$$ B\delta=A,\qquad t=1-\sum_M|\delta_M|\ge0,\qquad
y_e=x_e-\sum_{M\ni e}|\delta_M|\ge0,\qquad
\sum_{e\subset S} y_e\le t\quad(|S|=3).$$

Indeed average a successful law over the three cyclic color rotations. For each matching the three even labelings now have a common mass, and the three odd labelings have a common mass. Let $m_M$ be the total and $\delta_M$ the mass of all even minus all odd labelings. Then $m_M\ge|\delta_M|$ and the triple vector is $B\delta$. Subtracting $|\delta_M|$ from each $m_M$ gives a nonnegative residual matching measure, proving necessity. Conversely, the displayed conditions and the degree equations supply a residual matching measure by Edmonds' theorem, also when $t=0$. Add back $|\delta_M|$ to get $m_M$. Give each labeling mass $(m_M+\mathrm{sgn}(\sigma)\delta_M)/6$. This recovers $D_j(e)=x_e/3$ and $A$, hence every target histogram by Section 1. [Matching inequalities](../paired-matching-boundary/README.md#1-forgetting-the-labels-works-in-every-even-dimension)

For a general unitary, finding one compatible law of colored edges and triple correlations remains the nonnegativity problem. The frame gives a signed correction in closed form. Its sign constraints carry genuine information.

Run `python3 -I -B research/marcus-de-oliveira/five-alternating-triples/check.py`. It checks the rank statements, proves sufficiency of the five named coordinates by rational elimination, verifies the frame and sign conventions on all 90 colorings, and reconstructs every grouped row of a cyclic law with unequal even and odd masses.
