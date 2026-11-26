# One common parent law for every split-row orientation

This extends the [exact split-pair criterion](../split-one-eigenvalue-pair/README.md)
from PR99. On its dense complex Gaussian-rational fixture $U_0$, **one fixed
21-table parent law lifts to a positive 33-table law for every orientation of
the last row plane**. Every child mass can be at least $1/20000$. The same
statement, with floor $1/40000$, holds for every unitary within entrywise complex
distance $1/53760000$ of $U_0$, uniformly over all orientations.

The spectral pattern is

$$
A=\operatorname{diag}(x_0,x_0,x_1,x_1,\eta+\delta,\eta-\delta),\qquad
B=U_R\operatorname{diag}(b_0,b_0,b_1,b_1,b_2,b_2)U_R^*,
$$

where $U_R=\operatorname{diag}(I_4,R)U_0$ and $R\in U(2)$ is arbitrary.
All seven displayed spectral parameters are arbitrary complex numbers. The
weights depend on the orientation and the coarse unitary, **not** those
spectral parameters. This is a full continuous orientation family and an open
unitary subclass, not a proof for all coarse unitaries or all $6\times6$ spectra.

## Exact replay and evaluation

From the repository root, using Python 3.9+ and only its standard library:

```sh
python3 -B research/marcus-de-oliveira/split-row-orientations/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/split-row-orientations -p 'test_*.py' -v
(cd research/marcus-de-oliveira/split-row-orientations && shasum -a 256 -c SHA256SUMS)
```

The checker regenerates the prior 33 tables, 192 monomial coefficient rows,
six even kernel moves, six odd equations and all 31 unoriented facet normals.
It expands all 720 grouped Leibniz terms symbolically in the three Pauli
coordinates, checks every rational certificate entry, and verifies the sphere
support inequalities after shrinking capacities to enforce the child floor.
It then replays three genuinely complex rational $SU(2)$ rotations, constructs
their weights and checks all 192 coefficient identities. The tests also cover
the three axes, the south pole, a nonzero coarse-plane perturbation, tampering,
and rejection of floating-point or non-unit input.

Import `check.py` and call `weights_at((Fraction(-16,45), Fraction(172,225),
Fraction(121,225)))` for an exact 33-entry probability vector. The ordering is
the lexicographic table order inherited from PR99's `certificate.json`.
`nearby_weights(U,R)` accepts exact Gaussian-rational matrices in the stated
unitary neighborhood and returns the corrected law. The theorem includes
irrational orientations and matrices; the executable evaluator intentionally
accepts exact rational inputs only.

## Pauli response, including its gauge

Write $R^*\sigma_zR=n_x\sigma_x+n_y\sigma_y+n_z\sigma_z$, where
$|n|=1$ and $\sigma_y=\left[\begin{smallmatrix}0&-i\\i&0\end{smallmatrix}\right]$.
Thus $n=(\Re H_{01},-\Im H_{01},H_{00})$ for $H=R^*\sigma_zR$.
The overall phase of $R$ cancels from each grouped determinant ratio, so choose
$R=\left[\begin{smallmatrix}u&v\\-\bar v&\bar u\end{smallmatrix}\right]\in SU(2)$.
For entries $a_j=(U_0)_{4j}$ and $b_j=(U_0)_{5j}$, the last two factors of a
Leibniz term with columns $j,k$ are exactly

$$
\frac{a_jb_k-b_ja_k}{2}
+\frac{n_x(-a_ja_k+b_jb_k)}2
-\frac{i n_y(a_ja_k+b_jb_k)}2
+\frac{n_z(a_jb_k+b_ja_k)}2.
$$

Multiply by the first four row factors and permutation sign, divide by
$\det U_0$, group by the 33 tables, and take real parts **after complex
division**. This gives a real affine response $q(n)=c+En$ for all canonical
table masses. The checker derives this entire $33\times4$ coefficient array;
the three test rotations are additional checks, not an interpolation premise.

Merging the two singleton row colors gives 21 parents, nine with one child and
twelve with two. Each parent sum is independent of $n$: the split-row sum is
an exterior two-row determinant, invariant under $SU(2)$. Each double-child
difference changes sign under exchange of the split rows, so its constant
term vanishes. Consequently

$$d^{\rm can}(n)=Dn,\qquad t(n)=BDn=Tn,$$

with exact matrices $D\in\mathbb Q^{12\times3}$ and
$T\in\mathbb Q^{6\times3}$ stored in `certificate.json`. The six-row matrix
$B$ and the parent order are exactly PR99's conventions. For a rational unit
$n$ with $n_z\ne-1$, the unnormalized rational frame

$$
\widetilde R(n)=\begin{pmatrix}1+n_z&n_x-in_y\\-n_x-in_y&1+n_z\end{pmatrix}
$$

has $\widetilde R^*\widetilde R=2(1+n_z)I$. Its grouped determinant ratios
equal those of the normalized unitary. The south pole uses
$\left[\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right]$. This removes
square roots from the exact frame checks without pretending the frame itself
is unitary.

## Necessary and sufficient common-parent test

Let $q^0$ be the fixed canonical 21-parent vector and let $K$ be PR99's six
even kernel moves. Fix a repair $r\in\mathbb R^6$ and put
$m=q^0+K^Tr\ge0$. Denote the parent of double fiber $j$ by $p(j)$.
A lift at orientation $n$ exists exactly when

$$Bd=Tn,\qquad |d_j|\le m_{p(j)}.$$

The complete 31-normal description of this zonotope was derived and replayed
in PR99. Thus **this one parent law works for every orientation if and only if**

$$
\|T^Th\|_2\le \sum_{j=1}^{12}m_{p(j)}|h\cdot B_j|
\quad\text{for every one of the 31 normals }h.
$$

Indeed, maximize $|h\cdot Tn|$ on the unit sphere; its exact maximum is
$\|T^Th\|_2$. This exchanges two universal quantifiers, not the existential
choice of $d$, which may depend on $n$. Nonnegative support values permit
squaring both sides, leaving rational arithmetic only. Existence of a common
parent is a stronger requirement than allowing a different parent at every
orientation. Failure of the stronger requirement would not be an MDO
counterexample.

Our six rational repairs give strictly positive parent masses. With
$f=1/20000$, replace every double capacity by $m_{p(j)}-2f$ and require each
single mass at least $f$. All 31 support inequalities still pass, with strictly
positive squared slack. The minimum reduced squared slack is

$$
\frac{105334129914171950851804901039811769}
{4427849400851811351275024414062500000000}>0.
$$

The child masses $(m_{p(j)}\pm d_j)/2$ therefore satisfy the stated floor.
Normalization and all spectral polynomial identities follow from the even
kernel and six odd equations, as in PR99.

For evaluation, the feasible box section $Bd=Tn$ is compact and nonempty.
A vertex has at most six coordinates not at bounds; their columns are
independent and can be extended to a six-column basis since $\operatorname{rank}B=6$.
Enumerating six free columns and the signs of the remaining six bounds must
therefore find a solution. `imbalances_at` checks at most
$\binom{12}{6}2^6=59136$ candidates, all exactly. The three archived mixed
fixtures require 326, 538 and 1153 candidates respectively. No solver is
needed to replay or evaluate the certificate.

## A uniform open neighborhood of the whole orbit

PR99 proves that each canonical table coefficient changes by at most
$672\epsilon$ when two $6\times6$ unitaries differ entrywise by at most
$\epsilon$. Left multiplication by the same two-row unitary amplifies the
entrywise difference by at most $2$ (in fact $\sqrt2$ suffices). Hence

$$|q_j(U_R)-q_j((U_0)_R)|\le1344\epsilon.$$

For any baseline orientation law $w_R$, define

$$w_{U,R}=w_R+q(U_R)-q((U_0)_R).$$

The correction preserves normalization and all 192 polynomial identities.
If $\epsilon\le1/53760000$, each resulting mass is at least $1/40000$.
Moreover its parent vector is $m+q^0(U)-q^0(U_0)$, independent of $R$.
Thus a common parent still exists for every nearby coarse unitary; the radius
is uniform over the entire compact orientation sphere. This bound is
conservative and is not asserted optimal.

## Provenance and remaining question

The fixture, table algebra, facet enumeration and Lipschitz bound come from
[PR99](https://github.com/g4nesh/math-commons/pull/99), immutable source
`f4bf89382a6b49e210a41e93134bd82756ae8f12`. The new work is the algebraic
orientation response, ellipsoid containment criterion, common-parent exact
certificate, uniform orbit neighborhood, and exact child-weight evaluator.
A floating linear program suggested the six repairs; every published claim
is checked independently of that solver using exact arithmetic. An optional
numerical attempt to force a single linear imbalance map $d(n)=Ln$ did not
find a positive solution. This establishes no impossibility and is not used.

The dispatchable remaining question is whether every coarse unitary admits
a common parent satisfying the 31 norm inequalities. If not, separate that
failure from feasibility with orientation-dependent parents, and provide an
exact certificate for the claimed obstruction. The all-unitary theorem and
general single-pair-split MDO problem remain open here. These same-owner
checks are submitted for independent-owner review, not labeled as such a
review or a new universal theorem.
