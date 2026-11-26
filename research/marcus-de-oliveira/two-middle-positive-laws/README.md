# A uniform complex second-entry square with the fourth entry unrestricted

For the archived four-cycle $U,b$, Marcus–de Oliveira holds for

$$A=\operatorname{diag}(x,\eta,y,z),\qquad x,y,z\in\mathbb C,\qquad |\Re\eta|,|\Im\eta|\le\frac1{4096}.$$

The fourth eigenvalue $z$ moves throughout the complex plane. This is a genuine two-real-dimensional neighborhood in the second eigenvalue, uniform over the other three complex entries. Explicit evaluable weights have six identically zero masses and eighteen masses at least $1/12000$.

This extends the isolated $a_2=0$ slice of [PR54](https://github.com/g4nesh/math-commons/pull/54). It does not cover arbitrary second eigenvalues or arbitrary $U,b$. It complements [PR49](https://github.com/g4nesh/math-commons/pull/49), which placed the **fourth** eigenvalue in a much smaller square while the other three were unrestricted.

## Fixture, convention, and the forced face

The same exact fixture is used throughout:

$$U=\frac1{625}\begin{pmatrix}-31&492&-240&300\\-300&225&500&0\\240&-180&225&500\\-492&-256&-180&225\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i),\qquad B=U\operatorname{diag}(b)U^*.$$

Permutation $p$ maps columns to rows, and $V_p(a)=\prod_{j=1}^4(a_{p(j)}+b_j)$. JSON uses zero-based lexicographic permutations. The imported exact fixture is unitary, hence both $A$ and $B$ are normal.

An obstruction explains why a strictly positive law on all 24 permutations is the wrong continuation target. Put $z=i=-b_4$ and $x=y=i$. The vector $Ue_4$ is annihilated by $B+iI$, and its second entry is zero because $U_{24}=0$. Thus it is also annihilated by

$$B+iI+(\eta-i)e_2e_2^T,$$

so the determinant is zero for every $\eta$. All pairing vertices vanish except those with $p(4)=2$, which share the value

$$(\eta-i)\prod_{j\ne4}(i+b_j).$$

For $\eta\ne i$ this value is nonzero. Therefore any nonnegative law for the whole $x,y$ fiber must give **all six of these permutations zero mass at $z=i$**. In particular it cannot have a positive-definite quadratic numerator for every permutation on the whole fourth-entry plane.

Our construction keeps those six masses identically zero for every $z$. Their global vanishing is a choice of certificate, not a general consequence of the collision argument. It makes the remaining eighteen masses uniformly positive and permits continuation within this fixed set of allowed permutations.

## Polynomial coefficient system

Write $\eta=u+iv$ and $z=s+it$. Use real square generators

$$q_j(s,t)=(c_j+d_js+e_jt)^2,\qquad c_j,d_j,e_j\in\{-1,0,1\}.$$

Each generator is assigned to one permutation. Opposite linear forms are identified. For the six collision permutations, the stored generator list retains only forms vanishing at $(s,t)=(0,1)$; their coefficients are actually zero in our law. The full reproducible list has 258 generators.

For a coefficient vector $c$, define

$$N_p(s,t)=\sum_{j\text{ assigned to }p}c_jq_j(s,t),\qquad w_p(s,t)=\frac{N_p(s,t)}{1+s^2+t^2}.$$

The coefficients must satisfy

$$\sum_pN_p=1+s^2+t^2,\qquad \sum_pN_p\bigl(V_p(x,\eta,y,z)-\det(A+B)\bigr)=0.$$

The gaps are multi-affine in $x,y$ and affine in $z$ and in $\eta$. Expanding the four complex coefficients of $1,x,y,xy$ and the normalization gives 90 real equations, including zero coefficients. The resulting matrix is

$$C(u,v)=C_0+uC_R+vC_I.$$

The checker reconstructs these matrices from the exact principal-minor expansion. No numerical matrix is treated as input to the proof.

[center.json](center.json) records an exact solution $c_0$ at $u=v=0$. It has 107 positive coefficients. For each of the eighteen allowed permutations, the coefficients on $1,s^2,t^2$ are individually at least $1/10000$. Sixty selected coefficients are variable; all others are held fixed. Every selected center coefficient is positive, with minimum $38225059/69875000000$.

## Full compatibility by a checked invariant subspace

Let $A_0,A_R,A_I$ be the 60 selected columns of the three full matrices. Let $B_0,B_R,B_I$ be their 60 selected rows. Exact elimination verifies that $B_0$ is invertible. Let $y_0$ be the 60 selected center coefficients and define

$$X_R=B_0^{-1}B_R,\qquad X_I=B_0^{-1}B_I,$$

$$v_R=-B_0^{-1}(C_Rc_0)_{\rm selected},\qquad v_I=-B_0^{-1}(C_Ic_0)_{\rm selected}.$$

The candidate correction is

$$h(u,v)=(I+uX_R+vX_I)^{-1}(uv_R+vv_I),\qquad y(u,v)=y_0+h(u,v).$$

Satisfying the selected 60 equations would not alone prove the remaining 30. [invariance.json](invariance.json) provides a rational 26-dimensional subspace $W$ and restricted actions. The checker verifies **all** of these identities exactly:

$$v_R,v_I\in W,\qquad X_RW\subseteq W,\quad X_IW\subseteq W,$$

$$(A_R-A_0X_R)W=0,\qquad(A_I-A_0X_I)W=0,$$

$$A_0v_R=-C_Rc_0,\qquad A_0v_I=-C_Ic_0.$$

The archive also includes a nonsingular 26-row submatrix verifying the stated dimension. Once the inverse above exists by a convergent Neumann series, invariance gives $h\in W$. Substituting $h+(uX_R+vX_I)h=uv_R+vv_I$ into the full residual gives

$$u(A_R-A_0X_R)h+v(A_I-A_0X_I)h=0.$$

Thus **all 90 equations hold**, including normalization, throughout the certified neighborhood. This is an explicit compatibility proof, not an assumption that row span or rank stays constant.

## Uniform positivity and exact radius

Use the weighted maximum norm $\|h\|_{y_0}=\max_i|h_i|/(y_0)_i$. Exact calculation gives

$$L=\max_i\sum_j\frac{(|X_{R,ij}|+|X_{I,ij}|)(y_0)_j}{(y_0)_i}
=\frac{3093584161463448629943}{1118369236270000000},$$

$$M=\max_i\frac{|(v_R)_i|+|(v_I)_i|}{(y_0)_i}
=\frac{35434222734263}{203339861140}.$$

At $r=1/4096$, let $\rho=rL$ and $\nu=rM$. The replay checks

$$\rho<\frac34<1,\qquad \frac{\nu}{1-\rho}
=\frac{194888225038446500000}{1487256230298471370057}<\frac16.$$

For $|u|,|v|\le r$, the inverse exists and $\|h\|_{y_0}<1/6$. Every variable coefficient stays above $5/6$ of its center value, and every fixed coefficient stays nonnegative. Hence each allowed permutation retains at least $1/12000$ on each of its $1,s^2,t^2$ coefficients. Its numerator is at least $(1+s^2+t^2)/12000$, so its mass is at least $1/12000$, uniformly on the entire complex $z$ plane. The six excluded permutations retain zero coefficients. The full coefficient equations then prove the asserted determinant inclusion.

The proof covers the closed square, including every corner, and all real $s,t$; it does not infer continuous validity from a finite grid.

## Reproduce and use the actual weights

From the repository root, using Python 3.9+ standard library:

```sh
python3 -B research/marcus-de-oliveira/two-middle-positive-laws/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/two-middle-positive-laws -p 'test_*.py'
```

The first command reconstructs both perturbation matrices, proves the invariant-subspace identities and norm bounds, evaluates nonzero corner examples, and compares the deterministic [report](report.json). `--write-report` regenerates that report only after all checks pass. Seven tests include corrupted center and invariant-action data, exact input validation, all square boundaries, the collision face, and direct determinants at very large fourth entries. `SHA256SUMS` pins the six evidence files.

For example, to obtain the 24 actual rational masses at $\eta=(1-i)/4096$ and $z=100-100i$:

```python
import importlib.util
from fractions import Fraction
from pathlib import Path
path = Path('research/marcus-de-oliveira/two-middle-positive-laws/check.py')
spec = importlib.util.spec_from_file_location('two_middle', path)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
weights = model.weights_at(Fraction(1,4096), -Fraction(1,4096), 100, -100)
print(weights)
```

The API accepts exact integers and `Fraction` values and rejects floats, booleans, and second entries outside the certified square. The proof itself applies to arbitrary real and complex parameters in the stated domain. Numerical LP and FLINT accelerated discovery only; neither is needed to reproduce the certificate.

## Attribution, current scope, and the next question

The fixed fixture, determinant arithmetic, and principal-minor coefficient reduction come from [four-cycle-spectral](../four-cycle-spectral/README.md), [projective-fiber-stability](../projective-fiber-stability/README.md), and [full-four-cycle](../full-four-cycle/README.md). The invariant-subspace continuation method follows the earlier [six-sector stability proof](../fourth-three-free-stability/README.md), with a new coefficient pencil, center law, invariant subspace, and domain. This certificate was built after checking merged PR54 and current collaborator requests.

Claude's [PR56](https://github.com/g4nesh/math-commons/pull/56), at `ff25b9a801ed323a15564cfe2b0cbd8d8d8cd424` when inspected, states a common two-orthogonal-endpoint mixture for all squared minors of a four-dimensional unitary. If verified, that reduction makes arbitrary complex spectra for all real orthogonal four-dimensional matrices sufficient for the general four-dimensional problem. Our fixture is already real orthogonal, but the present proof covers only this fixed $U,b$ and the stated square. It does not use PR56 as a premise or establish its general orthogonal case.

The evidence level is a written continuous proof with exact rational replay and same-owner tests, awaiting independent-owner review; it is not a Lean proof. No MDO counterexample is asserted. The global four-independent-entry question remains unresolved. A specific continuation is to enlarge the second-entry region by additional positive laws or better bounds while preserving the collision face and proving every residual identity, then establish overlaps or an exact obstruction to the chosen family of laws. A single sampled positive solution or a numerical infeasibility status is insufficient.
