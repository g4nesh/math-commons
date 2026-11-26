# A moving fourth eigenvalue on six angular sectors

Fix the archived four-cycle $U$, $b=(3-2i,2,1+2i,-i)$ and $B=U\operatorname{diag}(b)U^*$. This certificate proves a genuine complex perturbation of the fourth eigenvalue while $U$ and $b$ remain fixed.

Write

$$A=\operatorname{diag}(\delta+x,\delta+\omega(k)t,\delta+y,\delta),\qquad x,y\in\mathbb C,\quad t\in\mathbb R,$$

where $\omega(k)=1+ik$ in chart `h`, or $\omega(k)=k+i$ in chart `v`. **For both charts and every $k\in[-1,1/2]$, the original determinant inclusion holds whenever**

$$|\Re\delta|\le 1/4194304,\qquad |\Im\delta|\le 1/4194304.$$

The proof gives computable nonnegative permutation weights, uniform over the two arbitrary complex variables $x,y$. It includes the closed endpoints, zero and unbounded radial parameters, and coincident eigenvalues. Each patch has its own larger sufficient radius:

| Source index | Chart | Closed slopes | Fourth-eigenvalue square radius |
|---|---|---|---|
| 0 | h | $[-1,-1/2]$ | $1/8192$ |
| 1 | h | $[-1/2,0]$ | $1/16384$ |
| 2 | h | $[0,1/2]$ | $1/65536$ |
| 4 | v | $[-1,-1/2]$ | $1/32768$ |
| 5 | v | $[-1/2,0]$ | $1/32768$ |
| 6 | v | $[0,1/2]$ | $1/4194304$ |

This is a partial extension of the [full three-variable theorem at fourth entry zero](https://github.com/g4nesh/math-commons/blob/0638abd5464ce9122d96bf801675102029e22c82/research/marcus-de-oliveira/complex-direction-cover/README.md). The two quartic angular sectors $k\in(1/2,1]$ are **not certified by this artifact when $\delta\ne0$**. The coordinator addresses them through a separate projective-fiber and compactness argument; that argument is not a premise of this explicit six-sector certificate. Neither a failed invariant-kernel test there nor this omission is a counterexample to pointwise containment. The radii here are conservative sufficient bounds, not a sharp stability boundary.

[Exact bounds](certificate.json) · [Verifier and weights](check.py) · [Invariant-space derivation](neighborhood.py) · [Translated determinant coefficients](compatibility.py) · [Tests](test_stability.py) · [Replay output](verification.json)

## Translation preserves the correct problem

Put $A'=\operatorname{diag}(x,\omega(k)t,y,0)$ and $B'=B+\delta I$. Then $A+B=A'+B'$ and the eigenvalues of $B'$ are $b_j+\delta$. Pairing products transform by exactly the same equality. Thus calculating translated gaps with $b\mapsto b+\delta$ and fourth diagonal entry zero proves a statement about the **original fixed $B$** and the moving fourth entry of $A$.

At $x=y=t=0$, $A=\delta I$ and every pairing product equals the determinant. The translated gap therefore has zero scalar term for every $\delta$. The implementation checks this cancellation. It avoids introducing a spurious constant-coefficient constraint at the scalar matrix.

For fixed diagonal variables, both the determinant and pairing products are monic degree-four polynomials in $\delta$. Their degree-three coefficients are the same trace. Hence their difference is holomorphic of degree at most two in $\delta$. With $\delta=u+iv$, all real coefficient systems have the exact form

$$C(u,v)=C_0+uC_u+vC_v+u^2C_{uu}+uvC_{uv}+v^2C_{vv}.$$

The code constructs these matrices from exact determinant and pairing polynomials. The degree bound justifies the coefficient extraction; an additional rational-point regression checks it away from the extraction points.

## Why full column-space stability would be wrong

Each archived quadratic patch supplies 64 strictly positive term coefficients $w_0$. Its 92 real equations include every mixed coefficient of the angular/radial gap polynomial and the two Bernstein-control normalization rows. Let their right side be $b_0$. The archived identity is $C_0w_0=b_0$.

Choose 64 independent rows $I$. Exact linear solves define

$$K_j=(C_0)_I^{-1}(C_j)_I,\qquad R_j=C_0K_j-C_j.$$

The matrices $R_j$ are **not zero**: the full perturbation columns do not remain in the base column space. A continuity argument that assumed this would leave a real compatibility gap. Direct first- and second-order tests on the archived weight vector passed, but those finite Taylor tests alone would also be insufficient.

Instead the replay constructs a matrix $H$ with 18 independent rows and proves, over the rationals,

$$\operatorname{row}(R_j)\subseteq\operatorname{row}(H),\qquad Hw_0=0,\qquad HK_j=L_jH$$

for all five perturbation matrices and suitable exact $L_j$. Each identity is checked in every entry. The selected pivot columns of $H$ form the identity matrix, so its rank is exactly 18. Consequently its kernel $W$ has dimension 46, contains $w_0$, and is invariant under every $K_j$. Modular arithmetic chooses candidate pivots only; exact solves and complete residual checks establish every assertion.

Let $K(u,v)=uK_u+vK_v+u^2K_{uu}+uvK_{uv}+v^2K_{vv}$. Whenever the inverse below exists, set

$$w(u,v)=(I+K(u,v))^{-1}w_0.$$

The inverse preserves $W$ by its convergent Neumann series. Thus $R_jw(u,v)=0$ for every $j$, and

$$C(u,v)w(u,v)=C_0(I+K(u,v))w(u,v)-\sum_jm_j(u,v)R_jw(u,v)=C_0w_0=b_0.$$

This proves exact compatibility of the full moving system on a neighborhood without assuming that its entire column space is constant.

## Explicit positivity and the infinite radial range

Set $S=\operatorname{diag}(w_0)$; every diagonal entry is positive. The replay computes exact infinity norms

$$M_j=\|S^{-1}K_jS\|_\infty,\qquad A_j=\|S^{-1}K_jw_0\|_\infty.$$

For a closed square of radius $r$, put

$$\rho=r(M_u+M_v)+r^2(M_{uu}+M_{uv}+M_{vv}),$$
$$\nu=r(A_u+A_v)+r^2(A_{uu}+A_{uv}+A_{vv}).$$

At each recorded radius, the exact rational tests give $\rho<1$ and $\nu/(1-\rho)\le1/2$. The matrix inverse exists and

$$\|S^{-1}w(u,v)-\mathbf1\|_\infty\le\nu/(1-\rho)\le1/2.$$

Every deformed coefficient is therefore at least half its archived positive value. All normalization equations remain exact. The support of each Bernstein control retains precisely the same radial square factors as before; their root sets have empty intersection. Thus their sum remains positive for every real $t$. At least one Bernstein factor is positive at every slope, including both endpoints. Dividing the resulting nonnegative permutation numerators by their positive sum gives the required convex weights.

The proof covers irrational complex $\delta,x,y$ as well as rational inputs. Executable rational examples are regression evidence; the invariant-subspace identities and norm inequalities prove the continuous region.

## Reproduction and next step

From the repository root, Python 3.9 or later and the standard library suffice:

```sh
python3 -B research/marcus-de-oliveira/fourth-three-free-stability/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/fourth-three-free-stability -p 'test_*.py'
```

`coefficients_at(index, real, imag)` solves the deformed patch and checks all 92 rows and the positive lower bound. `weights_at(index, real, imag, slope, radius)` returns the 24 exact permutation masses. All parameters must be integers or `Fraction` values. The output order is lexicographic column-to-row permutations. The source patch archive is SHA-256-bound in the certificate; its original 878-term theorem is cited above.

`neighborhood.py --index N` derives one patch; `--all` derives the six quadratic patches. Both are read-only unless an explicit `--report PATH` is supplied. Exact linear solves use fraction-free integer elimination, integer back substitution and full rational residual checks. No optimizer or floating-point linear algebra is used in replay.

The next request is to produce one explicit numerical fourth-eigenvalue radius across every complex middle-eigenvalue direction. Combine the separate qualitative compactness argument with quantitative collision neighborhoods, uniform rank bounds on the remaining compact projective region, and the explicit six-sector bounds here. Its deliverable should include a complete finite cover, rational lower bounds for the chosen minors and probability margins, an exact replay and stated dependence on $U,b$. The current quartic supports have noninvariant initial residual kernels; that sufficient test alone neither disproves their continuation nor pointwise containment. The present six-sector proof remains valid independently of the all-direction extension.

This is same-owner mathematical and executable evidence, awaiting independent-owner review. It does not resolve the general Marcus–de Oliveira conjecture.
