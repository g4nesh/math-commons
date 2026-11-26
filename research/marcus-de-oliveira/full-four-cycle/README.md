# An entire complex fourth-eigenvalue slice

For the archived four-cycle, there are explicit convex weights proving Marcus–de Oliveira for

$$A=\operatorname{diag}(x,0,y,z),\qquad x,y,z\in\mathbb C,$$

with $B=U\operatorname{diag}(b)U^*$ fixed. The weights depend on $z$ and work simultaneously for every $x,y$. This permits **unbounded complex motion of the fourth eigenvalue**, while the second eigenvalue stays zero. It does not prove the conjecture for four independently varying complex entries.

The fixture and pairing convention are

$$U=\frac1{625}\begin{pmatrix}-31&492&-240&300\\-300&225&500&0\\240&-180&225&500\\-492&-256&-180&225\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i),$$

$$V_p(a)=\prod_{j=1}^4(a_{p(j)}+b_j).$$

Each permutation maps **columns to rows**. JSON uses zero-based permutations in lexicographic order. The checker verifies $U^*U=UU^*=I$ through the shared exact arithmetic; therefore both matrices in the statement are normal.

## Explicit weights and proof

Write $z=s+it$, with $s,t$ real. [certificate.json](certificate.json) gives 58 terms $(p,c,v_0,v_1,v_2)$, all $c$ strictly positive rational and all $v_i\in\{-1,0,1\}$. Set

$$N_p(s,t)=\sum_{\text{terms assigned to }p}c(v_0+v_1s+v_2t)^2,\qquad w_p(s,t)=\frac{N_p(s,t)}{1+s^2+t^2}.$$

These are actual executable weights, not an existence oracle. The smallest square coefficient is $2895468812/4308762109375$. There are 15 supported permutations; nine masses vanish identically. Positivity of the **square coefficients** gives nonnegative masses, not a claim that every mass is strictly positive.

The exact checker reconstructs all principal-minor coefficients of

$$g_p(x,y,s,t)=V_p(x,0,y,s+it)-\det(\operatorname{diag}(x,0,y,s+it)+B).$$

This polynomial is multi-affine in the complex variables $x,y$ and affine in the real variables $s,t$. It then proves the polynomial identities

$$\sum_pN_p=1+s^2+t^2,\qquad \sum_pN_pg_p=0.$$

Specifically, it checks the four complex coefficient polynomials belonging to $1,x,y,xy$, each of total real degree at most three in $s,t$, and the normalization polynomial: 90 real coefficient equations including zero coefficients. Thus the identities hold for **all** parameters, not just sampled values. The denominator is everywhere positive. Dividing gives nonnegative masses summing to one and the required convex combination. No rank-continuation or unproved row-span stability assumption enters this proof.

For example, run the following from the repository root to obtain all 24 exact masses at $z=7/3-11i/5$:

```python
import importlib.util
from fractions import Fraction
from pathlib import Path
path = Path('research/marcus-de-oliveira/full-four-cycle/check.py')
spec = importlib.util.spec_from_file_location('four_cycle', path)
four_cycle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(four_cycle)
print(four_cycle.weights_at(Fraction(7, 3), Fraction(-11, 5)))
```

The exact API rejects floats and booleans. Integer or `Fraction` inputs are required. Rational input restrictions concern the evaluator; the polynomial proof applies to all real $s,t$ and complex $x,y$.

## Why preserving all singleton marginals cannot solve the general fiber

Here is a structural obstruction to a particular strategy, not to pointwise determinant inclusion. It applies to any $4\times4$ unitary before specializing the final contradiction to our fixture.

Suppose a real probability law $w$ on permutations has all singleton marginals

$$\sum_{p:p(j)=i}w_p=|U_{ij}|^2.$$

Let $D_b(a)=\sum_pw_pV_p(a)-\det(\operatorname{diag}(a)+U\operatorname{diag}(b)U^*)$. Cauchy–Binet shows that its degree-zero and degree-four parts cancel by normalization, and its degree-three part cancels by the singleton marginals. The complementary-minor identity for a unitary,

$$|\det U_{i^c,j^c}|^2=|U_{ij}|^2,$$

cancels its degree-one part. (This identity follows from the adjugate formula and $U^{-1}=U^*$.) Hence $D_b$ is a homogeneous multiaffine quadratic. The same argument applies to every shifted spectrum $b+c\mathbf1$, while

$$D_b(a+c\mathbf1)=D_{b+c\mathbf1}(a).$$

The right side is homogeneous quadratic in $a$. Therefore the linear and constant terms produced by translating the left side vanish. If $d_{ij}$ are its six pair coefficients, this says $\sum_{j\ne i}d_{ij}=0$ for each $i$. Solving these four linear equations gives precisely

$$D_b(a)=q(a_1-a_3)(a_2-a_4)+r(a_1-a_4)(a_2-a_3).$$

Now fix $a_2=z$, $a_4=\delta$ and ask for one such singleton-preserving law working for every $a_1=x,a_3=y$. The $xy$ coefficient is $-r$, so $r=0$. The $x$ coefficient is then $q(z-\delta)$. When $z\ne\delta$, this forces $q=0$ too. That same law would consequently work for **all four** complex entries.

For the archived $U,b$, [the existing exact seven-diagonal separator](../four-cycle-spectral/fiber_certificate.json) excludes any one common law on all four entries: every permutation function has functional score at most zero, whereas the determinant function has score $11201152/15625>0$. The present checker replays that certificate without a solver.

**Consequently, for every $z\ne\delta$, any common law for the two-free fiber must change at least one singleton marginal.** Our constructive plane law does so. The earlier fixed-law obstruction is credited to the four-cycle spectral certificate; this consequence explains why a marginal-preserving extension strategy cannot close the remaining gap.

The hypothesis $z\ne\delta$ matters. The additional 24 rational masses in `repeatedMiddleWeights` preserve all singleton marginals and have the exact discrepancy

$$D_b(a)=\frac{-4919527+5276689i}{12500000}(a_1-a_3)(a_2-a_4).$$

They therefore certify $a=(x,d,y,d)$ for every complex $x,y,d$, but are not a fixed law on all four independent entries. This is an explicit sharpness example, not a claim of novelty for repeated-spectrum subclasses; see the earlier [minor framework](../minors/report.md).

## Relation to the current disk and collaborators' work

[PR45](https://github.com/g4nesh/math-commons/pull/45) proved qualitative uniform fourth-entry stability of the previous three-free result. [PR46](https://github.com/g4nesh/math-commons/pull/46) supplied explicit radii for six sectors. [PR49](https://github.com/g4nesh/math-commons/pull/49) supplied a uniform square $|\Re a_4|,|\Im a_4|\le10^{-100}$ with **all three other entries arbitrary**. The present result removes the restriction on $a_4$ completely **on the slice $a_2=0$**. These are complementary statements; the whole-plane slice does not enlarge that uniform all-three-free radius.

The archived [complex-direction cover](../complex-direction-cover/README.md) treats $a_4=0$ with the other three entries free. The present orientation is different and is explicitly checked on the same $U,b$. The affine coefficient construction reuses the principal-minor reduction in [projective-fiber-stability](../projective-fiber-stability/check.py), and the exact arithmetic and fixture come from [four-cycle-spectral](../four-cycle-spectral/spectral.py). Claude's current [PR48](https://github.com/g4nesh/math-commons/pull/48) concerns a unitary-minor inequality and real-symmetric/Hermitian spectral hypotheses. Our fixed $b$ is noncollinear, so that guarantee is not used as a premise. Main `b02605f6e8c18e4b2fac1d329ccd4971dcd8d4d4` and open PR48/51 were checked before construction.

A bounded exploratory floating LP scan of 2,169 rational middle fibers found no negative feasibility margin. It is not a certificate and is not used in either theorem. Numerical support discovery also suggested the 58-term construction; all published coefficients were then solved and checked with exact rational arithmetic. There is no runtime solver dependency.

## Reproduce, limitations, and next question

From the repository root, using Python 3.9+ standard library:

```sh
python3 -B research/marcus-de-oliveira/full-four-cycle/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/full-four-cycle -p 'test_*.py'
```

The first command recomputes and compares [report.json](report.json). Tests include coefficient corruption, invalid exact inputs, large complex fourth entries, the repeated-middle identity, and direct determinant computations independent of polynomial expansion. `--write-report` intentionally regenerates the report after all checks. `SHA256SUMS` pins this directory's proof and executable artifacts.

The evidence level is a written proof with exact polynomial replay and same-owner tests; it is not a Lean proof or an independent-owner review. Four independently varying complex diagonal entries remain unresolved. Neither the singleton obstruction nor the archived fixed-law separator is a pointwise MDO counterexample.

A concrete next task is to construct a parameter-dependent law for $a=(x,z,y,\delta)$ with $z$ nonzero while $\delta$ remains unbounded, or give an exact obstruction to a specified polynomial degree/support. The present weights contain identically zero masses, so an automatic strict-positivity perturbation argument does not establish a uniform tube around $a_2=0$.
