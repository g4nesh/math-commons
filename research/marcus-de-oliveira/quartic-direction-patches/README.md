# Two exact quartic radial patches

Fix the archived four-cycle unitary and noncollinear spectrum

$$U=\frac1{625}\begin{pmatrix}-31&492&-240&300\\-300&225&500&0\\240&-180&225&500\\-492&-256&-180&225\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i),\qquad B=U\operatorname{diag}(b)U^*.$$

We construct explicit convex pairing weights for

$$A=\operatorname{diag}(x,z,y,0),\qquad x,y\in\mathbb C,$$

on the two separate closed angular intervals

$$z=(1+ik)t,\qquad k\in[1/2,9/16]\;\cup\;[15/16,1],\quad t\in\mathbb R.$$

For each fixed $z$, these weights work simultaneously for every complex $x,y$. They vary with the real radial coordinate $t$ and the angular coordinate $k$. In particular, this is a constructive family result for the original determinant inclusion. It does not supply one distribution valid for all three free complex eigenvalues, and it makes no claim for arbitrary $U$ or $b$.

These are independent Cramer-basis checks of selected regions in the coordinator's separate global angular construction. The standalone claim here is exactly the two intervals above, including their endpoints. This artifact does not cover the gap between them or the second angular chart. Indices are one-based: $a_1=x,a_2=z,a_3=y,a_4=0$. The different [two-free-eigenvalue result](https://github.com/g4nesh/math-commons/blob/cb253d9c3c7746fd2278db9d41a5eb1668ed52b4/research/marcus-de-oliveira/two-free-rank-one/README.md) fixes $a_3=2-i/4$ and uses a quadratic identity to prove pointwise inclusion without common two-variable weights. Its conclusion is not a premise here.

[Certificate](certificate.json) · [Exact replay](patches.py) · [Replay output](verification.json) · [Tests](test_patches.py) · [Exact reconstruction](derive.py)

## Explicit nonnegative weights

The certificate contains two rational patches, with centers $17/32$ and $1$. Each has 47 positive terms, degree-12 angular numerators and a degree-13 common denominator. On a patch with center $k_0$, a term consists of a column-to-row permutation $p$, zero, one or two integer centers chosen from $\{-2,-1,0,1,2\}$, and a polynomial $n(k-k_0)$ with rational coefficients. Its radial factor is respectively

$$1,\qquad (t-r)^2,\qquad (t-r)^2(t-s)^2.$$

Thus every radial factor $f(t)$ is nonnegative for every real $t$, including negative and arbitrarily large values. Every angular numerator is strictly positive throughout its assigned closed interval, as certified by its Bernstein coefficients.

Collect the terms attached to permutation $p$ as

$$N_p(k,t)=\sum_{j:p_j=p}n_j(k-k_0)f_j(t),\qquad Z(k,t)=\sum_pN_p(k,t).$$

The desired weights are the explicit real rational functions

$$w_p(k,t)=\frac{N_p(k,t)}{Z(k,t)}.$$

They are nonnegative and sum to one. The verifier establishes a strictly positive lower bound for $Z$ on each patch and for every real $t$. A positive constant term already suffices; otherwise it retains two positive square terms with distinct roots. For their uniform angular lower bounds $A,B>0$,

$$A(t-r)^2+B(t-s)^2\ge\frac{AB(r-s)^2}{A+B}>0.$$

The recorded bound may be conservative. Patch boundaries are included. If two patches meet, either formula is valid at the common endpoint; continuity between different valid formulas is unnecessary. At $z=0$ any certified chart with $t=0$ works.

## Why the weights match determinants for all free eigenvalues

For each permutation define the pairing gap

$$D_p(x,y,z)=\prod_{j=1}^4(a_{p(j)}+b_j)-\det(A+B).$$

It is affine in each of $x,y,z$ separately. In particular, the four evaluations $(x,y)\in\{0,1\}^2$ determine the entire bi-affine polynomial in the two arbitrary complex variables. In a given angular chart, each such gap is affine in $t$ and in $k$.

Multiplying by a degree-at-most-four radial factor gives degree at most five in $t$. The exact replay verifies

$$\sum_p N_p(k,t)D_p(x,y,z)=0$$

as a polynomial identity. It checks the real and imaginary parts of all six radial coefficients at the four Boolean pairs: 48 real coefficient rows. Each row is itself checked coefficient by coefficient in $k-k_0$. Multi-affinity then extends the equality to every complex $x,y$. Dividing by the strictly positive $Z$ gives the required determinant convex combination. No numerical tolerance, finite sample in $x,y$, or finite sample in $t$ enters this proof.

The same argument includes repeated eigenvalues, zero determinants and degenerate pairing hulls. Unitarity of the fixed $U$ is checked by the shared exact archived arithmetic, so $B$ is normal as required.

## How the angular certificates were obtained

At a rational center $k_0$, the initial scratch search considered 21 radial primitives and 24 permutations. A numerical linear program selects a 47-column support. It then solves all coefficients over the rational numbers, rejecting any support with a negative mass or an exact residual.

Let $C(k)=C_0+(k-k_0)C_1$ denote the resulting 49-by-47 system, including the 48 gap rows and a final normalization row. At each accepted support the script verifies an exact identity

$$C_1=C_0K.$$

Thus a center solution $v_0$ extends through $(I+(k-k_0)K)^{-1}v_0$. An exact Krylov relation produces a common denominator $q$ and numerator vector $n$ satisfying

$$(I+\delta K)n(\delta)=q(\delta)v_0,\qquad\delta=k-k_0.$$

The export stores the resulting rational polynomial coefficients. The standalone verifier does **not** trust the numerical rank, support discovery, matrix inverse or Krylov construction: it reconstructs all gap columns directly and checks every coefficient of $C(k)n(k-k_0)=q(k-k_0)e_{49}$ exactly. The final row also checks $\sum_jn_j=q$.

For positivity, the script converts $q$ and every $n_j$ to the Bernstein basis on the entire closed angular patch. Strict positivity of every coefficient proves strict positivity of the polynomial there. In the exploratory search, failed sign checks led to interval subdivision and a fresh support at the new midpoint. The replay reports the two covered intervals explicitly; it makes no gap-free angular-cover assertion.

## Reproduction and evidence level

Python 3.9 or later, standard library only, suffices for the archived proof and tests:

```sh
python3 -B research/marcus-de-oliveira/quartic-direction-patches/patches.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/quartic-direction-patches -p 'test_*.py'
```

`weights_at(patch, k, t)` returns the 24 exact masses in lexicographic column-to-row permutation order. Pass integers or `Fraction` values, and choose a patch containing $k$. The finite examples in the tests check the formula against direct physical determinants at both endpoints and off-grid eigenvalues; the polynomial and positivity identities prove the continuous regions.

`derive.py` reconstructs both archived continuations from their selected support alone, using only rational linear algebra and modular pivots that are followed by exact residual checks. Run it with `python3 -B research/marcus-de-oliveira/quartic-direction-patches/derive.py`. It checks that every regenerated coefficient matches the archive; no numerical solver is required. The original support discovery used SciPy, random objective seed 37, and the five centers above. Discovery is not part of the proof.

The broader adaptive support search reached its finite subdivision limit near $k=0.654$. This was a limitation of that particular support and Bernstein-sign workflow, not an infeasibility certificate. At the coordinator's request, the task then narrowed to these compact independent patches after a different global Bernstein construction succeeded. The two explicit certificates and their continuous proofs are the completed deliverable.

These artifacts are same-owner mathematical and computational evidence, submitted for independent-owner review. They are not a general proof of the Marcus–de Oliveira conjecture.
