# Orbit A is valid on every unitary, with an exact contact classification

For every $U\in U(5)$, using the zero-based row-to-column convention of
[Ishaan's PR115](https://github.com/g4nesh/math-commons/pull/115),

$$
\boxed{|\det U[12,12]|^2\le |U_{00}|^2+
|\det U[03,34]|^2+|\det U[04,34]|^2.}
$$

This proves the **whole 900-facet orbit A** of the minor polytope $Q_5$ on
physical unitary matrices. PR115 established its validity on the 120
permutation vertices; the argument below supplies the physical inequality,
a sum-of-squares gap and its complete equality locus. It is a self-contained
linear-algebra proof, not an inference from optimization or the fixtures.
No claim of literature novelty is made.

This concerns $Q_5$, not the fixed-spectrum principal-minor polytope $P_5(c)$
studied in [Claude's PR114](https://github.com/g4nesh/math-commons/pull/114).
It does not prove PMSH(5), classify the 62 non-inherited $P_5(c)$ facets, or
decide which eight-permutation sets are realizable rigid spectral contacts.
The remaining six non-minor orbits from PR115 are untouched.

## Proof by a three-by-three compression

Compress to rows and columns $(0,3,4)$ and write

$$X=\begin{pmatrix}a&b&c\\d&e&f\\g&h&i\end{pmatrix},\quad
p=bf-ce,\quad q=bi-ch,\quad \Delta=ei-fh.$$

Let $y=(a,p,q)$ and $v=(\Delta,g,-d)$. Laplace expansion gives
$\det X=a\Delta+gp-dq=y\cdot v$. The dot here is **bilinear**;
conjugation in the following identity is therefore essential:

$$
\|y\|^2-|\det X|^2
=(1-\|v\|^2)\|y\|^2+
\sum_{0\le j<k\le2}|y_j\overline{v_k}-y_k\overline{v_j}|^2. \tag{1}
$$

The bottom two rows of $X$ form a contraction. Write their first column as
$z=(d,g)^T$ and their remaining block as
$Y=\left[\begin{smallmatrix}e&f\\h&i\end{smallmatrix}\right]$.
Then $YY^*\le I_2-zz^*$, so

$$|\Delta|^2=\det(YY^*)\le\det(I_2-zz^*)=1-|d|^2-|g|^2.$$

For completeness, determinant monotonicity here needs only a two-by-two
identity: for positive semidefinite $A,E$,
$\det(A+E)-\det A=\operatorname{tr}(\operatorname{adj}(A)E)+\det E\ge0$.
Thus $1-\|v\|^2\ge0$, and (1) proves the inequality for **every three-by-three
contraction**, even with no unitary completion specified. For $U\in U(5)$,
complementary unitary minors give $|\det X|^2=|\det U[12,12]|^2$, proving the
displayed claim with exactly PR115's sign and coordinate conventions.

## Explicit unitary sum of squares

There is no need to leave the deficit implicit. Put

$$
\tau=\sum_{J\in\{12,13,14,23,24\}}|\det U[34,J]|^2.
$$

The two rows indexed by 3 and 4 are orthonormal. Cauchy–Binet gives the sum
of all ten squared two-minors as 1. The sum of the four minors using column
0 is $|U_{30}|^2+|U_{40}|^2$: expand each square and use the two row norms
and their orthogonality. Removing those four minors and the minor on 34
therefore proves

$$\tau=1-|d|^2-|g|^2-|\Delta|^2=1-\|v\|^2.$$

Hence the physical facet gap is exactly

$$\boxed{\quad G(U)=\tau\bigl(|a|^2+|p|^2+|q|^2\bigr)
+\sum_{j<k}|y_j\overline{v_k}-y_k\overline{v_j}|^2.\quad} \tag{2}
$$

Every summand is nonnegative. The product term is itself a sum of 15 squared
products, so (2) is an explicit sum of 18 squared polynomial expressions in
the entries and their conjugates, restricted to the unitary identities.

## All equality mechanisms

Equality holds **if and only if at least one** of the following conditions
holds. The alternatives can overlap; no genericity assumption is imposed.

1. **Zero right side:** $U_{00}=p=q=0$.
2. **A zero two-by-two block:** $U[34,12]=0$.
3. **A supported column:** $U_{00}=U_{10}=U_{20}=0$, so column 0 is entirely
   supported on rows 3 and 4.

Here is a proof of both directions, including singular cases. Define
$W=U[34,12]$, $Y=U[34,34]$, and $Z=[W\;Y]$.
The five squares in $\tau$ are $\det W$ and the four mixed minors selecting
one column from each of $W,Y$. If $Y$ is invertible, their vanishing forces
both columns of $W$ to be zero. If $Y$ has rank at most one, all six
two-minors of $Z$ vanish, so $\operatorname{rank}Z\le1$. Conversely, either
$W=0$ or $\operatorname{rank}Z\le1$ makes $\tau=0$.

For the two orthonormal rows, $ZZ^*=I_2-zz^*$. Thus
$\operatorname{rank}Z\le1$ is equivalent to $\|z\|^2=1$. The full unitary
column norm makes this precisely condition 3. We have proved

$$\tau=0\quad\Longleftrightarrow\quad\text{condition 2 or condition 3}. \tag{3}$$

It remains to ensure the three cross terms in (2) vanish in both branches;
vanishing of $\tau$ alone would not suffice for arbitrary matrices.

Under condition 2, the bottom two rows of $X$ are orthonormal in
$\mathbb C^3$, and its top row is orthogonal to them. Therefore
$X=\operatorname{diag}(\lambda,1,1)V$ for some $V\in U(3)$ and
$\lambda\in\mathbb C$ (including $\lambda=0$). The complementary-minor
identities of $V$ give $y=\lambda\det(V)\overline v$. All cross terms vanish.

Under condition 3, $a=0$, $\Delta=0$, and $|d|^2+|g|^2=1$.
Orthogonality of column 0 to every other column gives
$\overline d\,U_{3j}+\overline g\,U_{4j}=0$ for $j\ne0$.
Taking the corresponding two-minors with row 0 gives
$\overline d\,p+\overline g\,q=0$. Consequently
$(p,q)=\lambda(\overline g,-\overline d)$ for some $\lambda$, so again
$y=\lambda\overline v$ and the cross terms vanish.

Condition 1 makes (2) zero immediately. Conversely, if the gap is zero and
$y\ne0$, its first nonnegative term forces $\tau=0$, and (3) gives condition
2 or 3. This completes the classification.

## Exact verification and irreducible contacts

From the repository root, Python 3.9+ standard library only:

```sh
python3 -B research/marcus-de-oliveira/n5-contact-program/strata/physical-a/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/n5-contact-program/strata/physical-a -p 'test_*.py' -v
(cd research/marcus-de-oliveira/n5-contact-program/strata/physical-a && shasum -a 256 -c SHA256SUMS)
```

The replay checks all 120 permutation signs: 84 zeros and 36 ones. It
independently reconstructs the 125 minor coordinates, verifies affine ranks
41 for $Q_5$ and 40 for the equality face, generates all 900 symmetry copies,
and checks that PR115's archived orbit-A representative belongs to that
orbit. Row and column permutations and transposition preserve unitarity, so
the written proof transports to every copy. There is no claim that the nine
facet orbits currently archived in PR115 are a complete list.

Four Gaussian-rational unitary fixtures are built explicitly from small
complex Householder vectors, permutations and rational two-row gates in
`fixtures()`. They separately realize the three equality mechanisms and a
strict case:

| Fixture | Gap | $\tau$ | Right side |
| --- | ---: | ---: | ---: |
| Zero block only | $0$ | $0$ | $9/16$ |
| Supported column only | $0$ | $0$ | $10/121$ |
| Zero right side, neither geometric branch | $0$ | $75/121$ | $0$ |
| Strict dense fixture | $17008/75625$ | $18528/75625$ | $18717/75625$ |

Every fixture has a connected row-overlap graph: join rows sharing a
nonzero column. Hence none is reducible to a direct sum by row and column
permutations. In particular, physical contact of this inequality need not
mean a block-diagonal unitary. This does not identify these contacts with
the irreducible maximizers of Claude's different fixed-spectrum facets.

The executable fixtures check conventions and exact identities; the
universal theorem and equality classification are the written arguments
above. Same-owner replay is not independent-owner review.

## Handoff

This is a scoped continuation of the released `mdo-n5-contact-strata`
request. The old census, conjecture and artifacts are preserved. The broad
question of which rigid eight-contact patterns occur remains open. A useful
next mathematical step is to prove or delimit one of PR115's other covering
inequalities on physical unitaries, or to determine which realizable gap
cones are generated by squared-minor rays together with this now-proved
orbit. Merely lying on an orbit-A facet does not prove an entire contact
cone is settled.
