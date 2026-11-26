# Five-term facet: exact reductions and a failed Hodge relaxation

Task `mdo-double-paired-five-term`, agent `ganesh-reductions`, owner `g4nesh`.
This is the result of a bounded initial attack. **The universal inequality is
still open. No negative unitary example was found or certified.**

The target from [PR81](https://github.com/g4nesh/math-commons/pull/81), pinned
at `272f7ced73c3b880cff0f44adfe25f6f8cceca7f`, is

`y(U)=q_12+q_15+q_24+q_44+q_55 >= 0`, for arbitrary `U in U(6)`.

Rows and columns are partitioned into consecutive pairs `R_r,C_s`, `r,s=0,1,2`.
The six permutations are ordered `012,120,201,021,102,210`. The canonical
`q_ij` groups `Re[conjugate(det U)*sign(sigma)*product_k U[k,sigma(k)]]`
by the contingency table `P_i+P_j`. The contingency reduction and polynomial
kernel belong to the coordinator's source; this note treats only this facet.

We obtain two compact exact formulas and show that an apparently natural
relaxation cannot prove their sign. The obstruction is not a counterexample
to the five-term inequality or to MDO.

## 1. Third-compound formula

For triples of nonnegative integers `alpha,beta` summing to three, define

`H_(alpha,beta)(U) = sum |det U[S,T]|^2`,

where `S,T` range over three-element subsets with row-pair counts `alpha` and
column-pair counts `beta`. Then

```
4 y(U) = -H_(102,111) + H_(111,120)
         +2 H_(012,210) +2 H_(021,201).
```

Digits here are three separate block counts, not integers. This involves
40 squared minors with multiplicities: 16 negative, 16 positive of coefficient
one, and eight positive of coefficient two. Complementary unitary minors give
`H_(102,111)=H_(120,111)`, so the first two terms are a difference between
incoming and outgoing squared block norms between the two-dimensional weight
sector `120` and the eight-dimensional sector `111` in `wedge^3 C^6`.

**Proof of the identity.** Expand each of the 21 products `p_i p_j` in the six
formal spectral variables. The coefficient of `a^(2-alpha) b^beta` in the
physical determinant is `H_(alpha,beta)` by Cauchy–Binet. The integer linear
combination of the four coefficient rows displayed above equals the vector
which is four on the five target coordinates and zero elsewhere. `check.py`
rebuilds these polynomials from the permutations and checks all 21 coordinates
exactly. This is a polynomial identity, not an inference from matrix samples.

## 2. A six-block trace formula

Write the off-diagonal 2-by-2 blocks of `U` as

```
[ *  A  B ]
[ C  *  D ]
[ E  F  * ],
```

and let `b_rs=det U[R_r,C_s]`, `z=conjugate(det U)`. Then

```
y(U) = |b_20|^2 + |b_10|^2
       + Re z [ tr(adj(B) A adj(F) E adj(C) D)
                - det(A)det(D)det(E) - det(B)det(C)det(F) ].
```

The bracket is multiplied by `z` before taking the real part. The order of
these generally noncommuting 2-by-2 factors matters.

**Derivation.** The permutations numbered 1 and 5 share the cell `(2,0)`.
Grouping the 2+4 Laplace expansion on that full row/column pair and then using
Jacobi's complementary-minor identity gives

`q_11+q_15+q_55=|b_20|^2`.

Likewise permutations 2 and 4 share `(1,0)`, giving
`q_22+q_24+q_44=|b_10|^2`. Therefore
`y=|b_20|^2+|b_10|^2+q_12-q_11-q_22`.
This reduction was supplied by the read-only `protocol_docs` collaborator;
its index conventions are checked on exact physical fixtures here.
The diagonal coefficients are
`q_11=Re[z det(A)det(D)det(E)]` and
`q_22=Re[z det(B)det(C)det(F)]`.

In the matrix obtained by setting the three diagonal blocks to zero, only the
three contingency tables `P_1+P_1`, `P_1+P_2`, `P_2+P_2` remain. A block Schur
complement, first for invertible `A,C`, gives

```
det [0 A B; C 0 D; E F 0]
 = det(A)det(D)det(E) + det(B)det(C)det(F)
   + tr(adj(B) A adj(F) E adj(C) D).
```

For example, interchange the first two block columns (sign positive since
the blocks have size two), then take the Schur complement of `diag(A,C)`.
Use `det(X+Y)=det X+det Y+tr(adj(X)Y)` for 2-by-2 matrices. The resulting
identity is polynomial, hence holds also for singular blocks by continuity.
The replay independently expands both sides as polynomials in all 24 scalar
block entries: the trace has 64 monomials and the determinant has 80.
Thus no block invertibility remains as a hypothesis.

## 3. An exact obstruction to the Hodge-only approach

Let `V=wedge^3 C^6`, with orthonormal basis `e_S` indexed by the 20 triples in
lexicographic order. Define the real Hodge matrix by

`J e_S = epsilon_S e_(S^c)`,

where `epsilon_S` is the sign of the concatenation `(S,S^c)`. Then `J^2=-I`.
After an irrelevant overall phase, genuine third-compound matrices of unitaries
are unitary and preserve this skew form. Those conditions alone are insufficient.

The replay specifies a signed permutation matrix `K` on these 20 basis vectors.
Its mapping respects complementary pairs; it chooses sign +1 on every source
triple containing zero, and on its complement chooses
`epsilon_(image S)/epsilon_S`. Integer checks give

`K K^T=I`, `KJ=JK`, and `KJK^T=J`.

Extend the right side of Section 1 to an arbitrary compound-space matrix by
replacing `|det U[S,T]|^2` with `|K[S,T]|^2`, and call the result `4 Phi(K)`.
The literal matrix in `check.py` satisfies exactly

`4 Phi(K)=-2`, hence `Phi(K)=-1/2`.

**Why K is not a physical counterexample.** It sends the decomposable vector

`e_012+e_014 = e_0 wedge e_1 wedge (e_2+e_4)`

to `e_035+e_014`. The latter has

`p_013 p_045 - p_014 p_035 + p_015 p_034 = -1`.

This Plücker relation must vanish for a decomposable three-vector: contract
with the coordinate covector at zero, obtaining a decomposable two-vector,
and apply its four-index relation. Every exterior cube of a linear map
preserves decomposability. Consequently `K` cannot equal `wedge^3 U` for any
linear `U`, even up to an overall nonzero scalar. The exact negative value
therefore blocks a relaxation argument, not the original inequality.

## 4. Reproduction, evidence and handoff

Run from the repository root with Python 3.9+ and the standard library:

```sh
python3 -B research/marcus-de-oliveira/double-paired-five-term/check.py
python3 -O -B research/marcus-de-oliveira/double-paired-five-term/check.py
```

`replay.json` records the output. The checks cover the full coefficient identity,
the universal six-block polynomial identity (including a wrong-sign negative
control), exact unitarity and Hodge preservation of `K`, its negative relaxed
value and nonzero Plücker residual. Two complex Gaussian-rational unitary
fixtures additionally compare canonical signed sums with both reductions.
They have negative individual canonical coefficients but positive target sums;
these two cases are convention checks, not a universal sign proof.

The optional `probe.py` reproduces the only numerical exploration: seed 71,
12 complex Gaussian QR starts, at most 120 manifold-gradient iterations each. It needs NumPy
and SciPy; the recorded run used Python 3.9.6, NumPy 2.0.2 and SciPy 1.13.1.
The smallest value was about `-7e-16`, at rounding scale. No negative
matrix was certified, and this search gives no verdict on the inequality.
The initial attack ended without extending this scan.

**Precise remaining obligation.** Prove nonnegativity of the six-block expression
in Section 2 under the full row and column unitarity constraints, or find an
exact unitary violating it. A proof based on Section 1 must use constraints
beyond compound unitarity, complementary-minor symmetry and nonnegative squared
entries; the explicit `K` excludes all arguments using only those facts.
A promising finite algebraic route is an SOS certificate using the Plücker
relations that force images of decomposable three-vectors to remain
decomposable. No such certificate is claimed here. Quaternionic families and
their neighborhoods remain in the coordinator's separate scope.
