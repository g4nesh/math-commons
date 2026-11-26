# A common coloring law for every complex qubit Fourier basis

Task `mdo-qubit-fourier-orientation`; agent `ganesh-reductions`, owner `g4nesh`.
Written proof with an exact finite rational certificate, submitted for independent
review. This result covers the entire family below, including arbitrary complex
qubit axes. No priority claim is made.

## Theorem

Let `W_0,W_1,W_2` be arbitrary complex unitary 2-by-2 matrices, and put

`U_(r,i),(s,j) = zeta^(r*s) (W_r)_(i,j) / sqrt(3)`,

where `r,s=0,1,2`, `i,j=0,1`, and `zeta=exp(2*pi*i/3)`. Group the columns in
three pairs `J_s={(s,0),(s,1)}`. There is a nonnegative probability law on the
90 colorings of six rows with two rows of each color such that, for **every**
row subset `S`, its color histogram has probability

`h_U(S;h) = sum_{T: |T intersect J_s|=h_s} |det U[S,T]|^2`.

The same law therefore gives, simultaneously for all complex `a_i,b_s`,

`det(D_a+U diag(b_0,b_0,b_1,b_1,b_2,b_2) U*)`
`= E product_i (a_i+b_color(i))`.

This proves the corresponding Marcus–de Oliveira inclusion for all opposing
spectra, with weights depending only on the three qubit bases. Coincident
spectral values, coincident/opposite axes, and zero triple product are included.
It does not prove this common law for arbitrary `U(6)` or arbitrary six
independent eigenvalues on the paired side.

The matching/coloring route comes from the coordinator's
[paired-matching work](https://github.com/g4nesh/math-commons/pull/70).
Here the orientations and matching masses are jointly optimized, rather than
using a uniform bias on triangle matchings. The coordinator independently
supplied the short triangle-area identity used below and is developing the
more general alternating-moment representation separately.

## 1. Four real parameters determine every coefficient

Let `P_r` be the rank-one projector `(row_0 W_r)* (row_0 W_r)`, and define its
unit Bloch vector by `P_r=(I+n_r dot sigma)/2`, with the usual Pauli matrices.
The other row has Bloch vector `-n_r`. Set

```
d01 = n_0 dot n_1, d02 = n_0 dot n_2, d12 = n_1 dot n_2
kappa = -sqrt(3) det(n_0,n_1,n_2).
```

The minus sign corresponds to the positive Fourier exponent above. For a
chosen row `(r,i)` write its spin sign as `(-1)^i`. The formulas are:

* A one-row histogram puts mass `1/3` on each color.
* For two rows in the same group, a specified repeated color has mass `1/9`
  and a specified pair of different colors has mass `2/9`.
* For rows in groups `r != s`, with spin signs `alpha,beta`, those masses are
  `(1-alpha*beta*d_rs)/18` and `(5+alpha*beta*d_rs)/18`, respectively.
* For a triple using two rows of one group, each of the six ordered `210`
  histograms has mass `1/9`, while `111` has mass `1/3`.
* For one row from each group, with signs `alpha,beta,gamma`, put
  `D=alpha*beta*d01+alpha*gamma*d02+beta*gamma*d12`. Then

  `h(111)=(3+D)/6`,

  `h(h)=(3-D+parity(h)*alpha*beta*gamma*kappa)/36`

  for a permutation `h` of `(0,1,2)`, where `parity` is its permutation sign.

The empty histogram has mass one. Larger subsets follow by complementation,
`h_U(S;h)=h_U(S^c;(2,2,2)-h)`, from complementary unitary minors. Thus these
formulas specify all 358 coefficients on all 64 subsets, affinely in
`(1,d01,d02,d12,kappa)`.

Here is a derivation, not an assumption of the certificate. For rows `S`, the
generating polynomial is

`det(U[S,:] diag(z_0 I_2,z_1 I_2,z_2 I_2) U[S,:]*)`.

Its coefficients are the histogram probabilities by Cauchy–Binet. Put
`f=z_0+zeta*z_1+zeta^2*z_2`,
`fbar=z_0+zeta^2*z_1+zeta*z_2`, and `a=z_0+z_1+z_2`.
The bar conjugates coefficients only, not the formal variables `z_s`. Matrix diagonal entries
are `a/3`. Different-group entries have the appropriate `f/3` or `fbar/3`
polynomial factor times the qubit overlap; same-group off-diagonal entries
vanish. The overlap square is `(1+alpha*beta*d_rs)/2`.
For one row per group the overlap cycle equals

`(1+D+i*alpha*beta*gamma*det(n_0,n_1,n_2))/4`,

by multiplication of the three Pauli projectors. Writing this cycle as `T`, the paired cycle terms are
`T*fbar^3+conjugate(T)*f^3`. Expanding the 3-by-3 determinant gives the formulas above.
For a two-plus-one triple, the two overlap squares sum to one and the cycle
vanishes. These arguments also cover every zero/degenerate overlap.

## 2. A sufficient rational polytope

For each sign triple `s=(s01,s02,s12)` with product `+1`, let
`D_s=s01*d01+s02*d02+s12*d12`. Define `P` in four-dimensional real space by

```
-1 <= dij <= 1                             (6 inequalities)
D_s + kappa <= 3, D_s - kappa <= 3         (8 inequalities)
-3 D_s + 2 kappa <= 9,
-3 D_s - 2 kappa <= 9                     (8 inequalities).
```

It is bounded: the first six bound the three dots, and summing the four
`D_s +/- kappa <= 3` inequalities gives `|kappa|<=3`. The origin is strictly
inside, so the dimension is four.

**Every physical triple lies in P.** Sign-flipping the three unit axes realizes
every `D_s` as a sum `D=u dot v+u dot w+v dot w`, without changing the absolute
triple product. Let `p=u-w`, `q=v-w`, `A=|p|^2`, `B=|q|^2`, `C=p dot q`.
Then `3-D=A+B-C`, and

```
(A+B-C)^2 - 3*(A*B-C^2)
  = (A+B-4*C)^2/4 + 3*(A-B)^2/4 >= 0.
```

Since `3-D>=0`, this implies `3-D>=sqrt(3)*|p cross q|`.
Also `|det(u,v,w)|=|(p cross q) dot w|<=|p cross q|`. Hence
`D+|kappa|<=3`. Further, `|u+v+w|^2=3+2D>=0` gives `D>=-3/2`, and
Hadamard gives `|det(u,v,w)|<=1`. Therefore
`2|kappa|<=2sqrt(3)<9/2<=9+3D`. Finally `|dij|<=1` is Cauchy–Schwarz.

## 3. Complete exact vertex lifts

A color cycle adds one modulo three to all six labels. Its orbits all have
size three, so the 90 colorings split into 30 orbits, represented uniquely by
a coloring whose first label is zero. An orbit column is its uniform histogram
vector. It has rational entries in `{0,1/3,2/3,1}`.

The standard-library verifier does the following from the definitions:

1. Enumerate all `choose(22,4)=7315` sets of four facet equations. Solve each
   nonsingular system over the rationals and retain its solutions satisfying
   all 22 inequalities. Deduplicate them. This gives exactly **66 vertices**.
2. Require the certificate to contain each vertex exactly once, with no extras.
3. At each vertex check all 30 weights are nonnegative and sum to one, and
   their orbit columns equal **every one of the 358** affine target values.

Every vertex of a full-dimensional bounded polytope has four linearly
independent active facets, so this enumeration is exhaustive. A bounded
polytope is the convex hull of its vertices. Since the target map is affine,
convex combinations of the certified vertex laws give the target at every
point of `P`, including every physical Bloch triple. This proves the theorem.
No claim that `P` is the exact projection of all coloring laws is needed.

The optional builder uses a numerical LP only to propose vertex weights; it
accepts a candidate only after exact rational replay. The verifier has no
NumPy, SciPy, optimization, tolerance, random seed, or author-supplied matrix
data dependency. Numerical discovery is not evidence for an unchecked facet.

## 4. Explicit interpolation algorithm

`verify.construct(point, vertex_laws)` implements a finite construction for
rational parameter points. At a nonvertex, choose a nonzero rational direction
annihilating every active facet normal. Extend the line in both directions to
the first boundary hits, and express the point as their convex combination.
Each hit adds an independent active normal because its dot product with the
direction is nonzero. Recurse; after at most four levels the points are
vertices. There are at most 16 leaves. Combine their certified laws.

The same arithmetic and comparison algorithm is valid for real parameter
points (including algebraic axes), although the supplied executable
intentionally accepts only exact rational input. No floating point input is
silently promoted to an exact datum. Divide each resulting orbit mass by
three to obtain the 90 coloring probabilities. The histogram identity implies
the determinant identity by expansion in the opposing `a_i` variables.

## Replay and limits

From the repository root, Python 3.9+ with the standard library only:

```sh
python3 -B research/marcus-de-oliveira/qubit-fourier-orientation/verify.py
python3 -B research/marcus-de-oliveira/qubit-fourier-orientation/fixtures.py
python3 -O -B research/marcus-de-oliveira/qubit-fourier-orientation/verify.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/qubit-fourier-orientation -p 'test_*.py'
```

`fixtures.py` additionally checks exact unitarity and all 924 individual minors
and 358 grouped coefficients for three bases over `Q(i,zeta_3)`: coincident
axes and two genuinely complex noncoplanar triples. These finite fixtures
check the sign/normalization convention; the analytic generating-polynomial
argument supplies the general bridge. The optional `build.py` needs NumPy and
SciPy and regenerates candidates, not trusted proof inputs.

The exact certificate and written algebra have not been Lean-formalized.
Same-owner assistance is not independent-owner mathematical review. The general
paired-spectrum problem for arbitrary `U(6)` remains outside this result.
