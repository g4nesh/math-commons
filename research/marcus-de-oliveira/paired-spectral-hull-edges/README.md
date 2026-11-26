# An open all-unitary family beyond triangular half-product hulls

**Theorem.** There is an open neighborhood in `C^3 × C^3` of

```
a = (-5+4i, -3-5i, 0),       b = (-3-2i, 1, -i)
```

such that, after repeating each of the three entries on each side twice, the
Marcus–de Oliveira determinant inclusion holds for **every** relative unitary
`U ∈ U(6)`. At the center, and throughout a sufficiently small neighborhood,
the six half-products have **five** extreme points and the 21 genuine pairing
products have **eight** extreme points. Both spectral triples are simple and
noncollinear. The neighborhood is independent of `U`; no explicit numerical
radius is claimed.

This uses the 45 established grouped Cauchy–Binet inequalities. It does not
claim the four still-open universal signs of [PR81](https://github.com/g4nesh/math-commons/pull/81),
a common law for all spectra, or a general six-dimensional theorem. It goes
beyond the triangular half-product hypothesis of
[PR87](https://github.com/g4nesh/math-commons/pull/87). The fixed-spectrum edge
method is compatible with the obstruction in
[PR85](https://github.com/g4nesh/math-commons/pull/85): the 45 inequalities do
not suffice for every possible spectrum.

## 1. Coordinates and the universally positive generators

Use the parity order

`π = (012,120,201,021,102,210)`, `s=(1,1,1,-1,-1,-1)`,

and write `p_i=∏_(r=0)^2(a_r+b_(π_i(r)))`. All indices below are zero based.
The 21 unordered pairs `(i,j)`, `i≤j`, are in lexicographic order. Set
`v_(ij)=p_i p_j`. These are exactly the genuine six-coordinate pairing
products: each 3-by-3 contingency table with row and column sums two is one
`P_i+P_j`. The replay checks all tables arising from all 720 permutations.

For each even–odd edge `e={i,j}`, let `g_e` be its six-coordinate indicator.
For the 45 unordered pairs of these nine edges define the 21-coordinate row

`R_(ef)(i,j)=g_e(i)g_f(j)+g_f(i)g_e(j)`.

There is **no implicit factor two** on off-diagonal coordinates. In particular,
when `e=f` the displayed formula includes a factor two on the whole generator.

For every physical relative unitary, group the signed Leibniz weights
`conjugate(det U) sign(σ) ∏_k U_(k,σ(k))` by contingency and take their real
parts. The resulting real vector `q` satisfies

`∑ q_(ij)=1`,  `det(D_a+UD_bU*)=∑ q_(ij) v_(ij)`,  `q·R_(ef)≥0`.

Here the diagonal matrices repeat each entry of `a,b` twice. The determinant
identity first comes from `det[(D_a U+U D_b)U*]`. Taking real parts of the
weights preserves its **formal polynomial**, since its coefficients are real
by the principal-minor/Cauchy–Binet expansion

`∑_(S,T:|S|=|T|) |det U[S,T]|² ∏_(k∉S) a_(floor(k/2)) ∏_(l∈T) b_(floor(l/2))`.

The coefficients are sums of squared minors, hence nonnegative. Each `R_(ef)`
is a positive multiple of one grouped coefficient row: `check.py` reconstructs
all 141 rows of the 21 polynomials and checks this for all 45 generators,
including their factors. Normalization follows by putting all `a_r=1,b_s=0`.
This recaps the established PR81/PR87 positivity mechanism; no new universal
unitary inequality is assumed.

## 2. Exact certificates for the eight genuine hull edges

At the center, the six half-products are

```
(36-26i, -84-108i, 51+17i, 36+42i, 4-52i, -37-107i).
```

Their strict counterclockwise hull has index order `(1,5,0,2,3)`. The minimum
inward cross product of a nonendpoint half-product with one of these edges is
694, so this is genuinely a pentagon, not a degenerate triangle.

The strict counterclockwise hull of the 21 products has coordinate indices

`(20,1,7,8,17,11,6,10)`

and corresponding pair labels

`((5,5),(0,1),(1,2),(1,3),(3,5),(2,2),(1,1),(1,5))`.

Its complex vertices, in the same order, are

```
(-10080+7918i, -5832-1704i, -2448-6936i, 1512-7416i,
 3162-5406i, 2312+1734i, -4608+18144i, -8448+12984i).
```

For a directed hull edge `(l,r)`, let

`G_(l,r)(j)=cross(v_r-v_l, v_j-v_l)`,

where `cross(x,y)=Re(x)Im(y)-Im(x)Re(y)`. The hull is counterclockwise, so
these are inward gaps. Each edge has exactly two zero coordinates, `l,r`,
and every other coordinate is strictly positive. The archived certificate gives

`G_(l,r)=∑_h λ_(l,r,h) R_h`.

The coefficients are zero exactly when the generator fails to vanish at both
exposed coordinates; they are strictly positive on **every** allowed generator.
All weights are integers, except half-integers on the seventh edge.

| Edge coordinates | Allowed generators | Face rank | Smallest positive weight |
| --- | ---: | ---: | ---: |
| 20 → 1 | 31 | 13 | 645957 |
| 1 → 7 | 27 | 13 | 187272 |
| 7 → 8 | 30 | 13 | 322028 |
| 8 → 17 | 30 | 13 | 383592 |
| 17 → 11 | 31 | 13 | 31547 |
| 11 → 6 | 33 | 13 | 37050 |
| 6 → 10 | 33 | 13 | 919359 |
| 10 → 20 | 33 | 13 | 919360 |

[certificate.json](certificate.json) lists all 45 rational coefficients for
each edge. The standard-library checker verifies all 168 coordinates of the
eight identities exactly; it does not trust floating-point LP feasibility.

For the physical determinant `d=∑q_j v_j`, normalization gives

`cross(v_r-v_l,d-v_l)=q·G_(l,r)=∑_h λ_(l,r,h)(q·R_h)≥0`.

Thus `d` belongs to all eight inward closed half-planes. Their intersection
is exactly the displayed convex polygon. This proves the fixed-spectrum
statement simultaneously for every `U∈U(6)`.

## 3. Why the result holds on an open spectral neighborhood

The following finite-dimensional argument explains the rank check and gives a
reusable constructive criterion, without treating numerical strict feasibility
as proof.

The identity `∑s_i p_i=0` is the determinant of the rank-at-most-two matrix
`(a_r+b_c)`. Multiplying by `s_v p_v` gives six real coefficient relations
`k_v·(p_i p_j)=0`, where `k_v` has diagonal entry one at `(v,v)` and, for
`i<j`, entry `s_i s_j` when `v∈{i,j}`, zero otherwise. These six vectors are
independent because of their distinct diagonal entries. They also annihilate
the constant vector `1`: `1+s_v∑_(j≠v)s_j=0`.

Let `V={x∈R^21 : k_v·x=0 for all six v}`. Then `dim V=15`. The checker verifies
that every generator is in `V` and their span has rank 15. For a fixed edge
with endpoint coordinate indices `l,r`, put

`L_(l,r)=V ∩ {x_l=x_r=0}`.

Every moving edge gap belongs to this **fixed** real space: it is a real linear
combination of `Re(p_i p_j)`, `Im(p_i p_j)` and `1`, and vanishes at its own
two endpoint indices. Exact elimination shows `dim L_(l,r)=13` for each of
the eight edges. The generators vanishing at those two coordinates have rank
13, so they span all of `L_(l,r)`, not merely a proper subspace. A particular
13-column basis for each edge is recorded in [report.json](report.json).

At the center write the gap using the strictly positive weights in Section 2.
Hold all nonbasis coefficients fixed. Express the difference between the new
and old gap in the fixed 13-column basis, and add those 13 coordinates to the
basis coefficients. This is a linear right inverse on `L_(l,r)`, hence is
continuous. The difference tends to zero as the six spectral values approach
the center. All finitely many center coefficients are strictly positive, so
they remain positive in a neighborhood. Coefficients outside this face stay
zero. No rank of a varying matrix is assumed: the generators and the chosen
basis are constant.

Strict hull orientations and all nonendpoint gaps persist by continuity.
Consequently the same eight product edges remain the entire hull, while the
half-product hull still has five vertices. Noncollinearity of both original
triples persists too (their center oriented areas are 37 and -2). Intersect
the finitely many neighborhoods from these geometric conditions and the eight
basis corrections. Throughout the resulting open set the proof in Section 2
holds for **every** relative unitary. This neighborhood depends only on the
spectra, not on `U`; no compactness argument over unitaries is needed.

This proves existence of a positive neighborhood radius. Computing a convenient
explicit radius would require bounds for the eight fixed right inverses and
spectral derivatives; no numerical radius is asserted here. The certificate
criterion is sufficient and does not cover all paired spectra.

## 4. Replay and evidence level

Python 3.9+ standard library, from the repository root:

```sh
python3 -B research/marcus-de-oliveira/paired-spectral-hull-edges/check.py
python3 -O -B research/marcus-de-oliveira/paired-spectral-hull-edges/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/paired-spectral-hull-edges -p 'test_*.py'
```

[report.json](report.json) is compared against the rebuilt report. Tests cover
changed certificate weights, reversed edge orientation, missing edges and a
false positive face coefficient, as well as exact hull boundary behavior.
The replay derives the physical contingency list, all 45 positivity rows,
all eight rational identities and the dimensions needed for the continuity
proof. The proof uses standard determinant identities and finite-dimensional
continuity; it is not a Lean proof or an independent-owner review.

Discovery used Python 3.9.6, NumPy 2.0.2 and SciPy 1.13.1, seed 20261126.
The first draw of two 3-by-2 integer arrays in `[-5,5]` produced this center;
no further spectrum scan was needed. SciPy proposed supports, then exact
Fraction elimination reconstructed the archived weights. Discovery output is
not needed by the verifier and is not the proof. The conceptual face-span
continuity criterion was suggested by the coordinating agent.
