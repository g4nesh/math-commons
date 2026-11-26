# Exact dual rays for the split pattern

For the `(2,2,1,1) × (2,2,2)` pattern, the cone

\[
 C=\operatorname{row}(A)\cap\mathbb R_{\ge0}^{33}
\]

has **935 extreme rays**, in **62 orbits under the physical group of order 24**. Of these, 65 rays in ten orbits are positive multiples of grouped coefficient rows. The remaining 870 rays form 52 orbits. Their primitive integer coordinates range from zero to four; supports have sizes three through thirteen.

This is an exact polyhedral classification. It **does not establish nonnegativity of the 52 remaining physical sign types for every unitary**. The ordinary coefficient rows are identified algebraically; any analytic positivity theorem must still use the physical quantities and hypotheses in the [split-pattern reduction](../split-one-eigenvalue-pair/README.md). No unitary counterexample is claimed. Independent-owner mathematical review is pending.

The complete human-readable [orbit table](ORBITS.md) gives all 62 representatives and all 33 coordinate tables. [orbits.json](orbits.json) and [rays.json](rays.json) contain their exact integer coordinates. The indexing agrees with PR99: zero-based lexicographic order of flattened `4×3` tables, row sums `(2,2,1,1)`, column sums `(2,2,2)`. For example the first non-coefficient orbit, ID 3, is

\[
 e_{23}+e_{25}+e_{26}+e_{27}+e_{28}+e_{31}+e_{32}.
\]

Here `e_i` is a coordinate unit vector; a ray is a linear functional on the 33-table signed representation, not itself a probability law.

## Exact construction and proof obligations

Let `ROW=(0,0,1,1,2,3)` and `COL=(0,0,1,1,2,2)`. For every permutation `p` of six symbols, form `T_rs=#{i: ROW_i=r, COL_p(i)=s}`. There are exactly 33 distinct tables. Expand

\[
 P_T(x,y)=\prod_{r=0}^{3}\prod_{s=0}^{2}(x_r+y_s)^{T_{rs}}.
\]

The matrix `A` consists of the coefficients of these polynomials, ordered by the 192 lexicographic seven-exponent monomials. Fraction Gaussian elimination gives rank 21 and a twelve-dimensional nullspace. Put a primitive integer basis of this nullspace in the rows of `N`. Then `row(A)=ker(N)`, so

\[
 C=\{z\in\mathbb R^{33}:Nz=0,\ z\ge0\}.
\]

Every nonzero ray has a unique representative on `sum(z)=1`; this slice is bounded, and its vertices are exactly the extreme rays. Thus enumeration of the normalized polytope is equivalent to enumeration of the cone rays, with no direction at infinity omitted. [bounded.ine](bounded.ine) has the twelve equations `Nz=0`, the normalization equation, and all 33 coordinate inequalities. Its 10,000-output guard is a finite resource bound, not a truncation of this result: the run terminates normally with 935 vertices and zero unbounded rays, after visiting 199,641 bases.

For each reported primitive nonnegative vector `z`, the standalone checker verifies `Nz=0` and

\[
 |S|-\operatorname{rank}(N_{:,S})=1,\qquad S=\{i:z_i>0\}.
\]

This is an exact extremality certificate. If this kernel has dimension greater than one, choose a nonproportional vector supported on `S` in it and perturb `z` both ways by a sufficiently small real multiple; positivity on `S` is preserved, so the ray is not extreme. Conversely, any decomposition of `z` into two cone vectors has both summands supported on `S`; when this kernel is one-dimensional, both are nonnegative multiples of `z`.

**Completeness has a separate evidence boundary:** it is supplied by the completed exact lrslib reverse search, not inferred from testing the 935 candidates. The Python checker independently rebuilds the matrix, validates the archived enumeration input/output and every representative's extremality, and checks that all reported rays partition into full physical orbits. Re-running the pinned exact enumerator reproduces the global enumeration. No floating-point hull, LP tolerance, sampled directions, or presumed genericity is used.

For any real candidate signed table vector `q`, finite-dimensional cone duality now says that a nonnegative `m` with `Am=Aq` exists exactly when every listed ray `z` satisfies `z·q>=0`. Indeed `(ker A + R^33_+)^*=row(A)∩R^33_+`, and this finitely generated polyhedral sum is closed. The constant normalization row is already in `row(A)`. This is a reduction of the physical problem, not its resolution.

## Physical symmetry and concrete invariants

The action is precisely

\[
 S_2^{\mathrm{paired\ rows}}\times S_2^{\mathrm{singleton\ rows}}\times S_3^{\mathrm{columns}},
\]

with 24 elements. Paired row groups may swap, singleton row groups may swap, and column groups may permute. The checker verifies these actions preserve the actual coefficient rows. We do not enlarge the group to arbitrary four-row permutations or transpose the rectangular table.

A ray is made primitive with positive coordinates and represented by the lexicographically least element of its orbit. Each record includes support size, primitive coordinate sum, largest coordinate, orbit and stabilizer sizes, plus two invariants tied to the split geometry:

- **Equal-singleton mass:** sum of the ray coordinates on tables where the two singleton rows occupy the same column group.
- **Singleton imbalance:** the L1 distance from the ray to its image under swapping those two singleton rows.

The equal-column relation is preserved by the physical group; the singleton swap is central. Thus these quantities are invariant. They help organize the physical chart, but are **not claimed to be a complete separating set of invariants**. The explicit group action determines the orbit partition.

There are 25 orbits of size 24, 21 of size 12, 12 of size 6, three of size 3, and one of size 2, totaling 935 rays. Exactly ten of these orbits contain a primitive coefficient row. “Other” does not assert that no different known theorem can prove them.

## Replay

From the repository root, only Python's standard library is needed to check the archived exact certificates:

```sh
python3 -I -B research/marcus-de-oliveira/split-dual-rays/check.py
python3 -I -O -B research/marcus-de-oliveira/split-dual-rays/check.py
python3 -I -B research/marcus-de-oliveira/split-dual-rays/test_check.py
```

Six tests include altered sign/scaling, a vector outside the row space, a non-extreme positive sum, Boolean coordinate rejection, a missing output vertex and a missing completion footer. All validation uses explicit exceptions and remains active under `-O`.

To independently re-enumerate, download the official [lrslib 7.3a source archive](https://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/lrslib-073a.tar.gz). Its SHA256 is `7fb16ad5c2b7380ef8deaecbfaf53599eef074fd225012aa431b76ee0bded6a6`. Verify that hash before building. In the extracted `lrslib-073a` directory:

```sh
clang -O2 -I lrsarith-011 -DMGMP -DGMP \
  -o lrsmgmp lrs.c lrslib.c lrsarith-011/lrsgmp.c \
  lrsdriver.c lrsarith-011/mini-gmp.c
./lrsmgmp /absolute/path/to/split-dual-rays/bounded.ine > /tmp/split-rays-new.out
```

After building, the bounded fresh replay checks the complete set against the archive:

```sh
python3 -I -B research/marcus-de-oliveira/split-dual-rays/reenumerate.py \
  --lrs /absolute/path/to/lrslib-073a/lrsmgmp
```

This command uses a 60-second timeout and requires a successful process exit, empty stderr, normal completion totals and exact equality of all 935 rays. [enumeration-replay.json](enumeration-replay.json) records a successful second enumeration, including the locally built binary hash. It does not download or install anything.

The original run used Apple clang 17.0.0 on arm64 Darwin and took 15 elapsed seconds. The tool identifies itself as `lrslib_v.7.3_2024.5.30(64bit,lrsgmp.h,mini-gmp)`; `mini-gmp` supplies arbitrary-precision integer arithmetic, while `64bit` is the host word size. The stdout, with trailing whitespace removed, is retained in [lrs-output.txt](lrs-output.txt); the original process exited zero with empty stderr. Paths and timing in a fresh transcript may differ. Its exact primitive ray set must equal [rays.json](rays.json); use `read_lrs` in [algebra.py](algebra.py) or point `validate_output` in [check.py](check.py) at the fresh output and compare the sorted rays.

The [official lrslib manual](https://cgm.cs.mcgill.ca/~avis/C/lrslib/USERGUIDE.html) documents H/V representation and exact reverse search. The original split matrix and physical problem are from [PR99](https://github.com/g4nesh/math-commons/pull/99); this directory reconstructs them without importing that source. The result is an exact computational enumeration with individual rational certificates, not a Lean formalization or independent-owner proof review.
