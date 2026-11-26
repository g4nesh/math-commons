# What solved coarser spectra prove about the split rays

**Result.** For `(2,2,1,1) × (2,2,2)`, solved proper coarsenings establish 17 of the 52 non-coefficient dual-ray orbits from the [exact split-cone enumeration](../split-dual-rays/README.md). These are 186 individual rays. Together with the ten coefficient orbits, only **35 non-coefficient orbits** remain for a physical positivity argument.

The reductions use the independently reviewed [Theorem D](../minors/theorem-d.md), and, for merging the singleton rows, the written [universal both-paired theorem](../double-paired-five-term/GRASSMANN.md), whose independent-owner review remains pending. The new exact algebra audit likewise awaits independent-owner review. It does not prove a negative physical example or resolve the entire split family.

A second result describes the limit of this route: **no nonnegative linear combination of inequalities from any proper spectral coarsenings, even combined with the fine coefficient inequalities, can establish any of those 35 remaining extreme-ray inequalities purely as a coefficient-vector identity.** This is an extremal-cone obstruction, not an obstruction to nonlinear identities, unitary geometry or other analytic proofs.

## The nine elementary identifications

Keep the fine 33-table indexing and seven spectral variables `(x0,x1,x2,x3,y0,y1,y2)` from PR104. Identifying two variables on the same side adds their exponents. The coarsened coefficient row for an exponent is the sum of all fine rows that become that exponent. Denote the resulting matrix, still with 33 fine columns, by `A_c`.

| Coarsening | Multiplicities after merging | Rank of `A_c` | Coefficient rows | Established law |
| --- | --- | ---: | ---: | --- |
| `row01` | `(4,1,1) × (2,2,2)` | 8 | 96 | Theorem D, row side |
| `row02`, `row03`, `row12`, `row13` | a reordering of `(3,2,1) × (2,2,2)` | 12 | 122 | Theorem D, row side |
| `row23` | `(2,2,2) × (2,2,2)` | 15 | 141 | Universal both-paired theorem |
| `column01`, `column02`, `column12` | `(2,2,1,1) × (4,2)` | 8 | 96 | Theorem D, column side |

Theorem D supplies a law depending only on the relative unitary and the multiplicity pattern, valid identically in all spectral parameters. This stronger quantifier is essential. A separate pointwise hull inclusion without common coefficients would not justify a row-space argument.

For fixed physical `U`, let `q` be the real canonical signed vector on the 33 fine tables. A solved coarser law gives nonnegative masses on coarser tables. Each such table is realized by a permutation, so choose one such permutation and place its mass on the corresponding fine table. This gives `m >= 0` with

`A_c m = A_c q`.

If a fine dual ray `z` belongs to `row(A_c)`, then

`z·q = z·m >= 0`,

because `z` is coordinatewise nonnegative. The row-space expression for `z` may have negative coefficients; that does not matter, since it is the equality with a nonnegative *law* that proves positivity. Merely checking that `z` is constant on coarse table fibers would be insufficient: the exact row-space condition is what is certified here.

## An explicit reduction for orbit 3

Orbit 3 is

`z = e23 + e25 + e26 + e27 + e28 + e31 + e32`.

Merge row colors 1 and 2, meaning the paired block of zero-based physical rows 2,3 and the singleton row 4. The resulting row groups have sizes `(2,3,1)`. Theorem D therefore applies for every complex unitary, independently of the opposing spectrum.

Write the coarsened row variables as `(X0,X1,X2)=(x0,x1=x2,x3)`. If `[k]` denotes the coefficient vector of the monomial with exponent tuple `k` in `(X0,X1,X2,y0,y1,y2)`, the checker verifies the exact 33-coordinate identity

```text
z =  5[0,0,0,2,2,2] - 5[0,0,1,1,2,2] - 3[0,0,1,2,1,2]
   - 3[0,1,0,1,2,2] - 2[0,1,0,2,1,2] + 3[0,1,1,0,2,2]
   + 2[0,1,1,1,1,2] +  [0,2,0,0,2,2] +  [0,2,0,1,1,2]
   -  [0,2,1,0,1,2].
```

Equivalently, this ray is one on exactly three of the 15 coarsened tables and zero on the others. The three tables, written as row triples, and their fine fibers are:

| Coarse table | Fine indices |
| --- | --- |
| `(1,1,0) / (1,0,2) / (0,1,0)` | 23, 26 |
| `(1,1,0) / (1,1,1) / (0,0,1)` | 25, 27, 28 |
| `(2,0,0) / (0,2,1) / (0,0,1)` | 31, 32 |

Thus its physical value is the sum of these three probabilities under the coarser law. The exact row-space identity ensures that this sum is independent of which valid coarser law is chosen. Boundary or zero values cause no problem.

## Exact coverage and why further mergers add nothing

The covered non-coefficient orbit IDs are

`3, 8, 9, 13, 14, 17, 21, 23, 24, 25, 27, 29, 32, 47, 48, 49, 50`.

| Orbit | An elementary reduction | Law |
| ---: | --- | --- |
| 3 | `row12` | Theorem D |
| 8 | `row23` | Both-paired theorem |
| 9 | `row23` | Both-paired theorem |
| 13 | `row12` | Theorem D |
| 14 | `row23` | Both-paired theorem |
| 17 | `row23` | Both-paired theorem |
| 21 | `row23` | Both-paired theorem |
| 23 | `row23` | Both-paired theorem |
| 24 | `row23` | Both-paired theorem |
| 25 | `row23` | Both-paired theorem |
| 27 | `row23` | Both-paired theorem |
| 29 | `row23` | Both-paired theorem |
| 32 | `row02` | Theorem D |
| 47 | `row23` | Both-paired theorem |
| 48 | `row23` | Both-paired theorem |
| 49 | `row23` | Both-paired theorem |
| 50 | `row23` | Both-paired theorem |

The remaining non-coefficient IDs are

`10, 11, 12, 15, 16, 19, 20, 22, 26, 28, 30, 31, 33, 34, 35, 36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 56, 57, 58, 59, 60, 61`.

There are 684 individual rays in these remaining orbits. Proper coarsening covers 26 total orbits and 245 individual rays. The one coefficient orbit it misses is ID 38, containing six rays; that orbit is already nonnegative by the ordinary squared-minor/Cauchy–Binet argument. Hence the solved union contains 27 orbits and 251 rays, leaving the 35 listed above.

Every proper merging of color blocks within either spectrum factors through at least one of the nine elementary mergers. Its pulled-back row space is therefore contained in an elementary row space. This covers all further row mergers, column mergers and combinations of the two. The collection is invariant under the physical group; the checker also tests coverage on all 935 group images.

Put `C = row(A) ∩ R^33_+` and `C_c = row(A_c) ∩ R^33_+`. Every `C_c` is a subcone of `C`. Suppose an extreme ray `z` of `C` were a sum of vectors from several `C_c`. Extremality forces every nonzero summand to be proportional to `z`; hence `z` must itself lie in one of those `C_c`. Consequently their sum covers exactly the fine extreme rays covered by their union. Adding the fine coefficient generators changes this conclusion only for extreme rays proportional to a coefficient row.

For every uncovered ray and every elementary coarsening, the certificate gives a vector `k` with `A_c k=0` but `z·k != 0`. This proves nonmembership in `row(A_c)` exactly. For covered cases, it gives an explicit rational combination of actual coarsened coefficient rows equal to `z`. The certificate contains all nine decisions for every one of the 62 representatives, not just successful examples.

## Replay and evidence limits

From the repository root:

```sh
python3 -I -B research/marcus-de-oliveira/split-dual-coarsening/check.py
python3 -I -O -B research/marcus-de-oliveira/split-dual-coarsening/check.py
python3 -I -B research/marcus-de-oliveira/split-dual-coarsening/test_check.py
```

Only the standard library is needed. The checker reconstructs the fine matrix via the exact enumeration's [algebra.py](../split-dual-rays/algebra.py), re-verifies extremality of every representative, constructs all nine coarsenings, checks each membership/separation certificate, and verifies physical-orbit coverage. No floating-point LP, numerical unitary samples or tolerance is used. `--write` regenerates the local certificate and report deterministically; it does not write to the shared ledger.

The [certificate](certificate.json) and [report](report.json) use the immutable PR104 ray IDs. Their completeness as a classification of *all* fine extreme rays inherits PR104's pinned exact enumeration; the coarsening and extremality checks on each listed vector are independent exact rational calculations. The analytic transfer inherits the stated status and hypotheses of Theorem D and the both-paired theorem. This is same-owner new research, not an independent review of either source.
