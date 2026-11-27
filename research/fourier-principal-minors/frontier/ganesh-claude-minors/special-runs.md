# Two runs outside the scan: (35, 3) and (55, 2)

These complete the picture for orders 105 and 110. Both were run on 2026-11-27 with the table build of `pmscan.c` (called `pmscan` below). Subsets are written as hexadecimal bit masks: bit i is set when i belongs to the subset.

## (35, 3): the third presentation of 105

ord₃₅(3) = 12, so Φ₃₅ has two factors of degree 12 modulo 3 and the residue field has 3¹² = 531 441 elements. All 2³⁵ subsets were enumerated in both charts (under ten minutes per chart).

```sh
pmscan 35 3 12 1 0 1 0 2 1 0 1 2 1 2 2 list > chart_a.txt
pmscan 35 3 12 1 2 2 1 2 1 0 1 2 0 1 0 list > chart_b.txt
```

Last lines of the two outputs:

```
RESULT M=35 p=3 f=12 h=1,0,1,0,2,1,0,1,2,1,2,2 subsets=34359738368 zeros=99330 byrank=0,0,0,0,0,0,210,0,0,0,840,2100,1960,4760,4480,8120,9100,13160,11270,14140,12320,10780,2870,1330,1470,0,420,0,0,0,0,0,0,0,0,0
RESULT M=35 p=3 f=12 h=1,2,2,1,2,1,0,1,2,0,1,0 subsets=34359738368 zeros=99330 byrank=0,0,0,0,0,0,0,0,0,420,0,1470,1330,2870,10780,12320,14140,11270,13160,9100,8120,4480,4760,1960,2100,840,0,0,0,210,0,0,0,0,0,0
```

- Each chart has 99 330 vanishing principal minors. The ordinary test fails.
- The zero set of the second chart is the set of complements of the zero set of the first, as Jacobi's identity predicts when the charts are exchanged by ζ ↦ ζ⁻¹.
- 9 660 subsets are singular in both charts, at ranks 12 to 23 (280, 560, 700, 980, 910, 1400, 1400, 910, 980, 700, 560, 280). Any one of them is a covering chain of length one, so the selective test fails.
- An explicit one: {0, 1, 4, 5, 6, 8, 13, 14, 16, 19, 23, 25}, mask `2896173`. Direct check:

```
pmscan 35 3 12 1 0 1 0 2 1 0 1 2 1 2 2 check 2896173   ->  CHECK 2896173 singular
pmscan 35 3 12 1 2 2 1 2 1 0 1 2 0 1 0 check 2896173   ->  CHECK 2896173 singular
```

The two zero lists (about 1.6 MB each) are not committed; the commands above regenerate them.

## (55, 2): the third presentation of 110

ord₅₅(2) = 20, so Φ₅₅ has two factors of degree 20 modulo 2 and the residue field has 2²⁰ elements. Enumerating 2⁵⁵ subsets is out of reach, so a covering chain was searched for directly.

Every subset of size at most 8 containing 0 was tested by one elimination each (206 432 776 subsets per chart, under a minute):

```sh
pmscan 55 2 20 1 0 1 1 0 1 1 0 1 0 0 1 0 1 1 1 0 0 1 1 list lowrank 8     # 640 singular subsets, all of size 8
pmscan 55 2 20 1 1 1 0 0 1 1 1 0 1 0 0 1 0 1 1 0 1 1 0 list lowrank 8     # none
```

Take S = {0, 1, 2, 4, 33, 36, 42, 49}, mask `2041200000017`, from the first list. Its translate S + 5 = {5, 6, 7, 9, 38, 41, 47, 54} is disjoint from S, so S is contained in T = ℤ₅₅ ∖ (S + 5), mask `3f7dbffffffd1f`. Direct checks:

```
chart 1, S: CHECK 2041200000017 singular        chart 2, S: CHECK 2041200000017 nonsingular
chart 1, T: CHECK 3f7dbffffffd1f nonsingular    chart 2, T: CHECK 3f7dbffffffd1f singular
```

So S ⊂ T is a chain of two nested sets, S singular in the first chart and T in the second. Its length is min(2, 55) = 2, so it is a covering chain and the selective test fails at (2, 55). The ordinary test fails as well, since S is singular.

Why T had to be singular: translation preserves singularity in a chart, so S + 5 is singular in chart 1. Jacobi's identity gives det F[K, K](ζ) = 0 exactly when det F[Kᶜ, Kᶜ](ζ⁻¹) = 0, and ζ ↦ ζ⁻¹ exchanges the two charts because −1 is not a power of 2 modulo 55.
