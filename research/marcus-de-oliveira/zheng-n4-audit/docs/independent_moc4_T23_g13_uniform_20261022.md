# Audit: independent_moc4_T23_g13_uniform_20261022 (seven-contact set T_24 + {13})

Evidence level: EXACT for all rational identities (4)-(13), (10) and the degenerate-kernel identities; PROOF-READ for the exhaustiveness of the chart and the sign elimination (also NUMERIC: no counterexample in 400000 draws); CITED for the extra-tie cases, the one-circle theorem and App. B.

PDF: proofs/manin_transfer_20261021/independent_moc4_T23_g13_uniform_20261022.pdf (4 pages). Replay: `replay/independent_moc4_T23_g13_uniform_20261022__check.py`.

## Deviations and limits

- Chart sign convention as in the T_24 audit: the point at infinity contributes the factor 1 to each F_sigma.
- The kernel matrices M_j were computed by our own code from the reduced rational functions H_j. We did not use the document's general kernel formula as an input; it agrees with ours (it is the identity N(X)D(Y) - N(Y)D(X) = (X - Y)(X,1)M(Y,1)^T).
- The (15) identity for g5/g9 was checked in floating point at 3 points (with exact (6) behind it).
- The companion JSON group witnesses were not read; we computed our own (find_g).

## A. Statements

Theorem 1 (p. 1). For every x in X with P c(x) >= 0 vanishing on T_24 ∪ {13} = {2,3,4,6,7,8,13} (Table 5, seven-set 6; extra contacts allowed), Phi_{c(x)}(O) >= 0 for all O in SO(4). All spectra are covered: collisions (App. B), real A-spectrum (one-circle theorem), extra ties 0, 1, 10, 22 (cited cases), 16 and 17 allowed throughout. In the eight-distinct remaining layer, c(x) lies in the minor cone (delta <= 0).

## B. Dependencies

| input | cited | own resolution | status |
|---|---|---|---|
| tie 0 | supplement T3 | six-subsets of T+0 through 0 lie in the orbits of paper T_3 and T_7; also f54(0) > 0, so a tie at 0 forces delta <= 0 directly | CITED / EXACT (see note) |
| tie 1 | T14 | orbits of T_3, T_5, T_15, T_42, T_43, T_50 (T_15 = {1,2,3,4,6,7}) | CITED |
| tie 10 | T10 | orbits of T_1, T_10, T_11, T_18, T_20, T_47 | CITED |
| tie 22 | T11 | orbits of T_10, T_11, T_12, T_39, T_52, T_67 | CITED |
| one-circle theorem | paper Thm C.1 | real x | CITED |
| App. B | paper | collisions | CITED |

Note: under the usual "historical label = paper index - 1" rule, the label T3 for tie 0 would mean paper T_4. T_4 does not fit in T+0. The label must mean paper T_3, or the step is simply superseded by f54(0) > 0. Harmless.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | sources of T: 0,1,2,3,4,5,6,53,54,67, only 54 non-minor; delta = g0-g1-g16 = g5-g16 = g9-g17; decomposition (1) | EXACT (shared with T_24 audit) |
| 2 | delta <= 0 gives the minor cone | EXACT (replacement identities f54+f1, f54+f4) |
| 3 | reduction to (2): eight distinct labels, g0, g1, g10, g22 > 0, g13 = 0, delta > 0 | PROOF-READ + CITED |
| 4 | chart exhaustive: s, t, r, u nonzero and real; d, h, ell nonzero; X nonreal | PROOF-READ (identities w = ch, z-v = th(z-c), r = -d(z-1)/X EXACT) |
| 5 | (4) H9, H1; (5) H10; (6) H0 - H1 | EXACT |
| 6 | (7): -A1 q/u > 0, -tr/q > 0 | PROOF-READ: g1/g9 = -A1 q/(u abs(L1)^2), g10 = -(tr/q) g1 |
| 7 | (8), (9) as gap identities (difference is a real constant, so imaginary parts agree) | EXACT |
| 8 | (10) g22 = e g12 + lambda g13, with lambda = r(rst - su + tu)/(tu(r-u)) | EXACT |
| 9 | kernel representation g_j/g9 = Hermitian form / abs(J)^2 | EXACT (plus a float sanity check) |
| 10 | (11) det M13; (12), (13) | EXACT |
| 11 | M13 indefinite and det != 0, hence usrR1 > 0 | PROOF-READ (semidefinite real form vanishing at nonreal (X,1) is 0) |
| 12 | sign elimination to regions (15), (16) | PROOF-READ + NUMERIC (400000 draws, 9954 satisfying the hypotheses, 0 outside (15)/(16)) |
| 13 | (15) impossible by the g5/g9 identity | EXACT (identity (6)) + float check + PROOF-READ signs |
| 14 | (16) impossible: M12 - (sr/R1) M13 negative definite | EXACT (12), (13) + PROOF-READ signs |
| 15 | degenerate q+t-1 = 0: M13 = -t(s-1)(t-1)B [[1,-1],[-1,1]]; B = 0 gives a = u/q, g0 = q g1, g5 = -t g1; e/a = (t-u)/u < 0 | EXACT |
| 16 | t = 1 gives delta = 0 | EXACT (from (8) with g13 = 0) |

## D. Falsification

The sign-elimination implication was sampled with random (s,t,r,u): 400000 draws (seed 13), 9954 met (7), a > 0, e > 0 and det M13 < 0, and all of them fell in (15) or (16). The theorem's eight-distinct layer was not sampled directly on the seven-contact curve.

## E. Verdict

confirmed-modulo-citations. Citations: paper cases for ties 1, 10, 22 (tie 0 is also covered by f54(0) > 0), paper Thm C.1 and App. B. Presentation notes: the "(X,1)M(X,1)^T" notation means the Hermitian form; the label "T3" for tie 0 does not follow the usual index shift.
