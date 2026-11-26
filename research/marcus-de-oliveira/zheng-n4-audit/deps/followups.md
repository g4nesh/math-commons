Evidence level: EXACT for every generator list, new-contact list, orbit cover, group element, printed
action and graph-closure statement below (own code calling `lib/moc4lib.py`; the supplement's
`cases/proof-order.json` loaded as data only), and for the symbolic identities of item 1 (own sympy
code). PROOF-READ for what each proof document says it uses and where (text extracts, pages of the
compiled PDFs). CITED where a statement rests on a first-pass report in `docs/` that I did not redo.

# Follow-ups of the dependency map

## Deviations and limits

- Two trivial reads of our own `dependency-map.json` (no computation, under a second) were run with
  the system Python outside the job limiter before I noticed. Every later command, including all
  computations below, went through the limiter.
- "Dangerous new contact" j of a case T means: j is not in T and the contact cone of T ∪ {j} still
  has a non-minor generator (`generators(T | {j})` contains an index with `F_IS_MINOR` false). This
  is the coarsest criterion. A document that splits the cone into branches can use a finer one: in
  a branch where only one non-minor coefficient d is present, a contact where that profile is
  positive forces d = 0. Where a report uses the finer criterion I say so (T_85).
- "Earlier cover" of (T, j) means: some 6-subset of T ∪ {j} other than T lies in the G-orbit of a
  case of a strictly lower Table 6 stage (`orbit_rep`), or T ∪ {j} lies in the orbit of a Table 5
  set. Stages are taken from `dependency-map.json` (they equal Table 6).
- Item 4 compares my computed lists with what the first-pass reports record (their dependency
  lists, step tables and witness tables). I did not re-run each report's group search. I did check
  that every case named in a report for (k, j) appears in my exact list of earlier covers of (k, j).
- Scripts and raw outputs are in the worker's private scratch, not in this tree. Every check is a
  one-line call of a library function, named below, so it can be replayed directly.
- Nothing from the supplement was executed. No mathematics inside the case proofs was re-audited
  here, except the boundary16 interval argument (item 1), which I re-derived.

## 1. Possible two-way dependence T_8 versus T^(7)_1

What was checked. I read `root_moc4_T7_uniform` (T_8, 4 pp.), `root_moc4_T7_boundary16` (1 p.) and
`root_moc4_T7_boundary23` (1 p.) in full. I re-derived with sympy every formula the boundary16 note
takes from the T_8 note.

Set identities (EXACT). `find_g(T_8 ∪ {16}, AUX7[1])` = (rho, rho', eps) = ((3,2,1,0), (3,2,0,1), 1).
`find_g(T_8 ∪ {23}, AUX7[2])` = ((1,2,3,0), (3,0,2,1), 0). So the two notes prove T^(7)_1 and T^(7)_2.
The dangerous new contacts of T_8 are exactly {4,10,16,18,19,22,23}. This is the list on T_8 p. 3.

What boundary16 imports. All of it is on pp. 1-3 of the T_8 note, before the paragraph "All
seven-contact boundaries satisfy the stated support inequality" (p. 3):
- the canonical profile (1), and the finite replacements f46 + f0 and f46 + f1 that settle d <= 0
  (p. 1);
- the rational chart A = (∞, v, xv, z), C = (r, 1, 0, w), with the inverse (4) and its nonvanishing
  denominators (pp. 1-2);
- the identities (2): g16 = (1 - s) g12, g9 = -r b, a = -x d (p. 2), and with (1), g23 = g9 - d;
- the signs (3): x < 0, r < 0, 0 < s <= 1 (p. 2);
- the trace formulas (5) for d/b and g3/b, and the g10 identity (6), which has no division by rx - 1
  (p. 2);
- the sign result (7): t > 0, L > 0 (p. 3);
- the treatment of the collision and circle strata: Appendix B and the supplement's six-point circle
  theorem.

Where these are proved. Each is derived from three things only: the six contacts of T_8,
nonnegativity of the other gaps (g3 >= d, g9 >= d, g16 >= 0, g18 >= 0, a >= 0), and the chart. None
uses Theorem 1 of the T_8 note or either seven-contact statement. They hold on every supported
source with contacts containing T_8 and d > 0, so in particular on the face T_8 ∪ {16}. The first-pass
report on the T_8 note replayed (1), (2), (4), (5) and (6) EXACT (CITED).

My replay (EXACT, sympy). At s = 1, r = -R, x = -X:
- L = R(1 - RX) + (R + 1)t, J = -R - t, K = 0, c0 = rx - t, and |DZ|^2 = |Z|^2 |w|^2.
- (5) reduces to boundary16's d = tR^2/|w|^2. Its g3 denominator equals |(R+1)w + R^2(1+X)|^2.
- Substituting |w|^2 = tR^2/d: g3 >= d is equivalent to 2(1+X) Re w <= LG. Using t > 0, the g10
  identity (6) gives g10 >= 0 if and only if 2(1+X) Re w >= LV. Both with exactly boundary16's LG
  and LV.
- LG - LV = (X - R)(1 + d)L / ((R + 1)d) holds identically.

So the boundary16 argument is the s = 1 specialisation of the T_8 note's chart identities and nothing
more. (g16 = 0 and g12 > 0 force s = 1, because d = s g12 - g18 > 0.)

Graph (EXACT, `proof-order.json`). T_8 is not in the closure of T^(7)_1. T^(7)_1 is in the closure of
T_8. A theorem edge T_8 -> T^(7)_1 would create a cycle, but no such edge exists. The boundary23 note
uses only the cubic, Sing C, the facet profiles and the old SOC atom. It imports nothing from T_8.

Related observation. T_9 ∪ {9} contains an image of T_8 (stage 4, later than T_9 at stage 3). The T_9
polar-reduction note correctly routes this face to T^(7)_1 (stage 1). So results at stages 1 and 3
rest on chart identities that are written down only in a stage-4 document. They are unconditional
algebra, so the order is still sound.

Conclusion: RESOLVED (issue: presentation). This is lemma-level reuse of identities, not a logical
circle. The boundary16 note is not self-contained: it cites "adjacent original-product certificates"
instead of naming the T_8 note and its equations (2), (3), (5), (6), (7). The imported identities are
proved in the T_8 note itself only as "direct real traces" checked by a script. Our first-pass replay
and the derivation above supply them.

## 2. New contacts without a recorded witness: T_38 (10, 15, 20) and T_65 (15)

### T_38 = {1,2,3,8,11,12}, document `root_moc4_T37_uniform` (stage 5)

Library (EXACT). `generators(T_38)` = [0,1,2,9,10,15,53,56], with f56 the only non-minor profile.
This equals the document's source (1). The values of (f1, f2, f15, f9, f10, f56, f0, f53) at the
contacts (17, 15, 4, 10, 9, 6, 20, 14) form the identity matrix, as p. 1 states. So the coefficients
(q, p, j, c, e, d, b, a) are those gaps.

The dangerous new contacts are {4,9,10,14,15,17,20}. Contacts 10, 15 and 20 are the zero walls of the
coefficients c, p and b. For each, the cone of T ∪ {j} keeps the non-minor generator f56, so each needs
an argument.

How the document handles them (PROOF-READ). T_38 is proved by a whole-cone argument, not a curve
argument. The radial necessity G9 = 0 is combined with seventeen positivity charts over the closed
coefficient region.
- The four zero walls q, j, e, a (contacts 17, 4, 9, 14) are sent to earlier cases by the table on
  p. 1.
- p, c, b "remain arbitrary nonnegative real numbers" (p. 1). The strictness paragraph (p. 4) says
  that every strictly positive witness monomial avoids the p-fraction parameter, the free c and the
  b-fraction parameter: "p = c = b = 0, separately or together, causes no loss of strictness".

So contacts 15, 10 and 20 are handled inside the theorem, by the same charts.

The p. 1 table, replayed (EXACT, library convention rho ∘ sigma ∘ rho', inverse first if eps = 1).
Every printed image matches and lies in T ∪ {tie}:

| zero | tie | case | stage | (row, col, tr) | image |
|---|---|---|---|---|---|
| q | 17 | T_46 (printed T45) | 2 | ((1,3,0,2), (1,2,0,3), 1) | {1,2,8,11,12,17} |
| j | 4 | T_4 (printed T4) | 1 | ((0,1,2,3), (0,3,1,2), 0) | {1,2,3,4,8,11} |
| e | 9 | T_8 (printed T7) | 4 | ((0,1,3,2), (0,3,2,1), 0) | {1,2,3,8,9,11} |
| a | 14 | T_11 (printed T10) | 3 | ((3,1,2,0), (3,2,1,0), 1) | {1,2,3,8,12,14} |

Earlier covers would also exist for the internally handled walls: T_39 or T_44 for 10, T_29 for 15,
T_36 or T_40 for 20 (all stage 1). A citation-based repair is therefore available.

Caveat. The first-pass report (could-not-check) did not replay the 17 chart certificates (UNCHECKED).
The strictness on these three walls depends on them.

Conclusion: RESOLVED (no issue at the level of dependency structure). Contacts 10, 15 and 20 are zero
coefficients covered by the theorem's own closed charts. The supplement graph correctly records no
prerequisite for them. Whether the charts really stay strict on those walls is part of the open
T_38 certificate check.

### T_65 = {0,4,5,6,9,12}, document `root_moc4_T47_whole_assembly` and the T47 component notes (stage 6)

Library (EXACT). `generators(T_65)` = [0,1,2,5,7,9,35,45], with f35 the only non-minor profile.
This equals the source (1). The values of (f2, f9, f1, f45, f5, f35, f0, f7) at (15, 7, 17, 1, 11, 2,
18, 14) form the identity matrix. So q is the gap at contact 15.

The dangerous new contacts are {1,7,11,14,15,17,18}. The zero walls p, c, e, a, j, b (contacts 7, 1,
11, 14, 17, 18) go to table (2) on p. 1. With the historical labels read as T43 = T_44, T4 = T_4,
T7 = T_8, T37 = T_38, T41 = T_42, I found an element of G with these images (EXACT, `find_g` on
6-subsets):

| zero | tie | case (stage) | (row, col, tr) | image |
|---|---|---|---|---|
| p | 7 | T_44 (1) | ((0,1,2,3), (0,1,3,2), 0) | {0,4,5,6,7,9} |
| c | 1 | T_4 (1) | ((0,1,2,3), (0,1,3,2), 0) | {0,1,4,5,6,9} |
| e | 11 | T_8 (4) | ((0,1,3,2), (0,2,3,1), 0) | {0,4,5,6,9,11} |
| a | 14 | T_38 (5) | ((0,1,2,3), (0,2,1,3), 0) | {0,4,5,6,9,14} |
| j | 17 | T_42 (1) | ((2,1,0,3), (1,3,0,2), 0) | {0,4,5,6,12,17} |
| b | 18 | T_42 (1) | ((2,1,3,0), (2,3,0,1), 0) | {0,4,5,6,12,18} |

All are earlier than stage 6. Read with paper numbers instead (T_43, T_4, T_7, T_37, T_41, T_41),
only T_4 covers its face. The labels are historical (EXACT).

Contact 15 (q = 0). The only cases with an image in T_65 ∪ {15} are T_65 itself and T_66, T_67, T_69,
T_73, T_75 and T_86, all later. No Table 5 set matches. So no earlier case can supply this wall, and
it must be handled inside the proof.

The assembly (p. 1) keeps "q ≥ 0, without removing the cases covered by q = 0" and says (p. 3) "q = 0
is stated explicitly in every component". I checked all ten component notes (PROOF-READ, text
extracts). Each states q ≥ 0 in its hypotheses and names q = 0 among its covered walls or strictness
markers: I/II (p. 1-2, "includes q = 0"), IIIa-minus, IIIa-plus ("strictness even at p = q = 0"),
IIIb-minus-swapped, IIIb-plus-swapped, IIIb-plus-unswapped, high-j ("including q = 0"), low-j-p-below,
low-j-low-e (q = 0 endpoint theorem) and low-j-high-e.

Seven of the ten have first-pass reports that replayed strictness over all allowed zero patterns
(CITED). Three have no first-pass report: `independent_moc4_T47_I_II_uniform`,
`independent_moc4_T47_IIIb_plus_unswapped_uniform` and `independent_moc4_T47_high_j_uniform`.

Conclusion: RESOLVED (no issue in the dependency structure). The q = 0 wall (contact 15) is part of
every component's closed region, and no earlier witness exists or is needed. UNRESOLVED remainder: q = 0
in the three unreported components rests on their unreplayed certificates. The whole assembly has no
first-pass report.

## 3. Missing edge T^(7)_27 (S67) -> T_67

Use confirmed (PROOF-READ). `root_inventory66_native_whole` p. 1 states:
- the source on T_67 = {2,4,5,7,9,12} is [0,1,2,5,7,12,45,67,70], with f70 the only non-minor
  profile (`generators(T_67)` agrees, EXACT);
- the non-minor profile's extra zeros are 3, 8, 10, 13, 14, 16, 18 (this equals my dangerous list,
  EXACT);
- "T ∪ {14} is an exact row/column/transpose image of {1,4,5,7,9,12,14}", proved in
  `root_inventory67_g14_uniform`. The closing paragraph (p. 4) repeats that the proof "uses the
  separately proved seven-contact boundary theorem directly".

Transport (EXACT). `find_g(S67, T_67 ∪ {14})` = ((0,2,1,3), (0,1,3,2), 0). This g sends 14 to 9 and
T_68 = {1,4,5,7,9,12} to {2,4,5,7,12,14}. So T_68 (stage 6) also covers T_67 ∪ {14} directly. Citing
T_68 would be an equivalent repair, and T_68 itself rests on S67.

The other ties have earlier covers that match the graph prerequisites: 3 (T_39), 8 (T_26), 10 (T_35),
13 (T_47), 16 (T_75), 18 (T_36).

Independence of S67 from T_67:
- PROOF-READ. None of `root_inventory67_g14_uniform`, `root_inventory67_whole_assembly` or
  `independent_inventory67_native_harmonic_conditional` mentions inventory 66, the set
  {2,4,5,7,9,12} or T_67. The g14 note (p. 1) sends its six zero-coefficient walls to "literal old42
  closed six-tuples". The conditional note and the assembly use "the fixed preceding 42-family list".
- EXACT. Each zero wall S67 ∪ {t}, t ∈ {18,17,15,11,6,0}, still has the non-minor generator f35. Each
  has a cover of stage at most 4 that contains t: 18 (T_8, T_9, T_15, T_36, T_40), 17 (T_29, T_31,
  T_35, T_47), 15 (T_5, T_10, T_12, T_26, ...), 11 (T_42, T_46), 6 (T_29, T_31, T_35, T_47),
  0 (T_5, T_10, ...). This agrees with S67's graph prerequisites T_5, T_8, T_35, T_46.
- EXACT, graph. T_67 is not in the transitive closure of S67 (15 nodes). S67 is at stage 5 and T_67
  at stage 8. Adding S67 -> T_67, and even T_68 -> T_67 as well, leaves the graph acyclic (cycle
  search on all edges of `proof-order.json` plus the added ones).

Conclusion: RESOLVED (issue: presentation). The edge S67 -> T_67 is used and missing from
`proof-order.json` and Table 6's dependency record. The transport is correct, S67's proof does not
depend on T_67, and the order stays acyclic, with the stages unchanged since 5 < 8.

## 4. Completeness of new-contact coverage for the cases with first-pass reports

Method (EXACT). For each case whose primary or own document has a report in `docs/` (53 cases), I
computed the dangerous new contacts and the earlier covers of each. I then classified the argument
type and compared with what the report records.

A. Pointwise or whole-cone proofs. The argument holds on the whole closed contact cone or on every
eight-distinct configuration, so new contacts are just points of that region and need no witness:
- T_1, T_2, T_3, T_4, T_5, T_6, T_7, T_15, T_18, T_19, T_20, T_21, T_22, T_23 (stage 1): no earlier
  cover exists for any of their dangerous contacts, which is consistent with stage 1. The reports
  record minor decompositions or SOC capacity on the eight-distinct layer, plus Appendix B, C or the
  six-point circle theorem.
- T_10, T_11: pointwise chart and capacity branches. The T_11 A-branch sends the zero wall b = 0 to
  T_7.
- T_25, T_26, T_27, T_28, T_29, T_31, T_34, T_35, T_36, T_37, T_39, T_40, T_42, T_43, T_44: radial
  necessity on the closed cone. Every report records that zero coefficients and walls are included
  (continuity plus a strict constant term).
- T_38, T_65: whole-cone proofs with some zero walls sent to earlier cases (item 2).
- T_47: only one component has a report. Its recorded transfers (historical T28, T35, T0, T38, T41,
  T6 and the face T ∪ {14} = S26) match my earlier-cover lists for the dangerous contacts
  {0, 9, 14, 18, 19, 20, 22}. The whole assembly is unreported.

B. Curve or harmonic arguments. These need an argument for each dangerous j:

| case (stage) | dangerous j (computed) | recorded in the report | missing |
|---|---|---|---|
| T_8 (4) | 4,10,16,18,19,22,23 | 4,10,18,19,22 earlier cases; 16, 23 by T^(7)_1, T^(7)_2 | none |
| T_9 (3) | 6,7,9,18,20,21,23 | 6,7,21,23 by T_4..T_7; 9 by T^(7)_1; 18, 20 by Thms 2-3 (T^(7)_3, T^(7)_4) | none |
| T_12 (3) | 6,7,11,13,14,15,16,17,18,19,21,23 | all 12 | none |
| T_13 (4) | 6,7,9,10,11,13,14,15,16,17,18,19,20,22 | all 14 | none |
| T_14 (5) | 6,7,9,13,14,16,17,20,21,22,23 | all 11 | none |
| T_16 (5) | frame {0,2,5,6,8,18}: 4,7,10,16,19,22,23 (EXACT in that frame) | all 7 | none |
| T_17 (4) | 6,7,11,14,15,17,18 | all 7 | none |
| T_24 (5) | 1,10,12,13,16,17,22 | all 7 (T_3, T_1, T_10, T^(7)_6, T_25 or T_26, T^(7)_5) | none |
| T_30 (9) | 0,3,4,7,8,10,15,16,18,19,21,23 | all 12 | none |
| T_32 (11) | 1,5,6,9,13,16,17,18,20,21,23 | all 11 | none |
| T_33 (5) | 1,2,5,9,10,13,14,16,17,18,20,21,23 | all 13 | none |
| T_41 (9) | 2,5,7,9,10,13,14,16,18,20,21,23 | all 12 | none |
| T_45 (7) | 0,6,9,18,20,21,23 | all 7 | none |
| T_50 (8) | 0,7,12,13,14,16,20 | all 7 | none |
| T_51 (7) | 0,6,7,13,14,16,20 | all 7 | none |
| T_52 (6) | 0,6,7,12,13,16,20 | all 7 | none |
| T_83 (9) | 0,1,7,8,19,20,22 | all 7 | none |
| T_84 (9) | 0,3,4,6,7,8,15,16,17,19,23 | all 11 | none |
| T_85 (10) | 1,3,4,6,8,11,12,15,17,19,20,22 | 1,4,6,11,12,15,17 (zeros of f35, Table 3) by T_40, T_55, T_81; 3,8,19,20,22 are zeros of f28 only, and the report records that f35 > 0 there forces d = 0 in the d > 0 branch, while the d < 0 branch is a minor sum | none |

Every case the reports name for (k, j) is in my exact earlier-cover list for (k, j).

Pairs (k, j) with no recorded argument: none.

Conclusion: RESOLVED (no issue). The only coverage that rests on unreplayed material:
- T_38 contacts 10, 15, 20 (chart certificates UNCHECKED);
- T_65 contact 15 in three unreported components;
- the stage-1 whole-cone cases T_39 and T_40, whose reports are could-not-check for their
  certificates. These are certificate questions, not missing arguments.

## 5. Other inconsistencies met

- T_9 ∪ {9} contains an image of T_8, which is stage 4, later than T_9. The document correctly cites
  the seven-contact face T^(7)_1 instead. That face's proof uses identities from the T_8 note
  (item 1). There is no cycle. Presentation.
- Our own first-pass report on T_8 lists T_13 (same stage 4) as an alternative witness for j = 19,
  and T_16 (stage 5) for j = 18. The valid earlier witnesses T_2 and T_11 are listed too. Not a proof
  issue, but those alternatives should be dropped from the report, since they are not earlier than
  T_8.
- T_65 table (2) uses historical labels. The paper-number reading fails for five of six entries
  (EXACT). Presentation, as elsewhere.
- T_67 could cite T_68 instead of S67 for tie 14 (item 3). Equivalent; no issue.
- No printed symmetry transport that I replayed (the 4 actions of T_38 p. 1, the S67 transport of
  T_67, the T^(7)_1 and T^(7)_2 identifications) fails to map the cited set into the needed set.
  Every witness case met in items 1-4 is earlier in Table 6 than the citing case.
