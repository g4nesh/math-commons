Evidence level: EXACT for all set, orbit, group-action, graph, stage and coverage statements
(own integer code `build_dependency_map.py`, output `dependency-map.json`); PROOF-READ only for the
reading of theorem statements, tie tables and dependency paragraphs named below; keyword tags are
UNCHECKED heuristics (keyword search only). No mathematical step of
any case proof is audited here.

# Dependency map of the ninety contact cases

Reproduce: `MOC4_SUPP=<supplement root> python -I -B build_dependency_map.py <text-extract dir> [outdir]`
(or pass the supplement root as first argument). The text-extract directory holds `moc4.txt` and
`proofs/<dir>__<name>.txt` made with `pdftotext -layout`. Run time about 10 s of CPU.

## Deviations and limits

- Citation extraction from PDF text is regex-based (contact sets in braces, historical labels
  "T NN" near keywords, "inventoryNN", core letters, document file names, "T ∪ {k}" forms). It
  misses citations phrased only in words and produces some noise; every item in the
  missing-from-graph list below was read in context. Edges without textual support are reported,
  not treated as errors.
- "Edge support by boundary" means: some image under G of the source contact set lies in
  target ∪ {t} for a tie t at which a non-minor generator of the target cone also vanishes
  (a "relevant tie", computed from the 82 generators). This is the structural reason an earlier
  case can be used on a boundary face; whether the proof actually uses it there is read only
  where the text prints a tie table.
- Document homes (which document proves which node) come from file-name rules plus statement
  sets; joint documents are attached to every node they state. Pages are compiled-PDF pages
  counted from form feeds in the extracts.
- Nothing from the supplement was executed. Their checker scripts and JSON results were not
  used as evidence; where a witness exists only in checker data, the edge is listed as
  "no text support".
- Task 4 tags are coarse: a tag means the keyword family appears in the case's own or
  catalogued documents, not that the input is used or needed.

## Counts (EXACT)

- Cases 90, auxiliary seven-sets 27, branch nodes 4, graph nodes 121, distinct edges 398
  (`cases/proof-order.json`). Documents 208, 779 pages, all with a text extract.
- 82 generators (34 minor, 48 non-minor) recomputed; 214 orbits of six-sets, 124 minor-only;
  the 90 others are exactly the 90 Appendix E sets (one per orbit). The 27 auxiliary sets lie in
  27 distinct orbits and equal Table 5.
- `cases/index.tsv`: literal sets equal Appendix E; canonical representative and printed
  action correct for all 90; method groups equal Appendix G.

## Task 1. Graph checks and three numbering systems

- Acyclic: no strongly connected component of size > 1 (Tarjan). Every edge goes from a
  lower to a higher stage; stage = longest-path layer for every node.
- Table 6 versus graph: 0 discrepancies; every graph node appears in Table 6 and every Table 6
  entry is a graph node.
- Numbering: each case is identified by its Appendix E contact set. Inventory index = paper
  index - 1 for all 90; historical proof labels are paper index - 1 except
  T0 -> T_1, T1 -> T_2, T2 -> T_6, T3 -> T_3, T4 -> T_4, T5 -> T_5, the index.tsv labels
  I29 -> T_30, I47 -> T_48, I83 -> T_84, and the label "T47", which in
  inventory-52 and later notes denotes T_65 (the notes say so). For 19 printed tie -> witness
  tables (118 entries) the historical reading covers the face in all 118 entries, and the
  paper-number reading fails in most, so the historical reading is the intended one (EXACT).
- Every non-branch edge (392 of 398) is a one-tie boundary witness: an image of the source
  lies in target ∪ {t} for a relevant tie t of the target. The 6 branch edges are
  T_10^(1,2) -> T_10, T_11^(1,2) -> T_11, and T_7 -> T_11^(1,2).
- Edges with no mention in the PDF text of the target's documents: 254 of 398. All are
  structurally valid (previous item); their use is recorded only in unexecuted checker data
  and the graph file. Full list: `edges_without_text_support`.
- Printed witness tables: 118 entries, 118 cover their face, 117 are earlier in Table 6. The
  exception is in the T_49 allocation notes (below).

### Text citations missing from proof-order.json (47 extracted; reviewed)

Substantive (a node is used in a proof but has no graph edge):

1. T_67 <- T^(7)_27 (S67). `root_inventory66_native_whole` p.1 uses the proved seven-set
   {1,4,5,7,9,12,14} for tie 14 of T_67. No graph prerequisite of T_67 covers tie 14 (EXACT).
   Order is fine (stage 5 < 8): missing edge, not a cycle.
2. Within-stage lemma-level reuse (stage 1), contradicting the paper's statement that proofs
   inside a stage are independent, at the level of pointwise capacity lemmas, not contact
   theorems: the "T24 core" coherent capacity proved in T_25's document is used by T_28, T_37,
   T_39, T_42, T_44 (all stage 1) and T_65 (stage 6, fine); T_39 also uses the "T27 allocation
   lemma" (T_28), the "T28 resource lemma" (T_29); T_34's bad70 note cites a "T34" lemma
   (T_35). These form an acyclic order (T_25 -> T_28 -> T_39, T_29 -> T_39, T_35 -> T_34);
   no cycle results, but the stage-1 independence claim is false at document level.
3. Same-stage use in stage 6: `independent_inventory77_native_harmonic_uniform` (T_78) p.3
   reuses an "earlier inventory75 nonreal quadratic-coordinate" argument (T_76, also stage 6).
   Lemma-level; no cycle (T_76 does not cite T_78).
4. T^(7)_1 (`root_moc4_T7_boundary16` p.1) says it uses "the T7 source chart and spectral
   identities proved by the adjacent original" note, i.e. identities from T_8's document,
   while T_8 depends on T^(7)_1. Reverse direction at lemma level; not a theorem cycle if those
   identities are proved without the boundary note. UNCHECKED: needs a reader of the T_8 note.
5. T_65's low_j_high_e note uses Corollary 3 of `final_referee_T46_six_ray_coherent_capacity`
   (a T_47 document; stage 3 < 6, fine). T_53's conditional note reuses coordinates from the
   T_49 harmonic theorem (stage 4 < 8, fine). coreL and coreH cite primitive inequalities
   "recorded in the coreD package" (common pointwise atoms; all stage 1).
6. Declared-input lists: T_17's note declares hist T15 = T_16 (stage 5, after T_17 at stage 4);
   the tie it would serve (11) is covered by the earlier graph prerequisite T_15, so the
   declaration is stale, not a cycle. T_8's note (stage 4) declares T0..T6, T9, T10, i.e. up to
   T_11 (stage 3): fine. T_12 and T_13 declare "T0..T10, T14" (T_15 is stage 1): fine.
7. T_49 allocation notes I-IV (catalogued branches of T_49, stage 4) send a zero face to
   "T23" = T_24 (stage 5). The T_49 assembly says it does not use these notes. Stage inversion
   in an unused document only.

Not dependencies (reviewed): REVERSE hits where a seven-set or branch document names its own
base six-set (T^(7)_6, T^(7)_19, T^(7)_23, T^(7)_26, T^(7)_27, T_10^(1,2), T_11^(1,2)); the
mutual mention between the two T_11 branch notes ("the same chart as" the other); witness
lists already in the transitive closure (T_53, T_48, T_66, T_45, T_60, T_49); T_66's printed
witnesses T_10, T_24 are not graph edges, but graph prerequisites T_28, T_16 cover the same
ties (EXACT).

Unresolved references (8): seven-sets outside the 27 auxiliary ones, all proved inside the
citing case's own documents: T_27 uses T_27+{10}, +{13}, +{16} (`final_referee_T26_boundary10`,
`_boundary16`), and T_16's working-frame set {4,7,10,16,19,22,23} is a tie list, not a set.

## Task 3. Strongly connected components, dangling citations, homes

- Conditional hypotheses as edges (hypothesis node -> user case): all eleven (ten theorems,
  C3 has two) are already graph edges, and in no case does the hypothesis node's transitive
  closure contain the user (EXACT). Adding the text citations of items 1, 2, 3, 5, 6 and the
  allocation-note edge of item 7 (16 edges) and recomputing Tarjan SCCs: none of size > 1.
  Adding item 4 as a theorem edge T_8 -> T^(7)_1 would close the 2-cycle {T_8, T^(7)_1}; it is
  a lemma-level reuse of identities, flagged for the T_8 auditor (UNCHECKED whether those
  identities depend on the boundary note).
- Dangling: none. Every auxiliary and branch node has at least one proving document; every
  edge source is a node.
- Relevant ties with no graph witness (58 nodes) are mostly stage-1/2 whole-cone proofs and
  the cores, which treat those ties inside their own argument. Exceptions worth a reader:
  T_53 ties 5,7,11,14,15,18 (covered only transitively through T^(7)_19, as the text does);
  T_67 tie 14 (item 1); T_38 ties 10,15,20 and T_65 tie 15 (not examined).
- Every core document is tied to its contact set by the theorem statement set (EXACT orbit
  match): coreA..coreO and the S-sets, see `core_documents` in the JSON and the table below.

## Auxiliary seven-sets and branch statements: proving documents

Stage = Table 6. "Used by" = graph successors. T^(7)_3 and T^(7)_4 are proved as Theorems 3 and 2
of `independent_moc4_T8_uniform` (pp. 6 and 4); the polar-reduction note only mentions them.

| Node | Contacts | Stage | Proving documents (page, pages) | Direct prerequisites | Used by |
|---|---|---|---|---|---|
| T^(7)_1[S_0,1,2,3,6,15,17] | {0,1,2,3,6,15,17} | 1 | root_moc4_T7_boundary16_20261022 p.1 (1 pp) | - | T_8, T_9 |
| T^(7)_10[coreG] | {0,1,2,7,9,10,23} | 1 | independent_shared_coreG_uniform_20261023 p.1 (5 pp) | - | T_50, T_52 |
| T^(7)_11[coreC] | {0,1,2,7,9,16,17} | 1 | blind_moc4_coreC_uniform_20261023 p.1 (5 pp) | - | T_51, T_52 |
| T^(7)_12[coreD] | {0,1,2,7,9,16,23} | 1 | blind_moc4_coreD_low_pair_20261023 p.1 (3 pp); blind_moc4_coreD_uniform_20261023 p.1 (4 pp) | - | T_51 |
| T^(7)_13[coreM] | {0,1,2,7,10,16,17} | 1 | independent_shared_coreM_uniform_20261023 p.1 (5 pp) | - | T_55, T_63, T_80 |
| T^(7)_14[coreL] | {0,1,2,7,10,16,23} | 1 | root_shared_coreL_uniform_20261023 p.1 (4 pp) | - | T_55, T_63 |
| T^(7)_15[coreO] | {0,1,2,7,14,16,23} | 1 | root_shared_coreO_uniform_20261023 p.1 (3 pp) | - | T_61, T_64 |
| T^(7)_16[coreK] | {0,1,2,9,10,14,17} | 1 | final_referee_shared_coreK_uniform_20261023 p.1 (3 pp) | - | T_74 |
| T^(7)_17[coreE] | {0,1,2,9,10,16,17} | 1 | independent_coreE_four_chain_capacity_20261023 p.1 (2 pp); independent_coreE_uniform_20261023 p.1 (5 pp) | - | T_73 |
| T^(7)_2[S_0,1,2,3,6,15,22] | {0,1,2,3,6,15,22} | 1 | root_moc4_T7_boundary23_20261022 p.1 (1 pp) | - | T_8 |
| T^(7)_21[coreH] | {0,1,6,9,16,17,22} | 1 | root_shared_coreH_uniform_20261023 p.1 (3 pp) | - | T_81, T_83 |
| T^(7)_22[coreJ] | {0,2,3,7,11,12,23} | 1 | blind_shared_coreJ_uniform_20261023 p.1 (5 pp) | - | T_32, T_33 |
| T^(7)_24[coreN] | {0,3,7,11,12,16,20} | 1 | blind_shared_coreN_uniform_20261023 p.1 (4 pp) | - | T_89, T_90 |
| T^(7)_5[S_0,1,2,6,9,11,20] | {0,1,2,6,9,11,20} | 1 | independent_moc4_T23_g17_uniform_20261022 p.1 (4 pp) | - | T_24 |
| T^(7)_7[coreF] | {0,1,2,7,9,10,14} | 1 | blind_shared_coreF_uniform_20261023 p.1 (4 pp) | - | T_50 |
| T^(7)_8[coreB] | {0,1,2,7,9,10,16} | 1 | blind_moc4_coreB_uniform_20261023 p.1 (4 pp) | - | T_50, T_51 |
| T^(7)_9[coreA] | {0,1,2,7,9,10,17} | 1 | root_shared_seven_A_uniform_20261023 p.1 (4 pp) | - | T_50, T_52 |
| T_10^(1) | branch | 1 | final_referee_T9_first_branch_20261022 p.1 (4 pp) | - | T_10 |
| T_10^(2) | branch | 1 | final_referee_T9_second_branch_20261022 p.1 (3 pp) | - | T_10 |
| T^(7)_23[S48c0] | {0,3,5,6,9,12,17} | 2 | blind_moc4_T48_c0_all_j_chain_20261023 p.1 (2 pp); blind_moc4_T48_c0_native_interlacing_20261023 p.1 (3 pp); blind_moc4_T48_c0_upper_Heron_chain_20261023 p.1 (2 pp); independent_true_T48_c0_native_Heron_20261023 p.1 (3 pp); independent_true_T48_c0_native_Heron_20261023 p.1 (3 pp) | T_1, T_2, T_25, T_34, T_44 | T_49 |
| T^(7)_26[S46a0] | {1,3,5,6,7,8,14} | 2 | final_referee_T46_a0_uniform_20261023 p.1 (3 pp) | T_1, T_29, T_36, T_39, T_42, T_7 | T_47 |
| T^(7)_3[S_0,1,2,3,6,17,19] | {0,1,2,3,6,17,19} | 2 | independent_moc4_T8_polar_reduction_20261022 p.1 (4 pp); independent_moc4_T8_uniform_20261022 p.6 (8 pp) | T_4 | T_9 |
| T^(7)_4[S_0,1,2,3,6,17,22] | {0,1,2,3,6,17,22} | 2 | independent_moc4_T8_polar_reduction_20261022 p.1 (4 pp); independent_moc4_T8_uniform_20261022 p.4 (8 pp) | T_4 | T_9 |
| T_11^(1) | branch | 2 | final_referee_T10_A_branch_20261022 p.1 (4 pp) | T_7 | T_11 |
| T_11^(2) | branch | 2 | root_moc4_T10_B_branch_20261022 p.1 (4 pp) | T_7 | T_11 |
| T^(7)_18[S75] | {0,1,3,10,11,12,16} | 4 | final_referee_inventory75_boundary16_uniform_20261023 p.1 (3 pp) | T_11, T_29, T_31, T_40, T_44 | T_76, T_78 |
| T^(7)_6[S_0,1,2,6,9,11,23] | {0,1,2,6,9,11,23} | 4 | independent_moc4_T23_g13_uniform_20261022 p.1 (4 pp) | T_11, T_12, T_15, T_3 | T_24 |
| T^(7)_20[S56] | {0,1,5,7,9,12,17} | 5 | final_referee_inventory56_boundary17_large_product_20261023 p.1 (3 pp); final_referee_inventory56_boundary17_uniform_20261023 p.1 (3 pp) | T_17, T_26, T_31, T_34, T_43, T_44 | T_57, T_58, T_78 |
| T^(7)_27[S67] | {1,4,5,7,9,12,14} | 5 | root_inventory67_g14_uniform_20261023 p.1 (3 pp) | T_35, T_46, T_5, T_8 | T_68, T_71, T_76 |
| T^(7)_19[S52] | {0,1,4,6,9,12,17} | 7 | blind_moc4_inventory52_plus17_fourth_products_20261023 p.1 (3 pp); blind_moc4_inventory52_plus17_low_products_20261023 p.1 (3 pp); blind_moc4_inventory52_plus17_middle_products_20261023 p.1 (2 pp); blind_moc4_inventory52_plus17_third_products_20261023 p.1 (3 pp) | T_18, T_26, T_4, T_65 | T_53 |
| T^(7)_25[S78] | {1,3,4,8,11,12,15} | 7 | final_referee_inventory78_boundary15_uniform_20261023 p.1 (4 pp) | T_13, T_21, T_40, T_70, T_71, T_76 | T_79 |

## Conditional theorems (details in conditional-tracing.md)


| Id | Case | Hypothesis (paper name) | Discharging document, page | Status |
|---|---|---|---|---|
| C1 | T_51 | D = T^(7)_12 (coreD) | blind_moc4_coreD_uniform p.1 | DISCHARGED |
| C2 | T_52 | G = T^(7)_10 (coreG) | independent_shared_coreG_uniform p.1 | DISCHARGED |
| C3 | T_50 | F = T^(7)_7, G = T^(7)_10 | blind_shared_coreF_uniform p.1; coreG p.1 | DISCHARGED (assembly does not name the documents) |
| C4 | T_53 | T^(7)_19 (S52) | blind_moc4_inventory52_whole_assembly pp.1-3 | DISCHARGED |
| C5 | T_57 | T^(7)_20 (S56) | final_referee_inventory56_boundary17_uniform p.1 | DISCHARGED |
| C6 | T_58 | T^(7)_20 (S56) | same | DISCHARGED |
| C7 | T_68 | T^(7)_27 (S67) | root_inventory67_g14_uniform p.1 | DISCHARGED |
| C8 | T_76 | T^(7)_18 (S75) | final_referee_inventory75_boundary16_uniform p.1 | DISCHARGED |
| C9 | T_49 | T^(7)_23 (S48c0) | root_true_T48_whole_assembly p.3 | DISCHARGED |
| C10 | T_24 | T+16 (via T_25), T+17 in orbit of T^(7)_5 | root_moc4_T24_uniform p.1 (case T_25); independent_moc4_T23_g17_uniform p.1 | DISCHARGED |

No conditional theorem is DISCHARGED-BUT-CIRCULAR, PARTIALLY DISCHARGED or OPEN at the level of
statements and dependency structure. What remains open is the mathematics of each discharging
proof, which this map does not audit.

## Task 4. Common-input tags (UNCHECKED heuristic)

Counts over the 90 cases: two real responses / Theorem 1.2 (2R) 66, matching cubic / Sing F (MC)
65, Lemma 3.5 harmonic / maximum principle (H) 61, Appendix B shared zero / cross collision (Bs)
53, Rouché / wall perturbation (W) 52, Gram / PSD / SOS / capacity (G) 48, Appendix C Hessian rank
(Ch) 42, Appendix B coincident eigenvalues (Bc) 35, Appendix C one circle (Co) 32, resultant /
discriminant (R) 30, Birkhoff / height kernel (Bk) 25, single-collision compensation (B1) 8.
Common documents in the catalogue: `alternate_onecircle_20261022_elementary` (81 cases),
`principal_minors` (63), `root_unitary_two_orthogonals_20261022` (50),
`final_referee_shared_zero_birkhoff_20261022` (46), `one_circle` (39),
`final_referee_bounded_bad_harmonic_lemma_20261023` (33), `root_complete_negative_mass_20261022`
(26), `root_moc4_global_real_Hessian_inertia_20261022` (24). Known over-tagging: T_1 and T_2 have
the long joint note `independent_moc4_six_contact_reduction` as primary, so they carry nearly every tag. No
case-by-case reading of whether a tagged input is actually used was done (UNCHECKED); the tags
only point case auditors to where to look.

## Per-case table

Contacts = Appendix E literal set. Stage = Table 6 (equal to graph stage for all 90). Own docs =
documents homed to the case by file name; catalog branches = documents listed for the case in
`cases/catalog.json` (pages total). Graph prerequisites use paper numbering; S7.k = T^(7)_k.
Tags as in Task 4.

| Case | Contacts (App. E) | Grp | Stage | Primary document, page | Own docs / catalog branches (pages) | Graph prerequisites (direct) | Conditional input | Tags |
|---|---|---|---|---|---|---|---|---|
| T1 | {0,1,2,3,6,7} | I | 1 | independent_moc4_six_contact_reduction_20261022, p.6 | 1 / 0 (0 pp) | - | - | Bc Bs Bk Ch Co G H 2R MC |
| T2 | {0,1,2,3,6,8} | II | 1 | independent_moc4_six_contact_reduction_20261022, p.7 | 1 / 0 (0 pp) | - | - | Bc Bs Bk Ch Co G H 2R MC |
| T3 | {0,1,2,3,7,8} | III | 1 | independent_moc4_T3_uniform_20261022, p.1 | 1 / 14 (71 pp) | - | - | Bc Bs G 2R |
| T4 | {0,1,4,5,7,8} | III | 1 | root_moc4_T4_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Bc Bs G 2R MC |
| T5 | {1,2,3,4,6,8} | III | 1 | final_referee_moc4_T5_uniform_20261022, p.1 | 2 / 1 (3 pp) | - | - | Bc G 2R MC |
| T6 | {0,1,2,4,6,8} | III | 1 | independent_moc4_T2_uniform_20261022, p.1 | 2 / 1 (2 pp) | - | - | Bc Bs G 2R |
| T7 | {0,1,2,4,7,8} | II | 1 | root_moc4_T6_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Bs |
| T8 | {0,2,5,6,7,8} | VI | 4 | root_moc4_T7_uniform_20261022, p.1 | 1 / 2 (2 pp) | S7.1, S7.2, T2, T3, T6, T11 | - | Bs G H R |
| T9 | {0,1,2,3,11,12} | VI | 3 | independent_moc4_T8_uniform_20261022, p.1 | 3 / 2 (7 pp) | S7.1, T4, T5, T6, S7.3, S7.4 | - | Bc Bs G H 2R MC |
| T10 | {0,1,3,5,7,8} | III | 2 | final_referee_moc4_T9_uniform_20261022, p.1 | 1 / 2 (7 pp) | T10^(1), T10^(2) | - | Bc G 2R MC |
| T11 | {0,1,3,5,6,8} | III | 3 | final_referee_moc4_T10_uniform_20261022, p.1 | 1 / 2 (8 pp) | T11^(1), T11^(2) | - | G |
| T12 | {0,2,3,5,9,12} | VI | 3 | final_referee_moc4_T11_uniform_20261022, p.1 | 1 / 0 (0 pp) | T1, T2, T3, T4, T5, T6, T7, T10 | - | Bc Bs Co G H 2R |
| T13 | {1,2,3,4,8,12} | VI | 4 | final_referee_moc4_T12_uniform_20261022, p.1 | 1 / 0 (0 pp) | T1, T2, T3, T4, T5, T6, T7, T10, T11 | - | Co G H 2R |
| T14 | {0,2,3,5,11,12} | VI | 5 | independent_moc4_T13_uniform_20261022, p.1 | 1 / 0 (0 pp) | T1, T2, T4, T6, T7, T11, T9, T8 | - | Bc Bs H 2R MC |
| T15 | {1,2,3,4,6,7} | III | 1 | final_referee_moc4_T14_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Bc Bs Co G |
| T16 | {0,2,5,6,9,12} | VI | 5 | root_moc4_T15_uniform_20261022, p.1 | 1 / 0 (0 pp) | T15, T2, T3, T8 | - | Co H 2R |
| T17 | {0,1,4,5,9,12} | VI | 4 | root_moc4_T16_uniform_20261022, p.1 | 1 / 0 (0 pp) | T15, T4, T5, T9 | - | Bs H |
| T18 | {0,1,3,6,7,8} | II | 1 | root_moc4_T17_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Bs Co H 2R |
| T19 | {0,1,4,6,7,8} | II | 1 | root_moc4_T18_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Bs Co H 2R |
| T20 | {0,1,5,6,7,8} | II | 1 | root_moc4_T19_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | 2R |
| T21 | {0,3,5,9,11,12} | IV | 1 | root_moc4_T20_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | H 2R MC |
| T22 | {2,3,4,8,10,12} | IV | 1 | root_moc4_T21_T22_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | H 2R MC |
| T23 | {2,3,5,9,11,12} | IV | 1 | root_moc4_T21_T22_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | H 2R MC |
| T24 | {2,3,4,6,7,8} | VI | 5 | independent_moc4_T23_uniform_20261022, p.1 | 2 / 8 (30 pp) | S7.5, T15, T20, T25, T12, S7.6, T13 | C10: hyp T+17 = T^(7)_5 -> DISCHARGED | Bs Ch Co G H 2R |
| T25 | {3,4,7,8,10,12} | V | 1 | root_moc4_T24_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G MC |
| T26 | {2,3,5,7,8,12} | V | 2 | final_referee_T25_uniform_20261022, p.1 | 1 / 0 (0 pp) | T25 | - | Ch G MC |
| T27 | {2,3,4,7,8,12} | V | 1 | final_referee_T26_uniform_20261022, p.1 | 3 / 2 (5 pp) | - | - | Ch G MC |
| T28 | {2,3,5,6,9,12} | V | 1 | blind_moc4_T27_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G 2R MC |
| T29 | {0,4,5,6,7,8} | V | 1 | blind_moc4_T28_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G MC |
| T30 | {2,5,9,12,13,14} | VI | 9 | final_referee_inventory29_native_uniform_20261023, p.1 | 1 / 8 (44 pp) | T22, T29, T4, T40, T7, T9, T71, T79 | - | Bs Bk Ch H W 2R MC |
| T31 | {0,1,5,6,9,12} | V | 1 | root_moc4_T30_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G 2R MC |
| T32 | {0,2,3,7,11,12} | VI | 11 | blind_inventory31_native_whole_20261023, p.1 | 1 / 7 (31 pp) | S7.22, T21, T29, T31, T34, T9, T14, T78, T79, T77, T82 | - | Bc Bs Bk Co H W MC R |
| T33 | {0,3,6,7,11,12} | VI | 5 | blind_inventory32_native_whole_20261023, p.1 | 1 / 7 (31 pp) | S7.22, T19, T22, T29, T31, T4, T40, T7, T13 | - | Bc Bs Bk Co H W 2R MC R |
| T34 | {2,3,5,9,10,12} | V | 1 | blind_moc4_bad70_full_safe_disjunction_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G MC |
| T35 | {2,3,4,6,9,12} | V | 1 | final_referee_T34_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G W MC |
| T36 | {2,3,4,6,11,12} | V | 1 | final_referee_T35_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G W MC |
| T37 | {0,5,6,7,9,12} | V | 1 | independent_moc4_T36_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G 2R MC |
| T38 | {1,2,3,8,11,12} | V | 5 | root_moc4_T37_uniform_20261022, p.1 | 1 / 10 (33 pp) | T4, T46, T11, T8 | - | Ch G MC |
| T39 | {2,3,5,7,9,12} | V | 1 | blind_moc4_T38_uniform_20261022, p.1 | 1 / 1 (5 pp) | - | - | Ch G MC |
| T40 | {0,1,5,9,14,18} | V | 1 | independent_moc4_T39_uniform_20261022, p.1 | 4 / 5 (25 pp) | - | - | Ch G W 2R MC R |
| T41 | {0,1,3,6,11,12} | VI | 9 | independent_inventory40_native_uniform_20261023, p.1 | 1 / 9 (47 pp) | T19, T37, T40, T43, T5, T9, T57, T58, T71, T76, T79 | - | Bs B1 H W 2R |
| T42 | {1,2,3,6,7,8} | V | 1 | final_referee_T41_uniform_20261022, p.1 | 5 / 5 (19 pp) | - | - | Ch G H W 2R MC R |
| T43 | {1,2,4,6,7,8} | V | 1 | blind_moc4_T42_uniform_20261022, p.1 | 1 / 1 (5 pp) | - | - | Ch G MC |
| T44 | {1,4,5,6,7,8} | V | 1 | independent_moc4_T43_uniform_20261022, p.1 | 1 / 0 (0 pp) | - | - | Ch G 2R MC |
| T45 | {1,2,3,7,11,12} | VI | 7 | independent_moc4_T44_native_harmonic_uniform_20261023, p.1 | 7 / 9 (34 pp) | T43, T46, T9, T8, T38, T65 | - | Bc Bs Bk Ch Co G H W 2R MC R |
| T46 | {3,4,5,9,10,12} | V | 2 | independent_moc4_T45_uniform_20261022, p.1 | 1 / 2 (6 pp) | T15, T34, T7 | - | Ch G 2R MC |
| T47 | {1,3,5,6,7,8} | V | 3 | root_moc4_T46_whole_assembly_20261023, p.1 | 10 / 35 (101 pp) | T1, T29, T36, T39, T42, T7, S7.26 | - | Ch G W MC R |
| T48 | {0,2,5,9,12,13} | VI | 8 | root_true_inventory47_native_uniform_20261023, p.1 | 1 / 28 (93 pp) | T25, T36, T37, T4, T46, T12, T47, T49, T16, T70, T71, T69 | - | Bk G H W 2R MC R |
| T49 | {0,3,5,6,9,12} | VI | 4 | root_true_T48_whole_assembly_20261023, p.1 | 8 / 15 (54 pp) | T21, T34, T37, T42, T44, S7.23, T12 | C9: hyp S48c0 = T^(7)_23 -> DISCHARGED | Bc Bs Bk Ch Co G H W 2R MC R |
| T50 | {1,3,5,6,10,11} | VI | 8 | blind_inventory49_whole_20261023, p.1 | 2 / 21 (114 pp) | S7.10, S7.7, S7.8, S7.9, T43, T10, T66 | C3: hyp F = T^(7)_7 -> DISCHARGED; C3: hyp G = T^(7)_10 -> DISCHARGED | Bc Bs Bk Co H W 2R MC |
| T51 | {1,3,5,10,11,12} | VI | 7 | root_inventory50_bounded_native_whole_20261023, p.1 | 2 / 23 (114 pp) | S7.11, S7.12, S7.8, T65, T76 | C1: hyp D = T^(7)_12 -> DISCHARGED | Bs G H W 2R |
| T52 | {1,3,5,10,11,14} | VI | 6 | root_inventory51_bounded_native_whole_20261023, p.1 | 2 / 26 (126 pp) | S7.10, S7.11, S7.9, T39, T46, T24 | C2: hyp G = T^(7)_10 -> DISCHARGED | Bs H W 2R MC R |
| T53 | {0,1,4,6,9,12} | VI | 8 | blind_moc4_inventory52_whole_assembly_20261023, p.1 | 2 / 7 (21 pp) | S7.19 | C4: hyp S52 = T^(7)_19 -> DISCHARGED | Bk Ch Co G H W 2R MC R |
| T54 | {2,3,5,6,11,12} | VI | 9 | root_inventory53_native_uniform_20261023, p.1 | 1 / 10 (51 pp) | T18, T2, T28, T35, T42, T11, T47, T52, T51, T50, T53 | - | Bs Co H W 2R MC R |
| T55 | {1,2,3,6,11,12} | VI | 9 | independent_inventory54_native_uniform_20261023, p.1 | 1 / 30 (139 pp) | S7.13, S7.14, T9, T38, T45, T50, T73 | - | Bs B1 H W 2R MC |
| T56 | {3,4,5,8,9,12} | VI | 10 | independent_inventory55_native_uniform_20261023, p.1 | 1 / 10 (51 pp) | T19, T37, T43, T44, T6, T50, T41, T55 | - | Bc Bs B1 H W 2R |
| T57 | {0,1,5,7,9,12} | VI | 6 | final_referee_inventory56_57_whole_assembly_20261023, p.1 | 3 / 10 (52 pp) | T31, T39, T43, T44, T26, T17, S7.20 | C5: hyp S56 = T^(7)_20 -> DISCHARGED | Bc Bs Bk Ch G H W 2R MC R |
| T58 | {0,3,5,7,9,12} | VI | 6 | final_referee_inventory56_57_whole_assembly_20261023, p.1 | 2 / 10 (52 pp) | T19, T29, T34, T43, T5, S7.20 | C6: hyp S56 = T^(7)_20 -> DISCHARGED | Bc Bs Bk Ch Co G H W 2R MC R |
| T59 | {3,4,5,6,9,12} | VI | 8 | independent_inventory58_native_harmonic_uniform_20261023, p.1 | 1 / 4 (31 pp) | T28, T34, T36, T47, T49, T57, T69 | - | Bc Bs Bk Co H W 2R MC R |
| T60 | {2,3,5,7,11,12} | VI | 9 | blind_moc4_inventory59_native_harmonic_uniform_20261023, p.1 | 1 / 8 (24 pp) | T20, T23, T36, T4, T47, T53 | - | Bs Bk H W 2R MC R |
| T61 | {1,2,3,9,11,12} | VI | 9 | root_inventory60_native_uniform_20261023, p.1 | 1 / 20 (99 pp) | S7.15, T15, T8, T50, T62 | - | Bs B1 Co H W 2R R |
| T62 | {3,5,7,8,9,12} | VI | 8 | final_referee_inventory61_native_uniform_20261023, p.1 | 1 / 11 (62 pp) | T26, T8, T52, T51, T75 | - | Bs H W 2R |
| T63 | {1,2,3,10,11,12} | VI | 9 | independent_inventory62_native_uniform_20261023, p.1 | 1 / 30 (139 pp) | S7.13, S7.14, T38, T71, T51, T81, T73 | - | Bs B1 H W 2R MC |
| T64 | {0,3,5,9,14,18} | VI | 10 | independent_inventory63_native_uniform_20261023, p.1 | 1 / 17 (89 pp) | S7.15, T36, T14, T54, T55 | - | Bs Co H W 2R R |
| T65 | {0,4,5,6,9,12} | V | 6 | root_moc4_T47_whole_assembly_20261023, p.1 | 22 / 31 (100 pp) | T4, T42, T44, T8, T38 | - | Ch G W 2R MC R |
| T66 | {2,4,5,6,9,12} | VI | 7 | root_moc4_T65_native_harmonic_uniform_20261023, p.1 | 1 / 4 (23 pp) | T25, T28, T46, T16, T65 | - | Bc Bs Bk Ch Co H W 2R MC |
| T67 | {2,4,5,7,9,12} | VI | 8 | root_inventory66_native_whole_20261023, p.1 | 1 / 7 (40 pp) | T35, T36, T39, T26, T47, T75 | - | Bc Bs Bk Co H W 2R MC R |
| T68 | {1,4,5,7,9,12} | VI | 6 | root_inventory67_whole_assembly_20261023, p.1 | 2 / 6 (38 pp) | T29, T35, T42, T26, T17, S7.27, T38 | C7: hyp S67 = T^(7)_27 -> DISCHARGED | Bc Bs Bk Ch Co G H W 2R MC R |
| T69 | {0,2,5,7,9,12} | VI | 7 | final_referee_inventory68_native_harmonic_20261023, p.1 | 1 / 4 (23 pp) | T44, T46, T12, T16, T38, T65 | - | Bc Bs Bk Co G H W 2R MC R |
| T70 | {1,3,4,8,9,12} | VI | 5 | final_referee_inventory69_native_harmonic_20261023, p.1 | 1 / 4 (23 pp) | T29, T4, T42, T44, T13 | - | Bc Bs Bk Co H W 2R MC |
| T71 | {1,3,4,10,11,12} | VI | 6 | independent_inventory70_native_harmonic_uniform_20261023, p.1 | 1 / 5 (34 pp) | T29, T31, T37, T40, T42, T5, S7.27 | - | Bc Bk Co H W MC R |
| T72 | {1,3,4,7,8,12} | V | 1 | final_referee_inventory71_uniform_20261023, p.1 | 1 / 28 (132 pp) | - | - | Ch G W 2R MC |
| T73 | {0,4,5,7,9,12} | VI | 8 | root_inventory72_native_uniform_20261023, p.1 | 1 / 29 (134 pp) | S7.17, T37, T17, T52, T51, T69 | - | Bs B1 Co H W 2R MC |
| T74 | {3,4,5,7,8,12} | VI | 9 | final_referee_inventory73_native_uniform_20261023, p.1 | 1 / 13 (69 pp) | S7.16, T25, T27, T62, T73 | - | Bs Ch G H W 2R MC |
| T75 | {3,5,7,8,10,12} | VI | 7 | independent_inventory74_native_harmonic_uniform_20261023, p.1 | 1 / 4 (31 pp) | T25, T28, T26, T17, T65 | - | Bc Bs Ch Co H W 2R R |
| T76 | {0,1,3,10,11,12} | VI | 6 | final_referee_inventory75_whole_assembly_20261023, p.1 | 2 / 7 (41 pp) | T31, T40, T44, T46, T11, S7.18, S7.27 | C8: hyp S75 = T^(7)_18 -> DISCHARGED | Bc Bs Bk Ch Co G H W 2R R |
| T77 | {0,1,3,10,11,16} | VI | 9 | root_inventory76_bounded_masks_uniform_20261023, p.1 | 1 / 7 (40 pp) | T36, T11, T47, T70, T76, T78, T59, T67 | - | G H W 2R MC R |
| T78 | {1,4,5,6,9,12} | VI | 6 | independent_inventory77_native_harmonic_uniform_20261023, p.1 | 1 / 9 (46 pp) | T29, T34, T36, T11, S7.18, S7.20 | - | Bc Bs Bk Co H W 2R MC R |
| T79 | {1,3,4,8,11,12} | VI | 8 | final_referee_inventory78_native_uniform_20261023, p.1 | 1 / 7 (39 pp) | T40, T13, T70, T58, T71, T78, S7.25 | - | Bk Ch H W 2R MC R |
| T80 | {1,3,4,6,8,12} | VI | 9 | final_referee_inventory79_native_uniform_20261023, p.1 | 1 / 28 (132 pp) | S7.13, T44, T5, T69, T50, T73 | - | Bc Bs Ch Co H W 2R MC |
| T81 | {1,2,10,11,12,13} | VI | 7 | root_inventory80_native_uniform_20261023, p.1 | 1 / 28 (130 pp) | S7.21, T34, T37, T46, T76 | - | Bs B1 Co H W 2R MC |
| T82 | {2,5,9,10,12,13} | VI | 10 | final_referee_inventory81_native_uniform_20261023, p.1 | 1 / 4 (21 pp) | T31, T34, T79, T60 | - | Bs Bk Ch G H W 2R |
| T83 | {3,5,6,9,14,18} | VI | 9 | final_referee_inventory82_native_uniform_20261023, p.1 | 1 / 28 (130 pp) | S7.21, T40, T76, T45, T79 | - | Bc Bs Ch H W 2R |
| T84 | {2,5,9,11,12,13} | VI | 9 | final_referee_inventory83_native_uniform_20261023, p.1 | 1 / 8 (44 pp) | T23, T25, T28, T49, T59 | - | Bs Ch H W 2R R |
| T85 | {0,5,7,9,14,18} | VI | 10 | blind_inventory84_rank_drop_uniform_20261023, p.1 | 1 / 37 (165 pp) | T40, T81, T55 | - | Bc H W MC |
| T86 | {4,5,8,9,12,13} | VI | 8 | root_inventory85_native_uniform_corrected_20261023, p.1 | 2 / 10 (54 pp) | T37, T57, T58, T76, T75 | - | Bs B1 Bk H W MC R |
| T87 | {4,5,7,9,12,14} | VI | 10 | blind_inventory86_matching_section_uniform_20261023, p.1 | 1 / 31 (138 pp) | T68, T51, T61 | - | Bc Ch H W MC |
| T88 | {1,5,6,10,13,14} | V | 1 | final_referee_inventory87_uniform_20261023, p.1 | 1 / 13 (69 pp) | - | - | Bs Ch G W MC |
| T89 | {2,5,9,10,13,14} | VI | 12 | blind_inventory88_matching_curve_uniform_20261023, p.1 | 1 / 15 (81 pp) | S7.24, T79, T32 | - | Bc Ch H MC |
| T90 | {1,5,6,9,14,18} | VI | 2 | blind_inventory89_matching_curve_uniform_20261023, p.1 | 1 / 11 (45 pp) | S7.24, T40 | - | Bc Bs H W MC |
