Evidence level: mixed. Dependency structure, set identities, orbit memberships, printed group actions,
tie lists and witness coverage are EXACT (own integer code, `build_dependency_map.py`, results in
`dependency-map.json` under `conditional_theorems`, `printed_witness_table_checks`). Whether each
hypothesis is used, and whether the discharging text proves the needed statement, is PROOF-READ at
the level of theorem statements, hypothesis lines, boundary tables and dependency paragraphs only.
The mathematics inside the discharging proofs (Gram tables, radial polynomials, curve arguments) is
CITED / UNCHECKED here; it is assigned to the case auditors.

# Conditional theorems in the contact-case proofs: hypothesis and discharge

## Deviations and limits

- "Discharged" below means: a document in the supplement states and claims to prove exactly the
  hypothesis (same contact set up to the group G, verified by my code), the using case is not in that
  document's dependency closure (graph plus text citations), and the assembly note applies it with
  no remaining hypothesis. It does not mean the discharging proof is correct.
- Page numbers are pages of the compiled PDFs as counted from the text extracts (form feeds).
- Search for conditionals: all 208 text extracts were searched (patterns: theorem statements
  beginning "Assume"/"Suppose"/"Under the hypothesis", "conditional", "provided that", "explicit
  hypothesis", "open input", "is a hypothesis", "not proved in this note"). 55 documents match. Apart
  from the ten conditional theorems below, the hits are (a) coefficient-region lemmas ("Suppose
  b > 1, 0 < a < 1/b, ...", T46/T47/T48 allocation notes, T44 branches, T9/T10 branch notes),
  which are partial theorems on a region of one cone and are combined inside the same case, not
  hypotheses on another contact set; (b) "Assume eight distinct projective points" inside proofs;
  (c) the general reduction `final_referee_moc4_six_contact_reduction` (Theorem 3: if every point
  with at least six contacts is safe, all are), whose input is the full set of six-contact cases and
  which is applied only after all ninety are proved (this is the paper's Theorem 3.6 route, no
  circularity with any single case); (d) the two-endpoint reduction for T_9 in
  `independent_moc4_T8_polar_reduction` / `independent_moc4_T8_uniform` (Theorem 1 there reduces
  T_9 to two seven-contact faces; `independent_moc4_T8_uniform` proves them as Theorem 3 (p.6,
  T^(7)_3) and Theorem 2 (p.4, T^(7)_4), with literal sets {0,1,2,3,11,12,18} and
  {0,1,2,3,11,12,20}, whose orbits my code matches to T^(7)_3 and T^(7)_4).

## Files with "conditional" in the name (40 files)

| File (all in proofs/manin_transfer_20261021/) | Role |
|---|---|
| blind_inventory49_bounded_native_conditional_20261023.pdf | PDF proof, conditional theorem C3 (T_50) |
| blind_inventory49_bounded_native_conditional_frozen_manifest_20261023.json | manifest |
| blind_moc4_inventory52_native_harmonic_conditional_20261023.pdf | PDF proof, C4 (T_53) |
| blind_moc4_inventory52_native_harmonic_conditional_frozen_manifest_20261023.json | manifest |
| final_referee_inventory57_conditional_hostile_audit_20261023.txt | audit note (C6) |
| final_referee_inventory57_conditional_hostile_exact_20261023.{py,json} | exact-check script (not executed) and its JSON result |
| independent_inventory67_native_harmonic_conditional_20261023.pdf | PDF proof, C7 (T_68) |
| independent_inventory67_native_harmonic_conditional_exact_20261023.{py,json} | exact-check script / JSON result |
| independent_inventory67_native_harmonic_conditional_frozen_manifest_20261023.json | manifest |
| independent_inventory75_native_harmonic_conditional_20261023.pdf | PDF proof, C8 (T_76) |
| independent_inventory75_native_harmonic_conditional_exact_20261023.{py,json} | exact-check script / JSON result |
| independent_inventory75_native_harmonic_conditional_frozen_manifest_20261023.json | manifest |
| independent_true_T48_native_harmonic_conditional_20261023.pdf | PDF proof, C9 (T_49) |
| independent_true_T48_native_harmonic_conditional_exact_20261023.{py,json} | exact-check script / JSON result |
| root_inventory50_bounded_native_conditional_20261023.pdf | PDF proof, C1 (T_51) |
| root_inventory51_bounded_native_conditional_20261023.pdf | PDF proof, C2 (T_52) |
| root_inventory51_bounded_native_conditional_exact_20261023.{py,json} | exact-check script / JSON result |
| root_inventory51_bounded_native_conditional_frozen_manifest_20261023.json | manifest |
| root_inventory52_conditional_native_hostile_audit_20261023.txt | audit note (C4) |
| root_inventory52_conditional_native_hostile_exact_20261023.{py,json} | exact-check script / JSON result |
| root_inventory56_conditional_native_hostile_audit_20261023.txt | audit note (C5) |
| root_inventory56_conditional_native_hostile_exact_20261023.{py,json} | exact-check script / JSON result |
| root_inventory57_native_conditional_20261023.pdf | PDF proof, C6 (T_58) |
| root_inventory57_native_conditional_exact_20261023.{py,json} | exact-check script / JSON result |
| root_inventory57_native_conditional_frozen_manifest_20261023.json | manifest |
| root_inventory67_conditional_native_hostile_audit_20261023.txt | audit note (C7) |
| root_inventory67_conditional_native_hostile_exact_20261023.{py,json} | exact-check script / JSON result |
| root_inventory75_conditional_native_approved_frozen_manifest_20261023.json | manifest (approval record) |
| root_inventory75_conditional_native_hostile_audit_20261023.txt | audit note (C8) |
| root_inventory75_conditional_native_hostile_exact_20261023.{py,json} | exact-check script / JSON result |

(40 files: 8 PDFs, 10 exact-check scripts, 10 JSON results, 7 manifests, 5 audit notes; the
exact list with roles is `conditional_files` in `dependency-map.json`. None of the scripts was
executed, per the audit rules; only the PDF text was read.)
Two conditional theorems have no "conditional" in their file name:
`final_referee_inventory56_native_harmonic_20261023.pdf` (C5, T_57) and
`independent_moc4_T23_harmonic_reduction_20261022.pdf` (C10, T_24).

## Common pattern

Each conditional theorem concerns one paper case T. Its non-minor generator analysis leaves seven
"relevant" extra ties t (ties at which a non-minor generator vanishing on T also vanishes; Lemma
4.2). Six are closed by earlier cases (printed witnesses or an unnamed "list of N preceding
cases"); the seventh face T+{t} is the explicit hypothesis, proved in a separate seven-contact
document. My code recomputed the relevant ties of every one of the ten cases from the generator
definitions; in all ten they agree exactly with the tie list printed in the document (in the
document's own frame for T_58 and T_60). This is an EXACT confirmation that the conditional
statements isolate precisely one boundary face each.

## C1. T_51 = {1,3,5,10,11,12} (root_inventory50_bounded_native_conditional, p.1, Theorem 1)

- Hypothesis: support inequality on the seven-contact orbit of D = {0,1,2,7,9,16,23} = T^(7)_12
  (coreD), literal match with Table 5. Used at tie 20: T_51+{20} = {1,3,5,10,11,12,20} is in the
  G-orbit of D (EXACT).
- Other boundary inputs stated: ties 6,13 -> B = coreB (T^(7)_8), 14 -> C = coreC (T^(7)_11), 0,7,16
  -> "strictly prior closed 60" (unnamed). EXACT: coreB covers T+{6} and T+{13}, coreC covers
  T+{14}; for ties 0,7,16 earlier candidates exist (e.g. T_9, T_11, T_42, T_65), all before stage 7.
- Discharge: `blind_moc4_coreD_uniform_20261023.pdf`, p.1, Theorem 1 ("seven tied permutations
  include {0,1,2,7,9,16,23}"), with companion `blind_moc4_coreD_low_pair` (lowest-two-coefficient
  region). Applied in `root_inventory50_bounded_native_whole_20261023.pdf` p.1, which names the
  coreD document explicitly.
- Circularity: coreD is stage 1 with no graph prerequisites; its text cites only pointwise
  sharp-SOC/old-H atoms and the global Hessian/principal-minor facts; coreD_low_pair says it "must
  not be quoted as entire coreD or entire inventory50 safety". T_51 is not in its closure.
- Classification: DISCHARGED.

## C2. T_52 = {1,3,5,10,11,14} (root_inventory51_bounded_native_conditional, p.1, Theorem 1)

- Hypothesis: G = {0,1,2,7,9,10,23} = T^(7)_10 (coreG), literal. Used at tie 20; T_52+{20} in the
  orbit of G (EXACT). The note's index remark ("index 51 with representative {3,4,5,7,9,12}") is
  the inventory representative of the same orbit.
- Other inputs: ties 6,13 -> coreA, 12 -> coreC (EXACT coverage), 0,7,16 -> "full families"
  (earlier candidates exist: T_10, T_11, T_12, T_42, T_46).
- Discharge: `independent_shared_coreG_uniform_20261023.pdf`, p.1, Theorem 1, literal G. Applied by
  `root_inventory51_bounded_native_whole_20261023.pdf` p.1, naming the coreG document.
- Circularity: coreG (stage 1, no graph prerequisites) states on p.5 that it uses "neither whole49
  nor whole51"; its inputs are the sharp-SOC (f1,f9;f22) and old-H pointwise atoms and the
  matching-cubic/Hessian facts. No path to T_52.
- Classification: DISCHARGED.

## C3. T_50 = {1,3,5,6,10,11} (blind_inventory49_bounded_native_conditional, p.1, Theorem 1)

- Hypotheses: F = {0,1,2,7,9,10,14} = T^(7)_7 (coreF) and G = T^(7)_10 (coreG), literal. Used at
  ties 13 and 20; T_50+{13} in orbit(F), T_50+{20} in orbit(G) (EXACT).
- Witness table (whole note p.1): 0 -> T9 {0,1,3,5,7,8} = paper T_10; 7 -> T42 = T_43; 16 -> T65 =
  T_66; 12 -> B; 14 -> A; 13 -> F; 20 -> G. All seven printed (row, column, transpose) actions were
  replayed with the note's stated convention (rho(i) = column[pi(row[i])], transpose first): every
  image lies in T_50 + tie (EXACT). Note the labels are historical: reading "T9" as paper T_9 would
  be wrong (paper T_9 does not cover the face). All witnesses are in earlier stages (T_66: 7 < 8).
- Discharge: `blind_shared_coreF_uniform_20261023.pdf` p.1 (literal F) and coreG as in C2. The
  assembly `blind_inventory49_whole_20261023.pdf` (p.3) says "The established full F and G theorems
  discharge exactly the two explicit hypotheses" but does not name the two documents; they are
  identified only through the displayed sets. coreF p.4: "Neither whole49, whole51, nor core G was
  used."
- Classification: DISCHARGED (minor citation gap: discharging documents not named in the assembly).

## C4. T_53 = {0,1,4,6,9,12} (blind_moc4_inventory52_native_harmonic_conditional, p.1, Theorem 1)

- Hypothesis: S52 = {0,1,4,6,9,12,17} = T^(7)_19 = T_53+{17} literally.
- Witnesses (conditional p.1; assembly table p.1-2): ties 5,7,11,14,15,18 -> labels T4, T17, T25,
  T4, T25, T47 = paper T_4, T_18, T_26, T_4, T_26, T_65 (historical labels; the assembly warns that
  "T47" is not inventory 47). EXACT: each covers T_53+tie, all earlier than stage 8. The same labels
  are reused, with the same ties, for the zero-coefficient faces of the seven-contact face S52
  (8-element sets), also EXACT-covered.
- Discharge: the seven-contact theorem is proved inside the T_53 primary
  `blind_moc4_inventory52_whole_assembly_20261023.pdf`, pp.1-3 ("The strictly positive
  seven-contact cone", region table (2), four component notes
  `blind_moc4_inventory52_plus17_{low,middle,third,fourth}_products`). The assembly p.3 applies the
  conditional theorem.
- Circularity: the S52 argument uses only the six earlier cases above, the global real-Hessian
  theorem and the plus17 notes; none cites T_53's six-contact statement. Graph: T^(7)_19 stage 7 <
  T_53 stage 8.
- Remarks for the case auditor: the assembly itself flags an "incidental phrase 'including b = c'"
  in the third-products note as vacuous; region 1 relies on the "Hessian real-branch criterion"
  (Appendix C input).
- Graph remark: the graph gives T_53 no six-contact prerequisite covering ties 5,7,11,14,15,18;
  they are reached only transitively through T^(7)_19.
- Classification: DISCHARGED.

## C5. T_57 = {0,1,5,7,9,12} (final_referee_inventory56_native_harmonic, conditional theorem)

- Hypothesis (as summarised by the joint assembly p.1, "sole additional hypothesis the support
  inequality on S"): S = {0,1,5,7,9,12,17} = T^(7)_20 (S56) = T_57+{17} literally.
- Discharge: `final_referee_inventory56_boundary17_uniform_20261023.pdf` p.1 Theorem 1, literal S;
  its zero-coefficient witnesses (p.1) T43, T25, T42, T33, T30, T16 = paper T_44, T_26, T_43, T_34,
  T_31, T_17 for gaps 18, 15, 11, 14, 6, 4 are EXACT-covered and all earlier (stages 1-4 < 5); the
  note states the 52-family list contains neither inventory 56 nor 57. Applied in
  `final_referee_inventory56_57_whole_assembly_20261023.pdf` p.1.
- Classification: DISCHARGED.

## C6. T_58 = {0,3,5,7,9,12} (root_inventory57_native_conditional, p.1, Theorem 1)

- Hypothesis: the same S = {0,1,5,7,9,12,17}. The note works in the frame T* = {0,1,5,7,9,17}
  (image of T_58 under sigma -> sigma o (1 3)), where S = T*+{12}. EXACT: T_58+{17} is in the orbit
  of S, and T* is in the orbit of T_58.
- Other ties (4,6,11,14,15,18 in the T* frame; 2,6,11,13,16,23 literally, recomputed EXACT) ->
  "immutable 49-family list" (unnamed).
- Discharge: as C5, same document; applied in the same joint assembly p.1-2.
- Classification: DISCHARGED.

## C7. T_68 = {1,4,5,7,9,12} (independent_inventory67_native_harmonic_conditional, p.1)

- Hypothesis: {1,4,5,7,9,12,14} = T^(7)_27 (S67) = T_68+{14} literally.
- Other ties 0,6,11,15,17,18 -> "fixed preceding 42-family list" (unnamed); EXACT: each has earlier
  candidate covers (e.g. T_5, T_29, T_42, T_8).
- Discharge: `root_inventory67_g14_uniform_20261023.pdf` p.1 Theorem 1 (literal S67), summarised
  again in `root_inventory67_whole_assembly_20261023.pdf` p.1 ("The seven-contact boundary
  inequality"), which then applies the conditional theorem (p.2).
- Circularity: the g14 note's walls use "literal old42 closed six-tuples"; S67 is stage 5 < T_68
  stage 6. No path to T_68.
- Classification: DISCHARGED.

## C8. T_76 = {0,1,3,10,11,12} (independent_inventory75_native_harmonic_conditional, p.1)

- Hypothesis: {0,1,3,10,11,12,16} = T^(7)_18 (S75) = T_76+{16} literally.
- Second seven-face used (not a hypothesis): T_76+{14} is claimed to be an image of the proved
  {1,4,5,7,9,12,14} (S67). EXACT: true. S67 is stage 5 < 6.
- Ties 5,6,7,13,20 -> "immutable 42-family list" (unnamed); earlier candidates exist for each.
- Discharge: `final_referee_inventory75_boundary16_uniform_20261023.pdf` p.1 Theorem 1 (literal
  S75); its zero-coefficient witnesses T39, T28, T28, T43, T30, T10 = paper T_40, T_29, T_29, T_44,
  T_31, T_11 (gaps 20,14,13,6,7,5) are EXACT-covered and earlier; "No circular inventory75 input is
  used." Applied in `final_referee_inventory75_whole_assembly_20261023.pdf` p.1.
- Classification: DISCHARGED.

## C9. T_49 = {0,3,5,6,9,12} (independent_true_T48_native_harmonic_conditional, p.1)

- Hypothesis: {0,3,5,6,9,12,17} = T^(7)_23 (S48c0) = T_49+{17} literally (the "c = 0 face").
- Witnesses (p.2): ties 2,7,11,13,16,23 -> T11, T36, T20, T41, T33, T43 = paper T_12, T_37, T_21,
  T_42, T_34, T_44: EXACT-covered, all earlier (stages 1-3 < 4).
- Discharge: proved in the T_49 primary `root_true_T48_whole_assembly_20261023.pdf` p.3 ("The
  seven-contact frontier satisfies the stated support inequality"), using the component notes
  `independent_true_T48_c0_native_Heron`, `blind_moc4_T48_c0_native_interlacing` (orientation
  theorem only) and `blind_moc4_T48_c0_all_j_chain`; applied on p.4.
- Circularity: none found. S48c0 is stage 2, T_49 stage 4. The zero-coefficient faces of the
  frontier are said to use "an recorded earlier 42-family list" (unnamed).
- Side finding: the older allocation notes `blind_moc4_T48_allocation_{I,II,III,IV}` (listed as
  branches of T_49 in the catalogue, but explicitly not used by the assembly) send one zero face to
  "T23" = paper T_24, which is stage 5, after T_49 (stage 4). Harmless only because the assembly
  does not use those notes; an auditor should not rely on them.
- Classification: DISCHARGED.

## C10. T_24 = {2,3,4,6,7,8} (independent_moc4_T23_harmonic_reduction, p.1, Theorem 1)

- Hypotheses: T+{16} and T+{17} (two seven-contact families). T+{13} is handled separately.
- Discharge (stated in the primary `independent_moc4_T23_uniform_20261022.pdf`, p.1 table and
  p.4): T+{17} by `independent_moc4_T23_g17_uniform` (literal {2,3,4,6,7,8,17}; its orbit is
  T^(7)_5); T+{16} by the six-contact case "T24" = paper T_25 through a printed action (row 0312,
  column 1203, transpose); T+{13} by `independent_moc4_T23_g13_uniform` (T^(7)_6).
- EXACT checks (`conditional_theorems[C10]` in the JSON): the relevant ties of T_24 computed from
  the generators are 1,10,12,13,16,17,22, as printed; every printed witness (T14 = T_15,
  T19 = T_20, T12 = T_13, T11 = T_12, T24 = T_25, T+13 in orbit of T^(7)_6, T+17 in orbit of
  T^(7)_5) covers its face and is earlier than T_24. The literal hypothesis {2,3,4,6,7,8,17} is in
  the orbit of T^(7)_5 but is not the Table 5 representative (expected: the documents work in the
  T_24 frame). The printed action for tie 16 (row 0312, column 1203, transpose) maps T_25 into
  T_24+{16} under the docs/notation.md convention sigma -> rho o sigma o chi, but not under the
  convention printed in the inventory-49 note; the two families of notes use different action
  conventions. Coverage itself does not depend on the convention.
- Stages: T_25, T_15, T_20 stage 1; T_12 stage 3; T_13 stage 4; T^(7)_5 stage 1; T^(7)_6 stage 4;
  T_24 stage 5. Closure of T^(7)_5 does not contain T_24. Acyclic.
- Note: the primary says the conditional harmonic theorem "was separately proved in
  root_moc4_T23_harmonic_independent_exact_20261022.py", i.e. it points to a program rather than
  to the PDF reduction note for the proof of the reduction; the PDF note exists and states it.
- Classification: DISCHARGED.

## Summary

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
