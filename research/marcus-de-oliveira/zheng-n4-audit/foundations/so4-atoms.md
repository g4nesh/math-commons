Evidence level: EXACT (97 identity and finite checks by own code in `so4_atoms.py`, 0 failures) plus PROOF-READ arguments for the four core atom families SOC, old-H, chain and four-coefficient. NUMERIC only for the rank-two family and the face-specific Gram families. No atom was found to fail.

# Universal SO(4) inequalities ("atoms") used by the contact-case proofs

> **Update, 26 November 2026.** The rank-two family and the face-specific Gram families listed below as not proved were settled by two later passes: [so4-atoms-part2a.md](so4-atoms-part2a.md) (rank-two; the faces of historical T24, T25, T41 and inventory 56) and [so4-atoms-part2b.md](so4-atoms-part2b.md) (the faces of historical T44 to T47 and the allocation lemmas). The four core families were re-derived independently in [so4-atoms-verification.md](so4-atoms-verification.md).

Label: `atoms`. Files: `so4_atoms.py` (checks and attacks), `so4_atoms_output.txt` (captured run, 114 s), `so4_atoms.json` (machine-readable inventory).

## Deviations and limits

- **Proved here: four families.** SOC, old-H, chain and four-coefficient are proved completely, one representative per orbit, and transferred to the whole orbit by an exact covariance check. The rank-two family and the face-specific "coherent / parallel / double-SOC / partial-four / six-ray" Gram families are NOT proved here; they are only tested numerically (section 4). Each rests on a degree-six quaternion Gram identity whose cubic basis appears only in supplement checker data, not in any PDF.
- **Chain certificate.** The chain proof uses certificate DATA from the supplement: the cubic features, the sign maps and the congruence matrix T of the bad70 face. They are read with `json.load` from `MOC4_SUPP`. Nothing is executed. Every identity is verified exactly by my own code, so the proof does not depend on where the polynomials came from. No PDF of the supplement contains this proof: bad70 p.1 only says the identities are in its verifier.
- **Old-H argument.** The increment argument for old-H is the one written by the sibling gapcheck (`docs2/blind_moc4_bad70_full_safe_disjunction_20261022__gapcheck.md`). I re-read it line by line and re-verified its identities with independent code (A4).
- **Document inventory.** Two sources:
  - The "documents naming the family" lists come from a keyword scan of all 208 text extracts (page = form-feed count). They over-include. A document listed there names the family, but may not use it.
  - The role tuples in the cited-tuple table were read in context from the text. They cover about 60 uses, not every use. Many documents apply the atoms through "all 1152 literal actions" recorded only in checker data, and those were not enumerated.
- **Scope check.** A use in a form stronger than proved was checked only for those cited tuples (orbit membership with the printed roles, all pass). The coefficient forms printed with them were read and match the pointwise forms below. The face-specific Gram claims were not checked against proved statements because none is proved here.
- I did not read the supplement scripts beyond the docstrings and the data layout of the bad70 Gram check. I did not read the paper's Appendix C.7 beyond the passage on quaternion coordinates. That passage concerns (C.9)/(C.10), not the atoms. Table 1 group III says only that "Gram inequalities are established in their respective proofs".

## 1. Model, minors, covariance (EXACT, A0-A1)

Every O in SO(4) is x -> p x conj(q), and equally x -> p x q, for unit quaternions p, q: the double cover S^3 x S^3 -> SO(4), with q -> conj(q) a bijection of S^3. The matrix O(p,q) has entries bilinear in (p, q). Checked exactly in both conventions: O^T O = |p|^2|q|^2 I and det O = (|p|^2|q|^2)^2.

For all 34 minor generators the Leibniz sum is Phi_f = x_f^2 for a 2x2 minor x_f, and |p|^2|q|^2 x_f^2 for an entry; this is (4.2) in the model. Write x_f for the minor (or entry) with Phi_f = x_f^2 on S^3 x S^3.

Covariance:
- Checked with a generic symbolic O (16 independent entries) for all 1152 g: w_{g.sigma}(g.O) = eps(g) w_sigma(O), where eps(g) = det(g.O)/det(O) = +-1.
- Hence Phi^{O(4)}_{g.c}(g.O) = Phi^{O(4)}_c(O), with Phi^{O(4)} = Phi/det.
- A column sign flip D gives Phi^{O(4)}_c(OD) = Phi_c(O) and leaves every |x_f| unchanged. So a pointwise atom that holds on SO(4) also holds on O(4) in the Phi^{O(4)} form.
- G permutes the 82 f_j and preserves minor/non-minor (1152 distinct permutations).

Every atom below is a statement in the Phi_f and |x_f| alone. So an atom proved for one representative holds for every member g.(representative) of its orbit (PROOF-READ, using the exact identity above).

## 2. Inventory

Throughout, d > 0 is the non-minor coefficient; d = 0 is trivial because all other generators are minors. Minimising over the free coefficients turns each coefficient statement into an equivalent pointwise form; this is PROOF-READ, by AM-GM or Cauchy-Schwarz with the infimum over each ratio.

| family | coefficient statement (as used) | equivalent pointwise form on SO(4) | representative (paper f_j) | orbit under G (1152) | where stated | where the Leibniz form is proved |
|---|---|---|---|---|---|---|
| SOC ("sharp SOC", "old SOC", "144-orbit SOC", "pair-SOC") | a f_i + b f_j + d f_k safe if ab >= d^2; equivalently tau Phi_i + Phi_k + Phi_j/tau >= 0 for all tau > 0 | Phi_k + 2\|x_i x_j\| >= 0 | (f1 = n12,34, f9 = n14,23; f22 = h3412,+) | 144 (pair unordered) | T2 scalar note p.1 (scalar only); T3 p.1 eq. (1) | **here** (A3 + proof 3.1). Supplement: scalar only, the bridge is explicitly not supplied (T2 note p.2); the bridge exists only in its scripts |
| old-H ("old-H", "one-pair master", "H capacity") | t f_t + c f_c + d f_k safe if t(c - d) >= d^2 | Phi_c + Phi_k + 2\|x_t x_c\| >= 0 | (aux f3 = n14,34, main f15 = n12,14; f26 = h4132,-) | 288 (ordered) | coreB/coreF/coreH/coreL (reference), bad70 (H) endpoints | **here** (A4 + proof 3.2). Supplement: scalar form in alternate_onecircle_elementary (Theorem 1, as an SO(3) x SO(3) statement); the Leibniz bridge exists only in scripts |
| chain ("ordinary / three-minor chain", bad70 (PB), (CQ)) | p f_l + b f_m + a f_c + d f_k safe if a >= d^2/p + d^2/b | Phi_k + 2\|x_c\|(\|x_l\| + \|x_m\|) >= 0 | (centre f2 = n13,34; leaves f0 = n1,4, f12 = n14,24; f70 = h4231,+) | 576 | bad70 p.1 (PB), (CQ) as Gram forms | **here** (A5 + proof 3.3), via the supplement's certificate data. No PDF proof |
| four-coefficient ("four atom", "four-budget") | a f_a + b f_b + c f_c + d f_d safe if ab(c - d) >= d^2(a + b + c + d) | Phi_d + Phi_c + 2(\|x_a x_b\| + \|x_a x_c\| + \|x_b x_c\|) >= 0 | (f0 = n1,4, f6 = n4,1; c: f10 = n13,24; f54 = h4231,-) | 288 (a, b unordered) | final_referee_T5_four_budget_capacity p.1-2 | **supplement PDF, complete**; replayed EXACT here (A6) and in the first pass |
| rank-two | b f12 + j f2 + a f29 + d f35 safe if b >= d/4 + d^2/(4 min(a, j)) | Phi_35 + Phi_12/4 + \|x_12\| sqrt(Phi_2 + Phi_29) >= 0 | (leaves f2, f29; centre f12; f35 = h4312,-) | 144 | inventory56_rank_two p.1 | not proved here (Gram data only) |
| face Gram families ("coherent", "parallel", "double-SOC", "partial-four", "core coherent", "six-ray", "coherent-J stationary", "merged-C", "PCBE", "rational capacity") | explicit PSD constructions in "the ten completion relations" of a face's degree-six Gram | no closed pointwise form; each is a fixed-face coefficient region | T24 face (f0,f2,f59,f6,f53,f54,f55,f67); T25 face; T44 face; T46 face; T45, T47 faces | face-specific | T24 p.2, T25 p.1-2, T44 coherent p.1-2, T46 core p.1-2, T46 six-ray p.1, T41 bad17 p.2, T45 p.2, T47 notes | not proved here (Gram data only) |
| derived: old-H sum / Hsum / bad70 (H); H/chain and four/chain allocations; coreE four-chain | conic sums of the above with explicit allocations | - | - | - | bad70, T28, T38, T42, T36, T48 allocation notes, coreE four-chain | reduce to proved atoms; the allocation algebra is not audited here (coreE and bad70 (H) were in the first and second pass) |

### Cited tuples, checked against the orbits with the printed roles (EXACT, A2; all pass)

- **SOC.**
  - (1,11;17): T3 p.1.
  - (1,9;22): reference for T2, coreB/F/G/K/M/O, inventory87, coreE.
  - (2,12;26), (2,15;26), (12,15;26): T10 A-branch p.1, T10 p.2, coreB/C/D/F/G.
  - (5,2;54): T42 p.1. (2,5;70): T38 p.4. (1,9;18): T8 p.4 and T8 double boundary p.1. (2,39;46): T7 boundary23 p.1. (2,8;69): T35 p.2. (2,15;56): T37 p.2. (1,12;35): inventory67 g14 p.1.
- **old-H** (aux, main; k):
  - (3,15;26), (23,15;26): coreB, coreH, coreL, coreM, inventory71.
  - (12,2;70), (39,2;70): T38 p.4, bad70.
  - (3,5;54), (23,5;54): T42 p.1, T41 bad54 p.1.
  - (9,2;69), (10,8;69), (53,8;69): T35 p.2.
  - (29,12;35), (2,12;35): T39 p.2, inventory56.
  - (9,5;68): T34 p.2. (11,12;26): coreF p.2. (15,10;30), (8,10;30): inventory75 boundary16 p.1. (33,2;46), (5,2;46): T28 p.4. (9,2;56), (10,15;56), (53,15;56): T37 p.2.
  - In every case the reversed role order is NOT in the orbit. So the documents' "aux, main" order is meaningful, and they print it correctly.
- **chain** (centre; leaves; k):
  - (2;0,12;70), (2;39,67;70): bad70, coreB, coreC, T38. (15;0,3;26): coreB, coreE.
  - (5;0,3;54), (5;23,6;54): T42. (5;7,9;68): T34. (2;0,15;54): T41 bad54. (15;0,10;56): T37. (2;4,33;46), (2;5,6;46): T28. (59;53,55;54), (59;53,67;54): T24.
- **four** (a,b; c; d):
  - (0,6;10;54): T5. (6,0;5;34): T4 p.1. (0,13;3;26): coreE.
  - (0,7;12;40): T38. (7,0;1;68), (7,67;53;68): T34, T35. (0,67;2;40), (67,41;39;40): T27.
  - (7,21;9;35): inventory52 middle products. (4,6;2;28): T9 second branch. (4,7;3;22): T10, T10 B-branch.
  - (0,6;1;70), (6,67;39;70): T25 eq. (5). (67,55;53;54): T24, K >= 0.

### Documents naming each family (keyword scan; page numbers; over-inclusive)

- **SOC (64):**
  - blind_moc4: T27 (1,2), T38 (2,4,5,6), T42 (1,3), T44_coherent_chain_capacity (2), T44_conic_branch_exclusion (1), coreB (1,2), coreC (1,2), coreD (2,4).
  - blind_shared: coreF (2).
  - final_referee:
    - T10_A_branch (1,4), T25 (1,2,3,6), T26_boundary16 (1,2), T26 (1), T2_scalar_compensation (1,2), T35 (2).
    - T41_bad17_large_product (1), T41_bad17 (2), T41_bad54 (1), T41 (2), T5_four_budget (3), T9_first_branch (1,4).
    - inventory56_rank_two (3), inventory71 (1), inventory87 (2,4).
    - moc4_T10 (2), moc4_T11 (1), moc4_T12 (1), moc4_T14 (1,3), moc4_T5 (1), moc4_T9 (1).
    - shared_coreK (1,3).
  - independent:
    - coreE (1,2), inventory67_native_harmonic_conditional (4).
    - moc4: T23_g17 (2,4), T23_harmonic_reduction (1), T23 (2), T2 (1), T36 (1,4,5,6), T39_high_j (1), T39_high_p (1), T39 (1,2,3), T3 (1,4), T43 (1,2,4), T44_a_Heron_branch (3), T44_c_below_j (1,2), T45 (2), T47_IIIb_plus_unswapped (1,3), T47_I_II (1), T8_double_boundary (1,3), T8 (4,5,6,8).
    - shared: coreG (2,3,4,5), coreM (2,3,5).
  - root:
    - inventory67_g14 (1).
    - moc4: T24 (2,4), T37 (1,2), T44_e_below_p (1), T47_IIIa_plus (1,2), T7_boundary16 (1), T7_boundary23 (1), T7 (3).
    - shared: coreH (1,3), coreL (1,4), coreO (1,3), seven_A (1,4).
- **old-H (50):**
  - blind_moc4: T27 (1), T28 (2,3,5), T38 (2,4), T42 (1), T44_conic_branch_exclusion (1), T48_allocation_III_low_j (1), T48_allocation_II (1), T48_allocation_IV (1), bad70 (1,2), coreB (1,3), coreC (1,2,3), coreD (2,4).
  - blind_shared: coreF (2).
  - final_referee: T25 (1), T34 (1-4), T35 (1-3), T41_bad54 (1), T41 (2), inventory56_rank_two (3), inventory71 (2,4), inventory75_boundary16 (1), inventory75_whole_assembly (1).
  - fresh_strong_principal_minors (8).
  - independent:
    - moc4: T23_g17 (2,4), T36 (1,4), T39_high_j (1), T39_high_p (1), T39 (1-3), T43 (1,2), T44_a_Heron_branch (2,3), T44_c_below_j (1,2), T45 (1-4,6).
    - moc4_T47: IIIb_plus_large_pc (1), IIIb_plus_unswapped (1-3), I_II (2), high_j (1,3), merged_C_high_j (2).
    - shared: coreG (2,5), coreM (2,5).
  - root:
    - moc4: T24 (2,4), T30 (1), T37 (1,3), T44_e_below_p (1), T47_IIIa_plus (1), T47_whole_assembly (2).
    - shared: coreH (1), coreL (2), seven_A (1,4).
  - principal_minors (8/9, two copies).
- **chain (43):**
  - blind_moc4:
    - T28 (2,3,5), T38 (4), T42 (1,2), T44_coherent_chain_capacity (1,2).
    - T48: allocation_III_low_j (1,2), allocation_II (1), allocation_IV (1), allocation_I_low_j (1,2), allocation_I_small_lambda_low_j (1,2), c0_all_j_chain (1,2), c0_native_interlacing (3), c0_upper_Heron_chain (1).
    - coreB (1,3), coreC (1,2,4).
  - final_referee:
    - T34 (1-4), T35 (2,3), T41_bad54 (1), T41 (2).
    - inventory56_57_whole_assembly (1), inventory56_boundary17_large_product (1), inventory56_boundary17 (1), inventory56_rank_two (3), inventory71 (1-4), inventory75_boundary16 (1), inventory75_whole_assembly (1).
  - independent:
    - coreE_four_chain (1,2), coreE (1,2,4), inventory67_native_harmonic_conditional (4).
    - moc4: T23_g17 (4), T36 (1,4), T43 (4), T45 (2,4), T47_I_II (1).
    - shared_coreM (2-5), true_T48_c0_native_Heron (3).
  - root:
    - complete_negative_mass (2), inventory67_g14 (1,3).
    - moc4: T24 (4), T30 (1,2,4), T37 (1,2), T47_whole_assembly (2).
    - shared_seven_A (1,3,4), true_T48_whole_assembly (2,3).
- **four-coefficient (35):**
  - blind_moc4: T27 (1,2,4), T28 (2,5), T38 (2-4), T44_coherent_chain_capacity (2), inventory52_plus17_fourth (1,3) / middle (1,2) / third (1,3), inventory52_whole_assembly (2).
  - final_referee:
    - T25 (1,2), T26_boundary10 (1), T26 (1), T34 (1-4), T35 (2), T46_positive_a_large_product (2), T46_positive_a_small_or_balanced (1), T5_four_budget (1,3), T9_second_branch (1,2).
    - inventory56_57_whole_assembly (1), inventory56_boundary17_large_product (1), inventory56_boundary17 (2), inventory56_rank_two (3).
    - moc4_T10 (2), moc4_T11 (1), moc4_T12 (1), moc4_T9 (1).
  - independent: coreE_four_chain (1,2), coreE (1), moc4_T23_g17 (2,4), moc4_T23 (1), moc4_T43 (1-3), moc4_T44_c_below_j (1).
  - root: moc4_T10_B_branch (1,3), moc4_T24 (2,4), moc4_T30 (1,2), moc4_T4 (1,2).
- **face Gram families (46):**
  - blind_moc4: T27 (1,2,4), T38 (2,3), T44_coherent_chain_capacity (1,2), T48_negative_Phi (3), bad70 (1,2).
  - final_referee:
    - T25 (1,2,6), T34 (1), T41_bad17_large_product (1,3), T41_bad17 (2).
    - T46: a0_balanced_leaf (1), a0 (2), core_coherent_capacity (1,2,5), positive_a_large_product (1), positive_a_small_or_balanced (2), six_ray_coherent_capacity (1,4).
    - inventory29 (4), inventory56_57_whole_assembly (1), inventory56_boundary17 (2,3), inventory56_rank_two (1-3), inventory75_boundary16 (3), inventory78 (4), inventory81 (4), inventory87 (1,4).
    - single_collision_divisor_compensation (2).
  - independent:
    - coreE (5).
    - moc4: T23_g17 (2,4), T23 (1,2), T36 (1,2,4), T43 (1,2,4), T44_c_below_j (1,2), T45 (2).
    - moc4_T47: IIIb_plus_unswapped (1,3), I_II (1,2), high_j_small_b (1), high_j (1,3), merged_C_high_j (1).
  - root:
    - moc4: T24 (2,4), T44_e_below_p (1,2), T46_positive_a_mixed (1), T46_whole_assembly (2).
    - moc4_T47: IIIb_minus_swapped (1), low_j_high_e (1), low_j_p_below (1), zero_q_rational_capacity (1).
    - true_inventory47 (4).
  - one_circle (6).

## 3. Proofs

Notation for 3.1 and 3.2: on S^3 x S^3 put a_i = P_i^2 and A_i = Q_i^2, so that sum a = sum A = 1, and v = a - A, so that sum v = 0. An identity of bihomogeneous polynomials of equal bidegree that holds on S^3 x S^3 is a polynomial identity, by scaling. So "reduce modulo the sphere relations" means: homogenise both sides with |p|^2 and |q|^2 and compare exactly (A3, A4). For the plain-convention certificates (A5, A6) I reduce modulo N = |p|^2 - |q|^2 instead, by repeatedly replacing p0^2 with |q|^2 - p1^2 - p2^2 - p3^2. Since N is monic in p0, f is in (N) iff the reduced polynomial is 0. On unit p, q we have N = 0.

### 3.1 SOC: Phi_22 + 2|x_1 x_9| >= 0 (representative (f1, f9; f22))

**Bridge (EXACT, A3)**, model x -> p x conj(q), P = p, Q = q:

- x_1 = v_0 + v_1 and x_9 = v_0 + v_3;
- Phi_22 = 2(a_0A_2 + a_2A_0 + a_1A_3 + a_3A_1) - 4 p_0p_1p_2p_3 - 4 q_0q_1q_2q_3.

**Scalar theorem** (Theorem 2 of the T2 scalar-compensation note; re-derived here, with its identities EXACT in A3). Let a, A >= 0 be quadruples with equal sums, and put

- R = 2(a_0A_2 + a_1A_3 + a_2A_0 + a_3A_1) - 4 sqrt(a_0a_1a_2a_3) - 4 sqrt(A_0A_1A_2A_3),
- s = v_2 + v_3 and t = v_1 + v_2.

Then R + 2|st| >= 0.

Proof:
1. Put delta = sqrt(a_0a_2) - sqrt(a_1a_3) and Delta = sqrt(A_0A_2) - sqrt(A_1A_3). Then R = 2delta^2 + 2Delta^2 - 2(v_0v_2 + v_1v_3) identically, and st = v_1v_3 - v_0v_2 on sum v = 0. Hence R + 2|st| = 2delta^2 + 2Delta^2 - 4 min(v_0v_2, v_1v_3).
2. If the minimum is <= 0 we are done.
3. Otherwise (v_0, v_2) and (v_1, v_3) have opposite strict signs, since sum v = 0. Exchanging a and A swaps delta and Delta and fixes R and |st|, so we may assume v_0 = -alpha < 0, v_2 = -gamma < 0, v_1 = M > 0 and v_3 = N > 0.
4. The increment lemma sqrt((r+u)(s+w)) - sqrt(rs) >= sqrt(uw) (EXACT identity in A3) gives Delta - delta >= sqrt(alpha gamma) + sqrt(MN).
5. Hence 2delta^2 + 2Delta^2 >= (Delta - delta)^2 >= 4 sqrt(alpha gamma MN) >= 4 min(alpha gamma, MN).

No division occurs, so the zero-coordinate boundaries are included.

**Conclusion.** On S^3 x S^3:
- -4 p_0p_1p_2p_3 >= -4 sqrt(a_0a_1a_2a_3), and the same for q, so Phi_22 >= R.
- x_1 x_9 = (v_0+v_1)(v_0+v_3) = st on sum v = 0 (EXACT), since s = -(v_0+v_1) and t = -(v_0+v_3).

Hence Phi_22 + 2|x_1x_9| >= R + 2|st| >= 0. Minimising tau Phi_1 + Phi_22 + Phi_9/tau over tau > 0 gives exactly Phi_22 + 2|x_1x_9|, so the tau form holds for every tau > 0, unrestricted. Then a f1 + b f9 + d f22 >= 0 whenever ab >= d^2. By section 1 the statement holds on all 144 orbit members. **Status: proved, EXACT + PROOF-READ.**

### 3.2 old-H: Phi_2 + Phi_70 + 2|x_12 x_2| >= 0 (representative (aux f12, main f2; f70), in the orbit of (f3, f15; f26))

**Bridge (EXACT, A4).** Model x -> p x conj(q), Hadamard coordinates P = Wp and Q = Wq, where W = (1/2)[[1,1,1,1],[1,-1,-1,1],[1,-1,1,-1],[1,1,-1,-1]]. W is orthogonal with W^2 = I, so (P, Q) ranges over S^3 x S^3. Then:

- x_2 = v_0 + v_1 and x_12 = v_1 + v_2;
- Phi_2 + Phi_70 = 2[delta^2 + Delta^2 - (v_0v_2 + v_1v_3)], with delta = P_0Q_3 + P_1Q_2 and Delta = P_2Q_1 + P_3Q_0;
- x_2 x_12 = v_0v_2 - v_1v_3 on sum v = 0.

So the left side is 2[delta^2 + Delta^2 - 2 min(v_0v_2, v_1v_3)]. If the minimum is <= 0 we are done. Otherwise there are two cases.

**Case A: v_0, v_2 < 0 < v_1, v_3.** Write v_0 = -alpha, v_2 = -gamma, v_1 = M and v_3 = N.
1. By the triangle inequality, |delta| >= |X| with X = sqrt(a_1A_2) - sqrt(a_0A_3), and |Delta| >= |Y| with Y = sqrt(a_3A_0) - sqrt(a_2A_1).
2. Regrouping, X + Y = [sqrt((A_1+M)(a_2+gamma)) - sqrt(A_1a_2)] + [sqrt((A_3+N)(a_0+alpha)) - sqrt(A_3a_0)] >= sqrt(M gamma) + sqrt(N alpha).
3. So delta^2 + Delta^2 >= (X+Y)^2/2 >= 2 sqrt(alpha gamma MN) >= 2 min(...).

**Case B** (signs reversed) is symmetric, with X' = sqrt(a_0A_3) - sqrt(a_1A_2) and Y' = sqrt(a_2A_1) - sqrt(a_3A_0).

This is the argument of the docs2 gapcheck. I re-read it and re-verified its identities independently.

**Conclusion.** Minimising over the split of c gives: t f12 + c f2 + d f70 >= 0 whenever t(c - d) >= d^2. Write c = d + c' with tc' >= d^2. Then t Phi_12 + c' Phi_2 >= 2 sqrt(tc') |x_12x_2| >= 2d|x_12x_2|, and the remainder is d(Phi_2 + Phi_70 + 2|x_12 x_2|) >= 0. All 288 orbit members follow by section 1, including (f3, f15; f26). **Status: proved, EXACT + PROOF-READ.**

### 3.3 chain: Phi_70 + 2|x_2|(|x_0| + |x_12|) >= 0 (representative (f2; f0, f12; f70))

**Gram identity (EXACT, A5).** Model x -> p x q. Let F be the eight cubic features in each of the two parities given in the supplement's certificate data for the face {f0, f1, f2, f12, f39, f67, f70, f72}; A5 also checks that these are the generators of the contact set {2,3,5,9,10,12}. Put g = T^{-1} F, using the congruence T of the data. For the four signed coordinate maps S_j (quaternion multiplications) put u_k(y) = (g_k(S_j y))_{parity, j} in R^8. Then, modulo N = |p|^2 - |q|^2:

- 2 Phi_k = |p|^2 |u_k|^2 for each minor source k in {0, 1, 2, 12, 39, 67, 72};
- 2 Phi_70 = |p|^2 (|u_70|^2 - 2<u_0, u_2> - 2<u_2, u_12>) (the (PB) Gram of bad70 p.1);
- 2 Phi_70 = |p|^2 (|u_70|^2 + 2<u_2, u_39> - 2<u_2, u_67>) (the (CQ) Gram).

**Proof.** On unit p, q we have N = 0 and |x_k| = |u_k|/sqrt 2. By Cauchy-Schwarz,

2Phi_70 + 4|x_2|(|x_0| + |x_12|) = |u_70|^2 - 2<u_0,u_2> - 2<u_2,u_12> + 2|u_0||u_2| + 2|u_2||u_12| >= |u_70|^2 >= 0.

So the pointwise chain holds. Minimising p Phi_0 + b Phi_12 + (d^2/p + d^2/b) Phi_2 over p, b > 0 gives d·2|x_2|(|x_0|+|x_12|). So the coefficient form a >= d^2/p + d^2/b follows, and larger a only adds a minor. All 576 orbit members follow by section 1.

The sibling gapcheck's side observation was that the naive fixed-sign scalar reading of (PB) fails. That is consistent with this proof: the Gram vectors u_k are 8-dimensional, not scalars. **Status: proved, EXACT identity + PROOF-READ. The proof exists in the supplement only as checker data.**

### 3.4 four-coefficient: (f0, f6; f10; f54)

Replayed EXACT (A6), in the model x -> p x q:
- the bridge 4 Phi_f = n^2 Q_f mod N, for Q = (A-B)^2, (A+B)^2, D^2 and H;
- the certificate (4) of the T5 four-budget note, n[a(A-B)^2 + b(A+B)^2 + cD^2 + dH] = sum_j l(S_j y)^T Z l(S_j y) + m(S_j y)^T Z m(S_j y) mod N, coefficient by coefficient;
- det Z = 4d[ab(c-d) - d^2(a+b+c+d)].

Z is PSD iff ab(c - d) >= d^2(a + b + c + d), for a, b, c, d >= 0. This is the Schur complement argument of the PDF, re-derived: for d > 0 the condition forces c > d. Minimising over the capacity boundary by Cauchy-Schwarz gives the pointwise form Phi_d + Phi_c + 2(|x_ax_b| + |x_ax_c| + |x_bx_c|) >= 0 (PROOF-READ). All 288 orbit members follow by section 1. **Status: proved in the supplement and checked (EXACT + PROOF-READ).**

## 4. Attack (NUMERIC)

Method: multistart Nelder-Mead over (p, q) in R^8, normalised. Deterministic seeds: orbit sample 20261125, starts 1000+i. Each family's representative got 24 starts. Each other atom got 4 starts: 24 random orbit members per family, plus every cited tuple. The worst point was re-evaluated with mpmath at 50 digits.

| family | atoms tested (orbit) | min of pointwise form | mpmath at worst point | largest per-atom min |
|---|---|---|---|---|
| SOC | 35 (144) | 3.7e-32 | 4.2e-32 | 1.2e-12 |
| old-H | 43 (288) | 3.5e-31 | 3.5e-31 | 1.1e-12 |
| chain | 36 (576) | 1.8e-31 | 1.8e-31 | 1.5e-15 |
| four | 37 (288) | 1.0e-28 | 1.0e-28 | 9.1e-6 (not fully converged; positive) |
| rank-two | 25 (144) | 2.2e-26 | 2.2e-26 | 6.8e-14 |

- **SOC with explicit tau** (f1, f9; f22), 16 starts over (O, log tau): minimum 7.6e-30.
- **Sensitivity.** Multiplying the compensating term by 0.9 gives a negative minimum for SOC (-2.2e-5) and old-H (-8.9e-6), so the test can detect a violation. For chain it still gives 7e-27 >= 0. A separate scratch run found:
  - the chain stays nonnegative with factor 0.5 and turns negative at 0.3 (-9.4e-3);
  - scaling only the f0 leaf, the form is negative at 0 (-4.5e-2) and nonnegative at 0.5.
- **Face Gram families.** Boundary coefficient tests, d = 1, 6 random parameter sets on each stated boundary, 6 starts each. All minima are about 0 and none is negative:

| family | min Phi |
|---|---|
| T24 Lemma 1 (5) | 6e-31 |
| T24 Lemma 1 (6) (parallel) | 2e-30 |
| T24 double-SOC endpoint | 5e-31 |
| T25 Lemma 1 (9) | 1e-32 |
| T44 Prop. 1 | 3e-29 |
| T44 Prop. 2 | 7e-32 |
| T46 core coherent Thm 1 | 6e-31 |
| T46 six-ray Thm 1 | 5e-29 |
| inventory56 rank-two (1) | 1e-31 |
| coreE four-chain simple section | 6e-30 |

No atom has a negative minimum. Nothing in this section is a proof.

## 5. Sharpness and scope

**Sharpness:**
- SOC, old-H and four-coefficient are sharp: minima reach 0, and a 10% weaker compensation fails.
  - SOC equality: the T2 note's checkerboard family, reached at nonzero minors.
  - old-H equality: as in the gapcheck.
- Chain: **corrected after the independent re-verification** ([so4-atoms-verification.md](so4-atoms-verification.md)). An earlier version of this note said that the chain inequality survives halving its compensating term and that equality occurs only at degenerate points. Both statements were wrong. At half compensation the left-hand side is negative (-3.270e-4) at an exact rational point of SO(4); the threshold factor is about 0.78. Equality in the chain inequality is attained at a nondegenerate point with x_2 = 0.3175. The proof of the inequality itself is unaffected.
  - The coefficient form a >= d^2/p + d^2/b is therefore sufficient, with some slack pointwise.
  - The supplement claims necessity only for the fixed degree-six Gram cone (bad70 p.1, inventory56 p.1).
- Rank-two and the face Gram families reach 0 on their stated boundaries numerically.

**Stronger use than proved:** none found among the cited tuples:
- every tuple's roles (pair, aux/main order, centre/leaves, the c slot of the four atom) are in the proved orbit with the printed roles;
- the printed coefficient conditions are exactly the pointwise-equivalent forms above (tau unrestricted, d^2 on the correct side, absolute values present);
- the bad70 (H) "continuity" at b = 0 or c = 0 is unnecessary: the endpoint is itself an old-H atom.

**Not checked:**
- the face Gram constructions against any proved statement (none is proved here);
- the uses carried only by "all 1152 literal actions" in checker data;
- the allocation lemmas built from the atoms (H/chain resource allocation in T28, T36, T38 and T42; the T48 allocations).

**Attribution remark (not an error):** the T2 scalar note maps the old triple [11, 39, 17] to [9, 1, 22] by a stated row and column permutation. I checked that the two triples lie in one orbit (EXACT), not that particular group element.

## 6. Status per family

- **SOC (144)**: proved here, EXACT + PROOF-READ. The scalar theorem is from the T2 note; the bridge is mine.
- **old-H (288)**: proved here, EXACT + PROOF-READ. Bridge and increment argument as in the docs2 gapcheck, re-verified.
- **chain (576), = bad70 (PB)/(CQ)**: proved here, EXACT + PROOF-READ. The Gram identity was verified by my code from the supplement's certificate data; no PDF contains the proof.
- **four-coefficient (288)**: proved in the supplement (T5 four-budget note) and checked, EXACT + PROOF-READ.
- **derived sums** (old-H sum/Hsum, bad70 (H), coreE four-chain, H/chain and four/chain allocations): reduce to the four proved families. The allocation algebra is not audited here, except coreE and bad70 (H) by the first and second passes.
- **rank-two (144)**: NOT proved. Its Gram completion is stated in the PDF, but the cubic basis exists only in supplement checker data. NUMERIC only.
- **face Gram families** (T24 coherent/parallel/double-SOC, T25 partial-four, T44 Prop. 1-2, T46 core and six-ray, T41 bad17 coherent bound, T45, T47 merged-C/PCBE/IIIb/coherent-J, rational capacities): NOT proved, for the same reason. The 9 listed constructions are NUMERIC only; the rest are untested. Next step: reconstruct each face's coefficient-basis congruence, as done for bad70 in A5, and verify the stated PSD constructions against it.
