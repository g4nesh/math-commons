# Findings of the audit so far

Dated 26 November 2026. Task `mdo-zheng-n4-audit`. Evidence levels are those of [METHOD.md](METHOD.md).

The architecture of the proof is Zheng's: the reduction of unitary to orthogonal matrices, the contact reduction on the eight-point matching cone, the finite classification, and the case-by-case treatment of ninety contact sets. This audit checks that architecture and its supplement. Where the public material leaves a step unproved, the audit says so and supplies a proof when it can.

## Bottom line

**No error has been found in any part examined.** The conjecture for $n=4$ is **not yet established by this audit as a whole**: 98 of the 208 supplementary proof documents still await a first-pass report, five documents have large polynomial certificates that were not replayed, and no case has yet received its final verdict.

| Part of the proof | Result | Evidence |
| --- | --- | --- |
| Contact reduction (paper Sections 3, 4.2–4.3, 6) | Holds. The ninety contact inequalities, Corollary 2.2 and the classification imply the conjecture for $n=4$. | PROOF-READ twice, EXACT finite identities by two independent implementations. Rests on classical citations. [Report](reduction/contact-reduction.md) |
| 82 facets, 214 orbits, 124/90 split, Appendix E | Hold. | EXACT by three independent methods. [Report](reduction/finite-classification.md) |
| Theorem 1.2, unitary to two orthogonal matrices | Holds. | PROOF-READ by two referees who did not see each other's work, one with an own proof by theta functions. NUMERIC to 40 digits. [First](foundations/theorem-1-2.md) · [second](foundations/theorem-1-2-second.md). A different, elementary proof by another owner is in [`../unitary-factorization/`](../unitary-factorization/README.md) (Theorem F and Corollary F1 there; not yet independently reviewed). |
| Appendix B, coincident eigenvalues | Holds as written. | PROOF-READ, EXACT replays. [Report](foundations/appendix-b.md) |
| Appendix C, one spectrum on a circle or line | No error found. One step was only sketched in the paper; the audit supplies a proof. See finding 3. | [Part 1](foundations/appendix-c-part1.md) · [part 2](foundations/appendix-c-part2.md) · [open items](foundations/appendix-c-open-items.md) |
| Universal inequalities on SO(4) used by the cases | All that were examined hold. Several are not proved in any public document; the audit supplies proofs. See finding 2. | [Inventory and proofs](foundations/so4-atoms.md) · [independent re-verification](foundations/so4-atoms-verification.md) · [part 2a](foundations/so4-atoms-part2a.md) · [part 2b](foundations/so4-atoms-part2b.md) |
| Three-plus-three circle theorem | Holds. Its proof note is absent from the public supplement; the audit supplies a proof. See finding 1. | [Proof](foundations/three-plus-three-circle.md) · [independent re-verification](foundations/three-plus-three-circle-verification.md) |
| The 208 supplementary proof documents | 110 first-passed here: 101 confirmed or confirmed modulo citations, 4 with a step flagged and then settled (finding 4), 5 with certificates not replayed (finding 5). 92 are with cloud workers. 6 parallel write-ups are deferred. | [Reports](docs/) · [second pass](docs2/) · [manifest](manifest.json) |

## Finding 1. A theorem cited by twelve cases has no public proof; the audit supplies one

Twelve case proofs ($T_4$, $T_8$, $T_9$, $T_{10}$, $T_{12}$, $T_{13}$, $T_{15}$–$T_{20}$) invoke, as previously proved, a three-plus-three generalized-circle theorem: if three eigenvalues of one matrix and three of the negatives of the eigenvalues of the other lie on one circle or line, the support inequality holds whatever the remaining two eigenvalues are. This is not the paper's Theorem C.1 and it is not stated in the paper. The supplement's archive manifest lists a proof note for it (`final_referee_six_circle_profile_20261022.tex`). That note is not among the 208 public PDFs. A check program and two data files are public.

The audit reconstructed the statement from the public program and data, read but not executed, and proved it. After a Möbius normalization the gap vector lies in a real family of dimension at most five. Every nonnegative member of that family is a nonnegative combination of the 34 minor vectors and twelve further vectors, each a symmetry image of the inequality of finding 2(b). The finite part was replayed with own exact code: an own double-description enumeration (72 dual vertices with minors alone, 86 with the twelve vectors) and all 158 polynomial decompositions. Degenerate configurations follow by a closure argument. A second auditor who did not write the proof re-derived it and certified it, with 3,900 exact rational configurations and no failure. The theorem as proved covers all twelve uses.

## Finding 2. Universal inequalities on SO(4) that the cases cite

The case proofs end with inequalities on $\mathrm{SO}(4)$ that involve no spectra. They are stronger than the conjecture restricted to a case, so they could be false, and many case documents use them without saying where they are proved.

| Family | Orbit under the group of order 1152 | Where the public material proves it | Status here |
| --- | --- | --- | --- |
| (a) SOC: $\Phi_k+2\lvert x_ix_j\rvert\geq0$ | 144 | Scalar form only. The bridge to Leibniz sums is not in any PDF. | Proved here. Re-derived independently. |
| (b) old-H: $\Phi_c+\Phi_k+2\lvert x_tx_c\rvert\geq0$ | 288 | Scalar form only. | Proved here. Re-derived independently. |
| (c) chain: $\Phi_k+2\lvert x_c\rvert(\lvert x_l\rvert+\lvert x_m\rvert)\geq0$ | 576 | No PDF proof. A certificate exists in the supplement's data. | Proved here from that certificate, re-verified exactly with own code. An independent proof that avoids the certificate is not complete. |
| (d) four-coefficient | 288 | Proved in the supplement's four-budget note. | Checked. Re-derived independently by a different argument. |
| (e) rank-two | 144 | Stated; its Gram basis is only in the supplement's data. | Proved here. |
| (f) Gram constructions on the faces of historical $T_{24}$, $T_{25}$, $T_{44}$, $T_{46}$, $T_{47}$ | face-specific | As (e). | Proved here or checked in the document. Three sub-statements for $T_{47}$ are UNCHECKED. |
| (g) Allocation lemmas and the remaining capacities | — | — | Reduce exactly to the families above. |

Here $\Phi_k$ is the Leibniz sum of a non-minor facet vector and $x_j$ the $2\times2$ minor whose square is the Leibniz sum of a minor vector. Numerical minimization over $\mathrm{SO}(4)$ found no negative value for any family, and the same test detects deliberately weakened variants. About 60 cited uses were compared with the proved orbits; none uses a stronger form than was proved.

Two remarks in the first inventory were wrong and are corrected in place: the chain inequality does not survive halving its compensating term, and equality in it is attained at a nondegenerate point. Neither affects a proof.

## Finding 3. Appendix C

Theorem C.1 is the conjecture for $n=4$ when one spectrum lies on a circle or line. Two referees found no error. Proposition C.4 (the dual variety of the Gale quotient is its secant variety) was only sketched in the paper. The audit supplies a complete proof, and with it Corollary C.5. After this work Theorem C.1 rests on Theorem 1.2, one rank statement at stable fixed points that is implied by a cited theorem of Howard, Millson, Snowden and Vakil, projective biduality, and standard facts of geometric invariant theory. The use of Lam's double-dimer theorem was replaced by explicit weights.

Theorem C.1 is indispensable only for $T_7$, $T_{11}$, $T_{24}$ and one auxiliary seven-contact set, on which 51 cases depend. Repeated eigenvalues elsewhere are covered by Appendix B. Another owner's result in [`../unitary-factorization/`](../unitary-factorization/README.md) (Theorem P4 and Corollary P4′, not yet independently reviewed) reaches the same conclusion by a different route.

## Finding 4. Steps flagged in case documents and settled by the second pass

| Document | Step | Outcome |
| --- | --- | --- |
| `blind_moc4_bad70_full_safe_disjunction` ($T_{34}$, the paper's representative of group V) | Sufficient condition (H), p. 1 | One step not proved in the document; proof supplied by the audit; statement holds. [Second pass](docs2/blind_moc4_bad70_full_safe_disjunction_20261022__gapcheck.md) |
| `blind_moc4_T44_coherent_chain_capacity` | Gram representation in Propositions 1 and 2 | The document cites only a companion program. Proved by the audit with exact code. |
| `blind_moc4_T44_conic_branch_exclusion` | Two polynomial identities, pp. 1–2 | The polynomials are not printed. Rebuilt from the paper's definitions and verified exactly. |
| `root_moc4_T10_B_branch` | Reality of four chart parameters, p. 1 | Presentation only. Supported numerically; two of the four verified exactly. |

Seven documents were audited a second time from scratch by auditors who did not see the first report. All seven second verdicts agree with the first.

## Finding 5. What is not yet checked

- **Large certificates not replayed.** Five documents for $T_{38}$, $T_{39}$ and $T_{40}$ (`root_moc4_T37_uniform`, `blind_moc4_T38_uniform`, `independent_moc4_T39_high_j_uniform`, `…_high_p_uniform`, `…_low_pj_uniform`) rest on coefficient tables and Sturm sequences that the first pass tested only at sampled points. Their logic was read and found sound.
- **92 documents with cloud workers.** No report yet. See [CLOUD-WORKER.md](CLOUD-WORKER.md).
- **Six parallel write-ups** of the paper's Section 2 and Appendix C are deferred.
- **Dependency graph.** The supplement's graph is acyclic and agrees with the paper's Table 6, and all ten conditional theorems are discharged at the level of statements. Seven citations found in the texts are missing from the supplement's graph, and one possible two-way dependence between $T_8$ and a seven-contact set is not yet resolved. [Map](deps/dependency-map.md) · [conditional theorems](deps/conditional-tracing.md)
- **Independent review.** Nothing here has been reviewed by another owner. Nine requests for independent replays are open in the ledger.
- **Machine checking.** Nothing here is formalized.

## What a confirmed proof would still need

Every case document confirmed with its certificates replayed; the cited inputs named above either proved or accepted as citations; an assembled verdict per case along the dependency order; and an independent review by another owner. A formal verification would in addition need the algebraic geometry of Section 2 or the elementary factorization that replaces it, the potential theory of Lemma 3.5, and the polynomial certificates in a proof assistant.
