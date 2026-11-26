Evidence level: EXACT for the source, absorption, zero faces, a = 0 face, height formulas and the full radial identity (own code, 18 checks, 0 failures). PROOF-READ for the partition and the assembly logic. CITED for the three positive-a sector theorems (not audited here), the a = 0 face note (audited separately in this batch), the paper cases of (4), the four-coefficient atom and radial necessity. NUMERIC as support only.

# root_moc4_T46_whole_assembly_20261023 ("The support inequality on the T46 six-contact cone", 3 pp.)

Primary proof of paper case **T_47 = {1, 3, 5, 6, 7, 8}**. The supplement's label is historical T46. This is an assembly document, and most of its content is in the cited sector notes.

Replay: `replay/root_moc4_T46_whole_assembly_20261023__check.py`, output in `__check_output.txt`.

## Deviations and limits

1. The three sector theorems of table (6) are not in my batch. I read only their theorem statements, to confirm that they match (6). Their certificates (the 39 tables and the 1774 square terms) are UNCHECKED by me.
2. I recomputed G9 myself as the squared factor of the full 14x14 Hessian determinant of (A.2) on the 8-ray face. It has 1556 terms, which agrees with the document. I did not use the supplement's polynomial.

## A. Statements

- **Theorem 1 (whole, p. 1).** For x in X with P c(x) >= 0 vanishing on T = {1,3,5,6,7,8}, with extra contacts allowed: Phi_{c(x)}(O) >= 0 for all O in SO(4). The contact set is literally the paper's T_47. The statement is as strong as (*). Spectra are unrestricted, because the argument lives in coefficient space and uses only radial necessity, which holds for every x in X.

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| generators(T) = {f0..f6, f27, f28, f59} (Lemma 4.2) | core | EXACT (lib) | matches (1) |
| q=0, p=0, j=0, c=0, e=0, b=0 faces lead to historical T28, T35, T0, T38, T41, T6 = paper T_29, T_36, T_1, T_39, T_42, T_7 | paper cases | CITED (all Stage 1) | EXACT: each zero face of the canonical 8-ray form (with a present) contains a G-image of the paper set (group elements in the output) |
| a = 0 face = T ∪ {14} = Table 5 set T^(7)_26 | seven-contact | CITED: final_referee_T46_a0_uniform_20261023 (Stage 2; audited in this batch, confirmed-modulo-citations) | EXACT: the face equals T^(7)_26 |
| final_referee_T46_positive_a_small_or_balanced_20261023 (q,b <= 1 or q = b, any a >= 0) | branch | CITED, not audited | statement read; matches row 1 of (6) |
| root_moc4_T46_positive_a_mixed_uniform_20261023 (qb <= 1, max > 1) | branch | CITED, not audited | statement read; matches row 2 |
| final_referee_T46_positive_a_large_product_20261023 (qb > 1, q != b) | branch | CITED, not audited | statement read; matches row 3 |
| (5) ab(e - 1) >= a + b + e + 1 is safe | common lemma | CITED: four(4,6;2;28) in foundations so4-atoms-part2b G3 | roles (a: f4, b: f6, c: f2, d: f28) match the four-atom form ab(c - d) >= d^2(a + b + c + d) |
| six-ray / core coherent capacities | common lemma | CITED foundations G3 (used only inside the sector notes) | - |
| radial necessity | common lemma | CITED (audit item 5) | - |

T_47 is in Stage 3. Its inputs are in Stage 1, and T^(7)_26 is in Stage 2, so the order is acyclic.

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | rays (1); only f28 is non-minor | identity | EXACT | lib |
| 2 | (2a) f28 + f0 = f1 + f2 + f5 + f27 and (2b) f28 + f3 = f4 + f5 + f27 on all 24 values | identity | EXACT | - |
| 3 | absorption bookkeeping (q, c += u+v; p, e += u; a += v; d -= u+v) | identity | EXACT | symbolic, all 24 values |
| 4 | if f0 + f3 >= d, the cost is a minor sum; otherwise it has form (3) with d > 0 and the same cost vector | logic | PROOF-READ | the 8 rays are independent (EXACT), so (3) is unique |
| 5 | zero faces (4) | citation | EXACT + CITED | see B |
| 6 | a = 0 face is T ∪ {14} | identity | EXACT | - |
| 7 | C0 = -(q+1)L - jM < 0 and Q as printed | identity | EXACT | flint, paper cubic |
| 8 | det Hess C(g - (C0/Q) kappa) = -4 G9^2/Q^2, G9 with 1556 terms | identity | EXACT | Bareiss (230 s), then factorization: Q to the power 1 and G9 squared |
| 9 | Hessian rank at most twelve at a realizable lift | claim | GAP (presentation) | not justified in the document, and not used: only det = 0 is needed, and that follows from Hess C(u) v = 0 (radial necessity) |
| 10 | partition (6) of positive (q, b) | case split | PROOF-READ | exhaustive |
| 11 | each sector gives the support inequality | citation | CITED | not audited |

## D. Falsification

- NUMERIC necessity check, seed 7. At 10 matching vectors with gaps vanishing on T, |G9|/sum|terms| <= 2.7e-17. The random-coefficient control is >= 0.44.
- NUMERIC main statement. At 4 lib stratum points (seeds 0-3), the worst min Phi / max gap is -6.4e-17, which is rounding noise. No violation.

## E. Verdict

**confirmed-modulo-citations.** The assembly is correct, and every identity it states holds exactly. The proof of T_47 rests on three sector notes that I did not audit, on the a = 0 note (audited in this batch), on six Stage-1 cases, and on the foundations atoms.

Gaps:

- p. 2, "Hessian rank at most twelve": unjustified as stated, and unused. Only det Hess = 0 is needed, and that holds. Severity: presentation. I believe it is harmless.
