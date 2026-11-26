Evidence level: EXACT for every displayed identity, factorization, resultant, Newton-polygon and series claim, the normal cone and the boundary-tie reductions; PROOF-READ for the puncture analysis, the maximum principle, the constant exception and the Rouché perturbation; CITED for T_37, T_45, T_57, T_65, Appendix B (B.1, B.2, B.3), Lemma 4.2 / Prop. 4.1, (3.11) and Birkhoff; NUMERIC only for a brief falsification search.

# Audit of `final_referee_inventory85_native_uniform_20261023` (T_86, companion version)

PDF: `proofs/manin_transfer_20261021/final_referee_inventory85_native_uniform_20261023.pdf`, 7 pages. Title: "Two spectral charts and compensation at one collision divisor". Auditor label C43-T86.

Replay: `replay/final_referee_inventory85_native_uniform_20261023__check.py`. It runs the shared own-code replay `replay/root_inventory85_native_uniform_corrected_20261023__check.py` unchanged. Output: `..._check_output.txt`, 58 EXACT checks, all pass. The numeric section is in the root document's output.

## Deviations and limits

- The document is not the catalogue's primary proof of T_86. It is an earlier or alternative text of the same proof, listed as used for T_86. I compared it line by line with `root_inventory85_native_uniform_corrected_20261023` using the text extracts and the page images of p. 4. Every display, constant and equation number is identical except the one difference noted under E. The full step table is therefore the one in `docs/root_inventory85_native_uniform_corrected_20261023.md` (C1-C23), and every check there applies verbatim to this document.
- Same limits as the root report: the cited cases and Appendix B were not re-audited, and the supplement's verifier was not used.

## A. Statements

**Theorem 1 (p. 1).** For T = {4,5,8,9,12,13} = paper T_86 (check A1, identity group element): every x ∈ X with P c(x) ≥ 0 and (P c(x))_π = 0 for π ∈ T satisfies Φ_{c(x)}(O) ≥ 0 for every real orthogonal O. The statement covers all complex spectra; degenerate ones are handled at the outset by the within-spectrum and shared-zero laws. This is the same statement as the root version, at least as strong as (*) for T_86. The last line says "The conclusion concerns exactly the contact set (1)".

## B. Dependencies

Identical to the root report, table B1-B9. The text here speaks of "the list of 60 previously established contact-case inequalities". My code finds that every boundary tie embeds an earlier-stage case: j = 2, 3 → T_37; j = 6, 7, 16, 17 → T_65; j = 10, 11 → T_57; j = 14, 15, 18, 19 → T_45. All are in stages 1-7 of Table 6, and T_86 is in stage 8. The single-collision lemma use matches Corollary B.3 in form, lifts and cofactor (root report B5).

## C. Steps

The steps are root report C1-C23, all EXACT, PROOF-READ or CITED, with one exception: in C18, display (19) of this version is wrong as printed (below).

## D. Falsification

The numeric tests are shared with the root report. Over 8 stratum points, the worst min over SO(4) of Φ is −4.3e-13 (round-off). The chart round trip on 400 points has relative error ≤ 7.3e-12.

## E. Verdict

**gap-found** (one presentation gap only; the proof is otherwise identical to the root version, which is confirmed modulo citations).

- G1 (presentation), p. 4, display (19). The rescaled C-configuration is printed as C' = (∞, N/y, N, M). Multiplying C = (∞, v, x, t) by L gives (∞, vL, xL, tL) = (∞, N, N/y, M) by (8), so C1' and C2' are swapped (check L5; the root version prints the correct order). The denominator (20) and the conclusion are unaffected. In the correct order, [A0', C2'] = N with the lift (N, y), and [C1', C3'] = M − N = K. These give exactly the factors N and K of (20), so (20) was computed with the correct order. In the printed order, the second factor would be M − N/y instead. The claim is true and the repair is to swap the two entries.

Citations as in the root report: T_37, T_45, T_57, T_65 (or the alternatives listed); Appendix B (established); Lemma 4.2 / Prop. 4.1 and (3.11) (paper, lib-replayed); Birkhoff.
