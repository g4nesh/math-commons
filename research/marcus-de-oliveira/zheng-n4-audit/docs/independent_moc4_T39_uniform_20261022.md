# Audit of independent_moc4_T39_uniform_20261022 (primary proof of paper case T_40)

Evidence level: EXACT for every own step: the source cone, the three identities, the Appendix A completeness formulas and both rewritings, the transpose map, the capacity algebra and the radial C/Q formulas. PROOF-READ for the reduction logic and the four-region cover. CITED for the three chamber theorems, the SOC and old-H atoms, and Theorem 1.2. NUMERIC for the stratum test.

Auditor label C15-T38-T40. Document: proofs/manin_transfer_20261021/independent_moc4_T39_uniform_20261022.pdf, 4 pages. Replay: `replay/independent_moc4_T39_uniform_20261022__check.py` (0 failed checks).

## Deviations and limits

- This document is an assembly. Its Section 3 summarizes three chamber theorems proved in separate notes: low-pj, high-j and high-p. The high-j note is audited in this batch as independent_moc4_T39_high_j_uniform_20261022 (verdict could-not-check, no gap). The low-pj and high-p notes are supplement PDFs that I did not audit; they are CITED.
- The "mechanism summary" of Section 3 was read but not used as evidence.

## A. Statements

Theorem 1 (p. 1). For every x in X with P c(x) >= 0, contact set containing T_40 = {0,1,5,9,14,18}, and every O in O(4) (and so every unitary), Phi_{c(x)}(O) >= 0. Repeated or coincident spectral values, zero resources and extra contacts are included. The contact set equals paper T_40 literally (EXACT). Uniform disjunction (Section 4): every realizable canonical cost with d > 0 has pj >= d^2, or a+e > 0 and j >= d + d^2/(a+e).

## B. Dependencies

| kind | ref | used for | status |
|---|---|---|---|
| branch | independent_moc4_T39_high_j_uniform_20261022 | region j >= 1, pj < 1 | audited in this batch: could-not-check, no gap |
| branch | independent_moc4_T39_low_pj_uniform_20261022 | region p < 1, j < 1 | CITED, not audited |
| branch | independent_moc4_T39_high_p_uniform_20261022 | region p >= 1, pj < 1 | CITED, not audited |
| core-doc | SOC family (3) and old-H atoms on (f29,f12,f35), (f2,f12,f35) | capacities | CITED; NUMERIC spot check |
| common-lemma | Theorem 1.2 | U(4) | CITED |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | contact nullspace has dimension 8; generators {1,2,5,9,11,12,28,29,33,35}; non-minor 28, 35 | identity | EXACT |
| 2 | identities (1): f28+f12 = f1+f2+f11+f29; f35+f11 = f1+f5+f12+f33; f28+f35 = 2f1+f2+f5+f29+f33 | identity | EXACT |
| 3 | the cancellation procedure gives source g35 or g28 with nonnegative coefficients, or a minor sum | logic | PROOF-READ |
| 4 | App. A: (q,p,j,c,e,d,b,a) = (g11,g17,g6,g7,g15,g3-g12,g12,g4) on the nullspace; g3 = b+d, g20 = q+d, g23 = p+d, g8 = a+j+d, g19 = e+j+d, g22 = p+j+2d | identity | EXACT |
| 5 | App. A: for j >= m the minor-sum rewriting, and for j < m the second-source rewriting with coefficients (e-δ, p-m-δ, j, c, q-m, δ, a-δ, b-m); nonnegativity follows from step 4 | identity + inequality | EXACT + PROOF-READ |
| 6 | the inverse-permutation action (swap the spectra) preserves T_40 and sends the g28 source to g35 term by term; realizability is preserved | identity | EXACT (lib: G maps matching vectors to matching vectors) |
| 7 | capacity (4): splitting d in the ratio a : e gives threshold d + d^2/(a+e) | identity | EXACT (sympy) |
| 8 | C = -(ABc + Nj) and the -Q formula (d = 1) | identity | EXACT |
| 9 | radial necessity (F(w) = 0, det Hess = 0, real case via Euler) | analytic | PROOF-READ |
| 10 | four-region cover exhaustive; a+e = 0 is impossible in region 2 by the high-j theorem | case split | PROOF-READ |
| 11 | chamber exclusions | citation | CITED |

## D. Falsification

- NUMERIC: 12 stratum points with contact set exactly T_40 (seeds 0-11, margin 1e-2). 7 decompose to minor sums. The other 5 satisfy both pj >= 1 and the H-sum capacity. This is weak evidence, because the sampler does not reach the hard chambers.
- NUMERIC: min Phi >= -1e-12 at the SOC threshold p = j = 1, and at the H-sum threshold a = e = 1/2, j = 2 (30 restarts).

## E. Verdict

confirmed-modulo-citations. Every step of this document is EXACT or PROOF-READ. The proof rests on three CITED chamber theorems: low-pj and high-p not audited, and high-j audited here without a replay of its certificates. It also rests on the SOC and old-H atoms and on Theorem 1.2.
