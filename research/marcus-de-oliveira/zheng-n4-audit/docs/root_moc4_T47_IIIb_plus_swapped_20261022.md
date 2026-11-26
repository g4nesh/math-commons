# Audit of `root_moc4_T47_IIIb_plus_swapped_20261022`

Evidence level: EXACT for HG = alpha Phi + U0 + cU1, the chart algebra and all nine tables with term counts identical to the document (17604, 181025, 233360; 2373, 30345, 35745; 4386, 48680, 64615) and strictness; PROOF-READ for coverage; CITED for radial necessity and negative-Phi.

PDF: `proofs/manin_transfer_20261021/root_moc4_T47_IIIb_plus_swapped_20261022.pdf` (2 pages). Auditor label C24-T47.

## Deviations and limits

- G (the radial polynomial) was rebuilt from the paper's definitions only: w(g,h) = sum of the source coefficients times the Appendix D vectors f_j plus h*kappa, C(w) = C0 + hQ, G^2 := -Q^2 det Hess C(w(g,-C0/Q))/4 (d = 1). The 14x14 determinant was computed by Kronecker substitution (degree bounds found by univariate interpolation) and every G was cross-checked against an independent exact Fraction determinant at random rational points. G is only defined up to a global sign; the sign was fixed once (positive c^2-coefficient at the all-ones point) and the document's sign claims were then reproduced.
- The chart clearings are my own (each variable cleared by its own degree, or by a common degree where the document says so). Where my term counts equal the document's this is noted; where they differ, the sign conclusion was re-proved for my clearing.
- Nothing in the supplement was executed; its JSON certificates were not used.
- Cited inputs (radial necessity G = 0 for spectrally realizable sources, the negative-Phi exclusion, and the pointwise capacities) were not audited here and are CITED.
- Falsification: the exact certificates prove the radial claims outright. A floating-point realizability test was run only on a few points (budget): on 2 random realizable T_65-stratum points, |G| relative to the absolute-term sum was <= 6e-12 (consistent with G = 0 being necessary), Phi > 0 and H6 > 0. No random realizable point fell in any hard region (too few samples to mean much).

## A. Statements

**Theorem 1** (contact set [[0, 4, 5, 6, 9, 12]]). No root of G with Phi>=0 and a>1, 0<=b<1/a, p>0, q>=0, ea/(e+a)<j<=e, j<1/p+1/a, 0<=c<1/j. The wall j=ea/(e+a), a=1 and ab=1 are not claimed.

## B. Dependencies

| kind | reference | used for | match checked |
|---|---|---|---|
| core-doc | root_moc4_T47_negative_Phi_20261022 / root_moc4_T47_full_radial (radial necessity: realizable source => C(w)=0 and det Hess C = 0 at w(g,-C0/Q), hence G = 0; rests on the secant-variety rank bound) | G = 0 at every spectrally realizable source | False |
| core-doc | root_moc4_T47_negative_Phi_20261022 (negative-Heron theorem: Phi < 0, d,e > 0 => no realizable lift) | exclude Phi < 0 | False |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Contact set {0,4,5,6,9,12} is paper T_65 (identity element of G); its 8 generators are f2,f9,f1,f45,f5,f35,f0,f7 | identity | EXACT | lib generators/CONTACT_SETS |
| 2 | C(w(g,h)) = C0 + hQ with C0 = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} < 0 and Q as stated in the negative-Phi note | identity | EXACT | own flint computation from (A.2) and App. D |
| 3 | G: 965-term octic with det Hess C(w(g,-C0/Q)) = -4 d^2 G^2/Q^2 | identity | EXACT | own Kronecker determinant + 5-point independent Fraction det check |
| 4 | Realizable source => G = 0 (radial necessity), Q != 0 since C0 < 0 | citation | CITED | secant/Hessian-rank argument of the core notes, not audited here |
| 5 | three charts (p>=p*; p<p*, j<a; p<p*, j>=a) cover the region | case split | PROOF-READ | j + j^2/(a-j) = aj/(a-j); m + 1/p = a + 1/(p*P) for p = p*P/(1+P) (uses a - 1/a = 1/p*, EXACT) |
| 6 | alpha, U0, jU0+U1 tables positive with strict markers | inequality | EXACT | own replay, identical counts |

## D. Falsification

- shared 2-point realizability test

## E. Verdict

**confirmed-modulo-citations**. Radial exclusion on the IIIb-plus region replayed exactly (nine tables, identical term counts); rests on cited radial necessity and negative-Phi theorem.


Unchecked:

- radial necessity, negative-Phi

Replay: `replay/root_moc4_T47_IIIb_plus_swapped_20261022__check.py` (run `python -I -B root_moc4_T47_IIIb_plus_swapped_20261022__check.py ../lib` from `replay/`), output in `replay/root_moc4_T47_IIIb_plus_swapped_20261022__check_output.txt`.

