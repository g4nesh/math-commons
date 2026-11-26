# Audit of `root_moc4_T47_IIIa_plus_20261022`

Evidence level: EXACT for all sixteen tables (eight with term counts identical to the document; the j<1 low chart cleared differently, signs re-proved) and strictness; PROOF-READ for chart coverage; CITED for radial necessity, the negative-Phi theorem and Lemma 1 (replayed in the IIIa-minus audit).

PDF: `proofs/manin_transfer_20261021/root_moc4_T47_IIIa_plus_20261022.pdf` (2 pages). Auditor label C24-T47.

## Deviations and limits

- For the j<1 chart I cleared e = j + T^2E/((1+T)(1+E)) as a single fraction, which gives different term counts (alpha 1350 vs 1200; gamma 5292 vs 4508/4496; v0 26613/25794 vs 23322/22593; ESOC 24750 vs 22000). All my coefficients are positive, so the conclusion is re-proved independently. j>=1 counts are identical (1770, 5320/5275, 27618/27258, 30241).
- G (the radial polynomial) was rebuilt from the paper's definitions only: w(g,h) = sum of the source coefficients times the Appendix D vectors f_j plus h*kappa, C(w) = C0 + hQ, G^2 := -Q^2 det Hess C(w(g,-C0/Q))/4 (d = 1). The 14x14 determinant was computed by Kronecker substitution (degree bounds found by univariate interpolation) and every G was cross-checked against an independent exact Fraction determinant at random rational points. G is only defined up to a global sign; the sign was fixed once (positive c^2-coefficient at the all-ones point) and the document's sign claims were then reproduced.
- The chart clearings are my own (each variable cleared by its own degree, or by a common degree where the document says so). Where my term counts equal the document's this is noted; where they differ, the sign conclusion was re-proved for my clearing.
- Nothing in the supplement was executed; its JSON certificates were not used.
- Cited inputs (radial necessity G = 0 for spectrally realizable sources, the negative-Phi exclusion, and the pointwise capacities) were not audited here and are CITED.
- Falsification: the exact certificates prove the radial claims outright. A floating-point realizability test was run only on a few points (budget): on 2 random realizable T_65-stratum points, |G| relative to the absolute-term sum was <= 6e-12 (consistent with G = 0 being necessary), Phi > 0 and H6 > 0. No random realizable point fell in any hard region (too few samples to mean much).

## A. Statements

**Theorem 1** (contact set [[0, 4, 5, 6, 9, 12]]). Source g on T_65, d=1: if 0<a<=1, 0<=b<=1, p>0, q>=0, 0<j<=e, j>e/(1+e) (the equality j=e/(1+e) is delegated to IIIa-minus although the displayed hypothesis reads j>=e/(1+e)), j<1+1/p, 0<=c<1/j and Phi(c)>=0 then G(c)>0. p=0 and a=0 excluded (delegated).

## B. Dependencies

| kind | reference | used for | match checked |
|---|---|---|---|
| core-doc | root_moc4_T47_negative_Phi_20261022 / root_moc4_T47_full_radial (radial necessity: realizable source => C(w)=0 and det Hess C = 0 at w(g,-C0/Q), hence G = 0; rests on the secant-variety rank bound) | G = 0 at every spectrally realizable source | False |
| core-doc | root_moc4_T47_negative_Phi_20261022 (negative-Heron theorem: Phi < 0, d,e > 0 => no realizable lift) | exclude Phi < 0 | False |
| other-doc | C0 = -d{...} < 0 (own replay: formula EXACT) | Q != 0 at realizable points | True |
| other-doc | root_moc4_T47_IIIa_minus_20261022 Lemma 1 and the threshold j=e/(1+e) | comparison lemma; threshold | True |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Contact set {0,4,5,6,9,12} is paper T_65 (identity element of G); its 8 generators are f2,f9,f1,f45,f5,f35,f0,f7 | identity | EXACT | lib generators/CONTACT_SETS |
| 2 | C(w(g,h)) = C0 + hQ with C0 = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} < 0 and Q as stated in the negative-Phi note | identity | EXACT | own flint computation from (A.2) and App. D |
| 3 | G: 965-term octic with det Hess C(w(g,-C0/Q)) = -4 d^2 G^2/Q^2 | identity | EXACT | own Kronecker determinant + 5-point independent Fraction det check |
| 4 | Realizable source => G = 0 (radial necessity), Q != 0 since C0 < 0 | citation | CITED | secant/Hessian-rank argument of the core notes, not audited here |
| 5 | ESOC = j u0 + u1 is j times the affine remainder at c=1/j | identity | PROOF-READ | - |
| 6 | low chart covers j<=e<j/(1-j), j<1; high chart covers j in [1,1+1/p), e>=j | case split | EXACT | sympy identity e-j = (j^2/(1-j)) E/(1+E); limit of j |
| 7 | alpha>=0, gamma>0, v0>=0, ESOC>0 on all four charts | inequality | EXACT | own clearing; all coefficients positive; strictness over zero patterns |

## D. Falsification

- shared 2-point realizability test (G~0)

## E. Verdict

**confirmed-modulo-citations**. G>0 on the IIIa-plus sector replayed exactly on all four charts (own clearing); one presentation slip (j>=e/(1+e) displayed, equality actually delegated to IIIa-minus); conclusion rests on cited radial necessity and negative-Phi theorem.

Gaps:

- [presentation] p.1, Theorem 1 hypotheses: Hypothesis written j >= e/(1+e) but the proof covers only j > e/(1+e). Why: The equality case is deferred to IIIa-minus (stated in the next sentence); the theorem as displayed is stronger than what this note proves on its own. Believed true: True. Repair: State j > e/(1+e); IIIa-minus covers equality (its c-range 1+1/e equals 1/j there).

Unchecked:

- cited necessity/negative-Phi

Replay: `replay/root_moc4_T47_IIIa_plus_20261022__check.py` (run `python -I -B root_moc4_T47_IIIa_plus_20261022__check.py ../lib` from `replay/`), output in `replay/root_moc4_T47_IIIa_plus_20261022__check_output.txt`.

