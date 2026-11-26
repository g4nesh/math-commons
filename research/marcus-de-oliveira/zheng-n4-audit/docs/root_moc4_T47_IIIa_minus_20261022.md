# Audit of `root_moc4_T47_IIIa_minus_20261022`

Evidence level: EXACT for Lemma 1, the Phi decomposition, the identities HG = alpha Phi + u0 + c u1 and KG = gamma Phi + c(v0 - c u0), and all seven coefficient tables (term counts identical to the document) with strictness on every allowed zero pattern; PROOF-READ for Lemma 1's case analysis and chart coverage; CITED for radial necessity and the negative-Phi theorem.

PDF: `proofs/manin_transfer_20261021/root_moc4_T47_IIIa_minus_20261022.pdf` (2 pages). Auditor label C24-T47.

## Deviations and limits

- G (the radial polynomial) was rebuilt from the paper's definitions only: w(g,h) = sum of the source coefficients times the Appendix D vectors f_j plus h*kappa, C(w) = C0 + hQ, G^2 := -Q^2 det Hess C(w(g,-C0/Q))/4 (d = 1). The 14x14 determinant was computed by Kronecker substitution (degree bounds found by univariate interpolation) and every G was cross-checked against an independent exact Fraction determinant at random rational points. G is only defined up to a global sign; the sign was fixed once (positive c^2-coefficient at the all-ones point) and the document's sign claims were then reproduced.
- The chart clearings are my own (each variable cleared by its own degree, or by a common degree where the document says so). Where my term counts equal the document's this is noted; where they differ, the sign conclusion was re-proved for my clearing.
- Nothing in the supplement was executed; its JSON certificates were not used.
- Cited inputs (radial necessity G = 0 for spectrally realizable sources, the negative-Phi exclusion, and the pointwise capacities) were not audited here and are CITED.
- Falsification: the exact certificates prove the radial claims outright. A floating-point realizability test was run only on a few points (budget): on 2 random realizable T_65-stratum points, |G| relative to the absolute-term sum was <= 6e-12 (consistent with G = 0 being necessary), Phi > 0 and H6 > 0. No random realizable point fell in any hard region (too few samples to mean much).

## A. Statements

**Lemma 1** (contact set []). g = alpha t^2+beta t+gamma, f = H t^2+J t+K, H>0, K>=0, l>0, u0=H gamma-alpha K, u1=H beta-alpha J, v0=K beta-gamma J. If alpha>=0 and u0+l u1>0 and either u0>0 or (gamma>0, v0>=0), then g(t)>0 for all 0<=t<l with f(t)>=0.

**Theorem 1** (contact set [[0, 4, 5, 6, 9, 12]]). Source g on T_65 (d=1): if 0<a<=1, 0<=b<=1, e>0, 0<=j<=e/(1+e), p,q>=0, 0<=c<1+1/e and Phi(c)>=0 then G(c)>0 (G up to a global sign). With radial necessity and the negative-Phi theorem: no realizable source in this sector. a=0 not claimed. Pure polynomial statement, so no spectral degeneracy issue.

## B. Dependencies

| kind | reference | used for | match checked |
|---|---|---|---|
| core-doc | root_moc4_T47_negative_Phi_20261022 / root_moc4_T47_full_radial (radial necessity: realizable source => C(w)=0 and det Hess C = 0 at w(g,-C0/Q), hence G = 0; rests on the secant-variety rank bound) | G = 0 at every spectrally realizable source | False |
| core-doc | root_moc4_T47_negative_Phi_20261022 (negative-Heron theorem: Phi < 0, d,e > 0 => no realizable lift) | exclude Phi < 0 | False |
| other-doc | C0 = -d{...} < 0 (own replay: formula EXACT) | Q != 0 at realizable points | True |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Contact set {0,4,5,6,9,12} is paper T_65 (identity element of G); its 8 generators are f2,f9,f1,f45,f5,f35,f0,f7 | identity | EXACT | lib generators/CONTACT_SETS |
| 2 | C(w(g,h)) = C0 + hQ with C0 = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} < 0 and Q as stated in the negative-Phi note | identity | EXACT | own flint computation from (A.2) and App. D |
| 3 | G: 965-term octic with det Hess C(w(g,-C0/Q)) = -4 d^2 G^2/Q^2 | identity | EXACT | own Kronecker determinant + 5-point independent Fraction det check |
| 4 | Realizable source => G = 0 (radial necessity), Q != 0 since C0 < 0 | citation | CITED | secant/Hessian-rank argument of the core notes, not audited here |
| 5 | Lemma 1 identities | identity | EXACT | sympy |
| 6 | Lemma 1 case analysis (endpoint, u0<0 forces K>0) | inequality | PROOF-READ | read line by line |
| 7 | Phi = Hc^2+Jc+K with H0, K0 as stated | identity | EXACT | flint |
| 8 | Eend = e u0 + (e+1) u1 is the endpoint condition at l=(e+1)/e | identity | PROOF-READ | e>0 |
| 9 | Charts (a,b)=(L/S,(L+M)/S) or swapped, j = EU/((1+E)(U+W)) cover the sector incl. walls | case split | PROOF-READ | explicit inverses L=a, M=b-a, N=1-b etc. |
| 10 | Tables: alpha 1439, gamma 3635, v0 33089, Eend 31893 (a<=b); alpha 1439, u0 27768, Eend 31893 (a>=b), all coefficients positive, strict markers | inequality | EXACT | own clearing; counts identical; strictness checked over all zero patterns allowed by E>0, L>0 (resp. L+M>0), U+W>0 |

## D. Falsification

- 2 realizable T_65 points: G ~ 0 (rel 6e-12), Phi>0

## E. Verdict

**confirmed-modulo-citations**. Polynomial claim G>0 on the IIIa-minus sector (with Phi>=0) replayed exactly (all 7 tables, identical term counts); the no-lift conclusion rests on the cited radial necessity and negative-Phi theorem.


Unchecked:

- radial necessity and negative-Phi theorem (cited)

Replay: `replay/root_moc4_T47_IIIa_minus_20261022__check.py` (run `python -I -B root_moc4_T47_IIIa_minus_20261022__check.py ../lib` from `replay/`), output in `replay/root_moc4_T47_IIIa_minus_20261022__check_output.txt`.

