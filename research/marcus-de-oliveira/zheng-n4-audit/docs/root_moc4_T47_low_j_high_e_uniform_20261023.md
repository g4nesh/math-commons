# Audit of `root_moc4_T47_low_j_high_e_uniform_20261023`

Evidence level: EXACT for disc_p G = H6*H8 (123 and 260 terms), the full 1,289,451-term identity Z = R + 1280 P^6 V^16 S^2 Q^2 (28U^3+182U^2+399U+291)(T-1)^2 with exactly the four stated negative coefficients, R >= 0 (1,289,447 terms) with the monomial 81285120 T^2, and positive denominators; PROOF-READ for the homogenisation j = T/(T+Wj) (argument supplied); CITED for radial necessity, negative-Phi, the H6 >= 0 corollary and the Crat capacity.

PDF: `proofs/manin_transfer_20261021/root_moc4_T47_low_j_high_e_uniform_20261023.pdf` (2 pages). Auditor label C24-T47.

## Deviations and limits

- G (the radial polynomial) was rebuilt from the paper's definitions only: w(g,h) = sum of the source coefficients times the Appendix D vectors f_j plus h*kappa, C(w) = C0 + hQ, G^2 := -Q^2 det Hess C(w(g,-C0/Q))/4 (d = 1). The 14x14 determinant was computed by Kronecker substitution (degree bounds found by univariate interpolation) and every G was cross-checked against an independent exact Fraction determinant at random rational points. G is only defined up to a global sign; the sign was fixed once (positive c^2-coefficient at the all-ones point) and the document's sign claims were then reproduced.
- The chart clearings are my own (each variable cleared by its own degree, or by a common degree where the document says so). Where my term counts equal the document's this is noted; where they differ, the sign conclusion was re-proved for my clearing.
- Nothing in the supplement was executed; its JSON certificates were not used.
- Cited inputs (radial necessity G = 0 for spectrally realizable sources, the negative-Phi exclusion, and the pointwise capacities) were not audited here and are CITED.
- Falsification: the exact certificates prove the radial claims outright. A floating-point realizability test was run only on a few points (budget): on 2 random realizable T_65-stratum points, |G| relative to the absolute-term sum was <= 6e-12 (consistent with G = 0 being necessary), Phi > 0 and H6 > 0. No random realizable point fell in any hard region (too few samples to mean much).

## A. Statements

**Theorem 1** (contact set [[0, 4, 5, 6, 9, 12]]). Every realizable source on T_65 (d=1) with b>1, 0<a<1/b, p>b, q>=0, e>=b(b+1), 0<j<=1, c>0 satisfies Phi_c(O)>=0 for all O in SO(4). Proof: c>=Crat by cited capacity; for c<Crat, G - (p-b)H6 - M Phi > 0 with M>=0 contradicts G=0 given the cited Phi>=0, H6>=0.

## B. Dependencies

| kind | reference | used for | match checked |
|---|---|---|---|
| core-doc | root_moc4_T47_negative_Phi_20261022 / root_moc4_T47_full_radial (radial necessity: realizable source => C(w)=0 and det Hess C = 0 at w(g,-C0/Q), hence G = 0; rests on the secant-variety rank bound) | G = 0 at every spectrally realizable source | False |
| core-doc | root_moc4_T47_negative_Phi_20261022 (negative-Heron theorem: Phi < 0, d,e > 0 => no realizable lift) | exclude Phi < 0 | False |
| other-doc | final_referee_T47_negative_p_Heron_20261023 corollary: H6 >= 0 at realizable sources with a,c,d,e>0, p>=b | H6 >= 0 | False |
| common-lemma | final_referee_T46_six_ray_coherent_capacity_20261023 Corollary 3 (c >= Crat => support inequality) | capacity | False |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Contact set {0,4,5,6,9,12} is paper T_65 (identity element of G); its 8 generators are f2,f9,f1,f45,f5,f35,f0,f7 | identity | EXACT | lib generators/CONTACT_SETS |
| 2 | C(w(g,h)) = C0 + hQ with C0 = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} < 0 and Q as stated in the negative-Phi note | identity | EXACT | own flint computation from (A.2) and App. D |
| 3 | G: 965-term octic with det Hess C(w(g,-C0/Q)) = -4 d^2 G^2/Q^2 | identity | EXACT | own Kronecker determinant + 5-point independent Fraction det check |
| 4 | Realizable source => G = 0 (radial necessity), Q != 0 since C0 < 0 | citation | CITED | secant/Hessian-rank argument of the core notes, not audited here |
| 5 | H6 is the 123-term sextic factor of disc_p G = H6 H8 | identity | EXACT | flint factor |
| 6 | Crat algebra: N > D0^2 > 0 | inequality | EXACT | sympy + inspection |
| 7 | Z identity, four negative coefficients, R >= 0, strict T^2 monomial | identity | EXACT | own replay; counts identical (Z 1289451, R 1289447) |
| 8 | sign convention of H6 | citation | CITED | identity holds only for one sign of H6 (the sign returned by flint's factorisation relative to my normalised G); the document fixes the sign only by reference to companion arrays |
| 9 | j=1 wall by homogenisation | analytic-topological | PROOF-READ | argument supplied by auditor: Z has degree 2 in T, j enters only through n=T, d=1+T; replacing d by T+Wj is the standard homogenisation W^2 Z(T/W), which keeps coefficients |
| 10 | M >= 0, p-b > 0, H6>=0, Phi>=0 => G > 0, contradiction | inequality | PROOF-READ | - |

## D. Falsification

- 2 realizable points: H6 > 0 with the sign that makes the identity work (consistent with the cited corollary)

## E. Verdict

**confirmed-modulo-citations**. The 1.29M-term certificate is replayed exactly, including the four negative coefficients and the square; the H6 sign convention is not stated in the note (presentation); rests on the cited H6>=0, Phi>=0, radial necessity and Crat capacity.

Gaps:

- [presentation] p.1, definition of H6: H6 is defined only as 'the primitive sextic in disc_p G = H6 H8, with the sign fixed by the companion exact arrays'. Why: A polynomial factor is defined only up to sign; the identity needs one specific sign and the cited corollary must use the same one. The note does not state it. Believed true: True. Repair: State H6 explicitly or fix its sign by one coefficient; check the cited corollary uses that sign.

Unchecked:

- cited H6>=0 corollary, Crat capacity, necessity, negative-Phi

Replay: `replay/root_moc4_T47_low_j_high_e_uniform_20261023__check.py` (run `python -I -B root_moc4_T47_low_j_high_e_uniform_20261023__check.py ../lib` from `replay/`), output in `replay/root_moc4_T47_low_j_high_e_uniform_20261023__check_output.txt`.

