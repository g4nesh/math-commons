# Audit of `root_moc4_T47_low_j_low_e_uniform_20261023`

Evidence level: CITED for every substantive step (radial necessity, negative-Phi, Crat capacity, PCBE capacity, PCBE q-monotonicity note, q=0 endpoint theorem); EXACT only for the algebra turning PCBE failure into the chart conditions, and the setup (contact set, G, C0).

PDF: `proofs/manin_transfer_20261021/root_moc4_T47_low_j_low_e_uniform_20261023.pdf` (2 pages). Auditor label C24-T47.

## Deviations and limits

- G (the radial polynomial) was rebuilt from the paper's definitions only: w(g,h) = sum of the source coefficients times the Appendix D vectors f_j plus h*kappa, C(w) = C0 + hQ, G^2 := -Q^2 det Hess C(w(g,-C0/Q))/4 (d = 1). The 14x14 determinant was computed by Kronecker substitution (degree bounds found by univariate interpolation) and every G was cross-checked against an independent exact Fraction determinant at random rational points. G is only defined up to a global sign; the sign was fixed once (positive c^2-coefficient at the all-ones point) and the document's sign claims were then reproduced.
- The chart clearings are my own (each variable cleared by its own degree, or by a common degree where the document says so). Where my term counts equal the document's this is noted; where they differ, the sign conclusion was re-proved for my clearing.
- Nothing in the supplement was executed; its JSON certificates were not used.
- Cited inputs (radial necessity G = 0 for spectrally realizable sources, the negative-Phi exclusion, and the pointwise capacities) were not audited here and are CITED.
- Falsification: the exact certificates prove the radial claims outright. A floating-point realizability test was run only on a few points (budget): on 2 random realizable T_65-stratum points, |G| relative to the absolute-term sum was <= 6e-12 (consistent with G = 0 being necessary), Phi > 0 and H6 > 0. No random realizable point fell in any hard region (too few samples to mean much).

## A. Statements

**Theorem 1** (contact set [[0, 4, 5, 6, 9, 12]]). Every realizable source on T_65 (d=1) with b>1, 0<a<1/b, p>=b, q>=0, 0<j<=1, j<=eb/(e+b), 0<e<b(b+1), c>=0 satisfies the support inequality.

## B. Dependencies

| kind | reference | used for | match checked |
|---|---|---|---|
| core-doc | root_moc4_T47_negative_Phi_20261022 / root_moc4_T47_full_radial (radial necessity: realizable source => C(w)=0 and det Hess C = 0 at w(g,-C0/Q), hence G = 0; rests on the secant-variety rank bound) | G = 0 at every spectrally realizable source | False |
| core-doc | root_moc4_T47_negative_Phi_20261022 (negative-Heron theorem: Phi < 0, d,e > 0 => no realizable lift) | exclude Phi < 0 | False |
| other-doc | root_moc4_T47_zero_q_rational_capacity_20261023 (Crat) | c >= Crat | False |
| other-doc | PCBE capacity with L=b | capacity | False |
| other-doc | independent_moc4_T47_PCBE_low_e_q_monotone_20261023 | G(q) >= G(0) on the chart with Phi>=0 | False |
| other-doc | q=0 endpoint theorem (not named by file) | G(0) > 0 when c < Crat | False |
| other-doc | closed-source minor boundary theorem (c=0) | c=0 | False |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Contact set {0,4,5,6,9,12} is paper T_65 (identity element of G); its 8 generators are f2,f9,f1,f45,f5,f35,f0,f7 | identity | EXACT | lib generators/CONTACT_SETS |
| 2 | C(w(g,h)) = C0 + hQ with C0 = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} < 0 and Q as stated in the negative-Phi note | identity | EXACT | own flint computation from (A.2) and App. D |
| 3 | G: 965-term octic with det Hess C(w(g,-C0/Q)) = -4 d^2 G^2/Q^2 | identity | EXACT | own Kronecker determinant + 5-point independent Fraction det check |
| 4 | Realizable source => G = 0 (radial necessity), Q != 0 since C0 < 0 | citation | CITED | secant/Hessian-rank argument of the core notes, not audited here |
| 5 | PCBE failure <=> 0<=s<1 and j < (b+1)(1-s)^2/(b+2s+(b+1)s^2/e), s = sqrt(c/B) | identity | EXACT | sympy identity + monotonicity (PROOF-READ) |
| 6 | G increasing in q on the chart when Phi>=0 | citation | CITED | - |
| 7 | G(0) > 0 under c<Crat | citation | CITED | endpoint theorem not identified by file name |
| 8 | logic: G(q) >= G(0) > 0 contradicts G = 0 | inequality | PROOF-READ | valid provided Crat does not depend on q (true for the Crat formula used in the high-e note; the zero-q note's Crat not read) |

## D. Falsification

- shared 2-point realizability test

## E. Verdict

**confirmed-modulo-citations**. An assembly of cited results (Crat, PCBE capacity, PCBE q-monotonicity, an unnamed q=0 endpoint theorem); the only in-document algebra (PCBE failure conditions) checks exactly; the endpoint dependency is unnamed (presentation).

Gaps:

- [presentation] p.1-2, 'The independent algebraic endpoint theorem': G(0) > 0 under c < Crat is cited without naming the document. Why: An unnamed dependency cannot be resolved; I could not check that its region equals the PCBE-chart region used here or that Crat is q-independent. Believed true: True. Repair: Name the endpoint document and state its hypotheses.

Unchecked:

- all cited inputs; the endpoint and monotonicity certificates were not replayed (budget)

Replay: `replay/root_moc4_T47_low_j_low_e_uniform_20261023__check.py` (run `python -I -B root_moc4_T47_low_j_low_e_uniform_20261023__check.py ../lib` from `replay/`), output in `replay/root_moc4_T47_low_j_low_e_uniform_20261023__check_output.txt`.

