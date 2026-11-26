# Audit of independent_moc4_T44_native_harmonic_uniform_20261023

Evidence level: EXACT for the source, the seventh-tie resolution and every displayed chart identity, resultant and puncture projection; PROOF-READ for Lemmas 2-3 and the harmonic completion; CITED for the earlier contact cases, Appendix B / one-circle collisions, the chart's coverage and Birkhoff.

Auditor batch: C18-T43-T45. Source: `proofs/manin_transfer_20261021/independent_moc4_T44_native_harmonic_uniform_20261023.pdf` (4 pages). Replay: `replay/independent_moc4_T44_native_harmonic_uniform_20261023__check.py` (output in `replay/independent_moc4_T44_native_harmonic_uniform_20261023__check_output.txt`).

## Deviations and limits

- Batch budget about one hour for five documents; depth was put on the three primary proofs.
- Nothing from the supplement was executed or imported; no supplement JSON was loaded. All polynomials were rebuilt from the paper's definitions with the shared library (cubic (A.2), P, f_j of Appendix D).
- The large stored certificates (1000+-term radial factors, Bernstein tables) were NOT re-expanded symbolically. Instead the radial factor was reconstructed EXACTLY at sample points: G^2 = -Q^2 det Hess C(g - (C/Q) kappa)/4 by exact 14x14 rational determinants, interpolated along lines and square-rooted exactly; sign fixed by the d^9 coefficient = +1. Polynomial identities checked this way at finitely many random rational points are evidence (exact arithmetic, probabilistic), not proofs; they are labelled NUMERIC (exact-at-points).
- Library value check: the library cubic gives C, Q and the base normalisation G(d=1, rest 0) = 1 exactly as the documents state (constant factor 1); no disagreement with the library was found.

## A. Statements

- **main (Theorem 1)** (p. 1) - paper case(s) T_45: For every x in X with P c(x) >= 0 and contacts containing {1,2,3,7,11,12}, Phi_{c(x)}(O) >= 0 for all O in O(4). Degenerate spectra (collisions) and extra contacts are reduced to Appendix B / one-circle theorem and to earlier six-sets; the curve argument is run on eight-distinct, exactly-six-contact configurations.
- **Lemma 2** (p. 2): If m(tau) = f(tau) kappa + v(tau) in Sing C with f divergent and v holomorphic at 0, then H_kappa v(0) = 0; if -eta Im P v(0) >= 0 then its Leibniz sum is >= 0 on O(4).
- **Lemma 3** (p. 3): An eight-distinct counterexample with exactly these six contacts can be perturbed (Rouche on H11 - h) to one with generic real (r,s,u,h), keeping strict gaps and Phi_O < 0.

The contact set equals the paper's Appendix E set literally (no group element needed), and the generator list equals lib `generators(T)`. The statement is non-strict, covers all O in SO(4) (O(4) equivalently) and all spectra; it is at least as strong as (*).

## B. Dependencies

- [paper-case] T_7 (stage 1) - used for: seventh tie 0. Statement match checked: yes.
- [paper-case] T_36 (stage 1) - used for: seventh ties 6 and 9. Statement match checked: yes.
- [paper-case] T_42 (stage 1) - used for: seventh tie 18. Statement match checked: yes.
- [paper-case] T_40 (stage 1) - used for: seventh tie 20. Statement match checked: yes.
- [paper-case] T_46 (stage 2) - used for: seventh tie 21. Statement match checked: yes.
- [paper-case] T_29 (stage 1) - used for: seventh tie 23. Statement match checked: yes.
- [common-lemma] Appendix B (Thm B.1, Lemma B.2) and one-circle theorem (App. C) - used for: collision points on the boundary of Omega. Statement match checked: no (CITED).
- [other-doc] earlier 'original chart' document: every eight-distinct source with these contacts is in the chart with real nonzero r,s,u - used for: chart coverage. Statement match checked: no (CITED).
- [external] Birkhoff-von Neumann - used for: Lemma 2 (singleton span). Statement match checked: yes.
- [core-doc] paper (3.10)-(3.11): ker H_kappa = E (library-replayed) - used for: Lemma 2. Statement match checked: yes.
- [core-doc] Theorem 1.2 - used for: unitary extension. Statement match checked: no (CITED).

## C. Steps

| step | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | generators {0,1,2,7,8,12,15,18,53}; f18+f0=f2+f1+f7+f15; form (2) | identity | EXACT | lib; the 8 generators are independent |
| S2 | coefficients of (2) are gaps: q=g21,p=g23,j=g0,c=g6,e=g9,b=g20,a=g18; d=g4-g18; q+d,p+d,b+d,a+d are gaps g13,g17,g10,g4 (so d<=0 is a minor sum) | identity | EXACT | sympy on profiles |
| S3 | extra zeros of f18 are {0,6,9,18,20,21,23}; a new contact where f18>0 forces d=0 | identity | EXACT | zero set |
| S4 | each of the seven 7-sets contains a six-set in the orbit of an earlier case | citation | EXACT (resolution) / CITED (cases) | orbit_rep against Appendix E; stage from Table 6 |
| S5 | spinor chart: H1=0,H2=u(r-1),H3=-u,H7=1,H12=u(r-1)+1/s,H0=ruZ,H4=-uWY,H6=beta+r(u+1)Z+r/W; D=H4-H18; H21-H11+D identity; E=-sA2WKL; F1/E=uZW/K | identity | EXACT | 24 spinor products, bracket a_i+b_j |
| S6 | denominators of all 24 H_sigma lie in W*A2*L | identity | EXACT | factored denominators |
| S7 | curve (3): H11 L = N polynomial; quadratic in W; A(Z) as stated; C(Z) linear; basepoint (0,1) | identity | EXACT | sympy |
| S8 | Res_W(C,L), Res_W(C,A2), Res_W(C,K) = -Z Q2; Q2 at (2,3,4,5) | identity | EXACT | sympy resultants |
| S9 | finite poles real; D polynomial so d=-eta Im D -> 0 there | analytic | EXACT + PROOF-READ | all factors linear in Z and W over R; 20 random parameter points |
| S10 | W=inf: three real roots of A; (W-1)/A2 -> 1/Y; at Y=0, WY -> real limit | analytic | EXACT + PROOF-READ | limit (-hr+h-ru+r+u)/(ru-u-1) computed |
| S11 | Z=inf: (5) and no other branch; no (inf,inf) | analytic | EXACT | Newton polygon: all monomials have i-k<=1; Z^3W^2 coefficient -r^2 s(u+1) |
| S12 | K=0 points: Lemma 2 | analytic | PROOF-READ | 0 = grad C(v) + f H_kappa v; Birkhoff |
| S13 | spinors vanish only on W A2 K L = 0 | identity | EXACT | sympy solve |
| S14 | Lemma 3 (Rouche, genericity) and (2,3,4,5) generic | analytic | PROOF-READ + EXACT | simple roots of A and (5), gcd(A,C)=1, Q2 nonreal at (2,3,4,5) |
| S15 | harmonic completion | analytic | PROOF-READ | see report; boundary types: L_O=0, tie, collision, puncture |
| S16 | chart covers every eight-distinct source | citation | CITED + NUMERIC | 9 realizable points reproduced (error 1e-9) |

## D. Falsification

- NUMERIC: finite-pole and Z=inf projections real at 20 random rational parameter points
- NUMERIC: chart formulas reproduce 9 realizable T_45 points (max error 1.3e-9), r,s,u,h real
- NUMERIC: min Phi over SO(4) at those points ~ 2e-13

Exact replays:
- EXACT: T_45 inventory and identity (1)
- EXACT: gap expressions of all eight coefficients
- EXACT: seventh-tie resolution to T_7,T_36,T_42,T_40,T_46,T_29
- EXACT: 24 spinor products: H_0,H_2,H_3,H_4,H_6,H_7,H_12, D, H21-H11+D, E, F1/E
- EXACT: denominator divisors of all H_sigma
- EXACT: curve (3), A(Z), deg C(Z)=1, basepoint
- EXACT: three resultants and Q2 at (2,3,4,5)
- EXACT: equation (5) and the Newton polygon
- EXACT: spinor zeros lie on the divisors

## E. Verdict

**confirmed-modulo-citations**. Primary proof of paper T_45={1,2,3,7,11,12} (literal match): all chart identities (H_sigma, E, F1/E, D, curve (3), resultants, (5), Newton polygon) re-derived exactly from the stated spinors; every puncture is listed and d->0 or Lemma 2 applies at each; seventh ties reduce to T_7,T_36,T_42,T_40,T_46,T_29 (stages 1-2); no gap found.

All displayed algebra was re-derived from the spinors in the document. The puncture atlas is complete: finite poles (all real, D polynomial), W=inf (three real roots of A), Z=inf (two real V-roots; Newton polygon excludes other branches), K=0 (Lemma 2), basepoint (0,1) (D polynomial). Note: the bracket convention is a_i + b_j (B spinors carry b-values), not a_i - c_j.


Unchecked:
- the cited earlier six-set cases and Appendix B / one-circle theorem
- the earlier proof that the chart's denominators are nonzero for eight-distinct sources (supported numerically)

## Notes on the analytic steps (auditor)

- Curve inside X, contacts kept: the five values H1,H2,H3,H7,H12 are real constants of the chart and H11=h is the defining equation, so the six gaps -eta Im H_sigma vanish identically on the curve; lambda = i eta/E with eta fixed. One complex dimension remains, and H11 is nonconstant (quadratic W-coefficient A != 0).
- Lower bound on Omega: on Omega the cost is negative at O, so it is not a minor sum; by S2 this forces d > 0 with the seven minor coefficients equal to gaps (>= 0), hence L_O >= d Phi_{f18}(O). Where d -> 0, liminf L_O >= 0.
- Puncture list (complete in P^1 x P^1): finite points on W A2 L = 0 (projections are roots of linear factors with real coefficients, W solved linearly, so the points are real; D = W(Y-1/s)-beta is a polynomial with real coefficients, so d = -eta Im D -> 0, also at the basepoint (0,1) where L and K vanish together); W = inf at the three real roots of A ((W-1)/A2 -> 1/Y real gives q+d -> 0; at Y=0 the real-analytic branch gives WY -> (-hr+h-ru+r+u)/(ru-u-1), so a+d -> 0); Z = inf at the two real roots of (5) (D -> rV - beta); no point with Z and W both infinite (the Z^3 W^2 monomial dominates every other monomial). K = 0 points: H_sigma holomorphic and F1/E = uZW/K diverges, so Lemma 2 applies (derivation 0 = grad C(v) + f H_kappa v re-done; ker H_kappa = E is library-replayed; Birkhoff closes it). No reality is needed there.
- Several denominators vanishing at once is harmless here: the bound d -> 0 uses only that D is a polynomial with real values at real points; it is not a single-collision (Corollary B.3) argument. Corollary B.3 is not invoked.
- Boundary of Omega at finite non-puncture points: L_O = 0, a new tie (seventh contact; resolved in S4 to stage-1/2 cases), or a collision (spinors are nonzero off the divisors, so the configuration exists; covered by Appendix B / one-circle theorem, CITED). The minimum principle then gives a contradiction as written; the compact normalisation of the component reaches Z = inf because gcd(A,C)=1.
- Lemma 3: Rouche on a complex line through the counterexample gives nearby solutions of H11 = h' for nearby real (r,s,u,h'); strict gaps, distinctness and Phi_O < 0 persist by continuity. The generic conditions hold at (2,3,4,5) (checked: simple roots of A and of (5), gcd(A,C)=1, Q2 = -4(4Z^2-3Z+1) nonreal-rooted, so disjoint from the real pole projections), hence off a proper algebraic set.
- Convention note: the chart's bracket is a_i + b_j with B = (t,0,1,inf) the b-spectrum (paper convention); with a_i - c_j the stated formulas do not hold.
