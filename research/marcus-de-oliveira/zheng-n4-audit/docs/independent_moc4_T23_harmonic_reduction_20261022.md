Evidence level: EXACT for every algebraic identity of the document (decomposition, chart, poles, residues, exhaustiveness of the chart); PROOF-READ for the harmonic minimum argument; CITED for the boundary theorems it invokes; NUMERIC only for the cited SOC inequality.

# Audit of `independent_moc4_T23_harmonic_reduction_20261022` (3 pages)

Document: proofs/manin_transfer_20261021/independent_moc4_T23_harmonic_reduction_20261022.pdf. It is a branch (reduction) document used for paper case T_24. Auditor batch C10-T24-T26. Replay script: `replay/independent_moc4_T23_harmonic_reduction_20261022__check.py`. Its output is in `..._check_output.txt`, with 41 EXACT checks, all passing.

## Deviations and limits

- The supplement writes gaps as `g_j = -mu Im H_j`, with `F_sigma = prod (a_i - c_sigma(i))` and the point `c_4 = infinity` normalized away. I used the same convention, dropping the factor at infinity. This changes every F by one common nonzero constant, which the normalization `H_3 = 1` absorbs.
- Not checked by me: the "norm-kernel theorem" for T ∪ {13}, the "old SOC" inequality (tested numerically only), and the earlier paper cases used for the extra contacts 1, 10, 12 and 22. These are all CITED.
- I did not run a numerical search for violations of the main statement, because the statement is conditional and is itself a special case of the conjecture.

## A. Statements

Theorem 1 (p. 1). Let T = {2,3,4,6,7,8}. This is the paper's T_24 up to the group G (checked with `find_g`). Assume (*) holds for the seven-contact families T ∪ {16} and T ∪ {17}. Then (*) holds for T: for every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for all pi in T, Phi_{c(x)}(O) >= 0 for all O in SO(4). Contacts are non-strict and spectra are arbitrary. Degenerate configurations (coincident points, a_i + b_j = 0, a real spectrum) are sent to the cited Appendix B and one-circle theorems. The case T ∪ {13} is sent to a separate "norm-kernel theorem".

Contact-set bookkeeping (EXACT, lib):
- generators(T) = {0,1,2,3,4,5,6,53,54,67}. That is 10 generators of rank 8, with f54 the only nonminor one. The document's 8-term decomposition omits f1 and f4. This is legitimate: f54 + f1 = f0 + f6 + f53 and f54 + f4 = f2 + f3 + f6 + f53 (both replayed). A cone vector whose 8-basis coefficient delta is at most 0 therefore has a nonnegative minor decomposition.
- T ∪ {16} has no Table-5 label, but it contains G-images of T_25 (drop 6) and of T_26 (drop 7). T ∪ {17} is Table 5 #5 and also contains a T_26-image (drop 6). The two hypotheses are therefore discharged by paper cases T_25 and T_26, which are the other two documents of this batch.
- T ∪ {13} is Table 5 #6. It also contains a T_28-image (drop 6).

## B. Dependencies

| dependency | used for | status |
|---|---|---|
| Theorem 1.2 (two-orthogonal) | O(4) or SO(4) is enough | CITED (paper) |
| Lemma 4.2 / generators of N_T | the decomposition has nonnegative coefficients | EXACT for the generator list; Lemma 4.2 CITED |
| "old SOC": delta^2 <= g10 g13 implies safety for g10 f5 + g13 f2 + delta f54 (+ minors) | bounding g13 near a nonreal pole | CITED; NUMERIC support (6 samples at the boundary, worst min Phi 3.9e-14; a control at 3x the bound goes negative, -4.0e-2) |
| T ∪ {16}, T ∪ {17} | hypotheses | discharged by T_25 / T_26 via the group (EXACT inclusion) |
| T ∪ {13}: "norm-kernel theorem" | boundary contact 13 | CITED (Table 5 #6; also contains a T_28-image) |
| contacts j in {1,10,12,22} | "individually certified registry cases" | resolved by me: j=1 ⊇ T_3-image (drop 6); j=10 ⊇ T_1 (drop 3); j=12 ⊇ T_3 (drop 3); j=22 ⊇ T_10 (drop 6). CITED |
| contacts j with f54(j) > 0 (0,5,9,11,14,15,18,19,20,21,23) | force delta = 0, so a minor decomposition | EXACT (f54 values) |
| Appendix B (Theorem B.1, Lemma B.2) | collision or shared-zero boundary points, and the chart's excluded configurations | CITED (project-audited). Corollary B.3 is NOT used |
| Theorem C.1 (one spectrum on a circle or line) | real x, which makes the A spectrum real | CITED (paper) |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | 8-term decomposition of g with delta = g5 - g16 = g9 - g17 | identity | EXACT | reproduces all 10 generators; the span has rank 8 |
| 2 | delta <= 0 implies a minor decomposition; delta > 0 implies (1): L_O >= -M_O delta with 0 < delta <= min(g5,g9) | identity / inequality | EXACT + PROOF-READ (argument supplied by auditor via f54+f1, f54+f4) | |
| 3 | extra contacts reduce to listed cases | citation | EXACT inclusions; cases CITED | lib group search |
| 4 | the SOC gives the strict reverse inequality delta^2 > g10 g13 on the negative region | citation | CITED (NUMERIC support) | |
| 5 | chart A=(z,0,1,x), C=(c,v,w,inf) with H_2..H_8 = (0,1,s,r-u,r+t-u,r), H_9 = r+t-uX | identity | EXACT | sympy, rational in (x,s,t,r,u) |
| 6 | s,t,r,u real and nonzero | claim ("They are nonzero") | PROOF-READ, argument supplied by auditor | s = x(v-1)/(w(x-1)), r = (c-v)(z-1)/((c-z)(x-1)), u = (w-c)(v-z)/(w(c-z)(x-1)), t = c(v-z)/(w(c-z)). Each is nonzero for eight distinct points with no a_i = c_j |
| 7 | the chart is exhaustive | claim | EXACT (auditor) | the ratios of any configuration give back z, w, c, v identically. The only chart singularity, l = n = 0, forces real x (l = tx + sX); then z is real, so the A spectrum is real, which is the one-circle case |
| 8 | J relation for z - c; z - v = th(z-c); H1 formula; g10 = -(tr/q) g1; H5 -> 1 at infinity | identity | EXACT | H10 + (tr/q)H1 is the real constant (r^2+rst+rt-ru-su)/(r-u) |
| 9 | every finite nonreal pole of every H_j (j = 0..23) is a zero of J | identity | EXACT | denominator factors are only x, q(x-1)-t and J |
| 10 | at a nonreal zero of J (t ≠ 1): a simple pole, z = c = v nonreal, res H13/res H16 = z/(z-1) is not real, both residues are nonzero | identity + argument | EXACT (J divides the denominators of H13 and H16 exactly once; the ratio identity holds mod J; raw numerators (1-z)(x-z)(z-w) and z(w-z)(x-z)) + PROOF-READ (z is a real Möbius image of X with determinant -(1-t)r; z = x iff X = 0 or d = r) | |
| 11 | no negative-response sequence tends to a nonreal pole | analytic | PROOF-READ | H5, H9, H1, H10 have no nonreal pole, so g5, g9 and delta stay bounded. g1 > 0 at the limit (A1 ≠ 0, q ≠ 0), so g10 > 0. Hence g16 and g13 stay bounded. Two simple poles with an R-independent residue pair cannot both have bounded imaginary part |
| 12 | qA1 = 0 gives g1 ≡ 0, the extra contact 1 | identity | EXACT (from the H1 formula) | covered by the T_3 reduction |
| 13 | harmonic minimum on the x-sphere: real punctures (g9 = -mu u Im x -> 0, so liminf L_O >= 0 by (1)); infinity (g5 -> 0); regular boundary points (L_O = 0, a new tie, or a collision); nonreal punctures excluded | analytic | PROOF-READ | the puncture list is complete by step 9. L_O = -mu Im R(x) with R a real rational function, so L_O is harmonic off the poles |
| 14 | t = 1: a strict negative point persists under a small real change of t | perturbation | PROOF-READ | the chart is continuous at a non-pole point, the conditions are open, and the ratios stay real. Only strict six-contact points are perturbed |

## D. Falsification

The document makes no new inequality claim. The cited SOC (delta^2 <= g10 g13) was tested at the boundary on 6 random (g10, g13) pairs, seed 11, with 10 BFGS restarts each: worst min Phi = 3.9e-14. At delta = 3 sqrt(g10 g13) the minimum is -4.0e-2, which shows the test can detect a violation. All exact identities passed. No counterexample was found.

## E. Verdict

**confirmed-modulo-citations.** Every step is EXACT or PROOF-READ, apart from the cited SOC inequality, the norm-kernel theorem for T ∪ {13}, the paper cases T_1, T_3, T_10 (extra contacts) and T_25, T_26 (the hypotheses), Appendix B and Theorem C.1. The three Corollary B.3 conditions do not apply, because B.3 is not invoked. The nonreal puncture is handled by a direct two-residue argument instead.

Presentation-level gaps (the argument is recoverable; I supplied it):
1. p. 1: "They are nonzero" is not justified. Repair: the factorizations in step 6.
2. p. 2-3: "exhaustive rational chart" is asserted without proof. Repair: step 7 (the inverse map, plus the real-x exception).
3. p. 1: "the canonical decomposition immediately forces delta <= 0, hence nonnegativity" needs the relations f54+f1 and f54+f4 = nonnegative minor combinations, which are not stated. Repair: step 2.
