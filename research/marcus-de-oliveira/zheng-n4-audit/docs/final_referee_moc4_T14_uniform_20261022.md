# Audit: final_referee_moc4_T14_uniform_20261022 ("A six-contact inequality from a reciprocal half-plane obstruction")

Evidence level: EXACT for every linear-face identity, spectral-chart identity, the reciprocal decomposition (9) and the positive-weight identity; PROOF-READ for the sign and case logic; CITED for two auxiliary inputs that the paper does not list (the SOC (1) and the six-point circle theorem), each with NUMERIC support only; CITED for Appendix B and Theorem 1.2.

PDF: `proofs/manin_transfer_20261021/final_referee_moc4_T14_uniform_20261022.pdf` (3 pages). Primary proof of paper case **T_15 = {1, 2, 3, 4, 6, 7}** (supplement historical label T14, reference permutation pi_1).
Replay: `replay/final_referee_moc4_T14_uniform_20261022__check.py`. Output: `..._check_output.txt` (25 EXACT checks pass, 4 NUMERIC pass, exit 0).

## Deviations and limits

- The SOC (1) and the "three-A/three-C common-circle theorem" are used as established inputs. I did not audit their sources, and I could not find either among the inputs the paper lists (Appendix F: B.1, B.2, B.3, C.5, C.6, plus Theorem C.1 for a whole spectrum on a circle). The supplement's `docs/verification-scope.md` names the six-point circle theorem as a separate input, checked by its own script. For the SOC no source document is named; the note calls it the "old SOC". Both are tested numerically below.
- The supplement's companion checker was not read or run.

## A. Statements

**Theorem 1 (p. 1).** Let x in X, with g_sigma = Re(F_sigma - F_{pi_1}) >= 0 and g_i = 0 for i in T_15 = {1,2,3,4,6,7}, further contacts allowed. Then Phi_{c(x)}(O) >= 0 for all O in SO(4). All spectra are included: repeated points and A_i = C_j go to Appendix B, the case "v real" goes to the six-point circle theorem, and every other eight-distinct configuration is handled directly. The contact set is literally the paper's T_15, the contacts are non-strict, and every O is covered, so the statement is (*) for T_15.

Main intermediate claim: on every eight-distinct realizable source with Im v != 0 and d = g0 - g16 > 0, we have T >= X, and so ab >= (T/X) d^2 >= d^2 with a = g13, b = g10.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Lemma 4.2 / Appendix D | core-doc | the source cone {0..8, 53, 54} | generators EXACT |
| Thm B.1, Lemma B.2 | common-lemma | collisions | CITED (audited by the project) |
| "old safe SOC": a f2 + b f5 + d f54 is safe if a, b, d >= 0 and ab >= d^2 | common-lemma (unresolved source) | final step | CITED; not proved in this note; no document named; NUMERIC support: min over SO(4) of Phi_f54 + 2\|det O[13,34] det O[12,24]\| = -3.7e-30 (80 Nelder-Mead starts). Together with Phi_f2 = det O[13,34]^2 and Phi_f5 = det O[12,24]^2 this is exactly (1) |
| six-point circle theorem (A2, A3, A4 and C2, C3, C4 on R u {inf}) | common-lemma (supplement only) | the stratum v real | CITED; not in the paper's list; NUMERIC: 40 sampled realizable sources on that stratum, worst min Phi = 1.8e-16 |
| Theorem 1.2 | core-doc | unitary form | CITED (audited) |
| other paper cases | - | none used | - |

Stage order: T_15 is at stage 1 of Table 6 and uses no other contact case, which is consistent. It does use the two auxiliary inputs above, which Appendix G does not list. This is a traceability issue for the paper's dependency structure, not a flaw in this note's argument.

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | generators {0,...,8,53,54}; only f54 non-minor | finite | EXACT |
| 2 | g16 - g0 is 1 on f1, f4, f7, -1 on f54, 0 on the rest; the three replacement identities f54 + f1 = ..., f54 + f4 = ..., f54 + f7 = ... | finite | EXACT |
| 3 | d <= 0 implies a minor decomposition. For any nonnegative representation, lambda_54 <= lambda_1 + lambda_4 + lambda_7, and the identities trade it away | argument | PROOF-READ |
| 4 | (2) canonical decomposition for d > 0, with all coefficients gaps or d | identity | EXACT |
| 5 | (3) g11 = b + r + d, g19 = a + s + d | identity | EXACT |
| 6 | chart A = (z,1,0,x), C = (w,v,inf,tv) is exhaustive for eight distinct points (Moebius sends A2, A3, C3 to 1, 0, inf) | argument | PROOF-READ |
| 7 | (4) R2 = x - t, R3 = t(x-1), R4 = 1 - t; (P3 - P1)/(P2 - P4) = t, so t is real, kappa = i mu/(v(z-w)) and x is real | identity + argument | EXACT + PROOF-READ |
| 8 | J = x - t - R6 = (x - tv)[z(v+1-w) - v]/(v(z-w)) and R7 + t = tv/(x - tv) J; with Im v != 0 this forces J = 0, i.e. (5) | identity + argument | EXACT + PROOF-READ |
| 9 | (6) g11 = (1-x) b, g19 = ((t-1)/t) a | identity | EXACT (R11 - (1-x)R10 = x(1-t), R19 - ((t-1)/t)R13 = x - 1) |
| 10 | signs: -xb = r + d > 0 and -a/t = s + d > 0, so a, b > 0, x = -X, t = -T < 0, and (7) d <= Xb, d <= a/T | argument | PROOF-READ |
| 11 | mu > 0 without loss (conjugation) | argument | PROOF-READ (Phi is invariant under conjugating both spectra and kappa, for real O) |
| 12 | p, q, K = R10, z = p/(p-1), z - w, R13 = -t(L+1), U = -p(T+q), V = p(Xq+1)/(pq-p+1) | identity | EXACT |
| 13 | Im U = b/mu > 0, Im V = a/(T mu) > 0; Im D1 = CB > 0, Im D2 = -HA < 0, so D1, D2 != 0 | identity | EXACT |
| 14 | (8) R12 = (T-X) - Z/D1, R8 = (T-X) + Z/D2 | identity | EXACT |
| 15 | (9) decomposition of Z, and the positive-weight identity T C^2\|D1\|^2 m + X H^2\|D2\|^2 n = (T-X)(T C A\|D1\|^2 + X H B\|D2\|^2) | identity | EXACT (formal complex U, V) |
| 16 | m, n >= 0 and the bracket > 0 give T >= X; then ab >= (T/X) d^2 >= d^2 | inequality | PROOF-READ |
| 17 | SOC (1) and (2) give the support inequality | citation | CITED (see B) |
| 18 | v real: the six-point circle theorem | citation | CITED (see B) |

The case split is exhaustive: d <= 0 (step 3); d > 0 with a collision (App. B); eight distinct with v real (step 18); eight distinct with v nonreal (steps 8-17). There is no curve argument, so there are no punctures to check. Corollary B.3 is not used.

## D. Falsification

- NUMERIC (seed 7): 40 realizable sources from the chart (x, t < 0, z from (5)). Worst min Phi = 1.6e-14. On the 3 sampled sources with d > 0, the smallest (ab - d^2)/d^2 was 268 and the smallest T - X was 1.56.
- NUMERIC: the SOC bound above (worst -3.7e-30). Neither signed version Phi_f54 >= -+2 m1 m2 holds alone (both reach -0.5), so the absolute value is essential.
- NUMERIC: the v-real stratum (40 sources, worst 1.8e-16). No sampled source on that stratum had d > 0.

## E. Verdict

**confirmed-modulo-citations.** The note's own argument (steps 1-16) is complete and exactly replayed. Its conclusion rests on two cited inputs, the SOC (1) and the six-point circle theorem, plus Appendix B and Theorem 1.2. The project has not audited the first two, and the paper's Appendix F/G does not list them. Both are numerically supported. Until their sources are located and audited, the proof of T_15 is complete only modulo these two lemmas.

Presentation remark: on p. 1, "Thus the safe inequality family used below is already known in this coefficient inequality" gives no reference for the SOC.
