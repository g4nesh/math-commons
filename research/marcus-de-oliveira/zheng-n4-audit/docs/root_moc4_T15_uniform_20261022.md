# Audit: root_moc4_T15_uniform_20261022 ("A six-contact support inequality from a quadratic spectral curve")

Evidence level: EXACT for the contact-set identification, source cone, decomposition (1), coefficient signs, the seven seventh-contact witnesses, all chart identities (3), the curve (4), (5), the pole locus, the full-height formula and both branches over v = 0; PROOF-READ for the sign logic and the harmonic minimum principle; CITED for the paper cases T_2, T_3, T_8, T_15, Appendix B and Theorem 1.2; NUMERIC for one brief sampling test.

PDF: `proofs/manin_transfer_20261021/root_moc4_T15_uniform_20261022.pdf` (4 pages). Primary proof of paper case **T_16 = {0, 2, 5, 6, 9, 12}** (supplement historical label T15). The note works with T = {0, 2, 5, 6, 8, 18} and reference pi_0.
Replay: `replay/root_moc4_T15_uniform_20261022__check.py`. Output: `..._check_output.txt` (35 EXACT checks pass, 1 NUMERIC pass, exit 0).

## Deviations and limits

- The sign convention for the factor at infinity was fixed from the note's own display F0/D = -x(v-1): F_pi = -prod_{i>=1}(A_i - C_pi(i)). With the opposite sign the replay of R2, R6, R8 fails, as it should. With this sign every display matches.
- The note cites the six-point circle theorem for v real. I supply an argument that removes this dependency (see C, step 9). The supplement checker was not used.

## A. Statements

**Theorem 1 (p. 1).** Let x in X, with g = Re(F - F_{pi_0}) >= 0 and g_j = 0 for j in T = {0,2,5,6,8,18}, further contacts allowed. Then Phi_{c(x)}(O) >= 0 for all O in SO(4), for all spectra. T is exactly {sigma^{-1} : sigma in T_16}, the transpose, which is the group element (rho, rho', eps) = (2103, 2103, 0) of the library convention (EXACT). Since G preserves X and Phi, this is (*) for T_16: non-strict contacts, every O, collisions included.

Intermediate claims:
- (1)-(2): L_O(g) >= -C_O d with 0 < d <= g3 - g7 <= g3, where d = g12 - g16.
- d = -mu s r Im v/|v|^2 on the eight-distinct stratum.
- The source curve is the quadratic (4) in E over the v-line. Every R_j is holomorphic over the upper half-plane, and g3 -> 0 at both points over v = 0.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Lemma 4.2 / App. D | core-doc | source cone {1,2,4,5,6,9,39,45,46,47} | generators EXACT |
| paper T_2, T_8, T_3, T_15, T_2, T_3, T_8 for j = 4, 7, 10, 16, 19, 22, 23 | paper-case | seventh contacts not killed by f46 | CITED. The witness T' contained in T u {j} is EXACT for each j. T_2 and T_3 are at stage 1, T_15 at stage 1 and T_8 at stage 4, all before T_16 at stage 5. The note's labels T1, T7, T3, T14 are supplement historical labels for paper T_2, T_8, T_3, T_15 |
| Thm B.1, Lemma B.2 | common-lemma | collisions | CITED (audited) |
| paper Sect. 3.1 / proof of Prop. 4.3: every nonzero x in X has finite spectra and alpha != 0 | core-doc | full-height continuation at finite boundary points | CITED (paper) |
| six-point circle theorem | common-lemma | v real | cited by the note, but **not needed**: see step 9 |
| Theorem 1.2 | core-doc | unitary form | CITED (audited) |
| Corollary B.3 | - | not used | - |

## C. Steps

| # | claim | kind | status |
|---|---|---|---|
| 1 | generators; only f46 non-minor; f46 + f1 = f4 + f6 + f45 and f46 + f9 = f4 + f5 + f47; d equals [f46] - [f1] - [f9] on the generators | finite | EXACT |
| 2 | d <= 0 gives a nonnegative minor decomposition | argument | PROOF-READ |
| 3 | (1) canonical decomposition, in 24 coordinates | identity | EXACT |
| 4 | the seven coefficients other than d are >= 0 on every generator, hence on N_T; so (2) holds when d > 0 | finite | EXACT. The note's "all eight coefficients are nonnegative on the full source cone" is true only for d > 0 (presentation) |
| 5 | an extra contact j with f46(j) > 0 forces d <= 0; the zeros of f46 off T are exactly {4,7,10,16,19,22,23} | finite | EXACT |
| 6 | witnesses for the seven remaining j | citation | EXACT witnesses, CITED cases |
| 7 | chart A = (inf, v, xv, z), C = (r, 1, 0, w): R2 = 1 - x, R6 = x(r-1), R8 = r - x; hence x real, kappa = i mu/D, r real; R5 = -sx, so s real | identity + argument | EXACT + PROOF-READ |
| 8 | (3): R1 - R3 = E, R9 - R7 = -rE; (R12 - R16) + sr/v = rx - r + sx - x + 1; R16 - (1-s)R12 = const; the last combination equals x(x-1)(r+s-1) | identity | EXACT |
| 9 | v real: the note cites the six-point circle theorem. **Auditor's argument:** the identity of step 8 holds on the whole chart, so for real v, R12 - R16 is real and d = g12 - g16 = 0, and step 2 gives L_O >= 0 | argument | PROOF-READ (argument supplied by auditor; the citation becomes unnecessary) |
| 10 | gap relations g9 - g7 = -rb, a = -x(d+c), d = s g12, g16 = (1-s) g12; signs b > 0, r < 0, x < 0, 0 < s <= 1; s = 1 is the extra contact 16 | argument | PROOF-READ |
| 11 | normalization mu = 1 by conjugation and rescaling; then Im v > 0 | argument | PROOF-READ |
| 12 | w = v(E - sx)/(vE - s), z = (w(s+v) - v)/(s - 1 + w) invert (s, E); (4) Q = vD(R18 - u) = v(v-r)E^2 + B(v)E + C(v) with the printed B, C | identity | EXACT |
| 13 | v is nonconstant on every component (the leading coefficient v(v-r) is nonzero over nonreal v) | argument | PROOF-READ |
| 14 | the denominators of all 24 R_j are products of v, vx - 1, D; (5) Q(v, E_D) = -s(r-v)(vx-1)^2(s+v-1)/(v-1)^2, so D vanishes on the curve only over real v; E = inf only over v in {0, r, inf} | identity | EXACT |
| 15 | full height: F0/D = -x(v-1) != 0 over the upper half-plane, so every finite boundary point there is a nonzero point of X, hence realizable | identity + citation | EXACT + CITED (Sect. 3.1) |
| 16 | escape bound: d = -sr Im v/\|v\|^2 -> 0 as v tends to a nonzero real or infinity, on every branch | argument | PROOF-READ |
| 17 | over v = 0: Q(0, E) = rs(E - sx + 1); the finite point (0, sx - 1) has D = 1 - s != 0 and R3 = 1 - sx - x; the point E = inf (zeta = 1/E) is smooth with zeta = v/s + O(v^2), and R3 is regular there with value 1 - s - x; so g3 -> 0 and, by (2), d -> 0 | identity + analytic | EXACT + PROOF-READ |
| 18 | harmonic minimum principle on Omega (v in the upper half-plane, nontied gaps > 0, d > 0); boundary points are d = 0, an extra contact, or an escape; a negative constant is impossible | analytic | PROOF-READ |

Puncture inventory (point iii). The punctures of the compact normalization are the points over real v (including the poles over v = 0, 1/x, r, 1 - s, and E = inf) and the points over v = inf. At all of these except v = 0, the barrier d -> 0 is immediate from the exact formula. Over v = 0 the two points are listed exactly and handled by g3 -> 0 with d <= g3. No point of the upper half-plane is a pole.

## D. Falsification

NUMERIC (seed 31): 41 realizable sources were sampled on the curve (4), with real parameters x, r < 0, 0 < s < 1, u, v in the upper half-plane, and both roots E. The minimum of Phi over SO(4) was found by 12 BFGS restarts each. The worst value is 1.2e-18, so no violation. Every intermediate inequality ((2), the escape bound, the v = 0 barrier) follows from exact identities and needs no separate test.

## E. Verdict

**confirmed-modulo-citations.** The citations are paper cases T_2, T_3, T_8 and T_15 (stage order consistent), Theorem B.1, Lemma B.2, Theorem 1.2, and the surjectivity statement of Sect. 3.1. The cited six-point circle theorem is not needed; step 9 gives the replacement argument.

Presentation remarks:
- p. 1: the text has the typo "responsesty".
- p. 1: the wording of the coefficient-sign claim, see step 4.
- p. 1: the phrase "the transpose of the inventory family T15" is correct, as replayed.
