# Audit: final_referee_T41_bad54_uniform_20261022

Evidence level: EXACT for the cubic restriction, the nodal identity (3) (1011 terms reproduced), the Heron identities (4), the discriminants, the DT identity, the six Lemma 2 certificate tables and the symmetry group element; EXACT (randomized identity test: exact rational arithmetic at 4 to 8 random rational points) for det H_F = -4(q+1)^2 G^2/Q^2 and the Lemma 3 principal-minor, resultant and congruence identities; PROOF-READ for the logic of Lemmas 2 and 3; CITED for the capacities (1) and for the "complex rank bound"; NUMERIC for falsification of (5) and for the rank at real roots.

PDF: proofs/manin_transfer_20261021/final_referee_T41_bad54_uniform_20261022.pdf (5 pages). Auditor label C17-T42. Replay: `replay/final_referee_T41_bad54_uniform_20261022__check.py` (27 checks, all pass, about 42 s).

## Deviations and limits

- The Lemma 3 identities were checked at random rational base points by univariate exact interpolation in c, not by symbolic expansion. That covers the height-free principal block (c+1)^2 P(c), Res_c(G,P) = W^2, the S± square classes and P S± ≡ (U1 ± W)^2 mod G. The term counts 2176, 437 and 339 were not reproduced.
- Table row 6 of (5): my reduced numerator has 2842 terms; the note says 4135. The other five rows match exactly (213, 473, 372, 957, 1161). The difference is probably the note's convention of clearing by "canonical degree" rather than reducing the fraction. All coefficients are nonnegative either way. I also checked directly that R - R00(1 - delta_a a)^2 (1 - delta_b b)^2 has a one-signed nonnegative numerator in every row.
- The capacities (1) (SOC, old-H and chain inequalities) were not audited. They are CITED.
- I did not test the main statement numerically. The main statement is a special case of the conjecture.

## A. Statement

**Theorem 1 (p.1).** Let x be in X, and suppose c(x), up to a positive scale, equals g = q f6 + p f2 + j f5 + c f53 + e f0 + d f54 + b f3 + a f15 with all coefficients >= 0. Then Phi_{c(x)}(O) >= 0 for all O in SO(4).

The generators all lie in F_{T_42}, where T_42 = {1,2,3,6,7,8}. f54 = h_{4231,-} is the only non-minor, and the eight are linearly independent (EXACT). Every spectrum and every coefficient zero is covered, because the argument uses only grad C(x) = 0. This note is a branch; the rest of N_{T_42} is handled by the primary T41 note, which is not audited here.

**Lemma 2.** If all applicable capacities in (1) fail and Fa, Fb >= 0, then G > 0 for all q, c >= 0.

**Lemma 3.** No spectrally realizable source has (Fa < 0, Fb < 0), (Fa = 0, Fb < 0) or (Fb = 0, Fa < 0).

## B. Dependencies

- Capacities (1): pj >= 1; a(p-1) >= 1; b(j-1) >= 1; pae >= a + e; jbe >= b + e. They come from the SOC, the old-H family [3,5,54] (a_H(b_H - d) >= d^2) and the chains [2;0,15;54] and [5;0,3;54]. **CITED.** The stabilizer that exchanges (p,j) and (a,b) exists: (3210, 3210, transpose), which fixes T_42 (EXACT).
- The "complex rank bound": rank H_F(w) <= 12 whenever F = det H_F = 0, attributed to HMSV and root_moc4_global_real_Hessian_inertia. **CITED (external).** The auditor checked that this citation can be avoided: the bordered determinants on I ∪ {0} and I ∪ {1} (I = {2..13}) have numerators divisible by G for this source as well (EXACT, randomized, 4 points). So the Schur argument of the negative-Heron note gives rank <= 12 wherever det H_I != 0, which is a dense set of each root sheet, and closure gives the rest. NUMERIC: the rank is exactly 12 at 5 real roots.
- The paper: X ⊂ Sing C, cubic (A.2) (lib, EXACT).

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | d = 0 gives a minor sum; d = 1; (1) gives the support inequality (p.1) | citation | CITED (capacities); EXACT (stabilizer) | |
| 2 | C = -(q+1)(ab+ac+a+bc+b+c+ce+1+e) < 0 and Q as displayed; h = -C/Q (p.1) | identity | EXACT | |
| 3 | (2): det H_F = -4(q+1)^2 G^2/Q^2 (p.1) | identity | EXACT (randomized, 8 pts) | |
| 4 | (3): G = (c+1)^2 {V L^2 - B0 L - pjeC0}, degree 8, 1011 terms; L >= T (p.2) | identity | EXACT | |
| 5 | (4): B0^2 + 4VpjeC0 = Fa Fb; Disc_q G = (c+1)^2 W^2 Fa Fb; Disc_c G = a^2 b^2 (q+1)^2 Fa Fb; Disc_a Fa = 16epj^2 V; Disc_b Fb = 16ejp^2 V (p.2) | identity | EXACT | |
| 6 | capacity failure implies delta_a a < 1 and delta_b b < 1, with kappa = min(1, 1/e) (p.2) | inequality | PROOF-READ | e <= 1 uses a(p-1) < 1; e > 1 uses pae < a + e |
| 7 | (5) R >= R00 (1 - delta_a a)^2 (1 - delta_b b)^2 > 0: six charts cover the domain up to the interchange (p,a) <-> (j,b) and boundaries; positive tables (p.2) | certificate + case split | EXACT (tables); PROOF-READ (coverage) | coverage re-derived for e < 1 and e > 1 with pj < 1 |
| 8 | boundaries by continuity (p.3) | analytic | PROOF-READ | |
| 9 | F = bFa + aFb + R > 0; DT = VT^2 + F + pjeC0; C0 >= -T(pa + jb); V - p^2 je > 0; Phi(T) > 0 and Phi'(T) > 0, so G > 0 (p.3) | identity + inequality | EXACT + PROOF-READ | |
| 10 | Lemma 3: a negative factor with a nonpositive partner forces p, j, e, a, b > 0 (p.3) | inequality | PROOF-READ (auditor: Fa at j = 0 equals a^2(p-e)^2 >= 0, so Fa < 0 forces j > 0) | |
| 11 | height-free principal determinant (c+1)^2 P(c), P quadratic in c; Res_c = W^2; S± = Fa b^2 (q+1)^2 L7^2 and Fb a^2 (q+1)^2 M7^2; P S± ≡ (U1 ± W)^2 mod G (p.3) | identity | EXACT (randomized, 4 pts) | |
| 12 | P < 0 off W L7 = 0 (U1 = alpha P on G = 0); rank exactly 12 there by the complex rank bound (p.4) | inequality + citation | EXACT (U1 = alpha P - a2 G) + CITED (rank bound; auditor alternative above) | |
| 13 | equal weak sign of maximal principal minors; continuity; second derivative -8(q+1)^2 G_c^2/Q^2 < 0, so rank exactly 12 and some 12-minor is strictly negative (p.4) | analytic | PROOF-READ | |
| 14 | one-zero boundary Fa = 0, Fb < 0: alpha_q = V W^2 > 0; complete square, so G_q = 0; 4 alpha_q G_a = -(c+1)^2 W^2 (Fa)_a Fb != 0; a(q) with a' = 0; d²/dq² Fa(a(q)) = 8V^2 W^2/((c+1)^2 Fb) < 0; rank 12 by the a-curve (p.4) | identity + analytic | EXACT (alpha_q, Disc_q) + PROOF-READ (auditor re-derived the second derivative from a'' = -G_qq/G_a) | |
| 15 | final contradiction: delta in ker H(w), Q_v(delta) = 0, Q_v definite on K (p.4) | analytic | PROOF-READ | same argument as the negative-Heron note |
| 16 | the sign cases are exhaustive (p.5) | case split | PROOF-READ | opposite strict signs make Disc_q < 0 with alpha_q > 0 |

## D. Falsification

- (5) under failure of all capacities: 1,231,647 samples (seed 54, log-uniform in [1e-4, 1e4], with zeros). The minimum of (R - lower)/scale was -3.5e-15, which is roundoff, since the bound is attained in the limit. No violation.
- Rank: at 5 real c-roots of G (mpmath with 60 digits), the smallest singular values were about 1e-51, 1e-50 and 0.2-4.3. So the rank is exactly 12.

## E. Verdict

**confirmed-modulo-citations.** Citations: the capacities (1), and the complex rank bound (HMSV / root_moc4_global_real_Hessian_inertia), which the auditor's bordered-minor check could replace. No GAP was found. A presentation discrepancy: table row 6 has 2842 reduced numerator terms, not the 4135 stated.
