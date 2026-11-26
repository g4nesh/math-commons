Evidence level: EXACT for the chart's contact set, the four inverse ratios, (2)-(5), the forced-sign identities, (7) and the merged-stratum identities, (8), (9), the reduced denominator of L_V, (11), the SOS identity and the auxiliary-segment endpoint identities (43 exact checks in `replay/final_referee_T10_A_branch_20261022__check.py`). PROOF-READ for the sign, interval and convexity arguments. CITED for the derivation of the chart (B-branch note), T_7, Appendix B, Theorem C.1 and SOC. NUMERIC for the final bound.

# final_referee_T10_A_branch_20261022: branch A < 0 of {0,1,2,4,7,9} (paper T_11)

PDF: `proofs/manin_transfer_20261021/final_referee_T10_A_branch_20261022.pdf`, 4 pages. This is a branch document.

## Deviations and limits

- This note takes its chart from `root_moc4_T10_B_branch_20261022` (not in this batch). That note derives the chart from real ratios of contact differences, and I did not audit the derivation. What I checked exactly:
  - the chart reproduces exactly the contact set T for generic parameters;
  - the four inverse ratios h = z(v-1)/(v(z-1)), t = (z-c)(v-r)/((r-c)(v-z)), s = (w-v)(1-c)(z-r)/((z-c)(v-r)) and l = (w-r)z/((w-v)(z-r)) are identities.

  So the chart is the inverse of those ratios. That the ratios are real on every distinct supporting configuration, and that the chart's denominators are nonzero there, is CITED.
- The sign convention was fixed by matching (4): A_0 = ∞ contributes the factor +1 and kappa = i mu/((z-c)(v-r)). With rho = -1 every identity (2)-(5) holds exactly.
- The SOC on (f2, f12, f26) is CITED, source unresolved. Its orbit equivalence with the T9 SOC is EXACT; see the T10 main check.

## A. Statements

**Theorem 1 (p. 1).** Let T = {0,1,2,4,7,9} (G-equivalent to the paper's T_11). Every realizable supported gap vector with A = g16 - g3 < 0 is safe. On the eight-distinct, nonreal-v chart with a = g15 > 0, b = g8 > 0 and nonzero inverse denominators, g15 g8 >= (g3 - g16)^2. Additional ties are allowed, and the merged coordinate C = 0 is included. Collisions, the circle stratum and b = 0 go to separate results:
- b = 0 adds contact 8, and T ∪ {8} contains the paper's T_7 (EXACT);
- real v makes the A spectrum real (z = hv/(1+Hv) is real), so Theorem C.1 applies;
- collisions are covered by Appendix B.

## B. Dependencies

| dependency | status |
|---|---|
| chart derivation and nonvanishing of its denominators (B-branch note) | CITED. The inverse identities and the contact set are EXACT |
| paper case T_7 (b = 0 boundary) | CITED. The containment is EXACT |
| Appendix B (collisions), Theorem C.1 (real v) | CITED |
| SOC on (f2, f12, f26) | CITED (source unresolved). NUMERIC: no violation |
| decomposition (1) | EXACT (T10 main check) |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | contacts of the chart are exactly {0,1,2,4,7,9}; the inverse ratios h, t, s, l hold | EXACT |
| 2 | (2)-(5) | EXACT |
| 3 | forced signs: t < 1; t = (b-f)/M; g6 = b/beta gives 0 < theta < 1, theta >= e/(b+e) and beta e <= b theta; h(1-f) + e + f = 0; f != 1; D has the sign of h and H; l > 0; then s > 0, rho < 0, N0 < 0 | EXACT (identities) / PROOF-READ (signs) |
| 4 | b = 0 reduces to T_7 | EXACT containment, CITED case |
| 5 | merged C = 0: theta = 1 - k; (7); for t < 0 the third identity contradicts e > 0 and V >= 0 | EXACT / PROOF-READ |
| 6 | merged C = 0 with t > 0: p > 1, 0 < b < 1, s formula; G, 2 theta Re w and theta q formulas; V >= 0 rewritten as \|w - (p-1)\|^2 >= pb(p-1)/theta (I checked it equals p b theta \|v-1\|^2 >= p - 1 via (7)); then e^2 + e <= p(1-b)(e+b) <= e + b | EXACT / PROOF-READ |
| 7 | generic C != 0: (8), R = \|xi\|^2 = t beta W/(theta b), (9) | EXACT |
| 8 | (10) follows algebraically from (9); at k = 0 it holds as a rational-function identity because L_V reduces to denominator 2 M b theta (1+e)(f-1) | EXACT (reduction) / PROOF-READ |
| 9 | t > 0: W, Delta > 0; gamma has the sign of D (with the f < 1 formula for gamma, EXACT); (11) and L - (e^2 - b) = {theta(b+e) - e}M + eq >= 0; alpha != 0 because C != 0; so e^2 <= b | EXACT / PROOF-READ |
| 10 | t < 0: f < 1 is impossible (W > 0); for f > 1, delta < H via W(H) = theta h b/f and dW/d delta > 0 for delta >= H | EXACT / PROOF-READ |
| 11 | auxiliary segment to delta0 = (1-b)H: E(delta0) = (f-b)H, W(delta0) = -beta b H; WE is convex in delta, so W < 0 and R > 0 on the segment; SOS L_V^2 - R = alpha^2 P/(4 b^2 theta^2 (1+e)^2) with P > 0; alpha(delta0) = 0, R(delta0) = beta^2, L_V(delta0) = beta; so L_V > 0 on the segment and L_V > sqrt(R) at the actual point, contradicting (10) | EXACT / PROOF-READ (the segment varies only the auxiliary algebraic variable, and L_V is continuous because M = E/theta > 0) |
| 12 | return: e^2 <= b becomes g15 g8 >= (g3 - g16)^2, then SOC | PROOF-READ + CITED (SOC) |

## D. Falsification

2 x 400000 random chart points (seed 2027, both signs of mu). Only 6 were supporting with A < 0 and b > 0, all with t > 0. The minimum of g15 g8 - e^2 was 1.65, so no violation. The region is thin, so this test is weak.

## E. Verdict

**confirmed-modulo-citations.** I found no gap in the document's own argument. The things it relies on and does not prove are the derivation and exhaustiveness of the chart (B-branch note), T_7, Appendix B, Theorem C.1 and the SOC inequality, whose source is unresolved.
