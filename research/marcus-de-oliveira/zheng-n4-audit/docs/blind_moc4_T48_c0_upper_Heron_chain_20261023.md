Evidence level: EXACT for G, the discriminant factorization, the chart identity, the 46874-term certificate and its four walls; PROOF-READ for the root-existence step; CITED: only the chain atom, proved in foundations; NUMERIC for falsification only.

# blind_moc4_T48_c0_upper_Heron_chain_20261023 (2 pp.): chain capacity on the upper Heron branch, c = 0 face

Auditor label C26-T49. Script: `replay/blind_moc4_T48_c0_upper_Heron_chain_20261023__check.py` (0 EXACT failures).

## Deviations and limits

- G is rebuilt from the definitions; my G = -(supplement's G). Theorem 1 is a purely algebraic implication: it takes G = 0 and the Heron signs as hypotheses. Whether realizable c = 0 sources satisfy F6 >= 0 and L >= 0 is NOT claimed here; the document says so itself.

## A. Statements

- Display (not used in the proof). With B0 = (b+d)(q+d), t = B0 + dj, L = at + (d-j)B0 and F6 = L^2 - 4d^2 j^2 B0, one has Disc_p(G|c=0) = (a+d+e)^2 F6 F8 for a homogeneous octic F8. EXACT (F8 has 321 terms).
- Theorem 1. Suppose all eight coefficients are >= 0, c = 0, d, e > 0, j >= d, G = 0, F6 >= 0 and L >= 0. Then a > 0 and p >= d^2(1/a + 1/e). Hence (a, e) chain capacity, so Phi >= 0 pointwise. The source lies on the face of contact set T_49 u {17} = Table 5 set #23.

## B. Dependencies

Chain (f5; f0, f2; f40), in the proved orbit (EXACT). Nothing else: G = 0 is a hypothesis.

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | d = 1; F6 >= 0 and L >= 0 give L > 0 (L = 0 would make F6 = -4j^2B0 < 0), so X^2 - (L/j)X + B0 = 0 has a positive root X (both roots positive, double root if F6 = 0) | PROOF-READ | |
| 2 | a = ((j-1)B0 X + j(B0 + X^2))/(X(B0+j)) > 0 since j >= 1 | EXACT | a t + (1-j)B0 = j(B0/X + X) checked as a polynomial identity |
| 3 | p < (a+e)/(ae) is p = (a+e)/(ae) * P/(1+P), P >= 0 | PROOF-READ | |
| 4 | G(chart) times AD, AN, E, 1+P powers: 46874 terms, all of one sign; E^2 coefficient 4 (= -4 in the document's sign) | EXACT | flint; no further common factor |
| 5 | walls J = 0, P = 0, Q = 0, B = 0 keep the sign and the E^2 term | EXACT | |
| 6 | G < 0 contradicts G = 0; equality p = d^2(1/a+1/e) is the capacity | PROOF-READ | |

## D. Falsification

NUMERIC: 10 random c = 0 points at p = 1/a + 1/e (seed 17), min Phi (8 restarts): worst -1.4e-16, so the chain is sharp. The theorem itself is exact, so no sampling of it is needed.

## E. Verdict

**confirmed-modulo-citations** (the chain atom only, which is proved in foundations/so4-atoms.md). No gap.
