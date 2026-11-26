# Audit: final_referee_T41_negative_Heron_region_20261022

Evidence level: EXACT for the cubic restriction, the discriminant identity (1) and gamma; EXACT (randomized identity test: exact rational arithmetic at 4 to 6 random rational points) for the Hessian determinant identity, the principal and bordered minor identities, the resultant and the two square-class congruences (2); PROOF-READ for the rank and definiteness argument; NUMERIC for the rank-12 check at real roots.

PDF: proofs/manin_transfer_20261021/final_referee_T41_negative_Heron_region_20261022.pdf (3 pages). Auditor label C17-T42. Replay: `replay/final_referee_T41_negative_Heron_region_20261022__check.py` (17 checks, all pass, about 12 s).

## Deviations and limits

- Several polynomial identities are large: P has about 2000 terms, and L7, M7 and W14 have hundreds each. I did not expand them symbolically. I checked them at random rational base points instead, in exact arithmetic. At each point every quantity becomes a univariate polynomial in j, rebuilt by exact interpolation of exact 12x12 and 13x13 determinants. A false polynomial identity passes such a test only with negligible probability. Still, this is a randomized identity test, not a symbolic expansion. The term counts (393, 333, 9868, 6541, 1675, 531) were not reproduced.
- G is not written out in this note. I used the radial factor of the companion note final_referee_T41_bad17_uniform, built from its displayed formulas (G = V D^2 + R0 j D + A0 j^2). I then checked, independently, that it satisfies det H_F = -4(c+1)^2 G^2/Q^2 for the paper's cubic (A.2) on this source.
- The supplement's scripts were not run or read as evidence. The matching coordinates are the paper's (A.2). The cubic restriction C + hQ comes out exactly as displayed, with scale factor 1, so the coordinates and the f_j agree with the audit library.
- I did not test the main statement numerically. It is an impossibility statement about spectrally realizable points.

## A. Statement

**Proposition 1 (p.1).** Take the source g = q f3 + p f1 + j f53 + c f6 + e f15 + d f17 + b f16 + a f4, with all coefficients >= 0 and d > 0. All eight generators are generators of paper T_42 = {1,2,3,6,7,8}, the generator set is F_T = {0,1,2,3,4,5,6,15,16,17,53,54}, and f17 = h_{3421,+} is the only non-minor (EXACT). Assume Fe < 0 and Fq < 0, where

Fe = ((e+d+b)(p-a)+eb)^2 - 4bp(e-a)(e+d+b),  Fq = ((q+d+a)(p-b)+qa)^2 - 4ap(q-b)(q+d+a).

Then there is no x in X (X is the paper's matching variety) with Re x = g + h kappa for some real h. Equivalently, no c(x) with these coefficients is spectrally realizable. The statement covers every spectrum, because the argument uses only grad C(x) = 0. The strict signs force p, a, b, e, q > 0 and j > 0. Cases with Fe = 0 or Fq = 0, or with both quartics positive, are excluded from the statement.

## B. Dependencies

- The paper: X ⊂ Sing C, with the cubic (A.2) (replayed in lib, EXACT). This is the only property of X that is used.
- Companion note final_referee_T41_bad17_uniform: the definition of G only, which I re-derived and checked here. Nothing else is used.
- External: none. The note correctly avoids the HMSV rank theorem.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | normalise d=1; F(g+h kappa) = C + hQ with C = -(c+1)(ab+a+ae+b+bq+1+e+j+q+ej+eq+jq) < 0, so Q != 0 and h = -C/Q (p.1) | identity | EXACT | own cubic restriction |
| 2 | x = w + i delta in Sing C gives grad F(w) = grad F(delta), H(w) delta = 0 and F(w) = 0; delta != 0 because Q(w) != 0 (pp.1-2) | analytic | PROOF-READ | expansion of a cubic at w + i delta (auditor) |
| 3 | det H_F(w) = -4(c+1)^2 G^2/Q^2 (p.1) | identity | EXACT (randomized, 6 pts) | 14x14 exact determinants |
| 4 | gamma = V (e+1+b)^2 (c+1)^2 (q+1+a)^2 > 0, so j > 0 at every root (p.1) | identity | EXACT | symbolic |
| 5 | (1): beta^2 - 4 alpha gamma = (e+1+b)^2 (c+1)^2 (q+1+a)^2 Fe Fq > 0, so G_j != 0 at roots, including a degree drop (p.1) | identity | EXACT | symbolic |
| 6 | the principal block on I = {2..13} does not depend on h; its determinant is j^2 P(j) with P quadratic (p.2) | identity | EXACT (randomized, 4 pts); height-independence EXACT | symbolic Hessian; interpolation |
| 7 | Res_j(G,P) = W^2; T + 2 alpha W = Fe (c+1)^2 (q+1+a)^2 L7^2; T - 2 alpha W = Fq (c+1)^2 (e+1+b)^2 M7^2; P(T ± 2 alpha W) ≡ (U ± W)^2 mod G (p.2, eq. (2)) | identity | EXACT (randomized, 4 pts) | rational-square tests and exact polynomial remainders |
| 8 | U = alpha P - a2 G, so P = 0 on G = 0 forces W = 0; where W L7 != 0, P < 0 (p.2) | identity + inequality | EXACT (randomized) + PROOF-READ | |
| 9 | the bordered determinants on I ∪ {0} and I ∪ {1} have numerators divisible by G after h = -C/Q; Schur gives rank exactly 12 on the dense set (p.2) | identity + linear algebra | EXACT (randomized, 4 pts) + PROOF-READ | interpolation in j of Q^2 det |
| 10 | nonpositivity of all 12x12 principal minors and rank <= 12 extend by continuity to the whole local root sheet; W and L7 are nonzero polynomials in the base (p.2) | analytic | PROOF-READ | equal-sign property of the rank-r principal minors of a symmetric matrix of rank r (classical) |
| 11 | d²/dj² det H(w(j)) = -8(c+1)^2 G_j^2/Q^2 < 0, so the corank is at most 2 and the rank is exactly 12 everywhere (p.2) | analytic | PROOF-READ | follows from step 3 and G = 0 |
| 12 | Schur expansion: second derivative = 8 det A' det(Q_v restricted to K), with Q_v(x) = v · grad F(x) and v = w_j; this forces Q_v to be definite on K; but Q_v(delta) = v · grad F(w) = d/dj F(w(j)) = 0, a contradiction (p.2) | analytic | PROOF-READ | auditor re-derived it: x^T H(v) x = 2 Q_v(x); the Schur complement derivative is K^T H(v) K |

## D. Falsification

At five random real roots of G (real j, random rational base points, mpmath with 60 digits), the three smallest singular values of H_F were about 1e-50, 1e-50 and 0.02-5. So the rank is exactly 12, as claimed. NUMERIC, seed 993.

## E. Verdict

**confirmed.** Every step is EXACT (some by randomized exact identity test, as stated above) or PROOF-READ. There is no external dependency. No gap was found. The scope is limited: Fe = 0, Fq = 0 and both quartics positive are outside the statement, as the note says.
