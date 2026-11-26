# root_moc4_T7_boundary23_20261022 (seven-contact face T_8 ∪ {23})

Evidence level: EXACT for the lift, C, Q and every sign claim on G. EXACT-RANDOMIZED for the Hessian factorization: exact rational arithmetic at 8 random points, which is a probabilistic identity test and not a symbolic expansion. PROOF-READ for the singular-locus argument. CITED for X ⊂ Sing C (paper; lib-replayed) and the SOC atom.

PDF: proofs/manin_transfer_20261021/root_moc4_T7_boundary23_20261022.pdf (1 page). Replay: `replay/root_moc4_T7_boundary23_20261022__check.py` and its output (needs the environment variable MOC4_SUPP).

## Deviations and limits

- **G is taken from the supplement.** The PDF does not display the quintic G ("saved in the certificate"). It was loaded as data (json.load) from the supplement file root_moc4_T7_boundary23_radial_20261022.json. Every property of G was then checked with my own code. The JSON's C and Q equal the document's.
- **The factorization is checked at points only.** The identity Q^4 det H(c0 - (C/Q)kappa) = -4(d+v)^2(d+w)^2 Q^2 (c+d+h)^4 G^2 was checked by computing det Hess of the paper's cubic C exactly at 8 random rational points. The ratio was 1 every time. It was not expanded symbolically (the 14x14 determinant in 7 parameters was too costly for the budget).
- **Numeric test.** G > 0 at 2000 random points with ab < d^2 (min 1.8e-4).

## A. Statements

- **Theorem 1 (p. 1).** Let x be in X with g >= 0 and contacts ⊇ {0,2,5,6,7,8,23}. Generators: {2,3,4,5,39,45,46}, checked. Then Phi >= 0 for all O. In the nonminor case (d > 0), ab >= d^2. The argument uses X directly, so it covers every spectral configuration, collisions and repeated labels included.

## B. Dependencies

- X ⊂ Sing C, i.e. grad C vanishes on X: paper Appendix A, lib-replayed. CITED.
- SOC atom (F03) on (f2, f39; f46): CITED.

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | the generators and the canonical expansion with 7 rays | EXACT |
| 2 | C(c0 + z kappa) = C + zQ with the document's C and Q (constant 1, using the paper's cubic) | EXACT |
| 3 | from grad C(c + i delta) = 0: C(c) = 0 and det H(c) = 0 | PROOF-READ (polarization and Euler's identity, checked line by line) |
| 4 | C > 0 forces Q ≠ 0 and z = -C/Q | PROOF-READ |
| 5 | Q^4 det H = -4(d+v)^2(d+w)^2 Q^2 (c+d+h)^4 G^2 | EXACT-RANDOMIZED (8 points) |
| 6 | a = 0 or b = 0: all coefficients of G are positive, with d^5 present | EXACT |
| 7 | a^2(1+u)^2 G(..., b = d^2/(a(1+u)), ...) has 66 terms with positive integer coefficients; G(a,c,0,0,b,0,d) = -d(a+d)(b+d)(ab-d^2) | EXACT |
| 8 | conclusion: ab >= d^2, so the SOC atom applies | PROOF-READ + CITED |

## D. Falsification

G > 0 numerically on 2000 random points with ab < d^2 (seed 5; minimum 1.8e-4). No violation was found.

## E. Verdict

**confirmed-modulo-citations.** There is no gap. The one step not replayed symbolically (step 5) was checked by exact arithmetic at random points, and I believe it holds. A symbolic expansion would raise it to EXACT.
