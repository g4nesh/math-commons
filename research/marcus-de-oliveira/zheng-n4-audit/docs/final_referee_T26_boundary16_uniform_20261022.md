# Audit: final_referee_T26_boundary16_uniform_20261022 (face T_27 u {16})

Evidence level: EXACT for the generator list, the identities, the symmetry, C, Q, the factorization (3), the specialisation (4), both 136-coefficient endpoints and the 525-coefficient certificate (5), all by own code from the paper's definitions; PROOF-READ for the logic; CITED for the SOC capacity and paper Lemma 4.2; NUMERIC test of the cited SOC.

PDF: proofs/manin_transfer_20261021/final_referee_T26_boundary16_uniform_20261022.pdf (3 pages). Auditor label C11-T27-T28.

## Deviations and limits

- G6 was not taken from the supplement's JSON. It was obtained as the exact square root of -Q^4 det H / (4 Q^2 C^2) from my own determinant. Its sign is fixed by matching (4), since G6 enters only through G6^2.
- The cited SOC is tested numerically only.
- This document is not used by the primary T_27 proof. Its statement is the primary document's Lemma 1 at e = 0: there K = -abd - d^2(a+b+d) < 0, so pj >= d^2. With the names p -> a, j -> b, b -> h, a -> w this is exactly ab >= d^2 here. I checked by hand that the two forms of C agree.

## A. Statements

**Theorem 1 (p. 1).** For T' = T_27 u {16} = {2,3,4,7,8,12,16}: for every x in X with P c(x) >= 0 vanishing on T', and every O in SO(4), Phi_{c(x)}(O) >= 0. Contacts are non-strict and all spectra are included.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| paper Lemma 4.2 / Prop. 4.1 | core paper | CITED |
| X contained in Sing C | core paper | EXACT in lib |
| G invariance (4.7) | core paper | EXACT in lib |
| "old SOC": a f2 + b f5 + d f54 safe when ab >= d^2, for arbitrary nonnegative coefficients (the same as input (1) of the primary document) | other-doc | CITED; NUMERIC: no violation |

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | F_{T'} = {0,2,5,6,54,55,67,70}; non-minor 54, 70 | EXACT | lib |
| 2 | f54 + f70 = f0 + f6 + f55 + f67 | EXACT | |
| 3 | a symmetry fixing T' acts 54<->70, 2<->5, 55<->67, fixing 0 and 6 | EXACT | g = ((0,2,1,3),(0,2,1,3),0), the document's row/column permutation |
| 4 | reduction to (2) with d > 0 | PROOF-READ | |
| 5 | radial necessity C(Re x) = 0, det H(Re x) = 0 | PROOF-READ | as in the primary document |
| 6 | Q = -ab - d^2 < 0; C as displayed, < 0 | EXACT | |
| 7 | (3): Q^4 det H(y0 - (C/Q) kappa) = -4 Q^2 C^2 G6^2, G6 of 277 monomials, degree 6 | EXACT | determinant 7298 terms; the quotient is a perfect square |
| 8 | (4): G6 at q = c = h = w = 0 equals -d^2(a+d)(b+d)(ab-d^2) | EXACT | |
| 9 | G6 at a = 0 and at b = 0: 136 coefficients each, all strictly positive | EXACT | |
| 10 | (5): a^2(1+u)^2 G6(q, a, d^2/(a(1+u)), c, d, h, w) has 525 monomials, all positive (deg_b G6 = 2) | EXACT | |
| 11 | hence G6 > 0 whenever ab < d^2 (a, b > 0 via 10, a = 0 or b = 0 via 9, using d > 0) | PROOF-READ | every case of ab < d^2 is covered |
| 12 | apply the SOC, add the 4 minor terms, undo the symmetry | CITED | |

## D. Falsification

NUMERIC: the cited SOC on ab = d^2, 30 samples, worst min Phi / max gap = 1.0e-16 (no violation). The sign claims are exact coefficient certificates.

## E. Verdict

**confirmed-modulo-citations** (the cited SOC; paper Lemma 4.2). No gaps. This document is not needed for the T_27 proof.
