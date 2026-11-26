# Audit: final_referee_T26_boundary10_uniform_20261022 (faces T_27 u {10} and T_27 u {13})

Evidence level: EXACT for all identities, the factorization (5)-(6) and the sign identity (own code from the paper's definitions); PROOF-READ for the logic; CITED for the four-coefficient capacity (1) and paper Lemma 4.2; NUMERIC test of the cited capacity.

PDF: proofs/manin_transfer_20261021/final_referee_T26_boundary10_uniform_20261022.pdf (2 pages). Auditor label C11-T27-T28.

## Deviations and limits

- The supplement's generator, checker and JSON were not used. The 14 x 14 Hessian determinant was recomputed exactly from (A.2) and the Appendix D generators.
- The cited capacity (1) is tested numerically only.
- This document is not needed for paper case T_27. The primary proof (final_referee_T26_uniform_20261022) proves all of T_27 directly and says no seven-contact theorem is required. Its Lemma 1 at j = 0 contains the present statement: with j = 0 the alternative pj >= d^2 is impossible, so K >= 0. The supplement's catalogue still lists this document under T_27.

## A. Statements

**Theorem 1 (p. 1).** For T' = T_27 u {10} = {2,3,4,7,8,10,12} and for T'' = T_27 u {13}: for every x in X with P c(x) >= 0 vanishing on T' (respectively T''), and every O in SO(4), Phi_{c(x)}(O) >= 0. Contacts are non-strict and all spectra are included. Neither 7-set is needed by the paper's case list, so this is a stand-alone strengthening.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| paper Lemma 4.2 / Prop. 4.1 (generators of the contact cone) | core paper | CITED |
| X contained in Sing C (paper App. A.1) | core paper | EXACT in lib |
| group G invariance (4.7): symmetry and transposition | core paper | EXACT in lib |
| (1) four-coefficient capacity a f67 + b f55 + e f53 + d f54, ab(e-d) >= d^2(a+b+e+d), for arbitrary nonnegative coefficients | other-doc | CITED; NUMERIC on the boundary: no violation |

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | F_{T'} = {0,1,2,6,53,54,55,67,70}; non-minor exactly 54, 70 | EXACT | lib |
| 2 | identities (2) and (3) | EXACT | 24 values |
| 3 | a symmetry fixing T' acts 0<->67, 1<->53, 6<->55, 54<->70, fixing 2 | EXACT | g = ((1,3,0,2),(1,3,0,2),1) |
| 4 | transposition sigma -> sigma^{-1} maps T' onto T'' | EXACT | lib `g_act_set` |
| 5 | reduction to (4) (= the primary document's (4) with j = 0) | PROOF-READ | as in the primary document |
| 6 | radial necessity: C(Re x) = 0, det H(Re x) = 0, Re x != 0 | PROOF-READ | as in the primary document, step 5 |
| 7 | C = -(A0 + eT), Q = -d(d+e) | EXACT | own expansion |
| 8 | (5)-(6): Q^4 det H(y0 - (C/Q) kappa) = -4 Q^2 C^2 G6^2, with G6 = A0{A0 - e(ab - T)} + p(a+c+d)(b+d+q)(A0 + eT), 206 monomials | EXACT | full determinant (3319 terms, degree 3 in z) |
| 9 | sign identity A0 - e(ab - T) = -K + (c+q)(a+d)(b+d) + cq(a+b+d) + e{d(c+q) + cq}; A0 >= d(a+d)(b+d) > 0; second summand >= 0 | EXACT | polynomial identity and coefficient signs |
| 10 | K < 0 implies G6 > 0, which contradicts G6 = 0; so K >= 0 | PROOF-READ | |
| 11 | apply (1) and add the minor terms; transposition gives T'' | CITED / PROOF-READ | |

## D. Falsification

NUMERIC: the cited capacity (1) on its boundary, 30 samples, worst min Phi / max gap = 1.6e-16 (no violation). The sign claim in step 10 is an exact identity with visibly positive terms, so no sampling is needed.

## E. Verdict

**confirmed-modulo-citations** (the cited four-coefficient capacity; paper Lemma 4.2). No gaps. Note again that the paper's proof of T_27 does not depend on this document.
