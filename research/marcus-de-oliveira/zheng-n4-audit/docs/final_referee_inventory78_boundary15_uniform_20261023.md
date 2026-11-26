Evidence level: EXACT for the source, the boundary pullbacks, F(w0 + z kappa) = C + zQ, and the rational point (3) (Q, det A, the restored vector, inertia (7,5,2), rank 12). For the claims about the unprinted 584-term polynomial G (identity (1), property 1, factorization (2), the bordered identities of Lemma 3), the evidence is exact arithmetic on random slices only. Those claims are NUMERIC in the sense of the brief. Lemmas 2-4 and the theorem are PROOF-READ. Prop. C.6 and radial necessity are CITED.

# final_referee_inventory78_boundary15_uniform_20261023 (seven-set T^(7)_25 = {1,3,4,8,11,12,15}, used for T_79 and T_89)

PDF: proofs/manin_transfer_20261021/final_referee_inventory78_boundary15_uniform_20261023.pdf (4 pp.).
Replay: `replay/final_referee_inventory78_boundary15_uniform_20261023__check.py` and its output (all checks pass).

## Deviations and limits

- **G was not reconstructed as a polynomial.** The PDF does not print G. I computed det H exactly at rational points. On 25 random positive rational choices of (j,c,e,d,b,a), I computed det H * Q^2 exactly as a polynomial in q, by exact interpolation plus an extra point. On every slice:
  - -det H * Q^2 / 4 is a perfect square G(q)^2 with deg_q G = 3;
  - every coefficient of q^k, k >= 1, is positive;
  - G(0) = d C(0,...) P(d) with the correct sign;
  - the number of positive q-roots is 1 if P(d) < 0 and 0 otherwise.

  On 6 further slices with P(d) < 0:
  - G divides Q^2 det H_{0^0^} and Q^2 det H_{1^1^};
  - the Hessian at the positive root has inertia (7,5,2) (60-digit mpmath eigenvalues; zero eigenvalues below 1e-60).

  This is strong evidence but not an exact proof of the global coefficient claims.
- My coordinates are the library's (paper Appendix A, cubic (A.2)). The restored vector at (3) and det A come out identical to the PDF's numbers. The supplement's "original cubic" therefore coincides with the paper's C in these coordinates (EXACT: C(w0 + z kappa) = C_doc + z Q_doc with the same sign).

## A. Statement

Theorem 1. For every x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in S = {1,3,4,8,11,12,15} (= Table 5 T^(7)_25, literally), and every O in SO(4): Phi_{c(x)}(O) >= 0. All spectra are included. In fact the proof shows that no x in X has all seven source coefficients positive. Every other x is in a sub-face already covered.

Lemmas 2-4 concern the real algebraic set R = {(q,j,c,e,d,b,a) > 0 : G = 0}. Lemma 2: R is a connected analytic q-graph over {0 < d < d*}. Lemma 3: rank H = 12 on R. Lemma 4: inertia (7,5,2) on R.

## B. Dependencies

| dependency | check |
|---|---|
| source: f0, f1, f9, f10, f56, f53, f59 with exclusive gaps 20, 17, 10, 9, 6, 14, 2; only f56 non-minor | EXACT (lib) |
| zero coefficient q, j, c, e, b, a (added 20, 17, 10, 9, 14, 2) -> document labels T39, T20, T70, T69, T75, T12 = paper T40, T21, T71, T70, T76, T13 | EXACT: each S ∪ {j} contains a G-image of the cited case, all at stage <= 6 (S is stage 7) |
| other added ties: all generators minor | EXACT |
| radial necessity (task item 5): C(Re x) = 0 and H(Re x) singular | CITED. The case Im(alpha x) = 0 is supplied by the auditor: H(x)x = 2 grad C(x) = 0 with x != 0 |
| Prop. C.6 (n+, n- <= 6 for H(Re x), x in X) | CITED (proved in foundations/appendix-c-part1.md, modulo standard GIT surjectivity). Inertia is invariant under congruence and a sign change of C, so the coordinate convention does not matter |
| c(x) = Re x - Re F_{pi*}(x) kappa, so Re x = w0 + z kappa | PROOF-READ (definition (3.3)) |

## C. Steps

| # | claim (page) | kind | status |
|---|---|---|---|
| 1 | source, boundary pullbacks (p.1) | combinatorial | EXACT |
| 2 | F(w0 + z kappa) = C + zQ with the displayed C, Q (p.2) | identity | EXACT. Linear in z, as expected from grad C(kappa) = 0 |
| 3 | spectral realizability gives z = -C/Q | analytic | PROOF-READ (radial necessity; Q > 0 on the orthant) |
| 4 | (1) det H = -4G^2/Q^2 with G of degree 9 | identity | NUMERIC (exact values at 12 random points are -4 times rational squares, plus 25 exact q-slices) |
| 5 | property 1: G_q > 0 for q >= 0, the others positive | coefficient sign | NUMERIC (25 exact slices) |
| 6 | property 2: G(0,...) = d C P(d) | identity | NUMERIC (25 exact slices, sign included) |
| 7 | Lemma 2 from properties 1-2 | analytic | PROOF-READ. P(0) < 0, P' > 0, so d* is unique; G(0) < 0 iff d < d*; G strictly increasing in q, so the positive root is unique and simple; the base is connected |
| 8 | Lemma 3, upper bound: Q^2 det H_{i^i^} = G R_i (i = 0,1); Schur-complement argument; density of {det A != 0} on connected analytic R | identity + analytic | identities NUMERIC (G divides both bordered minors on 6 exact slices); the argument is PROOF-READ and correct |
| 9 | point (3): G = 0 (q = 0 and P(1/2) = 0 there, i.e. d = d*; det H = 0 exactly), Q = 101/10, det A = -4239165411/12005000 | finite | EXACT (Q, det A, det H = 0) |
| 10 | Lemma 3, lower bound: d^2/dq^2 det H = -8G_q^2/Q^2 at roots; rank <= 11 would force det = O(eps^3) | analytic | PROOF-READ (correct; it rests on (1)) |
| 11 | Lemma 4: at (3), the restored vector as printed; 12 nonzero pivots (7 positive, 5 negative) and a zero 2 x 2 remainder | finite | EXACT (own rational elimination; the vector agrees entry by entry) |
| 12 | inertia constant on R (continuity at the rank-12 point (3), constant rank, connectedness) | analytic | PROOF-READ |
| 13 | Theorem 1: a positive source would give inertia (7,5,2) at H(Re x), contradicting n+ <= 6 | citation | CITED (Prop. C.6) |

## D. Falsification

- Inertia on the actual sheet R: 6 random sheet points, all (7,5,2) (NUMERIC, 60 digits).
- No positive-source point in X was sought. The theorem says there is none.
- The coefficient claims survived 25 exact slices.

## E. Verdict

gap-found, presentation level only. The one gap: the PDF states identity (1) and the coefficient-table properties 1-2 of an unprinted 584-term polynomial and refers to a companion script for the proof. Under the audit rules a script claim is not evidence. My own exact-arithmetic slices (steps 4-6, 8) agree in every case, and I believe the claims are true. Repair: print G (or R_0, R_1), or replay G exactly by interpolation (about 2000 monomials in total over its four q-coefficients, signs fixed by property 2). Everything else is EXACT or PROOF-READ, modulo Prop. C.6 and radial necessity (CITED).
