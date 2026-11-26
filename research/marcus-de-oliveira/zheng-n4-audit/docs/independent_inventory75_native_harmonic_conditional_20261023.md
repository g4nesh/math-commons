Evidence level: EXACT for every displayed identity, inverse, resultant, leading coefficient, branch expansion and combinatorial claim (own sympy code, 0 failures); PROOF-READ for the analytic steps; CITED for the seven-contact hypothesis (supplied by the boundary16 document), the Table 5 set j=27, the earlier six-cases, Appendix B and Theorem 1.2; NUMERIC only for the falsification run.

# independent_inventory75_native_harmonic_conditional_20261023 (conditional proof of T_76)

PDF: `proofs/manin_transfer_20261021/independent_inventory75_native_harmonic_conditional_20261023.pdf`, 4 pages.
Replay: `replay/independent_inventory75_native_harmonic_conditional_20261023__check.py` (with `--numeric`).

## Deviations and limits

- No deviation. All formulas were rebuilt from the definitions with the document's conventions (A = (y,z,w,x), B = (t,0,1,∞), homogeneous brackets, reference permutation 0, E = F1 − F0).
- The curve C is −1 times my reduced N12 − h·D12. The resultants agree with ratio 1 after that normalization.
- The J = ∞ branch expansions were computed as truncated power series to order τ⁴, with exact symbolic coefficients. The limit of H14 + J uses only the coefficients up to the leading order. The code checks that these coefficients do not involve unsolved series terms.
- The document does not list its generic exclusions (p.2 "generic real parameters"; p.4 "proper finite algebraic conditions"). I reconstructed the list the argument actually uses (E below).

## A. Statements

**Theorem 1 (p. 1), conditional.** Assume (H): every x ∈ X whose cost is supported and vanishes on S16 = {0,1,3,10,11,12,16} has Φ ≥ 0 on O(4). Then every x ∈ X with P c(x) ≥ 0 and (Pc)_σ = 0 for σ ∈ T = {0,1,3,10,11,12} has Φ_{c(x)}(O) ≥ 0 for all O ∈ O(4), and hence for U(4). Contacts are non-strict. All complex spectra are covered: the non-eight-distinct ones go to the collision theorems (p.1). T is literally paper T_76 (EXACT). S16 is literally Table 5, j = 18.

## B. Dependencies

| kind | ref | used for | match checked |
|---|---|---|---|
| seven-contact | T(7)_18 = S16 (hypothesis H; proved in final_referee_inventory75_boundary16_uniform_20261023, stage 4) | seventh tie 16 | yes, literal |
| seven-contact | T(7)_27 = {1,4,5,7,9,12,14} (root_inventory67_g14 document, stage 5) | seventh tie 14 | yes: T∪{14} = g·T(7)_27 with g = (rho=0213, rho'=0312, eps=0) in 0-based one-line form; no six-case of stage ≤ 5 closes it |
| paper-case | T_9, T_11, T_17, T_38 | seventh tie 5 | yes (stage ≤ 5) |
| paper-case | T_39, T_44, T_38 | seventh tie 6 | yes |
| paper-case | T_29, T_31 | seventh tie 7 | yes |
| paper-case | T_46 | seventh tie 13 | yes |
| paper-case | T_36, T_40, T_38 | seventh tie 20 | yes |
| common-lemma | Appendix B (Theorem B.1, Lemma B.2), one-circle theorem | collisions | yes (they give MOC at every such spectrum) |
| core-doc | paper (3.8), (3.11), (3.12) (H_κ rank 4, ker H_κ = E, grad C(κ) = 0) and Birkhoff | height-only boundary K = 0 | yes: lib replays (3.10)-(3.12); rank 4 re-checked |
| core-doc | Lemma 4.2 | complete source | yes |
| external | Theorem 1.2 | unitary | yes |

T_76 is at stage 6 of Table 6. All of the above are at stage ≤ 5.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | source [0,1,2,8,10,12,15,30], only f30 non-minor; coefficients are the exclusive gaps g20,g16,g13,g14,g7,g6,g5,g2 = d (1) (p.1) | finite | EXACT | |
| 2 | zeros of f30 off T are 5,6,7,13,14,16,20, closed as in B; any other extra tie forces d = 0 (p.1) | finite + citation | EXACT + CITED | |
| 3 | E = −z(y−t)(x−w); r, s (2); the inverse (3); the t numerator minus its denominator is s(x−1)(w−1); coverage (pp.1-2) | identity + argument | EXACT + PROOF-READ | |
| 4 | (4) H0 = 0, H1 = 1, H3 = r+1, H10 = u, H11 = s+u (p.2) | identity | EXACT | |
| 5 | coordinates J, X; homogeneous inverse; (5) E = s(X−1)(J−r)LKA_N, F0/E = (r+1−JX)/K, H2 = J, H5 = r+1−JX (p.2) | identity | EXACT | the inverse also reproduces every original configuration exactly (auditor check) |
| 6 | every H_σ denominator divides (J−r)(X−1)LA_N (p.2) | identity | EXACT | all 24 |
| 7 | curve: 51 terms, bidegree (3,2); [J³]C, [X²]C, C(J,1), the resultants with L, A_N, K (6), P_H (pp.2-3) | identity | EXACT | |
| 8 | J − r, A_N punctures have real J; X = 1: H5 + J → r+1 (p.2) | identity + argument | EXACT + PROOF-READ | the A_N resultant roots are 0, u, r+1, r−s |
| 9 | L = 0 with nonreal J: (7) the limits of H6, H13, H20 (none of their denominators contains L); \|J−r\|² = rs; H13 = (r+1)h/J; (8); contradiction (p.3) | identity + argument | EXACT + PROOF-READ | |
| 10 | X = ∞ with finite J: J = 0 or f∞(J) = 0; H12 − (H6 + J) → 0 (p.3) | identity | EXACT | in fact the limit is 0 for every finite J |
| 11 | J = ∞: two simple branches X → 0 (slope (V−r)/V) and X → V/s; there is no (∞,∞) point since [J³X²]C = s; H14 + J → h − s, real, on both branches (9)-(11) (p.3) | identity | EXACT | truncated series |
| 12 | height divisor K = 0: all H are finite there; F0/E diverges except at J = r+1; H_κ v0 = 0, so v0 ∈ E, and Birkhoff gives Φ ≥ 0 at the limit (p.3) | identity + argument | EXACT + PROOF-READ | cubic expansion ∇C(fκ+v) = f·H_κ v + ∇C(v), using ∇C(κ) = 0 |
| 13 | H12 is nonconstant for every r, s ≠ 0 (numerator J-degree 3, denominator J-degree 2); Rouché restoration (p.4) | identity + analytic | EXACT + PROOF-READ | r, s ≠ 0 at eight-distinct points (x ≠ z; y ≠ 0; t ≠ 1) |
| 14 | harmonic completion (p.4) | analytic | PROOF-READ | as in the T_75 document |

## D. Falsification

- **Main statement (brief).** 120000 random draws of (r,s,u,h) and a complex J, solving C = 0 for X: 146 strictly supported points. On 40 of them, BFGS over SO(4) with 6 restarts: worst min Φ / max gap = 2.1e-15 (rounding; the minimum is 0). Seed 20261125.
- **Intermediate claims.** These are all exact identities here. The one sign argument, (8), is elementary and was replayed symbolically.

## E. Verdict

**confirmed-modulo-citations**, as a conditional theorem. Its hypothesis is discharged by final_referee_inventory75_boundary16_uniform_20261023, audited in the same batch.

Presentation gap (p.2 and p.4): the generic exclusions are not listed. The steps actually need u ≠ 1, s + u ≠ 0 and V = s + u − 1 ≠ 0. These are the factors of the resultants and of the J = ∞ leading coefficients. r, s ≠ 0 hold automatically at eight-distinct points. All of these are nonzero at (2,3,4,5), so the exceptional set is proper and the Rouché restoration applies. A component with constant J needs no exclusion: the X = ∞ barrier holds for every finite J. A component with constant X is excluded once V ≠ 0, because then the two J = ∞ branches are simple and not vertical.

Also p.2: "Res_X(C, J − r) = (J − r)²" is a resultant with an X-free polynomial and carries no information. The point it is meant to support (J = r is real) is immediate.
