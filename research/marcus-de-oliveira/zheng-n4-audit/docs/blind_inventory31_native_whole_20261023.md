# Audit: blind_inventory31_native_whole_20261023 (primary proof of paper case T_32)

Evidence level: EXACT for every algebraic identity, cone computation and group action the proof uses (own code, `replay/blind_inventory31_native_whole_20261023__check.py`, all checks pass); PROOF-READ for the analytic steps (curve, punctures, irreducibility, Rouché restoration); CITED for the boundary contact cases, Appendix B/C theorems, the bounded-negative-part harmonic lemma and Theorem 1.2; NUMERIC falsification on 26 stratum points.

PDF: `proofs/manin_transfer_20261021/blind_inventory31_native_whole_20261023.pdf`, 5 pages. Auditor label C13-T32-T34.

## Deviations and limits

- The document does not name which earlier theorem covers each of the ten added contacts in its list (6); it defers to its certificate. I resolved every one of them myself with the symmetry group (section B). I did not read or run the certificate.
- Boundary cases, the one-circle theorem (Appendix C, Thm C.1), the shared-zero theorem (Appendix B, Thm B.1/Lemma B.2), the harmonic lemma and Theorem 1.2 were not audited here. They are CITED.
- Numeric falsification used the library's stratum sampler (`random_point_on_contact_stratum`, margin 0.01, so the contact set is exactly T) and a 25-restart BFGS minimiser over SO(4). It gives an upper bound for the minimum and is evidence, not proof.
- Paper notation is used throughout. The document's "plus brackets" mean sigma-points F_sigma = prod_i (A_i + B_sigma(i)); a point at infinity contributes the factor 1 (homogeneous coordinates). I checked that this reading reproduces every displayed factorisation exactly.

## A. Statements

**Theorem 1 (p. 1).** Let T = {0, 2, 3, 7, 11, 12} (one-line 1234, 1324, 1342, 2143, 2431, 3124). For every x in X and every reference permutation, if P c(x) >= 0 and (P c(x))_pi = 0 for all pi in T, then Phi_{c(x)}(O) >= 0 for all O in SO(4) (and, via Theorem 1.2, the unitary statement). Contacts are non-strict ("contain T"); all spectra are included: repeated eigenvalues, a_i + b_j = 0, extra contacts, zero coefficients.

- T equals the paper's T_32 literally (group element = identity; EXACT).
- The statement is at least as strong as (*) for T_32. Degenerate spectra are dispatched inside the proof to the one-circle and shared-zero theorems (CITED). Extra contacts are dispatched to earlier cases (CITED, resolved in B).
- The boundary input J = T ∪ {23} = {0,2,3,7,11,12,23} is literally the paper's T^(7)_22 (Table 5).

## B. Dependencies

| dependency | kind | resolution (own code) | status |
|---|---|---|---|
| added contact 1 | paper case | T ∪ {1} ⊇ g(T_4) (and T_7, T_9; T_4, T_7 stage 1) | CITED |
| added contact 5 | paper case | ⊇ g(T_4) (stage 1) | CITED |
| added contact 6 | paper case | ⊇ g(T_29), g(T_7), g(T_36) (stage 1) | CITED |
| added contact 9 | paper case | ⊇ g(T_29), g(T_36) (stage 1) | CITED |
| added contact 13 | paper case | ⊇ g(T_19), g(T_21) (stage 1) | CITED |
| added contact 16 | aux 7-set / case | T ∪ {16} is in the orbit of T^(7)_25 (stage 7); also ⊇ g(T_79) (stage 8) | CITED |
| added contact 17 | paper case | ⊇ g(T_31), g(T_36) (stage 1) | CITED |
| added contact 18 | paper case | ⊇ g(T_29) (stage 1) | CITED |
| added contact 20 | paper case | ⊇ g(T_40) (stage 1) | CITED |
| added contact 21 | paper case | ⊇ g(T_34) (stage 1) | CITED |
| added contact 23 | aux 7-set | T ∪ {23} = T^(7)_22 literally (stage 1) | CITED |
| added contacts 4, 8, 10, 14, 15, 19, 22 | none | both non-minor profiles f18, f40 are > 0 there (EXACT), so a zero gap leaves a pure minor sum | EXACT |
| repeated eigenvalue in one spectrum | Appendix C Thm C.1 (one-circle) | at most three distinct points, which always lie on a generalised circle (PROOF-READ) | CITED |
| a_i + b_j = 0 | Appendix B Thm B.1 (shared zero) | as quoted by the document | CITED |
| harmonic lemma | common lemma (final_referee_bounded_bad_harmonic_lemma; paper Lemma 3.5) | hypotheses checked in step 13 | CITED |
| Theorem 1.2 | paper | only needed for the unitary restatement | CITED |

Every case used is at a stage below 11 in the paper's Table 6 (T_32 is stage 11), so the order is non-circular. In each row the group element maps a 6-subset of T ∪ {j} other than T onto the cited case, and the statement for that case covers costs whose contact set contains the image, so it applies.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | T = paper T_32 | finite | EXACT | `find_g` returns the identity |
| 2 | Contact rows have rank 6; f0,f1,f2,f7,f8,f12,f15,f39 vanish on T and span the 8-dim kernel; the gap formulas (1) | identity | EXACT | integer profiles from the library |
| 3 | The identities (2): f18+f0 = f1+f2+f7+f15 and f40+f1 = f0+f7+f39 | identity | EXACT | profile vectors |
| 4 | Every supported gap vector is a nonnegative combination of the ten-list (3); this list equals `generators(T)` | case split | PROOF-READ + EXACT | x0, x1 cannot both be negative (g23 = x0+x1); the substituted coefficients are the gaps g23, g21, g20, g18 (resp. g23, g16, g17) |
| 5 | f18 and f40 take value 1 at indices 4 and 10, so d18+d40 <= min(g4, g10) and -L_O <= M min(g4, g10) | inequality | EXACT + PROOF-READ | profile values; the minor terms are >= 0 by (4.2) |
| 6 | Extra contacts: pure minor at {4,8,10,14,15,19,22}; the others map to earlier cases | finite + citation | EXACT (group) / CITED (cases) | table in B |
| 7 | Normalisation A=(a,z,x,w), B=(t,0,1,inf); E, rE, sE, uE factorisations | identity | EXACT | sympy |
| 8 | Inverse parametrisation (7), both directions; Dt, t+z, a·Da = x(1+uW), Da(u=-1/W) = (1+sW)(w+1) | identity | EXACT | rational function field |
| 9 | All divisions in (7) are nonzero at distinct points | argument | PROOF-READ | uses the step 8 identities and distinctness |
| 10 | H2,H0,H3,H12,H7 = 0,1,r,s,u; contacts ⇔ H_i real; g_i = -η Im H_i | identity + argument | EXACT + PROOF-READ | α = iη/E is forced by the contact at 0 |
| 11 | H11 - h = C/D: C has 43 terms and bidegree (2,2); D as in (8); a - z = -WxD/(Da(1-xW)), so D ≠ 0 at distinct points; H4, Δ = H11 - H10, F2/E = 1/W | identity | EXACT | field arithmetic |
| 12 | Disc_x C at (2,3,4,7) = 3(1452W^4+…+35), squarefree; (10) ≠ 0 there; squarefree degree-4 discriminant implies C absolutely irreducible | identity + argument | EXACT + PROOF-READ | content and factorisation both force a square factor in the discriminant |
| 13 | Puncture atlas: (11), C_W(-1,-1) = -r(h-1)(s-1); second root W1 ≠ -1, -1/s under (10); [x^2]C, Res = -rs(r-u)^2(s-u)(u-1)(r+s-u+1); [W^2x^2]C = -ru(s+1); [W^2]C, x*, D/W^2 → rs x*, Δ → -u; weighted equation, V0, D/W → su(r+1)/(h+u); the infinity equation has a simple root at x=0 | identity | EXACT | all displayed formulas reproduced |
| 14 | The atlas covers every point of the compact normalisation | analytic | PROOF-READ | every point lies over (W,x) in P^1×P^1: finite with x ≠ -1 (H4 bounded), x = -1 (H10 or the corner branch), x = ∞, or W = ∞ (finite x); the joint corner is absent |
| 15 | H4 is nonconstant on the curve | argument | PROOF-READ | the numerator of H4 - k has W-degree 1 and Wx^2 coefficient -r ≠ 0 |
| 16 | Harmonic lemma excludes a generic negative component | citation | CITED | hypotheses: H_i meromorphic and rational in (W,x); L_O is a real combination of gaps; boundary contacts by step 6; d bounded by step 14; nonconstant by step 15 |
| 17 | C ≡ 0 iff r = s = u = -1, h = 1 (given rs ≠ 0); on (13), H10 ≡ 0, so g10 = 0 and the cost is a pure minor sum | identity | EXACT | Groebner basis [λ-1, r+1, s+1, u+1, h-1]; H10 ≡ 0 |
| 18 | Nongeneric strict counterexamples perturb to generic ones (Rouché) | analytic | PROOF-READ (details supplied by auditor) | C(·;p0) ≢ 0 off (13); restricted to a suitable complex line the zero is isolated, so by Hurwitz it persists for real p near p0; generic parameters are dense (complement of a proper real algebraic set, nonempty by the witness); the inverse (7), strict gaps, distinctness and Φ_O < 0 are open conditions; H11 = h holds by construction |
| 19 | Unitary restatement | citation | CITED | Theorem 1.2 |

## D. Falsification

NUMERIC: seeds 100–126 of `random_point_on_contact_stratum(T, margin=0.01)` (26 points found, 1 failure). On each point I minimised Phi over SO(4) with 25 BFGS restarts. Worst min Phi / sum(gaps) = -1.0e-17 (zero up to rounding; Phi vanishes at signed permutation matrices of contacts). No violation. The log is appended to the check output.

## E. Verdict

**confirmed-modulo-citations.** Every algebraic identity in the document replays exactly. I supplied two pieces it omits: the resolution of the boundary contacts to named cases, and the details of the Rouché/density step. Citations: paper cases T_4, T_7, T_19, T_29, T_31, T_34, T_36, T_40 (and T_79 as an alternative for 16), T^(7)_22, T^(7)_25; Thm C.1; Thm B.1; the harmonic lemma; Thm 1.2.

Presentation notes (not gaps): (i) the boundary targets of (6) are named only in the certificate; (ii) the Rouché paragraph does not say why generic parameters are dense near p0. Both are filled in above.
