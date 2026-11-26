Evidence level: EXACT for every displayed identity, parametrization, resultant, coefficient, branch expansion and combinatorial claim (own sympy code, 56 checks, 0 failures); PROOF-READ for the analytic steps (chart coverage, Rouché perturbation, pole atlas logic, harmonic minimum principle); CITED for the earlier contact cases, Appendix B, the one-circle theorem and Theorem 1.2; NUMERIC for the falsification run.

# independent_inventory74_native_harmonic_uniform_20261023 (paper T_75)

PDF: `proofs/manin_transfer_20261021/independent_inventory74_native_harmonic_uniform_20261023.pdf`, 5 pages. Primary proof of paper case T_75 = {3, 5, 7, 8, 10, 12}.
Replay: `replay/independent_inventory74_native_harmonic_uniform_20261023__check.py` (run with `--numeric`), output in `..._check_output.txt`.

## Deviations and limits

- No deviation from the assigned scope.
- Every formula I replayed was rebuilt from the definitions (homogeneous brackets, the 24 products, gaps = −η Im H_σ). I did not use or read the supplement's companion script or JSON. The text extract was enough: every displayed formula I used matched my own derivation, so no page image was needed.
- Three formulas are only correct up to a constant factor, which is harmless: the curve C (my reduced N12 − h·D12 is −1 times the document's normalization), and the four resultants (8) (ratio exactly 1 after that normalization).
- The analytic part (compact normalization, minimum principle, Rouché) is read, not machine-checked. The earlier contact cases, Appendix B and the one-circle theorem are cited, not re-audited here.
- The numerical test of the main statement is brief, as the brief asks: 120000 random curve points, 269 strictly supported, of which 40 were minimized over SO(4).

## A. Statements

**Theorem 1 (p. 1).** Let x ∈ X, with π* any element of T and c = c(x) = Re(x − F_{π*}(x)κ), such that P c ≥ 0 and (Pc)_σ = 0 for every σ ∈ T = {3, 5, 7, 8, 10, 12}. Then Φ_c(O) ≥ 0 for every O ∈ O(4), and hence also for every U ∈ U(4) via Theorem 1.2. Contacts are non-strict ("include T"). Spectra are arbitrary complex: the non-eight-distinct spectra (a repetition inside one spectrum, or a_i + b_j = 0) are sent to the collision theorems (p. 1, p. 2 "Ordinary finite collisions"). This is exactly (*) for paper T_75. The tuple is literally the paper's T_75, with g = identity (EXACT).

**Lemma 2 (p. 3).** At every eight-distinct point with exactly the contacts T, H12 has a nonconstant holomorphic germ in the chart (Y, Z), for the fixed real parameters (D, k, u). Hence any such counterexample can be moved, keeping strict support off T, the contacts T, η and the sign of the response at a fixed O, to one whose real parameters (D, k, u, h) avoid any given proper real algebraic set.

**Pole atlas (pp. 3-5).** Assume generic parameters: D, k, u, u−1, r, s, Dk+k+u−h−1 are all nonzero, and the Y² and Y⁰ coefficients of C are coprime. Then at every puncture of the compactified normalized curve (6), every approach through strictly supported points has d = g6 → 0.

## B. Dependencies

| kind | ref | used for | match checked |
|---|---|---|---|
| paper-case | T_28 or T_39 (stage 1), also T_26, T_47 | seventh ties 2 and 13 | yes: T∪{j} contains a 6-subset in the orbit (code) |
| paper-case | T_25 (stage 1), also T_24, T_26 | seventh ties 4 and 16 | yes |
| paper-case | T_15 (stage 1), also T_16, T_17, T_46 | seventh tie 9 | yes |
| paper-case | T_65 (stage 6, the only closer) | seventh ties 14 and 18 | yes; stage 6 < stage 7 of T_75 |
| core-doc | Lemma 4.2 / Prop. 4.1 (complete source of T) | the decomposition with coefficients equal to exclusive gaps | yes: generators(T) = [0,1,2,6,39,59,67,70] |
| common-lemma | Appendix B, Theorem B.1 (cross collisions) and Lemma B.2 (repetitions), established; one-circle theorem (Appendix C), established | non-eight-distinct spectra; collision boundary points of Ω | yes: B.1 and B.2 give MOC(4) at every such spectrum, which is stronger than needed |
| external | Theorem 1.2 (audited in foundations) | from O(4) to U(4) | yes |
| external | maximum principle for harmonic functions; Rouché; normalization of curves | analytic steps | classical |

Paper Table 6 puts T_75 at stage 7. Every closer used is at stage ≤ 6, so there is no circularity. Extra zeros σ with f70(σ) > 0 force d = 0, which is the trivial minor case.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | complete source on T is f0,f1,f2,f6,f39,f59,f67,f70; only f70 is non-minor; the coefficients are the exclusive gaps g18,g16,g13,g9,g4,g2,g14,g6 (p.1) | finite | EXACT | lib generators; exclusive positions with value 1 |
| 2 | (1) g0 ≥ d and g1 ≥ d (p.1) | finite | EXACT | f70(π0) = f70(π1) = 1, all coefficients ≥ 0 |
| 3 | (2) the zeros of f70 off T are 2,4,9,13,14,16,18; each T∪{j} contains an earlier closed case (p.1) | citation | EXACT + CITED | closers listed in B |
| 4 | reduction to eight-distinct spectra with exactly T (p.1) | citation | CITED | Appendix B / one-circle |
| 5 | E = F5 − F3 = x(y−t)(z−w); D and k formulas (4) (p.2) | identity | EXACT | sympy |
| 6 | solving for y, t, w; the raw equation E(H10 − u) = x(x−1)t/((x−z)(x(D+1)−1))·[(w−z)L − (x−z)(z−1)]; L ≠ 0; coverage (p.2) | identity + argument | EXACT + PROOF-READ | identity exact; each denominator is nonzero at an eight-distinct point (y, x, t ≠ 0; x ≠ 1, z; the right-hand side Dxt ≠ 0) |
| 7 | homogeneous inverse in (Y, Z); E = −Z(Z−1)L(Y+1)K; H3,H5,H7,H8,H10; F3/E = −L/(Z−1) (5) (p.2) | identity | EXACT | sympy |
| 8 | the new inverse sends the coordinates of any configuration back to its image under ξ ↦ (ξ−1)/(x−1) (onto chart) | identity | EXACT | auditor's own check, stronger than the text |
| 9 | every H_σ denominator divides Z·L·(Y+1)·K; no zero spinor off E = 0 (p.2) | identity + argument | EXACT + PROOF-READ | all 24 checked |
| 10 | (7): H0, H1, H2, H4, H9, H6 = N6(Z)/L; leading term of H6 (p.3) | identity | EXACT | |
| 11 | curve (6) has bidegree (2,4) and 160 terms (p.3) | identity | EXACT | 160 terms in (Y,Z,D,k,u,h) |
| 12 | Lemma 2: C ≢ 0 at an eight-distinct point (p.3) | identity + argument | EXACT + PROOF-READ | coefficient chain replayed; independently, the full solution set of "C ≡ 0" is {u=1, k=−1/D, h=D+1}, i.e. u = 1 and s = 0, which forces A2 = A3 |
| 13 | Lemma 2: Rouché on a complex line keeps a nearby real-parameter solution (p.3) | analytic | PROOF-READ | the germ is holomorphic at an eight-distinct point because E ≠ 0 there; the strict inequalities are open |
| 14 | generic exclusions are proper; the specialization (2,3,4,5) satisfies them; coprimality excludes a constant-Z component (p.3) | finite + argument | EXACT + PROOF-READ | values [2,3,4,3,6,10,7]; gcd of the Y², Y⁰ coefficients is constant. No constant-Y component either: Res_Y(C, Y+1) ≢ 0 (auditor) |
| 15 | resultants (8) (p.4) | identity | EXACT | ratio 1 |
| 16 | Z = 0: H6(0) = s, real; Y = −1 and K = 0 give real finite Z, so H6 is real unless L = 0 (p.4) | identity + argument | EXACT + PROOF-READ | |
| 17 | L = 0, Y = −1 (Z = s/u): C_Y = −Dk²r²s²/u⁴ and C_Z = krs/u, both nonzero; (Y+1)/L → u²/(Dkrs); H1 → −Dkr/u, real; then (1) gives d → 0 (p.4) | identity + argument | EXACT + PROOF-READ | L = 0 meets the finite curve only at Y = −1, since Res_Z(C, L) ∝ (Y+1) |
| 18 | Z = 1 (height): L = −Dk ≠ 0, so H6 is finite and real (p.4) | identity | EXACT | |
| 19 | Y = ∞ with Z finite: H0 → 1 (p.4) | identity | EXACT | |
| 20 | Z = ∞: [Z⁴]C = −D(u−1)²(Y+1); [Z³]C at Y = −1 is −u(u−1); branch (9) with ℓ = −u/(D(u−1)); (10); the second branch has H0 → 1+u; together they exhaust the two Y-roots (p.4) | identity + argument | EXACT + PROOF-READ | C is quadratic in Y, so there are 2 roots over Z = ∞, both found |
| 21 | barriers on branch (9): u > 0 gives H6 + H9/u finite and real; u < 0 gives H4 + (1−u)H6 finite (the 1/τ coefficients (1−u)/ℓ and −1/ℓ cancel) (pp.4-5) | identity | EXACT | series to O(τ²) |
| 22 | harmonic completion: L_O ≥ d·Φ_{f70}(O) on Ω; liminf ≥ 0 at every puncture; boundary points (L_O = 0, an added tie, a collision) are nonnegative; minimum principle; connectedness (p.5) | analytic | PROOF-READ | the standard argument, complete given 1-21 and the citations |

## D. Falsification

- **Main statement (brief).** I sampled 120000 random parameter tuples (D,k,u,h) ~ 2·N(0,1) and a complex Z, solved C = 0 for Y, and built the configuration from the inverse. I kept the points with both signs of η where all 18 off-T gaps are ≥ 1e-9 and d ≥ 0: 269 strictly supported points. On 40 of them I minimized Φ over SO(4) by BFGS with 6 restarts. Worst min Φ / max gap = −3.8e-15. That is rounding: the minimum is 0, attained at the permutation matrices of the contacts. Seed 20261125.
- **Intermediate claims.** Every inequality-type intermediate claim is exact here: (1), the reality of the limits, and the barriers. I did not find a counterexample to any of them. Nothing in the document is stronger than an exact identity plus a sign argument.

## E. Verdict

**confirmed-modulo-citations.** Citations: the earlier cases T_15, T_25, T_28 (or T_39) and T_65 for the seventh ties; Appendix B (Theorem B.1, Lemma B.2) and the one-circle theorem for collisions; Lemma 4.2; Theorem 1.2.

Presentation remarks (not gaps):
- p.1: the PDF does not name the earlier six-tuples that close the seventh ties. I supplied them (table B).
- p.2: the "remaining raw contact equation" is E·(H10 − u) = 0, written with raw products, not H10 − u. Replayed in that form.
