# Audit: final_referee_T2_scalar_compensation_20261022

Evidence level: EXACT (every algebraic identity, and the orbit claim) and PROOF-READ (inequality steps). NUMERIC for falsification.

PDF: `proofs/manin_transfer_20261021/final_referee_T2_scalar_compensation_20261022.pdf`, 2 pages. This is not the primary proof of a paper case. The supplement's catalogue lists it as used for paper T_3. Batch C01-T1-T3.

## Deviations and limits

- The note proves scalar inequalities only. It says so itself on p. 2: identifying S², T², H with Leibniz sums of f9, f1, f22 "requires that atom's separate exact coordinate bridge". That bridge was not looked for or checked here. So this note does **not** by itself establish the SOC support inequality that the T_3 note cites as (1). We tested that inequality only numerically (see the T_3 audit and section D).
- We did not check the phrase "earlier sharp family d − m4 + τ m6 + τ^{-1} m8" (p. 1), because its coordinates are not defined in the note. UNCHECKED.

## A. Statements

- **Lemma 1.** For r, s, u, v ≥ 0: √((r+u)(s+v)) − √(rs) ≥ √(uv).
- **Theorem 2.** Let (a,b,c,d) and (A,B,C,D) be nonnegative with equal sums. Put s = c + d − C − D and t = b + c − B − C. Then R + 2|st| ≥ 0, and hence τs² + R + τ⁻¹t² ≥ 0 for all τ > 0. Zero coordinates are included.
- **Theorem 3.** For p, q ∈ R⁴ with |p| = |q|: H + 2|ST| ≥ 0, and τS² + H + τ⁻¹T² ≥ 0 for all τ > 0. The page image confirms H = 2Σ_i p_i² q_{i+2 mod 4}² − 4p0p1p2p3 + 4q0q1q2q3.
- **Orbit claim (p. 1).** Facet triple [11, 39, 17] is mapped to [9, 1, 22] by row permutation (0,2,3,1) and column permutation (3,1,2,0), without transposition, and the orbit has 144 elements.

## B. Dependencies

None, apart from elementary algebra. The orbit claim uses the paper's facet numbering (identical to the supplement's) and the group G (lib).

## C. Steps

| # | claim | status | how |
|---|---|---|---|
| 1 | Lemma 1: the difference of squares is (√(rv) − √(su))² | EXACT | polynomial identity in square roots |
| 2 | R = 2δ² + 2Δ² − 2(v_a v_c + v_b v_d) | EXACT | |
| 3 | st = v_b v_d − v_a v_c, using Σv = 0 | EXACT | |
| 4 | (3): R + 2\|st\| = 2δ² + 2Δ² − 4 min(v_a v_c, v_b v_d) | PROOF-READ | −2(P+Q) + 2\|Q−P\| = −4 min(P,Q) |
| 5 | If both products are > 0, the pairs have opposite signs; WLOG v_b, v_d > 0 by exchanging the quadruples (R and \|st\| are invariant) | PROOF-READ | |
| 6 | Δ − δ ≥ √(αγ) + √(MN) by Lemma 1 twice; then 2δ² + 2Δ² ≥ (Δ−δ)² ≥ 4 min | PROOF-READ | the requirement Δ − δ ≥ 0 before squaring holds by the bound itself |
| 7 | τ-completion | EXACT | |
| 8 | Theorem 3: H ≥ R, because the signed products are ≥ −4\|·\| | EXACT + PROOF-READ | |
| 9 | Checkerboard identity F = 2δ² + 2Δ² + (v_a − v_c)² − 4 v_b v_d (equal masses) and the equality family | EXACT | the "exactly when" uniqueness was not checked (remark only) |
| 10 | Orbit claim | EXACT | The claim holds when the row permutation (0,2,3,1) acts through its inverse in lib's convention, σ ↦ ρ σ ρ'. The literal reading gives (8, 2, 69). Exactly 8 group elements fix the triple as a set with h ↦ h, so the orbit size is 1152/8 = 144. The T_3 triple (1, 11, 17) lies in the same orbit. |

## D. Falsification (NUMERIC)

- 50000 random pairs with equal mass, seed 5, about 25% zero coordinates: the minimum of (R + 2|st|)/mass² was 0 (it is attained at zero boundaries); the minimum of (H + 2|ST|)/|p|⁴ was 7.3e-4.
- The Leibniz version of the [9,1,22] family: the minimum over SO(4) × τ of Φ(τ f9 + f1/τ + f22) was 2.3e-17 (80 restarts). Tight and not violated.

## E. Verdict

**confirmed** for what the note states (Lemma 1, Theorems 2 and 3, the orbit claim). Every step is EXACT or PROOF-READ.

Scope caveat for its use in T_3: the bridge from (S², T², H) to the Leibniz sums of (f9, f1, f22) is explicitly not supplied. So this note does not prove the SOC input (1) of the T_3 note. That input must come from another document or be audited separately.

Presentation: the row-permutation convention is the inverse of the one used in our library.
