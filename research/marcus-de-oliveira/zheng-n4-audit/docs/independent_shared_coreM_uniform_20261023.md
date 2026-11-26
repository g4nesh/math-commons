# Audit: independent_shared_coreM_uniform_20261023 (seven-contact set M)

Evidence level: EXACT for every identity, the generator decomposition, the group relabellings, the Hessian identity (6), the face formula (11) and all 24 coefficient tables. PROOF-READ for the coefficient-bound allocations, the real-part bridge, C > 0 and Q > 0, the affine endpoint argument (23) and the chart coverage. CITED for the sharp pair-SOC, old-H and three-minor-chain atoms and for the paper's Theorem 1.2. NUMERIC only for falsification.

PDF: `proofs/manin_transfer_20261021/independent_shared_coreM_uniform_20261023.pdf`, 5 pages. Replay: `replay/independent_shared_coreM_uniform_20261023__check.py`. Its output is in `__check_output.txt`: 52 EXACT checks, all pass.

## Deviations and limits

- No deviation from the audit plan. Nothing from the supplement was executed or loaded. The degree-8 polynomial G (273 terms) was re-derived from the paper's cubic (A.2) by my own fraction-free determinant, without the supplement's radial JSON. All 24 tables were rebuilt from that polynomial.
- I did not re-prove the three pointwise atoms. Their group relabellings are checked EXACT, and they were tested numerically.
- Time: about 25 minutes, shared with the L document.

## A. Statements

The contact set is M = {0, 1, 2, 7, 10, 16, 17}, which is the paper's Table 5 set #13 (EXACT). M contains G-copies of T_55 (remove 2 or 17), T_56 (remove 0 or 16), T_63 (remove 1 or 10) and T_80 (remove 7). It is T_55 ∪ {21}, T_63 ∪ {17}, and T_80 ∪ {17} or T_80 ∪ {22}, up to G. The supplement's catalogue lists this document as used for T_72 and T_85, but M is not G-equivalent to T_72 ∪ {j} or T_85 ∪ {j} for any j. Those uses must therefore be indirect. This is information for the T_72 and T_85 auditors, not a defect here.

**Theorem 1 (p. 1).** Let x ∈ X and pi* ∈ M, and put c = c(x). Suppose P c >= 0 and (P c)_s = 0 for all s ∈ M; extra contacts are allowed. Then Phi_c(O) >= 0 for all O ∈ SO(4), and the unitary statement follows by the paper's Theorem 1.2. All spectra are covered: the proof uses only x ∈ X ⊂ Sing C, and every zero coefficient is handled.

Internal claims:
- (1) The basis of the contact kernel is f0, f2, f3, f12, f15, f23, f26 with coefficients q, p, t, b, c, v, d. Here (q, p, t, b, c, v) = (g23, g15, g14, g8, g4, g9) and d = g6 - g8.
- (2) g13 = p + d, g12 = t + d, g6 = b + d, g3 = v + d, and f8 + f26 = f2 + f3 + f12 + f23.
- (2'), (3), (4) Each of the following conditions alone suffices: p(b+c) >= d^2, b(p+c) >= d^2, c(p+b) >= d^2, S(c-d) >= d^2 with S = t + v, or (when q, t > 0) c >= d^2/q + d^2/t.
- (5) F(w(z)) = C + Q z, together with Q > 0 and C > 0 (7) at any hypothetical negative point.
- (6) Q^4 det H(w(-C/Q)) = -4 (d+v)^2 Q^2 G^2. G is linear in v, G = B + vA, with B of 137 terms and A of 136 terms (9). G is symmetric under p <-> b (10). The c = 0 face satisfies (11).
- (15)-(24) G > 0 on the remaining region, shown by 24 cleared tables with positive coefficients.

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| sharp pair-SOC (f1, f9; f22) | common lemma | CITED | EXACT: G maps it to (f2, f12; f26), (f2, f15; f26) and (f12, f15; f26). Phi is G-invariant (lib (4.7)). NUMERIC: no violation |
| old-H (f3, f15; f26) and its image (f23, f15; f26) | common lemma | CITED | the image is EXACT. NUMERIC: no violation |
| three-minor chain (centre f2; leaves f0, f12; f70) | common lemma | CITED | EXACT: G maps it to (centre f15; leaves f0, f3; f26), 2 elements. NUMERIC at c = 1/q + 1/t over a 3x3 grid: min Phi >= -1e-15 |
| paper (4.2), minor vectors | paper | EXACT in lib | d <= 0 side and minor remainders |
| paper App. A.1, X ⊂ Sing C | paper | EXACT in lib | the bridge |
| paper Theorem 1.2 | paper | CITED | unitary sentence only |

The proof does not use Appendix B, Corollary B.3, collision compensation or any other contact case.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | generators(M) = [0,2,3,8,12,15,23,26] (the document's "zero-active list"); rank 7; basis (1); the coefficients are the stated gaps | identity | EXACT | lib + symbolic |
| 2 | (2) gap formulas, the relation, and the d <= 0 minor rewrite | identity | EXACT | symbolic |
| 3 | SOC pair atoms, and (2'): p(b+c) >= d^2 etc. suffice | citation + inequality | CITED (relabel EXACT), allocation PROOF-READ | split p in the ratio b : c. The deltas are b·sqrt(p/(b+c)) and c·sqrt(p/(b+c)); their sum is sqrt(p(b+c)) |
| 4 | old-H split gives (3) S(c-d) >= d^2 | inequality | PROOF-READ (atom CITED) | the same algebra as in L |
| 5 | the chain gives (4) | citation | CITED (relabel EXACT) | |
| 6 | (5): C and Q | identity | EXACT | paper cubic, exact match |
| 7 | Q = pb + c(p+b) + [d^2 - S(c-d)] > 0 and the SC identity (7); C > 0 | identity + inequality | EXACT + PROOF-READ | |
| 8 | bridge: Re x satisfies F = 0 and det H = 0 | analytic | PROOF-READ (argument checked by auditor) | same as L step 9 |
| 9 | (6) Hessian identity; G has 273 terms, is homogeneous of degree 8 with G(0,...,0,1) = 1, and is linear in v with B and A of 137 and 136 terms; (10) p <-> b symmetry | identity | EXACT | own Bareiss determinant (819 divisions), exact division and square root |
| 10 | (11): G restricted to c = 0 equals d(d+v)(d+b)(d+p)[d(q+d) - pb](q+t+d)^2, so G > 0 when pb < d^2 | identity | EXACT | homogeneous check |
| 11 | when c > 0 and d = 1: (7) (t+v)(c-1) < 1; when ct > 1: (8) q < t/(ct-1) | inequality | PROOF-READ | |
| 12 | sector c < 1: chart (15) covers p, b < 1; chart (16) covers p >= 1 with p(b+c) < 1 (the case b >= 1 is reduced by the symmetry (10)); chain rows (17) cover ct < 1, ct > 1, ct = 1 | case split | PROOF-READ | inverses rederived. In (17), the ct > 1 q-range is exactly [0, t/(ct-1)). Two mains >= 1 would satisfy (2') |
| 13 | sector c > 1: chart (19) covers c(p+b) < 1; (20) gives v < V = 1/C - t with V > 0; (21) covers 1/c < t < 1/C with the q bound (8); G(v) = (1 - v/V)B + (v/V)G(V) | case split + inequality | PROOF-READ | the affine interpolation is exact because G is linear in v; V's numerator and denominator have positive coefficients (EXACT) |
| 14 | wall c = 1: p, b < 1; rows t < 1, t > 1, t = 1 | case split | PROOF-READ | |
| 15 | 24 cleared tables (B, A, G(V)) have all coefficients > 0, with the markers of (18), (22) and (24) and 31040 coefficients in total | finite certificate | EXACT | every count and marker coefficient is reproduced |
| 16 | conclusion, then the unitary form | logic + citation | PROOF-READ / CITED | |

## D. Falsification (NUMERIC)

- Chain atom at c = 1/q + 1/t for q, t ∈ {0.2, 1, 5}, 12 restarts each, seed 2: min Phi = 9e-16. The SOC and old-H atoms are as in the L report: no violation.
- Realizable points with contact set exactly M (margin 0.01, seeds 0 to 2): (d + v)·G / scale^9 is at most 5e-26 in absolute value, consistent with the necessary equation. min Phi over SO(4) >= -1e-13. Seed 0 has d > 0, and the coefficient bounds hold there.

## E. Verdict

**confirmed-modulo-citations.** Every algebraic step and every finite certificate was replayed exactly from the paper's definitions, and every table matches the document's counts and markers. The chart coverage was rederived and is complete: it includes the walls c = 0, c = 1, ct = 1, the zero leaves and q = 0. The bridge argument is sound and needs only X ⊂ Sing C. Citations: the SOC atom (f1, f9; f22), the old-H atom (f3, f15; f26) and the three-minor chain (f2; f0, f12; f70), all from the supplement, and paper Theorem 1.2.

Gaps: none. Presentation remarks:
- The equation numbering is duplicated or out of order: (1) and (2) are used twice, and (4) appears before (3).
- The text says "(2)" for the triple of pair conditions, but also uses (2) for the gap formulas.
