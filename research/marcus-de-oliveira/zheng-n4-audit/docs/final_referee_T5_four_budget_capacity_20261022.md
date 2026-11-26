Evidence level: EXACT (the polynomial certificate (4) mod N, the PSD criterion, the scalar-form algebra, the quaternion responses of f_0, f_6, f_10 and f_54, the T_5 generator list and profile identity, the gap form of the capacity); PROOF-READ (the Cauchy–Schwarz infimum, exhaustiveness of d = 0 / d > 0, the n = 0 case); NUMERIC (20000 random unit pairs).

# Audit: final_referee_T5_four_budget_capacity_20261022

PDF: `proofs/manin_transfer_20261021/final_referee_T5_four_budget_capacity_20261022.pdf` (3 pages). Batch F4-common-lemmas-b.
Replay: `replay/final_referee_T5_four_budget_capacity_20261022__check.py`, with output in `..._check_output.txt`. All 26 EXACT checks pass, and 2 NUMERIC checks pass.

## Deviations and limits

- The document does not fix a quaternion convention for (p, q) ↦ O. I searched 16 conventions (x ↦ p^{(*)} x q^{(*)}, the factors swapped, the matrix possibly transposed). Four of them reproduce all four responses: x ↦ p x q with the matrix [image of e_j as column j], and three equivalent variants. The symbolic identities were then proved for the first. Every O ∈ SO(4) is of the form x ↦ p x q / (|p||q|) (standard: SO(4) = S^3 × S^3 / ±1, and q ↦ q̄ is a bijection of S^3), so "for every O" is covered.
- The checker claims ("the Gram reduces to Z uniquely; rank ten on ten symmetric entries") and the claim that the new capacity covers an exact source in `final_referee_T5_old_SOC_obstruction_20261022.json` lying outside 144 SOC families are informational. Neither is used by Theorem 1. Both are UNCHECKED: the JSON was not loaded and the supplement's scripts were not run.
- An expectation in my first draft of the check was wrong, and I corrected it: generators(T_5) has 11 elements, not 8. See C.9. The document never claims 8 generators; it speaks of an 8-dimensional source span, which is correct.

## A. Statements

Notation: p, q ∈ R^4 with |p| = |q|; N = |p|^2 − |q|^2, n = |p|^2 + |q|^2. A, B, D, E, G and H = 2(A^2 + B^2) − EG are the displayed quadratics.

**Theorem 1 (p. 1).** For a, b, c, d ≥ 0 satisfying the capacity condition (1) ab(c − d) ≥ d^2(a + b + c + d), we have (2) a(A−B)^2 + b(A+B)^2 + cD^2 + dH ≥ 0 for all p, q with |p| = |q|. Equivalently (3): EG ≤ (|D| + 2 max(|A|, |B|))^2 for all such p, q. Equivalently: (3) implies (2) for every tuple satisfying (1).

**T_5 application (p. 2).** For x ↦ p x q (unit p, q), Φ_{f_0} = (A−B)^2, Φ_{f_6} = (A+B)^2, Φ_{f_10} = D^2 and Φ_{f_54} = H. Here f_0 = n_{1,4}, f_6 = n_{4,1}, f_10 = n_{13,24} and f_54 = h_{4231,−} in the paper's Appendix D numbering. For every c with (Pc)_π = 0 on T_5 = {1, 2, 3, 4, 6, 8} (the paper's T_5, 0-based lexicographic indices) and gaps g = Pc:

g = g22 f0 + g13 f2 + g12 f3 + g10 f5 + g17 f6 + g7 f10 + g16 f53 − A0 f54, where A0 = g16 − g0.

If A0 < 0 and g22 g17 (g7 + A0) ≥ A0^2 (g22 + g17 + g7 − A0), then Φ_c(O) ≥ 0 for all O ∈ SO(4). The document does NOT claim the T_5 case: its scope paragraph says that a separate argument is needed for all spectrally realizable configurations.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| Appendix D vectors f_j, generators(T_5), P, the Leibniz sum | paper / moc4lib | EXACT (moc4lib) |
| SO(4) = {x ↦ p x q : \|p\| = \|q\| = 1} | classical | PROOF-READ |
| Φ_{minor}(O) = det^2 ≥ 0 (paper (4.2)) | paper | EXACT (moc4lib) |
| `final_referee_T5_scalar_master_exact_20261022.py`, `final_referee_T5_four_budget_capacity_20261022.py` | supplement scripts | not used (by rule) |

No contact case, Appendix B or Corollary B.3 is invoked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | identity (4) modulo N | identity | EXACT | sympy: the difference reduces to 0 modulo the principal ideal (N), by division by the single generator N in lex order. The difference is not 0 before reduction. |
| 2 | (1) ⟺ Z ⪰ 0 for a, b, c, d ≥ 0 | inequality | EXACT + PROOF-READ | The Schur complement of d equals the displayed 3×3. The determinant of the 2×2 Schur complement equals 4[ab(c−d) − d^2(a+b+c+d)]/(c−d). For d = 0 the eigenvalues are 2a, 2b, c, 0. For d > 0, (1) forces c > d (and c = d makes the matrix not PSD). Also 400 random exact rational tuples: PSD by all principal minors ⟺ (1), 400/400. |
| 3 | (4) and Z ⪰ 0 give (2) | inequality | PROOF-READ | On N = 0 and n > 0 the right side is a sum of PSD forms. For n = 0, p = q = 0. |
| 4 | the choice d = 1, a = 1/α − 1, b = 1/β − 1, c = 1/γ makes (1) an equality on α+β+γ = 1 | identity | EXACT | |
| 5 | inf over the open simplex of u/α + v/β + w/γ = (√u + √v + √w)^2, zeros included | inequality | PROOF-READ | Cauchy–Schwarz. For zero entries, let the corresponding α tend to 0. |
| 6 | (2) gives (3) | algebra | EXACT (u + v = 2(A^2+B^2)) + PROOF-READ | |
| 7 | (3) gives (2) for all (1)-tuples | inequality | EXACT + PROOF-READ | (1) at d = 1 is the same polynomial as (a+1)(b+1)c[1 − 1/(a+1) − 1/(b+1) − 1/c]. Then Cauchy–Schwarz. The case d = 0 is trivial, and d > 0 is normalized by homogeneity. |
| 8 | quaternion responses of f_0, f_6, f_10, f_54 | identity | EXACT | 4 Φ_{f_j}(M(p,q)) − n^2 Q_j ∈ (N) symbolically for all four, and M M^T = (n/2)^2 I mod N |
| 9 | generators(T_5) and the 8-dimensional span | identity | EXACT | generators(T_5) = {f0, f1, f2, f3, f4, f5, f6, f9, f10, f53, f54}. All are minor except f54. The 8 profile vectors are independent, and all 11 span the 8-dimensional contact subspace (the 6 rows of P on T_5 are independent). |
| 10 | profile identity | identity | EXACT | Checked on all 11 generators, hence on the whole contact subspace |
| 11 | capacity in gap form, with a = g22, b = g17, c = g7, d = −A0 | identity | EXACT | Polynomial identity |
| 12 | the remaining summands are nonnegative minor costs | inequality | PROOF-READ | Coefficients g13, g12, g10 and g16 are gaps at noncontacts (≥ 0), and f2, f3, f5, f53 are minor |
| 13 | Gram uniqueness, rank 10; the SOC-obstruction source | informational | UNCHECKED | Not used |

## D. Falsification

- NUMERIC: 20000 random unit pairs (seed 2026). The minimum slack in (3) was 2.0e-3. The minimum of (2) with coefficients on the capacity boundary (d = 1, a = x, b = y, c = (x+1)(y+1)/(xy−1)) was 8.1e-3. No violation was found.
- Exactness makes Theorem 1 a proved statement. Falsification is relevant only to the application, which the document does not claim.

## E. Verdict

**confirmed.** Theorem 1, in all three equivalent forms, is proved by the exact certificate (4), which I replayed with my own code from the definitions in the PDF, together with the elementary PSD and Cauchy–Schwarz arguments. The T_5 identifications (quaternion responses, generators, profile identity, gap form of the capacity) are EXACT.

Scope remark, not a gap: the document proves safety of a T_5 cost only when A0 = g16 − g0 < 0 and the capacity inequality holds in gap variables. The regions A0 ≥ 0, and A0 < 0 with the capacity inequality failing, are not treated here. The document states this explicitly. Any T_5 proof (or any of the many cases listed as citing this document) must cover them by another argument, and its auditor must check that.
