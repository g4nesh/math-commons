Evidence level: EXACT for the a = 0 source, zero faces, height formulas, the radial identity with G8, the Lemma 1 algebra and all seven sign tables (own code, 23 checks, 0 failures). PROOF-READ for the case split. CITED for the singular-leaf Gram construction (audit foundations so4-atoms-part2b G3), radial necessity and the paper cases T_29, T_36, T_1, T_39, T_7. NUMERIC as support only.

# final_referee_T46_a0_all_small_product_uniform_20261023 ("A support inequality in the small-product region of a seven-contact T46 cone", 3 pp.)

This is a branch document, used for paper case T_47 through the seven-contact set T^(7)_26 = {1, 3, 5, 6, 7, 8, 14}.

Replay: `replay/final_referee_T46_a0_all_small_product_uniform_20261023__check.py`, output in `__check_output.txt`.

## Deviations and limits

1. G8 was rebuilt by my own code from (A.2) with the full Hessian. The document's 380-term G8 and degrees (2,2,2,2,2,8,3) agree. My term counts for all six tables also equal the document's: 815, 60, 210747, 196734, 3897, 4650.
2. The PSD construction behind Lemma 1 (singular leaf beta = -r) is not re-proved here. The foundations prove it, together with the rational section p >= k + D0^2/(jN) (part 2b, G3). I re-checked only the scalar algebra of the section, exactly.
3. The second-row strict monomial, which the document leaves to its JSON, is my own find: -A^25 V^24 J^7 C1^5.

## A. Statements

- **Lemma 1 (p. 1).** On the 8-ray T46 face, take d = 1, q, b, c, j > 0 and q != b. Put r = sqrt(qb). If p >= P* = k + D0^2/(jN), the cost q f27 + p f1 + j f59 + c f5 + e f2 + f28 + b f6 + a f4 has Phi >= 0 on SO(4), for any nonnegative additions on the minor rays. This is a universal SO(4) inequality: it does not need the cost to be realizable.
- **Theorem 2 (p. 2).** For x in X with c(x) on the a = 0 face (contact set containing T ∪ {14}) and qb <= d^2: Phi_{c(x)} >= 0 on SO(4). All zero coefficients and the wall qb = d^2 are included. Spectra are unrestricted (coefficient-space argument).

## B. Dependencies

| dependency | kind | status | check |
|---|---|---|---|
| singular-leaf Gram (Lemma 1's input "general six-ray Gram theorem") | common lemma | CITED, foundations G3 "singular leaf (q != b)" PROVED | p*, c* formulas and the section identity match (EXACT here too) |
| q, p, j, c, b = 0 lead to historical T28, T35, T0, T38, T6 = paper T_29, T_36, T_1, T_39, T_7 | paper cases | CITED (Stage 1) | EXACT containment on the a = 0 face |
| radial necessity | common lemma | CITED | - |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | Lemma 1: j - kappa t^2 = jN/(D0+A0)^2, consumed p = P*, c - consumed c > 0 via N - rho^2 A0/c = D0^2 + D0 A0 + A0^2, N >= D0^2 | identity | EXACT | sympy |
| 2 | Lemma 1: the PSD construction | inequality | CITED | foundations G3 |
| 3 | d = 0 is trivial; zero faces are cited | citation | EXACT + CITED | - |
| 4 | C0 < 0 and Q as printed, so Q != 0 and G8 = 0 is necessary | identity | EXACT | Q^3 det = -4 Q (q+d)^2 G8^2 |
| 5 | double-small table (homogenized in (A, Wq), (V, Wb)): 815 negative coefficients, strict monomial -32 A^2 V^3 | certificate | EXACT | covers 0 < q, b <= 1 |
| 6 | q = b = 1: 60 negative coefficients, constant -32 | certificate | EXACT | redundant with 5 |
| 7 | rows q > 1, qb < 1 / b > 1, qb < 1 / both qb = 1 walls, with p = P* S/(1+S): all coefficients negative, strict monomial free of S and E | certificate | EXACT | 210747 / 196734 / 3897 / 4650 terms; denominators certified positive |
| 8 | the regions exhaust qb <= 1: both <= 1, or exactly one > 1 (then q != b) | case split | PROOF-READ | - |

## D. Falsification

- NUMERIC: 12 samples of Lemma 1 exactly at p = P* (e = a = 0, seed 11, q, b, j, c in [0.05, 3], 20 BFGS restarts). Worst min Phi / max coef is 2.1e-16, so there is no violation. Control with p removed: -1.2e-2, so the test can detect failures.

## E. Verdict

**confirmed-modulo-citations.** Every table and identity holds exactly, and the split is exhaustive. The citations are the foundations singular-leaf Gram construction, radial necessity, and the paper cases T_29, T_36, T_1, T_39, T_7. No gaps.
