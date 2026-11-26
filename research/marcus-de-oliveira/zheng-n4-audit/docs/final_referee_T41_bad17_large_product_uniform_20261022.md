# Audit: final_referee_T41_bad17_large_product_uniform_20261022 (paper T_42, bad17 branch, region ab >= d^2)

Evidence level: EXACT for the restricted cubic, the compact Heron identities (Fe/Y^2 = f, Fq/X^2 = g, R0 = r XY, A0 polynomial, R0^2 - 4VA0 = Fe Fq), Lemma 1's algebra, all four coefficient tables of Lemma 2 (with matching term counts), the group transport of the T24 core and the T_42 stabilizer; EXACT pointwise (25 rational points) for the link between the document's G(c,j) and the full 14 x 14 Hessian; PROOF-READ for Proposition 1's deduction; CITED for the "T24 coherent capacity" and the "double-SOC Gram" (checked numerically at equality).

PDF: proofs/manin_transfer_20261021/final_referee_T41_bad17_large_product_uniform_20261022.pdf (3 pages). Auditor label C16-T40-T42.

## Deviations and limits

- The document does not display the identity linking G(c,j) to the Hessian ("its full fourteen-dimensional Hessian necessity is G(c,j) = 0"). I rebuilt G from the document's compact data alone: G = V D^2 + R0 j D + A0 j^2 with R0 = r(e/Y, q/X) XY and A0 = (R0^2 - Fe Fq)/(4V), and then compared with the Hessian of the paper's cubic (A.2) at 25 random rational points. Finding: -det H(w) Q^2/4 = (c+1)^2 G(c,j)^2 exactly at all 25 points (constant ratio 1). The extra factor (c+1)^2 is not mentioned; since c >= 0 it never vanishes, so "det H = 0 iff G = 0" holds on the region. This is a presentation point. The polynomial identity itself is EXACT only pointwise.
- The cited "T24 coherent capacity" (blind_moc4_T24_coherent_capacity note) and "double-SOC Gram" were not audited; I tested their use numerically.
- Paper case identification: the five T24 core rays f59, f53, f55, f67, f54 have common zero set {3,4,7,8,10,12,13,17,22}, which contains paper T_25 = {3,4,7,8,10,12}, so "T24" is the supplement's index for paper T_25 (stage 1, like T_42).

## A. Statements

Proposition 1 (p. 1). Let g = q f3 + p f1 + j f53 + c f6 + e f15 + d f17 + b f16 + a f4, all coefficients >= 0, d > 0, ab >= d^2. If g = c(x) for some x in X (any extra contacts, any spectra; the hypothesis is only on the vector), then Phi_g(O) >= 0 for all O in SO(4). More strongly, every real (c, j) with c >= 0, j > 0 satisfying the Hessian necessity G(c,j) = 0 obeys d^2/j <= p + d(ab - d^2)/((a+d)(b+d)) (eq. (1)).

Lemma 1. With D = XY(c+1) + Tcj: G(c,j) = V D^2 + R0 j D + A0 j^2; for a root with j > 0, c >= 0, J = jXY/D satisfies 0 < J <= j and G(0,J) = 0.

Lemma 2. For a, b, p >= 0, ab >= 1, x, y in [0,1]: with L = (ab-1)/((a+1)(b+1)) and h = 2V(p+L), h + r >= 0 and (h+r)^2 >= f g.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| X subset of {grad C = 0} (paper A.1), giving det H(w) = 0 for a realizable cost | core | CITED (library replays) |
| T24 coherent capacity: Phi >= 0 for the core d f17 + p f1 + j f53 + a f4 + b f16 when (1) holds and ab > d^2 (transported from [f59,f53,f55,f67,f54]) | other-doc | CITED; transport EXACT: exactly 2 elements of G send the ordered list [f59,f53,f55,f67,f54] to [f1,f53,f4,f16,f17], e.g. (rho, rho', eps) = ((0,3,2,1),(2,1,0,3),0) in the library's convention, which is the document's "row = (2,1,0,3), column = (0,3,2,1)" with the two roles named in the other order |
| double-SOC Gram: core safe when ab = d^2 and pj >= d^2 | other-doc | CITED |
| final_referee_T41_bad17_c_elimination / full_radial notes and scripts | supplement | not used; identity (2) rebuilt independently as described above |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | C = -(c+1){...} < 0, so Q != 0 for a realizable cost and h = -C/Q | identity | EXACT | sympy from (A.2) |
| 2 | Hessian necessity is G(c,j) = 0 | identity | EXACT at 25 points (global identity not proved here); extra factor (c+1)^2 noted | see Deviations |
| 3 | (2) form of G; R0 (48 terms) and A0 (101 terms) independent of j, c | identity | EXACT | built from the compact formulas; term counts match the document |
| 4 | D^2 G(0, jXY/D) = (XY)^2 G(c,j); D - XY = c(XY + Tj) >= 0; hence 0 < J <= j; at j = 0, G = V D^2 > 0 | identity + inequality | EXACT + PROOF-READ | |
| 5 | Heron coordinates: R0/XY = r, R0^2 - 4VA0 = Fe Fq, Fe/Y^2 = f, Fq/X^2 = g | identity | EXACT | |
| 6 | (4): t = 1/J solves V t^2 + r t + A0/(X^2Y^2) = 0, so t = (-r +- sqrt(fg))/(2V) with fg >= 0 | algebra | PROOF-READ | |
| 7 | Lemma 2 reduction a >= b by the swap (a,b)<->(b,a), x<->y | symmetry | EXACT | r is invariant and f <-> g under the swap of the compact formulas, so no spectral stabilizer is needed (the T_42 stabilizer does act this way: EXACT) |
| 8 | the two charts cover {a >= b, ab >= 1} | case split | PROOF-READ | b >= 1: b = 1+B, a = b + A; b <= 1: 1/b = 1+B, a = 1/b + A |
| 9 | Lemma 2 coefficient tables: H has 135 / 183 and E has 1454 / 1828 coefficients, all strictly positive | inequality certificate | EXACT | own expansion in (A,B,p,xi,eta); denominators (1+xi)(1+eta)(1+B)^k are positive; term counts equal the document's |
| 10 | 1/j <= 1/J <= (-r + sqrt(fg))/(2V) <= p + L, i.e. (1); p + L = 0 forces no root | deduction | PROOF-READ | |
| 11 | (1) and the cited capacity make the core safe; ab = d^2 by double-SOC Gram; q f3, e f15, c f6 are minors | citation | CITED | NUMERIC check at equality (D) |

## D. Falsification

- NUMERIC (seed 3, 400 samples with ab >= 1, about half of p, q, e, c set to 0): 279 positive roots j of the rebuilt G; the smallest slack p + L - 1/j was -3e-16 (rounding at equality), so (1) held everywhere.
- NUMERIC (seed 4, 12 samples, 12 BFGS restarts each): for p f1 + j f53 + a f4 + b f16 + f17 at equality in (1), the smallest Phi over SO(4) found was +3.7e-14. The cited capacity is therefore consistent with being sharp at (1), and no violation was found.

## E. Verdict

confirmed-modulo-citations. No gap. Citations: T24 (= paper T_25) coherent capacity and double-SOC Gram (both tested numerically at equality), and X subset of {grad C = 0}. Two presentation remarks: (i) the Hessian determinant equals -4 (c+1)^2 G^2/Q^2, not -4G^2/Q^2; harmless since c + 1 > 0; this was checked only at 25 exact points, not as a polynomial identity; (ii) the stabilizer argument for a >= b is not needed, because Lemma 2 is symmetric under the swap.

Replay: `replay/final_referee_T41_bad17_large_product_uniform_20261022__check.py`; output in the matching `__check_output.txt` (all EXACT checks pass).
