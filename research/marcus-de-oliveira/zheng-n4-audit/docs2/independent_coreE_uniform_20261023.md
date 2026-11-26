Evidence level: EXACT for the decomposition (1), the group transport of the cited inequalities, the algebra of (5), (6) and Lemma 1, the radial identities (8)-(10), the c = 0 factorization, all eleven positivity tables (15) and the inverse maps. PROOF-READ for the necessity argument, the case partitions and the assembly. CITED for the three reference Leibniz-sum inequalities behind (2), (3) and (4). NUMERIC support for the cited inequalities and for the main statement.

# Second-pass audit: independent_coreE_uniform_20261023

- Document: `proofs/manin_transfer_20261021/independent_coreE_uniform_20261023.pdf`, 5 pages, "A seven-contact inequality from a four-coefficient and chain allocation", dated 23 October 2026.
- Auditor label: K2-cores-EF-second. This is an independent second pass. The first auditor's report was not opened.
- Replay: `replay2/independent_coreE_uniform_20261023__check.py`, with output in `replay2/independent_coreE_uniform_20261023__check_output.txt`. It runs in 130 s, and all 43 EXACT checks pass.

## Deviations and limits

- All objects were rebuilt from the paper's definitions through the shared library `lib/moc4lib.py`: P, the f_j, the cubic C of (A.2), kappa = e_13 and the group G. No supplement script or data file was executed, imported or loaded. The document's "companion exact file" and its hashes were not used.
- The document says its 408-term polynomial G appears "in the checker output". I did not read that output. I derived G from scratch as the square root of -(Q^4 det H)/(4 d^2 Q^2), and its sign is fixed by the c = 0 factorization. Because all eleven term counts and all eleven lowest strict monomials of (15) match the document, my G is the same polynomial as the document's.
- The three reference inequalities (sharp SOC (f1,f9;f22), four-coefficient (f0,f6,f10;f54) and chain (f2;f0,f12;f70)) are proved in other supplement documents. I did not audit those documents. Here they are CITED. I checked exactly that G transports them onto the displayed images (2), (3), (4) and (4'), and I tested the images numerically. I attempted no proof of them. A quick probe showed that no fixed-sign identity of the form Phi_26 + 2s mu_i mu_j >= 0 holds.
- The final remark extending the result to unitary matrices rests on the paper's Theorem 1.2. That remark is not needed for the SO(4)/O(4) statement and was not audited.

## A. Statements

Notation: x in X, pi* in the contact set, c = c(x) = Re(x - F_pi*(x) kappa), and Phi_c(O) = sum_sigma sgn(sigma) (Pc)_sigma prod_i O_{i,sigma(i)}. Permutations are 0-based lexicographic, and f_j follows the numbering of Appendix D. The supplement uses the same numbering (docs/notation.md).

- **Theorem 1 (p. 1).** Let E = {0,1,2,9,10,16,17}. For every x in X with P c(x) >= 0 and (P c(x))_sigma = 0 for every sigma in E (non-strict: further contacts are allowed), Phi_{c(x)}(O) >= 0 for every O in O(4), equivalently SO(4). The statement covers all of X, so it includes repeated eigenvalues, a_i + b_j = 0 and all zero coefficients, because the proof uses only X ⊂ Sing C. E is exactly the paper's Table 5 seven-set no. 17 (EXACT). The document is not the primary proof of any paper case T_k.
- **Lemma 1 (p. 2).** Write g = q f0 + p f2 + t f3 + b f12 + c f15 + v f13 + d f26 with all coefficients >= 0. Suppose c, d > 0, v >= 0, q(c+v) > d^2 and t >= R_d(q,v,c) = d(v+d){q(c+v)+dv}/((c+v){q(c+v)-d^2}). Then Phi_g >= 0 on SO(4). The same holds with q and v interchanged. No realizability hypothesis is needed.
- **Intermediate claims.** (1): the unique decomposition of a supported cost. (5) and (6): sufficient coefficient conditions. (8): C(g + z kappa) = C_0 + Q z. (9): Q^4 det H(g + z kappa)|_{z=-C_0/Q} = -4 d^2 Q^2 G^2. (10): the symmetries of G. The c = 0 factorization. (15): eleven positivity tables.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| Sharp SOC: x f_i + y f_j + d f26 has Phi >= 0 when xy >= d^2, for {i,j} ⊂ {2,12,15}, as images of the reference (f1,f9;f22) | common lemma | CITED. Some g in G maps (f1,f9,f22) onto each of the three images (8 elements each, EXACT). NUMERIC: min_O [Phi_26 + 2 sqrt(Phi_i Phi_j)] = 0 up to 1e-29 |
| Four-coefficient (3): q f0 + v f13 + t f3 + d f26, with the reference (f0,f6,f10;f54) | common lemma | CITED. 4 elements of G map it onto the image, with the t-role f10 -> f3 (EXACT). NUMERIC: min >= 0 on the boundary t = (q+1)(v+1)/(qv-1) at 8 (q,v) pairs |
| Chain (4) and (4'): q f0 (or v f13) + t f3 + c f15 + d f26 with c >= d^2/q + d^2/t, from the reference (f2;f0,f12;f70) | common lemma | CITED. 2 elements of G map it onto each image, with c-role f2 -> f15 (EXACT). NUMERIC: min_O [Phi_26 + 2 sqrt(Phi_15)(sqrt(Phi_0 or 13) + sqrt(Phi_3))] >= 0 |
| X ⊂ Sing C, with C the cubic of (A.2) | core (paper App. A) | replayed EXACT by the shared library (grad C vanishes identically on matching vectors in 16 variables) |
| (4.2): Phi_{n_{I,J}}(O) = det O[I,J]^2 | core (paper Sect. 4.1) | EXACT on rational matrices plus a general Laplace/Jacobi proof (library README) |
| X is a cone stable under complex conjugation | paper Sect. 3 | PROOF-READ: m(conj a, conj b, conj alpha) = conj m |
| Theorem 1.2 (two-orthogonal reduction), used only for the unitary remark | paper | CITED and not needed |

The paper's Lemma 4.2 (cone generation) is not needed. The seven generators span the 7-dimensional space {c : (Pc)_E = 0}, and the coefficients are gap values, so they are nonnegative by hypothesis.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | generators(E) = {f0,f2,f3,f12,f13,f15,f26}. They are independent, rows [23,15,14,8,4,7,3] give I_7, and so q,p,t,b,c,v,d are the gaps at those rows and are >= 0 | identity | EXACT | library generators, rank, values |
| 2 | f0,f2,f3,f12,f15,f13 are minor vectors, so d = 0 gives Phi >= 0 | identity | EXACT | (4.2) via the library. F_IS_MINOR; f26 = h_{4132,-} |
| 3 | (2), (3), (4) and (4') hold | citation | CITED | group transport EXACT; images tested NUMERIC |
| 4 | (5) follows from the two chains by a proportional split; zero walls use one chain | identity | EXACT + PROOF-READ | symbolic identity |
| 5 | (6) follows from two SOC atoms with a shared central coefficient. s = b + c > 0 is forced when d > 0 | identity | EXACT + PROOF-READ | symbolic identity |
| 6 | Lemma 1: the allocation identities, positivity of X, W - X, q_C and t_F - delta, the choice delta = dv/(c+v) giving the total R_d, and the v = 0 wall equal to the chain | identity | EXACT + PROOF-READ | sympy: 7 identities |
| 7 | Realizable gives C(w) = 0 and det H(w) = 0 for w = Re m = g + Re(F_pi*) kappa | analytic | PROOF-READ | polarization with x = m/2, y = conj(m)/2 in Sing C. C(x+y) = 0, H(x+y)(x-y) = 0, and if x = y then H(w)w = 0 with w ≠ 0 |
| 8 | (8): C(g + z kappa) = C_0 + Q z, with C_0 and Q exactly as displayed (factor 1 against (A.2)), and C_0 > 0 when d > 0 | identity | EXACT | flint |
| 9 | (9): exact Bareiss on the 14 x 14 Hessian (all remainders zero). Q^4 det H(-C_0/Q) = -4 d^2 Q^2 G^2, where G has 408 terms and is homogeneous of degree 8. Hence d > 0 gives G = 0 | identity | EXACT | flint fmpq_mpoly, det with 7264 terms |
| 10 | (10): G is invariant under q <-> v and p <-> b | identity | EXACT | flint |
| 11 | c = 0: G = d A {(d+q)(d+v) - pb} B, plus the B(T_F) identity. The sign argument then uses pb < d^2 and t < T_F | identity + inequality | EXACT + PROOF-READ | flint and sympy; failure of (3) rewritten as (qv - d^2)(t - T_F) < 0 |
| 12 | The main partition into c >= 1 / row L / row P covers the whole failure set of (6) with c > 0 | case split | PROOF-READ + EXACT inverses | inverse maps checked symbolically |
| 13 | The auxiliary partition of (11) and (12) covers {t < R} or free t according to q(c+v) vs 1. The q <-> v orientation is legitimate because G and the failure set are symmetric | case split | PROOF-READ | Lemma 1 hypotheses hold in each row |
| 14 | The dominant-c partition by cw vs 1, with t < w/(cw-1) from the failure of (5). w = 0 and the zero walls are included | case split | PROOF-READ | |
| 15 | The eleven tables have all coefficients > 0, a strict monomial in the strictly positive parameters, and reduced denominators with positive coefficients and a strict monomial. The term counts 11763, 84446, 4181, 11409, 1091, 23837, 93285, 6427, 27531, 23799, 3955 and the markers 64C^2, 128C^2, C^3, C^2Q^2, 1, 2C^3, 2C^2Q^2, 2, 2, 16Q^2, 16 all match | inequality certificate | EXACT | flint, with m_i = deg_i G = (3,2,2,2,2,3) |
| 16 | Assembly: d = 0 is trivial; d > 0 gives G = 0; c = 0 is excluded; for c > 0 the failure point lies in a table domain where G > 0. Contradiction | logic | PROOF-READ | homogeneity of G and of all conditions justifies d = 1 |
| 17 | O(4) is equivalent to SO(4): Phi at O·diag(1,1,1,-1) equals -Phi(O) | identity | PROOF-READ | |
| 18 | Extension to unitary matrices | citation | CITED (Thm 1.2) | not needed |

## D. Falsification

- Cited inequalities, sharp forms: 30 multistart Nelder-Mead and BFGS runs over SO(4) for each SOC pair and each chain; 8 (q,v) pairs x 12 starts for (3) on its boundary. The worst value is -5.4e-30, which is roundoff. Separately, 50000 random SO(4) matrices (seed 20261125, scratch run) gave a worst value of +4.0e-4.
- Main statement: least squares in (a, b, alpha) found 20 realizable points with contact set exactly E and d >= 0.3 (seed 7, 50 starts). At these points the relative size of G is at most 3.9e-15, which independently confirms the necessity G = 0. The multistart minimum of Phi over SO(4), divided by the largest gap, is >= 5.5e-18. In a scratch run with 25 points, every point satisfied at least one of the sufficient conditions (2) to (7). Points with d = 0, found by the library sampler, are trivial.
- No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Every step is EXACT or PROOF-READ, except the reference inequalities (2), (3), (4) and (4'). Those are CITED to the supplement's SOC, four-coefficient and chain documents, which I did not audit; their group transport is EXACT and they are supported NUMERICALLY. The argument covers degenerate spectra and all zero walls, because it uses only the closed condition X ⊂ Sing C.

Presentation remarks, not gaps:
- Equation labels (10) and (12) are each used twice: the symmetries and the main-map table share (10), and the auxiliary table and the dominant-c map share (12). There is no (13).
- "The third factor" in the c = 0 paragraph means {(d+q)(d+v) - pb}.
