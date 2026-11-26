Evidence level: EXACT for the radial octic, both chart substitutions, the -2 coefficient structure, the decompositions G~ = R + (V^2-1)^2 W with R >= 0 and the strict markers, and the chart coverage identities; PROOF-READ for the case logic; CITED for Lemma 1 (foundations/so4-atoms-part2b.md G4) and the unitary transfer; NUMERIC check of Lemma 1 on its boundary.

# Audit of `independent_moc4_T47_IIIb_plus_large_pc_20261022`

Auditor: C22-T47. PDF: `proofs/manin_transfer_20261021/independent_moc4_T47_IIIb_plus_large_pc_20261022.pdf` (2 pages). Replay: `../replay/independent_moc4_T47_IIIb_plus_large_pc_20261022__check.py` (output `../replay/independent_moc4_T47_IIIb_plus_large_pc_20261022__check_output.txt`).

**Verdict: confirmed-modulo-citations.** Lemma 1 (Gram capacity pc >= a^2, j >= d^2/(a+c)) and Theorem 2 (support inequality on the high-j large-pc plus region of T_65). Both chart certificates replayed EXACT with identical term counts. Lemma 1 is CITED to so4-atoms-part2b G4, where it is proved.

## Deviations and limits

- Nothing from the supplement was executed or imported. The radial octic G, the principal minor P and every certificate polynomial were rebuilt from the paper's definitions (f_j of Appendix D via lib moc4lib, the cubic C of (A.2), the height chart w(g,h) = c(g) + h kappa) by the replay script.
- The source vector is c(g) = sum of the eight coefficients times the paper's f_j (0-based Appendix D indices, as the document prints them); d is set to 1 only where the document normalises d = 1.
- Chart polynomials are cleared by multiplying with prod_i den_i^(deg_i G) (each coordinate's denominator to the coordinate degree of G). Where my term counts equal the document's, the cleared polynomial is very likely the same; where they differ, my polynomial is the document's times extra positive denominators, which is enough for the sign argument.
- The sign of G is fixed by the square root; it agrees with the documents' convention, since the documents' positivity tables (low-j, large-pc, stationary) come out positive with my G, with identical term counts. The signs of H6 and H8 are then fixed by the displayed identities (4LG + H6 Phi = K5^2 and the S+- forms), so the hypotheses H6 < 0, H8 < 0 have a well-defined meaning.
- I did not re-derive the Gram of Lemma 1. so4-atoms-part2b.md G4 verifies it for exactly this face and these roles: q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7.

## A. Statements

- **Lemma 1** (p. 1; contact set [[0, 4, 5, 6, 9, 12]]). For the T_65 face cost c = q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7: if a, p, c, d > 0, pc >= a^2 and j >= d^2/(a+c), then Phi_c(O) >= 0 for all O in SO(4), for every q, e, b >= 0. This is a pointwise inequality and does not use realizability.
- **Theorem 2** (p. 1; contact set [[0, 4, 5, 6, 9, 12]]). d = 1, b > 1, 0 < a < 1/b, p > 0, q >= 0, j > 1, j <= e, j > eb/(e+b), j < 1 + 1/p, pc >= a^2. Every x in X with c(x) of this form has Phi_{c(x)} >= 0 on SO(4). Proof: Lemma 1 if j >= 1/(a+c); otherwise G > 0, which contradicts G = 0.

The contact set {0,4,5,6,9,12} is literally paper T_65 (lib CONTACT_SETS[65]; this is historical T47). Its cone N_T is simplicial and its generators are exactly the eight f_j used (lib `generators`). The statements concern realizable costs, that is, x in X with any spectra. No genericity of the spectra is used: radial necessity and Prop. C.6 hold for every x in X.

## B. Dependencies

- [common-lemma] radial necessity (paper p.27): for x = u + iv in X, C(u) = 0 and det Hess C(u) = 0: G = 0 at every spectrally realizable point, via det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (replayed EXACT) and C < 0 (EXACT), so Q != 0 and h = -C/Q. Match with the use: checked.
- [core-doc] Theorem 1.2 / Corollary 2.2 (audited in foundations): transfer of the orthogonal support inequality to unitary matrices. Not audited here (CITED).
- [common-lemma] foundations/so4-atoms-part2b.md G4 (IIIb_plus_large_pc Lemma 1: PROVED IN THE DOCUMENT AND CHECKED): Lemma 1. Match with the use: checked.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Lemma 1 Gram construction | citation | CITED | so4-atoms-part2b G4 (same face, same roles) |
| 2 | split: j >= 1/(a+c) uses Lemma 1; otherwise 0 < c < 1/j - a, a < 1/j | case split | PROOF-READ | exhaustive |
| 3 | C < 0 so Q != 0 and G = 0 is necessary | identity | EXACT | C and det identity replayed |
| 4 | chart 1 (1 < j < b) covers the hard region: j = 1+T, b = j+A, a = V/(b(1+V)), p = theta(P)/T, e = j + j^2 theta(U)/A, c = (1/j - a) theta(S) | identity | EXACT | sympy limits and identities; b(1+V) - jV = b + AV > 0 |
| 5 | chart 2 (j >= b): b = 1+A, j = b+T, a = theta(V)/j, p = theta(P)/(A+T), e = j+U, c = S/(j(1+V)(1+S)) | identity | EXACT | sympy |
| 6 | chart 1: cleared G~ has 45504 terms, exactly 10 negative coefficients, all -2, with V exponent >= 2; R = G~ - (V^2-1)^2 W has 45474 terms, all positive; a positive monomial free of U, Q | inequality | EXACT | flint, coefficient by coefficient |
| 7 | chart 2: 31086 terms, 12 negative (-2), R has 31050 positive terms; a positive monomial free of T, U, Q | inequality | EXACT | flint |
| 8 | G~ > 0 on both closed chart domains (walls e = j, q = 0, j = b included), so G > 0, contradicting G = 0 | inequality | PROOF-READ | positive denominators; strict markers |

## D. Falsification

Lemma 1 was tested on the tightest boundary of its region: no negative minimum. A main-statement search was not run, since it is a special case of the conjecture.
- NUMERIC: Lemma 1 on its boundary pc = a^2, j = d^2/(a+c), q = e = b = 0: 12 samples (seed 20261125, 8 BFGS restarts each), min Phi over SO(4) >= 2.9e-14 (no violation)

## E. Verdict and gaps

Verdict: **confirmed-modulo-citations**.

Gaps: none found.

Exact replays: G (965 terms) and det identity; chart 1 certificate (45504 / 10 x -2 / R 45474 >= 0 / markers); chart 2 certificate (31086 / 12 x -2 / R 31050 >= 0 / markers); chart coverage identities.

Unchecked: the high-j small-pc subregion, which the document says it does not prove.
