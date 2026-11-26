Evidence level: EXACT for the profile identities, Q, the positivity of (e+j)C - ejQ, the discriminant identity (3), the conic structure and determinant (4), the identity modulo W, gcd(G,G_j) = gcd(G,P) = 1; EXACT-SAMPLED (exact rational evaluation at 12 random points, a probabilistic identity test) for (2); NUMERIC for the rank bound, the signed-class consequence and the cited Hessian sign theorem; GAP for the signed-class congruences, whose certificate is absent.

# Audit: blind_moc4_T48_negative_Phi_uniform_20261023 (3 pp.)

"Excluding a negative discriminant factor on a six-contact cone". Batch C27-T49. A branch used for paper case T_49.

## Deviations and limits

- No deviation from DOC_AUDIT.md. The degree-nine polynomial G (1772 terms) was loaded as data from the supplement file `blind_moc4_T48_native_full_radial_20261023.json` (json.load only). Every property of G used below was recomputed by my own code. The library cubic and the companion `C_terms` agree exactly, so the normalizations are the same.
- The supplement contains no data or verifier for the polynomials R7, R~7, L+, L- of the signed-class congruences, nor for the three 13-minor quotients (872, 421, 613 terms). I searched the package for a companion file. These claims could only be tested numerically.
- I did not compute the full 2784-term polynomial P. gcd(G,P) = 1 was proved by restriction to a rational line (see step 8). The term count was not checked.
- (2) was checked by exact evaluation at 12 random positive rational points, not as a full polynomial identity.

## A. Statements

Contact set T = {0,3,5,6,9,12} = paper T_49 (EXACT). Rays in the fixed order (q,p,j,c,e,d,b,a) are f7, f5, f29, f39, f2, f40, f41, f0.

**Theorem 1 (p. 1).** Let x in X be a realizable source whose gap vector has the canonical form (1), g = q f7 + p f5 + j f29 + c f39 + e f2 + d f40 + b f41 + a f0, with q,p,j,e,b,a,d > 0 and c >= 0. Put V = a+d+e, T = qV - (a+d)(d+e) and Phi = {(q+d)Vb + pT}^2 + 4dp^2(a+d)(q-e)V. Then Phi >= 0. All complex spectra are allowed. The coefficient sign conditions exclude every extra contact (the zero cases are other contacts, see the companion harmonic note). The note proves nothing about the support inequality: it only shows that the region Phi < 0 contains no realizable source. The Phi >= 0 sector is explicitly left open (p. 3, "Remaining scope").

Intermediate claims: Q = (c+d)(d+e+j) - p(e+j), with C(m + h kappa) = C + hQ; det H(w) = -4G^2/Q^2 at h = -C/Q, eq. (2); (e+j)C - ejQ has 42 positive coefficients; rank H <= 12 on {G = 0, Q != 0}; gcd(G,P) = gcd(G,G_j) = 1; (3) B^2 - 4AC0 = (b+c+d)^2 (q+d)^2 Phi H6; P S+- = L+-^2 mod G; (4) det M = -2(e-b)pW(q+d)^2 Phi^2; e^2(a+d)^2 Phi = b^2V^2{...}^2 mod W.

## B. Dependencies

| kind | ref | used for | status |
|---|---|---|---|
| common-lemma | "stated global spectrally realizable Hessian theorem": every principal 12-minor of Hess C is >= 0 at a realizable source | the only sign input (p. 1) | CITED, source document not named in the note. NUMERIC: 300 random complex spectra and phases (seed 99), all 91 principal 12-minors of Hess C(Re x) >= 0 (worst relative value 0) |
| core-doc | radial necessity (paper p. 27): C(u) = 0, Hess C(u) v = 0, det Hess C(u) = 0 | forces h = -C/Q and G = 0 | established (brief, item 5). Use matches: w = gap vector + h kappa is the real part of a matching vector |
| other-doc | companion verifier for the congruences and the 13-minor quotients | steps 6, 9 | not present in the supplement |
| paper-case | T_47 with roles (27,1,59,5,28,6) | remark only (p. 3) | EXACT: exactly one g = (4213, 3421, eps 1) maps them in order to (7,5,29,39,40,41). Not used in Theorem 1 |

## C. Steps

| # | claim (location) | kind | status | check |
|---|---|---|---|---|
| 1 | profile identities, form (1) (p.1) | identity | EXACT | 24 values; rank 8 |
| 2 | C(m + h kappa) = C + hQ, Q formula | identity | EXACT | own cubic (A.2); equals companion C_terms |
| 3 | (2) det H(w) = -4G^2/Q^2 | identity | EXACT-SAMPLED | 12 exact rational points |
| 4 | (e+j)C - ejQ has 42 positive coefficients, so C = Q = 0 is impossible and h = -C/Q is forced | identity + argument | EXACT + PROOF-READ | if Q = 0 then C = 0 on the radial line, so both vanish |
| 5 | P = (0,5)-deleted principal block is height-free | identity | EXACT | H_kappa lives on (0,6), (1,5) |
| 6 | three 13-minors divisible by G; rank <= 12 on G = 0, Q != 0 | identity + closure | NUMERIC (exact claim not replayed; data absent) | 4604 real points of G = 0: numerical rank <= 12 at all, exactly 12 at 4599 |
| 7 | gcd(G, G_j) = 1 | identity | EXACT | flint gcd |
| 8 | gcd(G, P) = 1 | identity | EXACT | restriction to a rational line with G_top(v) != 0 gives coprime univariate polynomials. Any common factor F would satisfy F_top(v) != 0 (F_top divides G_top), so it would restrict to a nonconstant common factor |
| 9 | (3) discriminant identity, H6 has 159 terms and no j | identity | EXACT | flint exact division |
| 10 | real root with G_j != 0 gives disc > 0, hence Phi < 0 implies H6 < 0 | argument | PROOF-READ | |
| 11 | P S+- = L+-^2 mod G, R7, R~7 coprime to G, hence P < 0 and all 12-minors <= 0 on the dense generic part | identity | **GAP** (certificate absent) | NUMERIC consequence: at 568 real points of G = 0 with Phi < 0 no principal 12-minor is positive |
| 12 | rank factorization: every 12-minor is P times a square | argument | PROOF-READ | H = C^T H_I C |
| 13 | conic form of G, (4), identity mod W, W != 0 when Phi < 0 | identity | EXACT | flint and sympy |
| 14 | smoothness of G = 0 at Phi < 0, including e = b and c = 0 | argument | PROOF-READ | Euler relation for e != b; for e = b a simple zero of det M plus Jacobi's formula gives tr(adj M M_e) = kappa X^T M_e X != 0; (z,u) independent of e |
| 15 | density of generic points near a smooth real point, continuity | argument | PROOF-READ, conditional on step 11's coprimality claims | a smooth real point is Zariski dense in its component; needs every maintained factor coprime to G (G_j, P EXACT; H6, R7, R~7 not checked) |
| 16 | the cited Hessian theorem forces all 12-minors to vanish | citation | CITED | numeric test above |
| 17 | second derivative -8 dG(v)^2/Q^2 != 0 forces rank >= 12; a rank-12 symmetric matrix has a nonzero principal 12-minor; contradiction | argument | PROOF-READ | uses (2) along a real path |

Counts: 17 steps; 8 EXACT (1 sampled), 6 PROOF-READ, 1 NUMERIC, 1 CITED, 1 GAP.

## D. Falsification

- Main statement (NUMERIC, indirect): on the 568 real points of G = 0 with Phi < 0 found from 3000 random positive coefficient draws (seed 2026), the signed-class conclusion holds. The cited Hessian theorem was tested on 300 random realizable sources (seed 99) without a violation.
- The intermediate exact claims (3), (4) and the identity mod W were replayed exactly, so there was nothing to falsify numerically.

## E. Verdict

**gap-found** (one gap, substantive, believed true).

- Gap, p. 2 "Signed quadratic classes": the congruences P S+- = L+-^2 (mod G), the degree-seven factors R7, R~7 and L+- are neither printed nor present in the supplement, and neither is the coprimality of R7, R~7, H6 with G. The same is true of the three 13-minor quotients behind the rank-twelve chart. These finite polynomial facts carry the whole sign argument. Numerically their consequences hold: rank <= 12 at 4604 real points of G = 0, and every 12-minor <= 0 at 568 points with Phi < 0. Repair: publish the polynomials and an exact check; the claims are then routine to replay.
- The cited global Hessian sign theorem is not named. It has to be resolved to an audited document. Numeric evidence only here.
- Logic otherwise sound. The smoothness argument (step 14) and the rank argument (step 17) are correct as written.
