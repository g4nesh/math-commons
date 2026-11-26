# Audit: blind_shared_coreJ_uniform_20261023 ("A seven-contact cone and a connected radial branch")

Evidence level: every identity, generator decomposition, determinant, divisibility, gcd/irreducibility, sample value and parametrization in the document was replayed EXACT by our own code. The analytic and topological steps (rank bounds, constant inertia, connectedness) are PROOF-READ. The proof rests on one CITED theorem, the "global matching Hessian inertia theorem" (1). We have NUMERIC support for it but did not audit it.

PDF: `proofs/manin_transfer_20261021/blind_shared_coreJ_uniform_20261023.pdf`, 5 pages. Replay: `replay/blind_shared_coreJ_uniform_20261023__check.py`, output in `__check_output.txt`. On the last run 27 EXACT checks passed and none failed.

## Deviations and limits

- We did not audit the Hessian inertia theorem (1), the only geometric input. We tested it numerically: 200 random complex spectra with w = 2 Re x, where x is the matching vector. Every sample had numerical rank exactly 12, and the smallest normalized principal 12-minor was 1.1e-17, so no negative minor was found.
- The document's own checker and certificates were not run or read. Every polynomial (C, Q, A, B, G, the 14x14 determinant and the three borders) was rebuilt from the paper's cubic (A.2) and the paper's f_j. We used the shared library moc4lib, which was itself rebuilt from the paper.
- The unitary assertion ("simultaneous two-real-response theorem") is the paper's Theorem 1.2, CITED. It is not needed for statement (*) over SO(4).
- The coordinates are the paper's. The supplement's "original 14-dimensional lift" agrees with the paper's m_0..m_13 and kappa = e_13: the sample vector (19) is reproduced exactly from the paper's f_j, and the supplement's cubic F equals the paper's C exactly, sign included.

## A. Statement

Let J = {0, 2, 3, 7, 11, 12, 23}, the paper's Table 5 seven-set #22 (EXACT). Note that T_32 = {0,2,3,7,11,12} is a subset of J.

Theorem 1: let x be in X and pi* in J, and suppose P c(x) >= 0 with (P c(x))_pi = 0 for every pi in J; extra contacts are allowed. Then c(x) lies in cone{f2, f7, f8, f12, f15, f39}, which are all minor vectors. Hence Phi_{c(x)}(O) >= 0 for every O in SO(4), and for every O in O(4).

The theorem covers all spectra, degenerate ones included: the argument uses only that 2 Re x is a sum of two points of X, which is contained in Sing C, together with theorem (1). That holds for every matching vector, and (1) is a closed condition (see C). The statement is stronger than (*) for J.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| generators(J) = {2,7,8,12,15,18,39,40}, the contact kernel, and that the f_j are the paper's | paper App. D / lib | EXACT (lib, plus our checks) |
| X is contained in Sing C (grad C = 0 on X); X is conjugation-closed | paper Sect. 3 / A.1 | EXACT in lib (grad C = 0 on X); conjugation is PROOF-READ |
| Phi_{n_{I,J}} = det O[I,J]^2 >= 0 (4.2) | paper | EXACT/PROOF-READ in lib |
| the group action preserves X, contact sets and Phi (4.7) | paper / lib | EXACT in lib |
| (1): det H_C(w)[I,I] >= 0 for every 12-set I and every realizable real w = x + conj(x) | supplement "global matching Hessian inertia theorem" | CITED; NUMERIC support only |
| Theorem 1.2 (unitary to orthogonal) | paper | CITED (only for the unitary remark) |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | the 7 profiles [2,7,8,12,15,39,18] are independent, vanish on J, and rank P[J] = 7, so they span the contact kernel | identity | EXACT | rank computations |
| 2 | gap readings (3): (g21,g20,g9,g6,g18,g1) = (a,b,c,e,f,g); (g13,g16,g5,g17) = (a+rho, b+2rho, f+rho, g+rho) | identity | EXACT | all 24 gaps computed |
| 3 | the unique circuit (4) f18+f40 = f2+2f7+f15+f39; so for rho<0, (5) has nonnegative coefficients | identity | EXACT | 24 values; the 8 generators have rank 7 |
| 4 | the element (col,row,eps) = ((1,2,0,3),(1,2,0,3),1) follows the stated convention, stabilizes J, fixes the six minor vectors and swaps f18 with f40; it is the unique such element | identity | EXACT | lib group action |
| 5 | the reduction to excluding a realizable source (6) with d > 0 | logic | PROOF-READ | |
| 6 | if c is realizable then w = c + z kappa satisfies F(w) = 0 and det H(w) = 0 for some z: 2 Re x = x + conj(x) with both in Sing C; the trilinear form gives F(x+y) = 0 and H(x+y)(x-y) = 0; when x = y, Euler's identity applies; x = 0 is excluded by d > 0 | analytic | PROOF-READ | auditor re-derived the polarization |
| 7 | (7)-(8): C(w0 + z kappa) = C + zQ with the displayed C and Q; both have positive coefficients and contain d^3 and d^2 respectively | identity | EXACT | flint |
| 8 | B has 38 positive terms with d^6 coefficient 1; G has 64 terms whose negative part is -4acefg(b+d) | identity | EXACT | |
| 9 | (10): Q^4 det H(w0 - C kappa/Q) = -4 Q^2 C^2 G^2, so a realizable source has G = 0 | identity | EXACT | 14x14 Bareiss with 819 exact divisions over Z[a..d,z]; the z-degree is 3 |
| 10 | (11): Q^4 times each of the three 13x13 border determinants is divisible by G | identity | EXACT | 3 Bareiss runs (1950 divisions), then exact division by G |
| 11 | gcd(A,B) = 1; G is irreducible | identity | EXACT | flint gcd and factor |
| 12 | (12)/(19): sample values C = 320/3, Q = 40, z = -8/3, G = 0; the vector w as displayed; 3H integral; det(3H[I,I]) = -310542336000 | identity | EXACT | We also found rank H = 12 at the sample, with all 91 principal 12-minors <= 0 (87 nonzero) |
| 13 | (13): rank H <= 12 on G = 0 with AQ != 0, via the field Q(b,c,e,f,g,d) with a = -B/A; Sylvester's Schur-border identity; det H_I is a nonzero rational function because it is nonzero at the sample (where A, Q != 0) | algebraic | PROOF-READ | rests on 10 and 12 |
| 14 | (14)-(16): along alpha -> w(alpha), D(alpha) = -4C^2(alphaA+B)^2/Q^2 and D''(a) = -8C^2A^2/Q^2 < 0; a matrix of corank k makes the determinant vanish to order >= k; so rank >= 12 | analytic | PROOF-READ | uses that A and B do not depend on a, and that A < 0 |
| 15 | at a nonnegative root with d > 0: B > 0 gives a > 0 and A < 0, which forces c,e,f,g > 0 | inequality | PROOF-READ, with the A formula EXACT | |
| 16 | (17)-(18): A < 0 iff fg > kappa0(f+g+1); the parametrization gives A = -4ce(b+1)xi zeta; it is a bijection onto the root set with d = 1, so the root set is connected | identity/topology | EXACT (identities) + PROOF-READ (surjectivity: (f-k0)(g-k0) > k0(k0+1) forces f, g > k0) | sympy |
| 17 | constant rank 12 on a connected set gives constant inertia; the nonzero principal 12-minors of a rank-12 form share one sign; the sample gives sign -; so every branch point has some principal 12-minor < 0, contradicting (1) | linear algebra | PROOF-READ | (Cauchy-Binet: H = M^T D M) |
| 18 | conclusion: rho = 0, so c is a nonnegative sum of minor vectors | logic | PROOF-READ | |

Degenerate spectra and extra contacts: no genericity is used. (1) is needed at every realizable w, degenerate spectra included. If it holds on the dense set of generic spectra, it holds everywhere, because principal minors are continuous and ">= 0" is closed (argument supplied by auditor).

## D. Falsification

- NUMERIC, cited (1): 200 samples, seed 20261125, a third of them with a real spectrum a. The minimum normalized principal 12-minor was 1.1e-17, so no violation. Every sample had numerical rank 12.
- NUMERIC, main statement: `random_point_on_contact_stratum(J)` found realizable points on 8 of 8 seeds. At all of them |rho|/max gap <= 7.6e-11, consistent with rho = 0 for every realizable source.

## E. Verdict

**confirmed-modulo-citations.** Every computational claim was replayed exactly, including the counts the document states (819 divisions, 38 and 64 terms, det = -310542336000). The logic from these identities to the conclusion is sound. The single CITED premise is the Hessian inertia theorem (1); it passed a numerical test. No gaps were found.
