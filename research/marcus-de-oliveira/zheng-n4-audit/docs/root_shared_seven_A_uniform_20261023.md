# Audit: root_shared_seven_A_uniform_20261023 (seven-contact set {0,1,2,7,9,10,17})

Evidence level: all identities and the 15 chart certificates EXACT (own code from the definitions); chart coverage and case logic PROOF-READ; the sharp-SOC, old-H and chain atoms CITED and NUMERIC-tested. Verdict: **confirmed-modulo-citations**.

PDF: `proofs/manin_transfer_20261021/root_shared_seven_A_uniform_20261023.pdf` (4 pages). Replay: `replay/root_shared_seven_A_uniform_20261023__check.py` (output in `__check_output.txt`; exits 0).

## Deviations and limits

- G (degree 9, 486 terms) was re-derived by my own exact Hessian determinant. The certificates were rebuilt with my own substitution and clearing: each coordinate's denominator is raised to that coordinate's degree in G. All 15 term counts and markers coincide with table (9).
- The atoms are CITED. Their source documents are not named; the document gives only the reference profile triples [1,9,22], [11,1,34], [2,0,12,70].
- The numeric falsification is weak: the realizable stratum points found all have d of about 0.

## A. Statements

**Theorem 1 (p. 1).** Let T = {0,1,2,7,9,10,17}. This is literally the paper's Table 5 seven-set no. 9 (EXACT). Take x in X with P c(x) >= 0 and (P c(x))_pi = 0 on T (non-strict). Then Phi_{c(x)}(O) >= 0 for all O in SO(4).

- All spectra are covered: only Re x in Sing C + Sing C is used.
- Coordinates: g = q f0 + p f2 + t f3 + b f12 + c f15 + v f21 + d f26. Rows [23,15,14,8,4,16,11] form I_7, so all coefficients are nonnegative gaps (EXACT). The basis equals generators(T), and f26 = h_{4132,-} is the only non-minor.

## B. Dependencies

| dependency | kind | status | use |
|---|---|---|---|
| Sharp SOC (f1,f9;f22) -> pairs of {f2,f12,f15} with f26 | other-doc | CITED; group maps EXACT; NUMERIC 3 x 20, worst -5.8e-16 | H3 certificate |
| Old-H (f11,f1;f34) -> (f3,f15;f26): t f3 + c f15 + d f26 is safe if t(c - d) >= d^2 | other-doc | CITED; ordered group map EXACT; NUMERIC 20, worst 2.6e-17 | (2) |
| Chain (f2,f0,f12;f70) -> (f15,f0,f3;f26): q f0 + t f3 + c f15 + d f26 is safe if c >= d^2/q + d^2/t | other-doc | CITED; group map EXACT (main to c, the two auxiliaries to q and t; the condition is symmetric in q and t); NUMERIC 20 at c = 1/q + 1/t, worst -3.3e-16 | (2) |
| (4.2); X in Sing C; conj X = X | other-doc | EXACT in lib / PROOF-READ | d = 0 case, necessity |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | (1): the basis and the I_7 rows; the coefficients are gaps | identity | EXACT | own linear algebra |
| 2 | d = 0 is immediate | identity | EXACT | minor flags |
| 3 | H3 sufficiency from the three SOC atoms (balanced and dominant allocations; zero denominators) | inequality | PROOF-READ (atoms CITED) | capacities sum to s/2, or to sqrt(p(b+c)) in the dominant case |
| 4 | (2): failure gives H3 < 1, t(c-1) < 1, and cqt < q + t for q,t > 0; (3): a main >= 1 is uniquely dominant, and p(b+c) < 1 etc. | inequality | PROOF-READ | H3 >= sqrt(xy) for every pair; if M >= 1 and 2M <= s then s/2 >= 1 |
| 5 | F(w(z)) = C + zQ; C > 0, so Q != 0 | identity | EXACT | own expansion |
| 6 | realizable w gives F = 0 and det H = 0 | analytic | PROOF-READ | as for core G |
| 7 | (4): Q^4 det H(w(-C/Q)) = -4 Q^2 G^2; G has 486 terms, degree 9, coordinate degrees (3,2,2,2,2,3,9) | identity | EXACT | own Bareiss, 819 exact divisions |
| 8 | c = 0 factorization; positive when pb < d^2 | identity + inequality | EXACT / PROOF-READ | exact factor; d(d^2 - pb) > 0 |
| 9 | (5): main charts R0, Rp, Rb, Rc and the c = 1 chart cover every (p,b,c) allowed by (2)-(3) with c > 0 | case split + inverse | PROOF-READ (inverse EXACT) | R0 for all < 1; Rp, Rb for p or b >= 1 (P = 0 is the wall); Rc for c > 1; c = 1 forces p + b < 1 |
| 10 | (6), (7), (8): auxiliary rows cover all (q,t) allowed by old-H and chain failure, including q = 0, t = 0 and the interfaces cq = 1, q = 1 | case split + inverse | PROOF-READ (inverses EXACT) | for c < 1, old-H failure is automatic; cq <= 1 gives no chain bound; cq > 1 gives t < q/(cq-1); for c > 1, q <= 1 uses old-H t < 1/(c-1) |
| 11 | (9): 15 cleared numerators have all coefficients positive with the stated markers, and the denominators are positive | certificate | EXACT | all 15 counts and marker coefficients match; markers use only the strict parameters (R; R,U; P; U) |
| 12 | conclusion: G > 0 on the unsafe residual contradicts G = 0 | logic | PROOF-READ | |

## D. Falsification

- Atoms: 5 atoms x 20 samples at their boundaries, 8 restarts, seed 20261127. Worst minimum -5.8e-16 (roundoff).
- Theorem: 12 stratum points, 10 restarts. Worst min Phi / max gap -9.1e-12.
- Necessity sanity check at these points: C(Re x) is about 0 and H is singular. R != 0 occurs only where Q is about 0 or d = 0.

## E. Verdict

**confirmed-modulo-citations.** Citations: the sharp-SOC, old-H and chain atoms.

Presentation notes:
1. p. 2: "an uncost satisfying the support inequality must satisfy (2)" should read "a cost violating the support inequality must satisfy (2)".
2. The source documents of the three atoms are not named.
