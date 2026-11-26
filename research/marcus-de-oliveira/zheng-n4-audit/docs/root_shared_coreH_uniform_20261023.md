# Audit: root_shared_coreH_uniform_20261023 (seven-contact set H)

Evidence level: every identity and every positivity certificate EXACT (own code from the definitions); the case logic and the sign bound (6) PROOF-READ; the sharp-SOC and old-H atoms CITED and NUMERIC-tested. Verdict: **confirmed-modulo-citations**.

PDF: `proofs/manin_transfer_20261021/root_shared_coreH_uniform_20261023.pdf` (3 pages). Replay: `replay/root_shared_coreH_uniform_20261023__check.py` (output in `__check_output.txt`; exits 0).

## Deviations and limits

- G (the document's degree-8 polynomial) was re-derived by my own exact Bareiss determinant of the Hessian of (A.2) along the lift. Nothing from the supplement was run or read as data.
- The atoms are CITED. The document names the old-H source as the "coreD inequalities", and gives no named source for the SOC reference family. Group maps and numeric tests only.
- The numeric falsification is weak: all realizable floating points found on the stratum have d of about 0.

## A. Statements

**Theorem 1 (p. 1).** Let H = {0,1,6,9,16,17,22}. This is literally the paper's Table 5 seven-set no. 21 (EXACT). Take x in X with P c(x) >= 0 and (P c(x))_pi = 0 for pi in H (non-strict). Then Phi_{c(x)}(O) >= 0 for all O in SO(4) (equivalently O(4)).

- The proof uses only Re x in Sing C + Sing C, so all spectra are covered (degenerate included).
- The remark "applications to the six-contact cones with indices 80 and 82 require their own reductions" is outside this statement.
- Coordinates: g = a f15 + b f3 + q f2 + t f5 + u f29 + v f33 + d f31. All seven coefficients are gaps: rows [5,14,15,11,8,3,7] form I_7 (EXACT). f31 = h_{4123,+} is the only non-minor. The basis equals generators(H) (EXACT).

## B. Dependencies

| dependency | kind | status | use |
|---|---|---|---|
| Sharp SOC atom (f1,f9;f22) -> (f15,f3;f31): a f15 + b f3 + d f31 is safe if ab >= d^2 | other-doc | CITED; group map EXACT; NUMERIC 20 samples at ab = d^2, worst 9.1e-17 | (2) |
| Old-H atom (f3,f15;f26) -> (f2,f15;f31), (f29,f15;f31), (f5,f3;f31), (f33,f3;f31), as ordered (aux, main, non-minor) | core-doc (coreD per the document) | CITED; four ordered group maps EXACT; NUMERIC 4 x 20 samples, worst 2.6e-16 | (3),(4) |
| (4.2) minor vectors give squared minors; X in Sing C; conj X = X | other-doc | EXACT in shared lib / PROOF-READ | d = 0 case, remainders, necessity |

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | (1): canonical basis, rank 7, rows give I_7, coefficients are nonnegative gaps | identity | EXACT | own linear algebra |
| 2 | d = 0: sum of minors | identity | EXACT | minor flags + (4.2) |
| 3 | (2) and (3): SOC and four old-H atoms | citation | CITED | group maps EXACT, numeric |
| 4 | aggregation: xi(a - delta) >= delta^2 splits into the two old-H conditions in proportions q/xi, u/xi | inequality | PROOF-READ | q^2(a - delta)/xi >= (delta q/xi)^2 is equivalent to xi(a - delta) >= delta^2 |
| 5 | (4): lambda + mu >= d is sufficient; the resources are disjoint | inequality | PROOF-READ | cap identity EXACT (sympy); the old-H condition is monotone in delta |
| 6 | (5): F(w(z)) = C + Qz with C, Q as printed | identity | EXACT | own expansion |
| 7 | (6): at failure with d > 0, C <= -d xi psi - (3/4) d^2 (xi + psi) - d^3 < 0 | inequality | PROOF-READ (decomposition EXACT) | cap(X,m) < d implies X(m - d) < d^2 (trivially if m <= d); P <= xi^2/4 gives P(a-d) <= d^2 xi/4, also when a < d |
| 8 | realizable w gives F = 0 and det H = 0, so Q != 0 | analytic | PROOF-READ | as for core G (polarization; Euler identity when x is real; w != 0 since d > 0) |
| 9 | (7): Q^4 det H(w(-C/Q)) = -4 d^2 Q^2 G^2; G has 499 terms, degree 8, coefficient 1 at d^8 | identity | EXACT | own Bareiss, 819 exact divisions; R = dG |
| 10 | (8): three symmetries of G; the failure tests are invariant | identity | EXACT / PROOF-READ | polynomial equalities; lambda and mu depend only on (xi,a) and (psi,b) |
| 11 | (9): a,b < 1 (709 coefficients, constant 1); wall a = 1, b < 1 (235, constant 2); a = b = 1 is safe by (2) | certificate | EXACT | counts match |
| 12 | a > 1 implies ab < 1; (10)-(12) cover a > 1, b > 0; certificate has 476039 positive coefficients, marker 2Y^4; V_n and V_D have positive coefficients | certificate + inverse | EXACT | inverse map verified with sympy, including D = 1 + L + Y; Y > 0 because b > mu |
| 13 | (13): b = 0, a > 1; 2356 positive coefficients, constant 2; K^2/D0 = lambda^2/(a - lambda) | certificate + inverse | EXACT | sympy |
| 14 | orientation: put a nonzero auxiliary as q (resp. t), U = u/q, V = v/t; total zero gives ratio 0 | case split | PROOF-READ | q<->u and t<->v are independent exact symmetries of G and of the failure tests |
| 15 | exhaustiveness after d = 1 | case split | PROOF-READ | a,b <= 1, or one main > 1 (by symmetry a), with b > 0 or b = 0 |

## D. Falsification

- Atoms: 5 atoms x 20 samples at the boundary of each atom's own condition, BFGS with 8 restarts, seed 20261126. The worst minimum is 2.6e-16.
- Theorem: 12 stratum points (seeds 0..11), 10 restarts each. The worst min Phi / max gap is -1.2e-10 (roundoff).
- At these points C(Re x) is about 0 and H(Re x) is singular to 1e-18 relative (scratch check).

## E. Verdict

**confirmed-modulo-citations.** Citations: the sharp-SOC atom and the old-H atom (coreD). No gaps. The sign argument (6), which the document flags as needed because C is not sign-definite, is correct. All four certificate tables reproduce exactly, term counts included.
