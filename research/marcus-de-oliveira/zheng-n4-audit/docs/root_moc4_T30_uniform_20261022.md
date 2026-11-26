# Audit: root_moc4_T30_uniform_20261022 (primary proof of paper case T_31)

Evidence level: EXACT for the generators, the simplicial-cone claim, C and Q, the allocation algebra, the strict constant estimate (5), all seven domains with the domain-1 cut and the domain-3 univariate certificates, the zero-coefficient checks and the orbit witnesses of the three cited bounds. The radial identity Q^4 det H = -4 Q^2 G9^2 is EXACT at 10 random integer points (randomized identity test, error probability < 1e-42). PROOF-READ: the radial-necessity argument (supplied by the auditor), Lemma 1, domain coverage and the closure/limit arguments of Section 4. CITED: the old-H, chain and four-coefficient pointwise inequalities, backed by NUMERIC evidence. NUMERIC: a brief test of the main statement.

Document: proofs/manin_transfer_20261021/root_moc4_T30_uniform_20261022.pdf, 4 pages. Replay: `replay/root_moc4_T30_uniform_20261022__check.py`, 53 exact checks, all pass; output in `__check_output.txt`.

## Deviations and limits

- G9 (1183 terms) was loaded from `root_moc4_T30_full_radial_20261022.json` as data. It was not recomputed symbolically from the 14 x 14 Hessian determinant. My own C and Q agree with the JSON exactly. The identity Q^4 det Hess(g - (C/Q) kappa) = -4 Q^2 G9^2 was checked exactly at 10 random integer points (coordinates in [1, 10^6], d = 1). Both sides have degree <= 50 after multiplying by Q^10, so this is a Schwartz-Zippel test with failure probability < 1e-42, not a symbolic proof.
- My tables come from my own denominator clearing, so term counts differ from the document's. For domain 3 I get 109402/205803 terms against the stated 18189/33322, and 88 univariate groups instead of 21. The conclusions are identical.
- The three coefficient bounds (2)-(4) are cited as "previously proved uniform response inequalities". No source document is named. I checked exactly that the stated group witnesses exist, i.e. that some g in G maps each base profile tuple to the target tuple. I checked the transported inequalities only numerically.
- Main statement, brief test: 4 floating-point points on the T_31 stratum, degenerate and with many extra contacts, 10 BFGS restarts each. The worst min Phi was +2.6e-18. This was a scratch run and is not in the replay file.

## A. Statements

Theorem 1 (p. 1). Take T = {0,1,5,6,9,12}, which is literally the paper's T_31. If x in X has all 24 gaps >= 0 and zero gaps at T, then Phi_{c(x)}(O) >= 0 for every real orthogonal O. Contacts are non-strict, and extra contacts and zero coefficients are allowed. The proof uses only x in X and is algebraic on the closed contact cone, so it covers all spectra, degenerate ones included. The statement is exactly (*) for T_31.

Lemma 1 (p. 2). Let a, b, c > 0 and 0 <= t < min(ab, 1). The cost (1) with d = 1 satisfies the support inequality if e >= E(t) = t(a+1)(b+t)/(ab-t) and j >= J(t) = (1-t)/a + (1-t)^2/c. The lemma does not use spectral realizability.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| Prop. 4.1 / Lemma 4.2 | core-doc (paper) | generator list = lib `generators(T_31)` EXACT. The cone is moreover simplicial by the document's own identity-matrix argument (EXACT), so Lemma 4.2 is needed only for "these 8 span the kernel" (rank EXACT) |
| X subset Sing C, Euler | core-doc (paper App. A) | lib-replayed EXACT |
| four-coefficient inequality, base (0,6,10,54) -> (7,0,5,35): ab(e-1) >= a+b+e+1 | common-lemma | CITED, source unnamed; witness g EXACT; NUMERIC worst min Phi at the boundary 1.7e-17 |
| old-H, base (3,5,54) -> (9,1,35): j >= 1 + 1/c | common-lemma | CITED; witness EXACT; NUMERIC 8.3e-15 |
| chain, base (2,0,12,70) -> (1,7,9,35): j >= 1/a + 1/c | common-lemma | CITED; witness EXACT; NUMERIC 2.7e-15 |
| two-orthogonal theorem (paper Thm 1.2) | core-doc | CITED |

No other contact case is used, consistent with T_31 in Stage 1. Corollary B.3 is not used.

## C. Steps

| # | claim (location) | kind | status | how checked |
|---|---|---|---|---|
| 1 | contact equations rank 6; the 8 profiles {2,29,1,9,5,35,0,7} take identity-matrix values at permutations (15,4,17,7,11,3,18,14), so the coefficients are gaps >= 0 and (1) parametrizes the whole cone (p.1) | identity | EXACT | values; the 8 vectors are independent and lie in the 8-dim kernel |
| 2 | all profiles except f35 are minors (p.1) | identity | EXACT | lib flags |
| 3 | bounds (2)-(4) with stated orbit witnesses (p.1) | citation | CITED + EXACT witnesses + NUMERIC | |
| 4 | Lemma 1 and E'(t) formula; E increasing; Emax; failure of (2) iff e < Emax (p.2) | identity | EXACT (sympy) + PROOF-READ | the homogeneous four-coefficient bound solves to e >= E(t) |
| 5 | radial necessity: F(g + z kappa) = C + zQ = 0 and det Hess = 0 for real z (p.2) | analytic | PROOF-READ, argument supplied by auditor | as in the T_29 audit: g + z kappa = (y + ybar)/2, with polarization and Euler |
| 6 | C has only negative coefficients and contains -d^3, so Q != 0 (p.2) | identity | EXACT | own cubic = JSON |
| 7 | Q^4 det H(g - C/Q kappa) = -4 Q^2 G9^2 (p.2) | identity | EXACT at 10 random integer points (randomized) | own Hessian, flint determinant |
| 8 | G9 quadratic in j; k2 >= 0 coefficientwise; k0 at c=e=q=p=0 equals B; k0 at e=0 minus B >= 0 (p.2) | identity | EXACT | |
| 9 | (5) k0 >= B(1 - Le e)^2 in the four charts (a,b) = (zR, z/R) and swapped, z = T/(1+T) or 1+T (p.2-3) | inequality | EXACT | four tables, nonnegative; 1 - Le e = 1/(1+V) in the ab>1 chart |
| 10 | the seven domains cover all (a,b) >= 0 and each gives 0 <= j < J when not covered by a bound (p.3) | case split | PROOF-READ (auditor checked each parametrization's image) + EXACT (chart numerators/denominators nonnegative) | |
| 11 | Bernstein numerators nonnegative in domains 2-7 (p.3) | inequality | EXACT | the domain-3 endpoint by grouping in c: 88 groups, each primitive part checked >= 0 on c >= 0 by squarefree factorization and real-root isolation |
| 12 | domain 1: cut l = a - b; the [0,l] and [l,J] numerators are nonnegative; (7) G9(l) >= B (p.3) | inequality | EXACT | |
| 13 | strict positivity (6) on 0 <= j < J (p.3) | inequality | PROOF-READ from steps 9, 11, 12 | |
| 14 | closures: ab = 1, a = b, a = 1, b = 1; a = 0 (all coefficients >= 0, constant 1); c = 0 (k1, k2 >= 0 and (5)); b = 0 by the limit in domain 3; zero p, q, e (p.4) | analytic | EXACT (a=0, c=0 checks) + PROOF-READ (limits) | |
| 15 | conclusion: a realizable cost with negative response would need G9 = 0, contradicting G9 > 0 (p.4) | logic | PROOF-READ | |

## D. Falsification

- Cited bounds at their boundaries: 6 random samples per family, 12 BFGS restarts each (rng seed 7). The worst values were 8.3e-15 (old-H), 2.7e-15 (chain) and 1.7e-17 (four-coefficient). All three are tight, with minimum about 0, and none is violated.
- The intermediate positivity claims were replayed exactly, not sampled.
- Main statement: brief test, no violation found (see Deviations).

## E. Verdict

**confirmed-modulo-citations.** The citations are the three pointwise inequalities (old-H, chain, four-coefficient). The PDF does not name their source, so they are unresolved CITED dependencies. They are numerically consistent and tight, and their orbit witnesses were checked exactly. The paper's Lemma 4.2, App. A and Thm 1.2 are also cited. Every certificate in the document was replayed exactly. The one non-symbolic link is the randomized check of the radial identity, at 10 points with failure probability < 1e-42. The radial-necessity argument, which the document states in two sentences, was verified in full by the auditor.
