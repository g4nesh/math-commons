Evidence level: EXACT for the radial factor G (rebuilt from the paper's definitions), the allocation algebra, the chart identities, the full 149933-term sign certificate with its 18 square groups, the closed-wall homogenization and every zero-face embedding; PROOF-READ for the logic; CITED for the six zero-face contact cases; NUMERIC only for a falsification run.

# blind_moc4_T48_allocation_II_uniform_20261023 (2 pp.): mixed chain / old-H allocation cell of T_49

Auditor label C26-T49. Script: `replay/blind_moc4_T48_allocation_II_uniform_20261023__check.py` (output beside it; 0 EXACT failures). This script computes G from scratch when run with --fresh (144-208 s; the captured output is such a run); the other C26-T49 scripts reuse it through a hash-checked cache or recompute it with `--fresh`.

## Deviations and limits

- **G rebuilt independently.** Let c = q f7 + p f5 + j f29 + c f39 + e f2 + d f40 + b f41 + a f0 in V_R, let C be the cubic (A.2) and kappa = e_13. My code finds C(c + h kappa) = C0 + h Q0, computes det Hess C(c + h kappa) exactly by Bareiss elimination (degree 3 in h), and checks -Q0^3 det Hess C(c - (C0/Q0) kappa) = 4 G^2 Q0 with G a polynomial. So det Hess = -4G^2/Q0^2 exactly as the document says. G has 1772 terms, is homogeneous of degree 9 and is primitive. Only afterwards did I load the supplement's stored radial factor as data: it equals **-G**. All signs below are in my convention, so "all coefficients positive" here means "all negative" in the document.
- **Scope of the document.** generators(T_49) has ten elements, f0, f1, f2, f5, f7, f9, f29, f39, f40, f41. The document's eight roles omit f1 and f9. EXACT: f1 = f7 + f39 + f0 - f40 and f9 = f7 + f2 + f41 - f40, so every source has a unique eight-role form, in which d becomes d - x1 - x9. The theorem covers only sources whose eight-role coordinates are all >= 0. The other sources (the "full-source absorptions") are handled elsewhere. Observation (NUMERIC, seeds 0-11 in the script, and 40 more in a scratch run): every realizable point sampled with contact set exactly T_49 (12 + 40) had a negative f40 coordinate, so none of them was in this document's scope.
- The cell is one piece of a cover of the eight-coefficient face. I did not check that the cells of the T_49 notes cover the face; that is the whole-assembly document's job.
- I did not read the supplement's verifier scripts. Their published hashes are of their own canonical form; mine differ only by the overall sign convention and are recorded in the output.

## A. Statement

Theorem 1. Fix the eight-role source with d = 1 (homogeneity) and a, e, j > 0. Put Q = q + b, u = (1 + 2/j) / (2(1/a + 1/e + 1/j)) and v = 1 - u. Let x be any point of X (finite or not, collisions allowed) whose gap vector c(x) has these coordinates, all >= 0, and suppose 0 < u < 1, 2u <= a and Q <= 2v. Then Phi_{c(x)}(O) >= 0 for every O in SO(4). The equality walls a = 2u and Q = 2v are included, as is c = 0. The coefficient zeros q, p, j, e, b, a = 0 are delegated to other cases, and d = 0 is a sum of minors.

Not the primary proof of a paper case; it is a branch of T_49 = {0,3,5,6,9,12}, which is paper T_49 exactly (literal tuple).

## B. Dependencies

| dependency | use | check |
|---|---|---|
| chain atom, (centre f5; leaves f0, f2; f40), (f5; f7, f29; f40), (f5; f29, f41; f40) | pointwise allocation | EXACT: all three are in the proved 576-member chain orbit (foundations/so4-atoms.md, 3.3), with the roles as printed |
| old-H atom (aux f29, main f5; f40) | j-group current h | EXACT: in the proved 288-member orbit with aux = f29, main = f5; the reversed order is not in the orbit. The coefficient form t(c - d) >= d^2 matches the use (t = j h / v, d = h) |
| radial necessity (paper p. 27; audit item 5) | every realizable source has det Hess C(u) = 0 at u = c + h kappa with C(u) = 0 | established; combined with the EXACT identity above it gives G = 0 when Q0 != 0 |
| 42-term identity (e+j)C - ejQ0 | excludes C = Q0 = 0 | EXACT: 42 positive coefficients including d^3 e and d^3 j |
| zero faces q, p, j, e, b, a = 0 -> historical T38, T0, T3, T23, T42, T43 | resource-zero faces | EXACT: the contact set of each face contains an image under G of paper T_39, T_1, T_3, T_24, T_43, T_44 respectively (group elements printed). The proofs of those cases are CITED |
| supplement files "retained atoms" and "zero faces" JSON | identifications | not used; I re-derived both |

## C. Steps

| # | claim | kind | status | how |
|---|---|---|---|---|
| 1 | roles = eight of the ten generators of T_49; they vanish together exactly on T_49; f40 is the only non-minor | identity | EXACT | lib `generators`, F_VALUES |
| 2 | (a,e) chain with current u, j split h : q/2 : b/2 with old-H and two chains; totals are exact (e, a, q, b full, j shared exactly, non-minor currents sum to 1) | identity | EXACT | sympy |
| 3 | requirement = u^2(1/a+1/e) + v - Q/4 + v^2/j | identity | EXACT | sympy |
| 4 | u(1/a+1/e) - v/j = 1/2 and p* = 1 - u/2 + v/j - Q/4 >= 1/2 + v/j | identity, inequality | EXACT / PROOF-READ | sympy; the bound uses Q <= 2v and u + v = 1 |
| 5 | p >= p* implies the support inequality (sum of proved atoms) | citation | CITED (foundations, proved) | orbit checks of step B |
| 6 | realizable implies G = 0 (Q0 != 0 by the 42-term identity) | identity + citation | EXACT + established radial necessity | own Bareiss computation |
| 7 | the chart (L, A, B, C1, S, P, c) is onto the strict region; e is forced by the definition of u; EN > 0 | identity | EXACT | sympy: u(a,e,j) = L/(1+L), inverses of A and B; the other inverses are immediate |
| 8 | D G with D = 4S^2(1+L)^9(1+B)^6(1+C1)^6(1+P)^2 EN^2 is a polynomial with 149933 terms; all coefficients of one sign except 54, which lie in 18 C1-groups, each exactly (C1-1)^2 R(C1) with R >= 0; marker 432 B^6 S^2 with a nonnegative group | certificate | EXACT | flint, every term |
| 9 | the wall Q = 2v: homogenize B to degree 6 in (B, WB), set WB = 0 | certificate | EXACT | direct substitution with Q = 2vB/(B+WB); same sign pattern; marker survives WB = 0 |
| 10 | the wall a = 2u is A = 0; c = 0 | certificate | EXACT | marker independent of A and c, all variables >= 0 |
| 11 | p = p*: the capacity itself; p = 0 and the other zeros: zero faces | case split | PROOF-READ + CITED | step B |

Exhaustiveness: the hypotheses give q, p, j, e, b, a > 0, c >= 0 and p < p* (strict failure), or one of the delegated faces. The chart is onto that set (step 7), so the cases are exhaustive.

## D. Falsification

- NUMERIC: 12 random cell points with p = p* exactly (seed 11; not necessarily realizable), minimum of Phi over SO(4) by multistart BFGS (8 restarts): worst value -2.3e-16 relative to the coefficient scale. So the allocation is sharp, and no violation was found.
- NUMERIC: 12 realizable points on the T_49 stratum (seeds 0-11, margin 2e-2) satisfy G = 0 to 6e-17 relative accuracy. None lay in this cell (all had d < 0 in the eight-role form), so the main statement could not be sampled directly.

## E. Verdict

**confirmed-modulo-citations.** Every identity and certificate replays exactly, and the term counts match. The only citations are the six zero-face cases (T_39, T_1, T_3, T_24, T_43, T_44), whose embeddings are verified, and the proved atoms and radial necessity, which this audit has established. Presentation: "The face c = 0 is not covered and was stated throughout" contradicts the theorem's "This includes ... c = 0". Mathematically c = 0 is included, because the marker is independent of c. The phrase "the list of 42 previously established contact-case inequalities" is unexplained, but only the six cited faces are needed.
