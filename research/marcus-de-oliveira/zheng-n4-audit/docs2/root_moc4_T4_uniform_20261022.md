Evidence level: EXACT for every identity and the four-coefficient bound (own replay, 44/44 exact checks pass); PROOF-READ for the sign logic and the case split; CITED for paper Theorems B.1, C.1, 1.2 / Corollary 2.2, eq. (4.2) and the classical quaternion description of SO(4); NUMERIC falsification only as stated in D.

# Second-pass audit: root_moc4_T4_uniform_20261022 (paper case T_4)

Supplement PDF `proofs/manin_transfer_20261021/root_moc4_T4_uniform_20261022.pdf`, 5 pages: "The contact family T4: Hermitian forms and a four-coefficient inequality". Auditor label C02-T4-T5-second. This audit is independent: I did not open the first auditor's report or any other file under `docs/`, `docs2/` or `replay/`.

Replay: `replay2/root_moc4_T4_uniform_20261022__check.py`, with its output in `replay2/root_moc4_T4_uniform_20261022__check_output.txt`. It uses only `lib/moc4lib.py` and sympy/numpy/scipy. Nothing from the supplement was executed, imported or loaded.

## Deviations and limits

- To prove the imported four-coefficient bound (1), I replayed the polynomial certificate (4) of the companion note `final_referee_T5_four_budget_capacity_20261022` (p. 1). I re-typed it from that PDF's page image and checked it with my own code. I also had to establish its link to the Leibniz sums myself (check C9), because neither note states it in usable form. That note is not otherwise audited here.
- The quaternion parametrization O(p,q): x -> p x q of SO(4) is surjective. This classical fact is cited, not re-proved. Check C8 verifies exactly that O(p,q) is orthogonal for unit p, q.
- The document's own "arrowhead pencil" argument for (1) does not work as written (gap P1 below). I did not try to repair it. The replayed companion certificate proves (1) instead.
- The numeric falsification is floating point (BFGS over SO(4)) and gives evidence only.

## A. Statements

**Theorem 1 (p. 1).** Write F_sigma = prod_i (a_i - c_sigma(i)) with c = -b. Take any complex kappa and put g_sigma = Re kappa(F_sigma - F_0), with reference permutation 0. Suppose g_sigma >= 0 for all sigma and g_sigma = 0 for sigma in T4 = {0,1,4,5,7,8} (one-line 1234, 1243, 1423, 1432, 2143, 2314). Then the support inequality holds for every U in U(4) and all complex spectra, including repeated values and a_i = c_j.

- In paper notation: for every x in X with c(x) in N_{T_4}, Phi_{c(x)}(O) >= 0 for all O in SO(4). This is (*) for T_4.
- The contacts are non-strict ("contains T4").
- The doc's kappa equals the paper's alpha. The common Moebius factor is absorbed into it.
- The contact set is literally the paper's T_4 (check A1), so no group element is needed.
- The statement covers all spectra. The degenerate strata are handled in the text by Theorems C.1 and B.1, and the real-v stratum by an argument I supplied (step S11).

**Embedded lemma (1) (p. 2).** For a, b, c, d >= 0 with ab(c-d) >= d^2(a+b+c+d):

Phi_{a f6 + b f0 + c f5 + d f34}(O) >= 0 for all O in SO(4).

Here f0 = n_{1,4}, f5 = n_{12,24}, f6 = n_{4,1}, and f34 = h_{4321,+} = f0 - f2 + f6 + f33 (paper numbering, Appendix D).

## B. Dependencies

| dependency | used for | status |
|---|---|---|
| Paper Thm C.1 / Prop. 1.4(i) (one spectrum on a generalized circle) | a repeated value within a spectrum (at most 3 distinct values, so on a circle) | CITED. MOC for those spectra implies Phi >= 0 through (3.5) (PROOF-READ) |
| Paper Thm B.1 / Prop. 1.4(ii) (a_i + b_j = 0) | a shared value a_i = c_j | CITED |
| Paper Thm 1.2 / Cor. 2.2 (two-orthogonal reduction) | passing from SO(4) to U(4) | CITED (not needed for Prop. 4.3 itself) |
| Paper eq. (4.2), Phi_{n_{I,J}} = det O[I,J]^2 | nonnegativity of the minor terms | CITED (the lib README gives a proof) |
| Companion `final_referee_T5_four_budget_capacity_20261022`, identity (4) and Z >= 0 | inequality (1), moved from the T5 rays {f0,f6,f10,f54} to the T4 rays {f0,f6,f5,f34} | **replayed EXACT** (C7-C11). The transport g = (id, (0 2 1 3), 0), i.e. sigma -> sigma o (2 3) in 1-based terms, was found and checked by code. The sign argument for det(g.O) = -1 is PROOF-READ |
| "previously proved three-plus-three generalized-circle theorem" (p. 3) | the case v real | not a paper statement and not identified in the document. **Replaced by an auditor argument** (S11), so it is no longer a dependency |
| Classical: SO(4) = {x -> p x q : p, q unit quaternions} | surjectivity in the proof of (1) | external, standard |
| Lemma 4.2 / Prop. 4.1 | not needed: decomposition (2) is checked directly on the whole 8-dimensional space {c : (Pc)_T4 = 0} | none |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | The doc's gap convention (c = -b, gaps relative to F_0, kappa) equals the paper's c(x) with pi* = 0 | identity | PROOF-READ | definitions |
| S2 | The generators of the T4 face are {f0..f6, f11, f33, f34}, and the doc's minor table matches Appendix D | finite | EXACT | A2, A3 |
| S3 | f34 = f0 - f2 + f6 + f33; relations f0-f1-f2-f3+f4 = 0 and -f0+f2-f6-f33+f34 = 0 | identity | EXACT | A4-A6 |
| S4 | Decomposition (2): Pc = g19 f0 + g16 f1 + g14 f3 + g10 f5 + g15 f6 + g6 f11 + g13 f33 + d f34, with d = g2 - g13 | identity | EXACT | C1-C2, on a basis of the 8-dim space W |
| S5 | d <= 0: the replaced coefficients are (b-A+tau, g16-tau, A-tau, g14-tau, tau, g9, g2), all >= 0 | identity + inequality | EXACT + PROOF-READ | C3-C6 (g22 = g16+g19+d, g18 = g14+g19+d, g15 = g9-d hold linearly on W); tau = max(0, A-b) gives the signs |
| S6 | Inequality (1) | inequality with certificate | EXACT | C7 (transport), C8 (O(p,q) orthogonal), C9 (on unit p,q: Phi_f0 = (A-B)^2, Phi_f6 = (A+B)^2, Phi_f10 = D^2, Phi_f54 = H), C10 (certificate (4) mod N), C11 (Z >= 0 under the capacity condition by a Schur chain). With n = 2 and N = 0 on the unit spheres, 2 Phi = sum of l^T Z l + m^T Z m >= 0 |
| S7 | Repeated values give C.1; shared values give B.1; otherwise all eight points are distinct | case split | CITED (exhaustive: PROOF-READ) | |
| S8 | Moebius normalization A = (z,0,1,x), C = (w,v,inf,tv) | analytic | PROOF-READ (argument supplied by auditor) | The map sends c_2 to infinity, and every F_sigma picks up the same nonzero factor [prod_i(gamma a_i + delta) prod_{j != 2}(gamma c_j + delta)] / (gamma det^3). NUMERIC check: relative error 1e-13 (D2b) |
| S9 | Contact ratios on T4 are 0, 1-x, t-x, x(t-1), r-x+1, u. The contacts force t, x real and kappa v(z-w) = i mu with mu != 0. They force r, u real | identity + argument | EXACT + PROOF-READ | C13, C14. Note: the ratio of the differences is -t, not t (cosmetic). x != 1, t != 0, 1 and kappa != 0 by distinctness or nontriviality |
| S10 | The z, w formulas invert (z,w) -> (r,u). Their denominators are nonzero | identity | EXACT | C15-C17. The denominators equal -(v-w)(x-tv)(tv-1)/(z-w) and -(z-v)(x-tv)(tv-1)/(v(z-w)), which are nonzero for distinct points. The doc's stated reason ("z or w would equal infinity") is loose |
| S11 | v real reduces to a known case | case | PROOF-READ (argument supplied by auditor) | With v, x, t real, the six points 0, 1, x, v, inf, tv are real. For sigma(0) = 0, F_sigma = (z-w) times a real number (C37). Contact 1 gives Re kappa(z-w) = 0, since F_1 - F_0 = kappa(z-w) v(1-x) with v(1-x) real and nonzero. So all of {0..5} are contacts. The cone of T4 u {0..5} is generated by minor vectors only (C38), so Phi >= 0 |
| S12 | Scaled gaps (rho = -1): g2 = x/abs(v)^2, g6 = rt(x-1)/abs(tv-1)^2, c = -tux/abs(tv-x)^2, g9 = (x-1)c, a = g9 - d. Hence c > 0, x > 1, rho < 0, rt > 0, tu < 0 | identity + signs | EXACT + PROOF-READ | C19-C22. Scaled gap = (R(v) - R(vbar))/(v - vbar). r != 0 because r = 0 forces z = 0, w = v or tv = 1 |
| S13 | (3) k = hc - t(b+d) and Ee = tyHb - Jd; (4). t < 0 is impossible, so 0 < t < 1, r > 0, u < 0 | identity + signs | EXACT + PROOF-READ | C23-C25 |
| S14 | (5) K = tyb - xhd > 0 when E >= 0 | inequality | PROOF-READ | J > 0, so H, b > 0, and J > xhH |
| S15 | (6), (7): D = d abs(Q(v))^2 and P = K abs(Q(v))^2 are real Hermitian forms of bidegree <= (1,1), with det M_D = -r^3 t^2 uxyLE, D(x/t) = ruxy(hr-Hy), det M_P = -r^2 t^2 u^2 x^2 hyLE, P(x/t) = ru^2 x^2 y(y-h) | identity | EXACT | C27-C31 (the text extract misprints P(x/t) as y^2; the page image reads y). **Q(v) != 0** is needed for the L > 0 contradiction but is not stated. C18 gives Q(v) = -(z-v)(v-w)(x-tv)(tv-1)/(v(z-w)), which is nonzero |
| S16 | E < 0: hr - Hy > 0 (under L <= 0), then L > 0, then y > h, M_P > 0, so K > 0 | inequality | PROOF-READ + EXACT | C32, C33. The 2x2 semidefiniteness logic was checked by hand. The "First L > 0. Otherwise, ..." wording is confusing but correct |
| S17 | ybc - xd(b+d) >= (b+d)K/h > 0 and ab(c-d) - d^2(a+b+c+d) = c{...} > 0, so (1) applies | identity + inequality | EXACT | C34, C35 |
| S18 | SO(4) to U(4) | citation | CITED | Thm 1.2 / Cor. 2.2 |

## D. Falsification (NUMERIC)

All runs use seed 20261125 unless stated.

- **D1.** Inequality (1) on 40 random tuples on the boundary ab(c-d) = d^2(a+b+c+d) (12 BFGS restarts each): the worst normalized min Phi is 2.7e-18. It is >= 0 up to rounding, and 0 is attained at the contact permutation matrices.
- **D2.** Sampled the chart: 200000 points, 2033 admissible (all gaps >= 0, contacts at T4), none with d > 0. Minimizing Phi over SO(4) at 15 of them gave a worst min Phi / sum(gaps) of -2e-15.
- **D3, D4.** Independent least-squares points from the lib (seeds 0..39, margin 1e-3, contact set exactly T4) gave 4 points with d > 0. At each, mapped into the chart:
  - the chart gap vector is a positive multiple of the lib gap vector (relative spread 1.4e-10), which confirms rho < 0;
  - x > 1, 0 < t < 1, r > 0, u < 0 and K > 0 hold;
  - the capacity condition holds strictly;
  - min Phi / sum(gaps) is at least 5e-16.

  So the d > 0 branch is nonempty but thin.
- No violation was found.

## E. Verdict: confirmed-modulo-citations

Every step is EXACT or PROOF-READ, except the citations of paper Theorems C.1 and B.1 (degenerate spectra), Theorem 1.2 / Corollary 2.2 (unitary form), eq. (4.2), and the classical surjectivity of x -> p x q onto SO(4).

The issues below are presentation-level. I repaired each one, so none is left as a GAP step:

- **P1 (p. 2, the arrowhead pencil Z).** The note says (1) "follows independently" from Z >= 0, but gives no vector w(O) with Phi = w^T Z w. As written, that argument does not prove (1). The bound rests on the companion certificate, which I replayed exactly (C7-C11, with the link to Phi in C9).
- **P2 (p. 3, the case v real).** It cites an unidentified "three-plus-three generalized-circle theorem", which is not a paper statement. I supplied a short argument: contact 1 forces all of {0..5}, which is a minor-only face (C36-C38).
- **P3 (p. 4, L > 0 step).** Q(v) != 0 is used tacitly. It holds by the factorization in C18.
- **P4 (p. 3).** The justification that the z, w denominators are nonzero is loose. The factorizations in C17 settle it.
- **P5 (p. 3).** The ratio of the contact differences is -t, not t. This is cosmetic.
