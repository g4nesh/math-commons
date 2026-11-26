Evidence level: EXACT for every identity, generator list, decomposition, parametrization, Hermitian-form determinant and capacity certificate listed below (own sympy code, `replay/root_moc4_T4_uniform_20261022__check.py`); PROOF-READ for the sign and case logic; CITED for Theorem 1.2 / Corollary 2.2, Proposition 4.1, Theorem B.1 and Theorem C.1 of the paper; NUMERIC for the falsification runs.

# Audit of `root_moc4_T4_uniform_20261022` (primary proof of paper case T_4)

Document: `proofs/manin_transfer_20261021/root_moc4_T4_uniform_20261022.pdf`, 5 pages, "The contact family T4: Hermitian forms and a four-coefficient inequality". Auditor label C02-T4-T5.

## Deviations and limits

- The text extract was used for reading. The formulas I rely on were re-derived from the definitions, not copied, and every displayed identity I needed was replayed by my own code. The page images were not opened.
- The document justifies its key inequality (1) in two ways: by citing a "companion exact proof", and by saying that (1) "follows independently" from positive semidefiniteness of a 4x4 arrowhead matrix Z. It never says how the Leibniz sum is a quadratic form in Z. I did not reconstruct the T4 arrowhead link. A numerical search over simple minor-based vectors w(O) found none. Instead I proved (1) independently. I replayed exactly the quaternion Gram certificate of the companion note `final_referee_T5_four_budget_capacity_20261022`. I derived its link to the Leibniz sums myself (O(p,q) = L(p)R(q)). I then transported it to the T4 generators with an explicit group element. Step S4 is EXACT through this auditor-supplied route. The document's own route stays unexplained.
- The supplement's checker scripts were not run or used.
- Numerical sampling found no supported chart point with d > 0: 200,000 unbiased draws gave 63 supported points, and 200,000 draws in the region forced by the proof gave 3,067. So the d > 0 branch (the only branch where inequality (1) is needed) could not be falsification-tested on spectra. It may be empty or very thin. That does not affect the proof.

## A. Statements

Conventions. F_sigma = kappa * prod_i (a_i - c_{sigma(i)}) with c = -b. In the paper's notation this is alpha * prod_i (a_i + b_{sigma(i)}), with alpha = kappa. The gaps are g_sigma = Re kappa (F_sigma - F_0), so the reference permutation is pi* = pi_0 = identity, which lies in T4. Permutation indices are 0-based in lexicographic order. The facet numbering is the paper's Appendix D numbering. I checked this for every generator: f0 = n_{1,4}, f1 = n_{12,34}, f2 = n_{13,34}, f3 = n_{14,34}, f4 = n_{1,3}, f5 = n_{12,24}, f6 = n_{4,1}, f11 = n_{13,23}, f33 = n_{12,13}, f34 = h_{4321,+}.

**Theorem 1 (p. 1).** Let x be in X with c = c(x) for reference pi_0. Assume (Pc)_sigma >= 0 for all sigma, and (Pc)_sigma = 0 for sigma in T4 = {0, 1, 4, 5, 7, 8} (non-strict: extra contacts allowed). Then Phi_c(O) >= 0 for every O in SO(4), and so the support inequality holds for every U in U(4). The spectra may be arbitrary complex, including repeated values and a_i + b_j = 0.

- T4 equals the paper's T_4 exactly (the identity element of G; EXACT).
- The statement is (*) for T_4: contacts are non-strict, every O is covered, and all spectra are covered (degenerate strata go to cited theorems, see B).

Auxiliary claims used as lemmas:

- (1) For a, b, c, d >= 0 with ab(c - d) >= d^2 (a + b + c + d), Phi_{a f6 + b f0 + c f5 + d f34}(O) >= 0 for all O in SO(4).
- (2) On the 8-dimensional face space {Pc : (Pc)_T4 = 0}, g = g19 f0 + g16 f1 + g14 f3 + g10 f5 + g15 f6 + g6 f11 + g13 f33 + d f34, where d = g2 - g13.

## B. Dependencies

| dependency | kind | used for | status |
|---|---|---|---|
| Theorem 1.2 / Corollary 2.2 (two-real-orthogonal reduction) | core-doc (paper) | U(4) -> SO(4) | CITED |
| Proposition 4.1 (82 facets) and the f_j numbering | core-doc (paper) | generator list | CITED; the generator list itself is EXACT via lib `generators` |
| Theorem B.1 (a_i + b_j = 0) | common-lemma (paper App. B, audited by the project) | shared-value stratum | CITED. It gives MOC outright for such spectra (paper p. 41), so it covers the stratum exactly |
| Theorem C.1 via Proposition 1.4 (a spectrum with at most 3 distinct values lies on a generalized circle) | common-lemma (paper App. C) | repeated-value stratum | CITED. The hypotheses match: one whole spectrum is on a generalized circle |
| "three-plus-three generalized-circle theorem" | other-doc | the case v real | Not a result of the paper (C.1 needs a WHOLE spectrum on a circle). The auditor replaces it, see S9 |
| four-coefficient inequality (1) ("companion exact proof", transported by a symmetry) | other-doc: final_referee_T5_four_budget_capacity_20261022 | the d > 0 branch | Re-proved EXACTLY by the auditor (S4) |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| S1 | The generators of N_T4 are {f0, f1, f2, f3, f4, f5, f6, f11, f33, f34}; only f34 is non-minor; the minor labels are as stated | identity | EXACT | lib `generators`, `F_TABLE` |
| S2 | f34 = f0 - f2 + f6 + f33; f0 - f1 - f2 - f3 + f4 = 0; -f0 + f2 - f6 - f33 + f34 = 0 | identity | EXACT | 24-value vectors |
| S3 | Decomposition (2) holds on the face space (dim 8 = rank of generators) | identity | EXACT | checked on a basis |
| S4 | Inequality (1) | inequality | EXACT (route supplied by auditor) | (a) O(p,q) = L(p)R(q) parametrizes SO(4) (O^T O = \|p\|^2\|q\|^2 I replayed; surjectivity is classical). (b) On S^3 x S^3: Phi_f0 = (A-B)^2, Phi_f6 = (A+B)^2, Phi_f10 = D^2, Phi_f54 = H, as exact polynomial identities modulo \|p\|^2 = \|q\|^2 = 1. (c) The companion identity n(a(A-B)^2 + b(A+B)^2 + cD^2 + dH) = sum_j [l(S_j y)^T Z l(S_j y) + m(S_j y)^T Z m(S_j y)] mod N, replayed exactly. (d) Z >= 0 iff the capacity condition, by an exact Schur-complement chain. (e) Transport: g = ((3,2,1,0),(3,1,2,0),0) in G sends (f6, f0, f5, f34) to (f0, f6, f10, f54). g.O has det -1, and Phi_{c4}(O) = Phi_{g.c4}((g.O)J) with J = diag(1,1,1,-1) and (g.O)J in SO(4) (replayed on 8 exact rational O; the general fact is the paper's (4.7) plus the column-flip remark). The document's own arrowhead pencil: its Schur complement formula is EXACT, but its link to Phi is not given (see E, gap G1) |
| S5 | If d <= 0: adding A·(rel 2) + tau·(rel 1), with A = -d and tau = max(0, A - b), gives coefficients (b-A+tau, g16-tau, A-tau, g14-tau, tau) on f0..f4, g9 on f6, g2 on f33, all >= 0 by g22 = g16 + b - A and g18 = g14 + b - A | identity + inequality | EXACT (identities) + PROOF-READ (signs) | |
| S6 | Collision strata: repeated values -> C.1; a_i = c_j -> B.1 | citation | CITED | hypotheses match |
| S7 | Normalization A = (z, 0, 1, x), C = (w, v, inf, tv) for 8 distinct points; F_sigma = -prod_{sigma(i) != 2}(a_i - c_sigma(i)) | analytic | PROOF-READ | PGL_2 invariance of brackets (X is the cone over (P^1)^8 // SL_2) |
| S8 | Contact ratios at T4 are 0, 1-x, t-x, x(t-1), r-x+1, u; x, t real; kappa = i mu/[v(z-w)] | identity + argument | EXACT (ratios) + PROOF-READ (reality: from alpha(1-x), alpha(t-x), alpha x(t-1) in iR with t != 1 and x != 1) | |
| S9 | Solving the r, u equations gives the displayed z, w; the denominators are nonzero | identity | EXACT | The resultant in z is (z - v)·(a linear factor). The other solution of the bilinear system is z = w = v, a collision, so the chart is exhaustive. Den_z = 0 forces r v (tv - x) = 0, a collision |
| S9' | v real -> "three-plus-three theorem" | citation | GAP (presentation), repaired | Auditor's argument: v real gives z real (z·Den_z = r v (tv - x) with real coefficients and Den_z != 0). Then w is real, because r = z(1-tv)/v · (w-v)/(z-w) is a nondegenerate real Möbius function of w (z != v). So all eight points are real, and C.1 applies (indeed kappa is purely imaginary, so every gap vanishes and Phi = 0) |
| S10 | g2 = -rho x/\|v\|^2, g6 = -rho r t (x-1)/\|tv-1\|^2, c = g10 = rho t u x/\|tv-x\|^2, g9 = (x-1)c, a = g9 - d | identity | EXACT | gap = -rho times the divided difference (R(v) - R(vbar))/(v - vbar) |
| S11 | Signs c > 0, x > 1, rho < 0, rt > 0, tu < 0 (r = 0 or t = 0 is a collision) | inequality | PROOF-READ | |
| S12 | (3) k = hc - t(b+d), Ee = tyHb - Jd; (4) J = E + hy(r+1), J - xhH = -uy | identity | EXACT | |
| S13 | t < 0 impossible; so t > 0, r > 0, u < 0, and 0 < t < 1 | inequality | PROOF-READ | |
| S14 | E >= 0 implies K = tyb - xhd > 0 | inequality | PROOF-READ | |
| S15 | d\|Q(v)\|^2 and K\|Q(v)\|^2 are real Hermitian of bidegree (1,1); det M_D = -r^3 t^2 u x y L E, D(x/t,x/t) = r u x y (hr - Hy); det M_P = -r^2 t^2 u^2 x^2 h y L E, P(x/t,x/t) = r u^2 x^2 y (y - h); E = (h-y)(H+hy) + hyL | identity | EXACT | Also checked: Q(v) divides the numerator of z - w, so Q(v) != 0 on the chart. The document uses this implicitly |
| S16 | E < 0: L > 0, then y > h, then M_P is positive definite, so K > 0 | inequality | PROOF-READ | 2x2 symmetric matrix: det >= 0 and one negative value give NSD; det > 0 and one positive value give PD |
| S17 | ab(c-d) - d^2(a+b+c+d) = c{ybc - xd(b+d)} when a = yc - d; this is >= (b+d)cK/h > 0 | identity + inequality | EXACT + PROOF-READ | |
| S18 | Conclusion for U(4) | citation | CITED (Theorem 1.2) | |

## D. Falsification (NUMERIC)

- Chart sampling with fixed seeds 0 and 1, 200,000 draws each, spectra built from the document's chart and gaps computed directly from the sigma-points. Unbiased ranges gave 63 supported points; ranges restricted to x > 1, 0 < t < 1, r > 0, u < 0 gave 3,067. No supported point with d > 0 was found in either run, so the sign claims and the capacity could not be exercised on spectra. There were no contact or sign violations.
- Capacity inequality, on the boundary ab(c-d) = d^2(a+b+c+d) in the T5 form, which is equivalent by S4(e): 3 coefficient sets x 4 BFGS restarts. Minimum Phi = 5.7e-14 (no violation). As a sanity check, a = b = c = d = 1, which violates the capacity condition, gives min Phi = -0.0084, so the test can detect a failure. Output: `replay/final_referee_moc4_T5_uniform_20261022__check_output.txt`.
- All exact checks: 44 PASS lines and 0 FAIL in `replay/root_moc4_T4_uniform_20261022__check_output.txt` (run with `--numeric`, about 100 s).

## E. Verdict

**confirmed-modulo-citations.** Every identity was replayed exactly. Every inequality step was read and holds. The citations are Theorem 1.2 / Corollary 2.2, Proposition 4.1, Theorem B.1 and Theorem C.1 (paper).

Gaps (both presentation; both repaired by the auditor):

- G1 (p. 2, eq. (1) and the arrowhead Z): the claim that (1) "follows independently" from Z >= 0 is unsupported in the document, because no vector w(O) with Phi = w^T Z w is given. The inequality itself is true: the auditor proved it exactly through the companion quaternion certificate and an explicit group transport (S4).
- G2 (p. 2, "If v is real ... three-plus-three generalized-circle theorem"): no such theorem is in the paper, and C.1 needs a whole spectrum on a circle. Repair (S9'): when v is real, the chart forces z and w real, so all eight points are real and C.1 applies (all gaps vanish).
- Minor: Q(v) != 0, needed in the contradiction "D(v, vbar) = d\|Q(v)\|^2 > 0", is not stated. It holds because Q(v) divides the numerator of z - w (EXACT).
