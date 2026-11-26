# Audit: independent_coreE_uniform_20261023 ("core E", seven-contact set E)

Evidence level: EXACT for every identity, generator decomposition, group action, the determinant identity (9), the polynomial G and all eleven positivity tables (own code, rebuilt from the paper's definitions); PROOF-READ for the realizability necessity and the case coverage; CITED for the three reference inequalities (sharp SOC, four-coefficient, chain); NUMERIC for falsification tests.

PDF: proofs/manin_transfer_20261021/independent_coreE_uniform_20261023.pdf (5 pages). Replay: `replay/independent_coreE_uniform_20261023__check.py` (output in `__check_output.txt`, 71 PASS, 0 FAIL).

## Deviations and limits

- The determinant identity (9) was replayed at d = 1 (seven variables q,p,t,b,c,v,z) and re-homogenized. This is legitimate because both sides are homogeneous of degree 22 in (q,p,t,b,c,v,d) once z = -C/Q is substituted, and dehomogenization is injective on homogeneous polynomials of a fixed degree. The derived G is then homogeneous of degree 8.
- G was derived independently, as the exact square root of Q^4 det H(w(-C/Q)) / (-4 Q^2). It coincides with the document's 408 coefficients, which were read from its JSON as data for the comparison only. All proof checks use the derived G.
- The chart maps were transcribed from the PDF text, not from the supplement's data. The cleared tables have exactly the term counts the document states (11763, 84446, 4181, 11409, 1091, 23837, 93285, 6427, 27531, 23799, 3955), which corroborates the transcription.
- The three reference inequalities are CITED. Their source documents were not audited. Only the group actions that carry them to the literal families used here were replayed.
- The one-sentence unitary remark in the statement (via the two-real-orthogonal theorem) is not needed for (*). It is CITED to the paper's Theorem 1.2 and was not checked here.

## A. Statements

Notation: E = {0,1,2,9,10,16,17}. This is literally seven-set no. 17 of the paper's Table 5. Its generators are generators(E) = {f0,f2,f3,f12,f13,f15,f26} with f0 = n_{1,4}, f2 = n_{13,34}, f3 = n_{14,34}, f12 = n_{14,24}, f13 = n_{2,1}, f15 = n_{12,14} (minor vectors) and f26 = h_{4132,-} (non-minor).

- **Theorem 1 (p. 1).** For every x in X (every spectra a, b in C^4 including repeated eigenvalues and a_i + b_j = 0, every alpha) and every reference pi*: if P c(x) >= 0 and (P c(x))_s = 0 for all s in E, then Phi_{c(x)}(O) >= 0 for every O in SO(4). Contacts are non-strict, so additional contacts are included. This is (*) for the seven-set E. The argument is pointwise and algebraic (no curve, perturbation or density argument), so all spectra are covered.
- **Eq. (1).** Every such c(x) equals q f0 + p f2 + t f3 + b f12 + c f15 + v f13 + d f26, uniquely. The seven coefficients are the gaps at permutations [23,15,14,8,4,7,3], so all are >= 0.
- **Eqs. (2)-(6) (sufficient tests, pointwise, no realizability).**
  - (2) x f_i + y f_j + d f26 has Phi >= 0 on SO(4) when xy >= d^2, for {i,j} in {2,12,15}.
  - (3) q f0 + v f13 + t f3 + d f26 has Phi >= 0 when qv(t-d) >= d^2(q+v+t+d).
  - (4) q f0 + t f3 + c f15 + d f26 has Phi >= 0 when c >= d^2/q + d^2/t, and likewise with v f13 in place of q f0.
  - (5) c >= d^2/(q+v) + d^2/t is sufficient.
  - (6) p(b+c) >= d^2, b(p+c) >= d^2 and c(p+b) >= d^2 are each sufficient.
- **Lemma 1 (p. 2).** Let c, d > 0, v >= 0, q(c+v) > d^2 and p, b >= 0. If t >= R_d(q,v,c) = d(v+d){q(c+v)+dv} / ((c+v){q(c+v)-d^2}), then the full source (1) has Phi >= 0 on SO(4). The same holds with q and v interchanged. The lemma is pointwise and does not use realizability.
- **Eq. (8)-(9).**
  - Along w(z) = c + z kappa, C(w(z)) = Cst + Q z, where Cst = d{(d+q)(d+v)+dt} + cbp and Q = d^2 + dt + cb - ct + cp + bp.
  - Q^4 det H(w(-C/Q)) = -4 d^2 Q^2 G^2, where G has degree 8 and 408 terms.
  - Hence every realizable source with d > 0 has G = 0.
- **Eq. (10).** G is invariant under q <-> v and under p <-> b.

## B. Dependencies

| dependency | kind | status | match |
|---|---|---|---|
| X ⊂ Sing C, the cubic (A.2) | paper | EXACT in moc4lib (grad C vanishes identically on matching vectors in 16 variables) | exactly what the polarization argument needs |
| sharp SOC reference (f1, f9; f22): x f1 + y f9 + d f22 safe if xy >= d^2 | external doc (supplement "sharp SOC", foundation F03) | CITED | ordered group actions to (f2,f12;f26), (f2,f15;f26), (f12,f15;f26) EXACT |
| four-coefficient reference (f0, f6, f10; f54), with f10 in the t role | external doc | CITED | ordered action to (f0, f13, f3; f26) EXACT |
| chain reference (f2; f0, f12; f70), with f2 central | external doc | CITED | ordered actions to (f15; f0, f3; f26) and (f15; f13, f3; f26) EXACT |
| independent_coreE_four_chain_capacity_20261023 | other doc (same batch) | Lemma 1 is reproved inside this document. The companion's content is also audited in this batch | n/a |
| Theorem 1.2 (two real orthogonal) | paper | CITED, only for the unitary remark | not needed for (*) |

No paper case, seven-set, branch, Appendix B or Appendix C result is invoked. The proof has no puncture or curve argument.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | generators(E) are the 7 listed, f26 the only non-minor; rank of P at E is 7; values at rows [23,15,14,8,4,7,3] form I_7, so the coefficients are gaps and >= 0 | identity | EXACT | replay |
| 2 | d = 0 gives a sum of minor vectors, so Phi >= 0 | identity | EXACT | (4.2) of the paper (lib-replayed), step 1 |
| 3 | (2)-(4) hold | citation | CITED | source documents not audited. Ordered group actions EXACT. Numeric boundary tests pass and detect violations at 2d |
| 4 | (5) and (6) follow from (2) and (4) by splitting the non-minor weight | identity | EXACT | sympy |
| 5 | Lemma 1: the allocation identities, W > 0, X > 0, W - X = WB/K, delta = dv/(c+v) in (0,d), total equals R_d, and the v = 0 case equals the chain | identity, inequality | EXACT, with the positivity remarks PROOF-READ | sympy |
| 6 | A realizable w = Re x satisfies C(w) = 0 and det H(w) = 0. Write w = x/2 + x̄/2 with both halves in Sing C. Cubic polarization gives C(w) = 0 and H(w) Im x = 0. If Im x = 0 then H(w)w = 2 grad C(w) = 0. Either way H(w) has a kernel vector (w is nonzero when d > 0) | analytic-algebraic | PROOF-READ | argument checked line by line. It uses only X ⊂ Sing C and conjugation invariance of X |
| 7 | (8) C(w(z)) = Cst + Q z; Cst > 0 when d > 0, so Q != 0 and z = -Cst/Q | identity | EXACT | replay with the paper's cubic in paper coordinates. The document's Cst and Q agree exactly |
| 8 | (9) Q^4 det H(w(-C/Q)) = -4 d^2 Q^2 G^2 with the 408-term G | identity | EXACT | own fraction-free Bareiss over Z[q,p,t,b,c,v,z] (819 exact divisions); exact division and exact square root |
| 9 | (10) symmetries of G | identity | EXACT | replay |
| 10 | Case c = 0. G = d A {(d+q)(d+v) - pb} B (EXACT). Failure of (2) gives pb < d^2, so the middle factor is > 0. If qv > d^2, failure of (3) is t < T_F (EXACT), and (qv - d^2) B(T_F) = ... >= 0 (EXACT), so B > 0 and G > 0, a contradiction | identity plus case split | EXACT and PROOF-READ | sympy and flint |
| 11 | Main partition: c >= 1 goes to dominant-c; c < 1 with p, b < 1 goes to row L; otherwise p >= 1 after the p <-> b swap goes to row P. Row P is onto {p >= 1, p(b+c) < 1, c > 0} | case split | PROOF-READ, inverse map EXACT | the swap is legitimate because G and the set of failed tests are both p <-> b invariant |
| 12 | Auxiliary sectors. If q, v >= 1 (after q <-> v), chart (11) with t < R from Lemma 1. Otherwise v < 1 and three rows by the sign of q(c+v) - 1, with the middle row from Lemma 1 (v = 0 allowed) | case split | PROOF-READ | q(c+v) > 1 holds in (11) because c > 0. All denominators are positive |
| 13 | Dominant c. Map (13) is onto {c >= 1, c(p+b) < 1}; (q,v) -> (w,V); rows of (14) by the sign of cw - 1, with the middle row from failure of (5) (t = 0 is T = 0) | case split | PROOF-READ, inverse map EXACT | replay |
| 14 | The eleven cleared tables have all coefficients > 0 and a strictly positive marker; reduced denominators have nonnegative coefficients and a positive monomial in the strict parameters | inequality certificate | EXACT | own flint computation. Term counts equal the document's |
| 15 | Conclusion: an unsafe realizable source with d > 0 and c > 0 lies in some chart, where G > 0, contradicting step 8 | logic | PROOF-READ | |

## D. Falsification (NUMERIC)

Seed 2026, multistart BFGS over SO(4) with 8 restarts. Values are min Phi / max|Pc|.

- At equality of each reference inequality, in its literal image here, with 40 draws per family:
  - SOC (3 families): worst 1.7e-16.
  - Four-coefficient: worst -2.0e-18.
  - Chain (both): worst 2.6e-16.
  - Lemma 1 at t = R_1 with random p, b in [0, 0.3]: worst -7.1e-17.
- Sensitivity check: with d doubled, the same tests give clearly negative minima (SOC -1.5e-3, four -9.6e-4, chain -3.4e-5). So the test can detect a violation.
- Realizable points with contact set exactly E and d > 0 (stratum search with margin 0.01, 2 seeds): G / sum|terms| = 2.5e-14 and 3.2e-15, so G = 0 as step 8 predicts. Min Phi / max|Pc| = 1.6e-12 and 7.9e-12.
- No violation was found.

## E. Verdict

**confirmed-modulo-citations.** Every step is EXACT or PROOF-READ except the three reference inequalities (sharp SOC (f1,f9;f22), four-coefficient (f0,f6,f10;f54), chain (f2;f0,f12;f70)). These are CITED. Their ordered literal images were verified exactly and they pass numeric boundary tests. No gap was found. Presentation notes, not gaps: equation numbers (10) and (12) are each used twice (the symmetry display and the main-row table; the auxiliary maps and the dominant-c map). The necessity argument of step 6 is terse but complete.
