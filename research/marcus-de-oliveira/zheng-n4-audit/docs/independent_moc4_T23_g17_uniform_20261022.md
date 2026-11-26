# Audit: independent_moc4_T23_g17_uniform_20261022 (seven-contact set T_24 + {17})

Evidence level: EXACT for the source, the coefficient formulas, identity (2), identity (7), the p = e = 0 positivity, the capacity-transfer group element, and all 28 Bernstein numerators. EXACT-AT-POINTS (a randomized exact identity test) for the Hessian factorization (3). PROOF-READ for the necessity lemma, the region cover and the boundary extension. CITED for the capacities (4)-(6) and the double-SOC at ab = 1 (NUMERIC support).

PDF: proofs/manin_transfer_20261021/independent_moc4_T23_g17_uniform_20261022.pdf (4 pages). Replay: `replay/independent_moc4_T23_g17_uniform_20261022__check.py` (needs MOC4_SUPP; options --numeric, --bernstein).

## Deviations and limits

- The degree-8 polynomial G was loaded as DATA from the supplement file independent_moc4_T23_g17_full_radial_20261022.json. The variable correspondence was identified by matching its stored cubic and Q terms against -dW and Q, which we re-derived ourselves (exact match). G was then validated through identity (3) with our own 14x14 Hessian of the paper's cubic C, in exact rational arithmetic at 25 random rational points. This is a probabilistic identity test, not a symbolic proof of (3). We did not expand the full symbolic Hessian determinant.
- The Bernstein replay uses our own chart substitution and our own clearing of positive denominators. Our monomial counts therefore differ from the document's table in regions I+, II, IIIa+ and IIIb+, and agree in I-, IIIa- and IIIb-. All coefficients are strictly positive in all 28 numerators.
- Capacities (5), (6), SOC, old-H, K4 and the double-SOC at ab = 1 are statements about Leibniz sums proved elsewhere (coherent Gram certificates). They are CITED and only tested numerically at their thresholds.

## A. Statements

Theorem 1 (p. 1). For every x in X with P c(x) >= 0 vanishing on T_24 ∪ {17} = {2,3,4,6,7,8,17} (Table 5, seven-set 5; extra contacts allowed; all coefficients may vanish), Phi_{c(x)}(O) >= 0 for all O in SO(4). All spectra are included: the argument uses only the necessary equations C(y) = 0 and det Hess C(y) = 0 for the real part y, which hold on all of X, including degenerate strata.

## B. Dependencies

| dependency | used for | status |
|---|---|---|
| X subset Sing C (paper A.1) | C(Re x) = 0 and det Hess C(Re x) = 0 | replayed in lib; re-tested exactly on 8 Gaussian-rational points here |
| SOC pj >= d^2, old-H e(j-d) >= d^2, K4 >= 0 | thresholds 1/p, 1+1/e; K4 >= 0 is safe | CITED (supplement lemmas); NUMERIC at threshold |
| capacities (5) (ab > d^2) and (6) (ab < d^2, b > sqrt(ab)), transferred from the paper-T_25 source | thresholds 1/e - 1/L and 1/e + b(1-w)^2/(b-w)^2 | CITED (blind_moc4_T24_coherent_capacity and the T_25 documents); the transfer element g = ((0,3,1,2),(1,2,0,3),1) exists and maps facets (59,53,55,67,54,2) to (5,3,67,0,54,2) (EXACT). It maps T_25 to {2,3,4,7,8,16}, not into T_24+17, so the transfer is valid only because the capacities are realizability-free Leibniz-sum inequalities, invariant under (4.7) |
| double-SOC ej >= d^2 at ab = d^2 | ab = 1 layer | CITED |

## C. Steps

| # | claim | status |
|---|---|---|
| 1 | sources {0,2,3,5,53,54,67}, independent (rank 7), only f54 non-minor | EXACT |
| 2 | coefficients a=g1, b=g22, e=g12, p=g13, k=g16, j=g10, d=g9 = g5-g16 = g0-g1-g16 | EXACT |
| 3 | d = 0 gives a minor combination | PROOF-READ |
| 4 | C(Re x) = 0 and det Hess C(Re x) = 0 on X | PROOF-READ (polarization; the real case uses Euler H(y)y = 2 grad C(y) = 0) + EXACT tests |
| 5 | (2): C(g + h kappa) = -dW + hQ | EXACT |
| 6 | dW > 0, so Q != 0 and h = dW/Q | PROOF-READ |
| 7 | (3): Q^4 det Hess C(g + (dW/Q) kappa) = -4 d^2 Q^2 G^2, so G = 0 is necessary | EXACT-AT-POINTS (25 random rational points, exact arithmetic) |
| 8 | (7): k0 = d(b+d)W{(a+d+k)(-K4) + (a+d)(pW + aek)} | EXACT (symbolic, from the JSON G) |
| 9 | p = e = 0: all coefficients of G positive (coefficient of d^8 is 1) | EXACT |
| 10 | Bernstein criterion: k0 > 0 and nonnegative middle and endpoint coefficients give G > 0 on [0, J*) | PROOF-READ |
| 11 | 7 charts x 2 leaf orientations: middle and endpoint numerators have positive coefficients | EXACT (28 tables) |
| 12 | the charts cover every nonnegative (a,b,e,p) with K4 < 0 (ab != 1), and each J* is the threshold of a cited inequality | PROOF-READ (cuts e/(1+e) and eb/(e+b) checked; K4 < 0 iff e < L when ab > 1; K4 < 0 automatically when ab < 1) |
| 13 | boundary and limit cases (a = 0, e = 0, ab = 1, b = 1) by continuity of the Bernstein coefficients, with k0 > 0 from (7) | PROOF-READ (only nonstrict limits are needed besides k0 > 0) |
| 14 | thresholds imply safety | CITED |

## D. Falsification

- Capacities at their thresholds (d = 1, k = 0, 12 random draws, seed 5, 8 BFGS restarts): worst min Phi is 1.8e-14 (SOC), 2.7e-14 (old-H), 2.8e-14 (5), 1.1e-13 (6). No violation, and all four look tight.
- The Bernstein tables are the falsifiable core; they passed exactly.

## E. Verdict

confirmed-modulo-citations. Citations: the coherent Gram capacities (5) and (6), SOC, old-H, K4 and the double-SOC from supplement documents not audited here; paper A.1. Limit: identity (3) is checked by exact evaluation at 25 random rational points, not symbolically. It ties G (taken from supplement data) to the Hessian, and the whole Bernstein replay depends on it. A full symbolic expansion of the 14x14 Hessian determinant would close this. Presentation: the Gram sketches for (5) and (6) on p. 2 do not explain how a Gram configuration certifies the Leibniz-sum inequality, so they cannot be read stand-alone.
