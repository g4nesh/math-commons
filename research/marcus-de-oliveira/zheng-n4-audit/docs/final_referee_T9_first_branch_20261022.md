Evidence level: EXACT for every displayed identity: (1), (2), (3), the zeta forms, L, B0, the affine-in-e claim and the polynomial N, F(e*), R(e*), (6), the P(c*) expansion and the Theorem 2 identity (22 exact checks in `replay/final_referee_T9_first_branch_20261022__check.py`). PROOF-READ for the sign and convexity arguments. NUMERIC for the conclusions on region I.

# final_referee_T9_first_branch_20261022: region I of {0,1,3,5,7,8} (paper T_10)

PDF: `proofs/manin_transfer_20261021/final_referee_T9_first_branch_20261022.pdf`, 4 pages. This is a branch document, not a primary proof of a paper case.

## Deviations and limits

- The text extract drops the complex-conjugation bars. I read every formula I rely on from the page images. On p. 2 the "for completeness" forms carry conj(zeta), and with that reading they replay exactly.
- The chart and the main-document setup (A = (z0, 0, 1, x), C = (w0, ∞, v, 1 + t(v-1)), kappa) come from `final_referee_moc4_T9_uniform_20261022`. They are replayed in that document's check.
- The SOC inequality used at the end of Theorem 2 is CITED. Its source is unresolved; see the T9 main report.

## A. Statements

Setting: T = {0,1,3,5,7,8}, eight distinct finite points except C_1 = ∞, v nonreal, and region I: x = 1+X, t = -T, u = -U with X, T, r, U > 0. Every gap is >= 0 and d = g2 - g18 > 0. Gaps are normalized by g2 = 1.

- **Lemma 1 (p. 2).** Under (1), B0 < 0.
- **Theorem 1 (p. 3).** A0 h g19 - (TX - 1) d >= 0, that is x(1-t) g19 - (-t(x-1) - 1) d >= 0.
- **Theorem 2 (p. 4).** g6 g22 >= d^2. The SOC inequality on (f1, f11, f28) then gives Phi >= 0 for the T_10 cost in region I.

Coverage: this is the region-I part of the closed family. Collisions and real v are handled in the main document.

## B. Dependencies

- Main-document chart and identities: replayed EXACT in the T9 main check.
- SOC inequality (s f1 + b f11 + d f28 is safe when sb >= d^2): CITED, source unresolved. NUMERIC: no violation (T9 main report, section D).
- The decomposition (3) of the main document: EXACT there.

## C. Steps

| # | claim | status | note |
|---|---|---|---|
| 1 | (1): g14 = hn - d, g19 = G - n, g9 = A0(H - G), g18 = 1 - d, g6 = hG/T, g22 = [d + A0(H - n)]/X | EXACT | rebuilt from the 24 products; the normalized gap is the ratio of difference quotients in zeta |
| 2 | (2) formula for g20 | EXACT | |
| 3 | sign consequences: n >= d/h > 0, G >= n, H >= G, 0 < d <= 1; G, H > 0 | PROOF-READ | |
| 4 | (3): NL = M and 2 xi B0 = Ec - GN(Tn + H); also the conj(zeta) forms and K = -d + n conj(zeta) | EXACT | |
| 5 | L = Gd + Gc + e(G + Xn) > 0 and B0 = e(kn - c) - G(c + hn) | EXACT / PROOF-READ | |
| 6 | (4): F = -L B0 \|zeta - h\|^2 | PROOF-READ | Direct from (3): expand \|zeta-h\|^2 = N - 2h xi + h^2 |
| 7 | F is affine in e; dF/de = (T/h) N with the displayed N; F(e*) as displayed | EXACT | N is a sum of visibly nonnegative terms when 0 <= d <= 1, and its last term X d^2 A0 > 0 |
| 8 | Lemma 1 conclusion: B0 >= 0 would force e >= e*, so F(e) >= F(e*) > 0, contradicting (4) | PROOF-READ | |
| 9 | Theorem 1: the case k <= 0 is trivial; otherwise c < c* gives q > 0 and e < e*; R(0) > 0 and R(e*) > 0 | EXACT / PROOF-READ | |
| 10 | (6): L B0 F20 = -T R P(c) | EXACT | polynomial identity in (G, n, d, e, X, T) after N -> M/L and 2 xi L B0 -> EL |
| 11 | p2 > 0, P(0) = p0 < 0, P(c*) <= 0 by the displayed expansion with Q0, Q1, Q2 >= 0; convexity gives P < 0 on [0, c*) | EXACT / PROOF-READ | The expansion was replayed with T = (k+1)/X and W = Xh. P does not involve G, so fixing e, n, d is legitimate |
| 12 | Contradiction F20 < 0, so g20 < 0 | PROOF-READ | |
| 13 | Theorem 2 chain, including (d + hc)(d + A0 c) = TX d^2 + (c+d){A0 hc - (TX-1)d} | EXACT / PROOF-READ | hG >= d + hc uses a >= 0; the second factor uses H >= G |

## D. Falsification

400000 samples of region I (seed 2026, log-uniform parameters, nonreal v); 14098 were supporting with d > 0. The minimum of the left side of (5) was 1.2e-2 and the minimum of g6 g22 - d^2 was 2.7e-2. No violation.

## E. Verdict

**confirmed-modulo-citations.** I found no gap. The only external input is the SOC inequality, whose source is unresolved.
