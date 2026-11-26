# Audit: final_referee_single_collision_divisor_compensation_20261023

Evidence level: EXACT for the Theorem 2 flow law (6 exact unitaries) and for the divisibility of Delta_O - sum mu Pi by every one of the 28 collision brackets in homogeneous coordinates (one rational O in SO(4)); PROOF-READ for the general arguments, with one argument supplied by the auditor; CITED: the shared-zero Birkhoff theorem (audited in this batch, confirmed-modulo-citations) and the max-flow min-cut theorem.

Document: `proofs/manin_transfer_20261021/final_referee_single_collision_divisor_compensation_20261023.pdf`, 3 pages, "Uniform collision identities along a single denominator divisor".
Replay: `replay/final_referee_single_collision_divisor_compensation_20261023__check.py`, output in `..._check_output.txt` (all pass). It reuses the own-code law construction in `replay/final_referee_shared_zero_birkhoff_20261022__check.py`.

## Deviations and limits

- Content matches the paper. Lemma 1 together with Corollary 3 is the supplement's version of the paper's Corollary B.3 (pp. 41-42). Theorem 2 is the paper's Lemma B.2 (pp. 40-41).
- The general-U statements are PROOF-READ. The exact replays are instances: 4 Gaussian-rational unitaries and 2 rational SO(4) matrices for Theorem 2, and one rational O for the 28 homogeneous divisibility checks.
- The "subharmonic pasting argument" that Lemma 1's conclusion feeds into is not part of this document. It is not audited.

## A. Statements

- **Lemma 1 (p.1).** Fix the matrix (U = O in SO(4) in the application). Suppose a reduced collision divisor {d = 0} carries a probability law p, independent of all other spectral coordinates, with D_U = sum p_sigma F_sigma identically on {d = 0}. Let a holomorphic spectral curve have normalizing denominator E = d e, with e a holomorphic unit at the puncture. Then (D_U - sum p F)/E is holomorphic there, whatever the vanishing order of d(t). Consequently, on the supported branch (all gaps g_sigma >= 0), the normalized Leibniz sum L is bounded below near the puncture. The supplement writes g = -eta Im(.) with eta real; this is the paper's Re(eta' .) with eta' = i eta.
- **Theorem 2 (p.2)** = paper Lemma B.2. For every U in U(4) there is a probability law on S_4, depending only on U (that is, on the eigencolumns p, q), that represents all 16 principal minors of U diag(alpha, beta, gamma, gamma) U^* for all complex alpha, beta, gamma. Hence det(diag(x) - U diag(alpha, beta, gamma, gamma) U^*) = sum p_sigma prod (x_i - e_sigma(i)) for all x. The same holds for any repeated pair (by relabeling) and for a repeated pair in the other spectrum (interchange plus inverse permutations).
- **Corollary 3 (p.2).** For a denominator that is a product of collision brackets, every puncture at which exactly one factor vanishes has locally bounded negative part on the supported branch, at any vanishing order. Internal pairs use Theorem 2; cross pairs use the shared-zero theorem. The document explicitly makes no claim when two factors vanish together.

No contact sets appear.

## B. Dependencies

| ref | kind | used for | match checked |
|---|---|---|---|
| final_referee_shared_zero_birkhoff_20261022, Theorem 1 | common lemma (other doc, in this batch) | uniform law on cross divisors [A_i, C_j] = 0 | yes: law independent of the remaining labels, any pair via relabeling; checked EXACTLY that Delta_O - sum mu Pi vanishes on all 16 cross divisors |
| max-flow min-cut (classical) | external | existence of theta | hypotheses checked; the cut computation is written out below |
| Rueckert Nullstellensatz / factor theorem | external (classical) | divisibility by the reduced bracket | see C, step 2 |

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | Lemma 1: D_U - sum p F = d R with R holomorphic (polynomial) | identity / division | PROOF-READ | In an affine chart where d = x - y is linear, divide by the monic linear factor. The remainder is the function evaluated at x = y, which is 0 by hypothesis for all other coordinates. |
| 2 | the same holds in homogeneous (primitive-spinor) coordinates, including punctures where a colliding point goes to infinity | identity | PROOF-READ (argument supplied by auditor) + EXACT instances | The document gives one sentence. Supplied argument, as in paper B.3: Delta_O - sum mu Pi is multihomogeneous and vanishes on the dense open part of the irreducible hypersurface {[v,w] = 0} where the chart applies, so it is divisible by the bracket. Replayed EXACTLY for one rational O on all 16 cross, 6 [C,C] and 6 [A,A] brackets, by substituting v = t w and getting 0. The affine-chart identity Delta_O = det(diag(a) - O diag(e) O^t) is also checked. |
| 3 | pull back and divide by E = d e: get R/e, holomorphic at any order of d(t) | analytic | PROOF-READ | needs e(p) != 0, which is a hypothesis |
| 4 | L - sum p g = -eta Im(R/e), using sum p = 1 | identity | EXACT (formal) | |
| 5 | on the supported branch g >= 0, so L >= -C locally | inequality | PROOF-READ | |
| 6 | Thm 2: total edge mass 1, incident mass a_i + b_i | identity | PROOF-READ + EXACT on 6 U | |
| 7 | for every S, sum_{i<j in S} w_ij = \|p_S\|^2\|q_S\|^2 - \|<p_S,q_S>\|^2 <= \|p_S\|^2 | Lagrange identity | PROOF-READ + EXACT (all 15 S, 6 U) | |
| 8 | flow criterion: theta exists | max-flow min-cut | PROOF-READ | Written out by the auditor: a cut with edge set E', vertex set S' has finite capacity only if E' lies inside S'; then the capacity is 1 - w(E') + a(S') >= 1. So max flow = 1 = total demand. An exact augmenting-path max flow returned value 1 on all instances. |
| 9 | incoming mass at i is b_i | identity | EXACT on instances, PROOF-READ | |
| 10 | rank-two expansion of the principal minors of gamma I + (alpha-gamma) P + (beta-gamma) Q | identity | PROOF-READ + EXACT (all 16 I, 6 U) | e_2 of the rank-two part equals (alpha-gamma)(beta-gamma) sum w |
| 11 | expectation of the site product under theta gives the same expansion; splitting ordered-pair mass gives a law on S_4 | identity | EXACT on instances, PROOF-READ | |
| 12 | expanding det(diag(x) - M) in principal minors gives the full determinant law | identity | EXACT (polynomial identity, 6 U) | |
| 13 | interchange of spectra with inverse permutations | identity | EXACT (6 [A,A] brackets via homogeneous check) | det(D_a - U D_c U^*) = det(D_c - U^* D_a U) in order 4 |
| 14 | Cor 3: case analysis "internal pair / cross pair" covers every bracket type | case split | PROOF-READ | brackets are [A,C], [A,A], [C,C] |
| 15 | "ambient identities, so they hold when other points collide, provided no second factor vanishes" | logic | PROOF-READ | |

## D. Falsification

- With mpmath at 50 digits, along curves a_0(t) = e_0 + t^k (k = 1, 2, 3) with random fixed other labels and rational O: \|(D - sum p F)/t^k\| is constant to 6 digits at t = 1e-4, 1e-8, 1e-12. Without the law, \|(D - F_tau)/t^k\| grows like t^{-k}. NUMERIC; it agrees with Lemma 1.

## E. Verdict

**confirmed-modulo-citations.** Citations: the shared-zero Birkhoff theorem (audited in this batch; no gap found) and the classical max-flow min-cut theorem. One step (2, the homogeneous/primitive-spinor extension) is stated in a single sentence. I supplied the argument (the same one the paper gives in Corollary B.3) and replayed it exactly on instances, so it is marked PROOF-READ rather than GAP. Lemma 1's hypotheses are stated loosely ("original response", "local spectral chart"). Read with the paper's B.3 conventions, they are exactly what the proof uses. The document correctly disclaims the two-factor boundary, which must be handled elsewhere.
