# Optimizer search for counterexamples to the Marcus–de Oliveira conjecture, n = 5 to 9

**Result in one line.** About 8,500 local optimizations at n = 5 to 9, plus 5,300 calibration restarts at n ≤ 4, found no configuration with positive margin that survived a 50- or 60-digit recheck. Every configuration found with margin 0 is one of these: a generalized permutation U, a block-reducible U, or a degenerate spectrum. The degenerate spectra are concyclic a and −b, an eigenvalue collision, or a zero factor a_i + b_j = 0 together with a concyclic sub-configuration. No candidate reached the exact-certificate stage.

**Caveats, read first.**

- This is numerical evidence. It is not a proof, and it does not decide the conjecture for any n ≥ 5.
- Coverage is uneven. Most compute went to n = 5. n = 6 and 7 got about 1,800 restarts between them. n = 8 and 9 were covered only by a few fixed structured unitaries (one to three starts per U).
- The shared library `lib/moc_core.py` has a known defect, which was not fixed. When many z_σ nearly coincide, its hull routine can return a wrong hull, and then `margin` reports false positive values or the wrong sign. The lanes that worked near such configurations used their own robust hull (`lifted/lift_lib.py`, `hermitian/robust.py`). The tables distinguish floating-point values from the selected 50- and 60-digit rechecks using the full n! point set.
- This synthesis disagrees with two interpretations in the lane reports. (1) The hermitian lane's "edge contacts away from the trivial sets" are approaches to degenerate spectra (Section 6.3). (2) The structured lane's "possible new structured degeneracy" for the Fourier matrix F_5 is the known concyclic family. Both re-readings come from `synthesis/degeneracy_check.py`, which only reads saved run files.

## 1. The question

Let A and B be normal n × n complex matrices with eigenvalues a_1, …, a_n and b_1, …, b_n. Write D_a = diag(a) and D_b = diag(b). For U unitary, set

  d(U) = det(D_a + U D_b U*),  z_σ = ∏_i (a_i + b_σ(i)) for σ ∈ S_n,  H = conv{z_σ}.

The Marcus–de Oliveira conjecture (1973) states that d(U) ∈ H for every U. It is known for n ≤ 3, for Hermitian pairs (Fiedler), and for several structured families. A proof for n = 4 has been claimed (Zheng, October 2026) and is being audited in a separate lane of this repository. For n ≥ 5 it is open.

By Cauchy–Binet,

  d(U) = Σ_{k} Σ_{|S|=|T|=k} |det U[S,T]|² ∏_{i∉S} a_i ∏_{j∈T} b_j,

so d is linear in the vector of squared minors of U. For n ≥ 4 that vector can leave the convex hull of its values at permutation matrices (Drury 1992). An explicit witness is the rational 4 × 4 orthogonal matrix O*, at which one facet functional ψ of that polytope takes the value −158/2401. Over U(4), the minimum of ψ is about −0.0674422. The conjecture survives only if no functional of the product form above, coming from some (a, b), detects this excess. The searches here look for such an (a, b, U) directly.

**Margin.** The normalized signed margin is m = (signed distance from d to H) / diam(H), positive outside H. It is computed as m = max_θ [Re(e^{−iθ} d) − max_σ Re(e^{−iθ} z_σ)] / diam(H). m does not change under the invariances

- (a, b) → (a + t, b − t),
- (a, b) → (c a, c b),
- complex conjugation,
- swapping a with b together with U → U*,
- U → L U R with L, R diagonal unitaries,
- simultaneous permutation of indices.

**What would count as a violation.** All three of the following:

1. Float m > 1e-7.
2. The same sign and size after rebuilding U as an exact Cayley transform and recomputing at 50 or more digits.
3. An exact rational certificate: Gaussian-rational spectra, an exactly unitary U over Q(i), all n! points z_σ, and a rational direction with a strictly positive gap. The certificate is checked by `lib/check_cert.py`, which shares no code with the rest of the library, and, for n ≤ 6, by the repository verifier `certificates/verify.py`.

**Known zero-margin sets.** m = 0 happens in all of the following cases, none of which is a violation:

- U is a generalized permutation; then d = z_σ.
- U is block-reducible after permutations; then MOC reduces to smaller n, because d = d₁d₂ and conv(Z₁)·conv(Z₂) ⊂ conv(Z₁Z₂) by bilinearity.
- Some eigenvalues are repeated, or a or b is constant.
- H is a segment.

Near these sets, floating point produces spurious m up to about 1e-13. It also produces much larger fakes, listed in Section 5.

## 2. Method

**Shared toolkit** (`lib/`, documented in `lib/TOOLKIT.md`; 45 self-tests pass, see `lib/test_output.txt`).

- *Smooth objective.* `moc_opt.smooth_maximize` maximizes [Re(e^{−iθ}d) − LSE_τ(Re(e^{−iθ}z_σ))]/s over (a, b, U, θ) with L-BFGS and an analytic gradient. Here s is the rms spread of the z. The direction θ is a free variable, so the maximum over directions is exact; only the maximum over permutations is smoothed. τ is annealed from 1e-1 to 1e-4 relative to s. U = U₀·cayley(K) is recentred at each stage.
- *Penalties.* Optional penalties keep the eigenvalue gaps (`sep`), the distance to generalized permutations (`dperm`) and the distance to block form (`dblock`) above given bounds.
- *Polish.* `moc_opt.polish` is an active-set SLSQP step that removes the smoothing bias.
- *High precision.* `moc_hp.hp_margin` recomputes m at 60 digits on an exactly unitary Cayley(K).
- *Exact path.* `moc_exact.build_certificate` builds the exact certificate and `check_cert.py` checks it independently.

**Search lanes.** Open the lane report or search driver for the family you want to inspect.

| lane | n | idea |
|---|---|---|
| [joint5](joint5/README.md) | 5 | Joint ascent over (a, b, U) from eight start families: Haar, O* ⊕ 1 rotated out of block form, circle, line, cluster, roots of unity, nearly Hermitian, mixed. First unconstrained, then with sep = dperm = dblock up to 0.5 and a bound keeping d away from every z_σ. |
| [joint67](joint67/joint67.py) | 6, 7 | Joint ascent from six start families: Haar, O*/7 ⊕ Haar, Fourier with phases, U₂ ⊗ U₃, rotated direct sums, points on transposition edges. Free runs and runs with dperm = dblock = 0.15. A second-order test along transposition-edge curves. |
| [lifted](lifted/README.md) | 4 → 5 | Lifts of n = 4 boundary configurations U₄ ⊕ 1. Exact second variation in the off-block directions, adversarial maximization of its top eigenvalue over the spectra, then direct ascent from the lifts. Also a search for nontrivial n = 4 contacts at fixed nondegenerate spectra. |
| [structured](structured/README.md) | 5–9 | U held fixed and structured; spectra and θ optimized. U ranges over Fourier F₅…F₉, circulants, twisted complex Hadamards, Sylvester H₈, O* ⊕ {I, F, O*} rotated, and plane-rotation cycles. Runs are free or guarded (eigenvalue gap, non-collinearity, bounded common offset). |
| [hermitian](hermitian/hsearch.py) | 5, 6, 7 | Joint ascent restricted to spectral classes: a real with b having k = 4 to 7 nonreal entries; line/circle; two circles; unitary pairs; three-valued a; one outlier. Each class unconstrained, with `sep`, and with dperm = dblock = 0.15. A τ → 1e-6 continuation of the constrained optima. |
| [synthesis](synthesis/degeneracy_check.py) | 5–7 | No optimization. Degeneracy diagnostics of the saved near-contacts, and a numerical check of the reduction lemma in Section 6.2. |

All runs used `nice -n 10`, one process at a time, each run at most 10 minutes. Restart r of a run with seed S uses `numpy.random.default_rng([S, r])`. The structured lane uses `default_rng([seed, family_id, n, instance, restart, guard])`.

## 3. Calibration (`lib/calibration/`)

| check | result |
|---|---|
| n = 2, 3: 2000 restarts each (seeds 101, 102) | max float m 2.8e-13 and 7.0e-13; none of 4000 above 1e-7. Top 5 at 60 digits: n = 2, \|m\| ≤ 8e-59 (segment hull); n = 3, −1.1e-16 to −1e-24 |
| n = 4: 1000 restarts (seed 104) | max float m 4.2e-14; top 5 at 60 digits ≤ −8e-21 |
| n = 4 constrained, dperm = dblock = 0.15: 300 restarts (seed 105) | max float m 8.6e-15; top 5 negative at 60 digits. This run predates the fix that rejects polished points violating the penalties, so some of its records are trivial |
| ψ(O*) | −158/2401 exactly; ψ ∈ {0, 1, 2} on the 24 permutation matrices; min over U(4) = −0.06744224886812313 (200 restarts, seed 2026, all within 1e-6) |
| planted fakes | (1) A non-unitary U reads m = 0.57; it is flagged `ok = False` and the checker rejects it. (2) Dropping 6 permutations turns m = −0.138 into +0.124; the completeness guard raises and the forged certificate is rejected. (3) A permutation perturbed by 1e-9 reads +1.7e-16 in float and −5.2e-18 at 60 digits |
| repo verifier | on one n = 4 inside instance, the exact d agrees with `certificates/verify.py`, and both say "inside" |

The exact certificate path has never produced an "outside" certificate, because no such configuration was ever found. Its accept branch has been exercised only through arithmetic consistency.

## 4. Results by dimension and strategy

"Best 50/60-digit m" is the largest margin after the high-precision recheck. "Held off trivial U" means the best margin with U kept at distance at least 0.15 (sometimes 0.05 to 0.5) from generalized permutations and block forms; Section 6.3 explains why those values sit near 0.

| n | lane | optimizations | objective evals | max float m | best 50/60-digit m | held off trivial U |
|---|---|---|---|---|---|---|
| 5 | joint5, unconstrained (seed 501) | 4000 (8 × 500) | 1.97e6 | 9.13e-5, **fake** (−1.65e-11 at 60 digits); otherwise ≤ 2.4e-11 | all rechecked values between −1.6e-11 and −4e-17 | — |
| 5 | joint5, constrained (seeds 502–504) | 369 | 1.35e6 | 9.1e-16 | ≈ −1e-16 (contact) | m ≈ 0 at every bound up to 0.5 (with dv 0.2), reached only through degenerate spectra (Section 6) |
| 5 | lifted, direct ascent (seeds 503, 509) | 550 | 272,615 | 2.5e-13 | −1.9e-21, at a generalized permutation | — |
| 5 | lifted, second variation (seeds 501, 502, 508) | 36,808 lifts, 3,366 + 26 on ∂H | — | — | λ_max/diam ≤ −7.0e-11 (50 digits) | — |
| 5 | hermitian (seeds 501–532) | 1508 | 1.91e6 | +9.8e-6, **fake** (sign flip; −9.8e-6 at 60 digits) | trivial contacts only | see 6.3 |
| 5 | structured, guarded | part of 266 | — | — | F₅ −2.0e-3; cycle −6.1e-6 (float); O* ⊕ I −3.4e-4 | guard active at every optimum |
| 6 | joint67 free (seed 601) | 948 | 435,663 | 4.18e-15 | −8.6e-22, at a generalized permutation | — |
| 6 | joint67 constrained (seed 602) | 72 | 327,339 | 3.5e-16 | ≈ −1e-16 (contact) | 10 contacts, all degenerate; others ≈ −1.3e-5 (smoothing floor) |
| 6 | hermitian (seeds 601–631) | 606 | 0.75e6 | 3.2e-15 | negative | −4.8e-6 to −2.3e-5 at τ = 1e-4 |
| 6 | structured, guarded | part of 266 | — | — | F₆ −3.689e-5 (50 digits); circulant −1.1e-4; twisted Hadamard −1.8e-4; O* ⊕ F₂ −7.9e-5 | guard active |
| 7 | joint67 free (seed 701) | 85 | 41,985 | 1.65e-15 | −1.4e-20 to −7.5e-18, at generalized permutations | — |
| 7 | joint67 constrained (seed 702) | 15 | 70,209 | −5.5e-16 | −8.6e-16 (concyclic contact) | best other −9.243e-6 (60 digits), drifting toward concyclic (residual 5e-4) |
| 7 | hermitian (seeds 701–731) | 105 | 0.15e6 | 2.5e-15 | negative | ≈ −1.3e-5 to −2.1e-5 |
| 7 | structured, guarded | 1–3 starts per U | — | — | F₇ −7.7e-4 | guard active |
| 8 | structured | 1–2 starts per U | — | fakes up to 4.39 (cancellation) | F₈ −2.2e-4; Sylvester H₈ −8.4e-2; O* ⊕ O* −0.10 | guard active |
| 9 | structured | 1 counted guarded start, plus fakes | — | fake 25.5 (−0.279 at 50 digits) | F₉ −0.127184539061465825685981087722 (50 digits) | — |

Totals:

| lane | n | optimizations | objective evals |
|---|---|---|---|
| joint5 | 5 | 4369 | 3.32e6 |
| joint67 | 6, 7 | 1120 | 0.88e6 |
| lifted | 5 | 550 direct ascents, plus 36,808 analysed lifts and 337 fixed-spectra n = 4 contact searches (3 starts each) | 0.27e6 |
| structured | 5–9 | 266 | not tallied |
| hermitian | 5, 6, 7 | 2219 | 2.80e6 |

Most unconstrained runs ended at a generalized permutation. The counts by lane were

- joint5: 3970 of 4000 runs; 27 others end at repeated eigenvalues and 3 in the interior.
- joint67: all 948 runs at n = 6 and all 85 at n = 7.
- lifted: all 550 direct ascents.

This includes the starts at O* ⊕ 1, at Fourier matrices and at tensor products, where the squared-minor vector lies outside the permutation polytope.

**Fixed nondegenerate spectra show a strict gap.** In the lifted lane at n = 4, the spectra were held fixed and nondegenerate and U was kept at distance at least δ from generalized permutations and block forms. The best m was −6.7e-4 at δ = 0.05 (155 spectra) and −1.6e-3 at δ = 0.15 (160 spectra), with the block-distance constraint active every time. In the structured lane, every guarded optimum sits on a guard. Both say the same thing: when the spectra cannot degenerate, the margin returns to 0 only as U approaches the trivial sets.

## 5. Fakes that were rejected

Every float m > 1e-7 in the whole effort was a fake. All of them were caught by the high-precision recheck.

| float m | n | where | cause | true m |
|---|---|---|---|---|
| 25.5, 7.65 | 9 | structured, F₉ | a drifted to a common offset of about 1e9, so the float determinant cancels | −0.279, −0.304 |
| 4.39, 1.94, 1.40 | 8 | structured, H₈, O* ⊕ O*, F₈ | same cancellation (offsets 4e7 to 1e9) | −0.113, −0.277, −0.351 |
| 9.13e-5 | 5 | joint5 u_cluster, restart 142 | a and b nearly constant, diam H = 7.6e-27 | −1.65e-11 |
| 5.2e-5 | 5 | lifted, `out/fake_hullbug_n5.json` | wrong hull from `moc_core.hull_indices` when 24 z_σ coincide near 0 | −3.92e-3 |
| 9.84e-6, 9.00e-6 | 5 | hermitian a3c, a3r (repeated eigenvalues) | `moc_core.margin` flips the sign when z points nearly coincide | −9.84e-6, −9.00e-6 |
| up to 2e-9 | 5 | hermitian, no separation penalty | eigenvalues merging to 1e-9 | −9.3e-11 to −5e-13 |

Two lessons follow. First, the shared hull routine must not be trusted near coincident z_σ; replace it with the robust hull in `lifted/lift_lib.py`. Second, the spectra should be re-normalized (`canonical_spectra`) at every evaluation, not only at the start.

## 6. Where the margin reaches 0: contact families

### 6.1 Taxonomy

Every configuration found with m = 0 (to the precision of its spectra) falls into one of the following families. C1, C2 and C4 are the known trivial sets. C3 is a known true case. C5 is identified in this synthesis.

- **C1. Generalized permutation U.** Here d = z_σ, a hull vertex. This is the endpoint of essentially every unconstrained run. For U = P_σ e^{εK}, d = z_σ + ε² Σ_{p<q} |K_pq|² (z_τ − z_σ) + O(ε³), where τ = τ(p, q) runs over the transposition neighbours of σ. This was checked numerically to O(ε). So when every transposition neighbour is strictly behind the support line, the second variation is negative definite, and the maximum is locally isolated modulo diagonal phases. The "vertex contacts away from permutations" in `lib/contact_n4.json` are not counterexamples to this. The lifted lane showed they are degenerate spectra: entry 0 has a₀ + b₃ = 0, and entries 1 to 3 have all 24 z_σ collinear to 1e-12.

- **C2. Block-reducible U.** For example U = R_ij(φ)P_σ (a plane rotation times a permutation). Here d runs along the transposition edge [z_σ, z_{σ∘(ij)}].
  - joint67 checked 4,713 edge points at n = 6 (165 non-concyclic spectra, seed 611) and 1,281 at n = 7 (38 spectra, seed 711). The gradient is at most 1.5e-15, and the top Hessian eigenvalue is at most 9.3e-11. That value is finite-difference noise: it scales like 1/h.
  - The lifted lane computed the exact second variation off the block at n = 5 for 3,366 lifted edge contacts and the collision contact. It is strictly negative, with λ_max/diam between −0.565 and −7.0e-11, and the extreme value was confirmed at 50 digits.
  - λ_max → 0 only at *walls* in spectrum space. A wall is a codimension-one set where a further permutation that moves the new index reaches the supporting line. There, a flat 2 × 2 mixing direction opens, and λ_max is roughly proportional to that permutation's gap.

- **C3. Segment hulls: a and −b concyclic, or collinear.** The 2n points {a_i} ∪ {−b_j} lie on one circle or line. Then all z_σ have the same argument mod π, H is a segment, and d lies in it. This is the Hermitian case transported by a Möbius map. With φ(x) = (x − i)/(x + i) and A = φ(H₁), −B = φ(H₂), one has φ(H₁) − φ(H₂) = 2i(H₁ + i)^{−1}(H₁ − H₂)(H₂ + i)^{−1}. So d and every z_σ are the same complex constant times their Hermitian counterparts. C3 is the main sink of the constrained searches:
  - joint5: 87 of the 148 saved nontrivial-looking contacts. Their median concyclicity residual is 3.6e-16.
  - joint67: 9 contacts at n = 6 and 1 at n = 7, with residuals of order 1e-16.
  - structured: the Fourier F₅ "segment hulls". The unpolished optima are already within 3e-5 to 6e-5 of concyclic.
  - hermitian: the unitary-class optimum, with r = 1 + 7e-8 and both circles centred at 0, so concyclic to 3.6e-8. Also the real-a optima where Im b → 0, which approach the Hermitian case itself.

- **C4. Eigenvalue collisions and nearly constant spectra.** Many z_σ coincide, so the fibre of U over a hull point grows. This produced most hermitian-lane constrained optima (min gap 6e-8 to 7e-5) and the 27 repeated-eigenvalue endpoints of joint5.

- **C5. Nested zero-factor contacts (new here).** Some a_i + b_j = 0, so (n−1)! of the z_σ sit at 0, and 0 is a hull vertex. In addition, a complementary sub-configuration is concyclic. Then d lies on a hull edge through 0 at a U far from block form: in joint5's example, dist_genperm is 0.49 and dist_blockdiag is 0.38. Section 6.2 gives the structure and a proof for the simplest type.

### 6.2 The zero-factor family reduces to C3 one dimension down

joint5 reported 61 contacts with a_i + b_j ≈ 0 at U far from block form, and left open whether they are exact. joint67 found one more at n = 6. `synthesis/degeneracy_check.py` re-reads these configurations. In each, it takes the nonzero z_σ that lie on the supporting line and looks for a row set I and column set K such that σ(I) = K for all of them. Results for the 52 joint5 contacts with a single zero pair:

| type | count | nonzero z on the line | structure |
|---|---|---|---|
| 1 × 1 | 29 | 24 | the nonzero on-line z_σ are exactly {σ : σ(i) = k}; the complementary spectra {a_l : l ≠ i} ∪ {−b_l : l ≠ k} are concyclic to 1e-16; U has three zeros, either in row i (18 cases) or in column j (11 cases) |
| 2 × 2 | 22 | 12 | the nonzero on-line z_σ are exactly {σ : σ(I) = K} with \|I\| = \|K\| = 2; the complementary 3 + 3 spectra are concyclic to 1e-16; U has 5 zeros |
| unexplained | 1 | 14 | — |

The 9 contacts with two zero pairs include 2 more of the 2 × 2 type and 7 not explained by this test. The n = 6 contact (restart 37 of `joint67/runs/n6_constr_s602.json`) is of the 2 × 2 type: 48 nonzero z on the line, complement concyclic to 2.1e-16.

**Lemma (zero-factor reduction).** Suppose a_i + b_j = 0 and row i of U vanishes outside columns j and k. Then

  d(U) = |U_ik|² · (a_i + b_k) · d′,  d′ = det(D_{a′} + W D_{b′} W*),

where a′ = (a_l)_{l≠i}, b′ = (b_l)_{l≠k}, and W ∈ U(n−1) is determined by U. Hence, whenever MOC holds for (a′, b′) in dimension n − 1,

  d(U) ∈ conv({0} ∪ {z_σ : σ(i) = k}) ⊂ H.

*Proof.*

1. Shift (a, b) → (a − a_i, b + a_i), which changes neither d nor any z_σ. Now a_i = 0 = b_j.
2. With u_l = U e_l, write U D_b U* = b_k u_k u_k* + B₀, where B₀ = Σ_{l≠j,k} b_l u_l u_l*. The l = j term is absent because b_j = 0.
3. Since U_il = 0 for l ∉ {j, k}, the matrix X = D_a + B₀ has zero row i and zero column i. Hence adj X = det(X^{(i)}) e_i e_iᵀ, where X^{(i)} is X with row and column i deleted.
4. The matrix determinant lemma gives d = det X + b_k u_k* (adj X) u_k = b_k |U_ik|² det X^{(i)}.
5. After coordinate i is deleted, the vectors u_l (l ≠ j, k) are still orthonormal in C^{n−1}. Complete them with one unit vector placed in slot j, where the eigenvalue is b_j = 0. This gives X^{(i)} = D_{a′} + W D_{b′} W*.
6. Finally, b_k ∏_{l≠i}(a_l + b′_{τ(l)}) = z_σ, where σ(i) = k and σ = τ elsewhere; and z_σ = 0 whenever σ(i) = j. Since d = |U_ik|²·(b_k d′) + (1 − |U_ik|²)·0, the convex-hull statement follows. ∎

The column version follows by swapping a with b and U with U*.

If (a′, −b′) is concyclic, d′ lies on the line of its z′ by C3. Then every such U puts d on the line through 0 and the points {z_σ : σ(i) = k}. When that line supports H, this gives an exact contact set of positive dimension, far from block form, that is still covered by MOC in dimension n − 1.

Numerically, the identity holds on all 29 joint5 contacts of 1 × 1 type: the median relative error is 8.7e-12 and the maximum is 1.4e-8. The maximum is limited by treating |U_il|² < 1e-12 as exactly zero.

This also explains joint5's 80-digit Newton result. With a₃ + b₀ made exact, the maximum stalled at −5.7e-17·diam. That run left the complementary concyclicity at its float residual (2.6e-16 for that example, which is of 2 × 2 type). The 2 × 2 type is consistent with a two-step version of the lemma, but that is **not proved here**.

### 6.3 The constrained optima at −1e-5 are not a gap, and the "edge contacts" are degenerate

Several lanes report that, with U held at distance 0.15 from the trivial sets, the best m is about −1e-5 and the block distance sits on its bound. Two facts change how this should be read.

- **−1e-5 is the smoothing floor.** At τ = 1e-4 × (rms spread), a log-sum-exp over two equal values exceeds their maximum by τ log 2. The smoothed optimum therefore stops at m ≈ −τ log 2 · s/diam, which is of order 1e-5. Polish could remove this bias, but polish ignores the penalties. When the polished point leaves the constraint set it is discarded, and the smoothed point is kept. The hermitian lane's continuation to τ = 1e-6 shows m shrinking linearly in τ, for example −3.17e-7, −9.31e-8, …, −3.10e-9. This is the signature of the bias, not of a gap.
- **The points being approached are degenerate.** `synthesis/degeneracy_check.json` measures the 12 continuation endpoints of the hermitian lane. Every one is within 1e-8 to 1e-3 of a degenerate locus:
  - The unitary-class "edge contact" (seed 532, restart 10) has a and −b concyclic to 3.6e-8 and hull width/diam 6.6e-8. This is C3.
  - The real-a contacts (seed 503, restart 27, and seed 531, restart 7) have max |Im b| equal to 1.7e-6 and 1.4e-4. These are nearly Hermitian.
  - The others have eigenvalue gaps of 6e-8 to 7e-5 (C4) or min |a_i + b_j| of 2e-6 to 1e-3 (C5).
  - At n = 7, joint67's best constrained point (m = −9.243e-6) and the next three are drifting toward concyclic, with residuals 5e-4 to 1e-2.

So the hermitian lane's conclusion that "the contact set seems to reach out from block forms at least that far" is not supported. What reaches out is the set of degenerate spectra, where U is unconstrained anyway. The cleanest test of a genuine gap is the fixed-spectra one in Section 4, and it shows a gap that grows linearly with the distance from the trivial sets.

### 6.4 What the contact structure suggests

- **For a proof.** In every configuration examined here, the margin reaches 0 away from block-reducible U only on degenerate spectra. The degenerate loci are collisions, zero factors, and concyclic or nested-concyclic configurations, all of real codimension at least 2 when n ≥ 3. This suggests the following rigidity statement, which is a stronger, testable conjecture:

  > For spectra off a degenerate locus, d(U) ∈ ∂H only if U is block-reducible after permutations.

  Together with two further ingredients, this would give MOC by induction:
  1. block reduction: d = d₁d₂, and conv·conv ⊂ conv of products;
  2. a continuity argument in the spectra: the complement of a codimension-2 set is connected, and contacts at block U are strict local maxima off the block.

  The obstacle is the codimension-one walls of C2, where the second variation off the block vanishes. There, third- and fourth-order terms decide. A related analytic handle: at a critical point of Re(e^{−iθ}d) on U(n), the matrix e^{−iθ}[B, adj(A + B)] must be Hermitian, where B = U D_b U*. When A + B is invertible, this is the same as e^{i(arg d − θ)}(A + B)^{−1}[A, B](A + B)^{−1} being Hermitian. A contact at a non-block U would have to solve this equation with d on the support line.

- **For a disproof.** In this finite search, generic joint ascent did not produce a counterexample. Its observed near-contacts fell into C1 to C5. Two places are left uncovered:
  1. *The walls of C2.* These are spectra where three or more z_σ lie on a supporting edge and λ_max = 0, so the second variation does not decide.
  2. *High-contact strata.* The contact reduction (Theorem 3.6) in Zheng's n = 4 paper moves a counterexample at a fixed orthogonal matrix to one with at least six (= 2n − 2) permutation products on the supporting line. The repository's audit of that paper reads the argument as extending verbatim to n = 5, with at least 8 contacts, and to a fixed unitary U. That reading is proof-read only; it has not been verified. The searches here reached such strata only in degenerate form: segment hulls have all n! points on the line, and C5 has (n−1)! + 12 or 24.

## 7. Verification outcome

No configuration had m > 1e-7 after the high-precision recheck. So no exact certificate was built, `check_cert.py` was not run on an outside candidate, and the repository verifier was not invoked on any candidate. **certified violation: none.**

## 8. What this shows and what it does not

It shows the following:

- Across about 8,500 local optimizations at n = 5 to 9, many seeded at unitaries with squared-minor excess (O* embeddings, Fourier, Hadamard, tensor products), the margin never exceeded 0 beyond rounding.
- Every zero-margin configuration found belongs to a family where MOC holds for a known or elementary reason: C1 to C4, or C5 via the lemma above.
- With spectra fixed and nondegenerate, the margin stays strictly negative away from the trivial U.
- Block contacts are strict local maxima in the off-block directions, except on codimension-one walls.

It does not show the following:

- That MOC holds for any n ≥ 5. Local ascent from finitely many starts can miss a small basin.
- Anything about the walls of C2 or the high-contact strata, which were not targeted.
- Anything substantial at n = 8 and 9 (a handful of fixed U only) or at n = 7 (about 200 joint restarts and few constrained ones).
- Anything uniform in the size of the excess, or about orthogonal-versus-unitary questions.
- That the 2 × 2 type of C5 is a contact for an exact reason. That part is a numerical pattern.

## 9. Suggested next steps

1. **Fix the measuring instrument.**
   - Replace `moc_core.hull_indices` with the robust hull.
   - Re-normalize the spectra at every evaluation.
   - Add hard non-degeneracy constraints, so that no unguarded escape to C3 to C5 remains: hull width/diam ≥ w, concyclicity residual of every sub-configuration ≥ ρ, min |a_i + b_j| ≥ ζ, gaps ≥ δ.
   - Report the best m as a function of (w, ρ, ζ, δ) at τ → 0, with a constrained polish that moves the spectra.
2. **Search the high-contact strata directly.** At n = 5, parametrize spectra with at least 8 z_σ on a common supporting line and all others strictly behind it. Minimize the gap functional over U(5). Use an exact LP to discard contact sets whose gap vector lies in the cone of squared-minor vectors, as for the 124 easy orbits at n = 4.
3. **Analyse the walls of C2 to third and fourth order.** Maximize the resulting coefficient over the spectra.
4. **Proof side.**
   - Prove the 2 × 2 version of the zero-factor lemma and classify the degenerate loci.
   - Test the rigidity conjecture of Section 6.4 at n = 4. There it can be compared with the 90 contact cases of Zheng's paper.
   - Study the critical-point equation of Section 6.4 at non-block U.

## 10. Reproduction

Use Python 3 with NumPy, SciPy, SymPy and mpmath. For a bounded replay, run `python -I -B test_toolkit.py` from [lib/](lib/) for 45 checks. From [joint5/](joint5/), run `python -I -B hpcheck.py runs/cA_0.4.json 4 60`. All four archived positive floating-point margins become negative at 60 digits.

The longer examples below run from their named folders. Use a scratch checkout because several commands overwrite archived outputs. The joint search files store the command and seed in their summaries. Structured runs record their arguments in `argv`; the [structured report](structured/README.md) documents historical code differences. The Hermitian runs record seeds and options in their summaries, with batch arguments in the [list files](hermitian/). Time-limited searches can complete different numbers of restarts on different machines.

```
# toolkit (lib/)
python test_toolkit.py                          # 45 checks
calibration_dir=$(mktemp -d)
python -I -B calibrate.py fakes "$calibration_dir"
python -I -B calibrate.py verifier "$calibration_dir" "$calibration_dir" ../../certificates/verify.py
# Other stages include n2, n3, n4, n4c, and psi. The hp stage reads saved calibration runs.

# joint5/  (seed 501 unconstrained; 502-504 constrained)
python joint5.py --family haar --restarts 500 --seed 501 --time 140 --out runs/u_haar.json
python joint5.py --family mixed --restarts 300 --seed 502 --time 100 --sep 0.4 --dperm 0.4 --dblock 0.4 --out runs/cA_0.4.json
python joint5.py --family ostar --restarts 300 --seed 503 --time 75 --sep 0.3 --dperm 0.3 --dblock 0.3 --dv 0.1 --out runs/cB_ostar_0.3_dv0.1.json
python hpcheck.py runs/cA_0.4.json 4 60
python refine_hp.py runs/cA_0.4.json 3 80 5      # SNAP=1 makes a_3 + b_0 = 0 exact
python aggregate.py; python classify_contacts.py

# joint67/
python joint67.py 6 generic,ostar,dft,tensor,dsum,edge 3000 601 runs/n6_free_s601.json --sep 0.02 --time 420
python joint67.py 6 generic,ostar,dft,tensor,dsum,edge 3000 602 runs/n6_constr_s602.json --dperm 0.15 --dblock 0.15 --sep 0.05 --mu 10000 --time 300
python joint67.py 7 generic,ostar,dft,tensor,dsum,edge 3000 701 runs/n7_free_s701.json --sep 0.02 --time 330
python joint67.py 7 generic,ostar,dft,tensor,dsum,edge 3000 702 runs/n7_constr_s702.json --dperm 0.15 --dblock 0.15 --sep 0.05 --mu 10000 --time 360
python edgehess.py 6 1000 611 runs/n6_edgehess_s611.json --time 200
python edgehess.py 7 1000 711 runs/n7_edgehess_s711.json --time 220
python hpcheck.py hp_checked_configs.json runs/n6_free_s601.json runs/n6_constr_s602.json runs/n7_free_s701.json runs/n7_constr_s702.json

# lifted/
python run_lift.py check
python run_lift.py edge 501 300
python run_lift.py vertex 502
python n4_contacts.py 504 40 0.15; python n4_contacts.py 505 300 0.15; python n4_contacts.py 506 300 0.05
python run_lift.py adv 508 1e-12
python run_lift.py direct 503 150
python run_lift.py direct 509 400 adv_s508.json,edge_s501.json,vertex_s502.json
python hp_lift.py 50

# structured/  (full argument lists in structured/README.md)
MOC_FREE_EVERY=0 python run.py runs/r5_n9_fourier.json 505 150 fourier:9:1:2
python hpcheck.py runs/hpcheck.json runs/r2_guard.json runs/r3_n9_ostar_cycle.json

# hermitian/
sh batch.sh python list_n5.txt     # likewise list_n5b.txt, list_n6.txt, list_n7.txt, list_cs.txt
python hprobust.py runs/n5_herm_k4_c.json 2
python recheck.py
python cont.py runs/n5_herm_k4_c.json 2 runs/cont_n5_herm_k4_c.json

# synthesis/  (reads saved files only, a few seconds)
python degeneracy_check.py           # writes degeneracy_check.json
```

## 11. Files

- `lib/`: shared toolkit, calibration outputs, and `contact_n4.json`. Its four numerical "vertex contacts" are degenerate spectra; see C1.
- `joint5/`, `joint67/`, `lifted/`, `structured/`, `hermitian/`: lane code, run files (summaries, per-restart rows, best configurations at full precision) and high-precision checks.
- `synthesis/degeneracy_check.py`, `synthesis/degeneracy_check.json`: the degeneracy diagnostics of Sections 6.2 and 6.3.
