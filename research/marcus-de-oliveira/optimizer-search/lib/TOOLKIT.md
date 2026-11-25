# MOC search toolkit

Shared library for optimizer-based searches for a counterexample to the Marcus–de Oliveira conjecture. Python 3 with numpy, scipy and mpmath; the exact parts use `fractions` only. Run everything from this folder. Tests: `python -I -B test_toolkit.py` (45 checks, all pass; see `test_output.txt`).

Conventions. Permutations are 0-based tuples in lexicographic order, with `sigma[i]` the image of `i`. `perms(n)[k]` is permutation `k`. `perm_matrix(sigma)` has `P[i, sigma[i]] = 1`, so `d(P) = z_sigma`. The margin is `m = signed distance(d, H) / diam(H)`, positive outside `H`.

## API, one example each

| module.function | what it does | example |
|---|---|---|
| `moc_core.perms(n)` | all n! permutations (n ≤ 8), int8 array | `perms(4).shape == (24, 4)` |
| `moc_core.z_all(a, b)` | all z_sigma, vectorized (n ≤ 8) | `z = z_all(a, b)` |
| `moc_core.z_chunks(a, b)` | lean enumeration for n = 8, 9: one chunk per sigma(0), (n−1)! rows each | `for S, z in z_chunks(a, b): ...` |
| `moc_core.hull_lean(a, b)` | hull vertices and their permutations without storing all n! points | `V, S = hull_lean(a9, b9)` |
| `moc_core.det_d(a, b, U)` | d(U) | `det_d(a, b, perm_matrix(s))` equals `z_sigma` |
| `moc_core.margin(a, b, U)` | dict: `m`, `theta`, `active` (permutation indices within 1e-9·diam of the support line), `d`, `diam`, `unitarity_defect`, `ok` (False if U is not unitary to 1e-10 or the hull has fewer than 3 vertices) | `margin(a, b, U)["m"]` |
| `moc_core.margin_points(d, z, n)` | margin from a point list; raises unless `len(z) == n!` | `margin_points(d, z_all(a, b), 4)` |
| `moc_core.dist_genperm(U)` / `dist_blockdiag(U)` | distance from generalized permutations and from block-diagonal form up to permutations (0 on those sets) | `dist_genperm(U)[0]` |
| `moc_core.min_separation(a, b)` | minimum eigenvalue gaps, in canonical units | `min(min_separation(a, b))` |
| `moc_param.cayley(K)`, `inverse_cayley`, `expm_unitary`, `log_unitary`, `skew_from_real(X)` | U from skew-Hermitian K; `skew_from_real` maps any real n×n matrix onto the skew-Hermitian K | `U = cayley(skew_from_real(X))` |
| `moc_param.cayley_coords(U)` | `(K, R, gap)` with `U·diag(R) = cayley(K)`, choosing R so that eigenvalues stay away from −1; d is unchanged | used by the high-precision and exact paths |
| `moc_param.canonical_spectra(a, b)` | fixes the shift (mean b = 0), the scale (joint rms 1) and the phase (a0 − a1 > 0); m is unchanged | `a, b = canonical_spectra(a, b)` |
| `moc_param.embed(a, b, U, a_x, b_x, U_x=None)` | direct sum, with U ⊕ I by default | lift an n = 4 configuration to n = 5 |
| `moc_param.random_unitary`, `random_spectra(n, rng, kind)` | Haar U; `kind` is gaussian, circle or mixed | |
| `moc_opt.smooth_maximize(a, b, U, **opts)` | annealed L-BFGS on `[Re(e^{-iθ}d) − LSE_τ(Re(e^{-iθ}z))]/s` with analytic gradient. θ is a free variable, so the max over directions is exact; τ runs from 0.1 to 1e-4 (relative to the rms spread s). U = U0·cayley(K), and U0 is recentred at each stage | `a, b, U, th, info = smooth_maximize(a, b, U, dperm=0.15)` |
| options `sep`, `dperm`, `dblock`, `mu` | quadratic penalties when the eigenvalue gap, `dist_genperm` or `dist_blockdiag` drops below the limit | |
| `moc_opt.polish(a, b, U, theta)` | active-set SLSQP: maximize t subject to `Re(e^{-iθ}(d − z_σ)) ≥ t·s0` on a growing working set, with spectra moving under s = s0 and mean b fixed. It is **unconstrained** by `opts`; always re-evaluate with `margin` | |
| `moc_opt.run_restarts(n, R, seed, out_json, **opts)` | restart r uses `default_rng([seed, r])`, then smooth, then polish (kept only if the constraints still hold at 0.9 of each limit). Writes a summary and the 20 best configurations as JSON | `run_restarts(5, 100, 1, "out.json", dperm=0.1)` |
| `moc_hp.hp_margin(a, b, U, dps=60)` | recomputes U as Cayley(K) at `dps` digits (exactly unitary), then d, all z, the hull and m in mpmath (n ≤ 8) | `hp_margin(a, b, U)["m_str"]` |
| `moc_exact.build_certificate(a, b, U, theta, den, dir_den, out_json)` | rational K over Q(i), exact Cayley U, Gaussian-rational spectra, exact d, all n! z, exact membership, and a rational direction (p, q) with its exact gap. The status is `outside-certified` only if gap > 0 and U is exactly unitary | |
| `check_cert.py cert.json` | independent checker (standard library, no shared code): exact unitarity, Bareiss determinant over Z[i], its own permutation enumeration (count = n!), strict gap equal to the claimed one. Exit 0 iff certified | `python3 check_cert.py c.json` |
| `moc_exact.witness_json` | input for the repo verifier `certificates/verify.py` (n ≤ 6, integers of at most 64 digits: use small `den`) | |

Pipeline for a candidate: run `margin` and require `m > 1e-7` and `ok`, then `hp_margin` at 60 digits must give the same m, then `build_certificate` and finally `check_cert.py`.

## Calibration (`calibrate.py STAGE`, outputs in `calibration/`)

- **n2, n3.** 2000 restarts each (seeds 101 and 102; smooth + polish). Max float m is 2.8e-13 and 7.0e-13; 0 of 4000 runs exceed 1e-7. The 5 best of each, recomputed at 60 digits: n = 2 gives |m| ≤ 8e-59 (the hull is a segment, so m = 0 exactly), and n = 3 gives values from −1.1e-16 to −1e-24.
- **n4.** 1000 restarts (seed 104). Max float m is 4.2e-14, and the 5 best are ≤ −8e-21 at 60 digits. The constrained run n4c (300 restarts, seed 105, dperm = dblock = 0.15) has max m 8.6e-15, and its 5 best are all negative at 60 digits.
- **psi.** ψ(O*) = −158/2401, exact. ψ takes the values {0, 1, 2} on permutation matrices. min ψ over U(4) is −0.06744224886812313 (200 BFGS restarts, seed 2026, all within 1e-6). At the minimizer, |U|² has the 0.3239/0.0283 pattern.
- **fakes.** (1) Non-unitary "U": float m = 0.57 but `ok` = False with defect 0.69; the exact checker rejects it ("U is not unitary"). (2) Wrong hull: with the 6 top permutations dropped, a true m of −0.138 reads as +0.124; the completeness guard raises, and a forged certificate is rejected by the checker (negative gap over all 24 permutations). (3) U is a hull-vertex permutation perturbed by 1e-9: float m = +1.7e-16, and at 60 digits m = −5.2e-18.
- **verifier.** On one n = 4 instance (rational Cayley U and spectra, denominators ≤ 12), the toolkit and `certificates/verify.py` give the same exact d, and both say `inside`.

## Pitfalls

- **Trivial contact is everywhere.** m = 0 occurs at permutation U and at block-diagonal U. It also occurs on the whole fibre {U : d(U) = z_σ} over a hull vertex, which has positive dimension: the constrained run found vertex contacts at `dist_genperm` ≈ 0.15–0.22 (see `contact_n4.json`). Other zero-margin cases: repeated eigenvalues, constant a or b, and segment hulls (n = 2 or collinear spectra). Near these, float m reaches about +1e-13. Treat any m < 1e-7 as zero.
- **Conditioning.** Cayley fails near the eigenvalue −1; use `cayley_coords`, which chooses the phases. The adjugate is computed by SVD, so a singular D_a + B is fine. Normalize spectra (`canonical_spectra`), because d and z scale as cⁿ.
- **Permutation completeness.** Never compute a hull from a subset. `margin_points` and the checker enforce n! points.
- **polish ignores the penalties.** Check the distances afterwards.
