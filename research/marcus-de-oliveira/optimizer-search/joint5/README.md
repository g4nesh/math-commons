# joint5: joint optimization of the MOC margin at n = 5

Maximizes the normalized signed margin m over (a, b, U) at n = 5 with the shared toolkit
(`../lib`, unchanged): `smooth_maximize` (annealed L-BFGS) then `polish`, from many starts.
`joint5.py` adds start families and one extra penalty (`--dv`: |d - z_sigma| >= dv * rms spread of
the z, for every sigma). All runs used `nice -n 10`, one process at a time, each under 10 minutes.

## Result

No violation. 4369 restarts, 3.32 million objective evaluations. No float m above 1e-7 except one
fake (below). Every m that the float code put above 0 is at most 2.4e-11 (u_circle) and is a
trivial contact.

| run (seed) | bounds sep = dperm = dblock, dv | restarts | max float m | where the maximizers went |
|---|---|---|---|---|
| u_haar, u_ostar, u_circle, u_line, u_cluster, u_roots, u_herm, u_mixed (501) | none | 8 x 500 | 2.4e-11 (excluding the fake) | 3970 at generalized permutations (m = 0, d = z_sigma), 27 with repeated eigenvalues, 3 interior |
| cA_0.05 ... cA_0.4 (502) | 0.05, 0.1, 0.2, 0.3, 0.4; dv 0 | 60, 58, 61, 64, 65 | 2e-16 to 9e-16 | sliver hulls and zero-factor edge contacts (see below) |
| cB_ostar (503) | 0.3, dv 0.1 | 33 | 5.7e-16 | same |
| cB_mixed (504) | 0.5, dv 0.2 | 28 | 8.4e-16 | same |

Best m as a function of the bounds: 0 (float 1e-16) at every bound tried, up to sep = dperm = dblock = 0.5
with dv = 0.2. These bounds do not keep the optimizer away from m = 0, because the m = 0 set also
contains the two families below.

**Fake rejected.** u_cluster, restart 142: float m = 9.13e-5. a and b are each nearly constant (all
entries equal to 1e-9) and diam H = 7.6e-27. hp_margin at 60 digits gives m = -1.65e-11. This is the
repeated-eigenvalue pitfall.

## m = 0 configurations that are not permutation or block type (`contacts/classification.json`)

148 saved configurations have float m within 1e-15 of 0, U far from generalized permutations and from
block-diagonal form (both distances at least 0.38 at bound 0.4), separated spectra, and d not at any z_sigma.
All 148 fall into two classes:

1. **Sliver hulls (87).** All 120 z_sigma lie within 1e-10 * diam of one line: the spectra are nearly
   collinear after the shift invariance, which is the degenerate segment-hull case. The toolkit's
   `ok` flag does not catch this because the floating-point hull still has more than 2 vertices.
2. **Zero-factor edge contacts (61).** Some a_i + b_j is about 0 (one pair in 52 cases, two pairs in 9),
   so 24 (or 42) of the z_sigma lie at 0. Then 0 is a hull vertex and d lies inside a hull edge
   [0, z_rho], with U not block-reducible. Example: cA_0.4, best index 3 (seed 502, restart 15):
   d sits at t = 0.32 along the edge, dist_genperm 0.49, dist_blockdiag 0.38.
   `refine_hp.py` holds the spectra fixed and maximizes Re(e^{-i theta} d) over U by Newton's method at 80 digits.
   With the float spectra the maximum converges to F = -7.45e-17 * diam (`contacts/edge_cA_0.4_3.json`).
   After setting a_3 + b_0 = 0 exactly (`SNAP=1`), it converges to F = -5.70e-17 * diam
   (`contacts/edge_cA_0.4_3_snap.json`, U to 80 digits). The gradient stalls at about 1.7e-16
   along directions where the Hessian is close to singular. So exact contact is **not** established:
   either a further exact spectral condition is needed or this is a near miss at the 1e-17 level.
   This class is not on the toolkit's list of trivial sets.

A joint (U and spectra) version, `refine_joint_hp.py`, was also run on this example. Newton
diverged because the maximum is not smooth in the spectra: splitting a_i + b_j away from 0 moves
one of the coincident z outward in both directions. Its output was discarded.

## Files

- `joint5.py`: driver with start families haar, ostar (O* plus a 1 x 1 block, rotated by expm(eps K) with eps from 0.05 to 0.8), circle, line, cluster, roots, herm and mixed. Restart r uses `default_rng([seed, r])`. Each `runs/*.json` holds a summary (with the exact command), one row per restart, and the 40 best configurations at full double precision.
- `hpcheck.py RUN.json K DPS`: runs hp_margin on the K best configurations.
- `refine_hp.py`, `refine_joint_hp.py`: high-precision Newton tests described above.
- `aggregate.py` writes `summary.json`. `classify_contacts.py` writes `contacts/classification.json`.

## Not covered

Constrained runs are slow (1.5 to 2.5 s per restart), so only 369 constrained restarts were done. No bound was placed on
hull width or on min |a_i + b_j|. That bound is the natural next constraint, and with it the best m
as a function of the bounds becomes a meaningful quantity. No exact certificate was attempted,
because no candidate exceeded 1e-7 after the high-precision check. n = 6 and 7 were not run.
