# Lifted contacts, n = 4 to n = 5

Local attack on MOC at n = 5 around block configurations U0 = U4 (+) 1 whose determinant d(U0) lies on the boundary of H = conv{z_sigma}.

## Deviations

- **No nontrivial n = 4 contact was found.** The four "vertex contacts" in `../lib/contact_n4.json` are edge contacts on degenerate spectra. Entry 0 has a_0 + b_3 = 0 (to 7e-16), and d sits on the edge from z_0123 to the zero vertex. Entries 1 to 3 have all 24 z_sigma collinear to about 1e-12 relative, so H is a segment. In both cases `moc_core.margin`'s `theta` and `active` come from float noise, which is where the "vertex" label came from. With the spectra held fixed and nondegenerate (`n4_contacts.py`), the largest m over U at distance at least delta from generalized permutations and block forms is -1.6e-3 (delta = 0.15, 160 spectra) and -6.7e-4 (delta = 0.05, 155 spectra). In every case the block-distance constraint is active. So the lifts analysed here are of trivial contacts (transposition-block edges, U4 = R_ij(phi) P_sigma) and of the degenerate collision contact.
- **Toolkit hull bug.** When many z_sigma nearly coincide (for example a_i + b_j = 0), `lib/moc_core.hull_indices` returns a wrong hull, and `margin` reports false positive m. It gave +5.2e-5 at n = 5; see `out/fake_hullbug_n5.json`, where the 60-digit value is -3.9e-3. Everything here uses the robust hull in `lift_lib.py`, which merges near-duplicate points and treats near-collinear ones as collinear. Every reported value was rechecked with `moc_hp.hp_margin` or the 50-digit `hp_lift.py`.

## Method

Take an off-block perturbation U = U0 expm(eps K_x) with K_x[i,4] = x_i and x in C^4. The map d is even in x, because J = diag(1,1,1,1,-1) commutes with D_a, D_b and U0. So the first variation is 0, the mixed second derivatives with in-block and spectral directions vanish, and d = d0 + eps^2 w(x) + O(eps^4). Here w(x) = v^T W v is a complex quadratic form on v = (Re x, Im x) in R^8.

The z_sigma do not depend on U. So, at second order, a lifted contact gives a violation exactly when lambda_max(Re(e^{-i theta} W)) > 0 for an extreme outward normal theta at d0. That normal is the edge normal at an edge contact, or one of the two adjacent edge normals at a vertex. Block-preserving directions keep m <= 0 by MOC(4). W is computed analytically (`lift_lib.W_matrix`), checked against finite differences (relative error 1.9e-8 to 4.9e-8 on 5 random cases), and checked at U4 = I against the closed form w = sum_i |x_i|^2 (z_{(i 4)} - z_id) (relative error 4e-16).

## Runs (all with `nice -n 10`, one process at a time, each under 10 minutes)

| command | what | result |
|---|---|---|
| `run_lift.py check` | analytic W vs finite differences | `out/check.json` |
| `run_lift.py edge 501 300` | 300 random n = 4 spectra, every transposition hull edge, c in {0.2, 0.5, 0.8}, 8 random (a5, b5) each: 35,208 lifts, 3,366 on the boundary of H5 | lambda_max/diam in [-0.565, -1.94e-4], none > 0 |
| `run_lift.py vertex 502` | the 4 entries of `contact_n4.json` x 400 (a5, b5) | collision contact: 26 boundary lifts, lambda_max <= -7.5e-4; collinear ones: no boundary lift |
| `n4_contacts.py 504 40 0.15`, `505 300 0.15`, `506 300 0.05` | nontrivial n = 4 contacts, fixed spectra, 3 L-BFGS starts each | max m = -3.7e-3, -1.6e-3, -6.7e-4 |
| `run_lift.py adv 508 1e-12` | Nelder-Mead (1500 evaluations) maximizing lambda_max over all 5 spectral pairs and c (12 block-edge seeds; nondegeneracy guard 0.05 x rms) and over (a5, b5) for the collision contact | sup lambda_max = -7.0e-11, reached only as a permutation that moves index 5 approaches the supporting line |
| `run_lift.py direct 503 150`, `direct 509 400 adv_s508.json,edge_s501.json,vertex_s502.json` | U0 expm(eps K) with eps in {0.02, 0.1, 0.3}, then toolkit `smooth_maximize` and `polish`, rescored with the robust hull | 550 runs, 272,615 objective evaluations, max float m 2.5e-13; the top 5 of each run are negative at 60 digits (down to -1.9e-21); all end at generalized permutations |
| `hp_lift.py 50` | 50-digit rebuild of 4 representative lifts | d is on the edge line to 1e-52; lambda_max = -1.94e-4, -1.1e-7, -1.98e-9, -7.0e-11; the nearest other permutation gap is -3.0e-5, -5.2e-8, -8.4e-10, -2.9e-13 |

`adv_s507_tol1e-9.json` is a superseded run with a loose boundary tolerance. Its "positive" lambda values (up to 1e-7) all sit at lifts that are strictly inside by exactly the tolerance, and they stay negative at 50 digits.

## Finding

At every lifted contact that was examined, the off-block second variation is strictly negative. It tends to 0 only at transitions where a permutation that moves index 5 reaches the supporting line, which opens a flat 2 x 2 mixing direction. Near such a transition lambda_max is roughly proportional to that permutation's gap. There is no evidence of a violation.
