# Fixed structured U, optimized spectra

Hold a structured unitary U fixed and maximize the normalized signed margin m of
d(U) = det(D_a + U D_b U*) against H = conv{z_sigma} over the spectra (a, b) and
the direction theta. Uses the shared toolkit in `../lib` (not modified).

Files
- `fixedu.py`: spectra-only smooth objective (log-sum-exp over permutations,
  analytic gradient, U frozen), annealed L-BFGS (`maximize_spectra`), a fixed-U
  active-set SLSQP polish (`polish_fixed`), guard penalties, hull width.
- `unitaries.py`: Fourier F_n, circulant F diag(phases) F*, twisted complex
  Hadamards (F_p x I) D (I x F_q), Sylvester H_8, O* (+) V rotated out of block
  form by plane rotations (V = I, F_{n-4}, O*), O* (x) R(t), products of plane
  rotations along a random cycle.
- `run.py`: driver; `hpcheck.py`: 50-digit re-evaluation (`moc_hp.hp_margin`,
  all n! points) of every saved record with float m > 1e-7, plus the best
  non-degenerate guarded record of each file.
- `runs/`: JSON outputs (full-precision a, b, U, theta, seeds) and logs.

Modes. `free`: no penalties. `guard`: relative eigenvalue gap >= 0.1, each of a
and b not collinear (sqrt(lambda2/lambda1) >= 0.2), z cloud not collinear
(>= 0.05), and from run r5 on, common offset |mean a + mean b| <= 5 rms. A record
counts as non-degenerate when the minimum width of H is >= 1e-3 diam(H).

Commands (from this folder, `python -I -B`, each under `nice -n 10` and a time limit)

    r1: run.py runs/r1_n56.json 501 540 fourier:5:1:30 circulant:5:4:10 cycle_big:5:6:8 cycle_small:5:4:8 cycle_short:5:3:8 ostar_I:5:4:8 ostar_cplx:5:4:8 fourier:6:1:30 circulant:6:3:8 twisted:6:3:8 cycle_big:6:4:8 ostar_F:6:4:8 ostar_I:6:3:8
    r2: run.py runs/r2_guard.json 502 540 fourier:5:1:6 fourier:6:1:5 fourier:7:1:3 fourier:8:1:2 fourier:9:1:2 ostar_I:5:2:3 ostar_cplx:5:2:3 ostar_F:6:2:2 ostar_O:8:1:2 ostar_kron:8:1:2 cycle_big:5:3:3 cycle_short:6:2:3 cycle_big:7:2:2 circulant:6:2:3 twisted:6:2:2 sylvester8:8:1:2
    r3: run.py runs/r3_n9_ostar_cycle.json 503 480 fourier:9:1:1 ostar_I:5:2:2 ostar_cplx:5:2:2 cycle_big:5:2:2 cycle_short:6:1:2 ostar_F:6:1:2 circulant:6:1:2 twisted:6:1:2 sylvester8:8:1:1 ostar_O:8:1:1
    hp: hpcheck.py runs/hpcheck.json runs/r2_guard.json runs/r3_n9_ostar_cycle.json
    r5: MOC_FREE_EVERY=0 run.py runs/r5_n9_fourier.json 505 150 fourier:9:1:2 ; hpcheck.py runs/hpcheck_r5.json runs/r5_n9_fourier.json

Code versions. r1 ran an earlier version of `fixedu.py`/`run.py`: no polish, no
z-cloud or offset guard, 7 annealing stages with maxiter 250, a free and a
guarded start at every restart; its margins (about -1e-5) are the smoothing floor,
and the polished versions of its top records are segment hulls (see below).
r2 and r3 ran without the offset guard. A run r4 (seed 504, Fourier n = 9, 3
guarded starts) was killed by its outer time limit before writing JSON; its
printed best non-degenerate float margin was -9.0e-2. It is not counted below.

Findings
- No violation. Every float m > 1e-7 (n = 8 Fourier, n = 8 Sylvester, n = 8
  O* (+) O*, n = 9 Fourier: float m from 1.4 to 25.5) is a cancellation fake: the
  optimizer moved a to a common offset of 1e7 to 1e9 (the constant-spectrum
  limit), and at 50 digits the same configurations give m between -0.35 and -0.11
  (`runs/hpcheck.json`).
- Free mode converges to segment hulls: spectra with all n! points z_sigma
  collinear to about 1e-9 although a and b are well separated and not collinear.
  There the float margin is +-1e-16 and the 50-digit value is about -1e-16 (the
  rounding of the float spectra). These are degenerate (m = 0) configurations,
  not violations.
- Guarded, non-degenerate best margins (float; the n = 6 value confirmed at 50
  digits): Fourier n = 5: -2.0e-3, n = 6: -3.7e-5, n = 7: -7.7e-4, n = 8: -2.2e-4;
  cycle_big n = 5: -6.1e-6; circulant n = 6: -1.1e-4; twisted Hadamard n = 6:
  -1.8e-4; O* (+) I n = 5: -3.4e-4; O* (+) F_2 n = 6: -7.9e-5; Sylvester n = 8:
  -8.4e-2; O* (+) O* n = 8: -0.10. At these optima a guard is active (relative
  gap 0.10 in a or b, and often fatness 0.20), so within each family the
  supremum is approached at the trivial sets (collisions, collinear spectra).
- Fourier n = 9 (r5, offset guard on, one guarded start of 231 s, full 362880-point
  hull): m = -0.1272 in float and -0.127184539061465825685981087722 at 50 digits.
- cycle_short n = 6 (a rotation cycle on 4 of the 6 indices, so U is block
  diagonal) reaches m = -1.1e-16 with 121 active permutations: the known
  block-diagonal zero-margin set.
