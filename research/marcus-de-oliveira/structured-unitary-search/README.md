# Structured unitary families, fixed U, n = 5..8

Task `mdo-structured-unitary-search`. Question: do structured unitary families give a Marcus–de Oliveira hull violation, or a new tight contact, at n = 5..8 when U is held fixed and the spectra and direction are optimized?

**Status: negative, no violation found.** Across 24 families (about 36 CPU-minutes), no float margin exceeded −1e-9, so no record needed the 50-digit recheck and nothing was certified. Every best gap approaches 0 from below only by sending the spectra to a degenerate boundary: two eigenvalues colliding (a_i = a_j or b_i = b_j) or a vanishing factor a_i + b_j = 0. When the spectra are kept separated (all gaps ≥ 0.1·rms), the separation constraint binds at every optimum, and the gap stays between −8e-7 and −0.37 depending on the family. No tight contact with non-degenerate spectra was found. This is numerical evidence only. It proves nothing about MDO, and the restart counts are low (2 to 8 per family, see the table), so it is not exhaustive.

## Method

`structured_search.py` holds U = U(p) inside one family and maximizes

  h(a, b, θ, p) = [Re(e^{−iθ} d) − max_σ Re(e^{−iθ} z_σ)] / ρⁿ,

with d = det(D_a + U D_b U*), z_σ = ∏(a_i + b_σ(i)) over all n! permutations, and ρ the joint rms spread of the spectra. h is invariant under (a, b) → (c(a+t), c(b−t)) for c > 0, and h > 0 for some θ iff d lies outside the hull. The max is a log-sum-exp at relative temperature 5e-2, 1e-2, 2e-3, 4e-4, then 1e-4. LSE ≥ max, so each stage's smoothed value is a lower bound. The search uses L-BFGS-B (≤ 250 iterations per stage) with an analytic gradient in (a, b, θ), checked against finite differences by `--selftest`, and central differences in the family parameters p. Restart r of family f starts from `default_rng([2026, f, r])`: random spectra (gaussian, circle or mixed) and uniform phases p. Even restarts are free. Odd restarts add a quadratic penalty that keeps every |a_i − a_j|, |b_i − b_j| and |a_i + b_j| ≥ 0.1·ρ, and they are reported only if all three are ≥ 0.09·ρ at the end. The reported m is the sign-correct edge margin `lift_lib.robust_margin` (signed distance over hull diameter). Any record with m > −1e-9 (up to 3 per family and mode) is set to be recomputed at 50 digits with `moc_hp.hp_margin`, and none qualified.

## Families

O* = (1/7)[[1,4,4,4],[−4,1,4,−4],[−4,−4,1,4],[−4,4,−4,1]] (Drury/Zheng). V is the 5×5 unitary of `../orthogonal-hull-five/exact_ustar.py`, which lies outside conv r(O(5)). F_n is the Fourier matrix. R2(φ, ψ) = [[cos φ, −sin φ e^{iψ}], [sin φ e^{−iψ}, cos φ]]. DitaAxB(p) = (F_A ⊗ I_B)·blockdiag(diag(e^{ip_k}) F_B), the Diţă construction. It contains F_A ⊗ F_B and the affine Fourier families of F_6 and F_8. Diagonal phases on either side of U do not change d, so the entrywise phase perturbations that keep U unitary are these Diţă phases.

| family | n | params | restarts | best m, free spectra | contact (free) | best m, separated spectra | contact (separated) |
|---|---|---|---|---|---|---|---|
| F5 | 5 | 0 | 8 | -1.15e-05 | 2 active; min sep 9.5e-06 (a+b) | -5.44e-06 | 2 active; min sep 1.0e-01 (a-a) |
| Ostar+1 | 5 | 0 | 4 | -1.37e-06 | 2 active; min sep 9.3e-04 (a-a) | -1.76e-06 | 1 active; min sep 1.0e-01 (b-b) |
| V | 5 | 0 | 4 | -5.52e-06 | 1 active; min sep 7.7e-06 (a-a) | -1.61e-05 | 1 active; min sep 1.0e-01 (a-a) |
| F6 | 6 | 0 | 4 | -3.24e-07 | 20 active; min sep 2.9e-07 (b-b) | -1.63e-01 | 1 active; min sep 1.0e-01 (a-a) |
| Dita2x3 | 6 | 6 | 2 | -2.05e-05 | 4 active; min sep 3.4e-07 (a-a) | -2.20e-01 | 1 active; min sep 1.0e-01 (b-b) |
| Dita3x2 | 6 | 6 | 4 | -2.05e-07 | 70 active; min sep 1.6e-08 (b-b) | -9.83e-02 | 1 active; min sep 1.0e-01 (a+b) |
| Ostar+I2 | 6 | 0 | 4 | -3.33e-06 | 1 active; min sep 3.6e-05 (a-a) | -4.13e-05 | 1 active; min sep 1.0e-01 (b-b) |
| Ostar+F2 | 6 | 0 | 6 | -1.62e-06 | 4 active; min sep 3.6e-08 (b-b) | -8.29e-07 | 3 active; min sep 1.0e-01 (a-a) |
| Ostar+R2 | 6 | 2 | 4 | -5.62e-06 | 1 active; min sep 3.1e-05 (b-b) | -1.32e-06 | 2 active; min sep 1.0e-01 (b-b) |
| V+1 | 6 | 0 | 4 | -1.66e-06 | 2 active; min sep 2.7e-02 (a+b) | -2.40e-05 | 1 active; min sep 1.0e-01 (a-a) |
| F7 | 7 | 0 | 2 | -2.72e-01 | 8 active; min sep 3.1e-10 (b-b) | -2.47e-01 | 1 active; min sep 1.0e-01 (a-a) |
| Ostar+F3 | 7 | 0 | 4 | -1.04e-06 | 6 active; min sep 6.2e-08 (a-a) | -8.81e-03 | 1 active; min sep 1.0e-01 (b-b) |
| V+F2 | 7 | 0 | 4 | -2.12e-06 | 14 active; min sep 2.5e-07 (b-b) | -1.07e-01 | 1 active; min sep 1.0e-01 (b-b) |
| V+I2 | 7 | 0 | 6 | -2.90e-07 | 16 active; min sep 4.5e-04 (a+b) | -1.38e-06 | 1 active; min sep 1.0e-01 (a-a) |
| F8 | 8 | 0 | 2 | -8.47e-07 | 576 active; min sep 9.0e-11 (b-b) | -3.15e-01 | 1 active; min sep 1.0e-01 (a-a) |
| Dita2x4 | 8 | 8 | 2 | -3.71e-06 | 2 active; min sep 4.0e-03 (a+b) | -2.42e-01 | 1 active; min sep 1.0e-01 (b-b) |
| Dita4x2 | 8 | 8 | 2 | -3.64e-05 | 120 active; min sep 4.2e-07 (a-a) | -3.68e-01 | 1 active; min sep 1.0e-01 (a+b) |
| Ostar+Ostar | 8 | 0 | 2 | -2.75e-05 | 1 active; min sep 7.0e-06 (b-b) | -2.30e-02 | 1 active; min sep 1.0e-01 (a-a) |
| Ostar+F4 | 8 | 0 | 2 | -2.65e-06 | 8 active; min sep 3.1e-08 (a-a) | -1.42e-01 | 1 active; min sep 1.0e-01 (a-a) |
| Ostar+I4 | 8 | 0 | 2 | -1.02e-06 | 3 active; min sep 4.5e-05 (a-a) | -1.82e-05 | 1 active; min sep 1.0e-01 (a-a) |
| OstarxR2 | 8 | 2 | 2 | -4.43e-06 | 1 active; min sep 3.6e-06 (a-a) | -5.52e-02 | 1 active; min sep 1.0e-01 (b-b) |
| R2xOstar | 8 | 2 | 2 | -8.08e-05 | 1 active; min sep 4.2e-05 (b-b) | -8.45e-02 | 1 active; min sep 9.9e-02 (a-a) |
| V+F3 | 8 | 0 | 2 | -1.33e-06 | 6 active; min sep 2.2e-09 (a-a) | -2.17e-01 | 1 active; min sep 1.0e-01 (a+b) |
| V+I3 | 8 | 0 | 4 | -2.56e-07 | 36 active; min sep 1.1e-08 (b-b) | -2.51e-02 | 1 active; min sep 1.0e-01 (a-a) |

"Contact" counts the permutations within 1e-7·diam of the support line at the optimizer's θ. "min sep" is the smallest of |a_i − a_j|, |b_i − b_j| and |a_i + b_j| over ρ, labelled by which one it is. The separated runs all end on the 0.1 bound, so the penalty is active.

## Reading

- With free spectra, every family except F7 converges toward m = 0⁻ (best −2e-7 to −8e-5) while a separation shrinks: to 1e-5 to 1e-11 in 19 of the 23, and to 3e-2 to 4e-4 in V⊕1, Dita2x4, O*⊕1 and V⊕I2, where it is still falling. Repeated eigenvalues and a_i + b_j = 0 are known zero-margin degenerations. These contacts are therefore trivial and not specific to the structure.
- With separated spectra, the families split. Within 5e-5 of zero: the n = 5 fixed matrices (F5 −5.4e-6, O*⊕1 −1.8e-6, V −1.6e-5) and the direct sums with a small identity or 2×2 block (O*⊕I2, O*⊕F2, O*⊕R2, O*⊕I4, V⊕1, V⊕I2: −8e-7 to −4e-5). All of these end on the 0.1 separation bound, so the gap is still closing toward a degenerate limit. Far inside: the full Hadamard families (F6, F7, F8 and every Diţă family, −0.10 to −0.37), the Kronecker products O*⊗R2 and R2⊗O* (−0.06, −0.08), O*⊕O* (−0.023), O*⊕F3 (−0.009), O*⊕F4 (−0.14), V⊕F2 (−0.11), V⊕F3 (−0.22) and V⊕I3 (−0.025).
- F7 got only 2 restarts, and one of them was cut by the time limit, so F7 is under-explored. Every n = 8 family got only 2 to 4 restarts.

## Reproduce

Python 3.11 with numpy, scipy and mpmath. Run this from this folder with `OMP_NUM_THREADS=1`. The commands are the ones that produced `results/`. Wall time is per family, so the restart counts depend on the machine.

```
python -I -B structured_search.py --selftest
python -I -B structured_search.py 2026 40  results/n5.json  F5 Ostar+1 V
python -I -B structured_search.py 2026 60  results/n6.json  F6 Dita2x3 Dita3x2 Ostar+I2 Ostar+F2 Ostar+R2 V+1
python -I -B structured_search.py 2026 80  results/n7.json  F7 Ostar+F3 V+F2 V+I2
python -I -B structured_search.py 2026 120 results/n8a.json F8 Dita2x4 Dita4x2 Ostar+Ostar Ostar+F4
python -I -B structured_search.py 2026 120 results/n8b.json Ostar+I4 OstarxR2 R2xOstar V+F3 V+I3
python -I -B summarize.py
```

Each `results/n*.json` holds, per family, the restart count, the m quantiles of the free runs, and the 2 best free and 2 best separated records: a, b, θ, p, m, the active permutations and the separations. `SHA256SUMS` covers the committed files (LF).

## Next

- A larger budget per family (≥ 50 restarts at n = 7, 8), especially for F7.
- Joint optimization of U near these families. Fixed U cannot reach the lifted O* contacts.
- Families not run here: Haagerup, Björck circulant, Petrescu P7, Butson H(q, n) with q ≠ n.
