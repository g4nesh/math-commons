# The hard facets of P₅(c): what they reduce to

Task `mdo-pmsh-5-orbits` · agent `ganesh-claude-pmsh` · owner g4nesh · 26 November 2026. Continues [`../README.md`](../README.md) (task `mdo-pmsh-5`).

## Status

- **PMSH(5) is still not proved.** The facet harvest is partial, and 16 of the 34 orbits of hard facets found need inequalities that are only numerically true.
- **Proved (hand proof plus exact computation; not independently reviewed):**
  - Lemma A, an inequality on every contraction;
  - from it, 74 of the 92 hard facets found at c = (0, 1/4, 1/2, 3/4, 1), all members of 17 orbits that hold 1,380 facets of P₅(c).
- **Exact (rational arithmetic):** the orbits, the decomposition of every hard facet into facets of the minor polytope Q₅, nine orbits of facets of Q₅ that are new relative to [`../../n5-contact-program/strata/`](../../n5-contact-program/strata/README.md) (PR 115), and R1 at n = 5.
- **Numerical:**
  - which facets of Q₅ fail on U(5);
  - the spectral families of the hard facets;
  - two closed-form, spectrum-free (★)-type inequalities, for orbits 20 and 21;
  - a family test: no hard facet is crossed at 791 facet–spectrum pairs.
- **Correction to [`../README.md`](../README.md), N3.** "62 of 92 hard facets do not follow linearly from PMSH(4)" is wrong; the right count is **at most 33**. See R2.

## Setting

- Notation is that of [`../README.md`](../README.md): P(c) = conv{v_π(c)} with v_π(c)_T = ∏_{t∈T} c_{π(t)}.
- Q₅ = conv{x(σ) : σ ∈ S₅} is the minor polytope of the strata README.
  - Coordinates: x_{K→J}(σ) = 1[σ(K) = J] for 1 ≤ |K| = |J| ≤ 2. The 3×3 and 4×4 coordinates equal those of the complements.
  - For a unitary U, D(U) has entries x_{K→J} = |det U[K,J]|². With H = U·diag(c)·U*, Cauchy–Binet gives m(H) = M_c·D(U) + const, and m(diag(c_π)) = v_π(c).
  - So P(c) is the image of Q₅ under an affine map that keeps the vertex labels.
- **Gap vector.** For a facet ℓ·x ≤ h of P(c), its gap vector is f(σ) = h − ℓ·v_σ(c) on the 120 permutations. It is nonnegative and vanishes exactly on the facet's vertex set Z.
  - f is the list of values on the vertices of an affine function on aff(Q₅), namely h − ℓ·(M_c x + const).
  - So f lies in the gap cone of Q₅, whose extreme rays are the facets of Q₅ (Lemma of the strata README).
- **Hard facets.** A facet is *hard* if f is not a nonnegative combination of squared-minor facets x_{K→J} ≥ 0 (as in [`../README.md`](../README.md)).
- **The spectrum.** As in `../facet_classify.py`, c + 1 is used; it is equivalent to (4, 5, 6, 7, 8), so all vertex coordinates are integers.

## Results

### R1. The 92 hard facets form 34 orbits (exact)

- **Symmetry group.** Positions (π ↦ π∘σ, from H ↦ PHPᵀ) and the reversal (π ↦ r∘π with r(j) = 4 − j, from H ↦ I − H) generate a group of order 240 acting on P(c). m(I − H) is an integer-affine image of m(H).
- **Orbits.** The 92 hard facets of `../hard_n5_seed41_c0_0.25_0.5_0.75_1.npz` fall into 34 orbits. Completing them gives 4,020 hard facets of P₅(c), in orbits of size 60, 120 or 240.
- **Certification.** For each representative ([`orbits.py`](orbits.py)), in integer arithmetic:
  - the vertex set has affine rank 25 = dim P₅(c) − 1;
  - the affine function vanishing on it is one-signed on all 120 vertices.
- **Coverage.** The harvest reached most orbits through a single member. P₅(c) certainly has more hard orbits than these 34.

### R2. Correction: at most 33 of the 92 hard facets fail to follow linearly from PMSH(4)

- **Why status is an orbit property.** Over all shifts γ, the constraints that PMSH(4) imposes (README N3 (a), (b)) are invariant under the symmetry group:
  - positions permute the deleted index t;
  - with the spectrum c + 1, the reversal H ↦ 3 − H maps constraint (b) at γ to (b) at −3 − γ.
  So following linearly from PMSH(4) is a property of the whole orbit.
- **The mixed orbits.** Orbits 1, 2, 6 and 8 (53 of the 92 facets) contain members that the first N3 run cut to ≤ 1 and members it did not cut. Example: a member of orbit 1 with relaxation 1.131579 was cut to 0.99999999932 by the 15-shift grid.
- **What went wrong.** The finite shift grids were not closed under γ ↦ −3 − γ, so the cut members' certificates were never applied to their images. With orbit 5 (6 facets, all cut), 59 of the 92 follow linearly from PMSH(4); at most **33 (29 orbits)** do not.
- **Rerun with a reversal-closed 82-shift grid** ([`inherit_orbits.py`](inherit_orbits.py), `logs/inherit_orbits_sym82.txt`). This grid is the 40 log-spaced shifts, their images, 0 and −3.
  - Orbits 0, 9, 18, 19 and 26 are still exceeded: 1.104, 1.043, 1.194, 1.081, 1.054.
  - Orbit 1's representative is also not cut (1.099). The shifts that cut it, −1 and −1.5, lie inside the range of c + 1 and are not on this grid. A union grid confirms the cut: `logs/inherit_union_*.txt`; see the table note.
- **This also shrinks README N3's main claim.** N3 said that 62 hard facets need "new five-dimensional inequalities". By R4, 12 of the 29 remaining orbits follow from one such inequality, Lemma A, which is proved.

### R3. Every hard facet is an exact nonnegative combination of facets of Q₅; nine facet orbits of Q₅ are new

- **Method** ([`decompose.py`](decompose.py)). For every orbit representative, an LP writes f as Σ μ_F g_F + Σ ν_{K,J} x_{K→J} with μ, ν ≥ 0. The g_F run over the 57,780 non-minor facets of orbits A–G in the strata README, and the LP minimizes the total non-minor weight.
- **When the known list falls short.** For seven representatives the LP over the known facets is infeasible. An exact Carathéodory peel then decomposes f into extreme rays of the gap cone:
  - subtract squared-minor rays as far as possible;
  - walk from the remainder to an extreme ray of its face by exact random pushes inside the face;
  - subtract that ray, and repeat.
- **Nine new orbits.** Every non-minor ray found that is not among A–G is certified exactly as a facet of Q₅ (affine rank 40 and one-signed), twice, by two code paths ([`decompose.py`](decompose.py), [`exact_checks.py`](exact_checks.py)). Their orbits under the group of order 28,800 give nine orbits, N1–N9, of non-minor facets of Q₅ that are not in `facets_n5.json`.
- **Exactness.** All 34 decompositions are recomputed in rational arithmetic with positive rational weights (`exact_checks.py`; `logs/exact_checks.txt`).
- **Answers an open question of PR 115.** The strata README asks whether its seven orbits are complete. They are not: Q₅ has at least 16 orbits of non-minor facets.

Minimum-ℓ₁ integer forms ([`facet_forms.py`](facet_forms.py); x_{KJ} = x_{K→J}). The forms of A, B, D, E, F, G coincide with the strata README.

| Orbit | Vertices | Orbit size | Inequality on the vertices of Q₅ | Minimum over U(5) |
| --- | --- | --- | --- | --- |
| A | 84 | 900 | x12>12 ≤ x0>0 + x03>34 + x04>34 | **0 (Lemma A)** |
| N8 | 72 | 3,600 | x14>02 ≤ x4>0 + x04>12 + x12>01 + x13>01 | +3·10⁻¹¹ (numerically valid) |
| N9 | 74 | 3,600 | x01>24 ≤ x1>4 + x03>34 + x04>34 + x34>01 | −0.0674 |
| N7 | 68 | 3,600 | x34>02 ≤ x01>34 + x03>23 + x04>24 + x12>34 + x23>04 + x24>03 | −0.0674 |
| N5 | 66 | 14,400 | x04>14 ≤ x01>04 + x01>13 + x14>04 + x23>03 + x24>01 + x34>01 | −0.0674 |
| N6 | 66 | 28,800 | x24>23 ≤ x01>04 + x02>34 + x04>12 + x04>13 + x13>14 + x23>02 | −0.0674 |
| N4 | 64 | 7,200 | x01>01 ≤ x04>02 + x04>13 + x14>02 + x14>04 + x14>13 + x14>14 | −0.0674 |
| N3 | 59 | 28,800 | x02>34 + x03>23 ≤ x0>3 + x4>0 + x01>04 + x04>24 + x12>04 + x34>01 | −0.0674 |
| N2 | 58 | 14,400 | x01>34 + x23>23 ≤ x0>3 + x4>0 + x01>04 + x04>24 + x14>04 + x14>14 + x23>02 | −0.0523 |
| N1 | 54 | 28,800 | x01>14 + x13>04 ≤ x1>4 + x4>0 + x02>23 + x03>34 + x12>34 + x13>01 + x14>02 + x14>03 + x23>02 | −0.1349 |

Orbits B–G of the strata README (inequalities there) have minima over U(5) of: B, E, F, G −0.0674; C −0.0894; D −0.0759.

### R4. Lemma A: the block inequality holds on every contraction (proved)

> **Lemma A.** Let W be a 3×3 complex matrix with rows (a, b, c) and columns (a′, b′, c′) that is a submatrix of a contraction (in particular, of a unitary matrix). Then
> |det W|² ≤ |w_{aa′}|² + |det W[{a,b},{b′,c′}]|² + |det W[{a,c},{b′,c′}]|².

*Proof.*
1. Let p_i be the 2×2 minor of W on the columns (b′, c′) with row i deleted. Expanding along the column a′:
   det W = w_{aa′} p_a − w_{ba′} p_b + w_{ca′} p_c.
2. By Cauchy–Schwarz,
   |det W|² ≤ (|w_{aa′}|² + |p_b|² + |p_c|²) · (|p_a|² + |w_{ba′}|² + |w_{ca′}|²).
3. The second factor is at most 1. The rows b, c of W form a 2×3 contraction R = [v | N], with v = (w_{ba′}, w_{ca′})ᵀ and N the 2×2 block on (b′, c′).
   - R R* ≤ I gives N N* ≤ I − v v*.
   - Determinants are monotone on positive semidefinite matrices, so |det N|² ≤ det(I − v v*) = 1 − ‖v‖².
   - Since p_a = det N, the second factor is |p_a|² + ‖v‖² ≤ 1.
4. Finally, p_b = ±det W[{a,c},{b′,c′}] and p_c = ±det W[{a,b},{b′,c′}]. ∎

- **In Q₅ coordinates.** For rows S = {a,b,c} and columns S′ = {a′,b′,c′} of a 5×5 unitary, Jacobi's identity gives |det U[S,S′]| = |det U[Sᶜ,S′ᶜ]|. So Lemma A is the facet inequality of orbit A of Q₅ evaluated at D(U). The orbit is closed under row and column permutations and transposition, which preserve contractions.
- **Equality.** For W unitary (U = W ⊕ V), the cofactors give equality, so A is tight on a large set.
- **Numerical checks.** The minimum over U(5) found is +4·10⁻¹¹. On 200,000 random 3×3 contractions the gap is positive.

**R1 at n = 5 (exact; [`exact_checks.py`](exact_checks.py)).** D(U) ∈ aff(Q₅) for every 5×5 unitary U.
- aff(Q₅) has 84 independent affine relations, and they are spanned (exact rank) by identities valid for every unitary:
  - rows and columns of each compound matrix are unit vectors;
  - Σ_{T∋j} |det U[S,T]|² = Σ_{s∈S} |u_{sj}|² for |S| = 2. Proof: Cauchy–Binet for R = U[S,·] gives Σ_{T∌j} |det R_T|² = det(I − r_j r_j*) = 1 − ‖r_j‖², where r_j is column j of R. The transposed identity holds too.

> **Corollary A.** At c = (0, 1/4, 1/2, 3/4, 1), every hard facet in the 17 orbits 0, 1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 14, 15, 18, 24, 26, 31 holds for every Hermitian matrix with spectrum c. These orbits hold 1,380 facets, including 74 of the 92 hard facets found and every orbit that follows from PMSH(4).

*Proof.*
1. f = Σ μ_F g_F + Σ ν x_{K→J} on aff(Q₅) with rational μ, ν > 0 and every g_F in orbit A (R3, exact).
2. Evaluate at D(U), which lies in aff(Q₅) by R1.
3. Each g_F(D(U)) ≥ 0 by Lemma A, and every squared minor is ≥ 0. ∎

- **The same holds across spectra (numerical).** By R5 these decompositions keep the same support with positive weights over the whole spectral region tested, so the conclusion extends there.
- **Orbit 27** decomposes into a facet of N8 and squared minors. It follows from N8 ≥ 0 on U(5), which is numerically true (minimum +3·10⁻¹¹) but not proved.

### R5. Spectral families, and the remaining obstacles

- **Möbius invariance (proved).** Let φ(x) = (αx + β)/(γx + δ) with γc_i + δ ≠ 0 for all i, and let K = γH + δI.
  - Jacobi's identity gives det(K⁻¹)[T,T] = det K[Tᶜ,Tᶜ]/det K, and det K = ∏(γc_i + δ) is the same for every H with spectrum c.
  - Principal minors of aI + bL and of γH + δI are affine in those of L and H.
  - Hence m(φ(H)) = A(m(H)) for an affine bijection A that depends only on (φ, c). It maps v_π(c) to v_π(φ(c)), so P(φ(c)) = A(P(c)) with labels kept.
  - So PMSH(n) at c is equivalent to PMSH(n) at φ(c) (this extends the affine invariance of Proposition H2(b) of `../../minors/pmsh.md`), and the vertex sets and gap vectors of facets depend only on the Möbius class of c.
- **Coordinates for spectra.** Classes of five ordered points are written (0, 1, s, t, ∞) with 1 < s < t, with finite representative c(s, t) = (0, 1/2, s/(s+1), t/(t+1), 1). The arithmetic spectrum is (s, t) = (3, 9).
- **Families.** For each orbit, [`fit_families.py`](fit_families.py) keeps Z, recomputes the facet at 492 spectra, re-solves the decomposition on the same support, and fits each weight (non-minor weight normalized to 1) by a rational function of (s, t).
  - Results are in `logs/fit_families.txt`.
  - Single-facet orbits give clean formulas, mostly monomials in s, t, s − 1, t − 1, t − s.
  - For some orbits with several non-minor facets (for example orbit 4), weights are not fitted at degree ≤ 4: the log shows relative errors up to 10⁻². Those entries are not used anywhere.
  - Orbit 19's support at (3, 9) gives the facet at none of the grid points.
  - The families of orbits 15, 20 and 21 are checked exactly at three random rational spectra each ([`exact_families.py`](exact_families.py)).
  - Example, orbit 15 (proved by Lemma A), with a = (t − s)/s and b = (t − s)/(s − 1), for all 0 < a < b:
    g_A + a·x01>13 + a(1 + b)·x01>34 + b·(x12>04 + x13>04 + x14>04 + x24>12 + x34>12) + a(1 + b)/(b − a)·x23>01 ≥ 0.
- **Which facets of Q₅ fail on U(5)** ([`unitary_valid.py`](unitary_valid.py), 40 BFGS starts per orbit). Only A (proved) and N8 stay ≥ 0.
  - Ten orbits fail with minimum exactly −0.06744: B, E, F, G, N3–N7, N9. This is the minimum of the exotic facet g0 of Q₄ over U(4) (`../../pmsh/equivalence.md`), so these are most likely dimension-four obstructions lifted to n = 5.
  - C, D, N1 and N2 fail with other values: −0.0894, −0.0759, −0.1349, −0.0523.
  - Every hard facet orbit outside R4 uses at least one failing facet. These 16 orbits are the genuine obstacles: there the squared-minor terms must compensate, as in (★) for n = 4.
- **Closed forms for two obstacle families (numerical, family formulas exact-checked).**
  - Orbits 20 and 21 use one facet of orbit E. With w = t − 1 and λ = (t − s)/(t − 1), every weight is nondecreasing in w. So the infimum over spectra is the limit w → 0, where three eigenvalues collide at a fixed internal ratio λ, followed by AM–GM in λ.
  - The orbit-20 family over all spectra is therefore **equivalent** to
    > (E20) g + x04>01 + x04>14 + 2√(x01>12·(x04>01 + x04>14)) ≥ 0 on U(5),
  - and orbit 21 to
    > (E21) g + x01>23 + 2√(x01>23·(x02>12 + x04>12)) ≥ 0,
  - with g the E-facet of each decomposition. These have the same shape as (★), g0 + D_a + 2√(D_a D_b) ≥ 0.
  - Minimum over U(5), 60 starts each ([`envelope_min.py`](envelope_min.py)): +9·10⁻¹² and +5·10⁻¹¹, while g alone reaches −0.0674. Not proved.
- **Family test** ([`family_min.py`](family_min.py)). Every one of the 34 orbits was maximized over U(5) at 25 spectra (s ∈ {1.1, 1.5, 3, 8, 30}, t/s ∈ {1.1, 1.5, 3, 8, 30}), 24 starts with exact gradients, wherever the facet exists (791 facet–spectrum pairs).
  - Largest value 1 + 2.9·10⁻¹⁰: nothing is crossed.
  - Unlike the four spectra of `../README.md` N1, this reaches near-collisions such as (s, t) = (1.1, 1.21).

## Orbit table

Generated by [`summary.py`](summary.py).

| orbit | size | found | vertices | overshoot | PMSH(4), README N3 grid | PMSH(4), reversal-closed grid | non-minor Q5 facets (weight) | squared minors | status | family max - 1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
TABLE_PLACEHOLDER

- "PMSH(4), README N3 grid": "cut for some members" means some member was cut to ≤ 1 by the first N3 run; by R2 the whole orbit then follows linearly from PMSH(4).
- "reversal-closed grid": maximum over R_inh with the 82-shift grid of R2, for representatives run before the runs were stopped.
- "overshoot": relaxation overshoot over {D ≥ 0}.
- "weight": the rational multiplier of each non-minor facet, relative to the primitive integer gap vectors.

## What this means for PMSH(5)

At this spectrum, among the hard facets found, the obligation for PMSH(5) is now:

1. **Solved by Lemma A:** 17 orbits.
2. **N8 ≥ 0 on U(5):** orbit 27 only.
3. **(★)-type inequalities:** the 16 obstacle orbits, each a family over the two Möbius parameters (s, t) of the form g_F + Σ w_a(s, t)·x_a ≥ 0 with explicit positive weights. For orbits 20 and 21 the family is equivalent to one closed-form inequality, (E20) or (E21).

Ten of the eleven non-minor facets used by obstacle orbits (all but D) fail on U(5) with exactly the n = 4 minimum. A natural conjecture: the obstacle inequalities are images of (★) on 4×4 subsystems. If so, PMSH(5) at this spectrum would reduce to (★) and Lemma A. That has **not** been checked. N3 says the obstacles are not *linear* consequences of PMSH(4), but (★) applied to blocks of U is not a linear constraint in that sense.

## Limits

- **Coverage.** The harvest is partial: 34 orbits reached through 92 facets. The other hard facets of P₅(c) are unknown, and the decomposition and the conclusions above cover only the orbits found.
- **One spectrum.** The exact results hold at the Möbius class of (0, 1/4, 1/2, 3/4, 1). Across spectra the evidence is numerical: fitted families checked exactly at three points each for three orbits, and the family test.
- **Numerical evidence only:** that facets fail or hold on U(5), (E20), (E21) and N8. Local optimization can miss a region.
- **PMSH(4) rerun incomplete.** The reversal-closed inheritance LP was run for six representatives only, then stopped in favour of R3–R5.
- **No independent review yet.**

## Files

| File | What it does |
| --- | --- |
| [`orbits.py`](orbits.py) | Orbits of the 92 hard facets, exact facet certificates, integer gap vectors → `orbits.json` |
| [`inherit_orbits.py`](inherit_orbits.py) | PMSH(4) inheritance LP on orbit representatives with a reversal-closed shift grid (sparse) |
| [`decompose.py`](decompose.py) | Decomposition into facets of Q₅; exact Carathéodory peel; new facet orbits → `decompose.json` |
| [`exact_checks.py`](exact_checks.py) | R1 at n = 5, independent certification of N1–N9, exact recheck of all decompositions |
| [`facet_forms.py`](facet_forms.py) | Minimum-ℓ₁ integer forms of the 16 non-minor orbits → `facet_forms.json` |
| [`unitary_valid.py`](unitary_valid.py) | Minimum of each non-minor facet over U(5) |
| [`multipliers.py`](multipliers.py), [`fit_families.py`](fit_families.py) | Weights of each decomposition along the Möbius coordinates (s, t); rational fits |
| [`exact_families.py`](exact_families.py) | Exact check of the families of orbits 15, 20, 21 |
| [`envelope_min.py`](envelope_min.py) | (E20), (E21) over U(5) |
| [`family_min.py`](family_min.py) | Maximum of every orbit's facet over U(5) at 25 spectra |
| [`summary.py`](summary.py) | The orbit table |

**Reproduce.** Run from this directory with numpy, scipy (HiGHS) and sympy; no cddlib is needed. Times are on one core of a laptop.

```
python3 -I -B orbits.py                 # 30 s
python3 -I -B decompose.py              # 10 min
python3 -I -B exact_checks.py           # 5 min
python3 -I -B facet_forms.py
python3 -I -B unitary_valid.py 40 0 A   # one orbit label per run; 2-5 min each
python3 -I -B fit_families.py 0,1,...,33
python3 -I -B exact_families.py 3
python3 -I -B envelope_min.py 60
python3 -I -B family_min.py 0,1,...,33 24 1   # split into four runs: about 5 min each
python3 -I -B summary.py
GRID=union python3 -I -B inherit_orbits.py 1:2   # about 10 min per orbit
```

Logs of every run reported here are in [`logs/`](logs/).
