# PMSH(5): every hard facet found is touched but not crossed; two thirds of them need new five-dimensional inequalities

- **Task:** `mdo-pmsh-5`. Agent `ganesh-claude-pmsh` (Claude Code, model `claude-opus-5-5`), GitHub owner `g4nesh`. November 26, 2026.
- **Status: PMSH(5) is NOT proved here.** Numerical evidence (N1, N2) and a structural reduction (N3): which hard facets follow linearly from PMSH(4), and which do not. Nothing here is independently reviewed. **N3's count is corrected in [`facets/`](facets/README.md):** 69 of the 92 hard facets follow linearly from PMSH(4) (not 30), and 74 follow from the proved Lemma A; every hard facet found is covered by one of the two or lies in 16 obstacle orbits.
- **Why PMSH(5) matters.** By Theorem H1 of [`../minors/pmsh.md`](../minors/pmsh.md), PMSH(5) gives MOC for every 5×5 pair with one Hermitian matrix. With the lifting Theorem L of [`../hermitian-lifting/README.md`](../hermitian-lifting/README.md) it then gives MOC for A Hermitian and B with at most **five** nonreal eigenvalues, in every dimension. Known today: three (Shitov 2025); four (Corollary L4, conditional on review of the n = 4 base case).

## Statement tested

For Hermitian H with spectrum c, let m(H) = (det H[T,T])_{∅≠T⊆[n]} and P(c) = conv{v_π(c) : π ∈ S_n} with v_π(c)_T = ∏_{t∈T} c_{π(t)}.

> **PMSH(n).** m(H) ∈ P(c) for every Hermitian H with spectrum c.

PMSH is invariant under c ↦ αc + γ (Proposition H2(b) of `../minors/pmsh.md`), so spectra are normalized to min 0, max 1, and all computations use the shifted spectrum c + 1 (products in [1, 2ⁿ]; this matters, see "Pitfalls").

## Results

**N1 (numerical). The hard facets of P₅(c) are touched by U(5) but not crossed.** This is the main evidence.
- **Easy and hard.** Write D(U) for the vector of all squared minors |det U[S,T]|², so that m(U D_c U*) = M_c·D(U) + const. A supporting inequality ℓ·(x − g) ≤ 1 of P(c) is *easy* if it holds on the linear relaxation {D ≥ 0} ∩ aff{D(U)}, i.e. it follows from nonnegativity of squared minors (Farkas), and *hard* otherwise. Only hard inequalities can be violated by a unitary.
- **Method (`facet_classify.py`).** Shoot rays from the centroid g toward m(U D_c U*), half of them for Haar U and half for slightly perturbed permuted block-diagonal U. The dual of each ray LP is a supporting inequality, kept if valid on all vertices and identified by its vertex set; it is a true facet when that vertex set spans dimension 25. Classify each by the relaxation LP. Maximize every hard one over U(5): L-BFGS from 12 random starts with exact gradients in a Cayley chart (`facet_search.maximize_over_U`, gradient checked against finite differences to 10⁻⁸).
- **Control (n = 4, c = (0,1,2,3); `logs/n4_hard_control.txt`).** 100 of the 184 facets and 16 of the 96 hard facets are reached (22 hard inequalities in all); every hard maximum over U(4) is 1 to 2·10⁻¹³.
- **n = 5** (`logs/n5_hard.txt`), 6,000 rays per spectrum. "Hard" means relaxation overshoot above 10⁻⁶; overshoots between 10⁻⁹ and 10⁻⁶ are within the LP tolerance (the normals come from LP duals) and are listed separately as noise.

  | spectrum c | facets reached | hard facets | hard inequalities | noise-level | relaxation overshoot | max over U(5) − 1 | min over U(5) − 1 |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | (0, 1/4, 1/2, 3/4, 1) | 342 | 102 | 458 | 189 | 0.54 | 1.6·10⁻¹¹ | −8.4·10⁻¹³ |
  | (0, 0.129, 0.499, 0.601, 1) | 393 | 140 | 573 | 188 | 0.62 | 2.5·10⁻¹¹ | −1.1·10⁻¹² |
  | (0, 0.36, 0.70, 0.83, 1) | 330 | 104 | 551 | 216 | 0.53 | 6.0·10⁻¹¹ | −2.5·10⁻¹² |
  | (0, 0.1, 0.2, 0.9, 1) | 441 | 170 | 609 | 180 | 0.66 | 1.2·10⁻¹⁰ | −6.6·10⁻¹² |

  So 516 hard facets and 2,191 hard inequalities in all. **Every one has maximum exactly 1 over U(5), to 1.2·10⁻¹⁰**: each is touched by some unitary and none is crossed. The linear relaxation overshoots them by up to 0.66, so these are genuine tests. (A first pass with threshold 10⁻⁹ counted 821 hard facets; the extra ones are the noise-level ones, and their maxima were also exactly 1.)
- **Where they are touched.** For 1,044 of the 2,191 hard inequalities, some maximizer found is irreducible: the graph linking rows that share a column with |U_ij| > 10⁻³ is connected. For the other 1,147 every maximizer found is block-reducible (4+1, 3+2, …). So, numerically, PMSH(5) is tight at genuinely five-dimensional unitaries, and its contact set is not just inherited from dimension four. (At n = 4 the hard facets are likewise touched at irreducible U, which is why (★) is needed there.)
- **What this does not show.** The harvest reaches about 330–440 facets per spectrum; P₅(c) has more, and the true facet count is unknown (qhull cannot handle dimension 26 and pycddlib could not be built here). Twelve starts per inequality can miss a narrow violating region. Four spectra are not all spectra.

**N2 (numerical). Broad search over supporting inequalities, n = 5 and n = 6.** `facet_search.py` is the first version of N1 without classification: Haar and near-permutation rays, 6 starts with finite-difference gradients, all supporting inequalities reached (deduplicated by rounded normals, so the counts include near-duplicates and non-facets). No confirmed candidate; apart from the near-repeated spectrum, no value above 1 + 1.1·10⁻¹⁰. The special spectra are from `run_spectra.sh` (`logs/n5_special.txt`; its last two spectra were not run to the end), the random ones from `logs/n5_random_s*.txt` and `logs/n6_random_s*.txt`:

  | n | spectrum c | inequalities | max − 1 |
  | --- | --- | --- | --- |
  | 5 | (0, 1/4, 1/2, 3/4, 1) | 633 | 8·10⁻¹² |
  | 5 | (0, 0.01, 0.02, 0.03, 1) | 1,391 | 5·10⁻¹⁴ |
  | 5 | (0, 0.129, 0.499, 0.601, 1) | 667 | 0 |
  | 5 | (0, 0.361, 0.703, 0.831, 1) | 664 | 4·10⁻¹² |
  | 5 | (0, 0.344, 0.693, 0.816, 1) | 650 | 1.5·10⁻¹¹ |
  | 5 | (0, 0.094, 0.431, 0.567, 1) | 674 | 0 |
  | 5 | (0, 0.811, 0.855, 0.865, 1) | 995 | 6.6·10⁻¹¹ |
  | 5 | (0, 0.189, 0.251, 0.947, 1) | 780 | 1.2·10⁻¹¹ |
  | 5 | (0, 0.49, 0.5, 0.51, 1) | 1,123 | 2.8·10⁻¹² |
  | 5 | (0, 0.001, 0.01, 0.1, 1) | 1,445 | 4.5·10⁻¹⁴ |
  | 5 | (0, 0.01, 0.5, 0.99, 1) | 828 | 1.2·10⁻¹² |
  | 5 | (0, 0.3, 0.3001, 0.3002, 1) | 1,478 | 4.2·10⁻⁸ (48 values above 1, all inside by the independent LP) |
  | 6 | (0, 0.089, 0.606, 0.710, 0.781, 1) | 851 | 3.6·10⁻¹¹ |
  | 6 | (0, 0.089, 0.199, 0.366, 0.653, 1) | 711 | −7.5·10⁻¹² |
  | 6 | (0, 0.114, 0.129, 0.641, 0.694, 1) | 919 | 1.0·10⁻¹⁰ |

  Haar rays mostly reach easy facets, which is why N1 was needed: at n = 4 the hard region is thin (the relaxation overshoots by 0.004–0.24 in qhull's unit-normal scaling), and 3,000 random vertices of the relaxation at each of three spectra never land outside P(c) (`relaxation_vertices.py`, `logs/relaxation_vertices_n4.txt`).

**N2′ (reproduces a known count).** With all 184 qhull facets at n = 4, the Farkas test of N1 gives exactly 96 hard facets at each of four spectra, matching R3–R4 of [`../pmsh/README.md`](../pmsh/README.md) (`inherited_n4.py`, first column of its output).

**N3 (structural; LP results in floating point). Which hard facets follow linearly from dimension n − 1?**
- **Constraints inherited from dimension n − 1.** For every position t and real γ, PMSH(n − 1) forces:
  - (a) the projection of m(H) onto subsets not containing t lies in the same projection of P(c). Proof: the compression deleting t has spectrum s interlacing c; PMSH(n − 1) puts its minors in P(s), and the coupling Lemma C of `../hermitian-lifting/` (which needs only interlacing, not positivity) writes every v(s) as a mixture of v_τ(c) over injections τ;
  - (b) the projection of m(H + γI) onto subsets containing t lies in the same projection of P(c + γ). Proof: Jacobi's identity det K[T,T] = det K·det(K⁻¹)[Tᶜ,Tᶜ] for K = H + γI, (a) for K⁻¹, and the fact that scaling a convex hull by the nonzero real det K preserves membership; by continuity also when K is singular.
  These are linear in the squared-minor vector D. Adding them to the nonnegativity relaxation gives a relaxation R_inh; a hard facet *follows linearly from PMSH(n − 1)* if its maximum over R_inh is at most 1. (This uses PMSH(n − 1) only through Lemma C's injection hull, which is weaker than P(s) at the actual compression spectrum s.)
- **n = 4: none do** (`inherited_n4.py`, `logs/inherited_n4.txt`). With (a) for every t and (b) for every t and eight values of γ, all 96 hard facets of P₄(c) are still exceeded, **by exactly the same amount**, at c = (1,2,3,4), (0.3,1.1,2,3.7), (1,1.01,5,9) and (1,1.27,1.64,2). Control: replacing the vertex sets of the constraints (a) by 12 of the 24 permutations removes all excess (max 1.2·10⁻¹¹), so the constraints are imposed; they just never bind. This is consistent with PMSH(4) needing the genuinely four-dimensional inequality (★) of [`../pmsh/equivalence.md`](../pmsh/equivalence.md).
- **n = 5: a third do, two thirds do not** (`inherited.py`, c = (0, 1/4, 1/2, 3/4, 1), the hard facets of `hard_n5_seed41_c0_0.25_0.5_0.75_1.npz`). Of the 165 facets there with relaxation overshoot above 10⁻⁹, 73 overshoot by less than 10⁻⁶; that is within the LP tolerance (the normals come from LP duals), so they are treated as easy. The **92 genuinely hard facets** (overshoot 0.026–0.54):

  | shifts γ used in (b) | cut to ≤ 1 | still exceeded | log |
  | --- | --- | --- | --- |
  | 0, ±0.5, ±2.5 | 0 | 92 | `logs/n5_inherited_rows.txt` |
  | 0 and 15 values in [−10, 10] | 30 | 62 | `logs/n5_inherited_dense_rows.txt` (3 LP failures, redone below) |
  | 0 and 40 log-spaced values in [−40, 40] | 30 | **62** | `logs/n5_inherited_dense40_rows.txt`, `logs/n5_inherited_nan15_rows.txt` |

  The count is stable from 15 to 40 shifts. The 62 survivors are still exceeded by 0.023 to 0.37; 5 of them are not reduced at all.
- **Correction ([`facets/`](facets/README.md), R2).** The counts above are wrong, for two reasons.
  - Whether a facet follows linearly from PMSH(4) is a property of its orbit under positions and the reversal, because the reversal maps the constraint at γ to the one at −3 − γ. The shift grids were not closed under that map, and four orbits had members on both sides of the line.
  - The decisive shifts lie inside [−2, −1]. The 40-shift grid has only two there, so its "stability" compared two grids that are both sparse where it matters.
  With reversal-closed grids (32 and 50 shifts), **69 of the 92 follow linearly from PMSH(4)** (13 of 34 orbits); at most 23 do not.
- **Consequence (as first written; see the correction above and [`facets/`](facets/README.md)).** About a third of the hard facets of P₅(c) found are linear consequences of PMSH(4), which is proved (computer-assisted, unreviewed), but only through shifted Schur complements with well-chosen γ. The remaining ones, 62 of 92 at this spectrum, are not: PMSH(5) needs new, genuinely five-dimensional nonlinear inequalities on U(5), as PMSH(4) needed (★). A proof of PMSH(n) by linear programming from PMSH(n − 1) alone is ruled out at n = 4 and, for this family of inherited constraints, at n = 5.

**N4 (numerical). Splitting along one row.** For H = U D_c U* and a row r, put ρ_j = |U_{rj}|². A *split* is a family of Hermitian H_j of size n − 1 with spectrum c ∖ c_j such that Σ_j ρ_j m(H_j) = m(H with row and column r deleted) and Σ_j ρ_j c_j m(H_j) = (det H[S ∪ r])_S. A split plus PMSH(n − 1) for each H_j gives a PMSH certificate for H, by conditioning on the eigenvalue placed at r (Proposition H3 forces the marginal ρ). The obvious canonical choices of H_j fail (`../minors/pmsh.md`, negative results); `split_test.py` instead solves for the H_j by least squares. Result: a split along the first row exists to residual below 10⁻¹¹ for each of 40 random H at n = 4 and 30 at n = 5 (`logs/split_n4.txt`, `logs/split_n5.txt`). This is a reformulation, not a shortcut: at n = 5 it is equivalent to PMSH(5) if 4×4 Hermitian matrices with a fixed spectrum fill P₄, and for large n the number of conditions (2ⁿ − 2) exceeds the number of unknowns (n(n − 1)(n − 2)) from n = 9 on.

## Pitfalls met along the way (kept so others avoid them)

- **Joint ascent stalls.** Maximizing ℓ·m(U) − max_π ℓ·v_π(c) jointly over (U, ℓ) converges to a permutation matrix, where the gap is exactly 0: all 48 runs at n = 5 did. A gap of 0 from such a search says nothing.
- **Ill-conditioned spectra.** With c = (0, 0.001, 0.01, 0.1, 1) unshifted, the ray LP returned "facets" with entries near 10²¹ that are invalid on the polytope's own vertices, and the search reported values up to 135. Re-checking the same four matrices after c ↦ c + 1 put them strictly inside (ray depth 1.03–1.04, L1 distance 0). The scripts now shift, validate facets, and re-check candidates.
- **Haar rays and random relaxation vertices miss hard facets** (N2). Rays aimed at slightly perturbed block-diagonal unitaries reach them (N1).
- **Gauge ascent over the relaxation stalls.** Alternating "exit facet of the ray, then maximize that facet over the relaxation" stops at once: maximizing an easy facet lands exactly on the boundary.
- **Finite-difference gradients are too slow and too weak** at n = 5 (one start took seconds, and some inequalities were never brought up to 1). The exact Cayley-chart gradient fixed both.

## Open questions, in order

1. *(Done in [`facets/`](facets/README.md): 34 orbits; Lemma A, proved, settles 17 of them; 16 obstacle orbits remain.)* **The 62 hard facets of P₅(c) at c = (0, 1/4, 1/2, 3/4, 1) that do not follow from PMSH(4).** Classify them under S₅ (positions) and the reversal c ↦ 1 − c, write each as a polynomial inequality on U(5), and look for a few spectrum-free generators, as Theorem EQ did for n = 4. Their normals are the rows of `hard_n5_seed41_c0_0.25_0.5_0.75_1_notinherited15.npz` and `..._nan15.npz` whose value in `logs/n5_inherited_dense40_rows.txt` and `logs/n5_inherited_nan15_rows.txt` exceeds 1 + 10⁻⁶ (column 3 is the relaxation value, column 4 the value with inherited constraints).
2. **A stronger inheritance.** Use PMSH(4) at the actual spectrum of each compression (a nonlinear constraint) instead of Lemma C's injection hull. Does that cut the 62?
3. **Completeness.** The harvest reaches 330–440 facets per spectrum and its coverage of the hard facets is unknown. A full facet enumeration of P₅(c) at one rational spectrum (needs cdd or lrs) would settle what the obligations are.
4. **Contacts.** More than half of the hard inequalities are touched at irreducible U(5) (N1). Describe the contact set of the 62 non-inherited facets: it is where any proof must be sharp.

## Reproduction

numpy 2.0.2, scipy 1.13.1 (HiGHS), Python 3.9. Run from this directory.

Set OPENBLAS_NUM_THREADS=1 when running several jobs at once. Data: `hard_n5_seed31_c*.npz` are the hard inequalities of the four N1 spectra; `hard_n5_seed41_c0_0.25_0.5_0.75_1*.npz` are those used in N3 (spectrum `c`, normals scaled so that ℓ·(x − g) ≤ 1 on P(c + 1), vertex counts, relaxation values, facet flags).

```sh
python3 facet_classify.py 4 0 4000 12 0,1,2,3     # control, about 10 s
python3 facet_classify.py 5 31 6000 12 0,0.25,0.5,0.75,1   # one row of the N1 table, 10-20 min
python3 facet_search.py 4 0 600 4                 # N2 control, about 20 s
python3 facet_search.py 5 11 1500 6               # one n = 5 spectrum, 10-30 min
PY=python3 ./run_spectra.sh                       # special n = 5 spectra
python3 facet_classify.py 5 41 6000 0 0,0.25,0.5,0.75,1    # harvest only; writes the .npz used by N3
python3 inherited.py hard_n5_seed41_c0_0.25_0.5_0.75_1.npz -2.5,-0.5,0.5,2.5 --rows rows.txt   # about 30 min; use --range a:b to split
python3 inherited_n4.py 1,2,3,4 -2.5,-1.5,-0.5,0.5,1.5,-4.5,7
python3 inherited_n4.py 1,2,3,4 none --control
python3 split_test.py 5 1 30 12
```
