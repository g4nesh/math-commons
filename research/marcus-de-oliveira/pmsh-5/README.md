# PMSH(5): a facet-driven test, and why a linear induction cannot work

- **Task:** `mdo-pmsh-5`. Agent `ganesh-claude-pmsh` (Claude Code, model `claude-opus-5-5`), GitHub owner `g4nesh`. November 26, 2026.
- **Status: PMSH(5) is NOT proved here.** This is numerical evidence plus one negative structural result. Nothing here is independently reviewed.
- **Why PMSH(5) matters.** By Theorem H1 of [`../minors/pmsh.md`](../minors/pmsh.md), PMSH(5) gives MOC for every 5×5 pair with one Hermitian matrix. With the lifting Theorem L of [`../hermitian-lifting/README.md`](../hermitian-lifting/README.md) it then gives MOC for A Hermitian and B with at most **five** nonreal eigenvalues, in every dimension. Known today: three (Shitov 2025); four (Corollary L4, conditional on review of the n = 4 base case).

## Statement tested

For Hermitian H with spectrum c, let m(H) = (det H[T,T])_{∅≠T⊆[n]} and P(c) = conv{v_π(c) : π ∈ S_n} with v_π(c)_T = ∏_{t∈T} c_{π(t)}.

> **PMSH(n).** m(H) ∈ P(c) for every Hermitian H with spectrum c.

PMSH is invariant under c ↦ αc + γ (Proposition H2(b) of `../minors/pmsh.md`), so spectra are normalized to min 0, max 1, and all computations use the shifted spectrum c + 1 (products in [1, 2ⁿ]; this matters, see "Pitfalls").

## Results

**N1 (numerical). No violation of PMSH(5) or PMSH(6) at the spectra tested.**
- **Method (`facet_search.py`).** For a fixed spectrum, shoot rays from the centroid g of P(c) toward m(U D_c U*) for random U (Haar, and near-permutation). The dual of each ray LP is a facet ℓ of P(c), scaled so that ℓ·(x − g) ≤ 1 on P(c). Facets are kept only if they are valid on all vertices. Each distinct facet is then maximized over U(n) by L-BFGS from several random starts. PMSH(n) at c fails iff some maximum exceeds 1. A value above 1 is re-checked with an independent LP (L1 distance from m(U) to P(c)).
- **Counting caveat.** `facet_search.py` deduplicates by rounded normals, so its counts are of distinct *supporting inequalities*: they include near-duplicates and hyperplanes through lower-dimensional faces (any valid inequality is still a legitimate test). `facet_classify.py` (N2′) identifies true facets by their vertex sets.
- **Control (n = 4).** At c = (0, 0.27, 0.64, 1): 180 distinct supporting inequalities, every maximum over U(4) equal to 1 to 4·10⁻¹⁴, consistent with the proved PMSH(4).
- **n = 5.** Results so far (`logs/n5_*.txt`), 1,500 rays and 6 starts per facet:

  | spectrum c | supporting inequalities | max over U(5) of ℓ·(m − g) | candidates above 1 |
  | --- | --- | --- | --- |
  | (0, 1/4, 1/2, 3/4, 1) | 633 | 1 + 8·10⁻¹² | 0 |
  | (0, 0.129, 0.499, 0.601, 1) | 667 | 1 + 0 | 0 |
  | (0, 0.361, 0.703, 0.831, 1) | 664 | 1 + 4·10⁻¹² | 0 |
  | (0, 0.344, 0.693, 0.816, 1) | 650 | 1 + 1.5·10⁻¹¹ | 0 |
  | (0, 0.094, 0.431, 0.567, 1) | 674 | 1 + 0 | 0 |

  Further spectra (clustered, geometric, near-repeated) and n = 6 are still running; this table will be completed.
- **What this does and does not show.** Every inequality the rays reach is touched (maximum exactly 1) and none is crossed. Haar rays mostly reach **easy** facets, those implied by nonnegativity of squared minors (see N2). The hard facets, where a violation would have to occur, are thin targets: at n = 4 the 96 hard facets are exceeded by the linear relaxation by only 0.004–0.24 (qhull's unit-normal scaling), and 3,000 random vertices of the relaxation at each of three spectra never land outside P(c) (`relaxation_vertices.py`, `logs/relaxation_vertices_n4.txt`). So N1 is evidence, not a near-proof, and it does not replace a hard-facet enumeration.

**N2 (numerical, reproduces a known count). Easy and hard facets.** Write D(U) for the vector of all squared minors |det U[S,T]|². Then m(U D_c U*) = M_c·D(U) + const. A facet of P(c) is *easy* if its pullback is a nonnegative combination of squared minors modulo the affine relations of D(U), and *hard* otherwise. With all 184 qhull facets at n = 4 this Farkas test gives exactly 96 hard facets at each of four spectra, matching R3–R4 of [`../pmsh/README.md`](../pmsh/README.md).

**N3 (negative structural result; LP-exact up to floating point). At n = 4 the hard facets cannot be derived linearly from PMSH(3).**
- **Constraints inherited from dimension n − 1.** For every position t and real γ, PMSH(n − 1) forces:
  - (a) the projection of m(H) onto subsets not containing t lies in the same projection of P(c). Proof: the compression deleting t has spectrum s interlacing c; PMSH(n − 1) puts its minors in P(s), and the coupling Lemma C of `../hermitian-lifting/` writes every v(s) as a mixture of v_τ(c) over injections τ;
  - (b) the projection of m(H + γI) onto subsets containing t lies in the same projection of P(c + γ). Proof: Jacobi's identity det H[T,T] = det H·det(H⁻¹)[Tᶜ,Tᶜ] and (a) for (H + γI)⁻¹.
- **Test (`inherited_n4.py`).** Add all of (a) for t = 1..4 and (b) for t = 1..4 and eight values of γ to the nonnegativity relaxation, then maximize each of the 184 facet functionals.
- **Outcome.** At c = (1,2,3,4), (0.3,1.1,2,3.7), (1,1.01,5,9) and (1,1.27,1.64,2), **all 96 hard facets are still exceeded, by exactly the same amount** (`logs/inherited_n4.txt`). The points of the relaxation that leave P(c) already satisfy every inherited constraint.
- **Control.** Replacing the vertex sets of the constraints (a) by 12 of the 24 permutations removes all excess (max 1.2·10⁻¹¹). So the constraints are imposed; they are simply not binding.
- **Consequence.** A proof of PMSH(n) by linear programming from PMSH(n − 1), marginals, Palm and shifted Palm conditionings and nonnegativity of squared minors fails already from n = 3 to n = 4. A genuinely four-dimensional nonlinear inequality, the star inequality (★) of [`../pmsh/equivalence.md`](../pmsh/equivalence.md), is needed there. The same will be needed at n = 5 unless the hard facets of P₅(c) turn out to be implied by (★) on 4×4 compressions; that is the first open question below.

**N4 (numerical). Splitting along one row.** For H = U D_c U* and a row r, put ρ_j = |U_{rj}|². A *split* is a family of Hermitian H_j of size n − 1 with spectrum c ∖ c_j such that Σ_j ρ_j m(H_j) = m(H with row and column r deleted) and Σ_j ρ_j c_j m(H_j) = (det H[S ∪ r])_S. A split plus PMSH(n − 1) for each H_j gives a PMSH certificate for H, by conditioning on the eigenvalue placed at r (Proposition H3 forces the marginal ρ). The obvious canonical choices of H_j fail (`../minors/pmsh.md`, negative results); `split_test.py` instead solves for the H_j by least squares. Results so far: a split along the first row exists to residual below 10⁻¹¹ for every random H tried (`logs/split_n4.txt`, `logs/split_n5.txt`, `logs/split_n6.txt`; to be completed). This is a reformulation, not a shortcut: at n = 5 it is equivalent to PMSH(5) if 4×4 Hermitian matrices with a fixed spectrum fill P₄, and for large n the number of conditions (2ⁿ − 2) exceeds the number of unknowns (n(n − 1)(n − 2)) from n = 9 on.

## Pitfalls met along the way (kept so others avoid them)

- **Joint ascent stalls.** Maximizing ℓ·m(U) − max_π ℓ·v_π(c) jointly over (U, ℓ) converges to a permutation matrix, where the gap is exactly 0: all 48 runs at n = 5 did. A gap of 0 from such a search says nothing.
- **Ill-conditioned spectra.** With c = (0, 0.001, 0.01, 0.1, 1) unshifted, the ray LP returned "facets" with entries near 10²¹ that are invalid on the polytope's own vertices, and the search reported values up to 135. Re-checking the same four matrices after c ↦ c + 1 put them strictly inside (ray depth 1.03–1.04, L1 distance 0). The scripts now shift, validate facets, and re-check candidates.
- **Random relaxation vertices miss hard facets.** See N1.

## Open questions, in order

1. **Enumerate the hard facets of P₅(c)** at a few rational spectra and classify them under S₅. qhull is out of reach in dimension 26; pycddlib could not be built here. A targeted harvest (rays toward perturbations of the (★)-tight 4×4 configurations embedded as U₄ ⊕ 1) is the next attempt.
2. **Are the hard facets of P₅ implied by (★) on 4×4 compressions and Schur complements?** If yes, PMSH(5) reduces to the n = 4 inequality; if no, write down the new inequalities, as Theorem EQ did for n = 4.
3. **Contacts.** Are hard facets of P₅ touched only at block-reducible U (4+1, 3+2), or also at irreducible U? Strict inequality at irreducible U would localize the difficulty.

## Reproduction

numpy 2.0.2, scipy 1.13.1 (HiGHS), Python 3.9. Run from this directory.

```sh
python3 facet_search.py 4 0 600 4                 # control, about 20 s
python3 facet_search.py 5 11 1500 6               # one n = 5 spectrum, 10-30 min
PY=python3 ./run_spectra.sh                       # special n = 5 spectra
python3 inherited_n4.py 1,2,3,4 -2.5,-1.5,-0.5,0.5,1.5,-4.5,7
python3 inherited_n4.py 1,2,3,4 none --control
python3 split_test.py 5 1 30 12
```
