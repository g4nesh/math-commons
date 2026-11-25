# PMSH(4): reduction to one inequality on 4×4 unitaries

- **Task:** `mdo-pmsh-4`. Agent `ganesh-claude-minors` (Claude Code, model `claude-opus-5-5`), GitHub owner `g4nesh`. November 25, 2026.
- **Status:** **PMSH(4) is NOT proved here.** This is a bounded contribution that narrows it to an explicit, finite obligation and preserves the negative results.
- **Evidence levels:** marked on each item. No independent review yet.

## Update, second pass (same day): the spectrum is eliminated

**PMSH(4) is still not proved.** It is now equivalent to a single spectrum-free inequality. Full statement, proof and evidence: [equivalence.md](equivalence.md). Scripts: [`star/`](star/). Run logs: [`star/logs/`](star/logs/).

- **Theorem EQ (computer-assisted, unreviewed).** PMSH(4) holds for all Hermitian matrices ⟺ for every U ∈ U(4)

  D[12,34] ≤ D[13,13] + D[14,14] + D[13,23] + (√D[12,12] + √D[13,14])²,   (★)

  where D[S,T] = |det U[S,T]|².
- **How it supersedes the first pass.**
  - R3's "four inequality families at a fixed spectrum" are, for every spectrum, AM–GM consequences of (★) and its 288 symmetric images. Each hard facet is (exotic facet) + α·D_a + α′·D_b with (α−1)α′ = 1.
  - R4 (same structure at every spectrum) is now established for **all** ratios by an exact enumeration on the 42 intervals between 41 critical ratios: always 184 facets = 88 implied by positivity + 96 hard.
- **Evidence that (★) is true (numerical only).**
  - Direct minimization over U(4) gives minimum 0.
  - Numerical sum-of-squares certificates: degree 3 covers the regime D[12,12] ≥ D[13,14] with one certificate; degree 4 extends further but not uniformly.
  - By-product: the degree-3 bound for the exotic facet equals Drury's constant −0.067442.
- **Negative results kept.** First-level relaxation over O(6) is far too weak; a naive parameter-free strengthening is false; the square-root cross term cannot be dropped.
- **Handoff.** The remaining obligation is (★) itself. See "What would finish it" in [equivalence.md](equivalence.md). The handoff list at the end of this file is superseded by it, except item 1 (unitary-to-orthogonal reduction), which is still open and still useful.

The first-pass results R1–R7 below are unchanged. They were written by the same agent under model `claude-opus-5-5`; the second pass used `claude-fable-5-2`.

**Background.** PMSH(n) is defined in [`../minors/pmsh.md`](../minors/pmsh.md). For Hermitian H with spectrum c, the principal-minor vector m(H) = (det H[T,T])_T lies in conv{v_π(c) = (∏_{t∈T} c_{π(t)})_T : π ∈ Sₙ}. By Theorem H1 there, PMSH(4) implies MOC for n = 4 whenever one matrix is Hermitian. With Shitov's Theorem 19 that extends to Hermitian A with ≤ 4 nonreal eigenvalues of B, in every dimension.

**Notation.** For U ∈ U(4):
- x(U) := (B, D), where B_{sj} = |u_{sj}|² (16 entries) and D_{S,T} = |det U[S,T]|² for |S| = |T| = 2 (36 entries).
- Q := conv{x(P_π) : π ∈ S₄}, the 82-facet compound polytope.
- Φ_a is the linear map x ↦ (Σ_{|S|=|T|} a^{[4]∖S} x_{S,T})_T. Its image is m(H′) up to a positive factor, for H′ = U*D_{1/a}U.
- A *PMSH certificate* at (U, a) exists iff Φ_a(x(U)) ∈ Φ_a(Q).

## Results

**R1 (rigorous). The compound vector of every unitary lies in aff(Q).**
- **The identities.** Every U ∈ U(4) satisfies:
  - (i) B is doubly stochastic;
  - (ii) Σ_T D_{S,T} = Σ_S D_{S,T} = 1 (Cauchy–Binet);
  - (iii) Σ_{T∋j} D_{S,T} = Σ_{s∈S} B_{sj} and Σ_{S∋i} D_{S,T} = Σ_{t∈T} B_{it} (Lemma 3 of [`../minors/report.md`](../minors/report.md), with |R| = 1);
  - (iv) D_{S,T} = D_{Sᶜ,Tᶜ} (Jacobi).
- **Exact rank.** In exact integer arithmetic (`exact_rank.py`), these 64 equations have rank 39, so their solution set is a 13-dimensional affine space. It contains all 24 vertices, and dim aff(Q) = 13, so the solution set **equals** aff(Q).

**R2 (rigorous). Matching reformulation.**
- **Forced marginal.** By Proposition H3 of `../minors/pmsh.md`, a certificate λ is a Birkhoff decomposition of B. Modulo that, only two scalar conditions remain.
- **The two conditions.** Let the three perfect matchings of positions be M = {st|uv}, and those of eigen-indices be N = {jk|lm}. Put p_N := c_jc_k + c_lc_m and N_U(M,N) := |det U[S,T]|² + |det U[S,Tᶜ]|² (S ∈ M, T ∈ N). The conditions are
  E_λ[p_{π(M)}] = Σ_N N_U(M,N) p_N = det H[st] + det H[uv] for every M.
- **Proof.** Cauchy–Binet plus Jacobi complementarity. Given B, the six 2×2 conditions satisfy four linear relations (the incidence matrix of K₄ has rank 4), leaving two.
- **One parameter.** With c sorted, p_{12|34} > p_{13|24} > p_{14|23}, because the differences are (c₁−c₄)(c₂−c₃) and (c₁−c₂)(c₃−c₄). Doubly stochastic N ignores constant shifts of p, so the spectrum enters only through the ratio **r = (p_{12|34} − p_{13|24})/(p_{13|24} − p_{14|23}) ∈ (0, ∞)**.

**R3 (exact facets; Farkas classification by floating LP). At a fixed spectrum, PMSH(4) is equivalent to four inequalities over U(4).**
- **Facets.** At c = (3,4,6,12), the projected polytope has **184 facets in 16 S₄-orbits**, with exact rational normals (`pmsh_facets.py`).
- **The easy 88.** These are the orbits of sizes 4 (18 vertices: exactly |u_{tj}|² ≥ 0), 6 (16 vertices) and 12 (14 vertices). Each pulls back to a **nonnegative combination of Q's positivity facets** plus the R1 equations (Farkas LP, `classify.py`), so it holds for every unitary.
- **The hard 96.** These are the four orbits of 24 **simplicial** facets (11 vertices each). Every Farkas representation needs an exotic facet of Q, with minimal exotic weight 0.2957, 0.1634, 0.1545 or 0.1473. So at this spectrum PMSH(4) ⟺ four explicit polynomial inequalities over U(4), the generalized Fiedler inequalities. One orbit representative:
  `576m₁ − 576m₃ − 88m₁₂ − 88m₁₃ − 232m₁₄ + 8m₂₃ + 248m₂₄ + 152m₃₄ + 41m₁₂₃ − 19m₁₂₄ + 29m₁₃₄ − 51m₂₃₄ ≤ 2112`.
  The other three are printed by `python3 pmsh_facets.py 3,4,6,12`.

**R4 (computational). The structure looks uniform in the spectrum.**
- **Counts.** The projected polytope has 184 facets for all 10 spectra of `qtriv.py`, for a = (1,2,3,4), and for all 10 spectra completed in the sweep (`sweep_*.log`).
- **Classification.** For a = (1,2,3,4) and for 9 of the 10 sweep spectra, the Farkas classification gives exactly the 96 simplicial facets as the hard ones.
- **Maxima.** Every hard inequality has numerical max over U(4) ≤ 0, attained at permutation matrices. The largest value is −3·10⁻¹⁴. This used only 4 BFGS restarts per facet in the sweep (25 for a = (1,2,3,4)); the 192-restart adversarial search in `../minors/pmsh.md` is the stronger test.
- **Caveats (numerical failures, not mathematical observations).**
  - The 10th sweep spectrum was near-degenerate: two a-values within ~10⁻³ of each other and of 0. There the classification reported zero hard facets and a 12-vertex facet type, but the LP tolerances are meaningless at that scale. In the exact degenerate limit, PMSH holds by Proposition H4 anyway.
  - One sweep worker crashed with a qhull precision error at a clustered spectrum.
  - The remaining sweep workers were stopped after these 10 spectra.

**R5 (negative results, computational).**
- **Positivity alone is insufficient.** Let Q_triv := aff(Q) ∩ {B ≥ 0, D ≥ 0}. Then Φ_a(Q_triv) exceeds Φ_a(Q) by up to 1.5 (`qtriv.py`).
- **Uniform relaxation fails.** Relaxing all 48 exotic facets uniformly by any δ > 0 already enlarges the projection: δ_crit = 0 for all 12 spectra in `qdelta.py`. The reason is that a uniform relaxation adds points near the facets' vertices, and unitaries never reach those points (the caps stay ≥ 0.37 from every vertex). A proof must use the shape of the caps, not only Drury's depth bound.

**R6 (rigorous identity, verified numerically). SO(4) structure.**
- **The identity.** For O ∈ SO(4), write Λ²ℝ⁴ = Λ²₊ ⊕ Λ²₋ in the bases f^±_M = (e_S ± ε_S e_{Sᶜ})/√2, indexed by matchings, where ε_S is the Hodge sign. Then C₂(O) = R₊ ⊕ R₋ with R± ∈ SO(3), and
  **N_O(M,N) = ½(R₊(M,N)² + R₋(M,N)²)**.
- **Why it holds.** C₂(O) commutes with the Hodge star, and det O[S,T] = ½(R₊ ± R₋)_{MN}.
- **Check.** Error 1·10⁻¹⁵ over 20 random O.
- **So:** the matching data PMSH constrains is an average of two 3×3 orthostochastic matrices.

**R7 (computational, potentially decisive). conv{x(U)} = conv{x(O)} numerically.**
- **Evidence.** For 30 random linear functionals ℓ on the compound coordinates (seeds 1–3, 10 each), max over U(4) of ℓ equals max over O(4) of ℓ; the largest gap max_U − max_O is −1.3·10⁻¹⁵, so U never exceeds O (`conv_ou.py`; 30 BFGS restarts per functional on U(4), 60 on O(4)). The deepest exotic violation is also the same for both (`../minors/drury-constant.md`).
- **If proved:** every unitary's compound vector would be a convex combination of real orthogonal ones. PMSH(4), and MOC(4) itself, would then reduce to O(4), a 6-dimensional problem in (q_L, q_R) ∈ S³ × S³. This is consistent with the "unitary-to-orthogonal reduction" that `../literature/` attributes to Zheng's Theorem 1.2. It is evidence for that step, **not** a verification of it.

## Handoff: what would finish PMSH(4)

1. **Prove R7** (or audit Zheng's Theorem 1.2) for n = 4. This reduces everything to real orthogonal matrices.
2. **Prove the four simplicial-facet inequality families** on O(4), uniformly in r ∈ (0,∞).
   - By R6, every quantity is a polynomial in two SO(3) rotations and the orthostochastic matrix O∘O.
   - Each inequality is tight at 11 permutation matrices, so a certified proof needs second-order analysis there plus a global bound elsewhere: SOS, or interval branch-and-bound on S³ × S³.
3. **Prove R4 uniformity.** Show the facet combinatorics of Φ(Q) do not change with r, so there really are only four families.

## Reproduction

Python 3.9.6, numpy 2.0.2, scipy 1.13.1 (HiGHS), sympy 1.14.0. Run from this directory:

```sh
python3 exact_rank.py                 # R1, exact
python3 pmsh_facets.py 3,4,6,12       # R3, exact facet normals by orbit
python3 classify.py 1,2,3,4           # R3, Farkas classification + numerical max over U(4)
python3 qtriv.py; python3 qdelta.py   # R5 negative results
python3 sweep.py 1 6 4                # R4: seeds 1..6 were started with these arguments; 10 spectra completed (sweep_*.log) before the workers were stopped; worker 4 crashed (qhull)
python3 conv_ou.py 1 10               # R7 (seeds 1..3)
```

Run logs are in [`logs/`](logs/) as `.txt` files (the repository ignores `*.log`, so until the second-pass commits they were missing from the repo). Solver warnings and local file paths were stripped from them. The sweep workers were stopped before writing JSON, so their logs are the record.

**Compute:** about 2 process-hours, local (up to 9 parallel processes), within the 3 process-hours budgeted in this task's first checkpoint. No cloud or paid compute.
