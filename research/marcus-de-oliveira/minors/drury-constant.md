# Drury's 1992 constant in closed form

- **Objective / task:** `marcus-de-oliveira` / `mdo-minors` (`ganesh-claude-minors`, owner `g4nesh`).
- **Evidence level:**
  - Exact symbolic computation on an explicit one-parameter family.
  - That this family attains the global extremum over U(4) is **computational evidence**, not a proof.
- **Exact artifact commit:** added in 8c1b1fb (branch `claude/mdo-minors-compound-certificates`), on main as 304f790 (PR 20). Items 1 and 5 reworded on 25 Nov 2026 after the independent review [`../reviews/mdo-minors.md`](../reviews/mdo-minors.md), finding 2.

## Statement

Let g be the exotic facet functional of the n = 4 compound polytope ([report](report.md) §5.1):

  g(U) = |u₁₁|² + |u₃₃|² + |det U[{1,4},{3,4}]|² − |det U[{1,2},{1,4}]|²,

which is ≥ 0 on all permutation matrices. Drury's functional is 1 − g. Let S be the skew conference matrix

  S = [[0,1,1,−1],[−1,0,−1,−1],[−1,1,0,1],[1,1,−1,0]],  with Sᵀ = −S and S² = −3I,

and let O(y) := y·I + √((1−y²)/3)·S for y ∈ [−1,1]. Then:

1. **Orthogonal family.** O(y) is orthogonal, and {O(y)} is the half φ ∈ [0, π] of the one-parameter subgroup exp(φS/√3) = cos φ·I + sin φ·S/√3 (the positive root is taken in the definition). The other half gives the same minimum (numerically, on a 2001-point grid by the reviewer).
2. **Closed form on the family.** g(O(y)) = 2y² − (4√3/9)·y·(1−y²)^{3/2}.
3. **Minimum on the family.** The minimum over y is attained where y*² is the root near 0.028 of 16u³ − 24u² + 36u − 1, so y* ≈ 0.168231070580. The minimum value g_min ≈ −0.0674422488681229 is the unique real root of 4G³ − 12G² + 14G + 1.
4. **Drury's constant.** D := 1 − g_min is the real root of **4D³ + 2D − 7 = 0**:

   D = ∛(7/8 + 11√33/72) + ∛(7/8 − 11√33/72) ≈ 1.0674422488681232.

   The radicand is clean because 49/64 + 1/216 = 1331/1728 = (11/12)³.
5. **Zheng's example lies on the family up to symmetry.** g(O(1/7)) = −158/2401. O(1/7) = (I + 4S)/7 is not Zheng's O* itself: O(1/7) = D₁ P_α V P_β D₂ with V ∈ {O*, O*ᵀ}, permutation matrices P and diagonal sign matrices D (24 such representations; exact check in `moc/zheng_on_family.py`). These maps preserve all squared minors up to relabelling, and −158/2401 is the minimum of g over the 1,152 symmetric images of O*.

## Proof of 1–5, and the computational part

1. O(y)ᵀO(y) = y²I + ((1−y²)/3)·SᵀS + y·√((1−y²)/3)·(S + Sᵀ) = I.
2. Items 2–5 are exact sympy computations in `moc/drury_exact.py`:
   - item 2 by expanding the minors;
   - item 3 by the derivative and a resultant;
   - item 4 by substituting G = 1 − D, then Cardano's formula.
3. **Global optimality (numerical only).** min_{U∈U(4)} g(U) and min_{O∈O(4)} g(O) both equal −0.0674422489 (`moc/orth_vs_unitary.py`, seed 2026). The two searches were 200 and 400 BFGS + Nelder–Mead restarts, and they agree to 10⁻¹⁰.
4. **Consistency with the literature.** This matches Drury's published 1.06744.

## Conclusion and limits

The extremal violation of the compound-Birkhoff property for n = 4 is attained, numerically, by real orthogonal matrices. So "unitary vs orthogonal depth" gives **no** obstruction to any orthogonal-endpoint representation of unitary compounds. This is relevant to auditing Zheng's Theorem 1.2, but nothing is claimed about that theorem here.

The closed form is exact on the family. Proving that the family is globally extremal over U(4) remains open.

## Independent review

None yet.
