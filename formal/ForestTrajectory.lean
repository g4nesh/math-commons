import ForestCut

/-!
Finite occupied subsets are Boolean vertex functions. Each indexed gate either
preserves a subset or exchanges two endpoints. The finite vertex enumeration is
split around each edge, so occupancy counts count each endpoint exactly once.
The general uniqueness result below does not identify labeled particles: when
both endpoints have the same occupancy, a swap is the same subset.

This file concerns the combinatorial bridge only. No complex matrices,
exterior powers, Cauchy–Binet, amplitudes, or convex weights are formalized here.
-/
namespace MathCommons.ForestTrajectory

variable {V : Type} [DecidableEq V]

/-- An edge with an explicit decomposition of the finite vertex enumeration.
For a genuine set enumeration the caller supplies a duplicate-free `vertices`.
The arithmetic proofs only need that the two endpoints occur exactly once. -/
structure Edge (vertices : List V) where
  left : V
  right : V
  distinct : left ≠ right
  rest : List V
  enumerate : vertices.Perm (left :: right :: rest)
  avoids : ∀ x ∈ rest, x ≠ left ∧ x ≠ right

abbrev Subset (V : Type) := V → Bool

def swapSubset {vertices : List V} (e : Edge vertices) (s : Subset V) : Subset V :=
  fun x => if x = e.left then s e.right else if x = e.right then s e.left else s x

def cutCount (vertices : List V) (cut s : Subset V) : Nat :=
  vertices.countP (fun x => cut x && s x)

def Step {vertices : List V} (e : Edge vertices) (s t : Subset V) : Prop :=
  t = s ∨ t = swapSubset e s

theorem swapAtLeft {vertices : List V} (e : Edge vertices) (s : Subset V) :
    swapSubset e s e.left = s e.right := by simp [swapSubset]

theorem swapAtRight {vertices : List V} (e : Edge vertices) (s : Subset V) :
    swapSubset e s e.right = s e.left := by
  simp [swapSubset, Ne.symm e.distinct]

/-- Covers both double occupation and double vacancy. -/
theorem swapCollision {vertices : List V} (e : Edge vertices) (s : Subset V)
    (same : s e.left = s e.right) : swapSubset e s = s := by
  funext x
  by_cases left : x = e.left
  · subst x; simpa [swapSubset] using same.symm
  · by_cases right : x = e.right
    · subst x; simpa [swapSubset, Ne.symm e.distinct] using same
    · simp [swapSubset, left, right]

theorem swapEmpty {vertices : List V} (e : Edge vertices) :
    swapSubset e (fun _ => false) = (fun _ => false) := swapCollision e _ rfl

theorem swapFull {vertices : List V} (e : Edge vertices) :
    swapSubset e (fun _ => true) = (fun _ => true) := swapCollision e _ rfl

omit [DecidableEq V] in
theorem countSplit {vertices : List V} (e : Edge vertices) (cut s : Subset V) :
    cutCount vertices cut s = cutCount e.rest cut s +
      (if cut e.right && s e.right then 1 else 0) +
      (if cut e.left && s e.left then 1 else 0) := by
  unfold cutCount
  rw [e.enumerate.countP_eq, List.countP_cons, List.countP_cons]

theorem restCountSwap {vertices : List V} (e : Edge vertices) (cut s : Subset V) :
    cutCount e.rest cut (swapSubset e s) = cutCount e.rest cut s := by
  apply List.countP_congr
  intro x hx
  simp [swapSubset, (e.avoids x hx).1, (e.avoids x hx).2]

/-- A gate with both endpoints on the same side preserves that cut's count. -/
theorem swapPreservesNoncrossingCut {vertices : List V} (e : Edge vertices)
    (cut s : Subset V) (sameSide : cut e.left = cut e.right) :
    cutCount vertices cut (swapSubset e s) = cutCount vertices cut s := by
  rw [countSplit e, countSplit e, restCountSwap, swapAtLeft, swapAtRight]
  cases hl : cut e.left <;> cases hr : cut e.right <;>
    cases sl : s e.left <;> cases sr : s e.right <;> simp_all

theorem stepPreservesNoncrossingCut {vertices : List V} (e : Edge vertices)
    (cut s t : Subset V) (sameSide : cut e.left = cut e.right) (step : Step e s t) :
    cutCount vertices cut t = cutCount vertices cut s := by
  rcases step with rfl | rfl
  · rfl
  · exact swapPreservesNoncrossingCut e cut s sameSide

/-- If the cut separates the edge endpoints, an invisible swap is a collision. -/
theorem equalCountSwapCollision {vertices : List V} (e : Edge vertices)
    (cut s : Subset V) (separates : cut e.left ≠ cut e.right)
    (sameCount : cutCount vertices cut s = cutCount vertices cut (swapSubset e s)) :
    s e.left = s e.right := by
  rw [countSplit e, countSplit e, restCountSwap, swapAtLeft, swapAtRight] at sameCount
  cases hl : cut e.left <;> cases hr : cut e.right <;>
    cases sl : s e.left <;> cases sr : s e.right <;> simp_all <;> omega

/-- At a fixed current subset, the next subset is determined by its cut count. -/
theorem nextSubsetUnique {vertices : List V} (e : Edge vertices)
    (cut s t u : Subset V) (separates : cut e.left ≠ cut e.right)
    (st : Step e s t) (su : Step e s u)
    (sameCount : cutCount vertices cut t = cutCount vertices cut u) : t = u := by
  rcases st with ht | ht <;> rcases su with hu | hu
  · exact ht.trans hu.symm
  · rw [ht, hu] at sameCount ⊢
    exact (swapCollision e s (equalCountSwapCollision e cut s separates sameCount)).symm
  · rw [ht, hu] at sameCount ⊢
    exact swapCollision e s (equalCountSwapCollision e cut s separates sameCount.symm)
  · exact ht.trans hu.symm

/-- Explicit edge-cut certificate: exactly the indexed edge crosses its cut. -/
structure CutCertificate {vertices : List V} (edges : Nat → Edge vertices) (n : Nat) where
  side : Nat → Subset V
  separates : ∀ j, j < n → side j (edges j).left ≠ side j (edges j).right
  other : ∀ j, j < n → ∀ k, k < n → k ≠ j →
    side j (edges k).left = side j (edges k).right

def Trajectory {vertices : List V} (edges : Nat → Edge vertices) (n : Nat)
    (history : Nat → Subset V) : Prop :=
  ∀ k, k < n → Step (edges k) (history k) (history (k + 1))

/-- The one-opportunity hypothesis is derived from actual subset transitions. -/
theorem subsetCountsSingleOpportunity {vertices : List V}
    (edges : Nat → Edge vertices) (n : Nat) (cuts : CutCertificate edges n)
    (history : Nat → Subset V) (path : Trajectory edges n history)
    (j : Nat) (hj : j < n) :
    ForestCut.SingleOpportunity (fun k => (cutCount vertices (cuts.side j) (history k) : Int)) n j := by
  refine ⟨hj, ?_⟩
  intro k hk hne
  exact congrArg (fun value : Nat => (value : Int))
    (stepPreservesNoncrossingCut (edges k) (cuts.side j) _ _
      (cuts.other j hj k hk hne) (path k hk))

/-- Equal endpoint subsets determine every intermediate subset, including
collision, empty, and full cases. This is conditional on explicit cut data. -/
theorem trajectoryUnique {vertices : List V}
    (edges : Nat → Edge vertices) (n : Nat) (cuts : CutCertificate edges n)
    (a b : Nat → Subset V) (pathA : Trajectory edges n a) (pathB : Trajectory edges n b)
    (start : a 0 = b 0) (finish : a n = b n) :
    ∀ k, k ≤ n → a k = b k := by
  intro k
  induction k with
  | zero => intro _; exact start
  | succ k ih =>
    intro within
    have hk : k < n := by omega
    have previous : a k = b k := ih (by omega)
    have counts := ForestCut.cutHistoryUnique
      (fun t => (cutCount vertices (cuts.side k) (a t) : Int))
      (fun t => (cutCount vertices (cuts.side k) (b t) : Int)) n k
      (subsetCountsSingleOpportunity edges n cuts a pathA k hk)
      (subsetCountsSingleOpportunity edges n cuts b pathB k hk)
      (by dsimp only; rw [start]) (by dsimp only; rw [finish]) (k + 1) within
    have naturalCounts : cutCount vertices (cuts.side k) (a (k + 1)) =
        cutCount vertices (cuts.side k) (b (k + 1)) := Int.ofNat.inj counts
    apply nextSubsetUnique (edges k) (cuts.side k) (a k) _ _ (cuts.separates k hk)
    · exact pathA k hk
    · rw [previous]; exact pathB k hk
    · exact naturalCounts

#print axioms swapAtLeft
#print axioms swapAtRight
#print axioms swapCollision
#print axioms swapEmpty
#print axioms swapFull
#print axioms countSplit
#print axioms restCountSwap
#print axioms swapPreservesNoncrossingCut
#print axioms stepPreservesNoncrossingCut
#print axioms equalCountSwapCollision
#print axioms nextSubsetUnique
#print axioms subsetCountsSingleOpportunity
#print axioms trajectoryUnique

end MathCommons.ForestTrajectory
