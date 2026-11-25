import Std

/-!
# The unique crossing arithmetic in the forest-rotation argument

This file proves statements for arbitrary finite histories of integer cut counts.
If the count is unchanged at every step except one specified opportunity, its
signed change at that opportunity is determined by its two endpoint counts.
Two such histories with equal endpoints have the same crossing decision.

Connection to the written argument in
`research/marcus-de-oliveira/baselines/README.md`, section 2:
take `f k` to be the number of occupied vertices on one side of a deleted forest
edge after the first `k` gates. The cut is crossed by only the gate on that edge.

**Not formalized here:** graph/forest definitions; deleting an edge produces a
cut crossed by no other edge; subset cardinalities give these integer histories;
the permitted gate transitions give `SingleOpportunity` and `Crossing`; all cut
decisions together determine a subset path; exterior powers, Cauchy–Binet,
complex unitary matrices, and the resulting convex combination. The hypotheses
below are explicit, not axioms standing in for those bridges. This file alone
does not formalize the forest-product theorem or the general conjecture.
-/

namespace MathCommons.ForestCut

/-- Sum of the first `n` signed changes of a history, defined recursively. -/
def totalChange (f : Nat → Int) : Nat → Int
  | 0 => 0
  | n + 1 => totalChange f n + (f (n + 1) - f n)

/-- Telescoping holds for every history and every finite length. -/
theorem telescoping (f : Nat → Int) (n : Nat) :
    totalChange f n = f n - f 0 := by
  induction n with
  | zero => simp [totalChange]
  | succ n ih =>
    simp only [totalChange, ih]
    omega

/-- A history with no changes on an interval has equal endpoint values.
This lemma does not require integer-valued histories. -/
theorem constantOnInterval {α : Type} (f : Nat → α) (a b : Nat)
    (hab : a ≤ b)
    (stable : ∀ k, a ≤ k → k < b → f (k + 1) = f k) :
    f b = f a := by
  induction b with
  | zero =>
    have ha : a = 0 := by omega
    subst a
    rfl
  | succ b ih =>
    by_cases before : a ≤ b
    · have last : f (b + 1) = f b := stable b before (by omega)
      have previous : f b = f a := ih before (fun k lo hi => stable k lo (by omega))
      exact last.trans previous
    · have ha : a = b + 1 := by omega
      subst a
      rfl

/-- Among the first `n` transitions, only transition `j` may change the count.
The transition is from time `j` to time `j+1`; it need not actually change. -/
def SingleOpportunity (f : Nat → Int) (n j : Nat) : Prop :=
  j < n ∧ ∀ k, k < n → k ≠ j → f (k + 1) = f k

/-- The unique possible jump equals the difference of the endpoint counts. -/
theorem singleCutChange (f : Nat → Int) (n j : Nat)
    (one : SingleOpportunity f n j) :
    f (j + 1) - f j = f n - f 0 := by
  have before : f j = f 0 :=
    constantOnInterval f 0 j (by omega)
      (fun k _ hi => one.2 k (by have := one.1; omega) (by omega))
  have after : f n = f (j + 1) :=
    constantOnInterval f (j + 1) n (by have := one.1; omega)
      (fun k lo hi => one.2 k hi (by omega))
  rw [before, after]

/-- In the one-opportunity case, the telescoping sum is that one jump. -/
theorem totalChangeAtOpportunity (f : Nat → Int) (n j : Nat)
    (one : SingleOpportunity f n j) :
    totalChange f n = f (j + 1) - f j := by
  exact (telescoping f n).trans (singleCutChange f n j one).symm

/-- Two histories with the same endpoint counts have the same jump, provided
each can change at only the specified common opportunity. -/
theorem sameEndpointsSameCrossing (f g : Nat → Int) (n j : Nat)
    (oneF : SingleOpportunity f n j) (oneG : SingleOpportunity g n j)
    (start : f 0 = g 0) (finish : f n = g n) :
    f (j + 1) - f j = g (j + 1) - g j := by
  rw [singleCutChange f n j oneF, singleCutChange g n j oneG, start, finish]

/-- A stronger cut-level conclusion: the complete count histories agree
through the final time. This does not assert equality of the occupied subsets. -/
theorem cutHistoryUnique (f g : Nat → Int) (n j : Nat)
    (oneF : SingleOpportunity f n j) (oneG : SingleOpportunity g n j)
    (start : f 0 = g 0) (finish : f n = g n)
    (k : Nat) (within : k ≤ n) : f k = g k := by
  by_cases before : k ≤ j
  · have fBefore : f k = f 0 :=
      constantOnInterval f 0 k (by omega)
        (fun t _ hi => oneF.2 t (by have := oneF.1; omega) (by omega))
    have gBefore : g k = g 0 :=
      constantOnInterval g 0 k (by omega)
        (fun t _ hi => oneG.2 t (by have := oneG.1; omega) (by omega))
    exact fBefore.trans (start.trans gBefore.symm)
  · have fAfter : f n = f k :=
      constantOnInterval f k n within
        (fun t lo hi => oneF.2 t hi (by omega))
    have gAfter : g n = g k :=
      constantOnInterval g k n within
        (fun t lo hi => oneG.2 t hi (by omega))
    exact fAfter.symm.trans (finish.trans gAfter)

/-- The three possible occupation-count decisions for a two-coordinate gate. -/
inductive Crossing where
  | stay
  | enter
  | leave
  deriving DecidableEq

/-- Signed count change on a chosen side of the cut. -/
def crossingDelta : Crossing → Int
  | .stay => 0
  | .enter => 1
  | .leave => -1

/-- The signed changes distinguish the three crossing decisions. -/
theorem crossingDeltaInjective (x y : Crossing)
    (same : crossingDelta x = crossingDelta y) : x = y := by
  cases x <;> cases y <;> simp_all [crossingDelta]

/-- Equal endpoints and one crossing opportunity force the same stay/enter/leave
decision. The hypotheses link the decisions to the actual history changes. -/
theorem crossingDecisionUnique (f g : Nat → Int) (n j : Nat)
    (oneF : SingleOpportunity f n j) (oneG : SingleOpportunity g n j)
    (start : f 0 = g 0) (finish : f n = g n)
    (x y : Crossing)
    (representsF : f (j + 1) - f j = crossingDelta x)
    (representsG : g (j + 1) - g j = crossingDelta y) : x = y := by
  apply crossingDeltaInjective
  exact representsF.symm.trans
    ((sameEndpointsSameCrossing f g n j oneF oneG start finish).trans representsG)

#print axioms telescoping
#print axioms constantOnInterval
#print axioms singleCutChange
#print axioms totalChangeAtOpportunity
#print axioms sameEndpointsSameCrossing
#print axioms cutHistoryUnique
#print axioms crossingDeltaInjective
#print axioms crossingDecisionUnique

end MathCommons.ForestCut
