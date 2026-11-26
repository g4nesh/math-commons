import ForestPathProbability

/-!
Finite normalization of the existing gate probabilities and path products.
Successors are deduplicated in the equal-occupancy sectors. Histories are
enumerated by extending their next state. No matrix-product or Cauchy–Binet
identity is assumed here; all scalars remain Gaussian rationals / rationals.
-/
namespace MathCommons.ForestPathNormalization

open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability

def sumWith {α : Type} (f : α → Rat) : List α → Rat
  | [] => 0
  | x :: xs => f x + sumWith f xs

theorem sumWithAppend {α : Type} (f : α → Rat) (xs ys : List α) :
    sumWith f (xs ++ ys) = sumWith f xs + sumWith f ys := by
  induction xs with
  | nil => simp only [List.nil_append, sumWith]; grind
  | cons x xs ih => simp only [List.cons_append, sumWith, ih]; grind

theorem sumWithMap {α β : Type} (f : β → Rat) (g : α → β) (xs : List α) :
    sumWith f (xs.map g) = sumWith (fun x => f (g x)) xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp [sumWith, ih]

theorem sumWithFlatMap {α β : Type} (f : β → Rat) (g : α → List β) (xs : List α) :
    sumWith f (xs.flatMap g) = sumWith (fun x => sumWith f (g x)) xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp [List.flatMap_cons, sumWithAppend, sumWith, ih]

theorem sumWithScale {α : Type} (f : α → Rat) (c : Rat) (xs : List α) :
    sumWith (fun x => c * f x) xs = c * sumWith f xs := by
  induction xs with
  | nil => simp [sumWith]
  | cons x xs ih => simp only [sumWith, ih]; grind

theorem sumWithCongr {α : Type} (f g : α → Rat) (xs : List α)
    (same : ∀ x, f x = g x) : sumWith f xs = sumWith g xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp [sumWith, same, ih]

section Paths
variable {V : Type} [DecidableEq V] {vertices : List V}

/-- Distinct possible next subsets, counting a collision stay only once. -/
def successors (e : Edge vertices) (s : Subset V) : List (Subset V) :=
  if s e.left = s e.right then [s] else [s, swapSubset e s]

theorem successorMass (e : Edge vertices) (b : Gaussian) (s : Subset V) :
    sumWith (stepProb e b s) (successors e s) = 1 := by
  have stay := spectatorsAgreeOfStep e s s (Or.inl rfl)
  have swap := spectatorsAgreeOfStep e s (swapSubset e s) (Or.inr rfl)
  cases hl : s e.left <;> cases hr : s e.right <;>
    simp [successors, sumWith, stepProb, stay, swap, swapAtLeft, swapAtRight,
      localProb, hl, hr] <;> grind

theorem stepProbabilityNonnegative (e : Edge vertices) (a b c d : Gaussian)
    (unitary : UnitaryBlock a b c d) (s t : Subset V) : 0 ≤ stepProb e b s t := by
  obtain ⟨lower, upper⟩ := unitaryProbabilityBounds a b c d unitary
  unfold stepProb
  split
  · cases s e.left <;> cases s e.right <;> cases t e.left <;> cases t e.right <;>
      simp_all [localProb] <;> grind
  · decide

theorem successorLegal (e : Edge vertices) (s t : Subset V)
    (member : t ∈ successors e s) : Step e s t := by
  unfold successors at member
  split at member <;> simp only [List.mem_cons, List.not_mem_nil, or_false] at member
  · exact Or.inl member
  · exact member

theorem everyLegalSuccessor (e : Edge vertices) (s t : Subset V)
    (step : Step e s t) : t ∈ successors e s := by
  unfold successors
  split
  · rename_i collision
    rcases step with rfl | rfl
    · simp
    · simp [swapCollision e s collision]
  · rcases step with rfl | rfl <;> simp

theorem successorsNodup (e : Edge vertices) (s : Subset V) : (successors e s).Nodup := by
  unfold successors
  split
  · simp
  · rename_i different
    have ne : s ≠ swapSubset e s := by
      intro eq
      have atLeft := congrFun eq e.left
      rw [swapAtLeft] at atLeft
      exact different atLeft
    simp [ne]

/-- Change exactly the newly appended time, leaving the existing prefix fixed. -/
def extendHistory (history : Nat → Subset V) (n : Nat) (next : Subset V) : Nat → Subset V :=
  fun k => if k = n + 1 then next else history k

theorem pathProbabilityPrefixCongr (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (x y : Nat → Subset V) (n : Nat) (same : ∀ k, k ≤ n → x k = y k) :
    pathProbability edges bs x n = pathProbability edges bs y n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [pathProbability]
    rw [ih (fun k hk => same k (by omega)), same n (by omega), same (n+1) (by omega)]

theorem extendedPathProbability (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (history : Nat → Subset V) (n : Nat) (next : Subset V) :
    pathProbability edges bs (extendHistory history n next) (n+1) =
      pathProbability edges bs history n * stepProb (edges n) (bs n) (history n) next := by
  simp only [pathProbability]
  have samePrefix : ∀ k, k ≤ n → extendHistory history n next k = history k := by
    intro k hk
    simp [extendHistory, show k ≠ n+1 by omega]
  rw [pathProbabilityPrefixCongr edges bs _ history n samePrefix]
  simp [extendHistory]

/-- All legal finite histories from `start`; values after the final time stay at `start`.
This enumerates subset histories rather than labeled stay/swap coin outcomes. -/
def histories (edges : Nat → Edge vertices) (start : Subset V) : Nat → List (Nat → Subset V)
  | 0 => [fun _ => start]
  | n+1 => (histories edges start n).flatMap fun history =>
      (successors (edges n) (history n)).map fun next => extendHistory history n next

/-- Every enumerated history starts at the requested subset and takes legal steps. -/
theorem historyLegal (edges : Nat → Edge vertices) (start : Subset V) (n : Nat)
    (history : Nat → Subset V) (member : history ∈ histories edges start n) :
    history 0 = start ∧ ∀ k, k < n → Step (edges k) (history k) (history (k+1)) := by
  induction n generalizing history with
  | zero =>
    simp only [histories, List.mem_cons, List.not_mem_nil, or_false] at member
    subst history
    exact ⟨rfl, fun k hk => by omega⟩
  | succ n ih =>
    obtain ⟨old, oldMember, outputMember⟩ := List.mem_flatMap.mp member
    obtain ⟨next, nextMember, rfl⟩ := List.mem_map.mp outputMember
    obtain ⟨initial, legal⟩ := ih old oldMember
    constructor
    · simpa [extendHistory] using initial
    · intro k hk
      by_cases last : k = n
      · subst k
        simpa [extendHistory] using successorLegal (edges n) (old n) next nextMember
      · have earlier : k < n := by omega
        simpa [extendHistory, last, show k ≠ n+1 by omega, show k+1 ≠ n+1 by omega]
          using legal k earlier

/-- Every legal prefix has a representative in the finite enumeration. Future
values of an arbitrary input history are deliberately irrelevant. -/
theorem everyLegalPrefix (edges : Nat → Edge vertices) (start : Subset V) (n : Nat)
    (history : Nat → Subset V) (initial : history 0 = start)
    (legal : ∀ k, k < n → Step (edges k) (history k) (history (k+1))) :
    ∃ generated ∈ histories edges start n, ∀ k, k ≤ n → generated k = history k := by
  induction n with
  | zero =>
    refine ⟨fun _ => start, by simp [histories], ?_⟩
    intro k hk
    have : k = 0 := by omega
    subst k
    exact initial.symm
  | succ n ih =>
    obtain ⟨old, oldMember, same⟩ := ih (fun k hk => legal k (by omega))
    have nextMember : history (n+1) ∈ successors (edges n) (old n) := by
      apply everyLegalSuccessor
      rw [same n (by omega)]
      exact legal n (by omega)
    refine ⟨extendHistory old n (history (n+1)), ?_, ?_⟩
    · apply List.mem_flatMap.mpr
      refine ⟨old, oldMember, ?_⟩
      exact List.mem_map.mpr ⟨history (n+1), nextMember, rfl⟩
    · intro k hk
      by_cases last : k = n+1
      · simp [extendHistory, last]
      · simp [extendHistory, last, same k (by omega)]

/-- Enumerated histories have a fixed, harmless tail, so equality of their finite
prefixes is equality of the represented functions. -/
theorem historyFuture (edges : Nat → Edge vertices) (start : Subset V) (n : Nat)
    (history : Nat → Subset V) (member : history ∈ histories edges start n)
    (k : Nat) (future : n < k) : history k = start := by
  induction n generalizing history with
  | zero =>
    simp only [histories, List.mem_cons, List.not_mem_nil, or_false] at member
    subst history
    rfl
  | succ n ih =>
    obtain ⟨old, oldMember, outputMember⟩ := List.mem_flatMap.mp member
    obtain ⟨next, _, rfl⟩ := List.mem_map.mp outputMember
    simp [extendHistory, show k ≠ n+1 by omega, ih old oldMember (by omega)]

theorem historyDeterminedByPrefix (edges : Nat → Edge vertices) (start : Subset V) (n : Nat)
    (x y : Nat → Subset V) (mx : x ∈ histories edges start n)
    (my : y ∈ histories edges start n) (same : ∀ k, k ≤ n → x k = y k) : x = y := by
  funext k
  by_cases earlier : k ≤ n
  · exact same k earlier
  · rw [historyFuture edges start n x mx k (by omega),
      historyFuture edges start n y my k (by omega)]

/-- No legal finite subset history is counted twice, even in collision sectors. -/
theorem historiesNodup (edges : Nat → Edge vertices) (start : Subset V) (n : Nat) :
    (histories edges start n).Nodup := by
  induction n with
  | zero => simp [histories]
  | succ n ih =>
    unfold histories
    apply List.pairwise_flatMap.mpr
    constructor
    · intro old _
      apply (successorsNodup (edges n) (old n)).map (fun next => extendHistory old n next)
      intro x y different equal
      have atNext := congrFun equal (n+1)
      simp only [extendHistory] at atNext
      exact different atNext
    · apply List.Pairwise.imp_of_mem _ ih
      intro x y mx my different first mf second ms equal
      obtain ⟨nx, _, rfl⟩ := List.mem_map.mp mf
      obtain ⟨ny, _, rfl⟩ := List.mem_map.mp ms
      apply different
      apply historyDeterminedByPrefix edges start n x y mx my
      intro k hk
      have atK := congrFun equal k
      simpa [extendHistory, show k ≠ n+1 by omega] using atK

theorem pathMassOne (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) :
    sumWith (fun history => pathProbability edges bs history n) (histories edges start n) = 1 := by
  induction n with
  | zero => simp only [histories, sumWith, pathProbability]; grind
  | succ n ih =>
    simp only [histories, sumWithFlatMap, sumWithMap, extendedPathProbability]
    have inner (history : Nat → Subset V) :
        sumWith (fun next => pathProbability edges bs history n *
          stepProb (edges n) (bs n) (history n) next) (successors (edges n) (history n)) =
          pathProbability edges bs history n := by
      rw [sumWithScale, successorMass]
      simp
    rw [sumWithCongr _ (fun history => pathProbability edges bs history n) _ inner]
    exact ih

theorem pathProbabilityNonnegative (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (history : Nat → Subset V) (n : Nat)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    0 ≤ pathProbability edges bs history n := by
  induction n with
  | zero => change (0 : Rat) ≤ 1; decide
  | succ n ih =>
    exact Rat.mul_nonneg (ih (fun k hk => unitary k (by omega)))
      (stepProbabilityNonnegative (edges n) _ _ _ _ (unitary n (by omega)) _ _)

theorem pathSquaredMagnitudeMassOne (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start : Subset V) (n : Nat)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    sumWith (fun history => normSq (pathAmplitude (vertices := vertices) gaussianArithmetic
      (gateSchedule edges as bs cs ds) history n)) (histories edges start n) = 1 := by
  rw [sumWithCongr _ (fun history => pathProbability edges bs history n) _
    (fun history => pathMagnitude edges as bs cs ds history n unitary)]
  exact pathMassOne edges bs start n

end Paths

#print axioms sumWithAppend
#print axioms sumWithMap
#print axioms sumWithFlatMap
#print axioms sumWithScale
#print axioms sumWithCongr
#print axioms successorMass
#print axioms stepProbabilityNonnegative
#print axioms successorLegal
#print axioms everyLegalSuccessor
#print axioms successorsNodup
#print axioms pathProbabilityPrefixCongr
#print axioms extendedPathProbability
#print axioms historyLegal
#print axioms everyLegalPrefix
#print axioms historyFuture
#print axioms historyDeterminedByPrefix
#print axioms historiesNodup
#print axioms pathMassOne
#print axioms pathProbabilityNonnegative
#print axioms pathSquaredMagnitudeMassOne

end MathCommons.ForestPathNormalization
