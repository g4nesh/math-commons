import ForestSubsetExchange

/-! A common finite permutation law for every actual forest gate-product minor.
Permutations carry explicit inverse proofs; weights do not depend on a subset.
Gaussian-rational actual minors and finite `Edge` enumerations, Std Lean 4.24.0. -/
namespace MathCommons.ForestPermutationLaw
open Classical
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory ForestSubsetExchange

/-- A genuine bijection with a specified two-sided inverse. -/
structure VertexPermutation (V : Type) where
  forward : V → V
  inverse : V → V
  inverse_forward : ∀ x, inverse (forward x) = x
  forward_inverse : ∀ x, forward (inverse x) = x

def identityPermutation {V : Type} : VertexPermutation V :=
  ⟨id, id, fun _ => rfl, fun _ => rfl⟩

def subsetImage {V : Type} (p : VertexPermutation V) (s : Subset V) : Subset V :=
  fun x => s (p.inverse x)

theorem permutationInjective {V : Type} (p : VertexPermutation V) (x y : V)
    (same : p.forward x = p.forward y) : x = y := by
  have h := congrArg p.inverse same
  simpa only [p.inverse_forward] using h

theorem permutationSurjective {V : Type} (p : VertexPermutation V) (y : V) :
    ∃ x, p.forward x = y := ⟨p.inverse y, p.forward_inverse y⟩

theorem subsetImageForward {V : Type} (p : VertexPermutation V) (s : Subset V) (x : V) :
    subsetImage p s (p.forward x) = s x := by
  simp only [subsetImage, p.inverse_forward]

section Law
variable {V : Type} [DecidableEq V] {vertices : List V}

def vertexSwap (e : Edge vertices) (x : V) : V :=
  if x = e.left then e.right else if x = e.right then e.left else x

theorem vertexSwapInvolution (e : Edge vertices) (x : V) :
    vertexSwap e (vertexSwap e x) = x := by
  by_cases left : x = e.left
  · subst x; simp [vertexSwap, Ne.symm e.distinct]
  · by_cases right : x = e.right
    · subst x; simp [vertexSwap, Ne.symm e.distinct]
    · simp [vertexSwap, left, right]

def swapAfter (e : Edge vertices) (p : VertexPermutation V) : VertexPermutation V where
  forward := fun x => vertexSwap e (p.forward x)
  inverse := fun x => p.inverse (vertexSwap e x)
  inverse_forward := by intro x; rw [vertexSwapInvolution, p.inverse_forward]
  forward_inverse := by intro x; rw [p.forward_inverse, vertexSwapInvolution]

omit [DecidableEq V] in
theorem subsetImageIdentity (s : Subset V) : subsetImage identityPermutation s = s := rfl

theorem subsetImageSwap (e : Edge vertices) (p : VertexPermutation V) (s : Subset V) :
    subsetImage (swapAfter e p) s = swapSubset e (subsetImage p s) := by
  funext x
  simp only [subsetImage, swapAfter, vertexSwap, swapSubset]
  split
  · rfl
  · split <;> rfl

/-- The complete independent stay/swap coin law. Repeated permutation images
are deliberately retained as separate nonnegative contributions. -/
def coinLaw (edges : Nat → Edge vertices) (bs : Nat → Gaussian) :
    Nat → List (VertexPermutation V × Rat)
  | 0 => [(identityPermutation, 1)]
  | n+1 => (coinLaw edges bs n).flatMap fun entry =>
      [(entry.1, entry.2 * (1 - normSq (bs n))),
       (swapAfter (edges n) entry.1, entry.2 * normSq (bs n))]

def lawExpectation (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) (f : Subset V → Rat) : Rat :=
  sumWith (fun entry => entry.2 * f (subsetImage entry.1 start)) (coinLaw edges bs n)

def historyExpectation (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) (f : Subset V → Rat) : Rat :=
  sumWith (fun h => pathProbability edges bs h n * f (h n)) (histories edges start n)

/-- Independent coin choices aggregate to the canonical successor rule, even
when both choices give the same occupied subset. No bound on the weight needed. -/
theorem successorExpectation (e : Edge vertices) (b : Gaussian)
    (s : Subset V) (f : Subset V → Rat) :
    sumWith (fun next => stepProb e b s next * f next) (successors e s) =
      (1 - normSq b) * f s + normSq b * f (swapSubset e s) := by
  have stay := spectatorsAgreeOfStep e s s (Or.inl rfl)
  have swap := spectatorsAgreeOfStep e s (swapSubset e s) (Or.inr rfl)
  cases hl : s e.left <;> cases hr : s e.right
  · have same := swapCollision e s (hl.trans hr.symm)
    simp [successors, sumWith, stepProb, stay, localProb, hl, hr, same]
    grind
  · simp [successors, sumWith, stepProb, stay, swap, localProb, hl, hr, swapAtLeft, swapAtRight, Rat.add_zero]
  · simp [successors, sumWith, stepProb, stay, swap, localProb, hl, hr, swapAtLeft, swapAtRight, Rat.add_zero]
  · have same := swapCollision e s (hl.trans hr.symm)
    simp [successors, sumWith, stepProb, stay, localProb, hl, hr, same]
    grind

theorem lawExpectationStep (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) (f : Subset V → Rat) :
    lawExpectation edges bs start (n+1) f =
      lawExpectation edges bs start n (fun s =>
        (1 - normSq (bs n)) * f s + normSq (bs n) * f (swapSubset (edges n) s)) := by
  simp only [lawExpectation, coinLaw, sumWithFlatMap, sumWith]
  apply sumWithCongr
  intro entry
  rw [subsetImageSwap]
  grind

theorem historyExpectationStep (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) (f : Subset V → Rat) :
    historyExpectation edges bs start (n+1) f =
      historyExpectation edges bs start n (fun s =>
        (1 - normSq (bs n)) * f s + normSq (bs n) * f (swapSubset (edges n) s)) := by
  simp only [historyExpectation, histories, sumWithFlatMap, sumWithMap, extendedPathProbability]
  apply sumWithCongr
  intro h
  simp only [extendHistory, ↓reduceIte]
  have factor : sumWith (fun next =>
      pathProbability edges bs h n * stepProb (edges n) (bs n) (h n) next * f next)
      (successors (edges n) (h n)) =
      pathProbability edges bs h n * sumWith
        (fun next => stepProb (edges n) (bs n) (h n) next * f next)
        (successors (edges n) (h n)) := by
    rw [← sumWithScale]
    apply sumWithCongr
    intro next
    grind
  rw [factor, successorExpectation]

/-- The same permutation law pushes forward correctly for every starting
subset and every rational-valued test function; no forest or unitarity needed. -/
theorem commonPushforward (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) (f : Subset V → Rat) :
    lawExpectation edges bs start n f = historyExpectation edges bs start n f := by
  induction n generalizing f with
  | zero => simp [lawExpectation, historyExpectation, coinLaw, histories, sumWith,
      pathProbability, subsetImageIdentity]
  | succ n ih => rw [lawExpectationStep, historyExpectationStep, ih]


/-- Algebraic normalization of the single common law; no subset is selected. -/
theorem coinLawMassOne (edges : Nat → Edge vertices) (bs : Nat → Gaussian) (n : Nat) :
    sumWith (fun entry => entry.2) (coinLaw edges bs n) = 1 := by
  have push := commonPushforward edges bs (fun _ => false) n (fun _ => 1)
  simp only [lawExpectation, historyExpectation, Rat.mul_one] at push
  rw [push]
  exact pathMassOne edges bs (fun _ => false) n

theorem coinWeightsNonnegative (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (n : Nat) (bounds : ∀ k, k < n → 0 ≤ normSq (bs k) ∧ normSq (bs k) ≤ 1)
    (entry : VertexPermutation V × Rat) (member : entry ∈ coinLaw edges bs n) :
    0 ≤ entry.2 := by
  induction n generalizing entry with
  | zero =>
    simp only [coinLaw, List.mem_cons, List.not_mem_nil, or_false] at member
    subst entry
    change (0 : Rat) ≤ 1
    decide
  | succ n ih =>
    obtain ⟨old, oldMember, next⟩ := List.mem_flatMap.mp member
    have oldNonnegative := ih (fun k hk => bounds k (by omega)) old oldMember
    have current := bounds n (by omega)
    simp only [List.mem_cons, List.not_mem_nil, or_false] at next
    rcases next with rfl | rfl
    · exact Rat.mul_nonneg oldNonnegative (by grind)
    · exact Rat.mul_nonneg oldNonnegative current.1

theorem unitaryCoinWeightsNonnegative (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k))
    (entry : VertexPermutation V × Rat) (member : entry ∈ coinLaw edges bs n) :
    0 ≤ entry.2 := by
  exact coinWeightsNonnegative edges bs n
    (fun k hk => unitaryProbabilityBounds _ _ _ _ (unitary k hk)) entry member

/-- The push-forward endpoint mass of the common permutation law equals the
canonical history probability, including collision aggregation. -/
theorem commonEndpointLaw (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start finish : Subset V) (n : Nat) :
    sumWith (fun entry => if subsetImage entry.1 start = finish then entry.2 else 0)
      (coinLaw edges bs n) = endpointProbability edges bs start finish n := by
  have push := commonPushforward edges bs start n (fun s => if s = finish then 1 else 0)
  have term (weight : Rat) (s : Subset V) :
      weight * (if s = finish then 1 else 0) = if s = finish then weight else 0 := by
    split <;> simp only [Rat.mul_one, Rat.mul_zero]
  simp only [lawExpectation, historyExpectation, term] at push
  rw [← endpointProbabilityAsSum] at push
  exact push

/-- One finite nonnegative normalized law of genuine permutations simultaneously
represents every squared actual forest gate-product minor. The law is outside
(and independent of) the universal quantifiers over both occupied subsets. -/
theorem commonSquaredMinorLaw (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    (∀ entry ∈ coinLaw edges bs n, 0 ≤ entry.2) ∧
    sumWith (fun entry => entry.2) (coinLaw edges bs n) = 1 ∧
    ∀ start finish : Subset V,
      normSq (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
        (vertices.filter finish) (vertices.filter start)) =
      sumWith (fun entry => if subsetImage entry.1 start = finish then entry.2 else 0)
        (coinLaw edges bs n) := by
  refine ⟨?_, coinLawMassOne edges bs n, ?_⟩
  · intro entry member
    exact unitaryCoinWeightsNonnegative edges as bs cs ds n unitary entry member
  · intro start finish
    rw [productMagnitude edges as bs cs ds start finish n bridges unitary,
      commonEndpointLaw]

end Law
#print axioms permutationInjective
#print axioms permutationSurjective
#print axioms subsetImageForward
#print axioms vertexSwapInvolution
#print axioms subsetImageIdentity
#print axioms subsetImageSwap
#print axioms successorExpectation
#print axioms lawExpectationStep
#print axioms historyExpectationStep
#print axioms commonPushforward
#print axioms coinLawMassOne
#print axioms coinWeightsNonnegative
#print axioms unitaryCoinWeightsNonnegative
#print axioms commonEndpointLaw
#print axioms commonSquaredMinorLaw

end MathCommons.ForestPermutationLaw
