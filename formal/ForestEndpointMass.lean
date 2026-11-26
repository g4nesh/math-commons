import ForestPathNormalization

/-!
For one fixed forest endpoint, the squared norm of the finite signed amplitude
sum equals the probability mass of its canonical legal histories. The support
and bridge hypotheses force at most one nonzero summand. No Cauchy–Binet or
matrix-product identification is assumed. Empty endpoint fibers are included.
-/
namespace MathCommons.ForestEndpointMass

open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization

theorem sumWithZero {α : Type} (f : α → Rat) (xs : List α)
    (zero : ∀ x ∈ xs, f x = 0) : sumWith f xs = 0 := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [sumWith]
    rw [zero x List.mem_cons_self, ih (fun y hy => zero y (List.mem_cons_of_mem x hy))]
    exact Rat.zero_add 0

/-- The finite no-interference identity when at most one distinct summand is nonzero. -/
theorem normSumUnique {α : Type} (xs : List α) (f : α → Gaussian) (nodup : xs.Nodup)
    (unique : ∀ x ∈ xs, ∀ y ∈ xs, f x ≠ gaussianArithmetic.zero →
      f y ≠ gaussianArithmetic.zero → x = y) :
    normSq (scalarSum gaussianArithmetic (xs.map f)) = sumWith (fun x => normSq (f x)) xs := by
  induction xs with
  | nil => exact normSqZero
  | cons x xs ih =>
    obtain ⟨notMember, tailNodup⟩ := List.nodup_cons.mp nodup
    by_cases headZero : f x = gaussianArithmetic.zero
    · simp only [List.map_cons, scalarSum, sumWith]
      have normZero : normSq gaussianArithmetic.zero = 0 := normSqZero
      rw [headZero, gaussianZeroAdd, normZero, Rat.zero_add]
      apply ih tailNodup
      intro a ha b hb na nb
      exact unique a (List.mem_cons_of_mem x ha) b (List.mem_cons_of_mem x hb) na nb
    · have tailZero : ∀ y ∈ xs, f y = gaussianArithmetic.zero := by
        intro y hy
        by_cases zero : f y = gaussianArithmetic.zero
        · exact zero
        · have equal := unique x List.mem_cons_self y (List.mem_cons_of_mem x hy) headZero zero
          exfalso
          apply notMember
          rw [equal]
          exact hy
      have amplitudeTail := scalarSumMapZero xs f tailZero
      have normTail : sumWith (fun y => normSq (f y)) xs = 0 := by
        apply sumWithZero
        intro y hy
        rw [tailZero y hy]
        exact normSqZero
      simp only [List.map_cons, scalarSum, sumWith]
      rw [amplitudeTail, gaussianAddZero, normTail, Rat.add_zero]

section Endpoint
variable {V : Type} [DecidableEq V] {vertices : List V}

/-- A computable endpoint test on the supplied complete finite vertex enumeration. -/
def endpointHistories (edges : Nat → Edge vertices) (start finish : Subset V) (n : Nat) :
    List (Nat → Subset V) :=
  (histories edges start n).filter (fun history => vertices.all (fun x => history n x == finish x))

theorem endpointMember (edges : Nat → Edge vertices) (start finish : Subset V) (n : Nat)
    (history : Nat → Subset V) :
    history ∈ endpointHistories edges start finish n ↔
      history ∈ histories edges start n ∧ history n = finish := by
  simp only [endpointHistories, List.mem_filter]
  constructor
  · intro ⟨member, all⟩
    refine ⟨member, ?_⟩
    funext x
    have same := List.all_eq_true.mp all x ((edges 0).complete x)
    simpa only [beq_iff_eq] using same
  · intro ⟨member, same⟩
    refine ⟨member, ?_⟩
    apply List.all_eq_true.mpr
    intro x _
    simp [same]

theorem endpointHistoriesNodup (edges : Nat → Edge vertices) (start finish : Subset V) (n : Nat) :
    (endpointHistories edges start finish n).Nodup := by
  exact (historiesNodup edges start n).filter _

/-- Canonical finite histories with the same endpoints cannot both carry nonzero
amplitude in a forest schedule. No unitarity is needed for this support assertion. -/
theorem endpointNonzeroUnique (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n) (x y : Nat → Subset V)
    (mx : x ∈ endpointHistories edges start finish n)
    (my : y ∈ endpointHistories edges start finish n)
    (nx : pathAmplitude (vertices := vertices) gaussianArithmetic
      (gateSchedule edges as bs cs ds) x n ≠ gaussianArithmetic.zero)
    (ny : pathAmplitude (vertices := vertices) gaussianArithmetic
      (gateSchedule edges as bs cs ds) y n ≠ gaussianArithmetic.zero) : x = y := by
  obtain ⟨hx, fx⟩ := (endpointMember edges start finish n x).mp mx
  obtain ⟨hy, fy⟩ := (endpointMember edges start finish n y).mp my
  have sx := (historyLegal edges start n x hx).1
  have sy := (historyLegal edges start n y hy).1
  apply historyDeterminedByPrefix edges start n x y hx hy
  exact nonzeroAmplitudesUnique gaussianArithmetic edges (gateSchedule edges as bs cs ds) n
    bridges (fun k _ => twoCoordinateGateSupport gaussianArithmetic (edges k) _ _ _ _)
    x y nx ny (sx.trans sy.symm) (fx.trans fy.symm)

def endpointAmplitude (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian)
    (start finish : Subset V) (n : Nat) : Gaussian :=
  scalarSum gaussianArithmetic ((endpointHistories edges start finish n).map fun history =>
    pathAmplitude (vertices := vertices) gaussianArithmetic
      (gateSchedule edges as bs cs ds) history n)

def endpointProbability (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start finish : Subset V) (n : Nat) : Rat :=
  sumWith (fun history => pathProbability edges bs history n)
    (endpointHistories edges start finish n)

/-- Main forest endpoint identity. The left side is an explicit finite path sum,
not an asserted entry of a compound matrix of the gate product. -/
theorem endpointMagnitude (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian)
    (start finish : Subset V) (n : Nat) (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    normSq (endpointAmplitude edges as bs cs ds start finish n) =
      endpointProbability edges bs start finish n := by
  unfold endpointAmplitude endpointProbability
  rw [normSumUnique _ _ (endpointHistoriesNodup edges start finish n)
    (fun x mx y my nx ny => endpointNonzeroUnique edges as bs cs ds start finish n bridges
      x y mx my nx ny)]
  exact sumWithCongr _ _ _ (fun history => pathMagnitude edges as bs cs ds history n unitary)

theorem emptyEndpointFiber (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian)
    (start finish : Subset V) (n : Nat)
    (empty : endpointHistories edges start finish n = []) :
    endpointAmplitude edges as bs cs ds start finish n = gaussianArithmetic.zero ∧
      endpointProbability edges bs start finish n = 0 := by
  simp [endpointAmplitude, endpointProbability, empty, scalarSum, sumWith]

end Endpoint

#print axioms sumWithZero
#print axioms normSumUnique
#print axioms endpointMember
#print axioms endpointHistoriesNodup
#print axioms endpointNonzeroUnique
#print axioms endpointMagnitude
#print axioms emptyEndpointFiber

end MathCommons.ForestEndpointMass
