import ForestGateMagnitude

/-!
# Squared probabilities along a forest path

`localProb` is the local probability table. `stepMagnitude` proves that, for a unitary
two-coordinate gate and subsets agreeing off the gate's endpoints, the squared minor
between `vertices.filter t` and `vertices.filter s` equals `localProb`, in either
ambient order of the endpoints (`stepMagnitudeLeftFirst`, `stepMagnitudeRightFirst`).
`stepMagnitudeGeneral` drops the agreement hypothesis: `stepProb` is 0 when the subsets
disagree off the endpoints, matching the vanishing minor (`gaussianGateMinorStep`).

`pathMagnitude` is the path statement: for any subset history along a schedule of unitary
gates, `normSq (pathAmplitude ...) = pathProbability ...`. Zero factors, the empty schedule
(`emptySchedule`, `zeroStepPath`) and collision sectors (`collisionStay`) are covered.
`nonzeroProbabilityUnique` combines it with `ForestGateSupport.nonzeroAmplitudesUnique`.

Declarations used: `ForestGateMagnitude.orderedGateMagnitudes`, `filterTwoEndpoints`,
`decomposeTwoEndpoints`, `reverseEdge`, `reverseGate`, `reverseUnitary`, `unitaryWeights`,
`normSqMul`, `normSqOne`, `normSqZero`; `ForestGateSupport.gaussianGateMinorStep`,
`twoCoordinateGateSupport`, `pathAmplitude`, `nonzeroAmplitudesUnique`;
`ForestTrajectory.swapAtLeft`, `swapAtRight`.

Cauchy–Binet is left explicit: identifying the sum over histories of `pathAmplitude` with
a compound-matrix entry of the gate product is NOT proved or assumed here.
All scalars are Gaussian rationals. No Cauchy–Binet identity is used or claimed.
-/
namespace MathCommons.ForestPathProbability

open ForestTrajectory ForestGateSupport ForestGateMagnitude

section Step
variable {V : Type} [DecidableEq V] {vertices : List V}

/-- Local squared-probability table keyed by the endpoint occupancy of the source
subset (`sl`, `sr`) and of the target subset (`tl`, `tr`), with `q = |b|^2`.
Every combination that is not a stay or an exchange has probability 0. -/
def localProb (q : Rat) : Bool → Bool → Bool → Bool → Rat
  | false, false, false, false => 1
  | true, true, true, true => 1
  | true, false, true, false => 1 - q
  | false, true, false, true => 1 - q
  | true, false, false, true => q
  | false, true, true, false => q
  | _, _, _, _ => 0

/-- Outside a stay or an exchange of the two endpoints the local probability is 0. -/
theorem localProbIllegal (q : Rat) (sl sr tl tr : Bool)
    (illegal : ¬((tl = sl ∧ tr = sr) ∨ (tl = sr ∧ tr = sl))) :
    localProb q sl sr tl tr = 0 := by
  cases sl <;> cases sr <;> cases tl <;> cases tr <;> simp_all [localProb]

/-- A nonzero minor between two subsets forces a legal stay-or-swap on the endpoints. -/
theorem nonzeroEndpointsLegal (e : Edge vertices) (a b c d : Gaussian) (s t : Subset V)
    (nonzero : minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter t) (vertices.filter s) ≠ (0, 0)) :
    (t e.left = s e.left ∧ t e.right = s e.right) ∨
      (t e.left = s e.right ∧ t e.right = s e.left) := by
  rcases gaussianGateMinorStep e a b c d s t nonzero with rfl | rfl
  · exact Or.inl ⟨rfl, rfl⟩
  · exact Or.inr ⟨swapAtLeft e s, swapAtRight e s⟩

theorem stepMagnitudeLeftFirst (e : Edge vertices) (a b c d : Gaussian)
    (h : UnitaryBlock a b c d) (before middle after : List V)
    (order : vertices = before ++ e.left :: (middle ++ e.right :: after))
    (s t : Subset V) (spectators : ∀ x, x ≠ e.left → x ≠ e.right → t x = s x) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter t) (vertices.filter s)) =
      localProb (normSq b) (s e.left) (s e.right) (t e.left) (t e.right) := by
  have nd : (before ++ e.left :: (middle ++ e.right :: after)).Nodup := order ▸ e.nodup
  have avoidsAll : ∀ x ∈ before ++ (middle ++ after), x ≠ e.left ∧ x ≠ e.right := by
    simp only [List.nodup_append, List.nodup_cons, List.mem_append, List.mem_cons] at nd
    intro x hx
    simp only [List.mem_append] at hx
    grind
  have same : ∀ xs : List V, (∀ x ∈ xs, x ≠ e.left ∧ x ≠ e.right) →
      xs.filter t = xs.filter s := by
    intro xs hxs
    exact List.filter_congr (fun x hx => spectators x (hxs x hx).1 (hxs x hx).2)
  have hb := same before (fun x hx => avoidsAll x (List.mem_append_left _ hx))
  have hm := same middle
    (fun x hx => avoidsAll x (List.mem_append_right _ (List.mem_append_left _ hx)))
  have ha := same after
    (fun x hx => avoidsAll x (List.mem_append_right _ (List.mem_append_right _ hx)))
  have fs := (congrArg (List.filter s) order).trans
    (filterTwoEndpoints s e.left e.right before middle after)
  have ft := (congrArg (List.filter t) order).trans
    (filterTwoEndpoints t e.left e.right before middle after)
  have sub : (before.filter s ++ e.left :: (middle.filter s ++ e.right :: after.filter s)).Sublist
      (before ++ e.left :: (middle ++ e.right :: after)) :=
    (List.filter_sublist).append (((List.filter_sublist).append
      ((List.filter_sublist).cons₂ e.right)).cons₂ e.left)
  have hu := sub.nodup nd
  have hav : ∀ x ∈ before.filter s ++ (middle.filter s ++ after.filter s),
      x ≠ e.left ∧ x ≠ e.right := fun x hx =>
    avoidsAll x (((List.filter_sublist).append
      ((List.filter_sublist).append (List.filter_sublist))).subset hx)
  obtain ⟨t0, tll, trr, tlr, trl, tb⟩ :=
    orderedGateMagnitudes e a b c d h (before.filter s) (middle.filter s) (after.filter s) hu hav
  have zero : ¬((t e.left = s e.left ∧ t e.right = s e.right) ∨
      (t e.left = s e.right ∧ t e.right = s e.left)) →
      normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
        (vertices.filter t) (vertices.filter s)) = 0 := by
    intro illegal
    refine Classical.byContradiction fun nonzeroNorm => illegal ?_
    apply nonzeroEndpointsLegal e a b c d s t
    intro z
    exact nonzeroNorm (by rw [z]; exact normSqZero)
  by_cases leg : (t e.left = s e.left ∧ t e.right = s e.right) ∨
      (t e.left = s e.right ∧ t e.right = s e.left)
  · rw [ft, fs, hb, hm, ha]
    cases hsl : s e.left <;> cases hsr : s e.right <;> cases htl : t e.left <;>
      cases htr : t e.right <;>
      simp only [List.append_assoc, List.cons_append, List.nil_append,
        List.append_nil, Bool.false_eq_true, ↓reduceIte, localProb]
    all_goals first
      | exact t0 | exact tll | exact trr | exact tlr | exact trl | exact tb
      | (exfalso; simp [hsl, hsr, htl, htr] at leg)
  · rw [zero leg]
    exact (localProbIllegal _ _ _ _ _ leg).symm

/-- Exchanging the roles of the two endpoints leaves the table unchanged. -/
theorem localProbSwap (q : Rat) (sl sr tl tr : Bool) :
    localProb q sr sl tr tl = localProb q sl sr tl tr := by
  cases sl <;> cases sr <;> cases tl <;> cases tr <;> rfl

/-- The opposite ambient order, reduced to `stepMagnitudeLeftFirst` through `reverseEdge`. -/
theorem stepMagnitudeRightFirst (e : Edge vertices) (a b c d : Gaussian)
    (h : UnitaryBlock a b c d) (before middle after : List V)
    (order : vertices = before ++ e.right :: (middle ++ e.left :: after))
    (s t : Subset V) (spectators : ∀ x, x ≠ e.left → x ≠ e.right → t x = s x) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter t) (vertices.filter s)) =
      localProb (normSq b) (s e.left) (s e.right) (t e.left) (t e.right) := by
  have reversed := stepMagnitudeLeftFirst (reverseEdge e) d c b a (reverseUnitary a b c d h)
    before middle after order s t (fun x xl xr => spectators x xr xl)
  rw [reverseGate] at reversed
  rw [reversed, (unitaryWeights a b c d h).1]
  exact localProbSwap _ _ _ _ _

/-- The single-step identity in either ambient order. -/
theorem stepMagnitude (e : Edge vertices) (a b c d : Gaussian) (h : UnitaryBlock a b c d)
    (s t : Subset V) (spectators : ∀ x, x ≠ e.left → x ≠ e.right → t x = s x) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter t) (vertices.filter s)) =
      localProb (normSq b) (s e.left) (s e.right) (t e.left) (t e.right) := by
  rcases decomposeTwoEndpoints vertices e.left e.right e.distinct (e.complete _) (e.complete _)
    with ⟨before, middle, after, order⟩ | ⟨before, middle, after, order⟩
  · exact stepMagnitudeLeftFirst e a b c d h before middle after order s t spectators
  · exact stepMagnitudeRightFirst e a b c d h before middle after order s t spectators

/-- Decidable test that two subsets agree at every vertex other than the gate's endpoints. -/
def spectatorsAgree (e : Edge vertices) (s t : Subset V) : Bool :=
  vertices.all (fun x => decide (x = e.left) || decide (x = e.right) || t x == s x)

theorem spectatorsAgreeSpec (e : Edge vertices) (s t : Subset V)
    (agree : spectatorsAgree e s t = true) :
    ∀ x, x ≠ e.left → x ≠ e.right → t x = s x := by
  intro x xl xr
  have hx := List.all_eq_true.mp agree x (e.complete x)
  simp only [Bool.or_eq_true, decide_eq_true_eq, beq_iff_eq] at hx
  rcases hx with (l | r) | same
  · exact absurd l xl
  · exact absurd r xr
  · exact same

theorem spectatorsAgreeOfStep (e : Edge vertices) (s t : Subset V) (step : Step e s t) :
    spectatorsAgree e s t = true := by
  apply List.all_eq_true.mpr
  intro x _
  by_cases xl : x = e.left
  · simp [xl]
  · by_cases xr : x = e.right
    · simp [xr]
    · rcases step with rfl | rfl
      · simp
      · simp [swapSubset, xl, xr]

/-- Probability of one gate step: the local table when the spectators agree, else 0. -/
def stepProb (e : Edge vertices) (b : Gaussian) (s t : Subset V) : Rat :=
  if spectatorsAgree e s t then
    localProb (normSq b) (s e.left) (s e.right) (t e.left) (t e.right)
  else 0

/-- The single-step identity for arbitrary subsets, including spectator disagreement. -/
theorem stepMagnitudeGeneral (e : Edge vertices) (a b c d : Gaussian)
    (h : UnitaryBlock a b c d) (s t : Subset V) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter t) (vertices.filter s)) = stepProb e b s t := by
  unfold stepProb
  cases agree : spectatorsAgree e s t
  · have zero : minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
        (vertices.filter t) (vertices.filter s) = (0, 0) :=
      Classical.byContradiction fun nonzero => by
        have := spectatorsAgreeOfStep e s t (gaussianGateMinorStep e a b c d s t nonzero)
        rw [agree] at this
        exact Bool.false_ne_true this
    rw [zero, normSqZero]
    rfl
  · exact stepMagnitude e a b c d h s t (spectatorsAgreeSpec e s t agree)

/-- Ordered product of step probabilities along a subset history; 1 on the empty schedule. -/
def pathProbability (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (history : Nat → Subset V) : Nat → Rat
  | 0 => 1
  | n + 1 => pathProbability edges bs history n * stepProb (edges n) (bs n) (history n) (history (n + 1))

/-- The gate matrices of a schedule of two-coordinate gates. -/
def gateSchedule (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian) :
    Nat → V → V → Gaussian :=
  fun k => twoCoordinateGate gaussianArithmetic (edges k) (as k) (bs k) (cs k) (ds k)

/-- Main result: along ANY subset history of a schedule of unitary two-coordinate gates,
the squared modulus of the ordered product of gate minors equals the product of the
local probabilities. Illegal steps and spectator disagreement contribute 0; the empty
schedule gives 1. No Cauchy–Binet identity is used: identifying a sum of such products
with a compound-matrix entry remains a separate, explicit step. -/
theorem pathMagnitude (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian)
    (history : Nat → Subset V) (n : Nat)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    normSq (pathAmplitude (vertices := vertices) gaussianArithmetic
      (gateSchedule edges as bs cs ds) history n) = pathProbability edges bs history n := by
  induction n with
  | zero => exact normSqOne
  | succ n ih =>
    simp only [pathAmplitude, pathProbability]
    rw [normSqMul, ih (fun k hk => unitary k (by omega))]
    rw [gateSchedule, stepMagnitudeGeneral (edges n) (as n) (bs n) (cs n) (ds n)
      (unitary n (by omega))]

/-- The empty schedule has probability 1. -/
theorem emptySchedule (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (history : Nat → Subset V) : pathProbability edges bs history 0 = 1 := rfl

/-- A single zero step makes the whole path probability 0. -/
theorem zeroStepPath (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (history : Nat → Subset V) (n k : Nat) (hk : k < n)
    (zero : stepProb (edges k) (bs k) (history k) (history (k + 1)) = 0) :
    pathProbability edges bs history n = 0 := by
  induction n with
  | zero => omega
  | succ n ih =>
    simp only [pathProbability]
    by_cases last : k = n
    · subst last; rw [zero, Rat.mul_zero]
    · rw [ih (by omega), Rat.zero_mul]

/-- In a collision sector (both endpoints occupied or both empty), staying has probability 1. -/
theorem collisionStay (e : Edge vertices) (b : Gaussian) (s : Subset V)
    (collision : s e.left = s e.right) : stepProb e b s s = 1 := by
  have agree : spectatorsAgree e s s = true := spectatorsAgreeOfStep e s s (Or.inl rfl)
  unfold stepProb
  rw [agree]
  simp only [↓reduceIte]
  rw [← collision]
  cases s e.left <;> rfl

/-- With `ForestGateSupport.nonzeroAmplitudesUnique`: along a forest schedule, at most one
history between fixed endpoint subsets has nonzero path probability. -/
theorem nonzeroProbabilityUnique (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian)
    (n : Nat) (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k))
    (x y : Nat → Subset V)
    (nonzeroX : pathProbability edges bs x n ≠ 0) (nonzeroY : pathProbability edges bs y n ≠ 0)
    (start : x 0 = y 0) (finish : x n = y n) : ∀ k, k ≤ n → x k = y k := by
  have amplitude : ∀ z : Nat → Subset V, pathProbability edges bs z n ≠ 0 →
      pathAmplitude (vertices := vertices) gaussianArithmetic
        (gateSchedule edges as bs cs ds) z n ≠ gaussianArithmetic.zero := by
    intro z nz zero
    apply nz
    rw [← pathMagnitude edges as bs cs ds z n unitary, zero]
    exact normSqZero
  exact nonzeroAmplitudesUnique gaussianArithmetic edges (gateSchedule edges as bs cs ds) n
    bridges (fun k _ => twoCoordinateGateSupport gaussianArithmetic (edges k) _ _ _ _)
    x y (amplitude x nonzeroX) (amplitude y nonzeroY) start finish

end Step

#print axioms localProbIllegal
#print axioms nonzeroEndpointsLegal
#print axioms stepMagnitudeLeftFirst
#print axioms localProbSwap
#print axioms stepMagnitudeRightFirst
#print axioms stepMagnitude
#print axioms spectatorsAgreeSpec
#print axioms spectatorsAgreeOfStep
#print axioms stepMagnitudeGeneral
#print axioms pathMagnitude
#print axioms emptySchedule
#print axioms zeroStepPath
#print axioms collisionStay
#print axioms nonzeroProbabilityUnique

end MathCommons.ForestPathProbability
