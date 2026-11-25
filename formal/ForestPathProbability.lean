import ForestGateMagnitude

/-!
# Squared probabilities along a forest path

Pass 1: the local probability table `localProb` and the single-step lemma
`stepMagnitudeLeftFirst` for an ambient vertex order with `e.left` before `e.right`.
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

end Step

#print axioms localProbIllegal
#print axioms nonzeroEndpointsLegal
#print axioms stepMagnitudeLeftFirst

end MathCommons.ForestPathProbability
