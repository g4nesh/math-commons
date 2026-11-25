import ForestGateSupport

/-!
Exact local magnitude results for the concrete signed minor of ForestGateSupport.
All scalars are Gaussian rationals represented by pairs. The statements below
separate arithmetic consequences of unitarity from list-ordering obligations.
No Cauchy–Binet or matrix-library determinant correspondence is assumed.
-/
namespace MathCommons.ForestGateMagnitude
open ForestTrajectory ForestGateSupport

abbrev Gaussian := Rat × Rat

def normSq (z : Gaussian) : Rat := z.1 * z.1 + z.2 * z.2

def conjugate (z : Gaussian) : Gaussian := (z.1, -z.2)

def blockDet (a b c d : Gaussian) : Gaussian :=
  gaussianArithmetic.add (gaussianArithmetic.mul a d)
    (gaussianArithmetic.neg (gaussianArithmetic.mul b c))

/-- The explicit column and row equations of a unitary two-by-two block.
Redundant row equations keep every use of normalization inspectable. -/
structure UnitaryBlock (a b c d : Gaussian) : Prop where
  rowLeft : normSq a + normSq b = 1
  rowRight : normSq c + normSq d = 1
  columnLeft : normSq a + normSq c = 1
  columnRight : normSq b + normSq d = 1
  orthogonal : gaussianArithmetic.add
    (gaussianArithmetic.mul (conjugate a) b)
    (gaussianArithmetic.mul (conjugate c) d) = (0, 0)

theorem normSqZero : normSq (0, 0) = 0 := by simp only [normSq]; grind

theorem normSqOne : normSq (1, 0) = 1 := by simp only [normSq]; grind

theorem normSqNeg (z : Gaussian) :
    normSq (gaussianArithmetic.neg z) = normSq z := by
  simp only [normSq, gaussianArithmetic]
  grind

theorem normSqMul (z w : Gaussian) :
    normSq (gaussianArithmetic.mul z w) = normSq z * normSq w := by
  simp only [normSq, gaussianArithmetic]
  grind

theorem normSqSigned (odd : Bool) (z : Gaussian) :
    normSq (signed gaussianArithmetic odd z) = normSq z := by
  cases odd <;> simp [signed, normSqNeg]

theorem unitaryWeights (a b c d : Gaussian) (h : UnitaryBlock a b c d) :
    normSq c = normSq b ∧ normSq d = normSq a ∧ normSq a = 1 - normSq b := by
  have hl := h.rowLeft
  have hr := h.rowRight
  have hc := h.columnLeft
  grind

theorem gaussianOneMul (z : Gaussian) : gaussianArithmetic.mul gaussianArithmetic.one z = z := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianMulOne (z : Gaussian) : gaussianArithmetic.mul z gaussianArithmetic.one = z := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianAddZero (z : Gaussian) : gaussianArithmetic.add z gaussianArithmetic.zero = z := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianZeroAdd (z : Gaussian) : gaussianArithmetic.add gaussianArithmetic.zero z = z := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem singletonMinor {V : Type} (M : V → V → Gaussian) (r c : V) :
    minorDet gaussianArithmetic M [r] [c] = M r c := by
  simp only [minorDet, columnChoices, scalarSum, signed, List.map_nil, List.map_cons,
    List.isEmpty_nil, Bool.false_eq_true, ↓reduceIte, gaussianMulOne, gaussianAddZero]

theorem gaussianTwoByTwo {V : Type} (M : V → V → Gaussian) (r₁ r₂ c₁ c₂ : V) :
    minorDet gaussianArithmetic M [r₁, r₂] [c₁, c₂] =
      blockDet (M r₁ c₁) (M r₁ c₂) (M r₂ c₁) (M r₂ c₂) := by
  simp only [minorDet, columnChoices, scalarSum, signed, blockDet, List.map_nil,
    List.map_cons, List.isEmpty_nil, Bool.false_eq_true, ↓reduceIte, Bool.not_false,
    gaussianMulOne, gaussianAddZero]

theorem scalarSumMapZero {X : Type} (xs : List X) (f : X → Gaussian)
    (zero : ∀ x ∈ xs, f x = gaussianArithmetic.zero) :
    scalarSum gaussianArithmetic (xs.map f) = gaussianArithmetic.zero := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [List.map_cons, scalarSum]
    rw [zero x List.mem_cons_self, ih (fun y hy => zero y (List.mem_cons_of_mem x hy))]
    exact gaussianAddZero _

/-- Laplace expansion with a single supported entry in its first row. -/
theorem minorFirstColumn {V : Type} (M : V → V → Gaussian)
    (r c : V) (rs cs : List V) (zero : ∀ y ∈ cs, M r y = gaussianArithmetic.zero) :
    minorDet gaussianArithmetic M (r :: rs) (c :: cs) =
      gaussianArithmetic.mul (M r c) (minorDet gaussianArithmetic M rs cs) := by
  simp only [minorDet, columnChoices, List.map_cons, List.map_map, Function.comp_def, scalarSum]
  have tail : scalarSum gaussianArithmetic
      ((columnChoices cs).map (fun q => signed gaussianArithmetic (!q.odd)
        (gaussianArithmetic.mul (M r q.column)
          (minorDet gaussianArithmetic M rs (c :: q.remaining))))) = gaussianArithmetic.zero := by
    apply scalarSumMapZero
    intro q hq
    have hc : q.column ∈ cs :=
      (choicePerm cs q hq).mem_iff.mp List.mem_cons_self
    rw [zero q.column hc, gaussianArithmetic.zeroMul]
    cases q.odd <;> simp only [signed, Bool.not_true, Bool.not_false,
      Bool.false_eq_true, ↓reduceIte, gaussianArithmetic.negZero]
  rw [tail]
  exact gaussianAddZero _

/-- Arbitrarily many common spectator coordinates before the active coordinates
can be removed from the actual signed minor, without any magnitude assumption. -/
theorem identityPrefixMinor {V : Type} [DecidableEq V] (M : V → V → Gaussian)
    (spectators rows cols : List V) (unique : (spectators ++ cols).Nodup)
    (identity : ∀ x ∈ spectators, ∀ y, M x y = if x = y then gaussianArithmetic.one else gaussianArithmetic.zero) :
    minorDet gaussianArithmetic M (spectators ++ rows) (spectators ++ cols) =
      minorDet gaussianArithmetic M rows cols := by
  induction spectators with
  | nil => rfl
  | cons x xs ih =>
    have hx := (List.nodup_cons.mp unique).1
    have ht := (List.nodup_cons.mp unique).2
    rw [List.cons_append, List.cons_append]
    rw [minorFirstColumn M x x (xs ++ rows) (xs ++ cols)]
    · rw [identity x List.mem_cons_self x]
      simp only [↓reduceIte, gaussianOneMul]
      exact ih ht (fun y hy => identity y (List.mem_cons_of_mem x hy))
    · intro y hy
      have hxy : x ≠ y := by intro eq; subst y; exact hx hy
      simp only [identity x List.mem_cons_self y, hxy, ↓reduceIte]

#print axioms normSqZero
#print axioms normSqOne
#print axioms normSqNeg
#print axioms normSqMul
#print axioms normSqSigned
#print axioms unitaryWeights
#print axioms gaussianOneMul
#print axioms gaussianMulOne
#print axioms gaussianAddZero
#print axioms gaussianZeroAdd
#print axioms singletonMinor
#print axioms gaussianTwoByTwo
#print axioms scalarSumMapZero
#print axioms minorFirstColumn
#print axioms identityPrefixMinor

end MathCommons.ForestGateMagnitude
