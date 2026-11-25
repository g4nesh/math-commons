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


def parity : List V → Bool
  | [] => false
  | _ :: xs => !(parity xs)

theorem lastChoiceMember {V : Type} (xs : List V) (r : V) :
    (⟨r, xs, parity xs⟩ : ColumnChoice V) ∈ columnChoices (xs ++ [r]) := by
  induction xs with
  | nil => simp [columnChoices, parity]
  | cons x xs ih =>
    simp only [List.cons_append, columnChoices, List.mem_cons, parity]
    right
    exact List.mem_map.mpr ⟨⟨r, xs, parity xs⟩, ih, rfl⟩

/-- A sum whose keys are distinct and whose only possible nonzero term has key q. -/
theorem scalarSumUniqueKey {X K : Type} (xs : List X) (key : X → K)
    (unique : (xs.map key).Nodup) (q : X) (member : q ∈ xs) (f : X → Gaussian)
    (zero : ∀ x ∈ xs, key x ≠ key q → f x = gaussianArithmetic.zero) :
    scalarSum gaussianArithmetic (xs.map f) = f q := by
  induction xs with
  | nil => simp at member
  | cons x xs ih =>
    have hu := List.nodup_cons.mp unique
    rcases List.mem_cons.mp member with rfl | member
    · simp only [List.map_cons, scalarSum]
      have tail : scalarSum gaussianArithmetic (xs.map f) = gaussianArithmetic.zero := by
        apply scalarSumMapZero
        intro y hy
        apply zero y (List.mem_cons_of_mem q hy)
        intro eq
        exact hu.1 (List.mem_map.mpr ⟨y, hy, eq⟩)
      rw [tail, gaussianAddZero]
    · have different : key x ≠ key q := by
        intro eq
        exact hu.1 (List.mem_map.mpr ⟨q, member, eq.symm⟩)
      simp only [List.map_cons, scalarSum]
      rw [zero x List.mem_cons_self different, gaussianZeroAdd]
      exact ih hu.2 member (fun y hy => zero y (List.mem_cons_of_mem x hy))

/-- Exact ordering sign of a single supported entry in the last column. -/
theorem minorLastColumn {V : Type} (M : V → V → Gaussian)
    (r c : V) (rows before : List V) (unique : (before ++ [c]).Nodup)
    (zero : ∀ y ∈ before, M r y = gaussianArithmetic.zero) :
    minorDet gaussianArithmetic M (r :: rows) (before ++ [c]) =
      signed gaussianArithmetic (parity before)
        (gaussianArithmetic.mul (M r c) (minorDet gaussianArithmetic M rows before)) := by
  apply scalarSumUniqueKey (columnChoices (before ++ [c])) ColumnChoice.column
    (by simpa only [columnChoicesColumns] using unique)
    ⟨c, before, parity before⟩ (lastChoiceMember before c)
  intro q hq different
  have member : q.column ∈ before ++ [c] :=
    (choicePerm (before ++ [c]) q hq).mem_iff.mp List.mem_cons_self
  have inside : q.column ∈ before := by
    rcases List.mem_append.mp member with h | h
    · exact h
    · have eq : q.column = c := by simpa using h
      exact False.elim (different eq)
  rw [zero q.column inside, gaussianArithmetic.zeroMul]
  cases q.odd <;> simp only [signed, Bool.false_eq_true, ↓reduceIte,
    gaussianArithmetic.negZero]


theorem selectedChoiceMember {V : Type} (before after : List V) (c : V) :
    (⟨c, before ++ after, parity before⟩ : ColumnChoice V) ∈
      columnChoices (before ++ c :: after) := by
  induction before with
  | nil => simp [columnChoices, parity]
  | cons x xs ih =>
    simp only [List.cons_append, columnChoices, List.mem_cons, parity]
    right
    exact List.mem_map.mpr ⟨⟨c, xs ++ after, parity xs⟩, ih, rfl⟩

/-- The signed Laplace value for an arbitrary supported column position. -/
theorem minorSelectedColumn {V : Type} (M : V → V → Gaussian)
    (r c : V) (rows before after : List V) (unique : (before ++ c :: after).Nodup)
    (zero : ∀ y ∈ before ++ after, M r y = gaussianArithmetic.zero) :
    minorDet gaussianArithmetic M (r :: rows) (before ++ c :: after) =
      signed gaussianArithmetic (parity before)
        (gaussianArithmetic.mul (M r c)
          (minorDet gaussianArithmetic M rows (before ++ after))) := by
  apply scalarSumUniqueKey (columnChoices (before ++ c :: after)) ColumnChoice.column
    (by simpa only [columnChoicesColumns] using unique)
    ⟨c, before ++ after, parity before⟩ (selectedChoiceMember before after c)
  intro q hq different
  have member : q.column ∈ before ++ c :: after :=
    (choicePerm (before ++ c :: after) q hq).mem_iff.mp List.mem_cons_self
  have inside : q.column ∈ before ++ after := by
    rcases List.mem_append.mp member with h | h
    · exact List.mem_append_left _ h
    · rcases List.mem_cons.mp h with eq | h
      · exact False.elim (different eq)
      · exact List.mem_append_right _ h
  rw [zero q.column inside, gaussianArithmetic.zeroMul]
  cases q.odd <;> simp only [signed, Bool.false_eq_true, ↓reduceIte,
    gaussianArithmetic.negZero]


theorem blockDetLagrange (a b c d : Gaussian) :
    normSq (blockDet a b c d) =
      (normSq a + normSq c) * (normSq b + normSq d) -
      normSq (gaussianArithmetic.add
        (gaussianArithmetic.mul (conjugate a) b)
        (gaussianArithmetic.mul (conjugate c) d)) := by
  simp only [normSq, blockDet, conjugate, gaussianArithmetic]
  grind

theorem unitaryDetNorm (a b c d : Gaussian) (h : UnitaryBlock a b c d) :
    normSq (blockDet a b c d) = 1 := by
  rw [blockDetLagrange, h.columnLeft, h.columnRight, h.orthogonal, normSqZero]
  grind


section Gate
variable {V : Type} [DecidableEq V] {vertices : List V}

theorem gateSpectatorRow (e : Edge vertices) (a b c d : Gaussian) (x : V)
    (xl : x ≠ e.left) (xr : x ≠ e.right) (y : V) :
    twoCoordinateGate gaussianArithmetic e a b c d x y =
      if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
  simp only [twoCoordinateGate, xl, xr, false_and, ↓reduceIte]

theorem gateEntries (e : Edge vertices) (a b c d : Gaussian) :
    twoCoordinateGate gaussianArithmetic e a b c d e.left e.left = a ∧
    twoCoordinateGate gaussianArithmetic e a b c d e.left e.right = b ∧
    twoCoordinateGate gaussianArithmetic e a b c d e.right e.left = c ∧
    twoCoordinateGate gaussianArithmetic e a b c d e.right e.right = d := by
  simp [twoCoordinateGate, e.distinct, Ne.symm e.distinct]

/-- The complete local table with arbitrary leading occupied spectators.
The list of spectators may be empty; both occupied is the full minor on the
listed coordinates and neither occupied is their common spectator minor. -/
theorem prefixGateMagnitudes (e : Edge vertices) (a b c d : Gaussian)
    (h : UnitaryBlock a b c d) (xs : List V)
    (unique : xs.Nodup) (avoids : ∀ x ∈ xs, x ≠ e.left ∧ x ≠ e.right) :
    let M := twoCoordinateGate gaussianArithmetic e a b c d
    normSq (minorDet gaussianArithmetic M xs xs) = 1 ∧
    normSq (minorDet gaussianArithmetic M (xs ++ [e.left]) (xs ++ [e.left])) =
      1 - normSq b ∧
    normSq (minorDet gaussianArithmetic M (xs ++ [e.right]) (xs ++ [e.right])) =
      1 - normSq b ∧
    normSq (minorDet gaussianArithmetic M (xs ++ [e.left]) (xs ++ [e.right])) = normSq b ∧
    normSq (minorDet gaussianArithmetic M (xs ++ [e.right]) (xs ++ [e.left])) = normSq b ∧
    normSq (minorDet gaussianArithmetic M (xs ++ [e.left, e.right])
      (xs ++ [e.left, e.right])) = 1 := by
  let M := twoCoordinateGate gaussianArithmetic e a b c d
  have identity : ∀ x ∈ xs, ∀ y, M x y =
      if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
    intro x hx y
    exact gateSpectatorRow e a b c d x (avoids x hx).1 (avoids x hx).2 y
  have uniqueL : (xs ++ [e.left]).Nodup := by
    apply List.nodup_append.mpr
    refine ⟨unique, by simp, ?_⟩
    intro x hx y hy eq
    have yl : y = e.left := by simpa using hy
    exact (avoids x hx).1 (eq.trans yl)
  have uniqueR : (xs ++ [e.right]).Nodup := by
    apply List.nodup_append.mpr
    refine ⟨unique, by simp, ?_⟩
    intro x hx y hy eq
    have yr : y = e.right := by simpa using hy
    exact (avoids x hx).2 (eq.trans yr)
  have uniqueLR : (xs ++ [e.left, e.right]).Nodup := by
    apply List.nodup_append.mpr
    refine ⟨unique, by simp [e.distinct], ?_⟩
    intro x hx y hy eq
    have endpoints : y = e.left ∨ y = e.right := by simpa using hy
    rcases endpoints with yl | yr
    · exact (avoids x hx).1 (eq.trans yl)
    · exact (avoids x hx).2 (eq.trans yr)
  have empty := identityPrefixMinor M xs [] [] (by simpa using unique) identity
  have ll := identityPrefixMinor M xs [e.left] [e.left] uniqueL identity
  have rr := identityPrefixMinor M xs [e.right] [e.right] uniqueR identity
  have lr := identityPrefixMinor M xs [e.left] [e.right] uniqueR identity
  have rl := identityPrefixMinor M xs [e.right] [e.left] uniqueL identity
  have both := identityPrefixMinor M xs [e.left, e.right] [e.left, e.right] uniqueLR identity
  have entries := gateEntries e a b c d
  have weights := unitaryWeights a b c d h
  simp only [singletonMinor, gaussianTwoByTwo] at ll rr lr rl both
  change _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _
  dsimp only [M] at empty ll rr lr rl both
  rw [entries.1] at ll
  rw [entries.2.2.2] at rr
  rw [entries.2.1] at lr
  rw [entries.2.2.1] at rl
  rw [entries.1, entries.2.1, entries.2.2.1, entries.2.2.2] at both
  simp only [List.append_nil, emptyMinor] at empty
  exact ⟨by rw [empty]; exact normSqOne,
    by rw [ll]; exact weights.2.2,
    by rw [rr, weights.2.1]; exact weights.2.2,
    by rw [lr], by rw [rl]; exact weights.1,
    by rw [both]; exact unitaryDetNorm a b c d h⟩

end Gate
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
#print axioms lastChoiceMember
#print axioms scalarSumUniqueKey
#print axioms minorLastColumn
#print axioms selectedChoiceMember
#print axioms minorSelectedColumn
#print axioms blockDetLagrange
#print axioms unitaryDetNorm
#print axioms gateSpectatorRow
#print axioms gateEntries
#print axioms prefixGateMagnitudes
end MathCommons.ForestGateMagnitude




