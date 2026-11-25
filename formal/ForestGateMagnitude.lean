import ForestGateSupport

/-!
# Exact local unitary-gate magnitudes

All scalars are Gaussian rationals. `orderedGateMagnitudes` proves the complete
six-entry local table for the ACTUAL signed `ForestGateSupport.minorDet`:
neither/both occupied give 1; either one-particle stay gives 1-q; either exchange
gives q, where q = normSq b. `unitaryProbabilityBounds` proves 0 ≤ q ≤ 1.
The explicit occupied lists allow arbitrarily many spectators before, between,
and after the endpoints. Their combined list is required to be duplicate-free;
spectators must avoid the two endpoints. Empty lists and zero entry magnitudes
are allowed. Unitarity is given by explicit row/column norm equations and complex
column orthogonality, and the determinant norm is DERIVED from a Lagrange identity.

`interveningExchangeSign` and `reverseExchangeSign` retain the exact ordering
sign before squaring. `bothOccupiedMinor` accounts for cancellation of those signs.
`reverseGate` and `reverseUnitary` formally handle the opposite endpoint order by
replacing (e,a,b,c,d) with (reverseEdge e,d,c,b,a). `decomposeTwoEndpoints` and
`filterTwoEndpoints` formally expose the two ambient orders and four occupancy
cases. Applying those decompositions to arbitrary Boolean subset histories is
not packaged as one exported theorem here; the main table takes the resulting
occupied lists and their uniqueness/avoidance hypotheses explicitly.
`emptyFilteredMagnitude` and `fullFilteredMagnitude` give actual filtered empty
and full cases, with the ambient order explicitly supplied in the latter.

No Cauchy–Binet formula, matrix-library determinant identification, full-complex
scalar instance, product-gate probability identity or universal MOC result is
claimed. The zero-dimensional case needs no gate; Edge requires two distinct
vertices. The two-column Laplace and list-decomposition helper proofs were
contributed by the sibling protocol agent under the same GitHub owner; this is
collaboration, not independent review.
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


theorem minorTwoSupportedColumns {V : Type} (M : V → V → Gaussian)
    (r c₁ c₂ : V) (rs middle after : List V)
    (unique : (c₁ :: (middle ++ c₂ :: after)).Nodup)
    (zero : ∀ y ∈ middle ++ after, M r y = gaussianArithmetic.zero) :
    minorDet gaussianArithmetic M (r :: rs) (c₁ :: (middle ++ c₂ :: after)) =
      gaussianArithmetic.add
        (gaussianArithmetic.mul (M r c₁)
          (minorDet gaussianArithmetic M rs (middle ++ c₂ :: after)))
        (signed gaussianArithmetic (parity (c₁ :: middle))
          (gaussianArithmetic.mul (M r c₂)
            (minorDet gaussianArithmetic M rs (c₁ :: (middle ++ after))))) := by
  simp only [minorDet, columnChoices, List.map_cons, List.map_map, Function.comp_def, scalarSum]
  have tail : scalarSum gaussianArithmetic
      ((columnChoices (middle ++ c₂ :: after)).map
        (fun q => signed gaussianArithmetic (!q.odd)
          (gaussianArithmetic.mul (M r q.column)
            (minorDet gaussianArithmetic M rs (c₁ :: q.remaining))))) =
      signed gaussianArithmetic (parity (c₁ :: middle))
        (gaussianArithmetic.mul (M r c₂)
          (minorDet gaussianArithmetic M rs (c₁ :: (middle ++ after)))) := by
    apply scalarSumUniqueKey (columnChoices (middle ++ c₂ :: after)) ColumnChoice.column
      (by simpa only [columnChoicesColumns] using (List.nodup_cons.mp unique).2)
      ⟨c₂, middle ++ after, parity middle⟩ (selectedChoiceMember middle after c₂)
    intro q hq different
    have member : q.column ∈ middle ++ c₂ :: after :=
      (choicePerm (middle ++ c₂ :: after) q hq).mem_iff.mp List.mem_cons_self
    have inside : q.column ∈ middle ++ after := by
      rcases List.mem_append.mp member with h | h
      · exact List.mem_append_left _ h
      · rcases List.mem_cons.mp h with eq | h
        · exact False.elim (different eq)
        · exact List.mem_append_right _ h
    rw [zero q.column inside, gaussianArithmetic.zeroMul]
    cases q.odd <;> simp only [signed, Bool.not_true, Bool.not_false,
      Bool.false_eq_true, ↓reduceIte, gaussianArithmetic.negZero]
  rw [tail]
  rfl


/-- Two distinct members occur in one order or the other. No Nodup hypothesis is needed. -/
theorem decomposeTwoEndpoints {V : Type} (xs : List V) (left right : V)
    (different : left ≠ right) (hasLeft : left ∈ xs) (hasRight : right ∈ xs) :
    (∃ before middle after, xs = before ++ left :: (middle ++ right :: after)) ∨
    (∃ before middle after, xs = before ++ right :: (middle ++ left :: after)) := by
  obtain ⟨before, tail, rfl⟩ := List.append_of_mem hasLeft
  rcases List.mem_append.mp hasRight with inBefore | inTail
  · obtain ⟨pre, middle, eqBefore⟩ := List.append_of_mem inBefore
    right
    refine ⟨pre, middle, tail, ?_⟩
    rw [eqBefore]
    simp only [List.append_assoc, List.cons_append]
  · rcases List.mem_cons.mp inTail with equal | inTail
    · exact False.elim (different equal.symm)
    · obtain ⟨middle, after, eqTail⟩ := List.append_of_mem inTail
      left
      exact ⟨before, middle, after, by rw [eqTail]⟩

/-- Filtering a two-endpoint decomposition separates the four occupancy cases. -/
theorem filterTwoEndpoints {V : Type} (p : V → Bool) (left right : V)
    (before middle after : List V) :
    (before ++ left :: (middle ++ right :: after)).filter p =
      before.filter p ++ (if p left then [left] else []) ++ middle.filter p ++
        (if p right then [right] else []) ++ after.filter p := by
  cases hl : p left <;> cases hr : p right <;>
    simp [List.filter_append, hl, hr, List.append_assoc]


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


theorem gateEndpointSpectatorZero (e : Edge vertices) (a b c d : Gaussian)
    (x : V) (xl : x ≠ e.left) (xr : x ≠ e.right) :
    twoCoordinateGate gaussianArithmetic e a b c d e.left x = gaussianArithmetic.zero ∧
    twoCoordinateGate gaussianArithmetic e a b c d e.right x = gaussianArithmetic.zero := by
  simp [twoCoordinateGate, xl, xr, Ne.symm xl, Ne.symm xr, e.distinct, Ne.symm e.distinct]

/-- Exchange from the right endpoint to the left crosses any number of occupied
spectators, with the exact parity sign in the ordered minor. Trailing occupied
spectators are arbitrary as well. -/
theorem interveningExchangeSign (e : Edge vertices) (a b c d : Gaussian)
    (middle after : List V) (unique : (middle ++ after).Nodup)
    (avoids : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right) :
    minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (e.left :: (middle ++ after)) (middle ++ e.right :: after) =
      signed gaussianArithmetic (parity middle) b := by
  let M := twoCoordinateGate gaussianArithmetic e a b c d
  have identity : ∀ x ∈ middle ++ after, ∀ y, M x y =
      if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
    intro x hx y
    exact gateSpectatorRow e a b c d x (avoids x hx).1 (avoids x hx).2 y
  have nr : e.right ∉ middle ++ after := by
    intro hx; exact (avoids e.right hx).2 rfl
  have inserted : (middle ++ e.right :: after).Nodup :=
    List.perm_middle.nodup_iff.mpr (List.nodup_cons.mpr ⟨nr, unique⟩)
  rw [minorSelectedColumn _ e.left e.right (middle ++ after) middle after inserted]
  · have spectators := identityPrefixMinor M (middle ++ after) [] []
      (by simpa using unique) identity
    simp only [List.append_nil, emptyMinor] at spectators
    change signed gaussianArithmetic (parity middle)
      (gaussianArithmetic.mul (M e.left e.right)
        (minorDet gaussianArithmetic M (middle ++ after) (middle ++ after))) = _
    rw [spectators, gaussianMulOne]
    exact congrArg (signed gaussianArithmetic (parity middle)) (gateEntries e a b c d).2.1
  · intro x hx
    exact (gateEndpointSpectatorZero e a b c d x (avoids x hx).1 (avoids x hx).2).1

theorem interveningExchangeNorm (e : Edge vertices) (a b c d : Gaussian)
    (middle after : List V) (unique : (middle ++ after).Nodup)
    (avoids : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (e.left :: (middle ++ after)) (middle ++ e.right :: after)) = normSq b := by
  rw [interveningExchangeSign e a b c d middle after unique avoids, normSqSigned]

theorem signedNot (odd : Bool) (z : Gaussian) :
    signed gaussianArithmetic (!odd) z =
      gaussianArithmetic.neg (signed gaussianArithmetic odd z) := by
  cases odd <;> simp only [signed, Bool.not_false, Bool.not_true,
    Bool.false_eq_true, ↓reduceIte]
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

/-- The reverse exchange has the exact parity sign, with arbitrarily many
intervening and trailing occupied spectators. No unitarity is used here. -/
theorem reverseExchangeSign (e : Edge vertices) (a b c d : Gaussian)
    (middle after : List V) (unique : (middle ++ after).Nodup)
    (avoids : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right) :
    minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (middle ++ e.right :: after) (e.left :: (middle ++ after)) =
      signed gaussianArithmetic (parity middle) c := by
  let M := twoCoordinateGate gaussianArithmetic e a b c d
  induction middle with
  | nil =>
    simp only [List.nil_append] at unique avoids ⊢
    change minorDet gaussianArithmetic M (e.right :: after) (e.left :: after) = c
    rw [minorFirstColumn M e.right e.left after after]
    · have identity : ∀ x ∈ after, ∀ y, M x y =
          if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
        intro x hx y
        exact gateSpectatorRow e a b c d x (avoids x hx).1 (avoids x hx).2 y
      have tail := identityPrefixMinor M after [] [] (by simpa using unique) identity
      simp only [List.append_nil, emptyMinor] at tail
      rw [tail, gaussianMulOne]
      exact (gateEntries e a b c d).2.2.1
    · intro x hx
      exact (gateEndpointSpectatorZero e a b c d x (avoids x hx).1 (avoids x hx).2).2
  | cons x xs ih =>
    have nu := List.nodup_cons.mp unique
    have ax := avoids x List.mem_cons_self
    have atail : ∀ y ∈ xs ++ after, y ≠ e.left ∧ y ≠ e.right := by
      intro y hy; exact avoids y (List.mem_cons_of_mem x hy)
    have nl : e.left ∉ x :: (xs ++ after) := by
      intro hx; exact (avoids e.left hx).1 rfl
    have insert : ([e.left] ++ x :: (xs ++ after)).Nodup :=
      List.nodup_cons.mpr ⟨nl, unique⟩
    change minorDet gaussianArithmetic M (x :: (xs ++ e.right :: after))
      ([e.left] ++ x :: (xs ++ after)) = _
    rw [minorSelectedColumn M x x (xs ++ e.right :: after) [e.left] (xs ++ after) insert]
    · have diag : M x x = gaussianArithmetic.one := by
        simpa using gateSpectatorRow e a b c d x ax.1 ax.2 x
      rw [diag]
      simp only [gaussianOneMul]
      change gaussianArithmetic.neg (minorDet gaussianArithmetic M
        (xs ++ e.right :: after) (e.left :: (xs ++ after))) =
        signed gaussianArithmetic (!(parity xs)) c
      rw [ih nu.2 atail, signedNot]
    · intro y hy
      have different : x ≠ y := by
        intro eq
        rcases List.mem_cons.mp hy with yl | ht
        · exact ax.1 (eq.trans yl)
        · exact nu.1 (eq ▸ ht)
      change twoCoordinateGate gaussianArithmetic e a b c d x y = _
      simp only [gateSpectatorRow e a b c d x ax.1 ax.2 y, different, ↓reduceIte]


theorem reverseExchangeNorm (e : Edge vertices) (a b c d : Gaussian)
    (middle after : List V) (unique : (middle ++ after).Nodup)
    (avoids : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (middle ++ e.right :: after) (e.left :: (middle ++ after))) = normSq c := by
  rw [reverseExchangeSign e a b c d middle after unique avoids, normSqSigned]


theorem gateSpectatorMinor (e : Edge vertices) (a b c d : Gaussian)
    (xs : List V) (unique : xs.Nodup)
    (avoids : ∀ x ∈ xs, x ≠ e.left ∧ x ≠ e.right) :
    minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d) xs xs =
      gaussianArithmetic.one := by
  have identity : ∀ x ∈ xs, ∀ y,
      twoCoordinateGate gaussianArithmetic e a b c d x y =
      if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
    intro x hx y
    exact gateSpectatorRow e a b c d x (avoids x hx).1 (avoids x hx).2 y
  simpa only [List.append_nil, emptyMinor] using
    identityPrefixMinor (twoCoordinateGate gaussianArithmetic e a b c d) xs [] []
      (by simpa using unique) identity

theorem leftStayMinor (e : Edge vertices) (a b c d : Gaussian)
    (xs : List V) (unique : xs.Nodup)
    (avoids : ∀ x ∈ xs, x ≠ e.left ∧ x ≠ e.right) :
    minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (e.left :: xs) (e.left :: xs) = a := by
  rw [minorFirstColumn]
  · rw [gateSpectatorMinor e a b c d xs unique avoids, gaussianMulOne]
    exact (gateEntries e a b c d).1
  · intro x hx
    exact (gateEndpointSpectatorZero e a b c d x (avoids x hx).1 (avoids x hx).2).1

theorem rightStayMinor (e : Edge vertices) (a b c d : Gaussian)
    (middle after : List V) (unique : (middle ++ after).Nodup)
    (avoids : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right) :
    minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (middle ++ e.right :: after) (middle ++ e.right :: after) = d := by
  let M := twoCoordinateGate gaussianArithmetic e a b c d
  have nr : e.right ∉ middle ++ after := by intro hx; exact (avoids e.right hx).2 rfl
  have inserted : (middle ++ e.right :: after).Nodup :=
    List.perm_middle.nodup_iff.mpr (List.nodup_cons.mpr ⟨nr, unique⟩)
  have identity : ∀ x ∈ middle, ∀ y, M x y =
      if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
    intro x hx y
    have ax := avoids x (List.mem_append_left _ hx)
    exact gateSpectatorRow e a b c d x ax.1 ax.2 y
  rw [identityPrefixMinor M middle (e.right :: after) (e.right :: after) inserted identity]
  rw [minorFirstColumn M e.right e.right after after]
  · have ua := (List.sublist_append_right middle after).nodup unique
    have aa : ∀ x ∈ after, x ≠ e.left ∧ x ≠ e.right := by
      intro x hx; exact avoids x (List.mem_append_right _ hx)
    rw [gateSpectatorMinor e a b c d after ua aa, gaussianMulOne]
    exact (gateEntries e a b c d).2.2.2
  · intro x hx
    have ax := avoids x (List.mem_append_right _ hx)
    exact (gateEndpointSpectatorZero e a b c d x ax.1 ax.2).2

theorem oppositeSignsProduct (odd : Bool) (b c : Gaussian) :
    signed gaussianArithmetic (!odd)
      (gaussianArithmetic.mul b (signed gaussianArithmetic odd c)) =
      gaussianArithmetic.neg (gaussianArithmetic.mul b c) := by
  cases odd <;> simp only [signed, Bool.not_false, Bool.not_true,
    Bool.false_eq_true, ↓reduceIte]
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

/-- Double occupation has the block determinant for arbitrary intervening and
trailing spectators: the two permutation signs cancel in the crossed term. -/
theorem bothOccupiedMinor (e : Edge vertices) (a b c d : Gaussian)
    (middle after : List V) (unique : (middle ++ after).Nodup)
    (avoids : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right) :
    minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (e.left :: (middle ++ e.right :: after))
      (e.left :: (middle ++ e.right :: after)) = blockDet a b c d := by
  have nr : e.right ∉ middle ++ after := by intro hx; exact (avoids e.right hx).2 rfl
  have inserted : (middle ++ e.right :: after).Nodup :=
    List.perm_middle.nodup_iff.mpr (List.nodup_cons.mpr ⟨nr, unique⟩)
  have nl : e.left ∉ middle ++ e.right :: after := by
    intro hx
    rcases List.mem_append.mp hx with hm | ha
    · exact (avoids e.left (List.mem_append_left _ hm)).1 rfl
    · rcases List.mem_cons.mp ha with eq | ha
      · exact e.distinct eq
      · exact (avoids e.left (List.mem_append_right _ ha)).1 rfl
  rw [minorTwoSupportedColumns _ e.left e.left e.right
    (middle ++ e.right :: after) middle after (List.nodup_cons.mpr ⟨nl, inserted⟩)]
  · rw [rightStayMinor e a b c d middle after unique avoids,
      reverseExchangeSign e a b c d middle after unique avoids,
      (gateEntries e a b c d).1, (gateEntries e a b c d).2.1]
    simp only [parity, oppositeSignsProduct, blockDet]
  · intro x hx
    exact (gateEndpointSpectatorZero e a b c d x (avoids x hx).1 (avoids x hx).2).1


/-- All six local squared probabilities, with arbitrary occupied spectators
before, between and after the two ordered endpoints. The hypothesis describes
the occupied lists explicitly; conversion from a filtered ambient enumeration
is a separate combinatorial wrapper, not an assumed determinant identity. -/
theorem orderedGateMagnitudes (e : Edge vertices) (a b c d : Gaussian)
    (h : UnitaryBlock a b c d) (before middle after : List V)
    (unique : (before ++ e.left :: (middle ++ e.right :: after)).Nodup)
    (avoids : ∀ x ∈ before ++ (middle ++ after), x ≠ e.left ∧ x ≠ e.right) :
    let M := twoCoordinateGate gaussianArithmetic e a b c d
    let empty := before ++ (middle ++ after)
    let left := before ++ e.left :: (middle ++ after)
    let right := before ++ (middle ++ e.right :: after)
    let both := before ++ e.left :: (middle ++ e.right :: after)
    normSq (minorDet gaussianArithmetic M empty empty) = 1 ∧
    normSq (minorDet gaussianArithmetic M left left) = 1 - normSq b ∧
    normSq (minorDet gaussianArithmetic M right right) = 1 - normSq b ∧
    normSq (minorDet gaussianArithmetic M left right) = normSq b ∧
    normSq (minorDet gaussianArithmetic M right left) = normSq b ∧
    normSq (minorDet gaussianArithmetic M both both) = 1 := by
  let M := twoCoordinateGate gaussianArithmetic e a b c d
  have removeR : (middle ++ after).Sublist (middle ++ e.right :: after) :=
    (List.Sublist.refl middle).append ((List.Sublist.refl after).cons e.right)
  have removeBoth := removeR.cons e.left
  have removeLeft := (List.Sublist.refl (middle ++ e.right :: after)).cons e.left
  have removeRight := removeR.cons₂ e.left
  have coreUnique := (List.sublist_append_right before
    (e.left :: (middle ++ e.right :: after))).nodup unique
  have u := removeBoth.nodup coreUnique
  have u0 := ((List.Sublist.refl before).append removeBoth).nodup unique
  have ul := ((List.Sublist.refl before).append removeRight).nodup unique
  have ur := ((List.Sublist.refl before).append removeLeft).nodup unique
  have ac : ∀ x ∈ middle ++ after, x ≠ e.left ∧ x ≠ e.right := by
    intro x hx; exact avoids x (List.mem_append_right _ hx)
  have identity : ∀ x ∈ before, ∀ y, M x y =
      if x = y then gaussianArithmetic.one else gaussianArithmetic.zero := by
    intro x hx y
    have ax := avoids x (List.mem_append_left _ hx)
    exact gateSpectatorRow e a b c d x ax.1 ax.2 y
  have none := identityPrefixMinor M before (middle ++ after) (middle ++ after) u0 identity
  have ll := identityPrefixMinor M before (e.left :: (middle ++ after))
    (e.left :: (middle ++ after)) ul identity
  have rr := identityPrefixMinor M before (middle ++ e.right :: after)
    (middle ++ e.right :: after) ur identity
  have lr := identityPrefixMinor M before (e.left :: (middle ++ after))
    (middle ++ e.right :: after) ur identity
  have rl := identityPrefixMinor M before (middle ++ e.right :: after)
    (e.left :: (middle ++ after)) ul identity
  have bb := identityPrefixMinor M before (e.left :: (middle ++ e.right :: after))
    (e.left :: (middle ++ e.right :: after)) unique identity
  have weights := unitaryWeights a b c d h
  change _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _
  constructor
  · rw [none, gateSpectatorMinor e a b c d (middle ++ after) u ac]
    exact normSqOne
  constructor
  · rw [ll, leftStayMinor e a b c d (middle ++ after) u ac]
    exact weights.2.2
  constructor
  · rw [rr, rightStayMinor e a b c d middle after u ac, weights.2.1]
    exact weights.2.2
  constructor
  · rw [lr]
    exact interveningExchangeNorm e a b c d middle after u ac
  constructor
  · rw [rl, reverseExchangeNorm e a b c d middle after u ac]
    exact weights.1
  · rw [bb, bothOccupiedMinor e a b c d middle after u ac]
    exact unitaryDetNorm a b c d h

theorem normSqNonnegative (z : Gaussian) : 0 ≤ normSq z := by
  have square : ∀ x : Rat, 0 ≤ x * x := by
    intro x
    by_cases hx : 0 ≤ x
    · exact Rat.mul_nonneg hx hx
    · have hn : 0 ≤ -x := by grind
      have hh := Rat.mul_nonneg hn hn
      grind
  exact Rat.add_nonneg (square z.1) (square z.2)

theorem unitaryProbabilityBounds (a b c d : Gaussian) (h : UnitaryBlock a b c d) :
    0 ≤ normSq b ∧ normSq b ≤ 1 := by
  have ha := normSqNonnegative a
  have hb := normSqNonnegative b
  have hr := h.rowLeft
  grind


/-- Reversing the endpoint names reverses both block row and column order. -/
def reverseEdge (e : Edge vertices) : Edge vertices where
  nodup := e.nodup
  complete := e.complete
  left := e.right
  right := e.left
  distinct := Ne.symm e.distinct
  rest := e.rest
  enumerate := e.enumerate.trans (List.Perm.swap e.right e.left e.rest)
  avoids := by intro x hx; exact ⟨(e.avoids x hx).2, (e.avoids x hx).1⟩

theorem reverseGate (e : Edge vertices) (a b c d : Gaussian) :
    twoCoordinateGate gaussianArithmetic (reverseEdge e) d c b a =
      twoCoordinateGate gaussianArithmetic e a b c d := by
  funext x y
  simp only [twoCoordinateGate, reverseEdge]
  have distinct := e.distinct
  grind

theorem reverseUnitary (a b c d : Gaussian) (h : UnitaryBlock a b c d) :
    UnitaryBlock d c b a := by
  have h1 := h.rowLeft
  have h2 := h.rowRight
  have h3 := h.columnLeft
  have h4 := h.columnRight
  have hr := congrArg Prod.fst h.orthogonal
  have hi := congrArg Prod.snd h.orthogonal
  constructor
  · grind
  · grind
  · grind
  · grind
  · apply Prod.ext <;> simp only [gaussianArithmetic, conjugate] at * <;> grind


theorem emptyFilteredMagnitude (e : Edge vertices) (a b c d : Gaussian) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter (fun _ => false)) (vertices.filter (fun _ => false))) = 1 := by
  have empty : vertices.filter (fun _ => false) = [] := by simp
  rw [empty, emptyMinor]
  exact normSqOne

/-- Full-subset corollary with the ambient order explicitly supplied. The general
existence and reversal tools are `decomposeTwoEndpoints` and `reverseGate`. -/
theorem fullFilteredMagnitude (e : Edge vertices) (a b c d : Gaussian)
    (h : UnitaryBlock a b c d) (before middle after : List V)
    (order : vertices = before ++ e.left :: (middle ++ e.right :: after))
    (avoids : ∀ x ∈ before ++ (middle ++ after), x ≠ e.left ∧ x ≠ e.right) :
    normSq (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
      (vertices.filter (fun _ => true)) (vertices.filter (fun _ => true))) = 1 := by
  have unique : (before ++ e.left :: (middle ++ e.right :: after)).Nodup :=
    order ▸ e.nodup
  have table := orderedGateMagnitudes e a b c d h before middle after unique avoids
  have full : vertices.filter (fun _ => true) = vertices :=
    List.filter_eq_self.mpr (by intro x hx; rfl)
  rw [full]
  exact (congrArg (fun xs : List V => normSq (minorDet gaussianArithmetic
    (twoCoordinateGate gaussianArithmetic e a b c d) xs xs)) order).trans table.2.2.2.2.2

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
#print axioms minorTwoSupportedColumns
#print axioms decomposeTwoEndpoints
#print axioms filterTwoEndpoints
#print axioms gateSpectatorRow
#print axioms gateEntries
#print axioms prefixGateMagnitudes
#print axioms gateEndpointSpectatorZero
#print axioms interveningExchangeSign
#print axioms interveningExchangeNorm
#print axioms signedNot
#print axioms reverseExchangeSign
#print axioms reverseExchangeNorm
#print axioms gateSpectatorMinor
#print axioms leftStayMinor
#print axioms rightStayMinor
#print axioms oppositeSignsProduct
#print axioms bothOccupiedMinor
#print axioms orderedGateMagnitudes
#print axioms normSqNonnegative
#print axioms unitaryProbabilityBounds
#print axioms reverseGate
#print axioms reverseUnitary
#print axioms emptyFilteredMagnitude
#print axioms fullFilteredMagnitude

end MathCommons.ForestGateMagnitude
