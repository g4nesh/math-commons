import ForestProductMinor

/-!
Actual Gaussian-rational ordered determinants are alternating in arbitrary row
positions. Two successive signed column selections are antisymmetric; this
proves adjacent exchange, separated exchange, duplicate-row cancellation and
the parity of moving one row across a spectator block. The actual both-occupied
gate multiplication formula now follows without unproved alternation premises.
All column lists, including rectangular and empty cases, are allowed. This
module does not group arbitrary-order product tuples into canonical subsets.
-/

namespace MathCommons.ForestMinorAlternation
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor

theorem gaussianNegNeg (a : Gaussian) : gaussianArithmetic.neg (gaussianArithmetic.neg a) = a := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianNegAdd (a b : Gaussian) :
    gaussianArithmetic.neg (gaussianArithmetic.add a b) =
      gaussianArithmetic.add (gaussianArithmetic.neg a) (gaussianArithmetic.neg b) := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianMulNeg (a b : Gaussian) :
    gaussianArithmetic.mul a (gaussianArithmetic.neg b) =
      gaussianArithmetic.neg (gaussianArithmetic.mul a b) := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem signedNeg (odd : Bool) (a : Gaussian) :
    signed gaussianArithmetic odd (gaussianArithmetic.neg a) =
      gaussianArithmetic.neg (signed gaussianArithmetic odd a) := by
  cases odd <;> rfl

theorem gaussianSumNeg {α : Type} (xs : List α) (f : α → Gaussian) :
    scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.neg (f x)) =
      gaussianArithmetic.neg (scalarSum gaussianArithmetic (xs.map f)) := by
  induction xs with
  | nil => exact gaussianArithmetic.negZero.symm
  | cons x xs ih => simp only [List.map_cons, scalarSum, ih, gaussianNegAdd]

def doubleChoiceSum {V : Type} (cols : List V) (f : V → V → List V → Gaussian) : Gaussian :=
  scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
    signed gaussianArithmetic q.odd (scalarSum gaussianArithmetic
      ((columnChoices q.remaining).map fun r => signed gaussianArithmetic r.odd
        (f q.column r.column r.remaining))))

theorem doubleChoiceCons {V : Type} (c : V) (cols : List V)
    (f : V → V → List V → Gaussian) :
    doubleChoiceSum (c :: cols) f = gaussianArithmetic.add
      (scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
        signed gaussianArithmetic q.odd (f c q.column q.remaining)))
      (gaussianArithmetic.add
        (gaussianArithmetic.neg (scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
          signed gaussianArithmetic q.odd (f q.column c q.remaining))))
        (doubleChoiceSum cols (fun x y rest => f x y (c :: rest)))) := by
  have noSign (z : Gaussian) : signed gaussianArithmetic false z = z := rfl
  simp only [doubleChoiceSum, columnChoices, List.map_cons, List.map_map,
    scalarSum, Function.comp_def, noSign, signedNot, signedAdd, gaussianSumNeg,
    signedNeg, gaussianNegNeg, gaussianSumAdd]


theorem doubleChoiceSwap {V : Type} (cols : List V) (f : V → V → List V → Gaussian) :
    doubleChoiceSum cols (fun x y rest => f y x rest) = gaussianArithmetic.neg (doubleChoiceSum cols f) := by
  induction cols generalizing f with
  | nil => exact gaussianArithmetic.negZero.symm
  | cons c cols ih =>
    rw [doubleChoiceCons, doubleChoiceCons, ih (fun x y rest => f x y (c :: rest))]
    apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem rowDetTwoChoices {V : Type} (f g : V → Gaussian) (rows : List (V → Gaussian)) (cols : List V) :
    rowDet (f :: g :: rows) cols = doubleChoiceSum cols (fun x y rest =>
      gaussianArithmetic.mul (f x) (gaussianArithmetic.mul (g y) (rowDet rows rest))) := by
  simp only [rowDet, doubleChoiceSum]
  apply gaussianSumCongr
  intro q _
  apply congrArg (signed gaussianArithmetic q.odd)
  rw [← gaussianSumScale]
  apply gaussianSumCongr
  intro r _
  exact (signedScale r.odd (f q.column) _).symm

/-- Swapping the first two row vectors changes the actual Laplace determinant's
sign, for every remaining row list and column list, including rectangular ones. -/
theorem rowDetHeadSwap {V : Type} (f g : V → Gaussian) (rows : List (V → Gaussian)) (cols : List V) :
    rowDet (g :: f :: rows) cols = gaussianArithmetic.neg (rowDet (f :: g :: rows) cols) := by
  rw [rowDetTwoChoices, rowDetTwoChoices]
  rw [← doubleChoiceSwap]
  unfold doubleChoiceSum
  apply gaussianSumCongr
  intro q _
  apply congrArg (signed gaussianArithmetic q.odd)
  apply gaussianSumCongr
  intro r _
  apply congrArg (signed gaussianArithmetic r.odd)
  dsimp only
  rw [← gaussianMulAssoc, gaussianMulComm (g q.column) (f r.column), gaussianMulAssoc]

/-- An adjacent exchange at any row position changes the determinant's sign. -/
theorem rowDetAdjacentSwap {V : Type} (before rows : List (V → Gaussian))
    (f g : V → Gaussian) (cols : List V) :
    rowDet (before ++ g :: f :: rows) cols =
      gaussianArithmetic.neg (rowDet (before ++ f :: g :: rows) cols) := by
  induction before generalizing cols with
  | nil => exact rowDetHeadSwap f g rows cols
  | cons h before ih =>
    simp only [List.cons_append, rowDet]
    rw [← gaussianSumNeg]
    apply gaussianSumCongr
    intro q _
    rw [ih q.remaining, gaussianMulNeg, signedNeg]

/-- Adjacent duplicate rows cancel without any nonzero or dimension hypothesis. -/
theorem rowDetAdjacentDuplicate {V : Type} (before rows : List (V → Gaussian))
    (f : V → Gaussian) (cols : List V) :
    rowDet (before ++ f :: f :: rows) cols = gaussianArithmetic.zero := by
  have sign := rowDetAdjacentSwap before rows f f cols
  apply Prod.ext
  · have h := congrArg Prod.fst sign
    simp only [gaussianArithmetic] at h ⊢
    grind
  · have h := congrArg Prod.snd sign
    simp only [gaussianArithmetic] at h ⊢
    grind


/-- A signed relation between row-list determinants persists under a common
prefix of arbitrary row vectors. -/
theorem rowDetPrefixNeg {V : Type} (before left right : List (V → Gaussian))
    (relation : ∀ cols, rowDet left cols = gaussianArithmetic.neg (rowDet right cols))
    (cols : List V) : rowDet (before ++ left) cols = gaussianArithmetic.neg (rowDet (before ++ right) cols) := by
  induction before generalizing cols with
  | nil => exact relation cols
  | cons f before ih =>
    simp only [List.cons_append, rowDet]
    rw [← gaussianSumNeg]
    apply gaussianSumCongr
    intro q _
    rw [ih q.remaining, gaussianMulNeg, signedNeg]

theorem rowDetSeparatedHeadSwap {V : Type} (middle after : List (V → Gaussian))
    (f g : V → Gaussian) (cols : List V) :
    rowDet (g :: (middle ++ f :: after)) cols =
      gaussianArithmetic.neg (rowDet (f :: (middle ++ g :: after)) cols) := by
  induction middle generalizing cols with
  | nil => exact rowDetHeadSwap f g after cols
  | cons h middle ih =>
    rw [List.cons_append, rowDetHeadSwap h g (middle ++ f :: after) cols]
    have inside := rowDetPrefixNeg [h] (g :: (middle ++ f :: after))
      (f :: (middle ++ g :: after)) ih cols
    simp only [List.cons_append, List.nil_append] at inside ⊢
    rw [inside, gaussianNegNeg]
    exact rowDetHeadSwap f h (middle ++ g :: after) cols

/-- Exchanging any two row positions has negative sign, regardless of their
separation. This preserves the supplied order of all intervening rows. -/
theorem rowDetExchange {V : Type} (before middle after : List (V → Gaussian))
    (f g : V → Gaussian) (cols : List V) :
    rowDet (before ++ g :: (middle ++ f :: after)) cols =
      gaussianArithmetic.neg (rowDet (before ++ f :: (middle ++ g :: after)) cols) :=
  rowDetPrefixNeg before _ _ (rowDetSeparatedHeadSwap middle after f g) cols

theorem rowDetDuplicate {V : Type} (before middle after : List (V → Gaussian))
    (f : V → Gaussian) (cols : List V) :
    rowDet (before ++ f :: (middle ++ f :: after)) cols = gaussianArithmetic.zero := by
  have sign := rowDetExchange before middle after f f cols
  apply Prod.ext
  · have h := congrArg Prod.fst sign
    simp only [gaussianArithmetic] at h ⊢
    grind
  · have h := congrArg Prod.snd sign
    simp only [gaussianArithmetic] at h ⊢
    grind

/-- Two distinct selected row indices with equal row vectors force the actual
minor to vanish, without requiring a square or nonsingular matrix. -/
theorem minorEqualRowsZero {V : Type} (M : V → V → Gaussian) (rows cols : List V)
    (r s : V) (different : r ≠ s) (mr : r ∈ rows) (ms : s ∈ rows) (same : M r = M s) :
    minorDet gaussianArithmetic M rows cols = gaussianArithmetic.zero := by
  rw [minorAsRows]
  rcases decomposeTwoEndpoints rows r s different mr ms with
    ⟨before, middle, after, eq⟩ | ⟨before, middle, after, eq⟩
  · rw [eq]
    simp only [List.map_append, List.map_cons, same]
    exact rowDetDuplicate _ _ _ _ _
  · rw [eq]
    simp only [List.map_append, List.map_cons, same]
    exact rowDetDuplicate _ _ _ _ _

theorem minorDuplicateRow {V : Type} [DecidableEq V] (M : V → V → Gaussian)
    (rows cols : List V) (r s : V) (different : r ≠ s) (mr : r ∈ rows) (ms : s ∈ rows) :
    minorDet gaussianArithmetic (replaceRow M s (M r)) rows cols = gaussianArithmetic.zero := by
  apply minorEqualRowsZero _ rows cols r s different mr ms
  funext j
  simp [replaceRow, different]


theorem twoEndpointAvoids {V : Type} (before middle after : List V) (r s : V)
    (unique : (before ++ r :: (middle ++ s :: after)).Nodup) :
    ∀ x ∈ before ++ middle ++ after, x ≠ r ∧ x ≠ s := by
  obtain ⟨_, tailUnique, cross⟩ := List.nodup_append.mp unique
  obtain ⟨rAbsent, restUnique⟩ := List.nodup_cons.mp tailUnique
  obtain ⟨_, lastUnique, middleCross⟩ := List.nodup_append.mp restUnique
  have sAbsent := (List.nodup_cons.mp lastUnique).1
  intro x hx
  simp only [List.mem_append] at hx
  rcases hx with (hb | hm) | ha
  · exact ⟨cross x hb r List.mem_cons_self,
      cross x hb s (by simp)⟩
  · constructor
    · intro eq
      subst x
      exact rAbsent (List.mem_append.mpr (Or.inl hm))
    · exact middleCross x hm s List.mem_cons_self
  · constructor
    · intro eq
      subst x
      exact rAbsent (by simp [ha])
    · intro eq
      subst x
      exact sAbsent ha

/-- Swapping two selected matrix rows has the actual negative determinant sign.
The duplicate-free row enumeration ensures no further copies are changed. -/
theorem minorExchangeRows {V : Type} [DecidableEq V] (M : V → V → Gaussian)
    (rows cols : List V) (unique : rows.Nodup)
    (r s : V) (different : r ≠ s) (mr : r ∈ rows) (ms : s ∈ rows) :
    minorDet gaussianArithmetic (replaceRow (replaceRow M r (M s)) s (M r)) rows cols =
      gaussianArithmetic.neg (minorDet gaussianArithmetic M rows cols) := by
  let N := replaceRow (replaceRow M r (M s)) s (M r)
  have atR : N r = M s := by funext j; simp [N, replaceRow, different]
  have atS : N s = M r := by funext j; simp [N, replaceRow]
  have elsewhere : ∀ x, x ≠ r → x ≠ s → N x = M x := by
    intro x hr hs
    funext j
    simp [N, replaceRow, hr, hs]
  change minorDet gaussianArithmetic N rows cols = _
  rw [minorAsRows, minorAsRows]
  rcases decomposeTwoEndpoints rows r s different mr ms with
    ⟨before, middle, after, eq⟩ | ⟨before, middle, after, eq⟩
  · have avoids := twoEndpointAvoids before middle after r s (eq ▸ unique)
    have mb : before.map N = before.map M := by
      apply List.map_congr_left
      intro x hx
      have a := avoids x (by simp [hx])
      exact elsewhere x a.1 a.2
    have mm : middle.map N = middle.map M := by
      apply List.map_congr_left
      intro x hx
      have a := avoids x (by simp [hx])
      exact elsewhere x a.1 a.2
    have ma : after.map N = after.map M := by
      apply List.map_congr_left
      intro x hx
      have a := avoids x (by simp [hx])
      exact elsewhere x a.1 a.2
    rw [eq]
    simp only [List.map_append, List.map_cons, mb, mm, ma, atR, atS]
    exact rowDetExchange _ _ _ _ _ _
  · have avoids := twoEndpointAvoids before middle after s r (eq ▸ unique)
    have mb : before.map N = before.map M := by
      apply List.map_congr_left
      intro x hx
      have a := avoids x (by simp [hx])
      exact elsewhere x a.2 a.1
    have mm : middle.map N = middle.map M := by
      apply List.map_congr_left
      intro x hx
      have a := avoids x (by simp [hx])
      exact elsewhere x a.2 a.1
    have ma : after.map N = after.map M := by
      apply List.map_congr_left
      intro x hx
      have a := avoids x (by simp [hx])
      exact elsewhere x a.2 a.1
    rw [eq]
    simp only [List.map_append, List.map_cons, mb, mm, ma, atR, atS]
    exact rowDetExchange _ _ _ _ _ _

/-- The both-occupied gate update now has no unproved alternation premises. -/
theorem gateMinorBoth {V : Type} [DecidableEq V] {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (rows cols : List V)
    (unique : rows.Nodup) (left : e.left ∈ rows) (right : e.right ∈ rows) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) rows cols =
      gaussianArithmetic.mul (blockDet a b c d) (minorDet gaussianArithmetic M rows cols) := by
  apply gateMinorBothOfAlternation e a b c d M rows cols unique left right
  · exact minorDuplicateRow M rows cols e.left e.right e.distinct left right
  · exact minorDuplicateRow M rows cols e.right e.left (Ne.symm e.distinct) right left
  · exact minorExchangeRows M rows cols unique e.left e.right e.distinct left right


/-- Moving one row across a block of rows has exactly that block's parity. -/
theorem rowDetMove {V : Type} (before middle after : List (V → Gaussian))
    (f : V → Gaussian) (cols : List V) :
    rowDet (before ++ middle ++ f :: after) cols =
      signed gaussianArithmetic (parity middle) (rowDet (before ++ f :: (middle ++ after)) cols) := by
  induction middle generalizing before with
  | nil => simp only [List.append_nil, List.nil_append, parity, signed, Bool.false_eq_true, ↓reduceIte]
  | cons g middle ih =>
    have move := ih (before ++ [g])
    simp only [List.append_assoc, List.cons_append, List.nil_append] at move
    rw [show before ++ (g :: middle) ++ f :: after = before ++ g :: (middle ++ f :: after) by
      simp only [List.append_assoc, List.cons_append]]
    rw [move, rowDetAdjacentSwap before (middle ++ after) f g cols]
    rw [parity, signedNot, signedNeg]
    rfl

theorem parityMap {α β : Type} (xs : List α) (f : α → β) : parity (xs.map f) = parity xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [List.map_cons, parity, ih]

/-- The actual minor sign when a selected row moves past any number of selected
spectators. This is a list-order identity, without adjacency of vertex labels. -/
theorem minorMoveRow {V : Type} (M : V → V → Gaussian) (before middle after : List V)
    (r : V) (cols : List V) :
    minorDet gaussianArithmetic M (before ++ middle ++ r :: after) cols =
      signed gaussianArithmetic (parity middle)
        (minorDet gaussianArithmetic M (before ++ r :: (middle ++ after)) cols) := by
  rw [minorAsRows, minorAsRows]
  simp only [List.map_append, List.map_cons]
  rw [rowDetMove, parityMap]

#print axioms gaussianNegNeg
#print axioms gaussianNegAdd
#print axioms gaussianMulNeg
#print axioms signedNeg
#print axioms gaussianSumNeg
#print axioms doubleChoiceCons
#print axioms doubleChoiceSwap
#print axioms rowDetTwoChoices
#print axioms rowDetHeadSwap
#print axioms rowDetAdjacentSwap
#print axioms rowDetAdjacentDuplicate
#print axioms rowDetPrefixNeg
#print axioms rowDetSeparatedHeadSwap
#print axioms rowDetExchange
#print axioms rowDetDuplicate
#print axioms minorEqualRowsZero
#print axioms minorDuplicateRow
#print axioms twoEndpointAvoids
#print axioms minorExchangeRows
#print axioms gateMinorBoth
#print axioms rowDetMove
#print axioms parityMap
#print axioms minorMoveRow

end MathCommons.ForestMinorAlternation
