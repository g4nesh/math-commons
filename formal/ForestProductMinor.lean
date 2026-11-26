import ForestEndpointMass

/-!
Actual Gaussian-rational matrix multiplication and product minors. This module
proves unrestricted ordered-tuple expansion and full Cauchy–Binet in minor
orders zero, one and two, including arbitrary iteration in order two. It does
not identify that pair evolution with the existing canonical subset histories,
and it does not prove arbitrary-order determinant alternation or Cauchy–Binet.
The single gate squared-minor identity is unconditional apart from unitarity.
One collision simplification explicitly retains its alternation premises.
-/

namespace MathCommons.ForestProductMinor
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass

theorem gaussianAddAssoc (a b c : Gaussian) :
    gaussianArithmetic.add (gaussianArithmetic.add a b) c =
      gaussianArithmetic.add a (gaussianArithmetic.add b c) := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianAddComm (a b : Gaussian) :
    gaussianArithmetic.add a b = gaussianArithmetic.add b a := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianMulAdd (a b c : Gaussian) :
    gaussianArithmetic.mul a (gaussianArithmetic.add b c) =
      gaussianArithmetic.add (gaussianArithmetic.mul a b) (gaussianArithmetic.mul a c) := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianAddMul (a b c : Gaussian) :
    gaussianArithmetic.mul (gaussianArithmetic.add a b) c =
      gaussianArithmetic.add (gaussianArithmetic.mul a c) (gaussianArithmetic.mul b c) := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianMulComm (a b : Gaussian) :
    gaussianArithmetic.mul a b = gaussianArithmetic.mul b a := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianMulAssoc (a b c : Gaussian) :
    gaussianArithmetic.mul (gaussianArithmetic.mul a b) c =
      gaussianArithmetic.mul a (gaussianArithmetic.mul b c) := by
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianSumCongr {α : Type} (xs : List α) (f g : α → Gaussian)
    (same : ∀ x ∈ xs, f x = g x) :
    scalarSum gaussianArithmetic (xs.map f) = scalarSum gaussianArithmetic (xs.map g) := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [List.map_cons, scalarSum]
    rw [same x List.mem_cons_self, ih (fun y hy => same y (List.mem_cons_of_mem x hy))]

theorem gaussianSumAdd {α : Type} (xs : List α) (f g : α → Gaussian) :
    scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.add (f x) (g x)) =
      gaussianArithmetic.add (scalarSum gaussianArithmetic (xs.map f))
        (scalarSum gaussianArithmetic (xs.map g)) := by
  induction xs with
  | nil => exact gaussianArithmetic.addZeroZero.symm
  | cons x xs ih =>
    simp only [List.map_cons, scalarSum, ih]
    apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem gaussianSumScale {α : Type} (xs : List α) (f : α → Gaussian) (a : Gaussian) :
    scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul a (f x)) =
      gaussianArithmetic.mul a (scalarSum gaussianArithmetic (xs.map f)) := by
  induction xs with
  | nil => exact (gaussianArithmetic.mulZero a).symm
  | cons x xs ih => simp only [List.map_cons, scalarSum, ih, gaussianMulAdd]

theorem gaussianSumDelta {V : Type} [DecidableEq V] (xs : List V)
    (nodup : xs.Nodup) (q : V) (member : q ∈ xs) (a : Gaussian) :
    scalarSum gaussianArithmetic (xs.map fun x => if x = q then a else gaussianArithmetic.zero) = a := by
  have one := scalarSumUniqueKey xs (fun x => x) (by simpa using nodup) q member
    (fun x => if x = q then a else gaussianArithmetic.zero)
    (by intro x _ different; simp [different])
  simpa using one

theorem gaussianSumAppend {α : Type} (xs ys : List α) (f : α → Gaussian) :
    scalarSum gaussianArithmetic ((xs ++ ys).map f) =
      gaussianArithmetic.add (scalarSum gaussianArithmetic (xs.map f))
        (scalarSum gaussianArithmetic (ys.map f)) := by
  induction xs with
  | nil => simp only [List.nil_append, List.map_nil, scalarSum, gaussianZeroAdd]
  | cons x xs ih => simp only [List.cons_append, List.map_cons, scalarSum, ih, gaussianAddAssoc]

theorem gaussianSumFlatMap {α β : Type} (xs : List α) (g : α → List β) (f : β → Gaussian) :
    scalarSum gaussianArithmetic ((xs.flatMap g).map f) =
      scalarSum gaussianArithmetic (xs.map fun x => scalarSum gaussianArithmetic ((g x).map f)) := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [List.flatMap_cons, gaussianSumAppend, List.map_cons, scalarSum, ih]

theorem gaussianSumSwap {α β : Type} (xs : List α) (ys : List β) (f : α → β → Gaussian) :
    scalarSum gaussianArithmetic (xs.map fun x => scalarSum gaussianArithmetic (ys.map (f x))) =
      scalarSum gaussianArithmetic (ys.map fun y => scalarSum gaussianArithmetic (xs.map fun x => f x y)) := by
  induction xs with
  | nil =>
    simp only [List.map_nil, scalarSum]
    exact (scalarSumMapZero ys (fun _ => gaussianArithmetic.zero) (by intros; rfl)).symm
  | cons x xs ih =>
    simp only [List.map_cons, scalarSum]
    rw [ih, ← gaussianSumAdd]

theorem gaussianSumScaleRight {α : Type} (xs : List α) (f : α → Gaussian) (a : Gaussian) :
    gaussianArithmetic.mul (scalarSum gaussianArithmetic (xs.map f)) a =
      scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (f x) a) := by
  rw [gaussianMulComm, ← gaussianSumScale]
  apply gaussianSumCongr
  intro x _
  exact gaussianMulComm _ _

section Matrix
variable {V : Type} [DecidableEq V]

/-- Actual finite multiplication; row index first, column index second. -/
def matrixMul (vertices : List V) (M N : V → V → Gaussian) : V → V → Gaussian :=
  fun i j => scalarSum gaussianArithmetic (vertices.map fun k =>
    gaussianArithmetic.mul (M i k) (N k j))

def identityMatrix : V → V → Gaussian :=
  fun i j => if i = j then gaussianArithmetic.one else gaussianArithmetic.zero

/-- Gates act on the left, in chronological order. -/
def matrixProduct (vertices : List V) (matrices : Nat → V → V → Gaussian) :
    Nat → V → V → Gaussian
  | 0 => identityMatrix
  | n+1 => matrixMul vertices (matrices n) (matrixProduct vertices matrices n)

theorem matrixMulIdentity (vertices : List V) (nodup : vertices.Nodup)
    (complete : ∀ x, x ∈ vertices) (M : V → V → Gaussian) :
    matrixMul vertices M identityMatrix = M := by
  funext i j
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k =>
        if k = j then M i j else gaussianArithmetic.zero) := by
      apply gaussianSumCongr
      intro k _
      by_cases equal : k = j
      · subst k; simp [identityMatrix, gaussianMulOne]
      · simp [identityMatrix, equal, gaussianArithmetic.mulZero]
    _ = _ := gaussianSumDelta vertices nodup j (complete j) (M i j)

theorem matrixIdentityMul (vertices : List V) (nodup : vertices.Nodup)
    (complete : ∀ x, x ∈ vertices) (M : V → V → Gaussian) :
    matrixMul vertices identityMatrix M = M := by
  funext i j
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k =>
        if k = i then M i j else gaussianArithmetic.zero) := by
      apply gaussianSumCongr
      intro k _
      by_cases equal : k = i
      · subst k; simp [identityMatrix, gaussianOneMul]
      · simp [identityMatrix, equal, Ne.symm equal, gaussianArithmetic.zeroMul]
    _ = _ := gaussianSumDelta vertices nodup i (complete i) (M i j)

theorem gateProductOne {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) :
    matrixProduct vertices (gateSchedule edges as bs cs ds) 1 =
      gateSchedule edges as bs cs ds 0 := by
  exact matrixMulIdentity vertices (edges 0).nodup (edges 0).complete _

theorem gateMulSpectator {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (i j : V)
    (left : i ≠ e.left) (right : i ≠ e.right) :
    matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M i j = M i j := by
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k =>
        if k = i then M i j else gaussianArithmetic.zero) := by
      apply gaussianSumCongr
      intro k _
      by_cases equal : k = i
      · subst k; simp [twoCoordinateGate, left, right, gaussianOneMul]
      · simp [twoCoordinateGate, left, right, Ne.symm equal, equal, gaussianArithmetic.zeroMul]
    _ = _ := gaussianSumDelta vertices e.nodup i (e.complete i) (M i j)

theorem gateMulLeft {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (j : V) :
    matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M e.left j =
      gaussianArithmetic.add (gaussianArithmetic.mul a (M e.left j))
        (gaussianArithmetic.mul b (M e.right j)) := by
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k => gaussianArithmetic.add
        (if k = e.left then gaussianArithmetic.mul a (M e.left j) else gaussianArithmetic.zero)
        (if k = e.right then gaussianArithmetic.mul b (M e.right j) else gaussianArithmetic.zero)) := by
      apply gaussianSumCongr
      intro k _
      by_cases left : k = e.left
      · subst k; simp [twoCoordinateGate, e.distinct, gaussianAddZero]
      · by_cases right : k = e.right
        · subst k; simp [twoCoordinateGate, e.distinct, Ne.symm e.distinct, gaussianZeroAdd]
        · simp [twoCoordinateGate, left, right, Ne.symm left, gaussianArithmetic.zeroMul,
            gaussianArithmetic.addZeroZero]
    _ = _ := by
      rw [gaussianSumAdd, gaussianSumDelta vertices e.nodup e.left (e.complete e.left),
        gaussianSumDelta vertices e.nodup e.right (e.complete e.right)]

theorem gateMulRight {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (j : V) :
    matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M e.right j =
      gaussianArithmetic.add (gaussianArithmetic.mul c (M e.left j))
        (gaussianArithmetic.mul d (M e.right j)) := by
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k => gaussianArithmetic.add
        (if k = e.left then gaussianArithmetic.mul c (M e.left j) else gaussianArithmetic.zero)
        (if k = e.right then gaussianArithmetic.mul d (M e.right j) else gaussianArithmetic.zero)) := by
      apply gaussianSumCongr
      intro k _
      by_cases left : k = e.left
      · subst k; simp [twoCoordinateGate, e.distinct, Ne.symm e.distinct, gaussianAddZero]
      · by_cases right : k = e.right
        · subst k; simp [twoCoordinateGate, Ne.symm e.distinct, gaussianZeroAdd]
        · simp [twoCoordinateGate, left, right, Ne.symm right, gaussianArithmetic.zeroMul,
            gaussianArithmetic.addZeroZero]
    _ = _ := by
      rw [gaussianSumAdd, gaussianSumDelta vertices e.nodup e.left (e.complete e.left),
        gaussianSumDelta vertices e.nodup e.right (e.complete e.right)]

def replaceRow (M : V → V → Gaussian) (r : V) (f : V → Gaussian) : V → V → Gaussian :=
  fun i j => if i = r then f j else M i j

omit [DecidableEq V] in
theorem minorCongrRows (M N : V → V → Gaussian) (rows cols : List V)
    (same : ∀ r ∈ rows, M r = N r) :
    minorDet gaussianArithmetic M rows cols = minorDet gaussianArithmetic N rows cols := by
  induction rows generalizing cols with
  | nil => rfl
  | cons r rows ih =>
    simp only [minorDet]
    apply gaussianSumCongr
    intro q _
    rw [same r List.mem_cons_self,
      ih q.remaining (fun s hs => same s (List.mem_cons_of_mem r hs))]

theorem minorReplaceAbsent (M : V → V → Gaussian) (rows cols : List V)
    (r : V) (f : V → Gaussian) (absent : r ∉ rows) :
    minorDet gaussianArithmetic (replaceRow M r f) rows cols =
      minorDet gaussianArithmetic M rows cols := by
  apply minorCongrRows
  intro s hs
  funext j
  have different : s ≠ r := by intro equal; subst s; exact absent hs
  simp [replaceRow, different]

theorem signedAdd (odd : Bool) (a b : Gaussian) :
    signed gaussianArithmetic odd (gaussianArithmetic.add a b) =
      gaussianArithmetic.add (signed gaussianArithmetic odd a) (signed gaussianArithmetic odd b) := by
  cases odd <;> simp only [signed, Bool.false_eq_true, ↓reduceIte]
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

theorem signedScale (odd : Bool) (a b : Gaussian) :
    signed gaussianArithmetic odd (gaussianArithmetic.mul a b) =
      gaussianArithmetic.mul a (signed gaussianArithmetic odd b) := by
  cases odd <;> simp only [signed, Bool.false_eq_true, ↓reduceIte]
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

/-- Row additivity at any occupied row, including arbitrary spectator prefixes. -/
theorem minorRowAdd (M : V → V → Gaussian) (r : V) (f g : V → Gaussian)
    (rows cols : List V) (unique : rows.Nodup) (member : r ∈ rows) :
    minorDet gaussianArithmetic (replaceRow M r (fun j => gaussianArithmetic.add (f j) (g j))) rows cols =
      gaussianArithmetic.add (minorDet gaussianArithmetic (replaceRow M r f) rows cols)
        (minorDet gaussianArithmetic (replaceRow M r g) rows cols) := by
  induction rows generalizing cols with
  | nil => simp at member
  | cons s rows ih =>
    obtain ⟨notMember, tailUnique⟩ := List.nodup_cons.mp unique
    simp only [minorDet]
    rw [← gaussianSumAdd]
    apply gaussianSumCongr
    intro q _
    by_cases same : s = r
    · subst s
      simp only [replaceRow, ↓reduceIte]
      rw [minorReplaceAbsent M rows q.remaining r _ notMember,
        minorReplaceAbsent M rows q.remaining r f notMember,
        minorReplaceAbsent M rows q.remaining r g notMember,
        gaussianAddMul, signedAdd]
    · have tailMember : r ∈ rows := (List.mem_cons.mp member).resolve_left (Ne.symm same)
      simp only [replaceRow, if_neg same]
      rw [ih q.remaining tailUnique tailMember, gaussianMulAdd, signedAdd]

theorem minorRowScale (M : V → V → Gaussian) (r : V) (f : V → Gaussian) (a : Gaussian)
    (rows cols : List V) (unique : rows.Nodup) (member : r ∈ rows) :
    minorDet gaussianArithmetic (replaceRow M r (fun j => gaussianArithmetic.mul a (f j))) rows cols =
      gaussianArithmetic.mul a (minorDet gaussianArithmetic (replaceRow M r f) rows cols) := by
  induction rows generalizing cols with
  | nil => simp at member
  | cons s rows ih =>
    obtain ⟨notMember, tailUnique⟩ := List.nodup_cons.mp unique
    simp only [minorDet]
    rw [← gaussianSumScale]
    apply gaussianSumCongr
    intro q _
    by_cases same : s = r
    · subst s
      simp only [replaceRow, ↓reduceIte]
      rw [minorReplaceAbsent M rows q.remaining r _ notMember,
        minorReplaceAbsent M rows q.remaining r f notMember,
        gaussianMulAssoc, signedScale]
    · have tailMember : r ∈ rows := (List.mem_cons.mp member).resolve_left (Ne.symm same)
      simp only [replaceRow, if_neg same]
      rw [ih q.remaining tailUnique tailMember]
      rw [← gaussianMulAssoc, gaussianMulComm (M s q.column) a, gaussianMulAssoc, signedScale]

theorem minorRowLinearCombination (M : V → V → Gaussian) (r : V) (f g : V → Gaussian)
    (a b : Gaussian) (rows cols : List V) (unique : rows.Nodup) (member : r ∈ rows) :
    minorDet gaussianArithmetic (replaceRow M r (fun j => gaussianArithmetic.add
      (gaussianArithmetic.mul a (f j)) (gaussianArithmetic.mul b (g j)))) rows cols =
    gaussianArithmetic.add
      (gaussianArithmetic.mul a (minorDet gaussianArithmetic (replaceRow M r f) rows cols))
      (gaussianArithmetic.mul b (minorDet gaussianArithmetic (replaceRow M r g) rows cols)) := by
  rw [minorRowAdd M r _ _ rows cols unique member,
    minorRowScale M r f a rows cols unique member, minorRowScale M r g b rows cols unique member]

theorem replaceOwnRow (M : V → V → Gaussian) (r : V) : replaceRow M r (M r) = M := by
  funext i j
  by_cases same : i = r
  · subst i; simp [replaceRow]
  · simp [replaceRow, same]

/-- If neither gate endpoint is an occupied row, its minor is unchanged. -/
theorem gateMinorNeither {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (rows cols : List V)
    (left : e.left ∉ rows) (right : e.right ∉ rows) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) rows cols =
      minorDet gaussianArithmetic M rows cols := by
  apply minorCongrRows
  intro r hr
  funext j
  apply gateMulSpectator
  · intro equal; subst r; exact left hr
  · intro equal; subst r; exact right hr

/-- Actual single-occupied left-row update. The replacement row remains in its
original list position: converting to sorted subset order is a separate sign lemma. -/
theorem gateMinorLeftOnly {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (rows cols : List V)
    (unique : rows.Nodup) (left : e.left ∈ rows) (right : e.right ∉ rows) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) rows cols =
      gaussianArithmetic.add (gaussianArithmetic.mul a (minorDet gaussianArithmetic M rows cols))
        (gaussianArithmetic.mul b
          (minorDet gaussianArithmetic (replaceRow M e.left (M e.right)) rows cols)) := by
  calc
    _ = minorDet gaussianArithmetic (replaceRow M e.left (fun j => gaussianArithmetic.add
        (gaussianArithmetic.mul a (M e.left j)) (gaussianArithmetic.mul b (M e.right j)))) rows cols := by
      apply minorCongrRows
      intro r hr
      funext j
      by_cases same : r = e.left
      · subst r; simpa only [replaceRow, ↓reduceIte] using gateMulLeft e a b c d M j
      · have other : r ≠ e.right := by intro equal; subst r; exact right hr
        simpa only [replaceRow, if_neg same] using gateMulSpectator e a b c d M r j same other
    _ = _ := by
      rw [minorRowLinearCombination M e.left (M e.left) (M e.right) a b rows cols unique left,
        replaceOwnRow]

theorem gateMinorRightOnly {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (rows cols : List V)
    (unique : rows.Nodup) (left : e.left ∉ rows) (right : e.right ∈ rows) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) rows cols =
      gaussianArithmetic.add
        (gaussianArithmetic.mul c
          (minorDet gaussianArithmetic (replaceRow M e.right (M e.left)) rows cols))
        (gaussianArithmetic.mul d (minorDet gaussianArithmetic M rows cols)) := by
  calc
    _ = minorDet gaussianArithmetic (replaceRow M e.right (fun j => gaussianArithmetic.add
        (gaussianArithmetic.mul c (M e.left j)) (gaussianArithmetic.mul d (M e.right j)))) rows cols := by
      apply minorCongrRows
      intro r hr
      funext j
      by_cases same : r = e.right
      · subst r; simpa only [replaceRow, ↓reduceIte] using gateMulRight e a b c d M j
      · have other : r ≠ e.left := by intro equal; subst r; exact left hr
        simpa only [replaceRow, if_neg same] using gateMulSpectator e a b c d M r j other same
    _ = _ := by
      rw [minorRowLinearCombination M e.right (M e.left) (M e.right) c d rows cols unique right,
        replaceOwnRow]

theorem productEmptyMinor (vertices : List V) (matrices : Nat → V → V → Gaussian) (n : Nat) :
    minorDet gaussianArithmetic (matrixProduct vertices matrices n) [] [] = gaussianArithmetic.one := rfl

theorem productSingletonMinor (vertices : List V) (matrices : Nat → V → V → Gaussian)
    (n : Nat) (r c : V) :
    minorDet gaussianArithmetic (matrixProduct vertices matrices n) [r] [c] =
      matrixProduct vertices matrices n r c := by
  exact singletonMinor _ _ _

/-- The existing gate probability theorem now applies to the actual matrix
product of length one, in every occupation sector and endpoint order. -/
theorem oneGateProductMagnitude {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (unitary : UnitaryBlock (as 0) (bs 0) (cs 0) (ds 0))
    (start finish : Subset V) :
    normSq (minorDet gaussianArithmetic
      (matrixProduct vertices (gateSchedule edges as bs cs ds) 1)
      (vertices.filter finish) (vertices.filter start)) = stepProb (edges 0) (bs 0) start finish := by
  rw [gateProductOne]
  exact stepMagnitudeGeneral (edges 0) (as 0) (bs 0) (cs 0) (ds 0) unitary start finish

theorem replaceRowsCommute (M : V → V → Gaussian) (r s : V) (different : r ≠ s)
    (f g : V → Gaussian) :
    replaceRow (replaceRow M r f) s g = replaceRow (replaceRow M s g) r f := by
  funext i j
  by_cases left : i = r
  · subst i; simp [replaceRow, different]
  · by_cases right : i = s
    · subst i; simp [replaceRow, Ne.symm different]
    · simp [replaceRow, left, right]

theorem replaceOtherOwnRow (M : V → V → Gaussian) (r s : V) (different : r ≠ s)
    (f : V → Gaussian) : replaceRow (replaceRow M r f) s (M s) = replaceRow M r f := by
  rw [replaceRowsCommute M r s different, replaceOwnRow]

theorem gateAsRowUpdates {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) :
    matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M =
      replaceRow (replaceRow M e.right (fun j => gaussianArithmetic.add
        (gaussianArithmetic.mul c (M e.left j)) (gaussianArithmetic.mul d (M e.right j))))
        e.left (fun j => gaussianArithmetic.add
          (gaussianArithmetic.mul a (M e.left j)) (gaussianArithmetic.mul b (M e.right j))) := by
  funext i j
  by_cases left : i = e.left
  · subst i; simpa only [replaceRow, ↓reduceIte] using gateMulLeft e a b c d M j
  · by_cases right : i = e.right
    · subst i
      simpa only [replaceRow, if_neg (Ne.symm e.distinct), ↓reduceIte] using gateMulRight e a b c d M j
    · simpa only [replaceRow, if_neg left, if_neg right] using gateMulSpectator e a b c d M i j left right

/-- Unconditional four-term expansion in the both-occupied sector. Duplicate
row terms and the exchanged-row term are displayed rather than silently canceled. -/
theorem gateMinorBothExpansion {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (rows cols : List V)
    (unique : rows.Nodup) (left : e.left ∈ rows) (right : e.right ∈ rows) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) rows cols =
      gaussianArithmetic.add
        (gaussianArithmetic.mul a (gaussianArithmetic.add
          (gaussianArithmetic.mul c
            (minorDet gaussianArithmetic (replaceRow M e.right (M e.left)) rows cols))
          (gaussianArithmetic.mul d (minorDet gaussianArithmetic M rows cols))))
        (gaussianArithmetic.mul b (gaussianArithmetic.add
          (gaussianArithmetic.mul c
            (minorDet gaussianArithmetic
              (replaceRow (replaceRow M e.left (M e.right)) e.right (M e.left)) rows cols))
          (gaussianArithmetic.mul d
            (minorDet gaussianArithmetic (replaceRow M e.left (M e.right)) rows cols)))) := by
  rw [gateAsRowUpdates]
  rw [minorRowLinearCombination _ e.left (M e.left) (M e.right) a b rows cols unique left]
  rw [replaceRowsCommute M e.right e.left (Ne.symm e.distinct), replaceOwnRow]
  rw [minorRowLinearCombination M e.right (M e.left) (M e.right) c d rows cols unique right,
    replaceOwnRow]
  rw [replaceRowsCommute M e.right e.left (Ne.symm e.distinct)]
  rw [minorRowLinearCombination _ e.right (M e.left) (M e.right) c d rows cols unique right,
    replaceOtherOwnRow M e.left e.right e.distinct]

/-- A precise reduction of the collision sector to three ordinary determinant
alternation facts. These premises have NOT been proved by this module. -/
theorem gateMinorBothOfAlternation {vertices : List V} (e : Edge vertices)
    (a b c d : Gaussian) (M : V → V → Gaussian) (rows cols : List V)
    (unique : rows.Nodup) (left : e.left ∈ rows) (right : e.right ∈ rows)
    (duplicateLeft : minorDet gaussianArithmetic (replaceRow M e.right (M e.left)) rows cols = gaussianArithmetic.zero)
    (duplicateRight : minorDet gaussianArithmetic (replaceRow M e.left (M e.right)) rows cols = gaussianArithmetic.zero)
    (exchanged : minorDet gaussianArithmetic
      (replaceRow (replaceRow M e.left (M e.right)) e.right (M e.left)) rows cols =
        gaussianArithmetic.neg (minorDet gaussianArithmetic M rows cols)) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) rows cols =
      gaussianArithmetic.mul (blockDet a b c d) (minorDet gaussianArithmetic M rows cols) := by
  rw [gateMinorBothExpansion e a b c d M rows cols unique left right,
    duplicateLeft, duplicateRight, exchanged]
  apply Prod.ext <;> simp only [gaussianArithmetic, blockDet] <;> grind

end Matrix

theorem signedSum {α : Type} (xs : List α) (f : α → Gaussian) (odd : Bool) :
    signed gaussianArithmetic odd (scalarSum gaussianArithmetic (xs.map f)) =
      scalarSum gaussianArithmetic (xs.map fun x => signed gaussianArithmetic odd (f x)) := by
  induction xs with
  | nil => cases odd <;> rfl
  | cons x xs ih => simp only [List.map_cons, scalarSum, signedAdd, ih]

theorem laplaceSum {α : Type} (xs : List α) (weights values : α → Gaussian)
    (odd : Bool) (d : Gaussian) :
    signed gaussianArithmetic odd (gaussianArithmetic.mul
      (scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x) (values x))) d) =
      scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x)
        (signed gaussianArithmetic odd (gaussianArithmetic.mul (values x) d))) := by
  rw [gaussianSumScaleRight, signedSum]
  apply gaussianSumCongr
  intro x _
  rw [gaussianMulAssoc, signedScale]

theorem laplaceTailSum {α : Type} (xs : List α) (weights values : α → Gaussian)
    (odd : Bool) (d : Gaussian) :
    signed gaussianArithmetic odd (gaussianArithmetic.mul d
      (scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x) (values x)))) =
      scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x)
        (signed gaussianArithmetic odd (gaussianArithmetic.mul d (values x)))) := by
  rw [gaussianMulComm d, laplaceSum]
  apply gaussianSumCongr
  intro x _
  rw [gaussianMulComm (values x) d]

/-- The same signed row expansion with row vectors made explicit. -/
def rowDet {V : Type} : List (V → Gaussian) → List V → Gaussian
  | [], cols => if cols.isEmpty then gaussianArithmetic.one else gaussianArithmetic.zero
  | f :: rows, cols => scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
      signed gaussianArithmetic q.odd (gaussianArithmetic.mul (f q.column) (rowDet rows q.remaining)))

theorem minorAsRows {V : Type} (M : V → V → Gaussian) (rows cols : List V) :
    minorDet gaussianArithmetic M rows cols = rowDet (rows.map M) cols := by
  induction rows generalizing cols with
  | nil => rfl
  | cons r rows ih =>
    simp only [minorDet, List.map_cons, rowDet]
    apply gaussianSumCongr
    intro q _
    rw [ih]

theorem rowDetHeadSum {α V : Type} (xs : List α) (weights : α → Gaussian)
    (vectors : α → V → Gaussian) (rows : List (V → Gaussian)) (cols : List V) :
    rowDet ((fun j => scalarSum gaussianArithmetic
      (xs.map fun x => gaussianArithmetic.mul (weights x) (vectors x j))) :: rows) cols =
      scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x)
        (rowDet (vectors x :: rows) cols)) := by
  simp only [rowDet]
  calc
    _ = scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
        scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x)
          (signed gaussianArithmetic q.odd (gaussianArithmetic.mul (vectors x q.column)
            (rowDet rows q.remaining))))) := by
      apply gaussianSumCongr
      intro q _
      exact laplaceSum xs weights (fun x => vectors x q.column) q.odd (rowDet rows q.remaining)
    _ = _ := by
      rw [gaussianSumSwap]
      apply gaussianSumCongr
      intro x _
      exact gaussianSumScale _ _ _

theorem rowDetTailSum {α V : Type} (xs : List α) (weights : α → Gaussian)
    (tails : α → List (V → Gaussian)) (rows : List (V → Gaussian))
    (decomposition : ∀ cols, rowDet rows cols =
      scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x) (rowDet (tails x) cols)))
    (head : V → Gaussian) (cols : List V) :
    rowDet (head :: rows) cols = scalarSum gaussianArithmetic
      (xs.map fun x => gaussianArithmetic.mul (weights x) (rowDet (head :: tails x) cols)) := by
  simp only [rowDet]
  calc
    _ = scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
        scalarSum gaussianArithmetic (xs.map fun x => gaussianArithmetic.mul (weights x)
          (signed gaussianArithmetic q.odd (gaussianArithmetic.mul (head q.column)
            (rowDet (tails x) q.remaining))))) := by
      apply gaussianSumCongr
      intro q _
      rw [decomposition]
      exact laplaceTailSum xs weights (fun x => rowDet (tails x) q.remaining) q.odd (head q.column)
    _ = _ := by
      rw [gaussianSumSwap]
      apply gaussianSumCongr
      intro x _
      exact gaussianSumScale _ _ _

/-- All intermediate ordered row tuples, without prematurely deleting repeats.
The first component is the product of the selected matrix entries. -/
def rowExpansion {V : Type} (vertices : List V) (M : V → V → Gaussian) :
    List V → List (Gaussian × List V)
  | [] => [(gaussianArithmetic.one, [])]
  | r :: rows => vertices.flatMap fun k =>
      (rowExpansion vertices M rows).map fun term =>
        (gaussianArithmetic.mul (M r k) term.1, k :: term.2)

/-- A full matrix-product minor expansion over ordered intermediate row tuples.
This is before the alternating cancellation that turns tuples into subsets. -/
theorem productMinorTupleExpansion {V : Type} (vertices : List V)
    (M N : V → V → Gaussian) (rows cols : List V) :
    minorDet gaussianArithmetic (matrixMul vertices M N) rows cols =
      scalarSum gaussianArithmetic ((rowExpansion vertices M rows).map fun term =>
        gaussianArithmetic.mul term.1 (minorDet gaussianArithmetic N term.2 cols)) := by
  rw [minorAsRows]
  induction rows generalizing cols with
  | nil => simp only [rowExpansion, List.map_nil, List.map_cons, scalarSum,
      gaussianOneMul, gaussianAddZero, minorDet, rowDet]
  | cons r rows ih =>
    have tail : ∀ columns, rowDet (rows.map (matrixMul vertices M N)) columns =
        scalarSum gaussianArithmetic ((rowExpansion vertices M rows).map fun term =>
          gaussianArithmetic.mul term.1 (rowDet (term.2.map N) columns)) := by
      intro columns
      rw [ih]
      apply gaussianSumCongr
      intro term _
      rw [minorAsRows]
    simp only [List.map_cons]
    change rowDet ((fun j => scalarSum gaussianArithmetic (vertices.map fun k =>
      gaussianArithmetic.mul (M r k) (N k j))) :: rows.map (matrixMul vertices M N)) cols = _
    rw [rowDetHeadSum vertices (M r) N]
    simp only [rowExpansion, gaussianSumFlatMap, List.map_map, minorAsRows]
    apply gaussianSumCongr
    intro k _
    rw [rowDetTailSum (rowExpansion vertices M rows) Prod.fst (fun term => term.2.map N) _ tail]
    rw [← gaussianSumScale]
    apply gaussianSumCongr
    intro term _
    exact (gaussianMulAssoc _ _ _).symm

def orderedPairs {V : Type} : List V → List (V × V)
  | [] => []
  | x :: xs => (xs.map fun y => (x, y)) ++ orderedPairs xs

theorem gaussianSumPairs {V : Type} (xs : List V) (f : V → V → Gaussian)
    (diagonal : ∀ x ∈ xs, f x x = gaussianArithmetic.zero) :
    scalarSum gaussianArithmetic (xs.map fun x => scalarSum gaussianArithmetic (xs.map (f x))) =
      scalarSum gaussianArithmetic ((orderedPairs xs).map fun p =>
        gaussianArithmetic.add (f p.1 p.2) (f p.2 p.1)) := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [List.map_cons, scalarSum, orderedPairs, gaussianSumAppend, List.map_map]
    rw [diagonal x List.mem_cons_self, gaussianZeroAdd]
    rw [gaussianSumAdd]
    rw [ih (fun y hy => diagonal y (List.mem_cons_of_mem x hy))]
    rw [← gaussianAddAssoc, ← gaussianSumAdd]
    rfl

theorem productMinorTwoRows {V : Type} (vertices : List V)
    (M N : V → V → Gaussian) (r s : V) (cols : List V) :
    minorDet gaussianArithmetic (matrixMul vertices M N) [r, s] cols =
      scalarSum gaussianArithmetic (vertices.map fun k =>
        scalarSum gaussianArithmetic (vertices.map fun l =>
          gaussianArithmetic.mul (gaussianArithmetic.mul (M r k) (M s l))
            (minorDet gaussianArithmetic N [k, l] cols))) := by
  rw [productMinorTupleExpansion]
  simp [rowExpansion, gaussianSumFlatMap, List.map_map, scalarSum, gaussianMulOne,
    gaussianAddZero]

/-- Genuine order-two Cauchy–Binet for arbitrary finite matrices. Each pair of
enumeration positions is listed once; a duplicate-free enumeration therefore
lists each unordered pair of distinct intermediate vertices once. -/
theorem cauchyBinetTwo {V : Type} (vertices : List V)
    (M N : V → V → Gaussian) (r s c d : V) :
    minorDet gaussianArithmetic (matrixMul vertices M N) [r, s] [c, d] =
      scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p =>
        gaussianArithmetic.mul (minorDet gaussianArithmetic M [r, s] [p.1, p.2])
          (minorDet gaussianArithmetic N [p.1, p.2] [c, d])) := by
  rw [productMinorTwoRows]
  rw [gaussianSumPairs vertices _ (by
    intro k _
    rw [gaussianTwoByTwo]
    apply Prod.ext <;> simp only [gaussianArithmetic, blockDet] <;> grind)]
  apply gaussianSumCongr
  intro p _
  simp only [gaussianTwoByTwo]
  apply Prod.ext <;> simp only [gaussianArithmetic, blockDet] <;> grind

/-- Entries of the order-two compound kernel, with explicit ordered pairs. -/
def twoMinorKernel {V : Type} (M : V → V → Gaussian) : (V × V) → (V × V) → Gaussian :=
  fun p q => minorDet gaussianArithmetic M [p.1, p.2] [q.1, q.2]

/-- Iteration over canonical intermediate pairs. The initial kernel is the
actual identity-matrix minor, retaining signs for arbitrary queried row orders. -/
def twoMinorEvolution {V : Type} [DecidableEq V] (vertices : List V)
    (matrices : Nat → V → V → Gaussian) : Nat → (V × V) → (V × V) → Gaussian
  | 0 => twoMinorKernel identityMatrix
  | n+1 => matrixMul (orderedPairs vertices) (twoMinorKernel (matrices n))
      (twoMinorEvolution vertices matrices n)

/-- Actual iterated Cauchy–Binet in order two, for any schedule length and
ambient dimension. This has no unitary, forest, or determinant-identity premise. -/
theorem iteratedCauchyBinetTwo {V : Type} [DecidableEq V] (vertices : List V)
    (matrices : Nat → V → V → Gaussian) (n : Nat) (r s c d : V) :
    minorDet gaussianArithmetic (matrixProduct vertices matrices n) [r, s] [c, d] =
      twoMinorEvolution vertices matrices n (r, s) (c, d) := by
  induction n generalizing r s c d with
  | zero => rfl
  | succ n ih =>
    simp only [matrixProduct]
    rw [cauchyBinetTwo]
    unfold twoMinorEvolution matrixMul
    apply gaussianSumCongr
    intro p _
    rw [ih]
    rfl

theorem gateProductTwoMinorRecurrence {V : Type} [DecidableEq V] {vertices : List V}
    (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian) (n : Nat) (r s c d : V) :
    minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) (n+1))
      [r, s] [c, d] =
    scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p => gaussianArithmetic.mul
      (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n) [r, s] [p.1, p.2])
      (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
        [p.1, p.2] [c, d])) := by
  exact cauchyBinetTwo vertices _ _ r s c d

theorem cauchyBinetOne {V : Type} (vertices : List V) (M N : V → V → Gaussian) (r c : V) :
    minorDet gaussianArithmetic (matrixMul vertices M N) [r] [c] =
      scalarSum gaussianArithmetic (vertices.map fun k => gaussianArithmetic.mul
        (minorDet gaussianArithmetic M [r] [k]) (minorDet gaussianArithmetic N [k] [c])) := by
  simp only [singletonMinor, matrixMul]

end MathCommons.ForestProductMinor

#print axioms MathCommons.ForestProductMinor.gaussianAddAssoc
#print axioms MathCommons.ForestProductMinor.gaussianAddComm
#print axioms MathCommons.ForestProductMinor.gaussianMulAdd
#print axioms MathCommons.ForestProductMinor.gaussianAddMul
#print axioms MathCommons.ForestProductMinor.gaussianMulComm
#print axioms MathCommons.ForestProductMinor.gaussianMulAssoc
#print axioms MathCommons.ForestProductMinor.gaussianSumCongr
#print axioms MathCommons.ForestProductMinor.gaussianSumAdd
#print axioms MathCommons.ForestProductMinor.gaussianSumScale
#print axioms MathCommons.ForestProductMinor.gaussianSumDelta
#print axioms MathCommons.ForestProductMinor.gaussianSumAppend
#print axioms MathCommons.ForestProductMinor.gaussianSumFlatMap
#print axioms MathCommons.ForestProductMinor.gaussianSumSwap
#print axioms MathCommons.ForestProductMinor.gaussianSumScaleRight
#print axioms MathCommons.ForestProductMinor.matrixMulIdentity
#print axioms MathCommons.ForestProductMinor.matrixIdentityMul
#print axioms MathCommons.ForestProductMinor.gateProductOne
#print axioms MathCommons.ForestProductMinor.gateMulSpectator
#print axioms MathCommons.ForestProductMinor.gateMulLeft
#print axioms MathCommons.ForestProductMinor.gateMulRight
#print axioms MathCommons.ForestProductMinor.minorCongrRows
#print axioms MathCommons.ForestProductMinor.minorReplaceAbsent
#print axioms MathCommons.ForestProductMinor.signedAdd
#print axioms MathCommons.ForestProductMinor.signedScale
#print axioms MathCommons.ForestProductMinor.minorRowAdd
#print axioms MathCommons.ForestProductMinor.minorRowScale
#print axioms MathCommons.ForestProductMinor.minorRowLinearCombination
#print axioms MathCommons.ForestProductMinor.replaceOwnRow
#print axioms MathCommons.ForestProductMinor.gateMinorNeither
#print axioms MathCommons.ForestProductMinor.gateMinorLeftOnly
#print axioms MathCommons.ForestProductMinor.gateMinorRightOnly
#print axioms MathCommons.ForestProductMinor.productEmptyMinor
#print axioms MathCommons.ForestProductMinor.productSingletonMinor
#print axioms MathCommons.ForestProductMinor.oneGateProductMagnitude
#print axioms MathCommons.ForestProductMinor.replaceRowsCommute
#print axioms MathCommons.ForestProductMinor.replaceOtherOwnRow
#print axioms MathCommons.ForestProductMinor.gateAsRowUpdates
#print axioms MathCommons.ForestProductMinor.gateMinorBothExpansion
#print axioms MathCommons.ForestProductMinor.gateMinorBothOfAlternation
#print axioms MathCommons.ForestProductMinor.signedSum
#print axioms MathCommons.ForestProductMinor.laplaceSum
#print axioms MathCommons.ForestProductMinor.laplaceTailSum
#print axioms MathCommons.ForestProductMinor.minorAsRows
#print axioms MathCommons.ForestProductMinor.rowDetHeadSum
#print axioms MathCommons.ForestProductMinor.rowDetTailSum
#print axioms MathCommons.ForestProductMinor.productMinorTupleExpansion
#print axioms MathCommons.ForestProductMinor.gaussianSumPairs
#print axioms MathCommons.ForestProductMinor.productMinorTwoRows
#print axioms MathCommons.ForestProductMinor.cauchyBinetTwo
#print axioms MathCommons.ForestProductMinor.iteratedCauchyBinetTwo
#print axioms MathCommons.ForestProductMinor.gateProductTwoMinorRecurrence
#print axioms MathCommons.ForestProductMinor.cauchyBinetOne
