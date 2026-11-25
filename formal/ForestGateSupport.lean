import ForestTrajectory

/-!
# Support of a single two-coordinate gate minor

The determinant below is a concrete alternating Laplace expansion, not a
postulated support oracle. The scalar interface contains only elementary zero
laws needed to infer nonzero factors from a nonzero product; it does not assume
that products of nonzero factors remain nonzero. A concrete rational instance
is provided. Matching to Mathlib's complex determinant, the unitary transition
magnitudes, and Cauchy–Binet composition remain separate obligations.
-/
namespace MathCommons.ForestGateSupport

open ForestTrajectory

/-- Precisely the scalar operations and zero identities used in the support proof. -/
structure ZeroArithmetic (A : Type) where
  zero : A
  one : A
  add : A → A → A
  mul : A → A → A
  neg : A → A
  addZeroZero : add zero zero = zero
  zeroMul : ∀ x, mul zero x = zero
  mulZero : ∀ x, mul x zero = zero
  negZero : neg zero = zero

def rationalArithmetic : ZeroArithmetic Rat where
  zero := 0
  one := 1
  add := (· + ·)
  mul := (· * ·)
  neg := Neg.neg
  addZeroZero := Rat.zero_add 0
  zeroMul := Rat.zero_mul
  mulZero := Rat.mul_zero
  negZero := by decide

def scalarSum {A : Type} (R : ZeroArithmetic A) : List A → A
  | [] => R.zero
  | x :: xs => R.add x (scalarSum R xs)

def signed {A : Type} (R : ZeroArithmetic A) (odd : Bool) (x : A) : A :=
  if odd then R.neg x else x

/-- Selecting one column, keeping the others in order, and recording its parity. -/
structure ColumnChoice (V : Type) where
  column : V
  remaining : List V
  odd : Bool

def columnChoices {V : Type} : List V → List (ColumnChoice V)
  | [] => []
  | c :: cs => ⟨c, cs, false⟩ ::
      (columnChoices cs).map (fun q => ⟨q.column, c :: q.remaining, !q.odd⟩)

theorem choicePerm {V : Type} (cols : List V) (q : ColumnChoice V)
    (member : q ∈ columnChoices cols) : (q.column :: q.remaining).Perm cols := by
  induction cols generalizing q with
  | nil => simp [columnChoices] at member
  | cons c cs ih =>
    simp only [columnChoices, List.mem_cons, List.mem_map] at member
    rcases member with first | ⟨p, hp, rfl⟩
    · subst q; exact .rfl
    · exact (List.Perm.swap c p.column p.remaining).trans ((ih p hp).cons c)

/-- Concrete signed row expansion. Rectangular inputs evaluate to zero;
the empty square minor evaluates to the scalar one. -/
def minorDet {A V : Type} (R : ZeroArithmetic A) (M : V → V → A) : List V → List V → A
  | [], cols => if cols.isEmpty then R.one else R.zero
  | r :: rs, cols => scalarSum R ((columnChoices cols).map
      (fun q => signed R q.odd (R.mul (M r q.column) (minorDet R M rs q.remaining))))

theorem nonzeroSumTerm {A : Type} (R : ZeroArithmetic A) (xs : List A)
    (nonzero : scalarSum R xs ≠ R.zero) : ∃ x ∈ xs, x ≠ R.zero := by
  induction xs with
  | nil => exact False.elim (nonzero rfl)
  | cons x xs ih =>
    by_cases hx : x = R.zero
    · have tail : scalarSum R xs ≠ R.zero := by
        intro h
        apply nonzero
        simp only [scalarSum, hx, h, R.addZeroZero]
      obtain ⟨y, hy, hn⟩ := ih tail
      exact ⟨y, List.mem_cons_of_mem x hy, hn⟩
    · exact ⟨x, List.mem_cons_self, hx⟩

theorem nonzeroSignedFactors {A : Type} (R : ZeroArithmetic A) (odd : Bool) (x y : A)
    (nonzero : signed R odd (R.mul x y) ≠ R.zero) : x ≠ R.zero ∧ y ≠ R.zero := by
  have product : R.mul x y ≠ R.zero := by
    intro h
    apply nonzero
    cases odd <;> simp [signed, h, R.negZero]
  constructor
  · intro h; apply product; rw [h, R.zeroMul]
  · intro h; apply product; rw [h, R.mulZero]

/-- A recursively witnessed bijective row/column matching with nonzero entries. -/
inductive Matching {A V : Type} (R : ZeroArithmetic A) (M : V → V → A) :
    List V → List V → Prop where
  | nil : Matching R M [] []
  | cons {r : V} {rs : List V} {c : V} {rest cols : List V}
      (permutation : (c :: rest).Perm cols) (entry : M r c ≠ R.zero)
      (tail : Matching R M rs rest) : Matching R M (r :: rs) cols

theorem nonzeroMinorMatching {A V : Type} (R : ZeroArithmetic A) (M : V → V → A)
    (rows cols : List V) (nonzero : minorDet R M rows cols ≠ R.zero) :
    Matching R M rows cols := by
  induction rows generalizing cols with
  | nil =>
    cases cols with
    | nil => exact .nil
    | cons c cs => exact False.elim (nonzero rfl)
  | cons r rs ih =>
    obtain ⟨value, member, hn⟩ := nonzeroSumTerm R _ nonzero
    obtain ⟨q, hq, rfl⟩ := List.mem_map.mp member
    have factors := nonzeroSignedFactors R q.odd (M r q.column) _ hn
    exact .cons (choicePerm cols q hq) factors.1 (ih q.remaining factors.2)

theorem matchingLength {A V : Type} {R : ZeroArithmetic A} {M : V → V → A}
    {rows cols : List V} (matching : Matching R M rows cols) : rows.length = cols.length := by
  induction matching with
  | nil => rfl
  | cons permutation _ _ ih =>
    have length := permutation.length_eq
    simp only [List.length_cons] at *
    omega

theorem matchingRowSupport {A V : Type} {R : ZeroArithmetic A} {M : V → V → A}
    {rows cols : List V} (matching : Matching R M rows cols) :
    ∀ x ∈ rows, ∃ y ∈ cols, M x y ≠ R.zero := by
  induction matching with
  | nil => simp
  | @cons r rs c rest cols permutation entry tail ih =>
    intro x hx
    rcases List.mem_cons.mp hx with rfl | hx
    · exact ⟨c, permutation.mem_iff.mp List.mem_cons_self, entry⟩
    · obtain ⟨y, hy, hxy⟩ := ih x hx
      exact ⟨y, permutation.mem_iff.mp (List.mem_cons_of_mem c hy), hxy⟩

theorem matchingColumnSupport {A V : Type} {R : ZeroArithmetic A} {M : V → V → A}
    {rows cols : List V} (matching : Matching R M rows cols) :
    ∀ y ∈ cols, ∃ x ∈ rows, M x y ≠ R.zero := by
  induction matching with
  | nil => simp
  | @cons r rs c rest cols permutation entry tail ih =>
    intro y hy
    rcases List.mem_cons.mp (permutation.mem_iff.mpr hy) with rfl | hy
    · exact ⟨r, List.mem_cons_self, entry⟩
    · obtain ⟨x, hx, hxy⟩ := ih y hy
      exact ⟨x, List.mem_cons_of_mem r hx, hxy⟩

section Gate
variable {A V : Type} [DecidableEq V] {vertices : List V}

/-- Matrix support of a two-coordinate gate: every other off-diagonal entry is zero. -/
def GateSupport (R : ZeroArithmetic A) (e : Edge vertices) (M : V → V → A) : Prop :=
  ∀ x y, x ≠ y → ((x ≠ e.left ∧ x ≠ e.right) ∨ (y ≠ e.left ∧ y ≠ e.right)) →
    M x y = R.zero

/-- An arbitrary two-by-two block, with the identity on every other coordinate. -/
def twoCoordinateGate (R : ZeroArithmetic A) (e : Edge vertices) (a b c d : A) : V → V → A :=
  fun x y => if x = y then
      if x = e.left then a else if x = e.right then d else R.one
    else if x = e.left ∧ y = e.right then b
    else if x = e.right ∧ y = e.left then c else R.zero

theorem twoCoordinateGateSupport (R : ZeroArithmetic A) (e : Edge vertices) (a b c d : A) :
    GateSupport R e (twoCoordinateGate R e a b c d) := by
  intro x y different forbidden
  rcases forbidden with ⟨xl, xr⟩ | ⟨yl, yr⟩
  · simp [twoCoordinateGate, different, xl, xr]
  · simp [twoCoordinateGate, different, yl, yr]

/-- Equal cardinality and agreement away from the edge force stay or swap. -/
theorem outsideAgreementStep (e : Edge vertices) (s t : Subset V)
    (outside : ∀ x, x ≠ e.left → x ≠ e.right → t x = s x)
    (size : cutCount vertices (fun _ => true) t = cutCount vertices (fun _ => true) s) :
    Step e s t := by
  have rest : cutCount e.rest (fun _ => true) t = cutCount e.rest (fun _ => true) s := by
    apply List.countP_congr
    intro x hx
    simp only [Bool.true_and]
    rw [outside x (e.avoids x hx).1 (e.avoids x hx).2]
  rw [countSplit e, countSplit e, rest] at size
  simp only [Bool.true_and] at size
  by_cases sameLeft : t e.left = s e.left
  · have sameRight : t e.right = s e.right := by
      clear outside rest
      cases sl : s e.left <;> cases sr : s e.right <;>
        cases tl : t e.left <;> cases tr : t e.right <;> simp_all <;> omega
    left
    funext x
    by_cases xl : x = e.left
    · subst x; exact sameLeft
    · by_cases xr : x = e.right
      · subst x; exact sameRight
      · exact outside x xl xr
  · have exchanged : t e.left = s e.right ∧ t e.right = s e.left := by
      clear outside rest
      cases sl : s e.left <;> cases sr : s e.right <;>
        cases tl : t e.left <;> cases tr : t e.right <;> simp_all <;> omega
    right
    funext x
    by_cases xl : x = e.left
    · subst x; simpa only [swapAtLeft] using exchanged.1
    · by_cases xr : x = e.right
      · subst x; simpa only [swapAtRight] using exchanged.2
      · simp only [swapSubset, xl, xr, ↓reduceIte]
        exact outside x xl xr

/-- Nonzero Laplace minor implies a legal subset transition. The row subset is
the output `t`, and the column subset is the input `s`. No unitarity is needed. -/
theorem nonzeroMinorStep (R : ZeroArithmetic A) (e : Edge vertices) (M : V → V → A)
    (support : GateSupport R e M) (s t : Subset V)
    (nonzero : minorDet R M (vertices.filter t) (vertices.filter s) ≠ R.zero) :
    Step e s t := by
  have matching := nonzeroMinorMatching R M _ _ nonzero
  apply outsideAgreementStep e s t
  · intro x xl xr
    have forward : t x = true → s x = true := by
      intro tx
      obtain ⟨y, hy, entry⟩ := matchingRowSupport matching x
        (List.mem_filter.mpr ⟨e.complete x, tx⟩)
      have xy : x = y := by
        by_cases same : x = y
        · exact same
        · exact False.elim (entry (support x y same (Or.inl ⟨xl, xr⟩)))
      subst y
      exact (List.mem_filter.mp hy).2
    have reverse : s x = true → t x = true := by
      intro sx
      obtain ⟨y, hy, entry⟩ := matchingColumnSupport matching x
        (List.mem_filter.mpr ⟨e.complete x, sx⟩)
      have yx : y = x := by
        by_cases same : y = x
        · exact same
        · exact False.elim (entry (support y x same (Or.inr ⟨xl, xr⟩)))
      subst y
      exact (List.mem_filter.mp hy).2
    cases ht : t x <;> cases hs : s x <;> simp_all
  · simpa only [cutCount, Bool.true_and, List.countP_eq_length_filter] using matchingLength matching

/-- Fully concrete single-gate support corollary: four arbitrary scalar entries
specify the block and the determinant is the signed expansion defined above. -/
theorem twoCoordinateNonzeroMinorStep (R : ZeroArithmetic A) (e : Edge vertices)
    (a b c d : A) (s t : Subset V)
    (nonzero : minorDet R (twoCoordinateGate R e a b c d)
      (vertices.filter t) (vertices.filter s) ≠ R.zero) : Step e s t :=
  nonzeroMinorStep R e _ (twoCoordinateGateSupport R e a b c d) s t nonzero

/-- In either collision sector, a nonzero minor can only stay at the same subset. -/
theorem nonzeroMinorCollision (R : ZeroArithmetic A) (e : Edge vertices) (M : V → V → A)
    (support : GateSupport R e M) (s t : Subset V) (collision : s e.left = s e.right)
    (nonzero : minorDet R M (vertices.filter t) (vertices.filter s) ≠ R.zero) : t = s := by
  rcases nonzeroMinorStep R e M support s t nonzero with stay | swap
  · exact stay
  · exact swap.trans (swapCollision e s collision)

end Gate

#print axioms choicePerm
#print axioms nonzeroSumTerm
#print axioms nonzeroSignedFactors
#print axioms nonzeroMinorMatching
#print axioms matchingLength
#print axioms matchingRowSupport
#print axioms matchingColumnSupport
#print axioms twoCoordinateGateSupport
#print axioms outsideAgreementStep
#print axioms nonzeroMinorStep
#print axioms twoCoordinateNonzeroMinorStep
#print axioms nonzeroMinorCollision

end MathCommons.ForestGateSupport
