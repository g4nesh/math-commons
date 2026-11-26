import ForestMatrixSumExpansion

/-! Concrete Gaussian conjugation/scaling and actual low-order contraction.
No arbitrary-order Cauchy--Binet or transpose identity is assumed. -/
namespace MathCommons.ForestAdjointContraction
open Classical
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory
open ForestMinorAlternation ForestSubsetExchange ForestPermutationLaw ForestSpectralExpansion
open ForestMatrixSumExpansion

variable {V : Type}

theorem conjugateAdd (z w : Gaussian) :
    conjugate (gaussianArithmetic.add z w) = gaussianArithmetic.add (conjugate z) (conjugate w) := by
  apply Prod.ext <;> simp only [conjugate, gaussianArithmetic] <;> grind

theorem conjugateMul (z w : Gaussian) :
    conjugate (gaussianArithmetic.mul z w) = gaussianArithmetic.mul (conjugate z) (conjugate w) := by
  apply Prod.ext <;> simp only [conjugate, gaussianArithmetic] <;> grind

theorem conjugateSigned (odd : Bool) (z : Gaussian) :
    conjugate (signed gaussianArithmetic odd z) = signed gaussianArithmetic odd (conjugate z) := by
  cases odd <;> apply Prod.ext <;> simp only [conjugate, gaussianArithmetic, signed,
    Bool.false_eq_true, ↓reduceIte] <;> grind

theorem conjugateSum {X : Type} (xs : List X) (f : X → Gaussian) :
    conjugate (scalarSum gaussianArithmetic (xs.map f)) =
      scalarSum gaussianArithmetic (xs.map fun x => conjugate (f x)) := by
  induction xs with
  | nil => apply Prod.ext <;> simp only [List.map_nil, scalarSum, conjugate, gaussianArithmetic] <;> grind
  | cons x xs ih => simp only [List.map_cons, scalarSum, conjugateAdd, ih]

/-- Entrywise conjugation commutes with the actual minor in every order,
including rectangular lists; no transpose statement is assumed. -/
theorem conjugateMinor (M : V → V → Gaussian) (rows cols : List V) :
    conjugate (minorDet gaussianArithmetic M rows cols) =
      minorDet gaussianArithmetic (fun i j => conjugate (M i j)) rows cols := by
  induction rows generalizing cols with
  | nil => cases cols <;> apply Prod.ext <;> simp only [minorDet, conjugate, gaussianArithmetic, List.isEmpty_nil, List.isEmpty_cons, Bool.false_eq_true, ↓reduceIte] <;> grind
  | cons r rows ih =>
    simp only [minorDet, conjugateSum]
    apply gaussianSumCongr
    intro q _
    rw [conjugateSigned, conjugateMul, ih]

def columnScale (M : V → V → Gaussian) (b : V → Gaussian) : V → V → Gaussian :=
  fun i j => gaussianArithmetic.mul (M i j) (b j)

/-- Multiplying each selected column scales the actual minor by the product
of its factors, in every order and without nonzero or distinctness assumptions. -/
theorem columnScaleMinor (M : V → V → Gaussian) (b : V → Gaussian) (rows cols : List V) :
    minorDet gaussianArithmetic (columnScale M b) rows cols =
      gaussianArithmetic.mul (product cols b) (minorDet gaussianArithmetic M rows cols) := by
  induction rows generalizing cols with
  | nil =>
    cases cols with
    | nil => exact (gaussianOneMul gaussianArithmetic.one).symm
    | cons c cs => exact (gaussianArithmetic.mulZero _).symm
  | cons r rows ih =>
    simp only [minorDet]
    rw [← gaussianSumScale]
    apply gaussianSumCongr
    intro q hq
    have prod : product cols b = gaussianArithmetic.mul (b q.column) (product q.remaining b) :=
      (productPerm _ _ b (choicePerm cols q hq)).symm
    rw [ih, prod, ← signedScale]
    congr 1
    apply Prod.ext <;> simp only [columnScale, gaussianArithmetic] <;> grind

/-- Only entries in the selected rectangle affect an actual minor. -/
theorem minorCongrEntries (M N : V → V → Gaussian) (rows cols : List V)
    (same : ∀ i ∈ rows, ∀ j ∈ cols, M i j = N i j) :
    minorDet gaussianArithmetic M rows cols = minorDet gaussianArithmetic N rows cols := by
  induction rows generalizing cols with
  | nil => rfl
  | cons r rows ih =>
    simp only [minorDet]
    apply gaussianSumCongr
    intro q hq
    have member : q.column ∈ cols := (choicePerm cols q hq).mem_iff.mp List.mem_cons_self
    rw [same r List.mem_cons_self q.column member]
    congr 2
    apply ih
    intro i hi j hj
    exact same i (List.mem_cons_of_mem r hi) j
      ((choicePerm cols q hq).mem_iff.mp (List.mem_cons_of_mem q.column hj))

/-- Actual right multiplication by a diagonal matrix scales every selected
minor, provided its columns belong to the finite multiplication universe. -/
theorem rightDiagonalMinor [DecidableEq V] (vertices rows cols : List V) (U : V → V → Gaussian) (b : V → Gaussian)
    (unique : vertices.Nodup) (inside : ∀ j ∈ cols, j ∈ vertices) :
    minorDet gaussianArithmetic (matrixMul vertices U (diagonal b)) rows cols =
      gaussianArithmetic.mul (product cols b) (minorDet gaussianArithmetic U rows cols) := by
  have eq : minorDet gaussianArithmetic (matrixMul vertices U (diagonal b)) rows cols =
      minorDet gaussianArithmetic (columnScale U b) rows cols := by
    apply minorCongrEntries
    intro i _ j hj
    exact matrixMulDiagonal vertices U b unique i j (inside j hj)
  rw [eq, columnScaleMinor]

/-- Direct adjoint-minor identity in order zero. -/
theorem adjointMinorZero (U : V → V → Gaussian) :
    minorDet gaussianArithmetic (adjoint U) [] [] =
      conjugate (minorDet gaussianArithmetic U [] []) := by
  apply Prod.ext <;> simp only [minorDet, conjugate, gaussianArithmetic, List.isEmpty_nil, ↓reduceIte] <;> grind

/-- Direct adjoint-minor identity in order one. -/
theorem adjointMinorOne (U : V → V → Gaussian) (r c : V) :
    minorDet gaussianArithmetic (adjoint U) [r] [c] =
      conjugate (minorDet gaussianArithmetic U [c] [r]) := by
  simp only [singletonMinor, adjoint]

/-- Direct adjoint-minor identity in order two, retaining both supplied orders. -/
theorem adjointMinorTwo (U : V → V → Gaussian) (r s c d : V) :
    minorDet gaussianArithmetic (adjoint U) [r,s] [c,d] =
      conjugate (minorDet gaussianArithmetic U [c,d] [r,s]) := by
  simp only [gaussianTwoByTwo]
  apply Prod.ext <;> simp only [adjoint, blockDet, conjugate, gaussianArithmetic] <;> grind

theorem squaredFactor (z w : Gaussian) :
    gaussianArithmetic.mul (gaussianArithmetic.mul w z) (conjugate z) = scale (normSq z) w := by
  apply Prod.ext <;> simp only [gaussianArithmetic, conjugate, scale, normSq] <;> grind

/-- The order-one principal-minor contraction for arbitrary U and spectra. -/
theorem conjugatedMinorOne [DecidableEq V] (vertices : List V) (U : V → V → Gaussian) (b : V → Gaussian)
    (unique : vertices.Nodup) (r : V) :
    minorDet gaussianArithmetic (conjugatedDiagonal vertices U b) [r] [r] =
      scalarSum gaussianArithmetic (vertices.map fun j =>
        scale (normSq (minorDet gaussianArithmetic U [r] [j])) (b j)) := by
  rw [singletonMinor, conjugatedDiagonalEntry vertices U b unique r r]
  apply gaussianSumCongr
  intro j _
  rw [singletonMinor, gaussianMulComm (U r j) (b j), squaredFactor]

/-- Mixed order-two contraction, including off-principal matrix minors. -/
theorem conjugatedMinorTwoMixed [DecidableEq V] (vertices : List V) (U : V → V → Gaussian)
    (b : V → Gaussian) (unique : vertices.Nodup) (r s c d : V) :
    minorDet gaussianArithmetic (conjugatedDiagonal vertices U b) [r,s] [c,d] =
      scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p =>
        gaussianArithmetic.mul
          (gaussianArithmetic.mul (gaussianArithmetic.mul (b p.1) (b p.2))
            (minorDet gaussianArithmetic U [r,s] [p.1,p.2]))
          (conjugate (minorDet gaussianArithmetic U [c,d] [p.1,p.2]))) := by
  unfold conjugatedDiagonal
  rw [cauchyBinetTwo]
  apply gaussianSumCongr
  intro p hp
  have sub := (pairMemberSublist vertices p).mp hp
  have inside : ∀ j ∈ [p.1,p.2], j ∈ vertices := fun j hj => sub.subset hj
  rw [rightDiagonalMinor vertices [r,s] [p.1,p.2] U b unique inside, adjointMinorTwo]
  simp only [product, gaussianMulOne]

/-- The order-two principal-minor contraction for arbitrary Gaussian U in
arbitrary ambient dimension. Every unordered intermediate pair occurs once. -/
theorem conjugatedMinorTwo [DecidableEq V] (vertices : List V) (U : V → V → Gaussian) (b : V → Gaussian)
    (unique : vertices.Nodup) (r s : V) :
    minorDet gaussianArithmetic (conjugatedDiagonal vertices U b) [r,s] [r,s] =
      scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p =>
        scale (normSq (minorDet gaussianArithmetic U [r,s] [p.1,p.2]))
          (gaussianArithmetic.mul (b p.1) (b p.2))) := by
  rw [conjugatedMinorTwoMixed vertices U b unique r s r s]
  apply gaussianSumCongr
  intro p _
  exact squaredFactor _ _

/-- Order-zero contraction, including an empty multiplication universe. -/
theorem conjugatedMinorZero [DecidableEq V] (vertices : List V) (U : V → V → Gaussian)
    (b : V → Gaussian) :
    minorDet gaussianArithmetic (conjugatedDiagonal vertices U b) [] [] = gaussianArithmetic.one := rfl

/-- The order-two result at a canonical retained Boolean subset, with the
selected row list supplied by the actual coordinate filter. -/
theorem conjugatedTwoSubset [DecidableEq V] (vertices : List V) (U : V → V → Gaussian)
    (b : V → Gaussian) (unique : vertices.Nodup) (finish : Subset V)
    (size : (vertices.filter finish).length = 2) :
    minorDet gaussianArithmetic (conjugatedDiagonal vertices U b)
        (vertices.filter finish) (vertices.filter finish) =
      scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p =>
        scale (normSq (minorDet gaussianArithmetic U (vertices.filter finish) [p.1,p.2]))
          (gaussianArithmetic.mul (b p.1) (b p.2))) := by
  cases filtered : vertices.filter finish with
  | nil => simp [filtered] at size
  | cons r rows =>
    cases rows with
    | nil => simp [filtered] at size
    | cons s rest =>
      have empty : rest = [] := by simpa [filtered] using size
      subst rest
      exact conjugatedMinorTwo vertices U b unique r s

/-- The literal order-two diagonal-plus-matrix determinant separates into
its zero-, one- and two-coordinate principal contributions. -/
theorem twoDiagonalSum [DecidableEq V] (M : V → V → Gaussian) (a : V → Gaussian)
    (r s : V) (different : r ≠ s) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) M) [r,s] [r,s] =
      gaussianArithmetic.add (gaussianArithmetic.mul (a r) (a s))
        (gaussianArithmetic.add
          (gaussianArithmetic.add (gaussianArithmetic.mul (a r) (M s s))
            (gaussianArithmetic.mul (a s) (M r r)))
          (minorDet gaussianArithmetic M [r,s] [r,s])) := by
  simp only [gaussianTwoByTwo, blockDet, matrixAdd, diagonal, different,
    Ne.symm different, ↓reduceIte, gaussianZeroAdd]
  apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

/-- A closed actual two-coordinate spectral determinant expansion in any
ambient dimension; every coefficient is a squared actual minor of U. -/
theorem twoSpectralPencil [DecidableEq V] (vertices : List V) (U : V → V → Gaussian)
    (a b : V → Gaussian) (unique : vertices.Nodup) (r s : V) (different : r ≠ s) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) (conjugatedDiagonal vertices U b))
        [r,s] [r,s] =
      gaussianArithmetic.add (gaussianArithmetic.mul (a r) (a s))
        (gaussianArithmetic.add
          (gaussianArithmetic.add
            (gaussianArithmetic.mul (a r)
              (scalarSum gaussianArithmetic (vertices.map fun j =>
                scale (normSq (minorDet gaussianArithmetic U [s] [j])) (b j))))
            (gaussianArithmetic.mul (a s)
              (scalarSum gaussianArithmetic (vertices.map fun j =>
                scale (normSq (minorDet gaussianArithmetic U [r] [j])) (b j)))))
          (scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p =>
            scale (normSq (minorDet gaussianArithmetic U [r,s] [p.1,p.2]))
              (gaussianArithmetic.mul (b p.1) (b p.2))))) := by
  rw [twoDiagonalSum _ a r s different, conjugatedMinorTwo vertices U b unique r s]
  have hr := conjugatedMinorOne vertices U b unique r
  have hs := conjugatedMinorOne vertices U b unique s
  simp only [singletonMinor] at hr hs
  rw [hr, hs]
  simp only [singletonMinor]

/-- The same normalized full-permutation forest law now gives the coefficients
of an ACTUAL conjugated principal minor at every two-element retained subset.
This is a low-order matrix identity, not a general determinant-of-sum theorem. -/
theorem forestTwoContraction [DecidableEq V] {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k))
    (b : V → Gaussian) :
    (∀ entry ∈ coinLaw edges bs n, 0 ≤ entry.2) ∧
    sumWith (fun entry => entry.2) (coinLaw edges bs n) = 1 ∧
    ∀ finish : Subset V, (vertices.filter finish).length = 2 →
      minorDet gaussianArithmetic
        (conjugatedDiagonal vertices (matrixProduct vertices (gateSchedule edges as bs cs ds) n) b)
        (vertices.filter finish) (vertices.filter finish) =
      scalarSum gaussianArithmetic ((orderedPairs vertices).map fun p =>
        scale (sumWith (fun entry => if subsetImage entry.1 (pairSubset p) = finish then entry.2 else 0)
          (coinLaw edges bs n)) (gaussianArithmetic.mul (b p.1) (b p.2))) := by
  have law := commonSquaredMinorLaw edges as bs cs ds n bridges unitary
  refine ⟨law.1, law.2.1, ?_⟩
  intro finish size
  rw [conjugatedTwoSubset vertices _ b (edges 0).nodup finish size]
  apply gaussianSumCongr
  intro p hp
  rw [← pairFilter vertices (edges 0).nodup p hp, law.2.2 (pairSubset p) finish]

#print axioms conjugateAdd
#print axioms conjugateMul
#print axioms conjugateSigned
#print axioms conjugateSum
#print axioms conjugateMinor
#print axioms columnScaleMinor
#print axioms minorCongrEntries
#print axioms rightDiagonalMinor
#print axioms adjointMinorZero
#print axioms adjointMinorOne
#print axioms adjointMinorTwo
#print axioms squaredFactor
#print axioms conjugatedMinorOne
#print axioms conjugatedMinorTwoMixed
#print axioms conjugatedMinorTwo
#print axioms conjugatedMinorZero
#print axioms conjugatedTwoSubset
#print axioms twoDiagonalSum
#print axioms twoSpectralPencil
#print axioms forestTwoContraction

end MathCommons.ForestAdjointContraction
