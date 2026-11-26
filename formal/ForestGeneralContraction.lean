import ForestAdjointContraction

/-! Actual arbitrary-order transpose/adjoint and the full rational forest
 determinant law. Contraction is derived for a left gate product against an
 arbitrary right matrix; generic arbitrary-matrix Cauchy--Binet and transfer
 to arbitrary complex scalars are not claimed. -/
namespace MathCommons.ForestGeneralContraction
open Classical
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory
open ForestMinorAlternation ForestSubsetExchange ForestPermutationLaw ForestSpectralExpansion
open ForestMatrixSumExpansion ForestAdjointContraction

variable {V : Type}

theorem signedMulSum {X : Type} (xs : List X) (f : X → Gaussian) (odd : Bool) (a : Gaussian) :
    signed gaussianArithmetic odd (gaussianArithmetic.mul a (scalarSum gaussianArithmetic (xs.map f))) =
      scalarSum gaussianArithmetic (xs.map fun x => signed gaussianArithmetic odd (gaussianArithmetic.mul a (f x))) := by
  rw [← gaussianSumScale, signedSum]

theorem signedCross (odd₁ odd₂ : Bool) (a b d : Gaussian) :
    signed gaussianArithmetic (!odd₁) (gaussianArithmetic.mul a
      (signed gaussianArithmetic odd₂ (gaussianArithmetic.mul b d))) =
    signed gaussianArithmetic (!odd₂) (gaussianArithmetic.mul b
      (signed gaussianArithmetic odd₁ (gaussianArithmetic.mul a d))) := by
  cases odd₁ <;> cases odd₂ <;> apply Prod.ext <;>
    simp only [signed, Bool.not_false, Bool.not_true, Bool.false_eq_true, ↓reduceIte,
      gaussianArithmetic] <;> grind

theorem minorHeadSplit (M : V → V → Gaussian) (r c : V) (rows cols : List V) :
    minorDet gaussianArithmetic M (r :: rows) (c :: cols) =
      gaussianArithmetic.add (gaussianArithmetic.mul (M r c) (minorDet gaussianArithmetic M rows cols))
        (scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
          signed gaussianArithmetic (!q.odd) (gaussianArithmetic.mul (M r q.column)
            (minorDet gaussianArithmetic M rows (c :: q.remaining))))) := by
  simp only [minorDet, columnChoices, List.map_cons, List.map_map, Function.comp_def,
    scalarSum, signed, Bool.false_eq_true, ↓reduceIte]

/-- Column Laplace expansion derived from the actual row-recursive determinant,
by interchanging the two finite choice sums. No square or uniqueness premise. -/
theorem firstColumnExpansion (M : V → V → Gaussian) (rows cols : List V) (c : V) :
    minorDet gaussianArithmetic M rows (c :: cols) =
      scalarSum gaussianArithmetic ((columnChoices rows).map fun q =>
        signed gaussianArithmetic q.odd (gaussianArithmetic.mul (M q.column c)
          (minorDet gaussianArithmetic M q.remaining cols))) := by
  induction rows generalizing cols with
  | nil => rfl
  | cons r rows ih =>
    rw [minorHeadSplit]
    simp only [columnChoices, List.map_cons, List.map_map, Function.comp_def, scalarSum]
    change gaussianArithmetic.add _ _ = gaussianArithmetic.add
      (gaussianArithmetic.mul (M r c) (minorDet gaussianArithmetic M rows cols))
      (scalarSum gaussianArithmetic ((columnChoices rows).map fun p =>
        signed gaussianArithmetic (!p.odd) (gaussianArithmetic.mul (M p.column c)
          (minorDet gaussianArithmetic M (r :: p.remaining) cols))))
    congr 1
    calc
      _ = scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
          scalarSum gaussianArithmetic ((columnChoices rows).map fun p =>
            signed gaussianArithmetic (!q.odd) (gaussianArithmetic.mul (M r q.column)
              (signed gaussianArithmetic p.odd (gaussianArithmetic.mul (M p.column c)
                (minorDet gaussianArithmetic M p.remaining q.remaining)))))) := by
        apply gaussianSumCongr
        intro q _
        rw [ih, signedMulSum]
      _ = scalarSum gaussianArithmetic ((columnChoices rows).map fun p =>
          scalarSum gaussianArithmetic ((columnChoices cols).map fun q =>
            signed gaussianArithmetic (!q.odd) (gaussianArithmetic.mul (M r q.column)
              (signed gaussianArithmetic p.odd (gaussianArithmetic.mul (M p.column c)
                (minorDet gaussianArithmetic M p.remaining q.remaining)))))) := gaussianSumSwap _ _ _
      _ = _ := by
        apply gaussianSumCongr
        intro p _
        rw [minorDet, signedMulSum]
        apply gaussianSumCongr
        intro q _
        exact signedCross q.odd p.odd _ _ _

def transpose (M : V → V → Gaussian) : V → V → Gaussian := fun i j => M j i

/-- Transposition exchanges the two supplied index lists of every actual minor.
This includes rectangular and duplicate-index lists, with no determinant axiom. -/
theorem transposeMinor (M : V → V → Gaussian) (rows cols : List V) :
    minorDet gaussianArithmetic (transpose M) rows cols = minorDet gaussianArithmetic M cols rows := by
  induction rows generalizing cols with
  | nil => cases cols <;> rfl
  | cons r rows ih =>
    rw [firstColumnExpansion]
    change scalarSum gaussianArithmetic ((columnChoices cols).map _) = _
    apply gaussianSumCongr
    intro q _
    rw [ih]
    rfl

/-- The actual adjoint-minor identity in arbitrary order, now derived using
entrywise conjugation and the proved transpose theorem. -/
theorem adjointMinor (M : V → V → Gaussian) (rows cols : List V) :
    minorDet gaussianArithmetic (adjoint M) rows cols =
      conjugate (minorDet gaussianArithmetic M cols rows) := by
  rw [← transposeMinor M rows cols, conjugateMinor]
  rfl



/-- Associativity of the actual finite matrix product, derived from finite sums. -/
theorem matrixMulAssoc (vertices : List V) (M N P : V → V → Gaussian) :
    matrixMul vertices (matrixMul vertices M N) P = matrixMul vertices M (matrixMul vertices N P) := by
  funext i j
  simp only [matrixMul, gaussianSumScaleRight]
  rw [gaussianSumSwap]
  apply gaussianSumCongr
  intro k _
  rw [← gaussianSumScale]
  apply gaussianSumCongr
  intro l _
  exact gaussianMulAssoc _ _ _

/-- Arbitrary-order Cauchy--Binet with an actual left product of two-coordinate
matrices and an unrestricted right matrix. Sparse gate recurrence replaces
unproved generic ordered-tuple grouping. No forest or unitary premise. -/
theorem gateProductContraction [DecidableEq V] {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat) (N : V → V → Gaussian)
    (finish : Subset V) (cols : List V) :
    minorDet gaussianArithmetic
      (matrixMul vertices (matrixProduct vertices (gateSchedule edges as bs cs ds) n) N)
      (vertices.filter finish) cols =
    scalarSum gaussianArithmetic ((subsets vertices).map fun start =>
      gaussianArithmetic.mul
        (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
          (vertices.filter finish) (vertices.filter start))
        (minorDet gaussianArithmetic N (vertices.filter start) cols)) := by
  induction n generalizing finish with
  | zero =>
    simp only [matrixProduct, matrixIdentityMul vertices (edges 0).nodup (edges 0).complete]
    have one := scalarSumUniqueKey (subsets vertices) (fun x => x)
      (by simpa using subsetsNodup vertices (edges 0).nodup)
      finish (subsetsComplete vertices (edges 0).complete finish)
      (fun start => gaussianArithmetic.mul
        (minorDet gaussianArithmetic identityMatrix (vertices.filter finish) (vertices.filter start))
        (minorDet gaussianArithmetic N (vertices.filter start) cols)) (by
          intro start _ different
          dsimp only
          rw [identitySubsetMinor vertices (edges 0).nodup (edges 0).complete]
          simp only [Ne.symm different, ↓reduceIte, gaussianArithmetic.zeroMul])
    simpa only [identitySubsetMinor vertices (edges 0).nodup (edges 0).complete,
      ↓reduceIte, gaussianOneMul] using one.symm
  | succ n ih =>
    have cb (M : V → V → Gaussian) (t : Subset V) (columns : List V) :
        minorDet gaussianArithmetic (matrixMul vertices (gateSchedule edges as bs cs ds n) M)
          (vertices.filter t) columns =
        scalarSum gaussianArithmetic ((successors (edges n) t).map fun pred =>
          gaussianArithmetic.mul
            (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
              (vertices.filter t) (vertices.filter pred))
            (minorDet gaussianArithmetic M (vertices.filter pred) columns)) :=
      gateMinorPredecessors (edges n) (as n) (bs n) (cs n) (ds n) M t columns
    simp only [matrixProduct]
    rw [matrixMulAssoc, cb]
    calc
      _ = scalarSum gaussianArithmetic ((successors (edges n) finish).map fun pred =>
          scalarSum gaussianArithmetic ((subsets vertices).map fun start =>
            gaussianArithmetic.mul
              (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
                (vertices.filter finish) (vertices.filter pred))
              (gaussianArithmetic.mul
                (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
                  (vertices.filter pred) (vertices.filter start))
                (minorDet gaussianArithmetic N (vertices.filter start) cols)))) := by
        apply gaussianSumCongr
        intro pred _
        rw [ih, ← gaussianSumScale]
      _ = scalarSum gaussianArithmetic ((subsets vertices).map fun start =>
          scalarSum gaussianArithmetic ((successors (edges n) finish).map fun pred =>
            gaussianArithmetic.mul
              (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
                (vertices.filter finish) (vertices.filter pred))
              (gaussianArithmetic.mul
                (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
                  (vertices.filter pred) (vertices.filter start))
                (minorDet gaussianArithmetic N (vertices.filter start) cols)))) := gaussianSumSwap _ _ _
      _ = _ := by
        apply gaussianSumCongr
        intro start _
        rw [cb, gaussianSumScaleRight]
        apply gaussianSumCongr
        intro pred _
        exact (gaussianMulAssoc _ _ _).symm



/-- Actual multiplication on the left by a diagonal matrix scales rows. -/
theorem matrixDiagonalMul [DecidableEq V] (vertices : List V) (M : V → V → Gaussian)
    (b : V → Gaussian) (unique : vertices.Nodup) (i j : V) (member : i ∈ vertices) :
    matrixMul vertices (diagonal b) M i j = gaussianArithmetic.mul (b i) (M i j) := by
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k =>
        if k = i then gaussianArithmetic.mul (b i) (M i j) else gaussianArithmetic.zero) := by
      apply gaussianSumCongr
      intro k _
      by_cases same : k = i
      · subst k; simp [diagonal]
      · simp [diagonal, same, Ne.symm same, gaussianArithmetic.zeroMul]
    _ = _ := gaussianSumDelta vertices unique i member _

/-- Left diagonal minor scaling at arbitrary order, now derived using transpose. -/
theorem leftDiagonalMinor [DecidableEq V] (vertices rows cols : List V) (M : V → V → Gaussian)
    (b : V → Gaussian) (unique : vertices.Nodup) (inside : ∀ i ∈ rows, i ∈ vertices) :
    minorDet gaussianArithmetic (matrixMul vertices (diagonal b) M) rows cols =
      gaussianArithmetic.mul (product rows b) (minorDet gaussianArithmetic M rows cols) := by
  rw [← transposeMinor (matrixMul vertices (diagonal b) M) cols rows]
  have eq : minorDet gaussianArithmetic (transpose (matrixMul vertices (diagonal b) M)) cols rows =
      minorDet gaussianArithmetic (columnScale (transpose M) b) cols rows := by
    apply minorCongrEntries
    intro i _ j hj
    change matrixMul vertices (diagonal b) M j i = gaussianArithmetic.mul (M j i) (b j)
    rw [matrixDiagonalMul vertices M b unique j i (inside j hj), gaussianMulComm]
  rw [eq, columnScaleMinor, transposeMinor]

/-- The filtered product is the canonical selected-coordinate spectral monomial. -/
theorem productFilter (vertices : List V) (s : Subset V) (b : V → Gaussian) :
    product (vertices.filter s) b =
      product vertices (fun x => if s x then b x else gaussianArithmetic.one) := by
  induction vertices with
  | nil => rfl
  | cons x xs ih =>
    cases hs : s x <;> simp only [List.filter_cons, hs, Bool.false_eq_true, ↓reduceIte,
      product, gaussianOneMul, ih]

/-- Every principal minor of U D_b U* for an actual gate product U expands
into squared actual U minors in every order. No forest or unitary premise. -/
theorem gateConjugatedMinor [DecidableEq V] {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat) (b : V → Gaussian) (finish : Subset V) :
    minorDet gaussianArithmetic
      (conjugatedDiagonal vertices (matrixProduct vertices (gateSchedule edges as bs cs ds) n) b)
      (vertices.filter finish) (vertices.filter finish) =
    scalarSum gaussianArithmetic ((subsets vertices).map fun start =>
      scale (normSq (minorDet gaussianArithmetic
        (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
        (vertices.filter finish) (vertices.filter start)))
        (product vertices (fun x => if start x then b x else gaussianArithmetic.one))) := by
  unfold conjugatedDiagonal
  rw [matrixMulAssoc, gateProductContraction]
  apply gaussianSumCongr
  intro start _
  rw [leftDiagonalMinor vertices (vertices.filter start) (vertices.filter finish) _ b
    (edges 0).nodup (fun _ h => (List.mem_filter.mp h).1), adjointMinor]
  rw [← gaussianMulAssoc, gaussianMulComm
    (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start)) (product (vertices.filter start) b),
    squaredFactor, productFilter]

/-- Multiplying a Gaussian value by a rationally scaled Gaussian value. -/
theorem mulScale (q : Rat) (z w : Gaussian) :
    gaussianArithmetic.mul z (scale q w) = scale q (gaussianArithmetic.mul z w) := by
  apply Prod.ext <;> simp only [scale, gaussianArithmetic] <;> grind

/-- The actual determinant-of-sum spectral expansion for every gate product,
with all subset sizes and the same start/finish order as the permutation law. -/
theorem gateSpectralDeterminant [DecidableEq V] {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat) (a b : V → Gaussian) :
    minorDet gaussianArithmetic
      (matrixAdd (diagonal a)
        (conjugatedDiagonal vertices (matrixProduct vertices (gateSchedule edges as bs cs ds) n) b))
      vertices vertices =
    scalarSum gaussianArithmetic ((subsets vertices).map fun start =>
      scalarSum gaussianArithmetic ((subsets vertices).map fun finish =>
        scale (normSq (minorDet gaussianArithmetic
          (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
          (vertices.filter finish) (vertices.filter start)))
          (spectralMonomial vertices start finish a b))) := by
  rw [spectralPrincipalExpansion vertices a b _ (edges 0).nodup]
  calc
    _ = scalarSum gaussianArithmetic ((subsets vertices).map fun finish =>
        scalarSum gaussianArithmetic ((subsets vertices).map fun start =>
          scale (normSq (minorDet gaussianArithmetic
            (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
            (vertices.filter finish) (vertices.filter start)))
            (spectralMonomial vertices start finish a b))) := by
      apply gaussianSumCongr
      intro finish _
      rw [gateConjugatedMinor, ← gaussianSumScale]
      apply gaussianSumCongr
      intro start _
      exact mulScale _ _ _
    _ = _ := gaussianSumSwap _ _ _

/-- The ACTUAL forest determinant is the convex weighted sum of spectral
pairing products from one normalized nonnegative full-permutation law.
No determinant expansion, transpose identity or contraction is a premise. -/
theorem forestDeterminantLaw [DecidableEq V] {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k))
    (a b : V → Gaussian) :
    (∀ entry ∈ coinLaw edges bs n, 0 ≤ entry.2) ∧
    sumWith (fun entry => entry.2) (coinLaw edges bs n) = 1 ∧
    minorDet gaussianArithmetic
      (matrixAdd (diagonal a)
        (conjugatedDiagonal vertices (matrixProduct vertices (gateSchedule edges as bs cs ds) n) b))
      vertices vertices = weightedPairings vertices (coinLaw edges bs n) a b := by
  have law := forestWeightedSubsetExpansion edges as bs cs ds n bridges unitary a b
  refine ⟨law.1, law.2.1, ?_⟩
  rw [gateSpectralDeterminant]
  exact law.2.2.symm

/-- A single explicit law works for ALL spectra simultaneously; the spectral
quantifiers follow nonnegativity and normalization in the conclusion. -/
theorem forestCommonDeterminantLaw [DecidableEq V] {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    (∀ entry ∈ coinLaw edges bs n, 0 ≤ entry.2) ∧
    sumWith (fun entry => entry.2) (coinLaw edges bs n) = 1 ∧
    ∀ a b : V → Gaussian,
      minorDet gaussianArithmetic
        (matrixAdd (diagonal a)
          (conjugatedDiagonal vertices (matrixProduct vertices (gateSchedule edges as bs cs ds) n) b))
        vertices vertices = weightedPairings vertices (coinLaw edges bs n) a b := by
  have law := commonSquaredMinorLaw edges as bs cs ds n bridges unitary
  exact ⟨law.1, law.2.1, fun a b => (forestDeterminantLaw edges as bs cs ds n bridges unitary a b).2.2⟩

#print axioms signedMulSum
#print axioms signedCross
#print axioms minorHeadSplit
#print axioms firstColumnExpansion
#print axioms transposeMinor
#print axioms adjointMinor
#print axioms matrixMulAssoc
#print axioms gateProductContraction
#print axioms matrixDiagonalMul
#print axioms leftDiagonalMinor
#print axioms productFilter
#print axioms gateConjugatedMinor
#print axioms mulScale
#print axioms gateSpectralDeterminant
#print axioms forestDeterminantLaw
#print axioms forestCommonDeterminantLaw

end MathCommons.ForestGeneralContraction
