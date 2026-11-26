import ForestSpectralExpansion

/-! Actual determinant expansion under arbitrary diagonal updates.
The remaining conjugated-matrix Cauchy--Binet step is not assumed. -/
namespace MathCommons.ForestMatrixSumExpansion
open Classical
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory
open ForestMinorAlternation ForestSubsetExchange ForestPermutationLaw ForestSpectralExpansion

variable {V : Type} [DecidableEq V]

def matrixAdd (M N : V → V → Gaussian) : V → V → Gaussian :=
  fun i j => gaussianArithmetic.add (M i j) (N i j)

def diagonal (a : V → Gaussian) : V → V → Gaussian :=
  fun i j => if i = j then a i else gaussianArithmetic.zero

def adjoint (M : V → V → Gaussian) : V → V → Gaussian :=
  fun i j => conjugate (M j i)

def diagonalUpdate (M : V → V → Gaussian) (x : V) (a : Gaussian) : V → V → Gaussian :=
  replaceRow M x (fun j => gaussianArithmetic.add (M x j)
    (if x = j then a else gaussianArithmetic.zero))

def partialDiagonal (updates : List V) (a : V → Gaussian) (M : V → V → Gaussian) : V → V → Gaussian :=
  fun i j => gaussianArithmetic.add (M i j)
    (if i ∈ updates ∧ i = j then a i else gaussianArithmetic.zero)

/-- Removing a distinguished element preserves the relative order of all others. -/
theorem filterWithout (before after : List V) (x : V)
    (unique : (before ++ x :: after).Nodup) :
    (before ++ x :: after).filter (fun y => y != x) = before ++ after := by
  have absent : x ∉ before ++ after := by
    have h := List.nodup_append.mp unique
    intro hx
    rcases List.mem_append.mp hx with hb | ha
    · exact h.2.2 x hb x List.mem_cons_self rfl
    · exact (List.nodup_cons.mp h.2.1).1 ha
  have bef : before.filter (fun y => y != x) = before := by
    apply List.filter_eq_self.mpr
    intro y hy
    have ne : y ≠ x := by intro eq; subst y; exact absent (List.mem_append_left _ hy)
    simp [ne]
  have aft : after.filter (fun y => y != x) = after := by
    apply List.filter_eq_self.mpr
    intro y hy
    have ne : y ≠ x := by intro eq; subst y; exact absent (List.mem_append_right _ hy)
    simp [ne]
  simp [List.filter_append, bef, aft]

/-- A diagonal replacement row has a positive principal cofactor: the row and
column position signs cancel, for an arbitrary position in the ambient list. -/
theorem diagonalRowMinor (M : V → V → Gaussian) (rows : List V) (x : V) (a : Gaussian)
    (unique : rows.Nodup) (member : x ∈ rows) :
    minorDet gaussianArithmetic (replaceRow M x (fun j => if x = j then a else gaussianArithmetic.zero)) rows rows =
      gaussianArithmetic.mul a (minorDet gaussianArithmetic M
        (rows.filter (fun y => y != x)) (rows.filter (fun y => y != x))) := by
  obtain ⟨before, after, eq⟩ := List.mem_iff_append.mp member
  subst rows
  let N := replaceRow M x (fun j => if x = j then a else gaussianArithmetic.zero)
  have absent : x ∉ before ++ after := by
    have h := List.nodup_append.mp unique
    intro hx
    rcases List.mem_append.mp hx with hb | ha
    · exact h.2.2 x hb x List.mem_cons_self rfl
    · exact (List.nodup_cons.mp h.2.1).1 ha
  change minorDet gaussianArithmetic N (before ++ x :: after) (before ++ x :: after) = _
  have move := minorMoveRow N [] before after x (before ++ x :: after)
  simp only [List.nil_append] at move
  rw [move, minorSelectedColumn N x x (before ++ after) before after unique]
  · rw [signedTwice, filterWithout before after x unique]
    have entry : N x x = a := by simp [N, replaceRow]
    rw [entry]
    congr 1
    exact minorReplaceAbsent M (before ++ after) (before ++ after) x _ absent
  · intro y hy
    have ne : x ≠ y := by intro eq; subst y; exact absent hy
    simp [N, replaceRow, ne]

/-- Exact affine dependence on one diagonal entry of the actual signed minor. -/
theorem singleDiagonalUpdate (M : V → V → Gaussian) (rows : List V) (x : V) (a : Gaussian)
    (unique : rows.Nodup) (member : x ∈ rows) :
    minorDet gaussianArithmetic (diagonalUpdate M x a) rows rows =
      gaussianArithmetic.add (minorDet gaussianArithmetic M rows rows)
        (gaussianArithmetic.mul a (minorDet gaussianArithmetic M
          (rows.filter (fun y => y != x)) (rows.filter (fun y => y != x)))) := by
  unfold diagonalUpdate
  rw [minorRowAdd M x (M x) _ rows rows unique member, replaceOwnRow,
    diagonalRowMinor M rows x a unique member]

theorem partialDiagonalNil (a : V → Gaussian) (M : V → V → Gaussian) :
    partialDiagonal [] a M = M := by
  funext i j
  simp only [partialDiagonal, List.not_mem_nil, false_and, ↓reduceIte, gaussianAddZero]

theorem partialDiagonalCons (updates : List V) (a : V → Gaussian) (M : V → V → Gaussian)
    (x : V) (absent : x ∉ updates) :
    partialDiagonal (x :: updates) a M = diagonalUpdate (partialDiagonal updates a M) x (a x) := by
  funext i j
  by_cases hit : i = x
  · subst i
    by_cases same : x = j
    · subst j; simp [partialDiagonal, diagonalUpdate, replaceRow, absent, gaussianAddZero]
    · simp [partialDiagonal, diagonalUpdate, replaceRow, absent, same, gaussianAddZero]
  · simp [partialDiagonal, diagonalUpdate, replaceRow, hit]


def diagonalWeight (updates : List V) (a : V → Gaussian) (s : Subset V) : Gaussian :=
  product updates (fun x => if s x then a x else gaussianArithmetic.one)

def remaining (rows : List V) (s : Subset V) : List V := rows.filter (fun x => !(s x))

theorem remainingInsert (rows : List V) (s : Subset V) (x : V) :
    remaining rows (ForestSpectralExpansion.insert x s) = remaining (rows.filter (fun y => y != x)) s := by
  simp only [remaining, List.filter_filter]
  apply List.filter_congr
  intro y _
  by_cases same : y = x <;> simp [ForestSpectralExpansion.insert, same]

theorem diagonalWeightExcluded (updates : List V) (a : V → Gaussian) (s : Subset V)
    (x : V) (absent : x ∉ updates) (member : s ∈ subsets updates) :
    diagonalWeight (x :: updates) a s = diagonalWeight updates a s := by
  simp only [diagonalWeight, product, subsetOutside updates s member x absent,
    Bool.false_eq_true, ↓reduceIte, gaussianOneMul]

theorem diagonalWeightIncluded (updates : List V) (a : V → Gaussian) (s : Subset V)
    (x : V) (absent : x ∉ updates) :
    diagonalWeight (x :: updates) a (ForestSpectralExpansion.insert x s) =
      gaussianArithmetic.mul (a x) (diagonalWeight updates a s) := by
  simp only [diagonalWeight, product, ForestSpectralExpansion.insert, ↓reduceIte]
  congr 1
  apply productCongr
  intro y hy
  have ne : y ≠ x := by intro eq; subst y; exact absent hy
  simp only [ne, ↓reduceIte]

/-- Every principal minor of an arbitrary Gaussian matrix under any finite
set of distinct diagonal updates expands into its smaller principal minors. -/
theorem partialDiagonalExpansion (updates rows : List V) (a : V → Gaussian) (M : V → V → Gaussian)
    (uniqueUpdates : updates.Nodup) (uniqueRows : rows.Nodup)
    (inside : ∀ x ∈ updates, x ∈ rows) :
    minorDet gaussianArithmetic (partialDiagonal updates a M) rows rows =
      scalarSum gaussianArithmetic ((subsets updates).map fun s =>
        gaussianArithmetic.mul (diagonalWeight updates a s)
          (minorDet gaussianArithmetic M (remaining rows s) (remaining rows s))) := by
  induction updates generalizing rows with
  | nil =>
    simp only [partialDiagonalNil, subsets, List.map_cons, List.map_nil, scalarSum,
      diagonalWeight, product, remaining, Bool.not_false,
      gaussianOneMul, gaussianAddZero]
    have filt : rows.filter (fun _ => true) = rows := List.filter_eq_self.mpr (by intros; rfl)
    rw [filt]
  | cons x updates ih =>
    obtain ⟨absent, uniqueTail⟩ := List.nodup_cons.mp uniqueUpdates
    have member := inside x List.mem_cons_self
    have insideTail : ∀ y ∈ updates, y ∈ rows := fun y hy => inside y (List.mem_cons_of_mem x hy)
    have insideReduced : ∀ y ∈ updates, y ∈ rows.filter (fun z => z != x) := by
      intro y hy
      have ne : y ≠ x := by intro eq; subst y; exact absent hy
      simp only [List.mem_filter, bne_iff_ne]
      exact ⟨insideTail y hy, ne⟩
    rw [partialDiagonalCons updates a M x absent,
      singleDiagonalUpdate (partialDiagonal updates a M) rows x (a x) uniqueRows member,
      ih rows uniqueTail uniqueRows insideTail,
      ih (rows.filter (fun y => y != x)) uniqueTail (List.Pairwise.filter _ uniqueRows) insideReduced]
    simp only [subsets, gaussianSumAppend, List.map_map, Function.comp_def]
    have excluded :
        scalarSum gaussianArithmetic ((subsets updates).map fun s =>
          gaussianArithmetic.mul (diagonalWeight (x :: updates) a s)
            (minorDet gaussianArithmetic M (remaining rows s) (remaining rows s))) =
        scalarSum gaussianArithmetic ((subsets updates).map fun s =>
          gaussianArithmetic.mul (diagonalWeight updates a s)
            (minorDet gaussianArithmetic M (remaining rows s) (remaining rows s))) := by
      apply gaussianSumCongr
      intro s hs
      rw [diagonalWeightExcluded updates a s x absent hs]
    have included :
        scalarSum gaussianArithmetic ((subsets updates).map fun s =>
          gaussianArithmetic.mul (diagonalWeight (x :: updates) a (ForestSpectralExpansion.insert x s))
            (minorDet gaussianArithmetic M (remaining rows (ForestSpectralExpansion.insert x s)) (remaining rows (ForestSpectralExpansion.insert x s)))) =
        gaussianArithmetic.mul (a x)
          (scalarSum gaussianArithmetic ((subsets updates).map fun s =>
            gaussianArithmetic.mul (diagonalWeight updates a s)
              (minorDet gaussianArithmetic M (remaining (rows.filter (fun y => y != x)) s)
                (remaining (rows.filter (fun y => y != x)) s)))) := by
      rw [← gaussianSumScale]
      apply gaussianSumCongr
      intro s _
      rw [diagonalWeightIncluded updates a s x absent, remainingInsert, gaussianMulAssoc]
    rw [excluded, included]

/-- Canonical finite-subset principal-minor expansion of the ACTUAL determinant
of a diagonal matrix plus an arbitrary Gaussian matrix. No invertibility needed. -/
theorem diagonalSumExpansion (rows : List V) (a : V → Gaussian) (M : V → V → Gaussian)
    (unique : rows.Nodup) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) M) rows rows =
      scalarSum gaussianArithmetic ((subsets rows).map fun s =>
        gaussianArithmetic.mul (diagonalWeight rows a s)
          (minorDet gaussianArithmetic M (remaining rows s) (remaining rows s))) := by
  have eq : minorDet gaussianArithmetic (matrixAdd (diagonal a) M) rows rows =
      minorDet gaussianArithmetic (partialDiagonal rows a M) rows rows := by
    apply minorCongrRows
    intro i hi
    funext j
    simp only [matrixAdd, diagonal, partialDiagonal, hi, true_and]
    exact gaussianAddComm _ _
  rw [eq]
  exact partialDiagonalExpansion rows rows a M unique unique (fun _ h => h)


/-- Complement inside the supplied finite universe, false outside it. -/
def complementSubset (rows : List V) (s : Subset V) : Subset V :=
  fun x => if x ∈ rows then !(s x) else false

theorem complementMember (rows : List V) (s : Subset V) :
    complementSubset rows s ∈ subsets rows := by
  apply subsetMember
  intro x outside
  simp [complementSubset, outside]

theorem complementTwice (rows : List V) (s : Subset V) (member : s ∈ subsets rows) :
    complementSubset rows (complementSubset rows s) = s := by
  funext x
  by_cases hx : x ∈ rows
  · simp [complementSubset, hx]
  · simp [complementSubset, hx, subsetOutside rows s member x hx]

/-- Complement permutes the canonical finite subset list without multiplicity. -/
theorem complementsPerm (rows : List V) (unique : rows.Nodup) :
    ((subsets rows).map (complementSubset rows)).Perm (subsets rows) := by
  have nd := subsetsNodup rows unique
  have mapped := nodupMapOn (subsets rows) (complementSubset rows) nd
    (fun s hs t ht eq => by
      have h := congrArg (complementSubset rows) eq
      simpa only [complementTwice rows s hs, complementTwice rows t ht] using h)
  apply List.perm_iff_count.mpr
  intro s
  rw [mapped.count, nd.count]
  have equiv : s ∈ (subsets rows).map (complementSubset rows) ↔ s ∈ subsets rows := by
    constructor
    · intro h
      obtain ⟨t, _, rfl⟩ := List.mem_map.mp h
      exact complementMember rows t
    · intro h
      exact List.mem_map.mpr ⟨complementSubset rows s, complementMember rows s,
        complementTwice rows s h⟩
  simp only [equiv]

theorem sumPerm {X : Type} (xs ys : List X) (f : X → Gaussian) (perm : xs.Perm ys) :
    scalarSum gaussianArithmetic (xs.map f) = scalarSum gaussianArithmetic (ys.map f) := by
  induction perm with
  | nil => rfl
  | cons _ _ ih => simp only [List.map_cons, scalarSum, ih]
  | swap x y xs =>
    simp only [List.map_cons, scalarSum]
    rw [← gaussianAddAssoc, gaussianAddComm (f y) (f x), gaussianAddAssoc]
  | trans _ _ ih₁ ih₂ => exact ih₁.trans ih₂

/-- Principal-minor expansion indexed by the RETAINED coordinates, matching
rows=finish in the squared-minor spectral convention. -/
theorem retainedDiagonalExpansion (rows : List V) (a : V → Gaussian) (M : V → V → Gaussian)
    (unique : rows.Nodup) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) M) rows rows =
      scalarSum gaussianArithmetic ((subsets rows).map fun s =>
        gaussianArithmetic.mul (product rows (fun x => if s x then gaussianArithmetic.one else a x))
          (minorDet gaussianArithmetic M (rows.filter s) (rows.filter s))) := by
  rw [diagonalSumExpansion rows a M unique]
  let term := fun s => gaussianArithmetic.mul (diagonalWeight rows a s)
    (minorDet gaussianArithmetic M (remaining rows s) (remaining rows s))
  have perm := complementsPerm rows unique
  have sumPerm : scalarSum gaussianArithmetic ((subsets rows).map term) =
      scalarSum gaussianArithmetic (((subsets rows).map (complementSubset rows)).map term) := by
    exact (sumPerm _ _ term perm).symm
  rw [sumPerm]
  simp only [List.map_map, Function.comp_def]
  apply gaussianSumCongr
  intro s _
  have weight : diagonalWeight rows a (complementSubset rows s) =
      product rows (fun x => if s x then gaussianArithmetic.one else a x) := by
    apply productCongr
    intro x hx
    simp only [complementSubset, hx, ↓reduceIte]
    cases s x <;> rfl
  have keep : remaining rows (complementSubset rows s) = rows.filter s := by
    apply List.filter_congr
    intro x hx
    simp [complementSubset, hx]
  change gaussianArithmetic.mul _ _ = _
  rw [weight, keep]

/-- The ordinary finite matrix product U D_b U*, with explicit actual entries. -/
def conjugatedDiagonal (vertices : List V) (U : V → V → Gaussian) (b : V → Gaussian) : V → V → Gaussian :=
  matrixMul vertices (matrixMul vertices U (diagonal b)) (adjoint U)

/-- The proved first matrix-side stage for the actual spectral pencil.
Its coefficients are principal minors of U D_b U*. Their further expansion
into squared minors of U is a separate theorem, NOT a premise here. -/
theorem spectralPrincipalExpansion (vertices : List V) (a b : V → Gaussian)
    (U : V → V → Gaussian) (unique : vertices.Nodup) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) (conjugatedDiagonal vertices U b)) vertices vertices =
      scalarSum gaussianArithmetic ((subsets vertices).map fun s =>
        gaussianArithmetic.mul (product vertices (fun x => if s x then gaussianArithmetic.one else a x))
          (minorDet gaussianArithmetic (conjugatedDiagonal vertices U b)
            (vertices.filter s) (vertices.filter s))) :=
  retainedDiagonalExpansion vertices a (conjugatedDiagonal vertices U b) unique


/-- Actual multiplication by a diagonal matrix scales each included column. -/
theorem matrixMulDiagonal (vertices : List V) (U : V → V → Gaussian) (b : V → Gaussian)
    (unique : vertices.Nodup) (i j : V) (member : j ∈ vertices) :
    matrixMul vertices U (diagonal b) i j = gaussianArithmetic.mul (U i j) (b j) := by
  unfold matrixMul
  calc
    _ = scalarSum gaussianArithmetic (vertices.map fun k =>
        if k = j then gaussianArithmetic.mul (U i j) (b j) else gaussianArithmetic.zero) := by
      apply gaussianSumCongr
      intro k _
      by_cases same : k = j
      · subst k; simp [diagonal]
      · simp [diagonal, same, gaussianArithmetic.mulZero]
    _ = _ := gaussianSumDelta vertices unique j member _

/-- Explicit entries of the actual conjugated diagonal matrix. -/
theorem conjugatedDiagonalEntry (vertices : List V) (U : V → V → Gaussian) (b : V → Gaussian)
    (unique : vertices.Nodup) (i j : V) :
    conjugatedDiagonal vertices U b i j =
      scalarSum gaussianArithmetic (vertices.map fun k =>
        gaussianArithmetic.mul (gaussianArithmetic.mul (U i k) (b k)) (conjugate (U j k))) := by
  unfold conjugatedDiagonal matrixMul
  apply gaussianSumCongr
  intro k hk
  change gaussianArithmetic.mul (matrixMul vertices U (diagonal b) i k) _ = _
  rw [matrixMulDiagonal vertices U b unique i k hk]
  rfl

/-- All diagonal spectra, including zeros, have the expected actual determinant. -/
theorem diagonalDeterminant (rows : List V) (a : V → Gaussian) (unique : rows.Nodup) :
    minorDet gaussianArithmetic (diagonal a) rows rows = product rows a := by
  induction rows with
  | nil => rfl
  | cons x rows ih =>
    obtain ⟨absent, tail⟩ := List.nodup_cons.mp unique
    rw [minorFirstColumn (diagonal a) x x rows rows]
    · simp only [diagonal, ↓reduceIte, product, ih tail]
    · intro y hy
      have ne : x ≠ y := by intro eq; subst y; exact absent hy
      simp [diagonal, ne]

/-- The empty actual spectral determinant is one, with no ambient Edge needed. -/
theorem emptySpectralDeterminant (a b : V → Gaussian) (U : V → V → Gaussian) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) (conjugatedDiagonal [] U b)) [] [] =
      gaussianArithmetic.one := rfl

/-- A direct order-one check of the square-modulus coefficient convention. -/
theorem singletonSpectralDeterminant (x : V) (a b : V → Gaussian) (U : V → V → Gaussian) :
    minorDet gaussianArithmetic (matrixAdd (diagonal a) (conjugatedDiagonal [x] U b)) [x] [x] =
      gaussianArithmetic.add (a x) (scale (normSq (U x x)) (b x)) := by
  rw [singletonMinor]
  change gaussianArithmetic.add (if x = x then a x else gaussianArithmetic.zero)
    (conjugatedDiagonal [x] U b x x) = _
  rw [conjugatedDiagonalEntry [x] U b (by simp) x x]
  simp only [↓reduceIte, List.map_cons, List.map_nil, scalarSum, gaussianAddZero]
  congr 1
  apply Prod.ext <;> simp only [gaussianArithmetic, conjugate, scale, normSq] <;> grind

#print axioms filterWithout
#print axioms diagonalRowMinor
#print axioms singleDiagonalUpdate
#print axioms partialDiagonalNil
#print axioms partialDiagonalCons
#print axioms remainingInsert
#print axioms diagonalWeightExcluded
#print axioms diagonalWeightIncluded
#print axioms partialDiagonalExpansion
#print axioms diagonalSumExpansion
#print axioms complementMember
#print axioms complementTwice
#print axioms complementsPerm
#print axioms sumPerm
#print axioms retainedDiagonalExpansion
#print axioms spectralPrincipalExpansion
#print axioms matrixMulDiagonal
#print axioms conjugatedDiagonalEntry
#print axioms diagonalDeterminant
#print axioms emptySpectralDeterminant
#print axioms singletonSpectralDeterminant

end MathCommons.ForestMatrixSumExpansion
