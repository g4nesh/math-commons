import ForestPermutationLaw

/-! The scalar finite-subset half of the forest determinant bridge.
No matrix determinant-of-sum identity is assumed or claimed here. -/
namespace MathCommons.ForestSpectralExpansion
open Classical
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory
open ForestSubsetExchange ForestPermutationLaw


def product {V : Type} (xs : List V) (f : V → Gaussian) : Gaussian :=
  match xs with
  | [] => gaussianArithmetic.one
  | x :: rest => gaussianArithmetic.mul (f x) (product rest f)

theorem productCongr {V : Type} (xs : List V) (f g : V → Gaussian)
    (same : ∀ x ∈ xs, f x = g x) : product xs f = product xs g := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [product]
    rw [same x List.mem_cons_self, ih (fun y hy => same y (List.mem_cons_of_mem x hy))]

theorem productMap {V W : Type} (xs : List V) (p : V → W) (f : W → Gaussian) :
    product (xs.map p) f = product xs (fun x => f (p x)) := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [List.map_cons, product, ih]

theorem productPerm {V : Type} (xs ys : List V) (f : V → Gaussian)
    (perm : xs.Perm ys) : product xs f = product ys f := by
  induction perm with
  | nil => rfl
  | cons _ _ ih => simp only [product, ih]
  | swap x y xs =>
    simp only [product]
    rw [← gaussianMulAssoc, gaussianMulComm (f y) (f x), gaussianMulAssoc]
  | trans _ _ ih₁ ih₂ => exact ih₁.trans ih₂

section Subsets
variable {V : Type} [DecidableEq V]

def insert (x : V) (s : Subset V) : Subset V := fun y => if y = x then true else s y

/-- All Boolean subsets supported in the supplied finite list. -/
def subsets : List V → List (Subset V)
  | [] => [fun _ => false]
  | x :: xs => subsets xs ++ (subsets xs).map (insert x)

theorem subsetOutside (xs : List V) (s : Subset V) (member : s ∈ subsets xs)
    (x : V) (outside : x ∉ xs) : s x = false := by
  induction xs generalizing s with
  | nil =>
    simp only [subsets, List.mem_cons, List.not_mem_nil, or_false] at member
    subst s; rfl
  | cons y ys ih =>
    have ne : x ≠ y := by intro eq; exact outside (eq ▸ List.mem_cons_self)
    have absent : x ∉ ys := by intro h; exact outside (List.mem_cons_of_mem y h)
    simp only [subsets, List.mem_append, List.mem_map] at member
    rcases member with old | ⟨t, ht, rfl⟩
    · exact ih s old absent
    · simp only [insert, ne, ↓reduceIte]
      exact ih t ht absent

theorem subsetsNodup (xs : List V) (unique : xs.Nodup) : (subsets xs).Nodup := by
  induction xs with
  | nil => simp [subsets]
  | cons x xs ih =>
    obtain ⟨absent, tail⟩ := List.nodup_cons.mp unique
    apply List.nodup_append.mpr
    refine ⟨ih tail, nodupMapOn _ _ (ih tail) ?_, ?_⟩
    · intro s hs t ht same
      funext y
      by_cases hit : y = x
      · subst y; rw [subsetOutside xs s hs x absent, subsetOutside xs t ht x absent]
      · have h := congrFun same y
        simpa only [insert, hit, ↓reduceIte] using h
    · intro s hs t ht same
      obtain ⟨r, hr, eq⟩ := List.mem_map.mp ht
      have h := congrFun (same.trans eq.symm) x
      have off := subsetOutside xs s hs x absent
      simp only [insert, ↓reduceIte, off] at h
      cases h

theorem subsetMember (xs : List V) (s : Subset V)
    (supported : ∀ x, x ∉ xs → s x = false) : s ∈ subsets xs := by
  induction xs generalizing s with
  | nil =>
    have eq : s = fun _ => false := by funext x; exact supported x (by simp)
    simp [subsets, eq]
  | cons x xs ih =>
    let t : Subset V := fun y => if y = x then false else s y
    have tm : t ∈ subsets xs := by
      apply ih
      intro y outside
      by_cases hit : y = x
      · simp [t, hit]
      · simp only [t, hit, ↓reduceIte]
        exact supported y (by simp [hit, outside])
    cases hx : s x
    · have eq : t = s := by funext y; by_cases hit : y = x <;> simp [t, hit, hx]
      exact List.mem_append_left _ (eq ▸ tm)
    · apply List.mem_append_right
      apply List.mem_map.mpr
      refine ⟨t, tm, ?_⟩
      funext y
      by_cases hit : y = x <;> simp [insert, t, hit, hx]

theorem subsetsComplete (xs : List V) (complete : ∀ x, x ∈ xs) (s : Subset V) :
    s ∈ subsets xs := subsetMember xs s (fun x outside => False.elim (outside (complete x)))

def selectionTerm (xs : List V) (s : Subset V) (a b : V → Gaussian) : Gaussian :=
  product xs (fun x => if s x then b x else a x)

/-- Distributivity over every finite subset, proved by list induction. -/
theorem productOfSums (xs : List V) (unique : xs.Nodup) (a b : V → Gaussian) :
    product xs (fun x => gaussianArithmetic.add (a x) (b x)) =
      scalarSum gaussianArithmetic ((subsets xs).map fun s => selectionTerm xs s a b) := by
  induction xs with
  | nil => simp only [product, subsets, List.map_cons, List.map_nil, scalarSum,
      selectionTerm, product, gaussianAddZero]
  | cons x xs ih =>
    obtain ⟨absent, tail⟩ := List.nodup_cons.mp unique
    rw [product, ih tail, gaussianAddMul]
    simp only [subsets, gaussianSumAppend, List.map_map, Function.comp_def]
    have excluded : scalarSum gaussianArithmetic ((subsets xs).map fun s => selectionTerm (x::xs) s a b) =
        gaussianArithmetic.mul (a x) (scalarSum gaussianArithmetic ((subsets xs).map fun s => selectionTerm xs s a b)) := by
      rw [← gaussianSumScale]
      apply gaussianSumCongr
      intro s hs
      simp only [selectionTerm, product, subsetOutside xs s hs x absent, Bool.false_eq_true, ↓reduceIte]
    have included : scalarSum gaussianArithmetic ((subsets xs).map fun s => selectionTerm (x::xs) (insert x s) a b) =
        gaussianArithmetic.mul (b x) (scalarSum gaussianArithmetic ((subsets xs).map fun s => selectionTerm xs s a b)) := by
      rw [← gaussianSumScale]
      apply gaussianSumCongr
      intro s _
      simp only [selectionTerm, product, insert, ↓reduceIte]
      congr 1
      apply productCongr
      intro y hy
      have ne : y ≠ x := by intro eq; exact absent (eq ▸ hy)
      simp only [ne, ↓reduceIte]
    rw [excluded, included]

omit [DecidableEq V] in
theorem productPointwise (xs : List V) (a b : V → Gaussian) :
    product xs (fun x => gaussianArithmetic.mul (a x) (b x)) =
      gaussianArithmetic.mul (product xs a) (product xs b) := by
  induction xs with
  | nil => exact (gaussianOneMul gaussianArithmetic.one).symm
  | cons x xs ih =>
    simp only [product, ih]
    apply Prod.ext <;> simp only [gaussianArithmetic] <;> grind

omit [DecidableEq V] in
theorem selectionSplit (xs : List V) (s : Subset V) (a b : V → Gaussian) :
    selectionTerm xs s a b = gaussianArithmetic.mul
      (product xs (fun x => if s x then gaussianArithmetic.one else a x))
      (product xs (fun x => if s x then b x else gaussianArithmetic.one)) := by
  rw [← productPointwise]
  apply productCongr
  intro x _
  cases s x <;> simp only [Bool.false_eq_true, ↓reduceIte, gaussianMulOne, gaussianOneMul]

theorem permutationList (xs : List V) (unique : xs.Nodup)
    (complete : ∀ x, x ∈ xs) (p : VertexPermutation V) : (xs.map p.forward).Perm xs := by
  have mapped := nodupMapOn xs p.forward unique
    (fun x _ y _ eq => permutationInjective p x y eq)
  apply List.perm_iff_count.mpr
  intro x
  have member : x ∈ xs.map p.forward :=
    List.mem_map.mpr ⟨p.inverse x, complete _, p.forward_inverse x⟩
  rw [mapped.count, unique.count]
  simp only [member, complete x, ↓reduceIte]

theorem productReindex (xs : List V) (unique : xs.Nodup)
    (complete : ∀ x, x ∈ xs) (p : VertexPermutation V) (a : V → Gaussian) :
    product xs (fun x => a (p.forward x)) = product xs a := by
  rw [← productMap]
  exact productPerm _ _ a (permutationList xs unique complete p)

/-- A spectral coefficient with columns `start` and rows `finish`. -/
def spectralMonomial (xs : List V) (start finish : Subset V) (a b : V → Gaussian) : Gaussian :=
  gaussianArithmetic.mul
    (product xs (fun x => if finish x then gaussianArithmetic.one else a x))
    (product xs (fun x => if start x then b x else gaussianArithmetic.one))

theorem selectionAsMonomial (xs : List V) (unique : xs.Nodup)
    (complete : ∀ x, x ∈ xs) (p : VertexPermutation V) (s : Subset V)
    (a b : V → Gaussian) :
    selectionTerm xs s (fun x => a (p.forward x)) b =
      spectralMonomial xs s (subsetImage p s) a b := by
  rw [selectionSplit]
  unfold spectralMonomial
  congr 1
  calc
    product xs (fun x => if s x then gaussianArithmetic.one else a (p.forward x)) =
        product xs (fun x => if subsetImage p s (p.forward x) then gaussianArithmetic.one else a (p.forward x)) := by
      apply productCongr
      intro x _
      rw [subsetImageForward]
    _ = _ := productReindex xs unique complete p (fun x => if subsetImage p s x then gaussianArithmetic.one else a x)

/-- Pairings use column-to-row permutations, matching the actual minor law. -/
def pairingProduct (xs : List V) (p : VertexPermutation V) (a b : V → Gaussian) : Gaussian :=
  product xs (fun column => gaussianArithmetic.add (a (p.forward column)) (b column))

theorem pairingProductAsRows (xs : List V) (unique : xs.Nodup)
    (complete : ∀ x, x ∈ xs) (p : VertexPermutation V) (a b : V → Gaussian) :
    pairingProduct xs p a b =
      product xs (fun row => gaussianArithmetic.add (a row) (b (p.inverse row))) := by
  rw [← productReindex xs unique complete p]
  apply productCongr
  intro x _
  rw [p.inverse_forward]

theorem pairingSubsetExpansion (xs : List V) (unique : xs.Nodup)
    (complete : ∀ x, x ∈ xs) (p : VertexPermutation V) (a b : V → Gaussian) :
    pairingProduct xs p a b = scalarSum gaussianArithmetic
      ((subsets xs).map fun s => spectralMonomial xs s (subsetImage p s) a b) := by
  rw [pairingProduct, productOfSums xs unique]
  apply gaussianSumCongr
  intro s _
  exact selectionAsMonomial xs unique complete p s a b

end Subsets

def scale (q : Rat) (z : Gaussian) : Gaussian := (q * z.1, q * z.2)

theorem scaleZero (z : Gaussian) : scale 0 z = gaussianArithmetic.zero := by
  apply Prod.ext <;> simp only [scale, gaussianArithmetic, Rat.zero_mul]

theorem scaleAdd (q : Rat) (z w : Gaussian) :
    scale q (gaussianArithmetic.add z w) =
      gaussianArithmetic.add (scale q z) (scale q w) := by
  apply Prod.ext <;> simp only [scale, gaussianArithmetic] <;> grind

theorem scaleSum {X : Type} (xs : List X) (q : Rat) (f : X → Gaussian) :
    scale q (scalarSum gaussianArithmetic (xs.map f)) =
      scalarSum gaussianArithmetic (xs.map fun x => scale q (f x)) := by
  induction xs with
  | nil => apply Prod.ext <;> simp only [List.map_nil, scalarSum, scale, gaussianArithmetic, Rat.mul_zero]
  | cons x xs ih => simp only [List.map_cons, scalarSum, scaleAdd, ih]

theorem sumScale {X : Type} (xs : List X) (weight : X → Rat) (z : Gaussian) :
    scalarSum gaussianArithmetic (xs.map fun x => scale (weight x) z) =
      scale (sumWith weight xs) z := by
  induction xs with
  | nil => exact (scaleZero z).symm
  | cons x xs ih =>
    simp only [List.map_cons, scalarSum, sumWith, ih]
    apply Prod.ext <;> simp only [scale, gaussianArithmetic] <;> grind

/-- Regroup a finite weighted sum by a complete, duplicate-free set of keys. -/
theorem weightedGrouping {X K : Type} (xs : List X) (keys : List K)
    (key : X → K) (weight : X → Rat) (f : K → Gaussian)
    (unique : keys.Nodup) (complete : ∀ x ∈ xs, key x ∈ keys) :
    scalarSum gaussianArithmetic (xs.map fun x => scale (weight x) (f (key x))) =
      scalarSum gaussianArithmetic (keys.map fun k =>
        scale (sumWith (fun x => if key x = k then weight x else 0) xs) (f k)) := by
  have expand : scalarSum gaussianArithmetic (xs.map fun x => scale (weight x) (f (key x))) =
      scalarSum gaussianArithmetic (xs.map fun x => scalarSum gaussianArithmetic
        (keys.map fun k => if key x = k then scale (weight x) (f k) else gaussianArithmetic.zero)) := by
    apply gaussianSumCongr
    intro x hx
    symm
    have delta := scalarSumUniqueKey keys (fun k => k) (by simpa using unique)
      (key x) (complete x hx)
      (fun k => if key x = k then scale (weight x) (f k) else gaussianArithmetic.zero)
      (by intro k _ ne; simp only [Ne.symm ne, ↓reduceIte])
    simpa only [↓reduceIte] using delta
  rw [expand, gaussianSumSwap]
  apply gaussianSumCongr
  intro k _
  rw [← sumScale]
  apply gaussianSumCongr
  intro x _
  split
  · rfl
  · exact (scaleZero (f k)).symm

section Weighted
variable {V : Type} [DecidableEq V]

def weightedPairings (xs : List V) (law : List (VertexPermutation V × Rat))
    (a b : V → Gaussian) : Gaussian :=
  scalarSum gaussianArithmetic (law.map fun entry => scale entry.2 (pairingProduct xs entry.1 a b))

/-- Finite scalar expansion of weighted permutation products. No determinant
formula or positivity assumption is used. -/
theorem weightedSubsetExpansion (xs : List V) (unique : xs.Nodup)
    (complete : ∀ x, x ∈ xs) (law : List (VertexPermutation V × Rat))
    (a b : V → Gaussian) :
    weightedPairings xs law a b = scalarSum gaussianArithmetic
      ((subsets xs).map fun start => scalarSum gaussianArithmetic
        ((subsets xs).map fun finish => scale
          (sumWith (fun entry => if subsetImage entry.1 start = finish then entry.2 else 0) law)
          (spectralMonomial xs start finish a b))) := by
  unfold weightedPairings
  have expand : scalarSum gaussianArithmetic (law.map fun entry => scale entry.2 (pairingProduct xs entry.1 a b)) =
      scalarSum gaussianArithmetic (law.map fun entry => scalarSum gaussianArithmetic
        ((subsets xs).map fun start => scale entry.2
          (spectralMonomial xs start (subsetImage entry.1 start) a b))) := by
    apply gaussianSumCongr
    intro entry _
    rw [pairingSubsetExpansion xs unique complete, scaleSum]
  rw [expand, gaussianSumSwap]
  apply gaussianSumCongr
  intro start _
  exact weightedGrouping law (subsets xs) (fun entry => subsetImage entry.1 start)
    (fun entry => entry.2) (fun finish => spectralMonomial xs start finish a b)
    (subsetsNodup xs unique) (fun _ _ => subsetsComplete xs complete _)

omit [DecidableEq V] in
theorem productOnes (xs : List V) :
    product xs (fun _ => gaussianArithmetic.one) = gaussianArithmetic.one := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [product, ih, gaussianOneMul]

omit [DecidableEq V] in
theorem emptySpectralMonomial (xs : List V) (a b : V → Gaussian) :
    spectralMonomial xs (fun _ => false) (fun _ => false) a b = product xs a := by
  simp only [spectralMonomial, Bool.false_eq_true, ↓reduceIte, productOnes, gaussianMulOne]

omit [DecidableEq V] in
theorem fullSpectralMonomial (xs : List V) (a b : V → Gaussian) :
    spectralMonomial xs (fun _ => true) (fun _ => true) a b = product xs b := by
  simp only [spectralMonomial, ↓reduceIte, productOnes, gaussianOneMul]

/-- Concrete forest-law scalar expansion. Coefficients are the squared norms
of actual gate-product minors, with rows=finish and columns=start. The same
nonnegative normalized permutation law works for all Gaussian spectra.
The determinant of a matrix sum does not occur in this theorem. -/
theorem forestWeightedSubsetExpansion {vertices : List V} (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k))
    (a b : V → Gaussian) :
    (∀ entry ∈ coinLaw edges bs n, 0 ≤ entry.2) ∧
    sumWith (fun entry => entry.2) (coinLaw edges bs n) = 1 ∧
    weightedPairings vertices (coinLaw edges bs n) a b = scalarSum gaussianArithmetic
      ((subsets vertices).map fun start => scalarSum gaussianArithmetic
        ((subsets vertices).map fun finish => scale
          (normSq (minorDet gaussianArithmetic
            (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
            (vertices.filter finish) (vertices.filter start)))
          (spectralMonomial vertices start finish a b))) := by
  have law := commonSquaredMinorLaw edges as bs cs ds n bridges unitary
  refine ⟨law.1, law.2.1, ?_⟩
  rw [weightedSubsetExpansion vertices (edges 0).nodup (edges 0).complete]
  apply gaussianSumCongr
  intro start _
  apply gaussianSumCongr
  intro finish _
  rw [law.2.2 start finish]

end Weighted

#print axioms productCongr
#print axioms productMap
#print axioms productPerm
#print axioms subsetOutside
#print axioms subsetsNodup
#print axioms subsetMember
#print axioms subsetsComplete
#print axioms productOfSums
#print axioms productPointwise
#print axioms selectionSplit
#print axioms permutationList
#print axioms productReindex
#print axioms selectionAsMonomial
#print axioms pairingProductAsRows
#print axioms pairingSubsetExpansion
#print axioms scaleZero
#print axioms scaleAdd
#print axioms scaleSum
#print axioms sumScale
#print axioms weightedGrouping
#print axioms weightedSubsetExpansion
#print axioms productOnes
#print axioms emptySpectralMonomial
#print axioms fullSpectralMonomial
#print axioms forestWeightedSubsetExpansion
end MathCommons.ForestSpectralExpansion
