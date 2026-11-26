import ForestProductMinor

/-!
Actual order-two minors of an arbitrary-length two-coordinate gate product are
identified with the canonical legal-history endpoint amplitudes. The finite
ordered-pair/subset correspondence, actual identity-minor base case, and all
sum exchanges are proved here. For a forest schedule of unitary Gaussian-
rational blocks, the actual squared minors equal canonical endpoint probabilities
and sum to one across all canonical output pairs. No desired transport equality
is assumed. General minor orders and standard complex matrices remain outside
this module's scope.
-/

namespace MathCommons.ForestPairHistory
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor

section Pairs
variable {V : Type} [DecidableEq V]

def pairSubset (p : V × V) : Subset V := fun x => decide (x = p.1 ∨ x = p.2)

omit [DecidableEq V] in
theorem pairMemberSublist (vertices : List V) (p : V × V) :
    p ∈ orderedPairs vertices ↔ List.Sublist [p.1, p.2] vertices := by
  induction vertices with
  | nil => simp [orderedPairs]
  | cons x xs ih =>
    simp only [orderedPairs, List.mem_append, List.mem_map]
    constructor
    · intro h
      rcases h with ⟨y, hy, eq⟩ | h
      · cases eq
        exact List.Sublist.cons₂ x (List.singleton_sublist.mpr hy)
      · exact List.Sublist.cons x (ih.mp h)
    · intro h
      cases h with
      | cons _ h => exact Or.inr (ih.mpr h)
      | cons₂ _ h =>
        exact Or.inl ⟨p.2, List.singleton_sublist.mp h, by cases p; rfl⟩

theorem filterSingle (vertices : List V) (unique : vertices.Nodup) (x : V)
    (member : x ∈ vertices) : vertices.filter (fun y => decide (y = x)) = [x] := by
  induction vertices with
  | nil => simp at member
  | cons y ys ih =>
    obtain ⟨absent, tailUnique⟩ := List.nodup_cons.mp unique
    by_cases equal : y = x
    · subst y
      have zero : ys.filter (fun y => decide (y = x)) = [] := by
        apply List.filter_eq_nil_iff.mpr
        intro z hz
        simp only [decide_eq_true_eq]
        intro eq
        subst z
        exact absent hz
      simp [zero]
    · have tail : x ∈ ys := by simpa [Ne.symm equal] using member
      simp [equal, ih tailUnique tail]

theorem pairFilter (vertices : List V) (unique : vertices.Nodup) (p : V × V)
    (member : p ∈ orderedPairs vertices) : vertices.filter (pairSubset p) = [p.1, p.2] := by
  induction vertices with
  | nil => simp [orderedPairs] at member
  | cons x xs ih =>
    obtain ⟨absent, tailUnique⟩ := List.nodup_cons.mp unique
    simp only [orderedPairs, List.mem_append, List.mem_map] at member
    rcases member with ⟨y, hy, eq⟩ | member
    · cases eq
      have tail : xs.filter (pairSubset (x,y)) = xs.filter (fun z => decide (z=y)) := by
        apply List.filter_congr
        intro z hz
        have ne : z ≠ x := by intro eq; subst z; exact absent hz
        simp [pairSubset, ne]
      simp [pairSubset, tail, filterSingle xs tailUnique y hy]
    · have sub := (pairMemberSublist xs p).mp member
      have first : p.1 ∈ xs := sub.subset (by simp)
      have second : p.2 ∈ xs := sub.subset (by simp)
      have ne1 : x ≠ p.1 := by intro eq; subst x; exact absent first
      have ne2 : x ≠ p.2 := by intro eq; subst x; exact absent second
      simpa [List.filter_cons, pairSubset, ne1, ne2] using ih tailUnique member


theorem subsetOfFilter (vertices : List V) (complete : ∀ x, x ∈ vertices)
    (s : Subset V) (p : V × V) (filtered : vertices.filter s = [p.1,p.2]) :
    s = pairSubset p := by
  funext x
  have member : s x = true ↔ x = p.1 ∨ x = p.2 := by
    have h : x ∈ vertices.filter s ↔ x ∈ [p.1,p.2] := by rw [filtered]
    simpa [List.mem_filter, complete x] using h
  simp only [pairSubset]
  cases h : s x <;> simp_all

theorem pairSubsetInjective (vertices : List V) (unique : vertices.Nodup)
    (p q : V × V) (mp : p ∈ orderedPairs vertices) (mq : q ∈ orderedPairs vertices)
    (same : pairSubset p = pairSubset q) : p = q := by
  have a := pairFilter vertices unique p mp
  have b := pairFilter vertices unique q mq
  rw [same, b] at a
  cases p; cases q
  simpa using a.symm

omit [DecidableEq V] in
theorem orderedPairsNodup (vertices : List V) (unique : vertices.Nodup) :
    (orderedPairs vertices).Nodup := by
  induction vertices with
  | nil => simp [orderedPairs]
  | cons x xs ih =>
    obtain ⟨absent, tailUnique⟩ := List.nodup_cons.mp unique
    apply List.nodup_append.mpr
    refine ⟨?_, ih tailUnique, ?_⟩
    · apply List.pairwise_map.mpr
      exact tailUnique.imp (fun h eq => h (congrArg Prod.snd eq))
    · intro p mp q mq equal
      obtain ⟨y, hy, peq⟩ := List.mem_map.mp mp
      have qfirst := ((pairMemberSublist xs q).mp mq).subset (by simp : q.1 ∈ [q.1,q.2])
      have xeq : x = q.1 := by rw [← equal, ← peq]
      exact absent (xeq ▸ qfirst)

omit [DecidableEq V] in
theorem nodupMapOn {α β : Type} (xs : List α) (f : α → β) (unique : xs.Nodup)
    (injective : ∀ x ∈ xs, ∀ y ∈ xs, f x = f y → x = y) : (xs.map f).Nodup := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    obtain ⟨absent, tailUnique⟩ := List.nodup_cons.mp unique
    apply List.nodup_cons.mpr
    constructor
    · intro member
      obtain ⟨y, hy, eq⟩ := List.mem_map.mp member
      have same := injective y (List.mem_cons_of_mem x hy) x List.mem_cons_self eq
      exact absent (same ▸ hy)
    · exact ih tailUnique (fun y hy z hz => injective y (List.mem_cons_of_mem x hy)
        z (List.mem_cons_of_mem x hz))

theorem pairSubsetsNodup (vertices : List V) (unique : vertices.Nodup) :
    ((orderedPairs vertices).map pairSubset).Nodup :=
  nodupMapOn _ _ (orderedPairsNodup vertices unique) (fun p mp q mq => pairSubsetInjective vertices unique p q mp mq)

theorem twoSubsetRepresentative (vertices : List V) (complete : ∀ x, x ∈ vertices)
    (s : Subset V) (size : (vertices.filter s).length = 2) :
    ∃ p ∈ orderedPairs vertices, vertices.filter s = [p.1,p.2] ∧ s = pairSubset p := by
  cases filtered : vertices.filter s with
  | nil => simp [filtered] at size
  | cons x xs =>
    cases xs with
    | nil => simp [filtered] at size
    | cons y ys =>
      have empty : ys = [] := by simpa [filtered] using size
      subst ys
      refine ⟨(x,y), ?_, rfl, subsetOfFilter vertices complete s (x,y) filtered⟩
      apply (pairMemberSublist vertices (x,y)).mpr
      have h : List.Sublist (vertices.filter s) vertices := List.filter_sublist
      simpa [filtered] using h


open Classical in
theorem identitySubsetMinor (vertices : List V) (unique : vertices.Nodup)
    (complete : ∀ x, x ∈ vertices) (s t : Subset V) :
    minorDet gaussianArithmetic identityMatrix (vertices.filter t) (vertices.filter s) =
      if t = s then gaussianArithmetic.one else gaussianArithmetic.zero := by
  by_cases same : t = s
  · subst t
    simp only [↓reduceIte]
    have h := identityPrefixMinor (identityMatrix (V := V)) (vertices.filter s) [] []
      (by simpa using unique.filter s) (fun x _ y => rfl)
    simpa [minorDet] using h
  · simp only [same, ↓reduceIte]
    by_cases zero : minorDet gaussianArithmetic identityMatrix
        (vertices.filter t) (vertices.filter s) = gaussianArithmetic.zero
    · exact zero
    · have matching := nonzeroMinorMatching gaussianArithmetic identityMatrix _ _ zero
      have rows := matchingRowSupport matching
      have cols := matchingColumnSupport matching
      have equal : t = s := by
        funext x
        have iff : t x = true ↔ s x = true := by
          constructor
          · intro ht
            obtain ⟨y, hy, entry⟩ := rows x (List.mem_filter.mpr ⟨complete x, ht⟩)
            have eq : x = y := by
              by_cases h : x = y
              · exact h
              · exact False.elim (entry (by simp [identityMatrix, h]))
            subst y
            exact (List.mem_filter.mp hy).2
          · intro hs
            obtain ⟨y, hy, entry⟩ := cols x (List.mem_filter.mpr ⟨complete x, hs⟩)
            have eq : y = x := by
              by_cases h : y = x
              · exact h
              · exact False.elim (entry (by simp [identityMatrix, h]))
            subst y
            exact (List.mem_filter.mp hy).2
        cases ht : t x <;> cases hs : s x <;> simp_all
      exact False.elim (same equal)

end Pairs

section Histories
variable {V : Type} [DecidableEq V] {vertices : List V}
open Classical

theorem gaussianSumFilter {α : Type} (xs : List α) (test : α → Bool) (f : α → Gaussian) :
    scalarSum gaussianArithmetic ((xs.filter test).map f) =
      scalarSum gaussianArithmetic (xs.map fun x => if test x then f x else gaussianArithmetic.zero) := by
  induction xs with
  | nil => rfl
  | cons x xs ih => cases h : test x <;> simp [h, ih, scalarSum, gaussianZeroAdd]

theorem endpointAsSum (edges : Nat → Edge vertices) (as bs cs ds : Nat → Gaussian)
    (start finish : Subset V) (n : Nat) :
    endpointAmplitude edges as bs cs ds start finish n =
      scalarSum gaussianArithmetic ((histories edges start n).map fun h =>
        if h n = finish then pathAmplitude (vertices := vertices) gaussianArithmetic
          (gateSchedule edges as bs cs ds) h n else gaussianArithmetic.zero) := by
  unfold endpointAmplitude endpointHistories
  rw [gaussianSumFilter]
  apply gaussianSumCongr
  intro h _
  have iff : vertices.all (fun x => h n x == finish x) = true ↔ h n = finish := by
    constructor
    · intro same
      funext x
      simpa only [beq_iff_eq] using List.all_eq_true.mp same x ((edges 0).complete x)
    · intro same
      simp [same]
  simp only [iff]

omit [DecidableEq V] in
theorem pathAmplitudePrefixCongr (matrices : Nat → V → V → Gaussian)
    (x y : Nat → Subset V) (n : Nat) (same : ∀ k, k ≤ n → x k = y k) :
    pathAmplitude (vertices := vertices) gaussianArithmetic matrices x n =
      pathAmplitude (vertices := vertices) gaussianArithmetic matrices y n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [pathAmplitude]
    rw [ih (fun k hk => same k (by omega)), same n (by omega), same (n+1) (by omega)]

omit [DecidableEq V] in
theorem extendedPathAmplitude (matrices : Nat → V → V → Gaussian)
    (h : Nat → Subset V) (n : Nat) (next : Subset V) :
    pathAmplitude (vertices := vertices) gaussianArithmetic matrices (extendHistory h n next) (n+1) =
      gaussianArithmetic.mul (pathAmplitude (vertices := vertices) gaussianArithmetic matrices h n)
        (minorDet gaussianArithmetic (matrices n) (vertices.filter next) (vertices.filter (h n))) := by
  simp only [pathAmplitude]
  have samePrefix : ∀ k, k ≤ n → extendHistory h n next k = h k := by
    intro k hk
    simp [extendHistory, show k ≠ n+1 by omega]
  rw [pathAmplitudePrefixCongr matrices _ h n samePrefix]
  simp [extendHistory]


theorem successorMinorSum (e : Edge vertices) (a b c d : Gaussian) (s t : Subset V)
    (weight : Gaussian) :
    scalarSum gaussianArithmetic ((successors e s).map fun next =>
      if next = t then gaussianArithmetic.mul weight
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter next) (vertices.filter s)) else gaussianArithmetic.zero) =
      gaussianArithmetic.mul weight
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter s)) := by
  by_cases legal : Step e s t
  · have member := everyLegalSuccessor e s t legal
    have one := scalarSumUniqueKey (successors e s) (fun x => x)
      (by simpa using successorsNodup e s) t member
      (fun next => if next = t then gaussianArithmetic.mul weight
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter next) (vertices.filter s)) else gaussianArithmetic.zero)
      (by intro next _ ne; simp [ne])
    simpa using one
  · have zero : minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
        (vertices.filter t) (vertices.filter s) = gaussianArithmetic.zero := by
      by_cases h : minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter s) = gaussianArithmetic.zero
      · exact h
      · exact False.elim (legal (gaussianGateMinorStep e a b c d s t h))
    rw [zero, gaussianArithmetic.mulZero]
    apply scalarSumMapZero
    intro next member
    have ne : next ≠ t := by
      intro eq
      subst next
      exact legal (successorLegal e s t member)
    simp [ne]

theorem endpointAmplitudeRecurrence (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat) :
    endpointAmplitude edges as bs cs ds start finish (n+1) =
      scalarSum gaussianArithmetic ((histories edges start n).map fun h =>
        gaussianArithmetic.mul
          (pathAmplitude (vertices := vertices) gaussianArithmetic (gateSchedule edges as bs cs ds) h n)
          (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
            (vertices.filter finish) (vertices.filter (h n)))) := by
  rw [endpointAsSum]
  simp only [histories, gaussianSumFlatMap, List.map_map]
  apply gaussianSumCongr
  intro h _
  simp only [Function.comp_def, extendedPathAmplitude]
  simp only [extendHistory, ↓reduceIte]
  exact successorMinorSum (edges n) (as n) (bs n) (cs n) (ds n) (h n) finish _

theorem historySizeTwo (edges : Nat → Edge vertices) (start : Subset V)
    (size : (vertices.filter start).length = 2) (n : Nat) (h : Nat → Subset V)
    (member : h ∈ histories edges start n) : (vertices.filter (h n)).length = 2 := by
  obtain ⟨initial, legal⟩ := historyLegal edges start n h member
  have count := trajectoryPreservesCardinality edges n h legal n (by omega)
  simp only [cutCount, Bool.true_and, List.countP_eq_length_filter] at count
  rw [initial, size] at count
  exact count


/-- Actual order-two gate-product minors equal the canonical endpoint amplitude.
No unitarity or forest hypothesis is needed for this signed algebraic identity. -/
theorem productPairEndpointAmplitude (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat)
    (startSize : (vertices.filter start).length = 2)
    (finishSize : (vertices.filter finish).length = 2) :
    minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start) =
        endpointAmplitude edges as bs cs ds start finish n := by
  induction n generalizing finish with
  | zero =>
    rw [matrixProduct, identitySubsetMinor vertices (edges 0).nodup (edges 0).complete]
    rw [endpointAsSum]
    simp only [histories, List.map_cons, List.map_nil, scalarSum, pathAmplitude, gaussianAddZero]
    by_cases eq : finish = start
    · simp [eq]
    · simp [eq, Ne.symm eq]
  | succ n ih =>
    obtain ⟨p, _, startFilter, _⟩ := twoSubsetRepresentative vertices (edges 0).complete start startSize
    obtain ⟨q, _, finishFilter, _⟩ := twoSubsetRepresentative vertices (edges 0).complete finish finishSize
    have cb := gateProductTwoMinorRecurrence edges as bs cs ds n q.1 q.2 p.1 p.2
    rw [← finishFilter, ← startFilter] at cb
    rw [cb]
    have replace : scalarSum gaussianArithmetic ((orderedPairs vertices).map fun pair =>
        gaussianArithmetic.mul
          (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
            (vertices.filter finish) [pair.1,pair.2])
          (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
            [pair.1,pair.2] (vertices.filter start))) =
        scalarSum gaussianArithmetic ((orderedPairs vertices).map fun pair =>
          scalarSum gaussianArithmetic ((histories edges start n).map fun h =>
            gaussianArithmetic.mul
              (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
                (vertices.filter finish) (vertices.filter (pairSubset pair)))
              (if h n = pairSubset pair then
                pathAmplitude (vertices := vertices) gaussianArithmetic (gateSchedule edges as bs cs ds) h n
               else gaussianArithmetic.zero))) := by
      apply gaussianSumCongr
      intro pair member
      have filtered := pairFilter vertices (edges 0).nodup pair member
      rw [← filtered, ih (pairSubset pair) (by rw [filtered]; rfl), endpointAsSum]
      exact (gaussianSumScale _ _ _).symm
    rw [replace, gaussianSumSwap]
    rw [endpointAmplitudeRecurrence]
    apply gaussianSumCongr
    intro h member
    obtain ⟨pair, pairMember, pairFiltered, same⟩ := twoSubsetRepresentative vertices (edges 0).complete
      (h n) (historySizeTwo edges start startSize n h member)
    have one := scalarSumUniqueKey (orderedPairs vertices) pairSubset
      (pairSubsetsNodup vertices (edges 0).nodup) pair pairMember
      (fun p => gaussianArithmetic.mul
        (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
          (vertices.filter finish) (vertices.filter (pairSubset p)))
        (if h n = pairSubset p then
          pathAmplitude (vertices := vertices) gaussianArithmetic (gateSchedule edges as bs cs ds) h n
         else gaussianArithmetic.zero)) ?_
    · rw [one]
      simp only [same, ↓reduceIte]
      exact gaussianMulComm _ _
    · intro other _ different
      have ne : h n ≠ pairSubset other := by
        intro eq
        exact different (eq.symm.trans same)
      simp only [ne, ↓reduceIte, gaussianArithmetic.mulZero]

/-- The squared ACTUAL order-two minor is the canonical probability mass for a
forest gate schedule; includes zero length and empty endpoint fibers. -/
theorem productPairMagnitude (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat)
    (startSize : (vertices.filter start).length = 2)
    (finishSize : (vertices.filter finish).length = 2)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    normSq (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start)) =
        endpointProbability edges bs start finish n := by
  rw [productPairEndpointAmplitude edges as bs cs ds start finish n startSize finishSize]
  exact endpointMagnitude edges as bs cs ds start finish n bridges unitary


/-- An empty canonical endpoint fiber gives an actual zero minor, without
unitarity or forest assumptions. -/
theorem productPairEmptyFiber (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat)
    (startSize : (vertices.filter start).length = 2)
    (finishSize : (vertices.filter finish).length = 2)
    (empty : endpointHistories edges start finish n = []) :
    minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start) = gaussianArithmetic.zero := by
  rw [productPairEndpointAmplitude edges as bs cs ds start finish n startSize finishSize]
  exact (emptyEndpointFiber edges as bs cs ds start finish n empty).1

end Histories

section Mass
open Classical

theorem sumWithCongrOn {α : Type} (xs : List α) (f g : α → Rat)
    (same : ∀ x ∈ xs, f x = g x) : sumWith f xs = sumWith g xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [sumWith]
    rw [same x List.mem_cons_self, ih (fun y hy => same y (List.mem_cons_of_mem x hy))]

theorem sumWithFilter {α : Type} (xs : List α) (test : α → Bool) (f : α → Rat) :
    sumWith f (xs.filter test) = sumWith (fun x => if test x then f x else 0) xs := by
  induction xs with
  | nil => rfl
  | cons x xs ih => cases h : test x <;> simp [h, ih, sumWith, Rat.zero_add]

theorem sumWithAdd {α : Type} (xs : List α) (f g : α → Rat) :
    sumWith (fun x => f x + g x) xs = sumWith f xs + sumWith g xs := by
  induction xs with
  | nil => simp [sumWith, Rat.zero_add]
  | cons x xs ih => simp only [sumWith, ih]; grind

theorem sumWithSwap {α β : Type} (xs : List α) (ys : List β) (f : α → β → Rat) :
    sumWith (fun x => sumWith (f x) ys) xs = sumWith (fun y => sumWith (fun x => f x y) xs) ys := by
  induction xs with
  | nil =>
    simp only [sumWith]
    exact (sumWithZero (fun _ => 0) ys (by intros; rfl)).symm
  | cons x xs ih => simp only [sumWith, ih, sumWithAdd]

theorem sumWithUniqueKey {α β : Type} (xs : List α) (key : α → β)
    (unique : (xs.map key).Nodup) (q : α) (member : q ∈ xs) (f : α → Rat)
    (zero : ∀ x ∈ xs, key x ≠ key q → f x = 0) : sumWith f xs = f q := by
  induction xs with
  | nil => simp at member
  | cons x xs ih =>
    obtain ⟨absent, tailUnique⟩ := List.nodup_cons.mp unique
    rcases List.mem_cons.mp member with rfl | member
    · have tailZero : sumWith f xs = 0 := by
        apply sumWithZero
        intro y hy
        apply zero y (List.mem_cons_of_mem q hy)
        intro eq
        apply absent
        exact List.mem_map.mpr ⟨y, hy, eq⟩
      simp only [sumWith, tailZero, Rat.add_zero]
    · have different : key x ≠ key q := by
        intro eq
        apply absent
        exact List.mem_map.mpr ⟨q, member, eq.symm⟩
      rw [sumWith, zero x List.mem_cons_self different,
        ih tailUnique member (fun y hy => zero y (List.mem_cons_of_mem x hy)), Rat.zero_add]

variable {V : Type} [DecidableEq V] {vertices : List V}

theorem endpointProbabilityAsSum (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start finish : Subset V) (n : Nat) :
    endpointProbability edges bs start finish n = sumWith
      (fun h => if h n = finish then pathProbability edges bs h n else 0) (histories edges start n) := by
  unfold endpointProbability endpointHistories
  rw [sumWithFilter]
  apply sumWithCongrOn
  intro h _
  have iff : vertices.all (fun x => h n x == finish x) = true ↔ h n = finish := by
    constructor
    · intro same
      funext x
      simpa only [beq_iff_eq] using List.all_eq_true.mp same x ((edges 0).complete x)
    · intro same
      simp [same]
  simp only [iff]

/-- Every canonical legal history ends at exactly one canonical ordered pair. -/
theorem endpointPairMassOne (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) (size : (vertices.filter start).length = 2) :
    sumWith (fun p => endpointProbability edges bs start (pairSubset p) n) (orderedPairs vertices) = 1 := by
  have expand := sumWithCongr _ _ (orderedPairs vertices)
    (fun p => endpointProbabilityAsSum edges bs start (pairSubset p) n)
  rw [expand, sumWithSwap]
  have reduce : sumWith (fun h => sumWith
      (fun p => if h n = pairSubset p then pathProbability edges bs h n else 0)
      (orderedPairs vertices)) (histories edges start n) =
      sumWith (fun h => pathProbability edges bs h n) (histories edges start n) := by
    apply sumWithCongrOn
    intro h member
    obtain ⟨p, pm, _, same⟩ := twoSubsetRepresentative vertices (edges 0).complete (h n)
      (historySizeTwo edges start size n h member)
    have one := sumWithUniqueKey (orderedPairs vertices) pairSubset
      (pairSubsetsNodup vertices (edges 0).nodup) p pm
      (fun q => if h n = pairSubset q then pathProbability edges bs h n else 0)
      (by
        intro q _ different
        have ne : h n ≠ pairSubset q := by intro eq; exact different (eq.symm.trans same)
        simp [ne])
    simpa [same] using one
  rw [reduce]
  exact pathMassOne edges bs start n

/-- The squared actual order-two product minors form a normalized finite law
across all canonical output pairs. -/
theorem productPairMassOne (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start : Subset V) (n : Nat)
    (size : (vertices.filter start).length = 2)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    sumWith (fun p => normSq (minorDet gaussianArithmetic
      (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      [p.1,p.2] (vertices.filter start))) (orderedPairs vertices) = 1 := by
  have replace : sumWith (fun p => normSq (minorDet gaussianArithmetic
        (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
        [p.1,p.2] (vertices.filter start))) (orderedPairs vertices) =
      sumWith (fun p => endpointProbability edges bs start (pairSubset p) n) (orderedPairs vertices) := by
    apply sumWithCongrOn
    intro p member
    have filtered := pairFilter vertices (edges 0).nodup p member
    rw [← filtered]
    exact productPairMagnitude edges as bs cs ds start (pairSubset p) n size
      (by rw [filtered]; rfl) bridges unitary
  rw [replace]
  exact endpointPairMassOne edges bs start n size

end Mass
#print axioms pairMemberSublist
#print axioms filterSingle
#print axioms pairFilter
#print axioms subsetOfFilter
#print axioms pairSubsetInjective
#print axioms orderedPairsNodup
#print axioms nodupMapOn
#print axioms pairSubsetsNodup
#print axioms twoSubsetRepresentative
#print axioms identitySubsetMinor
#print axioms gaussianSumFilter
#print axioms endpointAsSum
#print axioms pathAmplitudePrefixCongr
#print axioms extendedPathAmplitude
#print axioms successorMinorSum
#print axioms endpointAmplitudeRecurrence
#print axioms historySizeTwo
#print axioms productPairEndpointAmplitude
#print axioms productPairMagnitude
#print axioms productPairEmptyFiber
#print axioms sumWithCongrOn
#print axioms sumWithFilter
#print axioms sumWithAdd
#print axioms sumWithSwap
#print axioms sumWithUniqueKey
#print axioms endpointProbabilityAsSum
#print axioms endpointPairMassOne
#print axioms productPairMassOne

end MathCommons.ForestPairHistory
