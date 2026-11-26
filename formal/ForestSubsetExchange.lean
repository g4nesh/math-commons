/- Canonical single-exchange signs, sparse actual gate-minor contraction,
and the arbitrary-order squared forest product-minor probability law.
Scalar domain: Gaussian rationals, using Std Lean 4.24.0 only. -/
import ForestMinorAlternation
import ForestPairHistory

namespace MathCommons.ForestSubsetExchange
open Classical
open ForestTrajectory ForestGateSupport ForestGateMagnitude ForestPathProbability
open ForestPathNormalization ForestEndpointMass ForestProductMinor ForestPairHistory ForestMinorAlternation

section Exchange
variable {V : Type} [DecidableEq V] {vertices : List V}

theorem signedTwice (odd : Bool) (z : Gaussian) :
    signed gaussianArithmetic odd (signed gaussianArithmetic odd z) = z := by
  cases odd
  · rfl
  · exact gaussianNegNeg z

theorem swapInvolutive (e : Edge vertices) (s : Subset V) : swapSubset e (swapSubset e s) = s := by
  funext x
  by_cases hl : x = e.left
  · subst x; simp only [swapAtLeft, swapAtRight]
  · by_cases hr : x = e.right
    · subst x; simp only [swapAtRight, swapAtLeft]
    · simp [swapSubset, hl, hr]

theorem stepSymmetric (e : Edge vertices) (s t : Subset V) (step : Step e s t) : Step e t s := by
  rcases step with rfl | rfl
  · exact Or.inl rfl
  · exact Or.inr (swapInvolutive e s).symm

theorem reverseSwap (e : Edge vertices) (s : Subset V) : swapSubset (reverseEdge e) s = swapSubset e s := by
  funext x
  have distinct := e.distinct
  simp only [swapSubset, reverseEdge]
  grind

theorem filterSwapAvoids (e : Edge vertices) (s : Subset V) (xs : List V)
    (avoids : ∀ x ∈ xs, x ≠ e.left ∧ x ≠ e.right) :
    xs.filter (swapSubset e s) = xs.filter s := by
  apply List.filter_congr
  intro x hx
  simp [swapSubset, (avoids x hx).1, (avoids x hx).2]

theorem replaceRowList (M : V → V → Gaussian) (before after cols : List V) (r s : V)
    (absent : r ∉ before ++ after) :
    minorDet gaussianArithmetic (replaceRow M r (M s)) (before ++ r :: after) cols =
      minorDet gaussianArithmetic M (before ++ s :: after) cols := by
  rw [minorAsRows, minorAsRows]
  have mb : before.map (replaceRow M r (M s)) = before.map M := by
    apply List.map_congr_left
    intro x hx
    funext j
    have ne : x ≠ r := by intro eq; subst x; exact absent (by simp [hx])
    simp [replaceRow, ne]
  have ma : after.map (replaceRow M r (M s)) = after.map M := by
    apply List.map_congr_left
    intro x hx
    funext j
    have ne : x ≠ r := by intro eq; subst x; exact absent (by simp [hx])
    simp [replaceRow, ne]
  have atR : replaceRow M r (M s) r = M s := by funext j; simp [replaceRow]
  simp only [List.map_append, List.map_cons, mb, ma, atR]

/-- The canonical filter sign when the left endpoint occurs first in the vertex enumeration. -/
theorem leftReplacementForward (e : Edge vertices) (s : Subset V)
    (before middle after : List V)
    (order : vertices = before ++ e.left :: (middle ++ e.right :: after))
    (left : s e.left = true) (right : s e.right = false)
    (M : V → V → Gaussian) (cols : List V) :
    minorDet gaussianArithmetic (replaceRow M e.left (M e.right)) (vertices.filter s) cols =
      signed gaussianArithmetic (parity (middle.filter s))
        (minorDet gaussianArithmetic M (vertices.filter (swapSubset e s)) cols) := by
  have avoids := twoEndpointAvoids before middle after e.left e.right (order ▸ e.nodup)
  have fb := filterSwapAvoids e s before (fun x hx => avoids x (by simp [hx]))
  have fm := filterSwapAvoids e s middle (fun x hx => avoids x (by simp [hx]))
  have fa := filterSwapAvoids e s after (fun x hx => avoids x (by simp [hx]))
  have original : vertices.filter s = before.filter s ++ e.left :: (middle.filter s ++ after.filter s) := by
    simpa [List.filter_append, left, right] using congrArg (List.filter s) order
  have swapped : vertices.filter (swapSubset e s) =
      before.filter s ++ middle.filter s ++ e.right :: after.filter s := by
    simpa [List.filter_append, swapAtLeft, swapAtRight, left, right, fb, fm, fa, List.append_assoc]
      using congrArg (List.filter (swapSubset e s)) order
  rw [original, swapped]
  have absent : e.left ∉ before.filter s ++ (middle.filter s ++ after.filter s) := by
    intro h
    have member : e.left ∈ before ++ middle ++ after := by
      simp only [List.mem_append, List.mem_filter] at h ⊢
      grind
    exact (avoids e.left member).1 rfl
  rw [replaceRowList M (before.filter s) (middle.filter s ++ after.filter s) cols e.left e.right absent]
  rw [minorMoveRow, signedTwice]

/-- The same selected-spectator parity when the left endpoint occurs second. -/
theorem leftReplacementBackward (e : Edge vertices) (s : Subset V)
    (before middle after : List V)
    (order : vertices = before ++ e.right :: (middle ++ e.left :: after))
    (left : s e.left = true) (right : s e.right = false)
    (M : V → V → Gaussian) (cols : List V) :
    minorDet gaussianArithmetic (replaceRow M e.left (M e.right)) (vertices.filter s) cols =
      signed gaussianArithmetic (parity (middle.filter s))
        (minorDet gaussianArithmetic M (vertices.filter (swapSubset e s)) cols) := by
  have avoids := twoEndpointAvoids before middle after e.right e.left (order ▸ e.nodup)
  have avoid (x : V) (hx : x ∈ before ++ middle ++ after) : x ≠ e.left ∧ x ≠ e.right :=
    ⟨(avoids x hx).2, (avoids x hx).1⟩
  have fb := filterSwapAvoids e s before (fun x hx => avoid x (by simp [hx]))
  have fm := filterSwapAvoids e s middle (fun x hx => avoid x (by simp [hx]))
  have fa := filterSwapAvoids e s after (fun x hx => avoid x (by simp [hx]))
  have original : vertices.filter s = before.filter s ++ middle.filter s ++ e.left :: after.filter s := by
    simpa [List.filter_append, left, right, List.append_assoc] using congrArg (List.filter s) order
  have swapped : vertices.filter (swapSubset e s) =
      before.filter s ++ e.right :: (middle.filter s ++ after.filter s) := by
    simpa [List.filter_append, swapAtLeft, swapAtRight, left, right, fb, fm, fa]
      using congrArg (List.filter (swapSubset e s)) order
  rw [original, swapped]
  have absent : e.left ∉ (before.filter s ++ middle.filter s) ++ after.filter s := by
    intro h
    have member : e.left ∈ before ++ middle ++ after := by
      simp only [List.mem_append, List.mem_filter] at h ⊢
      grind
    exact (avoid e.left member).1 rfl
  rw [replaceRowList M (before.filter s ++ middle.filter s) (after.filter s) cols e.left e.right absent]
  exact minorMoveRow M _ _ _ _ _


/-- Right-occupied exchange when the original left endpoint appears first. -/
theorem rightReplacementForward (e : Edge vertices) (s : Subset V)
    (before middle after : List V)
    (order : vertices = before ++ e.left :: (middle ++ e.right :: after))
    (left : s e.left = false) (right : s e.right = true)
    (M : V → V → Gaussian) (cols : List V) :
    minorDet gaussianArithmetic (replaceRow M e.right (M e.left)) (vertices.filter s) cols =
      signed gaussianArithmetic (parity (middle.filter s))
        (minorDet gaussianArithmetic M (vertices.filter (swapSubset e s)) cols) := by
  simpa only [reverseSwap] using
    leftReplacementBackward (reverseEdge e) s before middle after order right left M cols

theorem rightReplacementBackward (e : Edge vertices) (s : Subset V)
    (before middle after : List V)
    (order : vertices = before ++ e.right :: (middle ++ e.left :: after))
    (left : s e.left = false) (right : s e.right = true)
    (M : V → V → Gaussian) (cols : List V) :
    minorDet gaussianArithmetic (replaceRow M e.right (M e.left)) (vertices.filter s) cols =
      signed gaussianArithmetic (parity (middle.filter s))
        (minorDet gaussianArithmetic M (vertices.filter (swapSubset e s)) cols) := by
  simpa only [reverseSwap] using
    leftReplacementForward (reverseEdge e) s before middle after order right left M cols

theorem leftReplacementSign (e : Edge vertices) (s : Subset V)
    (left : s e.left = true) (right : s e.right = false) :
    ∃ odd : Bool, ∀ (M : V → V → Gaussian) (cols : List V),
      minorDet gaussianArithmetic (replaceRow M e.left (M e.right)) (vertices.filter s) cols =
        signed gaussianArithmetic odd (minorDet gaussianArithmetic M (vertices.filter (swapSubset e s)) cols) := by
  rcases decomposeTwoEndpoints vertices e.left e.right e.distinct (e.complete _) (e.complete _) with
    ⟨before, middle, after, order⟩ | ⟨before, middle, after, order⟩
  · exact ⟨parity (middle.filter s), leftReplacementForward e s before middle after order left right⟩
  · exact ⟨parity (middle.filter s), leftReplacementBackward e s before middle after order left right⟩

theorem rightReplacementSign (e : Edge vertices) (s : Subset V)
    (left : s e.left = false) (right : s e.right = true) :
    ∃ odd : Bool, ∀ (M : V → V → Gaussian) (cols : List V),
      minorDet gaussianArithmetic (replaceRow M e.right (M e.left)) (vertices.filter s) cols =
        signed gaussianArithmetic odd (minorDet gaussianArithmetic M (vertices.filter (swapSubset e s)) cols) := by
  obtain ⟨odd, sign⟩ := leftReplacementSign (reverseEdge e) s right left
  refine ⟨odd, ?_⟩
  intro M cols
  simpa only [reverseSwap] using sign M cols

theorem signedMulLeft (odd : Bool) (a b : Gaussian) :
    gaussianArithmetic.mul a (signed gaussianArithmetic odd b) =
      gaussianArithmetic.mul (signed gaussianArithmetic odd a) b := by
  cases odd
  · rfl
  · apply Prod.ext <;> simp only [signed, ↓reduceIte, gaussianArithmetic] <;> grind

/-- Actual sparse row multiplication has two coefficient functions, independent
of the remaining matrix and its chosen column list. Collision sectors have a
zero exchange coefficient. The coefficients are derived, not assumed. -/
theorem gateTwoCoefficients (e : Edge vertices) (a b c d : Gaussian) (t : Subset V) :
    ∃ alpha beta : Gaussian, (t e.left = t e.right → beta = gaussianArithmetic.zero) ∧
      ∀ (M : V → V → Gaussian) (cols : List V),
        minorDet gaussianArithmetic
          (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M) (vertices.filter t) cols =
        gaussianArithmetic.add
          (gaussianArithmetic.mul alpha (minorDet gaussianArithmetic M (vertices.filter t) cols))
          (gaussianArithmetic.mul beta (minorDet gaussianArithmetic M (vertices.filter (swapSubset e t)) cols)) := by
  have unique := e.nodup.filter t
  cases hl : t e.left <;> cases hr : t e.right
  · refine ⟨gaussianArithmetic.one, gaussianArithmetic.zero, by intros; rfl, ?_⟩
    intro M cols
    rw [gateMinorNeither e a b c d M _ cols (by simp [hl]) (by simp [hr])]
    simp only [gaussianOneMul, gaussianArithmetic.zeroMul, gaussianAddZero]
  · obtain ⟨odd, sign⟩ := rightReplacementSign e t hl hr
    refine ⟨d, signed gaussianArithmetic odd c, ?_, ?_⟩
    · intro collision; contradiction
    · intro M cols
      rw [gateMinorRightOnly e a b c d M _ cols unique (by simp [hl])
        (List.mem_filter.mpr ⟨e.complete _, hr⟩), sign M cols, signedMulLeft]
      exact gaussianAddComm _ _
  · obtain ⟨odd, sign⟩ := leftReplacementSign e t hl hr
    refine ⟨a, signed gaussianArithmetic odd b, ?_, ?_⟩
    · intro collision; contradiction
    · intro M cols
      rw [gateMinorLeftOnly e a b c d M _ cols unique
        (List.mem_filter.mpr ⟨e.complete _, hl⟩) (by simp [hr]), sign M cols, signedMulLeft]
  · refine ⟨blockDet a b c d, gaussianArithmetic.zero, by intros; rfl, ?_⟩
    intro M cols
    rw [gateMinorBoth e a b c d M _ cols unique
      (List.mem_filter.mpr ⟨e.complete _, hl⟩) (List.mem_filter.mpr ⟨e.complete _, hr⟩)]
    simp only [gaussianArithmetic.zeroMul, gaussianAddZero]


/-- A sparse Cauchy--Binet identity for an actual two-coordinate left factor,
with the canonical Boolean-filter ordering and no restriction on minor order. -/
theorem gateMinorPredecessors (e : Edge vertices) (a b c d : Gaussian)
    (M : V → V → Gaussian) (t : Subset V) (cols : List V) :
    minorDet gaussianArithmetic
      (matrixMul vertices (twoCoordinateGate gaussianArithmetic e a b c d) M)
      (vertices.filter t) cols =
    scalarSum gaussianArithmetic ((successors e t).map fun s =>
      gaussianArithmetic.mul
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter s))
        (minorDet gaussianArithmetic M (vertices.filter s) cols)) := by
  classical
  obtain ⟨alpha, beta, collisionZero, law⟩ := gateTwoCoefficients e a b c d t
  have self := law identityMatrix (vertices.filter t)
  have exchange := law identityMatrix (vertices.filter (swapSubset e t))
  rw [matrixMulIdentity vertices e.nodup e.complete] at self exchange
  by_cases collision : t e.left = t e.right
  · have same := swapCollision e t collision
    have zero := collisionZero collision
    simp only [identitySubsetMinor vertices e.nodup e.complete, ↓reduceIte,
      zero, gaussianArithmetic.zeroMul, gaussianMulOne, gaussianAddZero] at self
    rw [law M cols, zero, gaussianArithmetic.zeroMul, gaussianAddZero]
    simp only [successors, collision, ↓reduceIte, List.map_cons, List.map_nil, scalarSum,
      self, gaussianAddZero]
  · have different : swapSubset e t ≠ t := by
      intro same
      have atLeft := congrFun same e.left
      rw [swapAtLeft] at atLeft
      exact collision atLeft.symm
    have reverseDifferent : t ≠ swapSubset e t := Ne.symm different
    simp only [identitySubsetMinor vertices e.nodup e.complete, different, reverseDifferent,
      ↓reduceIte, gaussianMulOne, gaussianArithmetic.mulZero, gaussianAddZero,
      gaussianZeroAdd] at self exchange
    rw [law M cols]
    simp only [successors, collision, ↓reduceIte, List.map_cons, List.map_nil, scalarSum,
      self, exchange, gaussianAddZero]


/-- The reverse-neighbor list also contracts an endpoint indicator against an
actual gate minor; non-neighbors have zero actual minor. -/
theorem predecessorMinorSum (e : Edge vertices) (a b c d : Gaussian)
    (s t : Subset V) (weight : Gaussian) :
    scalarSum gaussianArithmetic ((successors e t).map fun pred =>
      gaussianArithmetic.mul
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter pred))
        (if s = pred then weight else gaussianArithmetic.zero)) =
      gaussianArithmetic.mul weight
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter s)) := by
  classical
  by_cases member : s ∈ successors e t
  · have one := scalarSumUniqueKey (successors e t) (fun x => x)
      (by simpa using successorsNodup e t) s member
      (fun pred => gaussianArithmetic.mul
        (minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter pred))
        (if s = pred then weight else gaussianArithmetic.zero))
      (by intro pred _ ne; simp only [Ne.symm ne, ↓reduceIte, gaussianArithmetic.mulZero])
    rw [one]
    simp only [↓reduceIte]
    exact gaussianMulComm _ _
  · have zero : minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
        (vertices.filter t) (vertices.filter s) = gaussianArithmetic.zero := by
      by_cases h : minorDet gaussianArithmetic (twoCoordinateGate gaussianArithmetic e a b c d)
          (vertices.filter t) (vertices.filter s) = gaussianArithmetic.zero
      · exact h
      · exact False.elim (member (everyLegalSuccessor e t s
          (stepSymmetric e s t (gaussianGateMinorStep e a b c d s t h))))
    rw [zero, gaussianArithmetic.mulZero]
    apply scalarSumMapZero
    intro pred hp
    have ne : s ≠ pred := by intro equal; exact member (equal ▸ hp)
    simp only [ne, ↓reduceIte, gaussianArithmetic.mulZero]

/-- Actual product minors, of every order, equal the canonical signed endpoint
amplitude. No unitarity or forest assumption enters this algebraic bridge. -/
theorem productEndpointAmplitude (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat) :
    minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start) =
        endpointAmplitude edges as bs cs ds start finish n := by
  classical
  induction n generalizing finish with
  | zero =>
    rw [matrixProduct, identitySubsetMinor vertices (edges 0).nodup (edges 0).complete]
    rw [endpointAsSum]
    simp only [histories, List.map_cons, List.map_nil, scalarSum, pathAmplitude, gaussianAddZero]
    by_cases eq : finish = start
    · simp [eq]
    · simp [eq, Ne.symm eq]
  | succ n ih =>
    have cb := gateMinorPredecessors (edges n) (as n) (bs n) (cs n) (ds n)
      (matrixProduct vertices (gateSchedule edges as bs cs ds) n) finish (vertices.filter start)
    change minorDet gaussianArithmetic
      (matrixMul vertices (gateSchedule edges as bs cs ds n)
        (matrixProduct vertices (gateSchedule edges as bs cs ds) n))
      (vertices.filter finish) (vertices.filter start) = _ at cb
    rw [matrixProduct, cb]
    have replace : scalarSum gaussianArithmetic ((successors (edges n) finish).map fun pred =>
        gaussianArithmetic.mul
          (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
            (vertices.filter finish) (vertices.filter pred))
          (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
            (vertices.filter pred) (vertices.filter start))) =
        scalarSum gaussianArithmetic ((successors (edges n) finish).map fun pred =>
          scalarSum gaussianArithmetic ((histories edges start n).map fun h =>
            gaussianArithmetic.mul
              (minorDet gaussianArithmetic (gateSchedule edges as bs cs ds n)
                (vertices.filter finish) (vertices.filter pred))
              (if h n = pred then
                pathAmplitude (vertices := vertices) gaussianArithmetic (gateSchedule edges as bs cs ds) h n
               else gaussianArithmetic.zero))) := by
      apply gaussianSumCongr
      intro pred _
      rw [ih pred, endpointAsSum]
      exact (gaussianSumScale _ _ _).symm
    calc
      _ = _ := replace
      _ = _ := by
        rw [gaussianSumSwap, endpointAmplitudeRecurrence]
        apply gaussianSumCongr
        intro h _
        exact predecessorMinorSum (edges n) (as n) (bs n) (cs n) (ds n) (h n) finish _

/-- The squared actual minor of a unitary forest gate product is the canonical
endpoint probability, without a restriction on the chosen subset size. -/
theorem productMagnitude (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    normSq (minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start)) =
        endpointProbability edges bs start finish n := by
  rw [productEndpointAmplitude]
  exact endpointMagnitude edges as bs cs ds start finish n bridges unitary

/-- An empty canonical endpoint fiber forces an actual zero product minor,
including without unitarity or a forest hypothesis. -/
theorem productEmptyFiber (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start finish : Subset V) (n : Nat)
    (empty : endpointHistories edges start finish n = []) :
    minorDet gaussianArithmetic (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start) = gaussianArithmetic.zero := by
  rw [productEndpointAmplitude]
  exact (emptyEndpointFiber edges as bs cs ds start finish n empty).1


/-- Finite deduplication with explicit classical equality. -/
noncomputable def distinctList {α : Type} : List α → List α
  | [] => []
  | x :: xs => if x ∈ distinctList xs then distinctList xs else x :: distinctList xs

omit [DecidableEq V] in
theorem distinctListMember {α : Type} (xs : List α) (x : α) :
    x ∈ distinctList xs ↔ x ∈ xs := by
  induction xs generalizing x with
  | nil => simp [distinctList]
  | cons y ys ih =>
    simp only [distinctList]
    split
    · rename_i member
      have original := (ih y).mp member
      simp only [ih, List.mem_cons]
      constructor
      · exact Or.inr
      · intro h; rcases h with rfl | h
        · exact original
        · exact h
    · simp [List.mem_cons, ih]

omit [DecidableEq V] in
theorem distinctListNodup {α : Type} (xs : List α) : (distinctList xs).Nodup := by
  induction xs with
  | nil => simp [distinctList]
  | cons x xs ih =>
    simp only [distinctList]
    split
    · exact ih
    · rename_i absent; exact List.nodup_cons.mpr ⟨absent, ih⟩

/-- Each reachable endpoint appears once, regardless of the number of legal
histories reaching it. This finite list uses classical function equality. -/
noncomputable def reachableEndpoints (edges : Nat → Edge vertices)
    (start : Subset V) (n : Nat) : List (Subset V) :=
  distinctList ((histories edges start n).map fun h => h n)

theorem reachableEndpointsNodup (edges : Nat → Edge vertices)
    (start : Subset V) (n : Nat) : (reachableEndpoints edges start n).Nodup := by
  exact distinctListNodup _

theorem historyEndpointMember (edges : Nat → Edge vertices)
    (start : Subset V) (n : Nat) (h : Nat → Subset V)
    (member : h ∈ histories edges start n) : h n ∈ reachableEndpoints edges start n := by
  exact (distinctListMember _ _).mpr (List.mem_map.mpr ⟨h, member, rfl⟩)

/-- The endpoint masses sum to one over distinct reachable subsets, in all
orders. As with pathMassOne, algebraic normalization alone needs no bounds. -/
theorem endpointMassOne (edges : Nat → Edge vertices) (bs : Nat → Gaussian)
    (start : Subset V) (n : Nat) :
    sumWith (fun finish => endpointProbability edges bs start finish n)
      (reachableEndpoints edges start n) = 1 := by
  have expand := sumWithCongr _ _ (reachableEndpoints edges start n)
    (fun finish => endpointProbabilityAsSum edges bs start finish n)
  rw [expand, sumWithSwap]
  have reduce : sumWith (fun h => sumWith
      (fun finish => if h n = finish then pathProbability edges bs h n else 0)
      (reachableEndpoints edges start n)) (histories edges start n) =
      sumWith (fun h => pathProbability edges bs h n) (histories edges start n) := by
    apply sumWithCongrOn
    intro h member
    have one := sumWithUniqueKey (reachableEndpoints edges start n) (fun x => x)
      (by simpa using reachableEndpointsNodup edges start n)
      (h n) (historyEndpointMember edges start n h member)
      (fun finish => if h n = finish then pathProbability edges bs h n else 0)
      (by intro finish _ different; simp [Ne.symm different])
    simpa using one
  rw [reduce]
  exact pathMassOne edges bs start n

/-- Squared actual forest product minors form a normalized law on the distinct
reachable subsets; unreachable subsets have zero minor by productEmptyFiber. -/
theorem productMassOne (edges : Nat → Edge vertices)
    (as bs cs ds : Nat → Gaussian) (start : Subset V) (n : Nat)
    (bridges : EveryEdgeIsBridge edges n)
    (unitary : ∀ k, k < n → UnitaryBlock (as k) (bs k) (cs k) (ds k)) :
    sumWith (fun finish => normSq (minorDet gaussianArithmetic
      (matrixProduct vertices (gateSchedule edges as bs cs ds) n)
      (vertices.filter finish) (vertices.filter start))) (reachableEndpoints edges start n) = 1 := by
  have replace := sumWithCongr _ _ (reachableEndpoints edges start n)
    (fun finish => productMagnitude edges as bs cs ds start finish n bridges unitary)
  rw [replace]
  exact endpointMassOne edges bs start n

end Exchange
#print axioms signedTwice
#print axioms swapInvolutive
#print axioms stepSymmetric
#print axioms reverseSwap
#print axioms filterSwapAvoids
#print axioms replaceRowList
#print axioms leftReplacementForward
#print axioms leftReplacementBackward
#print axioms rightReplacementForward
#print axioms rightReplacementBackward
#print axioms leftReplacementSign
#print axioms rightReplacementSign
#print axioms signedMulLeft
#print axioms gateTwoCoefficients
#print axioms gateMinorPredecessors
#print axioms predecessorMinorSum
#print axioms productEndpointAmplitude
#print axioms productMagnitude
#print axioms productEmptyFiber
#print axioms distinctListMember
#print axioms distinctListNodup
#print axioms reachableEndpointsNodup
#print axioms historyEndpointMember
#print axioms endpointMassOne
#print axioms productMassOne

end MathCommons.ForestSubsetExchange
