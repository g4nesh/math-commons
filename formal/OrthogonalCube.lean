import Std

set_option maxRecDepth 4096

/-!
# Claude's unit-cube quadratic lemma over an arbitrary ordered field

The domain is a `Lean.Grind.Field` with a lawful linear order and compatible
strict ordered-ring operations. These are Lean's standard algebra/order laws,
not extra sign assumptions for the lemma. No completeness or square roots are
used. The matrix, SO(3), SO(4), Hodge and PMSH reductions are separate claims.

Source: Claude's Lemma L in PR48, commit
41ca275e7354043979278d72f8db897b9e2905a0. `lemmaL` and `lemmaLMax`
derive all four facts F1–F4 from the original hypotheses, including the
boundary of the cube and zero parameters.
-/

namespace MathCommons.OrthogonalCube
open Std Lean.Grind

variable {R : Type} [Field R] [LE R] [LT R]
  [LawfulOrderLT R] [IsLinearOrder R] [OrderedRing R]

/-- The two-dimensional symmetric quadratic form with off-diagonal entry `-B`. -/
def quadratic (A B C x y : R) : R := A*x*x - 2*B*x*y + C*y*y

/-- Its associated symmetric bilinear form. -/
def bilinear (A B C x y x' y' : R) : R := A*x*x' - B*(x*y'+y*x') + C*y*y'

/-- The exact quadratic form in Lemma L. -/
def cubeForm (p q r u a b c d : R) : R :=
  p*r*(a*a+c*c) + q*u*(b*b+d*d) - 2*p*q*a*b - 2*r*u*c*d

/-- A square is nonnegative, including zero. -/
theorem squareNonnegative (x : R) : 0 ≤ x*x := by
  have h := OrderedRing.sq_nonneg (a := x)
  grind

/-- The diagonal/determinant criterion proves positivity; it is not assumed. -/
theorem quadraticNonnegative (A B C x y : R) (hA : 0 ≤ A) (hC : 0 ≤ C)
    (hdet : 0 ≤ A*C-B*B) : 0 ≤ quadratic A B C x y := by
  have hs := squareNonnegative (A*x-B*y)
  have ht := OrderedRing.mul_nonneg hdet (squareNonnegative y)
  by_cases ha : A = 0
  · subst A
    have hb : B = 0 := by
      have hbb := squareNonnegative B
      have hz : B*B = 0 := by grind
      grind
    subst B
    have hcy := OrderedRing.mul_nonneg hC (squareNonnegative y)
    grind [quadratic]
  · have hap : 0 < A := by grind
    by_cases hn : 0 ≤ quadratic A B C x y
    · exact hn
    have hneg : A*x*x - 2*B*x*y+C*y*y < 0 := by grind [quadratic]
    have hm := OrderedRing.mul_lt_mul_of_pos_left hneg hap
    grind [quadratic]

/-- The two-dimensional Lagrange identity, with no sign conditions. -/
theorem lagrangeIdentity (A B C x y x' y' : R) :
    quadratic A B C x y * quadratic A B C x' y'
      - bilinear A B C x y x' y' * bilinear A B C x y x' y'
      = (A*C-B*B)*(x*y'-y*x')*(x*y'-y*x') := by
  grind [quadratic, bilinear]

/-- The Cauchy inequality follows from the determinant criterion. -/
theorem quadraticCauchy (A B C x y x' y' : R) (hdet : 0 ≤ A*C-B*B) :
    bilinear A B C x y x' y' * bilinear A B C x y x' y'
      ≤ quadratic A B C x y * quadratic A B C x' y' := by
  have identity := lagrangeIdentity A B C x y x' y'
  have positive := OrderedRing.mul_nonneg hdet (squareNonnegative (x*y'-y*x'))
  grind

/-- The complementary determinant sign gives the reverse squared inequality. -/
theorem quadraticReverseCauchy (A B C x y x' y' : R) (hdet : 0 ≤ B*B-A*C) :
    quadratic A B C x y * quadratic A B C x' y'
      ≤ bilinear A B C x y x' y' * bilinear A B C x y x' y' := by
  have identity := lagrangeIdentity A B C x y x' y'
  have positive := OrderedRing.mul_nonneg hdet (squareNonnegative (x*y'-y*x'))
  grind

/-- F1's positive-semidefinite assertion follows from the original parameters. -/
theorem firstPartNonnegative (p q r u a b : R)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hr : 0 ≤ r) (hu : 0 ≤ u)
    (ordered : p*q ≤ r*u) : 0 ≤ quadratic (p*r) (p*q) (q*u) a b := by
  have hpr := OrderedRing.mul_nonneg hp hr
  have hqu := OrderedRing.mul_nonneg hq hu
  have hpq := OrderedRing.mul_nonneg hp hq
  have gap : 0 ≤ r*u-p*q := by grind
  have product := OrderedRing.mul_nonneg hpq gap
  have det : 0 ≤ (p*r)*(q*u)-(p*q)*(p*q) := by grind
  exact quadraticNonnegative (p*r) (p*q) (q*u) a b hpr hqu det

/-- F1's squared bilinear inequality is derived for arbitrary coordinates. -/
theorem firstPartCauchy (p q r u a b a' b' : R)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (ordered : p*q ≤ r*u) :
    bilinear (p*r) (p*q) (q*u) a b a' b' * bilinear (p*r) (p*q) (q*u) a b a' b'
      ≤ quadratic (p*r) (p*q) (q*u) a b * quadratic (p*r) (p*q) (q*u) a' b' := by
  have product := OrderedRing.mul_nonneg (OrderedRing.mul_nonneg hp hq)
    (show 0 ≤ r*u-p*q by grind)
  apply quadraticCauchy
  grind

/-- F2's reverse squared inequality is derived, including zero parameters. -/
theorem secondPartReverseCauchy (p q r u a b a' b' : R)
    (hr : 0 ≤ r) (hu : 0 ≤ u) (ordered : p*q ≤ r*u) :
    quadratic (p*r) (r*u) (q*u) a b * quadratic (p*r) (r*u) (q*u) a' b'
      ≤ bilinear (p*r) (r*u) (q*u) a b a' b' * bilinear (p*r) (r*u) (q*u) a b a' b' := by
  have product := OrderedRing.mul_nonneg (OrderedRing.mul_nonneg hr hu)
    (show 0 ≤ r*u-p*q by grind)
  apply quadraticReverseCauchy
  grind

/-- F3's polarization identity for complementary cube points. -/
theorem complementPolarization (A B C a b : R) :
    quadratic A B C a b + quadratic A B C (1-a) (1-b)
      + 2*bilinear A B C a b (1-a) (1-b) = A-2*B+C := by
  grind [quadratic, bilinear]

/-- The balance hypothesis makes the total value at one a nonnegative square. -/
theorem totalAtOne (p q r u : R) (balance : p+q=r+u) :
    (p*r-2*(p*q)+q*u)+(p*r-2*(r*u)+q*u) = 2*(p-u)*(p-u) := by
  grind

/-- Comparing nonnegative squares does not require square roots. -/
theorem strictOrderOfSquares (T S : R) (hS : 0 ≤ S) (hsq : T*T < S*S) : T < S := by
  by_cases h : T < S
  · exact h
  have difference : 0 ≤ T-S := by grind
  have sum : 0 ≤ T+S := by grind
  have product := OrderedRing.mul_nonneg difference sum
  grind

/-- The weighted Cauchy inequality needed in F4, proved by its polynomial identity. -/
theorem weightedCauchy (A C x y x' y' : R) (hA : 0 ≤ A) (hC : 0 ≤ C) :
    (A*x*x'+C*y*y')*(A*x*x'+C*y*y')
      ≤ (A*x*x+C*y*y)*(A*x'*x'+C*y'*y') := by
  have positive := OrderedRing.mul_nonneg (OrderedRing.mul_nonneg hA hC)
    (squareNonnegative (x*y'-y*x'))
  grind

/-- F4: two negative vectors in the nonnegative quadrant have negative pairing. -/
theorem strictWeightedCross (A C H x y x' y' : R)
    (hA : 0 ≤ A) (hC : 0 ≤ C) (hH : 0 ≤ H)
    (hx : 0 ≤ x) (hy : 0 ≤ y) (hx' : 0 ≤ x') (hy' : 0 ≤ y')
    (first : A*x*x+C*y*y < 2*H*x*y)
    (second : A*x'*x'+C*y'*y' < 2*H*x'*y') :
    A*x*x'+C*y*y' < H*(x*y'+y*x') := by
  have hxx := OrderedRing.mul_nonneg hA (squareNonnegative x)
  have hyy := OrderedRing.mul_nonneg hC (squareNonnegative y)
  have hxx' := OrderedRing.mul_nonneg hA (squareNonnegative x')
  have hyy' := OrderedRing.mul_nonneg hC (squareNonnegative y')
  have diagonal : 0 ≤ A*x*x+C*y*y := by grind
  have diagonal' : 0 ≤ A*x'*x'+C*y'*y' := by grind
  have rightPositive : 0 < 2*H*x*y := by grind
  have firstWeak : A*x*x+C*y*y ≤ 2*H*x*y := by grind
  have step1 := OrderedRing.mul_le_mul_of_nonneg_right firstWeak diagonal'
  have step2 := OrderedRing.mul_lt_mul_of_pos_left second rightPositive
  have cauchy := weightedCauchy A C x y x' y' hA hC
  have crossSquare := OrderedRing.mul_nonneg (squareNonnegative H)
    (squareNonnegative (x*y'-y*x'))
  have strictSquares : (A*x*x'+C*y*y')*(A*x*x'+C*y*y')
      < (H*(x*y'+y*x'))*(H*(x*y'+y*x')) := by
    grind
  have cross1 := OrderedRing.mul_nonneg hx hy'
  have cross2 := OrderedRing.mul_nonneg hy hx'
  have cross : 0 ≤ x*y'+y*x' := by grind
  exact strictOrderOfSquares _ _ (OrderedRing.mul_nonneg hH cross) strictSquares

/-- The last strict chain in the argument, with every comparison explicit. -/
theorem contradictionFromParts (a a' b b' m n : R)
    (ha : 0 ≤ a) (ha' : 0 ≤ a') (hab : a < b) (hab' : a' < b')
    (first : m*m ≤ a*a') (second : b*b' ≤ n*n)
    (balance : 0 ≤ a+a'+2*m-(b+b'+2*n)) (hn : 0 < n) : False := by
  have hm : 0 < m-n := by grind
  have hsum : 0 < m+n := by grind
  have squares := OrderedRing.mul_pos hm hsum
  have hdiff : 0 ≤ b-a := by grind
  have hdiff' : 0 ≤ b'-a' := by grind
  have hb : 0 ≤ b := by grind
  have product1 := OrderedRing.mul_nonneg hdiff ha'
  have product2 := OrderedRing.mul_nonneg hb hdiff'
  grind

/-- Lemma L in the case where the first two-dimensional part is semidefinite. -/
theorem lemmaLOrdered (p q r u a b c d : R)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hr : 0 ≤ r) (hu : 0 ≤ u)
    (balance : p+q=r+u) (ordered : p*q ≤ r*u)
    (ha : 0 ≤ a ∧ a ≤ 1) (hb : 0 ≤ b ∧ b ≤ 1)
    (hc : 0 ≤ c ∧ c ≤ 1) (hd : 0 ≤ d ∧ d ≤ 1) :
    0 ≤ cubeForm p q r u a b c d ∨ 0 ≤ cubeForm p q r u (1-a) (1-b) (1-c) (1-d) := by
  by_cases first : 0 ≤ cubeForm p q r u a b c d
  · exact Or.inl first
  by_cases second : 0 ≤ cubeForm p q r u (1-a) (1-b) (1-c) (1-d)
  · exact Or.inr second
  have negative : cubeForm p q r u a b c d < 0 := by grind
  have negative' : cubeForm p q r u (1-a) (1-b) (1-c) (1-d) < 0 := by grind
  let aa := quadratic (p*r) (p*q) (q*u) a b
  let aa' := quadratic (p*r) (p*q) (q*u) (1-a) (1-b)
  let bb := -quadratic (p*r) (r*u) (q*u) c d
  let bb' := -quadratic (p*r) (r*u) (q*u) (1-c) (1-d)
  let m := bilinear (p*r) (p*q) (q*u) a b (1-a) (1-b)
  let n := -bilinear (p*r) (r*u) (q*u) c d (1-c) (1-d)
  have haa : 0 ≤ aa := firstPartNonnegative p q r u a b hp hq hr hu ordered
  have haa' : 0 ≤ aa' := firstPartNonnegative p q r u (1-a) (1-b) hp hq hr hu ordered
  have hab : aa < bb := by grind [quadratic, cubeForm]
  have hab' : aa' < bb' := by grind [quadratic, cubeForm]
  have hm : m*m ≤ aa*aa' := firstPartCauchy p q r u a b (1-a) (1-b) hp hq ordered
  have reverse := secondPartReverseCauchy p q r u c d (1-c) (1-d) hr hu ordered
  have hn2 : bb*bb' ≤ n*n := by grind
  have polar1 := complementPolarization (p*r) (p*q) (q*u) a b
  have polar2 := complementPolarization (p*r) (r*u) (q*u) c d
  have total := totalAtOne p q r u balance
  have square := squareNonnegative (p-u)
  have htotal : 0 ≤ aa+aa'+2*m-(bb+bb'+2*n) := by
    grind
  have hpr := OrderedRing.mul_nonneg hp hr
  have hqu := OrderedRing.mul_nonneg hq hu
  have hru := OrderedRing.mul_nonneg hr hu
  have low : (p*r)*c*c+(q*u)*d*d < 2*(r*u)*c*d := by grind [quadratic]
  have low' : (p*r)*(1-c)*(1-c)+(q*u)*(1-d)*(1-d) < 2*(r*u)*(1-c)*(1-d) := by
    grind [quadratic]
  have cross := strictWeightedCross (p*r) (q*u) (r*u) c d (1-c) (1-d)
    hpr hqu hru hc.1 hd.1 (by grind) (by grind) low low'
  have hn : 0 < n := by grind [bilinear]
  exact False.elim (contradictionFromParts aa aa' bb bb' m n haa haa' hab hab' hm hn2 htotal hn)

-- Swapping the two parameter pairs and the two coordinate pairs preserves Q.
omit [LE R] [LT R] [LawfulOrderLT R] [IsLinearOrder R] [OrderedRing R] in
theorem cubeFormSwap (p q r u a b c d : R) :
    cubeForm p q r u a b c d = cubeForm r u p q c d a b := by
  grind [cubeForm]

/-- Unconditional Lemma L: only the original parameter and full-cube hypotheses. -/
theorem lemmaL (p q r u a b c d : R)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hr : 0 ≤ r) (hu : 0 ≤ u)
    (balance : p+q=r+u)
    (ha : 0 ≤ a ∧ a ≤ 1) (hb : 0 ≤ b ∧ b ≤ 1)
    (hc : 0 ≤ c ∧ c ≤ 1) (hd : 0 ≤ d ∧ d ≤ 1) :
    0 ≤ cubeForm p q r u a b c d ∨ 0 ≤ cubeForm p q r u (1-a) (1-b) (1-c) (1-d) := by
  by_cases ordered : p*q ≤ r*u
  · exact lemmaLOrdered p q r u a b c d hp hq hr hu balance ordered ha hb hc hd
  · have other : r*u ≤ p*q := by grind
    have result := lemmaLOrdered r u p q c d a b hr hu hp hq balance.symm other hc hd ha hb
    have swap := cubeFormSwap p q r u a b c d
    have swap' := cubeFormSwap p q r u (1-a) (1-b) (1-c) (1-d)
    grind

/-- The ordinary maximum, defined from the lawful linear order. -/
noncomputable def orderedMax (a b : R) : R := by
  classical
  exact if a ≤ b then b else a

/-- The maximum formulation of Lemma L, with no derived sign facts as premises. -/
theorem lemmaLMax (p q r u a b c d : R)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hr : 0 ≤ r) (hu : 0 ≤ u)
    (balance : p+q=r+u)
    (ha : 0 ≤ a ∧ a ≤ 1) (hb : 0 ≤ b ∧ b ≤ 1)
    (hc : 0 ≤ c ∧ c ≤ 1) (hd : 0 ≤ d ∧ d ≤ 1) :
    0 ≤ orderedMax (cubeForm p q r u a b c d) (cubeForm p q r u (1-a) (1-b) (1-c) (1-d)) := by
  have result := lemmaL p q r u a b c d hp hq hr hu balance ha hb hc hd
  unfold orderedMax
  split <;> grind

end MathCommons.OrthogonalCube

#print axioms MathCommons.OrthogonalCube.squareNonnegative
#print axioms MathCommons.OrthogonalCube.quadraticNonnegative
#print axioms MathCommons.OrthogonalCube.lagrangeIdentity
#print axioms MathCommons.OrthogonalCube.quadraticCauchy
#print axioms MathCommons.OrthogonalCube.quadraticReverseCauchy
#print axioms MathCommons.OrthogonalCube.firstPartNonnegative
#print axioms MathCommons.OrthogonalCube.firstPartCauchy
#print axioms MathCommons.OrthogonalCube.secondPartReverseCauchy
#print axioms MathCommons.OrthogonalCube.complementPolarization
#print axioms MathCommons.OrthogonalCube.totalAtOne

#print axioms MathCommons.OrthogonalCube.strictOrderOfSquares
#print axioms MathCommons.OrthogonalCube.weightedCauchy
#print axioms MathCommons.OrthogonalCube.strictWeightedCross
#print axioms MathCommons.OrthogonalCube.contradictionFromParts

#print axioms MathCommons.OrthogonalCube.lemmaLOrdered
#print axioms MathCommons.OrthogonalCube.cubeFormSwap
#print axioms MathCommons.OrthogonalCube.lemmaL
#print axioms MathCommons.OrthogonalCube.lemmaLMax
