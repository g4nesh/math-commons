import Std

/-!
# Claude's unit-cube quadratic lemma over an arbitrary ordered field

The domain is a `Lean.Grind.Field` with a lawful linear order and compatible
strict ordered-ring operations. These are Lean's standard algebra/order laws,
not extra sign assumptions for the lemma. No completeness or square roots are
used. The matrix, SO(3), SO(4), Hodge and PMSH reductions are separate claims.
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
