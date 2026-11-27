namespace Fixtures.CustomAxiom

/-- `p` is prime: at least 2, and its only divisors are 1 and itself. -/
def IsPrime (p : Nat) : Prop := 2 ≤ p ∧ ∀ d : Nat, d ∣ p → d = 1 ∨ d = p

/-- Euler's polynomial read as a universal claim. -/
def Claim : Prop := ∀ n : Nat, IsPrime (n ^ 2 + n + 41)

axiom euler_fails : ¬ Claim

theorem main : ¬ Claim := euler_fails

end Fixtures.CustomAxiom
