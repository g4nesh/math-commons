namespace Fixtures.Weakened

/-- `p` is prime: at least 2, and its only divisors are 1 and itself. -/
def IsPrime (p : Nat) : Prop := 2 ≤ p ∧ ∀ d : Nat, d ∣ p → d = 1 ∨ d = p

/-- Euler's polynomial read as a universal claim. -/
def Claim : Prop := ∀ n : Nat, IsPrime (n ^ 2 + n + 41)

/-- Proves only the n = 40 instance, not the negation of the claim. -/
theorem main : ¬ IsPrime (40 ^ 2 + 40 + 41) := by
  intro h40
  have h' : 2 ≤ 40 ^ 2 + 40 + 41 ∧ ∀ d : Nat, d ∣ 40 ^ 2 + 40 + 41 → d = 1 ∨ d = 40 ^ 2 + 40 + 41 := h40
  exact (h'.2 41 ⟨41, by decide⟩).elim (fun h1 => absurd h1 (by decide)) (fun h2 => absurd h2 (by decide))

end Fixtures.Weakened
