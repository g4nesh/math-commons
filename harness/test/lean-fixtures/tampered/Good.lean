namespace Fixtures.Good

/-- `p` is prime: at least 2, and its only divisors are 1 and itself. -/
def IsPrime (p : Nat) : Prop := 2 ≤ p ∧ ∀ d : Nat, d ∣ p → d = 1 ∨ d = p

/-- Silently weakened: only n ≤ 100. -/
def Claim : Prop := ∀ n : Nat, n ≤ 100 → IsPrime (n ^ 2 + n + 41)

theorem value_at_forty : 40 ^ 2 + 40 + 41 = 41 * 41 := by decide

theorem main : ¬ Claim := by
  intro h
  have h40 : 2 ≤ 40 ^ 2 + 40 + 41 ∧ ∀ d : Nat, d ∣ 40 ^ 2 + 40 + 41 → d = 1 ∨ d = 40 ^ 2 + 40 + 41 := h 40 (by decide)
  exact (h40.2 41 ⟨41, value_at_forty⟩).elim (fun h1 => absurd h1 (by decide)) (fun h2 => absurd h2 (by decide))

end Fixtures.Good
