import Mathlib

/-! # Selecting a block at a fixed precision

The natural-number arithmetic behind the largest admissible choice of `q`.
This file does not change the precision as the matrix is padded to larger sizes.
-/
namespace Cholesky

/-- Largest q satisfying 3q + 8 ≤ p and p + q even (for p ≥ 20). -/
def largestQ (p : ℕ) : ℕ :=
  2 * ((p - 8 - 3 * (p % 2)) / 6) + p % 2

theorem largestQ_properties (p : ℕ) (hp : 20 ≤ p) :
    3 ≤ largestQ p ∧
    3 * largestQ p + 8 ≤ p ∧
    (p + largestQ p) % 2 = 0 ∧
    (p - 3 * largestQ p = 8 ∨ p - 3 * largestQ p = 10 ∨
      p - 3 * largestQ p = 12) := by
  unfold largestQ
  omega

theorem largestQ_maximal (p q : ℕ)
    (hsize : 3 * q + 8 ≤ p) (hparity : (p + q) % 2 = 0) :
    q ≤ largestQ p := by
  unfold largestQ
  omega

theorem precision_at_least_17 (p q : ℕ) (hq : 3 ≤ q)
    (hsize : 3 * q + 8 ≤ p) : 17 ≤ p := by omega

end Cholesky
