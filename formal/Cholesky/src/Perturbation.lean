import src.Bounds

/-!
# Why a small perturbation cannot remove a negative direction

This file proves the last analytic implication in the breakdown argument.
The backward-error estimate is an explicit hypothesis. It is not a theorem here
about a floating-point implementation.
-/
namespace Cholesky

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- A perturbation that makes a negative direction nonnegative must be large. -/
theorem negative_direction_requires_error (A E : H →L[ℝ] H) (v : H) (b : ℝ)
    (hv : ‖v‖ = 1)
    (hnegative : inner ℝ v (A v) ≤ -b / 18)
    (hpositive : 0 ≤ inner ℝ v ((A + E) v)) : b / 18 ≤ ‖E‖ := by
  have hCS := real_inner_le_norm v (E v)
  have hop := E.le_opNorm v
  rw [hv, one_mul] at hCS
  rw [hv, mul_one] at hop
  rw [ContinuousLinearMap.add_apply, inner_add_right] at hpositive
  linarith

/-- The perturbation estimate used in the notes rules out a positive-semidefinite
result whenever A has a unit direction below -‖A‖/18. -/
theorem no_positive_perturbation (A E : H →L[ℝ] H) (v : H) (s : ℝ)
    (hA : 0 < ‖A‖) (hv : ‖v‖ = 1)
    (hnegative : inner ℝ v (A v) ≤ -‖A‖ / 18)
    (hs : 0 ≤ s) (hsmax : s ≤ 129 / 131072)
    (hbackward : ‖E‖ ≤ 9 * (s / (1 - s)) * (‖A‖ + ‖E‖)) :
    ¬ (∀ w : H, 0 ≤ inner ℝ w ((A + E) w)) := by
  intro hpositive
  exact breakdown_scalar_contradiction ‖A‖ ‖E‖ s hA (norm_nonneg _) hs hsmax hbackward
    (negative_direction_requires_error A E v ‖A‖ hv hnegative (hpositive v))

/-- In particular, no real Gram factor can satisfy these error estimates. -/
theorem no_gram_factor [CompleteSpace H] (A E : H →L[ℝ] H) (v : H) (s : ℝ)
    (hA : 0 < ‖A‖) (hv : ‖v‖ = 1)
    (hnegative : inner ℝ v (A v) ≤ -‖A‖ / 18)
    (hs : 0 ≤ s) (hsmax : s ≤ 129 / 131072)
    (hbackward : ‖E‖ ≤ 9 * (s / (1 - s)) * (‖A‖ + ‖E‖)) :
    ¬ ∃ R : H →L[ℝ] H, A + E = R.adjoint.comp R := by
  rintro ⟨R, hR⟩
  apply no_positive_perturbation A E v s hA hv hnegative hs hsmax hbackward
  intro w
  rw [hR, ContinuousLinearMap.comp_apply, ContinuousLinearMap.adjoint_inner_right]
  exact real_inner_self_nonneg

end Cholesky
