import Mathlib

/-!
# Scalar bounds from the optional Cholesky notes

All statements in this file are unconditional real-arithmetic theorems under their
explicit scalar hypotheses. They do not identify these expressions with the
condition number or the output of a floating-point algorithm.
-/
namespace Cholesky

noncomputable def conditionUpper (x : ℝ) : ℝ :=
  (17 / 2 + 39 * x) * (1 + 17 / (72 * x))

theorem conditionUpper_scaled (x : ℝ) (hx : x ≠ 0) :
    conditionUpper x * (27 * x) =
      27 * (17 / 2 + 39 * x) * (x + 17 / 72) := by
  unfold conditionUpper
  field_simp

theorem conditionUpper_lt_57 (x : ℝ) (hx : 0 < x) (hxmax : x ≤ 1 / 256) :
    conditionUpper x < 57 / (27 * x) := by
  apply (lt_div_iff₀ (by positivity : 0 < 27 * x)).2
  rw [conditionUpper_scaled x (ne_of_gt hx)]
  have h₁ : 0 ≤ 27 * (17 / 2 + 39 * x) := by positivity
  have h₂ : 0 ≤ x + 17 / 72 := by positivity
  have h := mul_le_mul
    (show 27 * (17 / 2 + 39 * x) ≤ 27 * (17 / 2 + 39 * (1 / 256 : ℝ)) by linarith)
    (show x + 17 / 72 ≤ (1 / 256 : ℝ) + 17 / 72 by linarith)
    h₂ (by norm_num)
  norm_num at h
  linarith

theorem conditionUpper_fixed_precision (x : ℝ)
    (hx : x = 1 / 256 ∨ x = 1 / 1024 ∨ x = 1 / 4096) :
    conditionUpper x ≤ 303691615 / 36864 ∧ conditionUpper x < 8239 := by
  rcases hx with rfl | rfl | rfl <;> norm_num [conditionUpper]

theorem conditionUpper_binary64 :
    conditionUpper (1 / 256) = 1224895 / 2304 ∧
    conditionUpper (1 / 256) < 532 := by
  norm_num [conditionUpper]

theorem family_lower_bound (x : ℝ) (hx : 0 < x) (hxmax : x ≤ 1 / 256) :
    (1156 : ℝ) / 3 ≤ 289 / (192 * x) := by
  apply (le_div_iff₀ (by positivity : 0 < 192 * x)).2
  nlinarith

theorem negative_eigenvalue_margin (k : ℝ) (hk : 64 ≤ k) :
    (1 : ℝ) / 18 < (k - 9) / (15 * k + 9) := by
  apply (lt_div_iff₀ (by linarith : 0 < 15 * k + 9)).2
  linarith

theorem update_budget (k u : ℝ) (hku : k * u ≤ 1 / 2048)
    (hu : u ≤ 1 / 131072) : (2 * k + 1) * u ≤ 129 / 131072 := by
  nlinarith

theorem perturbation_fraction (s : ℝ) (_hs : 0 ≤ s)
    (hsmax : s ≤ 129 / 131072) :
    9 * s / (1 - 10 * s) ≤ 1161 / 129782 ∧
    9 * s / (1 - 10 * s) < 1 / 100 := by
  have hd : 0 < 1 - 10 * s := by linarith
  have h : 9 * s / (1 - 10 * s) ≤ 1161 / 129782 := by
    apply (div_le_iff₀ hd).2
    linarith
  exact ⟨h, lt_of_le_of_lt h (by norm_num)⟩

/-- The norm contradiction used after the sign-pattern estimate. -/
theorem breakdown_scalar_contradiction (b e s : ℝ) (hb : 0 < b)
    (_he : 0 ≤ e) (hs : 0 ≤ s) (hsmax : s ≤ 129 / 131072)
    (hbackward : e ≤ 9 * (s / (1 - s)) * (b + e))
    (hnecessary : b / 18 ≤ e) : False := by
  have hd : 0 < 1 - s := by linarith
  have hmul := (le_div_iff₀ hd).mp
    (show e ≤ (9 * s * (b + e)) / (1 - s) by
      convert hbackward using 1
      ring)
  have hd10 : 0 < 1 - 10 * s := by linarith
  have herr : e ≤ (9 * s / (1 - 10 * s)) * b := by
    rw [div_mul_eq_mul_div]
    apply (le_div_iff₀ hd10).2
    nlinarith
  have hsmall := (perturbation_fraction s hs hsmax).2
  nlinarith

end Cholesky
