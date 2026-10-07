import Mathlib

/-!
# The constants `γ_m` of rounding-error analysis

The notes write `γ_m = m u / (1 - m u)`. This file collects the elementary
inequalities used when rounding errors are accumulated.
-/
namespace Cholesky

/-- `gamma u m = m u / (1 - m u)`, written `γ_m` in the notes. -/
noncomputable def gamma (u : ℝ) (m : ℕ) : ℝ := m * u / (1 - m * u)

variable {u : ℝ}

lemma gamma_zero (u : ℝ) : gamma u 0 = 0 := by simp [gamma]

lemma gamma_nonneg (hu : 0 ≤ u) {m : ℕ} (hm : m * u < 1) : 0 ≤ gamma u m := by
  unfold gamma
  exact div_nonneg (by positivity) (by linarith)

lemma one_add_gamma {m : ℕ} (hm : m * u < 1) : 1 + gamma u m = 1 / (1 - m * u) := by
  unfold gamma
  have : (1 : ℝ) - m * u ≠ 0 := by linarith
  field_simp
  ring

lemma gamma_mono (hu : 0 ≤ u) {a b : ℕ} (hab : a ≤ b) (hb : b * u < 1) :
    gamma u a ≤ gamma u b := by
  have hab' : (a : ℝ) * u ≤ b * u := mul_le_mul_of_nonneg_right (by exact_mod_cast hab) hu
  have ha : (a : ℝ) * u < 1 := lt_of_le_of_lt hab' hb
  have h1 := one_add_gamma (u := u) ha
  have h2 := one_add_gamma (u := u) hb
  have : 1 / (1 - a * u) ≤ 1 / (1 - b * u) :=
    one_div_le_one_div_of_le (by linarith) (by linarith)
  linarith

/-- The basic accumulation inequality `γ_a + γ_b + γ_a γ_b ≤ γ_{a+b}`. -/
lemma gamma_add_le (hu : 0 ≤ u) {a b : ℕ} (hab : ((a + b : ℕ) : ℝ) * u < 1) :
    gamma u a + gamma u b + gamma u a * gamma u b ≤ gamma u (a + b) := by
  push_cast at hab
  have ha : (a : ℝ) * u < 1 := by nlinarith [mul_nonneg (Nat.cast_nonneg b : (0 : ℝ) ≤ b) hu]
  have hb : (b : ℝ) * u < 1 := by nlinarith [mul_nonneg (Nat.cast_nonneg a : (0 : ℝ) ≤ a) hu]
  have hab' : ((a + b : ℕ) : ℝ) * u < 1 := by push_cast; exact hab
  have h1 := one_add_gamma (u := u) ha
  have h2 := one_add_gamma (u := u) hb
  have h3 := one_add_gamma (u := u) hab'
  have hkey : (1 + gamma u a) * (1 + gamma u b) ≤ 1 + gamma u (a + b) := by
    rw [h1, h2, h3]
    push_cast
    rw [div_mul_div_comm, one_mul]
    apply one_div_le_one_div_of_le (by linarith)
    nlinarith [mul_nonneg (mul_nonneg (Nat.cast_nonneg a : (0 : ℝ) ≤ a) hu)
      (mul_nonneg (Nat.cast_nonneg b : (0 : ℝ) ≤ b) hu)]
  nlinarith

lemma le_gamma_one (hu : 0 ≤ u) (hu1 : u < 1) : u ≤ gamma u 1 := by
  unfold gamma
  simp only [Nat.cast_one, one_mul]
  rw [le_div_iff₀ (by linarith)]
  nlinarith

lemma le_gamma (hu : 0 ≤ u) {m : ℕ} (hm1 : 1 ≤ m) (hm : m * u < 1) : u ≤ gamma u m :=
  le_trans (le_gamma_one hu (by
    have : (1 : ℝ) * u ≤ m * u := mul_le_mul_of_nonneg_right (by exact_mod_cast hm1) hu
    linarith)) (gamma_mono hu hm1 hm)

/-- One relative error `δ` contributes `|δ / (1 + δ)| ≤ γ_1`. -/
lemma abs_div_one_add_le (hu1 : u < 1) {δ : ℝ} (hδ : |δ| ≤ u) :
    |δ / (1 + δ)| ≤ gamma u 1 := by
  have hb := abs_le.mp hδ
  have hpos : 0 < 1 + δ := by linarith
  unfold gamma
  simp only [Nat.cast_one, one_mul]
  rw [abs_div, abs_of_pos hpos]
  apply div_le_div₀ (by linarith [abs_nonneg δ]) hδ (by linarith) (by linarith)

/-- Inverting `1 + δ` costs at most `γ_1`. -/
lemma abs_inv_one_add_sub_one_le (hu1 : u < 1) {δ : ℝ} (hδ : |δ| ≤ u) :
    |1 / (1 + δ) - 1| ≤ gamma u 1 := by
  have hb := abs_le.mp hδ
  have hpos : 0 < 1 + δ := by linarith
  have : 1 / (1 + δ) - 1 = -(δ / (1 + δ)) := by field_simp; ring
  rw [this, abs_neg]
  exact abs_div_one_add_le hu1 hδ

/-- A rounded square root contributes `|1 - 1 / (1 + δ)^2| ≤ γ_2`. -/
lemma abs_one_sub_inv_sq_le (hu : 0 ≤ u) (hu2 : 2 * u < 1) {δ : ℝ} (hδ : |δ| ≤ u) :
    |1 - 1 / (1 + δ) ^ 2| ≤ gamma u 2 := by
  have hb := abs_le.mp hδ
  have hpos : 0 < 1 + δ := by linarith
  unfold gamma
  have hden : 0 < 1 - (2 : ℝ) * u := by linarith
  push_cast
  rw [abs_le]
  constructor
  · -- `1 - 1/(1+δ)^2 ≥ 1 - 1/(1-u)^2 ≥ -2u/(1-2u)`
    rw [show (1 : ℝ) - 1 / (1 + δ) ^ 2 = ((1 + δ) ^ 2 - 1) / (1 + δ) ^ 2 by field_simp]
    rw [le_div_iff₀ (by positivity)]
    have : -(2 * u / (1 - 2 * u)) * (1 + δ) ^ 2 = -(2 * u) * (1 + δ) ^ 2 / (1 - 2 * u) := by
      ring
    rw [this, div_le_iff₀ hden]
    nlinarith [mul_nonneg hu (mul_nonneg hu hu), sq_nonneg (u + δ), mul_nonneg hu
      (show 0 ≤ u + δ by linarith), mul_nonneg (mul_nonneg hu hu) (show 0 ≤ u + δ by linarith)]
  · rw [show (1 : ℝ) - 1 / (1 + δ) ^ 2 = ((1 + δ) ^ 2 - 1) / (1 + δ) ^ 2 by field_simp]
    rw [div_le_div_iff₀ (by positivity) hden]
    nlinarith [mul_nonneg hu (mul_nonneg hu hu), mul_nonneg hu (show 0 ≤ u - δ by linarith),
      mul_nonneg (mul_nonneg hu hu) (show 0 ≤ u - δ by linarith), sq_nonneg δ,
      mul_nonneg (sq_nonneg δ) hu]

end Cholesky
