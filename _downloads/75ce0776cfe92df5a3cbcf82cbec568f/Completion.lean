import Mathlib

/-!
# Meinguet's product bound and the constant 1.1

This file proves the real inequality needed for the corollary in the notes.
It does not assume or formalize a floating-point Cholesky implementation.
-/
namespace Cholesky
open Finset Real

lemma sqrt_step (j : ℕ) :
    3 * Real.sqrt j ≤
      2 * (((j : ℝ) + 1) * Real.sqrt ((j : ℝ) + 1) - j * Real.sqrt j) := by
  have ha := Real.sqrt_nonneg (j : ℝ)
  have hb := Real.sqrt_nonneg ((j : ℝ) + 1)
  have ha2 := Real.sq_sqrt (show 0 ≤ (j : ℝ) by positivity)
  have hb2 := Real.sq_sqrt (show 0 ≤ (j : ℝ) + 1 by positivity)
  have hab : Real.sqrt (j : ℝ) ≤ Real.sqrt ((j : ℝ) + 1) :=
    Real.sqrt_le_sqrt (by linarith)
  have hprod := mul_nonneg (sub_nonneg.mpr hab)
    (show 0 ≤ (2 * ((j : ℝ) + 1) * Real.sqrt ((j : ℝ) + 1) +
      (2 * (j : ℝ) + 3) * Real.sqrt j) by positivity)
  nlinarith [sq_nonneg (Real.sqrt ((j : ℝ) + 1) - Real.sqrt j)]

lemma sum_sqrt_le (n : ℕ) :
    (∑ j ∈ range n, Real.sqrt (j : ℝ)) ≤ (2 / 3 : ℝ) * n * Real.sqrt n := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ]
    push_cast
    have hstep := sqrt_step n
    linarith

noncomputable def meinguetMu (n : ℕ) (u : ℝ) : ℝ :=
  (∏ j ∈ range (n - 1), (1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 : ℕ) + 2))) - 1

lemma log_chord (κ : ℝ) (hκ : 1 ≤ κ) :
    Real.log 2 / κ ≤ Real.log (1 + 1 / κ) := by
  have hkpos : 0 < κ := by linarith
  have ht0 : 0 ≤ 1 / κ := by positivity
  have ht1 : 1 / κ ≤ 1 := (div_le_one hkpos).2 hκ
  have h := strictConcaveOn_log_Ioi.concaveOn.2
    (show (1 : ℝ) ∈ Set.Ioi 0 by norm_num)
    (show (2 : ℝ) ∈ Set.Ioi 0 by norm_num)
    (show 0 ≤ 1 - 1 / κ by linarith) ht0 (by ring : (1 - 1 / κ) + 1 / κ = 1)
  simp only [smul_eq_mul, Real.log_one, mul_zero, zero_add] at h
  rw [show (1 - 1 / κ) * 1 + 1 / κ * 2 = 1 + 1 / κ by ring] at h
  simpa [div_eq_mul_inv, mul_comm] using h

lemma sum_update_weights (n : ℕ) (hn : 576 ≤ n) :
    (∑ j ∈ range (n - 1), (Real.sqrt (j + 1 : ℕ) + 2)) <
      (3 / 4 : ℝ) * n * Real.sqrt n := by
  have hn1 : 1 ≤ n := by omega
  have hnR : (576 : ℝ) ≤ n := by exact_mod_cast hn
  have hsqrt : (24 : ℝ) ≤ Real.sqrt n := by
    have h := Real.sqrt_le_sqrt hnR
    norm_num at h
    exact h
  have hshift := sum_range_succ' (fun j : ℕ => Real.sqrt (j : ℝ)) (n - 1)
  rw [Nat.sub_add_cancel hn1] at hshift
  simp only [Nat.cast_zero, Real.sqrt_zero, add_zero] at hshift
  have hsum := sum_sqrt_le n
  rw [hshift] at hsum
  have hnsub : ((n - 1 : ℕ) : ℝ) = n - 1 := by simp [Nat.cast_sub hn1]
  simp only [sum_add_distrib, sum_const, card_range, nsmul_eq_mul]
  rw [hnsub]
  nlinarith [mul_nonneg (show 0 ≤ (n : ℝ) by positivity) (sub_nonneg.mpr hsqrt)]

/-- The finite-precision scalar implication behind the 1.1 corollary.
Here `n * sqrt n` is exactly the real number denoted by `n^(3/2)` in the notes. -/
theorem meinguet_condition (n : ℕ) (u κ : ℝ)
    (hn : 576 ≤ n) (hu : 0 < u) (humax : u ≤ 1 / 131072)
    (hκ : 1 ≤ κ)
    (hbudget : (11 / 10 : ℝ) * (n * Real.sqrt n) * u * κ ≤ 1) :
    κ * meinguetMu n u < 1 := by
  have hkpos : 0 < κ := by linarith
  have hd : 0 < 1 - 5 * u := by linarith
  let f : ℕ → ℝ := fun j => u / (1 - 5 * u) * (Real.sqrt (j + 1 : ℕ) + 2)
  have hf : ∀ j, 0 ≤ f j := by intro j; dsimp [f]; positivity
  have hp : 0 < ∏ j ∈ range (n - 1), (1 + f j) :=
    prod_pos (fun j _ => by linarith [hf j])
  have hlog : Real.log (∏ j ∈ range (n - 1), (1 + f j)) ≤
      ∑ j ∈ range (n - 1), f j := by
    rw [Real.log_prod _ _ (fun j _ => ne_of_gt (by linarith [hf j]))]
    apply sum_le_sum
    intro j _
    have h := Real.log_le_sub_one_of_pos (show 0 < 1 + f j by linarith [hf j])
    linarith
  have hsum : (∑ j ∈ range (n - 1), f j) <
      (u / (1 - 5 * u)) * ((3 / 4 : ℝ) * n * Real.sqrt n) := by
    simpa only [f, ← mul_sum] using
      mul_lt_mul_of_pos_left (sum_update_weights n hn) (div_pos hu hd)
  have hprodBound : κ * Real.log (∏ j ∈ range (n - 1), (1 + f j)) <
      3 / ((22 / 5 : ℝ) * (1 - 5 * u)) := by
    have h := mul_lt_mul_of_pos_left (lt_of_le_of_lt hlog hsum) hkpos
    have heq : κ * (u / (1 - 5 * u) * ((3 / 4 : ℝ) * n * Real.sqrt n)) =
        ((3 / 4 : ℝ) * (n * Real.sqrt n) * u * κ) / (1 - 5 * u) := by ring
    rw [heq] at h
    have hnum : (3 / 4 : ℝ) * (n * Real.sqrt n) * u * κ ≤ 15 / 22 := by
      nlinarith
    exact lt_of_lt_of_le h (by
      rw [show 3 / ((22 / 5 : ℝ) * (1 - 5 * u)) = (15 / 22 : ℝ) / (1 - 5 * u) by field_simp; ring]
      exact div_le_div_of_nonneg_right hnum (le_of_lt hd))
  have hconstant : 3 / ((22 / 5 : ℝ) * (1 - 5 * u)) < Real.log 2 := by
    apply (div_lt_iff₀ (by positivity : 0 < (22 / 5 : ℝ) * (1 - 5 * u))).2
    have hlog2 := mul_lt_mul_of_pos_right Real.log_two_gt_d9
      (show 0 < (22 / 5 : ℝ) * (1 - 5 * u) by positivity)
    nlinarith
  have hcompare : Real.log (∏ j ∈ range (n - 1), (1 + f j)) <
      Real.log (1 + 1 / κ) := by
    apply lt_of_lt_of_le _ (log_chord κ hκ)
    apply (lt_div_iff₀ hkpos).2
    nlinarith [lt_trans hprodBound hconstant]
  have hP := (Real.log_lt_log_iff hp (by positivity)).mp hcompare
  change κ * ((∏ j ∈ range (n - 1), (1 + f j)) - 1) < 1
  have h := mul_lt_mul_of_pos_left hP hkpos
  have hinv : κ * (1 / κ) = 1 := by field_simp
  nlinarith

end Cholesky
