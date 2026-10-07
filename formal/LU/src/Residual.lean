import src.Norms

/-!
# Backward error and the residual

For `A x̂ = b` with computed `x̂`, a pair `(E, e)` *makes `x̂` exact* if
`(A + E) x̂ = b + e`. Its size is `τ = max {‖E‖/‖A‖, ‖e‖/‖b‖}` (`tau`), and the relative
backward error is the smallest size, `η(x̂) = min τ` (`eta`). With the residual
`r = b - A x̂` and `d = ‖A‖ ‖x̂‖ + ‖b‖`:

* `tauSet_nonempty`: such pairs exist, e.g. `(0, A x̂ - b)`;
* `residual_lower`: every pair has `‖r‖ ≤ τ d` (lower bound);
* `attaining_pair`: `E = (‖A‖/d) r x̂ᵀ/‖x̂‖`, `e = -(‖b‖/d) r` make `x̂` exact with
  `‖E‖/‖A‖ = ‖e‖/‖b‖ = ‖r‖/d`; for `x̂ = 0` every admissible pair has `e = -b`
  (`eta_zero_pairs`) and `(0, -b)` has `τ = 1`;
* `eta_isLeast`, `eta_eq`: **`η(x̂) = ‖r‖ / (‖A‖ ‖x̂‖ + ‖b‖)`**, and the minimum is attained;
* `etaFormula_le_one`: the scaled residual is at most 1.
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {n : ℕ}

/-- The size of a perturbation pair: `τ = max {‖E‖/‖A‖, ‖e‖/‖b‖}`. -/
noncomputable def tau (A E : Matrix (Fin n) (Fin n) ℝ) (b e : Fin n → ℝ) : ℝ :=
  max (norm2 E / norm2 A) (vnorm e / vnorm b)

/-- The sizes of all pairs `(E, e)` with `(A + E) x = b + e`. -/
def tauSet (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) : Set ℝ :=
  {τ | ∃ (E : Matrix (Fin n) (Fin n) ℝ) (e : Fin n → ℝ), (A + E) *ᵥ x = b + e ∧ τ = tau A E b e}

/-- The relative backward error `η(x) = min τ`. -/
noncomputable def eta (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) : ℝ :=
  sInf (tauSet A b x)

/-- The residual formula `‖b - A x‖ / (‖A‖ ‖x‖ + ‖b‖)`. -/
noncomputable def etaFormula (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) : ℝ :=
  vnorm (b - A *ᵥ x) / (norm2 A * vnorm x + vnorm b)

lemma tau_nonneg (A E : Matrix (Fin n) (Fin n) ℝ) (b e : Fin n → ℝ) : 0 ≤ tau A E b e :=
  le_trans (div_nonneg (norm2_nonneg _) (norm2_nonneg _)) (le_max_left _ _)

/-- Keeping `A` and replacing `b` by `A x` makes `x` exact. -/
lemma tauSet_nonempty (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) :
    (tauSet A b x).Nonempty :=
  ⟨_, 0, A *ᵥ x - b, by rw [add_zero]; abel, rfl⟩

lemma tauSet_bddBelow (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) :
    BddBelow (tauSet A b x) :=
  ⟨0, fun _ ⟨E, e, _, hτ⟩ => hτ ▸ tau_nonneg A E b e⟩

/-- The backward error is at most the size of any pair that makes `x` exact. -/
lemma eta_le_tau {A E : Matrix (Fin n) (Fin n) ℝ} {b e x : Fin n → ℝ}
    (h : (A + E) *ᵥ x = b + e) : eta A b x ≤ tau A E b e :=
  csInf_le (tauSet_bddBelow A b x) ⟨E, e, h, rfl⟩

/-- **Lower bound.** The perturbed equation gives `r = E x - e`, hence `‖r‖ ≤ τ d`. -/
theorem residual_lower {A E : Matrix (Fin n) (Fin n) ℝ} {b e x : Fin n → ℝ}
    (hA : 0 < norm2 A) (hb : b ≠ 0) (h : (A + E) *ᵥ x = b + e) :
    vnorm (b - A *ᵥ x) ≤ tau A E b e * (norm2 A * vnorm x + vnorm b) := by
  have hbpos := vnorm_pos hb
  have hr : b - A *ᵥ x = E *ᵥ x - e := by
    have : b = (A + E) *ᵥ x - e := by rw [h]; abel
    rw [this, Matrix.add_mulVec]; abel
  have hE : norm2 E ≤ tau A E b e * norm2 A := by
    rw [← div_le_iff₀ hA]; exact le_max_left _ _
  have he : vnorm e ≤ tau A E b e * vnorm b := by
    rw [← div_le_iff₀ hbpos]; exact le_max_right _ _
  rw [hr]
  calc vnorm (E *ᵥ x - e) ≤ vnorm (E *ᵥ x) + vnorm e := vnorm_sub_le _ _
    _ ≤ norm2 E * vnorm x + vnorm e := add_le_add (vnorm_mulVec_le _ _) le_rfl
    _ ≤ tau A E b e * norm2 A * vnorm x + tau A E b e * vnorm b :=
        add_le_add (mul_le_mul_of_nonneg_right hE (vnorm_nonneg x)) he
    _ = tau A E b e * (norm2 A * vnorm x + vnorm b) := by ring

/-- **Attaining the bound.** For `x ≠ 0` the pair `E = (‖A‖/d) r xᵀ/‖x‖`, `e = -(‖b‖/d) r`
makes `x` exact and has `‖E‖/‖A‖ = ‖e‖/‖b‖ = ‖r‖/d`. -/
theorem attaining_pair (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) (hA : 0 < norm2 A)
    (hb : b ≠ 0) (hx : x ≠ 0) :
    let r := b - A *ᵥ x
    let d := norm2 A * vnorm x + vnorm b
    let E := (norm2 A / d / vnorm x) • vecMulVec r x
    let e := -(vnorm b / d) • r
    (A + E) *ᵥ x = b + e ∧ norm2 E / norm2 A = vnorm r / d ∧ vnorm e / vnorm b = vnorm r / d := by
  intro r d E e
  have hbpos := vnorm_pos hb
  have hxpos := vnorm_pos hx
  have hd : 0 < d := by
    have := mul_nonneg (norm2_nonneg A) (vnorm_nonneg x); simp only [d]; linarith
  refine ⟨?_, ?_, ?_⟩
  · have hb' : b = A *ᵥ x + r := by simp only [r]; abel
    have key : norm2 A / d / vnorm x * vnorm x ^ 2 = 1 - vnorm b / d := by
      field_simp
      simp only [d]
      ring
    simp only [E, e]
    rw [Matrix.add_mulVec, Matrix.smul_mulVec, vecMulVec_mulVec', ← vnorm_sq, smul_smul, key]
    conv_rhs => arg 1; rw [hb']
    rw [sub_smul, one_smul, neg_smul]
    abel
  · simp only [E]
    rw [norm2_smul, norm2_vecMulVec, abs_of_nonneg (by positivity)]
    field_simp
  · simp only [e]
    rw [vnorm_smul, abs_neg, abs_of_nonneg (by positivity)]
    field_simp

/-- For `x = 0`, every pair that makes `x` exact has `e = -b`. -/
lemma eta_zero_pairs {A E : Matrix (Fin n) (Fin n) ℝ} {b e : Fin n → ℝ}
    (h : (A + E) *ᵥ (0 : Fin n → ℝ) = b + e) : e = -b := by
  rw [Matrix.mulVec_zero] at h
  exact (neg_eq_of_add_eq_zero_right h.symm).symm

/-- **The residual formula for the backward error**, with the minimum attained. -/
theorem eta_isLeast (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) (hA : 0 < norm2 A)
    (hb : b ≠ 0) : IsLeast (tauSet A b x) (etaFormula A b x) := by
  have hbpos := vnorm_pos hb
  have hd : 0 < norm2 A * vnorm x + vnorm b := by
    have := mul_nonneg (norm2_nonneg A) (vnorm_nonneg x); linarith
  constructor
  · by_cases hx : x = 0
    · -- `E = 0`, `e = -b`: `τ = 1 = ‖r‖ / d`
      refine ⟨0, -b, by rw [hx, Matrix.mulVec_zero, add_neg_cancel], ?_⟩
      have h1 : etaFormula A b x = 1 := by
        rw [etaFormula, hx, Matrix.mulVec_zero, sub_zero, vnorm_zero, mul_zero, zero_add,
          div_self (ne_of_gt hbpos)]
      rw [h1, tau, norm2_zero, zero_div, vnorm_neg, div_self (ne_of_gt hbpos)]
      simp
    · obtain ⟨h1, h2, h3⟩ := attaining_pair A b x hA hb hx
      exact ⟨_, _, h1, by rw [tau, h2, h3, max_self]; rfl⟩
  · rintro τ ⟨E, e, h, rfl⟩
    rw [etaFormula, div_le_iff₀ hd]
    exact residual_lower hA hb h

theorem eta_eq (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) (hA : 0 < norm2 A)
    (hb : b ≠ 0) : eta A b x = etaFormula A b x :=
  (eta_isLeast A b x hA hb).csInf_eq

/-- `η(0) = 1`. -/
theorem eta_zero (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) (hA : 0 < norm2 A)
    (hb : b ≠ 0) : eta A b 0 = 1 := by
  rw [eta_eq A b 0 hA hb, etaFormula, Matrix.mulVec_zero, sub_zero, vnorm_zero, mul_zero,
    zero_add, div_self (ne_of_gt (vnorm_pos hb))]

/-- The scaled residual is at most 1. -/
theorem etaFormula_le_one (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) :
    etaFormula A b x ≤ 1 := by
  unfold etaFormula
  have hd : 0 ≤ norm2 A * vnorm x + vnorm b :=
    add_nonneg (mul_nonneg (norm2_nonneg A) (vnorm_nonneg x)) (vnorm_nonneg b)
  rcases eq_or_lt_of_le hd with h | h
  · rw [← h, div_zero]; exact zero_le_one
  · rw [div_le_one h]
    calc vnorm (b - A *ᵥ x) ≤ vnorm b + vnorm (A *ᵥ x) := vnorm_sub_le _ _
      _ ≤ vnorm b + norm2 A * vnorm x := add_le_add le_rfl (vnorm_mulVec_le _ _)
      _ = norm2 A * vnorm x + vnorm b := by ring

end LU
