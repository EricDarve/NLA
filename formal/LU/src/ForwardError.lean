import src.Residual

/-!
# From backward error to forward error

Throughout, `A` is nonsingular, `A x = b` with `b ≠ 0`, and `x̂` is a computed solution;
`κ(A) = ‖A‖ ‖A⁻¹‖` (`cond2`).

* `error_eq`: `δx = x̂ - x = -A⁻¹ r` with `r = b - A x̂`.
* `relErr_le`: `‖δx‖/‖x‖ ≤ κ η (1 + ‖x̂‖/‖x‖)` (every step of the first-order derivation is
  rigorous); `relErr_le_of_norm_le`: if `‖x̂‖ ≤ (1 + t) ‖x‖` the bound is `(2 + t) κ η`,
  which is the first-order estimate `2 κ η` when `t = 0`.
* `banach`: **Banach lemma.** In a complete normed ring with `‖1‖ = 1` (every induced
  matrix norm), `‖X‖ < 1` gives that `I + X` is invertible with
  `‖(I + X)⁻¹‖ ≤ 1/(1 - ‖X‖)`; the inverse is the Neumann series `∑ (-X)^j`
  (`neumann_partial`, `banach_inverse`). `banach_matrix` is the case of `‖·‖₂`.
* `perturbation_bound`: **theorem `thm:perturbation_bound`.**
* `forward_error_bound`: if `κ η < 1` then `‖x̂ - x‖/‖x‖ ≤ 2 κ η/(1 - κ η)`.
-/
open scoped Matrix.Norms.L2Operator
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {n : ℕ}

/-! ## The residual and the error -/

/-- `δx = -A⁻¹ r`. -/
theorem error_eq {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) {x b xh : Fin n → ℝ}
    (hx : A *ᵥ x = b) : xh - x = -(A⁻¹ *ᵥ (b - A *ᵥ xh)) := by
  rw [Matrix.mulVec_sub, ← hx, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec,
    Matrix.nonsing_inv_mul _ hA, Matrix.one_mulVec, Matrix.one_mulVec, neg_sub]

lemma vnorm_le_of_mulVec {A : Matrix (Fin n) (Fin n) ℝ} {x b : Fin n → ℝ} (hx : A *ᵥ x = b) :
    vnorm b ≤ norm2 A * vnorm x := by rw [← hx]; exact vnorm_mulVec_le A x

lemma ne_zero_of_mulVec {A : Matrix (Fin n) (Fin n) ℝ} {x b : Fin n → ℝ} (hx : A *ᵥ x = b)
    (hb : b ≠ 0) : x ≠ 0 := by
  rintro rfl; exact hb (by rw [← hx, Matrix.mulVec_zero])

/-- The rigorous form of the first-order derivation: `‖δx‖/‖x‖ ≤ κ η (1 + ‖x̂‖/‖x‖)`. -/
theorem relErr_le {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) {x b xh : Fin n → ℝ}
    (hx : A *ᵥ x = b) (hb : b ≠ 0) :
    vnorm (xh - x) / vnorm x ≤ cond2 A * etaFormula A b xh * (1 + vnorm xh / vnorm x) := by
  have hxpos := vnorm_pos (ne_zero_of_mulVec hx hb)
  have hbpos := vnorm_pos hb
  have hd : 0 < norm2 A * vnorm xh + vnorm b := by
    have := mul_nonneg (norm2_nonneg A) (vnorm_nonneg xh); linarith
  have hr : vnorm (b - A *ᵥ xh) = etaFormula A b xh * (norm2 A * vnorm xh + vnorm b) := by
    rw [etaFormula, div_mul_cancel₀ _ (ne_of_gt hd)]
  have h1 : vnorm (xh - x) ≤ norm2 A⁻¹ * vnorm (b - A *ᵥ xh) := by
    rw [error_eq hA hx, vnorm_neg]; exact vnorm_mulVec_le _ _
  have hbx := vnorm_le_of_mulVec hx
  have hη : 0 ≤ etaFormula A b xh := div_nonneg (vnorm_nonneg _) hd.le
  rw [div_le_iff₀ hxpos]
  calc vnorm (xh - x) ≤ norm2 A⁻¹ * (etaFormula A b xh * (norm2 A * vnorm xh + vnorm b)) := by
        rw [← hr]; exact h1
    _ ≤ norm2 A⁻¹ * (etaFormula A b xh * (norm2 A * vnorm xh + norm2 A * vnorm x)) := by
        gcongr
        exact norm2_nonneg _
    _ = cond2 A * etaFormula A b xh * (1 + vnorm xh / vnorm x) * vnorm x := by
        unfold cond2; field_simp; ring

/-- If `‖x̂‖ ≤ (1 + t) ‖x‖`, then `‖δx‖/‖x‖ ≤ (2 + t) κ η`; for `‖x̂‖ ≈ ‖x‖` this is the
first-order estimate `2 κ η`. -/
theorem relErr_le_of_norm_le {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    {x b xh : Fin n → ℝ} (hx : A *ᵥ x = b) (hb : b ≠ 0) {t : ℝ}
    (ht : vnorm xh ≤ (1 + t) * vnorm x) :
    vnorm (xh - x) / vnorm x ≤ (2 + t) * cond2 A * etaFormula A b xh := by
  have hxpos := vnorm_pos (ne_zero_of_mulVec hx hb)
  have hκη : 0 ≤ cond2 A * etaFormula A b xh :=
    mul_nonneg (cond2_nonneg A) (div_nonneg (vnorm_nonneg _)
      (add_nonneg (mul_nonneg (norm2_nonneg _) (vnorm_nonneg _)) (vnorm_nonneg _)))
  have : vnorm xh / vnorm x ≤ 1 + t := by rwa [div_le_iff₀ hxpos]
  calc vnorm (xh - x) / vnorm x ≤ cond2 A * etaFormula A b xh * (1 + vnorm xh / vnorm x) :=
        relErr_le hA hx hb
    _ ≤ cond2 A * etaFormula A b xh * (1 + (1 + t)) := by gcongr
    _ = (2 + t) * cond2 A * etaFormula A b xh := by ring

/-! ## The Banach lemma -/

section Banach
variable {R : Type*} [NormedRing R] [CompleteSpace R] [NormOneClass R]

/-- The partial sums of the Neumann series: `(I + X) ∑_{j ≤ m} (-X)^j = I - (-X)^{m+1}`. -/
lemma neumann_partial (X : R) (m : ℕ) :
    (1 + X) * ∑ j ∈ Finset.range (m + 1), (-X) ^ j = 1 - (-X) ^ (m + 1) := by
  have := mul_neg_geom_sum (-X) (m + 1)
  rwa [sub_neg_eq_add] at this

/-- The Neumann series `S = ∑ (-X)^j` is a two-sided inverse of `I + X`. -/
lemma banach_inverse {X : R} (hX : ‖X‖ < 1) :
    (1 + X) * ∑' j, (-X) ^ j = 1 ∧ (∑' j, (-X) ^ j) * (1 + X) = 1 := by
  have h : ‖-X‖ < 1 := by rwa [norm_neg]
  refine ⟨?_, ?_⟩
  · have := mul_neg_geom_series (-X) h
    rwa [sub_neg_eq_add] at this
  · have := geom_series_mul_neg (-X) h
    rwa [sub_neg_eq_add] at this

/-- **Banach lemma.** -/
theorem banach {X : R} (hX : ‖X‖ < 1) :
    IsUnit (1 + X) ∧ ‖Ring.inverse (1 + X)‖ ≤ 1 / (1 - ‖X‖) := by
  obtain ⟨h1, h2⟩ := banach_inverse hX
  have hu : IsUnit (1 + X) := ⟨⟨1 + X, ∑' j, (-X) ^ j, h1, h2⟩, rfl⟩
  refine ⟨hu, ?_⟩
  have hinv : Ring.inverse (1 + X) = ∑' j, (-X) ^ j := by
    rw [Ring.inverse_of_isUnit hu]
    exact (Units.inv_eq_of_mul_eq_one_right (u := hu.unit) h1)
  rw [hinv, one_div]
  refine tsum_of_norm_bounded (hasSum_geometric_of_lt_one (norm_nonneg X) hX) fun j => ?_
  calc ‖(-X) ^ j‖ ≤ ‖-X‖ ^ j := norm_pow_le _ j
    _ = ‖X‖ ^ j := by rw [norm_neg]

end Banach

instance instNormOneClass [NeZero n] : NormOneClass (Matrix (Fin n) (Fin n) ℝ) :=
  ⟨norm2_one⟩

/-- The Banach lemma for `‖·‖₂`. -/
theorem banach_matrix [NeZero n] {X : Matrix (Fin n) (Fin n) ℝ} (hX : norm2 X < 1) :
    IsUnit (1 + X).det ∧ norm2 (1 + X)⁻¹ ≤ 1 / (1 - norm2 X) := by
  obtain ⟨hu, hb⟩ := banach (R := Matrix (Fin n) (Fin n) ℝ) hX
  refine ⟨(Matrix.isUnit_iff_isUnit_det _).mp hu, ?_⟩
  rw [norm2, Matrix.nonsing_inv_eq_ringInverse]
  exact hb

/-! ## The perturbation theorem -/

lemma neZero_of_ne_zero {b : Fin n → ℝ} (hb : b ≠ 0) : NeZero n :=
  ⟨by rintro rfl; exact hb (Subsingleton.elim _ _)⟩

lemma norm2_pos_of_isUnit [NeZero n] {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) :
    0 < norm2 A := by
  apply norm2_pos_of_ne_zero
  rintro rfl
  rw [Matrix.det_zero (Fin.pos_iff_nonempty.mp (NeZero.pos n))] at hA
  exact not_isUnit_zero hA

/-- `‖(A + E)⁻¹‖ ≤ ‖A⁻¹‖ / (1 - ‖A⁻¹‖ ‖E‖)`. -/
theorem inv_add_bound [NeZero n] {A E : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    (hE : norm2 A⁻¹ * norm2 E < 1) :
    IsUnit (A + E).det ∧ norm2 (A + E)⁻¹ ≤ norm2 A⁻¹ / (1 - norm2 A⁻¹ * norm2 E) := by
  have hX : norm2 (A⁻¹ * E) < 1 := lt_of_le_of_lt (norm2_mul_le _ _) hE
  obtain ⟨hu, hbd⟩ := banach_matrix hX
  have hfac : A + E = A * (1 + A⁻¹ * E) := by
    rw [Matrix.mul_add, Matrix.mul_one, ← Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hA,
      Matrix.one_mul]
  have hdet : IsUnit (A + E).det := by rw [hfac, Matrix.det_mul]; exact hA.mul hu
  refine ⟨hdet, ?_⟩
  have hinv : (A + E)⁻¹ = (1 + A⁻¹ * E)⁻¹ * A⁻¹ := by rw [hfac, Matrix.mul_inv_rev]
  have hden : 0 < 1 - norm2 A⁻¹ * norm2 E := by linarith
  rw [hinv]
  calc norm2 ((1 + A⁻¹ * E)⁻¹ * A⁻¹) ≤ norm2 (1 + A⁻¹ * E)⁻¹ * norm2 A⁻¹ := norm2_mul_le _ _
    _ ≤ 1 / (1 - norm2 (A⁻¹ * E)) * norm2 A⁻¹ :=
        mul_le_mul_of_nonneg_right hbd (norm2_nonneg _)
    _ ≤ 1 / (1 - norm2 A⁻¹ * norm2 E) * norm2 A⁻¹ :=
        mul_le_mul_of_nonneg_right (one_div_le_one_div_of_le hden
          (by linarith [norm2_mul_le A⁻¹ E])) (norm2_nonneg _)
    _ = norm2 A⁻¹ / (1 - norm2 A⁻¹ * norm2 E) := by ring

/-- **Perturbation bound for linear systems** (`thm:perturbation_bound`). -/
theorem perturbation_bound {A E : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    {x b xh e : Fin n → ℝ} (hx : A *ᵥ x = b) (hb : b ≠ 0) (hE : norm2 A⁻¹ * norm2 E < 1)
    (hxh : (A + E) *ᵥ xh = b + e) :
    IsUnit (A + E).det ∧ vnorm (xh - x) / vnorm x ≤
      cond2 A / (1 - cond2 A * (norm2 E / norm2 A)) * (norm2 E / norm2 A + vnorm e / vnorm b) := by
  haveI := neZero_of_ne_zero hb
  obtain ⟨hdet, hinv⟩ := inv_add_bound hA hE
  refine ⟨hdet, ?_⟩
  have hApos := norm2_pos_of_isUnit hA
  have hxpos := vnorm_pos (ne_zero_of_mulVec hx hb)
  have hbpos := vnorm_pos hb
  have hden : 0 < 1 - norm2 A⁻¹ * norm2 E := by linarith
  -- `(A + E) δx = e - E x`
  have hδ : (A + E) *ᵥ (xh - x) = e - E *ᵥ x := by
    rw [Matrix.mulVec_sub, hxh, Matrix.add_mulVec, hx]; abel
  have hδx : xh - x = (A + E)⁻¹ *ᵥ (e - E *ᵥ x) := by
    rw [← hδ, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
  have h1 : vnorm (xh - x) ≤ norm2 A⁻¹ / (1 - norm2 A⁻¹ * norm2 E) *
      (vnorm e + norm2 E * vnorm x) := by
    rw [hδx]
    calc vnorm ((A + E)⁻¹ *ᵥ (e - E *ᵥ x)) ≤ norm2 (A + E)⁻¹ * vnorm (e - E *ᵥ x) :=
          vnorm_mulVec_le _ _
      _ ≤ norm2 A⁻¹ / (1 - norm2 A⁻¹ * norm2 E) * (vnorm e + norm2 E * vnorm x) := by
          apply mul_le_mul hinv _ (vnorm_nonneg _) (div_nonneg (norm2_nonneg _) hden.le)
          exact le_trans (vnorm_sub_le _ _) (add_le_add le_rfl (vnorm_mulVec_le _ _))
  have hbx := vnorm_le_of_mulVec hx
  have hκ : cond2 A * (norm2 E / norm2 A) = norm2 A⁻¹ * norm2 E := by
    unfold cond2; field_simp
  rw [hκ, div_le_iff₀ hxpos]
  have he : vnorm e ≤ norm2 A * vnorm x * (vnorm e / vnorm b) := by
    have h := mul_le_mul_of_nonneg_left hbx (div_nonneg (vnorm_nonneg e) hbpos.le)
    rw [div_mul_cancel₀ _ (ne_of_gt hbpos)] at h
    calc vnorm e ≤ vnorm e / vnorm b * (norm2 A * vnorm x) := h
      _ = norm2 A * vnorm x * (vnorm e / vnorm b) := by ring
  calc vnorm (xh - x) ≤ norm2 A⁻¹ / (1 - norm2 A⁻¹ * norm2 E) *
        (vnorm e + norm2 E * vnorm x) := h1
    _ ≤ norm2 A⁻¹ / (1 - norm2 A⁻¹ * norm2 E) *
        (norm2 A * vnorm x * (vnorm e / vnorm b) + norm2 E * vnorm x) :=
        mul_le_mul_of_nonneg_left (add_le_add_right he _) (div_nonneg (norm2_nonneg _) hden.le)
    _ = cond2 A / (1 - norm2 A⁻¹ * norm2 E) * (norm2 E / norm2 A + vnorm e / vnorm b) *
        vnorm x := by
        unfold cond2; field_simp; ring

/-- **From backward error to forward error.** If `κ(A) η(x̂) < 1` then
`‖x̂ - x‖/‖x‖ ≤ 2 κ(A) η(x̂) / (1 - κ(A) η(x̂))`. The proof applies the perturbation
theorem to a pair attaining the minimum in the definition of `η`. -/
theorem forward_error_bound {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    {x b xh : Fin n → ℝ} (hx : A *ᵥ x = b) (hb : b ≠ 0) (hκη : cond2 A * eta A b xh < 1) :
    vnorm (xh - x) / vnorm x ≤ 2 * cond2 A * eta A b xh / (1 - cond2 A * eta A b xh) := by
  haveI := neZero_of_ne_zero hb
  have hApos := norm2_pos_of_isUnit hA
  have hbpos := vnorm_pos hb
  have hmin := eta_isLeast A b xh hApos hb
  rw [← eta_eq A b xh hApos hb] at hmin
  obtain ⟨E, e, hxh, hτ⟩ := hmin.1
  set η := eta A b xh
  have hE : norm2 E / norm2 A ≤ η := by rw [hτ]; exact le_max_left _ _
  have he : vnorm e / vnorm b ≤ η := by rw [hτ]; exact le_max_right _ _
  have hκ0 := cond2_nonneg A
  have hκE : cond2 A * (norm2 E / norm2 A) ≤ cond2 A * η := mul_le_mul_of_nonneg_left hE hκ0
  have hκ : cond2 A * (norm2 E / norm2 A) = norm2 A⁻¹ * norm2 E := by
    unfold cond2; field_simp
  have hE' : norm2 A⁻¹ * norm2 E < 1 := by rw [← hκ]; linarith
  obtain ⟨-, hbound⟩ := perturbation_bound hA hx hb hE' hxh
  have hden : 0 < 1 - cond2 A * η := by linarith
  have hden' : 0 < 1 - cond2 A * (norm2 E / norm2 A) := by linarith
  calc vnorm (xh - x) / vnorm x
      ≤ cond2 A / (1 - cond2 A * (norm2 E / norm2 A)) * (norm2 E / norm2 A + vnorm e / vnorm b) :=
        hbound
    _ ≤ cond2 A / (1 - cond2 A * η) * (η + η) := by
        apply mul_le_mul _ (add_le_add hE he)
          (add_nonneg (div_nonneg (norm2_nonneg _) (norm2_nonneg _))
            (div_nonneg (vnorm_nonneg _) (vnorm_nonneg _))) (div_nonneg hκ0 hden.le)
        exact div_le_div_of_nonneg_left hκ0 hden (by linarith)
    _ = 2 * cond2 A * η / (1 - cond2 A * η) := by ring

end LU
