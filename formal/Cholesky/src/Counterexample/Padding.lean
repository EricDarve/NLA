import src.Counterexample.Condition

/-!
# Padding with a diagonal block

`B_N = diag(A, 4 I)`.

* The algorithm processes the leading block exactly as it processes `A` alone
  (`run_castAdd`), whatever the rest of the matrix is. Hence, if Cholesky fails on `A`, it
  fails on `B_N` within the original block (`fails_of_leading`).
* If `‖D‖₂ ≤ ‖A‖₂` and `‖D⁻¹‖₂ ≤ ‖A⁻¹‖₂`, then `κ₂(diag(A, D)) = κ₂(A)` (`cond2_blkdiag`). For
  `A_{9k}` and `D = 4I` this holds because `4 ≤ 17/2 ≤ ‖A‖₂` and `1/4 ≤ 17/(8δ) ≤ ‖A⁻¹‖₂`.
-/
namespace Cholesky
open Matrix

variable {a m : ℕ}

/-! ## The leading block is processed independently -/

lemma newCol_castAdd (o : Ops) (M : Matrix (Fin (a + m)) (Fin (a + m)) ℝ) (k i : Fin a) :
    newCol o M (Fin.castAdd m k) (Fin.castAdd m i) =
      newCol o (M.submatrix (Fin.castAdd m) (Fin.castAdd m)) k i := by
  simp only [newCol, Matrix.submatrix_apply, cA_inj, cA_lt]

lemma stepAt_castAdd (o : Ops) (M : Matrix (Fin (a + m)) (Fin (a + m)) ℝ) (k : Fin a) :
    (stepAt o (Fin.castAdd m k) M).submatrix (Fin.castAdd m) (Fin.castAdd m) =
      stepAt o k (M.submatrix (Fin.castAdd m) (Fin.castAdd m)) := by
  ext i j
  simp only [stepAt, Matrix.submatrix_apply, Matrix.of_apply, cA_inj, cA_lt, cA_le,
    newCol_castAdd]

lemma run_castAdd (o : Ops) (M : Matrix (Fin (a + m)) (Fin (a + m)) ℝ) :
    ∀ s ≤ a, (run o s M).submatrix (Fin.castAdd m) (Fin.castAdd m) =
      run o s (M.submatrix (Fin.castAdd m) (Fin.castAdd m))
  | 0, _ => rfl
  | s + 1, hs => by
    rw [run_succ, run_succ, stepN, stepN, dif_pos (by omega : s < a + m), dif_pos (by omega : s < a)]
    have : (⟨s, by omega⟩ : Fin (a + m)) = Fin.castAdd m ⟨s, by omega⟩ := rfl
    rw [this, stepAt_castAdd, run_castAdd o M s (by omega)]

lemma pivot_castAdd (o : Ops) (M : Matrix (Fin (a + m)) (Fin (a + m)) ℝ) (k : Fin a) :
    pivot o M (Fin.castAdd m k) = pivot o (M.submatrix (Fin.castAdd m) (Fin.castAdd m)) k := by
  unfold pivot
  have := congrFun (congrFun (run_castAdd o M k k.2.le) k) k
  simp only [Matrix.submatrix_apply] at this
  rw [← this]
  rfl

/-- **Cholesky fails within the original block**: if it fails on the leading block, it
fails on the whole matrix. -/
theorem fails_of_leading (o : Ops) (M : Matrix (Fin (a + m)) (Fin (a + m)) ℝ)
    (h : ¬ Succeeds o (M.submatrix (Fin.castAdd m) (Fin.castAdd m))) : ¬ Succeeds o M := by
  intro hM
  apply h
  intro k
  rw [← pivot_castAdd]
  exact hM _

/-! ## Block-diagonal matrices -/

lemma blkdiag_mulVec (A : Matrix (Fin a) (Fin a) ℝ) (D : Matrix (Fin m) (Fin m) ℝ)
    (x : Fin a → ℝ) (y : Fin m → ℝ) :
    blk A 0 0 D *ᵥ Fin.append x y = Fin.append (A *ᵥ x) (D *ᵥ y) := by
  rw [blk_mulVec_append]; simp

lemma norm2_blkdiag_le (A : Matrix (Fin a) (Fin a) ℝ) (D : Matrix (Fin m) (Fin m) ℝ)
    (hD : norm2 D ≤ norm2 A) : norm2 (blk A 0 0 D) ≤ norm2 A := by
  apply norm2_le_of (norm2_nonneg A)
  intro v
  rw [← append_vL_vR v, blkdiag_mulVec]
  have h1 := vnorm_mulVec_le A (vL v)
  have h2 := vnorm_mulVec_le D (vR v)
  have hsq : vnorm (Fin.append (A *ᵥ vL v) (D *ᵥ vR v)) ^ 2 ≤
      (norm2 A * vnorm (Fin.append (vL v) (vR v))) ^ 2 := by
    rw [vnorm_append_sq, mul_pow, vnorm_append_sq]
    have e1 := mul_self_le_mul_self (vnorm_nonneg _) h1
    have e2 := mul_self_le_mul_self (vnorm_nonneg _) h2
    have e3 := mul_self_le_mul_self (mul_nonneg (norm2_nonneg D) (vnorm_nonneg _))
      (mul_le_mul_of_nonneg_right hD (vnorm_nonneg (vR v)))
    nlinarith
  exact le_of_pow_le_pow_left₀ two_ne_zero (mul_nonneg (norm2_nonneg A) (vnorm_nonneg _)) hsq

lemma le_norm2_blkdiag (A : Matrix (Fin a) (Fin a) ℝ) (D : Matrix (Fin m) (Fin m) ℝ) :
    norm2 A ≤ norm2 (blk A 0 0 D) := by
  apply norm2_le_of (norm2_nonneg _)
  intro x
  have := vnorm_mulVec_le (blk A 0 0 D) (Fin.append x 0)
  rw [blkdiag_mulVec, Matrix.mulVec_zero] at this
  have e1 : vnorm (Fin.append (A *ᵥ x) (0 : Fin m → ℝ)) = vnorm (A *ᵥ x) := by
    have := vnorm_append_sq (A *ᵥ x) (0 : Fin m → ℝ)
    rw [vnorm_zero] at this
    nlinarith [vnorm_nonneg (Fin.append (A *ᵥ x) (0 : Fin m → ℝ)), vnorm_nonneg (A *ᵥ x)]
  have e2 : vnorm (Fin.append x (0 : Fin m → ℝ)) = vnorm x := by
    have := vnorm_append_sq x (0 : Fin m → ℝ)
    rw [vnorm_zero] at this
    nlinarith [vnorm_nonneg (Fin.append x (0 : Fin m → ℝ)), vnorm_nonneg x]
  rw [e1, e2] at this
  exact this

lemma norm2_blkdiag (A : Matrix (Fin a) (Fin a) ℝ) (D : Matrix (Fin m) (Fin m) ℝ)
    (hD : norm2 D ≤ norm2 A) : norm2 (blk A 0 0 D) = norm2 A :=
  le_antisymm (norm2_blkdiag_le A D hD) (le_norm2_blkdiag A D)

lemma blkdiag_inv (A : Matrix (Fin a) (Fin a) ℝ) (D : Matrix (Fin m) (Fin m) ℝ)
    (hA : IsUnit A.det) (hD : IsUnit D.det) :
    (blk A 0 0 D)⁻¹ = blk A⁻¹ 0 0 D⁻¹ := by
  apply Matrix.inv_eq_right_inv
  rw [blk_mul, Matrix.mul_nonsing_inv _ hA, Matrix.mul_nonsing_inv _ hD, one_eq_blk]
  simp

/-- `κ₂(diag(A, D)) = κ₂(A)` when the eigenvalues of `D` lie between those of `A`. -/
theorem cond2_blkdiag (A : Matrix (Fin a) (Fin a) ℝ) (D : Matrix (Fin m) (Fin m) ℝ)
    (hA : IsUnit A.det) (hDu : IsUnit D.det)
    (hD : norm2 D ≤ norm2 A) (hDinv : norm2 D⁻¹ ≤ norm2 A⁻¹) :
    cond2 (blk A 0 0 D) = cond2 A := by
  unfold cond2
  rw [blkdiag_inv A D hA hDu, norm2_blkdiag A D hD, norm2_blkdiag _ _ hDinv]

lemma norm2_one_le : norm2 (1 : Matrix (Fin m) (Fin m) ℝ) ≤ 1 :=
  norm2_le_of zero_le_one fun v => by rw [Matrix.one_mulVec, one_mul]

lemma norm2_smul_one_le (c : ℝ) : norm2 (c • (1 : Matrix (Fin m) (Fin m) ℝ)) ≤ |c| := by
  rw [norm2_smul]
  nlinarith [norm2_one_le (m := m), abs_nonneg c, norm2_nonneg (1 : Matrix (Fin m) (Fin m) ℝ)]

lemma smul_one_inv (c : ℝ) (hc : c ≠ 0) :
    (c • (1 : Matrix (Fin m) (Fin m) ℝ))⁻¹ = c⁻¹ • 1 := by
  apply Matrix.inv_eq_right_inv
  rw [smul_mul_smul_comm, Matrix.one_mul, mul_inv_cancel₀ hc, one_smul]

lemma det_smul_one_unit (c : ℝ) (hc : c ≠ 0) : IsUnit (c • (1 : Matrix (Fin m) (Fin m) ℝ)).det := by
  rw [Matrix.det_smul, Matrix.det_one, mul_one]
  exact isUnit_iff_ne_zero.mpr (pow_ne_zero _ hc)

/-! ## The padded counterexample -/

namespace CE
namespace Par
variable (P : Par) (m : ℕ)

/-- `B = diag(A_{9k}, 4 I_m)`. -/
noncomputable def B : Matrix (Fin ((6 * P.k + P.k) + (P.k + P.k) + m))
    (Fin ((6 * P.k + P.k) + (P.k + P.k) + m)) ℝ :=
  blk P.A 0 0 ((4 : ℝ) • 1)

lemma B_leading : (P.B m).submatrix (Fin.castAdd m) (Fin.castAdd m) = P.A := by
  ext i j; simp [B]

/-- **Cholesky fails on every padded matrix** `B_N = diag(A_{9k}, 4I)`, in the original block. -/
theorem B_fails {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) :
    ¬ Succeeds (Ops.rounded r) (P.B m) ∧ ¬ Succeeds (Ops.fused r) (P.B m) := by
  obtain ⟨h1, h2⟩ := P.A_fails hr
  exact ⟨fails_of_leading _ _ (by rw [B_leading]; exact h1),
    fails_of_leading _ _ (by rw [B_leading]; exact h2)⟩

/-- **Padding preserves the condition number.** -/
theorem cond2_B : cond2 (P.B m) = cond2 P.A := by
  have hδ := P.δ_pos
  have hdet : IsUnit P.A.det := (Matrix.isUnit_iff_isUnit_det _).mp P.A_posDef.isUnit
  apply cond2_blkdiag _ _ hdet (det_smul_one_unit 4 (by norm_num))
  · calc norm2 ((4 : ℝ) • (1 : Matrix (Fin m) (Fin m) ℝ)) ≤ |4| := norm2_smul_one_le 4
      _ = 4 := by norm_num
      _ ≤ 17 / 2 := by norm_num
      _ ≤ norm2 P.A := P.norm2_A_ge
  · rw [smul_one_inv 4 (by norm_num)]
    have hδle : P.δ ≤ 12 / 256 := by rw [P.δ_eq_x]; linarith [P.x_le]
    calc norm2 ((4 : ℝ)⁻¹ • (1 : Matrix (Fin m) (Fin m) ℝ)) ≤ |(4 : ℝ)⁻¹| := norm2_smul_one_le _
      _ = 1 / 4 := by norm_num
      _ ≤ 17 / (8 * P.δ) := by rw [div_le_div_iff₀ (by norm_num) (by positivity)]; linarith
      _ ≤ norm2 P.A⁻¹ := P.norm2_Ainv_ge

/-- `B_N` is SPD. -/
theorem B_posDef : (P.B m).PosDef := by
  rw [posDef_iff_real]
  have hA := (posDef_iff_real P.A).mp P.A_posDef
  refine ⟨by simp only [B, blk_transpose, hA.1, Matrix.transpose_zero, Matrix.transpose_smul,
    Matrix.transpose_one], fun v hv => ?_⟩
  rw [← append_vL_vR v, B, blkdiag_mulVec, append_dotProduct]
  have h2 : vR v ⬝ᵥ (((4 : ℝ) • (1 : Matrix (Fin m) (Fin m) ℝ)) *ᵥ vR v) = 4 * vnorm (vR v) ^ 2 := by
    rw [Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_smul, smul_eq_mul, vnorm_sq]
  rw [h2]
  by_cases hx : vL v = 0
  · have hy : vR v ≠ 0 := by
      intro hy; apply hv; rw [← append_vL_vR v, hx, hy, append_zero]
    rw [hx, Matrix.mulVec_zero, dotProduct_zero, zero_add]
    have := vnorm_pos hy
    positivity
  · have := hA.2 _ hx
    have := sq_nonneg (vnorm (vR v))
    linarith

end Par
end CE
end Cholesky
