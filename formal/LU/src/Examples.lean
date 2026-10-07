import src.Stability

/-!
# The examples of the notes (exact arithmetic)

* `swap_no_lu`: `[[0, 1], [1, 0]]` is invertible but has no factorization `A = L U`
  (LU without pivoting fails at the zero pivot).
* `Aeps`: `A_ε = [[ε, 1], [1, π]]`. Without pivoting its factors are
  `L = [[1, 0], [1/ε, 1]]`, `U = [[ε, 1], [0, π - 1/ε]]` (`Aeps_lu`), and `π` is
  reconstructed by the cancellation `π = 1/ε + (π - 1/ε)` (`Aeps_cancel`). The routine
  (for `0 < ε < 1`) interchanges the rows and computes `L = [[1, 0], [ε, 1]]`,
  `U = [[1, π], [0, 1 - επ]]` (`Aeps_final`, `Aeps_plu`).
* `Aeps_absProd_unpivoted`: the two terms reconstructing `a₂₂` have total magnitude
  `2/ε - π`; `Aeps_absProd_pivoted`: after the interchange `|L| |U| = |P A|`.
* `A3`: the 3 × 3 example. The routine swaps rows at both steps (`A3_piv0`, `A3_piv1`), the
  factors are as computed in the notes (`A3_final`, `A3_plu`), and the computed solution of
  `A x = (3, 7, 17)` is `(1, 0, 1)` (`A3_solve`).
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

/-! ## A zero pivot -/

/-- An invertible matrix with no LU factorization. -/
theorem swap_no_lu :
    (!![0, 1; 1, 0] : Matrix (Fin 2) (Fin 2) ℝ).det ≠ 0 ∧
      ¬ ∃ L U : Matrix (Fin 2) (Fin 2) ℝ, IsUnitLower L ∧ IsUpper U ∧ !![0, 1; 1, 0] = L * U := by
  refine ⟨by simp [Matrix.det_fin_two], ?_⟩
  rintro ⟨L, U, hL, hU, h⟩
  have h00 := congrFun (congrFun h 0) 0
  have h10 := congrFun (congrFun h 1) 0
  simp [Matrix.mul_apply, Fin.sum_univ_two, hL.1, hL.2 0 1 (by decide), hU 1 0 (by decide)]
    at h00 h10
  rw [← h00, mul_zero] at h10
  exact one_ne_zero h10

/-! ## `A_ε` -/

/-- `A_ε = [[ε, 1], [1, π]]`. -/
noncomputable def Aeps (ε : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![ε, 1; 1, Real.pi]

/-- The factors without pivoting. -/
theorem Aeps_lu {ε : ℝ} (hε : ε ≠ 0) :
    (!![1, 0; 1 / ε, 1] : Matrix (Fin 2) (Fin 2) ℝ) * !![ε, 1; 0, Real.pi - 1 / ε] = Aeps ε := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Aeps, Matrix.mul_apply, Fin.sum_univ_two]
  field_simp

/-- The cancellation that reconstructs `a₂₂`. -/
theorem Aeps_cancel (ε : ℝ) : Real.pi = 1 / ε + (Real.pi - 1 / ε) := by ring

lemma pi_lt_four : Real.pi < 4 := by linarith [Real.pi_lt_d2]

/-- The pivot row at the first step is the second row. -/
lemma Aeps_piv {ε : ℝ} (hε : 0 ≤ ε) (hε1 : ε < 1) : pivRow (Aeps ε) 0 = 1 := by
  apply firstMax_eq
  refine ⟨by decide, fun i _ => ?_, fun i _ hi => ?_⟩
  · fin_cases i <;> simp [Aeps, abs_of_nonneg hε, hε1.le]
  · fin_cases i
    · simp [Aeps, abs_of_nonneg hε, hε1]
    · exact absurd hi (by decide)

/-- The routine on `A_ε`: one interchange, `l₂₁ = ε`, `u₂₂ = 1 - επ`. -/
theorem Aeps_final {ε : ℝ} (hε : 0 ≤ ε) (hε1 : ε < 1) :
    final (Ops.exact ℝ) (Aeps ε) = ⟨!![1, Real.pi; ε, 1 - ε * Real.pi], Equiv.swap 0 1⟩ := by
  rw [final, run_pred]
  show stepAt (Ops.exact ℝ) 0 ⟨Aeps ε, 1⟩ = _
  rw [stepAt]
  simp only [Aeps_piv hε hε1, one_mul]
  congr 1
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [elim, swapRows, Aeps, Equiv.swap_apply_def]

/-- `P A_ε = L U` with the factors of the notes. -/
theorem Aeps_plu {ε : ℝ} (hε : 0 ≤ ε) (hε1 : ε < 1) :
    factorP (Ops.exact ℝ) (Aeps ε) * Aeps ε = !![1, Real.pi; ε, 1] ∧
      factorL (Ops.exact ℝ) (Aeps ε) = !![1, 0; ε, 1] ∧
      factorU (Ops.exact ℝ) (Aeps ε) = !![1, Real.pi; 0, 1 - ε * Real.pi] ∧
      (!![1, 0; ε, 1] : Matrix (Fin 2) (Fin 2) ℝ) * !![1, Real.pi; 0, 1 - ε * Real.pi] =
        !![1, Real.pi; ε, 1] := by
  have hf := Aeps_final hε hε1
  refine ⟨?_, ?_, ?_, ?_⟩
  · ext i j
    rw [factorP, factorPerm, hf, permMatrix_mul_apply]
    fin_cases i <;> fin_cases j <;> simp [Aeps]
  · ext i j
    rw [factorL, hf]
    fin_cases i <;> fin_cases j <;> simp [lowerOf]
  · ext i j
    rw [factorU, hf]
    fin_cases i <;> fin_cases j <;> simp [upperOf]
  · ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.mul_apply, Fin.sum_univ_two]

/-- Without pivoting, the terms that reconstruct `a₂₂ = π` have total magnitude `2/ε - π`. -/
theorem Aeps_absProd_unpivoted {ε : ℝ} (hε : 0 < ε) (hεπ : ε * Real.pi < 1) :
    (absMat !![1, 0; 1 / ε, 1] * absMat !![ε, 1; 0, Real.pi - 1 / ε]) 1 1 = 2 / ε - Real.pi := by
  have h : Real.pi - ε⁻¹ < 0 := by
    have : Real.pi < 1 / ε := (lt_div_iff₀ hε).mpr (by linarith)
    rw [one_div] at this; linarith
  simp [Matrix.mul_apply, Fin.sum_univ_two, abs_of_pos hε]
  rw [abs_of_neg h]
  ring

/-- After the interchange, `|L| |U| = |P A_ε|`: no cancellation. -/
theorem Aeps_absProd_pivoted {ε : ℝ} (hε : 0 ≤ ε) (hεπ : ε * Real.pi ≤ 1) :
    absMat (!![1, 0; ε, 1] : Matrix (Fin 2) (Fin 2) ℝ) * absMat !![1, Real.pi; 0, 1 - ε * Real.pi] =
      absMat !![1, Real.pi; ε, 1] := by
  have h : 0 ≤ 1 - ε * Real.pi := by linarith
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_two, abs_of_nonneg hε, abs_of_nonneg h,
      abs_of_pos Real.pi_pos]

/-! ## The 3 × 3 example -/

/-- `A = [[2, 1, 1], [4, 3, 3], [8, 7, 9]]`. -/
noncomputable def A3 : Matrix (Fin 3) (Fin 3) ℝ := !![2, 1, 1; 4, 3, 3; 8, 7, 9]

/-- `b = (3, 7, 17)`. -/
noncomputable def b3 : Fin 3 → ℝ := ![3, 7, 17]

/-- The array after the first step. -/
noncomputable def A3s1 : Matrix (Fin 3) (Fin 3) ℝ :=
  !![8, 7, 9; 1 / 2, -1 / 2, -3 / 2; 1 / 4, -3 / 4, -5 / 4]

/-- The packed factors after the second step. -/
noncomputable def A3s2 : Matrix (Fin 3) (Fin 3) ℝ :=
  !![8, 7, 9; 1 / 4, -3 / 4, -5 / 4; 1 / 2, 2 / 3, -2 / 3]

/-- Step 1 chooses row 3 (`|8|` is largest). -/
lemma A3_piv0 : pivRow A3 0 = 2 := by
  apply firstMax_eq
  refine ⟨by decide, fun i _ => ?_, fun i _ hi => ?_⟩
  · fin_cases i <;> simp [A3] <;> norm_num
  · fin_cases i <;> simp [A3] at hi ⊢ <;> norm_num

lemma A3_step0 : stepAt (Ops.exact ℝ) 0 ⟨A3, 1⟩ = ⟨A3s1, Equiv.swap 0 2⟩ := by
  rw [stepAt]
  simp only [A3_piv0, one_mul]
  congr 1
  ext i j
  fin_cases i <;> fin_cases j <;> simp [elim, swapRows, A3, A3s1, Equiv.swap_apply_def] <;>
    norm_num

/-- Step 2 chooses row 3 again (`|-3/4| > |-1/2|`). -/
lemma A3_piv1 : pivRow A3s1 1 = 2 := by
  apply firstMax_eq
  refine ⟨by decide, fun i hi => ?_, fun i hi hi2 => ?_⟩
  · fin_cases i
    · exact absurd hi (by decide)
    · simp [A3s1]; norm_num
    · simp [A3s1]
  · fin_cases i
    · exact absurd hi (by decide)
    · simp [A3s1]; norm_num
    · exact absurd hi2 (by decide)

lemma A3_step1 : stepAt (Ops.exact ℝ) 1 ⟨A3s1, Equiv.swap 0 2⟩ =
    ⟨A3s2, Equiv.swap 0 2 * Equiv.swap 1 2⟩ := by
  rw [stepAt]
  simp only [A3_piv1]
  congr 1
  ext i j
  fin_cases i <;> fin_cases j <;> simp [elim, swapRows, A3s1, A3s2, Equiv.swap_apply_def] <;>
    norm_num

/-- **The routine swaps rows at both steps**; the second interchange also moves the
multiplier stored at the first step. -/
theorem A3_final : final (Ops.exact ℝ) A3 = ⟨A3s2, Equiv.swap 0 2 * Equiv.swap 1 2⟩ := by
  rw [final, run_pred]
  show stepN (Ops.exact ℝ) 1 (stepN (Ops.exact ℝ) 0 ⟨A3, 1⟩) = _
  rw [stepN, dif_pos (by norm_num), stepN, dif_pos (by norm_num)]
  exact (congrArg (stepAt (Ops.exact ℝ) 1) A3_step0).trans A3_step1

/-- The factors: `P A = L U` with the values of the notes. -/
theorem A3_plu :
    factorP (Ops.exact ℝ) A3 * A3 = !![8, 7, 9; 2, 1, 1; 4, 3, 3] ∧
      factorL (Ops.exact ℝ) A3 = !![1, 0, 0; 1 / 4, 1, 0; 1 / 2, 2 / 3, 1] ∧
      factorU (Ops.exact ℝ) A3 = !![8, 7, 9; 0, -3 / 4, -5 / 4; 0, 0, -2 / 3] ∧
      (!![1, 0, 0; 1 / 4, 1, 0; 1 / 2, 2 / 3, 1] : Matrix (Fin 3) (Fin 3) ℝ) *
        !![8, 7, 9; 0, -3 / 4, -5 / 4; 0, 0, -2 / 3] = !![8, 7, 9; 2, 1, 1; 4, 3, 3] := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · ext i j
    rw [factorP, factorPerm, A3_final, permMatrix_mul_apply]
    fin_cases i <;> fin_cases j <;> simp [A3, Equiv.swap_apply_def]
  · ext i j
    rw [factorL, A3_final]
    fin_cases i <;> fin_cases j <;> simp [lowerOf, A3s2]
  · ext i j
    rw [factorU, A3_final]
    fin_cases i <;> fin_cases j <;> simp [upperOf, A3s2]
  · ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.mul_apply, Fin.sum_univ_three] <;> norm_num

/-- `L y = P b` with `y = (17, -5/4, -2/3)`, and `U x = y` with `x = (1, 0, 1)`. -/
theorem A3_triangular :
    factorL (Ops.exact ℝ) A3 *ᵥ ![17, -5 / 4, -2 / 3] = factorP (Ops.exact ℝ) A3 *ᵥ b3 ∧
      factorU (Ops.exact ℝ) A3 *ᵥ ![1, 0, 1] = ![17, -5 / 4, -2 / 3] := by
  obtain ⟨-, hL, hU, -⟩ := A3_plu
  refine ⟨?_, ?_⟩
  · rw [hL, factorP, factorPerm, A3_final, permMatrix_mulVec]
    ext i
    fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three, b3,
      Equiv.swap_apply_def] <;> norm_num
  · rw [hU]
    ext i
    fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three]
    norm_num

/-- **The computed solution is `x = (1, 0, 1)`.** -/
theorem A3_solve : luSolve (Ops.exact ℝ) A3 b3 = ![1, 0, 1] := by
  have hd : ∀ i, factorU (Ops.exact ℝ) A3 i i ≠ 0 := by
    rw [A3_plu.2.2.1]
    intro i; fin_cases i <;> simp
  have h1 := exact_luSolve A3 b3 hd
  have h2 : A3 *ᵥ ![1, 0, 1] = b3 := by
    ext i
    fin_cases i <;> simp [A3, b3, Matrix.mulVec, dotProduct, Fin.sum_univ_three] <;> norm_num
  have hdet : IsUnit A3 := by
    rw [Matrix.isUnit_iff_isUnit_det, isUnit_iff_ne_zero]
    simp [A3, Matrix.det_fin_three]
    norm_num
  exact (Matrix.mulVec_injective_iff_isUnit.mpr hdet) (h1.trans h2.symm)

end LU
