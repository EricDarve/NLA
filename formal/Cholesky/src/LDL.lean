import src.SPD

/-!
# Connection with LU and the `L D Lᵀ` factorization

For a Cholesky factor `L`, put `D = diag(lᵢᵢ²)` and `L₀ = L D^{-1/2}`. Then `L₀` is unit
lower triangular, `A = L₀ D L₀ᵀ = L₀ U` with `U = D L₀ᵀ` upper triangular, and the
diagonal of `U` (the LU pivots) is `lᵢᵢ² > 0`. Since the unpivoted LU factorization of
a nonsingular matrix is unique (`lu_unique`), this is *the* LU factorization of `A`.
The pivots examined by the exact Cholesky algorithm are the same numbers
(`exact_pivot_eq_sq`).
-/
namespace Cholesky
open Matrix

variable {n : ℕ}

/-- `U` is upper triangular. -/
def IsUpperTriangular (U : Matrix (Fin n) (Fin n) ℝ) : Prop := ∀ i j, j < i → U i j = 0

section LDL
variable (L : Matrix (Fin n) (Fin n) ℝ)

/-- `D = diag(l₁₁², …, lₙₙ²)`. -/
noncomputable def ldlD : Matrix (Fin n) (Fin n) ℝ := diagonal fun i => L i i ^ 2

/-- `L₀ = L D^{-1/2}`. -/
noncomputable def ldlL0 : Matrix (Fin n) (Fin n) ℝ := L * diagonal fun i => 1 / L i i

/-- `U = D L₀ᵀ`. -/
noncomputable def ldlU : Matrix (Fin n) (Fin n) ℝ := ldlD L * (ldlL0 L)ᵀ

lemma ldlL0_apply (i j : Fin n) : ldlL0 L i j = L i j / L j j := by
  simp [ldlL0, Matrix.mul_diagonal, div_eq_mul_inv]

variable {L}
variable (hL : IsLowerTriangular L) (hd : ∀ i, 0 < L i i)
include hL hd

theorem ldlL0_unit_lower : IsLowerTriangular (ldlL0 L) ∧ ∀ i, ldlL0 L i i = 1 := by
  refine ⟨fun i j hij => ?_, fun i => ?_⟩
  · rw [ldlL0_apply, hL i j hij, zero_div]
  · rw [ldlL0_apply, div_self (ne_of_gt (hd i))]

omit hL in
/-- `L Lᵀ = L₀ D L₀ᵀ`. -/
theorem ldl_eq : L * Lᵀ = ldlL0 L * ldlD L * (ldlL0 L)ᵀ := by
  ext i j
  simp only [Matrix.mul_apply, ldlL0_apply, ldlD, Matrix.transpose_apply]
  apply Finset.sum_congr rfl
  intro k _
  have hk : L k k ≠ 0 := ne_of_gt (hd k)
  rw [Finset.sum_eq_single k]
  · simp only [diagonal_apply_eq]
    field_simp
  · intro b _ hb
    simp [diagonal_apply_ne _ hb]
  · simp

omit hL in
/-- `A = L₀ U`. -/
theorem lu_eq : L * Lᵀ = ldlL0 L * ldlU L := by
  rw [ldl_eq hd, ldlU, Matrix.mul_assoc]

theorem ldlU_upper : IsUpperTriangular (ldlU L) := by
  intro i j hij
  simp only [ldlU, ldlD, Matrix.diagonal_mul, Matrix.transpose_apply,
    (ldlL0_unit_lower hL hd).1 j i hij, mul_zero]

/-- The LU pivots are `uᵢᵢ = lᵢᵢ² > 0`. -/
theorem ldlU_diag (i : Fin n) : ldlU L i i = L i i ^ 2 ∧ 0 < ldlU L i i := by
  have h : ldlU L i i = L i i ^ 2 := by
    simp only [ldlU, ldlD, Matrix.diagonal_mul, Matrix.transpose_apply,
      (ldlL0_unit_lower hL hd).2 i, mul_one]
  exact ⟨h, by rw [h]; exact pow_pos (hd i) 2⟩

end LDL

/-! ## Uniqueness of the unpivoted LU factorization -/

lemma lower_blockTriangular {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L) :
    L.BlockTriangular OrderDual.toDual := fun i j hij => hL i j hij

lemma upper_blockTriangular {U : Matrix (Fin n) (Fin n) ℝ} (hU : IsUpperTriangular U) :
    U.BlockTriangular id := fun i j hij => hU i j hij

lemma det_upper {U : Matrix (Fin n) (Fin n) ℝ} (hU : IsUpperTriangular U) :
    U.det = ∏ i, U i i :=
  Matrix.det_of_upperTriangular (upper_blockTriangular hU)

lemma det_lower {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L) :
    L.det = ∏ i, L i i :=
  Matrix.det_of_lowerTriangular L (lower_blockTriangular hL)

/-- A matrix that is both lower triangular with unit diagonal and upper triangular is
the identity. -/
lemma eq_one_of_lower_upper {M : Matrix (Fin n) (Fin n) ℝ} (hl : IsLowerTriangular M)
    (hu : IsUpperTriangular M) (hd : ∀ i, M i i = 1) : M = 1 := by
  ext i j
  rcases lt_trichotomy i j with h | h | h
  · rw [hl i j h, Matrix.one_apply_ne (ne_of_lt h)]
  · subst h; rw [hd, Matrix.one_apply_eq]
  · rw [hu i j h, Matrix.one_apply_ne (ne_of_gt h)]

lemma diag_mul_lower {M N : Matrix (Fin n) (Fin n) ℝ} (hM : IsLowerTriangular M)
    (hN : IsLowerTriangular N) (i : Fin n) : (M * N) i i = M i i * N i i := by
  rw [Matrix.mul_apply, Finset.sum_eq_single i]
  · intro k _ hk
    rcases lt_or_gt_of_ne hk with h | h
    · rw [hN k i h, mul_zero]
    · rw [hM i k h, zero_mul]
  · simp

/-- **Uniqueness of LU.** If a matrix with nonzero determinant is written as `L₁ U₁` and
`L₂ U₂`, with unit lower triangular `Lᵢ` and upper triangular `Uᵢ`, the factors agree. -/
theorem lu_unique {L₁ L₂ U₁ U₂ : Matrix (Fin n) (Fin n) ℝ}
    (hL₁ : IsLowerTriangular L₁) (hL₂ : IsLowerTriangular L₂)
    (hd₁ : ∀ i, L₁ i i = 1) (hd₂ : ∀ i, L₂ i i = 1)
    (hU₁ : IsUpperTriangular U₁) (hU₂ : IsUpperTriangular U₂)
    (hdet : (L₁ * U₁).det ≠ 0) (h : L₁ * U₁ = L₂ * U₂) : L₁ = L₂ ∧ U₁ = U₂ := by
  have hL₂det : L₂.det = 1 := by rw [det_lower hL₂]; simp [hd₂]
  have hL₁det : L₁.det = 1 := by rw [det_lower hL₁]; simp [hd₁]
  have hU₁det : U₁.det ≠ 0 := by rwa [Matrix.det_mul, hL₁det, one_mul] at hdet
  have hL₂u : IsUnit L₂.det := by rw [hL₂det]; exact isUnit_one
  have hU₁u : IsUnit U₁.det := isUnit_iff_ne_zero.mpr hU₁det
  letI : Invertible L₂ := Matrix.invertibleOfIsUnitDet L₂ hL₂u
  letI : Invertible U₁ := Matrix.invertibleOfIsUnitDet U₁ hU₁u
  -- `M = L₂⁻¹ L₁ = U₂ U₁⁻¹`
  have hM : L₂⁻¹ * L₁ = U₂ * U₁⁻¹ := by
    calc L₂⁻¹ * L₁ = L₂⁻¹ * L₁ * (U₁ * U₁⁻¹) := by
            rw [Matrix.mul_nonsing_inv _ hU₁u, Matrix.mul_one]
      _ = L₂⁻¹ * (L₁ * U₁) * U₁⁻¹ := by simp only [Matrix.mul_assoc]
      _ = L₂⁻¹ * (L₂ * U₂) * U₁⁻¹ := by rw [h]
      _ = U₂ * U₁⁻¹ := by
            rw [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hL₂u, Matrix.one_mul]
  have hinvL : IsLowerTriangular L₂⁻¹ := fun i j hij =>
    blockTriangular_inv_of_blockTriangular (lower_blockTriangular hL₂) hij
  have hinvU : IsUpperTriangular U₁⁻¹ := fun i j hij =>
    blockTriangular_inv_of_blockTriangular (upper_blockTriangular hU₁) hij
  have hlow : IsLowerTriangular (L₂⁻¹ * L₁) := fun i j hij =>
    (lower_blockTriangular hinvL).mul (lower_blockTriangular hL₁) hij
  have hup : IsUpperTriangular (L₂⁻¹ * L₁) := by
    rw [hM]
    exact fun i j hij => (upper_blockTriangular hU₂).mul (upper_blockTriangular hinvU) hij
  have hinvdiag : ∀ i, L₂⁻¹ i i = 1 := by
    intro i
    have := diag_mul_lower hinvL hL₂ i
    rw [Matrix.nonsing_inv_mul _ hL₂u, Matrix.one_apply_eq, hd₂, mul_one] at this
    exact this.symm
  have hone : L₂⁻¹ * L₁ = 1 :=
    eq_one_of_lower_upper hlow hup (fun i => by rw [diag_mul_lower hinvL hL₁, hinvdiag, hd₁,
      one_mul])
  have hL : L₁ = L₂ := by
    calc L₁ = L₂ * (L₂⁻¹ * L₁) := by
          rw [← Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hL₂u, Matrix.one_mul]
      _ = L₂ := by rw [hone, Matrix.mul_one]
  refine ⟨hL, ?_⟩
  have hL₁u : IsUnit L₁.det := by rw [hL₁det]; exact isUnit_one
  calc U₁ = L₁⁻¹ * (L₁ * U₁) := by
        rw [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hL₁u, Matrix.one_mul]
    _ = L₁⁻¹ * (L₁ * U₂) := by rw [h, hL]
    _ = U₂ := by rw [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hL₁u, Matrix.one_mul]

/-- For an SPD matrix, `(L₀, U)` is the unique unpivoted LU factorization. -/
theorem cholesky_is_lu (A : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosDef)
    (L : Matrix (Fin n) (Fin n) ℝ) (hL : IsLowerTriangular L) (hd : ∀ i, 0 < L i i)
    (hLA : A = L * Lᵀ) (L' U' : Matrix (Fin n) (Fin n) ℝ) (hL' : IsLowerTriangular L')
    (hd' : ∀ i, L' i i = 1) (hU' : IsUpperTriangular U') (hLU : A = L' * U') :
    L' = ldlL0 L ∧ U' = ldlU L := by
  have hdet : (L' * U').det ≠ 0 := by
    rw [← hLU]; exact ne_of_gt (Matrix.PosDef.det_pos hA)
  obtain ⟨h1, h2⟩ := lu_unique hL' (ldlL0_unit_lower hL hd).1 hd' (ldlL0_unit_lower hL hd).2
    hU' (ldlU_upper hL hd) hdet (by rw [← hLU, hLA, lu_eq hd])
  exact ⟨h1, h2⟩

/-! ## The exact Cholesky pivots are the LU pivots -/

/-- The diagonal of the computed factor is the computed square root of the pivots. -/
theorem factor_diag_eq_sqrt_pivot (o : Ops) :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (i : Fin n),
      factor o A i i = o.sqrt (pivot o A i)
  | 0, _, i => Fin.elim0 i
  | n + 1, A, i => by
    refine Fin.cases ?_ (fun i => ?_) i
    · rw [factor_zero_zero, pivot_zero]
    · rw [factor_succ_succ, pivot_succ, factor_diag_eq_sqrt_pivot]

theorem exact_pivot_eq_sq (A : Matrix (Fin n) (Fin n) ℝ) (h : Succeeds Ops.exact A)
    (i : Fin n) : pivot Ops.exact A i = factor Ops.exact A i i ^ 2 := by
  rw [factor_diag_eq_sqrt_pivot]
  exact (Real.sq_sqrt (h i).le).symm

end Cholesky
