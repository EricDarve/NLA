import src.Factorization

/-!
# Symmetric positive definite matrices

The notes define SPD matrices by `Aᵀ = A` and `xᵀ A x > 0` for `x ≠ 0`; mathlib's
`Matrix.PosDef` is the same notion (`posDef_iff_real`). This file proves the facts of
the section "Symmetric Positive Definite Matrices":

* `quad_eq_sum_eigenvalues`: with `A = Q Λ Qᵀ` and `z = Qᵀ x`, `xᵀ A x = ∑ λᵢ zᵢ²`;
* `posDef_iff_eigenvalues_pos`: a symmetric matrix is SPD iff all eigenvalues are
  positive (proved as in the notes);
* `posDef_mulVec_eq_zero`, `posDef_isUnit`: an SPD matrix is nonsingular;
* `posDef_diag_pos`: `aᵢᵢ = eᵢᵀ A eᵢ > 0`;
* the example `!![1, 2; 2, 1]`;
* a positive semidefinite matrix with a zero pivot.
-/
namespace Cholesky
open Matrix

variable {n : ℕ}

section Spectral
variable {A : Matrix (Fin n) (Fin n) ℝ} (hA : A.IsHermitian)

/-- The orthogonal matrix `Q` of eigenvectors. -/
noncomputable def eigQ : Matrix (Fin n) (Fin n) ℝ := (hA.eigenvectorUnitary : Matrix _ _ ℝ)

lemma eigQ_star : star (eigQ hA) = (eigQ hA)ᵀ := by
  rw [Matrix.star_eq_conjTranspose, Matrix.conjTranspose_eq_transpose_of_trivial]

lemma eigQ_mul_transpose : eigQ hA * (eigQ hA)ᵀ = 1 := by
  rw [← eigQ_star]
  exact Matrix.mem_unitaryGroup_iff.mp hA.eigenvectorUnitary.2

lemma eigQ_transpose_mul : (eigQ hA)ᵀ * eigQ hA = 1 := by
  rw [← eigQ_star]
  exact Matrix.mem_unitaryGroup_iff'.mp hA.eigenvectorUnitary.2

/-- `A = Q Λ Qᵀ`. -/
lemma spectral_real : A = eigQ hA * diagonal hA.eigenvalues * (eigQ hA)ᵀ := by
  conv_lhs => rw [hA.spectral_theorem]
  rw [← eigQ_star]
  rfl

/-- `xᵀ A x = zᵀ Λ z = ∑ λᵢ zᵢ²` for `z = Qᵀ x`. -/
theorem quad_eq_sum_eigenvalues (x : Fin n → ℝ) :
    x ⬝ᵥ (A *ᵥ x) = ∑ i, hA.eigenvalues i * ((eigQ hA)ᵀ *ᵥ x) i ^ 2 := by
  conv_lhs => rw [spectral_real hA]
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
    ← Matrix.mulVec_transpose]
  simp only [dotProduct, Matrix.mulVec_diagonal]
  apply Finset.sum_congr rfl
  intro i _
  ring

/-- The change of variables `z = Qᵀ x` preserves the Euclidean length. -/
lemma eigQ_dot (x : Fin n → ℝ) :
    ((eigQ hA)ᵀ *ᵥ x) ⬝ᵥ ((eigQ hA)ᵀ *ᵥ x) = x ⬝ᵥ x := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.transpose_transpose,
    Matrix.mulVec_mulVec, eigQ_mul_transpose, Matrix.one_mulVec]

lemma eigQ_mulVec_eq_zero {x : Fin n → ℝ} (h : (eigQ hA)ᵀ *ᵥ x = 0) : x = 0 := by
  have : x = eigQ hA *ᵥ ((eigQ hA)ᵀ *ᵥ x) := by
    rw [Matrix.mulVec_mulVec, eigQ_mul_transpose, Matrix.one_mulVec]
  rw [this, h, Matrix.mulVec_zero]

lemma eigenvector_ne_zero (i : Fin n) : (⇑(hA.eigenvectorBasis i) : Fin n → ℝ) ≠ 0 := by
  intro h
  apply hA.eigenvectorBasis.orthonormal.ne_zero i
  ext j
  exact congrFun h j

lemma quad_eigenvector (i : Fin n) :
    (⇑(hA.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ (A *ᵥ ⇑(hA.eigenvectorBasis i)) =
      hA.eigenvalues i * ((⇑(hA.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ
        ⇑(hA.eigenvectorBasis i)) := by
  rw [hA.mulVec_eigenvectorBasis i, dotProduct_smul, smul_eq_mul]

/-- A symmetric matrix is positive definite iff all its eigenvalues are positive,
proved as in the notes. -/
theorem posDef_iff_eigenvalues_pos : A.PosDef ↔ ∀ i, 0 < hA.eigenvalues i := by
  constructor
  · intro h i
    have hv := eigenvector_ne_zero hA i
    have hq := ((posDef_iff_real A).mp h).2 _ hv
    rw [quad_eigenvector hA i] at hq
    have hvv : 0 < (⇑(hA.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ ⇑(hA.eigenvectorBasis i) := by
      simpa using Matrix.dotProduct_star_self_pos_iff.mpr hv
    exact pos_of_mul_pos_left hq hvv.le
  · intro h
    rw [posDef_iff_real]
    refine ⟨by
      have := hA.eq
      rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at this, fun x hx => ?_⟩
    rw [quad_eq_sum_eigenvalues hA]
    have hz : (eigQ hA)ᵀ *ᵥ x ≠ 0 := fun hz => hx (eigQ_mulVec_eq_zero hA hz)
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hz
    apply Finset.sum_pos' (fun j _ => mul_nonneg (h j).le (sq_nonneg _))
    exact ⟨i, Finset.mem_univ _, mul_pos (h i) (lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hi)))⟩

end Spectral

/-! ## Nonsingularity and positive diagonal -/

lemma posDef_mulVec_eq_zero {A : Matrix (Fin n) (Fin n) ℝ} (hA : A.PosDef) {x : Fin n → ℝ}
    (hx : A *ᵥ x = 0) : x = 0 := by
  by_contra h
  have := ((posDef_iff_real A).mp hA).2 x h
  rw [hx, dotProduct_zero] at this
  exact lt_irrefl _ this

/-- An SPD matrix is nonsingular. -/
theorem posDef_isUnit {A : Matrix (Fin n) (Fin n) ℝ} (hA : A.PosDef) : IsUnit A := by
  rw [← Matrix.mulVec_injective_iff_isUnit]
  intro x y hxy
  have : A *ᵥ (x - y) = 0 := by rw [Matrix.mulVec_sub, hxy, sub_self]
  exact sub_eq_zero.mp (posDef_mulVec_eq_zero hA this)

/-- `aᵢᵢ = eᵢᵀ A eᵢ > 0`. -/
theorem posDef_diag_pos {A : Matrix (Fin n) (Fin n) ℝ} (hA : A.PosDef) (i : Fin n) :
    0 < A i i := by
  have := ((posDef_iff_real A).mp hA).2 (Pi.single i 1) (by simp)
  simpa [dotProduct, Matrix.mulVec, Pi.single_apply] using this

/-! ## Positive diagonal entries do not imply positive definiteness -/

/-- The matrix `!![1, 2; 2, 1]` of the notes. -/
def example12 : Matrix (Fin 2) (Fin 2) ℝ := !![1, 2; 2, 1]

lemma example12_diag_pos : ∀ i, 0 < example12 i i := by
  intro i; fin_cases i <;> simp [example12]

/-- Its eigenvalues are exactly `3` and `-1`. -/
theorem example12_eigenvalues (μ : ℝ) :
    (∃ v : Fin 2 → ℝ, v ≠ 0 ∧ example12 *ᵥ v = μ • v) ↔ μ = 3 ∨ μ = -1 := by
  constructor
  · rintro ⟨v, hv, h⟩
    have h0 := congrFun h 0
    have h1 := congrFun h 1
    simp [example12, Matrix.mulVec, dotProduct, Fin.sum_univ_two] at h0 h1
    have hprod : (μ - 3) * (μ + 1) * v 0 = 0 := by linear_combination (1 - μ) * h0 - 2 * h1
    have hprod' : (μ - 3) * (μ + 1) * v 1 = 0 := by linear_combination (1 - μ) * h1 - 2 * h0
    have : (μ - 3) * (μ + 1) = 0 := by
      by_contra hne
      apply hv
      funext i
      fin_cases i
      · exact (mul_eq_zero.mp hprod).resolve_left hne
      · exact (mul_eq_zero.mp hprod').resolve_left hne
    rcases mul_eq_zero.mp this with h | h
    · left; linarith
    · right; linarith
  · rintro (rfl | rfl)
    · refine ⟨![1, 1], by simp [funext_iff], ?_⟩
      ext i; fin_cases i <;> simp [example12, Matrix.mulVec, dotProduct] <;> norm_num
    · refine ⟨![1, -1], by simp [funext_iff], ?_⟩
      ext i; fin_cases i <;> simp [example12, Matrix.mulVec, dotProduct] <;> norm_num

theorem example12_not_posDef : ¬ example12.PosDef := by
  intro h
  have := ((posDef_iff_real _).mp h).2 ![1, -1] (by simp [funext_iff])
  simp [example12, Matrix.mulVec, dotProduct] at this

/-! ## Positive semidefinite matrices can have zero pivots -/

/-- A positive semidefinite matrix on which exact Cholesky stops at its first pivot. -/
def psdExample : Matrix (Fin 2) (Fin 2) ℝ := !![0, 0; 0, 1]

theorem psdExample_posSemidef : psdExample.PosSemidef := by
  rw [posSemidef_iff_real]
  refine ⟨by ext i j; fin_cases i <;> fin_cases j <;> rfl, fun x => ?_⟩
  simp [psdExample, Matrix.mulVec, dotProduct]
  nlinarith [sq_nonneg (x 1)]

theorem psdExample_zero_pivot : pivot Ops.exact psdExample 0 = 0 ∧
    ¬ Succeeds Ops.exact psdExample := by
  refine ⟨rfl, fun h => ?_⟩
  have := h 0
  rw [pivot_zero] at this
  simp [psdExample] at this

end Cholesky
