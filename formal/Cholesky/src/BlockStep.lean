import Mathlib

/-!
# One elimination step: blocks and the Schur complement

For a matrix indexed by `Fin (n + 1)`, the first row and column play the roles of
`a₁₁`, `c`, and `cᵀ` in the notes, and the trailing block `B` is indexed by
`Fin.succ`. This file proves the block formulas for `L Lᵀ`, the formulas that
determine the first column of a Cholesky factor, the completing-the-square
identity, and the positive definiteness of the Schur complement.
-/
namespace Cholesky
open Matrix

variable {n : ℕ}

/-- `L` is lower triangular. -/
def IsLowerTriangular {m : ℕ} (L : Matrix (Fin m) (Fin m) ℝ) : Prop :=
  ∀ i j, i < j → L i j = 0

/-- The Schur complement `S = B - c cᵀ / a₁₁` of the leading entry. -/
noncomputable def schur (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j => A i.succ j.succ - A i.succ 0 * A j.succ 0 / A 0 0

/-- The trailing block of a matrix. -/
def trailing (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  A.submatrix Fin.succ Fin.succ

lemma posDef_iff_real (A : Matrix (Fin n) (Fin n) ℝ) :
    A.PosDef ↔ Aᵀ = A ∧ ∀ x : Fin n → ℝ, x ≠ 0 → 0 < x ⬝ᵥ (A *ᵥ x) := by
  simp [Matrix.PosDef, Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial]

lemma posSemidef_iff_real (A : Matrix (Fin n) (Fin n) ℝ) :
    A.PosSemidef ↔ Aᵀ = A ∧ ∀ x : Fin n → ℝ, 0 ≤ x ⬝ᵥ (A *ᵥ x) := by
  simp [Matrix.PosSemidef, Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial]

lemma symm_apply_of_transpose {A : Matrix (Fin n) (Fin n) ℝ} (h : Aᵀ = A) (i j : Fin n) :
    A i j = A j i := by
  calc A i j = Aᵀ j i := rfl
    _ = A j i := by rw [h]

/-! ## Block formulas for `L Lᵀ` -/

section Blocks
variable (L : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)

lemma mul_transpose_succ_succ (i j : Fin n) :
    (L * Lᵀ) i.succ j.succ =
      L i.succ 0 * L j.succ 0 + (trailing L * (trailing L)ᵀ) i j := by
  simp [Matrix.mul_apply, Fin.sum_univ_succ, trailing]

lemma mul_transpose_succ_zero (hL : IsLowerTriangular L) (i : Fin (n)) :
    (L * Lᵀ) i.succ 0 = L i.succ 0 * L 0 0 := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_succ]
  have : ∀ k : Fin n, L 0 k.succ = 0 := fun k => hL 0 k.succ (Fin.succ_pos k)
  simp [this]

lemma mul_transpose_zero_zero (hL : IsLowerTriangular L) :
    (L * Lᵀ) 0 0 = L 0 0 ^ 2 := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_succ]
  have : ∀ k : Fin n, L 0 k.succ = 0 := fun k => hL 0 k.succ (Fin.succ_pos k)
  simp [this, sq]

lemma trailing_lower (hL : IsLowerTriangular L) : IsLowerTriangular (trailing L) := by
  intro i j hij
  exact hL _ _ (Fin.succ_lt_succ_iff.mpr hij)

/-- Matching the blocks of `A = L Lᵀ` determines the first column of `L` and
leaves the Schur complement for the trailing factor. -/
theorem block_equations (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (hL : IsLowerTriangular L) (hpos : 0 < L 0 0) (hA : A = L * Lᵀ) :
    L 0 0 = Real.sqrt (A 0 0) ∧
    (∀ i : Fin n, L i.succ 0 = A i.succ 0 / L 0 0) ∧
    trailing L * (trailing L)ᵀ = schur A := by
  have h00 : A 0 0 = L 0 0 ^ 2 := by rw [hA, mul_transpose_zero_zero L hL]
  have hc : ∀ i : Fin n, A i.succ 0 = L i.succ 0 * L 0 0 := by
    intro i; rw [hA, mul_transpose_succ_zero L hL]
  refine ⟨?_, ?_, ?_⟩
  · rw [h00, Real.sqrt_sq hpos.le]
  · intro i; rw [hc i]; field_simp
  · ext i j
    have := mul_transpose_succ_succ L i j
    have hB : A i.succ j.succ = (L * Lᵀ) i.succ j.succ := by rw [hA]
    simp only [schur, Matrix.of_apply, hc, h00, hB, this]
    have hne : L 0 0 ≠ 0 := ne_of_gt hpos
    field_simp
    ring

end Blocks

/-! ## Completing the square -/

lemma quad_cons (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) (hA : Aᵀ = A)
    (t : ℝ) (y : Fin n → ℝ) :
    (Fin.cons t y : Fin (n + 1) → ℝ) ⬝ᵥ (A *ᵥ Fin.cons t y) =
      A 0 0 * t ^ 2 + 2 * t * (∑ i, A i.succ 0 * y i) +
        y ⬝ᵥ (trailing A *ᵥ y) := by
  have hsymm : ∀ j : Fin n, A 0 j.succ = A j.succ 0 := fun j =>
    symm_apply_of_transpose hA 0 j.succ
  simp only [dotProduct, Matrix.mulVec, Fin.sum_univ_succ, Fin.cons_zero, Fin.cons_succ,
    trailing, Matrix.submatrix_apply, hsymm]
  have e : ∀ i : Fin n, y i * (A i.succ 0 * t + ∑ j, A i.succ j.succ * y j) =
      t * (A i.succ 0 * y i) + y i * ∑ j, A i.succ j.succ * y j := by intro i; ring
  simp only [e, Finset.sum_add_distrib, ← Finset.mul_sum]
  ring

/-- The completing-the-square identity
`xᵀ A x = a₁₁ (t + cᵀ y / a₁₁)² + yᵀ S y`. -/
theorem complete_square (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) (hA : Aᵀ = A)
    (ha : A 0 0 ≠ 0) (t : ℝ) (y : Fin n → ℝ) :
    (Fin.cons t y : Fin (n + 1) → ℝ) ⬝ᵥ (A *ᵥ Fin.cons t y) =
      A 0 0 * (t + (∑ i, A i.succ 0 * y i) / A 0 0) ^ 2 + y ⬝ᵥ (schur A *ᵥ y) := by
  rw [quad_cons A hA]
  have hS : y ⬝ᵥ (schur A *ᵥ y) =
      y ⬝ᵥ (trailing A *ᵥ y) - (∑ i, A i.succ 0 * y i) ^ 2 / A 0 0 := by
    simp only [schur, trailing, dotProduct, Matrix.mulVec, Matrix.of_apply,
      Matrix.submatrix_apply, sub_mul, Finset.sum_sub_distrib, mul_sub]
    rw [sq, Finset.sum_mul_sum, Finset.sum_div]
    congr 1
    apply Finset.sum_congr rfl; intro i _
    rw [Finset.mul_sum, Finset.sum_div]
    apply Finset.sum_congr rfl; intro j _
    field_simp
  rw [hS]
  field_simp
  ring

lemma schur_transpose (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) (hA : Aᵀ = A) :
    (schur A)ᵀ = schur A := by
  ext i j
  simp only [schur, Matrix.transpose_apply, Matrix.of_apply]
  rw [symm_apply_of_transpose hA j.succ i.succ]
  ring

/-- The Schur complement of an SPD matrix is SPD. -/
theorem schur_posDef (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) (hA : A.PosDef) :
    (schur A).PosDef := by
  rw [posDef_iff_real] at hA ⊢
  obtain ⟨hsymm, hpos⟩ := hA
  have ha : 0 < A 0 0 := by
    have := hpos (Pi.single 0 1) (by simp)
    simpa [dotProduct, Matrix.mulVec, Pi.single_apply] using this
  refine ⟨schur_transpose A hsymm, fun y hy => ?_⟩
  set t := -(∑ i, A i.succ 0 * y i) / A 0 0
  have hx : (Fin.cons t y : Fin (n + 1) → ℝ) ≠ 0 := by
    intro h
    apply hy
    funext i
    have := congrFun h i.succ
    simpa using this
  have h := hpos _ hx
  rw [complete_square A hsymm (ne_of_gt ha)] at h
  have : t + (∑ i, A i.succ 0 * y i) / A 0 0 = 0 := by
    simp only [t]; field_simp; ring
  rw [this] at h
  simpa using h

end Cholesky
