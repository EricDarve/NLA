import Mathlib

/-!
# Packed storage of the factors

After elimination the array holds `L` strictly below the diagonal (its unit diagonal is
implicit) and `U` on and above it. `lowerOf M` and `upperOf M` extract the two factors
from a packed array `M`. For a matrix of size `m + 1` the product `L U` has the block form

`(L U)₀ⱼ = m₀ⱼ`, `(L U)ᵢ₀ = mᵢ₀ m₀₀`, `(L U)ᵢⱼ = mᵢ₀ m₀ⱼ + (L' U')ᵢⱼ` for `i, j ≥ 1`,

where `L', U'` are the factors packed in the trailing block (`mul_zero_apply`,
`mul_succ_zero`, `mul_succ_succ`).
-/
namespace LU
open Matrix

variable {𝕜 : Type*} [Field 𝕜] {n m : ℕ}

/-- `L` is unit lower triangular. -/
def IsUnitLower (L : Matrix (Fin n) (Fin n) 𝕜) : Prop :=
  (∀ i, L i i = 1) ∧ ∀ i j, i < j → L i j = 0

/-- `U` is upper triangular. -/
def IsUpper (U : Matrix (Fin n) (Fin n) 𝕜) : Prop :=
  ∀ i j, j < i → U i j = 0

/-- The unit lower triangular factor packed strictly below the diagonal of `M`. -/
def lowerOf (M : Matrix (Fin n) (Fin n) 𝕜) : Matrix (Fin n) (Fin n) 𝕜 :=
  Matrix.of fun i j => if i = j then 1 else if j < i then M i j else 0

/-- The upper triangular factor packed on and above the diagonal of `M`. -/
def upperOf (M : Matrix (Fin n) (Fin n) 𝕜) : Matrix (Fin n) (Fin n) 𝕜 :=
  Matrix.of fun i j => if i ≤ j then M i j else 0

lemma lowerOf_isUnitLower (M : Matrix (Fin n) (Fin n) 𝕜) : IsUnitLower (lowerOf M) :=
  ⟨fun i => by simp [lowerOf], fun i j hij => by
    simp [lowerOf, ne_of_lt hij, not_lt.mpr hij.le]⟩

lemma upperOf_isUpper (M : Matrix (Fin n) (Fin n) 𝕜) : IsUpper (upperOf M) :=
  fun i j hij => by simp [upperOf, not_le.mpr hij]

lemma lowerOf_apply_of_lt (M : Matrix (Fin n) (Fin n) 𝕜) {i j : Fin n} (h : j < i) :
    lowerOf M i j = M i j := by
  simp [lowerOf, ne_of_gt h, h]

lemma upperOf_apply_of_le (M : Matrix (Fin n) (Fin n) 𝕜) {i j : Fin n} (h : i ≤ j) :
    upperOf M i j = M i j := by
  simp [upperOf, h]

lemma det_of_isUnitLower {L : Matrix (Fin n) (Fin n) 𝕜} (hL : IsUnitLower L) : L.det = 1 := by
  rw [Matrix.det_of_lowerTriangular L (fun i j h => hL.2 i j h)]
  simp [hL.1]

lemma det_of_isUpper {U : Matrix (Fin n) (Fin n) 𝕜} (hU : IsUpper U) :
    U.det = ∏ i, U i i :=
  Matrix.det_of_upperTriangular (fun i j h => hU i j h)

/-! ## Block structure -/

section Blocks
variable (M : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜)

@[simp] lemma lowerOf_zero_zero : lowerOf M 0 0 = 1 := by simp [lowerOf]

@[simp] lemma lowerOf_zero_succ (j : Fin m) : lowerOf M 0 j.succ = 0 := by
  simp [lowerOf, (Fin.succ_ne_zero j).symm, not_lt.mpr (Fin.zero_le _)]

@[simp] lemma lowerOf_succ_zero (i : Fin m) : lowerOf M i.succ 0 = M i.succ 0 :=
  lowerOf_apply_of_lt M (Fin.succ_pos i)

@[simp] lemma lowerOf_succ_succ (i j : Fin m) :
    lowerOf M i.succ j.succ = lowerOf (M.submatrix Fin.succ Fin.succ) i j := by
  simp [lowerOf, Fin.succ_inj, Fin.succ_lt_succ_iff]

@[simp] lemma upperOf_zero (j : Fin (m + 1)) : upperOf M 0 j = M 0 j :=
  upperOf_apply_of_le M (Fin.zero_le j)

@[simp] lemma upperOf_succ_zero (i : Fin m) : upperOf M i.succ 0 = 0 := by
  simp [upperOf, not_le.mpr (Fin.succ_pos i)]

@[simp] lemma upperOf_succ_succ (i j : Fin m) :
    upperOf M i.succ j.succ = upperOf (M.submatrix Fin.succ Fin.succ) i j := by
  simp [upperOf, Fin.succ_le_succ_iff]

lemma mul_zero_apply (j : Fin (m + 1)) : (lowerOf M * upperOf M) 0 j = M 0 j := by
  rw [Matrix.mul_apply, Fin.sum_univ_succ]
  simp

lemma mul_succ_zero (i : Fin m) : (lowerOf M * upperOf M) i.succ 0 = M i.succ 0 * M 0 0 := by
  rw [Matrix.mul_apply, Fin.sum_univ_succ]
  simp

lemma mul_succ_succ (i j : Fin m) :
    (lowerOf M * upperOf M) i.succ j.succ = M i.succ 0 * M 0 j.succ +
      (lowerOf (M.submatrix Fin.succ Fin.succ) * upperOf (M.submatrix Fin.succ Fin.succ)) i j := by
  rw [Matrix.mul_apply, Fin.sum_univ_succ]
  simp [Matrix.mul_apply]

end Blocks

/-! ## Entrywise absolute values -/

/-- The entrywise absolute value `|A|` of a real matrix. -/
def absMat {p q : Type*} (A : Matrix p q ℝ) : Matrix p q ℝ := Matrix.of fun i j => |A i j|

@[simp] lemma absMat_apply {p q : Type*} (A : Matrix p q ℝ) (i : p) (j : q) :
    absMat A i j = |A i j| := rfl

lemma absMat_nonneg {p q : Type*} (A : Matrix p q ℝ) (i : p) (j : q) : 0 ≤ absMat A i j :=
  abs_nonneg _

lemma absProd_nonneg {p q r : Type*} [Fintype q] (A : Matrix p q ℝ) (B : Matrix q r ℝ)
    (i : p) (j : r) : 0 ≤ (absMat A * absMat B) i j :=
  Finset.sum_nonneg fun _ _ => mul_nonneg (abs_nonneg _) (abs_nonneg _)

lemma absMat_lowerOf {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) :
    absMat (lowerOf M) = lowerOf (absMat M) := by
  ext i j
  simp only [absMat_apply, lowerOf, Matrix.of_apply]
  split_ifs <;> simp

lemma absMat_upperOf {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) :
    absMat (upperOf M) = upperOf (absMat M) := by
  ext i j
  simp only [absMat_apply, upperOf, Matrix.of_apply]
  split_ifs <;> simp

/-- `|(A B)ᵢⱼ| ≤ (|A| |B|)ᵢⱼ`. -/
lemma abs_mul_apply_le {p q r : Type*} [Fintype q] (A : Matrix p q ℝ) (B : Matrix q r ℝ)
    (i : p) (j : r) : |(A * B) i j| ≤ (absMat A * absMat B) i j := by
  rw [Matrix.mul_apply, Matrix.mul_apply]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (le_of_eq ?_)
  simp [abs_mul]

end LU
