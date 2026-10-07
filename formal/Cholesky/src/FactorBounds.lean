import src.BlockStep

/-!
# Why the factor entries stay controlled

For `A = L Lᵀ`: `aᵢᵢ = ∑ₖ lᵢₖ²`, hence `|lᵢₖ| ≤ √aᵢᵢ`, and by Cauchy–Schwarz
`(|L| |L|ᵀ)ᵢⱼ = ∑ₖ |lᵢₖ lⱼₖ| ≤ √(aᵢᵢ aⱼⱼ)`.
-/
namespace Cholesky
open Matrix

variable {n : ℕ}

/-- The entrywise absolute value `|A|`. -/
def absMat {m p : Type*} (A : Matrix m p ℝ) : Matrix m p ℝ := Matrix.of fun i j => |A i j|

@[simp] lemma absMat_apply {m p : Type*} (A : Matrix m p ℝ) (i : m) (j : p) :
    absMat A i j = |A i j| := rfl

lemma absMat_transpose {m p : Type*} (A : Matrix m p ℝ) : absMat Aᵀ = (absMat A)ᵀ := rfl

lemma absMat_nonneg {m p : Type*} (A : Matrix m p ℝ) (i : m) (j : p) : 0 ≤ absMat A i j :=
  abs_nonneg _

/-- `aᵢᵢ = ∑ₖ lᵢₖ²`. -/
theorem diag_eq_sum_sq (L : Matrix (Fin n) (Fin n) ℝ) (i : Fin n) :
    (L * Lᵀ) i i = ∑ k, L i k ^ 2 := by
  simp [Matrix.mul_apply, sq]

/-- Only the entries with `k ≤ i` contribute. -/
theorem diag_eq_sum_sq_le (L : Matrix (Fin n) (Fin n) ℝ) (hL : IsLowerTriangular L)
    (i : Fin n) : (L * Lᵀ) i i = ∑ k ∈ Finset.univ.filter (· ≤ i), L i k ^ 2 := by
  rw [diag_eq_sum_sq, Finset.sum_filter]
  apply Finset.sum_congr rfl
  intro k _
  split_ifs with hk
  · rfl
  · rw [hL i k (lt_of_not_ge hk)]; ring

/-- `|lᵢₖ| ≤ √aᵢᵢ`. -/
theorem abs_entry_le_sqrt_diag (L : Matrix (Fin n) (Fin n) ℝ) (i k : Fin n) :
    |L i k| ≤ Real.sqrt ((L * Lᵀ) i i) := by
  rw [diag_eq_sum_sq]
  apply Real.abs_le_sqrt
  exact Finset.single_le_sum (f := fun k => L i k ^ 2) (fun _ _ => sq_nonneg _)
    (Finset.mem_univ k)

/-- `(|L| |L|ᵀ)ᵢⱼ = ∑ₖ |lᵢₖ lⱼₖ| ≤ √(aᵢᵢ aⱼⱼ)`. -/
theorem absProd_le (L : Matrix (Fin n) (Fin n) ℝ) (i j : Fin n) :
    (absMat L * (absMat L)ᵀ) i j = ∑ k, |L i k * L j k| ∧
      (absMat L * (absMat L)ᵀ) i j ≤ Real.sqrt ((L * Lᵀ) i i * (L * Lᵀ) j j) := by
  have heq : (absMat L * (absMat L)ᵀ) i j = ∑ k, |L i k * L j k| := by
    simp [Matrix.mul_apply, abs_mul]
  refine ⟨heq, ?_⟩
  rw [heq, diag_eq_sum_sq, diag_eq_sum_sq]
  apply Real.le_sqrt_of_sq_le
  have h := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun k => |L i k|) (fun k => |L j k|)
  simp only [sq_abs] at h
  simpa [abs_mul] using h

end Cholesky
