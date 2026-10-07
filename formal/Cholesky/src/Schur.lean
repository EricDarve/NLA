import Mathlib

/-!
# The exact and rounded Schur complements

These are statements about explicit real matrices. The name `roundedSchur` denotes
 the formula in the notes, not the output of a formalized rounding algorithm.
-/
namespace Cholesky
open Matrix Finset

noncomputable def rankOneShift (m : ℕ) (d c : ℝ) : Matrix (Fin m) (Fin m) ℝ :=
  Matrix.diagonal (fun _ => d) - Matrix.of (fun _ _ => c)

lemma rankOneShift_quadratic (m : ℕ) (d c : ℝ) (v : Fin m → ℝ) :
    star v ⬝ᵥ (rankOneShift m d c *ᵥ v) =
      d * (∑ i, v i ^ 2) - c * (∑ i, v i) ^ 2 := by
  rw [rankOneShift, Matrix.sub_mulVec, dotProduct_sub]
  have hdiag : Matrix.diagonal (fun _ : Fin m => d) *ᵥ v = fun i => d * v i := by
    ext i
    exact Matrix.mulVec_diagonal _ _ _
  rw [hdiag]
  simp only [dotProduct, star_trivial, Matrix.mulVec, Matrix.of_apply]
  simp_rw [← Finset.mul_sum]
  rw [← Finset.sum_mul]
  simp_rw [show ∀ i, v i * (d * v i) = d * v i ^ 2 by intro i; ring]
  rw [← Finset.mul_sum]
  ring

lemma sum_square_bound (m : ℕ) (v : Fin m → ℝ) :
    (∑ i, v i) ^ 2 ≤ (m : ℝ) * ∑ i, v i ^ 2 := by
  simpa using Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset (Fin m))
    (fun _ => (1 : ℝ)) v

lemma rankOneShift_posDef (m : ℕ) (d c : ℝ) (hc : 0 ≤ c)
    (hgap : c * m < d) : (rankOneShift m d c).PosDef := by
  constructor
  · ext i j
    simp [rankOneShift, Matrix.conjTranspose_apply, Matrix.diagonal_apply, eq_comm]
  · intro v hv
    rw [rankOneShift_quadratic]
    have hsq : 0 < ∑ i, v i ^ 2 := by
      have h := Matrix.dotProduct_star_self_pos_iff.mpr hv
      simpa [dotProduct, pow_two] using h
    have hCS := mul_le_mul_of_nonneg_left (sum_square_bound m v) hc
    nlinarith [mul_pos (sub_pos.mpr hgap) hsq]

/-- The exact Schur complement in the 9k-by-9k construction is positive definite. -/
theorem exactSchur_posDef (k : ℕ) (η : ℝ) (hk : 0 < k) (hη : 0 < η) :
    (rankOneShift (2 * k) (k * η) (η / 8)).PosDef := by
  apply rankOneShift_posDef _ _ _ (by positivity)
  have hkR : (0 : ℝ) < k := by exact_mod_cast hk
  push_cast
  nlinarith [mul_pos hkR hη]

noncomputable def roundedSchur (k : ℕ) (η : ℝ) :
    Matrix (Fin k ⊕ Fin k) (Fin k ⊕ Fin k) ℝ :=
  Matrix.fromBlocks
    (rankOneShift k (((k : ℝ) + 9 / 8) * η) (η / 8))
    (Matrix.of (fun _ _ => -η))
    (Matrix.of (fun _ _ => -η))
    (rankOneShift k (((k : ℝ) + 9 / 8) * η) (η / 8))

lemma rankOneShift_ones (m : ℕ) (d c : ℝ) :
    rankOneShift m d c *ᵥ (fun _ => 1) = fun _ => d - c * m := by
  rw [rankOneShift, Matrix.sub_mulVec]
  ext i
  simp only [Pi.sub_apply, Matrix.mulVec_diagonal]
  simp [Matrix.mulVec, dotProduct, mul_comm]

/-- The all-ones vector is an eigenvector for the negative eigenvalue in the notes. -/
theorem roundedSchur_ones (k : ℕ) (η : ℝ) :
    roundedSchur k η *ᵥ (fun _ => 1) = fun _ => -((k : ℝ) - 9) * η / 8 := by
  rw [roundedSchur, Matrix.fromBlocks_mulVec]
  ext i
  cases i <;> simp [Matrix.mulVec, dotProduct,
    rankOneShift, Matrix.diagonal_apply, Finset.sum_sub_distrib] <;> ring

/-- For k > 9, this explicit real matrix cannot be a Gram matrix. -/
theorem roundedSchur_not_posSemidef (k : ℕ) (η : ℝ) (hk : 9 < k) (hη : 0 < η) :
    ¬ (roundedSchur k η).PosSemidef := by
  intro h
  have hq := h.2 (fun _ => 1)
  rw [roundedSchur_ones] at hq
  simp only [dotProduct, Pi.star_apply, star_one, one_mul,
    Finset.sum_const, Finset.card_univ, Fintype.card_sum, Fintype.card_fin,
    nsmul_eq_mul] at hq
  have hkR : (9 : ℝ) < k := by exact_mod_cast hk
  push_cast at hq
  have hn : -((k : ℝ) - 9) * η / 8 < 0 := by nlinarith [mul_pos (sub_pos.mpr hkR) hη]
  have hp : (0 : ℝ) < k + k := by linarith
  have := mul_neg_of_pos_of_neg hp hn
  linarith

/-- The block matrix whose exact Schur complement is `S`. -/
noncomputable def liftedMatrix {m n : ℕ} (T : Matrix (Fin m) (Fin n) ℝ)
    (S : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin m ⊕ Fin n) (Fin m ⊕ Fin n) ℝ :=
  Matrix.fromBlocks 4 ((2 : ℝ) • T) ((2 : ℝ) • Tᵀ) (Tᵀ * T + S)

lemma liftedMatrix_quadratic {m n : ℕ} (T : Matrix (Fin m) (Fin n) ℝ)
    (S : Matrix (Fin n) (Fin n) ℝ) (x : Fin m → ℝ) (y : Fin n → ℝ) :
    star (Sum.elim x y) ⬝ᵥ (liftedMatrix T S *ᵥ Sum.elim x y) =
      ((2 : ℝ) • x + T *ᵥ y) ⬝ᵥ ((2 : ℝ) • x + T *ᵥ y) + star y ⬝ᵥ (S *ᵥ y) := by
  simp only [liftedMatrix, star_trivial, Matrix.fromBlocks_mulVec,
    Sum.elim_comp_inl, Sum.elim_comp_inr, sumElim_dotProduct_sumElim,
    Matrix.smul_mulVec, Matrix.add_mulVec, ← Matrix.mulVec_mulVec,
    dotProduct_add, add_dotProduct, dotProduct_smul, smul_dotProduct]
  rw [Matrix.dotProduct_mulVec y Tᵀ (T *ᵥ y), Matrix.vecMul_transpose,
    Matrix.dotProduct_mulVec y Tᵀ x, Matrix.vecMul_transpose]
  simp only [Matrix.ofNat_mulVec,
    dotProduct_smul, smul_eq_mul]
  rw [dotProduct_comm (T *ᵥ y) x]
  ring

/-- Positive definiteness of the full block matrix follows from that of `S`. -/
theorem liftedMatrix_posDef {m n : ℕ} (T : Matrix (Fin m) (Fin n) ℝ)
    (S : Matrix (Fin n) (Fin n) ℝ) (hS : S.PosDef) :
    (liftedMatrix T S).PosDef := by
  constructor
  · apply Matrix.IsHermitian.fromBlocks (Matrix.isHermitian_natCast 4)
    · simp
    · exact (Matrix.isHermitian_transpose_mul_self T).add hS.1
  · intro v hv
    let x := v ∘ Sum.inl
    let y := v ∘ Sum.inr
    have hvsplit : v = Sum.elim x y := (Sum.elim_comp_inl_inr v).symm
    rw [hvsplit, liftedMatrix_quadratic]
    have hnorm : 0 ≤ ((2 : ℝ) • x + T *ᵥ y) ⬝ᵥ ((2 : ℝ) • x + T *ᵥ y) := by
      simpa using dotProduct_star_self_nonneg ((2 : ℝ) • x + T *ᵥ y)
    by_cases hy : y = 0
    · have hx : x ≠ 0 := by
        intro hx
        apply hv
        rw [hvsplit, hx, hy]
        ext i
        cases i <;> rfl
      have hxx : 0 < x ⬝ᵥ x := by
        simpa using Matrix.dotProduct_star_self_pos_iff.mpr hx
      simp only [hy, Matrix.mulVec_zero, add_zero, star_zero, dotProduct_zero]
      rw [smul_dotProduct, dotProduct_smul]
      change 0 < 2 * (2 * (x ⬝ᵥ x))
      nlinarith
    · exact add_pos_of_nonneg_of_pos hnorm (hS.2 y hy)

/-- Any coupling block `T` gives an SPD matrix with the exact Schur complement
used in the counterexample. No representability claim is made here. -/
theorem constructedMatrix_posDef (k : ℕ) (η : ℝ)
    (T : Matrix (Fin (7 * k)) (Fin (2 * k)) ℝ) (hk : 0 < k) (hη : 0 < η) :
    (liftedMatrix T (rankOneShift (2 * k) (k * η) (η / 8))).PosDef :=
  liftedMatrix_posDef _ _ (exactSchur_posDef k η hk hη)

end Cholesky
