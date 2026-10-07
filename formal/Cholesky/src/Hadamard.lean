import Mathlib

/-!
# Walsh–Hadamard matrices

`hadamard m` is the Walsh–Hadamard matrix of size `2^m`, defined recursively as in the
notes: `H₁ = [1]`, `H_{2k} = [[H_k, H_k], [H_k, -H_k]]`. We prove that its entries are
`±1`, that it is symmetric, that `H Hᵀ = 2^m I`, that its first row is all ones, and that
every other row sums to zero. The scaled matrix `Q = H / √k` with `k = 4^q` is then
orthogonal with entries `±2^{-q}`.
-/
namespace Cholesky
open Matrix

lemma two_pow_succ_eq (m : ℕ) : 2 ^ m + 2 ^ m = 2 ^ (m + 1) := by ring

/-- `Fin (2^m) ⊕ Fin (2^m) ≃ Fin (2^(m+1))`, first block first. -/
def hadEquiv (m : ℕ) : Fin (2 ^ m) ⊕ Fin (2 ^ m) ≃ Fin (2 ^ (m + 1)) :=
  finSumFinEquiv.trans (finCongr (two_pow_succ_eq m))

/-- The Walsh–Hadamard matrix of size `2^m`. -/
def hadamard : (m : ℕ) → Matrix (Fin (2 ^ m)) (Fin (2 ^ m)) ℝ
  | 0 => 1
  | m + 1 => Matrix.reindex (hadEquiv m) (hadEquiv m)
      (Matrix.fromBlocks (hadamard m) (hadamard m) (hadamard m) (-hadamard m))

lemma hadamard_succ_apply (m : ℕ) (i j : Fin (2 ^ (m + 1))) :
    hadamard (m + 1) i j =
      Matrix.fromBlocks (hadamard m) (hadamard m) (hadamard m) (-hadamard m)
        ((hadEquiv m).symm i) ((hadEquiv m).symm j) := rfl

theorem hadamard_entry : ∀ (m : ℕ) (i j : Fin (2 ^ m)),
    hadamard m i j = 1 ∨ hadamard m i j = -1
  | 0, ⟨i, hi⟩, ⟨j, hj⟩ => by
    left
    have hi' : i = 0 := by simpa using hi
    have hj' : j = 0 := by simpa using hj
    subst hi' hj'
    rfl
  | m + 1, i, j => by
    rw [hadamard_succ_apply]
    rcases (hadEquiv m).symm i with a | a <;> rcases (hadEquiv m).symm j with b | b <;>
      simp only [Matrix.fromBlocks_apply₁₁, Matrix.fromBlocks_apply₁₂,
        Matrix.fromBlocks_apply₂₁, Matrix.fromBlocks_apply₂₂, Matrix.neg_apply]
    · exact hadamard_entry m a b
    · exact hadamard_entry m a b
    · exact hadamard_entry m a b
    · rcases hadamard_entry m a b with h | h <;> rw [h] <;> norm_num

theorem hadamard_transpose : ∀ m : ℕ, (hadamard m)ᵀ = hadamard m
  | 0 => Matrix.transpose_one
  | m + 1 => by
    rw [hadamard, Matrix.transpose_reindex, Matrix.fromBlocks_transpose,
      Matrix.transpose_neg, hadamard_transpose m]

theorem hadamard_mul_transpose : ∀ m : ℕ,
    hadamard m * (hadamard m)ᵀ = ((2 : ℝ) ^ m) • (1 : Matrix (Fin (2 ^ m)) (Fin (2 ^ m)) ℝ)
  | 0 => by simp [hadamard]
  | m + 1 => by
    have ih := hadamard_mul_transpose m
    rw [hadamard_transpose m] at ih
    rw [hadamard_transpose (m + 1)]
    rw [show hadamard (m + 1) = Matrix.reindex (hadEquiv m) (hadEquiv m)
      (Matrix.fromBlocks (hadamard m) (hadamard m) (hadamard m) (-hadamard m)) from rfl]
    simp only [Matrix.reindex_apply, Matrix.submatrix_mul_equiv,
      Matrix.fromBlocks_multiply, ih, Matrix.mul_neg, Matrix.neg_mul, neg_neg]
    rw [show (2 : ℝ) ^ m • (1 : Matrix (Fin (2 ^ m)) (Fin (2 ^ m)) ℝ) +
        (2 : ℝ) ^ m • (1 : Matrix (Fin (2 ^ m)) (Fin (2 ^ m)) ℝ) =
        (2 : ℝ) ^ (m + 1) • (1 : Matrix (Fin (2 ^ m)) (Fin (2 ^ m)) ℝ) by
      rw [← add_smul]; congr 1; ring]
    rw [add_neg_cancel,
      show (0 : Matrix (Fin (2 ^ m)) (Fin (2 ^ m)) ℝ) = (2 : ℝ) ^ (m + 1) • 0 by simp,
      ← Matrix.fromBlocks_smul, Matrix.fromBlocks_one]
    ext i j
    simp [Matrix.one_apply, (hadEquiv m).symm.injective.eq_iff]

lemma zero_lt_two_pow (m : ℕ) : 0 < 2 ^ m := pow_pos (by norm_num) m

/-- The index `0` of `Fin (2^m)`. -/
def had0 (m : ℕ) : Fin (2 ^ m) := ⟨0, zero_lt_two_pow m⟩

lemma hadEquiv_symm_had0 (m : ℕ) : (hadEquiv m).symm (had0 (m + 1)) = Sum.inl (had0 m) := by
  simp only [hadEquiv, Equiv.symm_trans_apply, finCongr_symm, finCongr_apply]
  have : (Fin.cast (two_pow_succ_eq m).symm (had0 (m + 1))) =
      Fin.castAdd (2 ^ m) (had0 m) := by
    ext; simp [had0]
  rw [this, finSumFinEquiv_symm_apply_castAdd]

/-- The first row of `H` is all ones. -/
theorem hadamard_row0 : ∀ (m : ℕ) (j : Fin (2 ^ m)), hadamard m (had0 m) j = 1
  | 0, ⟨j, hj⟩ => by
    have hj' : j = 0 := by simpa using hj
    subst hj'; rfl
  | m + 1, j => by
    rw [hadamard_succ_apply, hadEquiv_symm_had0]
    rcases (hadEquiv m).symm j with b | b
    · exact hadamard_row0 m b
    · exact hadamard_row0 m b

/-- Every row other than the first sums to zero. -/
theorem hadamard_row_sum (m : ℕ) (i : Fin (2 ^ m)) (hi : i ≠ had0 m) :
    ∑ j, hadamard m i j = 0 := by
  have h := congrFun (congrFun (hadamard_mul_transpose m) i) (had0 m)
  rw [Matrix.mul_apply, Matrix.smul_apply, Matrix.one_apply_ne hi, smul_zero] at h
  simpa [hadamard_row0] using h

/-! ## The scaled orthogonal matrix `Q = H / √k`, `k = 4^q` -/

/-- `Q = 2^{-q} H_{4^q}`. -/
noncomputable def hadQ (q : ℕ) : Matrix (Fin (2 ^ (2 * q))) (Fin (2 ^ (2 * q))) ℝ :=
  ((2 : ℝ) ^ (-(q : ℤ))) • hadamard (2 * q)

lemma hadQ_apply (q : ℕ) (i j : Fin (2 ^ (2 * q))) :
    hadQ q i j = (2 : ℝ) ^ (-(q : ℤ)) * hadamard (2 * q) i j := rfl

theorem hadQ_transpose (q : ℕ) : (hadQ q)ᵀ = hadQ q := by
  rw [hadQ, Matrix.transpose_smul, hadamard_transpose]

theorem hadQ_mul_transpose (q : ℕ) : hadQ q * (hadQ q)ᵀ = 1 := by
  rw [hadQ, Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, hadamard_mul_transpose,
    smul_smul, smul_smul]
  have : (2 : ℝ) ^ (-(q : ℤ)) * (2 : ℝ) ^ (-(q : ℤ)) * 2 ^ (2 * q) = 1 := by
    rw [← zpow_natCast, ← zpow_add₀ (by norm_num), ← zpow_add₀ (by norm_num)]
    push_cast
    ring_nf
    simp
  rw [this, one_smul]

theorem hadQ_transpose_mul (q : ℕ) : (hadQ q)ᵀ * hadQ q = 1 := by
  rw [hadQ_transpose]
  have := hadQ_mul_transpose q
  rwa [hadQ_transpose] at this

theorem hadQ_entry (q : ℕ) (i j : Fin (2 ^ (2 * q))) :
    hadQ q i j = (2 : ℝ) ^ (-(q : ℤ)) ∨ hadQ q i j = -(2 : ℝ) ^ (-(q : ℤ)) := by
  rw [hadQ_apply]
  rcases hadamard_entry (2 * q) i j with h | h <;> rw [h] <;> simp

lemma hadQ_sq (q : ℕ) (i j : Fin (2 ^ (2 * q))) :
    hadQ q i j ^ 2 = (2 : ℝ) ^ (-(2 * q : ℤ)) := by
  rcases hadQ_entry q i j with h | h <;> rw [h]
  · rw [← zpow_natCast, ← _root_.zpow_mul]; congr 1; push_cast; ring
  · rw [neg_sq, ← zpow_natCast, ← _root_.zpow_mul]; congr 1; push_cast; ring

theorem hadQ_row_sum (q : ℕ) (i : Fin (2 ^ (2 * q))) (hi : i ≠ had0 (2 * q)) :
    ∑ j, hadQ q i j = 0 := by
  simp only [hadQ_apply, ← Finset.mul_sum, hadamard_row_sum (2 * q) i hi, mul_zero]

end Cholesky
