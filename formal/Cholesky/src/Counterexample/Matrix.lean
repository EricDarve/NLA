import src.Counterexample.Params
import src.Counterexample.Blocks
import src.BlockStep

/-!
# The matrix `A_{9k}`

With `Q = H_k / √k`, `J` the all-ones matrix,

`T = [[t J_{6k,2k}], [(3/2) [I_k Q]]]`, `G = (9/4) [[I, Q], [Qᵀ, I]]`,
`A_{9k} = [[4 I_{7k}, 2T], [2Tᵀ, G + η J_{2k} + δ I_{2k}]]`.

We prove `Tᵀ T = G + (9/8) η J` (`T_gram`), so that the exact Schur complement of
`4 I_{7k}` is `S = δ I - (η/8) J` (`schur_eq`); its extreme eigenvalues are `(3/4) k η`
and `k η` (`S_quad_bounds`, `S_ones`, `S_perp`); and `A_{9k}` is SPD (`A_posDef`).
We also define the matrix `Ŝ` of the notes and prove its eigenvalue facts.
-/
namespace Cholesky
namespace CE
namespace Par
open Matrix

variable (P : Par)

/-- `Q = H_k / √k`. -/
noncomputable def Q : Matrix (Fin P.k) (Fin P.k) ℝ := hadQ P.q

lemma Q_mul_transpose : P.Q * P.Qᵀ = 1 := hadQ_mul_transpose P.q
lemma Q_transpose_mul : P.Qᵀ * P.Q = 1 := hadQ_transpose_mul P.q
lemma Q_transpose : P.Qᵀ = P.Q := hadQ_transpose P.q

lemma Q_entry (i j : Fin P.k) : P.Q i j = 1 / P.K ∨ P.Q i j = -(1 / P.K) := by
  rcases hadQ_entry P.q i j with h | h
  · left; rw [Q, h, P.zpow_neg_q]
  · right; rw [Q, h, P.zpow_neg_q]

lemma Q_sq (i j : Fin P.k) : P.Q i j ^ 2 = 1 / P.k := by
  rw [P.k_real]
  rcases P.Q_entry i j with h | h <;> rw [h] <;> field_simp

/-- The coupling block `T`: `6k` rows of `t`, then `(3/2) [I Q]`. -/
noncomputable def T : Matrix (Fin (6 * P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  blk (P.t • ones (6 * P.k) P.k) (P.t • ones (6 * P.k) P.k) ((3 / 2 : ℝ) • 1) ((3 / 2 : ℝ) • P.Q)

/-- `G = (9/4) [[I, Q], [Qᵀ, I]]`. -/
noncomputable def G : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  (9 / 4 : ℝ) • blk 1 P.Q P.Qᵀ 1

/-- `J_{2k}`. -/
def J2 : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ := ones _ _

/-- The trailing block `G + η J + δ I`. -/
noncomputable def C : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  P.G + P.η • P.J2 + P.δ • 1

/-- The exact Schur complement `S = δ I - (η/8) J`. -/
noncomputable def S : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  P.δ • 1 - (P.η / 8) • P.J2

/-- The stored matrix `A_{9k}`. -/
noncomputable def A : Matrix (Fin ((6 * P.k + P.k) + (P.k + P.k)))
    (Fin ((6 * P.k + P.k) + (P.k + P.k))) ℝ :=
  blk ((4 : ℝ) • 1) ((2 : ℝ) • P.T) ((2 : ℝ) • P.Tᵀ) P.C

/-- `a = (k + 9/8) η`. -/
noncomputable def a : ℝ := (P.k + 9 / 8) * P.η

/-- The trailing matrix `Ŝ` produced by the rounded leading elimination. -/
noncomputable def Shat : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  blk (P.a • 1 - (P.η / 8) • ones P.k P.k) (-P.η • ones P.k P.k) (-P.η • ones P.k P.k)
    (P.a • 1 - (P.η / 8) • ones P.k P.k)

/-! ## The Gram identity -/

/-- `Tᵀ T = G + (9/8) η J`. -/
theorem T_gram : P.Tᵀ * P.T = P.G + (9 / 8 * P.η) • P.J2 := by
  have h6 : 6 * (P.k : ℝ) * P.t ^ 2 = 9 / 8 * P.η := by
    rw [P.t_sq_eq, η]; ring
  rw [T, blk_transpose, blk_mul, G, J2, ones_eq_blk, smul_blk, smul_blk, blk_add]
  simp only [Matrix.transpose_smul, Matrix.transpose_one, Matrix.smul_mul, Matrix.mul_smul,
    Matrix.one_mul, Matrix.mul_one, P.Q_transpose_mul, smul_smul, ones_transpose_mul_ones]
  congr 1 <;> ext i j <;>
    simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, ones, Matrix.of_apply] <;>
    push_cast <;> linear_combination h6

lemma C_symm : P.Cᵀ = P.C := by
  have hG : P.Gᵀ = P.G := by
    rw [G, Matrix.transpose_smul, blk_transpose, Matrix.transpose_one, Matrix.transpose_transpose]
  have hJ : P.J2ᵀ = P.J2 := by ext i j; rfl
  simp only [C, Matrix.transpose_add, Matrix.transpose_smul, hG, hJ, Matrix.transpose_one]

/-- The trailing block is `Tᵀ T + S`, so `S` is the exact Schur complement of `4 I`. -/
theorem C_eq : P.C = P.Tᵀ * P.T + P.S := by
  rw [T_gram, C, S]
  ext i j
  simp only [Matrix.add_apply, Matrix.smul_apply, Matrix.sub_apply, smul_eq_mul]
  ring

/-- `S = C - (2Tᵀ)(4I)⁻¹(2T)`. -/
theorem schur_eq : P.C - ((2 : ℝ) • P.Tᵀ) *
    ((1 / 4 : ℝ) • (1 : Matrix (Fin (6 * P.k + P.k)) (Fin (6 * P.k + P.k)) ℝ)) *
    ((2 : ℝ) • P.T) = P.S := by
  have : ((2 : ℝ) • P.Tᵀ) *
      ((1 / 4 : ℝ) • (1 : Matrix (Fin (6 * P.k + P.k)) (Fin (6 * P.k + P.k)) ℝ)) *
      ((2 : ℝ) • P.T) = P.Tᵀ * P.T := by
    simp only [Matrix.smul_mul, Matrix.mul_smul, Matrix.mul_one, smul_smul]
    norm_num
  rw [this, C_eq, add_sub_cancel_left]

/-! ## The exact Schur complement -/

lemma S_quad (v : Fin (P.k + P.k) → ℝ) :
    v ⬝ᵥ (P.S *ᵥ v) = P.δ * vnorm v ^ 2 - P.η / 8 * (∑ i, v i) ^ 2 := by
  rw [S, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
    dotProduct_sub, dotProduct_smul, dotProduct_smul, J2, dot_ones_mulVec, vnorm_sq]
  simp [sq]

/-- `(3/4) k η ‖v‖² ≤ vᵀ S v ≤ k η ‖v‖²`. -/
theorem S_quad_bounds (v : Fin (P.k + P.k) → ℝ) :
    3 / 4 * P.k * P.η * vnorm v ^ 2 ≤ v ⬝ᵥ (P.S *ᵥ v) ∧
      v ⬝ᵥ (P.S *ᵥ v) ≤ P.k * P.η * vnorm v ^ 2 := by
  rw [S_quad]
  have h := sum_sq_le_card v
  have hη := P.η_pos
  push_cast at h
  have hδ : P.δ = P.k * P.η := rfl
  constructor
  · rw [hδ]; nlinarith [sq_nonneg (∑ i, v i)]
  · rw [hδ]; nlinarith [sq_nonneg (∑ i, v i)]

/-- `S 1 = (3/4) k η 1`: the smallest eigenvalue. -/
theorem S_ones : P.S *ᵥ (fun _ => 1) = (3 / 4 * P.k * P.η) • (fun _ => (1 : ℝ)) := by
  rw [S, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec, J2,
    ones_mulVec]
  ext i
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul, mul_one, δ]
  push_cast
  ring

/-- `S v = k η v` when `∑ v = 0`: the largest eigenvalue. -/
theorem S_perp (v : Fin (P.k + P.k) → ℝ) (hv : ∑ i, v i = 0) : P.S *ᵥ v = (P.k * P.η) • v := by
  ext i
  have : (P.J2 *ᵥ v) i = 0 := by rw [J2, ones_mulVec]; exact hv
  simp only [S, Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec, Pi.sub_apply,
    Pi.smul_apply, smul_eq_mul, this, mul_zero, sub_zero, δ]

lemma S_symm : P.Sᵀ = P.S := by
  ext i j; simp [S, J2, ones, Matrix.one_apply, eq_comm]

/-! ## Positive definiteness -/

/-- `vᵀ A v = ‖2x + T y‖² + yᵀ S y` for `v = (x, y)`. -/
theorem A_quad (x : Fin (6 * P.k + P.k) → ℝ) (y : Fin (P.k + P.k) → ℝ) :
    Fin.append x y ⬝ᵥ (P.A *ᵥ Fin.append x y) =
      ((2 : ℝ) • x + P.T *ᵥ y) ⬝ᵥ ((2 : ℝ) • x + P.T *ᵥ y) + y ⬝ᵥ (P.S *ᵥ y) := by
  rw [A, blk_mulVec_append, append_dotProduct, C_eq]
  simp only [Matrix.smul_mulVec, Matrix.one_mulVec, Matrix.add_mulVec, ← Matrix.mulVec_mulVec,
    dotProduct_add, add_dotProduct, dotProduct_smul, smul_dotProduct, smul_eq_mul]
  rw [Matrix.dotProduct_mulVec y P.Tᵀ x, Matrix.vecMul_transpose,
    Matrix.dotProduct_mulVec y P.Tᵀ (P.T *ᵥ y), Matrix.vecMul_transpose,
    dotProduct_comm (P.T *ᵥ y) x]
  ring

lemma A_symm : P.Aᵀ = P.A := by
  rw [A, blk_transpose, Matrix.transpose_smul, Matrix.transpose_one, Matrix.transpose_smul,
    Matrix.transpose_smul, Matrix.transpose_transpose, C_symm]

/-- `A_{9k}` is symmetric positive definite. -/
theorem A_posDef : P.A.PosDef := by
  rw [posDef_iff_real]
  refine ⟨P.A_symm, fun v hv => ?_⟩
  rw [← append_vL_vR v, A_quad]
  have hS := (P.S_quad_bounds (vR v)).1
  have hk : (0 : ℝ) < P.k := by have := P.k_ge; positivity
  have hη := P.η_pos
  have h1 : 0 ≤ ((2 : ℝ) • vL v + P.T *ᵥ vR v) ⬝ᵥ ((2 : ℝ) • vL v + P.T *ᵥ vR v) := by
    rw [← vnorm_sq]; positivity
  by_cases hy : vR v = 0
  · have hx : vL v ≠ 0 := by
      intro hx; apply hv; rw [← append_vL_vR v, hx, hy, append_zero]
    rw [hy, Matrix.mulVec_zero, add_zero, Matrix.mulVec_zero, dotProduct_zero, add_zero,
      smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul, ← vnorm_sq]
    have := vnorm_pos hx
    positivity
  · have := vnorm_pos hy
    have : 0 < 3 / 4 * P.k * P.η * vnorm (vR v) ^ 2 := by positivity
    linarith

/-! ## The matrix `Ŝ` -/

lemma Shat_symm : P.Shatᵀ = P.Shat := by
  rw [Shat, blk_transpose]
  congr 1 <;> ext i j <;> simp [ones, Matrix.one_apply, eq_comm]

/-- `vᵀ Ŝ v = a ‖v‖² - (η/8)(σ₁² + σ₂²) - 2η σ₁ σ₂` with `σᵢ` the block sums. -/
theorem Shat_quad (v : Fin (P.k + P.k) → ℝ) :
    v ⬝ᵥ (P.Shat *ᵥ v) = P.a * vnorm v ^ 2 - P.η / 8 * ((∑ i, vL v i) ^ 2 + (∑ i, vR v i) ^ 2)
      - 2 * P.η * (∑ i, vL v i) * (∑ i, vR v i) := by
  conv_lhs => rw [← append_vL_vR v]
  rw [Shat, blk_mulVec_append, append_dotProduct, vnorm_sq_split, vnorm_sq, vnorm_sq]
  simp only [Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_add,
    dotProduct_sub, dotProduct_smul, smul_eq_mul, dot_ones_mulVec]
  ring

/-- The all-ones vector is an eigenvector for `-(k - 9) η / 8`. -/
lemma const_eq_append (c : ℝ) :
    (fun _ : Fin (P.k + P.k) => c) = Fin.append (fun _ => c) (fun _ => c) := by
  ext i; refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    simp only [Fin.append_left, Fin.append_right]

theorem Shat_ones : P.Shat *ᵥ (fun _ => 1) = (-((P.k : ℝ) - 9) * P.η / 8) • (fun _ => (1 : ℝ)) := by
  rw [const_eq_append, Shat, blk_mulVec_append, append_smul]
  congr 1 <;> ext i <;>
    simp only [Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
      ones_mulVec, Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, Finset.sum_const,
      Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one, a] <;> ring

/-- `(1, -1)` is an eigenvector for `(15k + 9) η / 8`. -/
theorem Shat_alt : P.Shat *ᵥ Fin.append (fun _ => 1) (fun _ => -1) =
    ((15 * (P.k : ℝ) + 9) * P.η / 8) • Fin.append (fun _ => (1 : ℝ)) (fun _ => -1) := by
  rw [Shat, blk_mulVec_append, append_smul]
  congr 1 <;> ext i <;>
    simp only [Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
      ones_mulVec, Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, Finset.sum_const,
      Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one, a, mul_neg] <;> ring

/-- `-(k-9)η/8 ‖v‖² ≤ vᵀ Ŝ v ≤ (15k+9)η/8 ‖v‖²`. -/
theorem Shat_quad_bounds (v : Fin (P.k + P.k) → ℝ) :
    -((P.k : ℝ) - 9) * P.η / 8 * vnorm v ^ 2 ≤ v ⬝ᵥ (P.Shat *ᵥ v) ∧
      v ⬝ᵥ (P.Shat *ᵥ v) ≤ (15 * P.k + 9) * P.η / 8 * vnorm v ^ 2 := by
  rw [Shat_quad]
  set s1 := ∑ i, vL v i
  set s2 := ∑ i, vR v i
  have hv := vnorm_sq_split v
  -- `(s₁ ± s₂)² ≤ 2k ‖v‖²`
  have hp : (s1 + s2) ^ 2 ≤ (P.k + P.k : ℕ) * vnorm v ^ 2 := by
    have := sum_sq_le_card v
    rwa [Fin.sum_univ_add] at this
  have hm : (s1 - s2) ^ 2 ≤ (P.k + P.k : ℕ) * vnorm v ^ 2 := by
    have := sum_sq_le_card (Fin.append (vL v) (-vR v))
    rw [Fin.sum_univ_add, vnorm_append_sq, vnorm_neg] at this
    simp only [Fin.append_left, Fin.append_right, Pi.neg_apply, Finset.sum_neg_distrib] at this
    rw [hv]; linarith
  push_cast at hp hm
  have hη := P.η_pos
  have hdec : P.η / 8 * (s1 ^ 2 + s2 ^ 2) + 2 * P.η * s1 * s2 =
      9 * P.η / 16 * (s1 + s2) ^ 2 - 7 * P.η / 16 * (s1 - s2) ^ 2 := by ring
  constructor
  · rw [sub_sub, hdec, a]
    nlinarith [sq_nonneg (s1 - s2), mul_le_mul_of_nonneg_left hp (by positivity : (0 : ℝ) ≤ 9 * P.η / 16)]
  · rw [sub_sub, hdec, a]
    nlinarith [sq_nonneg (s1 + s2), mul_le_mul_of_nonneg_left hm (by positivity : (0 : ℝ) ≤ 7 * P.η / 16)]

/-- `‖Ŝ‖₂ ≤ (15k + 9) η / 8`. -/
theorem norm2_Shat_le : norm2 P.Shat ≤ (15 * P.k + 9) * P.η / 8 := by
  have hk : (64 : ℝ) ≤ P.k := by exact_mod_cast P.k_ge
  have hη := P.η_pos
  apply norm2_le_of_quad P.Shat_symm (by positivity)
  intro v
  obtain ⟨h1, h2⟩ := P.Shat_quad_bounds v
  rw [abs_le]
  constructor
  · nlinarith [sq_nonneg (vnorm v), mul_nonneg (by positivity : (0 : ℝ) ≤ 16 * P.k * P.η / 8)
      (sq_nonneg (vnorm v))]
  · exact h2

end Par
end CE
end Cholesky
