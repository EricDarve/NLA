import src.FactorBounds

/-!
# Vector and matrix norms

* `vnorm v` is the Euclidean norm `‖v‖₂`.
* `norm2 A` is the spectral norm `‖A‖₂`, i.e. mathlib's operator norm of `A` acting on
  Euclidean space (`Matrix.Norms.L2Operator`).
* `frob A` is the Frobenius norm `‖A‖_F`.
* `cond2 A = ‖A‖₂ ‖A⁻¹‖₂`.

Besides the basic properties we prove the facts used in the notes:
`‖A‖₂ ≤ ‖A‖_F`, `|E| ≤ B ⇒ ‖E‖₂ ≤ ‖B‖₂` for `B ≥ 0`, `|aᵢᵢ| ≤ ‖A‖₂`, `tr A ≤ n ‖A‖₂`,
and bounds on `‖A‖₂` and `‖A⁻¹‖₂` from bounds on the quadratic form `xᵀ A x`.
-/
open scoped Matrix.Norms.L2Operator

set_option linter.unusedSectionVars false

namespace Cholesky
open Matrix

variable {m n l : Type*} [Fintype m] [Fintype n] [Fintype l]

/-! ## Vectors -/

/-- The Euclidean norm. -/
noncomputable def vnorm (v : n → ℝ) : ℝ := ‖(WithLp.toLp 2 v : EuclideanSpace ℝ n)‖

lemma vnorm_eq (v : n → ℝ) : vnorm v = Real.sqrt (∑ i, v i ^ 2) := by
  simp [vnorm, EuclideanSpace.norm_eq]

lemma vnorm_nonneg (v : n → ℝ) : 0 ≤ vnorm v := norm_nonneg _

lemma vnorm_sq (v : n → ℝ) : vnorm v ^ 2 = v ⬝ᵥ v := by
  rw [vnorm_eq, Real.sq_sqrt (Finset.sum_nonneg fun i _ => sq_nonneg _)]
  simp [dotProduct, sq]

lemma vnorm_add_le (v w : n → ℝ) : vnorm (v + w) ≤ vnorm v + vnorm w := by
  unfold vnorm
  rw [WithLp.toLp_add]
  exact norm_add_le _ _

lemma vnorm_smul (c : ℝ) (v : n → ℝ) : vnorm (c • v) = |c| * vnorm v := by
  unfold vnorm
  rw [WithLp.toLp_smul, norm_smul, Real.norm_eq_abs]

lemma vnorm_neg (v : n → ℝ) : vnorm (-v) = vnorm v := by
  rw [show -v = (-1 : ℝ) • v by simp, vnorm_smul]; simp

lemma vnorm_sub_le (v w : n → ℝ) : vnorm (v - w) ≤ vnorm v + vnorm w := by
  rw [sub_eq_add_neg]; exact le_trans (vnorm_add_le _ _) (by rw [vnorm_neg])

lemma vnorm_sub_comm (v w : n → ℝ) : vnorm (v - w) = vnorm (w - v) := by
  rw [← vnorm_neg, neg_sub]

lemma vnorm_zero : vnorm (0 : n → ℝ) = 0 := by simp [vnorm]

lemma vnorm_eq_zero_iff (v : n → ℝ) : vnorm v = 0 ↔ v = 0 := by
  unfold vnorm
  rw [norm_eq_zero]
  constructor
  · intro h; exact (WithLp.toLp_eq_zero 2).mp h
  · intro h; rw [h]; simp

lemma vnorm_pos {v : n → ℝ} (hv : v ≠ 0) : 0 < vnorm v :=
  lt_of_le_of_ne (vnorm_nonneg v) (Ne.symm ((vnorm_eq_zero_iff v).not.mpr hv))

/-- Cauchy–Schwarz. -/
lemma abs_dot_le (v w : n → ℝ) : |v ⬝ᵥ w| ≤ vnorm v * vnorm w := by
  have h := abs_real_inner_le_norm (WithLp.toLp 2 v : EuclideanSpace ℝ n) (WithLp.toLp 2 w)
  rw [EuclideanSpace.inner_eq_star_dotProduct] at h
  simpa [vnorm, dotProduct_comm] using h

lemma dot_le (v w : n → ℝ) : v ⬝ᵥ w ≤ vnorm v * vnorm w :=
  le_trans (le_abs_self _) (abs_dot_le v w)

lemma vnorm_abs (v : n → ℝ) : vnorm (fun i => |v i|) = vnorm v := by
  simp [vnorm_eq, sq_abs]

lemma vnorm_le_of_sq {v w : n → ℝ} (h : v ⬝ᵥ v ≤ w ⬝ᵥ w) : vnorm v ≤ vnorm w := by
  rw [← vnorm_sq, ← vnorm_sq] at h
  nlinarith [vnorm_nonneg v, vnorm_nonneg w]

lemma vnorm_single [DecidableEq n] (i : n) : vnorm (Pi.single i (1 : ℝ)) = 1 := by
  rw [vnorm_eq]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hj; simp [hj]
  · simp

/-! ## The spectral norm -/

variable [DecidableEq m] [DecidableEq n] [DecidableEq l]

/-- The spectral norm `‖A‖₂`. -/
noncomputable def norm2 (A : Matrix m n ℝ) : ℝ := ‖A‖

lemma norm2_nonneg (A : Matrix m n ℝ) : 0 ≤ norm2 A := norm_nonneg _

lemma vnorm_mulVec_le (A : Matrix m n ℝ) (v : n → ℝ) :
    vnorm (A *ᵥ v) ≤ norm2 A * vnorm v := by
  have := Matrix.l2_opNorm_mulVec A (WithLp.toLp 2 v)
  simpa [vnorm, norm2] using this

/-- `‖A‖₂ ≤ c` as soon as `‖A x‖ ≤ c ‖x‖` for all `x`. -/
lemma norm2_le_of {A : Matrix m n ℝ} {c : ℝ} (hc : 0 ≤ c)
    (h : ∀ v, vnorm (A *ᵥ v) ≤ c * vnorm v) : norm2 A ≤ c := by
  unfold norm2
  rw [Matrix.l2_opNorm_def]
  apply ContinuousLinearMap.opNorm_le_bound _ hc
  intro x
  have := h (WithLp.ofLp x)
  simpa [vnorm] using this

lemma le_norm2_of {A : Matrix m n ℝ} {v : n → ℝ} (hv : v ≠ 0) {c : ℝ}
    (h : c * vnorm v ≤ vnorm (A *ᵥ v)) : c ≤ norm2 A := by
  have h2 := vnorm_mulVec_le A v
  have hp := vnorm_pos hv
  nlinarith

lemma norm2_add_le (A B : Matrix m n ℝ) : norm2 (A + B) ≤ norm2 A + norm2 B :=
  norm_add_le _ _

lemma norm2_sub_le (A B : Matrix m n ℝ) : norm2 (A - B) ≤ norm2 A + norm2 B :=
  norm_sub_le _ _

lemma norm2_smul (c : ℝ) (A : Matrix m n ℝ) : norm2 (c • A) = |c| * norm2 A := by
  unfold norm2; rw [norm_smul, Real.norm_eq_abs]

lemma norm2_neg (A : Matrix m n ℝ) : norm2 (-A) = norm2 A := norm_neg _

lemma norm2_mul_le (A : Matrix m n ℝ) (B : Matrix n l ℝ) :
    norm2 (A * B) ≤ norm2 A * norm2 B :=
  Matrix.l2_opNorm_mul A B

lemma norm2_transpose (A : Matrix m n ℝ) : norm2 Aᵀ = norm2 A := by
  have := Matrix.l2_opNorm_conjTranspose A
  rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at this

lemma norm2_transpose_mul_self (A : Matrix m n ℝ) : norm2 (Aᵀ * A) = norm2 A ^ 2 := by
  have := Matrix.l2_opNorm_conjTranspose_mul_self A
  rw [Matrix.conjTranspose_eq_transpose_of_trivial] at this
  rw [norm2, this, sq]; rfl

lemma norm2_mul_transpose_self (A : Matrix m n ℝ) : norm2 (A * Aᵀ) = norm2 A ^ 2 := by
  have := norm2_transpose_mul_self Aᵀ
  rwa [Matrix.transpose_transpose, norm2_transpose] at this

lemma norm2_zero : norm2 (0 : Matrix m n ℝ) = 0 := norm_zero

/-! ## The Frobenius norm -/

/-- The Frobenius norm `‖A‖_F`. -/
noncomputable def frob (A : Matrix m n ℝ) : ℝ := Real.sqrt (∑ i, ∑ j, A i j ^ 2)

lemma frob_nonneg (A : Matrix m n ℝ) : 0 ≤ frob A := Real.sqrt_nonneg _

lemma frob_sq (A : Matrix m n ℝ) : frob A ^ 2 = ∑ i, ∑ j, A i j ^ 2 :=
  Real.sq_sqrt (Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => sq_nonneg _)

lemma frob_absMat (A : Matrix m n ℝ) : frob (absMat A) = frob A := by
  simp [frob, sq_abs]

/-- `‖A‖₂ ≤ ‖A‖_F`. -/
theorem norm2_le_frob (A : Matrix m n ℝ) : norm2 A ≤ frob A := by
  apply norm2_le_of (frob_nonneg A)
  intro v
  have hsq : vnorm (A *ᵥ v) ^ 2 ≤ (frob A * vnorm v) ^ 2 := by
    rw [mul_pow, frob_sq, vnorm_sq, vnorm_sq, Finset.sum_mul]
    simp only [dotProduct, Matrix.mulVec]
    apply Finset.sum_le_sum
    intro i _
    have := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun j => A i j) v
    simpa [sq] using this
  exact le_of_pow_le_pow_left₀ two_ne_zero (mul_nonneg (frob_nonneg A) (vnorm_nonneg v)) hsq

/-! ## Entrywise bounds -/

/-- If `|E| ≤ B` entrywise then `‖E‖₂ ≤ ‖B‖₂`. -/
theorem norm2_le_of_abs_le {E B : Matrix m n ℝ} (h : ∀ i j, |E i j| ≤ B i j) :
    norm2 E ≤ norm2 B := by
  apply norm2_le_of (norm2_nonneg B)
  intro v
  have hsq : (E *ᵥ v) ⬝ᵥ (E *ᵥ v) ≤ (B *ᵥ fun j => |v j|) ⬝ᵥ (B *ᵥ fun j => |v j|) := by
    simp only [dotProduct, Matrix.mulVec]
    apply Finset.sum_le_sum
    intro i _
    have h1 : |∑ j, E i j * v j| ≤ ∑ j, B i j * |v j| := by
      refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum fun j _ => ?_)
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_right (h i j) (abs_nonneg _)
    have := mul_self_le_mul_self (abs_nonneg _) h1
    rwa [abs_mul_abs_self] at this
  calc vnorm (E *ᵥ v) ≤ vnorm (B *ᵥ fun j => |v j|) := vnorm_le_of_sq hsq
    _ ≤ norm2 B * vnorm (fun j => |v j|) := vnorm_mulVec_le _ _
    _ = norm2 B * vnorm v := by rw [vnorm_abs]

lemma norm2_le_norm2_absMat (A : Matrix m n ℝ) : norm2 A ≤ norm2 (absMat A) :=
  norm2_le_of_abs_le fun i j => le_of_eq (by simp)

/-- `|aᵢⱼ| ≤ ‖A‖₂`. -/
lemma abs_entry_le_norm2 (A : Matrix m n ℝ) (i : m) (j : n) : |A i j| ≤ norm2 A := by
  have h1 : |A i j| ≤ vnorm (A *ᵥ Pi.single j 1) := by
    rw [vnorm_eq]
    apply Real.abs_le_sqrt
    have : (A *ᵥ Pi.single j 1) i = A i j := by simp [Matrix.mulVec_single]
    rw [← this]
    exact Finset.single_le_sum (f := fun i => (A *ᵥ Pi.single j 1) i ^ 2)
      (fun _ _ => sq_nonneg _) (Finset.mem_univ i)
  have h2 := vnorm_mulVec_le A (Pi.single j 1)
  rw [vnorm_single, mul_one] at h2
  exact le_trans h1 h2

/-- `tr A ≤ n ‖A‖₂`. -/
theorem trace_le (A : Matrix n n ℝ) : A.trace ≤ Fintype.card n * norm2 A := by
  unfold Matrix.trace
  calc ∑ i, A.diag i ≤ ∑ _i : n, norm2 A :=
        Finset.sum_le_sum fun i _ => le_trans (le_abs_self _) (abs_entry_le_norm2 A i i)
    _ = Fintype.card n * norm2 A := by simp

/-! ## Quadratic forms -/

lemma quad_le_norm2 (A : Matrix n n ℝ) (v : n → ℝ) :
    v ⬝ᵥ (A *ᵥ v) ≤ norm2 A * vnorm v ^ 2 := by
  calc v ⬝ᵥ (A *ᵥ v) ≤ vnorm v * vnorm (A *ᵥ v) := dot_le _ _
    _ ≤ vnorm v * (norm2 A * vnorm v) :=
        mul_le_mul_of_nonneg_left (vnorm_mulVec_le A v) (vnorm_nonneg v)
    _ = norm2 A * vnorm v ^ 2 := by ring

lemma abs_quad_le_norm2 (A : Matrix n n ℝ) (v : n → ℝ) :
    |v ⬝ᵥ (A *ᵥ v)| ≤ norm2 A * vnorm v ^ 2 := by
  calc |v ⬝ᵥ (A *ᵥ v)| ≤ vnorm v * vnorm (A *ᵥ v) := abs_dot_le _ _
    _ ≤ vnorm v * (norm2 A * vnorm v) :=
        mul_le_mul_of_nonneg_left (vnorm_mulVec_le A v) (vnorm_nonneg v)
    _ = norm2 A * vnorm v ^ 2 := by ring

lemma dot_mulVec_symm {A : Matrix n n ℝ} (hA : Aᵀ = A) (v w : n → ℝ) :
    v ⬝ᵥ (A *ᵥ w) = (A *ᵥ v) ⬝ᵥ w := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hA]

/-- For a symmetric matrix, `|xᵀ A x| ≤ c ‖x‖²` for all `x` implies `‖A‖₂ ≤ c`. -/
theorem norm2_le_of_quad {A : Matrix n n ℝ} (hA : Aᵀ = A) {c : ℝ} (hc : 0 < c)
    (h : ∀ v, |v ⬝ᵥ (A *ᵥ v)| ≤ c * vnorm v ^ 2) : norm2 A ≤ c := by
  apply norm2_le_of hc.le
  intro v
  -- `2c (c² ‖v‖² - ‖Av‖²) = (Qv)ᵀ P (Qv) + (Pv)ᵀ Q (Pv)` with `P = c - A`, `Q = c + A`.
  set w := A *ᵥ v
  set a := v ⬝ᵥ v
  set b := v ⬝ᵥ w
  set d := w ⬝ᵥ w
  set e := w ⬝ᵥ (A *ᵥ w)
  have hvw : w ⬝ᵥ v = b := by rw [dotProduct_comm]
  have hAv : v ⬝ᵥ (A *ᵥ w) = d := by rw [dot_mulVec_symm hA]
  set w1 := c • v + w
  set w2 := c • v - w
  have hw1 : A *ᵥ w1 = c • w + A *ᵥ w := by
    simp only [w1, Matrix.mulVec_add, Matrix.mulVec_smul]; rfl
  have hw2 : A *ᵥ w2 = c • w - A *ᵥ w := by
    simp only [w2, Matrix.mulVec_sub, Matrix.mulVec_smul]; rfl
  have q1 : w1 ⬝ᵥ w1 = c ^ 2 * a + 2 * c * b + d := by
    simp only [w1, add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul,
      smul_eq_mul, hvw]; ring
  have q2 : w2 ⬝ᵥ w2 = c ^ 2 * a - 2 * c * b + d := by
    simp only [w2, sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul,
      smul_eq_mul, hvw]; ring
  have r1 : w1 ⬝ᵥ (A *ᵥ w1) = c ^ 2 * b + 2 * c * d + e := by
    rw [hw1]
    simp only [w1, add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul,
      smul_eq_mul, hAv]
    rw [show w ⬝ᵥ w = d from rfl]; ring
  have r2 : w2 ⬝ᵥ (A *ᵥ w2) = c ^ 2 * b - 2 * c * d + e := by
    rw [hw2]
    simp only [w2, sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul,
      smul_eq_mul, hAv]
    rw [show w ⬝ᵥ w = d from rfl]; ring
  have p1 := h w1
  have p2 := h w2
  rw [vnorm_sq, q1, r1] at p1
  rw [vnorm_sq, q2, r2] at p2
  have hp1 := (abs_le.mp p1).2
  have hp2 := (abs_le.mp p2).1
  have hd : d ≤ c ^ 2 * a := by nlinarith
  have hsq : vnorm w ^ 2 ≤ (c * vnorm v) ^ 2 := by
    rw [mul_pow, vnorm_sq, vnorm_sq]; exact hd
  exact le_of_pow_le_pow_left₀ two_ne_zero (mul_nonneg hc.le (vnorm_nonneg v)) hsq

/-- For a symmetric positive semidefinite matrix, `xᵀ A x ≤ c ‖x‖²` gives `‖A‖₂ ≤ c`. -/
theorem norm2_le_of_quad_psd {A : Matrix n n ℝ} (hA : Aᵀ = A) {c : ℝ} (hc : 0 < c)
    (hpsd : ∀ v, 0 ≤ v ⬝ᵥ (A *ᵥ v)) (h : ∀ v, v ⬝ᵥ (A *ᵥ v) ≤ c * vnorm v ^ 2) :
    norm2 A ≤ c :=
  norm2_le_of_quad hA hc fun v => abs_le.mpr ⟨by
    have := hpsd v; nlinarith [sq_nonneg (vnorm v), mul_nonneg hc.le (sq_nonneg (vnorm v))],
    h v⟩

/-- A lower bound on the quadratic form controls the inverse:
`μ ‖x‖² ≤ xᵀ A x` gives `‖A⁻¹ y‖ ≤ ‖y‖ / μ`. -/
theorem vnorm_le_of_quad_lower {A : Matrix n n ℝ} {μ : ℝ} (hμ : 0 < μ)
    (h : ∀ v, μ * vnorm v ^ 2 ≤ v ⬝ᵥ (A *ᵥ v)) (v : n → ℝ) :
    vnorm v ≤ vnorm (A *ᵥ v) / μ := by
  rw [le_div_iff₀ hμ]
  have h1 := h v
  have h2 := dot_le v (A *ᵥ v)
  by_cases hv : vnorm v = 0
  · rw [hv]; simp [vnorm_nonneg]
  have hvpos : 0 < vnorm v := lt_of_le_of_ne (vnorm_nonneg v) (Ne.symm hv)
  nlinarith

lemma isUnit_of_quad_lower {A : Matrix n n ℝ} {μ : ℝ} (hμ : 0 < μ)
    (h : ∀ v, μ * vnorm v ^ 2 ≤ v ⬝ᵥ (A *ᵥ v)) : IsUnit A := by
  rw [← Matrix.mulVec_injective_iff_isUnit]
  intro x y hxy
  have := vnorm_le_of_quad_lower hμ h (x - y)
  rw [Matrix.mulVec_sub, hxy, sub_self, vnorm_zero, zero_div] at this
  exact sub_eq_zero.mp ((vnorm_eq_zero_iff _).mp (le_antisymm this (vnorm_nonneg _)))

theorem norm2_inv_le_of_quad {A : Matrix n n ℝ} {μ : ℝ} (hμ : 0 < μ)
    (h : ∀ v, μ * vnorm v ^ 2 ≤ v ⬝ᵥ (A *ᵥ v)) : norm2 A⁻¹ ≤ 1 / μ := by
  have hu := isUnit_of_quad_lower hμ h
  have hdet : IsUnit A.det := (Matrix.isUnit_iff_isUnit_det A).mp hu
  apply norm2_le_of (by positivity)
  intro w
  have := vnorm_le_of_quad_lower hμ h (A⁻¹ *ᵥ w)
  rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec] at this
  rw [one_div_mul_eq_div]; exact this

/-! ## Condition number -/

/-- The condition number `κ₂(A) = ‖A‖₂ ‖A⁻¹‖₂`. -/
noncomputable def cond2 (A : Matrix n n ℝ) : ℝ := norm2 A * norm2 A⁻¹

lemma cond2_nonneg (A : Matrix n n ℝ) : 0 ≤ cond2 A :=
  mul_nonneg (norm2_nonneg _) (norm2_nonneg _)

/-- Orthogonal matrices preserve Euclidean length. -/
lemma vnorm_orthogonal {U : Matrix n n ℝ} (hU : Uᵀ * U = 1) (v : n → ℝ) :
    vnorm (U *ᵥ v) = vnorm v := by
  have : (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = v ⬝ᵥ v := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec, hU,
      Matrix.one_mulVec]
  rw [← sq_eq_sq₀ (vnorm_nonneg _) (vnorm_nonneg _), vnorm_sq, vnorm_sq, this]

lemma norm2_orthogonal_conj {U A : Matrix n n ℝ} (hU : Uᵀ * U = 1) (hU' : U * Uᵀ = 1) :
    norm2 (U * A * Uᵀ) = norm2 A := by
  apply le_antisymm
  · apply norm2_le_of (norm2_nonneg A)
    intro v
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, vnorm_orthogonal hU]
    calc vnorm (A *ᵥ (Uᵀ *ᵥ v)) ≤ norm2 A * vnorm (Uᵀ *ᵥ v) := vnorm_mulVec_le _ _
      _ = norm2 A * vnorm v := by
          rw [vnorm_orthogonal (by rw [Matrix.transpose_transpose]; exact hU')]
  · apply norm2_le_of (norm2_nonneg _)
    intro v
    have : A *ᵥ v = Uᵀ *ᵥ ((U * A * Uᵀ) *ᵥ (U *ᵥ v)) := by
      simp only [Matrix.mulVec_mulVec, ← Matrix.mul_assoc, hU, Matrix.one_mul]
      rw [Matrix.mul_assoc, hU, Matrix.mul_one]
    rw [this, vnorm_orthogonal (by rw [Matrix.transpose_transpose]; exact hU')]
    calc vnorm ((U * A * Uᵀ) *ᵥ (U *ᵥ v)) ≤ norm2 (U * A * Uᵀ) * vnorm (U *ᵥ v) :=
          vnorm_mulVec_le _ _
      _ = norm2 (U * A * Uᵀ) * vnorm v := by rw [vnorm_orthogonal hU]

/-- `κ₂(U A Uᵀ) = κ₂(A)` for orthogonal `U`, in particular for permutation matrices. -/
theorem cond2_orthogonal_conj {U A : Matrix n n ℝ} (hU : Uᵀ * U = 1) (hU' : U * Uᵀ = 1) :
    cond2 (U * A * Uᵀ) = cond2 A := by
  unfold cond2
  have hinv : (U * A * Uᵀ)⁻¹ = U * A⁻¹ * Uᵀ := by
    have hUinv : U⁻¹ = Uᵀ := Matrix.inv_eq_left_inv hU
    have hUTinv : Uᵀ⁻¹ = U := Matrix.inv_eq_left_inv hU'
    rw [Matrix.mul_inv_rev, Matrix.mul_inv_rev, hUinv, hUTinv, Matrix.mul_assoc]
  rw [hinv, norm2_orthogonal_conj hU hU', norm2_orthogonal_conj hU hU']

end Cholesky
