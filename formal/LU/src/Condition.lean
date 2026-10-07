import src.ForwardError

/-!
# The matrix condition number

* `one_le_norm_mul_norm_inv`, `one_le_cond2`: `κ(A) = ‖A‖ ‖A⁻¹‖ ≥ 1`, since
  `1 = ‖I‖ ≤ ‖A‖ ‖A⁻¹‖` (in any normed ring with `‖1‖ = 1`, i.e. any induced norm).
* Singular values: `σᵢ = √λᵢ(AᵀA)` (`singVal`). `norm2_eq_sigmaMax`: `‖A‖₂ = σ_max`;
  `norm2_inv_eq`: `‖A⁻¹‖₂ = 1/σ_min`; `cond2_eq_sigma`: **`κ₂(A) = σ_max/σ_min`**.
  Geometrically, `σ_max` and `σ_min` are the largest and smallest stretching of a unit vector
  (`sigmaMax_isGreatest`, `sigmaMin_isLeast`).
* `diagEps_*`: for `A = diag(1, ε)`, `b = (1, 0)` and `δb = (0, ε)`, the solution changes by
  `δx = (0, 1)`: a relative change `ε` in the data gives a relative change `1` in the
  solution, and `κ₂(A) = 1/ε`.
* `cond2_perm_mul`: `κ₂(P A) = κ₂(A)`; `gram_perm_mul`: `(PA)ᵀ(PA) = AᵀA`, so row permutations
  preserve the singular values (`singVal_perm_mul`).
-/
open scoped Matrix.Norms.L2Operator
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {n : ℕ}

/-! ## `κ ≥ 1` -/

theorem one_le_norm_mul_norm_inv {R : Type*} [NormedRing R] [NormOneClass R] (a : Rˣ) :
    1 ≤ ‖(a : R)‖ * ‖(↑a⁻¹ : R)‖ := by
  calc (1 : ℝ) = ‖(1 : R)‖ := norm_one.symm
    _ = ‖(a : R) * ↑a⁻¹‖ := by rw [Units.mul_inv]
    _ ≤ ‖(a : R)‖ * ‖(↑a⁻¹ : R)‖ := norm_mul_le _ _

/-- `κ₂(A) ≥ 1`. -/
theorem one_le_cond2 [NeZero n] {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) :
    1 ≤ cond2 A := by
  calc (1 : ℝ) = norm2 (1 : Matrix (Fin n) (Fin n) ℝ) := norm2_one.symm
    _ = norm2 (A * A⁻¹) := by rw [Matrix.mul_nonsing_inv _ hA]
    _ ≤ norm2 A * norm2 A⁻¹ := norm2_mul_le _ _

/-! ## Spectral facts -/

section Spectral
variable {M : Matrix (Fin n) (Fin n) ℝ} (hM : M.IsHermitian)

/-- The orthogonal matrix of eigenvectors. -/
noncomputable def eigQ : Matrix (Fin n) (Fin n) ℝ := (hM.eigenvectorUnitary : Matrix _ _ ℝ)

lemma eigQ_star : star (eigQ hM) = (eigQ hM)ᵀ := by
  rw [Matrix.star_eq_conjTranspose, Matrix.conjTranspose_eq_transpose_of_trivial]

lemma eigQ_mul_transpose : eigQ hM * (eigQ hM)ᵀ = 1 := by
  rw [← eigQ_star]
  exact Matrix.mem_unitaryGroup_iff.mp hM.eigenvectorUnitary.2

lemma spectral_real : M = eigQ hM * diagonal hM.eigenvalues * (eigQ hM)ᵀ := by
  conv_lhs => rw [hM.spectral_theorem]
  rw [← eigQ_star]
  rfl

theorem quad_eq_sum_eigenvalues (x : Fin n → ℝ) :
    x ⬝ᵥ (M *ᵥ x) = ∑ i, hM.eigenvalues i * ((eigQ hM)ᵀ *ᵥ x) i ^ 2 := by
  conv_lhs => rw [spectral_real hM]
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
    ← Matrix.mulVec_transpose]
  simp only [dotProduct, Matrix.mulVec_diagonal]
  apply Finset.sum_congr rfl
  intro i _
  ring

lemma eigQ_dot (x : Fin n → ℝ) :
    ((eigQ hM)ᵀ *ᵥ x) ⬝ᵥ ((eigQ hM)ᵀ *ᵥ x) = x ⬝ᵥ x := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.transpose_transpose,
    Matrix.mulVec_mulVec, eigQ_mul_transpose, Matrix.one_mulVec]

lemma eigenvector_ne_zero (i : Fin n) : (⇑(hM.eigenvectorBasis i) : Fin n → ℝ) ≠ 0 := by
  intro h
  apply hM.eigenvectorBasis.orthonormal.ne_zero i
  ext j
  exact congrFun h j

lemma quad_eigenvector (i : Fin n) :
    (⇑(hM.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ (M *ᵥ ⇑(hM.eigenvectorBasis i)) =
      hM.eigenvalues i * ((⇑(hM.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ
        ⇑(hM.eigenvectorBasis i)) := by
  rw [hM.mulVec_eigenvectorBasis i, dotProduct_smul, smul_eq_mul]

end Spectral

/-! ## Singular values -/

section Singular
variable (A : Matrix (Fin n) (Fin n) ℝ)

lemma gram_isHermitian : (Aᵀ * A).IsHermitian := by
  rw [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial, Matrix.transpose_mul,
    Matrix.transpose_transpose]

/-- `xᵀ (AᵀA) x = ‖A x‖²`. -/
lemma gram_quad (x : Fin n → ℝ) : x ⬝ᵥ ((Aᵀ * A) *ᵥ x) = vnorm (A *ᵥ x) ^ 2 := by
  rw [vnorm_sq, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
    Matrix.transpose_transpose, dotProduct_comm]

/-- The eigenvalues of `AᵀA`. -/
noncomputable def gramEig (i : Fin n) : ℝ := (gram_isHermitian A).eigenvalues i

/-- The singular values `σᵢ = √λᵢ(AᵀA)`. -/
noncomputable def singVal (i : Fin n) : ℝ := Real.sqrt (gramEig A i)

lemma gramEig_nonneg (i : Fin n) : 0 ≤ gramEig A i := by
  have hH := gram_isHermitian A
  have hv := eigenvector_ne_zero hH i
  have hq := quad_eigenvector hH i
  rw [gram_quad] at hq
  have hvv : 0 < (⇑(hH.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ ⇑(hH.eigenvectorBasis i) := by
    simpa using Matrix.dotProduct_star_self_pos_iff.mpr hv
  have : 0 ≤ hH.eigenvalues i * ((⇑(hH.eigenvectorBasis i) : Fin n → ℝ) ⬝ᵥ
      ⇑(hH.eigenvectorBasis i)) := by rw [← hq]; exact sq_nonneg _
  exact nonneg_of_mul_nonneg_left this hvv

lemma singVal_sq (i : Fin n) : singVal A i ^ 2 = gramEig A i :=
  Real.sq_sqrt (gramEig_nonneg A i)

/-- `‖A v‖ = σᵢ ‖v‖` for the `i`-th eigenvector `v` of `AᵀA`. -/
lemma vnorm_mul_eigenvector (i : Fin n) :
    vnorm (A *ᵥ ⇑((gram_isHermitian A).eigenvectorBasis i)) =
      singVal A i * vnorm (⇑((gram_isHermitian A).eigenvectorBasis i)) := by
  have hH := gram_isHermitian A
  have hq := quad_eigenvector hH i
  rw [gram_quad, ← vnorm_sq] at hq
  have : vnorm (A *ᵥ ⇑(hH.eigenvectorBasis i)) ^ 2 =
      (singVal A i * vnorm ⇑(hH.eigenvectorBasis i)) ^ 2 := by
    rw [mul_pow, singVal_sq]; exact hq
  exact (sq_eq_sq₀ (vnorm_nonneg _) (mul_nonneg (Real.sqrt_nonneg _) (vnorm_nonneg _))).mp this

/-- `‖A x‖² = ∑ σᵢ² zᵢ²` with `‖x‖² = ∑ zᵢ²`. -/
lemma vnorm_mul_sq (x : Fin n → ℝ) :
    vnorm (A *ᵥ x) ^ 2 = ∑ i, singVal A i ^ 2 *
      ((eigQ (gram_isHermitian A))ᵀ *ᵥ x) i ^ 2 := by
  rw [← gram_quad, quad_eq_sum_eigenvalues (gram_isHermitian A)]
  simp only [singVal_sq]
  rfl

lemma vnorm_sq_eq_sum (x : Fin n → ℝ) :
    vnorm x ^ 2 = ∑ i, ((eigQ (gram_isHermitian A))ᵀ *ᵥ x) i ^ 2 := by
  rw [vnorm_sq, ← eigQ_dot (gram_isHermitian A) x, dotProduct]
  simp [sq]

variable [NeZero n]

/-- The largest singular value. -/
noncomputable def sigmaMax : ℝ := Finset.univ.sup' Finset.univ_nonempty (singVal A)

/-- The smallest singular value. -/
noncomputable def sigmaMin : ℝ := Finset.univ.inf' Finset.univ_nonempty (singVal A)

lemma singVal_le_sigmaMax (i : Fin n) : singVal A i ≤ sigmaMax A :=
  Finset.le_sup' (singVal A) (Finset.mem_univ i)

lemma sigmaMin_le_singVal (i : Fin n) : sigmaMin A ≤ singVal A i :=
  Finset.inf'_le (singVal A) (Finset.mem_univ i)

lemma sigmaMin_nonneg : 0 ≤ sigmaMin A := by
  obtain ⟨i, -, hi⟩ := Finset.exists_mem_eq_inf' Finset.univ_nonempty (singVal A)
  rw [sigmaMin, hi]; exact Real.sqrt_nonneg _

/-- `σ_min ‖x‖ ≤ ‖A x‖ ≤ σ_max ‖x‖`. -/
theorem stretch_bounds (x : Fin n → ℝ) :
    sigmaMin A * vnorm x ≤ vnorm (A *ᵥ x) ∧ vnorm (A *ᵥ x) ≤ sigmaMax A * vnorm x := by
  have hmin := sigmaMin_nonneg A
  have hmax : 0 ≤ sigmaMax A := le_trans hmin (le_trans (sigmaMin_le_singVal A 0)
    (singVal_le_sigmaMax A 0))
  have hs : ∀ i, 0 ≤ singVal A i := fun i => Real.sqrt_nonneg _
  constructor
  · have : (sigmaMin A * vnorm x) ^ 2 ≤ vnorm (A *ᵥ x) ^ 2 := by
      rw [mul_pow, vnorm_sq_eq_sum A x, vnorm_mul_sq, Finset.mul_sum]
      exact Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_right
        (pow_le_pow_left₀ hmin (sigmaMin_le_singVal A i) 2) (sq_nonneg _)
    exact le_of_pow_le_pow_left₀ two_ne_zero (vnorm_nonneg _) this
  · have : vnorm (A *ᵥ x) ^ 2 ≤ (sigmaMax A * vnorm x) ^ 2 := by
      rw [mul_pow, vnorm_sq_eq_sum A x, vnorm_mul_sq, Finset.mul_sum]
      exact Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_right
        (pow_le_pow_left₀ (hs i) (singVal_le_sigmaMax A i) 2) (sq_nonneg _)
    exact le_of_pow_le_pow_left₀ two_ne_zero (mul_nonneg hmax (vnorm_nonneg _)) this

lemma exists_sigmaMax : ∃ v : Fin n → ℝ, v ≠ 0 ∧ vnorm (A *ᵥ v) = sigmaMax A * vnorm v := by
  obtain ⟨i, -, hi⟩ := Finset.exists_mem_eq_sup' Finset.univ_nonempty (singVal A)
  exact ⟨_, eigenvector_ne_zero (gram_isHermitian A) i, by
    rw [vnorm_mul_eigenvector, sigmaMax, hi]⟩

lemma exists_sigmaMin : ∃ v : Fin n → ℝ, v ≠ 0 ∧ vnorm (A *ᵥ v) = sigmaMin A * vnorm v := by
  obtain ⟨i, -, hi⟩ := Finset.exists_mem_eq_inf' Finset.univ_nonempty (singVal A)
  exact ⟨_, eigenvector_ne_zero (gram_isHermitian A) i, by
    rw [vnorm_mul_eigenvector, sigmaMin, hi]⟩

/-- `‖A‖₂ = σ_max`. -/
theorem norm2_eq_sigmaMax : norm2 A = sigmaMax A := by
  have hmax : 0 ≤ sigmaMax A := le_trans (sigmaMin_nonneg A) (le_trans
    (sigmaMin_le_singVal A 0) (singVal_le_sigmaMax A 0))
  apply le_antisymm (norm2_le_of hmax fun v => (stretch_bounds A v).2)
  obtain ⟨v, hv, h⟩ := exists_sigmaMax A
  exact le_norm2_of hv h.ge

/-- The geometric meaning: `σ_max` is the largest stretching of a unit vector. -/
theorem sigmaMax_isGreatest :
    IsGreatest {s | ∃ x : Fin n → ℝ, vnorm x = 1 ∧ s = vnorm (A *ᵥ x)} (sigmaMax A) := by
  obtain ⟨v, hv, h⟩ := exists_sigmaMax A
  have hvpos := vnorm_pos hv
  refine ⟨⟨(1 / vnorm v) • v, ?_, ?_⟩, ?_⟩
  · rw [vnorm_smul, abs_of_pos (by positivity)]; field_simp
  · rw [Matrix.mulVec_smul, vnorm_smul, abs_of_pos (by positivity), h]; field_simp
  · rintro s ⟨x, hx, rfl⟩
    have := (stretch_bounds A x).2
    rwa [hx, mul_one] at this

/-- The geometric meaning: `σ_min` is the smallest stretching of a unit vector. -/
theorem sigmaMin_isLeast :
    IsLeast {s | ∃ x : Fin n → ℝ, vnorm x = 1 ∧ s = vnorm (A *ᵥ x)} (sigmaMin A) := by
  obtain ⟨v, hv, h⟩ := exists_sigmaMin A
  have hvpos := vnorm_pos hv
  refine ⟨⟨(1 / vnorm v) • v, ?_, ?_⟩, ?_⟩
  · rw [vnorm_smul, abs_of_pos (by positivity)]; field_simp
  · rw [Matrix.mulVec_smul, vnorm_smul, abs_of_pos (by positivity), h]; field_simp
  · rintro s ⟨x, hx, rfl⟩
    have := (stretch_bounds A x).1
    rwa [hx, mul_one] at this

variable {A}

lemma sigmaMin_pos (hA : IsUnit A.det) : 0 < sigmaMin A := by
  obtain ⟨v, hv, h⟩ := exists_sigmaMin A
  rcases eq_or_lt_of_le (sigmaMin_nonneg A) with h0 | h0
  · exfalso
    rw [← h0, zero_mul, vnorm_eq_zero_iff] at h
    have := congrArg (fun w => A⁻¹ *ᵥ w) h
    simp only [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hA, Matrix.one_mulVec,
      Matrix.mulVec_zero] at this
    exact hv this
  · exact h0

/-- `‖A⁻¹‖₂ = 1/σ_min`. -/
theorem norm2_inv_eq (hA : IsUnit A.det) : norm2 A⁻¹ = 1 / sigmaMin A := by
  have hpos := sigmaMin_pos hA
  apply le_antisymm
  · apply norm2_le_of (by positivity)
    intro y
    have h := (stretch_bounds A (A⁻¹ *ᵥ y)).1
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hA, Matrix.one_mulVec] at h
    rw [one_div_mul_eq_div, le_div_iff₀ hpos, mul_comm]
    exact h
  · obtain ⟨v, hv, h⟩ := exists_sigmaMin A
    have hAv : A *ᵥ v ≠ 0 := by
      intro h0
      have := congrArg (fun w => A⁻¹ *ᵥ w) h0
      simp only [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hA, Matrix.one_mulVec,
        Matrix.mulVec_zero] at this
      exact hv this
    apply le_norm2_of hAv
    rw [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hA, Matrix.one_mulVec, h]
    field_simp
    rfl

/-- **`κ₂(A) = σ_max / σ_min`.** -/
theorem cond2_eq_sigma (hA : IsUnit A.det) : cond2 A = sigmaMax A / sigmaMin A := by
  rw [cond2, norm2_eq_sigmaMax, norm2_inv_eq hA]
  ring

end Singular

/-! ## Row permutations -/

lemma gram_perm_mul (σ : Equiv.Perm (Fin n)) (A : Matrix (Fin n) (Fin n) ℝ) :
    (σ.permMatrix ℝ * A)ᵀ * (σ.permMatrix ℝ * A) = Aᵀ * A := by
  rw [Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc (σ.permMatrix ℝ)ᵀ,
    permMatrix_transpose_mul, Matrix.one_mul]

lemma eigenvalues_congr {M N : Matrix (Fin n) (Fin n) ℝ} (h : M = N) (hM : M.IsHermitian)
    (hN : N.IsHermitian) : hM.eigenvalues = hN.eigenvalues := by
  subst h; rfl

/-- Row permutations preserve the singular values. -/
theorem singVal_perm_mul (σ : Equiv.Perm (Fin n)) (A : Matrix (Fin n) (Fin n) ℝ) :
    singVal (σ.permMatrix ℝ * A) = singVal A := by
  ext i
  simp only [singVal, gramEig]
  rw [eigenvalues_congr (gram_perm_mul σ A) (gram_isHermitian _) (gram_isHermitian A)]

/-- **`κ₂(P A) = κ₂(A)`.** -/
theorem cond2_perm_mul (σ : Equiv.Perm (Fin n)) (A : Matrix (Fin n) (Fin n) ℝ) :
    cond2 (σ.permMatrix ℝ * A) = cond2 A := by
  have hinv : (σ.permMatrix ℝ * A)⁻¹ = A⁻¹ * (σ.permMatrix ℝ)ᵀ := by
    rw [Matrix.mul_inv_rev, Matrix.inv_eq_left_inv (permMatrix_transpose_mul σ)]
  have h1 : (σ.permMatrix ℝ)ᵀᵀ * (σ.permMatrix ℝ)ᵀ = 1 := by
    rw [Matrix.transpose_transpose]; exact permMatrix_mul_transpose σ
  have h2 : (σ.permMatrix ℝ)ᵀ * (σ.permMatrix ℝ)ᵀᵀ = 1 := by
    rw [Matrix.transpose_transpose]; exact permMatrix_transpose_mul σ
  unfold cond2
  rw [norm2_perm_mul, hinv, norm2_mul_orthogonal h1 h2]

/-! ## The example `diag(1, ε)` -/

lemma vnorm_single_val (j : Fin n) (a : ℝ) : vnorm (Pi.single j a) = |a| := by
  rw [vnorm_eq, Finset.sum_eq_single j (fun i _ hi => by simp [hi]) (by simp)]
  simp [Real.sqrt_sq_eq_abs]

/-- `‖diag(a, c)‖₂ = max(|a|, |c|)`. -/
lemma norm2_diagonal_two (a c : ℝ) : norm2 (Matrix.diagonal ![a, c]) = max |a| |c| := by
  have hm : 0 ≤ max |a| |c| := le_trans (abs_nonneg a) (le_max_left _ _)
  apply le_antisymm
  · apply norm2_le_of hm
    intro v
    rw [show max |a| |c| * vnorm v = vnorm (max |a| |c| • v) by
      rw [vnorm_smul, abs_of_nonneg hm]]
    apply vnorm_le_of_sq
    have ha : a * a ≤ max |a| |c| * max |a| |c| := by
      rw [← abs_mul_abs_self a]; exact mul_self_le_mul_self (abs_nonneg a) (le_max_left _ _)
    have hc : c * c ≤ max |a| |c| * max |a| |c| := by
      rw [← abs_mul_abs_self c]; exact mul_self_le_mul_self (abs_nonneg c) (le_max_right _ _)
    simp only [Matrix.mulVec_diagonal, dotProduct, Fin.sum_univ_two, Pi.smul_apply, smul_eq_mul,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_fin_one]
    nlinarith [mul_le_mul_of_nonneg_right ha (mul_self_nonneg (v 0)),
      mul_le_mul_of_nonneg_right hc (mul_self_nonneg (v 1))]
  · have hv : ∀ j : Fin 2, (Pi.single j (1 : ℝ) : Fin 2 → ℝ) ≠ 0 := fun j h0 => by
      have := congrFun h0 j; simp at this
    apply max_le
    · apply le_norm2_of (hv 0)
      rw [vnorm_single, mul_one, Matrix.mulVec_single_one]
      rw [show (Matrix.diagonal ![a, c]).col 0 = Pi.single 0 a by
        ext i; fin_cases i <;> simp [Matrix.diagonal]]
      rw [vnorm_single_val]
    · apply le_norm2_of (hv 1)
      rw [vnorm_single, mul_one, Matrix.mulVec_single_one]
      rw [show (Matrix.diagonal ![a, c]).col 1 = Pi.single 1 c by
        ext i; fin_cases i <;> simp [Matrix.diagonal]]
      rw [vnorm_single_val]

/-- `A = diag(1, ε)`. -/
noncomputable def diagEps (ε : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := Matrix.diagonal ![1, ε]

/-- `x = (1, 0)` solves `A x = b = (1, 0)`, and `δx = (0, 1)` solves `A δx = δb = (0, ε)`. -/
theorem diagEps_solve (ε : ℝ) : diagEps ε *ᵥ ![1, 0] = ![1, 0] ∧
    diagEps ε *ᵥ (![1, 0] + ![0, 1]) = ![1, 0] + ![0, ε] := by
  constructor <;> ext i <;> fin_cases i <;> simp [diagEps, Matrix.mulVec_diagonal]

/-- A relative change `ε` in `b` gives a relative change `1` in `x`. -/
theorem diagEps_relative {ε : ℝ} (hε : 0 < ε) :
    vnorm (![0, ε] : Fin 2 → ℝ) / vnorm (![1, 0] : Fin 2 → ℝ) = ε ∧
      vnorm (![0, 1] : Fin 2 → ℝ) / vnorm (![1, 0] : Fin 2 → ℝ) = 1 := by
  have e1 : (![1, 0] : Fin 2 → ℝ) = Pi.single 0 1 := by ext i; fin_cases i <;> simp
  have e2 : (![0, ε] : Fin 2 → ℝ) = Pi.single 1 ε := by ext i; fin_cases i <;> simp
  have e3 : (![0, 1] : Fin 2 → ℝ) = Pi.single 1 1 := by ext i; fin_cases i <;> simp
  rw [e1, e2, e3, vnorm_single_val, vnorm_single_val, vnorm_single_val, abs_of_pos hε]
  simp

/-- **`κ₂(diag(1, ε)) = 1/ε`** for `0 < ε ≤ 1`. -/
theorem diagEps_cond {ε : ℝ} (hε : 0 < ε) (hε1 : ε ≤ 1) : cond2 (diagEps ε) = 1 / ε := by
  have hinv : (diagEps ε)⁻¹ = Matrix.diagonal ![1, 1 / ε] := by
    apply Matrix.inv_eq_left_inv
    rw [diagEps, Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
    congr 1
    ext i
    fin_cases i <;> simp [hε.ne']
  have h1 : 1 ≤ 1 / ε := by rw [le_div_iff₀ hε]; linarith
  unfold cond2
  rw [hinv, diagEps, norm2_diagonal_two, norm2_diagonal_two, abs_one, abs_of_pos hε,
    abs_of_pos (by positivity : (0 : ℝ) < 1 / ε), max_eq_left hε1, max_eq_right h1, one_mul]

end LU
