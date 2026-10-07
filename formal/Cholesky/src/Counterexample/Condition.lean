import src.Counterexample.Breakdown

/-!
# The condition number of `A_{9k}`

With `x = k^{3/2} u` (so `δ = 12x`, `λ_min(S) = 9x`, `n^{3/2} u = 27x`):

* `‖T‖₂² ≤ 9/2 + 27x` (`norm2_T_sq_le`), from `Tᵀ T = G + (9/8) η J`;
* `‖A‖₂ ≤ 17/2 + 39x` (`norm2_A_le`) and `‖A⁻¹‖₂ ≤ 1 + 17/(72x)` (`norm2_Ainv_le`);
* hence `κ₂(A) ≤ conditionUpper x < 57/(n^{3/2} u)` (`cond_A_lt`);
* `‖A‖₂ ≥ 17/2` and `‖A⁻¹‖₂ ≥ 17/(8δ)`, hence `κ₂(A) ≥ (2601/64)/(n^{3/2} u)` and
  `κ₂(A) ≥ 1156/3 > 385` (`cond_A_ge`).
-/
namespace Cholesky
namespace CE
namespace Par
open Matrix

variable (P : Par)

/-! ## Quadratic form of `G` and the norm of `T` -/

lemma vnorm_Q (w : Fin P.k → ℝ) : vnorm (P.Q *ᵥ w) = vnorm w :=
  vnorm_orthogonal P.Q_transpose_mul w

lemma vnorm_QT (w : Fin P.k → ℝ) : vnorm (P.Qᵀ *ᵥ w) = vnorm w := by
  rw [P.Q_transpose]; exact P.vnorm_Q w

/-- `yᵀ G y = (9/4)(‖y₁‖² + ‖y₂‖² + 2 y₁ᵀ Q y₂)`. -/
lemma G_quad (y : Fin (P.k + P.k) → ℝ) :
    y ⬝ᵥ (P.G *ᵥ y) = 9 / 4 * (vnorm (vL y) ^ 2 + vnorm (vR y) ^ 2 +
      2 * (vL y ⬝ᵥ (P.Q *ᵥ vR y))) := by
  conv_lhs => rw [← append_vL_vR y]
  rw [G, Matrix.smul_mulVec, blk_mulVec_append, dotProduct_smul, append_dotProduct,
    vnorm_sq, vnorm_sq]
  simp only [Matrix.one_mulVec, dotProduct_add, smul_eq_mul]
  rw [Matrix.dotProduct_mulVec (vR y) P.Qᵀ, Matrix.vecMul_transpose, dotProduct_comm (P.Q *ᵥ vR y)]
  ring

lemma G_quad_le (y : Fin (P.k + P.k) → ℝ) : y ⬝ᵥ (P.G *ᵥ y) ≤ 9 / 2 * vnorm y ^ 2 := by
  rw [P.G_quad, vnorm_sq_split y]
  have h := dot_le (vL y) (P.Q *ᵥ vR y)
  rw [P.vnorm_Q] at h
  nlinarith [sq_nonneg (vnorm (vL y) - vnorm (vR y))]

lemma T_quad (y : Fin (P.k + P.k) → ℝ) :
    vnorm (P.T *ᵥ y) ^ 2 = y ⬝ᵥ (P.G *ᵥ y) + 9 / 8 * P.η * (∑ i, y i) ^ 2 := by
  rw [vnorm_sq, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, ← dotProduct_comm,
    Matrix.mulVec_mulVec, P.T_gram, Matrix.add_mulVec, dotProduct_add, Matrix.smul_mulVec,
    dotProduct_smul, J2, dot_ones_mulVec, smul_eq_mul]
  ring

/-- `‖T y‖² ≤ (9/2 + 27x) ‖y‖²`. -/
theorem T_sq_le (y : Fin (P.k + P.k) → ℝ) :
    vnorm (P.T *ᵥ y) ^ 2 ≤ (9 / 2 + 27 * P.x) * vnorm y ^ 2 := by
  rw [P.T_quad]
  have h1 := P.G_quad_le y
  have h2 := sum_sq_le_card y
  have hη := P.η_pos
  have hkη : (9 / 8 * P.η) * ((P.k + P.k : ℕ) : ℝ) = 27 * P.x := by
    push_cast
    have := P.δ_eq_x; rw [δ] at this
    linarith
  nlinarith [mul_le_mul_of_nonneg_left h2 (by positivity : (0 : ℝ) ≤ 9 / 8 * P.η)]

/-- `‖T‖₂² ≤ 9/2 + 27x`. -/
theorem norm2_T_sq_le : norm2 P.T ^ 2 ≤ 9 / 2 + 27 * P.x := by
  have hx := P.x_pos
  have hτ : 0 ≤ Real.sqrt (9 / 2 + 27 * P.x) := Real.sqrt_nonneg _
  have hle : norm2 P.T ≤ Real.sqrt (9 / 2 + 27 * P.x) := by
    apply norm2_le_of hτ
    intro y
    have := P.T_sq_le y
    have hsq : vnorm (P.T *ᵥ y) ^ 2 ≤ (Real.sqrt (9 / 2 + 27 * P.x) * vnorm y) ^ 2 := by
      rw [mul_pow, Real.sq_sqrt (by positivity)]; exact this
    exact le_of_pow_le_pow_left₀ two_ne_zero (mul_nonneg hτ (vnorm_nonneg y)) hsq
  calc norm2 P.T ^ 2 ≤ Real.sqrt (9 / 2 + 27 * P.x) ^ 2 :=
        pow_le_pow_left₀ (norm2_nonneg _) hle 2
    _ = 9 / 2 + 27 * P.x := Real.sq_sqrt (by positivity)

/-! ## Upper bounds -/

lemma kη_eq : (P.k : ℝ) * P.η = 12 * P.x := by have := P.δ_eq_x; rw [δ] at this; exact this

/-- `‖A‖₂ ≤ 17/2 + 39x`. -/
theorem norm2_A_le : norm2 P.A ≤ 17 / 2 + 39 * P.x := by
  have hx := P.x_pos
  apply norm2_le_of_quad_psd P.A_symm (by positivity)
  · intro v
    have := ((posDef_iff_real P.A).mp P.A_posDef).2 v
    by_cases hv : v = 0
    · rw [hv]; simp
    · exact (this hv).le
  · intro v
    rw [← append_vL_vR v, P.A_quad, vnorm_append_sq]
    set a := vnorm (vL v)
    set b := vnorm (vR v)
    have hT := P.T_sq_le (vR v)
    have hS := (P.S_quad_bounds (vR v)).2
    rw [P.kη_eq] at hS
    have hlin : vnorm ((2 : ℝ) • vL v + P.T *ᵥ vR v) ≤ 2 * a + vnorm (P.T *ᵥ vR v) := by
      calc vnorm ((2 : ℝ) • vL v + P.T *ᵥ vR v) ≤ vnorm ((2 : ℝ) • vL v) + vnorm (P.T *ᵥ vR v) :=
            vnorm_add_le _ _
        _ = 2 * a + vnorm (P.T *ᵥ vR v) := by rw [vnorm_smul, abs_of_pos (by norm_num)]
    have hc : ((2 : ℝ) • vL v + P.T *ᵥ vR v) ⬝ᵥ ((2 : ℝ) • vL v + P.T *ᵥ vR v) ≤
        (2 * a + vnorm (P.T *ᵥ vR v)) ^ 2 := by
      rw [← vnorm_sq]; exact pow_le_pow_left₀ (vnorm_nonneg _) hlin 2
    have ha := vnorm_nonneg (vL v)
    have htn := vnorm_nonneg (P.T *ᵥ vR v)
    -- `(2a + t)² ≤ (4 + c)(a² + b²)` when `t² ≤ c b²`, `c = 9/2 + 27x`
    set c := 9 / 2 + 27 * P.x with hcdef
    set t := vnorm (P.T *ᵥ vR v)
    have hc : 0 < c := by positivity
    have hCS : (2 * a + t) ^ 2 ≤ (4 + c) * (a ^ 2 + b ^ 2) := by
      have h1 : 0 ≤ c * ((4 + c) * (a ^ 2 + b ^ 2) - (2 * a + t) ^ 2) := by
        nlinarith [sq_nonneg (c * a - 2 * t), mul_le_mul_of_nonneg_left hT hc.le]
      have h2 : 0 ≤ (4 + c) * (a ^ 2 + b ^ 2) - (2 * a + t) ^ 2 :=
        nonneg_of_mul_nonneg_right (by linarith) hc
      linarith
    have hb := sq_nonneg b
    nlinarith

lemma A_mulVec_append (x1 : Fin (6 * P.k + P.k) → ℝ) (y : Fin (P.k + P.k) → ℝ) :
    P.A *ᵥ Fin.append x1 y = Fin.append ((4 : ℝ) • x1 + (2 : ℝ) • (P.T *ᵥ y))
      ((2 : ℝ) • (P.Tᵀ *ᵥ x1) + P.C *ᵥ y) := by
  rw [A, blk_mulVec_append]
  simp only [Matrix.smul_mulVec, Matrix.one_mulVec]

lemma S_lower (w : Fin (P.k + P.k) → ℝ) : 9 * P.x * vnorm w ^ 2 ≤ w ⬝ᵥ (P.S *ᵥ w) := by
  have := (P.S_quad_bounds w).1
  have h := P.kη_eq
  nlinarith [sq_nonneg (vnorm w)]

lemma norm2_T_le : norm2 P.T ≤ Real.sqrt (9 / 2 + 27 * P.x) := by
  have := Real.abs_le_sqrt P.norm2_T_sq_le
  rwa [abs_of_nonneg (norm2_nonneg _)] at this

lemma cs2 (a b t : ℝ) : (b + t / 2 * a) ^ 2 ≤ (1 + t ^ 2 / 4) * (a ^ 2 + b ^ 2) := by
  nlinarith [sq_nonneg (t / 2 * b - a)]

lemma aux_sq_ineq (t w τ : ℝ) (h : t * t ≤ τ * w * (τ * w)) :
    (1 / 2 * t) ^ 2 + w ^ 2 ≤ (1 + τ ^ 2 / 4) * w ^ 2 := by nlinarith

lemma le_of_sq_le' {a c : ℝ} (hc : 0 ≤ c) (h : a ^ 2 ≤ c ^ 2) : a ≤ c :=
  le_of_pow_le_pow_left₀ two_ne_zero hc h

/-- `‖A⁻¹‖₂ ≤ 1 + 17/(72x)`, from the block solution of `A v = z`. -/
theorem norm2_Ainv_le : norm2 P.A⁻¹ ≤ 1 + 17 / (72 * P.x) := by
  have hx := P.x_pos
  have hdet : IsUnit P.A.det := (Matrix.isUnit_iff_isUnit_det _).mp P.A_posDef.isUnit
  -- `τ ≥ ‖T‖₂` and `ρ² = 1 + τ²/4`
  obtain ⟨τ, hτ0, hτ2, hT⟩ : ∃ τ : ℝ, 0 ≤ τ ∧ τ ^ 2 = 9 / 2 + 27 * P.x ∧ norm2 P.T ≤ τ :=
    ⟨_, Real.sqrt_nonneg _, Real.sq_sqrt (by positivity), P.norm2_T_le⟩
  obtain ⟨ρ, hρ0, hρ2⟩ : ∃ ρ : ℝ, 0 ≤ ρ ∧ ρ ^ 2 = 1 + τ ^ 2 / 4 :=
    ⟨_, Real.sqrt_nonneg _, Real.sq_sqrt (by positivity)⟩
  have hTv : ∀ y, vnorm (P.T *ᵥ y) ≤ τ * vnorm y := fun y =>
    (vnorm_mulVec_le P.T y).trans (mul_le_mul_of_nonneg_right hT (vnorm_nonneg y))
  have hTTv : ∀ w, vnorm (P.Tᵀ *ᵥ w) ≤ τ * vnorm w := fun w =>
    (vnorm_mulVec_le P.Tᵀ w).trans (mul_le_mul_of_nonneg_right
      (by rw [norm2_transpose]; exact hT) (vnorm_nonneg w))
  have hconst : 1 / 4 + ρ ^ 2 / (9 * P.x) = 1 + 17 / (72 * P.x) := by
    rw [hρ2, hτ2]; field_simp; ring
  apply norm2_le_of (by positivity)
  intro z
  set v := P.A⁻¹ *ᵥ z
  have hAv : P.A *ᵥ v = z := by
    simp only [v, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hdet, Matrix.one_mulVec]
  have hblock := P.A_mulVec_append (vL v) (vR v)
  rw [append_vL_vR, hAv] at hblock
  have h1 : (4 : ℝ) • vL v + (2 : ℝ) • (P.T *ᵥ vR v) = vL z := by
    have := congrArg vL hblock; rw [vL_append] at this; exact this.symm
  have h2 : (2 : ℝ) • (P.Tᵀ *ᵥ vL v) + P.C *ᵥ vR v = vR z := by
    have := congrArg vR hblock; rw [vR_append] at this; exact this.symm
  -- `S y = z₂ - Tᵀ z₁ / 2` and `x₁ = z₁/4 - T y/2`
  have hSy : P.S *ᵥ vR v = vR z - (1 / 2 : ℝ) • (P.Tᵀ *ᵥ vL z) := by
    rw [← h1, ← h2, P.C_eq, Matrix.add_mulVec, Matrix.mulVec_add, Matrix.mulVec_smul,
      Matrix.mulVec_smul, ← Matrix.mulVec_mulVec]
    ext i; simp only [Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]; ring
  have hx1 : vL v = (1 / 4 : ℝ) • vL z - (1 / 2 : ℝ) • (P.T *ᵥ vR v) := by
    ext i
    have := congrFun h1 i
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.sub_apply] at this ⊢
    linarith
  set a := vnorm (vL z)
  set b := vnorm (vR z)
  set w := vnorm (vR v)
  have ha : 0 ≤ a := vnorm_nonneg (vL z)
  have hb : 0 ≤ b := vnorm_nonneg (vR z)
  have hw := vnorm_nonneg (vR v)
  have hz : vnorm z ^ 2 = a ^ 2 + b ^ 2 := vnorm_sq_split z
  -- `‖y‖ ≤ (‖z₂‖ + τ ‖z₁‖ / 2) / (9x) ≤ ρ ‖z‖ / (9x)`
  have hy1 : w ≤ vnorm (P.S *ᵥ vR v) / (9 * P.x) :=
    vnorm_le_of_quad_lower (by positivity) P.S_lower (vR v)
  have hSy' : vnorm (P.S *ᵥ vR v) ≤ b + τ / 2 * a := by
    rw [hSy]
    calc vnorm (vR z - (1 / 2 : ℝ) • (P.Tᵀ *ᵥ vL z))
        ≤ vnorm (vR z) + vnorm ((1 / 2 : ℝ) • (P.Tᵀ *ᵥ vL z)) := vnorm_sub_le _ _
      _ = b + 1 / 2 * vnorm (P.Tᵀ *ᵥ vL z) := by rw [vnorm_smul, abs_of_pos (by norm_num)]
      _ ≤ b + τ / 2 * a := by have := hTTv (vL z); linarith
  have hCS : b + τ / 2 * a ≤ ρ * vnorm z :=
    le_of_sq_le' (mul_nonneg hρ0 (vnorm_nonneg z)) (by
      rw [mul_pow, hρ2, hz]; exact cs2 a b τ)
  have hy : w ≤ ρ * vnorm z / (9 * P.x) :=
    le_trans hy1 (div_le_div_of_nonneg_right (le_trans hSy' hCS) (by positivity))
  -- `v = (z₁/4, 0) + (-T y/2, y)`
  have hsplit : v = Fin.append ((1 / 4 : ℝ) • vL z) 0 +
      Fin.append (-((1 / 2 : ℝ) • (P.T *ᵥ vR v))) (vR v) := by
    rw [append_add, ← sub_eq_add_neg, ← hx1, zero_add, append_vL_vR]
  have hp1 : vnorm (Fin.append ((1 / 4 : ℝ) • vL z) (0 : Fin (P.k + P.k) → ℝ)) ≤ a / 4 :=
    le_of_sq_le' (div_nonneg ha (by norm_num)) (le_of_eq (by
      rw [vnorm_append_sq, vnorm_zero, vnorm_smul, abs_of_pos (by norm_num)]; ring))
  have hp2 : vnorm (Fin.append (-((1 / 2 : ℝ) • (P.T *ᵥ vR v))) (vR v)) ≤ ρ * w := by
    apply le_of_sq_le' (mul_nonneg hρ0 hw)
    have e : (ρ * vnorm (vR v)) ^ 2 = (1 + τ ^ 2 / 4) * vnorm (vR v) ^ 2 := by rw [mul_pow, hρ2]
    rw [vnorm_append_sq, vnorm_neg, vnorm_smul, abs_of_pos (by norm_num), e]
    have h3 := hTv (vR v)
    have ht := vnorm_nonneg (P.T *ᵥ vR v)
    exact aux_sq_ineq _ _ _ (mul_self_le_mul_self ht h3)
  have key : vnorm v ≤ a / 4 + ρ * w := by
    rw [hsplit]; exact le_trans (vnorm_add_le _ _) (add_le_add hp1 hp2)
  have haz : a ≤ vnorm z :=
    le_of_sq_le' (vnorm_nonneg z) (by rw [hz]; exact le_add_of_nonneg_right (sq_nonneg b))
  have h5 := mul_le_mul_of_nonneg_left hy hρ0
  have h6 : a / 4 ≤ vnorm z / 4 := div_le_div_of_nonneg_right haz (by norm_num)
  calc vnorm v ≤ vnorm z / 4 + ρ * (ρ * vnorm z / (9 * P.x)) := by linarith
    _ = (1 / 4 + ρ ^ 2 / (9 * P.x)) * vnorm z := by ring
    _ = (1 + 17 / (72 * P.x)) * vnorm z := by rw [hconst]

/-- **The condition number of `A_{9k}`:** `κ₂(A) ≤ (17/2 + 39x)(1 + 17/(72x))`. -/
theorem cond_A_le : cond2 P.A ≤ conditionUpper P.x := by
  unfold cond2 conditionUpper
  have hx := P.x_pos
  exact mul_le_mul P.norm2_A_le P.norm2_Ainv_le (norm2_nonneg _) (by positivity)

/-- `n^{3/2} u = 27 x` for `n = 9k` (with `n^{3/2} = n √n`). -/
lemma n32u : ((9 * P.k : ℕ) : ℝ) * Real.sqrt ((9 * P.k : ℕ) : ℝ) * P.u = 27 * P.x := by
  have hK := P.K_pos
  have hs : Real.sqrt ((9 * P.k : ℕ) : ℝ) = 3 * P.K := by
    rw [show ((9 * P.k : ℕ) : ℝ) = (3 * P.K) ^ 2 by push_cast; rw [P.k_real]; ring,
      Real.sqrt_sq (by positivity)]
  rw [hs, x]; push_cast; rw [P.k_real]; ring

/-- `κ₂(A_{9k}) < 57 / (n^{3/2} u)`. -/
theorem cond_A_lt :
    cond2 P.A < 57 / (((9 * P.k : ℕ) : ℝ) * Real.sqrt ((9 * P.k : ℕ) : ℝ) * P.u) := by
  rw [P.n32u]
  exact lt_of_le_of_lt P.cond_A_le (conditionUpper_lt_57 P.x P.x_pos P.x_le)

/-! ## Lower bounds -/

lemma aux_lower1 (a : ℝ) (ha : 9 ≤ a) :
    17 / 2 * ((4 / a) ^ 2 * a + 2) ≤ (2 * (4 / a) + 1) ^ 2 * a := by
  have ha0 : 0 < a := by linarith
  rw [← sub_nonneg]
  have : (2 * (4 / a) + 1) ^ 2 * a - 17 / 2 * ((4 / a) ^ 2 * a + 2) = (a + 8) * (a - 9) / a := by
    field_simp; ring
  rw [this]
  exact div_nonneg (mul_nonneg (by linarith) (by linarith)) ha0.le

def i0 : Fin P.k := ⟨0, by have := P.k_ge; omega⟩

lemma vnorm_single_k (i : Fin P.k) : vnorm (Pi.single i (1 : ℝ)) = 1 := vnorm_single i

/-- `‖A‖₂ ≥ 17/2`. -/
theorem norm2_A_ge : 17 / 2 ≤ norm2 P.A := by
  obtain ⟨e0, he0def⟩ : ∃ e0 : Fin P.k → ℝ, e0 = Pi.single P.i0 1 := ⟨_, rfl⟩
  obtain ⟨y0, hy0def⟩ : ∃ y0 : Fin (P.k + P.k) → ℝ, y0 = Fin.append e0 (P.Qᵀ *ᵥ e0) := ⟨_, rfl⟩
  have he0 : vnorm e0 = 1 := by rw [he0def]; exact vnorm_single _
  have hQe0 : vnorm (P.Qᵀ *ᵥ e0) = 1 := by rw [P.vnorm_QT, he0]
  have hy0 : vnorm y0 ^ 2 = 2 := by rw [hy0def, vnorm_append_sq, he0, hQe0]; norm_num
  have hG : y0 ⬝ᵥ (P.G *ᵥ y0) = 9 := by
    rw [P.G_quad, hy0def, vL_append, vR_append, he0, hQe0, Matrix.mulVec_mulVec,
      P.Q_mul_transpose, Matrix.one_mulVec, ← vnorm_sq, he0]
    norm_num
  obtain ⟨a, hadef⟩ : ∃ a : ℝ, a = vnorm (P.T *ᵥ y0) ^ 2 := ⟨_, rfl⟩
  have ha : 9 ≤ a := by
    have hT := P.T_quad y0
    have hη := P.η_pos
    have hpos : 0 ≤ 9 / 8 * P.η * (∑ i, y0 i) ^ 2 :=
      mul_nonneg (by linarith) (sq_nonneg _)
    rw [hadef]; linarith
  have ha0 : 0 < a := by linarith
  obtain ⟨c, hcdef⟩ : ∃ c : ℝ, c = 4 / a := ⟨_, rfl⟩
  obtain ⟨v, hvdef⟩ : ∃ v : Fin ((6 * P.k + P.k) + (P.k + P.k)) → ℝ,
    v = Fin.append (c • (P.T *ᵥ y0)) y0 := ⟨_, rfl⟩
  have hv2 : vnorm v ^ 2 = c ^ 2 * a + 2 := by
    rw [hvdef, vnorm_append_sq, vnorm_smul, mul_pow, sq_abs, hy0, hadef]
  have hquad : (2 * c + 1) ^ 2 * a ≤ v ⬝ᵥ (P.A *ᵥ v) := by
    have hcomb : (2 : ℝ) • (c • (P.T *ᵥ y0)) + P.T *ᵥ y0 = (2 * c + 1) • (P.T *ᵥ y0) := by
      rw [smul_smul, add_smul, one_smul]
    rw [hvdef, P.A_quad, hcomb]
    have hS := (P.S_quad_bounds y0).1
    have hk : (0 : ℝ) ≤ P.k := Nat.cast_nonneg _
    have hS0 : 0 ≤ y0 ⬝ᵥ (P.S *ᵥ y0) :=
      le_trans (mul_nonneg (mul_nonneg (mul_nonneg (by norm_num) hk) P.η_pos.le)
        (sq_nonneg _)) hS
    have : ((2 * c + 1) • (P.T *ᵥ y0)) ⬝ᵥ ((2 * c + 1) • (P.T *ᵥ y0)) = (2 * c + 1) ^ 2 * a := by
      rw [← vnorm_sq, vnorm_smul, mul_pow, sq_abs, hadef]
    rw [this]; linarith
  have hle := quad_le_norm2 P.A v
  have hv0 : 0 < vnorm v ^ 2 := by
    rw [hv2]; exact add_pos_of_nonneg_of_pos (mul_nonneg (sq_nonneg c) ha0.le) (by norm_num)
  have h1 := aux_lower1 a ha
  rw [← hcdef, ← hv2] at h1
  by_contra hlt
  push_neg at hlt
  have : norm2 P.A * vnorm v ^ 2 < 17 / 2 * vnorm v ^ 2 := mul_lt_mul_of_pos_right hlt hv0
  linarith

def i1 : Fin P.k := ⟨1, by have := P.k_ge; omega⟩
def i2 : Fin P.k := ⟨2, by have := P.k_ge; omega⟩

lemma i1_ne_i2 : P.i1 ≠ P.i2 := by intro h; have := congrArg Fin.val h; simp [i1, i2] at this

lemma i1_ne_0 : P.i1 ≠ had0 (2 * P.q) := by
  intro h; have := congrArg Fin.val h; simp [i1, had0] at this
lemma i2_ne_0 : P.i2 ≠ had0 (2 * P.q) := by
  intro h; have := congrArg Fin.val h; simp [i2, had0] at this

/-- The vector `w₀ = e₁ - e₂`. -/
noncomputable def w0 : Fin P.k → ℝ := Pi.single P.i1 1 - Pi.single P.i2 1

lemma w0_ne : P.w0 ≠ 0 := by
  intro h
  have := congrFun h P.i1
  simp [w0, P.i1_ne_i2] at this

lemma sum_w0 : ∑ i, P.w0 i = 0 := by
  simp [w0, Finset.sum_sub_distrib]

lemma sum_QT_w0 : ∑ j, (P.Qᵀ *ᵥ P.w0) j = 0 := by
  have h1 := hadQ_row_sum P.q P.i1 P.i1_ne_0
  have h2 := hadQ_row_sum P.q P.i2 P.i2_ne_0
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, w0, Pi.sub_apply,
    Pi.single_apply, mul_sub, mul_ite, mul_one, mul_zero, Finset.sum_sub_distrib,
    Finset.sum_ite_eq', Finset.mem_univ, if_true]
  have h1' : ∑ j : Fin P.k, P.Q P.i1 j = 0 := h1
  have h2' : ∑ j : Fin P.k, P.Q P.i2 j = 0 := h2
  rw [h1', h2', sub_self]

lemma δ_pos : 0 < P.δ := by rw [P.δ_eq_x]; have := P.x_pos; linarith

lemma aux_lower2 (δ Y : ℝ) (hδ : 0 < δ) (hY : 0 ≤ Y) :
    17 / (8 * δ) * ((1 / 2) ^ 2 * (9 / 2 * Y) + Y) ≤
      1 / 2 * (1 / 8 + 17 / (16 * δ)) * (9 / 2 * Y) + 17 / (8 * δ) * Y := by
  rw [← sub_nonneg]
  have : 1 / 2 * (1 / 8 + 17 / (16 * δ)) * (9 / 2 * Y) + 17 / (8 * δ) * Y -
      17 / (8 * δ) * ((1 / 2) ^ 2 * (9 / 2 * Y) + Y) = 9 / 32 * Y := by
    field_simp; ring
  rw [this]; positivity

/-- `‖A⁻¹‖₂ ≥ 17/(8δ)`. -/
theorem norm2_Ainv_ge : 17 / (8 * P.δ) ≤ norm2 P.A⁻¹ := by
  have hδ := P.δ_pos
  have hdet : IsUnit P.A.det := (Matrix.isUnit_iff_isUnit_det _).mp P.A_posDef.isUnit
  obtain ⟨y, hydef⟩ : ∃ y : Fin (P.k + P.k) → ℝ, y = Fin.append P.w0 (P.Qᵀ *ᵥ P.w0) :=
    ⟨_, rfl⟩
  have hsum : ∑ i, y i = 0 := by
    rw [hydef, Fin.sum_univ_add]
    simp only [Fin.append_left, Fin.append_right, P.sum_w0, P.sum_QT_w0, add_zero]
  have hGy : P.G *ᵥ y = (9 / 2 : ℝ) • y := by
    rw [hydef, G, Matrix.smul_mulVec, blk_mulVec_append, Matrix.one_mulVec, Matrix.one_mulVec,
      Matrix.mulVec_mulVec, P.Q_mul_transpose, Matrix.one_mulVec, append_smul, append_smul]
    congr 1 <;> ext i <;> simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul] <;> ring
  have hJy : P.J2 *ᵥ y = 0 := by rw [J2, ones_mulVec, hsum]; rfl
  have hTTy : P.Tᵀ *ᵥ (P.T *ᵥ y) = (9 / 2 : ℝ) • y := by
    rw [Matrix.mulVec_mulVec, P.T_gram, Matrix.add_mulVec, hGy, Matrix.smul_mulVec, hJy,
      smul_zero, add_zero]
  have hSy : P.S *ᵥ y = P.δ • y := P.S_perp y hsum
  have hCy : P.C *ᵥ y = (9 / 2 + P.δ) • y := by
    rw [P.C_eq, Matrix.add_mulVec, ← Matrix.mulVec_mulVec, hTTy, hSy, add_smul]
  obtain ⟨Y, hYdef⟩ : ∃ Y : ℝ, Y = vnorm y ^ 2 := ⟨_, rfl⟩
  have hTy2 : vnorm (P.T *ᵥ y) ^ 2 = 9 / 2 * Y := by
    rw [vnorm_sq, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, dotProduct_comm,
      hTTy, dotProduct_smul, ← vnorm_sq, hYdef, smul_eq_mul]
  have hy0 : y ≠ 0 := by
    rw [hydef]; intro h; exact P.w0_ne ((append_eq_zero_iff _ _).mp h).1
  have hYpos : 0 < Y := by rw [hYdef]; exact pow_pos (vnorm_pos hy0) 2
  set α : ℝ := 1 / 8 + 17 / (16 * P.δ) with hα
  set β : ℝ := 17 / (8 * P.δ) with hβ
  obtain ⟨z, hzdef⟩ : ∃ z : Fin ((6 * P.k + P.k) + (P.k + P.k)) → ℝ,
    z = Fin.append (-((1 / 2 : ℝ) • (P.T *ᵥ y))) y := ⟨_, rfl⟩
  obtain ⟨w', hwdef⟩ : ∃ w' : Fin ((6 * P.k + P.k) + (P.k + P.k)) → ℝ,
    w' = Fin.append (-(α • (P.T *ᵥ y))) (β • y) := ⟨_, rfl⟩
  have hAw : P.A *ᵥ w' = z := by
    rw [hwdef, hzdef, P.A_mulVec_append]
    simp only [Matrix.mulVec_neg, Matrix.mulVec_smul]
    rw [hTTy, hCy]
    have hδ0 : P.δ ≠ 0 := hδ.ne'
    have c1 : (4 : ℝ) • -(α • (P.T *ᵥ y)) + (2 : ℝ) • (β • (P.T *ᵥ y)) =
        -((1 / 2 : ℝ) • (P.T *ᵥ y)) := by
      rw [smul_neg, smul_smul, smul_smul, ← neg_smul, ← add_smul, ← neg_smul]
      congr 1; rw [hα, hβ]; field_simp; ring
    have c2 : (2 : ℝ) • -(α • ((9 / 2 : ℝ) • y)) + β • ((9 / 2 + P.δ) • y) = y := by
      rw [smul_neg, smul_smul, smul_smul, smul_smul, ← neg_smul, ← add_smul]
      conv_rhs => rw [← one_smul ℝ y]
      congr 1; rw [hα, hβ]; field_simp; ring
    rw [c1, c2]
  have hinv : P.A⁻¹ *ᵥ z = w' := by
    rw [← hAw, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
  have hzw : z ⬝ᵥ w' = 1 / 2 * α * (9 / 2 * Y) + β * Y := by
    rw [hzdef, hwdef, append_dotProduct, neg_dotProduct, dotProduct_neg, neg_neg, smul_dotProduct,
      dotProduct_smul, dotProduct_smul, ← vnorm_sq, ← vnorm_sq, hTy2, ← hYdef]
    simp only [smul_eq_mul]; ring
  have hz2 : vnorm z ^ 2 = (1 / 2) ^ 2 * (9 / 2 * Y) + Y := by
    rw [hzdef, vnorm_append_sq, vnorm_neg, vnorm_smul, abs_of_pos (by norm_num), mul_pow, hTy2,
      ← hYdef]
  have hle := quad_le_norm2 P.A⁻¹ z
  rw [hinv, hzw, hz2] at hle
  have hmain := aux_lower2 P.δ Y hδ hYpos.le
  rw [← hα, ← hβ] at hmain
  have hz2pos : 0 < (1 / 2) ^ 2 * (9 / 2 * Y) + Y := by positivity
  by_contra hlt
  push_neg at hlt
  have := mul_lt_mul_of_pos_right hlt hz2pos
  linarith

/-- **Lower bound.** `κ₂(A_{9k}) ≥ 289/(192x)`. -/
theorem cond_A_ge : 289 / (192 * P.x) ≤ cond2 P.A := by
  have hx := P.x_pos
  have h := mul_le_mul P.norm2_A_ge P.norm2_Ainv_ge
    (div_nonneg (by norm_num) (by linarith [P.δ_pos])) (norm2_nonneg _)
  have heq : 17 / 2 * (17 / (8 * P.δ)) = 289 / (192 * P.x) := by
    rw [P.δ_eq_x]; field_simp; ring
  rw [heq] at h
  exact h

/-- `κ₂(A_{9k}) ≥ (2601/64) / (n^{3/2} u)`. -/
theorem cond_A_ge' :
    2601 / 64 / (((9 * P.k : ℕ) : ℝ) * Real.sqrt ((9 * P.k : ℕ) : ℝ) * P.u) ≤ cond2 P.A := by
  rw [P.n32u]
  have hx := P.x_pos
  have : 2601 / 64 / (27 * P.x) = 289 / (192 * P.x) := by field_simp; ring
  rw [this]; exact P.cond_A_ge

/-- The family cannot have condition numbers near one: `κ₂(A_{9k}) ≥ 1156/3 > 385`. -/
theorem cond_A_ge_385 : (1156 : ℝ) / 3 ≤ cond2 P.A ∧ (385 : ℝ) < cond2 P.A := by
  have h1 := family_lower_bound P.x P.x_pos P.x_le
  have h2 := P.cond_A_ge
  constructor
  · linarith
  · linarith

end Par
end CE
end Cholesky