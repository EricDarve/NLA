import src.Counterexample.Sharpness
import src.SPD

/-!
# Small examples and numerical facts from the notes

* `nonmonotone_p3`, `nonmonotone_p4`: the SPD matrix `(1/32) [[4, 5], [5, 7]]`; Cholesky
  succeeds with `p = 3` (second pivot `1/32`) and fails with `p = 4` (second pivot `0`).
* `order_matters`: `fl(fl(a - b) - c) ≠ fl(a - fl(b + c))` in an example.
* The numbers in the table of numerical examples, the memory sizes, and the
  expansion `289/(144x) + 425/24 + 39x` of the condition-number bound.
* The intermediate inequalities of the `1.1` corollary.
* Remarks on scaling: `H = D⁻¹ A D⁻¹` has unit diagonal; `κ₂(P A Pᵀ) = κ₂(A)`.
-/
namespace Cholesky
open Matrix

/-! ## Failure is not monotone in precision -/

/-- `A = (1/32) [[4, 5], [5, 7]]`. -/
noncomputable def nmA : Matrix (Fin 2) (Fin 2) ℝ := !![4 / 32, 5 / 32; 5 / 32, 7 / 32]

theorem nmA_posDef : nmA.PosDef ∧ nmA.det = 3 / 1024 := by
  refine ⟨?_, by simp [nmA, Matrix.det_fin_two]; norm_num⟩
  rw [posDef_iff_real]
  refine ⟨by ext i j; fin_cases i <;> fin_cases j <;> rfl, fun x hx => ?_⟩
  simp only [nmA, dotProduct, Matrix.mulVec, Fin.sum_univ_two, Matrix.of_apply,
    Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.empty_val', Matrix.cons_val_fin_one]
  have : x 0 ≠ 0 ∨ x 1 ≠ 0 := by
    by_contra h; push_neg at h; apply hx; ext i; fin_cases i <;> simp [h.1, h.2]
  rcases this with h | h
  · nlinarith [sq_pos_of_ne_zero h, sq_nonneg (4 * x 0 + 5 * x 1), sq_nonneg (x 1)]
  · nlinarith [sq_pos_of_ne_zero h, sq_nonneg (4 * x 0 + 5 * x 1), sq_nonneg (x 1)]

lemma isFloat_mk {p : ℕ} (m : ℤ) (e : ℤ) (hm : |m| < 2 ^ p) : IsFloat p (m * (2 : ℝ) ^ e) :=
  ⟨m, e, hm, rfl⟩

/-- The pivots of a `2 × 2` run. -/
lemma pivot_two (o : Ops) (A : Matrix (Fin 2) (Fin 2) ℝ) :
    pivot o A 0 = A 0 0 ∧
      pivot o A 1 = o.upd (A 1 1) (o.div (A 1 0) (o.sqrt (A 0 0)))
        (o.div (A 1 0) (o.sqrt (A 0 0))) := by
  refine ⟨rfl, ?_⟩
  have := pivot_succ o A 0
  simp only [Fin.succ_zero_eq_one] at this
  rw [this, pivot_zero, trail_apply _ _ _ _ le_rfl]
  rfl

lemma sqrt_eighth_bounds' : 22 / 64 < Real.sqrt (1 / 8) ∧ Real.sqrt (1 / 8) < 23 / 64 := by
  constructor
  · rw [show (22 / 64 : ℝ) = Real.sqrt ((22 / 64) ^ 2) by rw [Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  · rw [show (23 / 64 : ℝ) = Real.sqrt ((23 / 64) ^ 2) by rw [Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)

lemma sqrt_eighth_bounds : 22 / 64 < (Real.sqrt 8)⁻¹ ∧ (Real.sqrt 8)⁻¹ < 23 / 64 := by
  have e : Real.sqrt (1 / 8) = (Real.sqrt 8)⁻¹ := by rw [one_div, Real.sqrt_inv]
  rw [← e]; exact sqrt_eighth_bounds'

/-- With `p = 3` the computed factor is `l₁₁ = 3/8`, `l₂₁ = 7/16`, and the second pivot is
`1/32 > 0`: the factorization succeeds. -/
theorem nonmonotone_p3 {r : ℝ → ℝ} (hr : IsRoundNearest 3 r) :
    (Ops.rounded r).sqrt (nmA 0 0) = 3 / 8 ∧
      (Ops.rounded r).div (nmA 1 0) (3 / 8) = 7 / 16 ∧
      pivot (Ops.rounded r) nmA 1 = 1 / 32 ∧ Succeeds (Ops.rounded r) nmA := by
  have hs := sqrt_eighth_bounds
  have h1 : (Ops.rounded r).sqrt (nmA 0 0) = 3 / 8 := by
    show r (Real.sqrt (4 / 32)) = 3 / 8
    rw [show (4 / 32 : ℝ) = 1 / 8 by norm_num]
    apply round_eq_of_near hr (-2)
    · exact ⟨3, -3, by norm_num, by norm_num⟩
    · exact ⟨6, by norm_num⟩
    · rw [abs_lt]; norm_num; constructor <;> linarith
    · rw [abs_of_pos (Real.sqrt_pos.mpr (by norm_num))]; norm_num; linarith
  have h2 : (Ops.rounded r).div (nmA 1 0) (3 / 8) = 7 / 16 := by
    show r (5 / 32 / (3 / 8)) = 7 / 16
    apply round_eq_of_near hr (-2)
    · exact ⟨7, -4, by norm_num, by norm_num⟩
    · exact ⟨7, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num
  have h3 : r (7 / 16 * (7 / 16)) = 3 / 16 := by
    apply round_eq_of_near hr (-3)
    · exact ⟨3, -4, by norm_num, by norm_num⟩
    · exact ⟨6, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num
  have h4 : r (7 / 32 - 3 / 16) = 1 / 32 := by
    rw [show (7 / 32 - 3 / 16 : ℝ) = 1 * 2 ^ (-5 : ℤ) by norm_num]
    have := round_eq_self hr (isFloat_mk 1 (-5) (by norm_num))
    push_cast at this; rw [this]; norm_num
  have hpiv : pivot (Ops.rounded r) nmA 1 = 1 / 32 := by
    rw [(pivot_two _ _).2, h1, h2]
    show r (nmA 1 1 - r (7 / 16 * (7 / 16))) = 1 / 32
    rw [h3]; exact h4
  refine ⟨h1, h2, hpiv, fun k => ?_⟩
  fin_cases k
  · show 0 < pivot _ nmA 0; rw [(pivot_two _ _).1]; simp [nmA]
  · show 0 < pivot _ nmA 1; rw [hpiv]; norm_num

/-- With `p = 4` the computed factor is `l₁₁ = 11/32`, `l₂₁ = 15/32`, and the second pivot
is `0`: the factorization fails. -/
theorem nonmonotone_p4 {r : ℝ → ℝ} (hr : IsRoundNearest 4 r) :
    (Ops.rounded r).sqrt (nmA 0 0) = 11 / 32 ∧
      (Ops.rounded r).div (nmA 1 0) (11 / 32) = 15 / 32 ∧
      pivot (Ops.rounded r) nmA 1 = 0 ∧ ¬ Succeeds (Ops.rounded r) nmA := by
  have hs := sqrt_eighth_bounds
  have h1 : (Ops.rounded r).sqrt (nmA 0 0) = 11 / 32 := by
    show r (Real.sqrt (4 / 32)) = 11 / 32
    rw [show (4 / 32 : ℝ) = 1 / 8 by norm_num]
    apply round_eq_of_near hr (-2)
    · exact ⟨11, -5, by norm_num, by norm_num⟩
    · exact ⟨11, by norm_num⟩
    · rw [abs_lt]; norm_num; constructor <;> linarith
    · rw [abs_of_pos (Real.sqrt_pos.mpr (by norm_num))]; norm_num; linarith
  have h2 : (Ops.rounded r).div (nmA 1 0) (11 / 32) = 15 / 32 := by
    show r (5 / 32 / (11 / 32)) = 15 / 32
    apply round_eq_of_near hr (-2)
    · exact ⟨15, -5, by norm_num, by norm_num⟩
    · exact ⟨15, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num
  have h3 : r (15 / 32 * (15 / 32)) = 7 / 32 := by
    apply round_eq_of_near hr (-3)
    · exact ⟨7, -5, by norm_num, by norm_num⟩
    · exact ⟨14, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num
  have hpiv : pivot (Ops.rounded r) nmA 1 = 0 := by
    rw [(pivot_two _ _).2, h1, h2]
    show r (nmA 1 1 - r (15 / 32 * (15 / 32))) = 0
    rw [h3]; simp [nmA, round_zero hr]
  refine ⟨h1, h2, hpiv, fun h => ?_⟩
  have := h 1
  rw [hpiv] at this
  exact lt_irrefl _ this

/-- All entries of the example are floats with `p = 3` (hence also with `p = 4`). -/
theorem nmA_float : ∀ i j, IsFloat 3 (nmA i j) := by
  intro i j
  fin_cases i <;> fin_cases j
  · exact ⟨1, -3, by norm_num, by simp [nmA]; norm_num⟩
  · exact ⟨5, -5, by norm_num, by simp [nmA]; norm_num⟩
  · exact ⟨5, -5, by norm_num, by simp [nmA]; norm_num⟩
  · exact ⟨7, -5, by norm_num, by simp [nmA]; norm_num⟩

/-! ## The order of operations matters -/

/-- With `p = 3`, `a = 1`, `b = c = 5/64`: `fl(fl(a - b) - c) = 3/4` but
`fl(a - fl(b + c)) = 7/8`. -/
theorem order_matters {r : ℝ → ℝ} (hr : IsRoundNearest 3 r) :
    r (r (1 - 5 / 64) - 5 / 64) = 3 / 4 ∧ r (1 - r (5 / 64 + 5 / 64)) = 7 / 8 := by
  have h1 : r (1 - 5 / 64) = 7 / 8 := by
    apply round_eq_of_near hr (-1)
    · exact ⟨7, -3, by norm_num, by norm_num⟩
    · exact ⟨7, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num
  have h2 : r (5 / 64 + 5 / 64) = 5 / 32 := by
    rw [show (5 / 64 + 5 / 64 : ℝ) = 5 * 2 ^ (-5 : ℤ) by norm_num]
    have := round_eq_self hr (isFloat_mk 5 (-5) (by norm_num))
    push_cast at this; rw [this]; norm_num
  refine ⟨?_, ?_⟩
  · rw [h1]
    apply round_eq_of_near hr (-1)
    · exact ⟨3, -2, by norm_num, by norm_num⟩
    · exact ⟨6, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num
  · rw [h2]
    apply round_eq_of_near hr (-1)
    · exact ⟨7, -3, by norm_num, by norm_num⟩
    · exact ⟨7, by norm_num⟩
    · rw [abs_lt]; norm_num
    · rw [abs_of_pos (by norm_num)]; norm_num

/-! ## The condition-number bound -/

/-- `(17/2 + 39x)(1 + 17/(72x)) = 289/(144x) + 425/24 + 39x`. -/
theorem conditionUpper_expand (x : ℝ) (hx : x ≠ 0) :
    conditionUpper x = 289 / (144 * x) + 425 / 24 + 39 * x := by
  unfold conditionUpper; field_simp; ring

/-- Among the three possible values of `x*`, the bound is largest at `x* = 2⁻¹²`. -/
theorem conditionUpper_largest :
    conditionUpper (1 / 256) ≤ conditionUpper (1 / 4096) ∧
      conditionUpper (1 / 1024) ≤ conditionUpper (1 / 4096) ∧
      conditionUpper (1 / 4096) = 303691615 / 36864 := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num [conditionUpper]

/-- The upper bounds in the table of numerical examples (rounded upward). -/
theorem table_bounds :
    conditionUpper ((2 : ℝ) ^ (-44 : ℤ)) ≤ 3.531e13 ∧
      conditionUpper ((2 : ℝ) ^ (-12 : ℤ)) < 8239 ∧
      conditionUpper ((2 : ℝ) ^ (-38 : ℤ)) ≤ 5.517e11 := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num [conditionUpper]

/-- The parameters of the three examples: `x = k^{3/2} u`. -/
theorem table_parameters :
    (CE.Par.ofPQ 53 3).x = (2 : ℝ) ^ (-44 : ℤ) ∧ 9 * (CE.Par.ofPQ 53 3).k = 576 ∧
      (CE.Par.ofPQ 24 4).x = (2 : ℝ) ^ (-12 : ℤ) ∧ 9 * (CE.Par.ofPQ 24 4).k = 2304 ∧
      (CE.Par.ofPQ 53 5).x = (2 : ℝ) ^ (-38 : ℤ) ∧ 9 * (CE.Par.ofPQ 53 5).k = 9216 := by
  have a1 : CE.Admissible 53 3 := ⟨by norm_num, by norm_num, by norm_num⟩
  have a2 : CE.Admissible 24 4 := ⟨by norm_num, by norm_num, by norm_num⟩
  have a3 : CE.Admissible 53 5 := ⟨by norm_num, by norm_num, by norm_num⟩
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [CE.Par.x_zpow, CE.Par.ofPQ_q a1, CE.Par.ofPQ_p a1]; norm_num
  · rw [CE.Par.ofPQ_k a1]; norm_num
  · rw [CE.Par.x_zpow, CE.Par.ofPQ_q a2, CE.Par.ofPQ_p a2]; norm_num
  · rw [CE.Par.ofPQ_k a2]; norm_num
  · rw [CE.Par.x_zpow, CE.Par.ofPQ_q a3, CE.Par.ofPQ_p a3]; norm_num
  · rw [CE.Par.ofPQ_k a3]; norm_num

/-- The larger binary64 bound is about `64` times smaller. -/
theorem table_ratio :
    63 < conditionUpper ((2 : ℝ) ^ (-44 : ℤ)) / conditionUpper ((2 : ℝ) ^ (-38 : ℤ)) ∧
      conditionUpper ((2 : ℝ) ^ (-44 : ℤ)) / conditionUpper ((2 : ℝ) ^ (-38 : ℤ)) < 65 := by
  constructor <;> norm_num [conditionUpper]

/-- The next admissible binary64 size is `147456`; a dense binary64 matrix of that size needs
`162 GiB`, and one of size `19456` needs about `2.82 GiB`. -/
theorem memory_sizes :
    9 * 4 ^ 7 = 147456 ∧ 147456 ^ 2 * 8 = 162 * 2 ^ 30 ∧
      (2.82 : ℝ) < 19456 ^ 2 * 8 / 2 ^ 30 ∧ (19456 : ℝ) ^ 2 * 8 / 2 ^ 30 < 2.83 := by
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num⟩

/-! ## The `1.1` corollary -/

/-- `∑_{j=1}^{n-1} (√j + 2) < (2/3) n^{3/2} + 2n ≤ (3/4) n^{3/2}` for `n ≥ 576`. -/
theorem sum_update_weights_steps (n : ℕ) (hn : 576 ≤ n) :
    (∑ j ∈ Finset.range (n - 1), (Real.sqrt (j + 1 : ℕ) + 2)) <
        (2 / 3 : ℝ) * n * Real.sqrt n + 2 * n ∧
      (2 / 3 : ℝ) * n * Real.sqrt n + 2 * n ≤ (3 / 4 : ℝ) * n * Real.sqrt n := by
  have hn1 : 1 ≤ n := by omega
  have hnR : (576 : ℝ) ≤ n := by exact_mod_cast hn
  have hsqrt : (24 : ℝ) ≤ Real.sqrt n := by
    have h := Real.sqrt_le_sqrt hnR
    rwa [show Real.sqrt 576 = 24 by
      rw [show (576 : ℝ) = 24 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]] at h
  constructor
  · have hshift := Finset.sum_range_succ' (fun j : ℕ => Real.sqrt (j : ℝ)) (n - 1)
    rw [Nat.sub_add_cancel hn1] at hshift
    simp only [Nat.cast_zero, Real.sqrt_zero, add_zero] at hshift
    have hsum := sum_sqrt_le n
    rw [hshift] at hsum
    have hlast : 0 < Real.sqrt ((n - 1 : ℕ) : ℝ) := Real.sqrt_pos.mpr (by
      have : 2 ≤ n := by omega
      exact_mod_cast (by omega : 0 < n - 1))
    have hnsub : ((n - 1 : ℕ) : ℝ) = n - 1 := by simp [Nat.cast_sub hn1]
    simp only [Finset.sum_add_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
    rw [hnsub]
    push_cast at hsum ⊢
    nlinarith
  · nlinarith [mul_nonneg (show 0 ≤ (n : ℝ) by positivity) (sub_nonneg.mpr hsqrt)]

/-- `3 / (4.4 (1 - 5 · 2^{-17})) < log 2`. -/
theorem meinguet_constant : 3 / (4.4 * (1 - 5 * (1 / 131072 : ℝ))) < Real.log 2 := by
  have := Real.log_two_gt_d9
  norm_num at this ⊢
  linarith

/-! ## Scaling remarks -/

/-- `H = D_A⁻¹ A D_A⁻¹` with `D_A = diag(√aᵢᵢ)` has unit diagonal. -/
theorem scaled_unit_diag {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (hA : ∀ i, 0 < A i i) (i : Fin n) :
    (Matrix.diagonal (fun i => 1 / Real.sqrt (A i i)) * A *
      Matrix.diagonal (fun i => 1 / Real.sqrt (A i i))) i i = 1 := by
  simp only [Matrix.diagonal_mul, Matrix.mul_diagonal]
  have h := Real.sqrt_pos.mpr (hA i)
  have h2 := Real.sq_sqrt (hA i).le
  field_simp
  rw [h2]

/-- For a positive diagonal matrix, `H = I`. -/
theorem scaled_diag_eq_one {n : ℕ} (d : Fin n → ℝ) (hd : ∀ i, 0 < d i) :
    Matrix.diagonal (fun i => 1 / Real.sqrt (d i)) * Matrix.diagonal d *
      Matrix.diagonal (fun i => 1 / Real.sqrt (d i)) = 1 := by
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
  congr 1; funext i
  have h := Real.sqrt_pos.mpr (hd i)
  have h2 := Real.sq_sqrt (hd i).le
  field_simp
  rw [h2]

lemma permMatrix_transpose_mul {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    (σ.permMatrix ℝ)ᵀ * σ.permMatrix ℝ = 1 := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Equiv.Perm.permMatrix,
    PEquiv.toMatrix_apply, Equiv.toPEquiv_apply, Option.mem_def, Option.some.injEq,
    Matrix.one_apply]
  rw [Finset.sum_eq_single (σ.symm i)]
  · simp
  · intro b _ hb
    have : σ b ≠ i := fun h => hb (by rw [← h]; simp)
    simp [this]
  · simp

lemma permMatrix_mul_transpose {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    σ.permMatrix ℝ * (σ.permMatrix ℝ)ᵀ = 1 := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Equiv.Perm.permMatrix,
    PEquiv.toMatrix_apply, Equiv.toPEquiv_apply, Option.mem_def, Option.some.injEq,
    Matrix.one_apply]
  rw [Finset.sum_eq_single (σ i)]
  · simp only [if_true]; split_ifs with h1 h2 h2 <;> simp_all [σ.injective.eq_iff]
  · intro b _ hb
    simp [Ne.symm hb]
  · simp

/-- `κ₂(P A Pᵀ) = κ₂(A)` for a permutation matrix `P`. -/
theorem cond2_perm {n : ℕ} (σ : Equiv.Perm (Fin n)) (A : Matrix (Fin n) (Fin n) ℝ) :
    cond2 (σ.permMatrix ℝ * A * (σ.permMatrix ℝ)ᵀ) = cond2 A :=
  cond2_orthogonal_conj (permMatrix_transpose_mul σ) (permMatrix_mul_transpose σ)

end Cholesky
