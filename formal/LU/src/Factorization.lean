import src.Algorithm

/-!
# Exact arithmetic: `P A = L U` for every square matrix

With exact arithmetic over `𝕜` (any normed field, in particular `ℝ` and `ℂ`):

* `exact_plu`: the routine computes `P A = L U`, where `A` is the original input,
  `P` the accumulated interchanges, `L` (unit lower triangular) and `U` (upper triangular)
  the packed factors. The induction follows the proof of the existence theorem:
  `exact_v_eq` is `v = α ℓ` (with `ℓ = 0` and `v = 0` when `α = 0`) and `exact_trail_eq` is
  `S = B - ℓ w`.
* `exists_lu_pivoting` (theorem `thm:existence_lu_pivoting`), for `ℝ` and `ℂ`.
* `norm_multiplier_le_one`: every stored multiplier has `|l_ik| ≤ 1`.
* `pivot_eq_diag`: the pivots are the diagonal entries of `U`;
  `pivot_eq_zero_col`: a zero pivot means the whole active column is zero.
* `pivot_ne_zero`: a nonsingular matrix has no zero pivots;
  `exists_diag_eq_zero`: a singular matrix has a zero diagonal entry in `U`, and then
  `Ax = b` is not uniquely solvable for every `b` (`not_unique_of_singular`).
* `solve_of_plu`: `L y = P b` and `U x = y` give `A x = b`.
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {𝕜 : Type*} [NormedField 𝕜] [DecidableEq 𝕜] {n m : ℕ}

/-! ## The first step in exact arithmetic -/

section ExactFirst
variable (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜)

/-- With `α = (P₁A)₀₀`, `v` the rest of the first column of `P₁ A`, and `ℓ` the stored first
column: `v = α ℓ`. -/
lemma exact_v_eq (i : Fin m) :
    swapped A i.succ 0 = swapped A 0 0 * (first (Ops.exact 𝕜) A).mat i.succ 0 := by
  rw [first_col]
  split_ifs with h
  · rw [h, zero_mul]; exact swapped_col_eq_zero A h _
  · rw [Ops.exact_div]; field_simp

/-- If `α = 0` then `ℓ = 0`. -/
lemma exact_mult_of_pivot_zero (h : swapped A 0 0 = 0) (i : Fin m) :
    (first (Ops.exact 𝕜) A).mat i.succ 0 = 0 := by
  rw [first_col, if_pos h]; exact swapped_col_eq_zero A h _

/-- The block the algorithm continues with is the Schur complement `S = B - ℓ w`. -/
lemma exact_trail_eq (i j : Fin m) :
    trail (Ops.exact 𝕜) A i j =
      swapped A i.succ j.succ - (first (Ops.exact 𝕜) A).mat i.succ 0 * swapped A 0 j.succ := by
  rw [trail_apply]
  split_ifs with h
  · rw [exact_mult_of_pivot_zero A h, zero_mul, sub_zero]
  · rw [Ops.exact_upd, first_col, if_neg h, Ops.exact_div]

end ExactFirst

/-! ## `P A = L U` -/

/-- **The routine computes `P A = L U`** (exact arithmetic). -/
theorem exact_plu : ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) 𝕜),
    factorP (Ops.exact 𝕜) A * A = factorL (Ops.exact 𝕜) A * factorU (Ops.exact 𝕜) A
  | 0, A => by ext i; exact Fin.elim0 i
  | m + 1, A => by
    have ih := exact_plu (trail (Ops.exact 𝕜) A)
    unfold factorP factorPerm factorL factorU at ih ⊢
    ext i j
    refine Fin.cases ?_ (fun i => ?_) i
    · rw [permMul_zero, mul_zero_apply, final_row]
    · refine Fin.cases ?_ (fun j => ?_) j
      · rw [permMul_succ, mul_succ_zero, final_col, final_row, exact_v_eq, mul_comm]
      · rw [permMul_succ, mul_succ_succ, final_col, final_row, final_trail, ← ih,
          permMatrix_mul_apply, exact_trail_eq]
        ring

/-- **Existence of LU factorization with row pivoting** (`thm:existence_lu_pivoting`). -/
theorem exists_lu_pivoting (A : Matrix (Fin n) (Fin n) 𝕜) :
    ∃ (σ : Equiv.Perm (Fin n)) (L U : Matrix (Fin n) (Fin n) 𝕜),
      IsUnitLower L ∧ IsUpper U ∧ σ.permMatrix 𝕜 * A = L * U :=
  ⟨_, _, _, lowerOf_isUnitLower _, upperOf_isUpper _, exact_plu A⟩

theorem exists_lu_pivoting_real (A : Matrix (Fin n) (Fin n) ℝ) :
    ∃ (σ : Equiv.Perm (Fin n)) (L U : Matrix (Fin n) (Fin n) ℝ),
      IsUnitLower L ∧ IsUpper U ∧ σ.permMatrix ℝ * A = L * U :=
  exists_lu_pivoting A

theorem exists_lu_pivoting_complex (A : Matrix (Fin n) (Fin n) ℂ) :
    ∃ (σ : Equiv.Perm (Fin n)) (L U : Matrix (Fin n) (Fin n) ℂ),
      IsUnitLower L ∧ IsUpper U ∧ σ.permMatrix ℂ * A = L * U :=
  exists_lu_pivoting A

/-! ## The multipliers -/

/-- **`|l_ik| ≤ 1`.** Every multiplier stored below the diagonal has magnitude at most one. -/
theorem norm_multiplier_le_one : ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) 𝕜) (i j : Fin n),
    j < i → ‖(final (Ops.exact 𝕜) A).mat i j‖ ≤ 1
  | 0, _, i, _, _ => Fin.elim0 i
  | m + 1, A, i, j, hij => by
    induction i using Fin.cases with
    | zero => exact absurd hij (not_lt.mpr (Fin.zero_le _))
    | succ i =>
      induction j using Fin.cases with
      | zero =>
        rw [final_col, first_col]
        split_ifs with h
        · rw [swapped_col_eq_zero A h, norm_zero]; exact zero_le_one
        · rw [Ops.exact_div, norm_div, div_le_one (norm_pos_iff.mpr h)]
          exact norm_swapped_le A _
      | succ j =>
        have := norm_multiplier_le_one (trail (Ops.exact 𝕜) A) i j (Fin.succ_lt_succ_iff.mp hij)
        rw [← final_trail] at this
        exact this

theorem norm_factorL_le_one (A : Matrix (Fin n) (Fin n) 𝕜) (i j : Fin n) :
    ‖factorL (Ops.exact 𝕜) A i j‖ ≤ 1 := by
  unfold factorL lowerOf
  simp only [Matrix.of_apply]
  split_ifs with h1 h2
  · rw [norm_one]
  · exact norm_multiplier_le_one A i j h2
  · rw [norm_zero]; exact zero_le_one

/-! ## The pivots -/

/-- The pivot selected at step `k`: the entry of largest magnitude in column `k` of the
active rows. -/
noncomputable def pivot (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) (k : Fin n) : 𝕜 :=
  (run o k ⟨A, 1⟩).mat (pivRow (run o k ⟨A, 1⟩).mat k) k

/-- **A zero pivot means that the entire active column is zero.** -/
theorem pivot_eq_zero_col (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) (k : Fin n)
    (h : pivot o A k = 0) (i : Fin n) (hi : k ≤ i) : (run o k ⟨A, 1⟩).mat i k = 0 := by
  have := le_at_firstMax (fun i => ‖(run o k ⟨A, 1⟩).mat i k‖) k hi
  simp only [] at this
  rw [show (run o k ⟨A, 1⟩).mat (firstMax (fun i => ‖(run o k ⟨A, 1⟩).mat i k‖) k) k = 0 from h,
    norm_zero] at this
  exact norm_eq_zero.mp (le_antisymm this (norm_nonneg _))

lemma rel_run (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) (t : ℕ) :
    Rel (run o (t + 1) ⟨A, 1⟩) (run o t ⟨trail o A, 1⟩) ((first o A).mat 0)
      (fun i => (first o A).mat i.succ 0) (first o A).perm := by
  have : run o (t + 1) ⟨A, 1⟩ = runFrom o 1 t (first o A) := by
    simp only [run, runFrom, LU.stepN, dif_pos (Nat.succ_pos m)]
    rfl
  rw [this]
  exact (rel_first o A).runFrom 0 t

lemma pivot_zero (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) :
    pivot o A 0 = swapped A 0 0 := by
  simp [pivot, run, runFrom, swapped, swapRows]

lemma pivot_succ (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) (k : Fin m) :
    pivot o A k.succ = pivot o (trail o A) k := by
  unfold pivot
  have h := rel_run o A k
  simp only [Fin.val_succ]
  rw [pivRow_succ, h.submatrix, h.trail]

/-- **The pivots are the diagonal entries of `U`** (in any arithmetic). -/
theorem pivot_eq_diag (o : Ops 𝕜) : ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) 𝕜) (k : Fin n),
    pivot o A k = factorU o A k k
  | 0, _, k => Fin.elim0 k
  | m + 1, A, k => by
    unfold factorU
    induction k using Fin.cases with
    | zero => rw [pivot_zero, upperOf_zero, final_row]
    | succ k =>
      rw [pivot_succ, pivot_eq_diag o (trail o A) k, upperOf_succ_succ, final_trail]
      rfl

lemma det_permMatrix_ne_zero (σ : Equiv.Perm (Fin n)) : (σ.permMatrix 𝕜).det ≠ 0 := by
  rw [Matrix.det_permutation]
  rcases Int.units_eq_one_or (Equiv.Perm.sign σ) with h | h <;> simp [h]

lemma factorL_isUnitLower (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) :
    IsUnitLower (factorL o A) := lowerOf_isUnitLower _

lemma factorU_isUpper (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) : IsUpper (factorU o A) :=
  upperOf_isUpper _

lemma det_factorU (A : Matrix (Fin n) (Fin n) 𝕜) :
    (factorP (Ops.exact 𝕜) A).det * A.det = ∏ i, factorU (Ops.exact 𝕜) A i i := by
  rw [← Matrix.det_mul, exact_plu, Matrix.det_mul, det_of_isUnitLower (factorL_isUnitLower _ _),
    one_mul, det_of_isUpper (factorU_isUpper _ _)]

/-- **A nonsingular matrix has no zero pivots** (exact arithmetic). -/
theorem pivot_ne_zero (A : Matrix (Fin n) (Fin n) 𝕜) (hA : A.det ≠ 0) (k : Fin n) :
    pivot (Ops.exact 𝕜) A k ≠ 0 := by
  rw [pivot_eq_diag]
  intro h0
  have := det_factorU A
  rw [Finset.prod_eq_zero (f := fun i => factorU (Ops.exact 𝕜) A i i) (Finset.mem_univ k) h0]
    at this
  exact mul_ne_zero (det_permMatrix_ne_zero _) hA this

/-- **For a singular matrix, `U` has a zero diagonal entry.** -/
theorem exists_diag_eq_zero (A : Matrix (Fin n) (Fin n) 𝕜) (hA : A.det = 0) :
    ∃ k, factorU (Ops.exact 𝕜) A k k = 0 := by
  have := det_factorU A
  rw [hA, mul_zero, eq_comm, Finset.prod_eq_zero_iff] at this
  obtain ⟨k, -, hk⟩ := this
  exact ⟨k, hk⟩

/-- For singular `A`, `A x = b` is not uniquely solvable for every right-hand side. -/
theorem not_unique_of_singular (A : Matrix (Fin n) (Fin n) 𝕜) (hA : A.det = 0) :
    ¬ ∀ b : Fin n → 𝕜, ∃! x, A *ᵥ x = b := by
  intro h
  obtain ⟨v, hv, hAv⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hA
  obtain ⟨x, -, hx⟩ := h 0
  have h1 := hx v hAv
  have h2 := hx 0 (Matrix.mulVec_zero A)
  exact hv (h1.trans h2.symm)

/-! ## Solving `A x = b` -/

lemma permMatrix_inv_mul (σ : Equiv.Perm (Fin n)) :
    (σ⁻¹).permMatrix 𝕜 * σ.permMatrix 𝕜 = 1 := by
  ext i j
  rw [permMatrix_mul_apply]
  simp [Equiv.Perm.permMatrix, PEquiv.toMatrix_apply, Matrix.one_apply]

/-- **Solving with the factors.** If `P A = L U`, `L y = P b` and `U x = y`, then `A x = b`:
no permutation of `x` has to be undone. -/
theorem solve_of_plu {σ : Equiv.Perm (Fin n)} {A L U : Matrix (Fin n) (Fin n) 𝕜}
    (h : σ.permMatrix 𝕜 * A = L * U) {b y x : Fin n → 𝕜}
    (hy : L *ᵥ y = σ.permMatrix 𝕜 *ᵥ b) (hx : U *ᵥ x = y) : A *ᵥ x = b := by
  have h1 : σ.permMatrix 𝕜 *ᵥ (A *ᵥ x) = σ.permMatrix 𝕜 *ᵥ b := by
    rw [Matrix.mulVec_mulVec, h, ← Matrix.mulVec_mulVec, hx, hy]
  have h2 := congrArg (fun v => (σ⁻¹).permMatrix 𝕜 *ᵥ v) h1
  simp only [Matrix.mulVec_mulVec, ← Matrix.mul_assoc, permMatrix_inv_mul, Matrix.one_mul,
    Matrix.one_mulVec] at h2
  exact h2

end LU
