import src.Ops
import src.Pivot
import src.Packed

/-!
# Gaussian elimination with partial pivoting, in place

This file models `lu_factorization_with_row_pivoting`. The state is the array together
with the accumulated row permutation (the routine's `P`). Step `k` (the notes count from 1)

1. chooses `p = pivRow M k`, the first row `p ≥ k` with `|m_pk|` maximal
   (`k + np.argmax(np.abs(A[k:, k]))`);
2. interchanges the entire rows `k` and `p` of the array and of `P` (`swapRows`);
3. if the pivot `m_kk` is nonzero, stores the multipliers `l_ik = m_ik / m_kk` below it and
   updates `m_ij ← m_ij - l_ik m_kj` for `i, j > k`; if it is zero, nothing is done
   (`elim`).

The arithmetic is a parameter `o : Ops 𝕜`. The routine runs `k = 0, …, n - 2`; we run
`k = 0, …, n - 1`, which is the same because the last step changes nothing
(`run_pred`).

The main structural fact is the reduction to the trailing block (`rel_final`): after the
first step, the remaining steps act on the trailing block exactly as the routine acts on a
matrix of size `n - 1`; the first row stays as it is, and the stored multipliers of the first
column are permuted by the later interchanges. This is the induction in the proof of the
existence theorem.
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

/-- The state of the in-place routine: the array and the row permutation. Row `i` of the
array descends from row `perm i` of the input. -/
structure State (𝕜 : Type*) (n : ℕ) where
  mat : Matrix (Fin n) (Fin n) 𝕜
  perm : Equiv.Perm (Fin n)

variable {𝕜 : Type*} [NormedField 𝕜] [DecidableEq 𝕜] {n m : ℕ}

/-- The pivot row chosen at step `k`. -/
noncomputable def pivRow (M : Matrix (Fin n) (Fin n) 𝕜) (k : Fin n) : Fin n :=
  firstMax (fun i => ‖M i k‖) k

/-- Interchange of the entire rows `k` and `p`. -/
def swapRows (M : Matrix (Fin n) (Fin n) 𝕜) (k p : Fin n) : Matrix (Fin n) (Fin n) 𝕜 :=
  M.submatrix (Equiv.swap k p) id

/-- Elimination with the pivot `m_kk` (after the interchange). -/
def elim (o : Ops 𝕜) (k : Fin n) (M : Matrix (Fin n) (Fin n) 𝕜) : Matrix (Fin n) (Fin n) 𝕜 :=
  if M k k = 0 then M else
    Matrix.of fun i j =>
      if k < i ∧ j = k then o.div (M i k) (M k k)
      else if k < i ∧ k < j then o.upd (M i j) (o.div (M i k) (M k k)) (M k j)
      else M i j

/-- Step `k` of the routine. -/
noncomputable def stepAt (o : Ops 𝕜) (k : Fin n) (S : State 𝕜 n) : State 𝕜 n :=
  ⟨elim o k (swapRows S.mat k (pivRow S.mat k)), S.perm * Equiv.swap k (pivRow S.mat k)⟩

/-- Step `s`, numbered by a natural number (a no-op when `s ≥ n`). -/
noncomputable def stepN (o : Ops 𝕜) (s : ℕ) (S : State 𝕜 n) : State 𝕜 n :=
  if h : s < n then stepAt o ⟨s, h⟩ S else S

/-- Steps `s, s + 1, …, s + t - 1`. -/
noncomputable def runFrom (o : Ops 𝕜) : ℕ → ℕ → State 𝕜 n → State 𝕜 n
  | _, 0, S => S
  | s, t + 1, S => runFrom o (s + 1) t (stepN o s S)

/-- The state after the first `t` steps. -/
noncomputable def run (o : Ops 𝕜) (t : ℕ) (S : State 𝕜 n) : State 𝕜 n := runFrom o 0 t S

/-- The final state, starting from the input `A` and `P = I`. -/
noncomputable def final (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) : State 𝕜 n :=
  run o n ⟨A, 1⟩

/-- The computed unit lower triangular factor. -/
noncomputable def factorL (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) :
    Matrix (Fin n) (Fin n) 𝕜 := lowerOf (final o A).mat

/-- The computed upper triangular factor. -/
noncomputable def factorU (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) :
    Matrix (Fin n) (Fin n) 𝕜 := upperOf (final o A).mat

/-- The computed row permutation. -/
noncomputable def factorPerm (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) : Equiv.Perm (Fin n) :=
  (final o A).perm

/-- The permutation matrix `P`: `(P A)ᵢⱼ = a_{perm i, j}`. -/
noncomputable def factorP (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) :
    Matrix (Fin n) (Fin n) 𝕜 := (factorPerm o A).permMatrix 𝕜

lemma permMatrix_mul_apply (σ : Equiv.Perm (Fin n)) (A : Matrix (Fin n) (Fin n) 𝕜) (i j : Fin n) :
    (σ.permMatrix 𝕜 * A) i j = A (σ i) j := by
  rw [PEquiv.toMatrix_toPEquiv_mul]; rfl

/-! ## The loop -/

lemma runFrom_succ' (o : Ops 𝕜) (s t : ℕ) (S : State 𝕜 n) :
    runFrom o s (t + 1) S = stepN o (s + t) (runFrom o s t S) := by
  induction t generalizing s S with
  | zero => simp [runFrom]
  | succ t ih =>
    rw [runFrom, ih, runFrom]
    congr 1
    omega

lemma run_succ (o : Ops 𝕜) (t : ℕ) (S : State 𝕜 n) :
    run o (t + 1) S = stepN o t (run o t S) := by
  simp [run, runFrom_succ']

lemma stepN_of_le (o : Ops 𝕜) {s : ℕ} (hs : n ≤ s) (S : State 𝕜 n) : stepN o s S = S := by
  simp [stepN, not_lt.mpr hs]

/-! ## Elimination -/

section Elim
variable {o : Ops 𝕜} {k i j : Fin n} {M : Matrix (Fin n) (Fin n) 𝕜}

lemma elim_of_pivot_eq_zero (h : M k k = 0) : elim o k M = M := by
  simp [elim, h]

lemma elim_apply_of_not_lt (hi : ¬ k < i) : elim o k M i j = M i j := by
  unfold elim
  split_ifs
  · rfl
  · simp [hi]

lemma elim_apply_of_lt_col (hj : j < k) : elim o k M i j = M i j := by
  unfold elim
  split_ifs
  · rfl
  · simp [ne_of_lt hj, not_lt.mpr hj.le]

lemma elim_apply_mult (h : M k k ≠ 0) (hi : k < i) :
    elim o k M i k = o.div (M i k) (M k k) := by
  simp [elim, h, hi]

lemma elim_apply_upd (h : M k k ≠ 0) (hi : k < i) (hj : k < j) :
    elim o k M i j = o.upd (M i j) (o.div (M i k) (M k k)) (M k j) := by
  simp [elim, h, hi, hj, ne_of_gt hj]

/-- The pivot row is not changed by elimination. -/
lemma elim_apply_self (j : Fin n) : elim o k M k j = M k j :=
  elim_apply_of_not_lt (lt_irrefl k)

end Elim

/-! ## The last step does nothing -/

lemma pivRow_last (M : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) :
    pivRow M (Fin.last m) = Fin.last m :=
  firstMax_last _

lemma stepAt_last (o : Ops 𝕜) (S : State 𝕜 (m + 1)) : stepAt o (Fin.last m) S = S := by
  have hsw : swapRows S.mat (Fin.last m) (Fin.last m) = S.mat := by
    ext i j; simp [swapRows]
  have hel : elim o (Fin.last m) S.mat = S.mat := by
    ext i j
    exact elim_apply_of_not_lt (not_lt.mpr (Fin.le_last i))
  simp only [stepAt, pivRow_last, hsw, hel, Equiv.swap_self, Equiv.Perm.mul_refl]

/-- The routine stops after step `n - 2`; running step `n - 1` as well changes nothing. -/
lemma run_pred (o : Ops 𝕜) (S : State 𝕜 n) : run o n S = run o (n - 1) S := by
  cases n with
  | zero => rfl
  | succ m =>
    rw [run_succ, Nat.add_sub_cancel]
    unfold stepN
    rw [dif_pos (Nat.lt_succ_self m)]
    exact stepAt_last o _

/-! ## Reduction to the trailing block -/

/-- The extension of a permutation of `{1, …, m}` that fixes `0`. -/
def liftPerm (σ : Equiv.Perm (Fin m)) : Equiv.Perm (Fin (m + 1)) where
  toFun i := Fin.cases 0 (fun j => (σ j).succ) i
  invFun i := Fin.cases 0 (fun j => (σ.symm j).succ) i
  left_inv i := by refine Fin.cases ?_ (fun j => ?_) i <;> simp
  right_inv i := by refine Fin.cases ?_ (fun j => ?_) i <;> simp

@[simp] lemma liftPerm_zero (σ : Equiv.Perm (Fin m)) : liftPerm σ 0 = 0 := rfl

@[simp] lemma liftPerm_succ (σ : Equiv.Perm (Fin m)) (i : Fin m) :
    liftPerm σ i.succ = (σ i).succ := rfl

lemma liftPerm_mul (σ τ : Equiv.Perm (Fin m)) :
    liftPerm (σ * τ) = liftPerm σ * liftPerm τ :=
  Equiv.ext fun i => by refine Fin.cases ?_ (fun j => ?_) i <;> simp

lemma liftPerm_one : liftPerm (1 : Equiv.Perm (Fin m)) = 1 :=
  Equiv.ext fun i => by refine Fin.cases ?_ (fun j => ?_) i <;> simp

lemma swap_succ_apply_zero (k p : Fin m) : Equiv.swap k.succ p.succ 0 = 0 :=
  Equiv.swap_apply_of_ne_of_ne (Fin.succ_ne_zero k).symm (Fin.succ_ne_zero p).symm

lemma swap_succ_apply_succ (k p i : Fin m) :
    Equiv.swap k.succ p.succ i.succ = (Equiv.swap k p i).succ := by
  simp only [Equiv.swap_apply_def, Fin.succ_inj]
  split_ifs <;> rfl

lemma swap_succ_eq_liftPerm (k p : Fin m) :
    Equiv.swap k.succ p.succ = liftPerm (Equiv.swap k p) :=
  Equiv.ext fun i => by
    refine Fin.cases ?_ (fun j => ?_) i
    · rw [swap_succ_apply_zero, liftPerm_zero]
    · rw [swap_succ_apply_succ, liftPerm_succ]

section Shift
variable (M : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜)

lemma pivRow_succ (k : Fin m) :
    pivRow M k.succ = (pivRow (M.submatrix Fin.succ Fin.succ) k).succ := by
  unfold pivRow
  rw [firstMax_succ]
  rfl

lemma swapRows_succ (k p : Fin m) :
    (swapRows M k.succ p.succ).submatrix Fin.succ Fin.succ =
      swapRows (M.submatrix Fin.succ Fin.succ) k p := by
  ext i j
  simp [swapRows, swap_succ_apply_succ]

lemma swapRows_succ_zero (k p : Fin m) (j : Fin (m + 1)) :
    swapRows M k.succ p.succ 0 j = M 0 j := by
  simp [swapRows, swap_succ_apply_zero]

lemma swapRows_succ_col (k p i : Fin m) :
    swapRows M k.succ p.succ i.succ 0 = M (Equiv.swap k p i).succ 0 := by
  simp [swapRows, swap_succ_apply_succ]

lemma elim_succ (o : Ops 𝕜) (k : Fin m) :
    (elim o k.succ M).submatrix Fin.succ Fin.succ = elim o k (M.submatrix Fin.succ Fin.succ) := by
  ext i j
  simp only [Matrix.submatrix_apply]
  by_cases h : M k.succ k.succ = 0
  · rw [elim_of_pivot_eq_zero h, elim_of_pivot_eq_zero (M := M.submatrix Fin.succ Fin.succ) h]
    rfl
  · simp only [elim, if_neg h, Matrix.of_apply, Matrix.submatrix_apply,
      Fin.succ_lt_succ_iff, Fin.succ_inj]

lemma elim_succ_zero (o : Ops 𝕜) (k : Fin m) (j : Fin (m + 1)) :
    elim o k.succ M 0 j = M 0 j :=
  elim_apply_of_not_lt (not_lt.mpr (Fin.zero_le _))

lemma elim_succ_col (o : Ops 𝕜) (k : Fin m) (i : Fin (m + 1)) :
    elim o k.succ M i 0 = M i 0 :=
  elim_apply_of_lt_col (Fin.succ_pos k)

end Shift

/-- `S` (of size `m + 1`) is described by the state `t` of size `m`: its first row is `row0`,
its first column below the diagonal holds the data `col0` permuted by `t.perm`, its
trailing block is `t.mat`, and its permutation is `π₀` followed by `t.perm`. -/
structure Rel (S : State 𝕜 (m + 1)) (t : State 𝕜 m) (row0 : Fin (m + 1) → 𝕜)
    (col0 : Fin m → 𝕜) (π₀ : Equiv.Perm (Fin (m + 1))) : Prop where
  row : ∀ j, S.mat 0 j = row0 j
  col : ∀ i, S.mat i.succ 0 = col0 (t.perm i)
  trail : ∀ i j, S.mat i.succ j.succ = t.mat i j
  perm : S.perm = π₀ * liftPerm t.perm

section Rel
variable {o : Ops 𝕜} {S : State 𝕜 (m + 1)} {t : State 𝕜 m} {row0 : Fin (m + 1) → 𝕜}
  {col0 : Fin m → 𝕜} {π₀ : Equiv.Perm (Fin (m + 1))}

lemma Rel.submatrix (h : Rel S t row0 col0 π₀) : S.mat.submatrix Fin.succ Fin.succ = t.mat := by
  ext i j; exact h.trail i j

/-- Step `k + 1` of the large matrix is step `k` of the trailing block. -/
lemma Rel.stepAt (h : Rel S t row0 col0 π₀) (k : Fin m) :
    Rel (stepAt o k.succ S) (stepAt o k t) row0 col0 π₀ := by
  have hp : pivRow S.mat k.succ = (pivRow t.mat k).succ := by rw [pivRow_succ, h.submatrix]
  refine ⟨fun j => ?_, fun i => ?_, fun i j => ?_, ?_⟩
  · simp only [LU.stepAt, hp]
    rw [elim_succ_zero, swapRows_succ_zero, h.row]
  · simp only [LU.stepAt, hp]
    rw [elim_succ_col, swapRows_succ_col, h.col]
    rfl
  · simp only [LU.stepAt, hp]
    have := congrFun (congrFun (elim_succ (swapRows S.mat k.succ (pivRow t.mat k).succ) o k) i) j
    simp only [Matrix.submatrix_apply] at this
    rw [this, swapRows_succ, h.submatrix]
  · simp only [LU.stepAt, hp]
    rw [h.perm, swap_succ_eq_liftPerm, liftPerm_mul, mul_assoc]

lemma Rel.stepN (h : Rel S t row0 col0 π₀) (s : ℕ) :
    Rel (stepN o (s + 1) S) (stepN o s t) row0 col0 π₀ := by
  unfold LU.stepN
  by_cases hs : s < m
  · rw [dif_pos (by omega : s + 1 < m + 1), dif_pos hs]
    exact h.stepAt ⟨s, hs⟩
  · rw [dif_neg (by omega : ¬ s + 1 < m + 1), dif_neg hs]
    exact h

lemma Rel.runFrom (h : Rel S t row0 col0 π₀) (s q : ℕ) :
    Rel (runFrom o (s + 1) q S) (runFrom o s q t) row0 col0 π₀ := by
  induction q generalizing s S t with
  | zero => exact h
  | succ q ih => exact ih (h.stepN s) (s + 1)

end Rel

/-- The state after the first step. -/
noncomputable def first (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) :
    State 𝕜 (m + 1) := stepAt o 0 ⟨A, 1⟩

/-- The trailing block after the first step: the matrix the routine continues with. -/
noncomputable def trail (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) :
    Matrix (Fin m) (Fin m) 𝕜 := (first o A).mat.submatrix Fin.succ Fin.succ

lemma rel_first (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) :
    Rel (first o A) ⟨trail o A, 1⟩ ((first o A).mat 0) (fun i => (first o A).mat i.succ 0)
      (first o A).perm :=
  ⟨fun _ => rfl, fun _ => rfl, fun _ _ => rfl, by rw [liftPerm_one, mul_one]⟩

/-- **Reduction to the trailing block.** The final state of a matrix of size `m + 1` is
described by the final state of `trail o A`. -/
lemma rel_final (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜) :
    Rel (final o A) (final o (trail o A)) ((first o A).mat 0)
      (fun i => (first o A).mat i.succ 0) (first o A).perm := by
  have : final o A = runFrom o 1 m (first o A) := by
    simp only [final, run, runFrom, LU.stepN, dif_pos (Nat.succ_pos m)]
    rfl
  rw [this]
  exact (rel_first o A).runFrom 0 m

/-! ## The first step -/

section First
variable (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜)

/-- The input after the first interchange, `P₁ A`. -/
noncomputable def swapped : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜 := swapRows A 0 (pivRow A 0)

lemma first_perm : (first o A).perm = Equiv.swap 0 (pivRow A 0) := by
  simp [first, stepAt]

lemma first_mat : (first o A).mat = elim o 0 (swapped A) := rfl

lemma first_row (j : Fin (m + 1)) : (first o A).mat 0 j = swapped A 0 j := by
  rw [first_mat, elim_apply_self]

lemma first_col (i : Fin m) : (first o A).mat i.succ 0 =
    if swapped A 0 0 = 0 then swapped A i.succ 0
    else o.div (swapped A i.succ 0) (swapped A 0 0) := by
  rw [first_mat]
  split_ifs with h
  · rw [elim_of_pivot_eq_zero h]
  · exact elim_apply_mult h (Fin.succ_pos i)

lemma trail_apply (i j : Fin m) : trail o A i j =
    if swapped A 0 0 = 0 then swapped A i.succ j.succ
    else o.upd (swapped A i.succ j.succ) (o.div (swapped A i.succ 0) (swapped A 0 0))
      (swapped A 0 j.succ) := by
  simp only [trail, Matrix.submatrix_apply, first_mat]
  split_ifs with h
  · rw [elim_of_pivot_eq_zero h]
  · exact elim_apply_upd h (Fin.succ_pos i) (Fin.succ_pos j)

/-- The pivot has the largest magnitude in the first column. -/
lemma norm_swapped_le (i : Fin (m + 1)) : ‖swapped A i 0‖ ≤ ‖swapped A 0 0‖ := by
  simp only [swapped, swapRows, Matrix.submatrix_apply, id, Equiv.swap_apply_left]
  exact le_at_firstMax (fun i => ‖A i 0‖) 0 (Fin.zero_le _)

/-- A zero pivot means that the whole column is zero. -/
lemma swapped_col_eq_zero (h : swapped A 0 0 = 0) (i : Fin (m + 1)) : swapped A i 0 = 0 := by
  have := norm_swapped_le A i
  rw [h, norm_zero] at this
  exact norm_eq_zero.mp (le_antisymm this (norm_nonneg _))

end First

/-! ## The final state of a matrix of size `m + 1` -/

section Final
variable (o : Ops 𝕜) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) 𝕜)

lemma final_row (j : Fin (m + 1)) : (final o A).mat 0 j = swapped A 0 j := by
  rw [(rel_final o A).row, first_row]

lemma final_col (i : Fin m) :
    (final o A).mat i.succ 0 = (first o A).mat ((final o (trail o A)).perm i).succ 0 :=
  (rel_final o A).col i

lemma final_trail : (final o A).mat.submatrix Fin.succ Fin.succ = (final o (trail o A)).mat :=
  (rel_final o A).submatrix

lemma final_perm : (final o A).perm =
    Equiv.swap 0 (pivRow A 0) * liftPerm (final o (trail o A)).perm := by
  rw [(rel_final o A).perm, first_perm]

/-- `(P A)ᵢⱼ` in terms of `P₁ A` and the trailing permutation. -/
lemma perm_apply_zero : (final o A).perm 0 = pivRow A 0 := by
  rw [final_perm, Equiv.Perm.mul_apply, liftPerm_zero, Equiv.swap_apply_left]

lemma permMul_zero (j : Fin (m + 1)) : ((final o A).perm.permMatrix 𝕜 * A) 0 j = swapped A 0 j := by
  rw [permMatrix_mul_apply, perm_apply_zero]
  simp [swapped, swapRows]

lemma permMul_succ (i : Fin m) (j : Fin (m + 1)) : ((final o A).perm.permMatrix 𝕜 * A) i.succ j =
    swapped A ((final o (trail o A)).perm i).succ j := by
  rw [permMatrix_mul_apply, final_perm, Equiv.Perm.mul_apply, liftPerm_succ]
  rfl

end Final

end LU
