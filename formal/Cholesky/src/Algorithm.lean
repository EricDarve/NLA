import src.BlockStep

/-!
# The in-place Cholesky algorithm

This file models the routine `cholesky_in_place` of the notes. For
`k = 0, …, n - 1` (the notes count from 1), step `k`

1. replaces `a_kk` by `sqrt a_kk`,
2. divides the entries below it by the new diagonal entry, and
3. updates the lower triangle of the trailing block,
   `a_ij ← a_ij - a_ik a_jk` for `k < j ≤ i`.

The arithmetic is a parameter `o : Ops`: exact real arithmetic, or rounded
arithmetic (`Ops.rounded`, one rounding for the product and one for the
subtraction, as in NumPy; `Ops.fused`, a fused multiply-add). The routine stops
when it sees a nonpositive pivot; `Succeeds` says that this never happens.

We prove the structural facts used later: the algorithm reads and writes only the
lower triangle, and after the first step it acts on the trailing block exactly
as the same algorithm acts on a matrix of size `n`.
-/
namespace Cholesky
open Matrix

/-- The elementary operations of the algorithm. `upd a b c` computes `a - b * c`. -/
structure Ops where
  sqrt : ℝ → ℝ
  div : ℝ → ℝ → ℝ
  upd : ℝ → ℝ → ℝ → ℝ

namespace Ops

/-- Exact real arithmetic. -/
noncomputable def exact : Ops :=
  ⟨Real.sqrt, fun a b => a / b, fun a b c => a - b * c⟩

/-- Rounded arithmetic with a rounding function `r`: the product and the
subtraction in an update are rounded separately. -/
noncomputable def rounded (r : ℝ → ℝ) : Ops :=
  ⟨fun a => r (Real.sqrt a), fun a b => r (a / b), fun a b c => r (a - r (b * c))⟩

/-- Rounded arithmetic in which each update is a fused multiply-add. -/
noncomputable def fused (r : ℝ → ℝ) : Ops :=
  ⟨fun a => r (Real.sqrt a), fun a b => r (a / b), fun a b c => r (a - b * c)⟩

/-- The standard model of floating-point arithmetic with unit roundoff `u`.
For updates it allows one relative error in the product and one in the
subtraction (a fused multiply-add is the case `ε = 0`). -/
structure StdModel (o : Ops) (u : ℝ) : Prop where
  sqrt : ∀ a, ∃ δ, |δ| ≤ u ∧ o.sqrt a = Real.sqrt a * (1 + δ)
  div : ∀ a b, ∃ δ, |δ| ≤ u ∧ o.div a b = a / b * (1 + δ)
  upd : ∀ a b c, ∃ δ ε, |δ| ≤ u ∧ |ε| ≤ u ∧ o.upd a b c = (a - b * c * (1 + ε)) * (1 + δ)

lemma exists_rel {u : ℝ} {r : ℝ → ℝ} (hr : ∀ x, |r x - x| ≤ u * |x|) (x : ℝ) :
    ∃ δ, |δ| ≤ u ∧ r x = x * (1 + δ) := by
  have hu : 0 ≤ u := by
    have := hr 1
    simp only [abs_one, mul_one] at this
    exact le_trans (abs_nonneg _) this
  by_cases hx : x = 0
  · refine ⟨0, by simpa using hu, ?_⟩
    have := hr 0
    simp only [hx, sub_zero, abs_zero, mul_zero] at this ⊢
    simpa using abs_nonpos_iff.mp this
  · refine ⟨(r x - x) / x, ?_, ?_⟩
    · rw [abs_div, div_le_iff₀ (abs_pos.mpr hx)]
      exact hr x
    · field_simp
      ring

theorem stdModel_exact {u : ℝ} (hu : 0 ≤ u) : exact.StdModel u where
  sqrt a := ⟨0, by simpa using hu, by simp [exact]⟩
  div a b := ⟨0, by simpa using hu, by simp [exact]⟩
  upd a b c := ⟨0, 0, by simpa using hu, by simpa using hu, by simp [exact]⟩

theorem stdModel_rounded {u : ℝ} {r : ℝ → ℝ} (hr : ∀ x, |r x - x| ≤ u * |x|) :
    (rounded r).StdModel u where
  sqrt a := exists_rel hr _
  div a b := exists_rel hr _
  upd a b c := by
    obtain ⟨ε, hε, he⟩ := exists_rel hr (b * c)
    obtain ⟨δ, hδ, hd⟩ := exists_rel hr (a - r (b * c))
    exact ⟨δ, ε, hδ, hε, by simp only [rounded]; rw [hd, he]⟩

theorem stdModel_fused {u : ℝ} {r : ℝ → ℝ} (hr : ∀ x, |r x - x| ≤ u * |x|) :
    (fused r).StdModel u where
  sqrt a := exists_rel hr _
  div a b := exists_rel hr _
  upd a b c := by
    have hu : 0 ≤ u := by
      have := hr 1
      simp only [abs_one, mul_one] at this
      exact le_trans (abs_nonneg _) this
    obtain ⟨δ, hδ, hd⟩ := exists_rel hr (a - b * c)
    exact ⟨δ, 0, hδ, by simpa using hu, by simp only [fused]; rw [hd]; ring⟩

end Ops

variable {n : ℕ}

/-- The lower triangle of a matrix (the stored factor). -/
def lowerPart (M : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j => if j ≤ i then M i j else 0

/-- The symmetric matrix defined by the lower triangle of `M`. -/
def symmOfLower (M : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j => if j ≤ i then M i j else M j i

/-- Two arrays with the same lower triangle. -/
def LowerEq (A B : Matrix (Fin n) (Fin n) ℝ) : Prop :=
  ∀ i j, j ≤ i → A i j = B i j

/-! ## One step and the full loop -/

/-- Column `k` after the square root and the divisions of step `k`. -/
def newCol (o : Ops) (A : Matrix (Fin n) (Fin n) ℝ) (k i : Fin n) : ℝ :=
  if i = k then o.sqrt (A k k) else if k < i then o.div (A i k) (o.sqrt (A k k)) else A i k

/-- Step `k` of the in-place algorithm. -/
def stepAt (o : Ops) (k : Fin n) (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.of fun i j =>
    if j = k then newCol o A k i
    else if k < j ∧ j ≤ i then o.upd (A i j) (newCol o A k i) (newCol o A k j)
    else A i j

/-- Step `s`, numbered by a natural number (a no-op when `s ≥ n`). -/
def stepN (o : Ops) (s : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  if h : s < n then stepAt o ⟨s, h⟩ A else A

/-- Steps `s, s + 1, …, s + m - 1`. -/
def runFrom (o : Ops) : ℕ → ℕ → Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin n) ℝ
  | _, 0, A => A
  | s, m + 1, A => runFrom o (s + 1) m (stepN o s A)

/-- The array after the first `m` steps. -/
def run (o : Ops) (m : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  runFrom o 0 m A

/-- The pivot examined at step `k`. -/
def pivot (o : Ops) (A : Matrix (Fin n) (Fin n) ℝ) (k : Fin n) : ℝ :=
  run o k A k k

/-- The routine completes: every pivot it examines is positive. -/
def Succeeds (o : Ops) (A : Matrix (Fin n) (Fin n) ℝ) : Prop :=
  ∀ k : Fin n, 0 < pivot o A k

/-- The computed factor: the lower triangle of the final array. -/
def factor (o : Ops) (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  lowerPart (run o n A)

/-! ## Basic facts -/

lemma lowerPart_lower (M : Matrix (Fin n) (Fin n) ℝ) : IsLowerTriangular (lowerPart M) := by
  intro i j hij
  simp [lowerPart, not_le.mpr hij]

lemma lowerEq_symmOfLower (A : Matrix (Fin n) (Fin n) ℝ) : LowerEq A (symmOfLower A) := by
  intro i j hij
  simp [symmOfLower, hij]

lemma symmOfLower_transpose (A : Matrix (Fin n) (Fin n) ℝ) :
    (symmOfLower A)ᵀ = symmOfLower A := by
  ext i j
  simp only [symmOfLower, Matrix.transpose_apply, Matrix.of_apply]
  rcases lt_trichotomy i j with h | h | h
  · simp [h.le, not_le.mpr h]
  · subst h; simp
  · simp [h.le, not_le.mpr h]

lemma symmOfLower_eq_self {A : Matrix (Fin n) (Fin n) ℝ} (hA : Aᵀ = A) : symmOfLower A = A := by
  ext i j
  simp only [symmOfLower, Matrix.of_apply]
  split_ifs
  · rfl
  · exact symm_apply_of_transpose hA j i

lemma runFrom_succ' (o : Ops) (s m : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) :
    runFrom o s (m + 1) A = stepN o (s + m) (runFrom o s m A) := by
  induction m generalizing s A with
  | zero => simp [runFrom]
  | succ m ih =>
    rw [runFrom, ih, runFrom]
    congr 1
    omega

lemma run_succ (o : Ops) (m : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) :
    run o (m + 1) A = stepN o m (run o m A) := by
  simp [run, runFrom_succ']

lemma runFrom_add (o : Ops) (s a b : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) :
    runFrom o s (a + b) A = runFrom o (s + a) b (runFrom o s a A) := by
  induction a generalizing s A with
  | zero => simp [runFrom]
  | succ a ih =>
    rw [show a + 1 + b = (a + b) + 1 by omega, runFrom, ih, runFrom]
    congr 1
    omega

/-! ## Only the lower triangle is read and written -/

lemma newCol_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B)
    (k i : Fin n) (hki : k ≤ i) : newCol o A k i = newCol o B k i := by
  unfold newCol
  rw [h k k le_rfl]
  split_ifs with h1 h2
  · rfl
  · rw [h i k hki]
  · exact absurd (lt_of_le_of_ne hki (Ne.symm h1)) h2

lemma stepAt_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B)
    (k : Fin n) : LowerEq (stepAt o k A) (stepAt o k B) := by
  intro i j hij
  simp only [stepAt, Matrix.of_apply]
  split_ifs with h1 h2
  · subst h1; exact newCol_lowerEq o h _ _ hij
  · rw [h i j hij, newCol_lowerEq o h k i (le_of_lt (lt_of_lt_of_le h2.1 h2.2)),
      newCol_lowerEq o h k j (le_of_lt h2.1)]
  · exact h i j hij

lemma stepN_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B) (s : ℕ) :
    LowerEq (stepN o s A) (stepN o s B) := by
  unfold stepN
  split_ifs
  · exact stepAt_lowerEq o h _
  · exact h

lemma runFrom_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B)
    (s m : ℕ) : LowerEq (runFrom o s m A) (runFrom o s m B) := by
  induction m generalizing s A B with
  | zero => exact h
  | succ m ih => exact ih (stepN_lowerEq o h s) (s + 1)

lemma pivot_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B) (k : Fin n) :
    pivot o A k = pivot o B k :=
  runFrom_lowerEq o h 0 k k k le_rfl

lemma succeeds_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B) :
    Succeeds o A ↔ Succeeds o B := by
  simp only [Succeeds, pivot_lowerEq o h]

lemma factor_lowerEq (o : Ops) {A B : Matrix (Fin n) (Fin n) ℝ} (h : LowerEq A B) :
    factor o A = factor o B := by
  ext i j
  simp only [factor, lowerPart, Matrix.of_apply]
  split_ifs with hij
  · exact runFrom_lowerEq o h 0 n i j hij
  · rfl

/-- The algorithm interprets its input as the symmetric matrix defined by the lower
triangle. -/
lemma factor_symmOfLower (o : Ops) (A : Matrix (Fin n) (Fin n) ℝ) :
    factor o (symmOfLower A) = factor o A :=
  (factor_lowerEq o (lowerEq_symmOfLower A)).symm

lemma succeeds_symmOfLower (o : Ops) (A : Matrix (Fin n) (Fin n) ℝ) :
    Succeeds o (symmOfLower A) ↔ Succeeds o A :=
  (succeeds_lowerEq o (lowerEq_symmOfLower A)).symm

lemma stepAt_upper (o : Ops) (k : Fin n) (A : Matrix (Fin n) (Fin n) ℝ) (i j : Fin n)
    (hij : i < j) : stepAt o k A i j = A i j := by
  simp only [stepAt, Matrix.of_apply]
  split_ifs with h1 h2
  · subst h1
    simp [newCol, ne_of_lt hij, not_lt.mpr hij.le]
  · exact absurd h2.2 (not_le.mpr hij)
  · rfl

/-- The strict upper triangle is never changed. -/
lemma run_upper (o : Ops) (m : ℕ) (A : Matrix (Fin n) (Fin n) ℝ) (i j : Fin n)
    (hij : i < j) : run o m A i j = A i j := by
  induction m with
  | zero => rfl
  | succ m ih =>
    rw [run_succ]
    unfold stepN
    split_ifs
    · rw [stepAt_upper o _ _ i j hij, ih]
    · exact ih

/-! ## Reduction to the trailing block -/

section Shift
variable {m : ℕ}

lemma newCol_succ (o : Ops) (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (k i : Fin m) :
    newCol o B k.succ i.succ = newCol o (B.submatrix Fin.succ Fin.succ) k i := by
  simp [newCol, Fin.succ_inj, Fin.succ_lt_succ_iff]

lemma stepAt_succ (o : Ops) (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (k : Fin m) :
    (stepAt o k.succ B).submatrix Fin.succ Fin.succ =
      stepAt o k (B.submatrix Fin.succ Fin.succ) := by
  ext i j
  simp only [stepAt, Matrix.submatrix_apply, Matrix.of_apply, Fin.succ_inj,
    Fin.succ_lt_succ_iff, Fin.succ_le_succ_iff, newCol_succ]

lemma stepN_succ (o : Ops) (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (s : ℕ) :
    (stepN o (s + 1) B).submatrix Fin.succ Fin.succ =
      stepN o s (B.submatrix Fin.succ Fin.succ) := by
  unfold stepN
  by_cases hs : s < m
  · have hs' : s + 1 < m + 1 := by omega
    rw [dif_pos hs', dif_pos hs]
    have : (⟨s + 1, hs'⟩ : Fin (m + 1)) = (⟨s, hs⟩ : Fin m).succ := rfl
    rw [this, stepAt_succ]
  · have hs' : ¬ s + 1 < m + 1 := by omega
    rw [dif_neg hs', dif_neg hs]

lemma stepN_succ_col0 (o : Ops) (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (s : ℕ)
    (i : Fin (m + 1)) : stepN o (s + 1) B i 0 = B i 0 := by
  unfold stepN
  split_ifs with hs
  · simp only [stepAt, Matrix.of_apply]
    have h1 : (0 : Fin (m + 1)) ≠ ⟨s + 1, hs⟩ := by
      intro h; have := congrArg Fin.val h; simp at this
    have h2 : ¬ ((⟨s + 1, hs⟩ : Fin (m + 1)) < 0) := by simp [Fin.lt_def]
    simp [h1, h2]
  · rfl

lemma runFrom_succ_submatrix (o : Ops) (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ)
    (s t : ℕ) :
    (runFrom o (s + 1) t B).submatrix Fin.succ Fin.succ =
      runFrom o s t (B.submatrix Fin.succ Fin.succ) := by
  induction t generalizing s B with
  | zero => rfl
  | succ t ih =>
    rw [runFrom, runFrom, ← stepN_succ, ih]

lemma runFrom_succ_col0 (o : Ops) (B : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ)
    (s t : ℕ) (i : Fin (m + 1)) : runFrom o (s + 1) t B i 0 = B i 0 := by
  induction t generalizing s B with
  | zero => rfl
  | succ t ih => rw [runFrom, ih, stepN_succ_col0]

/-- The array after the first step, restricted to the trailing block. -/
def trail (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) : Matrix (Fin m) (Fin m) ℝ :=
  (stepN o 0 A).submatrix Fin.succ Fin.succ

/-- The computed first column below the diagonal. -/
def firstCol (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (i : Fin m) : ℝ :=
  o.div (A i.succ 0) (o.sqrt (A 0 0))

lemma stepN_zero_eq (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) :
    stepN o 0 A = stepAt o 0 A := by
  simp [stepN]

lemma trail_apply (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (i j : Fin m)
    (hij : j ≤ i) :
    trail o A i j = o.upd (A i.succ j.succ) (firstCol o A i) (firstCol o A j) := by
  simp [trail, stepN_zero_eq, stepAt, newCol, firstCol, Fin.succ_ne_zero, hij,
    Fin.succ_pos]

lemma run_succ_submatrix (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (t : ℕ) :
    (run o (t + 1) A).submatrix Fin.succ Fin.succ = run o t (trail o A) := by
  simp only [run, runFrom]
  rw [runFrom_succ_submatrix]
  rfl

lemma run_succ_col0 (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (t : ℕ)
    (i : Fin (m + 1)) : run o (t + 1) A i 0 = stepAt o 0 A i 0 := by
  simp only [run, runFrom]
  rw [runFrom_succ_col0, stepN_zero_eq]

lemma pivot_zero (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) :
    pivot o A 0 = A 0 0 := rfl

lemma pivot_succ (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (k : Fin m) :
    pivot o A k.succ = pivot o (trail o A) k := by
  unfold pivot
  rw [Fin.val_succ, ← run_succ_submatrix]
  rfl

lemma succeeds_succ_iff (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) :
    Succeeds o A ↔ 0 < A 0 0 ∧ Succeeds o (trail o A) := by
  constructor
  · intro h
    exact ⟨h 0, fun k => by rw [← pivot_succ]; exact h k.succ⟩
  · rintro ⟨h0, h⟩ k
    refine Fin.cases ?_ (fun k => ?_) k
    · exact h0
    · rw [pivot_succ]; exact h k

/-- The computed factor of a matrix of size `m + 1`: the first column, and the factor
computed for the trailing block. -/
lemma factor_zero_zero (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) :
    factor o A 0 0 = o.sqrt (A 0 0) := by
  simp [factor, lowerPart, run_succ_col0, stepAt, newCol]

lemma factor_succ_zero (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (i : Fin m) :
    factor o A i.succ 0 = firstCol o A i := by
  simp [factor, lowerPart, run_succ_col0, stepAt, newCol, firstCol, Fin.succ_ne_zero,
    Fin.succ_pos]

lemma factor_zero_succ (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (j : Fin m) :
    factor o A 0 j.succ = 0 := by
  simp only [factor, lowerPart, Matrix.of_apply]
  rw [if_neg (not_le.mpr (Fin.succ_pos j))]

lemma factor_succ_succ (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (i j : Fin m) :
    factor o A i.succ j.succ = factor o (trail o A) i j := by
  have h := congrFun (congrFun (run_succ_submatrix o A m) i) j
  simp only [Matrix.submatrix_apply] at h
  simp only [factor, lowerPart, Matrix.of_apply, Fin.succ_le_succ_iff, h]

lemma trailing_factor (o : Ops) (A : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) :
    trailing (factor o A) = factor o (trail o A) := by
  ext i j
  simp [trailing, factor_succ_succ]

end Shift

end Cholesky
