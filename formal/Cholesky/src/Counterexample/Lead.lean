import src.Algorithm
import src.Counterexample.Blocks

/-!
# Eliminating a leading block `4 I`

Let `B = [[4 I_m, 2T], [2Tᵀ, C]]`. If the arithmetic computes `√4 = 2`, `0 / 2 = 0` and
`2T / 2 = T` exactly, and leaves an entry unchanged when it subtracts a product with a zero
factor, then the first `m` steps of the in-place algorithm

* see the pivots `4`;
* store `2` on the diagonal, `0` below it in the leading block, and `T` in the coupling block;
* apply to the trailing block the rounded rank-one updates by the rows of `T`, in order:
  `C⁽ˢ⁺¹⁾ᵢⱼ = upd(C⁽ˢ⁾ᵢⱼ, Tₛᵢ, Tₛⱼ)` (`Citer`).

After these `m` steps the algorithm continues on the trailing block exactly as it would on
the matrix `C⁽ᵐ⁾` (`run_natAdd`, `pivot_natAdd`).
-/
namespace Cholesky
open Matrix

variable {m N : ℕ}

/-- The trailing block after `s` rounded rank-one updates by the rows of `T`. -/
def Citer (o : Ops) (T : Matrix (Fin m) (Fin N) ℝ) (C : Matrix (Fin N) (Fin N) ℝ) :
    ℕ → Matrix (Fin N) (Fin N) ℝ
  | 0 => C
  | s + 1 => Matrix.of fun i j =>
      if h : s < m then
        (if j ≤ i then o.upd (Citer o T C s i j) (T ⟨s, h⟩ i) (T ⟨s, h⟩ j) else Citer o T C s i j)
      else Citer o T C s i j

/-- The stored matrix. -/
def leadMatrix (T : Matrix (Fin m) (Fin N) ℝ) (C : Matrix (Fin N) (Fin N) ℝ) :
    Matrix (Fin (m + N)) (Fin (m + N)) ℝ :=
  blk ((4 : ℝ) • 1) ((2 : ℝ) • T) ((2 : ℝ) • Tᵀ) C

/-- The array after `s ≤ m` steps (only its lower triangle matters). -/
def leadState (o : Ops) (T : Matrix (Fin m) (Fin N) ℝ) (C : Matrix (Fin N) (Fin N) ℝ)
    (s : ℕ) : Matrix (Fin (m + N)) (Fin (m + N)) ℝ :=
  blk (Matrix.of fun i j => if i = j then (if (j : ℕ) < s then 2 else 4) else 0)
    ((2 : ℝ) • T) (Matrix.of fun i j => if (j : ℕ) < s then T j i else 2 * T j i)
    (Citer o T C s)

/-- The exactness properties of the arithmetic needed for the leading block. -/
structure LeadOps (o : Ops) (T : Matrix (Fin m) (Fin N) ℝ) : Prop where
  sqrt4 : o.sqrt 4 = 2
  div0 : o.div 0 2 = 0
  divT : ∀ i j, o.div (2 * T i j) 2 = T i j
  upd0 : ∀ x b c, (x = 0 ∨ x = 4 ∨ ∃ i j, x = 2 * T i j) → (b = 0 ∨ c = 0) → o.upd x b c = x

section
variable {o : Ops} {T : Matrix (Fin m) (Fin N) ℝ} {C : Matrix (Fin N) (Fin N) ℝ}

lemma lowerEq_lead0 : LowerEq (leadMatrix T C) (leadState o T C 0) := by
  intro i j _
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;>
    simp [leadMatrix, leadState, Citer, Matrix.one_apply]

section Order
variable {i j : Fin m} {i' j' : Fin N}

lemma cA_inj : Fin.castAdd N i = Fin.castAdd N j ↔ i = j :=
  ⟨fun h => Fin.ext (by simpa using congrArg Fin.val h), fun h => by rw [h]⟩
lemma nA_inj : Fin.natAdd m i' = Fin.natAdd m j' ↔ i' = j' :=
  ⟨fun h => Fin.ext (by simpa using congrArg Fin.val h), fun h => by rw [h]⟩
lemma cA_lt : Fin.castAdd N i < Fin.castAdd N j ↔ i < j := by
  simp only [Fin.lt_iff_val_lt_val, Fin.coe_castAdd]
lemma cA_le : Fin.castAdd N i ≤ Fin.castAdd N j ↔ i ≤ j := by
  simp only [Fin.le_iff_val_le_val, Fin.coe_castAdd]
lemma nA_lt : Fin.natAdd m i' < Fin.natAdd m j' ↔ i' < j' := by
  simp only [Fin.lt_iff_val_lt_val, Fin.coe_natAdd, add_lt_add_iff_left]
lemma nA_le : Fin.natAdd m i' ≤ Fin.natAdd m j' ↔ i' ≤ j' := by
  simp only [Fin.le_iff_val_le_val, Fin.coe_natAdd, add_le_add_iff_left]
lemma cA_lt_nA : Fin.castAdd N i < Fin.natAdd m j' := by
  rw [Fin.lt_iff_val_lt_val, Fin.coe_castAdd, Fin.coe_natAdd]; have := i.2; omega
lemma nA_ne_cA : Fin.natAdd m j' ≠ Fin.castAdd N i := ne_of_gt cA_lt_nA
lemma not_nA_le_cA : ¬ Fin.natAdd m j' ≤ Fin.castAdd N i := not_le.mpr cA_lt_nA
lemma not_cA_lt_cA_of_le (h : j ≤ i) : ¬ Fin.castAdd N i < Fin.castAdd N j := by
  rw [cA_lt]; exact not_lt.mpr h

end Order

/-- One step of the leading elimination. -/
theorem leadState_step (ho : LeadOps o T) (s : ℕ) (hs : s < m) :
    LowerEq (stepN o s (leadState o T C s)) (leadState o T C (s + 1)) := by
  intro i j hij
  set s' : Fin m := ⟨s, hs⟩ with hs'
  have hk : (⟨s, by omega⟩ : Fin (m + N)) = Fin.castAdd N s' := rfl
  rw [stepN, dif_pos (by omega : s < m + N), hk]
  set M := leadState o T C s with hM
  have hlt : ∀ i1 : Fin m, (i1 : ℕ) < s ↔ i1 < s' := fun i1 => Iff.rfl
  have hlt1 : ∀ i1 : Fin m, (i1 : ℕ) < s + 1 ↔ i1 ≤ s' := fun i1 => by
    rw [Fin.le_def]; simp [s']; omega
  have hkk : M (Fin.castAdd N s') (Fin.castAdd N s') = 4 := by
    simp only [M, leadState, blk_ll, Matrix.of_apply, if_true, hlt, lt_irrefl, if_false]
  have hsq : o.sqrt (M (Fin.castAdd N s') (Fin.castAdd N s')) = 2 := by rw [hkk]; exact ho.sqrt4
  have hcol_lead : ∀ i1 : Fin m, s' < i1 → newCol o M (Fin.castAdd N s') (Fin.castAdd N i1) = 0 := by
    intro i1 hi1
    have h3 : M (Fin.castAdd N i1) (Fin.castAdd N s') = 0 := by
      simp only [M, leadState, blk_ll, Matrix.of_apply]; rw [if_neg (ne_of_gt hi1)]
    rw [newCol, if_neg (by rw [cA_inj]; exact ne_of_gt hi1), if_pos (cA_lt.mpr hi1), h3, hsq,
      ho.div0]
  have hcol_trail : ∀ i2 : Fin N, newCol o M (Fin.castAdd N s') (Fin.natAdd m i2) = T s' i2 := by
    intro i2
    have h3 : M (Fin.natAdd m i2) (Fin.castAdd N s') = 2 * T s' i2 := by
      simp only [M, leadState, blk_rl, Matrix.of_apply, hlt, lt_irrefl, if_false]
    rw [newCol, if_neg nA_ne_cA, if_pos cA_lt_nA, h3, hsq, ho.divT]
  have hcol_k : newCol o M (Fin.castAdd N s') (Fin.castAdd N s') = 2 := by
    rw [newCol, if_pos rfl, hsq]
  revert hij
  refine Fin.addCases (fun i1 => ?_) (fun i2 => ?_) i <;>
    refine Fin.addCases (fun j1 => ?_) (fun j2 => ?_) j <;> intro hij
  · -- leading block
    have hji : j1 ≤ i1 := cA_le.mp hij
    simp only [stepAt, Matrix.of_apply]
    rcases lt_trichotomy j1 s' with hj | rfl | hj
    · rw [if_neg (by rw [cA_inj]; exact ne_of_lt hj),
        if_neg (fun h => absurd (cA_lt.mp h.1) (not_lt.mpr hj.le))]
      simp only [M, leadState, blk_ll, Matrix.of_apply, hlt, hlt1, hj, hj.le, if_true]
    · rw [if_pos rfl]
      rcases eq_or_lt_of_le hji with hi | hi
      · rw [← hi, hcol_k]
        simp only [leadState, blk_ll, Matrix.of_apply, if_true, hlt1, le_refl]
      · rw [hcol_lead i1 hi]
        simp only [leadState, blk_ll, Matrix.of_apply]
        rw [if_neg (ne_of_gt hi)]
    · rw [if_neg (by rw [cA_inj]; exact ne_of_gt hj), if_pos ⟨cA_lt.mpr hj, hij⟩,
        hcol_lead i1 (lt_of_lt_of_le hj hji), hcol_lead j1 hj]
      have hnl : ¬ (j1 : ℕ) < s := fun h => absurd ((hlt j1).mp h) (not_lt.mpr hj.le)
      have hnl1 : ¬ (j1 : ℕ) < s + 1 := fun h => absurd ((hlt1 j1).mp h) (not_le.mpr hj)
      have hval : M (Fin.castAdd N i1) (Fin.castAdd N j1) = if i1 = j1 then 4 else 0 := by
        simp only [M, leadState, blk_ll, Matrix.of_apply, hnl, if_false]
      rw [hval, ho.upd0 _ _ _ (by split_ifs <;> simp) (Or.inl rfl)]
      simp only [leadState, blk_ll, Matrix.of_apply, hnl1, if_false]
  · exact absurd hij not_nA_le_cA
  · -- coupling block
    simp only [stepAt, Matrix.of_apply]
    rcases lt_trichotomy j1 s' with hj | rfl | hj
    · rw [if_neg (by rw [cA_inj]; exact ne_of_lt hj),
        if_neg (fun h => absurd (cA_lt.mp h.1) (not_lt.mpr hj.le))]
      simp only [M, leadState, blk_rl, Matrix.of_apply, hlt, hlt1, hj, hj.le, if_true]
    · rw [if_pos rfl, hcol_trail i2]
      simp only [leadState, blk_rl, Matrix.of_apply, hlt1, le_refl, if_true]
    · rw [if_neg (by rw [cA_inj]; exact ne_of_gt hj), if_pos ⟨cA_lt.mpr hj, hij⟩,
        hcol_trail i2, hcol_lead j1 hj]
      have hnl : ¬ (j1 : ℕ) < s := fun h => absurd ((hlt j1).mp h) (not_lt.mpr hj.le)
      have hnl1 : ¬ (j1 : ℕ) < s + 1 := fun h => absurd ((hlt1 j1).mp h) (not_le.mpr hj)
      have hval : M (Fin.natAdd m i2) (Fin.castAdd N j1) = 2 * T j1 i2 := by
        simp only [M, leadState, blk_rl, Matrix.of_apply, hnl, if_false]
      rw [hval, ho.upd0 _ _ _ (Or.inr (Or.inr ⟨j1, i2, rfl⟩)) (Or.inr rfl)]
      simp only [leadState, blk_rl, Matrix.of_apply, hnl1, if_false]
  · -- trailing block
    have hji : j2 ≤ i2 := nA_le.mp hij
    simp only [stepAt, Matrix.of_apply]
    rw [if_neg nA_ne_cA, if_pos ⟨cA_lt_nA, hij⟩, hcol_trail i2, hcol_trail j2]
    simp only [M, leadState, blk_rr, Citer, Matrix.of_apply]
    rw [dif_pos hs, if_pos hji]

/-- After `s ≤ m` steps the array is `leadState s` (on the lower triangle). -/
theorem lead_run (ho : LeadOps o T) : ∀ s ≤ m,
    LowerEq (run o s (leadMatrix T C)) (leadState o T C s)
  | 0, _ => lowerEq_lead0
  | s + 1, hs => by
    rw [run_succ]
    intro i j hij
    rw [stepN_lowerEq o (lead_run ho s (by omega)) s i j hij]
    exact leadState_step ho s (by omega) i j hij

/-- The pivots of the leading steps are `4`. -/
theorem lead_pivot (ho : LeadOps o T) (s : ℕ) (hs : s < m) :
    pivot o (leadMatrix T C) ⟨s, by omega⟩ = 4 := by
  unfold pivot
  rw [lead_run ho s hs.le _ _ le_rfl]
  have : (⟨s, by omega⟩ : Fin (m + N)) = Fin.castAdd N ⟨s, hs⟩ := rfl
  rw [this]
  simp only [leadState, blk_ll, Matrix.of_apply, if_true, lt_irrefl, if_false]

end

/-! ## Continuing on the trailing block -/

lemma stepN_natAdd (o : Ops) (B : Matrix (Fin (m + N)) (Fin (m + N)) ℝ) (s : ℕ) :
    (stepN o (m + s) B).submatrix (Fin.natAdd m) (Fin.natAdd m) =
      stepN o s (B.submatrix (Fin.natAdd m) (Fin.natAdd m)) := by
  unfold stepN
  by_cases hs : s < N
  · rw [dif_pos (by omega), dif_pos hs]
    have hk : (⟨m + s, by omega⟩ : Fin (m + N)) = Fin.natAdd m ⟨s, hs⟩ := rfl
    rw [hk]
    ext i j
    have e1 : ∀ a b : Fin N, Fin.natAdd m a = Fin.natAdd m b ↔ a = b := fun a b =>
      ⟨fun h => Fin.ext (by have := congrArg Fin.val h; simp at this; exact this),
        fun h => by rw [h]⟩
    have e2 : ∀ a b : Fin N, Fin.natAdd m a < Fin.natAdd m b ↔ a < b := fun _ _ => nA_lt
    have e3 : ∀ a b : Fin N, Fin.natAdd m a ≤ Fin.natAdd m b ↔ a ≤ b := fun _ _ => nA_le
    simp only [stepAt, newCol, Matrix.submatrix_apply, Matrix.of_apply, e1, e2, e3]
  · rw [dif_neg (by omega), dif_neg hs]

lemma runFrom_natAdd (o : Ops) (B : Matrix (Fin (m + N)) (Fin (m + N)) ℝ) (s t : ℕ) :
    (runFrom o (m + s) t B).submatrix (Fin.natAdd m) (Fin.natAdd m) =
      runFrom o s t (B.submatrix (Fin.natAdd m) (Fin.natAdd m)) := by
  induction t generalizing s B with
  | zero => rfl
  | succ t ih =>
    rw [runFrom, runFrom, show m + s + 1 = m + (s + 1) by ring, ih, stepN_natAdd]

/-- After the first `m` steps, the algorithm runs on the trailing block. -/
theorem run_natAdd (o : Ops) (B : Matrix (Fin (m + N)) (Fin (m + N)) ℝ) (t : ℕ) :
    (run o (m + t) B).submatrix (Fin.natAdd m) (Fin.natAdd m) =
      run o t ((run o m B).submatrix (Fin.natAdd m) (Fin.natAdd m)) := by
  rw [run, runFrom_add, zero_add]
  exact runFrom_natAdd o (run o m B) 0 t

theorem pivot_natAdd (o : Ops) (B : Matrix (Fin (m + N)) (Fin (m + N)) ℝ) (i : Fin N) :
    pivot o B (Fin.natAdd m i) =
      pivot o ((run o m B).submatrix (Fin.natAdd m) (Fin.natAdd m)) i := by
  unfold pivot
  have := congrFun (congrFun (run_natAdd o B i) i) i
  simp only [Matrix.submatrix_apply] at this
  rw [← this]
  rfl

/-- Success of the whole run implies success on the trailing block. -/
theorem succeeds_trailing (o : Ops) (B : Matrix (Fin (m + N)) (Fin (m + N)) ℝ)
    (h : Succeeds o B) : Succeeds o ((run o m B).submatrix (Fin.natAdd m) (Fin.natAdd m)) := by
  intro i
  rw [← pivot_natAdd]
  exact h _

lemma lowerEq_submatrix_natAdd {A B : Matrix (Fin (m + N)) (Fin (m + N)) ℝ} (h : LowerEq A B) :
    LowerEq (A.submatrix (Fin.natAdd m) (Fin.natAdd m)) (B.submatrix (Fin.natAdd m) (Fin.natAdd m)) :=
  fun _ _ hij => h _ _ (nA_le.mpr hij)

/-- For the stored matrix `[[4I, 2T], [2Tᵀ, C]]`: if the run succeeds, the run on `C⁽ᵐ⁾`
succeeds. -/
theorem lead_succeeds_Citer {o : Ops} {T : Matrix (Fin m) (Fin N) ℝ}
    {C : Matrix (Fin N) (Fin N) ℝ} (ho : LeadOps o T) (h : Succeeds o (leadMatrix T C)) :
    Succeeds o (Citer o T C m) := by
  have h1 := succeeds_trailing o _ h
  have h2 : LowerEq ((run o m (leadMatrix T C)).submatrix (Fin.natAdd m) (Fin.natAdd m))
      (Citer o T C m) := by
    intro i j hij
    have := lowerEq_submatrix_natAdd (lead_run (C := C) ho m le_rfl) i j hij
    rw [this]
    simp [leadState]
  exact (succeeds_lowerEq o h2).mp h1

end Cholesky
