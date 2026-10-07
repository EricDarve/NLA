import Mathlib

/-!
# Choosing the pivot row

Partial pivoting picks, at step `k`, a row `p ≥ k` whose entry in column `k` has the
largest magnitude. The routine of the notes computes `p = k + np.argmax(np.abs(A[k:, k]))`,
and NumPy returns the first maximal index. `IsFirstMax c k p` says that `p` is this index
for the magnitudes `c`; it exists and is unique (`exists_firstMax`, `IsFirstMax.unique`),
and `firstMax c k` denotes it.
-/
namespace LU

variable {n m : ℕ}

/-- `p` is the first index `≥ k` at which `c` attains its maximum over `{k, …, n - 1}`. -/
def IsFirstMax (c : Fin n → ℝ) (k p : Fin n) : Prop :=
  k ≤ p ∧ (∀ i, k ≤ i → c i ≤ c p) ∧ ∀ i, k ≤ i → i < p → c i < c p

lemma exists_firstMax (c : Fin n → ℝ) (k : Fin n) : ∃ p, IsFirstMax c k p := by
  classical
  set S := Finset.univ.filter (fun i : Fin n => k ≤ i) with hS
  have hSne : S.Nonempty := ⟨k, by simp [hS]⟩
  obtain ⟨q, hqS, hq⟩ := S.exists_max_image c hSne
  set T := S.filter (fun i => c i = c q) with hT
  have hTne : T.Nonempty := ⟨q, by simp [hT, hqS]⟩
  have hpT := T.min'_mem hTne
  simp only [hT, hS, Finset.mem_filter, Finset.mem_univ, true_and] at hpT
  refine ⟨T.min' hTne, hpT.1, fun i hi => ?_, fun i hi hip => ?_⟩
  · rw [hpT.2]
    exact hq i (by simp [hS, hi])
  · rw [hpT.2]
    rcases lt_or_eq_of_le (hq i (by simp [hS, hi])) with h | h
    · exact h
    · exfalso
      have : T.min' hTne ≤ i := T.min'_le i (by simp [hT, hS, hi, h])
      exact absurd hip (not_lt.mpr this)

lemma IsFirstMax.unique {c : Fin n → ℝ} {k p q : Fin n} (hp : IsFirstMax c k p)
    (hq : IsFirstMax c k q) : p = q := by
  rcases lt_trichotomy p q with h | h | h
  · exact absurd (hp.2.1 q hq.1) (not_le.mpr (hq.2.2 p hp.1 h))
  · exact h
  · exact absurd (hq.2.1 p hp.1) (not_le.mpr (hp.2.2 q hq.1 h))

/-- The first maximal index `≥ k`. -/
noncomputable def firstMax (c : Fin n → ℝ) (k : Fin n) : Fin n := (exists_firstMax c k).choose

lemma firstMax_spec (c : Fin n → ℝ) (k : Fin n) : IsFirstMax c k (firstMax c k) :=
  (exists_firstMax c k).choose_spec

lemma firstMax_eq {c : Fin n → ℝ} {k p : Fin n} (h : IsFirstMax c k p) : firstMax c k = p :=
  (firstMax_spec c k).unique h

lemma le_firstMax (c : Fin n → ℝ) (k : Fin n) : k ≤ firstMax c k := (firstMax_spec c k).1

lemma le_at_firstMax (c : Fin n → ℝ) (k : Fin n) {i : Fin n} (hi : k ≤ i) :
    c i ≤ c (firstMax c k) := (firstMax_spec c k).2.1 i hi

/-- The choice depends only on the candidates `c i`, `i ≥ k`. -/
lemma firstMax_congr {c c' : Fin n → ℝ} {k : Fin n} (h : ∀ i, k ≤ i → c i = c' i) :
    firstMax c k = firstMax c' k := by
  apply firstMax_eq
  obtain ⟨h1, h2, h3⟩ := firstMax_spec c' k
  refine ⟨h1, fun i hi => ?_, fun i hi hip => ?_⟩
  · rw [h i hi, h _ h1]; exact h2 i hi
  · rw [h i hi, h _ h1]; exact h3 i hi hip

/-- If no later candidate is larger, the current row is chosen (ties go to the first). -/
lemma firstMax_self {c : Fin n → ℝ} {k : Fin n} (h : ∀ i, k ≤ i → c i ≤ c k) :
    firstMax c k = k :=
  firstMax_eq ⟨le_rfl, h, fun _ hi hik => absurd hik (not_lt.mpr hi)⟩

/-- The last index has no competitors. -/
lemma firstMax_last (c : Fin (m + 1) → ℝ) : firstMax c (Fin.last m) = Fin.last m :=
  firstMax_self fun i hi => by rw [le_antisymm (Fin.le_last i) hi]

lemma IsFirstMax.succ {c : Fin (m + 1) → ℝ} {k p : Fin m}
    (h : IsFirstMax (fun i => c i.succ) k p) : IsFirstMax c k.succ p.succ := by
  obtain ⟨h1, h2, h3⟩ := h
  refine ⟨Fin.succ_le_succ_iff.mpr h1, fun i hi => ?_, fun i hi hip => ?_⟩
  · refine Fin.cases (fun hi => ?_) (fun i hi => ?_) i hi
    · exact absurd hi (not_le.mpr (Fin.succ_pos k))
    · exact h2 i (Fin.succ_le_succ_iff.mp hi)
  · refine Fin.cases (fun hi _ => ?_) (fun i hi hip => ?_) i hi hip
    · exact absurd hi (not_le.mpr (Fin.succ_pos k))
    · exact h3 i (Fin.succ_le_succ_iff.mp hi) (Fin.succ_lt_succ_iff.mp hip)

/-- Searching the rows `≥ k + 1` of a matrix of size `m + 1` is searching the rows `≥ k` of
its trailing part. -/
lemma firstMax_succ (c : Fin (m + 1) → ℝ) (k : Fin m) :
    firstMax c k.succ = (firstMax (fun i => c i.succ) k).succ :=
  firstMax_eq (firstMax_spec _ _).succ

end LU
