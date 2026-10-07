import src.Algorithm

/-!
# Operation counts

Step `k` (counting from 0) of the routine on an `n × n` array

* compares the `n - k` candidates of column `k` (`n - 1 - k` comparisons, `card_candidates`);
* interchanges at most two rows of `A` and two rows of `P` (`2n` entries each);
* computes `n - 1 - k` multipliers and updates `(n - 1 - k)²` entries, one multiplication and
  one subtraction each (`card_multipliers`, `card_updates`). Elimination changes no other
  entry (`elim_apply_of_not_mem`).

Hence the arithmetic costs `luFlops n = 2n³/3 - n²/2 - n/6` (`luFlops_real`), with leading
term `2n³/3` (`luFlops_ratio`), while the pivot searches and row interchanges cost
`n(n - 1)/2 + 4n²`, which is `O(n²)` (`overhead_le`). Forming the residual `b - A x̂`
costs `2n²` flops (`residualFlops`).
-/
namespace LU
open Finset Filter Topology

variable {𝕜 : Type*} [NormedField 𝕜] [DecidableEq 𝕜]

/-- Elimination at step `k` changes only the multipliers `(i, k)` and the trailing block
`(i, j)` with `i, j > k`. -/
lemma elim_apply_of_not_mem {n : ℕ} (o : Ops 𝕜) (k : Fin n) (M : Matrix (Fin n) (Fin n) 𝕜)
    {i j : Fin n} (h : ¬ (k < i ∧ k ≤ j)) : elim o k M i j = M i j := by
  by_cases hi : k < i
  · have hj : j < k := lt_of_not_ge fun hj => h ⟨hi, hj⟩
    exact elim_apply_of_lt_col hj
  · exact elim_apply_of_not_lt hi

lemma card_filter_lt_fin (n k : ℕ) :
    (univ.filter (fun i : Fin n => k < i.val)).card = n - 1 - k := by
  have : (univ.filter (fun i : Fin n => k < i.val)).map Fin.valEmbedding = Ioo k n := by
    ext x
    simp only [Finset.mem_map, Finset.mem_filter, Finset.mem_univ, true_and,
      Fin.valEmbedding_apply, Finset.mem_Ioo]
    constructor
    · rintro ⟨i, hi, rfl⟩; exact ⟨hi, i.2⟩
    · rintro ⟨h1, h2⟩; exact ⟨⟨x, h2⟩, h1, rfl⟩
  rw [← card_map, this, Nat.card_Ioo]
  omega

lemma card_filter_le_fin (n k : ℕ) :
    (univ.filter (fun i : Fin n => k ≤ i.val)).card = n - k := by
  have : (univ.filter (fun i : Fin n => k ≤ i.val)).map Fin.valEmbedding = Ico k n := by
    ext x
    simp only [Finset.mem_map, Finset.mem_filter, Finset.mem_univ, true_and,
      Fin.valEmbedding_apply, Finset.mem_Ico]
    constructor
    · rintro ⟨i, hi, rfl⟩; exact ⟨hi, i.2⟩
    · rintro ⟨h1, h2⟩; exact ⟨⟨x, h2⟩, h1, rfl⟩
  rw [← card_map, this, Nat.card_Ico]

/-- The pivot search at step `k` examines `n - k` candidates. -/
theorem card_candidates (n k : ℕ) :
    (univ.filter (fun i : Fin n => k ≤ i.val)).card = n - k := card_filter_le_fin n k

/-- Step `k` computes `n - 1 - k` multipliers. -/
theorem card_multipliers (n k : ℕ) :
    (univ.filter (fun i : Fin n => k < i.val)).card = n - 1 - k := card_filter_lt_fin n k

/-- Step `k` updates `(n - 1 - k)²` entries. -/
theorem card_updates (n k : ℕ) :
    (univ.filter (fun ij : Fin n × Fin n => k < ij.1.val ∧ k < ij.2.val)).card =
      (n - 1 - k) ^ 2 := by
  have : univ.filter (fun ij : Fin n × Fin n => k < ij.1.val ∧ k < ij.2.val) =
      (univ.filter (fun i : Fin n => k < i.val)) ×ˢ (univ.filter (fun i : Fin n => k < i.val)) := by
    ext ⟨i, j⟩; simp
  rw [this, card_product, card_filter_lt_fin, sq]

/-- Arithmetic of the factorization: at step `k`, `n - 1 - k` divisions and `(n - 1 - k)²`
updates of two flops each. -/
def luFlops (n : ℕ) : ℕ := ∑ k ∈ range n, ((n - 1 - k) + 2 * (n - 1 - k) ^ 2)

theorem luFlops_real (n : ℕ) :
    (luFlops n : ℝ) = 2 * (n : ℝ) ^ 3 / 3 - (n : ℝ) ^ 2 / 2 - n / 6 := by
  unfold luFlops
  rw [← Finset.sum_range_reflect]
  have : ∀ m ∈ range n, (n - 1 - (n - 1 - m)) + 2 * (n - 1 - (n - 1 - m)) ^ 2 =
      m + 2 * m ^ 2 := by
    intro m hm; rw [Finset.mem_range] at hm
    rw [show n - 1 - (n - 1 - m) = m by omega]
  rw [Finset.sum_congr rfl this]
  clear this
  push_cast
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; push_cast; ring

/-- **The leading cost is `2n³/3` flops.** -/
theorem luFlops_ratio :
    Tendsto (fun n : ℕ => (luFlops n : ℝ) / (n : ℝ) ^ 3) atTop (𝓝 (2 / 3)) := by
  have h1 : Tendsto (fun n : ℕ => (1 : ℝ) / n) atTop (𝓝 0) :=
    tendsto_const_div_atTop_nhds_zero_nat 1
  have h : Tendsto (fun n : ℕ => 2 / 3 - 1 / 2 * ((1 : ℝ) / n) - 1 / 6 * ((1 : ℝ) / n) ^ 2)
      atTop (𝓝 (2 / 3 - 1 / 2 * 0 - 1 / 6 * 0 ^ 2)) :=
    ((tendsto_const_nhds.sub (h1.const_mul _)).sub ((h1.pow 2).const_mul _))
  simp only [mul_zero, sub_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow] at h
  refine h.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with n hn
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  rw [luFlops_real]
  field_simp

/-- Comparisons in the pivot searches: `∑_k (n - 1 - k) = n(n - 1)/2`. -/
def pivotCompares (n : ℕ) : ℕ := ∑ k ∈ range n, (n - 1 - k)

theorem two_mul_pivotCompares (n : ℕ) : 2 * pivotCompares n = n * (n - 1) := by
  unfold pivotCompares
  rw [← Finset.sum_range_reflect]
  have : ∀ m ∈ range n, n - 1 - (n - 1 - m) = m := by
    intro m hm; rw [Finset.mem_range] at hm; omega
  rw [Finset.sum_congr rfl this, mul_comm, Finset.sum_range_id_mul_two]

/-- Entries moved by the interchanges: at most two rows of `A` and two of `P` per step. -/
def swapMoves (n : ℕ) : ℕ := ∑ _k ∈ range n, 4 * n

/-- **Pivoting adds `O(n²)` work.** -/
theorem overhead_le (n : ℕ) : pivotCompares n + swapMoves n ≤ 5 * n ^ 2 := by
  have h := two_mul_pivotCompares n
  have h2 : swapMoves n = 4 * n ^ 2 := by simp [swapMoves, sq]; ring
  have h3 : pivotCompares n ≤ n ^ 2 := by
    have : n * (n - 1) ≤ n * n := Nat.mul_le_mul_left n (Nat.sub_le n 1)
    nlinarith
  omega

/-- The residual `b - A x̂`: per row, `n` multiplications, `n - 1` additions, one
subtraction. -/
def residualFlops (n : ℕ) : ℕ := n * (n + (n - 1) + 1)

theorem residualFlops_eq (n : ℕ) : residualFlops n = 2 * n ^ 2 := by
  unfold residualFlops
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · rfl
  · rw [show n + (n - 1) + 1 = 2 * n by omega]; ring

end LU
