import src.Algorithm

/-!
# Arithmetic and storage

Operation counts for the in-place algorithm, read off from its loops.

* Step `k` (counting from 0) takes one square root, `n - 1 - k` divisions, and updates
  `n - j` entries of each column `j = k + 1, …, n - 1` (`updEntries`). These are exactly
  the entries changed by `stepAt` (`card_stepAt_updates`). With trailing size
  `m = n - 1 - k`, there are `m (m + 1) / 2` of them (`two_mul_updEntries`).
* Each update costs a multiplication and a subtraction, so the updates cost
  `∑_{m=1}^{n-1} m (m + 1) = (n³ - n) / 3` flops (`sum_update_flops`).
* The total `cholFlops n = n³/3 + n²/2 + n/6` has leading term `n³/3`
  (`cholFlops_ratio`), about half the `2n³/3` of LU (`luFlops`, `chol_lu_ratio`).
  Updating both triangles costs `~ 2n³/3` again (`fullUpdateFlops_ratio`).
* The lower triangle has `n (n + 1) / 2` entries (`card_lower_triangle`).
* A triangular solve costs `n²` flops, so a pair costs `2n²` (`triSolveFlops`).
-/
namespace Cholesky
open Finset Filter Topology

/-- Entries updated at step `k`: `n - j` entries in each column `j = k + 1, …, n - 1`,
as in `for j in range(k + 1, n): A[j:, j] -= ...`. -/
def updEntries (n k : ℕ) : ℕ := ∑ j ∈ Ico (k + 1) n, (n - j)

lemma two_mul_sum_range_succ (M : ℕ) : 2 * ∑ t ∈ range M, (t + 1) = M * (M + 1) := by
  induction M with
  | zero => simp
  | succ M ih => rw [Finset.sum_range_succ, mul_add, ih]; ring

lemma two_mul_updEntries (n k : ℕ) : 2 * updEntries n k = (n - 1 - k) * (n - k) := by
  unfold updEntries
  rw [Finset.sum_Ico_eq_sum_range]
  have : ∀ t ∈ range (n - (k + 1)),
      n - (k + 1 + t) = (fun t => t + 1) (n - (k + 1) - 1 - t) := by
    intro t ht; simp only [Finset.mem_range] at ht; simp only; omega
  rw [Finset.sum_congr rfl this, Finset.sum_range_reflect (fun t => t + 1) (n - (k + 1)),
    two_mul_sum_range_succ]
  by_cases hk : k < n
  · congr 1 <;> omega
  · rw [show n - (k + 1) = 0 by omega, show n - 1 - k = 0 by omega]; simp

/-- The entries changed by step `k` of `stepAt` are those with `k < j ≤ i`. Their number
is `updEntries n k`. -/
theorem card_stepAt_updates (n k : ℕ) :
    (univ.filter fun p : Fin n × Fin n => k < (p.2 : ℕ) ∧ p.2 ≤ p.1).card = updEntries n k := by
  rw [Finset.card_filter, Fintype.sum_prod_type_right]
  have h : ∀ j : Fin n, (∑ i : Fin n, if k < (j : ℕ) ∧ j ≤ i then 1 else 0) =
      if k < (j : ℕ) then n - j else 0 := by
    intro j
    by_cases hkj : k < (j : ℕ)
    · simp only [hkj, true_and, if_true]
      rw [← Finset.card_filter]
      have : (univ.filter fun i : Fin n => j ≤ i) = Finset.Ici j := by ext; simp
      rw [this, Fin.card_Ici]
    · simp [hkj]
  simp only [h]
  rw [Fin.sum_univ_eq_sum_range (fun j => if k < j then n - j else 0), ← Finset.sum_filter]
  unfold updEntries
  congr 1
  ext j; simp [Finset.mem_Ico]; omega

/-- The update flops: `∑_{m=1}^{n-1} m (m + 1) = (n³ - n) / 3`. -/
theorem sum_update_flops (n : ℕ) :
    3 * (∑ m ∈ range n, m * (m + 1)) = n ^ 3 - n := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    have h : n ≤ n ^ 3 := by
      rcases Nat.eq_zero_or_pos n with h | h
      · simp [h]
      · calc n = n ^ 1 := (pow_one n).symm
          _ ≤ n ^ 3 := Nat.pow_le_pow_right h (by norm_num)
    zify [h, show n + 1 ≤ (n + 1) ^ 3 by nlinarith]
    ring

lemma sum_update_flops_real (n : ℕ) :
    (∑ m ∈ range n, ((m : ℝ) * (m + 1))) = ((n : ℝ) ^ 3 - n) / 3 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; push_cast; ring

/-- Flops of step `k`: one square root, `n - 1 - k` divisions, and two flops per
updated entry. -/
def stepFlops (n k : ℕ) : ℕ := 1 + (n - 1 - k) + 2 * updEntries n k

/-- Total flops of the in-place algorithm. -/
def cholFlops (n : ℕ) : ℕ := ∑ k ∈ range n, stepFlops n k

lemma cholFlops_eq_sum (n : ℕ) :
    cholFlops n = ∑ m ∈ range n, (1 + m + m * (m + 1)) := by
  unfold cholFlops stepFlops
  rw [← Finset.sum_range_reflect]
  apply Finset.sum_congr rfl
  intro m hm
  rw [Finset.mem_range] at hm
  rw [two_mul_updEntries]
  have h1 : n - 1 - (n - 1 - m) = m := by omega
  have h2 : n - (n - 1 - m) = m + 1 := by omega
  rw [h1, h2]

/-- `cholFlops n = n³/3 + n²/2 + n/6`. -/
theorem cholFlops_real (n : ℕ) :
    (cholFlops n : ℝ) = (n : ℝ) ^ 3 / 3 + (n : ℝ) ^ 2 / 2 + n / 6 := by
  rw [cholFlops_eq_sum]
  push_cast
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; push_cast; ring

/-- The updates cost `(n³ - n) / 3`; the divisions `n (n - 1) / 2`; the square roots `n`. -/
theorem cholFlops_split (n : ℕ) :
    (cholFlops n : ℝ) = ((n : ℝ) ^ 3 - n) / 3 + (n : ℝ) * (n - 1) / 2 + n := by
  rw [cholFlops_real]; ring

lemma tendsto_inv_nat : Tendsto (fun n : ℕ => (1 : ℝ) / n) atTop (𝓝 0) :=
  tendsto_one_div_atTop_nhds_zero_nat

/-- The leading cost is `n³/3` flops. -/
theorem cholFlops_ratio :
    Tendsto (fun n : ℕ => (cholFlops n : ℝ) / ((n : ℝ) ^ 3 / 3)) atTop (𝓝 1) := by
  have h : Tendsto (fun n : ℕ => 1 + 3 / 2 * ((1 : ℝ) / n) + 1 / 2 * ((1 : ℝ) / n) ^ 2)
      atTop (𝓝 (1 + 3 / 2 * 0 + 1 / 2 * 0 ^ 2)) :=
    ((tendsto_const_nhds.add (tendsto_inv_nat.const_mul _)).add
      ((tendsto_inv_nat.pow 2).const_mul _))
  simp only [mul_zero, add_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true,
    zero_pow] at h
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with n hn
  have hn' : (n : ℝ) ≠ 0 := by positivity
  rw [cholFlops_real]
  field_simp
  ring

/-- Flops of unpivoted LU: at step `k`, `n - 1 - k` divisions and `(n - 1 - k)²` updates
of two flops each. -/
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

/-- Cholesky needs about half the arithmetic of LU. -/
theorem chol_lu_ratio :
    Tendsto (fun n : ℕ => (cholFlops n : ℝ) / (luFlops n : ℝ)) atTop (𝓝 (1 / 2)) := by
  have h : Tendsto (fun n : ℕ =>
      (1 / 3 + 1 / 2 * ((1 : ℝ) / n) + 1 / 6 * ((1 : ℝ) / n) ^ 2) /
        (2 / 3 - 1 / 2 * ((1 : ℝ) / n) - 1 / 6 * ((1 : ℝ) / n) ^ 2))
      atTop (𝓝 ((1 / 3 + 1 / 2 * 0 + 1 / 6 * 0 ^ 2) / (2 / 3 - 1 / 2 * 0 - 1 / 6 * 0 ^ 2))) := by
    apply Tendsto.div
    · exact (tendsto_const_nhds.add (tendsto_inv_nat.const_mul _)).add
        ((tendsto_inv_nat.pow 2).const_mul _)
    · exact (tendsto_const_nhds.sub (tendsto_inv_nat.const_mul _)).sub
        ((tendsto_inv_nat.pow 2).const_mul _)
    · norm_num
  norm_num at h
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop 2] with n hn
  have hn' : (n : ℝ) ≠ 0 := by positivity
  have hn1 : (2 : ℝ) ≤ n := by exact_mod_cast hn
  have hd : 4 * (n : ℝ) ^ 2 - 3 * n - 1 ≠ 0 := by nlinarith
  rw [cholFlops_real, luFlops_real]
  field_simp

/-- Updating the full trailing block (both triangles) costs `2 (n - 1 - k)²` flops at
step `k`, and `~ 2n³/3` in total. -/
def fullUpdateFlops (n : ℕ) : ℕ := ∑ k ∈ range n, 2 * (n - 1 - k) ^ 2

theorem fullUpdateFlops_real (n : ℕ) :
    (fullUpdateFlops n : ℝ) = 2 * (n : ℝ) ^ 3 / 3 - (n : ℝ) ^ 2 + n / 3 := by
  unfold fullUpdateFlops
  rw [← Finset.sum_range_reflect]
  have : ∀ m ∈ range n, 2 * (n - 1 - (n - 1 - m)) ^ 2 = 2 * m ^ 2 := by
    intro m hm; rw [Finset.mem_range] at hm
    rw [show n - 1 - (n - 1 - m) = m by omega]
  rw [Finset.sum_congr rfl this]
  clear this
  push_cast
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; push_cast; ring

theorem fullUpdateFlops_ratio :
    Tendsto (fun n : ℕ => (fullUpdateFlops n : ℝ) / ((n : ℝ) ^ 3 / 3)) atTop (𝓝 2) := by
  have h : Tendsto (fun n : ℕ => 2 - 3 * ((1 : ℝ) / n) + 1 * ((1 : ℝ) / n) ^ 2)
      atTop (𝓝 (2 - 3 * 0 + 1 * 0 ^ 2)) :=
    ((tendsto_const_nhds.sub (tendsto_inv_nat.const_mul _)).add
      ((tendsto_inv_nat.pow 2).const_mul _))
  norm_num at h
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with n hn
  have hn' : (n : ℝ) ≠ 0 := by positivity
  rw [fullUpdateFlops_real]
  field_simp

/-- The lower triangle of an `n × n` matrix has `n (n + 1) / 2` entries. -/
theorem card_lower_triangle (n : ℕ) :
    2 * (univ.filter fun p : Fin n × Fin n => p.2 ≤ p.1).card = n * (n + 1) := by
  rw [Finset.card_filter, Fintype.sum_prod_type]
  have h : ∀ i : Fin n, (∑ j : Fin n, if j ≤ i then 1 else 0) = (i : ℕ) + 1 := by
    intro i
    rw [← Finset.card_filter]
    have : (univ.filter fun j : Fin n => j ≤ i) = Finset.Iic i := by ext; simp
    rw [this, Fin.card_Iic]
  simp only [h]
  rw [Fin.sum_univ_eq_sum_range (fun i => i + 1), two_mul_sum_range_succ]

/-- Forward substitution: row `i` (counting from 0) takes `i` multiplications, `i`
subtractions, and one division. -/
def triSolveFlops (n : ℕ) : ℕ := ∑ i ∈ range n, (2 * i + 1)

theorem triSolveFlops_eq (n : ℕ) : triSolveFlops n = n ^ 2 := by
  unfold triSolveFlops
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; ring

/-- The two triangular solves cost `2n²` flops. -/
theorem two_solves_flops (n : ℕ) : 2 * triSolveFlops n = 2 * n ^ 2 := by
  rw [triSolveFlops_eq]

end Cholesky
