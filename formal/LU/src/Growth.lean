import src.Factorization
import src.Float

/-!
# Element growth

`stage o A k` is the array at the start of step `k` (counting from 0, so `stage o A 0 = A`;
the notes write `a^{(k+1)}`). The growth factor is

`ρ_n = max_k max_{i,j ≥ k} |a^{(k)}_ij| / max_{i,j} |a_ij|` (`growth`).

* `upd_norm_le`: `|a_ij - l_ik a_kj| ≤ |a_ij| + |a_kj|` when `|l_ik| ≤ 1`;
* `stage_bound`: in exact arithmetic `|a^{(k)}_ij| ≤ 2^k max |a_ij|` for `i, j ≥ k`, hence
  **`ρ_n ≤ 2^{n-1}`** (`growth_le`).

The matrix `wilk n` (1 on the diagonal, -1 below it, 0 above it in the first `n - 1`
columns, last column all ones) attains the bound. For any arithmetic that is exact on the
few values involved (exact arithmetic, `ExactOnW.exact`, and rounding to nearest,
`ExactOnW.rounded`):

* `run_wilk`: the array at step `k` is `wStage n k`; ties are resolved in favour of the first
  candidate, so no rows are swapped, every multiplier is `-1`, and the last column doubles;
* `wilk_factors`: `P = I`, `L` is unit lower triangular with `-1` below the diagonal, and `U`
  is the identity in the first `n - 1` columns with `u_{i,n} = 2^{i-1}`;
* `growth_wilk`: **`ρ_n = 2^{n-1}`**.
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {𝕜 : Type*} [NormedField 𝕜] [DecidableEq 𝕜] {n : ℕ}

/-! ## Maxima over finite index sets -/

lemma le_iSup_fin (f : Fin n → ℝ) (i : Fin n) : f i ≤ ⨆ i, f i :=
  le_ciSup (Set.finite_range f).bddAbove i

lemma iSup_fin_le [NeZero n] {f : Fin n → ℝ} {c : ℝ} (h : ∀ i, f i ≤ c) : ⨆ i, f i ≤ c :=
  haveI : Nonempty (Fin n) := ⟨0⟩
  ciSup_le h

/-! ## The growth factor -/

/-- The array at the start of step `k`. -/
noncomputable def stage (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) (k : ℕ) :
    Matrix (Fin n) (Fin n) 𝕜 := (run o k ⟨A, 1⟩).mat

/-- The largest active entry `max_{i,j ≥ k} |a^{(k)}_ij|`. -/
noncomputable def activeMax (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) (k : ℕ) : ℝ :=
  ⨆ i : Fin n, ⨆ j : Fin n, if k ≤ i.val ∧ k ≤ j.val then ‖stage o A k i j‖ else 0

/-- `max_{i,j} |a_ij|`. -/
noncomputable def maxAbs (A : Matrix (Fin n) (Fin n) 𝕜) : ℝ := ⨆ i : Fin n, ⨆ j : Fin n, ‖A i j‖

/-- The growth factor `ρ_n`. -/
noncomputable def growth (o : Ops 𝕜) (A : Matrix (Fin n) (Fin n) 𝕜) : ℝ :=
  (⨆ k : Fin n, activeMax o A k) / maxAbs A

lemma norm_le_maxAbs (A : Matrix (Fin n) (Fin n) 𝕜) (i j : Fin n) : ‖A i j‖ ≤ maxAbs A :=
  le_trans (le_iSup_fin (fun j => ‖A i j‖) j) (le_iSup_fin (fun i => ⨆ j, ‖A i j‖) i)

lemma maxAbs_nonneg (A : Matrix (Fin n) (Fin n) 𝕜) : 0 ≤ maxAbs A := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [maxAbs]
  · exact le_trans (norm_nonneg _) (norm_le_maxAbs A ⟨0, hn⟩ ⟨0, hn⟩)

/-- One update at most adds the two magnitudes when the multiplier is at most one. -/
lemma upd_norm_le {a l c : 𝕜} (hl : ‖l‖ ≤ 1) : ‖a - l * c‖ ≤ ‖a‖ + ‖c‖ :=
  calc ‖a - l * c‖ ≤ ‖a‖ + ‖l * c‖ := norm_sub_le _ _
    _ = ‖a‖ + ‖l‖ * ‖c‖ := by rw [norm_mul]
    _ ≤ ‖a‖ + 1 * ‖c‖ := by gcongr
    _ = ‖a‖ + ‖c‖ := by ring

lemma swap_ge {k p i : Fin n} (hp : k ≤ p) (hi : k ≤ i) : k ≤ Equiv.swap k p i := by
  rw [Equiv.swap_apply_def]
  split_ifs
  · exact hp
  · exact le_rfl
  · exact hi

/-- **Each step at most doubles the largest active entry.** -/
theorem stage_bound (A : Matrix (Fin n) (Fin n) 𝕜) : ∀ k : ℕ, ∀ i j : Fin n,
    k ≤ i.val → k ≤ j.val → ‖stage (Ops.exact 𝕜) A k i j‖ ≤ 2 ^ k * maxAbs A
  | 0, i, j, _, _ => by
    simp only [stage, run, runFrom, pow_zero, one_mul]
    exact norm_le_maxAbs A i j
  | k + 1, i, j, hi, hj => by
    have ih := stage_bound A k
    have hkn : k < n := by omega
    set kk : Fin n := ⟨k, hkn⟩
    have hki : kk < i := by rw [Fin.lt_iff_val_lt_val]; simp [kk]; omega
    have hkj : kk < j := by rw [Fin.lt_iff_val_lt_val]; simp [kk]; omega
    set M := stage (Ops.exact 𝕜) A k
    have hstep : stage (Ops.exact 𝕜) A (k + 1) = elim (Ops.exact 𝕜) kk
        (swapRows M kk (pivRow M kk)) := by
      simp only [stage, run_succ, stepN, dif_pos hkn]
      rfl
    have hp : kk ≤ pivRow M kk := le_firstMax _ _
    set M1 := swapRows M kk (pivRow M kk)
    have hM1 : ∀ i' j', kk ≤ i' → k ≤ j'.val → ‖M1 i' j'‖ ≤ 2 ^ k * maxAbs A := by
      intro i' j' hi' hj'
      exact ih _ _ (swap_ge hp hi') hj'
    have hpiv : ∀ i', kk ≤ i' → ‖M1 i' kk‖ ≤ ‖M1 kk kk‖ := by
      intro i' hi'
      simp only [M1, swapRows, Matrix.submatrix_apply, id, Equiv.swap_apply_left]
      exact le_at_firstMax (fun i => ‖M i kk‖) kk (swap_ge hp hi')
    have h2 : (2 : ℝ) ^ k * maxAbs A ≤ 2 ^ (k + 1) * maxAbs A := by
      gcongr
      · exact maxAbs_nonneg A
      · norm_num
      · omega
    rw [hstep]
    by_cases h0 : M1 kk kk = 0
    · rw [elim_of_pivot_eq_zero h0]
      exact le_trans (hM1 i j hki.le (by omega)) h2
    · rw [elim_apply_upd h0 hki hkj, Ops.exact_upd, Ops.exact_div]
      have hl : ‖M1 i kk / M1 kk kk‖ ≤ 1 := by
        rw [norm_div, div_le_one (norm_pos_iff.mpr h0)]; exact hpiv i hki.le
      calc ‖M1 i j - M1 i kk / M1 kk kk * M1 kk j‖ ≤ ‖M1 i j‖ + ‖M1 kk j‖ := upd_norm_le hl
        _ ≤ 2 ^ k * maxAbs A + 2 ^ k * maxAbs A :=
            add_le_add (hM1 i j hki.le (by omega)) (hM1 kk j le_rfl (by omega))
        _ = 2 ^ (k + 1) * maxAbs A := by ring

lemma activeMax_le (A : Matrix (Fin n) (Fin n) 𝕜) (k : ℕ) (hk : k < n) :
    activeMax (Ops.exact 𝕜) A k ≤ 2 ^ k * maxAbs A := by
  haveI : NeZero n := ⟨by omega⟩
  have h0 : 0 ≤ 2 ^ k * maxAbs A := mul_nonneg (by positivity) (maxAbs_nonneg A)
  apply iSup_fin_le
  intro i
  apply iSup_fin_le
  intro j
  split_ifs with h
  · exact stage_bound A k i j h.1 h.2
  · exact h0

/-- **`ρ_n ≤ 2^{n-1}`.** -/
theorem growth_le (A : Matrix (Fin n) (Fin n) 𝕜) : growth (Ops.exact 𝕜) A ≤ 2 ^ (n - 1) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [growth, maxAbs]
  haveI : NeZero n := ⟨by omega⟩
  have hm := maxAbs_nonneg A
  have hsup : (⨆ k : Fin n, activeMax (Ops.exact 𝕜) A k) ≤ 2 ^ (n - 1) * maxAbs A := by
    apply iSup_fin_le
    intro k
    calc activeMax (Ops.exact 𝕜) A k ≤ 2 ^ k.val * maxAbs A := activeMax_le A k k.2
      _ ≤ 2 ^ (n - 1) * maxAbs A := by
          gcongr
          · norm_num
          · omega
  rcases eq_or_lt_of_le hm with h | h
  · rw [growth, ← h, div_zero]; positivity
  · rw [growth, div_le_iff₀ h]; exact hsup

/-! ## A matrix that attains the bound -/

/-- `1` on the diagonal, `-1` below, `0` above in the first `n - 1` columns; last column all
ones. -/
def wilk (n : ℕ) : Matrix (Fin n) (Fin n) ℝ := Matrix.of fun i j =>
  if j.val + 1 = n then 1 else if i = j then 1 else if j < i then -1 else 0

/-- The array at the start of step `k`: the last column holds `2^{min(i, k)}`. -/
def wStage (n k : ℕ) : Matrix (Fin n) (Fin n) ℝ := Matrix.of fun i j =>
  if j.val + 1 = n then 2 ^ min i.val k else if i = j then 1 else if j < i then -1 else 0

/-- An arithmetic that is exact on the values occurring in the elimination of `wilk n`. -/
structure ExactOnW (o : Ops ℝ) : Prop where
  div : o.div (-1) 1 = -1
  upd0 : ∀ a : ℝ, (a = 0 ∨ a = 1 ∨ a = -1) → o.upd a (-1) 0 = a
  upd2 : ∀ k : ℕ, o.upd (2 ^ k) (-1) (2 ^ k) = 2 ^ (k + 1)

theorem ExactOnW.exact : ExactOnW (Ops.exact ℝ) where
  div := by norm_num
  upd0 a _ := by simp
  upd2 k := by simp [pow_succ]; ring

theorem ExactOnW.rounded {p : ℕ} (hp : 1 ≤ p) {r : ℝ → ℝ} (hr : IsRoundNearest p r) :
    ExactOnW (Ops.rounded r) where
  div := by
    show r (-1 / 1) = -1
    rw [div_one]; exact round_eq_self hr (isFloat_one hp).neg
  upd0 a ha := by
    show r (a - r (-1 * 0)) = a
    rw [mul_zero, round_zero hr, sub_zero]
    rcases ha with rfl | rfl | rfl
    · exact round_zero hr
    · exact round_eq_self hr (isFloat_one hp)
    · exact round_eq_self hr (isFloat_one hp).neg
  upd2 k := by
    show r (2 ^ k - r (-1 * 2 ^ k)) = 2 ^ (k + 1)
    rw [neg_one_mul, round_eq_self hr (isFloat_two_pow hp k).neg, sub_neg_eq_add, ← two_mul,
      ← pow_succ']
    exact round_eq_self hr (isFloat_two_pow hp (k + 1))

lemma wilk_eq_wStage (n : ℕ) : wilk n = wStage n 0 := by
  ext i j; simp [wilk, wStage]

lemma wStage_last {k : ℕ} {i j : Fin n} (hj : j.val + 1 = n) : wStage n k i j = 2 ^ min i.val k := by
  simp [wStage, hj]

lemma wStage_not_last {k : ℕ} {i j : Fin n} (hj : j.val + 1 ≠ n) :
    wStage n k i j = if i = j then 1 else if j < i then -1 else 0 := by
  simp [wStage, hj]

/-- One step of the elimination of `wilk n`. -/
lemma step_wilk (o : Ops ℝ) (ho : ExactOnW o) (k : ℕ) (hk : k + 2 ≤ n) :
    stepAt o ⟨k, by omega⟩ ⟨wStage n k, 1⟩ = ⟨wStage n (k + 1), 1⟩ := by
  set kk : Fin n := ⟨k, by omega⟩
  have hkl : kk.val + 1 ≠ n := by simp [kk]; omega
  have hkk : wStage n k kk kk = 1 := by rw [wStage_not_last hkl]; simp
  -- the pivot: all candidates have magnitude one, the first is chosen
  have hp : pivRow (wStage n k) kk = kk := by
    apply firstMax_self
    intro i hi
    simp only [wStage_not_last hkl]
    split_ifs <;> simp
  have hsw : swapRows (wStage n k) kk kk = wStage n k := by ext i j; simp [swapRows]
  simp only [stepAt, hp, hsw, Equiv.swap_self, Equiv.Perm.mul_refl]
  congr 1
  ext i j
  by_cases hi : kk < i
  · by_cases hj : j = kk
    · subst hj
      rw [elim_apply_mult (by rw [hkk]; exact one_ne_zero) hi, hkk, wStage_not_last hkl,
        wStage_not_last hkl, if_neg (ne_of_gt hi), if_pos hi]
      exact ho.div
    · by_cases hj2 : kk < j
      · rw [elim_apply_upd (by rw [hkk]; exact one_ne_zero) hi hj2, hkk]
        have hmult : o.div (wStage n k i kk) 1 = -1 := by
          rw [wStage_not_last hkl, if_neg (ne_of_gt hi), if_pos hi]; exact ho.div
        rw [hmult]
        by_cases hjl : j.val + 1 = n
        · rw [wStage_last hjl, wStage_last hjl, wStage_last hjl]
          have h1 : min i.val k = k := by
            have := Fin.lt_iff_val_lt_val.mp hi; simp [kk] at this; omega
          have h2 : min kk.val k = k := by simp [kk]
          have h3 : min i.val (k + 1) = k + 1 := by
            have := Fin.lt_iff_val_lt_val.mp hi; simp [kk] at this; omega
          rw [h1, h2, h3]
          exact ho.upd2 k
        · rw [wStage_not_last hjl, wStage_not_last hjl, wStage_not_last hjl,
            if_neg (ne_of_lt hj2), if_neg (not_lt.mpr hj2.le)]
          apply ho.upd0
          split_ifs <;> simp
      · have hjk : j < kk := lt_of_le_of_ne (not_lt.mp hj2) hj
        rw [elim_apply_of_lt_col hjk]
        have hjl : j.val + 1 ≠ n := by
          have := Fin.lt_iff_val_lt_val.mp hjk; simp [kk] at this; omega
        rw [wStage_not_last hjl, wStage_not_last hjl]
  · rw [elim_apply_of_not_lt hi]
    by_cases hjl : j.val + 1 = n
    · rw [wStage_last hjl, wStage_last hjl]
      have : i.val ≤ k := by
        have := not_lt.mp hi; rw [Fin.le_iff_val_le_val] at this; simpa [kk] using this
      rw [min_eq_left this, min_eq_left (by omega)]
    · rw [wStage_not_last hjl, wStage_not_last hjl]

/-- **The elimination of `wilk n`.** -/
theorem run_wilk (o : Ops ℝ) (ho : ExactOnW o) :
    ∀ k : ℕ, k + 1 ≤ n → run o k ⟨wilk n, 1⟩ = ⟨wStage n k, 1⟩
  | 0, _ => by rw [wilk_eq_wStage]; rfl
  | k + 1, hk => by
    rw [run_succ, run_wilk o ho k (by omega), stepN, dif_pos (by omega)]
    exact step_wilk o ho k (by omega)

theorem final_wilk (o : Ops ℝ) (ho : ExactOnW o) (hn : 1 ≤ n) :
    final o (wilk n) = ⟨wStage n (n - 1), 1⟩ := by
  rw [final, run_pred]
  exact run_wilk o ho (n - 1) (by omega)

/-- **The factors of `wilk n`**: no interchanges, `l_ij = -1` below the diagonal,
`u_{i,n} = 2^{i-1}` (here `2^i` for 0-based `i`), and `U` is the identity elsewhere. -/
theorem wilk_factors (o : Ops ℝ) (ho : ExactOnW o) (hn : 1 ≤ n) :
    factorPerm o (wilk n) = 1 ∧
      (∀ i j : Fin n, factorL o (wilk n) i j = if i = j then 1 else if j < i then -1 else 0) ∧
      (∀ i j : Fin n, factorU o (wilk n) i j =
        if j.val + 1 = n then 2 ^ i.val else if i = j then 1 else 0) := by
  have hf := final_wilk o ho hn
  refine ⟨by rw [factorPerm, hf], fun i j => ?_, fun i j => ?_⟩
  · rw [factorL, hf]
    simp only [lowerOf, Matrix.of_apply]
    split_ifs with h1 h2 <;> try rfl
    · have hjl : j.val + 1 ≠ n := by
        have := Fin.lt_iff_val_lt_val.mp h2; omega
      rw [wStage_not_last hjl, if_neg h1, if_pos h2]
  · rw [factorU, hf]
    simp only [upperOf, Matrix.of_apply]
    by_cases hij : i ≤ j
    · rw [if_pos hij]
      by_cases hjl : j.val + 1 = n
      · rw [wStage_last hjl, if_pos hjl, min_eq_left (by omega)]
      · rw [wStage_not_last hjl, if_neg hjl]
        split_ifs with h1 h2 <;> try rfl
        exact absurd h2 (not_lt.mpr hij)
    · rw [if_neg hij]
      have hji : j < i := lt_of_not_ge hij
      by_cases hjl : j.val + 1 = n
      · exfalso; have := Fin.lt_iff_val_lt_val.mp hji; omega
      · rw [if_neg hjl, if_neg (ne_of_gt hji)]

lemma maxAbs_wilk (hn : 1 ≤ n) : maxAbs (wilk n) = 1 := by
  haveI : NeZero n := ⟨by omega⟩
  apply le_antisymm
  · apply iSup_fin_le; intro i; apply iSup_fin_le; intro j
    simp only [wilk, Matrix.of_apply]
    split_ifs <;> simp
  · have := norm_le_maxAbs (wilk n) ⟨0, by omega⟩ ⟨n - 1, by omega⟩
    simpa [wilk, show n - 1 + 1 = n by omega] using this

lemma activeMax_wilk (k : ℕ) (hk : k + 1 ≤ n) : activeMax (Ops.exact ℝ) (wilk n) k = 2 ^ k := by
  haveI : NeZero n := ⟨by omega⟩
  have hst : stage (Ops.exact ℝ) (wilk n) k = wStage n k := by
    rw [stage, run_wilk _ ExactOnW.exact k hk]
  apply le_antisymm
  · apply iSup_fin_le; intro i; apply iSup_fin_le; intro j
    rw [hst]
    split_ifs with h
    · by_cases hjl : j.val + 1 = n
      · rw [wStage_last hjl, min_eq_right h.1]; simp
      · rw [wStage_not_last hjl]
        have : (1 : ℝ) ≤ 2 ^ k := one_le_pow₀ (by norm_num)
        split_ifs <;> simp [this]
    · positivity
  · have h1 := le_iSup_fin (fun j : Fin n => if k ≤ (⟨k, by omega⟩ : Fin n).val ∧ k ≤ j.val then
      ‖stage (Ops.exact ℝ) (wilk n) k ⟨k, by omega⟩ j‖ else 0) ⟨n - 1, by omega⟩
    have h2 := le_iSup_fin (fun i : Fin n => ⨆ j : Fin n, if k ≤ i.val ∧ k ≤ j.val then
      ‖stage (Ops.exact ℝ) (wilk n) k i j‖ else 0) ⟨k, by omega⟩
    refine le_trans ?_ (le_trans h1 h2)
    rw [if_pos ⟨le_rfl, by simp; omega⟩, hst, wStage_last (by simp; omega)]
    simp

/-- **The growth factor of `wilk n` is `2^{n-1}`.** -/
theorem growth_wilk (hn : 1 ≤ n) : growth (Ops.exact ℝ) (wilk n) = 2 ^ (n - 1) := by
  haveI : NeZero n := ⟨by omega⟩
  rw [growth, maxAbs_wilk hn, div_one]
  apply le_antisymm
  · apply iSup_fin_le
    intro k
    rw [activeMax_wilk k.val (by omega)]
    exact pow_le_pow_right₀ (by norm_num) (by omega)
  · have := le_iSup_fin (fun k : Fin n => activeMax (Ops.exact ℝ) (wilk n) k) ⟨n - 1, by omega⟩
    simp only [] at this
    rwa [activeMax_wilk (n - 1) (by omega)] at this

end LU
