import src.Counterexample.Rounding

/-!
# Following the rounded leading elimination

For rounding to nearest at precision `p` (either variant of the update), we compute the
trailing block after the first `7k` steps of the in-place algorithm on `A_{9k}`.

* The first `6k` rows of `T` subtract `t² = 9h/8` from every trailing entry. Diagonal
  entries do not change, within-group off-diagonal entries decrease exactly to `-η/8`,
  and cross-block entries decrease by `2h` per step, to `±9/(4√k) - η` (`phase1`).
* The next `k` rows cancel `G` exactly (`phase2`).
* The result is the matrix `Ŝ` of the notes (`trailing_eq_Shat`).
-/
namespace Cholesky
namespace CE
namespace Par
open Matrix

variable (P : Par)

/-! ## Entries of `T` -/

lemma T_ll (r : Fin (6 * P.k)) (c : Fin P.k) : P.T (Fin.castAdd P.k r) (Fin.castAdd P.k c) = P.t := by
  simp only [T, blk_ll, Matrix.smul_apply, ones, Matrix.of_apply, smul_eq_mul, mul_one]
lemma T_lr (r : Fin (6 * P.k)) (c : Fin P.k) : P.T (Fin.castAdd P.k r) (Fin.natAdd P.k c) = P.t := by
  simp only [T, blk_lr, Matrix.smul_apply, ones, Matrix.of_apply, smul_eq_mul, mul_one]
lemma T_rl (r c : Fin P.k) :
    P.T (Fin.natAdd (6 * P.k) r) (Fin.castAdd P.k c) = 3 / 2 * (if r = c then 1 else 0) := by
  simp only [T, blk_rl, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
lemma T_rr (r c : Fin P.k) : P.T (Fin.natAdd (6 * P.k) r) (Fin.natAdd P.k c) = 3 / 2 * P.Q r c := by
  simp only [T, blk_rr, Matrix.smul_apply, smul_eq_mul]

/-- The sign `±1` of the Hadamard entry. -/
noncomputable def hInt (i j : Fin P.k) : ℤ := if hadamard (2 * P.q) i j = 1 then 1 else -1

lemma hInt_cases (i j : Fin P.k) : P.hInt i j = 1 ∨ P.hInt i j = -1 := by
  unfold hInt; split_ifs <;> simp

lemma hadamard_eq_hInt (i j : Fin P.k) : hadamard (2 * P.q) i j = P.hInt i j := by
  unfold hInt
  rcases hadamard_entry (2 * P.q) i j with h | h <;> simp [h]

lemma Q_eq (i j : Fin P.k) : P.Q i j = (P.hInt i j : ℝ) / P.K := by
  rw [Q, hadQ_apply, P.hadamard_eq_hInt, P.zpow_neg_q]; ring

lemma hInt_sq (i j : Fin P.k) : (P.hInt i j : ℝ) ^ 2 = 1 := by
  rcases P.hInt_cases i j with h | h <;> rw [h] <;> norm_num

/-! ## Floats -/

lemma isFloat_three_mul (e : ℤ) : IsFloat P.p (3 * (2 : ℝ) ^ e) := by
  refine ⟨3, e, ?_, by norm_num⟩
  have := P.p_ge
  calc |(3 : ℤ)| = 3 := by norm_num
    _ < 2 ^ 2 := by norm_num
    _ ≤ 2 ^ P.p := pow_le_pow_right₀ (by norm_num) (by omega)

lemma isFloat_nine_mul (e : ℤ) : IsFloat P.p (9 * (2 : ℝ) ^ e) := by
  refine ⟨9, e, ?_, by norm_num⟩
  have := P.p_ge
  calc |(9 : ℤ)| = 9 := by norm_num
    _ < 2 ^ 4 := by norm_num
    _ ≤ 2 ^ P.p := pow_le_pow_right₀ (by norm_num) (by omega)

lemma isFloat_small_int (z : ℤ) (hz : |z| ≤ 4) : IsFloat P.p (z : ℝ) := by
  have := isFloat_of_abs_le P.p_pos z 0 (by
    have := P.p_ge
    calc |z| ≤ 4 := hz
      _ = 2 ^ 2 := by norm_num
      _ ≤ 2 ^ P.p := pow_le_pow_right₀ (by norm_num) (by omega))
  simpa using this

lemma t_float : IsFloat P.p P.t := P.isFloat_three_mul _

lemma two_t_float : IsFloat P.p (2 * P.t) := by
  rw [t, ← mul_assoc, mul_comm 2 3, mul_assoc,
    show (2 : ℝ) * 2 ^ (-(((P.p + P.q + 2) / 2 : ℕ) : ℤ)) =
      2 ^ (1 + -(((P.p + P.q + 2) / 2 : ℕ) : ℤ)) by rw [zpow_add₀ (by norm_num)]; norm_num]
  exact P.isFloat_three_mul _

lemma inv_K_zpow : 1 / P.K = (2 : ℝ) ^ (-(P.q : ℤ)) := by rw [P.zpow_neg_q]

/-- Every entry of `T` and of `2T` is a float. -/
theorem T_float (i : Fin (6 * P.k + P.k)) (j : Fin (P.k + P.k)) :
    IsFloat P.p (P.T i j) ∧ IsFloat P.p (2 * P.T i j) := by
  refine Fin.addCases (fun r => ?_) (fun r => ?_) i <;>
    refine Fin.addCases (fun c => ?_) (fun c => ?_) j
  · rw [T_ll]; exact ⟨P.t_float, P.two_t_float⟩
  · rw [T_lr]; exact ⟨P.t_float, P.two_t_float⟩
  · rw [T_rl]
    split_ifs
    · refine ⟨?_, ?_⟩
      · have := P.isFloat_three_mul (-1); norm_num at this ⊢; convert this using 1
      · have := P.isFloat_small_int 3 (by norm_num); norm_num at this ⊢; exact this
    · simp only [mul_zero]; exact ⟨isFloat_zero _, isFloat_zero _⟩
  · rw [T_rr, P.Q_eq]
    rcases P.hInt_cases r c with h | h <;> rw [h] <;> refine ⟨?_, ?_⟩
    · have := P.isFloat_three_mul (-1 + -(P.q : ℤ))
      rw [zpow_add₀ (by norm_num), ← P.inv_K_zpow] at this; convert this using 1; norm_num; ring
    · have := P.isFloat_three_mul (-(P.q : ℤ))
      rw [← P.inv_K_zpow] at this; convert this using 1; push_cast; ring
    · have := (P.isFloat_three_mul (-1 + -(P.q : ℤ))).neg
      rw [zpow_add₀ (by norm_num), ← P.inv_K_zpow] at this; convert this using 1; norm_num; ring
    · have := (P.isFloat_three_mul (-(P.q : ℤ))).neg
      rw [← P.inv_K_zpow] at this; convert this using 1; push_cast; ring

lemma sqrt_four : Real.sqrt 4 = 2 := by
  rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]

/-- The rounded arithmetic satisfies the hypotheses of the leading-block lemma. -/
theorem leadOps {o : Ops} {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ho : IsRN o r P.p) :
    LeadOps o P.T where
  sqrt4 := by rw [ho.sqrt, sqrt_four]; exact round_eq_self hr (by simpa using P.isFloat_small_int 2 (by norm_num))
  div0 := by rw [ho.div, zero_div, round_zero hr]
  divT i j := by
    rw [ho.div, mul_div_cancel_left₀ _ (by norm_num : (2 : ℝ) ≠ 0)]
    exact round_eq_self hr (P.T_float i j).1
  upd0 x b c hx hbc := by
    have h0 : b * c = 0 := by rcases hbc with rfl | rfl <;> simp
    rw [ho.upd x b c (by rw [h0]; exact isFloat_zero _), h0, sub_zero]
    apply round_eq_self hr
    rcases hx with rfl | rfl | ⟨i, j, rfl⟩
    · exact isFloat_zero _
    · simpa using P.isFloat_small_int 4 (by norm_num)
    · exact (P.T_float i j).2

/-! ## The identity matrix in blocks -/

lemma one_ll (i j : Fin P.k) : (1 : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ)
    (Fin.castAdd P.k i) (Fin.castAdd P.k j) = (1 : Matrix (Fin P.k) (Fin P.k) ℝ) i j := by
  rw [one_eq_blk, blk_ll]
lemma one_rl (i j : Fin P.k) : (1 : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ)
    (Fin.natAdd P.k i) (Fin.castAdd P.k j) = 0 := by
  rw [one_eq_blk, blk_rl]; rfl
lemma one_rr (i j : Fin P.k) : (1 : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ)
    (Fin.natAdd P.k i) (Fin.natAdd P.k j) = (1 : Matrix (Fin P.k) (Fin P.k) ℝ) i j := by
  rw [one_eq_blk, blk_rr]

lemma C_ll (i j : Fin P.k) : P.C (Fin.castAdd P.k i) (Fin.castAdd P.k j) =
    if i = j then P.dg else P.η := by
  simp only [C, G, J2, ones, Matrix.add_apply, Matrix.smul_apply, blk_ll, Matrix.of_apply,
    P.one_ll, Matrix.one_apply, smul_eq_mul, dg]
  split_ifs <;> ring
lemma C_rl (i j : Fin P.k) : P.C (Fin.natAdd P.k i) (Fin.castAdd P.k j) =
    9 / 4 * P.Q j i + P.η := by
  simp only [C, G, J2, ones, Matrix.add_apply, Matrix.smul_apply, blk_rl, Matrix.of_apply,
    P.one_rl, Matrix.transpose_apply, smul_eq_mul]
  ring
lemma C_rr (i j : Fin P.k) : P.C (Fin.natAdd P.k i) (Fin.natAdd P.k j) =
    if i = j then P.dg else P.η := by
  simp only [C, G, J2, ones, Matrix.add_apply, Matrix.smul_apply, blk_rr, Matrix.of_apply,
    P.one_rr, Matrix.one_apply, smul_eq_mul, dg]
  split_ifs <;> ring

/-! ## The first `6k` steps -/

/-- Within-group off-diagonal value after `s` small updates: `η - s t² = (48K² - 9s) E₀`. -/
noncomputable def offv (s : ℕ) : ℝ := ((48 * (P.Kn : ℤ) ^ 2 - 9 * s : ℤ) : ℝ) * P.E0

/-- The trailing block after `s ≤ 6k` steps (lower triangle). -/
noncomputable def P1 (s : ℕ) : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  blk (Matrix.of fun i j => if i = j then P.dg else P.offv s) 0
    (Matrix.of fun i j => (P.Ncross (P.hInt j i) s : ℝ) * P.E0)
    (Matrix.of fun i j => if i = j then P.dg else P.offv s)

lemma offv_zero : P.offv 0 = P.η := by
  rw [offv, P.η_eq, P.K_cast]; push_cast; ring

lemma cross_zero (i j : Fin P.k) : (P.Ncross (P.hInt j i) 0 : ℝ) * P.E0 = 9 / 4 * P.Q j i + P.η := by
  rw [← P.cross_eq, P.Q_eq]; push_cast; ring

lemma lowerEq_P1_zero {o : Ops} : LowerEq (Citer o P.T P.C 0) (P.P1 0) := by
  intro i j hij
  simp only [Citer]
  revert hij
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> intro hij
  · rw [C_ll]; simp only [P1, blk_ll, Matrix.of_apply, P.offv_zero]
  · exact absurd hij not_nA_le_cA
  · rw [C_rl]; simp only [P1, blk_rl, Matrix.of_apply, P.cross_zero]
  · rw [C_rr]; simp only [P1, blk_rr, Matrix.of_apply, P.offv_zero]

lemma citer_succ_apply {o : Ops} (s : ℕ) (h : s < 6 * P.k + P.k) (i j : Fin (P.k + P.k))
    (hij : j ≤ i) : Citer o P.T P.C (s + 1) i j =
      o.upd (Citer o P.T P.C s i j) (P.T ⟨s, h⟩ i) (P.T ⟨s, h⟩ j) := by
  simp only [Citer, Matrix.of_apply, dif_pos h, if_pos hij]

lemma tt_eq : P.t * P.t = 9 * P.E0 := by rw [← sq, P.t_sq]

lemma tt_float : IsFloat P.p (P.t * P.t) := by
  rw [tt_eq]
  have := P.isFloat_E0 9 (by
    rw [P.two_pow_p_nat]
    have hK : (8 : ℤ) ≤ P.Kn := by exact_mod_cast P.Kn_ge
    have hD : (1 : ℤ) ≤ P.Dn := by exact_mod_cast P.Dn_ge
    rw [abs_of_nonneg (by norm_num)]
    nlinarith [pow_pos (by linarith : (0 : ℤ) < P.Kn) 3,
      le_mul_of_one_le_right (by positivity : (0 : ℤ) ≤ (P.Kn : ℤ) ^ 3) hD])
  simpa using this

lemma round_offv {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (s : ℕ) (hs : s + 1 ≤ 6 * P.k) :
    r (P.offv s - P.t * P.t) = P.offv (s + 1) := by
  have hK : (8 : ℤ) ≤ P.Kn := by exact_mod_cast P.Kn_ge
  have hD : (1 : ℤ) ≤ P.Dn := by exact_mod_cast P.Dn_ge
  have hs' : (s : ℤ) + 1 ≤ 6 * (P.Kn : ℤ) ^ 2 := by
    have := P.k_nat; rw [this] at hs; exact_mod_cast hs
  have h : P.offv s - P.t * P.t = P.offv (s + 1) := by
    rw [offv, offv, P.tt_eq]; push_cast; ring
  rw [h]
  apply round_eq_self hr
  apply P.isFloat_E0
  rw [P.two_pow_p_nat, abs_le]
  have h3 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn :=
    le_mul_of_one_le_right (by positivity) hD
  have h32 : 8 * (P.Kn : ℤ) ^ 2 ≤ (P.Kn : ℤ) ^ 3 := by nlinarith [sq_nonneg (P.Kn : ℤ)]
  push_cast
  constructor <;> nlinarith

/-- **The first `6k` steps.** -/
theorem phase1 {o : Ops} {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ho : IsRN o r P.p) :
    ∀ s ≤ 6 * P.k, LowerEq (Citer o P.T P.C s) (P.P1 s)
  | 0, _ => P.lowerEq_P1_zero
  | s + 1, hs => by
    intro i j hij
    have ih := phase1 hr ho s (by omega)
    have hs' : s < 6 * P.k + P.k := by omega
    rw [P.citer_succ_apply s hs' i j hij, ih i j hij]
    have hrow : (⟨s, hs'⟩ : Fin (6 * P.k + P.k)) = Fin.castAdd P.k ⟨s, by omega⟩ := rfl
    rw [hrow]
    revert hij
    refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
      refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> intro hij
    · rw [T_ll, T_ll, ho.upd _ _ _ P.tt_float]
      simp only [P1, blk_ll, Matrix.of_apply]
      split_ifs
      · rw [← sq]; exact P.round_diag hr
      · exact P.round_offv hr s hs
    · exact absurd hij not_nA_le_cA
    · rw [T_lr, T_ll, ho.upd _ _ _ P.tt_float]
      simp only [P1, blk_rl, Matrix.of_apply]
      rw [← sq]
      exact P.round_cross hr _ (P.hInt_cases j i) s hs
    · rw [T_lr, T_lr, ho.upd _ _ _ P.tt_float]
      simp only [P1, blk_rr, Matrix.of_apply]
      split_ifs
      · rw [← sq]; exact P.round_diag hr
      · exact P.round_offv hr s hs


/-! ## The next `k` steps -/

/-- `∑_{l < r} Hₗᵢ Hₗⱼ`. -/
noncomputable def cSum (r : ℕ) (i j : Fin P.k) : ℤ :=
  ∑ l : Fin P.k, if (l : ℕ) < r then P.hInt l i * P.hInt l j else 0

lemma cSum_zero (i j : Fin P.k) : P.cSum 0 i j = 0 := by simp [cSum]

lemma cSum_succ (r : ℕ) (hr : r < P.k) (i j : Fin P.k) :
    P.cSum (r + 1) i j = P.cSum r i j + P.hInt ⟨r, hr⟩ i * P.hInt ⟨r, hr⟩ j := by
  unfold cSum
  rw [← Finset.sum_erase_add _ _ (Finset.mem_univ (⟨r, hr⟩ : Fin P.k)),
    ← Finset.sum_erase_add (s := Finset.univ) (a := (⟨r, hr⟩ : Fin P.k)) _ (Finset.mem_univ _)]
  simp only [lt_add_one, if_true, lt_irrefl, if_false, add_zero]
  congr 1
  apply Finset.sum_congr rfl
  intro l hl
  have hl' : (l : ℕ) ≠ r := fun h => (Finset.mem_erase.mp hl).1 (Fin.ext h)
  have : (l : ℕ) < r + 1 ↔ (l : ℕ) < r := by omega
  simp only [this]

lemma cSum_bound : ∀ (r : ℕ), r ≤ P.k → ∀ i j, |P.cSum r i j| ≤ r
  | 0, _, i, j => by simp [P.cSum_zero]
  | r + 1, hr, i, j => by
    rw [P.cSum_succ r (by omega)]
    have ih := cSum_bound r (by omega) i j
    have hp : |P.hInt ⟨r, by omega⟩ i * P.hInt ⟨r, by omega⟩ j| = 1 := by
      rcases P.hInt_cases ⟨r, by omega⟩ i with h | h <;>
        rcases P.hInt_cases ⟨r, by omega⟩ j with h' | h' <;> simp [h, h']
    calc |P.cSum r i j + P.hInt ⟨r, by omega⟩ i * P.hInt ⟨r, by omega⟩ j|
        ≤ |P.cSum r i j| + |P.hInt ⟨r, by omega⟩ i * P.hInt ⟨r, by omega⟩ j| := abs_add_le _ _
      _ ≤ r + 1 := by rw [hp]; linarith
      _ = ((r + 1 : ℕ) : ℤ) := by push_cast; ring

lemma cSum_full (i j : Fin P.k) (hij : i ≠ j) : P.cSum P.k i j = 0 := by
  have h := congrFun₂ (hadamard_mul_transpose (2 * P.q)) i j
  rw [Matrix.mul_apply, Matrix.smul_apply, Matrix.one_apply_ne hij, smul_zero] at h
  have hsymm := hadamard_transpose (2 * P.q)
  unfold cSum
  have : ∀ l : Fin P.k, (if (l : ℕ) < P.k then P.hInt l i * P.hInt l j else 0) =
      P.hInt l i * P.hInt l j := fun l => by rw [if_pos l.2]
  simp only [this]
  have hR : ((∑ l : Fin P.k, P.hInt l i * P.hInt l j : ℤ) : ℝ) = 0 := by
    push_cast
    rw [← h]
    apply Finset.sum_congr rfl
    intro l _
    rw [← P.hadamard_eq_hInt, ← P.hadamard_eq_hInt, Matrix.transpose_apply]
    rw [show hadamard (2 * P.q) l i = hadamard (2 * P.q) i l by
      rw [← Matrix.transpose_apply (hadamard (2 * P.q)) i l, hsymm]]
    congr 1
    rw [← Matrix.transpose_apply (hadamard (2 * P.q)) l j, hsymm]
  exact_mod_cast hR

/-- `η + δ` in units of `E₀`. -/
def Nηδ : ℤ := 48 * (P.Kn : ℤ) ^ 2 + 48 * (P.Kn : ℤ) ^ 4

/-- The trailing block after `6k + r` steps, `r ≤ k`, in units of `E₀`. -/
noncomputable def P2 (r : ℕ) : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ :=
  blk (Matrix.of fun i j => ((if i = j then (if (i : ℕ) < r then P.Nηδ else P.Nd)
        else -6 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) * P.E0) 0
    (Matrix.of fun i j => ((if (j : ℕ) < r then -48 * (P.Kn : ℤ) ^ 2
        else P.hInt j i * (2304 * (P.Kn : ℤ) ^ 3 * P.Dn) - 48 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) * P.E0)
    (Matrix.of fun i j => ((if i = j then P.Nd - 2304 * r * (P.Kn : ℤ) ^ 2 * P.Dn
        else -6 * (P.Kn : ℤ) ^ 2 - 2304 * P.cSum r i j * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) * P.E0)

lemma P1_eq_P2 : LowerEq (P.P1 (6 * P.k)) (P.P2 0) := by
  intro i j hij
  have hk : (6 * P.k : ℕ) = 6 * (P.Kn : ℤ) ^ 2 := by rw [P.k_nat]; push_cast; ring
  revert hij
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> intro hij
  · simp only [P1, P2, blk_ll, Matrix.of_apply, offv, P.dg_eq, Nat.not_lt_zero, if_false]
    split_ifs
    · rfl
    · push_cast [hk]; ring
  · exact absurd hij not_nA_le_cA
  · simp only [P1, P2, blk_rl, Matrix.of_apply, Ncross, Nat.not_lt_zero, if_false]
    push_cast [hk]; ring
  · simp only [P1, P2, blk_rr, Matrix.of_apply, offv, P.dg_eq, P.cSum_zero]
    split_ifs
    · push_cast; ring
    · push_cast [hk]; ring


/-! ### Floats on the grid `4u` -/

lemma isFloat_grid16 (N M : ℤ) (hN : N = 16 * P.Kn * M) (hM : |M| ≤ 2 ^ P.p) :
    IsFloat P.p ((N : ℝ) * P.E0) := by
  have := P.isFloat_grid4 M hM
  rw [hN, P.K_cast] at *; push_cast; convert this using 1; ring

section Bounds
variable {P}

lemma kbounds : (8 : ℤ) ≤ P.Kn ∧ (1 : ℤ) ≤ P.Dn ∧ (P.Kn : ℤ) = 8 * P.Ln ∧ (1 : ℤ) ≤ P.Ln := by
  refine ⟨by exact_mod_cast P.Kn_ge, by exact_mod_cast P.Dn_ge, by exact_mod_cast P.Kn_eq,
    by exact_mod_cast P.Ln_ge⟩

end Bounds

lemma F_ηδ : IsFloat P.p ((P.Nηδ : ℝ) * P.E0) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  apply P.isFloat_grid16 _ (3 * P.Kn + 3 * (P.Kn : ℤ) ^ 3) (by simp only [Nηδ]; ring)
  rw [P.two_pow_p_nat, abs_of_nonneg (by positivity)]
  have h1 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := le_mul_of_one_le_right (by positivity) hD
  have hK2 : (1 : ℤ) ≤ (P.Kn : ℤ) ^ 2 := by nlinarith
  nlinarith [mul_le_mul_of_nonneg_left hK2 (by linarith : (0 : ℤ) ≤ P.Kn)]

lemma F_d : IsFloat P.p ((P.Nd : ℝ) * P.E0) := by rw [← P.dg_eq]; exact P.dg_float

lemma F_m6 : IsFloat P.p (((-6 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) * P.E0) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  apply P.isFloat_grid16 _ (-3 * P.Ln) (by rw [hKL]; ring)
  rw [P.two_pow_p_nat, abs_le]
  have h1 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := le_mul_of_one_le_right (by positivity) hD
  have hK2 : (1 : ℤ) ≤ (P.Kn : ℤ) ^ 2 := by nlinarith
  have hK3 : (P.Kn : ℤ) ≤ (P.Kn : ℤ) ^ 3 := by
    nlinarith [mul_le_mul_of_nonneg_left hK2 (by linarith : (0 : ℤ) ≤ P.Kn)]
  constructor <;> nlinarith

lemma F_m48 : IsFloat P.p (((-48 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) * P.E0) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  apply P.isFloat_grid16 _ (-3 * P.Kn) (by ring)
  rw [P.two_pow_p_nat, abs_le]
  have h1 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := le_mul_of_one_le_right (by positivity) hD
  have hK2 : (1 : ℤ) ≤ (P.Kn : ℤ) ^ 2 := by nlinarith
  constructor <;> nlinarith [mul_le_mul_of_nonneg_left hK2 (by linarith : (0 : ℤ) ≤ P.Kn)]

lemma F_cross (ε : ℤ) (hε : ε = 1 ∨ ε = -1) :
    IsFloat P.p (((ε * (2304 * (P.Kn : ℤ) ^ 3 * P.Dn) - 48 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) * P.E0) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  apply P.isFloat_grid16 _ (ε * (144 * (P.Kn : ℤ) ^ 2 * P.Dn) - 3 * P.Kn) (by ring)
  rw [P.two_pow_p_nat, abs_le]
  have h1 : ((P.Kn : ℤ)) ^ 2 * P.Dn * 8 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := by
    nlinarith [mul_le_mul_of_nonneg_right hK (by positivity : (0 : ℤ) ≤ (P.Kn : ℤ) ^ 2 * P.Dn)]
  have hK2 : (8 : ℤ) * P.Kn ≤ (P.Kn : ℤ) ^ 2 := by nlinarith
  have h2 : (P.Kn : ℤ) ^ 2 ≤ (P.Kn : ℤ) ^ 2 * P.Dn := le_mul_of_one_le_right (by positivity) hD
  rcases hε with rfl | rfl <;> constructor <;> nlinarith

lemma F_g2d (r : ℕ) (hr : r ≤ P.k) :
    IsFloat P.p (((P.Nd - 2304 * r * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) * P.E0) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  have hr' : (r : ℤ) ≤ (P.Kn : ℤ) ^ 2 := by rw [P.k_nat] at hr; exact_mod_cast hr
  have hr0 : (0 : ℤ) ≤ r := by positivity
  apply P.isFloat_grid16 _ (P.Md - 144 * r * P.Kn * P.Dn) (by rw [P.Nd_eq, P.pow_q4]; ring)
  have hMd := P.Md_bound
  rw [P.two_pow_p_nat] at hMd ⊢
  rw [abs_le] at hMd ⊢
  have h1 : (r : ℤ) * P.Kn * P.Dn ≤ (P.Kn : ℤ) ^ 3 * P.Dn := by
    have := mul_le_mul_of_nonneg_right hr' (by positivity : (0 : ℤ) ≤ P.Kn * P.Dn); nlinarith
  have h0 : 0 ≤ (r : ℤ) * P.Kn * P.Dn := by positivity
  simp only [Md] at hMd ⊢
  constructor <;> nlinarith [pow_pos (by linarith : (0 : ℤ) < P.Kn) 3]

lemma F_g2o (c : ℤ) (hc : |c| ≤ P.k) :
    IsFloat P.p (((-6 * (P.Kn : ℤ) ^ 2 - 2304 * c * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) * P.E0) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  have hc' : |c| ≤ (P.Kn : ℤ) ^ 2 := by rw [P.k_nat] at hc; exact_mod_cast hc
  rw [abs_le] at hc'
  apply P.isFloat_grid16 _ (-3 * P.Ln - 144 * c * P.Kn * P.Dn) (by rw [hKL]; ring)
  rw [P.two_pow_p_nat, abs_le]
  have hKD : 0 ≤ (P.Kn : ℤ) * P.Dn := by positivity
  have h1 := mul_le_mul_of_nonneg_right hc'.2 hKD
  have h2 := mul_le_mul_of_nonneg_right hc'.1 hKD
  have h3 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := le_mul_of_one_le_right (by positivity) hD
  constructor <;> nlinarith [pow_pos (by linarith : (0 : ℤ) < P.Kn) 3]


/-! ### The products of the Hadamard rows -/

lemma nine_fourths_E0 : (9 / 4 : ℝ) = ((2304 * (P.Kn : ℤ) ^ 4 * P.Dn : ℤ) : ℝ) * P.E0 := by
  rw [P.nine_fourths, P.K_cast, P.D_cast]; push_cast; ring

lemma ninefourths_Q (i j : Fin P.k) :
    9 / 4 * P.Q i j = ((P.hInt i j * (2304 * (P.Kn : ℤ) ^ 3 * P.Dn) : ℤ) : ℝ) * P.E0 := by
  rw [P.Q_eq, div_eq_mul_one_div (P.hInt i j : ℝ), P.inv_K, P.K_cast, P.D_cast]; push_cast; ring

lemma ninefourths_QQ (r i j : Fin P.k) : 9 / 4 * (P.Q r i * P.Q r j) =
    ((P.hInt r i * P.hInt r j * (2304 * (P.Kn : ℤ) ^ 2 * P.Dn) : ℤ) : ℝ) * P.E0 := by
  have hK : P.K ≠ 0 := P.K_pos.ne'
  rw [P.Q_eq, P.Q_eq, P.nine_fourths]
  push_cast
  rw [← P.K_cast, ← P.D_cast]
  field_simp

lemma float_ninefourths : IsFloat P.p (9 / 4 : ℝ) := by
  have := P.isFloat_nine_mul (-2); norm_num at this ⊢; convert this using 1

lemma float_ninefourths_sign (ε : ℤ) (hε : ε = 1 ∨ ε = -1) (e : ℤ) :
    IsFloat P.p (9 / 4 * (ε * (2 : ℝ) ^ e)) := by
  have h := P.isFloat_nine_mul (-2 + e)
  rw [zpow_add₀ (by norm_num)] at h
  rcases hε with rfl | rfl
  · convert h using 1; norm_num; ring
  · convert h.neg using 1; norm_num; ring

lemma Q_zpow (i j : Fin P.k) : P.Q i j = P.hInt i j * (2 : ℝ) ^ (-(P.q : ℤ)) := by
  rw [P.Q_eq, P.zpow_neg_q]; ring

lemma prod_ll_float (r i j : Fin P.k) :
    IsFloat P.p (3 / 2 * (if r = i then 1 else 0) * (3 / 2 * (if r = j then 1 else 0))) := by
  split_ifs
  · convert P.float_ninefourths using 1; norm_num
  · simp only [mul_zero, mul_one]; exact isFloat_zero _
  · simp only [mul_zero, zero_mul, mul_one]; exact isFloat_zero _
  · simp only [mul_zero]; exact isFloat_zero _

lemma prod_rl_float (r i j : Fin P.k) :
    IsFloat P.p (3 / 2 * P.Q r i * (3 / 2 * (if r = j then 1 else 0))) := by
  split_ifs
  · rw [P.Q_zpow]
    convert P.float_ninefourths_sign _ (P.hInt_cases r i) (-(P.q : ℤ)) using 1; ring
  · simp only [mul_zero]; exact isFloat_zero _

lemma prod_rr_float (r i j : Fin P.k) :
    IsFloat P.p (3 / 2 * P.Q r i * (3 / 2 * P.Q r j)) := by
  rw [P.Q_zpow, P.Q_zpow]
  have hsign : P.hInt r i * P.hInt r j = 1 ∨ P.hInt r i * P.hInt r j = -1 := by
    rcases P.hInt_cases r i with h | h <;> rcases P.hInt_cases r j with h' | h' <;>
      simp [h, h']
  convert P.float_ninefourths_sign _ hsign (-(P.q : ℤ) + -(P.q : ℤ)) using 1
  rw [zpow_add₀ (by norm_num)]; push_cast; ring

lemma citer_lowerEq_trans {A B C : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ}
    (h1 : LowerEq A B) (h2 : LowerEq B C) : LowerEq A C :=
  fun i j hij => (h1 i j hij).trans (h2 i j hij)

/-- One Hadamard row: the step from `6k + r` to `6k + r + 1`. -/
theorem phase2_step {o : Ops} {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ho : IsRN o r P.p)
    (m : ℕ) (hm : m < P.k) (ih : LowerEq (Citer o P.T P.C (6 * P.k + m)) (P.P2 m)) :
    LowerEq (Citer o P.T P.C (6 * P.k + m + 1)) (P.P2 (m + 1)) := by
  intro i j hij
  have hs : 6 * P.k + m < 6 * P.k + P.k := by omega
  rw [P.citer_succ_apply _ hs i j hij, ih i j hij]
  set r' : Fin P.k := ⟨m, hm⟩
  have hrow : (⟨6 * P.k + m, hs⟩ : Fin (6 * P.k + P.k)) = Fin.natAdd (6 * P.k) r' := rfl
  rw [hrow]
  have hr'v : (r' : ℕ) = m := rfl
  have hlt : ∀ l : Fin P.k, ((l : ℕ) < m + 1 ↔ (l : ℕ) < m ∨ l = r') := fun l => by
    constructor
    · intro h; rcases Nat.lt_succ_iff_lt_or_eq.mp h with h | h
      · exact Or.inl h
      · exact Or.inr (Fin.ext h)
    · rintro (h | h)
      · omega
      · rw [h, hr'v]; omega
  revert hij
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> intro hij
  · -- the first group
    rw [T_rl, T_rl, ho.upd _ _ _ (P.prod_ll_float r' i j)]
    simp only [P2, blk_ll, Matrix.of_apply]
    by_cases hij' : i = j
    · subst hij'
      by_cases hi : r' = i
      · subst hi
        simp only [if_true, mul_one, hr'v, lt_irrefl, if_false, lt_add_one]
        rw [show (P.Nd : ℝ) * P.E0 - 3 / 2 * (3 / 2) = (P.Nηδ : ℝ) * P.E0 by
          rw [show (3 / 2 : ℝ) * (3 / 2) = 9 / 4 by norm_num, P.nine_fourths_E0, Nηδ, Nd]
          push_cast; ring]
        exact round_eq_self hr P.F_ηδ
      · simp only [if_true, hi, if_false, mul_zero, sub_zero]
        have hiff : (i : ℕ) < m + 1 ↔ (i : ℕ) < m := by
          rw [hlt]; exact ⟨fun h => h.resolve_right (Ne.symm hi), Or.inl⟩
        simp only [hiff]
        split_ifs
        · exact round_eq_self hr P.F_ηδ
        · exact round_eq_self hr P.F_d
    · have h0 : (3 / 2 : ℝ) * (if r' = i then 1 else 0) * (3 / 2 * (if r' = j then 1 else 0)) = 0 := by
        by_cases hi : r' = i
        · have : r' ≠ j := fun h => hij' (hi.symm.trans h)
          simp [this]
        · simp [hi]
      simp only [hij', if_false, h0, sub_zero]
      exact round_eq_self hr P.F_m6
  · exact absurd hij not_nA_le_cA
  · -- the cross block
    rw [T_rr, T_rl, ho.upd _ _ _ (P.prod_rl_float r' i j)]
    simp only [P2, blk_rl, Matrix.of_apply]
    by_cases hj : r' = j
    · subst hj
      simp only [if_true, mul_one, hr'v, lt_irrefl, if_false, lt_add_one]
      rw [show (((P.hInt r' i * (2304 * (P.Kn : ℤ) ^ 3 * P.Dn) - 48 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) *
          P.E0 - 3 / 2 * P.Q r' i * (3 / 2)) = ((-48 * (P.Kn : ℤ) ^ 2 : ℤ) : ℝ) * P.E0 by
        rw [show (3 / 2 : ℝ) * P.Q r' i * (3 / 2) = 9 / 4 * P.Q r' i by ring, P.ninefourths_Q]
        push_cast; ring]
      exact round_eq_self hr P.F_m48
    · simp only [hj, if_false, mul_zero, sub_zero]
      have hiff : (j : ℕ) < m + 1 ↔ (j : ℕ) < m := by
        rw [hlt]; exact ⟨fun h => h.resolve_right (Ne.symm hj), Or.inl⟩
      simp only [hiff]
      split_ifs
      · exact round_eq_self hr P.F_m48
      · exact round_eq_self hr (P.F_cross _ (P.hInt_cases j i))
  · -- the second group
    rw [T_rr, T_rr, ho.upd _ _ _ (P.prod_rr_float r' i j)]
    simp only [P2, blk_rr, Matrix.of_apply]
    have hprod : 3 / 2 * P.Q r' i * (3 / 2 * P.Q r' j) =
        ((P.hInt r' i * P.hInt r' j * (2304 * (P.Kn : ℤ) ^ 2 * P.Dn) : ℤ) : ℝ) * P.E0 := by
      rw [← P.ninefourths_QQ]; ring
    rw [hprod]
    by_cases hij' : i = j
    · subst hij'
      simp only [if_true]
      have hsq : P.hInt r' i * P.hInt r' i = 1 := by
        rcases P.hInt_cases r' i with h | h <;> rw [h] <;> norm_num
      rw [hsq, show (((P.Nd - 2304 * m * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) * P.E0 -
          ((1 * (2304 * (P.Kn : ℤ) ^ 2 * P.Dn) : ℤ) : ℝ) * P.E0) =
          ((P.Nd - 2304 * ((m + 1 : ℕ) : ℤ) * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) * P.E0 by
        push_cast; ring]
      exact round_eq_self hr (P.F_g2d (m + 1) hm)
    · simp only [hij', if_false]
      rw [show (((-6 * (P.Kn : ℤ) ^ 2 - 2304 * P.cSum m i j * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) *
          P.E0 - ((P.hInt r' i * P.hInt r' j * (2304 * (P.Kn : ℤ) ^ 2 * P.Dn) : ℤ) : ℝ) * P.E0) =
          ((-6 * (P.Kn : ℤ) ^ 2 - 2304 * P.cSum (m + 1) i j * (P.Kn : ℤ) ^ 2 * P.Dn : ℤ) : ℝ) * P.E0 by
        rw [P.cSum_succ m hm]; push_cast; ring]
      exact round_eq_self hr (P.F_g2o _ (le_trans (P.cSum_bound (m + 1) hm i j) (by
        exact_mod_cast hm)))

/-- **The next `k` steps.** -/
theorem phase2 {o : Ops} {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ho : IsRN o r P.p) :
    ∀ m ≤ P.k, LowerEq (Citer o P.T P.C (6 * P.k + m)) (P.P2 m)
  | 0, _ => P.citer_lowerEq_trans (P.phase1 hr ho (6 * P.k) le_rfl) P.P1_eq_P2
  | m + 1, hm => P.phase2_step hr ho m (by omega) (phase2 hr ho m (by omega))

/-- The rounded value after all `k` Hadamard rows is the matrix `Ŝ` of the notes. -/
lemma P2_eq_Shat : LowerEq (P.P2 P.k) P.Shat := by
  intro i j hij
  have hk : (P.k : ℤ) = (P.Kn : ℤ) ^ 2 := by rw [P.k_nat]; push_cast; ring
  have hη := P.η_eq
  have hδ := P.δ_eq
  rw [P.K_cast] at hη hδ
  revert hij
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> intro hij
  · simp only [P2, Shat, blk_ll, Matrix.of_apply, Matrix.sub_apply, Matrix.smul_apply,
      Matrix.one_apply, ones, smul_eq_mul, i.2, if_true]
    split_ifs
    · simp only [Nηδ, a, hη, P.k_real, P.K_cast]; push_cast; ring
    · simp only [hη]; push_cast; ring
  · exact absurd hij not_nA_le_cA
  · simp only [P2, Shat, blk_rl, Matrix.of_apply, Matrix.smul_apply, ones, smul_eq_mul, j.2,
      if_true, hη]
    push_cast; ring
  · simp only [P2, Shat, blk_rr, Matrix.of_apply, Matrix.sub_apply, Matrix.smul_apply,
      Matrix.one_apply, ones, smul_eq_mul]
    split_ifs with hij'
    · simp only [Nd, a, hη, P.k_real, P.K_cast]; push_cast [P.k_nat]; ring
    · rw [P.cSum_full i j hij']; simp only [hη]; push_cast; ring

/-- **The stored trailing matrix is `Ŝ`.** After the `7k` leading steps of the rounded
algorithm, the lower triangle of the trailing block equals that of `Ŝ`. -/
theorem trailing_eq_Shat {o : Ops} {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ho : IsRN o r P.p) :
    LowerEq (Citer o P.T P.C (6 * P.k + P.k)) P.Shat :=
  P.citer_lowerEq_trans (P.phase2 hr ho P.k le_rfl) P.P2_eq_Shat

end Par
end CE
end Cholesky
