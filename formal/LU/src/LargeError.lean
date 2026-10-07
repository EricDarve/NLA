import src.Growth
import src.Condition

/-!
# A right-hand side that produces a large backward error

Binary arithmetic with a `p`-bit significand (`u = 2^{-p}`), rounding to nearest `r`, and
`A = wilk n` with `n ≥ 2`. Let `δ = u 2^{n-3}` (`deltaN`) and `b = e₁ + δ eₙ`.

* The computed factors are exact (`wilk_factors` with `ExactOnW.rounded`): `P = I`,
  `l_ij = -1` below the diagonal, `U = I` except for the last column `u_{i,n} = 2^{i-1}`.
* `fwd_wilk`: the forward loop of the notes (`s = 0; s += L[i,j]*y[j]; y[i] = pb[i] - s`,
  `fwdLoop`) gives `ŷ₁ = 1`, `ŷᵢ = 2^{i-2}` for `2 ≤ i ≤ n`: all partial sums are powers of
  two, and `fl(2^{n-2} + δ) = 2^{n-2}` (`round_last`) loses `δ`.
* `bwd_wilk`: back substitution (`bwdLoop`, with the inner products accumulated left to right)
  gives `x̂₁ = x̂ₙ = 1/2` and `x̂ᵢ = 0` otherwise; these operations are exact.
* `wilk_mul_xhat`, `residual_wilk`: `A x̂ = e₁`, so `r = b - A x̂ = δ eₙ`;
  `residual_computed_exact`: the residual is also computed exactly.
* `eta_wilk`: **`η(x̂) = δ / (‖A‖₂/√2 + √(1 + δ²))`**.
* `norm2_wilk_bounds`: **`⌊n/2⌋ ≤ ‖A‖₂ ≤ ‖A‖_F ≤ n`**; the lower bound comes from the
  bottom-left block of `-1`s (`wilk_block`), whose 2-norm is `⌊n/2⌋` (`norm2_neg_ones`).
  `eta_wilk_bounds` gives the resulting two-sided bounds on `η`, and `η ≤ 1` always
  (`etaFormula_le_one`).

(Indices are 0-based in the code: `e₁` is `Pi.single 0`, `eₙ` is `Pi.single (n - 1)`.)
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

/-! ## The loops -/

/-- `s = 0; for t < m: s = fl(s + fl(l t * y t))`. -/
noncomputable def accum (r : ℝ → ℝ) (l y : ℕ → ℝ) : ℕ → ℝ
  | 0 => 0
  | m + 1 => r (accum r l y m + r (l m * y m))

/-- The forward loop of the notes for a unit lower triangular `L`: `y_i = fl(c_i - s_i)`
with `s_i = ∑_{j<i} l_ij y_j` accumulated in increasing order. `fwdLoop r L c i` holds
`y_0, …, y_{i-1}` (and zeros). -/
noncomputable def fwdLoop (r : ℝ → ℝ) (L : ℕ → ℕ → ℝ) (c : ℕ → ℝ) : ℕ → ℕ → ℝ
  | 0 => fun _ => 0
  | i + 1 => Function.update (fwdLoop r L c i) i (r (c i - accum r (L i) (fwdLoop r L c i) i))

/-- Back substitution `x_i = fl(fl(y_i - U[i, i+1:] @ x[i+1:]) / u_ii)` for
`i = n - 1, …, 0`, with the inner product accumulated left to right. `bwdLoop r U y n m`
holds `x_{n-m}, …, x_{n-1}` (and zeros). -/
noncomputable def bwdLoop (r : ℝ → ℝ) (U : ℕ → ℕ → ℝ) (y : ℕ → ℝ) (n : ℕ) : ℕ → ℕ → ℝ
  | 0 => fun _ => 0
  | m + 1 => Function.update (bwdLoop r U y n m) (n - 1 - m)
      (r (r (y (n - 1 - m) - accum r (fun t => U (n - 1 - m) (n - m + t))
        (fun t => bwdLoop r U y n m (n - m + t)) m) / U (n - 1 - m) (n - 1 - m)))

lemma accum_congr (r : ℝ → ℝ) {l l' y y' : ℕ → ℝ} {m : ℕ}
    (h : ∀ t, t < m → l t = l' t ∧ y t = y' t) : accum r l y m = accum r l' y' m := by
  induction m with
  | zero => rfl
  | succ m ih =>
    simp only [accum]
    rw [ih (fun t ht => h t (by omega)), (h m (by omega)).1, (h m (by omega)).2]

/-! ## The data -/

variable {p : ℕ} {r : ℝ → ℝ}

/-- `δ_n = u 2^{n-3} = 2^{n-3-p}`. -/
noncomputable def deltaN (p n : ℕ) : ℝ := (2 : ℝ) ^ ((n : ℤ) - 3 - p)

lemma deltaN_eq (p n : ℕ) : deltaN p n = unitRoundoff p * 2 ^ ((n : ℤ) - 3) := by
  rw [deltaN, unitRoundoff, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring_nf

lemma deltaN_pos (p n : ℕ) : 0 < deltaN p n := by unfold deltaN; positivity

/-- The entries of the computed `L` (below the diagonal) and `U`, as functions of `ℕ`. -/
def wL (i j : ℕ) : ℝ := if j < i then -1 else if i = j then 1 else 0

noncomputable def wU (n i j : ℕ) : ℝ := if j + 1 = n then 2 ^ i else if i = j then 1 else 0

/-- `P b = b = e₁ + δ eₙ`, as a function of `ℕ`. -/
noncomputable def cW (p n : ℕ) (t : ℕ) : ℝ :=
  (if t = 0 then 1 else 0) + (if t + 1 = n then deltaN p n else 0)

/-- The computed `ŷ`: `1, 1, 2, 4, …, 2^{n-2}`. -/
noncomputable def yv (t : ℕ) : ℝ := if t = 0 then 1 else 2 ^ t / 2

/-- The computed `x̂`: `1/2` in the first and last entries. -/
noncomputable def xv (n t : ℕ) : ℝ := if t = 0 ∨ t + 1 = n then 1 / 2 else 0

lemma two_pow_div_two (t : ℕ) : (2 : ℝ) ^ t / 2 = 2 ^ ((t : ℤ) - 1) := by
  rw [zpow_sub₀ (by norm_num), zpow_natCast, zpow_one]

lemma yv_zero : yv 0 = 1 := if_pos rfl

lemma yv_pos {t : ℕ} (ht : t ≠ 0) : yv t = 2 ^ t / 2 := if_neg ht

lemma cW_zero {n : ℕ} (hn : 2 ≤ n) : cW p n 0 = 1 := by
  unfold cW
  rw [if_pos rfl, if_neg (by omega), add_zero]

lemma cW_pos {n t : ℕ} (ht : t ≠ 0) : cW p n t = if t + 1 = n then deltaN p n else 0 := by
  simp only [cW, if_neg ht, zero_add]

lemma wL_lt {i j : ℕ} (h : j < i) : wL i j = -1 := if_pos h

/-! ## Forward substitution -/

section Forward
variable (hp : 1 ≤ p) (hr : IsRoundNearest p r)
include hp hr

lemma isFloat_half_pow (t : ℕ) : IsFloat p ((2 : ℝ) ^ t / 2) := by
  rw [two_pow_div_two]; exact isFloat_two_zpow hp _

/-- The partial sums of `∑_{j<k} l_ij y_j = -(1 + 1 + 2 + ⋯)` are `-2^{k-1}`. -/
lemma accum_fwd (i : ℕ) {y : ℕ → ℝ} (hy : ∀ t, t < i → y t = yv t) :
    ∀ k, 1 ≤ k → k ≤ i → accum r (wL i) y k = -((2 : ℝ) ^ k / 2)
  | 0, h, _ => absurd h (by norm_num)
  | 1, _, hk => by
    show r (0 + r (wL i 0 * y 0)) = _
    rw [wL_lt (by omega), hy 0 (by omega), yv_zero, mul_one,
      round_eq_self hr (isFloat_one hp).neg, zero_add, round_eq_self hr (isFloat_one hp).neg]
    norm_num
  | k + 2, _, hk => by
    show r (accum r (wL i) y (k + 1) + r (wL i (k + 1) * y (k + 1))) = _
    rw [accum_fwd i hy (k + 1) (by omega) (by omega), wL_lt (by omega), hy (k + 1) (by omega),
      yv_pos (Nat.succ_ne_zero k), neg_one_mul,
      round_eq_self hr (isFloat_half_pow hp hr (k + 1)).neg,
      show -((2 : ℝ) ^ (k + 1) / 2) + -((2 : ℝ) ^ (k + 1) / 2) = -((2 : ℝ) ^ (k + 2) / 2) by
        ring]
    exact round_eq_self hr (isFloat_half_pow hp hr (k + 2)).neg

/-- `fl(2^{n-2} + δ) = 2^{n-2}`: the next float above `2^{n-2}` is `2^{n-2} + 4δ`, and
`δ` is a quarter of the spacing. -/
lemma round_last {n : ℕ} (hn : 2 ≤ n) :
    r (deltaN p n + (2 : ℝ) ^ (n - 1) / 2) = (2 : ℝ) ^ (n - 1) / 2 := by
  have e1 : (2 : ℝ) ^ (n - 1) / 2 = 2 ^ ((n : ℤ) - 2) := by
    rw [two_pow_div_two, Nat.cast_sub (by omega)]; congr 1; push_cast; ring
  rw [e1]
  apply round_eq_of_near hp hr ((n : ℤ) - 2) (isFloat_two_zpow hp _)
  · refine ⟨2 ^ (p - 1), ?_⟩
    rw [Int.cast_pow, Int.cast_ofNat, ← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]
    congr 1
    rw [Nat.cast_sub hp]; push_cast; ring
  · rw [add_sub_cancel_right, abs_of_pos (deltaN_pos p n), deltaN]
    rw [show (2 : ℝ) ^ ((n : ℤ) - 2 + 1 - p) / 2 = 2 ^ ((n : ℤ) - 2 - p) by
      rw [show (n : ℤ) - 2 + 1 - p = ((n : ℤ) - 2 - p) + 1 by ring, zpow_add_one₀ (by norm_num)]
      ring]
    exact zpow_lt_zpow_right₀ (by norm_num) (by omega)
  · rw [abs_of_pos (by have := deltaN_pos p n; positivity)]
    have := deltaN_pos p n
    linarith

/-- **The computed forward substitution**: `ŷ = (1, 1, 2, …, 2^{n-2})`. -/
theorem fwd_wilk (n : ℕ) (hn : 2 ≤ n) :
    ∀ i, i ≤ n → ∀ t, t < i → fwdLoop r wL (cW p n) i t = yv t := by
  intro i
  induction i with
  | zero => intro _ t ht; exact absurd ht (by omega)
  | succ i ih =>
    intro hi t ht
    have ih' := ih (by omega)
    show Function.update (fwdLoop r wL (cW p n) i) i
      (r (cW p n i - accum r (wL i) (fwdLoop r wL (cW p n) i) i)) t = yv t
    by_cases hti : t = i
    · subst hti
      rw [Function.update_self]
      rcases Nat.eq_zero_or_pos t with h0 | hpos
      · subst h0
        show r (cW p n 0 - 0) = yv 0
        rw [sub_zero, cW_zero hn, yv_zero]
        exact round_eq_self hr (isFloat_one hp)
      · rw [accum_fwd hp hr t ih' t hpos le_rfl, sub_neg_eq_add,
          cW_pos (Nat.pos_iff_ne_zero.mp hpos), yv_pos (Nat.pos_iff_ne_zero.mp hpos)]
        by_cases hl : t + 1 = n
        · rw [if_pos hl, show t = n - 1 by omega]
          exact round_last hp hr hn
        · rw [if_neg hl, zero_add]
          exact round_eq_self hr (isFloat_half_pow hp hr t)
    · rw [Function.update_of_ne hti]
      exact ih' t (by omega)

end Forward

lemma wU_last {n i : ℕ} (hn : 1 ≤ n) : wU n i (n - 1) = 2 ^ i := by
  unfold wU; rw [if_pos (by omega)]

lemma wU_not_last {n i j : ℕ} (hj : j + 1 ≠ n) (hij : i ≠ j) : wU n i j = 0 := by
  unfold wU; rw [if_neg hj, if_neg hij]

lemma wU_diag {n i : ℕ} (hi : i + 1 ≠ n) : wU n i i = 1 := by
  unfold wU; rw [if_neg hi, if_pos rfl]

lemma xv_mid {n t : ℕ} (h0 : t ≠ 0) (hl : t + 1 ≠ n) : xv n t = 0 := by
  unfold xv; rw [if_neg (by omega)]

lemma xv_end {n t : ℕ} (h : t = 0 ∨ t + 1 = n) : xv n t = 1 / 2 := by
  unfold xv; rw [if_pos h]

/-! ## Back substitution -/

section Backward
variable (hp : 1 ≤ p) (hr : IsRoundNearest p r)
include hp hr

lemma isFloat_half : IsFloat p ((1 : ℝ) / 2) := by
  rw [one_div, ← _root_.zpow_neg_one]; exact isFloat_two_zpow hp (-1)

/-- The inner product `U[i, i+1:] @ x[i+1:]` of row `i = n - 1 - m` is `2^i · 1/2`, and
every partial sum is exact. -/
lemma accum_bwd {n m : ℕ} (hm : m + 1 ≤ n) {x : ℕ → ℝ}
    (hx : ∀ t, n - m ≤ t → t < n → x t = xv n t) (hm1 : 1 ≤ m) :
    accum r (fun t => wU n (n - 1 - m) (n - m + t)) (fun t => x (n - m + t)) m =
      (2 : ℝ) ^ (n - 1 - m) / 2 := by
  -- the first `m - 1` terms vanish
  have hzero : ∀ q, q < m → accum r (fun t => wU n (n - 1 - m) (n - m + t))
      (fun t => x (n - m + t)) q = 0 := by
    intro q hq
    induction q with
    | zero => rfl
    | succ q ih =>
      show r (accum r _ _ q + r (wU n (n - 1 - m) (n - m + q) * x (n - m + q))) = 0
      rw [ih (by omega), wU_not_last (by omega) (by omega), zero_mul, round_zero hr, add_zero,
        round_zero hr]
  obtain ⟨q, rfl⟩ : ∃ q, m = q + 1 := ⟨m - 1, by omega⟩
  show r (accum r _ _ q + r (wU n (n - 1 - (q + 1)) (n - (q + 1) + q) *
    x (n - (q + 1) + q))) = _
  rw [hzero q (by omega), zero_add, show n - (q + 1) + q = n - 1 by omega, wU_last (by omega),
    hx (n - 1) (by omega) (by omega), xv_end (Or.inr (by omega)), mul_one_div,
    round_eq_self hr (isFloat_half_pow hp hr _)]
  exact round_eq_self hr (isFloat_half_pow hp hr _)

/-- **The computed back substitution**: `x̂ = (1/2, 0, …, 0, 1/2)`. -/
theorem bwd_wilk {n : ℕ} (hn : 2 ≤ n) {y : ℕ → ℝ} (hy : ∀ t, t < n → y t = yv t) :
    ∀ m, m ≤ n → ∀ t, n - m ≤ t → t < n → bwdLoop r (wU n) y n m t = xv n t := by
  intro m
  induction m with
  | zero => intro _ t h1 h2; exact absurd h2 (by omega)
  | succ m ih =>
    intro hm t h1 h2
    have ih' := ih (by omega)
    show Function.update (bwdLoop r (wU n) y n m) (n - 1 - m) _ t = xv n t
    by_cases ht : t = n - 1 - m
    · subst ht
      rw [Function.update_self]
      rcases Nat.eq_zero_or_pos m with h0 | hpos
      · -- the last row: `x_n = fl(fl(2^{n-2}) / 2^{n-1}) = 1/2`
        subst h0
        show r (r (y (n - 1 - 0) - 0) / wU n (n - 1 - 0) (n - 1 - 0)) = xv n (n - 1 - 0)
        rw [Nat.sub_zero, sub_zero, hy (n - 1) (by omega), yv_pos (by omega),
          round_eq_self hr (isFloat_half_pow hp hr _), wU_last (by omega),
          xv_end (Or.inr (by omega))]
        rw [show (2 : ℝ) ^ (n - 1) / 2 / 2 ^ (n - 1) = 1 / 2 by
          field_simp]
        exact round_eq_self hr (isFloat_half hp hr)
      · rw [accum_bwd hp hr hm ih' hpos, hy _ (by omega)]
        by_cases hi0 : n - 1 - m = 0
        · -- the first row: `x_1 = fl(fl(1 - 1/2) / 1) = 1/2`
          rw [hi0, yv_zero, wU_diag (by omega), div_one, xv_end (Or.inl rfl)]
          norm_num
          rw [round_eq_self hr (isFloat_half hp hr)]
          exact round_eq_self hr (isFloat_half hp hr)
        · -- the middle rows: `x_i = fl(fl(2^{i-1} - 2^{i-1}) / 1) = 0`
          rw [yv_pos hi0, sub_self, round_zero hr, zero_div, round_zero hr,
            xv_mid hi0 (by omega)]
    · rw [Function.update_of_ne ht]
      exact ih' t (by omega) h2

end Backward

/-! ## Congruence: the loops read only the entries with indices below `n` -/

lemma fwdLoop_congr (r : ℝ → ℝ) {L L' : ℕ → ℕ → ℝ} {c c' : ℕ → ℝ} {n : ℕ}
    (hL : ∀ i t, t < i → i < n → L i t = L' i t) (hc : ∀ i, i < n → c i = c' i) :
    ∀ i, i ≤ n → fwdLoop r L c i = fwdLoop r L' c' i := by
  intro i
  induction i with
  | zero => intro _; rfl
  | succ i ih =>
    intro hi
    show Function.update _ _ _ = Function.update _ _ _
    rw [ih (by omega), hc i (by omega), accum_congr r (fun t ht => ⟨hL i t ht (by omega), rfl⟩)]

lemma bwdLoop_congr (r : ℝ → ℝ) {U U' : ℕ → ℕ → ℝ} {y y' : ℕ → ℝ} {n : ℕ}
    (hU : ∀ i j, i < n → j < n → U i j = U' i j) (hy : ∀ t, t < n → y t = y' t) :
    ∀ m, m ≤ n → bwdLoop r U y n m = bwdLoop r U' y' n m := by
  intro m
  induction m with
  | zero => intro _; rfl
  | succ m ih =>
    intro hm
    show Function.update _ _ _ = Function.update _ _ _
    rw [ih (by omega), hy _ (by omega), hU _ _ (by omega) (by omega),
      accum_congr r (fun t ht => ⟨hU _ _ (by omega) (by omega), rfl⟩)]

/-! ## The computed solution of the example -/

/-- A matrix as a function of natural-number indices (zero outside). -/
noncomputable def natMat {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (i j : ℕ) : ℝ :=
  if h : i < n ∧ j < n then M ⟨i, h.1⟩ ⟨j, h.2⟩ else 0

/-- A vector as a function of natural-number indices (zero outside). -/
noncomputable def natVec {n : ℕ} (v : Fin n → ℝ) (i : ℕ) : ℝ := if h : i < n then v ⟨i, h⟩ else 0

/-- `b = e₁ + δ eₙ`. -/
noncomputable def bW (p n : ℕ) : Fin n → ℝ := fun i => cW p n i.val

/-- `ŷ`: forward loop with the computed `L̂` and `P b`. -/
noncomputable def yhatW (r : ℝ → ℝ) (p n : ℕ) : ℕ → ℝ :=
  fwdLoop r (natMat (factorL (Ops.rounded r) (wilk n)))
    (natVec (factorP (Ops.rounded r) (wilk n) *ᵥ bW p n)) n

/-- `x̂`: back substitution with the computed `Û` and `ŷ`. -/
noncomputable def xhatW (r : ℝ → ℝ) (p n : ℕ) : Fin n → ℝ := fun i =>
  bwdLoop r (natMat (factorU (Ops.rounded r) (wilk n))) (yhatW r p n) n n i.val

section Solution
variable (hp : 1 ≤ p) (hr : IsRoundNearest p r)
include hp hr

/-- **The computed solution of the example**: the routine computes the exact factors of
`wilk n` (`P = I`), the forward loop gives `ŷ = (1, 1, 2, …, 2^{n-2})`, and back
substitution gives `x̂ = (1/2, 0, …, 0, 1/2)`. -/
theorem solution_wilk {n : ℕ} (hn : 2 ≤ n) :
    (∀ t, t < n → yhatW r p n t = yv t) ∧ xhatW r p n = fun i => xv n i.val := by
  obtain ⟨hP, hL, hU⟩ := wilk_factors (Ops.rounded r) (ExactOnW.rounded hp hr) (by omega : 1 ≤ n)
  have hLn : ∀ i t, t < i → i < n → natMat (factorL (Ops.rounded r) (wilk n)) i t = wL i t := by
    intro i t hti hi
    rw [natMat, dif_pos ⟨hi, by omega⟩, hL, wL_lt hti, if_neg (by
      intro h; have := congrArg Fin.val h; simp at this; omega), if_pos (by
      rw [Fin.lt_iff_val_lt_val]; exact hti)]
  have hc : ∀ i, i < n → natVec (factorP (Ops.rounded r) (wilk n) *ᵥ bW p n) i = cW p n i := by
    intro i hi
    rw [natVec, dif_pos hi, factorP, hP, Matrix.permMatrix_mulVec]
    rfl
  have hy : ∀ t, t < n → yhatW r p n t = yv t := by
    intro t ht
    rw [yhatW, fwdLoop_congr r hLn hc n le_rfl]
    exact fwd_wilk hp hr n hn n le_rfl t ht
  refine ⟨hy, ?_⟩
  have hUn : ∀ i j, i < n → j < n → natMat (factorU (Ops.rounded r) (wilk n)) i j = wU n i j := by
    intro i j hi hj
    rw [natMat, dif_pos ⟨hi, hj⟩, hU]
    simp only [wU, Fin.mk.injEq]
  ext i
  rw [xhatW, bwdLoop_congr r hUn (fun t _ => rfl) n le_rfl]
  exact bwd_wilk hp hr hn hy n le_rfl i.val (by omega) i.2

end Solution

/-! ## Residual and backward error -/

section Eta
variable {n : ℕ}

lemma xv_eq_single (hn : 2 ≤ n) : (fun i : Fin n => xv n i.val) =
    Pi.single ⟨0, by omega⟩ (1 / 2) + Pi.single ⟨n - 1, by omega⟩ (1 / 2) := by
  ext i
  simp only [Pi.add_apply, Pi.single_apply, Fin.ext_iff]
  by_cases h0 : i.val = 0
  · rw [xv_end (Or.inl h0), if_pos h0, if_neg (by omega), add_zero]
  · by_cases hl : i.val = n - 1
    · rw [xv_end (Or.inr (by omega)), if_neg h0, if_pos hl, zero_add]
    · rw [xv_mid h0 (by omega), if_neg h0, if_neg hl, add_zero]

lemma bW_eq_single (hn : 2 ≤ n) : bW p n =
    Pi.single ⟨0, by omega⟩ 1 + Pi.single ⟨n - 1, by omega⟩ (deltaN p n) := by
  ext i
  simp only [bW, cW, Pi.add_apply, Pi.single_apply, Fin.ext_iff]
  congr 1
  by_cases h : i.val + 1 = n
  · rw [if_pos h, if_pos (by omega)]
  · rw [if_neg h, if_neg (by omega)]

/-- **`A x̂ = e₁`**: half the sum of the first and last columns of `A`. -/
theorem wilk_mul_xhat (hn : 2 ≤ n) :
    wilk n *ᵥ (fun i : Fin n => xv n i.val) = Pi.single ⟨0, by omega⟩ 1 := by
  rw [xv_eq_single hn]
  ext i
  simp only [Matrix.mulVec, dotProduct_add, dotProduct_single]
  have hlast : wilk n i ⟨n - 1, by omega⟩ = 1 := by
    simp only [wilk, Matrix.of_apply]
    rw [if_pos (show n - 1 + 1 = n by omega)]
  rw [hlast]
  by_cases h0 : i.val = 0
  · have hi : i = ⟨0, by omega⟩ := Fin.ext h0
    rw [hi]
    have h00 : wilk n ⟨0, by omega⟩ ⟨0, by omega⟩ = 1 := by
      simp only [wilk, Matrix.of_apply]
      rw [if_neg (show ¬ (0 + 1 = n) by omega)]
      simp
    rw [h00, Pi.single_eq_same]; norm_num
  · have hne : i ≠ ⟨0, by omega⟩ := fun h => h0 (congrArg Fin.val h)
    have hi0 : wilk n i ⟨0, by omega⟩ = -1 := by
      simp only [wilk, Matrix.of_apply]
      rw [if_neg (show ¬ (0 + 1 = n) by omega), if_neg hne,
        if_pos (by rw [Fin.lt_iff_val_lt_val]; show 0 < i.val; omega)]
    rw [hi0, Pi.single_eq_of_ne hne]; norm_num

/-- **`r = b - A x̂ = δ eₙ`.** -/
theorem residual_wilk (hn : 2 ≤ n) :
    bW p n - wilk n *ᵥ (fun i : Fin n => xv n i.val) = Pi.single ⟨n - 1, by omega⟩ (deltaN p n) := by
  rw [wilk_mul_xhat hn, bW_eq_single hn, add_sub_cancel_left]

lemma vnorm_two_singles {i j : Fin n} (hij : i ≠ j) (a c : ℝ) :
    vnorm (Pi.single i a + Pi.single j c) = Real.sqrt (a ^ 2 + c ^ 2) := by
  rw [← Real.sqrt_sq (vnorm_nonneg _), vnorm_sq]
  congr 1
  simp only [dotProduct_add, dotProduct_single, Pi.add_apply, Pi.single_eq_same,
    Pi.single_eq_of_ne hij, Pi.single_eq_of_ne (Ne.symm hij)]
  ring

lemma wilk_ne_zero (hn : 1 ≤ n) : wilk n ≠ 0 := by
  intro h
  have := congrFun (congrFun h ⟨0, by omega⟩) ⟨0, by omega⟩
  simp [wilk] at this

/-- **The backward error of the computed solution**:
`η(x̂) = δ / (‖A‖₂/√2 + √(1 + δ²))`. -/
theorem eta_wilk (hp : 1 ≤ p) (hr : IsRoundNearest p r) (hn : 2 ≤ n) :
    eta (wilk n) (bW p n) (xhatW r p n) =
      deltaN p n / (norm2 (wilk n) / Real.sqrt 2 + Real.sqrt (1 + deltaN p n ^ 2)) := by
  have h01 : (⟨0, by omega⟩ : Fin n) ≠ ⟨n - 1, by omega⟩ := by
    intro h; have := congrArg Fin.val h; simp at this; omega
  have hb : bW p n ≠ 0 := by
    intro h
    have := congrFun h ⟨0, by omega⟩
    simp [bW, cW] at this
    split_ifs at this <;> linarith [deltaN_pos p n]
  have hA := norm2_pos_of_ne_zero (wilk_ne_zero (n := n) (by omega))
  rw [eta_eq _ _ _ hA hb, etaFormula, (solution_wilk hp hr hn).2, residual_wilk hn,
    vnorm_single_val, abs_of_pos (deltaN_pos p n), xv_eq_single hn, vnorm_two_singles h01,
    bW_eq_single hn, vnorm_two_singles h01]
  congr 2
  · rw [show ((1 : ℝ) / 2) ^ 2 + (1 / 2) ^ 2 = 1 / 2 by norm_num, Real.sqrt_div' 1
      (by norm_num : (0 : ℝ) ≤ 2), Real.sqrt_one]
    ring
  · rw [one_pow]

end Eta

/-! ## The norm of `A` -/

section Norm
variable {n : ℕ}

lemma sum_range_indicator_lt (m : ℕ) (c : ℝ) :
    ∀ k, m ≤ k → ∑ i ∈ Finset.range k, (if i < m then c else 0) = m * c := by
  intro k hk
  induction k, hk using Nat.le_induction with
  | base =>
    rw [Finset.sum_congr rfl (fun i hi => if_pos (Finset.mem_range.mp hi))]
    simp
  | succ k hk ih => rw [Finset.sum_range_succ, ih, if_neg (by omega), add_zero]

lemma sum_fin_indicator_lt {m : ℕ} (hm : m ≤ n) (c : ℝ) :
    ∑ i : Fin n, (if i.val < m then c else 0) = m * c := by
  rw [Fin.sum_univ_eq_sum_range (fun i => if i < m then c else 0)]
  exact sum_range_indicator_lt m c n hm

lemma sum_fin_indicator_ge {m : ℕ} (hm : m ≤ n) (c : ℝ) :
    ∑ i : Fin n, (if n - m ≤ i.val then c else 0) = m * c := by
  rw [Fin.sum_univ_eq_sum_range (fun i => if n - m ≤ i then c else 0),
    ← Finset.sum_range_reflect]
  rw [Finset.sum_congr rfl (g := fun i => if i < m then c else 0) (fun i hi => by
    have := Finset.mem_range.mp hi
    show (if n - m ≤ n - 1 - i then c else 0) = if i < m then c else 0
    by_cases h : i < m
    · rw [if_pos (by omega), if_pos h]
    · rw [if_neg (by omega), if_neg h])]
  exact sum_range_indicator_lt m c n hm

lemma wilk_sq_le_one (i j : Fin n) : wilk n i j ^ 2 ≤ 1 := by
  simp only [wilk, Matrix.of_apply]
  split_ifs <;> norm_num

/-- `‖A‖_F ≤ n`. -/
theorem frob_wilk_le : frob (wilk n) ≤ n := by
  rw [frob, Real.sqrt_le_left (Nat.cast_nonneg n)]
  calc ∑ i, ∑ j, wilk n i j ^ 2 ≤ ∑ _i : Fin n, ∑ _j : Fin n, (1 : ℝ) :=
        Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => wilk_sq_le_one i j
    _ = (n : ℝ) ^ 2 := by simp [sq]

/-- `⌊n/2⌋ ≤ ‖A‖₂`, from the block of `-1` entries in the bottom-left corner. -/
theorem half_le_norm2_wilk (hn : 2 ≤ n) : ((n / 2 : ℕ) : ℝ) ≤ norm2 (wilk n) := by
  set m := n / 2
  have hm1 : 1 ≤ m := by omega
  have hmn : m ≤ n := by omega
  set v : Fin n → ℝ := fun j => if j.val < m then 1 else 0
  have hv : v ≠ 0 := by
    intro h; have := congrFun h ⟨0, by omega⟩; simp [v] at this; omega
  have hv2 : vnorm v ^ 2 = m := by
    rw [vnorm_sq, dotProduct]
    simp only [v, mul_ite, mul_one, mul_zero]
    rw [Finset.sum_congr rfl (fun j _ => by rw [show (if j.val < m then
      (if j.val < m then (1 : ℝ) else 0) else 0) = if j.val < m then 1 else 0 by
        split_ifs <;> rfl])]
    rw [sum_fin_indicator_lt hmn, mul_one]
  have hrow : ∀ i : Fin n, n - m ≤ i.val → (wilk n *ᵥ v) i = -m := by
    intro i hi
    simp only [Matrix.mulVec, dotProduct, v, mul_ite, mul_one, mul_zero]
    rw [Finset.sum_congr rfl (g := fun j : Fin n => if j.val < m then (-1 : ℝ) else 0)
      (fun j _ => by
        show (if j.val < m then wilk n i j else 0) = if j.val < m then (-1 : ℝ) else 0
        by_cases hj : j.val < m
        · rw [if_pos hj, if_pos hj]
          simp only [wilk, Matrix.of_apply]
          rw [if_neg (by omega), if_neg (by intro h; rw [h] at hi; omega),
            if_pos (by rw [Fin.lt_iff_val_lt_val]; omega)]
        · rw [if_neg hj, if_neg hj])]
    rw [sum_fin_indicator_lt hmn]; ring
  have hAv : (m : ℝ) ^ 2 * vnorm v ^ 2 ≤ vnorm (wilk n *ᵥ v) ^ 2 := by
    rw [hv2, vnorm_sq, dotProduct]
    calc (m : ℝ) ^ 2 * m = ∑ i : Fin n, (if n - m ≤ i.val then (m : ℝ) ^ 2 else 0) := by
          rw [sum_fin_indicator_ge hmn]; ring
      _ ≤ ∑ i : Fin n, (wilk n *ᵥ v) i * (wilk n *ᵥ v) i := by
          apply Finset.sum_le_sum
          intro i _
          split_ifs with h
          · rw [hrow i h]; ring_nf; exact le_rfl
          · exact mul_self_nonneg _
  exact le_norm2_of hv
    (le_of_pow_le_pow_left₀ two_ne_zero (vnorm_nonneg _) (by rw [mul_pow]; exact hAv))

/-- **`⌊n/2⌋ ≤ ‖A‖₂ ≤ ‖A‖_F ≤ n`.** -/
theorem norm2_wilk_bounds (hn : 2 ≤ n) :
    ((n / 2 : ℕ) : ℝ) ≤ norm2 (wilk n) ∧ norm2 (wilk n) ≤ frob (wilk n) ∧ frob (wilk n) ≤ n :=
  ⟨half_le_norm2_wilk hn, norm2_le_frob _, frob_wilk_le⟩

/-- The resulting bounds on the backward error: it grows like `u 2^n / n` while
`u 2^n ≪ n`, and it is at most `1`. -/
theorem eta_wilk_bounds (hp : 1 ≤ p) (hr : IsRoundNearest p r) (hn : 2 ≤ n) :
    deltaN p n / (n / Real.sqrt 2 + Real.sqrt (1 + deltaN p n ^ 2)) ≤
        eta (wilk n) (bW p n) (xhatW r p n) ∧
      eta (wilk n) (bW p n) (xhatW r p n) ≤
        deltaN p n / (((n / 2 : ℕ) : ℝ) / Real.sqrt 2 + 1) ∧
      eta (wilk n) (bW p n) (xhatW r p n) ≤ 1 := by
  have hδ := deltaN_pos p n
  have hs2 : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  have hs1 : 1 ≤ Real.sqrt (1 + deltaN p n ^ 2) := by
    rw [Real.one_le_sqrt]; nlinarith
  obtain ⟨h1, h2, h3⟩ := norm2_wilk_bounds (n := n) hn
  have hpos : 0 < ((n / 2 : ℕ) : ℝ) / Real.sqrt 2 + 1 := by positivity
  have hD : 0 < norm2 (wilk n) / Real.sqrt 2 + Real.sqrt (1 + deltaN p n ^ 2) :=
    add_pos_of_nonneg_of_pos (div_nonneg (norm2_nonneg _) hs2.le) (by linarith)
  rw [eta_wilk hp hr hn]
  refine ⟨?_, ?_, ?_⟩
  · apply div_le_div_of_nonneg_left hδ.le hD
    gcongr
    linarith
  · apply div_le_div_of_nonneg_left hδ.le hpos
    gcongr
  · rw [div_le_one hD]
    have : deltaN p n ≤ Real.sqrt (1 + deltaN p n ^ 2) := by
      rw [Real.le_sqrt (by linarith) (by positivity)]; linarith
    have : 0 ≤ norm2 (wilk n) / Real.sqrt 2 := div_nonneg (norm2_nonneg _) hs2.le
    linarith

/-- The bottom-left `⌊n/2⌋ × ⌊n/2⌋` block of `A` consists of `-1` entries. -/
lemma wilk_block (hn : 2 ≤ n) :
    (wilk n).submatrix (fun i : Fin (n / 2) => (⟨n - n / 2 + i.val, by omega⟩ : Fin n))
      (fun j : Fin (n / 2) => (⟨j.val, by omega⟩ : Fin n)) = Matrix.of fun _ _ => -1 := by
  ext i j
  have hi := i.2
  have hj := j.2
  simp only [Matrix.submatrix_apply, wilk, Matrix.of_apply]
  rw [if_neg (show ¬ (j.val + 1 = n) by omega),
    if_neg (fun h => by have := congrArg Fin.val h; simp at this; omega),
    if_pos (by rw [Fin.lt_iff_val_lt_val]; show j.val < n - n / 2 + i.val; omega)]

/-- The `m × m` matrix of `-1` entries has 2-norm `m`. -/
theorem norm2_neg_ones (m : ℕ) :
    norm2 (Matrix.of fun _ _ => (-1 : ℝ) : Matrix (Fin m) (Fin m) ℝ) = m := by
  set J : Matrix (Fin m) (Fin m) ℝ := Matrix.of fun _ _ => -1
  have hJ : ∀ v : Fin m → ℝ, J *ᵥ v = fun _ => -(∑ j, v j) := by
    intro v; ext i; simp [J, Matrix.mulVec, dotProduct]
  have hJsq : ∀ v : Fin m → ℝ, vnorm (J *ᵥ v) ^ 2 = m * (∑ j, v j) ^ 2 := by
    intro v
    rw [hJ, vnorm_sq]
    simp [dotProduct, sq]
  apply le_antisymm
  · apply norm2_le_of (Nat.cast_nonneg m)
    intro v
    have hcs : (∑ j, v j) ^ 2 ≤ m * ∑ j, v j ^ 2 := by
      have := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun _ => (1 : ℝ)) v
      simpa using this
    have hsq : vnorm (J *ᵥ v) ^ 2 ≤ (m * vnorm v) ^ 2 := by
      rw [hJsq, mul_pow, vnorm_sq]
      have hv : v ⬝ᵥ v = ∑ j, v j ^ 2 := by simp [dotProduct, sq]
      rw [hv]
      have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
      nlinarith
    exact le_of_pow_le_pow_left₀ two_ne_zero (mul_nonneg (Nat.cast_nonneg m) (vnorm_nonneg v)) hsq
  · rcases Nat.eq_zero_or_pos m with h | h
    · subst h; simpa using norm2_nonneg J
    · have hv : (fun _ => (1 : ℝ)) ≠ (0 : Fin m → ℝ) := by
        intro h0; have := congrFun h0 ⟨0, h⟩; simp at this
      apply le_norm2_of hv
      have h1 : vnorm (fun _ : Fin m => (1 : ℝ)) ^ 2 = m := by
        rw [vnorm_sq]; simp [dotProduct]
      have h2 : vnorm (J *ᵥ fun _ => (1 : ℝ)) ^ 2 = (m * vnorm (fun _ : Fin m => (1 : ℝ))) ^ 2 := by
        rw [hJsq, mul_pow, h1]; simp; ring
      exact le_of_eq ((sq_eq_sq₀ (mul_nonneg (Nat.cast_nonneg m) (vnorm_nonneg _))
        (vnorm_nonneg _)).mp h2.symm)

end Norm

/-! ## The residual is computed exactly -/

section ResidualFP
variable {n : ℕ}

lemma natMat_wilk_zero (hn : 2 ≤ n) (i : Fin n) :
    natMat (wilk n) i.val 0 = if i.val = 0 then 1 else -1 := by
  rw [natMat, dif_pos ⟨i.2, by omega⟩]
  simp only [wilk, Matrix.of_apply]
  rw [if_neg (show ¬ (0 + 1 = n) by omega)]
  by_cases h0 : i.val = 0
  · rw [if_pos (Fin.ext h0), if_pos h0]
  · rw [if_neg (fun h => h0 (congrArg Fin.val h)), if_pos (by
      rw [Fin.lt_iff_val_lt_val]; show 0 < i.val; omega), if_neg h0]

lemma natMat_wilk_last (hn : 2 ≤ n) (i : Fin n) : natMat (wilk n) i.val (n - 1) = 1 := by
  rw [natMat, dif_pos ⟨i.2, by omega⟩]
  simp only [wilk, Matrix.of_apply]
  rw [if_pos (show n - 1 + 1 = n by omega)]

variable (hp : 1 ≤ p) (hr : IsRoundNearest p r)
include hp hr

/-- Both nonzero entries of `b` are floats. -/
lemma isFloat_b (hn : 2 ≤ n) (t : ℕ) : IsFloat p (cW p n t) := by
  unfold cW
  by_cases h0 : t = 0
  · rw [if_pos h0, if_neg (by omega), add_zero]; exact isFloat_one hp
  · rw [if_neg h0, zero_add]
    split_ifs
    · exact isFloat_two_zpow hp _
    · exact isFloat_zero p

/-- The partial sums of `∑_j a_ij x̂_j` before the last term are `±1/2`. -/
lemma accum_residual (hn : 2 ≤ n) (i : Fin n) :
    ∀ q, 1 ≤ q → q ≤ n - 1 →
      accum r (natMat (wilk n) i.val) (xv n) q = if i.val = 0 then 1 / 2 else -(1 / 2)
  | 0, h, _ => absurd h (by norm_num)
  | 1, _, _ => by
    show r (0 + r (natMat (wilk n) i.val 0 * xv n 0)) = _
    rw [xv_end (Or.inl rfl), natMat_wilk_zero hn, zero_add]
    split_ifs
    · rw [one_mul, round_eq_self hr (isFloat_half hp hr)]
      exact round_eq_self hr (isFloat_half hp hr)
    · rw [neg_one_mul, round_eq_self hr (isFloat_half hp hr).neg]
      exact round_eq_self hr (isFloat_half hp hr).neg
  | q + 2, _, hq => by
    show r (accum r (natMat (wilk n) i.val) (xv n) (q + 1) +
      r (natMat (wilk n) i.val (q + 1) * xv n (q + 1))) = _
    rw [accum_residual hn i (q + 1) (by omega) (by omega), xv_mid (by omega) (by omega),
      mul_zero, round_zero hr, add_zero]
    split_ifs
    · exact round_eq_self hr (isFloat_half hp hr)
    · exact round_eq_self hr (isFloat_half hp hr).neg

/-- **The residual is computed exactly**: with `A x̂` accumulated left to right,
`fl(bᵢ - fl(∑ⱼ aᵢⱼ x̂ⱼ)) = rᵢ`. -/
theorem residual_computed_exact (hn : 2 ≤ n) (i : Fin n) :
    r (cW p n i.val - accum r (natMat (wilk n) i.val) (xv n) n) =
      (bW p n - wilk n *ᵥ (fun i : Fin n => xv n i.val)) i := by
  have hlast' : ∀ q, q + 1 = n →
      accum r (natMat (wilk n) i.val) (xv n) (q + 1) = if i.val = 0 then 1 else 0 := by
    intro q hq
    show r (accum r (natMat (wilk n) i.val) (xv n) q +
      r (natMat (wilk n) i.val q * xv n q)) = _
    rw [accum_residual hp hr hn i q (by omega) (by omega), show q = n - 1 by omega,
      natMat_wilk_last hn, xv_end (Or.inr (by omega)), one_mul,
      round_eq_self hr (isFloat_half hp hr)]
    split_ifs
    · rw [show (1 : ℝ) / 2 + 1 / 2 = 1 by norm_num]; exact round_eq_self hr (isFloat_one hp)
    · rw [neg_add_cancel]; exact round_zero hr
  have hlast : accum r (natMat (wilk n) i.val) (xv n) n = if i.val = 0 then 1 else 0 := by
    have := hlast' (n - 1) (by omega)
    rwa [Nat.sub_add_cancel (by omega : 1 ≤ n)] at this
  have hexact : (bW p n - wilk n *ᵥ (fun i : Fin n => xv n i.val)) i =
      cW p n i.val - (if i.val = 0 then 1 else 0) := by
    rw [Pi.sub_apply, wilk_mul_xhat hn, Pi.single_apply]
    simp only [bW, Fin.ext_iff]
  rw [hlast, hexact]
  apply round_eq_self hr
  by_cases h0 : i.val = 0
  · rw [if_pos h0, h0, cW_zero hn, sub_self]; exact isFloat_zero p
  · rw [if_neg h0, sub_zero]; exact isFloat_b hp hr hn _

end ResidualFP

end LU
