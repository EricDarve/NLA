import src.Float
import src.Examples

/-!
# `A_ε` in double precision

Without pivoting, elimination on `A_ε = [[ε, 1], [1, π]]` computes `l₂₁ = fl(1/ε)` and
`u₂₂ = fl(π - fl(l₂₁ · 1))`. In binary64 (`p = 53`, `u = 2^{-53}`) with `ε = 10⁻¹⁸`, the
stored inputs are `ε̂ = fl(10⁻¹⁸)` and `π̂ = fl(π)`, and

* `Aeps_double`: `u₂₂ = -l₂₁`, exactly the value obtained with `π` replaced by `0`: forming
  `π - 1/ε` loses the contribution of `π`;
* `Aeps_double_product`: the computed factors reconstruct `(L̂ Û)₂₂ = l₂₁ + u₂₂ = 0` instead
  of `π̂`.

The proof locates `l₂₁` in the binade `[2^59, 2^60)`, where floats are spaced `2^7` apart,
and `|π̂| < 2^6`.
-/
namespace LU

/-- The multiplier `l₂₁ = fl(1 / fl(10⁻¹⁸))` is at least `2^59 + 68`. -/
lemma Aeps_double_mult (r : ℝ → ℝ) (hr : IsRoundNearest 53 r) :
    (2 : ℝ) ^ 59 + 68 ≤ (Ops.rounded r).div 1 (r (1 / 10 ^ 18)) := by
  have hp : 1 ≤ 53 := by norm_num
  have hu : unitRoundoff 53 = 1 / 2 ^ 53 := by
    rw [unitRoundoff, zpow_neg, zpow_natCast, one_div]
  have hrel := roundNearest_relErr hp hr
  have h1 := hrel (1 / 10 ^ 18)
  rw [hu, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 10 ^ 18)] at h1
  set εh := r (1 / 10 ^ 18)
  have hε_lo : 1 / 10 ^ 18 * (1 - 1 / 2 ^ 53) ≤ εh := by
    have := (abs_le.mp h1).1; linarith
  have hε_hi : εh ≤ 1 / 10 ^ 18 * (1 + 1 / 2 ^ 53) := by
    have := (abs_le.mp h1).2; linarith
  have hεpos : 0 < εh := lt_of_lt_of_le (by norm_num) hε_lo
  have hinv_lo : 10 ^ 18 / (1 + 1 / 2 ^ 53) ≤ 1 / εh := by
    rw [div_le_div_iff₀ (by norm_num) hεpos]
    nlinarith
  have h2 := hrel (1 / εh)
  rw [hu, abs_of_pos (by positivity : (0 : ℝ) < 1 / εh)] at h2
  have hl_lo : 1 / εh * (1 - 1 / 2 ^ 53) ≤ r (1 / εh) := by
    have := (abs_le.mp h2).1; linarith
  show (2 : ℝ) ^ 59 + 68 ≤ r (1 / εh)
  calc (2 : ℝ) ^ 59 + 68 ≤ 10 ^ 18 / (1 + 1 / 2 ^ 53) * (1 - 1 / 2 ^ 53) := by norm_num
    _ ≤ 1 / εh * (1 - 1 / 2 ^ 53) := mul_le_mul_of_nonneg_right hinv_lo (by norm_num)
    _ ≤ r (1 / εh) := hl_lo

/-- `0 ≤ π̂ ≤ 4`. -/
lemma round_pi_bounds (r : ℝ → ℝ) (hr : IsRoundNearest 53 r) :
    0 ≤ r Real.pi ∧ r Real.pi ≤ 4 := by
  have hu : unitRoundoff 53 = 1 / 2 ^ 53 := by
    rw [unitRoundoff, zpow_neg, zpow_natCast, one_div]
  have h := roundNearest_relErr (by norm_num) hr Real.pi
  rw [hu, abs_of_pos Real.pi_pos] at h
  have h4 := Real.pi_lt_d2
  have h0 := Real.pi_pos
  constructor
  · have := (abs_le.mp h).1; nlinarith
  · have := (abs_le.mp h).2; norm_num at this ⊢; linarith

/-- **Forming `π - 1/ε` loses `π`** (binary64, `ε = 10⁻¹⁸`). -/
theorem Aeps_double (r : ℝ → ℝ) (hr : IsRoundNearest 53 r) :
    (Ops.rounded r).upd (r Real.pi) ((Ops.rounded r).div 1 (r (1 / 10 ^ 18))) 1 =
        -(Ops.rounded r).div 1 (r (1 / 10 ^ 18)) ∧
      (Ops.rounded r).upd 0 ((Ops.rounded r).div 1 (r (1 / 10 ^ 18))) 1 =
        -(Ops.rounded r).div 1 (r (1 / 10 ^ 18)) := by
  have hp : 1 ≤ 53 := by norm_num
  have hbig := Aeps_double_mult r hr
  set l := (Ops.rounded r).div 1 (r (1 / 10 ^ 18))
  have hlf : IsFloat 53 l := (hr _).1
  have hl1 : r (l * 1) = l := by rw [mul_one]; exact round_eq_self hr hlf
  have hlpos : 0 < l := lt_of_lt_of_le (by positivity) hbig
  have h59 : (2 : ℝ) ^ (59 : ℤ) = 2 ^ 59 := zpow_natCast 2 59
  have hgrid : ∃ M : ℤ, -l = M * (2 : ℝ) ^ ((59 : ℤ) + 1 - (53 : ℕ)) :=
    hlf.neg.mem_grid 59 (by rw [abs_neg, abs_of_pos hlpos, h59]; linarith)
  obtain ⟨hπ0, hπ4⟩ := round_pi_bounds r hr
  have h7 : (2 : ℝ) ^ ((59 : ℤ) + 1 - (53 : ℕ)) = 128 := by norm_num
  constructor
  · show r (r Real.pi - r (l * 1)) = -l
    rw [hl1]
    apply round_eq_of_near hp hr 59 hlf.neg hgrid
    · rw [h7, show r Real.pi - l - -l = r Real.pi by ring, abs_of_nonneg hπ0]
      linarith
    · rw [h59, abs_of_neg (by linarith), neg_sub]
      linarith
  · show r (0 - r (l * 1)) = -l
    rw [hl1, zero_sub]
    exact round_eq_self hr hlf.neg

/-- The computed factors reconstruct `(L̂ Û)₂₂ = l₂₁ · 1 + 1 · u₂₂ = 0`, not `π̂`. -/
theorem Aeps_double_product (r : ℝ → ℝ) (hr : IsRoundNearest 53 r) :
    (Ops.rounded r).div 1 (r (1 / 10 ^ 18)) * 1 +
      1 * (Ops.rounded r).upd (r Real.pi) ((Ops.rounded r).div 1 (r (1 / 10 ^ 18))) 1 = 0 := by
  rw [(Aeps_double r hr).1]
  ring

end LU
