import src.Ops

/-!
# Binary floating-point arithmetic with rounding to nearest

Binary numbers with a `p`-bit significand and an unbounded exponent range, so that
underflow and overflow cannot occur:

`IsFloat p x ↔ x = m · 2^e` with integers `m, e` and `|m| < 2^p`.

`IsRoundNearest p r` says that `r x` is a float at minimal distance from `x`. Ties may
be broken in any way. We prove

* `exists_roundNearest`: such functions exist (the hypothesis is not vacuous);
* `roundNearest_relErr`: `|r x - x| ≤ u |x|` with `u = 2^{-p}`, so rounded arithmetic
  satisfies the standard model (`stdModel_of_roundNearest`);
* `round_eq_self`, `round_eq_of_near`: the tools used to follow rounded computations
  exactly in the counterexamples.
-/
namespace LU

/-- `x` is a binary floating-point number with a `p`-bit significand. -/
def IsFloat (p : ℕ) (x : ℝ) : Prop := ∃ m e : ℤ, |m| < 2 ^ p ∧ x = m * (2 : ℝ) ^ e

/-- `r` rounds to nearest: `r x` is a float closest to `x`. -/
def IsRoundNearest (p : ℕ) (r : ℝ → ℝ) : Prop :=
  ∀ x, IsFloat p (r x) ∧ ∀ f, IsFloat p f → |r x - x| ≤ |f - x|

/-- The unit roundoff `u = 2^{-p}`. -/
noncomputable def unitRoundoff (p : ℕ) : ℝ := (2 : ℝ) ^ (-(p : ℤ))

variable {p : ℕ}

lemma unitRoundoff_pos (p : ℕ) : 0 < unitRoundoff p := by
  unfold unitRoundoff; positivity

lemma isFloat_zero (p : ℕ) : IsFloat p 0 := ⟨0, 0, by simp, by simp⟩

lemma IsFloat.neg {x : ℝ} (h : IsFloat p x) : IsFloat p (-x) := by
  obtain ⟨m, e, hm, rfl⟩ := h
  exact ⟨-m, e, by simpa using hm, by push_cast; ring⟩

lemma isFloat_neg_iff {x : ℝ} : IsFloat p (-x) ↔ IsFloat p x :=
  ⟨fun h => by simpa using h.neg, fun h => h.neg⟩

/-- `m · 2^e` is a float when `|m| ≤ 2^p` (the case `|m| = 2^p` is `±2^{e+p}`). -/
lemma isFloat_of_abs_le (hp : 1 ≤ p) (m e : ℤ) (hm : |m| ≤ 2 ^ p) :
    IsFloat p (m * (2 : ℝ) ^ e) := by
  rcases lt_or_eq_of_le hm with h | h
  · exact ⟨m, e, h, rfl⟩
  · refine ⟨if 0 ≤ m then 1 else -1, e + p, ?_, ?_⟩
    · have : (1 : ℤ) < 2 ^ p := by
        calc (1 : ℤ) < 2 ^ 1 := by norm_num
          _ ≤ 2 ^ p := pow_le_pow_right₀ (by norm_num) hp
      split_ifs <;> simpa using this
    · have hm' : (m : ℝ) = if 0 ≤ m then (2 : ℝ) ^ p else -(2 : ℝ) ^ p := by
        split_ifs with hs
        · rw [abs_of_nonneg hs] at h; exact_mod_cast h
        · rw [abs_of_neg (not_le.mp hs)] at h
          have : m = -(2 ^ p) := by linarith
          rw [this]; push_cast; ring
      rw [hm', zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0), zpow_natCast]
      split_ifs <;> push_cast <;> ring

/-- A float of magnitude at least `2^b` is a multiple of `2^{b+1-p}`. -/
lemma IsFloat.mem_grid {f : ℝ} (hf : IsFloat p f) (b : ℤ) (hb : (2 : ℝ) ^ b ≤ |f|) :
    ∃ M : ℤ, f = M * (2 : ℝ) ^ (b + 1 - p) := by
  obtain ⟨m, e, hm, rfl⟩ := hf
  have he : b + 1 - p ≤ e := by
    by_contra hlt
    push_neg at hlt
    have h1 : e ≤ b - p := by omega
    have hm' : (|m| : ℝ) < 2 ^ p := by exact_mod_cast hm
    have : |(m : ℝ) * (2 : ℝ) ^ e| < (2 : ℝ) ^ b := by
      rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 2 ^ e)]
      calc |(m : ℝ)| * 2 ^ e < 2 ^ p * 2 ^ e := by
            exact mul_lt_mul_of_pos_right hm' (by positivity)
        _ ≤ 2 ^ p * 2 ^ (b - p) := by
            gcongr
            · norm_num
        _ = 2 ^ b := by
            rw [← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring_nf
    linarith
  refine ⟨m * 2 ^ (e - (b + 1 - p)).toNat, ?_⟩
  push_cast
  rw [mul_assoc, ← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]
  congr 2
  rw [Int.toNat_of_nonneg (by omega)]
  ring

lemma round_eq_self {r : ℝ → ℝ} (hr : IsRoundNearest p r) {x : ℝ} (hx : IsFloat p x) :
    r x = x := by
  have := (hr x).2 x hx
  rw [sub_self, abs_zero] at this
  exact sub_eq_zero.mp (abs_nonpos_iff.mp this)

lemma round_zero {r : ℝ → ℝ} (hr : IsRoundNearest p r) : r 0 = 0 :=
  round_eq_self hr (isFloat_zero p)

/-- A float of magnitude below `2^b` has magnitude at most `2^b - 2^{b-p}`: the largest
float below a power of two is one half-spacing below it. -/
lemma IsFloat.abs_le_of_lt (hp : 1 ≤ p) {g : ℝ} (hg : IsFloat p g) (b : ℤ)
    (hgb : |g| < (2 : ℝ) ^ b) : |g| ≤ (2 : ℝ) ^ b - (2 : ℝ) ^ (b - p) := by
  have h2 : (2 : ℝ) ^ b = 2 ^ p * 2 ^ (b - p) := by
    rw [← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring_nf
  by_cases h : (2 : ℝ) ^ (b - 1) ≤ |g|
  · obtain ⟨M, hM⟩ := hg.mem_grid (b - 1) h
    rw [show b - 1 + 1 - (p : ℤ) = b - p by ring] at hM
    have hpos : (0 : ℝ) < 2 ^ (b - p) := by positivity
    rw [hM, abs_mul, abs_of_pos hpos] at hgb ⊢
    rw [h2] at hgb ⊢
    have hM1 : |(M : ℝ)| < 2 ^ p := (mul_lt_mul_iff_of_pos_right hpos).mp hgb
    have hM2 : |M| < 2 ^ p := by exact_mod_cast hM1
    have hM3 : |M| ≤ 2 ^ p - 1 := by omega
    have hM4 : |(M : ℝ)| ≤ 2 ^ p - 1 := by exact_mod_cast hM3
    nlinarith
  · push_neg at h
    have h3 : (2 : ℝ) ^ (b - p) ≤ 2 ^ (b - 1) :=
      zpow_le_zpow_right₀ (by norm_num) (by omega)
    have h4 : (2 : ℝ) ^ b = 2 * 2 ^ (b - 1) := by
      rw [← zpow_one_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring_nf
    linarith

/-- The rounding tool for the examples: if `|x| ≥ 2^b` and `f` is a float on the grid
`2^{b+1-p} ℤ` within half a spacing of `x`, then every rounding to nearest returns `f`. -/
lemma round_eq_of_near (hp : 1 ≤ p) {r : ℝ → ℝ} (hr : IsRoundNearest p r) (b : ℤ) {x f : ℝ}
    (hf : IsFloat p f) (hfM : ∃ M : ℤ, f = M * (2 : ℝ) ^ (b + 1 - p))
    (hxf : |x - f| < (2 : ℝ) ^ (b + 1 - p) / 2)
    (hbx : (2 : ℝ) ^ b ≤ |x|) : r x = f := by
  set σ := (2 : ℝ) ^ (b + 1 - p)
  have hσ : 0 < σ := by positivity
  have hσ2 : σ / 2 = (2 : ℝ) ^ (b - p) := by
    simp only [σ]
    rw [show b + 1 - (p : ℤ) = 1 + (b - p) by ring, zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]
    ring
  obtain ⟨hg, hmin⟩ := hr x
  have hgx : |r x - x| < σ / 2 := lt_of_le_of_lt (hmin f hf) (by rwa [abs_sub_comm])
  have hgb : (2 : ℝ) ^ b ≤ |r x| := by
    by_contra hlt
    push_neg at hlt
    have h1 := hg.abs_le_of_lt hp b hlt
    have h2 := abs_sub_abs_le_abs_sub x (r x)
    rw [abs_sub_comm] at h2
    linarith
  obtain ⟨M, hM⟩ := hg.mem_grid b hgb
  obtain ⟨N, hN⟩ := hfM
  have hdiff : |r x - f| < σ := by
    calc |r x - f| ≤ |r x - x| + |x - f| := abs_sub_le _ _ _
      _ < σ / 2 + σ / 2 := add_lt_add hgx hxf
      _ = σ := by ring
  rw [hM, hN, ← sub_mul, abs_mul, abs_of_pos hσ] at hdiff
  have h1 : |((M - N : ℤ) : ℝ)| < 1 := by
    push_cast
    exact (mul_lt_iff_lt_one_left hσ).mp hdiff
  have h2 : M = N := by
    have : |M - N| < 1 := by exact_mod_cast h1
    have := abs_lt.mp this
    omega
  rw [hM, hN, h2]

/-- Powers of two are floats. -/
lemma isFloat_two_zpow (hp : 1 ≤ p) (e : ℤ) : IsFloat p ((2 : ℝ) ^ e) :=
  ⟨1, e, by
    have : (1 : ℤ) < 2 ^ p := by
      calc (1 : ℤ) < 2 ^ 1 := by norm_num
        _ ≤ 2 ^ p := pow_le_pow_right₀ (by norm_num) hp
    simpa using this, by simp⟩

lemma isFloat_two_pow (hp : 1 ≤ p) (k : ℕ) : IsFloat p ((2 : ℝ) ^ k) := by
  simpa using isFloat_two_zpow hp (k : ℤ)

lemma isFloat_one (hp : 1 ≤ p) : IsFloat p 1 := by simpa using isFloat_two_zpow hp 0

/-- Multiplying a float by a power of two gives a float. -/
lemma IsFloat.mul_two_zpow {x : ℝ} (hx : IsFloat p x) (k : ℤ) : IsFloat p (x * 2 ^ k) := by
  obtain ⟨m, e, hm, rfl⟩ := hx
  exact ⟨m, e + k, hm, by rw [zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring⟩

/-! ## The standard model -/

lemma binade (x : ℝ) (hx : x ≠ 0) :
    (2 : ℝ) ^ Int.log 2 |x| ≤ |x| ∧ |x| < (2 : ℝ) ^ (Int.log 2 |x| + 1) :=
  ⟨Int.zpow_log_le_self (by norm_num) (abs_pos.mpr hx),
    Int.lt_zpow_succ_log_self (by norm_num) _⟩

/-- Rounding to nearest has relative error at most `u = 2^{-p}`. -/
theorem roundNearest_relErr (hp : 1 ≤ p) {r : ℝ → ℝ} (hr : IsRoundNearest p r) (x : ℝ) :
    |r x - x| ≤ unitRoundoff p * |x| := by
  by_cases hx : x = 0
  · rw [hx, round_zero hr]; simp
  obtain ⟨hlo, hhi⟩ := binade x hx
  set b := Int.log 2 |x|
  set σ := (2 : ℝ) ^ (b + 1 - p)
  have hσ : 0 < σ := by positivity
  set N := round (x / σ)
  have hN : |x / σ - N| ≤ 1 / 2 := abs_sub_round _
  have hfx : |(N : ℝ) * σ - x| ≤ σ / 2 := by
    have : (N : ℝ) * σ - x = -((x / σ - N) * σ) := by field_simp; ring
    rw [this, abs_neg, abs_mul, abs_of_pos hσ]
    nlinarith
  have hσx : σ / 2 ≤ unitRoundoff p * |x| := by
    have : σ / 2 = (2 : ℝ) ^ b * unitRoundoff p := by
      simp only [σ, unitRoundoff]
      rw [show b + 1 - (p : ℤ) = b + 1 + (-(p : ℤ)) by ring, zpow_add₀ (by norm_num),
        zpow_add₀ (by norm_num)]
      ring
    rw [this, mul_comm]
    exact mul_le_mul_of_nonneg_left hlo (unitRoundoff_pos p).le
  have hfloat : IsFloat p ((N : ℝ) * σ) := by
    apply isFloat_of_abs_le hp
    have hxσ : |x / σ| < 2 ^ p := by
      rw [abs_div, abs_of_pos hσ, div_lt_iff₀ hσ]
      calc |x| < (2 : ℝ) ^ (b + 1) := hhi
        _ = 2 ^ p * σ := by
          simp only [σ]
          rw [← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring_nf
    have h1 : |(N : ℝ)| < 2 ^ p + 1 := by
      have := abs_sub_abs_le_abs_sub (N : ℝ) (x / σ)
      rw [abs_sub_comm] at this
      linarith
    have h2 : |N| < 2 ^ p + 1 := by exact_mod_cast h1
    omega
  calc |r x - x| ≤ |(N : ℝ) * σ - x| := (hr x).2 _ hfloat
    _ ≤ σ / 2 := hfx
    _ ≤ unitRoundoff p * |x| := hσx

theorem stdModel_of_roundNearest (hp : 1 ≤ p) {r : ℝ → ℝ} (hr : IsRoundNearest p r) :
    (Ops.rounded r).StdModel (unitRoundoff p) ∧ (Ops.fused r).StdModel (unitRoundoff p) :=
  ⟨Ops.stdModel_rounded (roundNearest_relErr hp hr),
    Ops.stdModel_fused (roundNearest_relErr hp hr)⟩

/-! ## Existence of rounding to nearest -/

lemma exists_nearest_pos (hp : 1 ≤ p) (x : ℝ) (hx : 0 < x) :
    ∃ f, IsFloat p f ∧ ∀ g, IsFloat p g → |f - x| ≤ |g - x| := by
  obtain ⟨hlo, hhi⟩ := binade x (ne_of_gt hx)
  rw [abs_of_pos hx] at hlo hhi
  set b := Int.log 2 x
  set σ := (2 : ℝ) ^ (b + 1 - p)
  have hσ : 0 < σ := by positivity
  have hσb : (2 : ℝ) ^ b = 2 ^ (p - 1 : ℕ) * σ := by
    simp only [σ]
    rw [← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]
    congr 1; push_cast [Nat.cast_sub hp]; ring
  have hσb1 : (2 : ℝ) ^ (b + 1) = 2 ^ p * σ := by
    simp only [σ]
    rw [← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; ring_nf
  set y := x / σ
  have hy : y * σ = x := by simp only [y]; field_simp
  have hylo : (2 : ℝ) ^ (p - 1 : ℕ) ≤ y := by
    rw [le_div_iff₀ hσ]; rw [← hσb]; exact hlo
  have hyhi : y < 2 ^ p := by
    rw [div_lt_iff₀ hσ, ← hσb1]; exact hhi
  set f1 := (⌊y⌋ : ℝ) * σ
  set f2 := (⌈y⌉ : ℝ) * σ
  have hfl : (2 ^ (p - 1 : ℕ) : ℤ) ≤ ⌊y⌋ := by
    apply Int.le_floor.mpr; push_cast; exact hylo
  have hce : ⌈y⌉ ≤ (2 ^ p : ℤ) := by
    apply Int.ceil_le.mpr; push_cast; exact hyhi.le
  have hfl0 : 0 ≤ ⌊y⌋ := le_trans (by positivity) hfl
  have hf1 : IsFloat p f1 := by
    apply isFloat_of_abs_le hp
    rw [abs_of_nonneg hfl0]
    exact le_trans (Int.floor_le_ceil y) hce
  have hf2 : IsFloat p f2 := by
    apply isFloat_of_abs_le hp
    rw [abs_of_nonneg (le_trans hfl0 (Int.floor_le_ceil y))]
    exact hce
  have hf1x : f1 ≤ x := by
    rw [← hy]; exact mul_le_mul_of_nonneg_right (Int.floor_le y) hσ.le
  have hf2x : x ≤ f2 := by
    rw [← hy]; exact mul_le_mul_of_nonneg_right (Int.le_ceil y) hσ.le
  have hf1b : (2 : ℝ) ^ b ≤ f1 := by
    rw [hσb]
    exact mul_le_mul_of_nonneg_right (by exact_mod_cast hfl) hσ.le
  have hf2b : f2 ≤ (2 : ℝ) ^ (b + 1) := by
    rw [hσb1]
    exact mul_le_mul_of_nonneg_right (by exact_mod_cast hce) hσ.le
  -- the closer of the two neighbours
  refine ⟨if x - f1 ≤ f2 - x then f1 else f2, by split_ifs <;> assumption, fun g hg => ?_⟩
  have hbest : ∀ g, IsFloat p g → min (x - f1) (f2 - x) ≤ |g - x| := by
    intro g hg
    by_cases hgb : (2 : ℝ) ^ b ≤ |g|
    · obtain ⟨M, hM⟩ := hg.mem_grid b hgb
      have : M ≤ ⌊y⌋ ∨ ⌈y⌉ ≤ M := by
        by_contra h
        push_neg at h
        have := Int.ceil_le_floor_add_one y
        omega
      rcases this with h | h
      · have : g ≤ f1 := by rw [hM]; exact mul_le_mul_of_nonneg_right (by exact_mod_cast h) hσ.le
        rw [abs_of_nonpos (by linarith)]
        exact le_trans (min_le_left _ _) (by linarith)
      · have : f2 ≤ g := by rw [hM]; exact mul_le_mul_of_nonneg_right (by exact_mod_cast h) hσ.le
        rw [abs_of_nonneg (by linarith)]
        exact le_trans (min_le_right _ _) (by linarith)
    · push_neg at hgb
      have : g < f1 := lt_of_le_of_lt (le_abs_self g) (lt_of_lt_of_le hgb hf1b)
      rw [abs_of_neg (by linarith)]
      exact le_trans (min_le_left _ _) (by linarith)
  have := hbest g hg
  split_ifs with h
  · rw [abs_of_nonpos (by linarith)]
    rw [min_eq_left h] at this; linarith
  · rw [abs_of_nonneg (by linarith)]
    rw [min_eq_right (le_of_lt (not_le.mp h))] at this; linarith

lemma exists_nearest (hp : 1 ≤ p) (x : ℝ) :
    ∃ f, IsFloat p f ∧ ∀ g, IsFloat p g → |f - x| ≤ |g - x| := by
  rcases lt_trichotomy x 0 with hx | hx | hx
  · obtain ⟨f, hf, hmin⟩ := exists_nearest_pos hp (-x) (by linarith)
    refine ⟨-f, hf.neg, fun g hg => ?_⟩
    have := hmin (-g) hg.neg
    rw [show -f - x = -(f - -x) by ring, abs_neg, show g - x = -(-g - -x) by ring, abs_neg]
    exact this
  · exact ⟨0, isFloat_zero p, fun g _ => by rw [hx]; simp⟩
  · exact exists_nearest_pos hp x hx

/-- Rounding to nearest exists (for any tie-breaking rule we may choose one). -/
theorem exists_roundNearest (hp : 1 ≤ p) : ∃ r : ℝ → ℝ, IsRoundNearest p r := by
  choose r hr using exists_nearest hp
  exact ⟨r, fun x => ⟨(hr x).1, fun f hf => (hr x).2 f hf⟩⟩

end LU
