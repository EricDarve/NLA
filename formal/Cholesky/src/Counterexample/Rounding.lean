import src.Counterexample.Lead
import src.Counterexample.Matrix

/-!
# Rounding facts for the counterexample

`IsRN o r p` says that `o` performs each operation with the rounding `r`, either with a
separately rounded product (`Ops.rounded`, as in NumPy) or with a fused multiply-add
(`Ops.fused`). When the product is a float, both variants compute `r (a - b c)`.

We prove that every number stored in `A_{9k}` is a float, that the arithmetic satisfies
the hypotheses of the leading-block lemma, and the two rounding steps of the notes:

* `round_cross`: subtracting `9h/8` from a cross-block entry, whose spacing is `2h`,
  decreases it by `2h`;
* `round_diag`: subtracting `9h/8` from the diagonal value `9/4 + η + δ`, whose spacing is
  `4u`, leaves it unchanged.
-/
namespace Cholesky
namespace CE

/-- The arithmetic rounds every operation with `r`. -/
structure IsRN (o : Ops) (r : ℝ → ℝ) (p : ℕ) : Prop where
  sqrt : ∀ a, o.sqrt a = r (Real.sqrt a)
  div : ∀ a b, o.div a b = r (a / b)
  upd : ∀ x b c, IsFloat p (b * c) → o.upd x b c = r (x - b * c)

lemma isRN_rounded {r : ℝ → ℝ} {p : ℕ} (hr : IsRoundNearest p r) :
    IsRN (Ops.rounded r) r p where
  sqrt _ := rfl
  div _ _ := rfl
  upd x b c hbc := by simp only [Ops.rounded]; rw [round_eq_self hr hbc]

lemma isRN_fused {r : ℝ → ℝ} {p : ℕ} : IsRN (Ops.fused r) r p where
  sqrt _ := rfl
  div _ _ := rfl
  upd _ _ _ _ := rfl

namespace Par
variable (P : Par)

lemma p_pos : 1 ≤ P.p := by have := P.p_ge; omega

/-! ## Grids -/

lemma E0_zpow : P.E0 = (2 : ℝ) ^ (-((P.p + P.q + 2 : ℕ) : ℤ)) := rfl

/-- `N E₀` is a float when `|N| ≤ 2^p`. -/
lemma isFloat_E0 (N : ℤ) (hN : |N| ≤ 2 ^ P.p) : IsFloat P.p (N * P.E0) :=
  isFloat_of_abs_le P.p_pos N _ hN

/-- `4u = 2^{2-p} = 16 K E₀`. -/
lemma four_u : (2 : ℝ) ^ (2 - (P.p : ℤ)) = 16 * P.K * P.E0 := by
  have h1 : (2 : ℝ) ^ (2 - (P.p : ℤ)) = 4 * P.u := by
    rw [u, unitRoundoff, show (2 : ℤ) - P.p = 2 + -(P.p : ℤ) by ring, zpow_add₀ (by norm_num)]
    norm_num
  rw [h1, P.u_eq]; ring

/-- `M · 16 K E₀ = M · 4u` is a float when `|M| ≤ 2^p`. -/
lemma isFloat_grid4 (M : ℤ) (hM : |M| ≤ 2 ^ P.p) : IsFloat P.p (M * (16 * P.K * P.E0)) := by
  rw [← P.four_u]; exact isFloat_of_abs_le P.p_pos M _ hM

/-- A real number written with the integer `Kn = 2^q`, `Dn = 4^e`. -/
def Kn : ℕ := 2 ^ P.q
def Dn : ℕ := 4 ^ P.e
def Ln : ℕ := 2 ^ P.q'

lemma Kn_eq : P.Kn = 8 * P.Ln := by simp only [Kn, Ln, q, pow_add]; norm_num; ring
lemma K_cast : P.K = (P.Kn : ℝ) := by simp [K, Kn]
lemma D_cast : P.D = (P.Dn : ℝ) := by simp [D, Dn]
lemma L_cast : P.L = (P.Ln : ℝ) := by simp [L, Ln]
lemma Ln_ge : 1 ≤ P.Ln := Nat.one_le_two_pow
lemma Dn_ge : 1 ≤ P.Dn := Nat.one_le_pow _ _ (by norm_num)
lemma two_pow_p_nat : (2 : ℤ) ^ P.p = 256 * (P.Kn : ℤ) ^ 3 * P.Dn := by
  have := P.two_pow_p
  rw [P.K_cast, P.D_cast] at this
  exact_mod_cast this
lemma k_nat : P.k = P.Kn ^ 2 := by simp only [k, Kn]; ring

/-! ## The rounding lemma in units of `E₀` -/

/-- Rounding in units of `E₀ = 2^c`: `x = X E₀` lies in the binade starting at `2^a E₀`,
whose spacing is `2^g E₀` with `g = a + 1 - p`, and `M 2^g` is the grid point within half a
spacing. -/
lemma round_units {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (a g : ℕ)
    (hg : (g : ℤ) = a + 1 - P.p) (X M : ℤ) (hM : |M| ≤ 2 ^ P.p)
    (hX : |X - M * 2 ^ g| * 2 < 2 ^ g) (hbin : 2 ^ (a + 1) + 2 ^ g ≤ 2 * |X|) :
    r (X * P.E0) = (M * 2 ^ g : ℤ) * P.E0 := by
  set c : ℤ := -((P.p + P.q + 2 : ℕ) : ℤ)
  have hE : P.E0 = (2 : ℝ) ^ c := rfl
  have hc : (2 : ℝ) ^ c > 0 := by positivity
  have hσ : (2 : ℝ) ^ ((a + c : ℤ) + 1 - P.p) = 2 ^ g * 2 ^ c := by
    rw [← zpow_natCast, ← zpow_add₀ (by norm_num)]; congr 1; omega
  have hfl : IsFloat P.p (((M * 2 ^ g : ℤ) : ℝ) * P.E0) := by
    have := isFloat_of_abs_le P.p_pos M ((g : ℤ) + c) hM
    rw [zpow_add₀ (by norm_num), zpow_natCast, ← mul_assoc, ← hE] at this
    push_cast; exact this
  apply round_eq_of_near hr (a + c) (by
      push_cast at hfl ⊢
      convert hfl using 2)
  · exact ⟨M, by rw [hσ, hE]; push_cast; ring⟩
  · rw [hσ, hE]
    have : (X : ℝ) * 2 ^ c - ((M * 2 ^ g : ℤ) : ℝ) * 2 ^ c = ((X - M * 2 ^ g : ℤ) : ℝ) * 2 ^ c := by
      push_cast; ring
    rw [this, abs_mul, abs_of_pos hc]
    have hX' : (|((X - M * 2 ^ g : ℤ) : ℝ)|) * 2 < 2 ^ g := by exact_mod_cast hX
    nlinarith
  · rw [hσ, zpow_add₀ (by norm_num), zpow_natCast, hE, abs_mul, abs_of_pos hc]
    have hbin' : (2 : ℝ) ^ (a + 1) + 2 ^ g ≤ 2 * |(X : ℝ)| := by exact_mod_cast hbin
    rw [pow_succ] at hbin'
    nlinarith

/-- `(M 2^g) E₀` is a float when `|M| ≤ 2^p`. -/
lemma isFloat_E0_grid (M : ℤ) (g : ℕ) (hM : |M| ≤ 2 ^ P.p) :
    IsFloat P.p (((M * 2 ^ g : ℤ) : ℝ) * P.E0) := by
  have := isFloat_of_abs_le P.p_pos M ((g : ℤ) + -((P.p + P.q + 2 : ℕ) : ℤ)) hM
  rw [zpow_add₀ (by norm_num), zpow_natCast, ← mul_assoc] at this
  push_cast; exact this

lemma Kn_ge : 8 ≤ P.Kn := by rw [P.Kn_eq]; have := P.Ln_ge; omega

lemma two_pow_q_nat : (2 : ℤ) ^ P.q = P.Kn := by simp [Kn]

/-! ## The numbers in units of `E₀` -/

/-- The diagonal value `9/4 + η + δ`. -/
noncomputable def dg : ℝ := 9 / 4 + P.η + P.δ

/-- The integer `Nd` with `9/4 + η + δ = Nd E₀`. -/
def Nd : ℤ := 2304 * (P.Kn : ℤ) ^ 4 * P.Dn + 48 * (P.Kn : ℤ) ^ 2 + 48 * (P.Kn : ℤ) ^ 4

lemma dg_eq : P.dg = (P.Nd : ℝ) * P.E0 := by
  rw [dg, P.nine_fourths, P.η_eq, P.δ_eq, Nd, P.K_cast, P.D_cast]; push_cast; ring

/-- `Nd = 16 K (144 K³ D + 3K + 3K³)`. -/
def Md : ℤ := 144 * (P.Kn : ℤ) ^ 3 * P.Dn + 3 * P.Kn + 3 * (P.Kn : ℤ) ^ 3

lemma Nd_eq : P.Nd = P.Md * 2 ^ (P.q + 4) := by
  rw [Nd, Md, pow_add, P.two_pow_q_nat]; ring

lemma Md_bound : |P.Md| ≤ 2 ^ P.p := by
  rw [P.two_pow_p_nat, Md]
  have hK : (8 : ℤ) ≤ P.Kn := by exact_mod_cast P.Kn_ge
  have hD : (1 : ℤ) ≤ P.Dn := by exact_mod_cast P.Dn_ge
  have h1 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn :=
    le_mul_of_one_le_right (by positivity) hD
  have hK2 : (1 : ℤ) ≤ (P.Kn : ℤ) ^ 2 := by nlinarith
  have h2 : (P.Kn : ℤ) ≤ (P.Kn : ℤ) ^ 3 := by
    nlinarith [mul_le_mul_of_nonneg_left hK2 (by linarith : (0 : ℤ) ≤ P.Kn)]
  rw [abs_of_nonneg (by positivity)]
  nlinarith

theorem dg_float : IsFloat P.p P.dg := by
  rw [P.dg_eq, P.Nd_eq]; exact P.isFloat_E0_grid _ _ P.Md_bound

lemma pow_q4 : (2 : ℤ) ^ (P.q + 4) = 16 * P.Kn := by rw [pow_add, P.two_pow_q_nat]; ring

/-- **The diagonal entries are unchanged**: the spacing near `9/4` is `4u`, and
`9h/8 < 2u`. -/
theorem round_diag {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) : r (P.dg - P.t ^ 2) = P.dg := by
  have hK : (8 : ℤ) ≤ P.Kn := by exact_mod_cast P.Kn_ge
  have hD : (1 : ℤ) ≤ P.Dn := by exact_mod_cast P.Dn_ge
  have hx : P.dg - P.t ^ 2 = ((P.Nd - 9 : ℤ) : ℝ) * P.E0 := by
    rw [P.dg_eq, P.t_sq]; push_cast; ring
  have e1 : (2 : ℤ) ^ (P.p + P.q + 3 + 1) = 4096 * (P.Kn : ℤ) ^ 4 * P.Dn := by
    rw [show P.p + P.q + 3 + 1 = P.p + (P.q + 4) by ring, pow_add, P.two_pow_p_nat, P.pow_q4]
    ring
  rw [hx, P.dg_eq, P.Nd_eq]
  apply P.round_units hr (P.p + P.q + 3) (P.q + 4) (by push_cast; ring) _ _ P.Md_bound
  · rw [← P.Nd_eq, show P.Nd - 9 - P.Nd = -9 by ring, P.pow_q4]
    norm_num; linarith
  · rw [e1, P.pow_q4]
    have h4 : ((P.Kn : ℤ)) ^ 4 ≤ (P.Kn : ℤ) ^ 4 * P.Dn :=
      le_mul_of_one_le_right (by positivity) hD
    have hK2 : (1 : ℤ) ≤ (P.Kn : ℤ) ^ 2 := by nlinarith
    have h5 : (P.Kn : ℤ) ≤ (P.Kn : ℤ) ^ 4 := by
      nlinarith [mul_le_mul_of_nonneg_left hK2 (by positivity : (0 : ℤ) ≤ (P.Kn : ℤ) ^ 2),
        mul_le_mul_of_nonneg_left hK2 (by linarith : (0 : ℤ) ≤ P.Kn)]
    have hpos : 0 ≤ P.Md * (16 * (P.Kn : ℤ)) - 9 := by
      simp only [Md]; nlinarith [pow_pos (by linarith : (0 : ℤ) < P.Kn) 4]
    rw [abs_of_nonneg hpos]
    simp only [Md]
    nlinarith [pow_pos (by linarith : (0 : ℤ) < P.Kn) 4, pow_pos (by linarith : (0 : ℤ) < P.Kn) 2]

/-- The cross-block value after `s` small updates, `±9/(4√k) + η - 2hs`, in units of `E₀`. -/
def Ncross (ε : ℤ) (s : ℕ) : ℤ := ε * (2304 * (P.Kn : ℤ) ^ 3 * P.Dn) + 48 * (P.Kn : ℤ) ^ 2 - 16 * s

lemma cross_eq (ε : ℤ) (s : ℕ) : (9 / 4 : ℝ) * (ε / P.K) + P.η - 2 * P.h * s =
    (P.Ncross ε s : ℝ) * P.E0 := by
  rw [div_eq_mul_one_div (ε : ℝ), P.inv_K, P.η_eq, P.h_eq, Ncross, P.K_cast, P.D_cast]
  push_cast
  ring

/-- **Each cross-block subtraction of `9h/8` rounds to a decrease of `2h`.** -/
theorem round_cross {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ε : ℤ) (hε : ε = 1 ∨ ε = -1)
    (s : ℕ) (hs : s + 1 ≤ 6 * P.k) :
    r ((P.Ncross ε s : ℝ) * P.E0 - P.t ^ 2) = (P.Ncross ε (s + 1) : ℝ) * P.E0 := by
  have hK : (8 : ℤ) ≤ P.Kn := by exact_mod_cast P.Kn_ge
  have hD : (1 : ℤ) ≤ P.Dn := by exact_mod_cast P.Dn_ge
  have hs' : (s : ℤ) + 1 ≤ 6 * (P.Kn : ℤ) ^ 2 := by
    have := P.k_nat; rw [this] at hs; exact_mod_cast hs
  have hs0 : (0 : ℤ) ≤ s := by positivity
  have hx : (P.Ncross ε s : ℝ) * P.E0 - P.t ^ 2 = ((P.Ncross ε s - 9 : ℤ) : ℝ) * P.E0 := by
    rw [P.t_sq]; push_cast; ring
  set M : ℤ := ε * (144 * (P.Kn : ℤ) ^ 3 * P.Dn) + 3 * (P.Kn : ℤ) ^ 2 - s - 1
  have hN : P.Ncross ε (s + 1) = M * 2 ^ 4 := by simp only [Ncross, M]; push_cast; ring
  have h3 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn :=
    le_mul_of_one_le_right (by positivity) hD
  have h32 : 8 * (P.Kn : ℤ) ^ 2 ≤ (P.Kn : ℤ) ^ 3 := by nlinarith [sq_nonneg (P.Kn : ℤ)]
  have hM : |M| ≤ 2 ^ P.p := by
    rw [P.two_pow_p_nat]
    rcases hε with rfl | rfl
    · rw [abs_le]; constructor <;> nlinarith
    · rw [abs_le]; constructor <;> nlinarith
  have e1 : (2 : ℤ) ^ (P.p + 3 + 1) = 4096 * (P.Kn : ℤ) ^ 3 * P.Dn := by
    rw [show P.p + 3 + 1 = P.p + 4 by ring, pow_add, P.two_pow_p_nat]; ring
  have hdiff : P.Ncross ε s - 9 - M * 2 ^ 4 = 7 := by
    rw [← hN]; simp only [Ncross]; push_cast; ring
  rw [hx, hN]
  apply P.round_units hr (P.p + 3) 4 (by push_cast; ring) _ _ hM
  · rw [hdiff]; norm_num
  · rw [e1]
    rcases hε with rfl | rfl
    · have hpos : 0 ≤ P.Ncross 1 s - 9 := by simp only [Ncross]; nlinarith
      rw [abs_of_nonneg hpos]; simp only [Ncross]; nlinarith
    · have hneg : P.Ncross (-1) s - 9 < 0 := by simp only [Ncross]; nlinarith
      rw [abs_of_neg hneg]; simp only [Ncross]; nlinarith

end Par
end CE
end Cholesky
