import src.Float
import src.Hadamard

/-!
# Parameters of the counterexample family

The notes fix integers `p, q` with `q ≥ 3`, `p + q` even, `3q + 8 ≤ p`, and put
`u = 2^{-p}`, `k = 4^q`, `n = 9k`,

`h = 2u/√k`, `η = 6kh = 12√k u`, `δ = kη`, `t = 3 · 2^{-(p+q+2)/2}`.

We write an admissible pair as `q = q' + 3`, `p = 3q + 8 + 2e` (`admissible_iff`). All
quantities are integer multiples of `E₀ = 2^{-(p+q+2)}`; with `K = 2^q` and `D = 4^e`
we have `2^p = 256 K³ D` and, in units of `E₀`,

`h = 8`, `η = 48 K²`, `δ = 48 K⁴`, `t² = 9`, `9/4 = 2304 K⁴ D`, `1/√k = 1024 K³ D`,
`u = 4K`.
-/
namespace Cholesky
namespace CE

/-- The admissibility conditions of the notes. -/
def Admissible (p q : ℕ) : Prop := 3 ≤ q ∧ (p + q) % 2 = 0 ∧ 3 * q + 8 ≤ p

theorem admissible_iff (p q : ℕ) :
    Admissible p q ↔ ∃ q' e : ℕ, q = q' + 3 ∧ p = 3 * q + 8 + 2 * e := by
  constructor
  · rintro ⟨hq, hpar, hp⟩
    refine ⟨q - 3, (p - 3 * q - 8) / 2, by omega, by omega⟩
  · rintro ⟨q', e, rfl, rfl⟩
    refine ⟨by omega, by omega, by omega⟩

/-- `3q + 8 ≤ p` is equivalent to `k^{3/2} u ≤ 1/256`, where `k^{3/2} = 2^{3q}`. -/
theorem size_iff (p q : ℕ) :
    3 * q + 8 ≤ p ↔ (2 : ℝ) ^ (3 * q) * (2 : ℝ) ^ (-(p : ℤ)) ≤ 1 / 256 := by
  rw [← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0),
    show (1 / 256 : ℝ) = (2 : ℝ) ^ (-8 : ℤ) by norm_num,
    zpow_le_zpow_iff_right₀ (by norm_num : (1 : ℝ) < 2)]
  push_cast
  omega

/-- The parameters, written as `q = q' + 3` and `p = 3q + 8 + 2e`. -/
structure Par where
  q' : ℕ
  e : ℕ

namespace Par
variable (P : Par)

def q : ℕ := P.q' + 3
def p : ℕ := 3 * P.q + 8 + 2 * P.e
/-- `k = 4^q`, the size of the Hadamard block. -/
def k : ℕ := 2 ^ (2 * P.q)

lemma admissible : Admissible P.p P.q := (admissible_iff _ _).mpr ⟨P.q', P.e, rfl, rfl⟩

/-- `K = 2^q = √k`. -/
noncomputable def K : ℝ := 2 ^ P.q
/-- `D = 4^e`. -/
noncomputable def D : ℝ := 4 ^ P.e
/-- `L = 2^{q'}`, so that `K = 8 L`. -/
noncomputable def L : ℝ := 2 ^ P.q'
/-- The unit `E₀ = 2^{-(p+q+2)}`. -/
noncomputable def E0 : ℝ := (2 : ℝ) ^ (-((P.p + P.q + 2 : ℕ) : ℤ))
/-- The unit roundoff `u = 2^{-p}`. -/
noncomputable def u : ℝ := unitRoundoff P.p

noncomputable def h : ℝ := 2 * P.u / P.K
noncomputable def η : ℝ := 6 * P.k * P.h
noncomputable def δ : ℝ := P.k * P.η
/-- `t = 3 · 2^{-(p+q+2)/2}`. -/
noncomputable def t : ℝ := 3 * (2 : ℝ) ^ (-(((P.p + P.q + 2) / 2 : ℕ) : ℤ))

lemma K_eq : P.K = 8 * P.L := by
  simp only [K, L, q, pow_add]; norm_num; ring

lemma L_ge : 1 ≤ P.L := one_le_pow₀ (by norm_num)
lemma D_ge : 1 ≤ P.D := one_le_pow₀ (by norm_num)
lemma K_ge : 8 ≤ P.K := by rw [K_eq]; linarith [P.L_ge]
lemma K_pos : 0 < P.K := by linarith [P.K_ge]
lemma D_pos : 0 < P.D := by linarith [P.D_ge]
lemma E0_pos : 0 < P.E0 := by unfold E0; positivity

lemma k_real : (P.k : ℝ) = P.K ^ 2 := by
  simp only [k, K]; push_cast; ring

lemma two_pow_p : (2 : ℝ) ^ P.p = 256 * P.K ^ 3 * P.D := by
  have h1 : (2 : ℝ) ^ (3 * P.q) = (2 ^ P.q) ^ 3 := by rw [← pow_mul, mul_comm]
  have h2 : (2 : ℝ) ^ (2 * P.e) = 4 ^ P.e := by rw [pow_mul]; norm_num
  simp only [p, K, D, pow_add, h1, h2]; norm_num; ring

/-- `E₀ · 1024 K⁴ D = 1`. -/
lemma E0_mul : P.E0 * (1024 * P.K ^ 4 * P.D) = 1 := by
  have h1 : (2 : ℝ) ^ (P.p + P.q + 2) = 1024 * P.K ^ 4 * P.D := by
    rw [pow_add, pow_add, P.two_pow_p]; simp only [K]; ring
  rw [← h1, E0, ← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0), neg_add_cancel,
    zpow_zero]

lemma E0_eq : P.E0 = 1 / (1024 * P.K ^ 4 * P.D) :=
  eq_one_div_of_mul_eq_one_left P.E0_mul

lemma u_eq : P.u = 4 * P.K * P.E0 := by
  have hu : P.u * (2 : ℝ) ^ P.p = 1 := by
    rw [u, unitRoundoff, ← zpow_natCast, ← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]; simp
  rw [P.two_pow_p] at hu
  have hK := P.K_pos; have hD := P.D_pos
  rw [P.E0_eq]
  field_simp
  nlinarith [hu]

lemma h_eq : P.h = 8 * P.E0 := by
  rw [h, P.u_eq]; field_simp [P.K_pos.ne']; ring

lemma η_eq : P.η = 48 * P.K ^ 2 * P.E0 := by
  rw [η, P.h_eq, P.k_real]; ring

lemma δ_eq : P.δ = 48 * P.K ^ 4 * P.E0 := by
  rw [δ, P.η_eq, P.k_real]; ring

lemma t_sq : P.t ^ 2 = 9 * P.E0 := by
  have heven : (P.p + P.q + 2) / 2 * 2 = P.p + P.q + 2 := by
    simp only [p, q]; omega
  set w := (P.p + P.q + 2) / 2
  calc P.t ^ 2 = 9 * ((2 : ℝ) ^ (-(w : ℤ)) * (2 : ℝ) ^ (-(w : ℤ))) := by rw [t]; ring
    _ = 9 * (2 : ℝ) ^ (-(w : ℤ) + -(w : ℤ)) := by rw [← zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0)]
    _ = 9 * P.E0 := by rw [E0]; congr 2; push_cast; omega

/-- `t² = (9/8) h`. -/
theorem t_sq_eq : P.t ^ 2 = 9 / 8 * P.h := by
  rw [P.t_sq, P.h_eq]; ring

lemma nine_fourths : (9 / 4 : ℝ) = 2304 * P.K ^ 4 * P.D * P.E0 := by
  have := P.E0_mul; linarith

/-- `1/√k = 2^{-q} = 1024 K³ D E₀`. -/
lemma inv_K : 1 / P.K = 1024 * P.K ^ 3 * P.D * P.E0 := by
  rw [P.E0_eq]; have := P.K_pos; have := P.D_pos; field_simp

lemma zpow_neg_q : (2 : ℝ) ^ (-(P.q : ℤ)) = 1 / P.K := by
  rw [K, zpow_neg, zpow_natCast, one_div]

/-- `η = 12 √k u`. -/
theorem η_eq_u : P.η = 12 * P.K * P.u := by
  rw [P.η_eq, P.u_eq]; ring

/-- `x = k^{3/2} u = 2^{3q} u`. -/
noncomputable def x : ℝ := P.K ^ 3 * P.u

lemma x_eq : P.x = 4 * P.K ^ 4 * P.E0 := by rw [x, P.u_eq]; ring

/-- `δ = 12 x`. -/
theorem δ_eq_x : P.δ = 12 * P.x := by rw [P.δ_eq, P.x_eq]; ring

/-- `x ≤ 1/256`, i.e. `k^{3/2} u ≤ 1/256`. -/
theorem x_le : P.x ≤ 1 / 256 := by
  rw [P.x_eq, P.E0_eq]
  have hK := P.K_pos; have hD := P.D_ge
  rw [show 4 * P.K ^ 4 * (1 / (1024 * P.K ^ 4 * P.D)) = 1 / (256 * P.D) by field_simp; ring]
  exact one_div_le_one_div_of_le (by norm_num) (by linarith)

lemma x_pos : 0 < P.x := by rw [P.x_eq]; have := P.K_pos; have := P.E0_pos; positivity

lemma u_pos : 0 < P.u := unitRoundoff_pos _

/-- `k u ≤ 1/2048`. -/
lemma ku_le : (P.k : ℝ) * P.u ≤ 1 / 2048 := by
  rw [P.k_real, P.u_eq, P.E0_eq]
  have hK := P.K_ge; have hD := P.D_ge
  rw [show P.K ^ 2 * (4 * P.K * (1 / (1024 * P.K ^ 4 * P.D))) = 1 / (256 * P.K * P.D) by
    field_simp; ring]
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]
  nlinarith

/-- `u ≤ 2^{-17}`. -/
lemma u_le : P.u ≤ 1 / 131072 := by
  rw [P.u_eq, P.E0_eq]
  have hK := P.K_ge; have hD := P.D_ge
  rw [show 4 * P.K * (1 / (1024 * P.K ^ 4 * P.D)) = 1 / (256 * P.K ^ 3 * P.D) by
    field_simp; ring]
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]
  nlinarith [pow_le_pow_left₀ (by norm_num) hK 3]

lemma p_ge : 17 ≤ P.p := by simp only [p, q]; omega

lemma k_ge : 64 ≤ P.k := by
  simp only [k, q]
  calc 64 = 2 ^ (2 * 3) := by norm_num
    _ ≤ 2 ^ (2 * (P.q' + 3)) := Nat.pow_le_pow_right (by norm_num) (by omega)

lemma η_pos : 0 < P.η := by rw [P.η_eq]; have := P.K_pos; have := P.E0_pos; positivity

end Par
end CE
end Cholesky
