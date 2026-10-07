import src.Gamma
import src.Factorization
import src.FactorBounds

/-!
# Backward error of the Cholesky factorization

**Theorem (`thm:backward_error_cholesky`).** In any arithmetic satisfying the standard
model with unit roundoff `u` (in particular rounding to nearest without underflow or
overflow, `Float.lean`), if the in-place algorithm completes on a matrix of size `n`
with `(n + 1) u < 1`, then the computed factor satisfies

`A + E = L̂ L̂ᵀ`, `Eᵀ = E`, `|E| ≤ γ_{n+1} |L̂| |L̂|ᵀ`.

`A` is the symmetric matrix defined by the lower triangle of the input. The proof only
uses that every pivot is positive; it does not need `A` to be SPD. It is by induction
on `n` along the block structure: the first column contributes `γ₂` (square root) and
`γ₁` (division), and the rounded Schur-complement update adds `γ₁ + γ_n + γ₁ γ_n ≤
γ_{n+1}` to the error of the trailing factor.
-/
namespace Cholesky
open Matrix

variable {u : ℝ}

/-! ## The three elementary error estimates -/

lemma sqrt_step_err (hu : 0 ≤ u) (hu2 : 2 * u < 1) {a d δ : ℝ} (ha : 0 < a)
    (hδ : |δ| ≤ u) (hd : d = Real.sqrt a * (1 + δ)) :
    |d ^ 2 - a| ≤ gamma u 2 * d ^ 2 := by
  have hb := abs_le.mp hδ
  have h1 : 0 < 1 + δ := by linarith
  have hd2 : d ^ 2 = a * (1 + δ) ^ 2 := by rw [hd, mul_pow, Real.sq_sqrt ha.le]
  have : d ^ 2 - a = d ^ 2 * (1 - 1 / (1 + δ) ^ 2) := by
    rw [hd2]; field_simp
  rw [this, abs_mul, abs_of_nonneg (sq_nonneg d), mul_comm]
  exact mul_le_mul_of_nonneg_right (abs_one_sub_inv_sq_le hu hu2 hδ) (sq_nonneg d)

lemma div_step_err (hu1 : u < 1) {a d l δ : ℝ} (hd : d ≠ 0) (hδ : |δ| ≤ u)
    (hl : l = a / d * (1 + δ)) : |l * d - a| ≤ gamma u 1 * |l * d| := by
  have hb := abs_le.mp hδ
  have h1 : 0 < 1 + δ := by linarith
  have : l * d - a = l * d * (1 - 1 / (1 + δ)) := by
    rw [hl]; field_simp
  rw [this, abs_mul, mul_comm]
  have h := abs_inv_one_add_sub_one_le hu1 hδ
  rw [abs_sub_comm] at h
  exact mul_le_mul_of_nonneg_right h (abs_nonneg _)

lemma upd_step_err {a b q s P g1 gN gM δ ε : ℝ} (hu1 : u < 1)
    (hδ : |δ| ≤ u) (hε : |ε| ≤ u) (hs : s = (a - b * (1 + ε)) * (1 + δ))
    (hP : 0 ≤ P) (hq : |q| ≤ P) (hF : |q - s| ≤ gN * P) (hgN : 0 ≤ gN)
    (hg1 : gamma u 1 = g1) (hgM1 : u ≤ gM) (hgM : g1 + gN + g1 * gN ≤ gM) :
    |b + q - a| ≤ gM * (|b| + P) := by
  have hb := abs_le.mp hδ
  have h1 : 0 < 1 + δ := by linarith
  have hg1' : |1 - 1 / (1 + δ)| ≤ g1 := by
    rw [← hg1, abs_sub_comm]; exact abs_inv_one_add_sub_one_le hu1 hδ
  have hg1nn : 0 ≤ g1 := le_trans (abs_nonneg _) hg1'
  have ha : a = s / (1 + δ) + b * (1 + ε) := by rw [hs]; field_simp; ring
  have hkey : b + q - a = -(b * ε) + s * (1 - 1 / (1 + δ)) + (q - s) := by
    rw [ha]; ring
  have hsabs : |s| ≤ P + gN * P := by
    calc |s| = |q - (q - s)| := by ring_nf
      _ ≤ |q| + |q - s| := abs_sub _ _
      _ ≤ P + gN * P := add_le_add hq hF
  rw [hkey]
  calc |-(b * ε) + s * (1 - 1 / (1 + δ)) + (q - s)|
      ≤ |-(b * ε)| + |s * (1 - 1 / (1 + δ))| + |q - s| := abs_add_three _ _ _
    _ = |b| * |ε| + |s| * |1 - 1 / (1 + δ)| + |q - s| := by
        rw [abs_neg, abs_mul, abs_mul]
    _ ≤ |b| * u + (P + gN * P) * g1 + gN * P := by
        gcongr
    _ ≤ gM * (|b| + P) := by
        nlinarith [abs_nonneg b, mul_nonneg hg1nn hP, mul_nonneg hgN hP,
          mul_nonneg (abs_nonneg b) (sub_nonneg.mpr hgM1),
          mul_nonneg hP (sub_nonneg.mpr hgM)]

/-! ## Block formulas for `|L| |L|ᵀ` -/

variable {n : ℕ}

lemma absMat_lower {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L) :
    IsLowerTriangular (absMat L) := by
  intro i j hij; simp [hL i j hij]

lemma trailing_absMat (L : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) :
    trailing (absMat L) = absMat (trailing L) := rfl

lemma abs_mul_transpose_le (L : Matrix (Fin n) (Fin n) ℝ) (i j : Fin n) :
    |(L * Lᵀ) i j| ≤ (absMat L * (absMat L)ᵀ) i j := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, absMat_apply]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (le_of_eq ?_)
  simp [abs_mul]

lemma absProd_nonneg (L : Matrix (Fin n) (Fin n) ℝ) (i j : Fin n) :
    0 ≤ (absMat L * (absMat L)ᵀ) i j := by
  rw [Matrix.mul_apply]
  exact Finset.sum_nonneg fun k _ => mul_nonneg (abs_nonneg _) (abs_nonneg _)

/-! ## The theorem -/

lemma sqrt_pos_of_stdModel {o : Ops} (ho : o.StdModel u) (hu1 : u < 1) {a : ℝ} (ha : 0 < a) :
    0 < o.sqrt a := by
  obtain ⟨δ, hδ, h⟩ := ho.sqrt a
  rw [h]
  have := (abs_le.mp hδ).1
  exact mul_pos (Real.sqrt_pos.mpr ha) (by linarith)

/-- The componentwise bound on the lower triangle. -/
theorem backward_error_lower (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u) :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ), ((n + 1 : ℕ) : ℝ) * u < 1 → Succeeds o A →
      ∀ i j, j ≤ i →
        |(factor o A * (factor o A)ᵀ) i j - A i j| ≤
          gamma u (n + 1) * (absMat (factor o A) * (absMat (factor o A))ᵀ) i j
  | 0, _, _, _, i, _, _ => Fin.elim0 i
  | n + 1, A, hnu, h, i, j, hij => by
    rw [succeeds_succ_iff] at h
    obtain ⟨h0, htrail⟩ := h
    have hnR : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have hnu2 : ((n : ℝ) + 2) * u < 1 := by push_cast at hnu; linarith
    have hnu' : ((n + 1 : ℕ) : ℝ) * u < 1 := by push_cast; nlinarith
    have hu1 : u < 1 := by nlinarith
    have hu2 : 2 * u < 1 := by nlinarith
    have ih := backward_error_lower o ho hu (trail o A) hnu' htrail
    have hL : IsLowerTriangular (factor o A) := lowerPart_lower _
    have hLa : IsLowerTriangular (absMat (factor o A)) := absMat_lower hL
    have hgM : u ≤ gamma u (n + 1 + 1) :=
      le_gamma hu (by omega) (by push_cast; linarith)
    -- the first diagonal entry
    obtain ⟨δ0, hδ0, hd⟩ := ho.sqrt (A 0 0)
    have hdpos : 0 < factor o A 0 0 := by
      rw [factor_zero_zero]; exact sqrt_pos_of_stdModel ho hu1 h0
    have hd' : factor o A 0 0 = Real.sqrt (A 0 0) * (1 + δ0) := by rw [factor_zero_zero, hd]
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨i, rfl⟩
    · -- (0, 0)
      have hj : j = 0 := le_antisymm hij (Fin.zero_le _)
      subst hj
      rw [mul_transpose_zero_zero _ hL, mul_transpose_zero_zero _ hLa, absMat_apply, sq_abs]
      calc |factor o A 0 0 ^ 2 - A 0 0| ≤ gamma u 2 * factor o A 0 0 ^ 2 :=
            sqrt_step_err hu hu2 h0 hδ0 hd'
        _ ≤ gamma u (n + 1 + 1) * factor o A 0 0 ^ 2 := by
            gcongr
            exact gamma_mono hu (by omega) (by push_cast; linarith)
    rcases Fin.eq_zero_or_eq_succ j with rfl | ⟨j, rfl⟩
    · -- (i + 1, 0)
      rw [mul_transpose_succ_zero _ hL, mul_transpose_succ_zero _ hLa, absMat_apply,
        absMat_apply, ← abs_mul]
      obtain ⟨δ, hδ, hl⟩ := ho.div (A i.succ 0) (o.sqrt (A 0 0))
      have hl' : factor o A i.succ 0 = A i.succ 0 / factor o A 0 0 * (1 + δ) := by
        rw [factor_succ_zero, factor_zero_zero]; exact hl
      calc |factor o A i.succ 0 * factor o A 0 0 - A i.succ 0|
          ≤ gamma u 1 * |factor o A i.succ 0 * factor o A 0 0| :=
            div_step_err hu1 (ne_of_gt hdpos) hδ hl'
        _ ≤ gamma u (n + 1 + 1) * |factor o A i.succ 0 * factor o A 0 0| := by
            gcongr
            exact gamma_mono hu (by omega) (by push_cast; linarith)
    · -- (i + 1, j + 1)
      have hij' : j ≤ i := Fin.succ_le_succ_iff.mp hij
      rw [mul_transpose_succ_succ _ i j, mul_transpose_succ_succ _ i j, trailing_absMat,
        trailing_factor, absMat_apply, absMat_apply, ← abs_mul]
      obtain ⟨δ, ε, hδ, hε, hs⟩ := ho.upd (A i.succ j.succ) (firstCol o A i) (firstCol o A j)
      have hs' : trail o A i j = (A i.succ j.succ - factor o A i.succ 0 * factor o A j.succ 0 *
          (1 + ε)) * (1 + δ) := by
        rw [trail_apply _ _ _ _ hij', factor_succ_zero, factor_succ_zero]; exact hs
      have hF := ih i j hij'
      have hgN : 0 ≤ gamma u (n + 1) := gamma_nonneg hu (by push_cast; linarith)
      have hsum := gamma_add_le hu (a := 1) (b := n + 1) (by push_cast; linarith)
      rw [show 1 + (n + 1) = n + 1 + 1 by ring] at hsum
      exact upd_step_err hu1 hδ hε hs' (absProd_nonneg _ i j)
        (abs_mul_transpose_le _ i j) hF hgN rfl hgM hsum

/-- **Backward error of Cholesky factorization.** With `A` the symmetric matrix defined by
the lower triangle of the input and `L̂` the computed factor, `E = L̂ L̂ᵀ - A` is symmetric
and `|E| ≤ γ_{n+1} |L̂| |L̂|ᵀ` entrywise. -/
theorem backward_error (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u)
    (A : Matrix (Fin n) (Fin n) ℝ) (hnu : ((n + 1 : ℕ) : ℝ) * u < 1) (h : Succeeds o A) :
    let L := factor o A
    let E := L * Lᵀ - symmOfLower A
    symmOfLower A + E = L * Lᵀ ∧ Eᵀ = E ∧ IsLowerTriangular L ∧ (∀ i, 0 < L i i) ∧
      ∀ i j, |E i j| ≤ gamma u (n + 1) * (absMat L * (absMat L)ᵀ) i j := by
  intro L E
  have hu1 : u < 1 := by
    have : (1 : ℝ) * u ≤ ((n + 1 : ℕ) : ℝ) * u :=
      mul_le_mul_of_nonneg_right (by norm_cast; omega) hu
    linarith
  have hsymmL : (L * Lᵀ)ᵀ = L * Lᵀ := by rw [Matrix.transpose_mul, Matrix.transpose_transpose]
  have hsymmP : (absMat L * (absMat L)ᵀ)ᵀ = absMat L * (absMat L)ᵀ := by
    rw [Matrix.transpose_mul, Matrix.transpose_transpose]
  have hE : Eᵀ = E := by
    simp only [E, Matrix.transpose_sub, hsymmL, symmOfLower_transpose]
  refine ⟨by simp [E], hE, lowerPart_lower _,
    factor_diag_pos o (fun a ha => sqrt_pos_of_stdModel ho hu1 ha) A h, ?_⟩
  have hlow : ∀ i j, j ≤ i → |E i j| ≤ gamma u (n + 1) * (absMat L * (absMat L)ᵀ) i j := by
    intro i j hij
    have := backward_error_lower o ho hu A hnu h i j hij
    simpa [E, symmOfLower, hij] using this
  intro i j
  rcases le_total j i with hij | hij
  · exact hlow i j hij
  · have h1 := hlow j i hij
    rw [← symm_apply_of_transpose hE, ← symm_apply_of_transpose hsymmP] at h1
    exact h1

end Cholesky
