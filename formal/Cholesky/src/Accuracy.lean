import src.Solve
import src.SPD

/-!
# Backward and forward error of the computed solution

* `eta A b x` is the normwise backward error: the smallest `ε ≥ 0` such that
  `(A + ΔA) x = b + Δb` with `‖ΔA‖₂ ≤ ε ‖A‖₂` and `‖Δb‖₂ ≤ ε ‖b‖₂`.
* `rigal_gaches`: for `b ≠ 0` it equals `‖b - A x‖₂ / (‖A‖₂ ‖x‖₂ + ‖b‖₂)` and is
  attained.
* `eta_le_of_perturbation`: the pair `(ΔA, 0)` gives `η(x̂) ≤ ‖ΔA‖₂ / ‖A‖₂`.
* `eta_solve_le`: for the computed Cholesky solution `η(x̂) ≤ β_n = n γ_{3n+1} / (1 - γ_{n+1})`;
  `beta_asymp`: `β_n ~ n (3n + 1) u`.
* `forward_error`: `‖x̂ - x‖ / ‖x‖ ≤ 2 κ β / (1 - κ β)` when `η(x̂) ≤ β` and `κ₂(A) β < 1`.
* `cond2_spd`: for SPD `A`, `κ₂(A) = λ_max / λ_min`.
-/
open Filter Topology

namespace Cholesky
open Matrix

variable {n : ℕ}

/-- The relative sizes `ε` of perturbations `(ΔA, Δb)` that make `x` an exact solution. -/
def bwdErrSet (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) : Set ℝ :=
  {ε | 0 ≤ ε ∧ ∃ (ΔA : Matrix (Fin n) (Fin n) ℝ) (Δb : Fin n → ℝ),
      (A + ΔA) *ᵥ x = b + Δb ∧ norm2 ΔA ≤ ε * norm2 A ∧ vnorm Δb ≤ ε * vnorm b}

/-- The normwise backward error `η(x)`. -/
noncomputable def eta (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) : ℝ :=
  sInf (bwdErrSet A b x)

/-- The residual formula for `η`. -/
noncomputable def etaFormula (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) : ℝ :=
  vnorm (b - A *ᵥ x) / (norm2 A * vnorm x + vnorm b)

lemma vecMulVec_mulVec' (a c v : Fin n → ℝ) : vecMulVec a c *ᵥ v = (c ⬝ᵥ v) • a := by
  ext i
  simp [vecMulVec, Matrix.mulVec, dotProduct, Finset.mul_sum, mul_comm, mul_left_comm]

lemma norm2_vecMulVec_le (a c : Fin n → ℝ) : norm2 (vecMulVec a c) ≤ vnorm a * vnorm c := by
  apply norm2_le_of (mul_nonneg (vnorm_nonneg a) (vnorm_nonneg c))
  intro v
  rw [vecMulVec_mulVec', vnorm_smul]
  calc |c ⬝ᵥ v| * vnorm a ≤ vnorm c * vnorm v * vnorm a :=
        mul_le_mul_of_nonneg_right (abs_dot_le c v) (vnorm_nonneg a)
    _ = vnorm a * vnorm c * vnorm v := by ring

/-- **Rigal–Gaches.** For `b ≠ 0` the backward error is given by the residual, and the
minimum is attained. -/
theorem rigal_gaches (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) (hb : b ≠ 0) :
    IsLeast (bwdErrSet A b x) (etaFormula A b x) := by
  unfold etaFormula
  have hbpos := vnorm_pos hb
  set r := b - A *ᵥ x with hr
  set D := norm2 A * vnorm x + vnorm b with hDdef
  have hD : 0 < D := by
    have := mul_nonneg (norm2_nonneg A) (vnorm_nonneg x); linarith
  have hb' : b = A *ᵥ x + r := by rw [hr]; abel
  constructor
  · refine ⟨div_nonneg (vnorm_nonneg _) hD.le, ?_⟩
    by_cases hx : x = 0
    · refine ⟨0, -r, ?_, ?_, ?_⟩
      · rw [hx, Matrix.mulVec_zero, hb', hx, Matrix.mulVec_zero, zero_add, add_neg_cancel]
      · rw [norm2_zero]
        exact mul_nonneg (div_nonneg (vnorm_nonneg _) hD.le) (norm2_nonneg A)
      · have : D = vnorm b := by rw [hDdef, hx, vnorm_zero, mul_zero, zero_add]
        rw [vnorm_neg, this, div_mul_cancel₀ _ (ne_of_gt hbpos)]
    · have hs := vnorm_pos hx
      have hDb : norm2 A * vnorm x = D - vnorm b := by rw [hDdef]; ring
      have key : norm2 A / (D * vnorm x) * vnorm x ^ 2 = 1 - vnorm b / D := by
        calc norm2 A / (D * vnorm x) * vnorm x ^ 2 = norm2 A * vnorm x / D := by
              field_simp
          _ = (D - vnorm b) / D := by rw [hDb]
          _ = 1 - vnorm b / D := by field_simp
      refine ⟨(norm2 A / (D * vnorm x)) • vecMulVec r x, -(vnorm b / D) • r, ?_, ?_, ?_⟩
      · rw [Matrix.add_mulVec, Matrix.smul_mulVec, vecMulVec_mulVec', ← vnorm_sq, smul_smul, key]
        conv_rhs => arg 1; rw [hb']
        rw [sub_smul, one_smul, neg_smul]
        abel
      · rw [norm2_smul, abs_of_nonneg (div_nonneg (norm2_nonneg A) (by positivity))]
        calc norm2 A / (D * vnorm x) * norm2 (vecMulVec r x)
            ≤ norm2 A / (D * vnorm x) * (vnorm r * vnorm x) :=
              mul_le_mul_of_nonneg_left (norm2_vecMulVec_le r x)
                (div_nonneg (norm2_nonneg A) (by positivity))
          _ = vnorm r / D * norm2 A := by field_simp
      · rw [vnorm_smul, abs_neg, abs_of_nonneg (div_nonneg (vnorm_nonneg b) hD.le)]
        exact le_of_eq (by ring)
  · rintro ε ⟨hε, ΔA, Δb, hx, hΔA, hΔb⟩
    have hrr : r = ΔA *ᵥ x - Δb := by
      have : b = (A + ΔA) *ᵥ x - Δb := by rw [hx]; abel
      rw [hr, this, Matrix.add_mulVec]; abel
    have : vnorm r ≤ ε * D := by
      rw [hrr]
      calc vnorm (ΔA *ᵥ x - Δb) ≤ vnorm (ΔA *ᵥ x) + vnorm Δb := vnorm_sub_le _ _
        _ ≤ norm2 ΔA * vnorm x + vnorm Δb := add_le_add (vnorm_mulVec_le _ _) le_rfl
        _ ≤ ε * norm2 A * vnorm x + ε * vnorm b :=
            add_le_add (mul_le_mul_of_nonneg_right hΔA (vnorm_nonneg x)) hΔb
        _ = ε * D := by rw [hDdef]; ring
    rwa [div_le_iff₀ hD]

theorem eta_eq (A : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ) (hb : b ≠ 0) :
    eta A b x = etaFormula A b x :=
  (rigal_gaches A b x hb).csInf_eq

/-- The perturbation `(ΔA, 0)` bounds the backward error. -/
theorem eta_le_of_perturbation (A ΔA : Matrix (Fin n) (Fin n) ℝ) (b x : Fin n → ℝ)
    (hb : b ≠ 0) (hA : 0 < norm2 A) (hΔ : (A + ΔA) *ᵥ x = b) :
    eta A b x ≤ norm2 ΔA / norm2 A := by
  apply csInf_le (rigal_gaches A b x hb).bddBelow
  refine ⟨div_nonneg (norm2_nonneg _) hA.le, ΔA, 0, by rw [hΔ, add_zero], ?_, ?_⟩
  · rw [div_mul_cancel₀ _ (ne_of_gt hA)]
  · rw [vnorm_zero]; exact mul_nonneg (div_nonneg (norm2_nonneg _) hA.le) (vnorm_nonneg _)

/-- `β_n = n γ_{3n+1} / (1 - γ_{n+1})`. -/
noncomputable def betaN (u : ℝ) (n : ℕ) : ℝ := n * gamma u (3 * n + 1) / (1 - gamma u (n + 1))

/-- **Backward error of a Cholesky solve.** `η(x̂) ≤ β_n`. -/
theorem eta_solve_le (o : Ops) {u : ℝ} (ho : o.StdModel u) (hu : 0 ≤ u)
    (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) (hb : b ≠ 0)
    (hA : 0 < norm2 (symmOfLower A)) (hn : 1 ≤ n) (hnu : ((3 * n + 1 : ℕ) : ℝ) * u < 1)
    (h : Succeeds o A) :
    eta (symmOfLower A) b (choleskySolve o A b) ≤ betaN u n := by
  have hnR : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hn1 : ((n + 1 : ℕ) : ℝ) * u < 1 := by push_cast at hnu ⊢; nlinarith
  have h2n1 : 2 * (((n + 1 : ℕ) : ℝ) * u) < 1 := by push_cast at hnu ⊢; nlinarith
  have hγ1 : gamma u (n + 1) < 1 := by
    unfold gamma
    rw [div_lt_one (by linarith)]; linarith
  have hγ3 : 0 ≤ gamma u (3 * n + 1) := gamma_nonneg hu hnu
  obtain ⟨ΔA, hΔ, hΔb⟩ := solve_backward_error o ho hu A b hnu h
  have hc := gamma_solve_le hu n hnu
  set L := factor o A
  have hdom : ∀ i j, |ΔA i j| ≤ (gamma u (3 * n + 1) • (absMat L * (absMat L)ᵀ)) i j := by
    intro i j
    rw [Matrix.smul_apply, smul_eq_mul]
    exact le_trans (hΔb i j) (mul_le_mul_of_nonneg_right hc (absProd_nonneg L i j))
  have h1 : norm2 ΔA ≤ gamma u (3 * n + 1) * frob L ^ 2 := by
    calc norm2 ΔA ≤ norm2 (gamma u (3 * n + 1) • (absMat L * (absMat L)ᵀ)) :=
          norm2_le_of_abs_le hdom
      _ = gamma u (3 * n + 1) * norm2 (absMat L * (absMat L)ᵀ) := by
          rw [norm2_smul, abs_of_nonneg hγ3]
      _ ≤ gamma u (3 * n + 1) * frob L ^ 2 :=
          mul_le_mul_of_nonneg_left (norm2_absProd_le_frob_sq L) hγ3
  have h2 := frob_sq_le o ho hu A hn1 h hγ1
  have h3 : (symmOfLower A).trace ≤ n * norm2 (symmOfLower A) := by
    simpa using trace_le (symmOfLower A)
  calc eta (symmOfLower A) b (choleskySolve o A b) ≤ norm2 ΔA / norm2 (symmOfLower A) :=
        eta_le_of_perturbation _ _ _ _ hb hA hΔ
    _ ≤ gamma u (3 * n + 1) * (n * norm2 (symmOfLower A) / (1 - gamma u (n + 1))) /
          norm2 (symmOfLower A) := by
        apply div_le_div_of_nonneg_right _ hA.le
        exact le_trans h1 (mul_le_mul_of_nonneg_left
          (le_trans h2 (div_le_div_of_nonneg_right h3 (by linarith))) hγ3)
    _ = betaN u n := by
        unfold betaN
        have : 1 - gamma u (n + 1) ≠ 0 := by linarith
        field_simp

/-- `β_n = n (3n + 1) u (1 + o(1))` as `u → 0`. -/
theorem beta_asymp (n : ℕ) (hn : 0 < n) :
    Tendsto (fun u : ℝ => betaN u n / (n * (3 * n + 1) * u)) (𝓝[>] 0) (𝓝 1) := by
  have hc : Tendsto (fun u : ℝ => (1 - (n + 1) * u) / ((1 - (3 * n + 1) * u) *
      (1 - 2 * (n + 1) * u))) (𝓝 0) (𝓝 ((1 - (n + 1) * 0) / ((1 - (3 * n + 1) * 0) *
      (1 - 2 * (n + 1) * 0)))) := by
    apply Tendsto.div
    · exact tendsto_const_nhds.sub (tendsto_id.const_mul _)
    · exact (tendsto_const_nhds.sub (tendsto_id.const_mul _)).mul
        (tendsto_const_nhds.sub (tendsto_id.const_mul _))
    · simp
  simp only [mul_zero, sub_zero, div_one, mul_one] at hc
  refine (hc.mono_left nhdsWithin_le_nhds).congr' ?_
  have hev : ∀ᶠ u in 𝓝[>] (0 : ℝ), 0 < u ∧ 4 * (3 * (n : ℝ) + 1) * u < 1 := by
    have h1 : ∀ᶠ u in 𝓝[>] (0 : ℝ), 0 < u := self_mem_nhdsWithin
    have h2 : ∀ᶠ u in 𝓝[>] (0 : ℝ), 4 * (3 * (n : ℝ) + 1) * u < 1 := by
      apply nhdsWithin_le_nhds
      have : Tendsto (fun u : ℝ => 4 * (3 * (n : ℝ) + 1) * u) (𝓝 0) (𝓝 0) := by
        simpa using ((tendsto_id : Tendsto (fun u : ℝ => u) (𝓝 0) (𝓝 0)).const_mul
          (4 * (3 * (n : ℝ) + 1)))
      exact this (Iio_mem_nhds (by norm_num))
    exact h1.and h2
  filter_upwards [hev] with u ⟨hu, hu4⟩
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have e1 : (1 : ℝ) - (n + 1) * u ≠ 0 := by nlinarith
  have e2 : (1 : ℝ) - (3 * n + 1) * u ≠ 0 := by nlinarith
  have e3 : (1 : ℝ) - 2 * (n + 1) * u ≠ 0 := by nlinarith
  have e4 : (1 : ℝ) - (n + 1) * u / (1 - (n + 1) * u) ≠ 0 := by
    rw [show (1 : ℝ) - (n + 1) * u / (1 - (n + 1) * u) = (1 - 2 * (n + 1) * u) /
      (1 - (n + 1) * u) by field_simp; ring]
    exact div_ne_zero e3 e1
  unfold betaN gamma
  push_cast
  field_simp
  ring

/-! ## Forward error -/

/-- **Forward error.** If `A x = b`, `η(x̂) ≤ β` and `κ₂(A) β < 1`, then
`‖x̂ - x‖ / ‖x‖ ≤ 2 κ₂(A) β / (1 - κ₂(A) β)`. -/
theorem forward_error (A : Matrix (Fin n) (Fin n) ℝ) (hA : IsUnit A.det)
    (x b xh : Fin n → ℝ) (hx : A *ᵥ x = b) (hb : b ≠ 0) (β : ℝ)
    (hη : eta A b xh ≤ β) (hκ : cond2 A * β < 1) :
    vnorm (xh - x) / vnorm x ≤ 2 * cond2 A * β / (1 - cond2 A * β) := by
  have hx0 : x ≠ 0 := by rintro rfl; apply hb; rw [← hx, Matrix.mulVec_zero]
  have hxpos := vnorm_pos hx0
  have hres : vnorm (b - A *ᵥ xh) ≤ β * (norm2 A * vnorm xh + vnorm b) := by
    rw [eta_eq A b xh hb] at hη
    have hD : 0 < norm2 A * vnorm xh + vnorm b := by
      have := vnorm_pos hb; have := mul_nonneg (norm2_nonneg A) (vnorm_nonneg xh); linarith
    rwa [etaFormula, div_le_iff₀ hD] at hη
  have hβ : 0 ≤ β := le_trans (rigal_gaches A b xh hb).1.1 (by rw [← eta_eq A b xh hb]; exact hη)
  have herr : xh - x = A⁻¹ *ᵥ (A *ᵥ xh - b) := by
    rw [Matrix.mulVec_sub, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hA, Matrix.one_mulVec,
      ← hx, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hA, Matrix.one_mulVec]
  have hbx : vnorm b ≤ norm2 A * vnorm x := by rw [← hx]; exact vnorm_mulVec_le A x
  have hxh : vnorm xh ≤ vnorm x + vnorm (xh - x) := by
    have := vnorm_add_le x (xh - x); rwa [add_sub_cancel] at this
  set e := vnorm (xh - x)
  have hκ0 : 0 ≤ cond2 A := cond2_nonneg A
  have hmain : e ≤ cond2 A * β * (2 * vnorm x + e) := by
    calc e = vnorm (A⁻¹ *ᵥ (A *ᵥ xh - b)) := by rw [show e = vnorm (xh - x) from rfl, herr]
      _ ≤ norm2 A⁻¹ * vnorm (A *ᵥ xh - b) := vnorm_mulVec_le _ _
      _ = norm2 A⁻¹ * vnorm (b - A *ᵥ xh) := by rw [vnorm_sub_comm]
      _ ≤ norm2 A⁻¹ * (β * (norm2 A * vnorm xh + vnorm b)) :=
          mul_le_mul_of_nonneg_left hres (norm2_nonneg _)
      _ ≤ norm2 A⁻¹ * (β * (norm2 A * (vnorm x + e) + norm2 A * vnorm x)) := by
          apply mul_le_mul_of_nonneg_left _ (norm2_nonneg _)
          apply mul_le_mul_of_nonneg_left _ hβ
          exact add_le_add (mul_le_mul_of_nonneg_left hxh (norm2_nonneg _)) hbx
      _ = cond2 A * β * (2 * vnorm x + e) := by unfold cond2; ring
  rw [div_le_div_iff₀ hxpos (by linarith)]
  nlinarith

/-! ## The condition number of an SPD matrix -/

section SPDcond
variable [NeZero n] (A : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosDef)

/-- The largest eigenvalue. -/
noncomputable def lamMax : ℝ := Finset.univ.sup' Finset.univ_nonempty hA.1.eigenvalues

/-- The smallest eigenvalue. -/
noncomputable def lamMin : ℝ := Finset.univ.inf' Finset.univ_nonempty hA.1.eigenvalues

include hA in
theorem cond2_spd :
    norm2 A = lamMax A hA ∧ norm2 A⁻¹ = 1 / lamMin A hA ∧
      cond2 A = lamMax A hA / lamMin A hA := by
  have hH := hA.1
  have hsymm : Aᵀ = A := ((posDef_iff_real A).mp hA).1
  have hpos := (posDef_iff_eigenvalues_pos hH).mp hA
  obtain ⟨imax, -, hmax⟩ := Finset.exists_mem_eq_sup' Finset.univ_nonempty hH.eigenvalues
  obtain ⟨imin, -, hmin⟩ := Finset.exists_mem_eq_inf' Finset.univ_nonempty hH.eigenvalues
  have hle : ∀ i, hH.eigenvalues i ≤ lamMax A hA := fun i =>
    Finset.le_sup' hH.eigenvalues (Finset.mem_univ i)
  have hge : ∀ i, lamMin A hA ≤ hH.eigenvalues i := fun i =>
    Finset.inf'_le hH.eigenvalues (Finset.mem_univ i)
  have hmaxeq : lamMax A hA = hH.eigenvalues imax := hmax
  have hmineq : lamMin A hA = hH.eigenvalues imin := hmin
  have hmaxpos : 0 < lamMax A hA := by rw [hmaxeq]; exact hpos imax
  have hminpos : 0 < lamMin A hA := by rw [hmineq]; exact hpos imin
  have hquad := quad_eq_sum_eigenvalues hH
  have hz := eigQ_dot hH
  have hupper : ∀ v, v ⬝ᵥ (A *ᵥ v) ≤ lamMax A hA * vnorm v ^ 2 := by
    intro v
    rw [hquad, vnorm_sq, ← hz, dotProduct, Finset.mul_sum]
    exact Finset.sum_le_sum fun i _ => by
      rw [← sq]; exact mul_le_mul_of_nonneg_right (hle i) (sq_nonneg _)
  have hlower : ∀ v, lamMin A hA * vnorm v ^ 2 ≤ v ⬝ᵥ (A *ᵥ v) := by
    intro v
    rw [hquad, vnorm_sq, ← hz, dotProduct, Finset.mul_sum]
    exact Finset.sum_le_sum fun i _ => by
      rw [← sq]; exact mul_le_mul_of_nonneg_right (hge i) (sq_nonneg _)
  have hpsd : ∀ v, 0 ≤ v ⬝ᵥ (A *ᵥ v) := fun v =>
    le_trans (mul_nonneg hminpos.le (sq_nonneg _)) (hlower v)
  have hdet : IsUnit A.det := (Matrix.isUnit_iff_isUnit_det A).mp (posDef_isUnit hA)
  have h1 : norm2 A = lamMax A hA := by
    apply le_antisymm (norm2_le_of_quad_psd hsymm hmaxpos hpsd hupper)
    set e : Fin n → ℝ := ⇑(hH.eigenvectorBasis imax)
    apply le_norm2_of (eigenvector_ne_zero hH imax)
    rw [show A *ᵥ e = hH.eigenvalues imax • e from hH.mulVec_eigenvectorBasis imax, vnorm_smul,
      abs_of_pos (hpos imax), hmaxeq]
  have h2 : norm2 A⁻¹ = 1 / lamMin A hA := by
    apply le_antisymm (norm2_inv_le_of_quad hminpos hlower)
    set e : Fin n → ℝ := ⇑(hH.eigenvectorBasis imin)
    have he : A *ᵥ e = lamMin A hA • e := by
      rw [hmineq]; exact hH.mulVec_eigenvectorBasis imin
    have hinv : A⁻¹ *ᵥ e = (1 / lamMin A hA) • e := by
      have : e = A *ᵥ ((1 / lamMin A hA) • e) := by
        rw [Matrix.mulVec_smul, he, smul_smul, one_div_mul_cancel (ne_of_gt hminpos), one_smul]
      conv_lhs => rw [this]
      rw [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hdet, Matrix.one_mulVec]
    apply le_norm2_of (eigenvector_ne_zero hH imin)
    rw [hinv, vnorm_smul, abs_of_pos (by positivity)]
  refine ⟨h1, h2, ?_⟩
  unfold cond2
  rw [h1, h2]
  ring

end SPDcond

end Cholesky
