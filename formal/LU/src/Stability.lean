import src.BackwardError
import src.ForwardError

/-!
# LU with partial pivoting is backward stable up to growth

* `eta_lu_le`: the pair `(E, 0)` of `solve_backward` gives
  `η(x̂) ≤ ‖E‖/‖A‖ ≤ γ_{3n} ‖ |L̂| |Û| ‖ / ‖A‖`, using `‖P E‖ = ‖E‖`.
* `gamma_asymp`: `γ_{3n} ≈ 3 n u`, i.e. `γ_{3n} / (3 n u) → 1` as `u → 0`.
-/
open Filter Topology
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {n : ℕ} {u : ℝ}

/-- **Backward error of the computed solution.** -/
theorem eta_lu_le (o : Ops ℝ) (ho : o.StdModel u) (hu : 0 ≤ u)
    (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) (hA : 0 < norm2 A)
    (hnu : ((3 * n : ℕ) : ℝ) * u < 1) (hd : ∀ i, factorU o A i i ≠ 0) :
    ∃ E : Matrix (Fin n) (Fin n) ℝ, (A + E) *ᵥ luSolve o A b = b ∧
      eta A b (luSolve o A b) ≤ norm2 E / norm2 A ∧
      norm2 E / norm2 A ≤
        gamma u (3 * n) * norm2 (absMat (factorL o A) * absMat (factorU o A)) / norm2 A := by
  obtain ⟨E, hE, hPE⟩ := solve_backward o ho hu A b hnu hd
  have hγ : 0 ≤ gamma u (3 * n) := gamma_nonneg hu hnu
  refine ⟨E, hE, ?_, ?_⟩
  · have h0 : (A + E) *ᵥ luSolve o A b = b + 0 := by rw [add_zero]; exact hE
    calc eta A b (luSolve o A b) ≤ tau A E b 0 := eta_le_tau h0
      _ = norm2 E / norm2 A := by
          rw [tau, vnorm_zero, zero_div]
          exact max_eq_left (div_nonneg (norm2_nonneg _) (norm2_nonneg _))
  · apply div_le_div_of_nonneg_right _ hA.le
    have hdom : ∀ i j, |(factorP o A * E) i j| ≤
        (gamma u (3 * n) • (absMat (factorL o A) * absMat (factorU o A))) i j := by
      intro i j
      rw [Matrix.smul_apply, smul_eq_mul]
      exact hPE i j
    calc norm2 E = norm2 (factorP o A * E) := (norm2_perm_mul _ E).symm
      _ ≤ norm2 (gamma u (3 * n) • (absMat (factorL o A) * absMat (factorU o A))) :=
          norm2_le_of_abs_le hdom
      _ = gamma u (3 * n) * norm2 (absMat (factorL o A) * absMat (factorU o A)) := by
          rw [norm2_smul, abs_of_nonneg hγ]

/-- `γ_{3n} ≈ 3 n u` for small `u`. -/
theorem gamma_asymp (n : ℕ) (hn : 0 < n) :
    Tendsto (fun u : ℝ => gamma u (3 * n) / (3 * n * u)) (𝓝[>] 0) (𝓝 1) := by
  have hc : Tendsto (fun u : ℝ => 1 / (1 - 3 * n * u)) (𝓝 0) (𝓝 (1 / (1 - 3 * n * 0))) := by
    apply Tendsto.div tendsto_const_nhds
    · exact tendsto_const_nhds.sub (tendsto_id.const_mul _)
    · simp
  simp only [mul_zero, sub_zero, div_one] at hc
  refine (hc.mono_left nhdsWithin_le_nhds).congr' ?_
  have hev : ∀ᶠ u in 𝓝[>] (0 : ℝ), 0 < u ∧ 3 * (n : ℝ) * u < 1 := by
    have h1 : ∀ᶠ u in 𝓝[>] (0 : ℝ), 0 < u := self_mem_nhdsWithin
    have h2 : ∀ᶠ u in 𝓝[>] (0 : ℝ), 3 * (n : ℝ) * u < 1 := by
      apply nhdsWithin_le_nhds
      have : Tendsto (fun u : ℝ => 3 * (n : ℝ) * u) (𝓝 0) (𝓝 0) := by
        simpa using ((tendsto_id : Tendsto (fun u : ℝ => u) (𝓝 0) (𝓝 0)).const_mul
          (3 * (n : ℝ)))
      exact this (Iio_mem_nhds (by norm_num))
    exact h1.and h2
  filter_upwards [hev] with u ⟨hu, hu3⟩
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have e1 : (1 : ℝ) - 3 * n * u ≠ 0 := by linarith
  unfold gamma
  push_cast
  field_simp

end LU
