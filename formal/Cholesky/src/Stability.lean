import src.BackwardError
import src.Norms

/-!
# Normwise backward stability of Cholesky

From the componentwise theorem `|E| ≤ γ |L̂| |L̂|ᵀ` (with `γ = γ_{n+1} < 1`):

* the diagonal equations `∑ₖ l̂ᵢₖ² = aᵢᵢ + eᵢᵢ ≤ aᵢᵢ + γ ∑ₖ l̂ᵢₖ²`;
* `‖L̂‖_F² ≤ tr(A) / (1 - γ)`;
* `‖E‖₂ ≤ γ ‖|L̂| |L̂|ᵀ‖₂ ≤ γ ‖L̂‖_F² ≤ n γ / (1 - γ) ‖A‖₂`;
* `A + E` is SPD and `L̂` is its (unique) Cholesky factor;
* `n γ_{n+1} / (1 - γ_{n+1}) = n (n + 1) u / (1 - 2 (n + 1) u) ~ n (n + 1) u`.

In exact arithmetic the same argument gives `‖|L| |L|ᵀ‖₂ ≤ ‖L‖_F² = tr(A) ≤ n ‖A‖₂`.
-/
open Filter Topology

namespace Cholesky
open Matrix

variable {n : ℕ} {u : ℝ}

/-- `‖|L| |L|ᵀ‖₂ ≤ ‖|L|‖₂² ≤ ‖L‖_F²`. -/
theorem norm2_absProd_le_frob_sq (L : Matrix (Fin n) (Fin n) ℝ) :
    norm2 (absMat L * (absMat L)ᵀ) ≤ frob L ^ 2 := by
  calc norm2 (absMat L * (absMat L)ᵀ) = norm2 (absMat L) ^ 2 := norm2_mul_transpose_self _
    _ ≤ frob (absMat L) ^ 2 :=
        pow_le_pow_left₀ (norm2_nonneg _) (norm2_le_frob _) 2
    _ = frob L ^ 2 := by rw [frob_absMat]

lemma frob_sq_eq_trace (L : Matrix (Fin n) (Fin n) ℝ) : frob L ^ 2 = (L * Lᵀ).trace := by
  rw [frob_sq, Matrix.trace]
  apply Finset.sum_congr rfl
  intro i _
  simp [Matrix.mul_apply, sq]

/-- In exact arithmetic `‖|L| |L|ᵀ‖₂ ≤ ‖L‖_F² = tr(A) ≤ n ‖A‖₂`. -/
theorem exact_absProd_bound (L : Matrix (Fin n) (Fin n) ℝ) :
    norm2 (absMat L * (absMat L)ᵀ) ≤ (L * Lᵀ).trace ∧ (L * Lᵀ).trace ≤ n * norm2 (L * Lᵀ) := by
  refine ⟨by rw [← frob_sq_eq_trace]; exact norm2_absProd_le_frob_sq L, ?_⟩
  simpa using trace_le (L * Lᵀ)

section Computed
variable (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u) (A : Matrix (Fin n) (Fin n) ℝ)
  (hnu : ((n + 1 : ℕ) : ℝ) * u < 1) (h : Succeeds o A)
include ho hu hnu h

/-- The diagonal equations: `∑ₖ l̂ᵢₖ² = aᵢᵢ + eᵢᵢ ≤ aᵢᵢ + γ ∑ₖ l̂ᵢₖ²`. -/
theorem diag_equation (i : Fin n) :
    ∑ k, factor o A i k ^ 2 ≤ symmOfLower A i i + gamma u (n + 1) * ∑ k, factor o A i k ^ 2 := by
  obtain ⟨hAE, -, -, -, hE⟩ := backward_error o ho hu A hnu h
  have h1 := hE i i
  have hdiag : (factor o A * (factor o A)ᵀ) i i = ∑ k, factor o A i k ^ 2 := by
    simp [Matrix.mul_apply, sq]
  have hP : (absMat (factor o A) * (absMat (factor o A))ᵀ) i i = ∑ k, factor o A i k ^ 2 := by
    simp [Matrix.mul_apply, ← sq, sq_abs]
  have hEii : (factor o A * (factor o A)ᵀ - symmOfLower A) i i =
      ∑ k, factor o A i k ^ 2 - symmOfLower A i i := by
    rw [Matrix.sub_apply, hdiag]
  rw [hEii, hP] at h1
  linarith [le_abs_self (∑ k, factor o A i k ^ 2 - symmOfLower A i i)]

/-- `‖L̂‖_F² ≤ tr(A) / (1 - γ)`. -/
theorem frob_sq_le (hγ : gamma u (n + 1) < 1) :
    frob (factor o A) ^ 2 ≤ (symmOfLower A).trace / (1 - gamma u (n + 1)) := by
  rw [le_div_iff₀ (by linarith), frob_sq, Matrix.trace, Finset.sum_mul]
  apply Finset.sum_le_sum
  intro i _
  have := diag_equation o ho hu A hnu h i
  simp only [Matrix.diag]
  linarith

/-- `‖E‖₂ ≤ γ ‖|L̂| |L̂|ᵀ‖₂ ≤ γ ‖L̂‖_F²`. -/
theorem norm2_E_le :
    norm2 (factor o A * (factor o A)ᵀ - symmOfLower A) ≤
        gamma u (n + 1) * norm2 (absMat (factor o A) * (absMat (factor o A))ᵀ) ∧
      gamma u (n + 1) * norm2 (absMat (factor o A) * (absMat (factor o A))ᵀ) ≤
        gamma u (n + 1) * frob (factor o A) ^ 2 := by
  obtain ⟨-, -, -, -, hE⟩ := backward_error o ho hu A hnu h
  have hγ : 0 ≤ gamma u (n + 1) := gamma_nonneg hu hnu
  refine ⟨?_, mul_le_mul_of_nonneg_left (norm2_absProd_le_frob_sq _) hγ⟩
  rw [← abs_of_nonneg hγ, ← norm2_smul]
  apply norm2_le_of_abs_le
  intro i j
  rw [Matrix.smul_apply, smul_eq_mul]
  exact hE i j

/-- **Backward stability.** `‖E‖₂ ≤ n γ / (1 - γ) ‖A‖₂` with `γ = γ_{n+1} < 1`. -/
theorem norm2_E_le_rel (hγ : gamma u (n + 1) < 1) :
    norm2 (factor o A * (factor o A)ᵀ - symmOfLower A) ≤
      n * gamma u (n + 1) / (1 - gamma u (n + 1)) * norm2 (symmOfLower A) := by
  have hγ0 : 0 ≤ gamma u (n + 1) := gamma_nonneg hu hnu
  have h1 := (norm2_E_le o ho hu A hnu h)
  have h2 := frob_sq_le o ho hu A hnu h hγ
  have h3 : (symmOfLower A).trace ≤ n * norm2 (symmOfLower A) := by
    simpa using trace_le (symmOfLower A)
  have h4 : (symmOfLower A).trace / (1 - gamma u (n + 1)) ≤
      n * norm2 (symmOfLower A) / (1 - gamma u (n + 1)) :=
    div_le_div_of_nonneg_right h3 (by linarith)
  calc norm2 (factor o A * (factor o A)ᵀ - symmOfLower A)
      ≤ gamma u (n + 1) * frob (factor o A) ^ 2 := le_trans h1.1 h1.2
    _ ≤ gamma u (n + 1) * (n * norm2 (symmOfLower A) / (1 - gamma u (n + 1))) :=
        mul_le_mul_of_nonneg_left (le_trans h2 h4) hγ0
    _ = n * gamma u (n + 1) / (1 - gamma u (n + 1)) * norm2 (symmOfLower A) := by ring

/-- A completed factorization is the exact Cholesky factor of the nearby SPD matrix `A + E`. -/
theorem nearby_spd :
    (factor o A * (factor o A)ᵀ).PosDef ∧
      ∀ M : Matrix (Fin n) (Fin n) ℝ, IsLowerTriangular M → (∀ i, 0 < M i i) →
        M * Mᵀ = factor o A * (factor o A)ᵀ → M = factor o A := by
  obtain ⟨-, -, hL, hd, -⟩ := backward_error o ho hu A hnu h
  exact ⟨posDef_of_factorization _ hL hd,
    fun M hM hMd hMM => cholesky_unique M (factor o A) hM hL hMd hd hMM⟩

end Computed

/-! ## The size of the bound -/

lemma gamma_ratio (m : ℕ) (hu : 2 * m * u < 1) :
    gamma u m / (1 - gamma u m) = m * u / (1 - 2 * m * u) := by
  unfold gamma
  have h1 : (1 : ℝ) - m * u ≠ 0 := by nlinarith
  have h2 : (1 : ℝ) - 2 * m * u ≠ 0 := by linarith
  have h3 : 1 - m * u / (1 - m * u) = (1 - 2 * m * u) / (1 - m * u) := by
    field_simp; ring
  rw [h3, div_div_div_cancel_right₀ h1]

/-- `n γ_{n+1} / (1 - γ_{n+1}) = n (n + 1) u / (1 - 2 (n + 1) u)`. -/
theorem stability_constant_eq (h : 2 * ((n + 1 : ℕ) : ℝ) * u < 1) :
    n * gamma u (n + 1) / (1 - gamma u (n + 1)) =
      n * (n + 1) * u / (1 - 2 * (n + 1) * u) := by
  rw [mul_div_assoc, gamma_ratio (n + 1) h]
  push_cast; ring

/-- For small `u` the constant is `n (n + 1) u (1 + o(1))`. -/
theorem stability_constant_asymp (n : ℕ) (hn : 0 < n) :
    Tendsto (fun u : ℝ => (n * gamma u (n + 1) / (1 - gamma u (n + 1))) / (n * (n + 1) * u))
      (𝓝[>] 0) (𝓝 1) := by
  have hlim : Tendsto (fun u : ℝ => 1 / (1 - 2 * (n + 1) * u)) (𝓝 0) (𝓝 1) := by
    have : Tendsto (fun u : ℝ => 1 / (1 - 2 * (n + 1) * u)) (𝓝 0)
        (𝓝 (1 / (1 - 2 * (n + 1) * 0))) := by
      apply Tendsto.div tendsto_const_nhds
      · exact tendsto_const_nhds.sub (tendsto_id.const_mul _)
      · simp
    simpa using this
  refine (hlim.mono_left nhdsWithin_le_nhds).congr' ?_
  have hev : ∀ᶠ u in 𝓝[>] (0 : ℝ), 0 < u ∧ 2 * ((n + 1 : ℕ) : ℝ) * u < 1 := by
    have h1 : ∀ᶠ u in 𝓝[>] (0 : ℝ), 0 < u := self_mem_nhdsWithin
    have h2 : ∀ᶠ u in 𝓝[>] (0 : ℝ), 2 * ((n + 1 : ℕ) : ℝ) * u < 1 := by
      apply nhdsWithin_le_nhds
      have : Tendsto (fun u : ℝ => 2 * ((n + 1 : ℕ) : ℝ) * u) (𝓝 0) (𝓝 0) := by
        simpa using ((tendsto_id : Tendsto (fun u : ℝ => u) (𝓝 0) (𝓝 0)).const_mul
          (2 * ((n + 1 : ℕ) : ℝ)))
      exact this (Iio_mem_nhds (by norm_num))
    exact h1.and h2
  filter_upwards [hev] with u ⟨hu, hu2⟩
  rw [stability_constant_eq hu2]
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have : (1 : ℝ) - 2 * (n + 1) * u ≠ 0 := by push_cast at hu2; linarith
  field_simp

end Cholesky
