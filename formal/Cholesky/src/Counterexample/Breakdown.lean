import src.Counterexample.Trace
import src.Stability
import src.Bounds

/-!
# Why a nonpositive pivot must follow

The proof "Breakdown from the sign pattern" of the notes, for any matrix `S` with
nonpositive off-diagonal entries and a direction `v` with `vᵀ S v ≤ -‖S‖₂ ‖v‖² / 18`:

* if the rounded algorithm completed on `S`, every off-diagonal entry of `L̂` would be
  nonpositive (`factor_offdiag_nonpos`), so `|L̂| = 2D - L̂` and `‖|L̂|‖₂ ≤ 3 ‖L̂‖₂`;
* the backward error bound (which does not need `S` to be SPD) gives
  `‖E‖₂ ≤ γ ‖|L̂|‖₂² ≤ 9 γ ‖S + E‖₂`, hence `‖E‖₂ < ‖S‖₂ / 100`;
* but `L̂ L̂ᵀ = S + E` is positive semidefinite, which forces `‖E‖₂ ≥ ‖S‖₂ / 18`.

Applied to `Ŝ` this proves that Cholesky fails on `A_{9k}` (`A_fails`), for both the
separately rounded and the fused update.
-/
namespace Cholesky
open Matrix

variable {u : ℝ}

/-! ## The sign pattern -/

/-- If the input has nonpositive off-diagonal entries, so has the computed factor. -/
theorem factor_offdiag_nonpos (o : Ops) (ho : o.StdModel u) (hu1 : u < 1) :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ), Succeeds o A → (∀ i j, j < i → A i j ≤ 0) →
      ∀ i j, j < i → factor o A i j ≤ 0
  | 0, _, _, _, i, _, _ => Fin.elim0 i
  | n + 1, A, h, hA, i, j, hij => by
    rw [succeeds_succ_iff] at h
    obtain ⟨h0, htrail⟩ := h
    have hd : 0 < o.sqrt (A 0 0) := by
      obtain ⟨δ, hδ, e⟩ := ho.sqrt (A 0 0)
      rw [e]; have := (abs_le.mp hδ).1
      exact mul_pos (Real.sqrt_pos.mpr h0) (by linarith)
    have hcol : ∀ i : Fin n, firstCol o A i ≤ 0 := by
      intro i
      obtain ⟨δ, hδ, e⟩ := ho.div (A i.succ 0) (o.sqrt (A 0 0))
      rw [firstCol, e]
      have := (abs_le.mp hδ).1
      exact mul_nonpos_of_nonpos_of_nonneg (div_nonpos_of_nonpos_of_nonneg (hA _ _ (Fin.succ_pos i))
        hd.le) (by linarith)
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨i, rfl⟩
    · exact absurd hij (Fin.not_lt_zero j)
    rcases Fin.eq_zero_or_eq_succ j with rfl | ⟨j, rfl⟩
    · rw [factor_succ_zero]; exact hcol i
    · rw [factor_succ_succ]
      apply factor_offdiag_nonpos o ho hu1 (trail o A) htrail _ i j (Fin.succ_lt_succ_iff.mp hij)
      intro a b hab
      rw [trail_apply _ _ _ _ hab.le]
      obtain ⟨δ, ε, hδ, hε, e⟩ := ho.upd (A a.succ b.succ) (firstCol o A a) (firstCol o A b)
      rw [e]
      have h1 := (abs_le.mp hδ).1
      have h2 := (abs_le.mp hε).1
      have hprod : 0 ≤ firstCol o A a * firstCol o A b * (1 + ε) :=
        mul_nonneg (mul_nonneg_of_nonpos_of_nonpos (hcol a) (hcol b)) (by linarith)
      have := hA a.succ b.succ (Fin.succ_lt_succ_iff.mpr hab)
      exact mul_nonpos_of_nonpos_of_nonneg (by linarith) (by linarith)

/-! ## Norm estimates -/

variable {n : ℕ}

lemma norm2_diagonal_le (d : Fin n → ℝ) {c : ℝ} (hc : 0 ≤ c) (hd : ∀ i, |d i| ≤ c) :
    norm2 (Matrix.diagonal d) ≤ c := by
  apply norm2_le_of hc
  intro v
  have hsq : (Matrix.diagonal d *ᵥ v) ⬝ᵥ (Matrix.diagonal d *ᵥ v) ≤ (c • v) ⬝ᵥ (c • v) := by
    simp only [Matrix.mulVec_diagonal, dotProduct, Pi.smul_apply, smul_eq_mul]
    apply Finset.sum_le_sum
    intro i _
    have := mul_self_le_mul_self (abs_nonneg _) (hd i)
    rw [abs_mul_abs_self] at this
    nlinarith [mul_self_nonneg (v i)]
  calc vnorm (Matrix.diagonal d *ᵥ v) ≤ vnorm (c • v) := vnorm_le_of_sq hsq
    _ = c * vnorm v := by rw [vnorm_smul, abs_of_nonneg hc]

/-- `|L| = 2D - L` when the diagonal is positive and the off-diagonal part nonpositive. -/
lemma absMat_eq_of_signs {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L)
    (hd : ∀ i, 0 < L i i) (hoff : ∀ i j, j < i → L i j ≤ 0) :
    absMat L = (2 : ℝ) • Matrix.diagonal (fun i => L i i) - L := by
  ext i j
  simp only [absMat_apply, Matrix.sub_apply, Matrix.smul_apply, Matrix.diagonal_apply,
    smul_eq_mul]
  rcases lt_trichotomy i j with h | rfl | h
  · rw [hL i j h, if_neg (ne_of_lt h)]; simp
  · rw [if_pos rfl, abs_of_pos (hd i)]; ring
  · rw [if_neg (ne_of_gt h), abs_of_nonpos (hoff i j h)]; ring

/-- `‖|L̂|‖₂ ≤ 3 ‖L̂‖₂`. -/
lemma norm2_absMat_le_three {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L)
    (hd : ∀ i, 0 < L i i) (hoff : ∀ i j, j < i → L i j ≤ 0) :
    norm2 (absMat L) ≤ 3 * norm2 L := by
  rw [absMat_eq_of_signs hL hd hoff]
  have hD : norm2 (Matrix.diagonal (fun i => L i i)) ≤ norm2 L :=
    norm2_diagonal_le _ (norm2_nonneg L) (fun i => abs_entry_le_norm2 L i i)
  calc norm2 ((2 : ℝ) • Matrix.diagonal (fun i => L i i) - L)
      ≤ norm2 ((2 : ℝ) • Matrix.diagonal (fun i => L i i)) + norm2 L := norm2_sub_le _ _
    _ = 2 * norm2 (Matrix.diagonal (fun i => L i i)) + norm2 L := by
        rw [norm2_smul, abs_of_pos (by norm_num)]
    _ ≤ 3 * norm2 L := by linarith

/-! ## The contradiction -/

/-- **Breakdown from the sign pattern.** -/
theorem no_success (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u) (S : Matrix (Fin n) (Fin n) ℝ)
    (hS : Sᵀ = S) (hoff : ∀ i j, j < i → S i j ≤ 0)
    (hnu : ((n + 1 : ℕ) : ℝ) * u ≤ 129 / 131072)
    (v : Fin n → ℝ) (hv : v ≠ 0) (hS0 : 0 < norm2 S)
    (hneg : v ⬝ᵥ (S *ᵥ v) ≤ -(norm2 S / 18) * vnorm v ^ 2) : ¬ Succeeds o S := by
  intro h
  have hnu1 : ((n + 1 : ℕ) : ℝ) * u < 1 := by linarith
  obtain ⟨hAE, -, hL, hd, hE⟩ := backward_error o ho hu S hnu1 h
  set L := factor o S
  rw [symmOfLower_eq_self hS] at hAE hE
  set E := L * Lᵀ - S with hEdef
  set s := ((n + 1 : ℕ) : ℝ) * u
  have hs0 : 0 ≤ s := by positivity
  have hu1 : u < 1 := by
    have : (1 : ℝ) * u ≤ s := mul_le_mul_of_nonneg_right (by norm_cast; omega) hu
    linarith
  have hγ : gamma u (n + 1) = s / (1 - s) := rfl
  have hγ0 : 0 ≤ gamma u (n + 1) := gamma_nonneg hu hnu1
  -- sign pattern
  have hoffL : ∀ i j, j < i → L i j ≤ 0 :=
    factor_offdiag_nonpos o ho hu1 S h (fun i j hij => by
      have := hoff i j hij
      simpa [symmOfLower, hij.le] using this)
  -- `‖E‖ ≤ 9 γ (‖S‖ + ‖E‖)`
  have hEle : norm2 E ≤ 9 * (s / (1 - s)) * (norm2 S + norm2 E) := by
    have h1 : norm2 E ≤ gamma u (n + 1) * norm2 (absMat L * (absMat L)ᵀ) := by
      rw [← abs_of_nonneg hγ0, ← norm2_smul]
      apply norm2_le_of_abs_le
      intro i j
      rw [Matrix.smul_apply, smul_eq_mul]
      exact hE i j
    have h2 : norm2 (absMat L * (absMat L)ᵀ) ≤ 9 * norm2 (L * Lᵀ) := by
      rw [norm2_mul_transpose_self, norm2_mul_transpose_self]
      have := norm2_absMat_le_three hL hd hoffL
      have := norm2_nonneg (absMat L)
      nlinarith
    have h3 : norm2 (L * Lᵀ) ≤ norm2 S + norm2 E := by
      rw [show L * Lᵀ = S + E by rw [hEdef]; abel]; exact norm2_add_le _ _
    calc norm2 E ≤ gamma u (n + 1) * norm2 (absMat L * (absMat L)ᵀ) := h1
      _ ≤ gamma u (n + 1) * (9 * (norm2 S + norm2 E)) :=
          mul_le_mul_of_nonneg_left (le_trans h2 (by linarith)) hγ0
      _ = 9 * (s / (1 - s)) * (norm2 S + norm2 E) := by rw [hγ]; ring
  -- positivity of `L Lᵀ = S + E` forces `‖E‖ ≥ ‖S‖ / 18`
  have hpsd : 0 ≤ v ⬝ᵥ ((L * Lᵀ) *ᵥ v) := by
    rw [quad_mul_transpose, ← vnorm_sq]; positivity
  have hsplit : v ⬝ᵥ ((L * Lᵀ) *ᵥ v) = v ⬝ᵥ (S *ᵥ v) + v ⬝ᵥ (E *ᵥ v) := by
    rw [show L * Lᵀ = S + E by rw [hEdef]; abel, Matrix.add_mulVec, dotProduct_add]
  have hEq := quad_le_norm2 E v
  have hvpos := vnorm_pos hv
  have hnec : norm2 S / 18 ≤ norm2 E := by
    have : norm2 S / 18 * vnorm v ^ 2 ≤ norm2 E * vnorm v ^ 2 := by nlinarith
    exact le_of_mul_le_mul_right this (by positivity)
  exact breakdown_scalar_contradiction (norm2 S) (norm2 E) s hS0 (norm2_nonneg _) hs0 hnu hEle hnec

/-! ## Cholesky fails on `A_{9k}` -/

namespace CE
namespace Par
variable (P : Par)

lemma Shat_offdiag : ∀ i j : Fin (P.k + P.k), j < i → P.Shat i j ≤ 0 := by
  have hη := P.η_pos
  intro i j hij
  revert hij
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j <;> intro hij
  · simp only [Shat, blk_ll, Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply, ones,
      Matrix.of_apply, smul_eq_mul]
    rw [if_neg (ne_of_gt (cA_lt.mp hij))]; linarith
  · exact absurd (le_of_lt hij) not_nA_le_cA
  · simp only [Shat, blk_rl, Matrix.smul_apply, ones, Matrix.of_apply, smul_eq_mul]; linarith
  · simp only [Shat, blk_rr, Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply, ones,
      Matrix.of_apply, smul_eq_mul]
    rw [if_neg (ne_of_gt (nA_lt.mp hij))]; linarith

lemma ones_ne_zero' : (fun _ : Fin (P.k + P.k) => (1 : ℝ)) ≠ 0 := by
  intro h
  have := congrFun h (Fin.castAdd P.k ⟨0, by have := P.k_ge; omega⟩)
  simp at this

/-- The negative direction: `1ᵀ Ŝ 1 ≤ -‖Ŝ‖₂ ‖1‖² / 18`. -/
lemma Shat_neg_direction :
    (fun _ => (1 : ℝ)) ⬝ᵥ (P.Shat *ᵥ fun _ => 1) ≤
      -(norm2 P.Shat / 18) * vnorm (fun _ : Fin (P.k + P.k) => (1 : ℝ)) ^ 2 := by
  have hk : (64 : ℝ) ≤ P.k := by exact_mod_cast P.k_ge
  have hη := P.η_pos
  have hnorm := P.norm2_Shat_le
  rw [P.Shat_ones, dotProduct_smul, smul_eq_mul, ← vnorm_sq]
  have hv : 0 ≤ vnorm (fun _ : Fin (P.k + P.k) => (1 : ℝ)) ^ 2 := sq_nonneg _
  have hmargin := negative_eigenvalue_margin P.k hk
  rw [lt_div_iff₀ (by linarith)] at hmargin
  nlinarith [mul_nonneg hv hη.le]

lemma norm2_Shat_pos : 0 < norm2 P.Shat := by
  have hk : (64 : ℝ) ≤ P.k := by exact_mod_cast P.k_ge
  have hη := P.η_pos
  have hv := vnorm_pos P.ones_ne_zero'
  have h := abs_quad_le_norm2 P.Shat (fun _ => 1)
  rw [P.Shat_ones, dotProduct_smul, smul_eq_mul, ← vnorm_sq, abs_mul,
    abs_of_neg (by nlinarith)] at h
  by_contra hle
  push_neg at hle
  have := norm2_nonneg P.Shat
  have h0 : norm2 P.Shat = 0 := le_antisymm hle this
  rw [h0, zero_mul, abs_of_nonneg (sq_nonneg _)] at h
  have hpos : 0 < -(-((P.k : ℝ) - 9) * P.η / 8) := by nlinarith
  nlinarith [mul_pos hpos (pow_pos hv 2)]

lemma budget : ((P.k + P.k + 1 : ℕ) : ℝ) * P.u ≤ 129 / 131072 := by
  have h := update_budget P.k P.u P.ku_le P.u_le
  push_cast at h ⊢; linarith

/-- `A_{9k}` is the matrix `leadMatrix T C` of the leading-block lemma. -/
lemma A_eq : P.A = leadMatrix P.T P.C := rfl

/-- **Cholesky fails on `A_{9k}`**, for any rounding to nearest and either update variant. -/
theorem A_fails_of {o : Ops} {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) (ho : IsRN o r P.p)
    (hstd : o.StdModel P.u) : ¬ Succeeds o P.A := by
  intro h
  rw [A_eq] at h
  have h1 := lead_succeeds_Citer (P.leadOps hr ho) h
  have h2 : Succeeds o P.Shat := (succeeds_lowerEq o (P.trailing_eq_Shat hr ho)).mp h1
  exact no_success o hstd P.u_pos.le P.Shat P.Shat_symm P.Shat_offdiag P.budget _
    P.ones_ne_zero' P.norm2_Shat_pos P.Shat_neg_direction h2

theorem A_fails {r : ℝ → ℝ} (hr : IsRoundNearest P.p r) :
    ¬ Succeeds (Ops.rounded r) P.A ∧ ¬ Succeeds (Ops.fused r) P.A := by
  obtain ⟨h1, h2⟩ := stdModel_of_roundNearest P.p_pos hr
  exact ⟨P.A_fails_of hr (isRN_rounded hr) h1, P.A_fails_of hr isRN_fused h2⟩

end Par
end CE
end Cholesky
