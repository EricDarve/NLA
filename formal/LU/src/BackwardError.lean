import src.Factorization
import src.Solve

/-!
# Rounding errors: backward error of LU with partial pivoting

Theorem `thm:backward_error_lu` for real arithmetic in the standard model with unit
roundoff `u` (rounding to nearest satisfies it, `stdModel_of_roundNearest`):

* `factor_backward`: the computed factors satisfy `P A + F = L̂ Û` with
  `|F| ≤ γ_n |L̂| |Û|`. Here `P` is the permutation that the computation itself chose.
  The pivoting does not enter the error bound: the induction is the one of the existence
  proof, with one rounding error per multiplier and two per update.
* `solve_backward`: if `Û` has no zero diagonal entry, the computed solution of
  `L̂ ŷ = P b`, `Û x̂ = ŷ` satisfies `(A + E) x̂ = b` with `|P E| ≤ γ_{3n} |L̂| |Û|`.

The routine always completes in this model: a zero pivot is skipped.
-/
set_option linter.unusedSectionVars false

namespace LU
open Matrix

variable {u : ℝ} {n m : ℕ}

/-! ## Two scalar inequalities -/

/-- A rounded multiplier `ℓ = (v / α)(1 + δ)` reconstructs `v` with error `γ₁ |ℓ| |α|`. -/
lemma mult_err (hu1 : u < 1) {v α ℓ δ : ℝ} (hα : α ≠ 0) (hδ : |δ| ≤ u)
    (hℓ : ℓ = v / α * (1 + δ)) : |ℓ * α - v| ≤ gamma u 1 * (|ℓ| * |α|) := by
  have h1 : 1 + δ ≠ 0 := by have := abs_le.mp hδ; linarith
  have e1 : ℓ * α - v = v * (1 + δ) * (δ / (1 + δ)) := by rw [hℓ]; field_simp; ring
  have e2 : |ℓ| * |α| = |v * (1 + δ)| := by rw [← abs_mul, hℓ]; congr 1; field_simp
  rw [e1, e2, abs_mul (v * (1 + δ)), mul_comm (gamma u 1)]
  exact mul_le_mul_of_nonneg_left (abs_div_one_add_le hu1 hδ) (abs_nonneg _)

/-- One update `t = (B - ℓ w (1 + ε))(1 + δ)` followed by a trailing reconstruction
`t + e` with `|e| ≤ γ_m M` and `|t + e| ≤ M`. -/
lemma upd_err (hu : 0 ≤ u) {k : ℕ} (hk : ((k + 1 : ℕ) : ℝ) * u < 1)
    {B ℓ w t e M δ ε : ℝ} (hδ : |δ| ≤ u) (hε : |ε| ≤ u)
    (ht : t = (B - ℓ * w * (1 + ε)) * (1 + δ)) (hM : |t + e| ≤ M) (he : |e| ≤ gamma u k * M) :
    |ℓ * w + (t + e) - B| ≤ gamma u (k + 1) * (|ℓ| * |w| + M) := by
  have hkR : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  push_cast at hk
  have hu1 : u < 1 := by nlinarith
  have h1 : 1 + δ ≠ 0 := by have := abs_le.mp hδ; linarith
  have hgk : 0 ≤ gamma u k := gamma_nonneg hu (by nlinarith)
  have hsum := gamma_add_le hu (a := 1) (b := k) (by push_cast; linarith)
  rw [show 1 + k = k + 1 by ring] at hsum
  have hug : u ≤ gamma u (k + 1) := le_gamma hu (by omega) (by push_cast; linarith)
  have hM0 : 0 ≤ M := le_trans (abs_nonneg _) hM
  have hB : B = t / (1 + δ) + ℓ * w * (1 + ε) := by rw [ht]; field_simp; ring
  have key : ℓ * w + (t + e) - B = t * (1 - 1 / (1 + δ)) + e - ℓ * w * ε := by
    rw [hB]; field_simp; ring
  have ht' : |t| ≤ (1 + gamma u k) * M := by
    have : |t| ≤ |t + e| + |e| := by
      calc |t| = |(t + e) - e| := by ring_nf
        _ ≤ |t + e| + |e| := abs_sub _ _
    linarith
  have hd : |1 - 1 / (1 + δ)| ≤ gamma u 1 := by
    rw [abs_sub_comm]; exact abs_inv_one_add_sub_one_le hu1 hδ
  have hg1 : 0 ≤ gamma u 1 := le_trans (abs_nonneg _) hd
  rw [key]
  calc |t * (1 - 1 / (1 + δ)) + e - ℓ * w * ε|
      ≤ |t| * |1 - 1 / (1 + δ)| + |e| + |ℓ| * |w| * |ε| := by
        rw [← abs_mul, ← abs_mul, ← abs_mul]
        exact le_trans (abs_sub _ _) (add_le_add (abs_add_le _ _) le_rfl)
    _ ≤ (1 + gamma u k) * M * gamma u 1 + gamma u k * M + |ℓ| * |w| * u := by
        gcongr
    _ = (gamma u 1 + gamma u k + gamma u 1 * gamma u k) * M + u * (|ℓ| * |w|) := by ring
    _ ≤ gamma u (k + 1) * M + gamma u (k + 1) * (|ℓ| * |w|) := by
        gcongr
    _ = gamma u (k + 1) * (|ℓ| * |w| + M) := by ring

/-! ## The factorization -/

lemma absProd_succ_zero (M : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (i : Fin m) :
    (absMat (lowerOf M) * absMat (upperOf M)) i.succ 0 = |M i.succ 0| * |M 0 0| := by
  rw [absMat_lowerOf, absMat_upperOf, mul_succ_zero]; rfl

lemma absProd_succ_succ (M : Matrix (Fin (m + 1)) (Fin (m + 1)) ℝ) (i j : Fin m) :
    (absMat (lowerOf M) * absMat (upperOf M)) i.succ j.succ = |M i.succ 0| * |M 0 j.succ| +
      (absMat (lowerOf (M.submatrix Fin.succ Fin.succ)) *
        absMat (upperOf (M.submatrix Fin.succ Fin.succ))) i j := by
  rw [absMat_lowerOf, absMat_upperOf, mul_succ_succ, absMat_lowerOf, absMat_upperOf]; rfl

/-- **Backward error of the factorization**: `|L̂ Û - P A| ≤ γ_n |L̂| |Û|`. -/
theorem factor_backward (o : Ops ℝ) (ho : o.StdModel u) (hu : 0 ≤ u) :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ), (n : ℝ) * u < 1 → ∀ i j,
      |(factorL o A * factorU o A - factorP o A * A) i j| ≤
        gamma u n * (absMat (factorL o A) * absMat (factorU o A)) i j
  | 0, _, _, i, _ => Fin.elim0 i
  | m + 1, A, hnu, i, j => by
    have hmR : (0 : ℝ) ≤ m := Nat.cast_nonneg m
    have hnu' : ((m + 1 : ℕ) : ℝ) * u < 1 := by exact_mod_cast hnu
    push_cast at hnu
    have hu1 : u < 1 := by nlinarith
    have hmu : (m : ℝ) * u < 1 := by nlinarith
    have ih := factor_backward o ho hu (trail o A) hmu
    have hγ1 : gamma u 1 ≤ gamma u (m + 1) := gamma_mono hu (by omega) hnu'
    have hγm : gamma u m ≤ gamma u (m + 1) := gamma_mono hu (by omega) hnu'
    have hgm : 0 ≤ gamma u m := gamma_nonneg hu hmu
    have hgn : 0 ≤ gamma u (m + 1) := gamma_nonneg hu hnu'
    unfold factorL factorU factorP factorPerm at ih ⊢
    rw [Matrix.sub_apply]
    refine Fin.cases ?_ (fun i => ?_) i
    · rw [permMul_zero, mul_zero_apply, final_row, sub_self, abs_zero]
      exact mul_nonneg hgn (absProd_nonneg _ _ _ _)
    · refine Fin.cases ?_ (fun j => ?_) j
      · rw [permMul_succ, mul_succ_zero, absProd_succ_zero, final_col, final_row, first_col]
        split_ifs with h
        · rw [h, swapped_col_eq_zero A h, mul_zero, sub_zero, abs_zero]
          exact mul_nonneg hgn (by positivity)
        · obtain ⟨δ, hδ, hℓ⟩ := ho.div (swapped A ((final o (trail o A)).perm i).succ 0)
            (swapped A 0 0)
          exact le_trans (mult_err hu1 h hδ hℓ)
            (mul_le_mul_of_nonneg_right hγ1 (mul_nonneg (abs_nonneg _) (abs_nonneg _)))
      · rw [permMul_succ, mul_succ_succ, absProd_succ_succ, final_col, final_row, final_trail]
        have hij := ih i j
        rw [Matrix.sub_apply, permMatrix_mul_apply, trail_apply] at hij
        rw [first_col]
        set i' := (final o (trail o A)).perm i
        set N := (final o (trail o A)).mat
        have hM : |(lowerOf N * upperOf N) i j| ≤ (absMat (lowerOf N) * absMat (upperOf N)) i j :=
          abs_mul_apply_le _ _ _ _
        split_ifs at hij ⊢ with h
        · rw [swapped_col_eq_zero A h, abs_zero, zero_mul, zero_mul, zero_add, zero_add]
          calc |(lowerOf N * upperOf N) i j - swapped A i'.succ j.succ|
              ≤ gamma u m * (absMat (lowerOf N) * absMat (upperOf N)) i j := hij
            _ ≤ gamma u (m + 1) * (absMat (lowerOf N) * absMat (upperOf N)) i j :=
                mul_le_mul_of_nonneg_right hγm (absProd_nonneg _ _ _ _)
        · obtain ⟨δ, ε, hδ, hε, ht⟩ := ho.upd (swapped A i'.succ j.succ)
            (o.div (swapped A i'.succ 0) (swapped A 0 0)) (swapped A 0 j.succ)
          have := upd_err hu hnu' hδ hε ht (e := (lowerOf N * upperOf N) i j -
            o.upd (swapped A i'.succ j.succ) (o.div (swapped A i'.succ 0) (swapped A 0 0))
              (swapped A 0 j.succ)) (by rw [add_sub_cancel]; exact hM) hij
          rw [add_sub_cancel] at this
          exact this

/-- `P A + F = L̂ Û` with `|F| ≤ γ_n |L̂| |Û|`. -/
theorem factor_backward' (o : Ops ℝ) (ho : o.StdModel u) (hu : 0 ≤ u)
    (A : Matrix (Fin n) (Fin n) ℝ) (hnu : (n : ℝ) * u < 1) :
    ∃ F : Matrix (Fin n) (Fin n) ℝ, factorP o A * A + F = factorL o A * factorU o A ∧
      ∀ i j, |F i j| ≤ gamma u n * (absMat (factorL o A) * absMat (factorU o A)) i j :=
  ⟨factorL o A * factorU o A - factorP o A * A, by abel, factor_backward o ho hu A hnu⟩

/-! ## The complete solve -/

/-- The computed solution of `A x = b`: `L̂ ŷ = P b` by forward substitution, then
`Û x̂ = ŷ` by back substitution. Forming `P b` only moves entries. -/
noncomputable def luSolve (o : Ops ℝ) (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) :
    Fin n → ℝ :=
  bwdSolve o (factorU o A) (fwdSolve o (factorL o A) (factorP o A *ᵥ b))

/-- In exact arithmetic the computed solution solves `A x = b` (when `U` has no zero
diagonal entry). -/
theorem exact_luSolve (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ)
    (hd : ∀ i, factorU (Ops.exact ℝ) A i i ≠ 0) : A *ᵥ luSolve (Ops.exact ℝ) A b = b := by
  have hL := factorL_isUnitLower (Ops.exact ℝ) A
  exact solve_of_plu (exact_plu A)
    (exact_fwdSolve _ _ hL.isLower fun i => by rw [hL.1 i]; exact one_ne_zero)
    (exact_bwdSolve _ _ (factorU_isUpper _ A) hd)

lemma permMatrix_mul_inv (σ : Equiv.Perm (Fin n)) :
    σ.permMatrix ℝ * (σ⁻¹).permMatrix ℝ = 1 := by
  have := permMatrix_inv_mul (𝕜 := ℝ) σ⁻¹
  rwa [inv_inv] at this

/-- **Backward error of the complete solve**: `(A + E) x̂ = b` with
`|P E| ≤ γ_{3n} |L̂| |Û|`. -/
theorem solve_backward (o : Ops ℝ) (ho : o.StdModel u) (hu : 0 ≤ u)
    (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) (hnu : ((3 * n : ℕ) : ℝ) * u < 1)
    (hd : ∀ i, factorU o A i i ≠ 0) :
    ∃ E : Matrix (Fin n) (Fin n) ℝ, (A + E) *ᵥ luSolve o A b = b ∧
      ∀ i j, |(factorP o A * E) i j| ≤
        gamma u (3 * n) * (absMat (factorL o A) * absMat (factorU o A)) i j := by
  have hnR : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hn0 : (n : ℝ) * u < 1 := by push_cast at hnu; nlinarith
  set L := factorL o A
  set U := factorU o A
  set σ := factorPerm o A
  have hP : factorP o A = σ.permMatrix ℝ := rfl
  have hL : IsUnitLower L := factorL_isUnitLower o A
  have hU : IsUpper U := factorU_isUpper o A
  have hLd : ∀ i, L i i ≠ 0 := fun i => by rw [hL.1 i]; exact one_ne_zero
  set y := fwdSolve o L (factorP o A *ᵥ b)
  obtain ⟨Δ1, h1, h1b⟩ := fwdSolve_backward o ho hu L (factorP o A *ᵥ b) hL.isLower hLd hn0
  obtain ⟨Δ2, h2, h2b⟩ := bwdSolve_backward o ho hu U y hU hd hn0
  have hF := factor_backward o ho hu A hn0
  set x := bwdSolve o U y
  set G := (L + Δ1) * (U + Δ2) - factorP o A * A
  refine ⟨(σ⁻¹).permMatrix ℝ * G, ?_, ?_⟩
  · have hx : luSolve o A b = x := rfl
    have hG : factorP o A * A + G = (L + Δ1) * (U + Δ2) := add_sub_cancel _ _
    have : A + (σ⁻¹).permMatrix ℝ * G = (σ⁻¹).permMatrix ℝ * ((L + Δ1) * (U + Δ2)) := by
      rw [← hG, Matrix.mul_add, hP, ← Matrix.mul_assoc, permMatrix_inv_mul, Matrix.one_mul]
    rw [hx, this, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, h2, h1, Matrix.mulVec_mulVec,
      hP, permMatrix_inv_mul, Matrix.one_mulVec]
  · intro i j
    rw [← Matrix.mul_assoc, hP, permMatrix_mul_inv, Matrix.one_mul]
    have hγ : 0 ≤ gamma u n := gamma_nonneg hu hn0
    have hexp : G = (L * U - factorP o A * A) + L * Δ2 + Δ1 * U + Δ1 * Δ2 := by
      simp only [G, Matrix.add_mul, Matrix.mul_add]; abel
    rw [hexp]
    have hLabs : ∀ i j, |L i j| ≤ 1 * |L i j| := fun i j => by rw [one_mul]
    have hUabs : ∀ i j, |U i j| ≤ 1 * |U i j| := fun i j => by rw [one_mul]
    have e1 := hF i j
    have e2 := abs_mul_le_of zero_le_one hLabs h2b i j
    have e3 := abs_mul_le_of hγ h1b hUabs i j
    have e4 := abs_mul_le_of hγ h1b h2b i j
    simp only [one_mul, mul_one] at e2 e3 e4
    have hc := gamma_three_le hu n hnu
    have hpos := absProd_nonneg L U i j
    calc |((L * U - factorP o A * A) + L * Δ2 + Δ1 * U + Δ1 * Δ2) i j|
        ≤ |(L * U - factorP o A * A) i j| + |(L * Δ2) i j| + |(Δ1 * U) i j| +
            |(Δ1 * Δ2) i j| := by
          simp only [Matrix.add_apply]
          exact le_trans (abs_add_le _ _) (add_le_add (le_trans (abs_add_le _ _)
            (add_le_add (abs_add_le _ _) le_rfl)) le_rfl)
      _ ≤ gamma u n * (absMat L * absMat U) i j + gamma u n * (absMat L * absMat U) i j +
            gamma u n * (absMat L * absMat U) i j +
            gamma u n * gamma u n * (absMat L * absMat U) i j := by
          gcongr
      _ = (gamma u n + 2 * gamma u n + gamma u n ^ 2) * (absMat L * absMat U) i j := by ring
      _ ≤ gamma u (3 * n) * (absMat L * absMat U) i j :=
          mul_le_mul_of_nonneg_right hc hpos

end LU
