import src.Stability
import src.LDL

/-!
# Solving `A x = b` with two triangular solves

`fwdSolve` solves `L y = b` by forward substitution and `bwdSolve` solves `U x = y` by
back substitution, both in the arithmetic `o`. Forward substitution computes
`yᵢ = ((bᵢ - lᵢ₀ y₀) - lᵢ₁ y₁ - ⋯) / lᵢᵢ`, subtracting in increasing order of `j`; back
substitution is forward substitution on the reversed system.

* `fwdSolve_backward`, `bwdSolve_backward`: `(T + ΔT) x̂ = b` with `|ΔT| ≤ γ_n |T|`.
* `exact_solve`: in exact arithmetic `L y = b` and `Lᵀ x = y` give `A x = b`.
* `solve_backward_error`: for the computed solution of `A x = b`,
  `(A + ΔA) x̂ = b` with `|ΔA| ≤ (γ_{n+1} + 2γ_n + γ_n²) |L̂| |L̂|ᵀ ≤ γ_{3n+1} |L̂| |L̂|ᵀ`.
-/
namespace Cholesky
open Matrix

variable {u : ℝ}

/-- Forward substitution for `L y = b`. -/
def fwdSolve (o : Ops) : {n : ℕ} → Matrix (Fin n) (Fin n) ℝ → (Fin n → ℝ) → Fin n → ℝ
  | 0, _, _ => 0
  | _ + 1, L, b =>
    Fin.cons (o.div (b 0) (L 0 0))
      (fwdSolve o (trailing L) fun i => o.upd (b i.succ) (L i.succ 0) (o.div (b 0) (L 0 0)))

/-- Back substitution for an upper triangular `U x = y`: forward substitution on the
reversed system. -/
def bwdSolve (o : Ops) {n : ℕ} (U : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ) : Fin n → ℝ :=
  fun i => fwdSolve o (U.submatrix Fin.rev Fin.rev) (fun j => y (Fin.rev j)) (Fin.rev i)

/-- The computed solution of `A x = b`: `L̂ ŷ = b`, then `L̂ᵀ x̂ = ŷ`. -/
def choleskySolve (o : Ops) {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) :
    Fin n → ℝ :=
  bwdSolve o (factor o A)ᵀ (fwdSolve o (factor o A) b)

lemma abs_inv_one_add_le (hu1 : u < 1) {δ : ℝ} (hδ : |δ| ≤ u) :
    |1 / (1 + δ)| ≤ 1 + gamma u 1 := by
  have h := abs_inv_one_add_sub_one_le hu1 hδ
  calc |1 / (1 + δ)| = |(1 / (1 + δ) - 1) + 1| := by ring_nf
    _ ≤ |1 / (1 + δ) - 1| + |1| := abs_add_le _ _
    _ ≤ gamma u 1 + 1 := by rw [abs_one]; linarith
    _ = 1 + gamma u 1 := by ring

/-- Backward error of forward substitution. -/
theorem fwdSolve_backward (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u) :
    ∀ {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ), IsLowerTriangular L →
      (∀ i, L i i ≠ 0) → (n : ℝ) * u < 1 →
      ∃ ΔL : Matrix (Fin n) (Fin n) ℝ, (L + ΔL) *ᵥ fwdSolve o L b = b ∧
        ∀ i j, |ΔL i j| ≤ gamma u n * |L i j|
  | 0, _, b, _, _, _ => ⟨0, by ext i; exact Fin.elim0 i, fun i => Fin.elim0 i⟩
  | n + 1, L, b, hL, hd, hnu => by
    push_cast at hnu
    have hnR : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have hu1 : u < 1 := by nlinarith
    have hnu' : (n : ℝ) * u < 1 := by nlinarith
    set y0 := o.div (b 0) (L 0 0) with hy0
    obtain ⟨δ0, hδ0, hy0'⟩ := ho.div (b 0) (L 0 0)
    choose δ ε hδ hε hb' using fun i : Fin n => ho.upd (b i.succ) (L i.succ 0) y0
    obtain ⟨ΔS, hS, hSb⟩ := fwdSolve_backward o ho hu (trailing L)
      (fun i => o.upd (b i.succ) (L i.succ 0) y0) (trailing_lower L hL)
      (fun i => hd i.succ) hnu'
    set y' := fwdSolve o (trailing L) (fun i => o.upd (b i.succ) (L i.succ 0) y0)
    have hy : fwdSolve o L b = Fin.cons y0 y' := rfl
    let ΔL : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ := Matrix.of fun i j =>
      Fin.cases (motive := fun _ => ℝ)
        (Fin.cases (motive := fun _ => ℝ) (L 0 0 * (1 / (1 + δ0) - 1)) (fun _ => 0) j)
        (fun i' => Fin.cases (motive := fun _ => ℝ) (L i'.succ 0 * ε i')
          (fun j' => (trailing L + ΔS) i' j' / (1 + δ i') - trailing L i' j') j) i
    have hδ0' := abs_le.mp hδ0
    have h1δ0 : 1 + δ0 ≠ 0 := by linarith
    have hγ1 : gamma u 1 ≤ gamma u (n + 1) := gamma_mono hu (by omega) (by push_cast; linarith)
    have hγn : gamma u n ≤ gamma u (n + 1) := gamma_mono hu (by omega) (by push_cast; linarith)
    have hgN : 0 ≤ gamma u n := gamma_nonneg hu hnu'
    have hsum := gamma_add_le hu (a := 1) (b := n) (by push_cast; linarith)
    rw [show 1 + n = n + 1 by ring] at hsum
    refine ⟨ΔL, ?_, ?_⟩
    · rw [hy]
      ext i
      refine Fin.cases ?_ (fun i => ?_) i
      · simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, Fin.cons_zero, Fin.cons_succ,
          Matrix.add_apply, ΔL, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ]
        have : ∀ j : Fin n, L 0 j.succ = 0 := fun j => hL 0 j.succ (Fin.succ_pos j)
        simp only [this, zero_mul, Finset.sum_const_zero, add_zero]
        rw [hy0, hy0']
        have h00 := hd 0
        field_simp
        ring
      · have hrow := congrFun hS i
        simp only [Matrix.mulVec, dotProduct, Matrix.add_apply] at hrow
        simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_succ, Fin.cons_zero, Fin.cons_succ,
          Matrix.add_apply, ΔL, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ]
        have hδi := abs_le.mp (hδ i)
        have h1δ : 1 + δ i ≠ 0 := by linarith
        have : ∀ j : Fin n, (L i.succ j.succ + ((trailing L i j + ΔS i j) / (1 + δ i) -
            trailing L i j)) * y' j = (trailing L i j + ΔS i j) * y' j / (1 + δ i) := by
          intro j; simp only [trailing, Matrix.submatrix_apply]; field_simp; ring
        simp only [this, ← Finset.sum_div]
        rw [hrow, hb' i]
        field_simp
        ring
    · intro i j
      refine Fin.cases ?_ (fun i => ?_) i <;> refine Fin.cases ?_ (fun j => ?_) j
      · simp only [ΔL, Matrix.of_apply, Fin.cases_zero]
        rw [abs_mul, mul_comm]
        exact mul_le_mul_of_nonneg_right
          (le_trans (by rw [abs_sub_comm, show (1 : ℝ) - 1 / (1 + δ0) = -(1 / (1 + δ0) - 1) by
            ring, abs_neg]; exact abs_inv_one_add_sub_one_le hu1 hδ0) hγ1) (abs_nonneg _)
      · simp only [ΔL, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ, abs_zero]
        exact mul_nonneg (le_trans hgN hγn) (abs_nonneg _)
      · simp only [ΔL, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ]
        rw [abs_mul, mul_comm]
        exact mul_le_mul_of_nonneg_right
          (le_trans (hε i) (le_gamma hu (by omega) (by push_cast; linarith))) (abs_nonneg _)
      · simp only [ΔL, Matrix.of_apply, Fin.cases_succ]
        have hδi := abs_le.mp (hδ i)
        have h1δ : (0 : ℝ) < 1 + δ i := by linarith
        have hb := hSb i j
        have e : (trailing L + ΔS) i j / (1 + δ i) - trailing L i j =
            trailing L i j * (1 / (1 + δ i) - 1) + ΔS i j * (1 / (1 + δ i)) := by
          simp only [Matrix.add_apply]; field_simp; ring
        rw [e]
        have hT : trailing L i j = L i.succ j.succ := rfl
        calc |trailing L i j * (1 / (1 + δ i) - 1) + ΔS i j * (1 / (1 + δ i))|
            ≤ |trailing L i j| * gamma u 1 + gamma u n * |trailing L i j| * (1 + gamma u 1) := by
              refine le_trans (abs_add_le _ _) (add_le_add ?_ ?_)
              · rw [abs_mul]
                exact mul_le_mul_of_nonneg_left (abs_inv_one_add_sub_one_le hu1 (hδ i))
                  (abs_nonneg _)
              · rw [abs_mul]
                exact mul_le_mul hb (abs_inv_one_add_le hu1 (hδ i)) (abs_nonneg _)
                  (mul_nonneg hgN (abs_nonneg _))
          _ = (gamma u 1 + gamma u n + gamma u 1 * gamma u n) * |trailing L i j| := by ring
          _ ≤ gamma u (n + 1) * |L i.succ j.succ| := by
              rw [hT]; exact mul_le_mul_of_nonneg_right hsum (abs_nonneg _)

lemma rev_lower {n : ℕ} {U : Matrix (Fin n) (Fin n) ℝ} (hU : IsUpperTriangular U) :
    IsLowerTriangular (U.submatrix Fin.rev Fin.rev) := by
  intro i j hij
  exact hU _ _ (Fin.rev_lt_rev.mpr hij)

/-- Backward error of back substitution. -/
theorem bwdSolve_backward (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u) {n : ℕ}
    (U : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ) (hU : IsUpperTriangular U)
    (hd : ∀ i, U i i ≠ 0) (hnu : (n : ℝ) * u < 1) :
    ∃ ΔU : Matrix (Fin n) (Fin n) ℝ, (U + ΔU) *ᵥ bwdSolve o U y = y ∧
      ∀ i j, |ΔU i j| ≤ gamma u n * |U i j| := by
  obtain ⟨ΔU', h1, h2⟩ := fwdSolve_backward o ho hu (U.submatrix Fin.rev Fin.rev)
    (fun j => y (Fin.rev j)) (rev_lower hU) (fun i => hd _) hnu
  refine ⟨ΔU'.submatrix Fin.rev Fin.rev, ?_, fun i j => by simpa using h2 (Fin.rev i) (Fin.rev j)⟩
  ext i
  have := congrFun h1 (Fin.rev i)
  simp only [Matrix.mulVec, dotProduct, Matrix.add_apply, Matrix.submatrix_apply,
    Fin.rev_rev] at this ⊢
  rw [← this]
  exact (Fintype.sum_equiv Fin.revPerm _ _ (fun j => by simp [bwdSolve, Fin.rev_rev])).symm

/-! ## Exact arithmetic -/

theorem exact_fwdSolve {n : ℕ} (L : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ)
    (hL : IsLowerTriangular L) (hd : ∀ i, L i i ≠ 0) : L *ᵥ fwdSolve Ops.exact L b = b := by
  obtain ⟨ΔL, h1, h2⟩ := fwdSolve_backward Ops.exact (Ops.stdModel_exact le_rfl) le_rfl L b hL hd
    (by simp)
  have : ΔL = 0 := by
    ext i j
    have := h2 i j
    simp [gamma] at this
    simpa using this
  simpa [this] using h1

theorem exact_bwdSolve {n : ℕ} (U : Matrix (Fin n) (Fin n) ℝ) (y : Fin n → ℝ)
    (hU : IsUpperTriangular U) (hd : ∀ i, U i i ≠ 0) : U *ᵥ bwdSolve Ops.exact U y = y := by
  obtain ⟨ΔU, h1, h2⟩ := bwdSolve_backward Ops.exact (Ops.stdModel_exact le_rfl) le_rfl U y hU hd
    (by simp)
  have : ΔU = 0 := by
    ext i j
    have := h2 i j
    simp [gamma] at this
    simpa using this
  simpa [this] using h1

lemma transpose_upper {n : ℕ} {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L) :
    IsUpperTriangular Lᵀ := fun i j hij => hL j i hij

/-- In exact arithmetic, `L y = b` and `Lᵀ x = y` give `A x = b`. -/
theorem exact_solve {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosDef) (b : Fin n → ℝ) :
    A *ᵥ choleskySolve Ops.exact A b = b := by
  have hsymm : Aᵀ = A := ((posDef_iff_real A).mp hA).1
  obtain ⟨hL, hd, hLL⟩ := exact_factor_spec A (by rwa [symmOfLower_eq_self hsymm])
  rw [symmOfLower_eq_self hsymm] at hLL
  have hd' : ∀ i, factor Ops.exact A i i ≠ 0 := fun i => ne_of_gt (hd i)
  calc A *ᵥ choleskySolve Ops.exact A b
      = (factor Ops.exact A * (factor Ops.exact A)ᵀ) *ᵥ choleskySolve Ops.exact A b := by
        rw [← hLL]
    _ = b := by
        unfold choleskySolve
        rw [← Matrix.mulVec_mulVec, exact_bwdSolve _ _ (transpose_upper hL) (fun i => hd' i),
          exact_fwdSolve _ _ hL hd']

/-! ## The computed solution -/

lemma abs_mul_le_of {n : ℕ} {M N X Y : Matrix (Fin n) (Fin n) ℝ} {a b : ℝ} (ha : 0 ≤ a)
    (hM : ∀ i j, |M i j| ≤ a * |X i j|) (hN : ∀ i j, |N i j| ≤ b * |Y i j|)
    (i j : Fin n) : |(M * N) i j| ≤ a * b * (absMat X * absMat Y) i j := by
  rw [Matrix.mul_apply, Matrix.mul_apply, Finset.mul_sum]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum fun k _ => ?_)
  rw [abs_mul]
  calc |M i k| * |N k j| ≤ (a * |X i k|) * (b * |Y k j|) :=
        mul_le_mul (hM i k) (hN k j) (abs_nonneg _) (mul_nonneg ha (abs_nonneg _))
    _ = a * b * (absMat X i k * absMat Y k j) := by simp; ring

/-- `γ_{n+1} + 2γ_n + γ_n² ≤ γ_{3n+1}`. -/
lemma gamma_solve_le (hu : 0 ≤ u) (n : ℕ) (h : ((3 * n + 1 : ℕ) : ℝ) * u < 1) :
    gamma u (n + 1) + 2 * gamma u n + gamma u n ^ 2 ≤ gamma u (3 * n + 1) := by
  have h2n := gamma_add_le hu (a := n) (b := n) (by push_cast at h ⊢; nlinarith)
  have h3n := gamma_add_le hu (a := n + 1) (b := n + n) (by push_cast at h ⊢; linarith)
  rw [show n + 1 + (n + n) = 3 * n + 1 by ring] at h3n
  have hp : 0 ≤ gamma u (n + 1) * gamma u (n + n) :=
    mul_nonneg (gamma_nonneg hu (by push_cast at h ⊢; nlinarith))
      (gamma_nonneg hu (by push_cast at h ⊢; nlinarith))
  nlinarith

/-- **Backward error of the computed solution.** -/
theorem solve_backward_error (o : Ops) (ho : o.StdModel u) (hu : 0 ≤ u) {n : ℕ}
    (A : Matrix (Fin n) (Fin n) ℝ) (b : Fin n → ℝ) (hnu : ((3 * n + 1 : ℕ) : ℝ) * u < 1)
    (h : Succeeds o A) :
    ∃ ΔA : Matrix (Fin n) (Fin n) ℝ, (symmOfLower A + ΔA) *ᵥ choleskySolve o A b = b ∧
      ∀ i j, |ΔA i j| ≤ (gamma u (n + 1) + 2 * gamma u n + gamma u n ^ 2) *
        (absMat (factor o A) * (absMat (factor o A))ᵀ) i j := by
  have hnR : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hn1 : ((n + 1 : ℕ) : ℝ) * u < 1 := by push_cast at hnu ⊢; nlinarith
  have hn0 : (n : ℝ) * u < 1 := by push_cast at hnu; nlinarith
  obtain ⟨hAE, -, hL, hd, hE⟩ := backward_error o ho hu A hn1 h
  set L := factor o A
  have hd' : ∀ i, L i i ≠ 0 := fun i => ne_of_gt (hd i)
  set y := fwdSolve o L b
  obtain ⟨Δ1, h1, h1b⟩ := fwdSolve_backward o ho hu L b hL hd' hn0
  obtain ⟨Δ2, h2, h2b⟩ := bwdSolve_backward o ho hu Lᵀ y (transpose_upper hL) (fun i => hd' i) hn0
  set x := bwdSolve o Lᵀ y
  refine ⟨(L + Δ1) * (Lᵀ + Δ2) - symmOfLower A, ?_, ?_⟩
  · have hx : choleskySolve o A b = x := rfl
    rw [hx, add_sub_cancel, ← Matrix.mulVec_mulVec, h2, h1]
  · intro i j
    have hγ : 0 ≤ gamma u n := gamma_nonneg hu hn0
    have hexp : (L + Δ1) * (Lᵀ + Δ2) - symmOfLower A =
        (L * Lᵀ - symmOfLower A) + L * Δ2 + Δ1 * Lᵀ + Δ1 * Δ2 := by
      simp only [Matrix.add_mul, Matrix.mul_add]; abel
    rw [hexp]
    have hLabs : ∀ i j, |L i j| ≤ 1 * |L i j| := fun i j => by rw [one_mul]
    have hLTabs : ∀ i j, |Lᵀ i j| ≤ 1 * |Lᵀ i j| := fun i j => by rw [one_mul]
    have e1 := hE i j
    have e2 := abs_mul_le_of zero_le_one hLabs h2b i j
    have e3 := abs_mul_le_of hγ h1b hLTabs i j
    have e4 := abs_mul_le_of hγ h1b h2b i j
    simp only [one_mul, mul_one, absMat_transpose] at e2 e3 e4
    calc |((L * Lᵀ - symmOfLower A) + L * Δ2 + Δ1 * Lᵀ + Δ1 * Δ2) i j|
        ≤ |(L * Lᵀ - symmOfLower A) i j| + |(L * Δ2) i j| + |(Δ1 * Lᵀ) i j| +
            |(Δ1 * Δ2) i j| := by
          simp only [Matrix.add_apply]
          exact le_trans (abs_add_le _ _) (add_le_add (le_trans (abs_add_le _ _)
            (add_le_add (abs_add_le _ _) le_rfl)) le_rfl)
      _ ≤ gamma u (n + 1) * (absMat L * (absMat L)ᵀ) i j +
            gamma u n * (absMat L * (absMat L)ᵀ) i j +
            gamma u n * (absMat L * (absMat L)ᵀ) i j +
            gamma u n * gamma u n * (absMat L * (absMat L)ᵀ) i j := by
          gcongr
      _ = (gamma u (n + 1) + 2 * gamma u n + gamma u n ^ 2) *
            (absMat L * (absMat L)ᵀ) i j := by ring

end Cholesky
