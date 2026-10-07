import src.Algorithm

/-!
# Existence and uniqueness of the Cholesky factorization

* `exact_factor_mul_transpose`: if the exact in-place algorithm completes, its output
  `L` satisfies `L Lᵀ = A` (for the symmetric matrix defined by the lower triangle).
* `exact_succeeds_of_posDef`: for an SPD matrix every exact pivot is positive.
* `posDef_of_factorization`: the converse direction of the theorem,
  `xᵀ L Lᵀ x = ‖Lᵀ x‖² > 0`.
* `cholesky_unique`: the factor with positive diagonal is unique.
* `cholesky_iff`: the theorem `thm:cholesky_existence` of the notes.
-/
namespace Cholesky
open Matrix

variable {n : ℕ}

/-- The Schur complement computed by one exact step. -/
lemma symmOfLower_trail_exact (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) (hA : 0 < A 0 0) :
    symmOfLower (trail Ops.exact A) = schur (symmOfLower A) := by
  have hs : Real.sqrt (A 0 0) ^ 2 = A 0 0 := Real.sq_sqrt hA.le
  have hs0 : Real.sqrt (A 0 0) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hA)
  ext i j
  simp only [symmOfLower, schur, Matrix.of_apply, Fin.succ_le_succ_iff]
  have h0i : (0 : Fin (n + 1)) ≤ i.succ := Fin.zero_le _
  have h0j : (0 : Fin (n + 1)) ≤ j.succ := Fin.zero_le _
  simp only [h0i, h0j, if_true, le_refl]
  split_ifs with hij
  · rw [trail_apply _ _ _ _ hij]
    simp only [Ops.exact, firstCol]
    field_simp
    rw [hs]; ring
  · rw [trail_apply _ _ _ _ (le_of_lt (not_le.mp hij))]
    simp only [Ops.exact, firstCol]
    field_simp
    rw [hs]
    ring

/-- If the exact algorithm completes, its output is a Cholesky factor. -/
theorem exact_factor_mul_transpose :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ), Succeeds Ops.exact A →
      factor Ops.exact A * (factor Ops.exact A)ᵀ = symmOfLower A
  | 0, A, _ => by ext i; exact Fin.elim0 i
  | n + 1, A, h => by
    rw [succeeds_succ_iff] at h
    obtain ⟨h0, htrail⟩ := h
    have ih := exact_factor_mul_transpose (trail Ops.exact A) htrail
    rw [symmOfLower_trail_exact A h0] at ih
    have hs : Real.sqrt (A 0 0) ^ 2 = A 0 0 := Real.sq_sqrt h0.le
    have hs0 : Real.sqrt (A 0 0) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr h0)
    have hL := lowerPart_lower (run Ops.exact (n + 1) A)
    set L := factor Ops.exact A with hLdef
    have hLl : IsLowerTriangular L := hL
    ext i j
    refine Fin.cases ?_ (fun i => ?_) i <;> refine Fin.cases ?_ (fun j => ?_) j
    · rw [mul_transpose_zero_zero L hLl, hLdef, factor_zero_zero]
      simp [symmOfLower, Ops.exact, hs]
    · rw [← Matrix.transpose_apply (L * Lᵀ), Matrix.transpose_mul, Matrix.transpose_transpose,
        mul_transpose_succ_zero L hLl, hLdef, factor_zero_zero, factor_succ_zero]
      simp only [symmOfLower, Matrix.of_apply, Ops.exact, firstCol]
      rw [if_neg (not_le.mpr (Fin.succ_pos j))]
      field_simp
    · rw [mul_transpose_succ_zero L hLl, hLdef, factor_zero_zero, factor_succ_zero]
      simp only [symmOfLower, Matrix.of_apply, Ops.exact, firstCol]
      rw [if_pos (Fin.zero_le _)]
      field_simp
    · rw [mul_transpose_succ_succ L i j, hLdef, trailing_factor, ih, factor_succ_zero,
        factor_succ_zero]
      simp only [symmOfLower, schur, Matrix.of_apply, Fin.succ_le_succ_iff, Ops.exact,
        firstCol, if_pos (Fin.zero_le _), le_refl, if_true]
      field_simp
      rw [hs]
      ring

/-- For an SPD matrix, every pivot of the exact algorithm is positive. -/
theorem exact_succeeds_of_posDef :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ), (symmOfLower A).PosDef → Succeeds Ops.exact A
  | 0, _, _ => fun k => Fin.elim0 k
  | n + 1, A, h => by
    have h0 : 0 < A 0 0 := by
      have := h.2 (Pi.single 0 1) (by simp)
      simpa [dotProduct, Matrix.mulVec, Pi.single_apply, symmOfLower] using this
    rw [succeeds_succ_iff]
    refine ⟨h0, exact_succeeds_of_posDef _ ?_⟩
    rw [symmOfLower_trail_exact A h0]
    exact schur_posDef _ h

/-- If the arithmetic returns a positive square root of a positive number, a completed
run produces a factor with positive diagonal. -/
theorem factor_diag_pos (o : Ops) (ho : ∀ a, 0 < a → 0 < o.sqrt a) :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ), Succeeds o A → ∀ i, 0 < factor o A i i
  | 0, _, _, i => Fin.elim0 i
  | n + 1, A, h, i => by
    rw [succeeds_succ_iff] at h
    refine Fin.cases ?_ (fun i => ?_) i
    · rw [factor_zero_zero]; exact ho _ h.1
    · rw [factor_succ_succ]; exact factor_diag_pos o ho _ h.2 i

lemma exact_sqrt_pos (a : ℝ) (ha : 0 < a) : 0 < Ops.exact.sqrt a := Real.sqrt_pos.mpr ha

/-! ## The theorem -/

lemma isLowerTriangular_blockTriangular {L : Matrix (Fin n) (Fin n) ℝ}
    (hL : IsLowerTriangular L) : L.BlockTriangular OrderDual.toDual := by
  intro i j hij
  exact hL i j hij

lemma det_ne_zero_of_lower {L : Matrix (Fin n) (Fin n) ℝ} (hL : IsLowerTriangular L)
    (hdiag : ∀ i, L i i ≠ 0) : L.det ≠ 0 := by
  rw [Matrix.det_of_lowerTriangular L (isLowerTriangular_blockTriangular hL)]
  exact Finset.prod_ne_zero_iff.mpr (fun i _ => hdiag i)

lemma quad_mul_transpose (L : Matrix (Fin n) (Fin n) ℝ) (x : Fin n → ℝ) :
    x ⬝ᵥ ((L * Lᵀ) *ᵥ x) = (Lᵀ *ᵥ x) ⬝ᵥ (Lᵀ *ᵥ x) := by
  rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose]

/-- The converse direction: `xᵀ L Lᵀ x = ‖Lᵀ x‖² > 0` because `L` is invertible. -/
theorem posDef_of_factorization (L : Matrix (Fin n) (Fin n) ℝ) (hL : IsLowerTriangular L)
    (hdiag : ∀ i, 0 < L i i) : (L * Lᵀ).PosDef := by
  rw [posDef_iff_real]
  refine ⟨by rw [Matrix.transpose_mul, Matrix.transpose_transpose], fun x hx => ?_⟩
  rw [quad_mul_transpose]
  have hdet : Lᵀ.det ≠ 0 := by
    rw [Matrix.det_transpose]
    exact det_ne_zero_of_lower hL (fun i => ne_of_gt (hdiag i))
  have hinj : Function.Injective (Lᵀ *ᵥ ·) := by
    rw [Matrix.mulVec_injective_iff_isUnit, Matrix.isUnit_iff_isUnit_det]
    exact isUnit_iff_ne_zero.mpr hdet
  have hy : Lᵀ *ᵥ x ≠ 0 := by
    intro h
    apply hx
    apply hinj
    simp only [h, Matrix.mulVec_zero]
  simpa using Matrix.dotProduct_star_self_pos_iff.mpr hy

/-- Uniqueness of the Cholesky factor. -/
theorem cholesky_unique :
    ∀ {n : ℕ} (L M : Matrix (Fin n) (Fin n) ℝ), IsLowerTriangular L → IsLowerTriangular M →
      (∀ i, 0 < L i i) → (∀ i, 0 < M i i) → L * Lᵀ = M * Mᵀ → L = M
  | 0, L, M, _, _, _, _, _ => by ext i; exact Fin.elim0 i
  | n + 1, L, M, hL, hM, hLd, hMd, h => by
    obtain ⟨l0, lc, ls⟩ := block_equations L (L * Lᵀ) hL (hLd 0) rfl
    obtain ⟨m0, mc, ms⟩ := block_equations M (L * Lᵀ) hM (hMd 0) h
    have h00 : L 0 0 = M 0 0 := by rw [l0, m0]
    have hc : ∀ i : Fin n, L i.succ 0 = M i.succ 0 := by
      intro i; rw [lc i, mc i, h00]
    have htr : trailing L = trailing M :=
      cholesky_unique (trailing L) (trailing M) (trailing_lower L hL) (trailing_lower M hM)
        (fun i => hLd i.succ) (fun i => hMd i.succ) (by rw [ls, ms])
    ext i j
    refine Fin.cases ?_ (fun i => ?_) i <;> refine Fin.cases ?_ (fun j => ?_) j
    · exact h00
    · rw [hL 0 j.succ (Fin.succ_pos j), hM 0 j.succ (Fin.succ_pos j)]
    · exact hc i
    · exact congrFun (congrFun htr i) j

/-- Existence of the Cholesky factorization, computed by the exact algorithm. -/
theorem cholesky_exists (A : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosDef) :
    ∃ L : Matrix (Fin n) (Fin n) ℝ, IsLowerTriangular L ∧ (∀ i, 0 < L i i) ∧ A = L * Lᵀ := by
  have hsymm : Aᵀ = A := ((posDef_iff_real A).mp hA).1
  have hA' : (symmOfLower A).PosDef := by rwa [symmOfLower_eq_self hsymm]
  have hsucc := exact_succeeds_of_posDef A hA'
  refine ⟨factor Ops.exact A, lowerPart_lower _,
    factor_diag_pos Ops.exact exact_sqrt_pos A hsucc, ?_⟩
  rw [exact_factor_mul_transpose A hsucc, symmOfLower_eq_self hsymm]

/-- **Theorem (existence and uniqueness).** A real symmetric matrix is positive definite
if and only if it has a factorization `A = L Lᵀ` with `L` lower triangular and
`lᵢᵢ > 0`; this factor is unique. -/
theorem cholesky_iff (A : Matrix (Fin n) (Fin n) ℝ) (hA : Aᵀ = A) :
    A.PosDef ↔
      ∃ L : Matrix (Fin n) (Fin n) ℝ, IsLowerTriangular L ∧ (∀ i, 0 < L i i) ∧ A = L * Lᵀ := by
  constructor
  · exact cholesky_exists A
  · rintro ⟨L, hL, hd, rfl⟩
    exact posDef_of_factorization L hL hd

theorem cholesky_existsUnique (A : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosDef) :
    ∃! L : Matrix (Fin n) (Fin n) ℝ, IsLowerTriangular L ∧ (∀ i, 0 < L i i) ∧ A = L * Lᵀ := by
  obtain ⟨L, hL⟩ := cholesky_exists A hA
  refine ⟨L, hL, fun M hM => ?_⟩
  exact cholesky_unique M L hM.1 hL.1 hM.2.1 hL.2.1 (by rw [← hM.2.2, ← hL.2.2])

/-! ## The exact in-place algorithm -/

/-- In exact arithmetic the routine completes if and only if the matrix defined by the
lower triangle is SPD. -/
theorem exact_succeeds_iff (A : Matrix (Fin n) (Fin n) ℝ) :
    Succeeds Ops.exact A ↔ (symmOfLower A).PosDef := by
  constructor
  · intro h
    rw [← exact_factor_mul_transpose A h]
    exact posDef_of_factorization _ (lowerPart_lower _) (factor_diag_pos _ exact_sqrt_pos A h)
  · exact exact_succeeds_of_posDef A

/-- At completion the lower triangle holds the Cholesky factor `L`: it is the unique lower
triangular factor with positive diagonal. -/
theorem exact_factor_spec (A : Matrix (Fin n) (Fin n) ℝ) (hA : (symmOfLower A).PosDef) :
    IsLowerTriangular (factor Ops.exact A) ∧ (∀ i, 0 < factor Ops.exact A i i) ∧
      symmOfLower A = factor Ops.exact A * (factor Ops.exact A)ᵀ := by
  have h := exact_succeeds_of_posDef A hA
  exact ⟨lowerPart_lower _, factor_diag_pos _ exact_sqrt_pos A h,
    (exact_factor_mul_transpose A h).symm⟩

end Cholesky
