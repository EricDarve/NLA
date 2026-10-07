import Mathlib

/-!
# The complex Hermitian case

A complex matrix is Hermitian positive definite if `Aᴴ = A` and `xᴴ A x > 0` for every
nonzero `x` (mathlib's `Matrix.PosDef` with the order `ComplexOrder`, in which `0 < z`
means that `z` is real and positive). The factorization is `A = L Lᴴ` with `L` lower
triangular and positive real diagonal, and it is unique (`cholesky_iff`,
`cholesky_unique`). The proof is the same block recursion as in the real case, with
conjugate transposes.
-/
open scoped ComplexOrder

namespace Cholesky
namespace Herm
open Matrix

variable {n : ℕ}

/-- `L` is lower triangular. -/
def IsLower {m : ℕ} (L : Matrix (Fin m) (Fin m) ℂ) : Prop := ∀ i j, i < j → L i j = 0

lemma pos_real {z : ℂ} (hz : 0 < z) : z = ((z.re : ℝ) : ℂ) ∧ 0 < z.re := by
  obtain ⟨h1, h2⟩ := Complex.pos_iff.mp hz
  exact ⟨Complex.ext (by simp) (by simp [← h2]), h1⟩

lemma star_of_pos {z : ℂ} (hz : 0 < z) : star z = z := by
  obtain ⟨h1, -⟩ := pos_real hz
  rw [h1, Complex.star_def, Complex.conj_ofReal]

/-- The trailing block. -/
def trail (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  A.submatrix Fin.succ Fin.succ

section Blocks
variable (L : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ)

lemma mulH_zero_zero (hL : IsLower L) : (L * Lᴴ) 0 0 = L 0 0 * star (L 0 0) := by
  simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, Fin.sum_univ_succ]
  have : ∀ k : Fin n, L 0 k.succ = 0 := fun k => hL 0 k.succ (Fin.succ_pos k)
  simp [this]

lemma mulH_succ_zero (hL : IsLower L) (i : Fin n) :
    (L * Lᴴ) i.succ 0 = L i.succ 0 * star (L 0 0) := by
  simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, Fin.sum_univ_succ]
  have : ∀ k : Fin n, L 0 k.succ = 0 := fun k => hL 0 k.succ (Fin.succ_pos k)
  simp [this]

lemma mulH_succ_succ (i j : Fin n) :
    (L * Lᴴ) i.succ j.succ = L i.succ 0 * star (L j.succ 0) + (trail L * (trail L)ᴴ) i j := by
  simp [Matrix.mul_apply, Fin.sum_univ_succ, trail]

lemma trail_lower (hL : IsLower L) : IsLower (trail L) := fun _ _ hij =>
  hL _ _ (Fin.succ_lt_succ_iff.mpr hij)

end Blocks

/-! ## The Schur complement -/

/-- `S = B - c cᴴ / a₁₁`. -/
noncomputable def schur (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  Matrix.of fun i j => A i.succ j.succ - A i.succ 0 * A 0 j.succ / A 0 0

lemma diag_pos {A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ} (hA : A.PosDef) : 0 < A 0 0 := by
  have := hA.2 (Pi.single 0 1) (by simp)
  simpa [dotProduct, Matrix.mulVec, Pi.single_apply] using this

/-- With `t = -(cᴴ y)/a₁₁` and `x = (t, y)`, `xᴴ A x = yᴴ S y`. -/
lemma quad_schur (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ) (ha : A 0 0 ≠ 0) (y : Fin n → ℂ) :
    star (Fin.cons (-(∑ j, A 0 j.succ * y j) / A 0 0) y : Fin (n + 1) → ℂ) ⬝ᵥ
        (A *ᵥ Fin.cons (-(∑ j, A 0 j.succ * y j) / A 0 0) y) = star y ⬝ᵥ (schur A *ᵥ y) := by
  set c := ∑ j, A 0 j.succ * y j
  set t := -c / A 0 0
  have ht : A 0 0 * t + c = 0 := by simp only [t]; field_simp; ring
  simp only [dotProduct, Matrix.mulVec, Fin.sum_univ_succ, Fin.cons_zero, Fin.cons_succ,
    Pi.star_apply, schur, Matrix.of_apply]
  rw [show A 0 0 * t + ∑ x, A 0 x.succ * y x = A 0 0 * t + c from rfl, ht, mul_zero, zero_add]
  have e : ∀ i : Fin n, star (y i) * (A i.succ 0 * t + ∑ j, A i.succ j.succ * y j) =
      t * (star (y i) * A i.succ 0) + star (y i) * ∑ j, A i.succ j.succ * y j := by
    intro i; ring
  simp only [e, Finset.sum_add_distrib, ← Finset.mul_sum]
  have e2 : ∀ i : Fin n, star (y i) * ∑ j, (A i.succ j.succ - A i.succ 0 * A 0 j.succ / A 0 0) * y j =
      star (y i) * ∑ j, A i.succ j.succ * y j - (star (y i) * A i.succ 0) * c / A 0 0 := by
    intro i
    have h1 : (star (y i) * A i.succ 0) * c / A 0 0 =
        ∑ j, star (y i) * (A i.succ 0 * A 0 j.succ / A 0 0 * y j) := by
      simp only [c, Finset.mul_sum, Finset.sum_div]
      apply Finset.sum_congr rfl; intro j _; ring
    rw [h1, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl; intro j _; ring
  simp only [e2, Finset.sum_sub_distrib, ← Finset.sum_div, ← Finset.sum_mul]
  simp only [t]
  field_simp
  ring

/-- The Schur complement of a Hermitian positive definite matrix is Hermitian positive
definite. -/
theorem schur_posDef (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ) (hA : A.PosDef) :
    (schur A).PosDef := by
  have ha := diag_pos hA
  have hstar := star_of_pos ha
  refine ⟨?_, fun y hy => ?_⟩
  · ext i j
    simp only [schur, Matrix.conjTranspose_apply, Matrix.of_apply, star_sub, star_div₀, star_mul',
      hstar]
    rw [hA.1.apply i.succ j.succ, hA.1.apply 0 j.succ, hA.1.apply i.succ 0]
    ring
  · have hx : (Fin.cons (-(∑ j, A 0 j.succ * y j) / A 0 0) y : Fin (n + 1) → ℂ) ≠ 0 := by
      intro h; apply hy; funext i; simpa using congrFun h i.succ
    have := hA.2 _ hx
    rwa [quad_schur A (ne_of_gt ha)] at this

/-! ## Existence, the converse, and uniqueness -/

theorem cholesky_exists : ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℂ), A.PosDef →
    ∃ L : Matrix (Fin n) (Fin n) ℂ, IsLower L ∧ (∀ i, 0 < L i i) ∧ A = L * Lᴴ
  | 0, A, _ => ⟨0, fun i => Fin.elim0 i, fun i => Fin.elim0 i, by ext i; exact Fin.elim0 i⟩
  | n + 1, A, hA => by
    have ha := diag_pos hA
    obtain ⟨har, hre⟩ := pos_real ha
    set d : ℂ := ((Real.sqrt (A 0 0).re : ℝ) : ℂ)
    have hd : 0 < d := Complex.zero_lt_real.mpr (Real.sqrt_pos.mpr hre)
    have hd0 : d ≠ 0 := ne_of_gt hd
    have hdd : d * star d = A 0 0 := by
      rw [star_of_pos hd, har]
      simp only [d, ← Complex.ofReal_mul, Real.mul_self_sqrt hre.le]
    obtain ⟨LS, hLS, hLSd, hSeq⟩ := cholesky_exists (schur A) (schur_posDef A hA)
    let L : Matrix (Fin (n + 1)) (Fin (n + 1)) ℂ := Matrix.of fun i j =>
      Fin.cases (motive := fun _ => ℂ) (Fin.cases (motive := fun _ => ℂ) d (fun _ => 0) j)
        (fun i' => Fin.cases (motive := fun _ => ℂ) (A i'.succ 0 / d) (fun j' => LS i' j') j) i
    have hL : IsLower L := by
      intro i j hij
      refine Fin.cases (fun hij => ?_) (fun i => fun hij => ?_) i hij
      · rcases Fin.eq_zero_or_eq_succ j with rfl | ⟨j, rfl⟩
        · exact absurd hij (lt_irrefl _)
        · simp [L]
      · rcases Fin.eq_zero_or_eq_succ j with rfl | ⟨j, rfl⟩
        · exact absurd hij (Fin.not_lt_zero _)
        · simp only [L, Matrix.of_apply, Fin.cases_succ]
          exact hLS i j (Fin.succ_lt_succ_iff.mp hij)
    refine ⟨L, hL, fun i => ?_, ?_⟩
    · refine Fin.cases ?_ (fun i => ?_) i
      · simpa [L] using hd
      · simpa [L] using hLSd i
    · have htr : trail L = LS := by ext i j; simp [trail, L]
      ext i j
      refine Fin.cases ?_ (fun i => ?_) i <;> refine Fin.cases ?_ (fun j => ?_) j
      · rw [mulH_zero_zero L hL]; simp only [L, Matrix.of_apply, Fin.cases_zero]; rw [hdd]
      · have hH : (L * Lᴴ)ᴴ = L * Lᴴ := by
          rw [Matrix.conjTranspose_mul, Matrix.conjTranspose_conjTranspose]
        rw [← congrFun (congrFun hH 0) j.succ, Matrix.conjTranspose_apply, mulH_succ_zero L hL]
        simp only [L, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ, star_mul', star_div₀,
          star_star]
        rw [star_of_pos hd, hA.1.apply 0 j.succ]
        field_simp
      · rw [mulH_succ_zero L hL]
        simp only [L, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ]
        rw [star_of_pos hd]; field_simp
      · rw [mulH_succ_succ L i j, htr, ← hSeq]
        simp only [schur, L, Matrix.of_apply, Fin.cases_zero, Fin.cases_succ, star_div₀]
        rw [star_of_pos hd, hA.1.apply 0 j.succ]
        have hdd' : A 0 0 = d * d := by rw [← hdd, star_of_pos hd]
        rw [hdd']
        field_simp
        ring

lemma det_ne_zero_of_lower {L : Matrix (Fin n) (Fin n) ℂ} (hL : IsLower L) (hd : ∀ i, L i i ≠ 0) :
    L.det ≠ 0 := by
  rw [Matrix.det_of_lowerTriangular L (fun i j hij => hL i j hij)]
  exact Finset.prod_ne_zero_iff.mpr (fun i _ => hd i)

/-- The converse: `L Lᴴ` is Hermitian positive definite. -/
theorem posDef_of_factorization (L : Matrix (Fin n) (Fin n) ℂ) (hL : IsLower L)
    (hd : ∀ i, 0 < L i i) : (L * Lᴴ).PosDef := by
  apply Matrix.PosDef.mul_conjTranspose_self
  rw [Matrix.vecMul_injective_iff_isUnit, Matrix.isUnit_iff_isUnit_det]
  exact isUnit_iff_ne_zero.mpr (det_ne_zero_of_lower hL (fun i => ne_of_gt (hd i)))

/-- Uniqueness of the factor with positive real diagonal. -/
theorem cholesky_unique : ∀ {n : ℕ} (L M : Matrix (Fin n) (Fin n) ℂ), IsLower L → IsLower M →
    (∀ i, 0 < L i i) → (∀ i, 0 < M i i) → L * Lᴴ = M * Mᴴ → L = M
  | 0, L, M, _, _, _, _, _ => by ext i; exact Fin.elim0 i
  | n + 1, L, M, hL, hM, hLd, hMd, h => by
    have h00 : L 0 0 = M 0 0 := by
      have e := congrFun (congrFun h 0) 0
      rw [mulH_zero_zero L hL, mulH_zero_zero M hM, star_of_pos (hLd 0), star_of_pos (hMd 0)] at e
      obtain ⟨a1, b1⟩ := pos_real (hLd 0)
      obtain ⟨a2, b2⟩ := pos_real (hMd 0)
      rw [a1, a2] at e ⊢
      rw [← Complex.ofReal_mul, ← Complex.ofReal_mul, Complex.ofReal_inj] at e
      congr 1; nlinarith
    have hne : L 0 0 ≠ 0 := ne_of_gt (hLd 0)
    have hc : ∀ i : Fin n, L i.succ 0 = M i.succ 0 := by
      intro i
      have e := congrFun (congrFun h i.succ) 0
      rw [mulH_succ_zero L hL, mulH_succ_zero M hM, ← h00, star_of_pos (hLd 0)] at e
      exact mul_right_cancel₀ hne e
    have htr : trail L * (trail L)ᴴ = trail M * (trail M)ᴴ := by
      ext i j
      have e := congrFun (congrFun h i.succ) j.succ
      rw [mulH_succ_succ L i j, mulH_succ_succ M i j, hc i, hc j] at e
      exact add_left_cancel e
    have ht := cholesky_unique (trail L) (trail M) (trail_lower L hL) (trail_lower M hM)
      (fun i => hLd i.succ) (fun i => hMd i.succ) htr
    ext i j
    refine Fin.cases ?_ (fun i => ?_) i <;> refine Fin.cases ?_ (fun j => ?_) j
    · exact h00
    · rw [hL 0 j.succ (Fin.succ_pos j), hM 0 j.succ (Fin.succ_pos j)]
    · exact hc i
    · exact congrFun (congrFun ht i) j

/-- **The complex case.** A Hermitian matrix is positive definite iff it has a factorization
`A = L Lᴴ` with `L` lower triangular with positive real diagonal; the factor is unique. -/
theorem cholesky_iff (A : Matrix (Fin n) (Fin n) ℂ) :
    A.PosDef ↔ ∃ L : Matrix (Fin n) (Fin n) ℂ, IsLower L ∧ (∀ i, 0 < L i i) ∧ A = L * Lᴴ :=
  ⟨cholesky_exists A, fun ⟨L, hL, hd, h⟩ => h ▸ posDef_of_factorization L hL hd⟩

end Herm
end Cholesky
