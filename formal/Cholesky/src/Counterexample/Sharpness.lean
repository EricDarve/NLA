import src.Counterexample.Padding
import src.Parameters
import src.Completion

/-!
# The smallest condition number permitting failure

For a fixed arithmetic `o` (rounding to nearest at precision `p`),

`c*_{n,u} = inf { u κ₂(A) : A a stored SPD n × n matrix on which the algorithm fails }`

(`cStar`). We prove:

* `succeeds_of_cond_lt`: `κ₂(A) < c*/u` guarantees success;
* `cStar_ge_u`: `c* ≥ u`, because `κ₂(A) ≥ 1`;
* `cStar_ge_meinguet`, `cStar_ge_meinguet_mu`: *assuming* Meinguet's completion theorem
  (`MeinguetGuarantee`, an explicit hypothesis, not an axiom),
  `c* ≥ u / μ_n(u)` and `c* ≥ 1 / (1.1 n^{3/2})`;
* `fixed_precision`: at a fixed precision `p ≥ 20`, the largest admissible block has
  `x* ∈ {2⁻⁸, 2⁻¹⁰, 2⁻¹²}` and `κ₂ < 8239`; `binary64`, `binary32`: the special cases;
* `sharp_family`: for every `N ≥ 2304` a stored SPD `N × N` matrix on which Cholesky fails,
  with `κ₂ < 8239 max{1, 1/(N^{3/2} u)}`;
* `cStar_lt`: `c*_{N,u} < 8239 max{u, 1/N^{3/2}}`.
-/
namespace Cholesky
open Matrix

/-! ## Reindexing by a cast of the size -/

/-- A matrix of size `n` viewed as a matrix of size `n'` when `n = n'`. -/
def castMat {n n' : ℕ} (h : n = n') (A : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n') (Fin n') ℝ :=
  A.submatrix (Fin.cast h.symm) (Fin.cast h.symm)

lemma castMat_rfl {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) : castMat rfl A = A := by
  ext i j; rfl

lemma castMat_props {n n' : ℕ} (h : n = n') (A : Matrix (Fin n) (Fin n) ℝ) (p : ℕ) (o : Ops) :
    ((∀ i j, IsFloat p (castMat h A i j)) ↔ ∀ i j, IsFloat p (A i j)) ∧
      ((castMat h A).PosDef ↔ A.PosDef) ∧ (Succeeds o (castMat h A) ↔ Succeeds o A) ∧
      cond2 (castMat h A) = cond2 A := by
  subst h
  rw [castMat_rfl]
  exact ⟨Iff.rfl, Iff.rfl, Iff.rfl, rfl⟩

/-! ## The threshold `c*` -/

/-- The values `u κ₂(A)` over stored SPD matrices on which `o` fails. -/
def failSet (p n : ℕ) (o : Ops) : Set ℝ :=
  {c | ∃ A : Matrix (Fin n) (Fin n) ℝ, (∀ i j, IsFloat p (A i j)) ∧ A.PosDef ∧
    ¬ Succeeds o A ∧ c = unitRoundoff p * cond2 A}

/-- `c*_{n,u}`. -/
noncomputable def cStar (p n : ℕ) (o : Ops) : ℝ := sInf (failSet p n o)

variable {p n : ℕ} {o : Ops}

lemma failSet_bdd : BddBelow (failSet p n o) :=
  ⟨0, by rintro c ⟨A, -, -, -, rfl⟩; exact mul_nonneg (unitRoundoff_pos p).le (cond2_nonneg A)⟩

lemma cStar_le (A : Matrix (Fin n) (Fin n) ℝ) (hfl : ∀ i j, IsFloat p (A i j)) (hpd : A.PosDef)
    (hfail : ¬ Succeeds o A) : cStar p n o ≤ unitRoundoff p * cond2 A :=
  csInf_le failSet_bdd ⟨A, hfl, hpd, hfail, rfl⟩

/-- `κ₂(A) < c*/u` guarantees success. -/
theorem succeeds_of_cond_lt (A : Matrix (Fin n) (Fin n) ℝ) (hfl : ∀ i j, IsFloat p (A i j))
    (hpd : A.PosDef) (h : cond2 A < cStar p n o / unitRoundoff p) : Succeeds o A := by
  by_contra hf
  have h1 := cStar_le A hfl hpd hf
  have hu := unitRoundoff_pos p
  rw [lt_div_iff₀ hu] at h
  linarith

lemma norm2_one_eq (hn : 1 ≤ n) : norm2 (1 : Matrix (Fin n) (Fin n) ℝ) = 1 := by
  apply le_antisymm norm2_one_le
  have i0 : Fin n := ⟨0, by omega⟩
  apply le_norm2_of (v := Pi.single i0 1) (by simp)
  rw [Matrix.one_mulVec, one_mul]

/-- `κ₂(A) ≥ 1`. -/
lemma cond2_ge_one (hn : 1 ≤ n) (A : Matrix (Fin n) (Fin n) ℝ) (hA : IsUnit A.det) :
    1 ≤ cond2 A := by
  have := norm2_mul_le A A⁻¹
  rwa [Matrix.mul_nonsing_inv _ hA, norm2_one_eq hn] at this

theorem cStar_ge_u (hn : 1 ≤ n) (hne : (failSet p n o).Nonempty) : unitRoundoff p ≤ cStar p n o := by
  apply le_csInf hne
  rintro c ⟨A, -, hpd, -, rfl⟩
  have hdet : IsUnit A.det := (Matrix.isUnit_iff_isUnit_det _).mp hpd.isUnit
  exact le_mul_of_one_le_right (unitRoundoff_pos p).le (cond2_ge_one hn A hdet)

/-- **Meinguet's completion theorem**, as a hypothesis: `κ₂(A) μ_n(u) < 1` guarantees success.
It is cited in the notes and not proved here. -/
def MeinguetGuarantee (p n : ℕ) (o : Ops) : Prop :=
  ∀ A : Matrix (Fin n) (Fin n) ℝ, (∀ i j, IsFloat p (A i j)) → A.PosDef →
    cond2 A * meinguetMu n (unitRoundoff p) < 1 → Succeeds o A

lemma one_le_prod_real (s : Finset ℕ) (f : ℕ → ℝ) (h : ∀ j ∈ s, 1 ≤ f j) :
    1 ≤ ∏ j ∈ s, f j := by
  have := Finset.prod_le_prod (s := s) (f := fun _ => (1 : ℝ)) (g := f)
    (fun _ _ => zero_le_one) h
  simpa using this

lemma meinguetMu_pos (hn : 2 ≤ n) {u : ℝ} (hu : 0 < u) (hu5 : u < 1 / 5) : 0 < meinguetMu n u := by
  unfold meinguetMu
  have hd : 0 < 1 - 5 * u := by linarith
  have hf : ∀ j ∈ Finset.range (n - 1), (1 : ℝ) ≤ 1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 : ℕ) + 2) :=
    fun j _ => by have := Real.sqrt_nonneg ((j + 1 : ℕ) : ℝ); have := div_pos hu hd; nlinarith
  have h0 : (1 : ℝ) < 1 + u / (1 - 5 * u) * (Real.sqrt ((0 + 1 : ℕ) : ℝ) + 2) := by
    have := Real.sqrt_nonneg ((0 + 1 : ℕ) : ℝ); have := div_pos hu hd; nlinarith
  have hprod : (1 : ℝ) < ∏ j ∈ Finset.range (n - 1),
      (1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 : ℕ) + 2)) := by
    rw [show n - 1 = (n - 2) + 1 by omega, Finset.prod_range_succ']
    calc (1 : ℝ) = 1 * 1 := by ring
      _ < (∏ j ∈ Finset.range (n - 2), (1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 + 1 : ℕ) + 2))) *
          (1 + u / (1 - 5 * u) * (Real.sqrt ((0 + 1 : ℕ) : ℝ) + 2)) := by
        have hrest := one_le_prod_real (Finset.range (n - 2))
          (fun j => 1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 + 1 : ℕ) + 2)) (fun j _ => by
            have := Real.sqrt_nonneg ((j + 1 + 1 : ℕ) : ℝ); have := div_pos hu hd; nlinarith)
        nlinarith
  linarith

/-- Assuming Meinguet's theorem, `c* ≥ u / μ_n(u)`. -/
theorem cStar_ge_meinguet_mu (hM : MeinguetGuarantee p n o) (hn : 2 ≤ n)
    (hu5 : unitRoundoff p < 1 / 5) (hne : (failSet p n o).Nonempty) :
    unitRoundoff p / meinguetMu n (unitRoundoff p) ≤ cStar p n o := by
  have hμ := meinguetMu_pos hn (unitRoundoff_pos p) hu5
  apply le_csInf hne
  rintro c ⟨A, hfl, hpd, hfail, rfl⟩
  have : 1 ≤ cond2 A * meinguetMu n (unitRoundoff p) := by
    by_contra h; exact hfail (hM A hfl hpd (lt_of_not_ge h))
  rw [div_le_iff₀ hμ]
  have hu := unitRoundoff_pos p
  nlinarith

/-- Assuming Meinguet's theorem, `c* ≥ 1 / (1.1 n^{3/2})` for `n ≥ 576`, `u ≤ 2^{-17}`. -/
theorem cStar_ge_meinguet (hM : MeinguetGuarantee p n o) (hn : 576 ≤ n)
    (hu : unitRoundoff p ≤ 1 / 131072) (hne : (failSet p n o).Nonempty) :
    1 / (11 / 10 * (n * Real.sqrt n)) ≤ cStar p n o := by
  apply le_csInf hne
  rintro c ⟨A, hfl, hpd, hfail, rfl⟩
  have hdet : IsUnit A.det := (Matrix.isUnit_iff_isUnit_det _).mp hpd.isUnit
  have hκ := cond2_ge_one (by omega) A hdet
  have hpos : 0 < 11 / 10 * ((n : ℝ) * Real.sqrt n) := by
    have : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
    positivity
  have hgt : 1 < 11 / 10 * ((n : ℝ) * Real.sqrt n) * unitRoundoff p * cond2 A := by
    by_contra h
    push_neg at h
    exact hfail (hM A hfl hpd (by
      have := meinguet_condition n (unitRoundoff p) (cond2 A) hn (unitRoundoff_pos p) hu hκ h
      linarith))
  rw [div_le_iff₀ hpos]
  linarith

/-! ## Every stored entry is a float -/

namespace CE
namespace Par
variable (P : Par)

lemma η_float : IsFloat P.p P.η := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  rw [P.η_eq, P.K_cast]
  have := P.isFloat_E0 (48 * (P.Kn : ℤ) ^ 2) (by
    rw [P.two_pow_p_nat, abs_of_nonneg (by positivity)]
    have h1 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := le_mul_of_one_le_right (by positivity) hD
    nlinarith [mul_le_mul_of_nonneg_right hK (by positivity : (0 : ℤ) ≤ (P.Kn : ℤ) ^ 2)])
  push_cast at this; exact this

lemma cross0_float (ε : ℤ) (hε : ε = 1 ∨ ε = -1) : IsFloat P.p ((9 / 4 : ℝ) * (ε / P.K) + P.η) := by
  obtain ⟨hK, hD, hKL, hL⟩ := kbounds (P := P)
  have h := P.cross_eq ε 0
  simp only [Nat.cast_zero, mul_zero, sub_zero] at h
  rw [h, show P.Ncross ε 0 = (ε * (144 * (P.Kn : ℤ) ^ 3 * P.Dn) + 3 * (P.Kn : ℤ) ^ 2) * 2 ^ 4 by
    simp only [Ncross]; push_cast; ring]
  apply P.isFloat_E0_grid
  rw [P.two_pow_p_nat, abs_le]
  have h3 : ((P.Kn : ℤ)) ^ 3 ≤ (P.Kn : ℤ) ^ 3 * P.Dn := le_mul_of_one_le_right (by positivity) hD
  have h32 : 8 * (P.Kn : ℤ) ^ 2 ≤ (P.Kn : ℤ) ^ 3 := by nlinarith [sq_nonneg (P.Kn : ℤ)]
  rcases hε with rfl | rfl <;> constructor <;> nlinarith

lemma Q_entry_cross (i j : Fin P.k) : IsFloat P.p (9 / 4 * P.Q i j + P.η) := by
  rw [P.Q_eq]; exact P.cross0_float _ (P.hInt_cases i j)

lemma one_lr (i j : Fin P.k) : (1 : Matrix (Fin (P.k + P.k)) (Fin (P.k + P.k)) ℝ)
    (Fin.castAdd P.k i) (Fin.natAdd P.k j) = 0 := by
  rw [one_eq_blk, blk_lr]; rfl

lemma C_lr (i j : Fin P.k) : P.C (Fin.castAdd P.k i) (Fin.natAdd P.k j) = 9 / 4 * P.Q i j + P.η := by
  simp only [C, G, J2, ones, Matrix.add_apply, Matrix.smul_apply, blk_lr, Matrix.of_apply,
    P.one_lr, smul_eq_mul]
  ring

lemma C_float (i j : Fin (P.k + P.k)) : IsFloat P.p (P.C i j) := by
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j
  · rw [C_ll]; split_ifs
    · exact P.dg_float
    · exact P.η_float
  · rw [C_lr]; exact P.Q_entry_cross i j
  · rw [C_rl]; exact P.Q_entry_cross j i
  · rw [C_rr]; split_ifs
    · exact P.dg_float
    · exact P.η_float

/-- **All entries of `A_{9k}` are exactly representable.** -/
theorem A_float (i j : Fin ((6 * P.k + P.k) + (P.k + P.k))) : IsFloat P.p (P.A i j) := by
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j
  · simp only [A, blk_ll, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
    split_ifs
    · simpa using P.isFloat_small_int 4 (by norm_num)
    · simpa using isFloat_zero P.p
  · simp only [A, blk_lr, Matrix.smul_apply, smul_eq_mul]; exact (P.T_float i j).2
  · simp only [A, blk_rl, Matrix.smul_apply, Matrix.transpose_apply, smul_eq_mul]
    exact (P.T_float j i).2
  · simp only [A, blk_rr]; exact P.C_float i j

theorem B_float (m : ℕ) (i j : Fin ((6 * P.k + P.k) + (P.k + P.k) + m)) :
    IsFloat P.p (P.B m i j) := by
  refine Fin.addCases (fun i => ?_) (fun i => ?_) i <;>
    refine Fin.addCases (fun j => ?_) (fun j => ?_) j
  · simp only [B, blk_ll]; exact P.A_float i j
  · simp only [B, blk_lr, Matrix.zero_apply]; exact isFloat_zero _
  · simp only [B, blk_rl, Matrix.zero_apply]; exact isFloat_zero _
  · simp only [B, blk_rr, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
    split_ifs
    · simpa using P.isFloat_small_int 4 (by norm_num)
    · simpa using isFloat_zero P.p

/-- `x = k^{3/2} u = 2^{3q - p}`. -/
lemma x_zpow : P.x = (2 : ℝ) ^ ((3 * P.q : ℕ) - (P.p : ℤ)) := by
  rw [x, u, unitRoundoff, K, ← pow_mul, ← zpow_natCast, ← zpow_add₀ (by norm_num)]
  congr 1; push_cast; ring

lemma nk : (6 * P.k + P.k) + (P.k + P.k) = 9 * P.k := by ring

end Par

/-- The parameters for an admissible pair `(p, q)`. -/
def Par.ofPQ (p q : ℕ) : Par := ⟨q - 3, (p - 3 * q - 8) / 2⟩

lemma Par.ofPQ_q {p q : ℕ} (h : Admissible p q) : (Par.ofPQ p q).q = q := by
  obtain ⟨hq, -, -⟩ := h; simp only [Par.ofPQ, Par.q]; omega

lemma Par.ofPQ_p {p q : ℕ} (h : Admissible p q) : (Par.ofPQ p q).p = p := by
  obtain ⟨hq, hpar, hp⟩ := h; simp only [Par.ofPQ, Par.p, Par.q]; omega

lemma Par.ofPQ_k {p q : ℕ} (h : Admissible p q) : (Par.ofPQ p q).k = 4 ^ q := by
  simp only [Par.k, Par.ofPQ_q h, pow_mul]; norm_num

/-! ## Fixed precision -/

/-- **The largest admissible block at a fixed precision `p ≥ 20`.** Its parameter is
`x* = 2^{-r}` with `r ∈ {8, 10, 12}`, so `κ₂ ≤ 303691615/36864 < 8239`. -/
theorem fixed_precision (p : ℕ) (hp : 20 ≤ p) :
    let P := Par.ofPQ p (largestQ p)
    (P.x = 1 / 256 ∨ P.x = 1 / 1024 ∨ P.x = 1 / 4096) ∧
      cond2 P.A ≤ 303691615 / 36864 ∧ cond2 P.A < 8239 := by
  intro P
  obtain ⟨h3, hsize, hpar, hr⟩ := largestQ_properties p hp
  have hadm : Admissible p (largestQ p) := ⟨h3, hpar, hsize⟩
  have hx : P.x = 1 / 256 ∨ P.x = 1 / 1024 ∨ P.x = 1 / 4096 := by
    rw [P.x_zpow, Par.ofPQ_q hadm, Par.ofPQ_p hadm]
    have e : ((3 * largestQ p : ℕ) : ℤ) - (p : ℤ) = -((p - 3 * largestQ p : ℕ) : ℤ) := by
      push_cast [Nat.cast_sub (by omega : 3 * largestQ p ≤ p)]; ring
    rw [e]
    rcases hr with h | h | h <;> rw [h] <;> norm_num
  have hb := conditionUpper_fixed_precision P.x hx
  exact ⟨hx, le_trans P.cond_A_le hb.1, lt_of_le_of_lt P.cond_A_le hb.2⟩

/-- **binary64**: `p = 53`, `q* = 15`, `x* = 2⁻⁸`, block size `9 · 2³⁰`, `κ₂ < 532`. -/
theorem binary64 :
    largestQ 53 = 15 ∧ (Par.ofPQ 53 15).k = 2 ^ 30 ∧ (Par.ofPQ 53 15).x = 1 / 256 ∧
      cond2 (Par.ofPQ 53 15).A ≤ 1224895 / 2304 ∧ cond2 (Par.ofPQ 53 15).A < 532 := by
  have hadm : Admissible 53 15 := ⟨by norm_num, by norm_num, by norm_num⟩
  have hx : (Par.ofPQ 53 15).x = 1 / 256 := by
    rw [Par.x_zpow, Par.ofPQ_q hadm, Par.ofPQ_p hadm]; norm_num
  have hb := conditionUpper_binary64
  refine ⟨by decide, by rw [Par.ofPQ_k hadm]; norm_num, hx, ?_, ?_⟩
  · rw [← hb.1, ← hx]; exact Par.cond_A_le _
  · rw [← hx] at hb; exact lt_of_le_of_lt (Par.cond_A_le _) hb.2

/-- **binary32**: `p = 24`; the largest admissible block has size `2304`, `k^{3/2} u = 1/4096`,
and `κ₂ < 8239`. -/
theorem binary32 :
    largestQ 24 = 4 ∧ 9 * (Par.ofPQ 24 4).k = 2304 ∧ (Par.ofPQ 24 4).x = 1 / 4096 ∧
      cond2 (Par.ofPQ 24 4).A < 8239 := by
  have hadm : Admissible 24 4 := ⟨by norm_num, by norm_num, by norm_num⟩
  have hx : (Par.ofPQ 24 4).x = 1 / 4096 := by
    rw [Par.x_zpow, Par.ofPQ_q hadm, Par.ofPQ_p hadm]; norm_num
  refine ⟨by decide, by rw [Par.ofPQ_k hadm]; norm_num, hx, ?_⟩
  have hb := conditionUpper_fixed_precision _ (Or.inr (Or.inr hx))
  exact lt_of_le_of_lt (Par.cond_A_le _) hb.2

/-! ## Every dimension `N ≥ 2304` -/

lemma pad_ineq (n0 N : ℕ) (hn0 : 0 < n0) (hN0 : 0 < N) (hN : N < 16 * n0) {u : ℝ} (hu : 0 < u) :
    57 / ((n0 : ℝ) * Real.sqrt n0 * u) < 3648 / ((N : ℝ) * Real.sqrt N * u) := by
  have hn0' : (0 : ℝ) < n0 := by exact_mod_cast hn0
  have hN0' : (0 : ℝ) < N := by exact_mod_cast hN0
  have hN' : (N : ℝ) < 16 * n0 := by exact_mod_cast hN
  have hs0 := Real.sqrt_pos.mpr hn0'
  have hsN := Real.sqrt_pos.mpr hN0'
  have hsq : Real.sqrt N < 4 * Real.sqrt n0 := by
    rw [show (4 : ℝ) * Real.sqrt n0 = Real.sqrt (16 * n0) by
      rw [Real.sqrt_mul (by norm_num), show Real.sqrt 16 = 4 by
        rw [show (16 : ℝ) = 4 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]]
    exact Real.sqrt_lt_sqrt hN0'.le hN'
  have key : (N : ℝ) * Real.sqrt N < 64 * ((n0 : ℝ) * Real.sqrt n0) := by
    have := mul_lt_mul hN' hsq.le hsN (by positivity)
    nlinarith
  rw [div_lt_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_pos hu (mul_pos hN0' hsN)]

/-- **A failing family in every dimension `N ≥ 2304`** at a fixed precision `p ≥ 20`, with
`κ₂ < 8239 max{1, 1/(N^{3/2} u)}`. -/
theorem sharp_family (p N : ℕ) (hp : 20 ≤ p) (hN : 2304 ≤ N) :
    ∃ B : Matrix (Fin N) (Fin N) ℝ, (∀ i j, IsFloat p (B i j)) ∧ B.PosDef ∧
      (∀ r, IsRoundNearest p r → ¬ Succeeds (Ops.rounded r) B ∧ ¬ Succeeds (Ops.fused r) B) ∧
      cond2 B < 8239 * max 1 (1 / ((N : ℝ) * Real.sqrt N * unitRoundoff p)) := by
  classical
  set pred : ℕ → Prop := fun q => 3 ≤ q ∧ (p + q) % 2 = 0 ∧ 3 * q + 8 ≤ p ∧ 9 * 4 ^ q ≤ N
  obtain ⟨h3s, hsizes, hpars, -⟩ := largestQ_properties p hp
  have hq0 : pred (if p % 2 = 1 then 3 else 4) := by
    split_ifs with h
    · refine ⟨le_refl _, by omega, by omega, by norm_num; omega⟩
    · refine ⟨by norm_num, by omega, by omega, by norm_num; omega⟩
  have hq0le : (if p % 2 = 1 then 3 else 4) ≤ largestQ p :=
    largestQ_maximal p _ hq0.2.2.1 hq0.2.1
  set q := Nat.findGreatest pred (largestQ p)
  have hq : pred q := Nat.findGreatest_spec hq0le hq0
  have hqle : q ≤ largestQ p := Nat.findGreatest_le _
  have hadm : Admissible p q := ⟨hq.1, hq.2.1, hq.2.2.1⟩
  set P := Par.ofPQ p q
  have hPp : P.p = p := Par.ofPQ_p hadm
  have hPk : P.k = 4 ^ q := Par.ofPQ_k hadm
  have hn0 : (6 * P.k + P.k) + (P.k + P.k) ≤ N := by rw [P.nk, hPk]; exact hq.2.2.2
  set m := N - ((6 * P.k + P.k) + (P.k + P.k))
  have hNm : (6 * P.k + P.k) + (P.k + P.k) + m = N := by omega
  obtain ⟨hfl, hpd, hsucc, hcond⟩ := castMat_props hNm (P.B m) p (Ops.rounded (fun x => x))
  refine ⟨castMat hNm (P.B m), ?_, ?_, ?_, ?_⟩
  · rw [hfl]; intro i j; rw [← hPp]; exact P.B_float m i j
  · rw [hpd]; exact P.B_posDef m
  · intro r hr
    rw [← hPp] at hr
    obtain ⟨h1, h2⟩ := P.B_fails m hr
    exact ⟨fun h => h1 ((castMat_props hNm (P.B m) p _).2.2.1.mp h),
      fun h => h2 ((castMat_props hNm (P.B m) p _).2.2.1.mp h)⟩
  · rw [hcond, P.cond2_B]
    have hu : 0 < unitRoundoff p := unitRoundoff_pos p
    by_cases hqe : q = largestQ p
    · have hfix := (fixed_precision p hp).2.2
      have hPeq : P = Par.ofPQ p (largestQ p) := by rw [← hqe]
      rw [hPeq]
      exact lt_of_lt_of_le hfix (le_mul_of_one_le_right (by norm_num) (le_max_left _ _))
    · have hlt : q < largestQ p := lt_of_le_of_ne hqle hqe
      have hq2 : q + 2 ≤ largestQ p := by
        have := hq.2.1; omega
      have hnot : ¬ pred (q + 2) := Nat.findGreatest_is_greatest (by omega) hq2
      have hbig : N < 16 * (9 * 4 ^ q) := by
        by_contra h
        apply hnot
        refine ⟨by omega, by have := hq.2.1; omega, by omega, ?_⟩
        rw [pow_add]; norm_num at h ⊢; omega
      have h57 := P.cond_A_lt
      rw [show P.u = unitRoundoff p by rw [Par.u, hPp]] at h57
      have hk9 : 9 * P.k = 9 * 4 ^ q := by rw [hPk]
      rw [hk9] at h57
      have hpad := pad_ineq (9 * 4 ^ q) N (by positivity) (by omega) hbig hu
      have hfinal : (3648 : ℝ) / ((N : ℝ) * Real.sqrt N * unitRoundoff p) ≤
          8239 * max 1 (1 / ((N : ℝ) * Real.sqrt N * unitRoundoff p)) := by
        have hpos : 0 < (N : ℝ) * Real.sqrt N * unitRoundoff p := by
          have : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
          have := Real.sqrt_pos.mpr this
          positivity
        calc (3648 : ℝ) / ((N : ℝ) * Real.sqrt N * unitRoundoff p) =
            3648 * (1 / ((N : ℝ) * Real.sqrt N * unitRoundoff p)) := by ring
          _ ≤ 8239 * (1 / ((N : ℝ) * Real.sqrt N * unitRoundoff p)) :=
              mul_le_mul_of_nonneg_right (by norm_num) (by positivity)
          _ ≤ 8239 * max 1 (1 / ((N : ℝ) * Real.sqrt N * unitRoundoff p)) :=
              mul_le_mul_of_nonneg_left (le_max_right _ _) (by norm_num)
      push_cast at h57 hpad
      exact lt_of_lt_of_le (lt_trans h57 hpad) hfinal

/-- **Upper bound on `c*`**: `c*_{N,u} < 8239 max{u, 1/N^{3/2}}` for `p ≥ 20`, `N ≥ 2304`. -/
theorem cStar_lt (p N : ℕ) (hp : 20 ≤ p) (hN : 2304 ≤ N) {r : ℝ → ℝ}
    (hr : IsRoundNearest p r) :
    cStar p N (Ops.rounded r) < 8239 * max (unitRoundoff p) (1 / ((N : ℝ) * Real.sqrt N)) ∧
      cStar p N (Ops.fused r) < 8239 * max (unitRoundoff p) (1 / ((N : ℝ) * Real.sqrt N)) := by
  obtain ⟨B, hfl, hpd, hfail, hcond⟩ := sharp_family p N hp hN
  have hu := unitRoundoff_pos p
  have hNpos : (0 : ℝ) < (N : ℝ) * Real.sqrt N := by
    have : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
    have := Real.sqrt_pos.mpr this
    positivity
  have hmax : unitRoundoff p * (8239 * max 1 (1 / ((N : ℝ) * Real.sqrt N * unitRoundoff p))) =
      8239 * max (unitRoundoff p) (1 / ((N : ℝ) * Real.sqrt N)) := by
    rw [mul_left_comm, mul_max_of_nonneg _ _ hu.le, mul_one]
    congr 2; field_simp
  have key : unitRoundoff p * cond2 B < 8239 * max (unitRoundoff p) (1 / ((N : ℝ) * Real.sqrt N)) := by
    rw [← hmax]; exact mul_lt_mul_of_pos_left hcond hu
  exact ⟨lt_of_le_of_lt (cStar_le B hfl hpd (hfail r hr).1) key,
    lt_of_le_of_lt (cStar_le B hfl hpd (hfail r hr).2) key⟩

/-- **Both bounds on `c*`.** The lower bound `1/(1.1 N^{3/2})` uses Meinguet's theorem,
passed as the hypothesis `hM`. -/
theorem cStar_bounds (p N : ℕ) (hp : 20 ≤ p) (hN : 2304 ≤ N) {r : ℝ → ℝ}
    (hr : IsRoundNearest p r) (hM : MeinguetGuarantee p N (Ops.rounded r)) :
    max (unitRoundoff p) (1 / (11 / 10 * ((N : ℝ) * Real.sqrt N))) ≤ cStar p N (Ops.rounded r) ∧
      cStar p N (Ops.rounded r) < 8239 * max (unitRoundoff p) (1 / ((N : ℝ) * Real.sqrt N)) := by
  obtain ⟨B, hfl, hpd, hfail, -⟩ := sharp_family p N hp hN
  have hne : (failSet p N (Ops.rounded r)).Nonempty := ⟨_, B, hfl, hpd, (hfail r hr).1, rfl⟩
  have hu17 : unitRoundoff p ≤ 1 / 131072 := by
    unfold unitRoundoff
    rw [show (1 / 131072 : ℝ) = (2 : ℝ) ^ (-(17 : ℤ)) by norm_num]
    exact zpow_le_zpow_right₀ (by norm_num) (by omega)
  refine ⟨max_le (cStar_ge_u (by omega) hne) (cStar_ge_meinguet hM (by omega) hu17 hne),
    (cStar_lt p N hp hN hr).1⟩

end CE

end Cholesky
