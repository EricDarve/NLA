import src.Counterexample.Sharpness

/-!
# Remaining facts stated in the optional section

* `meinguet_prod_le`: Meinguet's theorem replaces the product by the larger
  `(1 + f_{n-1})^{n-1}`.
* `spacing_facts`: `η` and `δ` are multiples of `4u`, the cross-block shift `η` is a multiple
  of `4u/√k`; `η√k = 12ku ≤ 3/512 < 1/4`; `48k < 2^{2q+6}` and `2q + 6 ≤ p`; `-η/8` and
  `9/(4k)` are multiples of `4u` (this uses `q ≥ 3`).
* `constant_57`: `27 (17/2 + 39/256)(1/256 + 17/72) = 3674685/65536 < 57`.
* `norm2_T_sq_ge`: `‖T‖₂² ≥ 9/2`.
* `cStar_admissible`: at an admissible size `n = 9k`, `c* < 57 / n^{3/2}`.
-/
namespace Cholesky
open Matrix

/-- Meinguet's theorem uses `(1 + f_{n-1})^{n-1}`, which is larger than the product. -/
theorem meinguet_prod_le (n : ℕ) (u : ℝ) (hu : 0 ≤ u) (hu5 : u < 1 / 5) :
    (∏ j ∈ Finset.range (n - 1), (1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 : ℕ) + 2))) ≤
      (1 + u / (1 - 5 * u) * (Real.sqrt (n - 1 : ℕ) + 2)) ^ (n - 1) := by
  have hc : 0 ≤ u / (1 - 5 * u) := div_nonneg hu (by linarith)
  calc (∏ j ∈ Finset.range (n - 1), (1 + u / (1 - 5 * u) * (Real.sqrt (j + 1 : ℕ) + 2)))
      ≤ ∏ _j ∈ Finset.range (n - 1), (1 + u / (1 - 5 * u) * (Real.sqrt (n - 1 : ℕ) + 2)) := by
        apply Finset.prod_le_prod
        · intro j _; have := Real.sqrt_nonneg ((j + 1 : ℕ) : ℝ); positivity
        · intro j hj
          rw [Finset.mem_range] at hj
          have : Real.sqrt ((j + 1 : ℕ) : ℝ) ≤ Real.sqrt ((n - 1 : ℕ) : ℝ) :=
            Real.sqrt_le_sqrt (by exact_mod_cast (by omega : j + 1 ≤ n - 1))
          nlinarith
    _ = (1 + u / (1 - 5 * u) * (Real.sqrt (n - 1 : ℕ) + 2)) ^ (n - 1) := by
        rw [Finset.prod_const, Finset.card_range]

/-- `27 (17/2 + 39/256)(1/256 + 17/72) = 3674685/65536 < 57`. -/
theorem constant_57 :
    (27 : ℝ) * (17 / 2 + 39 * (1 / 256)) * (1 / 256 + 17 / 72) = 3674685 / 65536 ∧
      (3674685 : ℝ) / 65536 < 57 := by
  constructor <;> norm_num

namespace CE
namespace Par
variable (P : Par)

/-- The spacing facts of the construction. -/
theorem spacing_facts :
    P.η = (3 * P.K) * (4 * P.u) ∧ P.δ = (3 * P.K ^ 3) * (4 * P.u) ∧
      P.η = (3 * P.k) * (4 * P.u / P.K) ∧
      P.η * P.K = 12 * P.k * P.u ∧ P.η * P.K ≤ 3 / 512 ∧ (3 / 512 : ℝ) < 1 / 4 ∧
      48 * P.k < 2 ^ (2 * P.q + 6) ∧ 2 * P.q + 6 ≤ P.p ∧
      P.η / 8 = (3 * P.L) * (4 * P.u) ∧
      9 / (4 * (P.k : ℝ)) = (9 * 2 ^ (P.p - 2 * P.q - 4 : ℕ)) * (4 * P.u) := by
  have hK := P.K_pos
  have hku := P.ku_le
  have hk := P.k_real
  have hKL := P.K_eq
  refine ⟨?_, ?_, ?_, ?_, ?_, by norm_num, ?_, ?_, ?_, ?_⟩
  · rw [P.η_eq_u]; ring
  · rw [δ, P.η_eq_u, hk]; ring
  · rw [P.η_eq_u, hk]; field_simp; ring
  · rw [P.η_eq_u, hk]; ring
  · rw [P.η_eq_u, show 12 * P.K * P.u * P.K = 12 * (P.K ^ 2 * P.u) by ring, ← hk]
    linarith
  · calc 48 * P.k < 64 * P.k := by have := P.k_ge; omega
      _ = 2 ^ (2 * P.q + 6) := by rw [pow_add]; simp only [k]; ring
  · simp only [p, q]; omega
  · rw [P.η_eq_u, hKL]; ring
  · have h1 : (2 : ℝ) ^ (P.p - 2 * P.q - 4 : ℕ) * (2 : ℝ) ^ (2 * P.q + 4) = 2 ^ P.p := by
      rw [← pow_add]; congr 1; simp only [p, q]; omega
    have hu : P.u * (2 : ℝ) ^ P.p = 1 := by
      rw [u, unitRoundoff, ← zpow_natCast, ← zpow_add₀ (by norm_num)]; simp
    rw [hk, show P.K ^ 2 = (2 : ℝ) ^ (2 * P.q) by rw [K, ← pow_mul, mul_comm]]
    rw [div_eq_iff (by positivity)]
    have : (2 : ℝ) ^ (2 * P.q + 4) = 16 * 2 ^ (2 * P.q) := by rw [pow_add]; norm_num; ring
    rw [show 9 * (2 : ℝ) ^ (P.p - 2 * P.q - 4 : ℕ) * (4 * P.u) * (4 * 2 ^ (2 * P.q)) =
      9 * ((2 : ℝ) ^ (P.p - 2 * P.q - 4 : ℕ) * 2 ^ (2 * P.q + 4)) * P.u by rw [this]; ring, h1]
    linear_combination (-9 : ℝ) * hu

/-- `‖T‖₂² ≥ 9/2`. -/
theorem norm2_T_sq_ge : 9 / 2 ≤ norm2 P.T ^ 2 := by
  obtain ⟨e0, he0def⟩ : ∃ e0 : Fin P.k → ℝ, e0 = Pi.single P.i0 1 := ⟨_, rfl⟩
  obtain ⟨y0, hy0def⟩ : ∃ y0 : Fin (P.k + P.k) → ℝ, y0 = Fin.append e0 (P.Qᵀ *ᵥ e0) := ⟨_, rfl⟩
  have he0 : vnorm e0 = 1 := by rw [he0def]; exact vnorm_single _
  have hQe0 : vnorm (P.Qᵀ *ᵥ e0) = 1 := by rw [P.vnorm_QT, he0]
  have hy0 : vnorm y0 ^ 2 = 2 := by rw [hy0def, vnorm_append_sq, he0, hQe0]; norm_num
  have hG : y0 ⬝ᵥ (P.G *ᵥ y0) = 9 := by
    rw [P.G_quad, hy0def, vL_append, vR_append, he0, hQe0, Matrix.mulVec_mulVec,
      P.Q_mul_transpose, Matrix.one_mulVec, ← vnorm_sq, he0]
    norm_num
  have hT := P.T_quad y0
  have hη := P.η_pos
  have hpos : 0 ≤ 9 / 8 * P.η * (∑ i, y0 i) ^ 2 := mul_nonneg (by linarith) (sq_nonneg _)
  have h9 : 9 ≤ vnorm (P.T *ᵥ y0) ^ 2 := by linarith
  have hle := vnorm_mulVec_le P.T y0
  have hsq := mul_self_le_mul_self (vnorm_nonneg _) hle
  have hy0' : vnorm y0 * vnorm y0 = 2 := by rw [← sq]; exact hy0
  nlinarith [norm2_nonneg P.T]

end Par

/-- **At an admissible size** `n = 9k`: `c*_{n,u} < 57 / n^{3/2}`. -/
theorem cStar_admissible (p q : ℕ) (h : Admissible p q) {r : ℝ → ℝ}
    (hr : IsRoundNearest p r) :
    cStar p ((6 * (Par.ofPQ p q).k + (Par.ofPQ p q).k) + ((Par.ofPQ p q).k + (Par.ofPQ p q).k))
        (Ops.rounded r) <
      57 / (((9 * (Par.ofPQ p q).k : ℕ) : ℝ) * Real.sqrt ((9 * (Par.ofPQ p q).k : ℕ) : ℝ)) := by
  set P := Par.ofPQ p q
  have hPp : P.p = p := Par.ofPQ_p h
  have hfl : ∀ i j, IsFloat p (P.A i j) := fun i j => by rw [← hPp]; exact P.A_float i j
  rw [← hPp] at hr
  have hle := cStar_le P.A hfl P.A_posDef (P.A_fails hr).1
  have hlt := P.cond_A_lt
  have hu : P.u = unitRoundoff p := by rw [Par.u, hPp]
  rw [hu] at hlt
  have hupos := unitRoundoff_pos p
  have hpos : 0 < ((9 * P.k : ℕ) : ℝ) * Real.sqrt ((9 * P.k : ℕ) : ℝ) := by
    have : (0 : ℝ) < ((9 * P.k : ℕ) : ℝ) := by
      have := P.k_ge; exact_mod_cast (by omega : 0 < 9 * P.k)
    have := Real.sqrt_pos.mpr this
    positivity
  calc cStar p _ (Ops.rounded r) ≤ unitRoundoff p * cond2 P.A := hle
    _ < unitRoundoff p * (57 / (((9 * P.k : ℕ) : ℝ) * Real.sqrt ((9 * P.k : ℕ) : ℝ) *
          unitRoundoff p)) := mul_lt_mul_of_pos_left hlt hupos
    _ = 57 / (((9 * P.k : ℕ) : ℝ) * Real.sqrt ((9 * P.k : ℕ) : ℝ)) := by
        field_simp

end CE
end Cholesky
