import src.Condition

/-!
# Sensitivity

For a problem `x = f(d)` the local absolute sensitivity is

`S_f(d) = lim_{ε → 0⁺} sup_{0 < ‖δd‖ ≤ ε} ‖f(d + δd) - f(d)‖ / ‖δd‖`.

`quotSup f d ε` is the supremum and `HasSensitivity f d S` says that the limit is `S`.

* `hasSensitivity_of_hasFDerivAt`: **`S_f(d) = ‖Df(d)‖`** for `f` differentiable at `d`
  (normed spaces over `ℝ`, `d` in a nontrivial space);
* `hasSensitivity_of_hasDerivAt`: for a scalar function, `S_f(d) = |f'(d)|`;
* `forward_le_sens`: `‖f(d + δd) - f(d)‖ ≤ S_f(d) ‖δd‖ + o(‖δd‖)`;
* `relative_le_cond`: the relative condition number `κ_f(d) = S_f(d) ‖d‖/‖f(d)‖` bounds the
  amplification of relative errors, up to `o(‖δd‖)`;
* `sens_linear`, `relCond_linear_le`: for `b ↦ A⁻¹ b` the sensitivity is `‖A⁻¹‖₂` and the
  relative condition number is at most `κ₂(A)`.
-/
open Filter Topology Asymptotics
open scoped Matrix.Norms.L2Operator

namespace LU

section General
variable {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [NormedAddCommGroup F]
  [NormedSpace ℝ F]

/-- `sup_{0 < ‖δd‖ ≤ ε} ‖f(d + δd) - f(d)‖ / ‖δd‖`. -/
noncomputable def quotSup (f : E → F) (d : E) (ε : ℝ) : ℝ :=
  sSup ((fun h => ‖f (d + h) - f d‖ / ‖h‖) '' {h : E | 0 < ‖h‖ ∧ ‖h‖ ≤ ε})

/-- `S_f(d) = S`: the suprema converge to `S` as `ε → 0⁺`. -/
def HasSensitivity (f : E → F) (d : E) (S : ℝ) : Prop :=
  Tendsto (quotSup f d) (𝓝[>] 0) (𝓝 S)

/-- **`S_f(d) = ‖Df(d)‖`.** -/
theorem hasSensitivity_of_hasFDerivAt [Nontrivial E] {f : E → F} {f' : E →L[ℝ] F} {d : E}
    (hf : HasFDerivAt f f' d) : HasSensitivity f d ‖f'‖ := by
  rw [HasSensitivity, Metric.tendsto_nhdsWithin_nhds]
  intro ε' hε'
  have hlo := hasFDerivAt_iff_isLittleO_nhds_zero.mp hf
  have hη : (0 : ℝ) < ε' / 3 := by positivity
  obtain ⟨δ, hδ, hball⟩ := Metric.eventually_nhds_iff.mp (hlo.def hη)
  refine ⟨δ, hδ, fun ε hε hεδ => ?_⟩
  simp only [Set.mem_Ioi] at hε
  rw [dist_zero_right, Real.norm_eq_abs, abs_of_pos hε] at hεδ
  set S := (fun h => ‖f (d + h) - f d‖ / ‖h‖) '' {h : E | 0 < ‖h‖ ∧ ‖h‖ ≤ ε}
  have key : ∀ h : E, 0 < ‖h‖ → ‖h‖ ≤ ε →
      |‖f (d + h) - f d‖ / ‖h‖ - ‖f' h‖ / ‖h‖| ≤ ε' / 3 := by
    intro h hh hhε
    have hb : ‖f (d + h) - f d - f' h‖ ≤ ε' / 3 * ‖h‖ := by
      have := hball (y := h) (by rw [dist_zero_right]; linarith)
      simpa [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg h)] using this
    rw [← sub_div, abs_div, abs_of_pos hh, div_le_iff₀ hh]
    calc |‖f (d + h) - f d‖ - ‖f' h‖| ≤ ‖f (d + h) - f d - f' h‖ := abs_norm_sub_norm_le _ _
      _ ≤ ε' / 3 * ‖h‖ := hb
  obtain ⟨h₀, hh₀⟩ := exists_ne (0 : E)
  have hh₀pos : 0 < ‖h₀‖ := norm_pos_iff.mpr hh₀
  have hscale : ∀ x : E, 0 < ‖x‖ → ‖(ε / ‖x‖) • x‖ = ε := fun x hx => by
    rw [norm_smul, Real.norm_eq_abs, abs_of_pos (by positivity)]; field_simp
  have hne : S.Nonempty :=
    ⟨_, (ε / ‖h₀‖) • h₀, ⟨by rw [hscale h₀ hh₀pos]; exact hε, by rw [hscale h₀ hh₀pos]⟩, rfl⟩
  have hup : ∀ s ∈ S, s ≤ ‖f'‖ + ε' / 3 := by
    rintro s ⟨h, ⟨hh, hhε⟩, rfl⟩
    have h1 := (abs_le.mp (key h hh hhε)).2
    have h2 : ‖f' h‖ / ‖h‖ ≤ ‖f'‖ := by rw [div_le_iff₀ hh]; exact f'.le_opNorm h
    linarith
  have hbdd : BddAbove S := ⟨_, hup⟩
  have hsup : sSup S ≤ ‖f'‖ + ε' / 3 := csSup_le hne hup
  have hlow : ‖f'‖ - 2 * (ε' / 3) ≤ sSup S := by
    by_cases hc : ‖f'‖ ≤ ε' / 3
    · obtain ⟨s, hs⟩ := hne
      have hs0 : 0 ≤ s := by obtain ⟨h, _, rfl⟩ := hs; positivity
      linarith [le_csSup hbdd hs]
    · push_neg at hc
      obtain ⟨x, hx1, hx⟩ := f'.exists_lt_apply_of_lt_opNorm (by linarith : ‖f'‖ - ε' / 3 < ‖f'‖)
      have hx0 : x ≠ 0 := by rintro rfl; simp at hx; linarith
      have hxpos : 0 < ‖x‖ := norm_pos_iff.mpr hx0
      set h := (ε / ‖x‖) • x with hdef
      have hh : ‖h‖ = ε := hscale x hxpos
      have hmem : ‖f (d + h) - f d‖ / ‖h‖ ∈ S :=
        ⟨h, ⟨by rw [hh]; exact hε, by rw [hh]⟩, rfl⟩
      have h1 := (abs_le.mp (key h (by rw [hh]; exact hε) (by rw [hh]))).1
      have h2 : ‖f' h‖ / ‖h‖ = ‖f' x‖ / ‖x‖ := by
        rw [hh, hdef, map_smul, norm_smul, Real.norm_eq_abs, abs_of_pos (by positivity)]
        field_simp
      have h3 : ‖f' x‖ ≤ ‖f' x‖ / ‖x‖ := by
        rw [le_div_iff₀ hxpos]
        nlinarith [norm_nonneg (f' x)]
      have := le_csSup hbdd hmem
      linarith
  show dist (sSup S) ‖f'‖ < ε'
  rw [Real.dist_eq, abs_lt]
  constructor <;> linarith

/-- For a scalar function, `S_f(d) = |f'(d)|`. -/
theorem hasSensitivity_of_hasDerivAt {f : ℝ → ℝ} {f' d : ℝ} (hf : HasDerivAt f f' d) :
    HasSensitivity f d |f'| := by
  have := hasSensitivity_of_hasFDerivAt hf.hasFDerivAt
  rwa [ContinuousLinearMap.norm_smulRight_apply, norm_one, one_mul, Real.norm_eq_abs] at this

/-- **`forward error ≤ S_f(d) ‖δd‖ + o(‖δd‖)`.** -/
theorem forward_le_sens {f : E → F} {f' : E →L[ℝ] F} {d : E} (hf : HasFDerivAt f f' d) :
    ∃ g : E → ℝ, g =o[𝓝 0] (fun h => ‖h‖) ∧ ∀ h, ‖f (d + h) - f d‖ ≤ ‖f'‖ * ‖h‖ + g h := by
  refine ⟨fun h => ‖f (d + h) - f d - f' h‖, ?_, fun h => ?_⟩
  · exact (hasFDerivAt_iff_isLittleO_nhds_zero.mp hf).norm_left.trans_isBigO
      (isBigO_refl _ _ |>.norm_right)
  · calc ‖f (d + h) - f d‖ = ‖f' h + (f (d + h) - f d - f' h)‖ := by congr 1; abel
      _ ≤ ‖f' h‖ + ‖f (d + h) - f d - f' h‖ := norm_add_le _ _
      _ ≤ ‖f'‖ * ‖h‖ + ‖f (d + h) - f d - f' h‖ := add_le_add (f'.le_opNorm h) le_rfl

/-- The relative condition number `κ_f(d) = S_f(d) ‖d‖ / ‖f(d)‖`. -/
noncomputable def relCond (S : ℝ) (f : E → F) (d : E) : ℝ := S * ‖d‖ / ‖f d‖

/-- `κ_f(d)` bounds the amplification of relative errors:
`‖f(d+δd) - f(d)‖/‖f(d)‖ ≤ κ_f(d) ‖δd‖/‖d‖ + o(‖δd‖)`. -/
theorem relative_le_cond {f : E → F} {f' : E →L[ℝ] F} {d : E} (hf : HasFDerivAt f f' d)
    (hd : d ≠ 0) (hfd : f d ≠ 0) :
    ∃ g : E → ℝ, g =o[𝓝 0] (fun h => ‖h‖) ∧ ∀ h,
      ‖f (d + h) - f d‖ / ‖f d‖ ≤ relCond ‖f'‖ f d * (‖h‖ / ‖d‖) + g h := by
  obtain ⟨g, hg, hle⟩ := forward_le_sens hf
  have hd' : 0 < ‖d‖ := norm_pos_iff.mpr hd
  have hfd' : 0 < ‖f d‖ := norm_pos_iff.mpr hfd
  refine ⟨fun h => 1 / ‖f d‖ * g h, hg.const_mul_left _, fun h => ?_⟩
  rw [relCond, div_le_iff₀ hfd']
  calc ‖f (d + h) - f d‖ ≤ ‖f'‖ * ‖h‖ + g h := hle h
    _ = (‖f'‖ * ‖d‖ / ‖f d‖ * (‖h‖ / ‖d‖) + 1 / ‖f d‖ * g h) * ‖f d‖ := by field_simp

end General

/-! ## Linear systems -/

section Linear
variable {n : ℕ}

/-- With `A` fixed, `b ↦ A⁻¹ b` has sensitivity `‖A⁻¹‖₂`. -/
theorem sens_linear [NeZero n] (A : Matrix (Fin n) (Fin n) ℝ) (b : EuclideanSpace ℝ (Fin n)) :
    HasSensitivity (Matrix.toEuclideanCLM (n := Fin n) (𝕜 := ℝ) A⁻¹) b (norm2 A⁻¹) := by
  haveI : Nontrivial (EuclideanSpace ℝ (Fin n)) := by
    refine ⟨⟨EuclideanSpace.single 0 1, 0, fun h => ?_⟩⟩
    have := congrArg (fun v => v 0) h
    simp at this
  have := hasSensitivity_of_hasFDerivAt
    ((Matrix.toEuclideanCLM (n := Fin n) (𝕜 := ℝ) A⁻¹).hasFDerivAt (x := b))
  rwa [← Matrix.cstar_norm_def] at this

/-- The relative condition number of `b ↦ A⁻¹ b` is at most `κ₂(A)`. -/
theorem relCond_linear_le {A : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    (b : EuclideanSpace ℝ (Fin n)) (hb : b ≠ 0) :
    relCond (norm2 A⁻¹) (Matrix.toEuclideanCLM (n := Fin n) (𝕜 := ℝ) A⁻¹) b ≤ cond2 A := by
  set T := Matrix.toEuclideanCLM (n := Fin n) (𝕜 := ℝ)
  have hb' : b = T A (T A⁻¹ b) := by
    rw [← ContinuousLinearMap.mul_apply, ← map_mul, Matrix.mul_nonsing_inv _ hA, map_one]
    rfl
  have hx : T A⁻¹ b ≠ 0 := by
    intro h; rw [h, map_zero] at hb'; exact hb hb'
  have hxpos : 0 < ‖T A⁻¹ b‖ := norm_pos_iff.mpr hx
  have hbx : ‖b‖ ≤ norm2 A * ‖T A⁻¹ b‖ := by
    conv_lhs => rw [hb']
    exact (T A).le_opNorm _
  rw [relCond, cond2, div_le_iff₀ hxpos]
  calc norm2 A⁻¹ * ‖b‖ ≤ norm2 A⁻¹ * (norm2 A * ‖T A⁻¹ b‖) :=
        mul_le_mul_of_nonneg_left hbx (norm2_nonneg _)
    _ = norm2 A * norm2 A⁻¹ * ‖T A⁻¹ b‖ := by ring

end Linear

end LU
