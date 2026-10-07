import Mathlib

/-!
# The arithmetic of the algorithms

Elimination uses two operations: the multiplier `div a b = a / b` and the update
`upd a b c = a - b * c`. The triangular solves use the same two operations. The
arithmetic is a parameter `o : Ops K`:

* `Ops.exact`: exact arithmetic in a field `K` (for example `ℝ` or `ℂ`);
* `Ops.rounded r`: real arithmetic rounded by `r`, with one rounding for the product and
  one for the subtraction, as in `A[k+1:, k+1:] -= np.outer(...)`;
* `Ops.fused r`: the update is a fused multiply-add.

`Ops.StdModel o u` is the standard model of floating-point arithmetic with unit
roundoff `u`.
-/
namespace LU

/-- The elementary operations. `upd a b c` computes `a - b * c`. -/
structure Ops (K : Type*) where
  div : K → K → K
  upd : K → K → K → K

namespace Ops

/-- Exact arithmetic in a field. -/
def exact (K : Type*) [Field K] : Ops K :=
  ⟨fun a b => a / b, fun a b c => a - b * c⟩

@[simp] lemma exact_div (K : Type*) [Field K] (a b : K) : (exact K).div a b = a / b := rfl
@[simp] lemma exact_upd (K : Type*) [Field K] (a b c : K) : (exact K).upd a b c = a - b * c := rfl

/-- Real arithmetic rounded by `r`: the product and the subtraction in an update are
rounded separately. -/
noncomputable def rounded (r : ℝ → ℝ) : Ops ℝ :=
  ⟨fun a b => r (a / b), fun a b c => r (a - r (b * c))⟩

/-- Rounded arithmetic in which each update is a fused multiply-add. -/
noncomputable def fused (r : ℝ → ℝ) : Ops ℝ :=
  ⟨fun a b => r (a / b), fun a b c => r (a - b * c)⟩

/-- The standard model of floating-point arithmetic with unit roundoff `u`. An update
has one relative error in the product and one in the subtraction (a fused
multiply-add is the case `ε = 0`). -/
structure StdModel (o : Ops ℝ) (u : ℝ) : Prop where
  div : ∀ a b, ∃ δ, |δ| ≤ u ∧ o.div a b = a / b * (1 + δ)
  upd : ∀ a b c, ∃ δ ε, |δ| ≤ u ∧ |ε| ≤ u ∧ o.upd a b c = (a - b * c * (1 + ε)) * (1 + δ)

lemma exists_rel {u : ℝ} {r : ℝ → ℝ} (hr : ∀ x, |r x - x| ≤ u * |x|) (x : ℝ) :
    ∃ δ, |δ| ≤ u ∧ r x = x * (1 + δ) := by
  have hu : 0 ≤ u := by
    have := hr 1
    simp only [abs_one, mul_one] at this
    exact le_trans (abs_nonneg _) this
  by_cases hx : x = 0
  · refine ⟨0, by simpa using hu, ?_⟩
    have := hr 0
    simp only [hx, sub_zero, abs_zero, mul_zero] at this ⊢
    simpa using abs_nonpos_iff.mp this
  · refine ⟨(r x - x) / x, ?_, ?_⟩
    · rw [abs_div, div_le_iff₀ (abs_pos.mpr hx)]
      exact hr x
    · field_simp
      ring

theorem stdModel_exact {u : ℝ} (hu : 0 ≤ u) : (exact ℝ).StdModel u where
  div a b := ⟨0, by simpa using hu, by simp⟩
  upd a b c := ⟨0, 0, by simpa using hu, by simpa using hu, by simp⟩

theorem stdModel_rounded {u : ℝ} {r : ℝ → ℝ} (hr : ∀ x, |r x - x| ≤ u * |x|) :
    (rounded r).StdModel u where
  div a b := exists_rel hr _
  upd a b c := by
    obtain ⟨ε, hε, he⟩ := exists_rel hr (b * c)
    obtain ⟨δ, hδ, hd⟩ := exists_rel hr (a - r (b * c))
    exact ⟨δ, ε, hδ, hε, by simp only [rounded]; rw [hd, he]⟩

theorem stdModel_fused {u : ℝ} {r : ℝ → ℝ} (hr : ∀ x, |r x - x| ≤ u * |x|) :
    (fused r).StdModel u where
  div a b := exists_rel hr _
  upd a b c := by
    have hu : 0 ≤ u := by
      have := hr 1
      simp only [abs_one, mul_one] at this
      exact le_trans (abs_nonneg _) this
    obtain ⟨δ, hδ, hd⟩ := exists_rel hr (a - b * c)
    exact ⟨δ, 0, hδ, by simpa using hu, by simp only [fused]; rw [hd]; ring⟩

end Ops

end LU
