# Lean verification of the Cholesky chapter

This project checks the mathematical statements and proofs of
[`content/cholesky.md`](../../content/cholesky.md) in Lean 4 with mathlib. It covers
the main chapter and the optional section on numerical breakdown, including the
counterexample family. The only exceptions are results the chapter cites without
proof and the reported numerical experiments; see
[What is not proved in Lean](#what-is-not-proved-in-lean).

## Build and check

The project pins Lean and mathlib to **v4.24.0**. With `elan` installed, run from this
directory:

```sh
lake exe cache get
lake build
lake env lean Audit.lean
```

`lake build` checks every proof with warnings treated as errors. `Audit.lean` prints
the axioms used by the main theorems. Only Lean's standard `propext`,
`Classical.choice`, and `Quot.sound` occur. The project contains no `sorry`, no added
axioms, and no `native_decide`.

The Lean sources are in [`src/`](src/). The root module [`Cholesky.lean`](Cholesky.lean)
imports all of them. Module names follow the folders (for example
`src.Counterexample.Matrix`), while every theorem lives in the namespace `Cholesky`.

The numerical experiments are cross-checked separately, in Python, by
[`numerics/check_table.py`](numerics/check_table.py) (requires NumPy).

## How the chapter is modelled

* **Matrices** are `Matrix (Fin n) (Fin n) ℝ`. Indices start at 0, while the notes
  count from 1. SPD is mathlib's `Matrix.PosDef`; `posDef_iff_real` shows that this is
  `Aᵀ = A` and `xᵀ A x > 0` for `x ≠ 0`. The quantity `n^{3/2}` is written `n * √n`.
* **The in-place algorithm** ([`Algorithm.lean`](src/Algorithm.lean)) follows `cholesky_in_place`. At step
  `k` it takes the square root of the pivot, divides the entries below it, and updates
  `a_ij ← a_ij - a_ik a_jk` for `k < j ≤ i`. `pivot o A k` is the pivot examined at step
  `k`, and `Succeeds o A` says that every pivot is positive. `factor o A` is the
  computed lower triangle. The routine reads only the lower triangle of its input
  array, and `symmOfLower A` is the symmetric matrix that triangle defines.
* **Arithmetic** is a parameter `o : Ops`:
  * `Ops.exact` is exact real arithmetic;
  * `Ops.rounded r` rounds the product and the subtraction of each update separately,
    as NumPy does;
  * `Ops.fused r` uses a fused multiply-add.

  `Ops.StdModel o u` is the standard model of floating-point arithmetic with unit
  roundoff `u`.
* **Floating point** ([`Float.lean`](src/Float.lean)): `IsFloat p x` means `x = m · 2^e` with `|m| < 2^p`,
  and the exponent range is unbounded, so there is no underflow or overflow.
  `IsRoundNearest p r` says that `r x` is a nearest float, with ties broken in any
  way. `exists_roundNearest` shows that such roundings exist.
  `roundNearest_relErr` and `stdModel_of_roundNearest` give the standard model with
  `u = 2^{-p}` (`unitRoundoff p`).
* **Norms** ([`Norms.lean`](src/Norms.lean)): `vnorm` is the Euclidean norm. `norm2` is the spectral
  norm, i.e. mathlib's ℓ² operator norm. `frob` is the Frobenius norm, and
  `cond2 A = ‖A‖₂ ‖A⁻¹‖₂`.
* **Triangular solves** ([`Solve.lean`](src/Solve.lean)): `fwdSolve` subtracts `lᵢⱼ yⱼ` in increasing
  order of `j`, and `bwdSolve` subtracts in decreasing order. The notes do not fix an
  order, and the proofs work for these ones.

## Where each statement is proved

The tables follow the sections of the notes. Theorem names are relative to the
namespace `Cholesky`; for example, `CE.Par.A_fails` is `Cholesky.CE.Par.A_fails`.

### Symmetric positive definite matrices — [`SPD.lean`](src/SPD.lean), [`Complex.lean`](src/Complex.lean)

| Statement | Lean |
| --- | --- |
| `xᵀAx = zᵀΛz = ∑ λᵢ zᵢ²` with `z = Qᵀx` | `quad_eq_sum_eigenvalues`, `eigQ_dot` |
| A symmetric matrix is SPD iff all eigenvalues are positive (proved as in the notes) | `posDef_iff_eigenvalues_pos` |
| An SPD matrix is nonsingular; `aᵢᵢ = eᵢᵀAeᵢ > 0` | `posDef_isUnit`, `posDef_diag_pos` |
| `!![1,2;2,1]` has positive diagonal, eigenvalues `3` and `-1`, and is not SPD | `example12_diag_pos`, `example12_eigenvalues`, `example12_not_posDef` |
| Complex case: Hermitian PD iff `A = LLᴴ` with positive real diagonal; unique | `Herm.cholesky_iff`, `Herm.cholesky_unique` |

### Deriving the factorization — [`BlockStep.lean`](src/BlockStep.lean)

| Statement | Lean |
| --- | --- |
| Block formula for `L Lᵀ` | `mul_transpose_zero_zero`, `mul_transpose_succ_zero`, `mul_transpose_succ_succ` |
| Matching blocks: `l₁₁ = √a₁₁`, `ℓ = c/l₁₁`, `L_S L_Sᵀ = B - ℓℓᵀ = S` | `block_equations` |
| `S` is symmetric | `schur_transpose` |
| Completing the square | `complete_square` |
| `S` is SPD | `schur_posDef` |

### Existence and uniqueness — [`Factorization.lean`](src/Factorization.lean)

| Statement | Lean |
| --- | --- |
| **Theorem** `thm:cholesky_existence` | `cholesky_iff`, `cholesky_existsUnique`, `cholesky_unique` |
| Converse: `xᵀLLᵀx = ‖Lᵀx‖² > 0` | `posDef_of_factorization`, `quad_mul_transpose` |
| Every pivot is positive in exact arithmetic | `exact_succeeds_of_posDef`, `exact_succeeds_iff` |
| Positive semidefinite matrices can have zero pivots ([`SPD.lean`](src/SPD.lean)) | `psdExample_posSemidef`, `psdExample_zero_pivot` |

### Connection with LU — [`LDL.lean`](src/LDL.lean)

| Statement | Lean |
| --- | --- |
| `L₀ = L D^{-1/2}` is unit lower triangular | `ldlL0_unit_lower` |
| `A = L₀ D L₀ᵀ = L₀ U`, `U = D L₀ᵀ` upper triangular | `ldl_eq`, `lu_eq`, `ldlU_upper` |
| The LU pivots are `lᵢᵢ² > 0` | `ldlU_diag`, `exact_pivot_eq_sq` |
| Cholesky and unpivoted LU give the same factorization | `lu_unique`, `cholesky_is_lu` |

### The in-place algorithm — [`Algorithm.lean`](src/Algorithm.lean), [`Factorization.lean`](src/Factorization.lean)

| Statement | Lean |
| --- | --- |
| At completion the lower triangle contains `L` | `exact_factor_spec`, `exact_factor_mul_transpose` |
| Only the lower triangle is read and overwritten; it is interpreted as a symmetric matrix | `factor_lowerEq`, `factor_symmOfLower`, `run_upper` |
| After the first step the algorithm acts on the trailing block | `run_succ_submatrix`, `pivot_succ`, `factor_succ_succ` |

### Arithmetic and storage — [`Cost.lean`](src/Cost.lean)

| Statement | Lean |
| --- | --- |
| Step `k` updates `m(m+1)/2` entries (`m = n-1-k`), exactly those changed by `stepAt` | `card_stepAt_updates`, `two_mul_updEntries` |
| `∑_{m=1}^{n-1} m(m+1) = (n³ - n)/3` | `sum_update_flops`, `sum_update_flops_real` |
| Total with divisions and square roots; leading cost `n³/3` | `cholFlops_split`, `cholFlops_real`, `cholFlops_ratio` |
| LU costs `2n³/3`; Cholesky costs about half | `luFlops_real`, `chol_lu_ratio` |
| Updating both triangles costs `~2n³/3` | `fullUpdateFlops_ratio` |
| `n(n+1)/2` stored entries | `card_lower_triangle` |

### Solving a linear system — [`Solve.lean`](src/Solve.lean), [`Cost.lean`](src/Cost.lean)

| Statement | Lean |
| --- | --- |
| `Ly = b`, `Lᵀx = y` solve `Ax = b` | `exact_solve`, `exact_fwdSolve`, `exact_bwdSolve` |
| The two solves cost `2n²` flops | `triSolveFlops_eq`, `two_solves_flops` |

### Why the factor entries stay controlled — [`FactorBounds.lean`](src/FactorBounds.lean), [`Stability.lean`](src/Stability.lean)

| Statement | Lean |
| --- | --- |
| `aᵢᵢ = ∑ₖ lᵢₖ²` and `|lᵢₖ| ≤ √aᵢᵢ` | `diag_eq_sum_sq_le`, `abs_entry_le_sqrt_diag` |
| `(|L||L|ᵀ)ᵢⱼ = ∑ₖ |lᵢₖ lⱼₖ| ≤ √(aᵢᵢ aⱼⱼ)` | `absProd_le` |
| In exact arithmetic `‖|L||L|ᵀ‖₂ ≤ tr A ≤ n‖A‖₂` | `exact_absProd_bound` |

### Backward error of the factorization — [`BackwardError.lean`](src/BackwardError.lean), [`Stability.lean`](src/Stability.lean), [`Norms.lean`](src/Norms.lean)

| Statement | Lean |
| --- | --- |
| **Theorem** `thm:backward_error_cholesky`: `A + E = L̂L̂ᵀ`, `E` symmetric, `|E| ≤ γ_{n+1}|L̂||L̂|ᵀ` (no SPD assumption needed) | `backward_error` |
| Diagonal equations `∑ₖ l̂ᵢₖ² ≤ aᵢᵢ + γ ∑ₖ l̂ᵢₖ²` | `diag_equation` |
| `‖L̂‖_F² ≤ tr(A)/(1-γ)` | `frob_sq_le` |
| `‖|L̂|‖₂ ≤ ‖L̂‖_F`, `tr A ≤ n‖A‖₂`, `|E| ≤ B ⇒ ‖E‖₂ ≤ ‖B‖₂` | `norm2_le_frob`, `frob_absMat`, `trace_le`, `norm2_le_of_abs_le` |
| `‖E‖₂ ≤ γ‖|L̂||L̂|ᵀ‖₂ ≤ γ‖L̂‖_F² ≤ nγ/(1-γ)‖A‖₂` | `norm2_E_le`, `norm2_absProd_le_frob_sq`, `norm2_E_le_rel` |
| A completed factorization is the Cholesky factor of a nearby SPD matrix | `nearby_spd` |
| `nγ_{n+1}/(1-γ_{n+1}) = n(n+1)u/(1-2(n+1)u) ≈ n(n+1)u` | `stability_constant_eq`, `stability_constant_asymp` |

### Backward and forward error of the computed solution — [`Solve.lean`](src/Solve.lean), [`Accuracy.lean`](src/Accuracy.lean)

| Statement | Lean |
| --- | --- |
| `η(x̂)` is the smallest relative perturbation; for `b ≠ 0` it equals the residual formula | `rigal_gaches`, `eta_eq` |
| `(A+ΔA)x̂ = b`, `|ΔA| ≤ (γ_{n+1} + 2γ_n + γ_n²)|L̂||L̂|ᵀ`, and `γ_{n+1} + 2γ_n + γ_n² ≤ γ_{3n+1}` | `solve_backward_error`, `gamma_solve_le`, `fwdSolve_backward`, `bwdSolve_backward` |
| `η(x̂) ≤ ‖ΔA‖₂/‖A‖₂` | `eta_le_of_perturbation` |
| `η(x̂) ≤ β_n = nγ_{3n+1}/(1-γ_{n+1})` and `β_n ≈ n(3n+1)u` | `eta_solve_le`, `beta_asymp` |
| `‖x̂ - x‖/‖x‖ ≤ 2κβ/(1-κβ)` | `forward_error` |
| For SPD matrices `κ₂ = λ_max/λ_min` | `cond2_spd` |

### Optional section: guarantees of successful completion — [`Completion.lean`](src/Completion.lean), [`Examples.lean`](src/Examples.lean), [`Counterexample/Remarks.lean`](src/Counterexample/Remarks.lean)

| Statement | Lean |
| --- | --- |
| Meinguet's theorem uses `(1+f_{n-1})^{n-1}`, larger than the product | `meinguet_prod_le` |
| `1.1 n^{3/2} uκ ≤ 1` implies `κ μ_n(u) < 1` for `n ≥ 576`, `u ≤ 2^{-17}` | `meinguet_condition` |
| `∑(√j + 2) < (2/3)n^{3/2} + 2n ≤ (3/4)n^{3/2}`, the log estimates, and concavity | `sum_update_weights_steps`, `meinguet_constant`, `log_chord` |

### Optional section: the counterexample — [`Counterexample/`](src/Counterexample/), [`Hadamard.lean`](src/Hadamard.lean)

| Statement | Lean |
| --- | --- |
| Definition of `c*`; `κ < c*/u` guarantees success | `cStar`, `succeeds_of_cond_lt` |
| Admissible `(p, q)`; `3q+8 ≤ p ⟺ k^{3/2}u ≤ 1/256` | `CE.admissible_iff`, `CE.size_iff` |
| Walsh–Hadamard matrices; `QᵀQ = I`; entries `±2^{-q}` | `hadamard_mul_transpose`, `hadamard_entry`, `hadQ_transpose_mul`, `hadQ_entry` |
| `t² = 9h/8`; `t`, `t²` and all stored entries are floats | `CE.Par.t_sq_eq`, `CE.Par.t_float`, `CE.Par.tt_float`, `CE.Par.T_float`, `CE.Par.A_float` |
| `η`, `δ` are multiples of `4u`; the cross shifts are multiples of `4u/√k`; `η√k = 12ku ≤ 3/512 < 1/4`; `48k < 2^{2q+6}`, `2q+6 ≤ p`; `-η/8` and `9/(4k)` are multiples of `4u` | `CE.Par.spacing_facts` |
| `TᵀT = G + (9/8)ηJ`; the exact Schur complement is `S = δI - (η/8)J` | `CE.Par.T_gram`, `CE.Par.C_eq`, `CE.Par.schur_eq` |
| `λ_min(S) = (3/4)kη`, `λ_max(S) = kη`; the stored matrix is SPD | `CE.Par.S_quad_bounds`, `CE.Par.S_ones`, `CE.Par.S_perp`, `CE.Par.A_posDef` |
| The leading block is factored exactly and the trailing block receives the rounded updates | `lead_run`, `lead_pivot`, `leadState_step` |
| Cross entries decrease by `2h` per step; diagonal entries are unchanged; off-diagonal updates are exact | `CE.Par.round_cross`, `CE.Par.round_diag`, `CE.Par.round_offv`, `CE.Par.phase1` |
| The Hadamard rows cancel `G` exactly | `CE.Par.phase2` |
| The stored trailing matrix is `Ŝ` | `CE.Par.trailing_eq_Shat` |
| `λ_min(Ŝ) = -(k-9)η/8`, `λ_max(Ŝ) = (15k+9)η/8`, `‖Ŝ‖₂ ≤ (15k+9)η/8`, ratio `> 1/18` | `CE.Par.Shat_ones`, `CE.Par.Shat_alt`, `CE.Par.Shat_quad_bounds`, `CE.Par.norm2_Shat_le`, `negative_eigenvalue_margin` |
| Breakdown proof: sign pattern, `|L̂| = 2D - L̂`, `‖|L̂|‖₂ ≤ 3‖L̂‖₂`, `‖E‖₂ < ‖Ŝ‖₂/100`, contradiction | `factor_offdiag_nonpos`, `absMat_eq_of_signs`, `norm2_absMat_le_three`, `perturbation_fraction`, `update_budget`, `no_success` |
| **Cholesky fails on `A_{9k}`**, with separate rounding and with fused multiply-add | `CE.Par.A_fails` |
| `‖T‖₂² ≤ 9/2 + 27x`, `‖A‖₂ ≤ 17/2 + 39x`, `‖A⁻¹‖₂ ≤ 1 + 17/(72x)` | `CE.Par.norm2_T_sq_le`, `CE.Par.norm2_A_le`, `CE.Par.norm2_Ainv_le` |
| `κ₂(A_{9k}) < 57/(n^{3/2}u)`; the constant `3674685/65536 < 57` | `CE.Par.cond_A_le`, `CE.Par.cond_A_lt`, `conditionUpper_lt_57`, `constant_57` |
| `‖T‖₂² ≥ 9/2`, `λ_max(A) ≥ 17/2`, `‖A⁻¹‖₂ ≥ 17/(8δ)`, `κ₂ ≥ (2601/64)/(n^{3/2}u)` | `CE.Par.norm2_T_sq_ge`, `CE.Par.norm2_A_ge`, `CE.Par.norm2_Ainv_ge`, `CE.Par.cond_A_ge`, `CE.Par.cond_A_ge'` |
| At admissible sizes, `c* < 57/n^{3/2}` (and `≥ 1/(1.1n^{3/2})` given Meinguet's theorem) | `CE.cStar_admissible`, `cStar_ge_meinguet` |
| Fixed precision `p ≥ 20`: largest block, `x* ∈ {2⁻⁸,2⁻¹⁰,2⁻¹²}`, `κ₂ ≤ 303691615/36864 < 8239` | `largestQ_properties`, `CE.fixed_precision` |
| The bound `289/(144x) + 425/24 + 39x` is largest at `x* = 2^{-12}` | `conditionUpper_expand`, `conditionUpper_largest` |
| Padding: failure in the original block, same condition number, SPD, stored | `fails_of_leading`, `cond2_blkdiag`, `CE.Par.B_fails`, `CE.Par.cond2_B`, `CE.Par.B_posDef`, `CE.Par.B_float` |
| binary64 (`q* = 15`, `x* = 2⁻⁸`, `κ₂ < 532`) and binary32 (block size `2304`, `x = 1/4096`, `κ₂ < 8239`) | `CE.binary64`, `CE.binary32` |
| For `p ≥ 20` and every `N ≥ 2304`, a failing stored SPD matrix with `κ₂ < 8239 max{1, 1/(N^{3/2}u)}` | `CE.sharp_family` |
| `max{u, 1/(1.1N^{3/2})} ≤ c* < 8239 max{u, 1/N^{3/2}}` | `CE.cStar_bounds`, `CE.cStar_lt`, `cStar_ge_u`, `cond2_ge_one` |
| The two evaluation orders can round differently | `order_matters` |
| The family has `κ₂ ≥ 1156/3 > 385` | `CE.Par.cond_A_ge_385`, `family_lower_bound` |
| Table parameters and upper bounds; ratio about 64; memory sizes `162 GiB` and `2.82 GiB` | `table_parameters`, `table_bounds`, `table_ratio`, `memory_sizes` |
| Not monotone in precision: `(1/32)[[4,5],[5,7]]` succeeds with `p = 3`, fails with `p = 4` | `nmA_posDef`, `nmA_float`, `nonmonotone_p3`, `nonmonotone_p4` |
| `H = D_A⁻¹AD_A⁻¹` has unit diagonal; `H = I` for a positive diagonal matrix; `κ₂(PAPᵀ) = κ₂(A)` | `scaled_unit_diag`, `scaled_diag_eq_one`, `cond2_perm` |

## What is not proved in Lean

* **Cited theorems.** Several results the chapter states without proof are not
  formalized: Wilkinson's condition `20 n^{3/2} uκ ≤ 1`, Meinguet's completion
  theorem, Kiełbasiński's theorem, Demmel's criterion, and the bounds of Sun and of
  Chang. Wherever a result depends on Meinguet's theorem, the theorem is an explicit
  hypothesis (`MeinguetGuarantee`), not an axiom. This affects only the lower bounds
  `c* ≥ u/μ_n(u)` and `c* ≥ 1/(1.1 n^{3/2})`. Bibliographic details such as page
  numbers and Wilkinson's enlarged rounding parameter are not checked.
* **Numerical experiments.** The failing steps and computed pivots in the table
  (565, 2252, 9002) are checked in Python by `numerics/check_table.py`, not in Lean.
  The script also checks that accumulating each inner product from zero and then
  subtracting it succeeds on these examples. Lean does prove the facts that make these
  experiments reliable:
  * the first `7k` steps produce exactly `Ŝ` (`trailing_eq_Shat`);
  * the remaining steps act on `Ŝ` alone (`run_natAdd`);
  * padding leaves every pivot of the leading block unchanged (`pivot_castAdd`).

  So the padded runs must fail at the same step with the same pivot.
* **Qualitative remarks**, such as "there is no need to form `A⁻¹`" or "the elapsed
  time also depends on the implementation", are not formal statements.

The files [`Bounds.lean`](src/Bounds.lean), [`Completion.lean`](src/Completion.lean), [`Parameters.lean`](src/Parameters.lean), [`Perturbation.lean`](src/Perturbation.lean)
and [`Schur.lean`](src/Schur.lean) contain the first, partial formalization. They are still used: for
example, `meinguet_condition` and the scalar bounds come from there.
