# Lean verification of the LU-with-pivoting chapter

This project checks the mathematical statements and proofs of
[`content/lu_pivoting.md`](../../content/lu_pivoting.md) in Lean 4 with mathlib. It covers
the main chapter and the optional section on a right-hand side with a large backward
error. The exceptions are results the chapter cites without proof, informal remarks,
and the printed numerical experiment; see
[What is not proved in Lean](#what-is-not-proved-in-lean).

## Build and check

The project pins Lean and mathlib to **v4.24.0**, the same versions as
[`../Cholesky`](../Cholesky). With `elan` installed, run from this directory:

```sh
lake exe cache get
lake build
lake env lean Audit.lean
```

`lake build` checks every proof with warnings treated as errors. `Audit.lean` prints the
axioms used by 135 main theorems. Only Lean's standard `propext`, `Classical.choice`, and
`Quot.sound` occur. The project contains no `sorry`, no added axioms, and no
`native_decide`.

The Lean sources are in [`src/`](src/). The root module [`LU.lean`](LU.lean) imports all
of them, and every theorem lives in the namespace `LU`.

The numerical claims are cross-checked separately, in Python, by
[`numerics/check_examples.py`](numerics/check_examples.py) (requires NumPy).

## How the chapter is modelled

* **Matrices** are `Matrix (Fin n) (Fin n) 𝕜`. Indices start at 0, while the notes count
  from 1. The algorithm and the existence theorem are proved over any normed field `𝕜`,
  in particular `ℝ` and `ℂ`. The error analysis is over `ℝ`, as in the notes.
* **The in-place routine** ([`Algorithm.lean`](src/Algorithm.lean)) follows
  `lu_factorization_with_row_pivoting`. Its state is the array and the row permutation
  (the routine's `P`). Step `k` does the following:
  * It chooses `pivRow M k`, the *first* row `p ≥ k` with `|m_pk|` maximal. This is
    `k + np.argmax(...)`, see `IsFirstMax` in [`Pivot.lean`](src/Pivot.lean).
  * It swaps the entire rows `k` and `p` (`swapRows`).
  * It runs `elim`: if the pivot is nonzero, it stores the multipliers and updates the
    trailing block; a zero pivot is skipped.

  The routine runs `k = 0, …, n - 2`. `final o A` runs one more step, which changes
  nothing (`run_pred`). The factors are read from the packed array: `factorL o A`
  (unit diagonal implicit) and `factorU o A`. `factorP o A` is the permutation matrix,
  with `(P A)ᵢⱼ = a_{perm i, j}`.
* **Arithmetic** is a parameter `o : Ops 𝕜` ([`Ops.lean`](src/Ops.lean)):
  * `Ops.exact 𝕜` is exact arithmetic.
  * `Ops.rounded r` rounds each multiplier, and rounds the product and the subtraction
    of each update separately, as NumPy does.
  * `Ops.fused r` uses a fused multiply-add.

  `Ops.StdModel o u` is the standard model with unit roundoff `u`. The pivot search
  compares the computed values exactly.
* **Floating point** ([`Float.lean`](src/Float.lean)): `IsFloat p x` means `x = m · 2^e`
  with `|m| < 2^p`. The exponent range is unbounded, so there is no underflow or
  overflow. `IsRoundNearest p r` says that `r x` is a nearest float, with ties broken in
  any way. Such roundings exist (`exists_roundNearest`). They satisfy the standard model
  with `u = 2^{-p}` (`roundNearest_relErr`, `stdModel_of_roundNearest`).
* **Norms** ([`Norms.lean`](src/Norms.lean)):
  * `vnorm` is the Euclidean norm.
  * `norm2` is the spectral norm, i.e. mathlib's ℓ² operator norm.
  * `frob` is the Frobenius norm.
  * `cond2 A = ‖A‖₂ ‖A⁻¹‖₂`.
* **Backward error** ([`Residual.lean`](src/Residual.lean)) follows the boxed definition
  of the notes:
  * `tau A E b e = max {‖E‖/‖A‖, ‖e‖/‖b‖}`;
  * `tauSet A b x` holds these sizes for all pairs with `(A + E) x = b + e`;
  * `eta A b x` is the infimum of `tauSet A b x`.

  `eta_isLeast` shows that the minimum is attained.
* **Triangular solves.** For the backward-error theorem, `fwdSolve` and `bwdSolve`
  ([`Solve.lean`](src/Solve.lean)) subtract `lᵢⱼ yⱼ` one term at a time. The bound holds for
  any order of the inner products, and we prove it for this one. In the optional example
  the order matters, so [`LargeError.lean`](src/LargeError.lean) follows the loops of the
  notes exactly:
  * `fwdLoop` accumulates `s` from zero, then sets `y[i] = pb[i] - s`;
  * `bwdLoop` accumulates the products of `U[i, i+1:] @ x[i+1:]` from left to right.

## Where each statement is proved

The tables follow the sections of the notes. Theorem names are relative to the namespace
`LU`.

### Why interchange rows? — [`Examples.lean`](src/Examples.lean), [`Cancellation.lean`](src/Cancellation.lean)

| Statement | Lean |
| --- | --- |
| LU can fail at a zero pivot of an invertible matrix | `swap_no_lu` |
| Unpivoted factors of `A_ε`: `L = [[1,0],[1/ε,1]]`, `U = [[ε,1],[0,π-1/ε]]` | `Aeps_lu` |
| `π = 1/ε + (π - 1/ε)` | `Aeps_cancel` |
| At `ε = 10⁻¹⁸`, binary64 loses `π` when forming `π - 1/ε`; the computed factors reconstruct `a₂₂` as `0` | `Aeps_double`, `Aeps_double_product` |
| `P A_ε = [[1,π],[ε,1]]`, `L = [[1,0],[ε,1]]`, `U = [[1,π],[0,1-επ]]`, computed by the routine | `Aeps_final`, `Aeps_plu` |

### Partial pivoting and `PA = LU` — [`Pivot.lean`](src/Pivot.lean), [`Algorithm.lean`](src/Algorithm.lean), [`Factorization.lean`](src/Factorization.lean)

| Statement | Lean |
| --- | --- |
| A row with `\|a_pk\| = max_{i ≥ k} \|a_ik\|` exists; ties go to the first candidate | `exists_firstMax`, `IsFirstMax.unique`, `firstMax_self` |
| Running step `n - 1` as well changes nothing | `stepAt_last`, `run_pred` |
| After the first step the routine acts on the trailing block; later interchanges reorder the stored multipliers (`Qℓ`) | `rel_final`, `final_row`, `final_col`, `final_trail`, `final_perm` |
| `\|l_ik\| ≤ 1` | `norm_multiplier_le_one`, `norm_factorL_le_one` |
| A zero pivot means the entire active column is zero | `pivot_eq_zero_col` |
| The pivots are the diagonal of `U` | `pivot_eq_diag` |
| A nonsingular matrix has no zero pivots in exact arithmetic | `pivot_ne_zero`, `det_factorU` |
| `PA = LU` with `A` the original input | `exact_plu` |
| `Ly = Pb`, `Ux = y` give `Ax = b`; no permutation of `x` | `solve_of_plu` |

### Why the factorization always exists — [`Factorization.lean`](src/Factorization.lean)

| Statement | Lean |
| --- | --- |
| **Theorem** `thm:existence_lu_pivoting` (every square matrix, `𝔽 = ℝ` or `ℂ`) | `exists_lu_pivoting`, `exists_lu_pivoting_real`, `exists_lu_pivoting_complex` |
| `α = 0` forces `v = 0` and `ℓ = 0`; `v = αℓ` in either case | `swapped_col_eq_zero`, `exact_mult_of_pivot_zero`, `exact_v_eq` |
| `S = B - ℓw` | `exact_trail_eq` |
| `LU = [[α, w], [Qv, QℓW + QS]] = PA` (the induction) | `exact_plu` with `mul_zero_apply`, `mul_succ_zero`, `mul_succ_succ` |
| For singular `A`, `U` has a zero diagonal entry, and `Ax = b` is not uniquely solvable for every `b` | `exists_diag_eq_zero`, `not_unique_of_singular` |

### In-place implementation — [`Examples.lean`](src/Examples.lean), [`Cost.lean`](src/Cost.lean)

| Statement | Lean |
| --- | --- |
| The 3 × 3 example swaps rows at both steps, and the second swap moves a stored multiplier | `A3_piv0`, `A3_step0`, `A3_piv1`, `A3_step1`, `A3_final` |
| Its factors satisfy `PA = LU`; `y = (17, -5/4, -2/3)`; the computed solution is `x = (1, 0, 1)` | `A3_plu`, `A3_triangular`, `A3_solve`, `exact_luSolve` |
| Step `k`: `n - k` pivot candidates, `n - 1 - k` multipliers, `(n - 1 - k)²` updated entries; nothing else changes | `card_candidates`, `card_multipliers`, `card_updates`, `elim_apply_of_not_mem` |
| Leading arithmetic cost `2n³/3` | `luFlops_real`, `luFlops_ratio` |
| Pivot searches and row swaps add `O(n²)` work | `two_mul_pivotCompares`, `overhead_le` |

### Forward error, backward error, and the residual — [`Residual.lean`](src/Residual.lean), [`Norms.lean`](src/Norms.lean), [`Cost.lean`](src/Cost.lean)

| Statement | Lean |
| --- | --- |
| Pairs making `x̂` exact always exist (`E = 0`, `e = Ax̂ - b`) | `tauSet_nonempty` |
| `η` is at most the size of any such pair | `eta_le_tau` |
| Lower bound: `r = Ex̂ - e` and `‖r‖ ≤ τ d` | `residual_lower` |
| The attaining pair; `‖r x̂ᴴ‖₂ = ‖r‖₂‖x̂‖₂`; `‖E‖/‖A‖ = ‖e‖/‖b‖ = ‖r‖/d` | `attaining_pair`, `norm2_vecMulVec` |
| For `x̂ = 0` every admissible pair has `e = -b`, and `η(0) = 1` | `eta_zero_pairs`, `eta_zero` |
| **Residual formula** `η(x̂) = ‖r‖/(‖A‖‖x̂‖ + ‖b‖)`, the minimum is attained | `eta_isLeast`, `eta_eq` |
| The scaled residual is at most 1 | `etaFormula_le_one` |
| Computing the residual costs `2n²` flops | `residualFlops_eq` |

### Sensitivity — [`Sensitivity.lean`](src/Sensitivity.lean)

| Statement | Lean |
| --- | --- |
| `S_f(d) = ‖Df(d)‖` for differentiable `f` | `hasSensitivity_of_hasFDerivAt` |
| Scalar case `S_f(d) = \|f'(d)\|` | `hasSensitivity_of_hasDerivAt` |
| Forward error `≤ S_f(d)‖δd‖ + o(‖δd‖)` | `forward_le_sens` |
| `κ_f(d) = S_f(d)‖d‖/‖f(d)‖` bounds the amplification of relative errors | `relCond`, `relative_le_cond` |
| With `A` fixed, `b ↦ A⁻¹b` has sensitivity `‖A⁻¹‖`; using `‖b‖ ≤ ‖A‖‖x‖`, its relative condition number is at most `κ(A)` | `sens_linear`, `relCond_linear_le` |

### The matrix condition number — [`Condition.lean`](src/Condition.lean)

| Statement | Lean |
| --- | --- |
| `κ(A) ≥ 1` since `1 = ‖I‖ ≤ ‖A‖‖A⁻¹‖` (any induced norm; the 2-norm) | `one_le_norm_mul_norm_inv`, `one_le_cond2`, `norm2_one` |
| **`κ₂(A) = σ_max/σ_min`**, with `σᵢ = √λᵢ(AᵀA)` | `norm2_eq_sigmaMax`, `norm2_inv_eq`, `cond2_eq_sigma` |
| `σ_max`, `σ_min` are the largest and smallest stretching of a unit vector | `stretch_bounds`, `sigmaMax_isGreatest`, `sigmaMin_isLeast` |
| `diag(1, ε)`: `x = (1,0)`, `δx = (0,1)` for `δb = (0,ε)`; relative changes `ε` and `1`; `κ₂ = 1/ε` | `diagEps_solve`, `diagEps_relative`, `diagEps_cond` |
| Row permutations preserve the singular values, and `κ₂(PA) = κ₂(A)` | `gram_perm_mul`, `singVal_perm_mul`, `norm2_perm_mul`, `cond2_perm_mul` |

### Forward-error estimates, Banach lemma, perturbation theorem — [`ForwardError.lean`](src/ForwardError.lean)

| Statement | Lean |
| --- | --- |
| `δx = -A⁻¹r` | `error_eq` |
| Rigorous steps 1–4: `‖δx‖/‖x‖ ≤ κη(1 + ‖x̂‖/‖x‖)` | `relErr_le` |
| First-order estimate: if `‖x̂‖ ≤ (1+t)‖x‖` then `‖δx‖/‖x‖ ≤ (2+t)κη` | `relErr_le_of_norm_le` |
| **Lemma** `lem:banach`: partial sums, the Neumann series as inverse, the bound `1/(1 - ‖X‖)` | `neumann_partial`, `banach_inverse`, `banach`, `banach_matrix` |
| `‖(A+E)⁻¹‖ ≤ ‖A⁻¹‖/(1 - ‖A⁻¹‖‖E‖)` | `inv_add_bound` |
| **Theorem** `thm:perturbation_bound` | `perturbation_bound` |
| **Rigorous bound**: `κη < 1 ⇒ ‖x̂ - x‖/‖x‖ ≤ 2κη/(1 - κη)`, via a pair attaining `η` | `forward_error_bound` |

### Backward error of LU and the role of cancellation — [`BackwardError.lean`](src/BackwardError.lean), [`Solve.lean`](src/Solve.lean), [`Stability.lean`](src/Stability.lean), [`Examples.lean`](src/Examples.lean)

| Statement | Lean |
| --- | --- |
| Rounding to nearest satisfies the standard model | `roundNearest_relErr`, `stdModel_of_roundNearest` |
| **Theorem** `thm:backward_error_lu`, factorization: `PA + F = L̂Û`, `\|F\| ≤ γ_n\|L̂\|\|Û\|` | `factor_backward`, `factor_backward'` (with `mult_err`, `upd_err`) |
| Triangular solves: `(T + ΔT)x̂ = b`, `\|ΔT\| ≤ γ_n\|T\|` | `fwdSolve_backward`, `bwdSolve_backward` |
| **Theorem** `thm:backward_error_lu`, complete solve: `(A+E)x̂ = b`, `\|PE\| ≤ γ_{3n}\|L̂\|\|Û\|` | `solve_backward`, `gamma_three_le` |
| `η(x̂) ≤ ‖E‖/‖A‖ ≤ γ_{3n}‖\|L̂\|\|Û\|‖/‖A‖`, using `‖PE‖ = ‖E‖` | `eta_lu_le`, `norm2_perm_mul`, `norm2_le_of_abs_le` |
| `γ_{3n} ≈ 3nu` | `gamma_asymp` |
| `A_ε` without pivoting: the terms reconstructing `a₂₂` total `2/ε - π`; after the swap `\|L\|\|U\| = \|PA\|` | `Aeps_absProd_unpivoted`, `Aeps_absProd_pivoted` |

### Element growth — [`Growth.lean`](src/Growth.lean)

| Statement | Lean |
| --- | --- |
| `\|a_ij - l_ik a_kj\| ≤ \|a_ij\| + \|a_kj\|` | `upd_norm_le` |
| Each step at most doubles the largest active entry; **`ρ_n ≤ 2^{n-1}`** | `stage_bound`, `activeMax_le`, `growth_le` |
| The matrix `wilk n`: ties go to the first candidate, so there are no swaps; every multiplier is `-1`; the last column doubles | `step_wilk`, `run_wilk`, `final_wilk` |
| Its factors: `L` with `-1` below the diagonal, `u_in = 2^{i-1}`, `P = I` | `wilk_factors` |
| **`ρ_n = 2^{n-1}`** | `maxAbs_wilk`, `activeMax_wilk`, `growth_wilk` |
| The same computation in rounded arithmetic | `ExactOnW.exact`, `ExactOnW.rounded` |

### Optional section: a right-hand side with a large backward error — [`LargeError.lean`](src/LargeError.lean)

| Statement | Lean |
| --- | --- |
| `δ_n = u 2^{n-3}`; both nonzero entries of `b` are floats; the computed factors are exact | `deltaN_eq`, `isFloat_b`, `wilk_factors`, `ExactOnW.rounded` |
| Forward loop: the partial sums are powers of two, and `ŷ₁ = 1`, `ŷᵢ = 2^{i-2}` | `accum_fwd`, `fwd_wilk` |
| `fl(2^{n-2} + δ_n) = 2^{n-2}` (`δ_n` is a quarter of the spacing) | `round_last`, `round_eq_of_near` |
| Back substitution is exact: `x̂₁ = x̂ₙ = 1/2`, `x̂ᵢ = 0` otherwise | `accum_bwd`, `bwd_wilk` |
| End to end: the routine's factors, `P b`, and the loops of the notes give `ŷ` and `x̂` | `solution_wilk`, `fwdLoop_congr`, `bwdLoop_congr` |
| `Ax̂ = e₁`, `r = δ_n eₙ`; the residual is also computed exactly | `wilk_mul_xhat`, `residual_wilk`, `residual_computed_exact` |
| **`η(x̂) = δ_n/(‖A‖₂/√2 + √(1+δ_n²))`** (`‖x̂‖ = 1/√2`, `‖b‖ = √(1+δ_n²)`) | `eta_wilk` |
| **`⌊n/2⌋ ≤ ‖A‖₂ ≤ ‖A‖_F ≤ n`** | `half_le_norm2_wilk`, `norm2_le_frob`, `frob_wilk_le`, `norm2_wilk_bounds` |
| The bottom-left block of size `⌊n/2⌋` consists of `-1`s, and its 2-norm is `⌊n/2⌋` | `wilk_block`, `norm2_neg_ones` |
| Two-sided bounds on `η`: it grows like `u2ⁿ/n` while `u2ⁿ ≪ n`, approaches 1 when `δ_n ≫ n`, and never exceeds 1 | `eta_wilk_bounds` |

## What is not proved in Lean

* **Cited results and literature.** The chapter refers to Higham's book, to Carson and
  Higham for the solve bound, and to Higham's discussion of growth factors. The
  backward-error theorem is proved here, so it does not depend on these references.
  Two statements are not formalized:
  * "large growth can occur in matrices from applications";
  * "complex arithmetic has analogous bounds with different constants".
* **Definitions and informal remarks**:
  * backward stability ("of order `u`, allowing for a modest factor");
  * "ill-conditioned" and "well-conditioned";
  * memory use;
  * "pass a copy";
  * "the zero-pivot test is not a test for numerical rank";
  * the summary list at the end of the chapter.

  The approximations `γ_{3n} ≈ 3nu` and `‖x̂‖ ≈ ‖x‖`, and the growth `u2ⁿ/n`, are
  replaced by precise statements: `gamma_asymp`, `relErr_le_of_norm_le`, and
  `eta_wilk_bounds`.
* **Python details.** The input checks, dtypes, and `np.allclose` are not modelled.
  Neither is the order in which NumPy evaluates `@`, which is unspecified. The Lean
  model accumulates inner products from left to right. In the optional example every
  inner product has at most two nonzero terms, and they are exact, so the order does not
  matter there. That is argued here, not proved for all orders.
* **Model of the routine.** The routine skips zero pivots, so in this model it always
  completes. The hypothesis "if LU with partial pivoting completes" is therefore
  automatic. The factorization bound needs only `nu < 1`. The complete solve uses
  `3nu < 1`.
* **The printed experiment.** The backward errors printed by the notes' code cell
  (`n = 10, …, 70`, `u = 2^{-53}`) are not computed in Lean. Lean proves the exact
  formula `eta_wilk` for every `n ≥ 2` and every precision `p ≥ 1`.
  [`numerics/check_examples.py`](numerics/check_examples.py) reruns the cell in binary64
  and checks the values against this formula. It also checks:
  * the factors, `ŷ`, `x̂` and `r`;
  * the norm bounds;
  * the 3 × 3 example and the `ε = 10⁻¹⁸` cancellation.
* **Real matrices only for the norm-based results.** The residual formula, the condition
  number results, and the perturbation theorem are proved for real matrices with the
  2-norm, as in [`../Cholesky`](../Cholesky). The notes' `x̂ᴴ` becomes `x̂ᵀ`. The Banach
  lemma is proved in any complete normed ring with `‖1‖ = 1`. This covers every induced
  norm, for real or complex matrices.
