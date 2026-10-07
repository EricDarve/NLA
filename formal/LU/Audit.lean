import LU

/-! Inspect the axioms used by the theorems of this project.
Only Lean's standard `propext`, `Classical.choice`, and `Quot.sound`
should occur. In particular, `sorryAx` must not occur. -/

-- Why interchange rows?
#print axioms LU.swap_no_lu
#print axioms LU.Aeps_lu
#print axioms LU.Aeps_cancel
#print axioms LU.Aeps_final
#print axioms LU.Aeps_plu
#print axioms LU.Aeps_double
#print axioms LU.Aeps_double_product
-- Partial pivoting, in place
#print axioms LU.exists_firstMax
#print axioms LU.IsFirstMax.unique
#print axioms LU.firstMax_self
#print axioms LU.run_pred
#print axioms LU.rel_final
#print axioms LU.final_col
#print axioms LU.norm_multiplier_le_one
#print axioms LU.norm_factorL_le_one
#print axioms LU.pivot_eq_zero_col
#print axioms LU.pivot_eq_diag
#print axioms LU.pivot_ne_zero
#print axioms LU.exact_plu
#print axioms LU.solve_of_plu
-- Existence of the factorization
#print axioms LU.exact_v_eq
#print axioms LU.exact_mult_of_pivot_zero
#print axioms LU.exact_trail_eq
#print axioms LU.exists_lu_pivoting
#print axioms LU.exists_lu_pivoting_real
#print axioms LU.exists_lu_pivoting_complex
#print axioms LU.exists_diag_eq_zero
#print axioms LU.not_unique_of_singular
-- The 3 × 3 example and the cost
#print axioms LU.A3_piv0
#print axioms LU.A3_piv1
#print axioms LU.A3_final
#print axioms LU.A3_plu
#print axioms LU.A3_triangular
#print axioms LU.A3_solve
#print axioms LU.elim_apply_of_not_mem
#print axioms LU.card_candidates
#print axioms LU.card_multipliers
#print axioms LU.card_updates
#print axioms LU.luFlops_real
#print axioms LU.luFlops_ratio
#print axioms LU.two_mul_pivotCompares
#print axioms LU.overhead_le
#print axioms LU.residualFlops_eq
-- Backward error and the residual
#print axioms LU.tauSet_nonempty
#print axioms LU.eta_le_tau
#print axioms LU.residual_lower
#print axioms LU.attaining_pair
#print axioms LU.norm2_vecMulVec
#print axioms LU.eta_zero_pairs
#print axioms LU.eta_isLeast
#print axioms LU.eta_eq
#print axioms LU.eta_zero
#print axioms LU.etaFormula_le_one
-- Sensitivity
#print axioms LU.hasSensitivity_of_hasFDerivAt
#print axioms LU.hasSensitivity_of_hasDerivAt
#print axioms LU.forward_le_sens
#print axioms LU.relative_le_cond
#print axioms LU.sens_linear
#print axioms LU.relCond_linear_le
-- The matrix condition number
#print axioms LU.one_le_norm_mul_norm_inv
#print axioms LU.one_le_cond2
#print axioms LU.stretch_bounds
#print axioms LU.norm2_eq_sigmaMax
#print axioms LU.norm2_inv_eq
#print axioms LU.cond2_eq_sigma
#print axioms LU.sigmaMax_isGreatest
#print axioms LU.sigmaMin_isLeast
#print axioms LU.diagEps_solve
#print axioms LU.diagEps_relative
#print axioms LU.diagEps_cond
#print axioms LU.norm2_perm_mul
#print axioms LU.gram_perm_mul
#print axioms LU.singVal_perm_mul
#print axioms LU.cond2_perm_mul
-- Forward error
#print axioms LU.error_eq
#print axioms LU.relErr_le
#print axioms LU.relErr_le_of_norm_le
#print axioms LU.neumann_partial
#print axioms LU.banach_inverse
#print axioms LU.banach
#print axioms LU.banach_matrix
#print axioms LU.inv_add_bound
#print axioms LU.perturbation_bound
#print axioms LU.forward_error_bound
-- Backward error of LU
#print axioms LU.roundNearest_relErr
#print axioms LU.stdModel_of_roundNearest
#print axioms LU.exists_roundNearest
#print axioms LU.gamma_add_le
#print axioms LU.gamma_three_le
#print axioms LU.mult_err
#print axioms LU.upd_err
#print axioms LU.factor_backward
#print axioms LU.factor_backward'
#print axioms LU.fwdSolve_backward
#print axioms LU.bwdSolve_backward
#print axioms LU.solve_backward
#print axioms LU.exact_luSolve
#print axioms LU.eta_lu_le
#print axioms LU.gamma_asymp
#print axioms LU.Aeps_absProd_unpivoted
#print axioms LU.Aeps_absProd_pivoted
-- Element growth
#print axioms LU.upd_norm_le
#print axioms LU.stage_bound
#print axioms LU.activeMax_le
#print axioms LU.growth_le
#print axioms LU.ExactOnW.exact
#print axioms LU.ExactOnW.rounded
#print axioms LU.step_wilk
#print axioms LU.run_wilk
#print axioms LU.final_wilk
#print axioms LU.wilk_factors
#print axioms LU.maxAbs_wilk
#print axioms LU.activeMax_wilk
#print axioms LU.growth_wilk
-- A right-hand side that produces a large backward error
#print axioms LU.round_eq_of_near
#print axioms LU.deltaN_eq
#print axioms LU.isFloat_b
#print axioms LU.accum_fwd
#print axioms LU.round_last
#print axioms LU.fwd_wilk
#print axioms LU.accum_bwd
#print axioms LU.bwd_wilk
#print axioms LU.fwdLoop_congr
#print axioms LU.bwdLoop_congr
#print axioms LU.solution_wilk
#print axioms LU.wilk_mul_xhat
#print axioms LU.residual_wilk
#print axioms LU.residual_computed_exact
#print axioms LU.eta_wilk
#print axioms LU.frob_wilk_le
#print axioms LU.half_le_norm2_wilk
#print axioms LU.wilk_block
#print axioms LU.norm2_neg_ones
#print axioms LU.norm2_wilk_bounds
#print axioms LU.eta_wilk_bounds
