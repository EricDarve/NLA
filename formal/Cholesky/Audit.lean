import Cholesky

/-! Inspect the axioms used by the theorems of this project.
Only Lean's standard `propext`, `Classical.choice`, and `Quot.sound`
should occur. In particular, `sorryAx` must not occur. -/

-- Symmetric positive definite matrices
#print axioms Cholesky.posDef_iff_real
#print axioms Cholesky.quad_eq_sum_eigenvalues
#print axioms Cholesky.posDef_iff_eigenvalues_pos
#print axioms Cholesky.posDef_isUnit
#print axioms Cholesky.posDef_diag_pos
#print axioms Cholesky.example12_diag_pos
#print axioms Cholesky.example12_eigenvalues
#print axioms Cholesky.example12_not_posDef
#print axioms Cholesky.Herm.cholesky_iff
#print axioms Cholesky.Herm.cholesky_unique
-- Deriving the factorization
#print axioms Cholesky.mul_transpose_succ_succ
#print axioms Cholesky.block_equations
#print axioms Cholesky.complete_square
#print axioms Cholesky.schur_transpose
#print axioms Cholesky.schur_posDef
-- Existence and uniqueness
#print axioms Cholesky.cholesky_iff
#print axioms Cholesky.cholesky_existsUnique
#print axioms Cholesky.cholesky_unique
#print axioms Cholesky.posDef_of_factorization
#print axioms Cholesky.psdExample_posSemidef
#print axioms Cholesky.psdExample_zero_pivot
-- Connection with LU
#print axioms Cholesky.ldlL0_unit_lower
#print axioms Cholesky.ldl_eq
#print axioms Cholesky.lu_eq
#print axioms Cholesky.ldlU_upper
#print axioms Cholesky.ldlU_diag
#print axioms Cholesky.lu_unique
#print axioms Cholesky.cholesky_is_lu
#print axioms Cholesky.exact_pivot_eq_sq
-- The in-place algorithm
#print axioms Cholesky.exact_succeeds_iff
#print axioms Cholesky.exact_factor_spec
#print axioms Cholesky.exact_factor_mul_transpose
#print axioms Cholesky.factor_symmOfLower
#print axioms Cholesky.factor_lowerEq
#print axioms Cholesky.run_upper
-- Arithmetic and storage
#print axioms Cholesky.card_stepAt_updates
#print axioms Cholesky.two_mul_updEntries
#print axioms Cholesky.sum_update_flops
#print axioms Cholesky.cholFlops_split
#print axioms Cholesky.cholFlops_ratio
#print axioms Cholesky.luFlops_real
#print axioms Cholesky.chol_lu_ratio
#print axioms Cholesky.fullUpdateFlops_ratio
#print axioms Cholesky.card_lower_triangle
-- Solving a linear system
#print axioms Cholesky.exact_solve
#print axioms Cholesky.triSolveFlops_eq
#print axioms Cholesky.two_solves_flops
-- Factor entries
#print axioms Cholesky.diag_eq_sum_sq_le
#print axioms Cholesky.abs_entry_le_sqrt_diag
#print axioms Cholesky.absProd_le
#print axioms Cholesky.exact_absProd_bound
-- Floating-point arithmetic
#print axioms Cholesky.exists_roundNearest
#print axioms Cholesky.roundNearest_relErr
#print axioms Cholesky.stdModel_of_roundNearest
-- Backward error of the factorization
#print axioms Cholesky.backward_error
#print axioms Cholesky.diag_equation
#print axioms Cholesky.frob_sq_le
#print axioms Cholesky.norm2_le_frob
#print axioms Cholesky.trace_le
#print axioms Cholesky.norm2_E_le
#print axioms Cholesky.norm2_E_le_rel
#print axioms Cholesky.nearby_spd
#print axioms Cholesky.stability_constant_eq
#print axioms Cholesky.stability_constant_asymp
-- Backward and forward error of the computed solution
#print axioms Cholesky.fwdSolve_backward
#print axioms Cholesky.bwdSolve_backward
#print axioms Cholesky.solve_backward_error
#print axioms Cholesky.gamma_solve_le
#print axioms Cholesky.rigal_gaches
#print axioms Cholesky.eta_eq
#print axioms Cholesky.eta_le_of_perturbation
#print axioms Cholesky.eta_solve_le
#print axioms Cholesky.beta_asymp
#print axioms Cholesky.forward_error
#print axioms Cholesky.cond2_spd
-- Guarantees of successful completion (the 1.1 corollary of Meinguet's condition)
#print axioms Cholesky.meinguet_prod_le
#print axioms Cholesky.sum_update_weights_steps
#print axioms Cholesky.sum_update_weights
#print axioms Cholesky.meinguet_constant
#print axioms Cholesky.log_chord
#print axioms Cholesky.meinguet_condition
-- The counterexample family
#print axioms Cholesky.CE.admissible_iff
#print axioms Cholesky.CE.size_iff
#print axioms Cholesky.hadamard_mul_transpose
#print axioms Cholesky.hadamard_entry
#print axioms Cholesky.hadQ_transpose_mul
#print axioms Cholesky.hadQ_entry
#print axioms Cholesky.CE.Par.t_sq_eq
#print axioms Cholesky.CE.Par.t_float
#print axioms Cholesky.CE.Par.tt_float
#print axioms Cholesky.CE.Par.T_float
#print axioms Cholesky.CE.Par.A_float
#print axioms Cholesky.CE.Par.spacing_facts
#print axioms Cholesky.CE.Par.T_gram
#print axioms Cholesky.CE.Par.C_eq
#print axioms Cholesky.CE.Par.schur_eq
#print axioms Cholesky.CE.Par.S_quad_bounds
#print axioms Cholesky.CE.Par.S_ones
#print axioms Cholesky.CE.Par.S_perp
#print axioms Cholesky.CE.Par.A_posDef
#print axioms Cholesky.lead_run
#print axioms Cholesky.lead_pivot
#print axioms Cholesky.run_natAdd
#print axioms Cholesky.CE.Par.round_diag
#print axioms Cholesky.CE.Par.round_cross
#print axioms Cholesky.CE.Par.phase1
#print axioms Cholesky.CE.Par.phase2
#print axioms Cholesky.CE.Par.trailing_eq_Shat
#print axioms Cholesky.CE.Par.Shat_ones
#print axioms Cholesky.CE.Par.Shat_alt
#print axioms Cholesky.CE.Par.Shat_quad_bounds
#print axioms Cholesky.CE.Par.norm2_Shat_le
#print axioms Cholesky.negative_eigenvalue_margin
#print axioms Cholesky.factor_offdiag_nonpos
#print axioms Cholesky.norm2_absMat_le_three
#print axioms Cholesky.perturbation_fraction
#print axioms Cholesky.update_budget
#print axioms Cholesky.no_success
#print axioms Cholesky.CE.Par.A_fails
-- Condition number and the remaining gap
#print axioms Cholesky.CE.Par.norm2_T_sq_le
#print axioms Cholesky.CE.Par.norm2_T_sq_ge
#print axioms Cholesky.CE.Par.norm2_A_le
#print axioms Cholesky.CE.Par.norm2_Ainv_le
#print axioms Cholesky.CE.Par.cond_A_le
#print axioms Cholesky.CE.Par.cond_A_lt
#print axioms Cholesky.conditionUpper_lt_57
#print axioms Cholesky.constant_57
#print axioms Cholesky.CE.Par.norm2_A_ge
#print axioms Cholesky.CE.Par.norm2_Ainv_ge
#print axioms Cholesky.CE.Par.cond_A_ge
#print axioms Cholesky.CE.Par.cond_A_ge'
#print axioms Cholesky.CE.cStar_admissible
-- Bounded condition numbers at fixed precision; padding
#print axioms Cholesky.largestQ_properties
#print axioms Cholesky.largestQ_maximal
#print axioms Cholesky.CE.fixed_precision
#print axioms Cholesky.conditionUpper_fixed_precision
#print axioms Cholesky.conditionUpper_expand
#print axioms Cholesky.conditionUpper_largest
#print axioms Cholesky.fails_of_leading
#print axioms Cholesky.cond2_blkdiag
#print axioms Cholesky.CE.Par.B_fails
#print axioms Cholesky.CE.Par.cond2_B
#print axioms Cholesky.CE.Par.B_posDef
#print axioms Cholesky.CE.Par.B_float
#print axioms Cholesky.CE.binary64
#print axioms Cholesky.CE.binary32
#print axioms Cholesky.conditionUpper_binary64
-- Sharpness
#print axioms Cholesky.succeeds_of_cond_lt
#print axioms Cholesky.cond2_ge_one
#print axioms Cholesky.cStar_ge_u
#print axioms Cholesky.cStar_ge_meinguet_mu
#print axioms Cholesky.cStar_ge_meinguet
#print axioms Cholesky.CE.sharp_family
#print axioms Cholesky.CE.cStar_lt
#print axioms Cholesky.CE.cStar_bounds
-- Operation order, condition numbers near one, numerical examples
#print axioms Cholesky.order_matters
#print axioms Cholesky.CE.Par.cond_A_ge_385
#print axioms Cholesky.family_lower_bound
#print axioms Cholesky.table_parameters
#print axioms Cholesky.table_bounds
#print axioms Cholesky.table_ratio
#print axioms Cholesky.memory_sizes
-- Failure is not monotone in precision
#print axioms Cholesky.nmA_posDef
#print axioms Cholesky.nmA_float
#print axioms Cholesky.nonmonotone_p3
#print axioms Cholesky.nonmonotone_p4
-- Further reading
#print axioms Cholesky.scaled_unit_diag
#print axioms Cholesky.scaled_diag_eq_one
#print axioms Cholesky.cond2_perm
-- Earlier lemmas (kept for reference)
#print axioms Cholesky.exactSchur_posDef
#print axioms Cholesky.roundedSchur_not_posSemidef
#print axioms Cholesky.liftedMatrix_posDef
#print axioms Cholesky.no_gram_factor
#print axioms Cholesky.breakdown_scalar_contradiction
