import Cholesky

/-! Inspect the axioms used by every theorem in this project.
Only Lean's standard `propext`, `Classical.choice`, and `Quot.sound`
should occur. In particular, `sorryAx` must not occur. -/

#print axioms Cholesky.conditionUpper_scaled
#print axioms Cholesky.conditionUpper_lt_57
#print axioms Cholesky.conditionUpper_fixed_precision
#print axioms Cholesky.conditionUpper_binary64
#print axioms Cholesky.family_lower_bound
#print axioms Cholesky.negative_eigenvalue_margin
#print axioms Cholesky.update_budget
#print axioms Cholesky.perturbation_fraction
#print axioms Cholesky.breakdown_scalar_contradiction
#print axioms Cholesky.sqrt_step
#print axioms Cholesky.sum_sqrt_le
#print axioms Cholesky.log_chord
#print axioms Cholesky.sum_update_weights
#print axioms Cholesky.meinguet_condition
#print axioms Cholesky.largestQ_properties
#print axioms Cholesky.largestQ_maximal
#print axioms Cholesky.precision_at_least_17
#print axioms Cholesky.negative_direction_requires_error
#print axioms Cholesky.no_positive_perturbation
#print axioms Cholesky.no_gram_factor
#print axioms Cholesky.rankOneShift_quadratic
#print axioms Cholesky.sum_square_bound
#print axioms Cholesky.rankOneShift_posDef
#print axioms Cholesky.exactSchur_posDef
#print axioms Cholesky.rankOneShift_ones
#print axioms Cholesky.roundedSchur_ones
#print axioms Cholesky.roundedSchur_not_posSemidef
#print axioms Cholesky.liftedMatrix_quadratic
#print axioms Cholesky.liftedMatrix_posDef
#print axioms Cholesky.constructedMatrix_posDef
