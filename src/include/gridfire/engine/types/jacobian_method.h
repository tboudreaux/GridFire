/**
 * @file jacobian_method.h
 * @brief Selects how GraphEngine computes the species Jacobian.
 */
#pragma once

namespace gridfire::engine {
    /**
     * @enum JacobianMethod
     * @brief How GraphEngine computes the species Jacobian dYdot_i/dY_j.
     *
     * Both methods differentiate the same right-hand side and agree to round-off, including the derivatives with
     * respect to species that are currently absent (Y_j = 0). They differ only in cost. Set through
     * GridFireConfig::engine.jacobian.
     */
    enum class JacobianMethod {
        /// Assembled directly from the precomputed reaction topology (reactant powers, symmetry factors,
        /// stoichiometry) and the screening model's analytic derivatives. Tens of times faster than the tape.
        /// Falls back to AUTOMATIC_DIFFERENTIATION for a screening model without analytic derivatives.
        ANALYTIC,
        /// Dense reverse-mode sweeps over the CppAD tape of calculateAllDerivatives. Kept as the reference
        /// implementation and as the fallback for screening models without analytic derivatives.
        AUTOMATIC_DIFFERENTIATION
    };
}
