#pragma once

#include "gridfire/reaction/reaction.h"

#include "fourdst/composition/atomicSpecies.h"

#include "cppad/cppad.hpp"

#include <vector>

namespace gridfire::screening {
    /**
     * @class ScreeningModel
     * @brief An abstract base class for plasma screening models.
     *
     * This class defines the interface for models that calculate the enhancement
     * factor for nuclear reaction rates due to the electrostatic screening of
     * interacting nuclei by the surrounding plasma. Concrete implementations of
     * this class will provide specific screening prescriptions (e.g., WEAK,
     * BARE, STRONG, etc.).
     *
     * The interface provides methods for calculating screening factors for both
     * standard double-precision inputs and for CppAD's automatic differentiation
     * types, allowing the screening contributions to be included in Jacobian
     * calculations.
     */
    class ScreeningModel {
    public:
        /// @brief Alias for CppAD Automatic Differentiation type for double precision.
        using ADDouble = CppAD::AD<double>;
        /**
         * @brief Virtual destructor.
         *
         * Ensures that derived class destructors are called correctly.
         */
        virtual ~ScreeningModel() = default;

        /**
         * @brief Calculates screening factors for a set of reactions.
         *
         * This is a pure virtual function that must be implemented by derived
         * classes. It computes the screening enhancement factor for each reaction
         * in the provided set based on the given plasma conditions.
         *
         * @param reactions The set of logical reactions in the network.
         * @param species A vector of all atomic species involved in the network.
         * @param Y A vector of the molar abundances (mol/g) for each species.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The plasma density in g/cm^3.
         * @return A vector of screening factors (dimensionless), one for each reaction
         *         in the `reactions` set, in the same order.
         *
         * @b Pre-conditions
         * - The size of the `Y` vector must match the size of the `species` vector.
         * - `T9` and `rho` must be positive.
         *
         * @b Post-conditions
         * - The returned vector will have the same size as the `reactions` set.
         * - Each element in the returned vector will be >= 1.0.
         *
         * @b Usage
         * @code
         * // Assume 'model' is a std::unique_ptr<ScreeningModel> to a concrete implementation
         * // and other parameters (reactions, species, Y, T9, rho) are initialized.
         * std::vector<double> screening_factors = model->calculateScreeningFactors(
         *     reactions, species, Y, T9, rho
         * );
         * for (size_t i = 0; i < reactions.size(); ++i) {
         *     // ... use screening_factors[i] ...
         * }
         * @endcode
         */
        virtual std::vector<double> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            const double T9,
            const double rho
            ) const = 0;

        /**
         * @brief Calculates screening factors using CppAD types for automatic differentiation.
         *
         * This is a pure virtual function that provides an overload of
         * `calculateScreeningFactors` for use with CppAD. It allows the derivatives
         * of the screening factors with respect to abundances, temperature, and
         * density to be computed automatically.
         *
         * @param reactions The set of logical reactions in the network.
         * @param species A vector of all atomic species involved in the network.
         * @param Y A vector of the molar abundances (mol/g) for each species, as AD types.
         * @param T9 The temperature in units of 10^9 K, as an AD type.
         * @param rho The plasma density in g/cm^3, as an AD type.
         * @return A vector of screening factors (dimensionless), as AD types.
         *
         * @b Note
         * This method is essential for including the effects of screening in the
         * Jacobian matrix of the reaction network.
         */
        virtual std::vector<ADDouble> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<ADDouble>& Y,
            const ADDouble T9,
            const ADDouble rho
        ) const = 0;
    };
}