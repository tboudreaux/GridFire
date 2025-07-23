#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "cppad/cppad.hpp"

namespace gridfire::screening {
    /**
     * @class BareScreeningModel
     * @brief A screening model that applies no screening effect.
     *
     * This class implements the `ScreeningModel` interface but returns a
     * screening factor of 1.0 for all reactions, regardless of the plasma
     * conditions. It represents the case of bare, unscreened nuclei and serves
     * as a baseline or can be used when screening effects are negligible or
     * intentionally ignored.
     *
     * @implements ScreeningModel
     */
    class BareScreeningModel final : public ScreeningModel {
    /// @brief Alias for CppAD Automatic Differentiation type for double precision.
    using ADDouble = CppAD::AD<double>;
    public:
        /**
         * @brief Calculates screening factors, which are always 1.0.
         *
         * This implementation returns a vector of screening factors where every
         * element is 1.0, effectively applying no screening correction to the
         * reaction rates.
         *
         * @param reactions The set of logical reactions in the network.
         * @param species A vector of all atomic species (unused).
         * @param Y A vector of the molar abundances (unused).
         * @param T9 The temperature (unused).
         * @param rho The plasma density (unused).
         * @return A vector of doubles, with each element being 1.0, of the same
         *         size as the `reactions` set.
         *
         * @b Algorithm
         * The function simply creates and returns a `std::vector<double>` of the
         * same size as the input `reactions` set, with all elements initialized to 1.0.
         *
         * @b Usage
         * @code
         * BareScreeningModel bare_model;
         * // ... (initialize reactions, species, Y, T9, rho)
         * std::vector<double> factors = bare_model.calculateScreeningFactors(
         *     reactions, species, Y, T9, rho
         * );
         * // 'factors' will contain [1.0, 1.0, ...]
         * @endcode
         */
        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const override;

        /**
         * @brief Calculates screening factors for AD types, which are always 1.0.
         *
         * This implementation returns a vector of AD-typed screening factors where
         * every element is 1.0. This is the automatic differentiation-compatible
         * version.
         *
         * @param reactions The set of logical reactions in the network.
         * @param species A vector of all atomic species (unused).
         * @param Y A vector of the molar abundances as AD types (unused).
         * @param T9 The temperature as an AD type (unused).
         * @param rho The plasma density as an AD type (unused).
         * @return A vector of ADDouble, with each element being 1.0, of the same
         *         size as the `reactions` set.
         */
        [[nodiscard]] std::vector<ADDouble> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<ADDouble>& Y,
            const ADDouble T9,
            const ADDouble rho
        ) const override;
    private:
        /**
         * @brief Template implementation for calculating screening factors.
         *
         * This private helper function contains the core logic for both the `double`
         * and `ADDouble` versions of `calculateScreeningFactors`. It is templated
         * to handle both numeric types seamlessly.
         *
         * @tparam T The numeric type, either `double` or `CppAD::AD<double>`.
         * @param reactions The set of reactions for which to calculate factors.
         * @param species A vector of all atomic species (unused).
         * @param Y A vector of molar abundances (unused).
         * @param T9 The temperature (unused).
         * @param rho The density (unused).
         * @return A vector of type `T` with all elements initialized to 1.0.
         */
        template <typename T>
        [[nodiscard]] std::vector<T> calculateFactors_impl(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            const T T9,
            const T rho
        ) const;
    };

    /**
     * @brief Template implementation for the bare screening model.
     *
     * This function provides the actual implementation for `calculateFactors_impl`.
     * It creates a vector of the appropriate numeric type (`T`) and size, and
     * initializes all its elements to 1.0, representing no screening.
     *
     * @tparam T The numeric type, either `double` or `CppAD::AD<double>`.
     * @param reactions The set of reactions, used to determine the size of the output vector.
     * @param species Unused parameter.
     * @param Y Unused parameter.
     * @param T9 Unused parameter.
     * @param rho Unused parameter.
     * @return A `std::vector<T>` of the same size as `reactions`, with all elements set to 1.0.
     */
    template<typename T>
    std::vector<T> BareScreeningModel::calculateFactors_impl(
        const reaction::LogicalReactionSet &reactions,
        const std::vector<fourdst::atomic::Species> &species,
        const std::vector<T> &Y,
        const T T9,
        const T rho
    ) const {
        return std::vector<T>(reactions.size(), T(1.0)); // Bare screening returns 1.0 for all reactions
    }
}
