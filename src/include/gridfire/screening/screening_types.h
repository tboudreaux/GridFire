#pragma once

#include "gridfire/screening/screening_abstract.h"

#include <memory>

namespace gridfire::screening {
    /**
     * @enum ScreeningType
     * @brief Enumerates the available plasma screening models.
     *
     * This enum provides a set of identifiers for the different screening
     * prescriptions that can be used in the reaction rate calculations.
     */
    enum class ScreeningType {
        BARE,    ///< No screening applied. The screening factor is always 1.0.
        /**
         * @brief Weak screening model (Salpeter, 1954).
         *
         * This model is suitable for non-degenerate, non-relativistic plasmas
         * where the electrostatic potential energy between ions is small compared
         * to their thermal kinetic energy. The screening enhancement factor is
         * calculated as `exp(H_12)`.
         *
         * @b Algorithm
         * 1. A composition-dependent term, `ζ = ∑(Z_i^2 + Z_i) * Y_i`, is calculated,
         *    where Z_i is the charge and Y_i is the molar abundance of each species.
         * 2. A prefactor is computed: `prefactor = 0.188 * sqrt(ρ / T₇³) * sqrt(ζ)`,
         *    where ρ is the density and T₇ is the temperature in 10^7 K.
         * 3. For a reaction between two nuclei with charges Z₁ and Z₂, the enhancement
         *    term is `H_12 = prefactor * Z₁ * Z₂`.
         * 4. The final screening factor is `exp(H_12)`.
         * A special calculation is performed for the triple-alpha reaction.
         */
        WEAK,   
    };

    /**
     * @brief A factory function to select and create a screening model.
     *
     * This function returns a `std::unique_ptr` to a concrete implementation of
     * the `ScreeningModel` abstract base class, based on the specified `ScreeningType`.
     * This allows for easy switching between different screening prescriptions at runtime.
     *
     * @param type The `ScreeningType` enum value specifying which model to create.
     * @return A `std::unique_ptr<ScreeningModel>` holding an instance of the
     *         requested screening model.
     *
     * @b Algorithm
     * The function uses a `switch` statement to determine which concrete model to
     * instantiate. If the provided `type` does not match a known case, it defaults
     * to creating a `BareScreeningModel` to ensure safe behavior.
     *
     * @b Post-conditions
     * - A non-null `std::unique_ptr<ScreeningModel>` is always returned.
     *
     * @b Usage
     * @code
     * // Select the weak screening model
     * auto screening_model = gridfire::screening::selectScreeningModel(gridfire::screening::ScreeningType::WEAK);
     *
     * // Use the model to calculate screening factors
     * // (assuming other parameters are initialized)
     * std::vector<double> factors = screening_model->calculateScreeningFactors(
     *     reactions, species, Y, T9, rho
     * );
     * @endcode
     */
    std::unique_ptr<ScreeningModel> selectScreeningModel(ScreeningType type);
}