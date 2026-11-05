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
        /**
         * @brief Chugunov & DeWitt (2009) screening, valid from weak through strong coupling (CD09ScreeningModel).
         *
         * Their eq. (A4) interpolation between the Salpeter weak limit (with the electrons' contribution weighted by
         * their degeneracy) and the strong-coupling fit, eqs. (23)-(25). Suited to dense or degenerate plasmas (low-mass
         * stars, red-giant and white-dwarf cores) where weak screening is not valid.
         */
        CD09,
        /**
         * @brief Salpeter weak screening blended into Chugunov & DeWitt (2009) as coupling grows
         *        (BlendedScreeningModel of WeakScreeningModel and CD09ScreeningModel).
         *
         * Pure Salpeter while its exponent \f$H\f$ is below ScreeningBlendWindow::weakExponent, where weak screening
         * is accurate (e.g. every pp-chain and CNO reaction in the Sun), pure CD09 above
         * ScreeningBlendWindow::strongExponent, and a smooth blend between.
         */
        BLENDED,
    };

    /**
     * @brief Where BlendedScreeningModel hands over from its weak to its strong prescription.
     *
     * The handover is measured by the weak prescription's own exponent \f$H = \ln f_{\rm weak}\f$ (for Salpeter,
     * \f$H = Z_1 Z_2 e^2 / \lambda_D k_B T\f$, the small parameter of the weak-screening expansion): the weak result
     * is used for \f$H \le\f$ weakExponent, the strong one for \f$H \ge\f$ strongExponent, and their exponents are
     * blended smoothly in \f$\ln H\f$ between. strongExponent must not exceed 2, where WeakScreeningModel caps
     * \f$H\f$, so the cap never enters the blend.
     */
    struct ScreeningBlendWindow {
        double weakExponent = 0.5;   ///< Below this weak-screening exponent the weak prescription is used alone
        double strongExponent = 1.5; ///< Above this weak-screening exponent the strong prescription is used alone
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
    std::unique_ptr<ScreeningModel> selectScreeningModel(ScreeningType type, const ScreeningBlendWindow& window = {});
}