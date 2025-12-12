/**
 * @file dynamic_engine_diagnostics.h
 * @brief Diagnostics utilities for DynamicEngine instances.
 *
 * These free functions provide diagnostic tools for analyzing and reporting
 * on the behavior of DynamicEngine instances in GridFire. They include
 * functions for identifying limiting species, inspecting species balance,
 * and evaluating Jacobian stiffness. Generally these functions are expensive to call
 * and should be avoided in performance-critical code paths or during large simulations.
 * These are primarily intended for debugging and analysis during development.
 *
 * @section Functions
 *  - report_limiting_species: Identifies species that are limiting the
 *    accuracy of the integration based on specified tolerances.
 *  - inspect_species_balance: Analyzes the production and destruction rates
 *    of a specified species to identify imbalances.
 *  - inspect_jacobian_stiffness: Evaluates the stiffness of the Jacobian
 *    matrix to identify potential numerical issues.
 *  - save_jacobian_to_file: Saves the Jacobian matrix to a file for
 *    external analysis.
 *
 * @note Each of these functions, aside from save_jacobian_to_file, may return results in JSON format for easy
 * integration with other tools and workflows. If the `json` parameter is set
 * then the output will be formatted as JSON; otherwise, it will be printed
 * to standard output and `std::nullopt` will be returned.
 */


#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/scratchpads/blob.h"

#include <vector>
#include <string>

#include "nlohmann/json.hpp"

namespace gridfire::engine::diagnostics {
    /**
     * @brief Report the species that are limiting the accuracy of the integration. This is useful for identifying
     * species that may be limiting the timestepping due to their high relative errors compared to the specified tolerances.
     * @param engine Constant reference to the DynamicEngine instance.
     * @param Y_full Vector of the current species molar abundances sorted in the same order as the dynamic engine species list.
     * @param E_full Vector of the current species molar abundance errors sorted in the same order as the dynamic engine species list.
     * @param relTol Relative tolerance for the integration.
     * @param absTol Absolute tolerance for the integration.
     * @param top_n The number of top limiting species to report. Default is 10.
     * @param json Flag indicating whether to return the results in JSON format. If false, results are printed to standard output.
     * @return std::optional<nlohmann::json> JSON object containing the limiting species report if `json` is true; otherwise, std::nullopt.
     */
    std::optional<nlohmann::json> report_limiting_species(
        scratch::StateBlob& ctx,
        const DynamicEngine &engine,
        const std::vector<double> &Y_full,
        const std::vector<double> &E_full,
        double relTol,
        double absTol,
        size_t top_n = 10,
        bool json = false
    );


    /**
     * @brief Inspect the production and destruction balance of a specific species in the reaction network. This function analyzes the reactions that create and destroy the specified species,
     * providing insights into potential imbalances that may affect the accuracy of the simulation.
     * @param engine Constant reference to the DynamicEngine instance.
     * @param species_name The name of the species to inspect.
     * @param comp The current composition of the system as a fourdst::composition::Composition object.
     * @param T9 The temperature in GK (10^9 K).
     * @param rho The density in g/cm^3.
     * @param json Flag indicating whether to return the results in JSON format. If false, results are printed to standard output.
     * @return std::optional<nlohmann::json> JSON object containing the species balance report if `json` is true; otherwise, std::nullopt.
     */
    std::optional<nlohmann::json> inspect_species_balance(
        scratch::StateBlob& ctx,
        const DynamicEngine& engine,
        const std::string& species_name,
        const fourdst::composition::Composition &comp,
        double T9,
        double rho,
        bool json
    );

    /**
     * @brief Inspect the stiffness of the Jacobian matrix for the reaction network at the given temperature and density. This function evaluates the diagonal and off-diagonal elements of the Jacobian to identify potential numerical issues related to stiffness.
     * @param engine Constant reference to the DynamicEngine instance.
     * @param comp The current composition of the system as a fourdst::composition::Composition object.
     * @param T9 The temperature in GK (10^9 K).
     * @param rho The density in g/cm^3.
     * @param json Flag indicating whether to return the results in JSON format. If false, results are printed to standard output.
     * @return std::optional<nlohmann::json> JSON object containing the Jacobian stiffness report if `json` is true; otherwise, std::nullopt.
     */
    std::optional<nlohmann::json> inspect_jacobian_stiffness(
        scratch::StateBlob& ctx,
        const DynamicEngine &engine,
        const fourdst::composition::Composition &comp,
        double T9,
        double rho,
        bool json
    );

    /**
     * @brief Save the Jacobian matrix to a file for external analysis.
     * @param jacobian Constant reference to the NetworkJacobian instance to be saved.
     * @param filename The name of the file where the Jacobian will be saved.
     */
    void save_jacobian_to_file(
        const NetworkJacobian& jacobian,
        const std::string& filename
    );

}
