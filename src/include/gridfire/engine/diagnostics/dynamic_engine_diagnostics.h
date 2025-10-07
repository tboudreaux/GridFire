#pragma once

#include "gridfire/engine/engine_abstract.h"

#include <vector>
#include <string>

namespace gridfire::diagnostics {
    void report_limiting_species(
        const DynamicEngine& engine,
        const std::vector<double>& Y_full,
        const std::vector<double>& E_full,
        const std::vector<double>& dydt_full,
        double relTol,
        double absTol,
        size_t top_n = 10
    );

    void inspect_species_balance(
        const DynamicEngine& engine,
        const std::string& species_name,
        const fourdst::composition::Composition &comp,
        double T9,
        double rho
    );

    void inspect_jacobian_stiffness(
        const DynamicEngine& engine,
        const fourdst::composition::Composition &comp,
        double T9,
        double rho
    );

}