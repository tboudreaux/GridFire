#include "gridfire/engine/diagnostics/dynamic_engine_diagnostics.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/utils/table_format.h"
#include "fourdst/atomic/species.h"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

namespace gridfire::diagnostics {
    void report_limiting_species(
        const DynamicEngine& engine,
        const std::vector<double>& Y_full,
        const std::vector<double>& E_full,
        const std::vector<double>& dydt_full,
        const double relTol,
        const double absTol,
        const size_t top_n
    ) {
        struct SpeciesError {
            std::string name;
            double ratio;
            double abundance;
            double dydt;
        };

        const auto& species_list = engine.getNetworkSpecies();
        std::vector<SpeciesError> errors;

        for (size_t i = 0; i < species_list.size(); ++i) {
            const double weight = relTol * std::abs(Y_full[i]) + absTol;
            if (weight > 1e-99) { // Avoid division by zero for zero-abundance species
                const double ratio = std::abs(E_full[i]) / weight;
                errors.push_back({
                    std::string(species_list[i].name()),
                    ratio,
                    Y_full[i],
                    dydt_full[i]
                });
            }
        }

        // Sort by error ratio in descending order
        std::ranges::sort(
            errors,
            [](const auto& a, const auto& b) {
                return a.ratio > b.ratio;
            }
        );

        std::vector<std::string> sorted_speciesNames;
        std::vector<double> sorted_err_ratios;
        std::vector<double> sorted_abundances;
        std::vector<double> sorted_dydt;

        for (size_t i = 0; i < std::min(top_n, errors.size()); ++i) {
            sorted_speciesNames.push_back(errors[i].name);
            sorted_err_ratios.push_back(errors[i].ratio);
            sorted_abundances.push_back(errors[i].abundance);
            sorted_dydt.push_back(errors[i].dydt);
        }

        std::vector<std::unique_ptr<utils::ColumnBase>> columns;
        columns.push_back(std::make_unique<utils::Column<std::string>>("Species", sorted_speciesNames));
        columns.push_back(std::make_unique<utils::Column<double>>("Error Ratio", sorted_err_ratios));
        columns.push_back(std::make_unique<utils::Column<double>>("Abundance", sorted_abundances));
        columns.push_back(std::make_unique<utils::Column<double>>("dY/dt", sorted_dydt));

        std::cout << utils::format_table("Timestep Limiting Species", columns) << std::endl;
    }

    void inspect_species_balance(
        const DynamicEngine& engine,
        const std::string& species_name,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) {
        const auto& species_obj = fourdst::atomic::species.at(species_name);

        std::vector<std::string> creation_ids, destruction_ids;
        std::vector<int> creation_stoichiometry, destruction_stoichiometry;
        std::vector<double> creation_flows, destruction_flows;
        double total_creation_flow = 0.0;
        double total_destruction_flow = 0.0;

        for (const auto& reaction : engine.getNetworkReactions()) {
            const int stoichiometry = reaction->stoichiometry(species_obj);
            if (stoichiometry == 0) continue;

            const double flow = engine.calculateMolarReactionFlow(*reaction, comp, T9, rho);

            if (stoichiometry > 0) {
                creation_ids.emplace_back(reaction->id());
                creation_stoichiometry.push_back(stoichiometry);
                creation_flows.push_back(flow);
                total_creation_flow += stoichiometry * flow;
            } else {
                destruction_ids.emplace_back(reaction->id());
                destruction_stoichiometry.push_back(stoichiometry);
                destruction_flows.push_back(flow);
                total_destruction_flow += std::abs(stoichiometry) * flow;
            }
        }

        {
            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<std::string>>("Reaction ID", creation_ids));
            columns.push_back(std::make_unique<utils::Column<int>>("Stoichiometry", creation_stoichiometry));
            columns.push_back(std::make_unique<utils::Column<double>>("Molar Flow", creation_flows));
            std::cout << utils::format_table("Creation Reactions for " + species_name, columns) << std::endl;
        }

        {
            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<std::string>>("Reaction ID", destruction_ids));
            columns.push_back(std::make_unique<utils::Column<int>>("Stoichiometry", destruction_stoichiometry));
            columns.push_back(std::make_unique<utils::Column<double>>("Molar Flow", destruction_flows));
            std::cout << utils::format_table("Destruction Reactions for " + species_name, columns) << std::endl;
        }

        std::cout << "--- Balance Summary for " << species_name << " ---" << std::endl;
        std::cout << "  Total Creation Rate:  " << std::scientific << total_creation_flow << " [mol/g/s]" << std::endl;
        std::cout << "  Total Destruction Rate: " << std::scientific << total_destruction_flow << " [mol/g/s]" << std::endl;
        std::cout << "  Net dY/dt:              " << std::scientific << (total_creation_flow - total_destruction_flow) << std::endl;
        std::cout << "-----------------------------------" << std::endl;
    }

    void inspect_jacobian_stiffness(
        const DynamicEngine& engine,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) {
        const NetworkJacobian jac = engine.generateJacobianMatrix(comp, T9, rho);
        const auto& species_list = engine.getNetworkSpecies();

        double max_diag = 0.0;
        double max_off_diag = 0.0;
        std::optional<fourdst::atomic::Species> max_diag_species = std::nullopt;
        std::optional<std::pair<fourdst::atomic::Species, fourdst::atomic::Species>> max_off_diag_species = std::nullopt;

        for (const auto& rowSpecies : species_list) {
            for (const auto& colSpecies : species_list) {
                const double val = std::abs(jac(rowSpecies, colSpecies));
                if (rowSpecies == colSpecies) {
                    if (val > max_diag) { max_diag = val; max_diag_species = colSpecies; }
                } else {
                    if (val > max_off_diag) { max_off_diag = val; max_off_diag_species = {rowSpecies, colSpecies};}
                }
            }
        }

        std::cout << "\n--- Jacobian Stiffness Report ---" << std::endl;
        if (max_diag_species.has_value()) {
            std::cout << "  Largest Diagonal Element (d(dYi/dt)/dYi): " << std::scientific << max_diag
                      << " for species " << max_diag_species->name() << std::endl;
        }
        if (max_off_diag_species.has_value()) {
            std::cout << "  Largest Off-Diagonal Element (d(dYi/dt)/dYj): " << std::scientific << max_off_diag
                      << " for d(" << max_off_diag_species->first.name()
                      << ")/d(" << max_off_diag_species->second.name() << ")" << std::endl;
        }
        std::cout << "---------------------------------" << std::endl;
    }

}