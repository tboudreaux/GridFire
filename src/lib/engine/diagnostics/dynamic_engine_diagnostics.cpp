#include "gridfire/engine/diagnostics/dynamic_engine_diagnostics.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/utils/table_format.h"
#include "fourdst/atomic/species.h"

#include <vector>
#include <string>
#include <algorithm>

namespace gridfire::engine::diagnostics {
    std::optional<nlohmann::json> report_limiting_species(
        const DynamicEngine &engine,
        const std::vector<double> &Y_full,
        const std::vector<double> &E_full,
        const double relTol,
        const double absTol,
        const size_t top_n,
        const bool json
    ) {
        struct SpeciesError {
            std::string name;
            double error;
            double ratio;
            double abundance;
        };

        const auto& species_list = engine.getNetworkSpecies();
        std::vector<SpeciesError> errors;

        for (size_t i = 0; i < species_list.size(); ++i) {
            const double weight = relTol * std::abs(Y_full[i]) + absTol;
            if (weight > 1e-99) { // Avoid division by zero for zero-abundance species
                const double ratio = std::abs(E_full[i]) / weight;
                errors.push_back({
                    std::string(species_list[i].name()),
                    E_full[i],
                    ratio,
                    Y_full[i],
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
        std::vector<double> sorted_errors;

        for (size_t i = 0; i < std::min(top_n, errors.size()); ++i) {
            sorted_speciesNames.push_back(errors[i].name);
            sorted_err_ratios.push_back(errors[i].ratio);
            sorted_abundances.push_back(errors[i].abundance);
            sorted_errors.push_back(errors[i].error);
        }

        std::vector<std::unique_ptr<utils::ColumnBase>> columns;
        columns.push_back(std::make_unique<utils::Column<std::string>>("Species", sorted_speciesNames));
        columns.push_back(std::make_unique<utils::Column<double>>("Error Ratio", sorted_err_ratios));
        columns.push_back(std::make_unique<utils::Column<double>>("Abundance", sorted_abundances));
        columns.push_back(std::make_unique<utils::Column<double>>("Error", sorted_errors));

        if (json) {
            return utils::to_json(columns);
        }

        utils::print_table("Timestep Limiting Species", columns);
        return std::nullopt;
    }

    std::optional<nlohmann::json> inspect_species_balance(
        const DynamicEngine& engine,
        const std::string& species_name,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho,
        bool json
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

        nlohmann::json j;
        {
            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<std::string>>("Reaction ID", creation_ids));
            columns.push_back(std::make_unique<utils::Column<int>>("Stoichiometry", creation_stoichiometry));
            columns.push_back(std::make_unique<utils::Column<double>>("Molar Flow", creation_flows));
            if (json) {
                j["Creation_Reactions_" + species_name] = utils::to_json(columns);
            }
            else {
                utils::print_table("Creation Reactions for " + species_name, columns);
            }
        }

        {
            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<std::string>>("Reaction ID", destruction_ids));
            columns.push_back(std::make_unique<utils::Column<int>>("Stoichiometry", destruction_stoichiometry));
            columns.push_back(std::make_unique<utils::Column<double>>("Molar Flow", destruction_flows));
            if (json) {
                j["Destruction_Reactions_" + species_name] = utils::to_json(columns);
            } else {
                utils::print_table("Destruction Reactions for " + species_name, columns);
            }
        }

        std::vector<std::unique_ptr<utils::ColumnBase>> summary_columns;
        summary_columns.push_back(std::make_unique<utils::Column<std::string>>("Metric", std::vector<std::string>{
            "Total Creation Rate [mol/g/s]",
            "Total Destruction Rate [mol/g/s]",
            "Net dY/dt [mol/g/s]"
        }));
        summary_columns.push_back(std::make_unique<utils::Column<double>>("Value", std::vector<double>{
            total_creation_flow,
            total_destruction_flow,
            total_creation_flow - total_destruction_flow
        }));

        if (json) {
            j["Species_Balance_Summary_" + species_name] = utils::to_json(summary_columns);
            return j;
        }

        utils::print_table("Species Balance Summary for " + species_name, summary_columns);
        return std::nullopt;
    }

    std::optional<nlohmann::json> inspect_jacobian_stiffness(
        const DynamicEngine &engine,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho,
        const bool json
    ) {
        NetworkJacobian jac = engine.generateJacobianMatrix(comp, T9, rho);

        jac = regularize_jacobian(jac, comp);

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

        std::vector<std::unique_ptr<utils::ColumnBase>> jacobian_columns;
        jacobian_columns.push_back(std::make_unique<utils::Column<std::string>>("Metric", std::vector<std::string>{
            "Largest Diagonal Element (d(dYi/dt)/dYi)",
            "Largest Off-Diagonal Element (d(dYi/dt)/dYj)"
        }));
        jacobian_columns.push_back(std::make_unique<utils::Column<double>>("Value", std::vector<double>{
            max_diag,
            max_off_diag
        }));
        jacobian_columns.push_back(std::make_unique<utils::Column<std::string>>("Species", std::vector<std::string>{
            max_diag_species.has_value() ? std::string(max_diag_species->name()) : "N/A",
            max_off_diag_species.has_value() ?
                ("d(" + std::string(max_off_diag_species->first.name()) + ")/d(" + std::string(max_off_diag_species->second.name()) + ")")
                : "N/A"
        }));

        if (json) {
            nlohmann::json j;
            j["Jacobian_Stiffness"] = utils::to_json(jacobian_columns);
            return j;
        }
        utils::print_table("Jacobian Stiffness Diagnostics", jacobian_columns);
        return std::nullopt;
    }
}