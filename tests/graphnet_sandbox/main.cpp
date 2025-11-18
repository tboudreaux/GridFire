#include <iostream>
#include <fstream>

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/solver/strategies/CVODE_solver_strategy.h"
#include "gridfire/policy/stellar_policy.h"
#include "gridfire/utils/table_format.h"

#include "gridfire/network.h"

#include "fourdst/composition/composition.h"

#include "fourdst/plugin/bundle/bundle.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"
#include "quill/Backend.h"

#include <clocale>
#include <functional>

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"
#include "gridfire/exceptions/error_solver.h"

#include <boost/json/src.hpp>

#include "gridfire/reaction/reaclib.h"
#include "gridfire/solver/strategies/CVODE_solver_strategy.h"


static std::terminate_handler g_previousHandler = nullptr;
boost::json::object g_reaction_contribution_history;
static std::vector<std::pair<double, std::unordered_map<std::string, std::pair<double, double>>>> g_callbackHistory;
static bool s_wrote_abundance_history = false;
static bool s_wrote_reaction_history = false;
void quill_terminate_handler();


inline std::unique_ptr<gridfire::partition::PartitionFunction> build_partition_function() {
    using gridfire::partition::BasePartitionType;
    const auto partitionFunction = gridfire::partition::CompositePartitionFunction({
        BasePartitionType::RauscherThielemann,
        BasePartitionType::GroundState
    });
    return std::make_unique<gridfire::partition::CompositePartitionFunction>(partitionFunction);
}

gridfire::NetIn init(const double temp) {
    std::setlocale(LC_ALL, "");
    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::TraceL2);

    using namespace gridfire;
    const std::vector<double> X            = {0.7081145999999999, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};


    const fourdst::composition::Composition composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);

    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = temp;
    netIn.density = 1.5e2;
    netIn.energy = 0;

    netIn.tMax = 3e17;
    netIn.dt0 = 1e-12;

    return netIn;
}

void log_results(const gridfire::NetOut& netOut, const gridfire::NetIn& netIn) {
    std::vector<fourdst::atomic::Species> logSpecies = {
        fourdst::atomic::H_1,
        fourdst::atomic::He_3,
        fourdst::atomic::He_4,
        fourdst::atomic::C_12,
        fourdst::atomic::N_14,
        fourdst::atomic::O_16,
        fourdst::atomic::Ne_20,
        fourdst::atomic::Mg_24
    };

    std::vector<double> initial;
    std::vector<double> final;
    std::vector<double> delta;
    std::vector<double> fractional;
    for (const auto& species : logSpecies) {
        double initial_X = netIn.composition.getMassFraction(species);
        double final_X = netOut.composition.getMassFraction(species);
        double delta_X = final_X - initial_X;
        double fractionalChange = (delta_X) / initial_X * 100.0;

        initial.push_back(initial_X);
        final.push_back(final_X);
        delta.push_back(delta_X);
        fractional.push_back(fractionalChange);
    }

    initial.push_back(0.0); // Placeholder for energy
    final.push_back(netOut.energy);
    delta.push_back(netOut.energy);
    fractional.push_back(0.0); // Placeholder for energy

    initial.push_back(0.0);
    final.push_back(netOut.dEps_dT);
    delta.push_back(netOut.dEps_dT);
    fractional.push_back(0.0);

    initial.push_back(0.0);
    final.push_back(netOut.dEps_dRho);
    delta.push_back(netOut.dEps_dRho);
    fractional.push_back(0.0);

    initial.push_back(netIn.composition.getMeanParticleMass());
    final.push_back(netOut.composition.getMeanParticleMass());
    delta.push_back(final.back() - initial.back());
    fractional.push_back((final.back() - initial.back()) / initial.back() * 100.0);

    std::vector<std::string> rowLabels = [&]() -> std::vector<std::string> {
        std::vector<std::string> labels;
        for (const auto& species : logSpecies) {
            labels.push_back(std::string(species.name()));
        }
        labels.push_back("ε");
        labels.push_back("dε/dT");
        labels.push_back("dε/dρ");
        labels.push_back("<μ>");
        return labels;
    }();


    gridfire::utils::Column<std::string> paramCol("Parameter", rowLabels);
    gridfire::utils::Column<double> initialCol("Initial", initial);
    gridfire::utils::Column<double> finalCol  ("Final", final);
    gridfire::utils::Column<double> deltaCol  ("δ", delta);
    gridfire::utils::Column<double> percentCol("% Change", fractional);

    std::vector<std::unique_ptr<gridfire::utils::ColumnBase>> columns;
    columns.push_back(std::make_unique<gridfire::utils::Column<std::string>>(paramCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(initialCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(finalCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(deltaCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(percentCol));


    gridfire::utils::print_table("Simulation Results", columns);
}


void record_abundance_history_callback(const gridfire::solver::CVODESolverStrategy::TimestepContext& ctx) {
    s_wrote_abundance_history = true;
    const auto& engine = ctx.engine;
    // std::unordered_map<std::string, std::pair<double, double>> abundances;
    std::vector<double> Y;
    for (const auto& species : engine.getNetworkSpecies()) {
        const size_t sid = engine.getSpeciesIndex(species);
        double y =  N_VGetArrayPointer(ctx.state)[sid];
        Y.push_back(y > 0.0 ? y : 0.0); // Regularize tiny negative abundances to zero
    }

    fourdst::composition::Composition comp(engine.getNetworkSpecies(), Y);


    std::unordered_map<std::string, std::pair<double, double>> abundances;
    for (const auto& sp : comp | std::views::keys) {
        abundances.emplace(std::string(sp.name()), std::make_pair(sp.mass(), comp.getMolarAbundance(sp)));
    }
    g_callbackHistory.emplace_back(ctx.t, abundances);
}

size_t g_iters = 0;
void record_contribution_callback(const gridfire::solver::CVODESolverStrategy::TimestepContext& ctx) {
    s_wrote_reaction_history = true;
    boost::json::object timestep;
    boost::json::object reaction_contribution;
    boost::json::object species_abundance;
    std::set<fourdst::atomic::Species> activeSpecies(ctx.engine.getNetworkSpecies().begin(), ctx.engine.getNetworkSpecies().end());
    for (const auto& [species, contributions] : ctx.reactionContributionMap) {
        boost::json::object species_obj;
        for (const auto& [reaction_id, contribution] : contributions) {
            species_obj[reaction_id] = contribution;
        }
        reaction_contribution[std::string(species.name())] = species_obj;

        double y;
        if (activeSpecies.contains(species)) {
            const size_t sid = ctx.engine.getSpeciesIndex(species);
            y =  N_VGetArrayPointer(ctx.state)[sid];
        } else {
            y = 0.0;
        }

        species_abundance[std::string(species.name())] = y;
    }
    timestep["t"] = ctx.t;
    timestep["dt"] = ctx.dt;
    timestep["reaction_contribution"] = reaction_contribution;
    timestep["species_abundance"] = species_abundance;
    g_reaction_contribution_history[std::to_string(g_iters)] = timestep;
    g_iters++;
}

void save_callback_data(const std::string_view filename) {
    std::set<std::string> unique_species;
    for (const auto &abundances: g_callbackHistory | std::views::values) {
        for (const auto &species_name: abundances | std::views::keys) {
            unique_species.insert(species_name);
        }
    }
    std::ofstream csvFile(filename.data(), std::ios::out);
    csvFile << "t,";

    size_t i = 0;
    for (const auto& species_name : unique_species) {
        csvFile << species_name;
        if (i < unique_species.size() - 1) {
            csvFile << ",";
        }
        i++;
    }

    csvFile << "\n";

    for (const auto& [time, data] : g_callbackHistory) {
        csvFile << time << ",";
        size_t j = 0;
        for (const auto& species_name : unique_species) {
            if (!data.contains(species_name)) {
                csvFile << "0.0";
            } else {
                csvFile << data.at(species_name).second;
            }
            if (j < unique_species.size() - 1) {
                csvFile << ",";
            }
            ++j;
        }
        csvFile << "\n";
    }

    csvFile.close();
}

void log_callback_data() {
    if (s_wrote_abundance_history) {
        std::cout << "Saving abundance history to abundance_history.csv" << std::endl;
        save_callback_data("abundance_history.csv");
    }

    if (s_wrote_reaction_history) {
        std::cout << "Saving reaction history to reaction_contribution_history.json" << std::endl;
        std::ofstream jsonFile("reaction_contribution_history.json", std::ios::out);
        jsonFile << boost::json::serialize(g_reaction_contribution_history);
        jsonFile.close();
    }
}

void quill_terminate_handler()
{
    log_callback_data();
    quill::Backend::stop();
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

void callback_main(const gridfire::solver::CVODESolverStrategy::TimestepContext& ctx) {
    record_abundance_history_callback(ctx);
    record_contribution_callback(ctx);
}

int main() {
    using namespace gridfire;
    const NetIn netIn = init(1.5e7);

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    stellarPolicy.construct();
    DynamicEngine& engine = stellarPolicy.construct();

    solver::CVODESolverStrategy solver(engine);
    solver.set_callback(solver::CVODESolverStrategy::TimestepCallback(callback_main));
    const NetOut netOut = solver.evaluate(netIn);
    log_results(netOut, netIn);
    log_callback_data();
}