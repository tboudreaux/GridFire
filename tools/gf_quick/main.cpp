// ReSharper disable CppUnusedIncludeDirective
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <format>

#include "gridfire/gridfire.h"
#include <cppad/utility/thread_alloc.hpp> // Required for parallel_setup

#include "fourdst/composition/composition.h"
#include "fourdst/logging/logging.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils/utils.h"

#include "quill/Logger.h"
#include "quill/Backend.h"
#include "CLI/CLI.hpp"

#include <clocale>

#include "gridfire/utils/gf_omp.h"

#include "nlohmann/json.hpp"

struct IntermediateResult {
    double time{};
    fourdst::composition::Composition comp;
    gridfire::reaction::ReactionSet reactions;
    std::vector<double> reaction_flows;
    double current_energy{};
    double current_neutrino_loss_rate{};

    gridfire::reaction::ReactionSet inactive_reactions;
    std::vector<double> inactive_reaction_flows;
};

static std::vector<IntermediateResult> g_callbackHistory;

gridfire::NetIn init(const double temp, const double rho, const double tMax) {
    std::setlocale(LC_ALL, "");
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::TraceL2);

    using namespace gridfire;
    const std::vector<double> X            = {0.7081145999999999, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};


    const fourdst::composition::Composition composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);

    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = temp;
    netIn.density = rho;
    netIn.energy = 0;

    netIn.tMax = tMax;
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

    initial.push_back(0.0);
    final.push_back(netOut.specific_neutrino_energy_loss);
    delta.push_back(netOut.specific_neutrino_energy_loss);
    fractional.push_back(0.0);

    initial.push_back(0.0);
    final.push_back(netOut.specific_neutrino_flux);
    delta.push_back(netOut.specific_neutrino_flux);
    fractional.push_back(0.0);

    initial.push_back(netIn.composition.getMeanParticleMass());
    final.push_back(netOut.composition.getMeanParticleMass());
    delta.push_back(final.back() - initial.back());
    fractional.push_back((final.back() - initial.back()) / initial.back() * 100.0);

    std::vector<std::string> rowLabels = [&]() -> std::vector<std::string> {
        std::vector<std::string> labels;
        for (const auto& species : logSpecies) {
            labels.emplace_back(species.name());
        }
        labels.emplace_back("ε");
        labels.emplace_back("dε/dT");
        labels.emplace_back("dε/dρ");
        labels.emplace_back("Eν");
        labels.emplace_back("Fν");
        labels.emplace_back("<μ>");
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


void record_abundance_history_callback(const gridfire::solver::PointSolverTimestepContext& ctx) {
    const auto& engine = ctx.engine;
    std::vector<double> Y;
    for (const auto& species : engine.getNetworkSpecies(ctx.state_ctx)) {
        const size_t sid = engine.getSpeciesIndex(ctx.state_ctx, species);
        const double y = ctx.abundance(sid);
        Y.push_back(y > 0.0 ? y : 0.0);
    }
    const fourdst::composition::Composition comp(engine.getNetworkSpecies(ctx.state_ctx), Y);
    IntermediateResult stepResult;
    stepResult.comp = comp;
    stepResult.time = ctx.t;
    stepResult.current_energy = ctx.current_total_energy;
    stepResult.current_neutrino_loss_rate = ctx.current_neutrino_energy_loss_rate;
    stepResult.reactions = engine.getNetworkReactions(ctx.state_ctx);

    for (const auto& reactionPtr : stepResult.reactions) {
        double flow = engine.calculateMolarReactionFlow(ctx.state_ctx, *reactionPtr, comp, ctx.T9, ctx.rho);
        stepResult.reaction_flows.push_back(flow);
    }

    stepResult.inactive_reactions = engine.getInactiveNetworkReactions(ctx.state_ctx);
    for (const auto& reactionPtr : stepResult.inactive_reactions) {
        double flow = engine.getInactiveReactionMolarReactionFlow(ctx.state_ctx, *reactionPtr, comp, ctx.T9, ctx.rho);
        stepResult.inactive_reaction_flows.push_back(flow);
    }
    g_callbackHistory.push_back(stepResult);
}



void callback_main(const gridfire::solver::PointSolverTimestepContext& ctx) {
    record_abundance_history_callback(ctx);
}

void save_callback(const std::string& filename) {
    // Save to JSON
    nlohmann::json j;
    for (const auto& record : g_callbackHistory) {
        nlohmann::json entry;
        entry["time"] = record.time;
        entry["current_energy"] = record.current_energy;
        entry["current_neutrino_loss_rate"] = record.current_neutrino_loss_rate;
        // make a sub-json for composition
        nlohmann::json comp_json;
        for (const auto& [species, abundance] : record.comp) {
            comp_json[species.name()] = abundance;
        }
        entry["composition"] = comp_json;
        entry["reactions"] = nlohmann::json::array();
        for (const auto& [reaction, flow] : std::views::zip(record.reactions, record.reaction_flows)) {
            nlohmann::json reaction_info;
            reaction_info["id"] = reaction->id();
            reaction_info["flow"] = flow;
            reaction_info["species"] = nlohmann::json::array();
            reaction_info["Q"] = reaction->qValue();
            for (const auto& sp : reaction->all_species()) {
                nlohmann::json species_info;
                species_info["name"] = sp.name();
                species_info["stoichiometry"] = reaction->stoichiometry(sp);
                reaction_info["species"].push_back(species_info);
            }
            entry["reactions"].push_back(reaction_info);
        }

        entry["inactive_reactions"] = nlohmann::json::array();
        for (const auto& [reaction, flow] : std::views::zip(record.inactive_reactions, record.inactive_reaction_flows)) {
            nlohmann::json reaction_info;
            reaction_info["id"] = reaction->id();
            reaction_info["flow"] = flow;
            reaction_info["species"] = nlohmann::json::array();
            reaction_info["Q"] = reaction->qValue();
            for (const auto& sp : reaction->all_species()) {
                nlohmann::json species_info;
                species_info["name"] = sp.name();
                species_info["stoichiometry"] = reaction->stoichiometry(sp);
                reaction_info["species"].push_back(species_info);
            }
            entry["inactive_reactions"].push_back(reaction_info);
        }
        j.push_back(entry);
    }
    std::ofstream ofs(filename);
    ofs << j.dump(4);
}

int main(int argc, char** argv) {
    GF_PAR_INIT();
    using namespace gridfire;

    double temp = 1.5e7;
    double rho = 1.5e2;
    double tMax = 3.1536e+16;
    bool save_intermediate_results = false;
    bool display_trigger = false;
    std::string output_filename = "abundance_history.json";


    CLI::App app("GridFire Quick CLI Test");
    app.add_option("--temp", temp, "Initial Temperature")->default_val(std::format("{:5.2E}", temp));
    app.add_option("--rho", rho, "Initial Density")->default_val(std::format("{:5.2E}", rho));
    app.add_option("--tmax", tMax, "Maximum Time")->default_val(std::format("{:5.2E}", tMax));
    app.add_option("--save_intermediate_results", save_intermediate_results, "Save Intermediate Results")->default_val("false");
    app.add_option("--output", output_filename, "Output filename for intermediate results")->default_val("abundance_history.json");
    app.add_option("--display_trigger_explanations", display_trigger, "Display trigger explanations during run")->default_val("false");

    CLI11_PARSE(app, argc, argv);
    NetIn netIn = init(temp, rho, tMax);

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    auto [engine, ctx_template] = stellarPolicy.construct();

    solver::PointSolverContext solver_context(*ctx_template);
    solver::PointSolver solver(engine);
    if (save_intermediate_results) {
        solver_context.callback = solver::TimestepCallback(callback_main);
    }

    NetOut result = solver.evaluate(solver_context, netIn, display_trigger);
    log_results(result, netIn);

    if (save_intermediate_results) {
        save_callback(output_filename);
    }
}
