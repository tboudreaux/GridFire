// ReSharper disable CppUnusedIncludeDirective
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <format>

#include "gridfire/gridfire.h"

#include "fourdst/composition/composition.h"
#include "fourdst/logging/logging.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"

#include "quill/Logger.h"
#include "quill/Backend.h"
#include "CLI/CLI.hpp"

#include <clocale>
#include <cmath>

#include "gridfire/utils/gf_omp.h"

#include "nlohmann/json.hpp"

struct IntermediateResult {
    double time;
    fourdst::composition::Composition comp;
    double current_energy;
    double current_neutrino_loss_rate;
};

static std::vector<IntermediateResult> g_callbackHistory;

gridfire::NetIn init(const double temp, const double rho, const double tMax) {
    std::setlocale(LC_ALL, "");
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::Info);

    using namespace gridfire;
    constexpr double XpXn = 7.17;
    constexpr double Xn = 1.0 / (1.0 + XpXn);
    constexpr double Xp = 1 - Xn;

    const std::vector<double> X            = {Xp, Xn};
    const std::vector<std::string> symbols = {"H-1", "n-1"};


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
        double y =  N_VGetArrayPointer(ctx.state)[sid];
        Y.push_back(y > 0.0 ? y : 0.0); // Regularize tiny negative abundances to zero
    }

    const fourdst::composition::Composition comp(engine.getNetworkSpecies(ctx.state_ctx), Y);
    IntermediateResult stepResult;
    stepResult.comp = comp;
    stepResult.time = ctx.t;
    stepResult.current_energy = ctx.current_total_energy;
    stepResult.current_neutrino_loss_rate = ctx.current_neutrino_energy_loss_rate;
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
        j.push_back(entry);
    }
    std::ofstream ofs(filename);
    ofs << j.dump(4);
}

double T9(const double age) {
    return 10.0/std::sqrt(age);
}

double density(const double age) {
    return 4e-5 * std::pow(T9(age), 3);
}

int main(int argc, char** argv) {
    GF_PAR_INIT();
    using namespace gridfire;

    double tMax = 3600;
    double h = 0.1;

    CLI::App app("GridFire Quick BBN Test");
    app.add_option("--tmax", tMax, "Maximum Time in seconds")->default_val(std::format("{:5.2E}", tMax));
    app.add_option("--s", h, "Geometric Timestep Scale Factor")->default_val(std::format("{:5.2f}", h));

    CLI11_PARSE(app, argc, argv);

    NetIn netIn = init(0, 0, tMax);
    const engine::GraphEngine engine(netIn.composition);
    auto blob = std::make_unique<engine::scratch::StateBlob>();
    blob->enroll<engine::scratch::GraphEngineScratchPad>();
    auto* graph_engine_state = engine::scratch::get_state<engine::scratch::GraphEngineScratchPad, false>(*blob);
    graph_engine_state->initialize(engine);


    solver::PointSolverContext solver_ctx(*blob);
    solver::PointSolver solver(engine);
    solver_ctx.stdout_logging=true;

    double current_time = 180;
    nlohmann::json j;
    nlohmann::json meta;

    meta["tMax"] = tMax;
    meta["tStart"] = current_time;
    meta["h"] = h;
    j.push_back(meta);

    nlohmann::json steps;

    while (current_time < tMax) {
        nlohmann::json entry;
        double current_dt = h * current_time;
        double next_time = current_time + current_dt;

        netIn.tMax = current_dt;
        double current_temp = T9(current_time) * 1e9;
        double next_temp = T9(next_time) * 1e9;
        double burn_temp = (current_temp + next_temp)/2.0;

        double current_density = density(current_time);
        double next_density = density(next_time);
        double burn_density = (current_density + next_density)/2.0;

        netIn.temperature = burn_temp;
        netIn.density = burn_density;


        fourdst::composition::Composition initial_comp = netIn.composition;

        NetOut result = solver.evaluate(solver_ctx, netIn);
        netIn.composition = result.composition;
        std::println("Time: {:5.2E} (+{:5.2E}), Burn Temp: {:5.2E}, Burn Density: {:5.2E}", current_time, current_dt, burn_temp, burn_density);

        const fourdst::composition::Composition& comp = result.composition;

        auto Xi = [&](const std::string& symbol) -> double {
            if (!initial_comp.contains(symbol)) {
                return 0;
            }
            return initial_comp.getMassFraction(symbol);
        };

        entry["t"] = current_time;
        entry["T"] = burn_temp;
        entry["D"] = burn_density;
        entry["mu"] = comp.getMeanParticleMass();
        auto lifetimes = engine.getSpeciesTimescales(*blob, comp, burn_temp, burn_density);
        if (lifetimes.has_value()) {
            entry["tau_n"] = lifetimes.value().at(fourdst::atomic::n_1);
        }

        for (const auto& [sp, finalY] : comp) {
            double initialY = Xi(std::string(sp.name()));
            entry[std::format("{}_f", sp.name())] = finalY;
            entry[std::format("{}_i", sp.name())] = initialY;
        }

        steps.push_back(entry);


        current_time += current_dt;
    }
    j.push_back(steps);
    std::ofstream out("BBNResults.json");
    out << j.dump(4);
}
