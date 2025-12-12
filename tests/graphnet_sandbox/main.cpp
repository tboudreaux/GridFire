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
#include "fourdst/composition/utils.h"

#include "quill/Logger.h"
#include "quill/Backend.h"
#include "CLI/CLI.hpp"

#include <clocale>

#include "gridfire/reaction/reaclib.h"


static std::terminate_handler g_previousHandler = nullptr;
static std::vector<std::pair<double, std::unordered_map<std::string, std::pair<double, double>>>> g_callbackHistory;
static bool s_wrote_abundance_history = false;
void quill_terminate_handler();

gridfire::NetIn init(const double temp, const double rho, const double tMax) {
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


void record_abundance_history_callback(const gridfire::solver::CVODESolverStrategy::TimestepContext& ctx) {
    s_wrote_abundance_history = true;
    const auto& engine = ctx.engine;
    // std::unordered_map<std::string, std::pair<double, double>> abundances;
    std::vector<double> Y;
    for (const auto& species : engine.getNetworkSpecies(ctx.state_ctx)) {
        const size_t sid = engine.getSpeciesIndex(ctx.state_ctx, species);
        double y =  N_VGetArrayPointer(ctx.state)[sid];
        Y.push_back(y > 0.0 ? y : 0.0); // Regularize tiny negative abundances to zero
    }

    fourdst::composition::Composition comp(engine.getNetworkSpecies(ctx.state_ctx), Y);


    std::unordered_map<std::string, std::pair<double, double>> abundances;
    for (const auto& sp : comp | std::views::keys) {
        abundances.emplace(std::string(sp.name()), std::make_pair(sp.mass(), comp.getMolarAbundance(sp)));
    }
    g_callbackHistory.emplace_back(ctx.t, abundances);
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

void log_callback_data(const double temp) {
    if (s_wrote_abundance_history) {
        std::cout << "Saving abundance history to abundance_history.csv" << std::endl;
        save_callback_data("abundance_history_" + std::to_string(temp) + ".csv");
    }

}

void quill_terminate_handler()
{
    log_callback_data(1.5e7);
    quill::Backend::stop();
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

void callback_main(const gridfire::solver::CVODESolverStrategy::TimestepContext& ctx) {
    record_abundance_history_callback(ctx);
}

int main() {
    using namespace gridfire;

    constexpr size_t breaks = 1;
    double temp = 1.5e7;
    double rho = 1.5e2;
    double tMax = 3.1536e+16/breaks;

    const NetIn netIn = init(temp, rho, tMax);

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    policy::ConstructionResults construct = stellarPolicy.construct();
    std::println("Sandbox Engine Stack: {}", stellarPolicy);
    std::println("Scratch Blob State: {}", *construct.scratch_blob);


    constexpr size_t runs = 1000;
    auto startTime = std::chrono::high_resolution_clock::now();

    // arrays to store timings
    std::array<std::chrono::duration<double>, runs> setup_times;
    std::array<std::chrono::duration<double>, runs> eval_times;
    std::array<NetOut, runs> serial_results;
    for (size_t i = 0; i < runs; ++i) {
        auto start_setup_time = std::chrono::high_resolution_clock::now();
        std::print("Run {}/{}\r", i + 1, runs);
        solver::CVODESolverStrategy solver(construct.engine, *construct.scratch_blob);
        // solver.set_callback(solver::CVODESolverStrategy::TimestepCallback(callback_main));
        solver.set_stdout_logging_enabled(false);
        auto end_setup_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> setup_elapsed = end_setup_time - start_setup_time;
        setup_times[i] = setup_elapsed;

        auto start_eval_time = std::chrono::high_resolution_clock::now();
        const NetOut netOut = solver.evaluate(netIn);
        auto end_eval_time = std::chrono::high_resolution_clock::now();
        serial_results[i] = netOut;
        std::chrono::duration<double> eval_elapsed = end_eval_time - start_eval_time;
        eval_times[i] = eval_elapsed;

        // log_results(netOut, netIn);
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;
    std::println("");

    // Summarize serial timings
    double total_setup_time = 0.0;
    double total_eval_time = 0.0;
    for (size_t i = 0; i < runs; ++i) {
        total_setup_time += setup_times[i].count();
        total_eval_time += eval_times[i].count();
    }
    std::println("Average Setup Time over {} runs: {:.6f} seconds", runs, total_setup_time / runs);
    std::println("Average Evaluation Time over {} runs: {:.6f} seconds", runs, total_eval_time / runs);
    std::println("Total Time for {} runs: {:.6f} seconds", runs, elapsed.count());
    std::println("Final H-1 Abundances Serial: {}", serial_results[0].composition.getMolarAbundance(fourdst::atomic::H_1));

    // OPTIONAL: Prevent CppAD from returning memory to the system
    // during execution to reduce overhead (can speed up tight loops)
    CppAD::thread_alloc::hold_memory(true);

    std::array<NetOut, runs> parallelResults;
    std::array<std::chrono::duration<double>, runs> setupTimes;
    std::array<std::chrono::duration<double>, runs> evalTimes;
    std::array<std::unique_ptr<gridfire::engine::scratch::StateBlob>, runs> workspaces;
    for (size_t i = 0; i < runs; ++i) {
        workspaces[i] = construct.scratch_blob->clone_structure();
    }


    // Parallel runs
    startTime = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < runs; ++i) {
        auto start_setup_time = std::chrono::high_resolution_clock::now();
        solver::CVODESolverStrategy solver(construct.engine, *workspaces[i]);
        solver.set_stdout_logging_enabled(false);
        auto end_setup_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> setup_elapsed = end_setup_time - start_setup_time;
        setupTimes[i] = setup_elapsed;
        auto start_eval_time = std::chrono::high_resolution_clock::now();
        parallelResults[i] = solver.evaluate(netIn);
        auto end_eval_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> eval_elapsed = end_eval_time - start_eval_time;
        evalTimes[i] = eval_elapsed;
    }
    endTime = std::chrono::high_resolution_clock::now();
    elapsed = endTime - startTime;
    std::println("");

    // Summarize parallel timings
    total_setup_time = 0.0;
    total_eval_time = 0.0;
    for (size_t i = 0; i < runs; ++i) {
        total_setup_time += setupTimes[i].count();
        total_eval_time += evalTimes[i].count();
    }

    std::println("Average Parallel Setup Time over {} runs: {:.6f} seconds", runs, total_setup_time / runs);
    std::println("Average Parallel Evaluation Time over {} runs: {:.6f} seconds", runs, total_eval_time / runs);
    std::println("Total Parallel Time for {} runs: {:.6f} seconds", runs, elapsed.count());

    std::println("Final H-1 Abundances Parallel: {}", utils::iterable_to_delimited_string(parallelResults, ",", [](const auto& result) {
        return result.composition.getMolarAbundance(fourdst::atomic::H_1);
    }));
}