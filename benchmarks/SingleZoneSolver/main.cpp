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

#include <clocale>

#include "gridfire/reaction/reaclib.h"
#include "gridfire/utils/gf_omp.h"

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


int main() {
    GF_PAR_INIT()
    using namespace gridfire;

    constexpr size_t breaks = 1;
    constexpr double temp = 1.5e7;
    constexpr double rho = 1.5e2;
    constexpr double tMax = 3.1536e+16/breaks;

    const NetIn netIn = init(temp, rho, tMax);

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    const policy::ConstructionResults construct = stellarPolicy.construct();
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
        solver::CVODESolverStrategy solver(construct.engine, *construct.scratch_blob);
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


    std::array<NetOut, runs> parallelResults;
    std::array<std::chrono::duration<double>, runs> setupTimes;
    std::array<std::chrono::duration<double>, runs> evalTimes;
    std::array<std::unique_ptr<gridfire::engine::scratch::StateBlob>, runs> workspaces;
    for (size_t i = 0; i < runs; ++i) {
        workspaces[i] = construct.scratch_blob->clone_structure();
    }


    // Parallel runs
    startTime = std::chrono::high_resolution_clock::now();

    GF_OMP(parallel for,)
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

    std::println("========== Summary ==========");
    std::println("Serial Runs:");
    std::println("  Average Setup Time: {:.6f} seconds", total_setup_time / runs);
    std::println("  Average Evaluation Time: {:.6f} seconds", total_eval_time / runs);
    std::println("Parallel Runs:");
    std::println("  Average Setup Time: {:.6f} seconds", total_setup_time / runs);
    std::println("  Average Evaluation Time: {:.6f} seconds", total_eval_time / runs);
    std::println("Difference:");
    std::println("  Setup Time Difference: {:.6f} seconds", (total_setup_time / runs) - (total_setup_time / runs));
    std::println("  Evaluation Time Difference: {:.6f} seconds", (total_eval_time / runs) - (total_eval_time / runs));
    std::println("  Setup Time Fractional Difference: {:.2f}%", ((total_setup_time / runs) - (total_setup_time / runs)) / (total_setup_time / runs) * 100.0);
    std::println("  Evaluation Time Fractional Difference: {:.2f}%", ((total_eval_time / runs) - (total_eval_time / runs)) / (total_eval_time / runs) * 100.0);
}