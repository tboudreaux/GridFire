// ReSharper disable CppUnusedIncludeDirective
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <format>

#include "fourdst/logging/logging.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"

#include "quill/Logger.h"
#include "quill/Backend.h"
#include "CLI/CLI.hpp"

#include <clocale>

#include "gridfire/gridfire.h"
#include "fourdst/composition/composition.h"
#include "gridfire/utils/gf_omp.h"

#include <atomic>
#include <new>
#include <cstdlib>

static std::terminate_handler g_previousHandler = nullptr;
void quill_terminate_handler();

gridfire::NetIn init(const double temp, const double rho, const double tMax) {
    std::setlocale(LC_ALL, "");
    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::Info);

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


void quill_terminate_handler()
{
    quill::Backend::stop();
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

int main(int argc, char* argv[]) {
    using namespace gridfire;

    double temp = 1.5e7;
    double rho = 1.6e2;
    double tMax = 3e17;

    std::string output_filename = "gridfire_timings.csv";

    CLI::App app("GridFire Timeing Benchmarks");
    app.add_option("--temperature", temp, "Temperature in degrees")->default_val(std::format("{:5.2E}", temp));
    app.add_option("--density", rho, "Density in Kg")->default_val(std::format("{:5.2E}", rho));
    app.add_option("--tmax", tMax, "Maximum time in seconds")->default_val(std::format("{:5.2E}", tMax));
    app.add_option("--output", output_filename, "Output filename for intermediate results")->default_val("gridfire_timings.csv");

    CLI11_PARSE(app, argc, argv);

    const NetIn netIn = init(temp, rho, tMax);

    std::unique_ptr<engine::GraphEngine> engine;


    struct TimingInfo  {
        double depth;
        int num_reactions;
        int num_species;
        double timing;
    };

    std::vector<TimingInfo> timings;

    int prev_reactions = 0;
    int prev_species = 0;
    engine = std::make_unique<engine::GraphEngine>(netIn.composition, 1);
    for (int depth = 1; depth <= 100; depth++) {
        engine = std::make_unique<engine::GraphEngine>(netIn.composition, depth);
        auto blob = engine->constructStateBlob();

        TimingInfo info;
        info.depth = depth;
        info.num_species = engine->getNetworkSpecies(*blob).size();
        info.num_reactions = engine->getNetworkReactions(*blob).size();

        if (prev_reactions == info.num_reactions && prev_species == info.num_species) {
            std::println("Found end of useful graph traversal at a depth of {}", depth);
            break;
        }
        const solver::PointSolver localSolver(*engine);
        solver::PointSolverContext solverCtx(*blob);


        solverCtx.stdout_logging = true;

        try {
            auto start = std::chrono::high_resolution_clock::now();
            auto result = localSolver.evaluate(solverCtx, netIn, false, false);
            auto end = std::chrono::high_resolution_clock::now();

            double ns = std::chrono::duration<double>(end - start).count();

            info.timing = ns;

        prev_reactions = info.num_reactions;
        prev_species = info.num_species;

            timings.push_back(info);
        } catch (gridfire::exceptions::CVODESolverFailureError& e) {
            continue;
        }
    }

    std::ofstream csvFile(output_filename, std::ios::out);
    csvFile << std::format("# Temperature (K): {}", temp);
    csvFile << std::format("# Density: {}", rho);
    csvFile << std::format("# TMax: {}", tMax);
    csvFile << "depth,reactions,species,time\n";
    for (const auto& [depth, numReactions, numSpecies, ns]: timings) {
        std::string line = std::format("{},{},{},{}\n", depth, numReactions, numSpecies, ns);
        csvFile << line;
    }
    csvFile.close();

    std::println("Timeing Benchmarks results written to {}", output_filename);
}
