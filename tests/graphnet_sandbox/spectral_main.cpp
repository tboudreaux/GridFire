#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#include "gridfire/gridfire.h"

#include "fourdst/composition/composition.h"
#include "fourdst/logging/logging.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"

#include "quill/Logger.h"
#include "quill/Backend.h"
#include "CLI/CLI.hpp"

#include <clocale>


static std::terminate_handler g_previousHandler = nullptr;
static std::vector<std::pair<double, std::unordered_map<std::string, std::pair<double, double>>>> g_callbackHistory;
static bool s_wrote_abundance_history = false;
void quill_terminate_handler();

std::vector<double> linspace(const double start, const double end, const size_t num) {
    std::vector<double> result;
    if (num == 0) return result;
    if (num == 1) {
        result.push_back(start);
        return result;
    }
    const double step = (end - start) / static_cast<double>(num - 1);
    for (size_t i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

std::vector<gridfire::NetIn> init(const double tMin, const double tMax, const double rhoMin, const double rhoMax, const double nShells, const double evolveTime) {
    std::setlocale(LC_ALL, "");
    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::TraceL2);
    LOG_INFO(logger, "Initializing GridFire Spectral Solver Sandbox...");

    using namespace gridfire;
    const std::vector<double> X            = {0.7081145999999999, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};


    const fourdst::composition::Composition composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);

    std::vector<NetIn> netIns;
    for (const auto& [T, ρ]: std::views::zip(linspace(tMin, tMax, nShells), linspace(rhoMax, rhoMin, nShells))) {
        NetIn netIn;
        netIn.composition = composition;
        netIn.temperature = T;
        netIn.density = ρ;
        netIn.energy = 0;

        netIn.tMax = evolveTime;
        netIn.dt0 = 1e-12;
        netIns.push_back(netIn);
    }

    return netIns;
}
void quill_terminate_handler()
{
    quill::Backend::stop();
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

int main(int argc, char** argv) {
    using namespace gridfire;

    CLI::App app{"GridFire Sandbox Application."};

    double tMin = 1.5e7;
    double tMax = 1.7e7;
    double rhoMin = 1.5e2;
    double rhoMax = 1.5e2;
    double nShells = 15;
    double evolveTime = 3.1536e+16;

    app.add_option("--tMin", tMin, "Minimum time in seconds");
    app.add_option("--tMax", tMax, "Maximum time in seconds");
    app.add_option("--rhoMin", rhoMin, "Minimum density in g/cm^3");
    app.add_option("--rhoMax", rhoMax, "Maximum density in g/cm^3");
    app.add_option("--nShells", nShells, "Number of shells");
    app.add_option("--evolveTime", evolveTime, "Maximum time in seconds");

    CLI11_PARSE(app, argc, argv);

    const std::vector<NetIn> netIns = init(tMin, tMax, rhoMin, rhoMax, nShells, evolveTime);

    policy::MainSequencePolicy stellarPolicy(netIns[0].composition);
    stellarPolicy.construct();
    policy::ConstructionResults construct = stellarPolicy.construct();

    solver::SpectralSolverStrategy solver(construct.engine);
    std::vector<double> mass_coords = linspace(1e-5, 1.0, nShells);

    std::vector<NetOut> results = solver.evaluate(netIns, mass_coords, *construct.scratch_blob);
}
