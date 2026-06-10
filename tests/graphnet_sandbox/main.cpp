// ReSharper disable CppUnusedIncludeDirective
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <format>

#include <cppad/utility/thread_alloc.hpp> // Required for parallel_setup

#include "fourdst/logging/logging.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils/utils.h"

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
    double coupling_ratio = 5.0;
    std::string output_filename = "coupling.dat";

    CLI::App app("GridFire Test Coupling");
    app.add_option("--temperature", temp, "Temperature in degrees")->default_val(std::format("{:5.2E}", temp));
    app.add_option("--density", rho, "Density in Kg")->default_val(std::format("{:5.2E}", rho));
    app.add_option("--tmax", tMax, "Maximum time in seconds")->default_val(std::format("{:5.2E}", tMax));
    app.add_option("--coupling_ratio", coupling_ratio, "Coupling ratio for multiscale partitioning")->default_val(std::format("{:.2f}", coupling_ratio));
    app.add_option("--output", output_filename, "Output filename for intermediate results")->default_val("coupling.dat");

    CLI11_PARSE(app, argc, argv);

    const NetIn netIn = init(temp, rho, tMax);

    auto base_engine = std::make_unique<engine::GraphEngine>(netIn.composition, 3);
    auto base_blob = base_engine->constructStateBlob();
    auto qse_engine = std::make_unique<engine::MultiscalePartitioningEngineView>(*base_engine);
    auto blob = qse_engine->constructStateBlob(base_blob.get());

    auto* state = engine::scratch::get_state<engine::scratch::MultiscalePartitioningEngineViewScratchPad, true>(*blob);

    const solver::PointSolver localSolver(*base_engine);
    solver::PointSolverContext solverCtx(*base_blob);

    auto result = localSolver.evaluate(solverCtx, netIn, false, false);

}