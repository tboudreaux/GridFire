#include <iostream>
#include <fstream>

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_approx8.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/partition/partition_types.h"
#include "gridfire/engine/views/engine_defined.h"
#include "gridfire/io/network_file.h"

#include "gridfire/solver/solver.h"

#include "gridfire/network.h"

#include "fourdst/composition/composition.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"
#include "quill/Backend.h"
#include "quill/Frontend.h"

#include <chrono>
#include <functional>

#include "gridfire/partition/composite/partition_composite.h"

// Keep a copy of the previous handler
static std::terminate_handler g_previousHandler = nullptr;

void measure_execution_time(const std::function<void()>& callback, const std::string& name)
{
    // variable names in camelCase
    const auto startTime   = std::chrono::steady_clock::now();
    callback();
    const auto endTime     = std::chrono::steady_clock::now();
    const auto duration    = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);
    std::cout << "Execution time for " << name << ": "
              << duration.count()/1e9 << " s\n";
}

void quill_terminate_handler()
{
    // 1. Stop the Quill backend (flushes all sinks and joins thread)
    quill::Backend::stop();

    // 2. Optionally, flush individual loggers:
    // for (auto logger : quill::get_all_loggers())
    //   logger->flush_log();

    // 3. Chain to the previous terminate handler (or abort if none)
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

int main() {
    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::TraceL2);
    LOG_DEBUG(logger, "Starting Adaptive Engine View Example...");

    using namespace gridfire;
    const std::vector<double> comp = {0.708, 0.0, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "H-2", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};

    fourdst::composition::Composition composition;
    composition.registerSymbol(symbols, true);
    composition.setMassFraction(symbols, comp);
    composition.finalize(true);


    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 1.5e7;
    netIn.density = 1e2;
    netIn.energy = 0.0;

    netIn.tMax = 3.1536e17;
    netIn.dt0 = 1e12;

    NetOut netOut;

    // netIn.dt0 = 1e12;
    // approx8::Approx8Network approx8Network;
    // measure_execution_time([&]() {
    //     netOut = approx8Network.evaluate(netIn);
    // }, "Approx8 Network Initialization");
    // std::cout << "Approx8 Network H-1: " << netOut.composition.getMassFraction("H-1") << " in " << netOut.num_steps << " steps." << std::endl;

    using partition::BasePartitionType;
    const auto partitionFunction = partition::CompositePartitionFunction({
        BasePartitionType::RauscherThielemann,
        BasePartitionType::GroundState
    });
    std::cout << "Partition Function for Mg-24: " << partitionFunction.evaluate(12, 24, 8) << std::endl;
    std::cout << "Partition Function for F-23: " << partitionFunction.evaluate(9, 23, 8) << std::endl;
    std::cout << "Partition Function for O-13: " << partitionFunction.evaluate(8, 13, 8) << std::endl;

    netIn.dt0 = 1e-15;

    GraphEngine ReaclibEngine(composition, partitionFunction);
    std::cout << ReaclibEngine.getPartitionFunction().type() << std::endl;
    // ReaclibEngine.setPrecomputation(true);
    // // AdaptiveEngineView adaptiveEngine(ReaclibEngine);
    // io::SimpleReactionListFileParser parser{};
    // FileDefinedEngineView approx8EngineView(ReaclibEngine, "approx8.net", parser);
    // approx8EngineView.setScreeningModel(screening::ScreeningType::WEAK);
    // solver::QSENetworkSolver solver(approx8EngineView);
    // netOut = solver.evaluate(netIn);

    // measure_execution_time([&]() {
    //     netOut = solver.evaluate(netIn);
    // }, "Approx8 Network Evaluation (Precomputation)");
    // ReaclibEngine.setPrecomputation(false);
    // std::cout << "Precomputation H-1: " << netOut.composition.getMassFraction("H-1") << " in " << netOut.num_steps << " steps." << std::endl;
}