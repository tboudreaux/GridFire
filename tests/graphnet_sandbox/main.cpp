#include <iostream>
#include <fstream>

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_approx8.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/partition/partition_types.h"
#include "gridfire/engine/views/engine_defined.h"
#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/engine/procedures/priming.h"
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

static std::terminate_handler g_previousHandler = nullptr;

void measure_execution_time(const std::function<void()>& callback, const std::string& name)
{
    const auto startTime   = std::chrono::steady_clock::now();
    callback();
    const auto endTime     = std::chrono::steady_clock::now();
    const auto duration    = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);
    std::cout << "Execution time for " << name << ": "
              << duration.count()/1e9 << " s\n";
}

void quill_terminate_handler()
{
    quill::Backend::stop();
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

int main() {
    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::Info);
    LOG_DEBUG(logger, "Starting Adaptive Engine View Example...");

    using namespace gridfire;
    const std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};


    fourdst::composition::Composition composition;
    composition.registerSymbol(symbols, true);
    composition.setMassFraction(symbols, comp);
    composition.finalize(true);
    using partition::BasePartitionType;
    const auto partitionFunction = partition::CompositePartitionFunction({
        BasePartitionType::RauscherThielemann,
        BasePartitionType::GroundState
    });

    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 5e9;
    netIn.density = 1.6e6;
    netIn.energy = 0;
    // netIn.tMax = 3.1536e17; // ~ 10Gyr
    netIn.tMax = 1e-14;
    netIn.dt0 = 1e-12;

    GraphEngine ReaclibEngine(composition, partitionFunction, NetworkBuildDepth::SecondOrder);
    ReaclibEngine.setUseReverseReactions(true);
    // ReaclibEngine.setScreeningModel(screening::ScreeningType::WEAK);
    //
    MultiscalePartitioningEngineView partitioningView(ReaclibEngine);
    AdaptiveEngineView adaptiveView(partitioningView);
    //
    solver::DirectNetworkSolver solver(adaptiveView);
    NetOut netOut;
    netOut = solver.evaluate(netIn);
    std::cout << "Initial H-1: " << netIn.composition.getMassFraction("H-1") << std::endl;
    std::cout << "NetOut H-1: " << netOut.composition.getMassFraction("H-1") << std::endl;
    std::cout << "Consumed " << (netIn.composition.getMassFraction("H-1") - netOut.composition.getMassFraction("H-1")) * 100 << " % H-1 by mass" << std::endl;
    // measure_execution_time([&](){netOut = solver.evaluate(netIn);}, "DirectNetworkSolver Evaluation");
    // std::cout << "DirectNetworkSolver completed in " << netOut.num_steps << " steps.\n";
    // std::cout << "Final composition:\n";
    // for (const auto& [symbol, entry] : netOut.composition) {
    //     std::cout << symbol << ": " << entry.mass_fraction() << "\n";
    // }

}