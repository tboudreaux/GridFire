#include <iostream>
#include <fstream>

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_approx8.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/partition/partition_types.h"
#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/solver/strategies/CVODE_solver_strategy.h"

#include "gridfire/network.h"

#include "fourdst/composition/composition.h"

#include "fourdst/plugin/bundle/bundle.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"
#include "quill/Backend.h"
#include "quill/Frontend.h"

#include <chrono>
#include <functional>

#include "gridfire/engine/views/engine_defined.h"
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

int main(int argc, char* argv[]){

    // Valid usages are either
    // ./graphnet_sandbox
    //or
    // ./graphnet_sandbox --plug <plugin_bundle_path>
    if (argc == 3 && std::string(argv[1]) == "--plug") {
        std::filesystem::path pluginBundlePath(argv[2]);
        if (!std::filesystem::exists(pluginBundlePath)) {
            std::cerr << "Error: Plugin bundle path does not exist: " << pluginBundlePath << "\n";
            std::cerr << "Usage: " << argv[0] << " [--plug <plugin_bundle_path>]\n";
            return 1;
        }
        std::cout << "Loading plugin bundle from: " << pluginBundlePath << "\n";
        fourdst::plugin::bundle::PluginBundle pluginBundle(pluginBundlePath);
    }
    if (argc == 2 && std::string(argv[1]) != "--plug") {
        std::cerr << "Invalid argument: " << argv[1] << "\n";
        std::cerr << "Usage: " << argv[0] << " [--plug <plugin_bundle_path>]\n";
        return 1;
    }
    if (argc == 2 && std::string(argv[1]) == "--plug") {
        std::cerr << "Error: No plugin bundle path provided.\n";
        std::cerr << "Usage: " << argv[0] << " [--plug <plugin_bundle_path>]\n";
        return 1;
    }
    if (argc > 3) {
        std::cerr << "Too many arguments provided.\n";
        std::cerr << "Usage: " << argv[0] << " [--plug <plugin_bundle_path>]\n";
        return 1;
    }


    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::TraceL3);
    LOG_INFO(logger, "Starting Adaptive Engine View Example...");

    using namespace gridfire;
    const std::vector<double> comp = {0.708, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};


    fourdst::composition::Composition composition;
    composition.registerSymbol(symbols, true);
    composition.setMassFraction(symbols, comp);
    bool didFinalize = composition.finalize(true);
    if (!didFinalize) {
        std::cerr << "Failed to finalize initial composition." << std::endl;
        return 1;
    }
    using partition::BasePartitionType;
    const auto partitionFunction = partition::CompositePartitionFunction({
        BasePartitionType::RauscherThielemann,
        BasePartitionType::GroundState
    });

    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 1.5e7;
    netIn.density = 1.6e2;
    netIn.energy = 0;
    netIn.tMax = 3e16;
    // netIn.tMax = 1e-14;
    netIn.dt0 = 1e-12;

    GraphEngine ReaclibEngine(composition, partitionFunction, NetworkBuildDepth::SecondOrder);
    ReaclibEngine.setUseReverseReactions(false);
    ReaclibEngine.setPrecomputation(false);
    // DefinedEngineView ppEngine({"p(p,e+)d", "d(p,g)he3", "he3(he3,2p)he4"}, ReaclibEngine);

    MultiscalePartitioningEngineView partitioningView(ReaclibEngine);
    AdaptiveEngineView adaptiveView(partitioningView);

    solver::CVODESolverStrategy solver(adaptiveView);
    NetOut netOut;


    netOut = solver.evaluate(netIn);
    // consumptionFile.close();
    std::cout << "Initial H-1: " << netIn.composition.getMassFraction("H-1") << std::endl;
    std::cout << "NetOut H-1: " << netOut.composition.getMassFraction("H-1") << std::endl;

    double initialHydrogen = netIn.composition.getMassFraction("H-1");
    double finalHydrogen = netOut.composition.getMassFraction("H-1");
    double fractionalConsumedHydrogen = (initialHydrogen - finalHydrogen) / initialHydrogen * 100.0;
    std::cout << "Fractional consumed hydrogen: " << fractionalConsumedHydrogen << "%" << std::endl;
    std::cout << netOut << std::endl;

}