// ReSharper disable CppUnusedIncludeDirective
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <format>

#include <cppad/utility/thread_alloc.hpp> // Required for parallel_setup

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

static std::atomic<size_t> g_allocated_bytes{0};

void* operator new(std::size_t size) {
    g_allocated_bytes += size;
    if (void* ptr = std::malloc(size)) {
        return ptr;
    }
    throw std::bad_alloc();
}

void operator delete(void* ptr, std::size_t size) {
    g_allocated_bytes -= size;
    std::free(ptr);
}

void operator delete(void* ptr) {
    std::free(ptr);
}

struct MemoryScopeTracker {
    size_t start_bytes;
    MemoryScopeTracker() : start_bytes(g_allocated_bytes.load()) {}
    size_t bytes_allocated() const {
        return g_allocated_bytes.load() - start_bytes;
    }

    void reset_tracking() {
        start_bytes = 0;
        g_allocated_bytes = 0;
    }
};


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
    std::string output_filename = "gf_mem.csv";

    CLI::App app("GridFire Memory Benchmarks");
    app.add_option("--temperature", temp, "Temperature in degrees")->default_val(std::format("{:5.2E}", temp));
    app.add_option("--density", rho, "Density in Kg")->default_val(std::format("{:5.2E}", rho));
    app.add_option("--tmax", tMax, "Maximum time in seconds")->default_val(std::format("{:5.2E}", tMax));
    app.add_option("--output", output_filename, "Output filename for intermediate results")->default_val("gf_mem.csv");

    CLI11_PARSE(app, argc, argv);

    const NetIn netIn = init(temp, rho, tMax);

    std::unique_ptr<engine::GraphEngine> engine;

    int prev_reactions = 0;
    int prev_species = 0;
    engine = std::make_unique<engine::GraphEngine>(netIn.composition, 1);

    MemoryScopeTracker tracker;
    std::ofstream mem_file(output_filename, std::ios::out);
    mem_file << "depth,species,reactions,engine_memory_bytes,solver_memory_bytes\n";
    for (int depth = 1; depth <= 100; depth++) {
        tracker.reset_tracking();
        engine = std::make_unique<engine::GraphEngine>(netIn.composition, depth);
        auto blob = engine->constructStateBlob();
        size_t engine_usage = tracker.bytes_allocated();

        size_t current_num_species = engine->getNetworkSpecies(*blob).size();
        size_t current_num_reactions = engine->getNetworkReactions(*blob).size();

        if (prev_reactions == current_num_reactions && prev_species == current_num_species) {
            std::println("Found end of useful graph traversal at a depth of {}", depth);
            break;
        }
        tracker.reset_tracking();
        const solver::PointSolver localSolver(*engine);
        solver::PointSolverContext solverCtx(*blob);
        size_t solver_usage = tracker.bytes_allocated();

        mem_file << std::format("{},{},{},{},{}\n", depth, current_num_species, current_num_reactions, engine_usage, solver_usage);

        prev_reactions = current_num_reactions;
        prev_species = current_num_species;
    }
    mem_file.close();
    std::println("Memory benchmarks results written to {}", output_filename);
}
