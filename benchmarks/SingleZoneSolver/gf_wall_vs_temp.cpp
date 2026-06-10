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
#include "fourdst/composition/utils/utils.h"

#include "quill/Logger.h"
#include "quill/Backend.h"

#include <clocale>

#include "gridfire/reaction/reaclib.h"
#include "gridfire/utils/gf_omp.h"

template <std::floating_point T>
[[nodiscard]] constexpr auto linspace(T start, T end, std::size_t num_points) -> std::vector<T> {
    if (num_points == 0) {
        return {};
    }
    
    if (num_points == 1) {
        return {start};
    }

    return std::views::iota(0uz, num_points) 
        | std::views::transform([=](std::size_t i) -> T {
            const T t = static_cast<T>(i) / static_cast<T>(num_points - 1);
            
            return std::lerp(start, end, t);
        }) 
        | std::ranges::to<std::vector<T>>(); 
}

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

    constexpr double temp_init = 1.5e7;
    constexpr double rho_init = 1.5e2;
    constexpr double tMax = 3.1536e+12;

    NetIn netIn = init(temp_init, rho_init, tMax);

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    const policy::ConstructionResults construct = stellarPolicy.construct();
    std::println("Sandbox Engine Stack: {}", stellarPolicy);
    std::println("Scratch Blob State: {}", *construct.scratch_blob);


    // arrays to store timings

    // Total number of interpolated data points
    constexpr size_t N = 20;
    std::array<double, N*N> eval_times{};
    auto density = linspace(10.0, 5.0e2, N);
    auto temperature = linspace(4e6,3e7, N);

    solver::PointSolverContext solverCtx(*construct.scratch_blob);
    solverCtx.set_stdout_logging(false);
    solver::PointSolver solver(construct.engine);

    auto startTime = std::chrono::high_resolution_clock::now();


    size_t i = 0;
    for (const auto temp : temperature) {
        for (const auto dens : density) {
            std::println("Evaluation {:3}/{:5} ({:3.0f}%): ρ = {:10.4E}, T = {:10.4E}", i + 1, N*N, 100.0*((static_cast<double>(i)+1.0)/(N*N)), dens, temp);
            netIn.temperature = temp;
            netIn.density = dens;
            try {
                auto start_eval_time = std::chrono::high_resolution_clock::now();
                const NetOut netOut = solver.evaluate(solverCtx, netIn);
                auto end_eval_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> eval_elapsed = end_eval_time - start_eval_time;
                eval_times[i] = eval_elapsed.count();
            } catch (const gridfire::exceptions::GridFireError& e) {
                std::cerr << "Error during evaluation " << (i + 1) << ": " << e.what() << std::endl;
                eval_times[i] = std::numeric_limits<double>::quiet_NaN();
            }
            i++;
        }
    }
    auto endTime = std::chrono::high_resolution_clock::now();

    std::println("Total time for {} evaluations: {} seconds", N, (endTime - startTime).count());
    for (size_t j = 0; j < static_cast<size_t>(N*N); ++j) {
        std::println("Evaluation {}: {} seconds", j + 1, eval_times[j]);
    }

    std::ofstream outfile("gf_wall_vs_temp_results.csv");
    outfile << "Evaluation,Density,Temperature,TimeSeconds\n";
    size_t j = 0;
    for (const auto temp: temperature) {
        for (const auto dens: density ) {
            outfile << (j + 1) << "," << dens << ","<< temp << "," << eval_times[j] << "\n";
            j++;
        }
    }
}