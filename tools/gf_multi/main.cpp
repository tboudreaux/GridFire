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

#include "gridfire/utils/gf_omp.h"

static std::terminate_handler g_previousHandler = nullptr;
static std::vector<std::pair<double, std::unordered_map<std::string, std::pair<double, double>>>> g_callbackHistory;
static bool s_wrote_abundance_history = false;
void quill_terminate_handler();

enum class ScalingTypes {
    LINEAR,
    LOG,
    GEOM
};

std::map<std::string, ScalingTypes> scaling_type_map = {
    {"linear", ScalingTypes::LINEAR},
    {"log",    ScalingTypes::LOG},
    {"geom",   ScalingTypes::GEOM}
};

template<typename T>
concept IsOrderableLinear = std::floating_point<T>;

template<typename T>
concept IsOrderableLog = std::floating_point<T>;

template<IsOrderableLinear T>
constexpr std::vector<T> linspace(T start, T end, size_t N) {
    if (N == 0) {
        return {};
    }
    if (N == 1) {
        return {start};
    }

    std::vector<T> result{};
    result.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        const T t = static_cast<T>(i) / static_cast<T>(N - 1);
        result[i] = std::lerp(start, end, t);
    }

    return result;
}

template<IsOrderableLog T>
std::vector<T> logspace(T start, T end, size_t N) {
    if (N == 0) {
        return {};
    }
    if (N == 1) {
        return {start};
    }

    std::vector<T> result{};
    result.resize(N);

    const T log_start = start;
    const T log_end = end;

    for (std::size_t i = 0; i < N; ++i) {
        const T t = static_cast<T>(i) / static_cast<T>(N - 1);
        const T exponent = std::lerp(log_start, log_end, t);
        result[i] = std::pow(static_cast<T>(10), exponent);
    }

    return result;
}

template<IsOrderableLog T>
std::vector<T> geomspace(T start, T end, size_t N) {
    if (start <= 0 || end <= 0) {
        throw std::domain_error("geomspace requires positive start/end values");
    }

    const T log_start = std::log10(start);
    const T log_end = std::log10(end);

    std::vector<T> result = logspace<T>(log_start, log_end, N);

    result[0] = start;
    result[N - 1] = end;

    return result;
}



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


void record_abundance_history_callback(const gridfire::solver::PointSolverTimestepContext& ctx) {
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

void callback_main(const gridfire::solver::PointSolverTimestepContext& ctx) {
    record_abundance_history_callback(ctx);
}

int main(int argc, char** argv) {
    GF_PAR_INIT();
    using namespace gridfire;

    double tMin = 1e7;
    double tMax = 3e7;
    ScalingTypes tempScaling = ScalingTypes::LINEAR;

    double rhoMin = 1e2;
    double rhoMax = 1e3;
    ScalingTypes rhoScaling = ScalingTypes::LOG;

    double evolveTime = 1e13;
    size_t shells = 10;


    CLI::App app("GridFire Quick CLI Test");
    // Add temp, rho, and tMax as options if desired
    app.add_option("--tMin", tMin, "Initial Temperature")->default_val(std::format("{:5.2E}", tMin));
    app.add_option("--tMax", tMax, "Maximum Temperature")->default_val(std::format("{:5.2E}", tMax));
    app.add_option("--tempScaling", tempScaling, "Temperature Scaling Type (linear, log, geom)")->default_val(ScalingTypes::LINEAR)->transform(CLI::CheckedTransformer(scaling_type_map, CLI::ignore_case));

    app.add_option("--rhoMin", rhoMin, "Initial Density")->default_val(std::format("{:5.2E}", rhoMin));
    app.add_option("--rhoMax", rhoMax, "Maximum Density")->default_val(std::format("{:5.2E}", rhoMax));
    app.add_option("--rhoScaling", rhoScaling, "Density Scaling Type (linear, log, geom)")->default_val(ScalingTypes::LOG)->transform(CLI::CheckedTransformer(scaling_type_map, CLI::ignore_case));

    app.add_option("--shell", shells, "Number of Shells")->default_val(shells);

    app.add_option("--evolveTime", evolveTime, "Evolution Time")->default_val(std::format("{:5.2E}", evolveTime));


    CLI11_PARSE(app, argc, argv);
    NetIn rootNetIn = init(tMin, tMax, evolveTime);

    std::vector<double> temps;
    std::vector<double> densities;

    switch (tempScaling) {
        case ScalingTypes::LINEAR:
            temps = linspace<double>(tMin, tMax, static_cast<size_t>(shells));
            break;
        case ScalingTypes::LOG:
            temps = logspace<double>(std::log10(tMin), std::log10(tMax), static_cast<size_t>(shells));
            break;
        case ScalingTypes::GEOM:
            temps = geomspace<double>(tMin, tMax, static_cast<size_t>(shells));
            break;
    }
    switch (rhoScaling) {
        case ScalingTypes::LINEAR:
            densities = linspace<double>(rhoMin, rhoMax, static_cast<size_t>(shells));
            break;
        case ScalingTypes::LOG:
            densities = logspace<double>(std::log10(rhoMin), std::log10(rhoMax), static_cast<size_t>(shells));
            break;
        case ScalingTypes::GEOM:
            densities = geomspace<double>(rhoMin, rhoMax, static_cast<size_t>(shells));
            break;
    }

    std::vector<NetIn> netIns;
    netIns.reserve(shells);
    for (size_t i = 0; i < static_cast<size_t>(shells); ++i) {
        NetIn netIn = rootNetIn;
        netIn.temperature = temps[i];
        netIn.density = densities[i];
        netIns.emplace_back(netIn);
    }

    policy::MainSequencePolicy stellarPolicy(rootNetIn.composition);
    auto [engine, ctx_template] = stellarPolicy.construct();

    solver::PointSolver point_solver(engine);
    solver::GridSolver grid_solver(engine, point_solver);
    solver::GridSolverContext solver_ctx(*ctx_template);


    std::vector<NetOut> results = grid_solver.evaluate(solver_ctx, netIns);
    for (size_t i = 0; i < results.size(); ++i) {
        std::cout << "=== Shell " << i << " ===" << std::endl;
        log_results(results[i], netIns[i]);
    }
}