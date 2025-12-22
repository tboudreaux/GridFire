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

using namespace fourdst::composition;
Composition rescale(const Composition& comp, double target_X, double target_Z) {
    // 1. Validate inputs
    if (target_X < 0.0 || target_Z < 0.0 || (target_X + target_Z) > 1.0 + 1e-14) {
        throw std::invalid_argument("Target mass fractions X and Z must be non-negative and sum to <= 1.0");
    }

    // Force high precision for the target Y to ensure X+Y+Z = 1.0 exactly in our logic
    long double ld_target_X = static_cast<long double>(target_X);
    long double ld_target_Z = static_cast<long double>(target_Z);
    long double ld_target_Y = 1.0L - ld_target_X - ld_target_Z;

    // Clamp Y to 0 if it dipped slightly below due to precision (e.g. X+Z=1.0000000001)
    if (ld_target_Y < 0.0L) ld_target_Y = 0.0L;

    // 2. Manually calculate current Mass Totals (bypass getCanonicalComposition to avoid crashes)
    long double total_mass_H = 0.0L;
    long double total_mass_He = 0.0L;
    long double total_mass_Z = 0.0L;

    // We need to iterate and identify species types manually
    // Standard definition: H (z=1), He (z=2), Metals (z>2)
    // Note: We use long double accumulators to prevent summation drift
    for (const auto& [spec, molar_abundance] : comp) {
        // Retrieve atomic properties.
        // Note: usage assumes fourdst::atomic::Species has .z() and .mass()
        // consistent with the provided composition.cpp
        int z = spec.z();
        double a = spec.mass();

        long double mass_contribution = static_cast<long double>(molar_abundance) * static_cast<long double>(a);

        if (z == 1) {
            total_mass_H += mass_contribution;
        } else if (z == 2) {
            total_mass_He += mass_contribution;
        } else {
            total_mass_Z += mass_contribution;
        }
    }

    long double total_mass_current = total_mass_H + total_mass_He + total_mass_Z;

    // Edge case: Empty composition
    if (total_mass_current <= 0.0L) {
         // Return empty or throw? If input was empty, return empty.
         if (comp.size() == 0) return comp;
         throw std::runtime_error("Input composition has zero total mass.");
    }

    // 3. Calculate Scaling Factors
    // Factor = (Target_Mass_Fraction / Old_Mass_Fraction)
    //        = (Target_Mass_Fraction) / (Old_Group_Mass / Total_Mass)
    //        = (Target_Mass_Fraction * Total_Mass) / Old_Group_Mass

    long double scale_H = 0.0L;
    long double scale_He = 0.0L;
    long double scale_Z = 0.0L;

    if (ld_target_X > 1e-16L) {
        if (total_mass_H <= 1e-19L) {
             throw std::runtime_error("Cannot rescale Hydrogen to " + std::to_string(target_X) +
                                     " because input has no Hydrogen.");
        }
        scale_H = (ld_target_X * total_mass_current) / total_mass_H;
    }

    if (ld_target_Y > 1e-16L) {
        if (total_mass_He <= 1e-19L) {
             throw std::runtime_error("Cannot rescale Helium to " + std::to_string((double)ld_target_Y) +
                                     " because input has no Helium.");
        }
        scale_He = (ld_target_Y * total_mass_current) / total_mass_He;
    }

    if (ld_target_Z > 1e-16L) {
        if (total_mass_Z <= 1e-19L) {
             throw std::runtime_error("Cannot rescale Metals to " + std::to_string(target_Z) +
                                     " because input has no Metals.");
        }
        scale_Z = (ld_target_Z * total_mass_current) / total_mass_Z;
    }

    // 4. Apply Scaling and Construct New Vectors
    std::vector<fourdst::atomic::Species> new_species;
    std::vector<double> new_abundances;
    new_species.reserve(comp.size());
    new_abundances.reserve(comp.size());

    for (const auto& [spec, abundance] : comp) {
        new_species.push_back(spec);

        long double factor = 0.0L;
        int z = spec.z();

        if (z == 1) {
            factor = scale_H;
        } else if (z == 2) {
            factor = scale_He;
        } else {
            factor = scale_Z;
        }

        // Calculate new abundance in long double then cast back
        long double new_val_ld = static_cast<long double>(abundance) * factor;
        new_abundances.push_back(static_cast<double>(new_val_ld));
    }

    return Composition(new_species, new_abundances);
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

    double temp = 1.5e7;
    double rho = 1.5e2;
    double tMax = 3.1536e+16;


    CLI::App app("GridFire Quick CLI Test");
    app.add_option("--temp", temp, "Initial Temperature")->default_val(std::format("{:5.2E}", temp));
    app.add_option("--rho", rho, "Initial Density")->default_val(std::format("{:5.2E}", rho));
    app.add_option("--tmax", tMax, "Maximum Time")->default_val(std::format("{:5.2E}", tMax));

    CLI11_PARSE(app, argc, argv);
    NetIn netIn = init(temp, rho, tMax);

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    auto [engine, ctx_template] = stellarPolicy.construct();

    solver::PointSolverContext solver_context(*ctx_template);
    solver::PointSolver solver(engine);

    NetOut result = solver.evaluate(solver_context, netIn);
    log_results(result, netIn);
}
