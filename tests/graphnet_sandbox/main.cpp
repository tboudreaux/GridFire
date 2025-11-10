#include <iostream>
#include <fstream>

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/solver/strategies/CVODE_solver_strategy.h"
#include "gridfire/policy/stellar_policy.h"
#include "gridfire/utils/table_format.h"

#include "gridfire/network.h"

#include "fourdst/composition/composition.h"

#include "fourdst/plugin/bundle/bundle.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"
#include "quill/Backend.h"

#include <clocale>
#include <functional>

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"


static std::terminate_handler g_previousHandler = nullptr;

void quill_terminate_handler()
{
    quill::Backend::stop();
    if (g_previousHandler)
        g_previousHandler();
    else
        std::abort();
}

gridfire::NetIn init() {
    std::setlocale(LC_ALL, "");
    g_previousHandler = std::set_terminate(quill_terminate_handler);
    quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
    logger->set_log_level(quill::LogLevel::TraceL1);

    using namespace gridfire;
    const std::vector<double> X            = {0.7081145999999999, 2.94e-5, 0.276, 0.003, 0.0011, 9.62e-3, 1.62e-3, 5.16e-4};
    const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};


    fourdst::composition::Composition composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);

    NetIn netIn;
    netIn.composition = composition;
    netIn.temperature = 1.5e7;
    netIn.density = 1.6e2;
    netIn.energy = 0;

    netIn.tMax = 3e16;
    netIn.dt0 = 1e-12;

    return netIn;
}

void log_results(const gridfire::NetOut& netOut, const gridfire::NetIn& netIn) {
    double initialH1_X = netIn.composition.getMassFraction("H-1");
    double finalH1_X = netOut.composition.getMassFraction("H-1");
    double deltaH1_X = finalH1_X - initialH1_X;
    double fractionalChangeH1 = (deltaH1_X) / initialH1_X * 100.0;

    double initialHe4_X = netIn.composition.getMassFraction("He-4");
    double finalHe4_X = netOut.composition.getMassFraction("He-4");
    double deltaHe4_X = finalHe4_X - initialHe4_X;
    double fractionalChangeHe4 = (deltaHe4_X) / initialHe4_X * 100.0;

    double initialBe7_X = 0.0;
    double finalBe7_X = netOut.composition.getMassFraction("Be-7");
    double deltaBe7_X = finalBe7_X - initialBe7_X;
    double fractionalChangeBe7 = (deltaBe7_X) / initialBe7_X * 100.0;

    double initialC12_X = netIn.composition.getMassFraction("C-12");
    double finalC12_X = netOut.composition.getMassFraction("C-12");
    double deltaC12_X = finalC12_X - initialC12_X;
    double fractionalChangeC12 = (deltaC12_X) / initialC12_X * 100.0;

    double finalEnergy = netOut.energy;
    double dEpsdT = netOut.dEps_dT;
    double dEpsdRho = netOut.dEps_dRho;

    double initialMeanMolecularMass = netIn.composition.getMeanParticleMass();
    double finalMeanMolecularMass = netOut.composition.getMeanParticleMass();

    std::vector<std::string> rowLabels = {"H-1", "He-4", "Be-7", "C-12", "ε", "dε/dT", "dε/dρ", "<μ>"};
    std::vector<std::string> header = {"Parameter", "Initial", "Final", "Δ", "% Change"};
    std::vector<double> H1Data = {initialH1_X, finalH1_X, deltaH1_X, fractionalChangeH1};
    std::vector<double> He4Data = {initialHe4_X, finalHe4_X, deltaHe4_X, fractionalChangeHe4};
    std::vector<double> Be7Data = {initialBe7_X, finalBe7_X, deltaBe7_X, fractionalChangeBe7};
    std::vector<double> C12Data = {initialC12_X, finalC12_X, deltaC12_X, fractionalChangeC12};
    std::vector<double> energyData = {0.0, finalEnergy, finalEnergy, 0.0};
    std::vector<double> dEpsdTData = {0.0, dEpsdT, dEpsdT, 0.0};
    std::vector<double> dEpsdRhoData = {0.0, dEpsdRho, dEpsdRho, 0.0};
    std::vector<double> meanMolecularMassData = {initialMeanMolecularMass, finalMeanMolecularMass, finalMeanMolecularMass - initialMeanMolecularMass,
                                                 (finalMeanMolecularMass - initialMeanMolecularMass) / initialMeanMolecularMass * 100.0};

    gridfire::utils::Column<std::string> paramCol(header[0], rowLabels);
    gridfire::utils::Column<double> initialCol(header[1], {H1Data[0], He4Data[0], Be7Data[0], C12Data[0], energyData[0], dEpsdTData[0], dEpsdRhoData[0], meanMolecularMassData[0]});
    gridfire::utils::Column<double> finalCol  (header[2], {H1Data[1], He4Data[1], Be7Data[1], C12Data[1], energyData[1], dEpsdTData[1], dEpsdRhoData[1], meanMolecularMassData[1]});
    gridfire::utils::Column<double> deltaCol  (header[3], {H1Data[2], He4Data[2], Be7Data[2], C12Data[2], energyData[2], dEpsdTData[2], dEpsdRhoData[2], meanMolecularMassData[2]});
    gridfire::utils::Column<double> percentCol(header[4], {H1Data[3], He4Data[3], Be7Data[3], C12Data[3], energyData[3], dEpsdTData[3], dEpsdRhoData[3], meanMolecularMassData[3]});

    std::vector<std::unique_ptr<gridfire::utils::ColumnBase>> columns;
    columns.push_back(std::make_unique<gridfire::utils::Column<std::string>>(paramCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(initialCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(finalCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(deltaCol));
    columns.push_back(std::make_unique<gridfire::utils::Column<double>>(percentCol));


    gridfire::utils::print_table("Simulation Results", columns);
}


static std::vector<std::pair<double, std::unordered_map<std::string, std::pair<double, double>>>> g_callbackHistory;
void record_abundance_history_callback(const gridfire::solver::CVODESolverStrategy::TimestepContext& ctx) {
    const auto& engine = ctx.engine;
    std::unordered_map<std::string, std::pair<double, double>> abundances;
    for (const auto& species : engine.getNetworkSpecies()) {
        const size_t sid = engine.getSpeciesIndex(species);
        abundances[std::string(species.name())] = std::make_pair(species.mass(), N_VGetArrayPointer(ctx.state)[sid]);
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

int main() {
    using namespace gridfire;
    const NetIn netIn = init();

    policy::MainSequencePolicy stellarPolicy(netIn.composition);
    DynamicEngine& engine = stellarPolicy.construct();

    solver::CVODESolverStrategy solver(engine);
    solver.set_callback(solver::CVODESolverStrategy::TimestepCallback(record_abundance_history_callback));

    const NetOut netOut = solver.evaluate(netIn, false);

    log_results(netOut, netIn);
    save_callback_data("abundance_history.csv");
}