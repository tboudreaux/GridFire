// Per-reaction energy generation of a validation case's GridFire network at one state; the GridFire counterpart of
// mesa_rates/ (MESA), for diagnosing GridFire/MESA differences reaction by reaction.
//
//   gf_mesa_rates <cases.json> <case name> <composition file> <T [K]> <rho [g/cm^3]> [excluded_reactions.txt]
//
// The composition file has one "<MESA iso name> <mass fraction>" per line (the same file mesa_rates reads). Output
// (CSV on stdout):
//   total,<eps_nuc>,<eps_neu>,<screening model>
//   reaction,<MESA name>,<GridFire id>,<forward rate>,<net molar flow [mol/g/s]>,<Q from masses [MeV]>,<E_nu [MeV]>,<eps [erg/g/s]>
#include "mesa_case.h"

#include "gridfire/io/generative/mesa.h"
#include "gridfire/reaction/weak/neutrino_energies.h"
#include "gridfire/screening/screening_types.h"

#include "fourdst/constants/const.h"
#include "fourdst/logging/logging.h"
#include "quill/Logger.h"

#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>

int main(int argc, char** argv) {
    if (argc != 6 && argc != 7) {
        std::cerr << "usage: gf_mesa_rates <cases.json> <case name> <composition file> <T [K]> <rho [g/cm^3]> [excluded_reactions.txt]\n";
        return 2;
    }
    fourdst::logging::LogManager::getInstance().getLogger("log")->set_log_level(quill::LogLevel::Warning);
    const std::string wanted = argv[2];
    const double T9 = std::stod(argv[4]) / 1e9;
    const double rho = std::stod(argv[5]);

    for (const auto& burn_case : gridfire::test::mesa::load_cases(argv[1])) {
        if (burn_case.name != wanted) continue;
        const auto excluded = argc == 7 ? gridfire::test::mesa::read_excluded_reactions(argv[6]) : std::vector<std::string>{};
        auto network = gridfire::test::mesa::build_network(burn_case, excluded);
        auto& ctx = *network.context;
        const auto& species = network.engine->getNetworkSpecies(ctx);

        std::map<std::string, double> X_by_iso;
        std::ifstream comp_in(argv[3]);
        for (std::string iso; comp_in >> iso;) {
            double x;
            comp_in >> x;
            X_by_iso[iso] = x;
        }
        std::vector<double> Y(species.size(), 0.0);
        std::unordered_map<size_t, fourdst::atomic::Species> index_map;
        double Ye = 0.0;
        for (size_t i = 0; i < species.size(); ++i) {
            index_map.emplace(i, species[i]);
            const auto it = X_by_iso.find(gridfire::io::generative::get_mesa_iso_name(species[i]));
            if (it != X_by_iso.end()) Y[i] = it->second / species[i].mass();
            Ye += species[i].z() * Y[i];
        }
        const fourdst::composition::Composition comp(species, Y);

        const auto rhs = network.engine->calculateRHSAndEnergy(ctx, comp, T9, rho, false);
        if (!rhs) {
            std::cerr << "RHS evaluation failed\n";
            return 1;
        }
        const auto screening = network.engine->getScreeningModel(ctx);
        std::cout << std::format("total,{:.16e},{:.16e},{}\n", rhs->nuclearEnergyGenerationRate, rhs->neutrinoEnergyLossRate,
                                 screening == gridfire::screening::ScreeningType::BARE ? "BARE" : "WEAK");

        const auto& k = fourdst::constant::Constants::getInstance();
        const double c = k.get("c").value, u = k.get("u").value, Na = k.get("N_a").value;
        constexpr double MeV = 1.602176634e-6;
        for (const auto& reaction : network.engine->getNetworkReactions(ctx)) {
            double dm = 0.0; // amu released per reaction
            for (const auto& sp : reaction->reactants()) dm += sp.mass();
            for (const auto& sp : reaction->products()) dm -= sp.mass();
            const double Q = dm * u * c * c / MeV;
            const double E_nu = gridfire::reaction::weak::tabulated_neutrino_energy_MeV(*reaction).value_or(std::nan(""));
            const double rate = reaction->calculate_rate(T9, rho, Ye, 0.0, Y, index_map);
            const double flow = network.engine->calculateMolarReactionFlow(ctx, *reaction, comp, T9, rho);
            const double eps = flow * Na * (Q - (std::isnan(E_nu) ? 0.0 : E_nu)) * MeV;
            std::cout << std::format("reaction,{},{},{:.16e},{:.16e},{:.10f},{:.6f},{:.16e}\n",
                                     gridfire::io::generative::get_mesa_reaction_name(*reaction), reaction->id(), rate,
                                     flow, Q, E_nu, eps);
        }
        return 0;
    }
    std::cerr << "no case named " << wanted << "\n";
    return 1;
}
