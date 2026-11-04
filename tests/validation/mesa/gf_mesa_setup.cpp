// Writes everything bbq needs to run a validation case on exactly GridFire's network.
//
//   gf_mesa_setup <cases.json> <case name> <output dir> [excluded_reactions.txt]
//
// Output:
//   network.net               MESA network file (GridFire's network, exported)
//   reactions.txt             GridFire reaction ids, sorted (identity of the network for the comparison test)
//   species_map.txt           "<GridFire species> <MESA iso name>" for every network species
//   initial_composition.txt   "<GridFire species> <mass fraction>" for every network species (0 if not seeded)
//   reaction_map.txt          "<GridFire reaction id>\t<MESA reaction name>\t<participants>" (tab separated, since
//                             reaction ids can contain spaces; participants: sorted MESA
//                             iso names of all reactants and products, comma separated; used to match MESA errors)
#include "mesa_case.h"

#include "gridfire/io/generative/mesa.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"

#include <cstdio>
#include <format>
#include <fstream>
#include <iostream>
#include <algorithm>

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "usage: gf_mesa_setup <cases.json> <case name> <output dir> [excluded_reactions.txt]\n";
        return 2;
    }
    fourdst::logging::LogManager::getInstance().getLogger("log")->set_log_level(quill::LogLevel::Warning);

    const auto cases = gridfire::test::mesa::load_cases(argv[1]);
    const std::string wanted = argv[2];
    const std::filesystem::path out_dir = argv[3];
    std::filesystem::create_directories(out_dir);

    for (const auto& burn_case : cases) {
        if (burn_case.name != wanted) continue;

        const auto excluded = argc == 5 ? gridfire::test::mesa::read_excluded_reactions(argv[4]) : std::vector<std::string>{};
        auto network = gridfire::test::mesa::build_network(burn_case, excluded);
        const auto& species = network.engine->getNetworkSpecies(*network.context);

        std::ofstream(out_dir / "network.net") << gridfire::io::generative::export_engine_to_mesa_net(*network.engine, *network.context, false);

        std::ofstream reactions(out_dir / "reactions.txt");
        for (const auto& id : gridfire::test::mesa::sorted_reaction_ids(network)) reactions << id << "\n";

        std::ofstream reaction_map(out_dir / "reaction_map.txt");
        for (const auto& reaction : network.engine->getNetworkReactions(*network.context)) {
            std::vector<std::string> participants;
            for (const auto& sp : reaction->reactants()) participants.push_back(gridfire::io::generative::get_mesa_iso_name(sp));
            for (const auto& sp : reaction->products()) participants.push_back(gridfire::io::generative::get_mesa_iso_name(sp));
            std::ranges::sort(participants);
            std::string joined;
            for (const auto& p : participants) joined += (joined.empty() ? "" : ",") + p;
            reaction_map << reaction->id() << "\t" << gridfire::io::generative::get_mesa_reaction_name(*reaction) << "\t" << joined << "\n";
        }

        std::ofstream species_map(out_dir / "species_map.txt");
        std::ofstream composition(out_dir / "initial_composition.txt");
        for (const auto& sp : species) {
            species_map << sp.name() << " " << gridfire::io::generative::get_mesa_iso_name(sp) << "\n";
            const double X = network.initial_composition.contains(sp) ? network.initial_composition.getMassFraction(sp) : 0.0;
            composition << sp.name() << " " << std::format("{:.17e}", X) << "\n";
        }

        std::cout << std::format("{}: {} species, {} reactions written to {}\n", burn_case.name, species.size(),
                                 network.engine->getNetworkReactions(*network.context).size(), out_dir.string());
        return 0;
    }
    std::cerr << "no case named " << wanted << " in " << argv[1] << "\n";
    return 1;
}
