#include "gridfire/reaction/weak/neutrino_energies.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <algorithm>
#include <cctype>
#include <string>
#include <string_view>
#include <vector>

namespace {
    struct NeutrinoEnergyEntry {
        std::string_view reactants; ///< Sorted MESA iso names, space separated
        std::string_view products; ///< Sorted MESA iso names, space separated
        std::string_view lepton; ///< "e+" (emitted positron), "e-" (captured/emitted electron), or "" for any
        double energy_MeV;
    };

    constexpr NeutrinoEnergyEntry k_neutrino_energies[] = {
#include "neutrino_energies_data.inc"
    };

    std::string iso_name(const fourdst::atomic::Species& species) {
        if (species.z() == 0 && species.a() == 1) return "neut";
        std::string element(species.el());
        std::ranges::transform(element, element.begin(), [](const unsigned char c) { return std::tolower(c); });
        return element + std::to_string(species.a());
    }

    std::string joined_sorted(const std::vector<fourdst::atomic::Species>& species) {
        std::vector<std::string> names;
        names.reserve(species.size());
        for (const auto& sp : species) names.push_back(iso_name(sp));
        std::ranges::sort(names);
        std::string joined;
        for (const auto& name : names) joined += (joined.empty() ? "" : " ") + name;
        return joined;
    }

    /// The lepton named in a REACLIB id, e.g. "e+" for p(p,e+)d and "e-" for be7(e-,)li7.
    std::string_view lepton_of(const std::string_view id) {
        if (id.find("e+") != std::string_view::npos) return "e+";
        if (id.find("e-") != std::string_view::npos) return "e-";
        return "";
    }
}

namespace gridfire::reaction::weak {
    std::optional<double> tabulated_neutrino_energy_MeV(const Reaction& reaction) {
        const std::string reactants = joined_sorted(reaction.reactants());
        const std::string products = joined_sorted(reaction.products());
        const std::string_view lepton = lepton_of(reaction.id());
        for (const auto& entry : k_neutrino_energies) {
            if (entry.reactants == reactants && entry.products == products &&
                (entry.lepton.empty() || entry.lepton == lepton)) {
                return entry.energy_MeV;
            }
        }
        return std::nullopt;
    }
}
