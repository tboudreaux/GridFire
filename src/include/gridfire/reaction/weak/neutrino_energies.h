#pragma once

#include "gridfire/reaction/reaction.h"

#include <optional>

namespace gridfire::reaction::weak {
    /**
     * @brief Average energy carried away by the neutrino in a weak reaction [MeV], if tabulated.
     *
     * REACLIB weak rates carry no neutrino information. This looks the reaction up in a table generated from MESA's
     * rate data (utils/neutrino/generate_neutrino_energies.py): pp, pep, hep, 7Be, 8B, 13N and 15O from Bahcall (2002),
     * 17F from Bahcall (1997), and most heavier parent -> daughter decays and captures from weaklib (evaluated by MESA
     * at logT = 7.5, logRho = 4). Using the same values as MESA keeps GridFire's energy accounting comparable.
     *
     * @return The average neutrino energy, or std::nullopt for reactions not in the table (callers fall back to an
     *         approximation).
     */
    std::optional<double> tabulated_neutrino_energy_MeV(const Reaction& reaction);
}
