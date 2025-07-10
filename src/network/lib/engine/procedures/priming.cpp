#include "gridfire/engine/procedures/priming.h"
#include "gridfire/engine/views/engine_priming.h"
#include "gridfire/solver/solver.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

namespace gridfire {
    fourdst::composition::Composition primeNetwork(const NetIn& netIn, DynamicEngine& engine) {
        using fourdst::composition::Composition;
        using fourdst::atomic::Species;


        NetIn currentConditions = netIn;
        std::vector<Species> speciesToPrime;
        for (const auto &entry: netIn.composition | std::views::values) {
            if (entry.mass_fraction() == 0.0) {
                speciesToPrime.push_back(entry.isotope());
            }
        }


        if (speciesToPrime.empty()) {
            return netIn.composition; // No priming needed
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature to T9
        const double rho = netIn.density; // Density in g/cm^3

        for (const auto& primingSpecies :speciesToPrime) {
            NetworkPrimingEngineView primer(primingSpecies, engine);

            const auto Y = primer.mapNetInToMolarAbundanceVector(currentConditions);
            const double destructionRateConstant = calculateDestructionRateConstant(primer, primingSpecies, Y, T9, rho);
            const double creationRate = calculateCreationRate(primer, primingSpecies, Y, T9, rho);
            const double equilibriumMolarAbundance = creationRate / destructionRateConstant;
            const double equilibriumMassFraction = equilibriumMolarAbundance * primingSpecies.mass();

            currentConditions.composition.setMassFraction(
                std::string(primingSpecies.name()),
                equilibriumMassFraction
            );
        }

        return currentConditions.composition;
    }

    double calculateDestructionRateConstant(
        const DynamicEngine& engine,
        const fourdst::atomic::Species& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) {
        const int speciesIndex = engine.getSpeciesIndex(species);
        std::vector<double> Y_scaled(Y.begin(), Y.end());
        Y_scaled[speciesIndex] = 1.0; // Set the abundance of the species to 1.0 for rate constant calculation
        double destructionRateConstant = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            if (reaction.contains_reactant(species)) {
                const int stoichiometry = reaction.stoichiometry(species);
                destructionRateConstant += std::abs(stoichiometry) * engine.calculateMolarReactionFlow(reaction, Y_scaled, T9, rho);
            }
        }
        return destructionRateConstant;
    }

    double calculateCreationRate(
        const DynamicEngine& engine,
        const fourdst::atomic::Species& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) {
        double creationRate = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            const int stoichiometry = reaction.stoichiometry(species);
            if (stoichiometry > 0) {
                creationRate += stoichiometry * engine.calculateMolarReactionFlow(reaction, Y, T9, rho);
            }
        }
        return creationRate;
    }
}