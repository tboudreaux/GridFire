#pragma once
#include "fourdst/composition/composition.h"
#include "fourdst/composition/atomicSpecies.h"

namespace gridfire::utils {
    inline double massFractionFromMolarAbundanceAndComposition (
        const fourdst::composition::Composition& composition,
        const fourdst::atomic::Species& species,
        const double Yi
    ) {
        double sum = 0;
        for (const auto& [symbol, entry] : composition) {
            if (entry.isotope() == species) {
                sum += species.mass() * Yi;
            } else {
                sum += entry.isotope().mass() * composition.getMolarAbundance(symbol);
            }
        }
        return (species.mass() * Yi) / sum;
    };

    /**
     * @brief Convert a vector of molar abundances into a vector of mass fractions
     * @param molarAbundances Vector of molar abundances
     * @param molarMasses Vector of molar masses
     *
     * @note The vectors molarAbundances and molarMasses must be parallel. This function does not provide any checks
     * to ensure that the correct molar mass is being used with the correct molar abundance.
     * @return A vector of molar masses such that each molar mass < 1 and the sum of all is = 1
     */
    inline std::vector<double> massFractionFromMolarAbundanceAndMolarMass (
        const std::vector<double>& molarAbundances,
        const std::vector<double>& molarMasses
    ) noexcept {
        assert(molarMasses.size() == molarAbundances.size());
        assert(!molarMasses.empty());

        double totalMass = 0;
        std::vector<double> masses;
        masses.reserve(molarMasses.size());
        for (const auto [m, Y] : std::views::zip(molarMasses, molarAbundances)) {
            const double mass = m * Y;
            totalMass += mass;
            masses.push_back(mass);
        }

        assert(totalMass > 0);

        std::vector<double> massFractions;
        massFractions.reserve(masses.size());
        std::ranges::transform(
            masses,
            std::back_inserter(massFractions),
            [&totalMass](const double speciesMass) {
                const double Xi = speciesMass / totalMass;
                if (std::abs(Xi) < 1e-16 && Xi < 0) {
                    return 0.0;
                }
                return Xi;
            });

        return massFractions;
    }

    inline std::vector<double> molarMassVectorFromComposition(
        const fourdst::composition::Composition& composition
    ) {
        std::vector<double> molarMassVector;
        molarMassVector.reserve(composition.getRegisteredSymbols().size());
        for (const auto &entry: composition | std::views::values) {
            molarMassVector.push_back(entry.isotope().mass());
        }
        return molarMassVector;
    }
}