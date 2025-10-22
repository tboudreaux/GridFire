#pragma once
#include "fourdst/composition/composition.h"
#include "fourdst/composition/atomicSpecies.h"

namespace gridfire::utils {
    inline double massFractionFromMolarAbundance (
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
}