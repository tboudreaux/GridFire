/**
 * @file energy_channels.h
 * @brief Automatic classification of reactions into energy generation channels (H burning, He burning, ...).
 */
#pragma once

#include "gridfire/reaction/reaction.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <string>
#include <vector>

namespace gridfire::engine {
    /**
     * @brief Partition of a reaction network into named energy generation channels.
     *
     * Every reaction is assigned to exactly one channel, so summing the energy released by each channel recovers the
     * total nuclear energy generation rate exactly.
     */
    struct EnergyChannels {
        std::vector<std::string> names; ///< Channel names, in output order. The last channel is always "other".
        std::vector<size_t> reaction_channel; ///< Channel index for each reaction, indexed like the input ReactionSet.
    };

    /**
     * @brief Classify reactions into energy channels from an ordered list of fuels.
     *
     * The user supplies fuels (e.g. H-1, He-4, C-12, Ne-20, O-16, Si-28), ordered by burning stage; no reactions are
     * listed by hand. The classification is purely topological and runs once per network:
     *
     * 1. Each fuel is labelled with its own stage. Every other species inherits the lowest stage of the reactions that
     *    produce it (a fixed point over the reaction graph), so intermediates such as He-3, Be-7, or N-13 are labelled
     *    as hydrogen burning products.
     * 2. Each reaction takes the label of its lightest reactant (the "projectile"). This places He-3 + He-4 in the pp
     *    chain while placing N-14 + He-4 in helium burning.
     * 3. Optionally, hydrogen burning is split into pp and CNO: a hydrogen-stage reaction involving any nucleus with
     *    Z >= 6 is CNO (this also captures the NeNa and MgAl cycles).
     *
     * Reactions whose lightest reactant has no label (no fuel lineage) are placed in "other".
     *
     * @note During advanced burning, alpha captures are labelled as helium burning even when the alphas come from
     *       photodisintegration (e.g. Ne-20(a,g)Mg-24 during neon burning).
     *
     * @param reactions Reactions to classify.
     * @param fuels Fuel species ordered by burning stage.
     * @param splitHydrogenPPCNO Split the H-1 stage into "H_pp" and "H_CNO".
     */
    EnergyChannels classify_energy_channels(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& fuels,
        bool splitHydrogenPPCNO
    );
}
