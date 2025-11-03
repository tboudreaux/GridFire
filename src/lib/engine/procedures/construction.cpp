#include "gridfire/engine/procedures/construction.h"
#include "gridfire/reaction/weak/weak_interpolator.h"
#include "gridfire/reaction/weak/weak.h"
#include "gridfire/reaction/weak/weak_types.h"

#include <ranges>
#include <stdexcept>
#include <memory>
#include <cmath>

#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/reaclib.h"

#include "fourdst/composition/composition.h"

#include "fourdst/logging/logging.h"

#include "quill/Logger.h"
#include "quill/LogMacros.h"
namespace {
    // Simple heuristic to check if a reaclib reaction is a strong or weak reaction
    bool reaclib_reaction_is_weak(const gridfire::reaction::Reaction& reaction) {
        const std::vector<fourdst::atomic::Species>& reactants = reaction.reactants();
        const std::vector<fourdst::atomic::Species>& products = reaction.products();

        if (reactants.size() != products.size()) {
            return false;
        }

        if (reactants.size() != 1 || products.size() != 1) {
            return false;
        }

        if (std::floor(reactants[0].a()) != std::floor(products[0].a())) {
            return false;
        }

        return true;
    }


    gridfire::reaction::ReactionSet register_weak_reactions(
        const std::optional<gridfire::rates::weak::WeakRateInterpolator> &weakInterpolator,
        const gridfire::NetworkConstructionFlags reactionTypes
    ) {
        gridfire::reaction::ReactionSet weak_reaction_pool;
        assert(weakInterpolator.has_value());
        if (!has_flag(reactionTypes, gridfire::NetworkConstructionFlags::WEAK)) {
            return weak_reaction_pool;
        }

        for (const auto& parent_species: weakInterpolator->available_isotopes()) {
            std::expected<fourdst::atomic::Species, fourdst::atomic::SpeciesErrorType> upProduct = fourdst::atomic::az_to_species(
                parent_species.a(),
                parent_species.z() + 1
            );
            std::expected<fourdst::atomic::Species, fourdst::atomic::SpeciesErrorType> downProduct = fourdst::atomic::az_to_species(
                parent_species.a(),
                parent_species.z() - 1
            );
            if (downProduct.has_value()) { // Only add the reaction if the Species map contains the product
                if (has_flag(reactionTypes, gridfire::NetworkConstructionFlags::BETA_PLUS)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::BETA_PLUS_DECAY,
                            *weakInterpolator
                        )
                    );
                }
                if (has_flag(reactionTypes, gridfire::NetworkConstructionFlags::ELECTRON_CAPTURE)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::ELECTRON_CAPTURE,
                            *weakInterpolator
                        )
                    );
                }
            }
            if (upProduct.has_value()) { // Only add the reaction if the Species map contains the product
                if (has_flag(reactionTypes, gridfire::NetworkConstructionFlags::BETA_MINUS)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::BETA_MINUS_DECAY,
                            *weakInterpolator
                        )
                    );
                }
                if (has_flag(reactionTypes, gridfire::NetworkConstructionFlags::POSITRON_CAPTURE)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::POSITRON_CAPTURE,
                            *weakInterpolator
                        )
                    );
                }
            }
        }

        return weak_reaction_pool;
    }

    gridfire::reaction::ReactionSet register_strong_reactions(
        const gridfire::NetworkConstructionFlags reaction_types
    ) {
        gridfire::reaction::ReactionSet strong_reaction_pool;
        if (has_flag(reaction_types, gridfire::NetworkConstructionFlags::STRONG)) {
            const auto& allReaclibReactions = gridfire::reaclib::get_all_reaclib_reactions();
            for (const auto& reaction : allReaclibReactions) {
                if (!reaction->is_reverse() && !reaclib_reaction_is_weak(*reaction)) { // Only add reactions of the correct direction and which are not weak. Weak reactions are handled from the WRL separately which provides much higher quality weak reactions than reaclib does
                    strong_reaction_pool.add_reaction(reaction->clone());
                }
            }
        }
        return strong_reaction_pool;
    }
}

namespace gridfire {
    using reaction::ReactionSet;
    using reaction::Reaction;
    using fourdst::composition::Composition;
    using fourdst::atomic::Species;

    ReactionSet build_nuclear_network(
        const Composition& composition,
        const std::optional<rates::weak::WeakRateInterpolator> &weakInterpolator,
        BuildDepthType maxLayers,
        NetworkConstructionFlags ReactionTypes
    ) {
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        LOG_INFO(logger, "Constructing network topology from reaction types : {}", NetworkConstructionFlagsToString(ReactionTypes));

        if (ReactionTypes == NetworkConstructionFlags::NONE) {
            LOG_ERROR(logger, "Reaction types is set to NONE. No reactions will be collected");
            throw std::logic_error("Reaction types is set to NONE. No reactions will be collected");
        }

        int depth;
        if (std::holds_alternative<NetworkBuildDepth>(maxLayers)) {
            depth = static_cast<int>(std::get<NetworkBuildDepth>(maxLayers));
        } else {
            depth = std::get<int>(maxLayers);
        }
        if (depth == 0) {
            LOG_ERROR(logger, "Network build depth is set to 0. No reactions will be collected.");
            throw std::logic_error("Network build depth is set to 0. No reactions will be collected.");
        }

        // --- Step 1: Create a single master pool that owns all possible reactions ---
        ReactionSet master_reaction_pool;

        // Clone all relevant REACLIB reactions into the master pool
        master_reaction_pool.extend(register_strong_reactions(ReactionTypes));

        // --- Clone all possible weak reactions into the master reaction pool if the construction function is told to use weak reactions ---
        master_reaction_pool.extend(register_weak_reactions(weakInterpolator, ReactionTypes));

        // --- Step 2: Use non-owning raw pointers for the fast build algorithm ---
        std::vector<Reaction*> remainingReactions;
        remainingReactions.reserve(master_reaction_pool.size());
        for(const auto& reaction : master_reaction_pool) {
            remainingReactions.push_back(reaction.get());
        }

        if (depth == static_cast<int>(NetworkBuildDepth::Full)) {
            LOG_INFO(logger, "Building full nuclear network with a total of {} reactions.", remainingReactions.size());
            return master_reaction_pool; // No need to build layer by layer if we want the full network
        }

        // --- Step 3: Execute the layered network build using observing pointers ---
        std::unordered_set<Species> availableSpecies;
        for (const auto& entry : composition | std::views::values) {
            if (entry.mass_fraction() > 0.0) {
                availableSpecies.insert(entry.isotope());
            }
        }

        std::vector<Reaction*> collectedReactionPtrs;
        LOG_INFO(logger, "Starting network construction with {} available species.", availableSpecies.size());

        for (int layer = 0; layer < depth && !remainingReactions.empty(); ++layer) {
            [[maybe_unused]] size_t collectedThisLayer = 0;
            [[maybe_unused]] size_t collectedStrong = 0;
            [[maybe_unused]] size_t collectedWeak = 0;
            LOG_TRACE_L1(logger, "Collecting reactions for layer {} with {} remaining reactions. Currently there are {} available species", layer, remainingReactions.size(), availableSpecies.size());
            std::vector<Reaction*> reactionsForNextPass;
            std::unordered_set<Species> newProductsThisLayer;
            bool newReactionsAdded = false;

            reactionsForNextPass.reserve(remainingReactions.size());

            for (Reaction* reaction : remainingReactions) {
                bool allReactantsAvailable = true;
                for (const auto& reactant : reaction->reactants()) {
                    if (!availableSpecies.contains(reactant)) {
                        allReactantsAvailable = false;
                        break;
                    }
                }

                if (allReactantsAvailable) {
                    collectedReactionPtrs.push_back(reaction);
                    if (reaction->type() == reaction::ReactionType::WEAK) {
                        collectedWeak++;
                    } else {
                        collectedStrong++;
                    }
                    collectedThisLayer++;
                    newReactionsAdded = true;

                    for (const auto& product : reaction->products()) {
                        newProductsThisLayer.insert(product);
                    }
                } else {
                    reactionsForNextPass.push_back(reaction);
                }
            }

            if (!newReactionsAdded) {
                LOG_INFO(logger, "No new reactions added in layer {}. Stopping network construction.", layer);
                break;
            }

            [[maybe_unused]] size_t oldProductCount = availableSpecies.size();
            availableSpecies.insert(newProductsThisLayer.begin(), newProductsThisLayer.end());
            [[maybe_unused]] size_t newProductCount = availableSpecies.size() - oldProductCount;
            LOG_TRACE_L1(
                logger,
                "Layer {}: Collected {} new reactions ({} strong, {} weak). New products this layer: {}",
                layer,
                collectedThisLayer,
                collectedStrong,
                collectedWeak,
                newProductCount
            );
            remainingReactions = std::move(reactionsForNextPass);
        }

        // --- Step 4: Construct the final ReactionSet by moving ownership ---
        LOG_INFO(logger, "Network construction completed. Assembling final set of {} reactions.", collectedReactionPtrs.size());
        ReactionSet finalReactionSet;

        std::unordered_set<Reaction*> collectedPtrSet(collectedReactionPtrs.begin(), collectedReactionPtrs.end());

        for (auto& reaction_ptr : master_reaction_pool) {
            if (reaction_ptr && collectedPtrSet.contains(reaction_ptr.get())) {
                finalReactionSet.add_reaction(std::move(reaction_ptr));
            }
        }

        return finalReactionSet;
    }


}