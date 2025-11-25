#include "gridfire/engine/procedures/construction.h"
#include "gridfire/reaction/weak/weak_interpolator.h"
#include "gridfire/reaction/weak/weak.h"
#include "gridfire/reaction/weak/weak_types.h"

#include <ranges>
#include <stdexcept>
#include <memory>

#include "fourdst/atomic/species.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/reaclib.h"

#include "fourdst/composition/composition_abstract.h"

#include "fourdst/logging/logging.h"

#include "quill/Logger.h"
#include "quill/LogMacros.h"
namespace {
    // Simple heuristic to check if a reaclib reaction is a strong or weak reaction
    /*  A weak reaction is defined here as one where:
        - The number of reactants is equal to the number of products
        - There is only one reactant and one product
        - The mass number (A) of the reactant is equal to the mass number (A) of the product
    */
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
        const gridfire::rates::weak::WeakRateInterpolator &weakInterpolator,
        const gridfire::engine::NetworkConstructionFlags reactionTypes
    ) {
        gridfire::reaction::ReactionSet weak_reaction_pool;
        if (!has_flag(reactionTypes, gridfire::engine::NetworkConstructionFlags::WRL_WEAK)) {
            return weak_reaction_pool;
        }

        for (const auto& parent_species: weakInterpolator.available_isotopes()) {
            std::expected<fourdst::atomic::Species, fourdst::atomic::SpeciesErrorType> upProduct = fourdst::atomic::az_to_species(
                parent_species.a(),
                parent_species.z() + 1
            );
            std::expected<fourdst::atomic::Species, fourdst::atomic::SpeciesErrorType> downProduct = fourdst::atomic::az_to_species(
                parent_species.a(),
                parent_species.z() - 1
            );
            if (downProduct.has_value()) { // Only add the reaction if the Species map contains the product
                if (has_flag(reactionTypes, gridfire::engine::NetworkConstructionFlags::WRL_BETA_PLUS)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::BETA_PLUS_DECAY,
                            weakInterpolator
                        )
                    );
                }
                if (has_flag(reactionTypes, gridfire::engine::NetworkConstructionFlags::WRL_ELECTRON_CAPTURE)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::ELECTRON_CAPTURE,
                            weakInterpolator
                        )
                    );
                }
            }
            if (upProduct.has_value()) { // Only add the reaction if the Species map contains the product
                if (has_flag(reactionTypes, gridfire::engine::NetworkConstructionFlags::WRL_BETA_MINUS)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::BETA_MINUS_DECAY,
                            weakInterpolator
                        )
                    );
                }
                if (has_flag(reactionTypes, gridfire::engine::NetworkConstructionFlags::WRL_POSITRON_CAPTURE)) {
                    weak_reaction_pool.add_reaction(
                        std::make_unique<gridfire::rates::weak::WeakReaction>(
                            parent_species,
                            gridfire::rates::weak::WeakReactionType::POSITRON_CAPTURE,
                            weakInterpolator
                        )
                    );
                }
            }
        }

        return weak_reaction_pool;
    }

    gridfire::reaction::ReactionSet register_strong_reactions(
        const gridfire::engine::NetworkConstructionFlags reaction_types
    ) {
        gridfire::reaction::ReactionSet strong_reaction_pool;
        if (has_flag(reaction_types, gridfire::engine::NetworkConstructionFlags::REACLIB_STRONG)) {
            const auto& allReaclibReactions = gridfire::reaclib::get_all_reaclib_reactions();
            for (const auto& reaction : allReaclibReactions) {
                const bool isWeakReaction = reaclib_reaction_is_weak(*reaction);
                const bool okayToUseReaclibWeakReaction = has_flag(reaction_types, gridfire::engine::NetworkConstructionFlags::REACLIB_WEAK);

                const bool reaclibWeakOkay = !isWeakReaction || okayToUseReaclibWeakReaction;
                if (!reaction->is_reverse() && reaclibWeakOkay) {
                    strong_reaction_pool.add_reaction(reaction->clone());
                }
            }
        }
        return strong_reaction_pool;
    }

    bool validate_unique_weak_set(gridfire::engine::NetworkConstructionFlags flag) {
        // This method ensures that weak reactions will only be fetched from either reaclib or the weak reaction library (WRL)
        // but not both
        const std::array<gridfire::engine::NetworkConstructionFlags, 4> WRL_Flags = {
            gridfire::engine::NetworkConstructionFlags::WRL_BETA_PLUS,
            gridfire::engine::NetworkConstructionFlags::WRL_ELECTRON_CAPTURE,
            gridfire::engine::NetworkConstructionFlags::WRL_POSITRON_CAPTURE,
            gridfire::engine::NetworkConstructionFlags::WRL_BETA_MINUS
        };

        if (!has_flag(flag, gridfire::engine::NetworkConstructionFlags::REACLIB_WEAK)) {
            return true;
        }
        for (const auto& WRLReactionType : WRL_Flags) {
            if (has_flag(flag, WRLReactionType)) {
                return false;
            }
        }
        return true;
    }
}

namespace gridfire::engine {
    using reaction::ReactionSet;
    using reaction::Reaction;
    using fourdst::atomic::Species;

    ReactionSet build_nuclear_network(
        const fourdst::composition::CompositionAbstract &composition,
        const rates::weak::WeakRateInterpolator &weakInterpolator,
        BuildDepthType maxLayers,
        NetworkConstructionFlags ReactionTypes
    ) {
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        if (!validate_unique_weak_set(ReactionTypes)) {
            std::string msg = "Cannot construct network since weak reactions from both reaclib and WRL were requested. Only one weak rate source may be used. In network construction this likely means that the flag NetworkConstructionFlags::REACLIB_WEAK was used along with one of the weak flags (BETA_PLUS, BETA_MINUS, ELECTRON_CAPTURE, POSITRON_CAPTURE). These flags are for WRL rates and may not be used in conjunction with reaclib weak rates.";
            LOG_ERROR(logger, "{}", msg);
            throw std::logic_error(msg);
        }
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
        for (const auto& sp : composition.getRegisteredSpecies()) {
            if (composition.getMolarAbundance(sp) > 0.0) {
                availableSpecies.insert(sp);
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
