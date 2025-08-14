#include "gridfire/reaction/weak/weak_rate_library.h"
#include "gridfire/reaction/weak/weak.h"

#include "fourdst/composition/species.h"

#include <array>
#include <ranges>
#include <unordered_map>
#include <expected>

#define GRIDFIRE_WEAK_REACTION_LIB_SENTINEL -60.0

namespace gridfire::rates::weak {
    WeakReactionMap::WeakReactionMap() {
        using namespace fourdst::atomic;


        for (const auto& weak_reaction_record : UNIFIED_WEAK_DATA) {
            Species species = az_to_species(weak_reaction_record.A, weak_reaction_record.Z);

            if (weak_reaction_record.log_beta_minus > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReaction{
                        WeakReactionType::BETA_MINUS_DECAY,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_beta_minus,
                        weak_reaction_record.log_antineutrino_loss_bd
                    }
                );
            }
            if (weak_reaction_record.log_beta_plus > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReaction{
                        WeakReactionType::BETA_PLUS_DECAY,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_beta_plus,
                        weak_reaction_record.log_neutrino_loss_ec
                    }
                );
            }
            if (weak_reaction_record.log_electron_capture > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReaction{
                        WeakReactionType::ELECTRON_CAPTURE,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_electron_capture,
                        weak_reaction_record.log_neutrino_loss_ec
                    }
                );
            }
            if (weak_reaction_record.log_positron_capture > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReaction{
                        WeakReactionType::POSITRON_CAPTURE,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_positron_capture,
                        weak_reaction_record.log_antineutrino_loss_bd
                    }
                );
            }
        }
    }

    std::vector<WeakReaction> WeakReactionMap::get_all_reactions() const {
        std::vector<WeakReaction> reactions;
        for (const auto &species_reactions: m_weak_network | std::views::values) {
            reactions.insert(reactions.end(), species_reactions.begin(), species_reactions.end());
        }
        return reactions;
    }

    std::expected<std::vector<WeakReaction>, bool> WeakReactionMap::get_species_reactions(const fourdst::atomic::Species &species) const {
        if (m_weak_network.contains(species)) {
            return m_weak_network.at(species);
        }
        return std::unexpected(false);
    }

    std::expected<std::vector<WeakReaction>, bool> WeakReactionMap::get_species_reactions(const std::string &species_name) const {
        fourdst::atomic::Species species = fourdst::atomic::species.at(species_name);
        if (m_weak_network.contains(species)) {
            return m_weak_network.at(species);
        }
        return std::unexpected(false);
    }
}
