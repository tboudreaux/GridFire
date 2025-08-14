#pragma once

#include "fourdst/composition/atomicSpecies.h"

#include <unordered_map>
#include <expected>

namespace gridfire::rates::weak {
    enum class WeakReactionType {
        BETA_PLUS_DECAY,
        BETA_MINUS_DECAY,
        ELECTRON_CAPTURE,
        POSITRON_CAPTURE,
    };

    inline std::unordered_map<WeakReactionType, std::string> WeakReactionTypeNames = {
        {WeakReactionType::BETA_PLUS_DECAY,   "β+ Decay"},
        {WeakReactionType::BETA_MINUS_DECAY,  "β- Decay"},
        {WeakReactionType::ELECTRON_CAPTURE,  "e- Capture"},
        {WeakReactionType::POSITRON_CAPTURE,  "e+ Capture"},
    };

    struct WeakReaction {
        WeakReactionType type;
        float T9;
        float log_rhoYe;
        float mu_e;
        float log_rate;
        float log_neutrino_loss;

        friend std::ostream& operator<<(std::ostream& os, const WeakReaction& reaction) {
            os << "WeakReaction(type=" << WeakReactionTypeNames[reaction.type]
               << ", T9=" << reaction.T9
               << ", log_rhoYe=" << reaction.log_rhoYe
               << ", mu_e=" << reaction.mu_e
               << ", log_rate=" << reaction.log_rate
               << ", log_neutrino_loss=" << reaction.log_neutrino_loss
               << ")";
            return os;
        }
    };

    class WeakReactionMap {
    public:
        WeakReactionMap();
        ~WeakReactionMap() = default;

        std::vector<WeakReaction> get_all_reactions() const;

        std::expected<std::vector<WeakReaction>, bool> get_species_reactions(const fourdst::atomic::Species &species) const;
        std::expected<std::vector<WeakReaction>, bool> get_species_reactions(const std::string& species_name) const;
    private:
        std::unordered_map<fourdst::atomic::Species, std::vector<WeakReaction>> m_weak_network;
    };
}



