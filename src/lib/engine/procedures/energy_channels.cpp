#include "gridfire/engine/procedures/energy_channels.h"

#include <algorithm>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
    constexpr size_t UNLABELED = std::numeric_limits<size_t>::max();

    const fourdst::atomic::Species& lightest_reactant(const gridfire::reaction::Reaction& reaction) {
        const auto& reactants = reaction.reactants();
        return *std::ranges::min_element(reactants, [](const auto& a, const auto& b) {
            return a.a() != b.a() ? a.a() < b.a() : a.z() < b.z();
        });
    }

    size_t label_of(
        const std::unordered_map<fourdst::atomic::Species, size_t>& labels,
        const fourdst::atomic::Species& species
    ) {
        const auto it = labels.find(species);
        return it == labels.end() ? UNLABELED : it->second;
    }

    bool involves_heavy_nucleus(const gridfire::reaction::Reaction& reaction) {
        auto heavy = [](const fourdst::atomic::Species& s) { return s.z() >= 6; };
        return std::ranges::any_of(reaction.reactants(), heavy) || std::ranges::any_of(reaction.products(), heavy);
    }
}

namespace gridfire::engine {
    EnergyChannels classify_energy_channels(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& fuels,
        const bool splitHydrogenPPCNO
    ) {
        EnergyChannels result;

        // --- Channel names: one per fuel stage (H-1 optionally split into pp and CNO), then "other" ---
        std::vector<size_t> stage_channel(fuels.size());
        size_t cno_channel = UNLABELED;
        for (size_t stage = 0; stage < fuels.size(); ++stage) {
            const auto& fuel = fuels[stage];
            const bool element_is_unique = std::ranges::count_if(fuels, [&](const auto& f) { return f.z() == fuel.z(); }) == 1;
            const std::string name = element_is_unique ? std::string(fuel.el()) : std::string(fuel.name());
            stage_channel[stage] = result.names.size();
            if (splitHydrogenPPCNO && fuel.z() == 1 && fuel.a() == 1) {
                result.names.push_back(name + "_pp");
                cno_channel = result.names.size();
                result.names.push_back(name + "_CNO");
            } else {
                result.names.push_back(name);
            }
        }
        const size_t other_channel = result.names.size();
        result.names.emplace_back("other");

        // --- Species lineage: fuels keep their stage, everything else inherits the lowest producing stage ---
        std::unordered_map<fourdst::atomic::Species, size_t> labels;
        for (size_t stage = 0; stage < fuels.size(); ++stage) {
            labels.emplace(fuels[stage], stage);
        }

        bool changed = true;
        while (changed) { // Labels only decrease and are bounded below, so this terminates
            changed = false;
            for (const auto& reaction : reactions) {
                const size_t stage = label_of(labels, lightest_reactant(*reaction));
                if (stage == UNLABELED) continue;
                for (const auto& product : reaction->products()) {
                    if (std::ranges::contains(fuels, product)) continue;
                    if (stage < label_of(labels, product)) {
                        labels[product] = stage;
                        changed = true;
                    }
                }
            }
        }

        // --- Reaction channels from the lightest reactant ---
        result.reaction_channel.reserve(reactions.size());
        for (const auto& reaction : reactions) {
            const size_t stage = label_of(labels, lightest_reactant(*reaction));
            if (stage == UNLABELED) {
                result.reaction_channel.push_back(other_channel);
            } else if (stage_channel[stage] + 1 == cno_channel && involves_heavy_nucleus(*reaction)) {
                result.reaction_channel.push_back(cno_channel);
            } else {
                result.reaction_channel.push_back(stage_channel[stage]);
            }
        }

        return result;
    }
}
