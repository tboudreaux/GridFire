#include "gridfire/reaction/reaction.h"

#include<string_view>
#include<string>
#include<vector>
#include<unordered_set>
#include<algorithm>
#include <ranges>

#include "quill/LogMacros.h"

#include "fourdst/composition/atomicSpecies.h"

#include "xxhash64.h"

namespace gridfire::reaction {
    using namespace fourdst::atomic;

    Reaction::Reaction(
        const std::string_view id,
        const std::string_view peName,
        const int chapter,
        const std::vector<Species>& reactants,
        const std::vector<Species>& products,
        const double qValue,
        const std::string_view label,
        const RateCoefficientSet& sets,
        const bool reverse) :
    m_id(id),
    m_peName(peName),
    m_chapter(chapter),
    m_qValue(qValue),
    m_reactants(reactants),
    m_products(products),
    m_sourceLabel(label),
    m_rateCoefficients(sets),
    m_reverse(reverse) {}

    double Reaction::calculate_rate(const double T9) const {
        return calculate_rate<double>(T9);
    }

    CppAD::AD<double> Reaction::calculate_rate(const CppAD::AD<double> T9) const {
        return calculate_rate<CppAD::AD<double>>(T9);
    }

    bool Reaction::contains(const Species &species) const {
        return contains_reactant(species) || contains_product(species);
    }


    bool Reaction::contains_reactant(const Species& species) const {
        for (const auto& reactant : m_reactants) {
            if (reactant == species) {
                return true;
            }
        }
        return false;
    }

    bool Reaction::contains_product(const Species& species) const {
        for (const auto& product : m_products) {
            if (product == species) {
                return true;
            }
        }
        return false;
    }

    std::unordered_set<Species> Reaction::all_species() const {
        auto rs = reactant_species();
        auto ps = product_species();
        rs.insert(ps.begin(), ps.end());
        return rs;
    }

    std::unordered_set<Species> Reaction::reactant_species() const {
        std::unordered_set<Species> reactantsSet;
        for (const auto& reactant : m_reactants) {
            reactantsSet.insert(reactant);
        }
        return reactantsSet;
    }

    std::unordered_set<Species> Reaction::product_species() const {
        std::unordered_set<Species> productsSet;
        for (const auto& product : m_products) {
            productsSet.insert(product);
        }
        return productsSet;
    }

    int Reaction::stoichiometry(const Species& species) const {
        int s = 0;
        for (const auto& reactant : m_reactants) {
            if (reactant == species) {
                s--;
            }
        }
        for (const auto& product : m_products) {
            if (product == species) {
                s++;
            }
        }
        return s;
    }

    size_t Reaction::num_species() const {
        return all_species().size();
    }

    std::unordered_map<Species, int> Reaction::stoichiometry() const {
        std::unordered_map<Species, int> stoichiometryMap;
        for (const auto& reactant : m_reactants) {
            stoichiometryMap[reactant]--;
        }
        for (const auto& product : m_products) {
            stoichiometryMap[product]++;
        }
        return stoichiometryMap;
    }

    double Reaction::excess_energy() const {
        double reactantMass = 0.0;
        double productMass = 0.0;
        constexpr double AMU2MeV = 931.494893; // Conversion factor from atomic mass unit to MeV
        for (const auto& reactant : m_reactants) {
            reactantMass += reactant.mass();
        }
        for (const auto& product : m_products) {
            productMass += product.mass();
        }
        return (reactantMass - productMass) * AMU2MeV;
    }

    uint64_t Reaction::hash(uint64_t seed) const {
        return XXHash64::hash(m_id.data(), m_id.size(), seed);
    }



    LogicalReaction::LogicalReaction(const std::vector<Reaction>& reactants) :
        Reaction(reactants.front().peName(),
                 reactants.front().peName(),
                 reactants.front().chapter(),
                 reactants.front().reactants(),
                 reactants.front().products(),
                 reactants.front().qValue(),
                 reactants.front().sourceLabel(),
                 reactants.front().rateCoefficients(),
                 reactants.front().is_reverse()) {

        m_sources.reserve(reactants.size());
        m_rates.reserve(reactants.size());
        for (const auto& reaction : reactants) {
            if (std::abs(std::abs(reaction.qValue()) - std::abs(m_qValue)) > 1e-6) {
                LOG_ERROR(
                    m_logger,
                    "LogicalReaction constructed with reactions having different Q-values. Expected {} got {}.",
                    m_qValue,
                    reaction.qValue()
                );
                m_logger -> flush_log();
                throw std::runtime_error("LogicalReaction constructed with reactions having different Q-values. Expected " + std::to_string(m_qValue) + " got " + std::to_string(reaction.qValue()) + " (difference : " + std::to_string(std::abs(reaction.qValue() - m_qValue)) + ").");
            }
            m_sources.push_back(std::string(reaction.sourceLabel()));
            m_rates.push_back(reaction.rateCoefficients());
        }
    }

    void LogicalReaction::add_reaction(const Reaction& reaction) {
        if (reaction.peName() != m_id) {
            LOG_ERROR(m_logger, "Cannot add reaction with different peName to LogicalReaction. Expected {} got {}.", m_id, reaction.peName());
            m_logger -> flush_log();
            throw std::runtime_error("Cannot add reaction with different peName to LogicalReaction. Expected " + std::string(m_id) + " got " + std::string(reaction.peName()) + ".");
        }
        for (const auto& source : m_sources) {
            if (source == reaction.sourceLabel()) {
                LOG_ERROR(m_logger, "Cannot add reaction with duplicate source label {} to LogicalReaction.", reaction.sourceLabel());
                m_logger -> flush_log();
                throw std::runtime_error("Cannot add reaction with duplicate source label " + std::string(reaction.sourceLabel()) + " to LogicalReaction.");
            }
        }
        if (std::abs(reaction.qValue() - m_qValue) > 1e-6) {
            LOG_ERROR(m_logger, "LogicalReaction constructed with reactions having different Q-values. Expected {} got {}.", m_qValue, reaction.qValue());
            m_logger -> flush_log();
            throw std::runtime_error("LogicalReaction constructed with reactions having different Q-values. Expected " + std::to_string(m_qValue) + " got " + std::to_string(reaction.qValue()) + ".");
        }
        m_sources.push_back(std::string(reaction.sourceLabel()));
        m_rates.push_back(reaction.rateCoefficients());
    }

    double LogicalReaction::calculate_rate(const double T9) const {
        return calculate_rate<double>(T9);
    }

    CppAD::AD<double> LogicalReaction::calculate_rate(const CppAD::AD<double> T9) const {
        return calculate_rate<CppAD::AD<double>>(T9);
    }

    LogicalReactionSet packReactionSetToLogicalReactionSet(const ReactionSet& reactionSet) {
        std::unordered_map<std::string_view, std::vector<Reaction>> groupedReactions;

        for (const auto& reaction: reactionSet) {
            groupedReactions[reaction.peName()].push_back(reaction);
        }

        std::vector<LogicalReaction> reactions;
        reactions.reserve(groupedReactions.size());

        for (const auto &reactionsGroup: groupedReactions | std::views::values) {
            LogicalReaction logicalReaction(reactionsGroup);
            reactions.push_back(logicalReaction);
        }
        return LogicalReactionSet(std::move(reactions));
    }
}

namespace std {
    template<>
    struct hash<gridfire::reaction::Reaction> {
        size_t operator()(const gridfire::reaction::Reaction& r) const noexcept {
            return r.hash(0);
        }
    };

    template<>
    struct hash<gridfire::reaction::ReactionSet> {
        size_t operator()(const gridfire::reaction::ReactionSet& s) const noexcept {
            return s.hash(0);
        }
    };

    template<>
    struct hash<gridfire::reaction::LogicalReactionSet> {
        size_t operator()(const gridfire::reaction::LogicalReactionSet& s) const noexcept {
            return s.hash(0);
        }
    };
} // namespace std
