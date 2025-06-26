#include "gridfire/reaction/reaction.h"

#include<string_view>
#include<string>
#include<vector>
#include<memory>
#include<unordered_set>
#include<algorithm>

#include "quill/LogMacros.h"

#include "fourdst/composition/atomicSpecies.h"

#include "xxhash64.h"

namespace gridfire::reaction {
    using namespace fourdst::atomic;

    Reaction::Reaction(
        const std::string_view id,
        const double qValue,
        const std::vector<Species>& reactants,
        const std::vector<Species>& products,
        const bool reverse) :
        m_id(id),
        m_qValue(qValue),
        m_reactants(std::move(reactants)),
        m_products(std::move(products)),
        m_reverse(reverse) {}

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

    ReactionSet::ReactionSet(
        std::vector<std::unique_ptr<Reaction>> reactions) :
        m_reactions(std::move(reactions)) {
        if (m_reactions.empty()) {
            return; // Case where the reactions will be added later.
        }
        m_reactionNameMap.reserve(reactions.size());
        for (const auto& reaction : m_reactions) {
            m_id += reaction->id();
            m_reactionNameMap.emplace(reaction->id(), reaction.get());
        }
    }

    ReactionSet::ReactionSet(const ReactionSet &other) {
        m_reactions.reserve(other.m_reactions.size());
        for (const auto& reaction_ptr: other.m_reactions) {
            m_reactions.push_back(reaction_ptr->clone());
        }

        m_reactionNameMap.reserve(other.m_reactionNameMap.size());
        for (const auto& reaction_ptr : m_reactions) {
            m_reactionNameMap.emplace(reaction_ptr->id(), reaction_ptr.get());
        }
    }

    ReactionSet & ReactionSet::operator=(const ReactionSet &other) {
        if (this != &other) {
            ReactionSet temp(other);
            std::swap(m_reactions, temp.m_reactions);
            std::swap(m_reactionNameMap, temp.m_reactionNameMap);
        }
        return *this;
    }

    void ReactionSet::add_reaction(std::unique_ptr<Reaction> reaction) {
        m_reactions.emplace_back(std::move(reaction));
        m_id += m_reactions.back()->id();
        m_reactionNameMap.emplace(m_reactions.back()->id(), m_reactions.back().get());
    }

    void ReactionSet::remove_reaction(const std::unique_ptr<Reaction>& reaction) {
        if (!m_reactionNameMap.contains(std::string(reaction->id()))) {
            // LOG_INFO(m_logger, "Attempted to remove reaction {} that does not exist in ReactionSet. Skipping...", reaction->id());
            return;
        }

        m_reactionNameMap.erase(std::string(reaction->id()));

        std::erase_if(m_reactions, [&reaction](const std::unique_ptr<Reaction>& r) {
            return *r == *reaction;
        });
    }

    bool ReactionSet::contains(const std::string_view& id) const {
        for (const auto& reaction : m_reactions) {
            if (reaction->id() == id) {
                return true;
            }
        }
        return false;
    }

    bool ReactionSet::contains(const Reaction& reaction) const {
        for (const auto& r : m_reactions) {
            if (*r == reaction) {
                return true;
            }
        }
        return false;
    }

    void ReactionSet::sort(double T9) {
        std::ranges::sort(m_reactions,
                          [&T9](const std::unique_ptr<Reaction>& r1, const std::unique_ptr<Reaction>& r2) {
                              return r1->calculate_rate(T9) < r2->calculate_rate(T9);
                          });
    }

    bool ReactionSet::contains_species(const Species& species) const {
        for (const auto& reaction : m_reactions) {
            if (reaction->contains(species)) {
                return true;
            }
        }
        return false;
    }

    bool ReactionSet::contains_reactant(const Species& species) const {
        for (const auto& r : m_reactions) {
            if (r->contains_reactant(species)) {
                return true;
            }
        }
        return false;
    }

    bool ReactionSet::contains_product(const Species& species) const {
        for (const auto& r : m_reactions) {
            if (r->contains_product(species)) {
                return true;
            }
        }
        return false;
    }

    const Reaction& ReactionSet::operator[](const size_t index) const {
        if (index >= m_reactions.size()) {
            throw std::out_of_range("Index" + std::to_string(index) + " out of range for ReactionSet of size " + std::to_string(m_reactions.size()) + ".");
        }
        return *m_reactions[index];
    }

    const Reaction& ReactionSet::operator[](const std::string_view& id) const {
        if (auto it = m_reactionNameMap.find(std::string(id)); it != m_reactionNameMap.end()) {
            return *it->second;
        }
        throw std::out_of_range("Species " + std::string(id) + " does not exist in ReactionSet.");
    }

    bool ReactionSet::operator==(const ReactionSet& other) const {
        if (size() != other.size()) {
            return false;
        }
        return hash() == other.hash();
    }

    bool ReactionSet::operator!=(const ReactionSet& other) const {
        return !(*this == other);
    }

    uint64_t ReactionSet::hash(uint64_t seed) const {
        if (m_reactions.empty()) {
            return XXHash64::hash(nullptr, 0, seed);
        }
        std::vector<uint64_t> individualReactionHashes;
        individualReactionHashes.reserve(m_reactions.size());
        for (const auto& reaction : m_reactions) {
            individualReactionHashes.push_back(reaction->hash(seed));
        }

        std::ranges::sort(individualReactionHashes);

        const void* data = static_cast<const void*>(individualReactionHashes.data());
        size_t sizeInBytes = individualReactionHashes.size() * sizeof(uint64_t);
        return XXHash64::hash(data, sizeInBytes, seed);
    }

    REACLIBReaction::REACLIBReaction(
        const std::string_view id,
        const std::string_view peName,
        const int chapter,
        const std::vector<Species> &reactants,
        const std::vector<Species> &products,
        const double qValue,
        const std::string_view label,
        const REACLIBRateCoefficientSet &sets,
        const bool reverse) :
        Reaction(id, qValue, reactants, products, reverse),
        m_peName(peName),
        m_chapter(chapter),
        m_sourceLabel(label),
        m_rateCoefficients(sets) {}

    std::unique_ptr<Reaction> REACLIBReaction::clone() const {
        return std::make_unique<REACLIBReaction>(*this);
    }


    double REACLIBReaction::calculate_rate(const double T9) const {
        return calculate_rate<double>(T9);
    }

    CppAD::AD<double> REACLIBReaction::calculate_rate(const CppAD::AD<double> T9) const {
        return calculate_rate<CppAD::AD<double>>(T9);
    }

    REACLIBReactionSet::REACLIBReactionSet(std::vector<REACLIBReaction> reactions) :
        ReactionSet(std::vector<std::unique_ptr<Reaction>>()) {
        // Convert REACLIBReaction to unique_ptr<Reaction> and store in m_reactions
        m_reactions.reserve(reactions.size());
        m_reactionNameMap.reserve(reactions.size());
        for (auto& reaction : reactions) {
            m_reactions.emplace_back(std::make_unique<REACLIBReaction>(std::move(reaction)));
            m_reactionNameMap.emplace(std::string(reaction.id()), m_reactions.back().get());
        }
    }

    std::unordered_set<std::string> REACLIBReactionSet::peNames() const {
        std::unordered_set<std::string> peNames;
        for (const auto& reactionPtr: m_reactions) {
            if (const auto* reaction = dynamic_cast<REACLIBReaction*>(reactionPtr.get())) {
                peNames.insert(std::string(reaction->peName()));
            } else {
                // LOG_ERROR(m_logger, "Failed to cast Reaction to REACLIBReaction in REACLIBReactionSet::peNames().");
                throw std::runtime_error("Failed to cast Reaction to REACLIBReaction in REACLIBReactionSet::peNames(). This should not happen, please check your reaction setup.");
            }
        }
        return peNames;
    }

    REACLIBLogicalReaction::REACLIBLogicalReaction(const std::vector<REACLIBReaction>& reactants) :
        Reaction(reactants.front().peName(),
                 reactants.front().qValue(),
                 reactants.front().reactants(),
                 reactants.front().products(),
                 reactants.front().is_reverse()),
        m_chapter(reactants.front().chapter()) {

        m_sources.reserve(reactants.size());
        m_rates.reserve(reactants.size());
        for (const auto& reaction : reactants) {
            if (std::abs(reaction.qValue() - m_qValue) > 1e-6) {
                LOG_ERROR(m_logger, "REACLIBLogicalReaction constructed with reactions having different Q-values. Expected {} got {}.", m_qValue, reaction.qValue());
                throw std::runtime_error("REACLIBLogicalReaction constructed with reactions having different Q-values. Expected " + std::to_string(m_qValue) + " got " + std::to_string(reaction.qValue()) + ".");
            }
            m_sources.push_back(std::string(reaction.sourceLabel()));
            m_rates.push_back(reaction.rateCoefficients());
        }
    }

    REACLIBLogicalReaction::REACLIBLogicalReaction(const REACLIBReaction& reaction) :
        Reaction(reaction.peName(),
                 reaction.qValue(),
                 reaction.reactants(),
                 reaction.products(),
                 reaction.is_reverse()),
        m_chapter(reaction.chapter()) {
        m_sources.reserve(1);
        m_rates.reserve(1);
        m_sources.push_back(std::string(reaction.sourceLabel()));
        m_rates.push_back(reaction.rateCoefficients());
    }

    void REACLIBLogicalReaction::add_reaction(const REACLIBReaction& reaction) {
        if (reaction.peName() != m_id) {
            LOG_ERROR(m_logger, "Cannot add reaction with different peName to REACLIBLogicalReaction. Expected {} got {}.", m_id, reaction.peName());
            throw std::runtime_error("Cannot add reaction with different peName to REACLIBLogicalReaction. Expected " + std::string(m_id) + " got " + std::string(reaction.peName()) + ".");
        }
        for (const auto& source : m_sources) {
            if (source == reaction.sourceLabel()) {
                LOG_ERROR(m_logger, "Cannot add reaction with duplicate source label {} to REACLIBLogicalReaction.", reaction.sourceLabel());
                throw std::runtime_error("Cannot add reaction with duplicate source label " + std::string(reaction.sourceLabel()) + " to REACLIBLogicalReaction.");
            }
        }
        if (std::abs(reaction.qValue() - m_qValue) > 1e-6) {
            LOG_ERROR(m_logger, "REACLIBLogicalReaction constructed with reactions having different Q-values. Expected {} got {}.", m_qValue, reaction.qValue());
            throw std::runtime_error("REACLIBLogicalReaction constructed with reactions having different Q-values. Expected " + std::to_string(m_qValue) + " got " + std::to_string(reaction.qValue()) + ".");
        }
        m_sources.push_back(std::string(reaction.sourceLabel()));
        m_rates.push_back(reaction.rateCoefficients());
    }

    std::unique_ptr<Reaction> REACLIBLogicalReaction::clone() const {
        return std::make_unique<REACLIBLogicalReaction>(*this);
    }

    double REACLIBLogicalReaction::calculate_rate(const double T9) const {
        return calculate_rate<double>(T9);
    }

    CppAD::AD<double> REACLIBLogicalReaction::calculate_rate(const CppAD::AD<double> T9) const {
        return calculate_rate<CppAD::AD<double>>(T9);
    }

    REACLIBLogicalReactionSet::REACLIBLogicalReactionSet(const REACLIBReactionSet &reactionSet) :
        ReactionSet(std::vector<std::unique_ptr<Reaction>>()) {

        std::unordered_map<std::string_view, std::vector<REACLIBReaction>> grouped_reactions;

        for (const auto& reaction_ptr : reactionSet) {
            if (const auto* reaclib_reaction = dynamic_cast<const REACLIBReaction*>(reaction_ptr.get())) {
                grouped_reactions[reaclib_reaction->peName()].push_back(*reaclib_reaction);
            }
        }
        m_reactions.reserve(grouped_reactions.size());
        m_reactionNameMap.reserve(grouped_reactions.size());
        for (const auto& [peName, reactions_for_peName] : grouped_reactions) {
            m_peNames.insert(std::string(peName));
            auto logical_reaction = std::make_unique<REACLIBLogicalReaction>(reactions_for_peName);
            m_reactionNameMap.emplace(logical_reaction->id(), logical_reaction.get());
            m_reactions.push_back(std::move(logical_reaction));
        }
    }

    std::unordered_set<std::string> REACLIBLogicalReactionSet::peNames() const {
        return m_peNames;
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
} // namespace std
