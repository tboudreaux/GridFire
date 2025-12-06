#include "gridfire/reaction/reaction.h"

#include<string_view>
#include<string>
#include <utility>
#include<vector>
#include<unordered_set>
#include<algorithm>
#include <ranges>

#include "quill/LogMacros.h"

#include "fourdst/atomic/atomicSpecies.h"

#include "xxhash64.h"
#include "gridfire/exceptions/exceptions.h"
#include "gridfire/reaction/reaclib.h"

namespace {
    std::string_view safe_check_reactant_id(const std::vector<std::unique_ptr<gridfire::reaction::ReaclibReaction>>& reactions) {
        if (reactions.empty()) {
            throw std::runtime_error("No reactions found in the REACLIB reaction set.");
        }
        return reactions.front()->peName();
    }
}

namespace gridfire::reaction {
    using namespace fourdst::atomic;

    ReaclibReaction::ReaclibReaction(
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
    m_sourceLabel(label),
    m_rateCoefficients(sets),
    m_reverse(reverse) {
        for (const auto& reactant : reactants) {
            m_reactants[reactant]++;
        }
        for (const auto& product : products) {
            m_products[product]++;
        }
    }

    double ReaclibReaction::calculate_rate(
        const double T9,
        const double rho,
        double Ye,
        double mue,
        const std::vector<double> &Y,
        const std::unordered_map<size_t, Species>& index_to_species_map
    ) const {
        return calculate_rate<double>(T9);
    }

    CppAD::AD<double> ReaclibReaction::calculate_rate(
        const CppAD::AD<double> T9,
        const CppAD::AD<double> rho,
        CppAD::AD<double> Ye,
        CppAD::AD<double> mue,
        const std::vector<CppAD::AD<double>>& Y,
        const std::unordered_map<size_t, Species>& index_to_species_map
    ) const {
        return calculate_rate<CppAD::AD<double>>(T9);
    }

    double ReaclibReaction::calculate_log_rate_partial_deriv_wrt_T9(
        const double T9,
        const double rho,
        double Ye,
        double mue,
        const fourdst::composition::Composition& comp
    ) const {
        constexpr double r_p13 = 1.0 / 3.0;
        constexpr double r_p53 = 5.0 / 3.0;
        constexpr double r_p23 = 2.0 / 3.0;
        constexpr double r_p43 = 4.0 / 3.0;

        const double T9_m1 = 1.0 / T9;
        const double T9_m23 = std::pow(T9, -r_p23);
        const double T9_m43 = std::pow(T9, -r_p43);

        const double d_log_k_fwd_dT9 =
            -m_rateCoefficients.a1 * T9_m1 * T9_m1
            - r_p13 * m_rateCoefficients.a2 * T9_m43
            + r_p13 * m_rateCoefficients.a3 * T9_m23
            + m_rateCoefficients.a4
            + r_p53 * m_rateCoefficients.a5 * std::pow(T9, r_p23)
            + m_rateCoefficients.a6 * T9_m1;

        return d_log_k_fwd_dT9; // Return the derivative of the log rate with respect to T9
    }

    std::optional<std::vector<RateCoefficientSet>> ReaclibReaction::getRateCoefficients() const {
        std::vector<RateCoefficientSet> rateCoefficients;
        rateCoefficients.push_back(m_rateCoefficients);
        return rateCoefficients;
    }

    bool ReaclibReaction::contains(const Species &species) const {
        return contains_reactant(species) || contains_product(species);
    }


    bool ReaclibReaction::contains_reactant(const Species& species) const {
        return m_reactants.contains(species);
    }

    bool ReaclibReaction::contains_product(const Species& species) const {
        return m_products.contains(species);
    }

    std::unordered_set<Species> ReaclibReaction::all_species() const {
        auto rs = reactant_species();
        auto ps = product_species();
        rs.insert(ps.begin(), ps.end());
        return rs;
    }

    std::unordered_set<Species> ReaclibReaction::reactant_species() const {
        std::unordered_set<Species> reactantsSet;
        for (const auto& reactant : m_reactants | std::views::keys) {
            reactantsSet.insert(reactant);
        }
        return reactantsSet;
    }

    std::unordered_set<Species> ReaclibReaction::product_species() const {
        std::unordered_set<Species> productsSet;
        for (const auto& product : m_products | std::views::keys) {
            productsSet.insert(product);
        }
        return productsSet;
    }

    int ReaclibReaction::stoichiometry(const Species& species) const {
        int s = 0;
        for (const auto& [reactant, count] : m_reactants) {
            if (reactant == species) {
                s -= count;
            }
        }
        for (const auto& [product, count] : m_products) {
            if (product == species) {
                s += count;
            }
        }
        return s;
    }

    size_t ReaclibReaction::num_species() const {
        return all_species().size();
    }

    std::unordered_map<Species, int> ReaclibReaction::stoichiometry() const {
        std::unordered_map<Species, int> stoichiometryMap;
        for (const auto& sp : all_species()) {
            stoichiometryMap[sp] = stoichiometry(sp);
        }
        return stoichiometryMap;
    }

    double ReaclibReaction::excess_energy() const {
        double reactantMass = 0.0;
        double productMass = 0.0;
        constexpr double AMU2MeV = 931.494893; // Conversion factor from atomic mass unit to MeV
        for (const auto& reactant : m_reactants | std::views::keys) {
            reactantMass += reactant.mass();
        }
        for (const auto& product : m_products | std::views::keys) {
            productMass += product.mass();
        }
        return (reactantMass - productMass) * AMU2MeV;
    }

    uint64_t ReaclibReaction::hash(const uint64_t seed) const {
        return XXHash64::hash(m_id.data(), m_id.size(), seed);
    }

    std::unique_ptr<Reaction> ReaclibReaction::clone() const {
        return std::make_unique<ReaclibReaction>(*this);
    }


    LogicalReaclibReaction::LogicalReaclibReaction(
        const std::vector<std::unique_ptr<ReaclibReaction>>& reactions
    ) : LogicalReaclibReaction(reactions, false) {}

    LogicalReaclibReaction::LogicalReaclibReaction(
        const std::vector<std::unique_ptr<ReaclibReaction>> &reactions,
        const bool reverse
    ) :
    ReaclibReaction(
        safe_check_reactant_id(reactions),
        reactions.front()->peName(),
        reactions.front()->chapter(),
        reactions.front()->reactants(),
        reactions.front()->products(),
        reactions.front()->qValue(),
        reactions.front()->sourceLabel(),
        reactions.front()->rateCoefficients(),
        reverse)
    {
        m_sources.reserve(reactions.size());
        m_rates.reserve(reactions.size());
        for (const auto& reaction : reactions) {
            if (std::abs(reaction->qValue() - m_qValue) > 1e-6) {
                LOG_ERROR(
                    m_logger,
                    "LogicalReaclibReaction constructed with reactions having different Q-values. Expected {} got {}.",
                    m_qValue,
                    reaction->qValue()
                );
                m_logger -> flush_log();
                throw std::runtime_error("LogicalReaclibReaction constructed with reactions having different Q-values. Expected " + std::to_string(m_qValue) + " got " + std::to_string(reaction->qValue()) + " (difference : " + std::to_string(std::abs(reaction->qValue() - m_qValue)) + ").");
            }
            m_sources.emplace_back(reaction->sourceLabel());
            m_rates.push_back(reaction->rateCoefficients());
        }

        std::unordered_set<bool> reaction_weak_types;
        for (const auto& reaction : reactions) {
            reaction_weak_types.insert(reaction::reaction_is_weak(*reaction));
        }

        if (reaction_weak_types.size() != 1) {
            LOG_ERROR(
                m_logger,
                "LogicalReaclibReaction constructed with reactions of mixed weak/strong types. Each LogicalReaclibReaction must contain only weak or only strong reactions."
            );

            throw exceptions::ReactionError(
                "LogicalReaclibReaction constructed with reactions of mixed weak/strong types. Each LogicalReaclibReaction must contain only weak or only strong reactions.",
                m_id
            );
        }

        m_weak = *reaction_weak_types.begin();
    }


    void LogicalReaclibReaction::add_reaction(const ReaclibReaction& reaction) {
        if (reaction.peName() != m_id) {
            LOG_ERROR(m_logger, "Cannot add reaction with different peName to LogicalReaclibReaction. Expected {} got {}.", m_id, reaction.peName());
            m_logger -> flush_log();
            throw std::runtime_error("Cannot add reaction with different peName to LogicalReaclibReaction. Expected " + std::string(m_id) + " got " + std::string(reaction.peName()) + ".");
        }
        for (const auto& source : m_sources) {
            if (source == reaction.sourceLabel()) {
                LOG_ERROR(m_logger, "Cannot add reaction with duplicate source label {} to LogicalReaclibReaction.", reaction.sourceLabel());
                m_logger -> flush_log();
                throw std::runtime_error("Cannot add reaction with duplicate source label " + std::string(reaction.sourceLabel()) + " to LogicalReaclibReaction.");
            }
        }
        if (std::abs(reaction.qValue() - m_qValue) > 1e-6) {
            LOG_ERROR(m_logger, "LogicalReaclibReaction constructed with reactions having different Q-values. Expected {} got {}.", m_qValue, reaction.qValue());
            m_logger -> flush_log();
            throw std::runtime_error("LogicalReaclibReaction constructed with reactions having different Q-values. Expected " + std::to_string(m_qValue) + " got " + std::to_string(reaction.qValue()) + ".");
        }
        m_sources.emplace_back(reaction.sourceLabel());
        m_rates.push_back(reaction.rateCoefficients());
    }

    double LogicalReaclibReaction::calculate_rate(
        const double T9,
        const double rho,
        double Ye,
        double mue, const std::vector<double> &Y, const std::unordered_map<size_t, Species>& index_to_species_map
    ) const {
        if (m_cached_rates.contains(T9)) {
            return m_cached_rates.at(T9);
        }
        const double rate = calculate_rate<double>(T9);
        m_cached_rates[T9] = rate;
        return rate;
    }

    double LogicalReaclibReaction::calculate_log_rate_partial_deriv_wrt_T9(
        const double T9,
        const double rho,
        double Ye, double mue, const fourdst::composition::Composition& comp
    ) const {
        constexpr double r_p13 = 1.0 / 3.0;
        constexpr double r_p53 = 5.0 / 3.0;
        constexpr double r_p23 = 2.0 / 3.0;
        constexpr double r_p43 = 4.0 / 3.0;

        double totalRate = 0.0;
        double totalRateDerivative = 0.0;


        const double T9_m1 = 1.0 / T9;
        const double T913 = std::pow(T9, r_p13);
        const double T953 = std::pow(T9, r_p53);
        const double logT9 = std::log(T9);

        const double T9_m1_sq = T9_m1 * T9_m1;
        const double T9_m23 = std::pow(T9, -r_p23);
        const double T9_m43 = std::pow(T9, -r_p43);
        const double T9_p23 = std::pow(T9, r_p23);


        // ReSharper disable once CppUseStructuredBinding
        for (const auto& coeffs : m_rates) {
            const double exponent = coeffs.a0 +
                                    coeffs.a1 * T9_m1 +
                                    coeffs.a2 / T913 +
                                    coeffs.a3 * T913 +
                                    coeffs.a4 * T9 +
                                    coeffs.a5 * T953 +
                                    coeffs.a6 * logT9;
            const double individualRate = std::exp(exponent);

            const double d_exponent_T9 =
                -coeffs.a1 * T9_m1_sq
                - r_p13 * coeffs.a2 * T9_m43
                + r_p13 * coeffs.a3 * T9_m23
                + coeffs.a4
                + r_p53 * coeffs.a5 * T9_p23
                + coeffs.a6 * T9_m1;

            const double individualRateDerivative = individualRate * d_exponent_T9;
            totalRate += individualRate;
            totalRateDerivative += individualRateDerivative;
        }

        if (totalRate == 0.0) {
            return 0.0; // Avoid division by zero
        }

        return totalRateDerivative / totalRate;
    }

    std::unique_ptr<Reaction> LogicalReaclibReaction::clone() const {
        return std::make_unique<LogicalReaclibReaction>(*this);
    }

    CppAD::AD<double> LogicalReaclibReaction::calculate_rate(
        const CppAD::AD<double> T9,
        const CppAD::AD<double> rho,
        CppAD::AD<double> Ye,
        CppAD::AD<double> mue, const std::vector<CppAD::AD<double>>& Y, const std::unordered_map<size_t, Species>& index_to_species_map
    ) const {
        return calculate_rate<CppAD::AD<double>>(T9);
    }

    std::optional<std::vector<RateCoefficientSet>> LogicalReaclibReaction::getRateCoefficients() const {
        std::vector<RateCoefficientSet> rates;
        for (const auto& rate : m_rates) {
            rates.push_back(rate);
        }
        return rates;
    }

    ReactionSet::ReactionSet(
        std::vector<std::unique_ptr<Reaction>>&& reactions
    ) :
    m_reactions(std::move(reactions)) {
        if (m_reactions.empty()) {
            return; // Case where the reactions will be added later.
        }
        m_reactionNameMap.reserve(m_reactions.size());
        m_reactionHashes.reserve(m_reactions.size());
        size_t i = 0;
        for (const auto& reaction : m_reactions) {
            m_id += reaction->id();
            m_reactionNameMap.emplace(std::string(reaction->id()), i);
            m_reactionHashes.insert(reaction->hash(0));
            i++;
        }
    }

    ReactionSet::ReactionSet(const std::vector<Reaction *> &reactions) {
        m_reactions.reserve(reactions.size());
        m_reactionNameMap.reserve(reactions.size());
        m_reactionHashes.reserve(reactions.size());
        size_t i = 0;
        for (const auto& reaction : reactions) {
            m_reactions.push_back(reaction->clone());
            m_id += reaction->id();
            m_reactionNameMap.emplace(std::string(reaction->id()), i);
            m_reactionHashes.insert(reaction->hash(0));
            i++;
        }
    }

    ReactionSet::ReactionSet() = default;

    ReactionSet::ReactionSet(const ReactionSet &other) {
        m_reactions.reserve(other.m_reactions.size());
        m_reactionHashes.reserve(other.m_reactions.size());
        m_reactionNameMap.reserve(other.m_reactionNameMap.size());
        size_t i = 0;
        for (const auto& reaction: other.m_reactions) {
            m_reactions.push_back(reaction->clone());
            m_reactionNameMap.emplace(std::string(reaction->id()), i);
            m_reactionHashes.insert(reaction->hash(0));
            i++;
        }
    }

    ReactionSet& ReactionSet::operator=(const ReactionSet &other) {
        if (this != &other) {
            ReactionSet temp(other);
            std::swap(m_reactions, temp.m_reactions);
            std::swap(m_reactionNameMap, temp.m_reactionNameMap);
        }
        return *this;
    }

    void ReactionSet::add_reaction(const Reaction& reaction) {
        const std::size_t new_index = m_reactions.size();

        auto reaction_id = std::string(reaction.id());

        m_reactions.emplace_back(reaction.clone());
        m_id += reaction_id;

        m_reactionNameMap.emplace(std::move(reaction_id), new_index);

        m_reactionHashes.insert(reaction.hash(0));
    }

    void ReactionSet::add_reaction(std::unique_ptr<Reaction>&& reaction) {
        const std::size_t new_index = m_reactionNameMap.size();

        auto reaction_id = std::string(reaction->id());
        size_t reaction_hash = reaction->hash(0);

        m_reactions.emplace_back(std::move(reaction));

        m_id += reaction_id;

        m_reactionNameMap.emplace(std::move(reaction_id), new_index);

        m_reactionHashes.insert(reaction_hash);
    }

    void ReactionSet::extend(const ReactionSet &other) {
        for (const auto& reaction : other.m_reactions) {
            add_reaction(*reaction);
        }
    }

    std::optional<std::unique_ptr<Reaction>> ReactionSet::get(const std::string_view &id) const {
        if (!contains(id)) {
            return std::nullopt;
        }
        return std::make_optional(m_reactions[m_reactionNameMap.at(std::string(id))]->clone());
    }

    std::unique_ptr<Reaction> ReactionSet::get(size_t index) const {
        return m_reactions.at(index)->clone();
    }

    void ReactionSet::remove_reaction(const Reaction& reaction) {
        const size_t rh = reaction.hash(0);
        if (!m_reactionHashes.contains(rh)) {
            return;
        }

        std::erase_if(m_reactions, [&rh](const auto& r_ptr) {
            return r_ptr->hash(0) == rh;
        });

        m_reactionNameMap.erase(std::string(reaction.id()));

        m_id.clear();
        for (const auto& r_ptr : m_reactions) {
            m_id += r_ptr->id();
        }

        m_reactionHashes.erase(rh);
    }

    bool ReactionSet::contains(const std::string_view& id) const {
        return m_reactionNameMap.contains(std::string(id));
    }

    bool ReactionSet::contains(const Reaction& reaction) const {
        const size_t rh = reaction.hash(0);
        return m_reactionHashes.contains(rh);
    }

    void ReactionSet::clear() {
        m_reactions.clear();
        m_reactionNameMap.clear();
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
            m_logger -> flush_log();
            throw std::out_of_range("Index" + std::to_string(index) + " out of range for ReactionSet of size " + std::to_string(m_reactions.size()) + ".");
        }
        return *m_reactions[index];
    }

    const Reaction& ReactionSet::operator[](const std::string_view& id) const {
        if (!contains(id)) {
            m_logger -> flush_log();
            throw std::out_of_range("Reaction " + std::string(id) + " does not exist in ReactionSet.");
        }
        const size_t idx = m_reactionNameMap.at(std::string(id));
        return *m_reactions[idx];
    }

    bool ReactionSet::operator==(const ReactionSet& other) const {
        if (size() != other.size()) {
            return false;
        }
        return hash(0) == other.hash(0);
    }

    bool ReactionSet::operator!=(const ReactionSet& other) const {
        return !(*this == other);
    }

    uint64_t ReactionSet::hash(const uint64_t seed) const {
        if (m_reactions.empty()) {
            return XXHash64::hash(nullptr, 0, seed);
        }
        std::vector<uint64_t> individualReactionHashes;
        individualReactionHashes.reserve(m_reactions.size());
        for (const auto& reaction : m_reactions) {
            individualReactionHashes.push_back(reaction->hash(seed));
        }

        std::ranges::sort(individualReactionHashes);

        const auto data = static_cast<const void*>(individualReactionHashes.data());
        const size_t sizeInBytes = individualReactionHashes.size() * sizeof(uint64_t);
        return XXHash64::hash(data, sizeInBytes, seed);
    }

    std::unordered_set<Species> ReactionSet::getReactionSetSpecies() const {
        std::unordered_set<Species> species;
        for (const auto& reaction : m_reactions) {
            const auto reactionSpecies = reaction->all_species();
            species.insert(reactionSpecies.begin(), reactionSpecies.end());
        }
        return species;
    }

    ReactionSet packReactionSet(const ReactionSet& reactionSet) {
        std::unordered_map<std::string, std::vector<std::unique_ptr<ReaclibReaction>>> groupedReaclibReactions;
        ReactionSet finalReactionSet;

        for (const auto& reaction_ptr : reactionSet) {
            switch (reaction_ptr->type()) {
                case ReactionType::REACLIB_WEAK:
                    [[fallthrough]];
                case ReactionType::REACLIB: {
                    const auto& reaclib_cast_reaction = static_cast<const ReaclibReaction&>(*reaction_ptr); // NOLINT(*-pro-type-static-cast-downcast)
                    std::string rid = std::format("{}{}", reaclib_cast_reaction.peName(), (reaclib_cast_reaction.is_reverse() ? "_r" : ""));
                    groupedReaclibReactions[rid].push_back(std::make_unique<ReaclibReaction>(reaclib_cast_reaction));
                    break;
                }
                case ReactionType::LOGICAL_REACLIB_WEAK:
                    [[fallthrough]];
                case ReactionType::LOGICAL_REACLIB: {
                    // It doesn't make sense to pack an already-packed reaction.
                    throw std::runtime_error("packReactionSet: Cannot pack a LogicalReaclibReaction.");
                }
                case ReactionType::WEAK: {
                    finalReactionSet.add_reaction(*reaction_ptr);
                    break;
                }
            }
        }

        // Now, process the grouped REACLIB reactions
        for (const auto &reactionsGroup: groupedReaclibReactions | std::views::values) {
            if (reactionsGroup.empty()) {
                continue;
            }
            if (reactionsGroup.size() == 1) {
                finalReactionSet.add_reaction(reactionsGroup.front()->clone());
            }
            else {
                const auto logicalReaction = std::make_unique<LogicalReaclibReaction>(reactionsGroup, reactionsGroup.front()->is_reverse());
                finalReactionSet.add_reaction(logicalReaction->clone());
            }
        }

        return finalReactionSet;
    }

    bool reaction_is_weak(const reaction::Reaction& reaction) {
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
