#pragma once

#include <string_view>

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include <unordered_map>
#include <vector>
#include <unordered_set>


#include "cppad/cppad.hpp"
#include "xxhash64.h"

/**
 * @file reaction.h
 * @brief Defines classes for representing and managing nuclear reactions.
 *
 * This file contains the core data structures for handling nuclear reactions,
 * including individual reactions from specific sources (`Reaction`), collections
 * of reactions (`ReactionSet`), and logical reactions that aggregate rates from
* multiple sources (`LogicalReaction`, `LogicalReactionSet`).
 */
namespace gridfire::reaction {
    /**
     * @struct RateCoefficientSet
     * @brief Holds the seven coefficients for the REACLIB rate equation.
     *
     * This struct stores the parameters (a0-a6) used to calculate reaction rates
     * as a function of temperature.
     */
    struct RateCoefficientSet {
        double a0; ///< Coefficient a0
        double a1; ///< Coefficient a1
        double a2; ///< Coefficient a2
        double a3; ///< Coefficient a3
        double a4; ///< Coefficient a4
        double a5; ///< Coefficient a5
        double a6; ///< Coefficient a6

        /**
         * @brief Overloads the stream insertion operator for easy printing.
         * @param os The output stream.
         * @param r The RateCoefficientSet to print.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const RateCoefficientSet& r) {
            os << "[" << r.a0 << ", " << r.a1 << ", " << r.a2 << ", "
               << r.a3 << ", " << r.a4 << ", " << r.a5 << ", " << r.a6 << "]";
            return os;
        }
    };

    /**
     * @class Reaction
     * @brief Represents a single nuclear reaction from a specific data source.
     *
     * This class encapsulates all properties of a single nuclear reaction as defined
     * in formats like REACLIB, including reactants, products, Q-value, and rate
     * coefficients from a particular evaluation (source).
     *
     * Example:
     * @code
     * // Assuming species and rate coefficients are defined
     * Reaction p_gamma_d(
     *     "H_1_H_1_to_H_2", "p(p,g)d", 1, {H_1, H_1}, {H_2}, 5.493, "st08", rate_coeffs
     * );
     * double rate = p_gamma_d.calculate_rate(0.1); // T9 = 0.1
     * @endcode
     */
    class Reaction {
    public:
        /**
         * @brief Virtual destructor.
         */
        virtual ~Reaction() = default;

        /**
         * @brief Constructs a Reaction object.
         * @param id A unique identifier for the reaction.
         * @param peName The name in (projectile, ejectile) notation (e.g., "p(p,g)d").
         * @param chapter The REACLIB chapter number, defining reaction structure.
         * @param reactants A vector of reactant species.
         * @param products A vector of product species.
         * @param qValue The Q-value of the reaction in MeV.
         * @param label The source label for the rate data (e.g., "wc12", "st08").
         * @param sets The set of rate coefficients.
         * @param reverse True if this is a reverse reaction rate.
         */
        Reaction(
            const std::string_view id,
            const std::string_view peName,
            const int chapter,
            const std::vector<fourdst::atomic::Species> &reactants,
            const std::vector<fourdst::atomic::Species> &products,
            const double qValue,
            const std::string_view label,
            const RateCoefficientSet &sets,
            const bool reverse = false);

        /**
         * @brief Calculates the reaction rate for a given temperature.
         * @param T9 The temperature in units of 10^9 K.
         * @return The calculated reaction rate.
         */
        [[nodiscard]] virtual double calculate_rate(const double T9) const;

        /**
         * @brief Calculates the reaction rate for a given temperature using CppAD types.
         * @param T9 The temperature in units of 10^9 K, as a CppAD::AD<double>.
         * @return The calculated reaction rate, as a CppAD::AD<double>.
         */
        [[nodiscard]] virtual CppAD::AD<double> calculate_rate(const CppAD::AD<double> T9) const;

        [[nodiscard]] virtual double calculate_forward_rate_log_derivative(const double T9) const;

        /**
         * @brief Gets the reaction name in (projectile, ejectile) notation.
         * @return The reaction name (e.g., "p(p,g)d").
         */
        [[nodiscard]] virtual std::string_view peName() const { return m_peName; }

        /**
         * @brief Gets the REACLIB chapter number.
         * @return The chapter number.
         */
        [[nodiscard]] int chapter() const { return m_chapter; }

        /**
         * @brief Gets the source label for the rate data.
         * @return The source label (e.g., "wc12w", "st08").
         */
        [[nodiscard]] std::string_view sourceLabel() const { return m_sourceLabel; }

        /**
         * @brief Gets the set of rate coefficients.
         * @return A const reference to the RateCoefficientSet.
         */
        [[nodiscard]] const RateCoefficientSet& rateCoefficients() const { return m_rateCoefficients; }

        /**
         * @brief Checks if the reaction involves a given species as a reactant or product.
         * @param species The species to check for.
         * @return True if the species is involved, false otherwise.
         */
        [[nodiscard]] bool contains(const fourdst::atomic::Species& species) const;

        /**
         * @brief Checks if the reaction involves a given species as a reactant.
         * @param species The species to check for.
         * @return True if the species is a reactant, false otherwise.
         */
        [[nodiscard]] bool contains_reactant(const fourdst::atomic::Species& species) const;

        /**
         * @brief Checks if the reaction involves a given species as a product.
         * @param species The species to check for.
         * @return True if the species is a product, false otherwise.
         */
        [[nodiscard]] bool contains_product(const fourdst::atomic::Species& species) const;

        /**
         * @brief Gets a set of all unique species involved in the reaction.
         * @return An unordered_set of all reactant and product species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> all_species() const;

        /**
         * @brief Gets a set of all unique reactant species.
         * @return An unordered_set of reactant species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> reactant_species() const;

        /**
         * @brief Gets a set of all unique product species.
         * @return An unordered_set of product species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> product_species() const;

        /**
         * @brief Gets the number of unique species involved in the reaction.
         * @return The count of unique species.
         */
        [[nodiscard]] size_t num_species() const;

        /**
         * @brief Calculates the stoichiometric coefficient for a given species.
         * @param species The species for which to find the coefficient.
         * @return The stoichiometric coefficient (negative for reactants, positive for products).
         */
        [[nodiscard]] int stoichiometry(const fourdst::atomic::Species& species) const;

        /**
         * @brief Gets a map of all species to their stoichiometric coefficients.
         * @return An unordered_map from species to their integer coefficients.
         */
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, int> stoichiometry() const;

        /**
         * @brief Gets the unique identifier of the reaction.
         * @return The reaction ID.
         */
        [[nodiscard]] std::string_view id() const { return m_id; }

        /**
         * @brief Gets the Q-value of the reaction.
         * @return The Q-value in whatever units the reaction was defined in (usually MeV).
         */
        [[nodiscard]] double qValue() const { return m_qValue; }

        /**
         * @brief Gets the vector of reactant species.
         * @return A const reference to the vector of reactants.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& reactants() const { return m_reactants; }

        /**
         * @brief Gets the vector of product species.
         * @return A const reference to the vector of products.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& products() const { return m_products; }

        /**
         * @brief Checks if this is a reverse reaction rate.
         * @return True if it is a reverse rate, false otherwise.
         */
        [[nodiscard]] bool is_reverse() const { return m_reverse; }

        /**
         * @brief Calculates the excess energy from the mass difference of reactants and products.
         * @return The excess energy in MeV.
         */
        [[nodiscard]] double excess_energy() const;

        /**
         * @brief Compares this reaction with another for equality based on their IDs.
         * @param other The other Reaction to compare with.
         * @return True if the reaction IDs are the same.
         */
        bool operator==(const Reaction& other) const { return m_id == other.m_id; }

        /**
         * @brief Compares this reaction with another for inequality.
         * @param other The other Reaction to compare with.
         * @return True if the reactions are not equal.
         */
        bool operator!=(const Reaction& other) const { return !(*this == other); }

        /**
         * @brief Computes a hash for the reaction based on its ID.
         * @param seed The seed for the hash function.
         * @return A 64-bit hash value.
         * @details Uses the XXHash64 algorithm on the reaction's ID string.
         */
        [[nodiscard]] uint64_t hash(uint64_t seed = 0) const;

        friend std::ostream& operator<<(std::ostream& os, const Reaction& r) {
            return os << "(Reaction:" << r.m_id << ")";
        }

    protected:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::string m_id; ///< Unique identifier for the reaction (e.g., "h1+h1=>h2+e+nu").
        std::string m_peName; ///< Name of the reaction in (projectile, ejectile) notation (e.g. "p(p,g)d").
        int m_chapter;    ///< Chapter number from the REACLIB database, defining the reaction structure.
        double m_qValue = 0.0; ///< Q-value of the reaction in MeV.
        std::vector<fourdst::atomic::Species> m_reactants; ///< Reactants of the reaction.
        std::vector<fourdst::atomic::Species> m_products; ///< Products of the reaction.
        std::string m_sourceLabel; ///< Source label for the rate data (e.g., "wc12w", "st08").
        RateCoefficientSet m_rateCoefficients; ///< The seven rate coefficients.
        bool m_reverse = false; ///< Flag indicating if this is a reverse reaction rate.
    private:
        /**
         * @brief Template implementation for calculating the reaction rate.
         * @tparam T The numeric type (double or CppAD::AD<double>).
         * @param T9 The temperature in units of 10^9 K.
         * @return The calculated reaction rate.
         * @details The rate is calculated using the standard REACLIB formula:
         * `rate = exp(a0 + a1/T9 + a2/T9^(1/3) + a3*T9^(1/3) + a4*T9 + a5*T9^(5/3) + a6*ln(T9))`
         */
        template <typename T>
        [[nodiscard]] T calculate_rate(const T T9) const {
            const T T913 = CppAD::pow(T9, 1.0/3.0);
            const T T953 = CppAD::pow(T9, 5.0/3.0);
            const T logT9 = CppAD::log(T9);
            const T exponent = m_rateCoefficients.a0 +
                   m_rateCoefficients.a1 / T9 +
                   m_rateCoefficients.a2 / T913 +
                   m_rateCoefficients.a3 * T913 +
                   m_rateCoefficients.a4 * T9 +
                   m_rateCoefficients.a5 * T953 +
                   m_rateCoefficients.a6 * logT9;
            return CppAD::exp(exponent);
        }
    };



    /**
     * @class LogicalReaction
     * @brief Represents a "logical" reaction that aggregates rates from multiple sources.
     *
     * A LogicalReaction shares the same reactants and products but combines rates
     * from different evaluations (e.g., "wc12" and "st08" for the same physical
     * reaction). The total rate is the sum of the individual rates.
     * It inherits from Reaction, using the properties of the first provided reaction
     * as its base properties (reactants, products, Q-value, etc.).
     */
    class LogicalReaction final : public Reaction {
    public:
        /**
         * @brief Constructs a LogicalReaction from a vector of `Reaction` objects.
         * @param reactions A vector of reactions that represent the same logical process.
         * @throws std::runtime_error if the provided reactions have inconsistent Q-values.
         */
        explicit LogicalReaction(const std::vector<Reaction> &reactions);

        /**
         * @brief Adds another `Reaction` source to this logical reaction.
         * @param reaction The reaction to add.
         * @throws std::runtime_error if the reaction has a different `peName`, a duplicate
         *         source label, or an inconsistent Q-value.
         */
        void add_reaction(const Reaction& reaction);

        /**
         * @brief Gets the number of source rates contributing to this logical reaction.
         * @return The number of aggregated rates.
         */
        [[nodiscard]] size_t size() const { return m_rates.size(); }

        /**
         * @brief Gets the list of source labels for the aggregated rates.
         * @return A vector of source label strings.
         */
        [[nodiscard]] std::vector<std::string> sources() const { return m_sources; }

        /**
         * @brief Calculates the total reaction rate by summing all source rates.
         * @param T9 The temperature in units of 10^9 K.
         * @return The total calculated reaction rate.
         */
        [[nodiscard]] double calculate_rate(const double T9) const override;

        [[nodiscard]] virtual double calculate_forward_rate_log_derivative(const double T9) const override;

        /**
         * @brief Calculates the total reaction rate using CppAD types.
         * @param T9 The temperature in units of 10^9 K, as a CppAD::AD<double>.
         * @return The total calculated reaction rate, as a CppAD::AD<double>.
         */
        [[nodiscard]] CppAD::AD<double> calculate_rate(const CppAD::AD<double> T9) const override;

        /** @name Iterators
         *  Provides iterators to loop over the rate coefficient sets.
         */
        ///@{
        auto begin() { return m_rates.begin(); }
        [[nodiscard]] auto begin() const { return m_rates.cbegin(); }
        auto end() { return m_rates.end(); }
        [[nodiscard]] auto end() const { return m_rates.cend(); }
        ///@}
        ///

        friend std::ostream& operator<<(std::ostream& os, const LogicalReaction& r) {
            os << "(LogicalReaction: " << r.id() << ", reverse: " << r.is_reverse() << ")";
            return os;
        }

    private:
        std::vector<std::string> m_sources; ///< List of source labels.
        std::vector<RateCoefficientSet> m_rates; ///< List of rate coefficient sets from each source.

    private:
        /**
         * @brief Template implementation for calculating the total reaction rate.
         * @tparam T The numeric type (double or CppAD::AD<double>).
         * @param T9 The temperature in units of 10^9 K.
         * @return The total calculated reaction rate.
         * @details This method iterates through all stored `RateCoefficientSet`s,
         * calculates the rate for each, and returns their sum.
         */
        template <typename T>
        [[nodiscard]] T calculate_rate(const T T9) const {
            T sum = static_cast<T>(0.0);
            const T T913 = CppAD::pow(T9, 1.0/3.0);
            const T T953 = CppAD::pow(T9, 5.0/3.0);
            const T logT9 = CppAD::log(T9);
            // ReSharper disable once CppUseStructuredBinding
            for (const auto& rate : m_rates) {
                const T exponent = rate.a0 +
                       rate.a1 / T9 +
                       rate.a2 / T913 +
                       rate.a3 * T913 +
                       rate.a4 * T9 +
                       rate.a5 * T953 +
                       rate.a6 * logT9;
                sum += CppAD::exp(exponent);
            }
            return sum;
        }
    };

    template <typename ReactionT>
    class TemplatedReactionSet final {
    public:
        /**
         * @brief Constructs a ReactionSet from a vector of reactions.
         * @param reactions The initial vector of Reaction objects.
         */
        explicit TemplatedReactionSet(std::vector<ReactionT> reactions);

        TemplatedReactionSet();

        /**
         * @brief Copy constructor.
         * @param other The ReactionSet to copy.
         */
        TemplatedReactionSet(const TemplatedReactionSet<ReactionT>& other);

        /**
         * @brief Copy assignment operator.
         * @param other The ReactionSet to assign from.
         * @return A reference to this ReactionSet.
         */
        TemplatedReactionSet<ReactionT>& operator=(const TemplatedReactionSet<ReactionT>& other);

        /**
         * @brief Adds a reaction to the set.
         * @param reaction The Reaction to add.
         */
        void add_reaction(ReactionT reaction);

        /**
         * @brief Removes a reaction from the set.
         * @param reaction The Reaction to remove.
         */
        void remove_reaction(const ReactionT& reaction);

        /**
         * @brief Checks if the set contains a reaction with the given ID.
         * @param id The ID of the reaction to find.
         * @return True if the reaction is in the set, false otherwise.
         */
        [[nodiscard]] bool contains(const std::string_view& id) const;

        /**
         * @brief Checks if the set contains the given reaction.
         * @param reaction The Reaction to find.
         * @return True if the reaction is in the set, false otherwise.
         */
        [[nodiscard]] bool contains(const Reaction& reaction) const;

        /**
         * @brief Gets the number of reactions in the set.
         * @return The size of the set.
         */
        [[nodiscard]] size_t size() const { return m_reactions.size(); }

        /**
         * @brief Removes all reactions from the set.
         */
        void clear();

        /**
         * @brief Checks if any reaction in the set involves the given species.
         * @param species The species to check for.
         * @return True if the species is involved in any reaction.
         */
        [[nodiscard]] bool contains_species(const fourdst::atomic::Species& species) const;

        /**
         * @brief Checks if any reaction in the set contains the given species as a reactant.
         * @param species The species to check for.
         * @return True if the species is a reactant in any reaction.
         */
        [[nodiscard]] bool contains_reactant(const fourdst::atomic::Species& species) const;

        /**
         * @brief Checks if any reaction in the set contains the given species as a product.
         * @param species The species to check for.
         * @return True if the species is a product in any reaction.
         */
        [[nodiscard]] bool contains_product(const fourdst::atomic::Species& species) const;

        /**
         * @brief Accesses a reaction by its index.
         * @param index The index of the reaction to access.
         * @return A const reference to the Reaction.
         * @throws std::out_of_range if the index is out of bounds.
         */
        [[nodiscard]] const ReactionT& operator[](size_t index) const;

        /**
         * @brief Accesses a reaction by its ID.
         * @param id The ID of the reaction to access.
         * @return A const reference to the Reaction.
         * @throws std::out_of_range if no reaction with the given ID exists.
         */
        [[nodiscard]] const ReactionT& operator[](const std::string_view& id) const;

        /**
         * @brief Compares this set with another for equality.
         * @param other The other ReactionSet to compare with.
         * @return True if the sets are equal (same size and hash).
         */
        bool operator==(const TemplatedReactionSet& other) const;

        /**
         * @brief Compares this set with another for inequality.
         * @param other The other ReactionSet to compare with.
         * @return True if the sets are not equal.
         */
        bool operator!=(const TemplatedReactionSet& other) const;

        /**
         * @brief Computes a hash for the entire set.
         * @param seed The seed for the hash function.
         * @return A 64-bit hash value.
         * @details The algorithm computes the hash of each individual reaction,
         * sorts the hashes, and then computes a final hash over the sorted list
         * of hashes. This ensures the hash is order-independent.
         */
        [[nodiscard]] uint64_t hash(uint64_t seed = 0) const;

        /** @name Iterators
         *  Provides iterators to loop over the reactions in the set.
         */
        ///@{
        auto begin() { return m_reactions.begin(); }
        [[nodiscard]] auto begin() const { return m_reactions.cbegin(); }
        auto end() { return m_reactions.end(); }
        [[nodiscard]] auto end() const { return m_reactions.cend(); }
        ///@}
        ///
        friend std::ostream& operator<<(std::ostream& os, const TemplatedReactionSet<ReactionT>& r) {
            os << "(ReactionSet: [";
            size_t counter = 0;
            for (const auto& reaction : r.m_reactions) {
                os << reaction;
                if (counter < r.m_reactions.size() - 2) {
                    os << ", ";
                } else if (counter == r.m_reactions.size() - 2) {
                    os << " and ";
                }
                ++counter;
            }
            os << "])";
            return os;
        }

        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> getReactionSetSpecies() const;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<ReactionT> m_reactions;
        std::string m_id;
        std::unordered_map<std::string, ReactionT> m_reactionNameMap; ///< Maps reaction IDs to Reaction objects for quick lookup.

    };

    using ReactionSet = TemplatedReactionSet<Reaction>; ///< A set of reactions, typically from a single source like REACLIB.
    using LogicalReactionSet = TemplatedReactionSet<LogicalReaction>; ///< A set of logical reactions.

    LogicalReactionSet packReactionSetToLogicalReactionSet(const ReactionSet& reactionSet);

    template <typename ReactionT>
    TemplatedReactionSet<ReactionT>::TemplatedReactionSet(
        std::vector<ReactionT> reactions
    ) :
    m_reactions(std::move(reactions)) {
        if (m_reactions.empty()) {
            return; // Case where the reactions will be added later.
        }
        m_reactionNameMap.reserve(reactions.size());
        for (const auto& reaction : m_reactions) {
            m_id += reaction.id();
            m_reactionNameMap.emplace(reaction.id(), reaction);
        }
    }

    template<typename ReactionT>
    TemplatedReactionSet<ReactionT>::TemplatedReactionSet() {}

    template <typename ReactionT>
    TemplatedReactionSet<ReactionT>::TemplatedReactionSet(const TemplatedReactionSet<ReactionT> &other) {
        m_reactions.reserve(other.m_reactions.size());
        for (const auto& reaction_ptr: other.m_reactions) {
            m_reactions.push_back(reaction_ptr);
        }

        m_reactionNameMap.reserve(other.m_reactionNameMap.size());
        for (const auto& reaction_ptr : m_reactions) {
            m_reactionNameMap.emplace(reaction_ptr.id(), reaction_ptr);
        }
    }

    template <typename ReactionT>
    TemplatedReactionSet<ReactionT>& TemplatedReactionSet<ReactionT>::operator=(const TemplatedReactionSet<ReactionT> &other) {
        if (this != &other) {
            TemplatedReactionSet temp(other);
            std::swap(m_reactions, temp.m_reactions);
            std::swap(m_reactionNameMap, temp.m_reactionNameMap);
        }
        return *this;
    }

    template <typename ReactionT>
    void TemplatedReactionSet<ReactionT>::add_reaction(ReactionT reaction) {
        m_reactions.emplace_back(reaction);
        m_id += m_reactions.back().id();
        m_reactionNameMap.emplace(m_reactions.back().id(), m_reactions.back());
    }

    template <typename ReactionT>
    void TemplatedReactionSet<ReactionT>::remove_reaction(const ReactionT& reaction) {
        if (!m_reactionNameMap.contains(std::string(reaction.id()))) {
            return;
        }

        m_reactionNameMap.erase(std::string(reaction.id()));

        std::erase_if(m_reactions, [&reaction](const Reaction& r) {
            return r == reaction;
        });
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::contains(const std::string_view& id) const {
        for (const auto& reaction : m_reactions) {
            if (reaction.id() == id) {
                return true;
            }
        }
        return false;
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::contains(const Reaction& reaction) const {
        for (const auto& r : m_reactions) {
            if (r == reaction) {
                return true;
            }
        }
        return false;
    }

    template <typename ReactionT>
    void TemplatedReactionSet<ReactionT>::clear() {
        m_reactions.clear();
        m_reactionNameMap.clear();
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::contains_species(const fourdst::atomic::Species& species) const {
        for (const auto& reaction : m_reactions) {
            if (reaction.contains(species)) {
                return true;
            }
        }
        return false;
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::contains_reactant(const fourdst::atomic::Species& species) const {
        for (const auto& r : m_reactions) {
            if (r.contains_reactant(species)) {
                return true;
            }
        }
        return false;
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::contains_product(const fourdst::atomic::Species& species) const {
        for (const auto& r : m_reactions) {
            if (r.contains_product(species)) {
                return true;
            }
        }
        return false;
    }

    template <typename ReactionT>
    const ReactionT& TemplatedReactionSet<ReactionT>::operator[](const size_t index) const {
        if (index >= m_reactions.size()) {
            m_logger -> flush_log();
            throw std::out_of_range("Index" + std::to_string(index) + " out of range for ReactionSet of size " + std::to_string(m_reactions.size()) + ".");
        }
        return m_reactions[index];
    }

    template <typename ReactionT>
    const ReactionT& TemplatedReactionSet<ReactionT>::operator[](const std::string_view& id) const {
        if (auto it = m_reactionNameMap.find(std::string(id)); it != m_reactionNameMap.end()) {
            return it->second;
        }
        m_logger -> flush_log();
        throw std::out_of_range("Species " + std::string(id) + " does not exist in ReactionSet.");
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::operator==(const TemplatedReactionSet<ReactionT>& other) const {
        if (size() != other.size()) {
            return false;
        }
        return hash() == other.hash();
    }

    template <typename ReactionT>
    bool TemplatedReactionSet<ReactionT>::operator!=(const TemplatedReactionSet<ReactionT>& other) const {
        return !(*this == other);
    }

    template <typename ReactionT>
    uint64_t TemplatedReactionSet<ReactionT>::hash(uint64_t seed) const {
        if (m_reactions.empty()) {
            return XXHash64::hash(nullptr, 0, seed);
        }
        std::vector<uint64_t> individualReactionHashes;
        individualReactionHashes.reserve(m_reactions.size());
        for (const auto& reaction : m_reactions) {
            individualReactionHashes.push_back(reaction.hash(seed));
        }

        std::ranges::sort(individualReactionHashes);

        const auto data = static_cast<const void*>(individualReactionHashes.data());
        const size_t sizeInBytes = individualReactionHashes.size() * sizeof(uint64_t);
        return XXHash64::hash(data, sizeInBytes, seed);
    }

    template<typename ReactionT>
    std::unordered_set<fourdst::atomic::Species> TemplatedReactionSet<ReactionT>::getReactionSetSpecies() const {
        std::unordered_set<fourdst::atomic::Species> species;
        for (const auto& reaction : m_reactions) {
            const auto reactionSpecies = reaction.all_species();
            species.insert(reactionSpecies.begin(), reactionSpecies.end());
        }
        return species;
    }
}

