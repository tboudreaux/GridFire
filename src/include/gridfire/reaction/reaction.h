#pragma once

#include <ranges>
#include <string_view>

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <optional>

#include "cppad/cppad.hpp"
#include "fourdst/composition/composition.h"

/**
 * @file reaction.h
 * @brief Defines classes for representing and managing nuclear reactions.
 *
 * This file contains the core data structures for handling nuclear reactions,
 * including individual reactions from specific sources (`Reaction`), collections
 * of reactions (`ReactionSet`), and logical reactions that aggregate rates from
* multiple sources (`LogicalReaclibReaction`, `LogicalReactionSet`).
 */
namespace gridfire::reaction {
    enum class ReactionType {
        WEAK,
        REACLIB,
        REACLIB_WEAK,
        LOGICAL_REACLIB,
        LOGICAL_REACLIB_WEAK,
    };

    static std::unordered_map<ReactionType, std::string> ReactionTypeNames = {
        {ReactionType::WEAK, "weak"},
        {ReactionType::REACLIB, "reaclib"},
        {ReactionType::REACLIB_WEAK, "reaclib_weak"},
        {ReactionType::LOGICAL_REACLIB, "logical_reaclib"},
        {ReactionType::LOGICAL_REACLIB_WEAK, "logical_reaclib_weak"},

    };

    static std::unordered_map<ReactionType, std::string> ReactionPhysicalTypeNames = {
        {ReactionType::WEAK, "Weak"},
        {ReactionType::REACLIB, "Strong"},
        {ReactionType::LOGICAL_REACLIB, "Strong"},
        {ReactionType::REACLIB_WEAK, "Weak"},
        {ReactionType::LOGICAL_REACLIB_WEAK, "Weak"},
    };

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
         * @brief Virtual destructor for correct polymorphic cleanup.
         */
        virtual ~Reaction() = default;

        /**
         * @brief Compute the temperature- and composition-dependent reaction rate.
         *
         * This is the primary interface used by the network to obtain the rate of a single
         * reaction at the given thermodynamic state and composition. The exact units and
         * normalization are defined by the concrete implementation (e.g., REACLIB typically
         * provides NA<sigma v> with units depending on the reaction order). Implementations may
         * use density/electron properties for weak processes or screening, and the composition
         * vector for multi-body reactions.
         *
         * @param T9 Temperature in GK (10^9 K).
         * @param rho Mass density (g cm^-3). May be unused for some reaction types.
         * @param Ye Electron fraction. May be unused depending on the reaction type.
         * @param mue Electron chemical potential. May be unused depending on the reaction type.
         * @param Y Composition vector (molar abundances or number fractions) indexed consistently
         *          with index_to_species_map.
         * @param index_to_species_map Mapping from state-vector index to Species, used to interpret Y.
         * @return The reaction rate for the forward direction, with units/normalization defined by the
         *         specific model (implementation must document its convention).
         */
        [[nodiscard]] virtual double calculate_rate(
            double T9,
            double rho,
            double Ye,
            double mue,
            const std::vector<double> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const = 0;

        /**
         * @brief AD-enabled reaction rate for algorithmic differentiation.
         *
         * This overload mirrors calculate_rate(double, ...) but operates on CppAD types to
         * enable derivative calculations w.r.t. its inputs.
         *
         * @param T9 Temperature in GK as CppAD::AD<double>.
         * @param rho Mass density as CppAD::AD<double>.
         * @param Ye Electron fraction as CppAD::AD<double>.
         * @param mue Electron chemical potential as CppAD::AD<double>.
         * @param Y Composition vector as CppAD::AD<double> values.
         * @param index_to_species_map Mapping from state-vector index to Species, used to interpret Y.
         * @return The reaction rate as a CppAD::AD<double> value.
         */
        [[nodiscard]] virtual CppAD::AD<double> calculate_rate(
            CppAD::AD<double> T9,
            CppAD::AD<double> rho,
            CppAD::AD<double> Ye,
            CppAD::AD<double> mue,
            const std::vector<CppAD::AD<double>>& Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const = 0;

        /**
         * @brief A stable, unique identifier for this reaction instance.
         * @return String view of the reaction ID (stable across runs and suitable for lookups).
         */
        [[nodiscard]] virtual std::string_view id() const = 0;

        /**
         * @brief Ordered list of reactant species.
         *
         * Multiplicity is represented by duplicates, e.g., (p, p) would list H1 twice.
         * @return Const reference to the vector of reactants.
         */
        [[nodiscard]] virtual const std::vector<fourdst::atomic::Species>& reactants() const = 0;

        /**
         * @brief Ordered list of product species.
         *
         * Multiplicity is represented by duplicates if applicable.
         * @return Const reference to the vector of products.
         */
        [[nodiscard]] virtual const std::vector<fourdst::atomic::Species>& products() const = 0;

        /**
         * @brief True if the species appears as a reactant or a product.
         * @param species Species to test.
         * @return Whether the species participates in the reaction (either side).
         */
        [[nodiscard]] virtual bool contains(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief True if the species appears among the reactants.
         * @param species Species to test.
         * @return Whether the species is a reactant.
         */
        [[nodiscard]] virtual bool contains_reactant(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief True if the species appears among the products.
         * @param species Species to test.
         * @return Whether the species is a product.
         */
        [[nodiscard]] virtual bool contains_product(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Whether this object represents a reverse (backward) rate.
         *
         * Implementations may pair forward/reverse rates for detailed balance. This flag indicates
         * that the parameterization corresponds to the reverse direction.
         * @return True for a reverse rate, false for a forward rate.
         */
        [[nodiscard]] virtual bool is_reverse() const = 0;

        /**
         * @brief Set of all unique species appearing in the reaction.
         * @return Unordered set of all reactants and products (no duplicates).
         */
        [[nodiscard]] virtual std::unordered_set<fourdst::atomic::Species> all_species() const = 0;

        /**
         * @brief Set of unique reactant species.
         * @return Unordered set of reactant species (no duplicates).
         */
        [[nodiscard]] virtual std::unordered_set<fourdst::atomic::Species> reactant_species() const = 0;

        /**
         * @brief Set of unique product species.
         * @return Unordered set of product species (no duplicates).
         */
        [[nodiscard]] virtual std::unordered_set<fourdst::atomic::Species> product_species() const = 0;

        [[nodiscard]] virtual size_t countReactantOccurrences(const fourdst::atomic::Species& species) const = 0;

        [[nodiscard]] virtual size_t countProductOccurrences(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Number of unique species involved in the reaction.
         * @return Count of distinct species across reactants and products.
         */
        [[nodiscard]] virtual size_t num_species() const = 0;

        /**
         * @brief Full stoichiometry map for this reaction.
         *
         * Coefficients are negative for reactants and positive for products; multiplicity is reflected
         * in the magnitude (e.g., 2H -> He gives H: -2, He: +1).
         * @return Map from Species to integer stoichiometric coefficient.
         */
        [[nodiscard]] virtual std::unordered_map<fourdst::atomic::Species, int> stoichiometry() const = 0;

        /**
         * @brief Stoichiometric coefficient for a particular species.
         * @param species Species for which to query the coefficient.
         * @return Negative for reactants, positive for products, zero if absent.
         */
        [[nodiscard]] virtual int stoichiometry(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Stable content-based hash for this reaction.
         *
         * Intended for use in caches, sets, and order-independent hashing of Reaction collections.
         * Implementations should produce the same value across processes for the same content and seed.
         * @param seed Seed value to initialize/mix into the hash.
         * @return 64-bit hash value.
         */
        [[nodiscard]] virtual uint64_t hash(uint64_t seed) const = 0;

        /**
         * @brief Q-value of the reaction (typically MeV), positive if exothermic.
         * @return Reaction Q-value used for energy accounting.
         */
        [[nodiscard]] virtual double qValue() const = 0;

        /**
         * @brief Convenience: energy generation rate from this reaction (double version).
         *
         * Default implementation multiplies the scalar rate by the reaction Q-value. Electron
         * quantities (Ye, mue) are ignored in this default, so override in derived classes if
         * needed. Sign convention follows qValue().
         *
         * @param T9 Temperature in GK (10^9 K).
         * @param rho Mass density (g cm^-3).
         * @param Ye Electron fraction (ignored by default implementation).
         * @param mue Electron chemical potential (ignored by default implementation).
         * @param Y Composition vector.
         * @param index_to_species_map Mapping from state-vector index to Species.
         * @return Energy generation rate, typically rate * qValue().
         */
        [[nodiscard]] virtual double calculate_energy_generation_rate(
            const double T9,
            const double rho,
            const double Ye,
            double mue,
            const std::vector<double>& Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const {
            return calculate_rate(T9, rho, 0, 0, Y, index_to_species_map) * qValue();
        }

        /**
         * @brief Convenience: AD-enabled energy generation rate (AD version).
         *
         * Default implementation multiplies the AD rate by the reaction Q-value. Electron
         * quantities (Ye, mue) are ignored in this default, so override if they contribute.
         *
         * @param T9 Temperature in GK as CppAD::AD<double>.
         * @param rho Mass density as CppAD::AD<double>.
         * @param Ye Electron fraction as CppAD::AD<double> (ignored by default).
         * @param mue Electron chemical potential as CppAD::AD<double> (ignored by default).
         * @param Y Composition vector as CppAD::AD<double> values.
         * @param index_to_species_map Mapping from state-vector index to Species.
         * @return Energy generation rate as CppAD::AD<double>.
         */
        [[nodiscard]] virtual CppAD::AD<double> calculate_energy_generation_rate(
            const CppAD::AD<double>& T9,
            const CppAD::AD<double>& rho,
            const CppAD::AD<double> &Ye,
            const CppAD::AD<double> &mue,
            const std::vector<CppAD::AD<double>>& Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const {
            return calculate_rate(T9, rho, {}, {}, Y, index_to_species_map) * qValue();
        }

        /**
         * @brief Logarithmic partial derivative of the rate with respect to temperature.
         *
         * Implementations return d(ln rate)/d(ln T9) or an equivalent measure (as documented by the
         * concrete class), evaluated at the provided state.
         *
         * @param T9 Temperature in GK (10^9 K).
         * @param rho Mass density (g cm^-3).
         * @param Ye Electron fraction.
         * @param mue Electron chemical potential.
         * @param comp Composition object providing composition in a convenient form.
         * @return The logarithmic temperature derivative of the rate.
         */
        [[nodiscard]] virtual double calculate_log_rate_partial_deriv_wrt_T9(
            double T9,
            double rho,
            double Ye,
            double mue,
            const fourdst::composition::Composition& comp
        ) const = 0;

        /**
         * @brief Category of this reaction (e.g., REACLIB, WEAK, LOGICAL_REACLIB).
         * @return Enumerated reaction type for runtime dispatch and filtering.
         */
        [[nodiscard]] virtual ReactionType type() const = 0;

        /**
         * @brief Polymorphic deep copy.
         * @return A std::unique_ptr owning a new Reaction equal to this one.
         */
        [[nodiscard]] virtual std::unique_ptr<Reaction> clone() const = 0;

        friend std::ostream& operator<<(std::ostream& os, const Reaction& r) {
            os << "Reaction(ID: " << r.id() << ")";
            return os;
        }

        [[nodiscard]] virtual std::optional<std::vector<RateCoefficientSet>> getRateCoefficients() const = 0;

    };

    // Simple heuristic to check if a reaclib reaction is a strong or weak reaction
    /*  A weak reaction is defined here as one where:
        - The number of reactants is equal to the number of products
        - There is only one reactant and one product
        - The mass number (A) of the reactant is equal to the mass number (A) of the product
    */
    bool reaction_is_weak(const Reaction& reaction);

    class ReaclibReaction : public Reaction {
    public:
        ~ReaclibReaction() override = default;

        /**
         * @brief Constructs a Reaction object.
         * @param id A unique identifier for the reaction.
         * @param peName The name in (projectile, ejectile) notation (e.g., "p(p,g)d").
         * @param chapter The REACLIB chapter number, defining reaction structure.
         * @param reactants A vector of reactant species.
         * @param products A vector of product species.
         * @param qValue The Q-value of the reaction in MeV.
         * @param label The sources label for the rate data (e.g., "wc12", "st08").
         * @param sets The set of rate coefficients.
         * @param reverse True if this is a reverse reaction rate.
         */
        ReaclibReaction(
            std::string_view id,
            std::string_view peName,
            int chapter,
            const std::vector<fourdst::atomic::Species> &reactants,
            const std::vector<fourdst::atomic::Species> &products,
            double qValue,
            std::string_view label,
            const RateCoefficientSet &sets,
            bool reverse = false);

        /**
         * @brief Calculates the reaction rate for a given temperature.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho Density [Not used in this implementation].
         * @param Ye
         * @param mue
         * @param Y
         * @param index_to_species_map
         * @return The calculated reaction rate.
         */
        [[nodiscard]] double calculate_rate(
            double T9,
            double rho,
            double Ye,
            double mue,
            const std::vector<double> &Y,
            const std::unordered_map<size_t,
            fourdst::atomic::Species>& index_to_species_map
        ) const override;

        /**
         * @brief Calculates the reaction rate for a given temperature using CppAD types.
         * @param T9 The temperature in units of 10^9 K, as a CppAD::AD<double>.
         * @param rho Density, as a CppAD::AD<double> [Not used in this implementation].
         * @param Ye
         * @param mue
         * @param Y Molar abundances of species, as a vector of CppAD::AD<double> [Not used in this implementation].
         * @param index_to_species_map
         * @return The calculated reaction rate, as a CppAD::AD<double>.
         */
        [[nodiscard]] CppAD::AD<double> calculate_rate(
            CppAD::AD<double> T9,
            CppAD::AD<double> rho,
            CppAD::AD<double> Ye,
            CppAD::AD<double> mue, const std::vector<CppAD::AD<double>>& Y, const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const override;

        [[nodiscard]] double calculate_log_rate_partial_deriv_wrt_T9(
            double T9,
            double rho,
            double Ye,
            double mue,
            const fourdst::composition::Composition& comp
        ) const override;

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

        [[nodiscard]] ReactionType type() const override {
            if (reaction::reaction_is_weak(*this)) {
                return ReactionType::REACLIB_WEAK;
            }
            return ReactionType::REACLIB;
        }

        /**
         * @brief Gets the set of rate coefficients.
         * @return A const reference to the RateCoefficientSet.
         */
        [[nodiscard]] const RateCoefficientSet& rateCoefficients() const { return m_rateCoefficients; }

        [[nodiscard]] std::optional<std::vector<RateCoefficientSet>> getRateCoefficients() const override;

        /**
         * @brief Checks if the reaction involves a given species as a reactant or product.
         * @param species The species to check for.
         * @return True if the species is involved, false otherwise.
         */
        [[nodiscard]] bool contains(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Checks if the reaction involves a given species as a reactant.
         * @param species The species to check for.
         * @return True if the species is a reactant, false otherwise.
         */
        [[nodiscard]] bool contains_reactant(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Checks if the reaction involves a given species as a product.
         * @param species The species to check for.
         * @return True if the species is a product, false otherwise.
         */
        [[nodiscard]] bool contains_product(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Gets a set of all unique species involved in the reaction.
         * @return An unordered_set of all reactant and product species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> all_species() const override;

        /**
         * @brief Gets a set of all unique reactant species.
         * @return An unordered_set of reactant species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> reactant_species() const override;

        /**
         * @brief Gets a set of all unique product species.
         * @return An unordered_set of product species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> product_species() const override;

        /**
         * @brief Gets the number of unique species involved in the reaction.
         * @return The count of unique species.
         */
        [[nodiscard]] size_t num_species() const override;

        /**
         * @brief Calculates the stoichiometric coefficient for a given species.
         * @param species The species for which to find the coefficient.
         * @return The stoichiometric coefficient (negative for reactants, positive for products).
         */
        [[nodiscard]] int stoichiometry(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Gets a map of all species to their stoichiometric coefficients.
         * @return An unordered_map from species to their integer coefficients.
         */
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, int> stoichiometry() const override;

        /**
         * @brief Gets the unique identifier of the reaction.
         * @return The reaction ID.
         */
        [[nodiscard]] std::string_view id() const override { return m_id; }

        /**
         * @brief Gets the Q-value of the reaction.
         * @return The Q-value in whatever units the reaction was defined in (usually MeV).
         */
        [[nodiscard]] double qValue() const override { return m_qValue; }

        /**
         * @brief Gets the vector of reactant species.
         * @return A const reference to the vector of reactants.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& reactants() const override {
            if (!m_reactantsVec) {
                m_reactantsVec.emplace(std::vector<fourdst::atomic::Species>());
                m_reactantsVec->reserve(m_reactants.size());
                for (const auto& [reactant, count] : m_reactants) {
                    for (size_t i = 0; i < count; ++i) {
                        m_reactantsVec->push_back(reactant);
                    }
                }
            }
            return m_reactantsVec.value();
        }

        /**
         * @brief Gets the vector of product species.
         * @return A const reference to the vector of products.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& products() const override {
            if (!m_productsVec) {
                m_productsVec.emplace(std::vector<fourdst::atomic::Species>());
                m_productsVec->reserve(m_products.size());
                for (const auto& [product, count] : m_products) {
                    for (size_t i = 0; i < count; ++i) {
                        m_productsVec->push_back(product);
                    }
                }
            }
            return m_productsVec.value();
        }

        /**
         * @brief Checks if this is a reverse reaction rate.
         * @return True if it is a reverse rate, false otherwise.
         */
        [[nodiscard]] bool is_reverse() const override { return m_reverse; }

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
        bool operator==(const ReaclibReaction& other) const { return m_id == other.m_id; }

        /**
         * @brief Compares this reaction with another for inequality.
         * @param other The other Reaction to compare with.
         * @return True if the reactions are not equal.
         */
        bool operator!=(const ReaclibReaction& other) const { return !(*this == other); }

        /**
         * @brief Computes a hash for the reaction based on its ID.
         * @param seed The seed for the hash function.
         * @return A 64-bit hash value.
         * @details Uses the XXHash64 algorithm on the reaction's ID string.
         */
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;

        [[nodiscard]] std::unique_ptr<Reaction> clone() const override;

        friend std::ostream& operator<<(std::ostream& os, const ReaclibReaction& r) {
            return os << "(ReaclibReaction:" << r.m_id << ")";
        }

        size_t countReactantOccurrences(const fourdst::atomic::Species& species) const override {
            if (!m_reactants.contains(species)) return 0;
            return m_reactants.at(species);
        }

        size_t countProductOccurrences(const fourdst::atomic::Species& species) const override {
            if (!m_products.contains(species)) return 0;
            return m_products.at(species);
        }

    protected:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::string m_id; ///< Unique identifier for the reaction (e.g., "h1+h1=>h2+e+nu").
        std::string m_peName; ///< Name of the reaction in (projectile, ejectile) notation (e.g. "p(p,g)d").
        int m_chapter;    ///< Chapter number from the REACLIB database, defining the reaction structure.
        double m_qValue = 0.0; ///< Q-value of the reaction in MeV.
        std::unordered_map<fourdst::atomic::Species, size_t> m_reactants; ///< Reactants of the reaction.
        std::unordered_map<fourdst::atomic::Species, size_t> m_products; ///< Products of the reaction.

        mutable std::optional<std::vector<fourdst::atomic::Species>> m_reactantsVec;
        mutable std::optional<std::vector<fourdst::atomic::Species>> m_productsVec;
        mutable std::optional<std::size_t> m_hashCache = std::nullopt; ///< Cached hash value for the reaction.

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

    class WeakReaclibReaction final : public ReaclibReaction {
    public:
        using ReaclibReaction::ReaclibReaction;

        [[nodiscard]] ReactionType type() const override { return ReactionType::REACLIB_WEAK; }

        [[nodiscard]] std::unique_ptr<Reaction> clone() const override {
            return std::make_unique<WeakReaclibReaction>(
                m_id,
                m_peName,
                m_chapter,
                reactants(),
                products(),
                m_qValue,
                m_sourceLabel,
                m_rateCoefficients,
                m_reverse
            );
        }
    };


    /**
     * @class LogicalReaclibReaction
     * @brief Represents a "logical" reaction that aggregates rates from multiple sources.
     *
     * A LogicalReaclibReaction shares the same reactants and products but combines rates
     * from different evaluations (e.g., "wc12" and "st08" for the same physical
     * reaction). The total rate is the sum of the individual rates.
     * It inherits from Reaction, using the properties of the first provided reaction
     * as its base properties (reactants, products, Q-value, etc.).
     */
    class LogicalReaclibReaction final : public ReaclibReaction {
    public:
        /**
         * @brief Constructs a LogicalReaction from a vector of `Reaction` objects. Implicitly assumes that the
         * logical reaction is for a forward (i.e. not reverse) reaction.
         * @param reactions A vector of reactions that represent the same logical process.
         * @throws std::runtime_error if the provided reactions have inconsistent Q-values.
         */
        explicit LogicalReaclibReaction(const std::vector<std::unique_ptr<ReaclibReaction>> &reactions);

        /**
         * @breif Constructs a LogicalReaction from a vector of `Reaction` objects and allows the user
         * to specify if the logical set is for a reverse reaction explicitly
         * @param reactions A vector of reactions that represent the same logical process
         * @param reverse A flag to control if this logical reaction is reverse or not
         * @returns std::runtime_error if the provided reactions have inconsistent Q-values.
         */
        explicit LogicalReaclibReaction(const std::vector<std::unique_ptr<ReaclibReaction>> &reactions, bool reverse);

        /**
         * @brief Adds another `Reaction` source to this logical reaction.
         * @param reaction The reaction to add.
         * @throws std::runtime_error if the reaction has a different `peName`, a duplicate
         *         source label, or an inconsistent Q-value.
         */
        void add_reaction(const ReaclibReaction& reaction);

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
         * @param rho
         * @param Ye
         * @param mue
         * @param Y
         * @param index_to_species_map
         * @return The total calculated reaction rate.
         */
        [[nodiscard]] double calculate_rate(
            double T9,
            double rho,
            double Ye,
            double mue, const std::vector<double> &Y, const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const override;

        [[nodiscard]] double calculate_log_rate_partial_deriv_wrt_T9(
            double T9,
            double rho,
            double Ye, double mue, const fourdst::composition::Composition& comp
        ) const override;

        [[nodiscard]] ReactionType type() const override {
            if (m_weak) {
                return ReactionType::LOGICAL_REACLIB_WEAK;
            }
            return ReactionType::LOGICAL_REACLIB;
        }

        [[nodiscard]] std::unique_ptr<Reaction> clone() const override;

        /**
         * @brief Calculates the total reaction rate using CppAD types.
         * @param T9 The temperature in units of 10^9 K, as a CppAD::AD<double>.
         * @param rho
         * @param Ye
         * @param mue
         * @param Y
         * @param index_to_species_map
         * @return The total calculated reaction rate, as a CppAD::AD<double>.
         */
        [[nodiscard]] CppAD::AD<double> calculate_rate(
            CppAD::AD<double> T9,
            CppAD::AD<double> rho,
            CppAD::AD<double> Ye,
            CppAD::AD<double> mue, const std::vector<CppAD::AD<double>>& Y, const std::unordered_map<size_t,fourdst::atomic::Species>& index_to_species_map
        ) const override;

        [[nodiscard]] std::optional<std::vector<RateCoefficientSet>> getRateCoefficients() const override;

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

        friend std::ostream& operator<<(std::ostream& os, const LogicalReaclibReaction& r) {
            os << "(LogicalReaclibReaction: " << r.id() << ", reverse: " << r.is_reverse() << ")";
            return os;
        }

    private:
        std::vector<std::string> m_sources; ///< List of source labels.
        std::vector<RateCoefficientSet> m_rates; ///< List of rate coefficient sets from each source.
        bool m_weak = false;

        mutable std::unordered_map<double, double> m_cached_rates;

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
            for (const auto& [rate, source] : std::views::zip(m_rates, m_sources)) {
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

    class ReactionSet final {
    public:
        /**
         * @brief Constructs a ReactionSet from a vector of reactions.
         * @param reactions The initial vector of Reaction objects.
         */
        explicit ReactionSet(std::vector<std::unique_ptr<Reaction>>&& reactions);

        explicit ReactionSet(const std::vector<Reaction*>& reactions);

        ReactionSet();

        /**
         * @brief Copy constructor.
         * @param other The ReactionSet to copy.
         */
        ReactionSet(const ReactionSet& other);

        /**
         * @brief Copy assignment operator.
         * @param other The ReactionSet to assign from.
         * @return A reference to this ReactionSet.
         */
        ReactionSet& operator=(const ReactionSet& other);

        /**
         * @brief Adds a reaction to the set.
         * @param reaction The Reaction to add.
         */
        void add_reaction(const Reaction& reaction);

        void add_reaction(std::unique_ptr<Reaction>&& reaction);

        void extend(const ReactionSet& other);

        [[nodiscard]] std::optional<std::unique_ptr<Reaction>> get(const std::string_view& id) const;

        [[nodiscard]] std::unique_ptr<Reaction> get(size_t index) const;

        /**
         * @brief Removes a reaction from the set.
         * @param reaction The Reaction to remove.
         */
        void remove_reaction(const Reaction& reaction);

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

        [[nodiscard]] bool empty() const {return m_reactions.empty(); }

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
        [[nodiscard]] const Reaction& operator[](size_t index) const;

        /**
         * @brief Accesses a reaction by its ID.
         * @param id The ID of the reaction to access.
         * @return A const reference to the Reaction.
         * @throws std::out_of_range if no reaction with the given ID exists.
         */
        [[nodiscard]] const Reaction& operator[](const std::string_view& id) const;

        /**
         * @brief Compares this set with another for equality.
         * @param other The other ReactionSet to compare with.
         * @return True if the sets are equal (same size and hash).
         */
        bool operator==(const ReactionSet& other) const;

        /**
         * @brief Compares this set with another for inequality.
         * @param other The other ReactionSet to compare with.
         * @return True if the sets are not equal.
         */
        bool operator!=(const ReactionSet& other) const;

        /**
         * @brief Computes a hash for the entire set.
         * @param seed The seed for the hash function.
         * @return A 64-bit hash value.
         * @details The algorithm computes the hash of each individual reaction,
         * sorts the hashes, and then computes a final hash over the sorted list
         * of hashes. This ensures the hash is order-independent.
         */
        [[nodiscard]] uint64_t hash(uint64_t seed) const;

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

        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> getReactionSetSpecies() const;

        friend std::ostream& operator<<(std::ostream& os, const ReactionSet& rs) {
            os << "(ReactionSet: {";
            size_t i = 0;
            for (const auto& reaction : rs.m_reactions) {
                os << *reaction;
                if (i < rs.m_reactions.size() - 1) {
                    os << ", ";
                }
                ++i;
            }
            os << "})";
            return os;
        }
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<std::unique_ptr<Reaction>> m_reactions;
        std::string m_id;
        std::unordered_map<std::string, size_t> m_reactionNameMap; ///< Maps reaction IDs to Reaction objects for quick lookup.
        std::unordered_set<size_t> m_reactionHashes;
        mutable std::optional<uint64_t> m_hashCache = std::nullopt;

    };

    ReactionSet packReactionSet(const ReactionSet& reactionSet);

}
