#pragma once

#include <string_view>

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include <unordered_map>
#include <vector>
#include <memory>
#include <unordered_set>
#include <cstdint>


#include "cppad/cppad.hpp"

namespace gridfire::reaction {
    /**
     * @struct REACLIBRateCoefficientSet
     * @brief Holds the seven fitting parameters for a single REACLIB rate set.
     * @details The thermonuclear reaction rate for a single set is calculated as:
     * rate = exp(a0 + a1/T9 + a2/T9^(-1/3) + a3*T9^(1/3) + a4*T9 + a5*T9^(5/3) + a6*ln(T9))
     * where T9 is the temperature in billions of Kelvin. The total rate for a
     * reaction is the sum of the rates from all its sets.
     */

    struct REACLIBRateCoefficientSet {
        const double a0;
        const double a1;
        const double a2;
        const double a3;
        const double a4;
        const double a5;
        const double a6;

        friend std::ostream& operator<<(std::ostream& os, const REACLIBRateCoefficientSet& r) {
            os << "[" << r.a0 << ", " << r.a1 << ", " << r.a2 << ", "
               << r.a3 << ", " << r.a4 << ", " << r.a5 << ", " << r.a6 << "]";
            return os;
        }
    };

    class Reaction {
    public:
        Reaction(
            const std::string_view id,
            const double qValue,
            const std::vector<fourdst::atomic::Species>& reactants,
            const std::vector<fourdst::atomic::Species>& products,
            const bool reverse = false
        );
        virtual ~Reaction() = default;

        virtual std::unique_ptr<Reaction> clone() const = 0;

        virtual double calculate_rate(double T9) const = 0;
        virtual CppAD::AD<double> calculate_rate(const CppAD::AD<double> T9) const = 0;

        virtual std::string_view peName() const { return ""; }

        [[nodiscard]] bool contains(const fourdst::atomic::Species& species) const;
        [[nodiscard]] bool contains_reactant(const fourdst::atomic::Species& species) const;
        [[nodiscard]] bool contains_product(const fourdst::atomic::Species& species) const;

        std::unordered_set<fourdst::atomic::Species> all_species() const;
        std::unordered_set<fourdst::atomic::Species> reactant_species() const;
        std::unordered_set<fourdst::atomic::Species> product_species() const;

        size_t num_species() const;

        int stoichiometry(const fourdst::atomic::Species& species) const;
        std::unordered_map<fourdst::atomic::Species, int> stoichiometry() const;

        std::string_view id() const { return m_id; }
        double qValue() const { return m_qValue; }
        const std::vector<fourdst::atomic::Species>& reactants() const { return m_reactants; }
        const std::vector<fourdst::atomic::Species>& products() const { return m_products; }
        bool is_reverse() const { return m_reverse; }

        double excess_energy() const;

        bool operator==(const Reaction& other) const { return m_id == other.m_id; }
        bool operator!=(const Reaction& other) const { return !(*this == other); }
        [[nodiscard]] uint64_t hash(uint64_t seed) const;

    protected:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::string m_id;
        double m_qValue = 0.0; ///< Q-value of the reaction
        std::vector<fourdst::atomic::Species> m_reactants; ///< Reactants of the reaction
        std::vector<fourdst::atomic::Species> m_products; ///< Products of the reaction
        bool m_reverse = false;
    };

    class ReactionSet {
    public:
        explicit ReactionSet(std::vector<std::unique_ptr<Reaction>> reactions);
        ReactionSet(const ReactionSet& other);
        ReactionSet& operator=(const ReactionSet& other);
        virtual ~ReactionSet() = default;

        virtual void add_reaction(std::unique_ptr<Reaction> reaction);
        virtual void remove_reaction(const std::unique_ptr<Reaction>& reaction);

        bool contains(const std::string_view& id) const;
        bool contains(const Reaction& reaction) const;

        size_t size() const { return m_reactions.size(); }

        void sort(double T9=1.0);

        bool contains_species(const fourdst::atomic::Species& species) const;
        bool contains_reactant(const fourdst::atomic::Species& species) const;
        bool contains_product(const fourdst::atomic::Species& species) const;

        [[nodiscard]] const Reaction& operator[](size_t index) const;
        [[nodiscard]] const Reaction& operator[](const std::string_view& id) const;

        bool operator==(const ReactionSet& other) const;
        bool operator!=(const ReactionSet& other) const;

        [[nodiscard]] uint64_t hash(uint64_t seed = 0) const;

        auto begin() { return m_reactions.begin(); }
        auto begin() const { return m_reactions.cbegin(); }
        auto end() { return m_reactions.end(); }
        auto end() const { return m_reactions.cend(); }
    protected:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<std::unique_ptr<Reaction>> m_reactions;
        std::string m_id;
        std::unordered_map<std::string, Reaction*> m_reactionNameMap; ///< Maps reaction IDs to Reaction objects for quick lookup

    };

    /**
     * @struct REACLIBReaction
     * @brief Represents a single nuclear reaction from the JINA REACLIB database.
     * @details This struct is designed to be constructed at compile time (constexpr) from
     * the data parsed by the Python generation script. It stores all necessary
     * information to identify a reaction and calculate its rate.
     */
    class REACLIBReaction final : public Reaction {
    public:
        /**
         * @brief Constructs a REACLIBReaction object at compile time.
         * @param id A unique string identifier generated by the Python script.
         * @param peName
         * @param chapter The REACLIB chapter number, defining the reaction structure.
         * @param reactants A vector of strings with the names of the reactant species.
         * @param products A vector of strings with the names of the product species.
         * @param qValue The Q-value of the reaction in MeV.
         * @param label The source label for the rate data (e.g., "wc12w", "st08").
         * @param sets A vector of RateFitSet, containing the fitting coefficients for the rate.
         * @param reverse A boolean indicating if the reaction is reversed (default is false).
         */
        REACLIBReaction(
            const std::string_view id,
            const std::string_view peName,
            const int chapter,
            const std::vector<fourdst::atomic::Species> &reactants,
            const std::vector<fourdst::atomic::Species> &products,
            const double qValue,
            const std::string_view label,
            const REACLIBRateCoefficientSet &sets,
            const bool reverse = false);

        [[nodiscard]] std::unique_ptr<Reaction> clone() const override;

        [[nodiscard]] double calculate_rate(const double T9) const override;
        [[nodiscard]] CppAD::AD<double> calculate_rate(const CppAD::AD<double> T9) const override;

        template <typename GeneralScalarType>
        [[nodiscard]] GeneralScalarType calculate_rate(const GeneralScalarType T9) const {
            const GeneralScalarType T913 = CppAD::pow(T9, 1.0/3.0);
            const GeneralScalarType rateExponent = m_rateCoefficients.a0 +
                m_rateCoefficients.a1 / T9 +
                m_rateCoefficients.a2 / T913 +
                m_rateCoefficients.a3 * T913 +
                m_rateCoefficients.a4 * T9 +
                m_rateCoefficients.a5 * CppAD::pow(T9, 5.0/3.0) +
                m_rateCoefficients.a6 * CppAD::log(T9);
            return CppAD::exp(rateExponent);
        }

        [[nodiscard]] std::string_view peName() const override { return m_peName; }

        [[nodiscard]] int chapter() const { return m_chapter; }

        [[nodiscard]] std::string_view sourceLabel() const { return m_sourceLabel; }

        [[nodiscard]] const REACLIBRateCoefficientSet& rateCoefficients() const { return m_rateCoefficients; }

        friend std::ostream& operator<<(std::ostream& os, const REACLIBReaction& reaction);
    private:
        std::string m_peName; ///< Name of the reaction in (projectile, ejectile)  notation (e.g. p(p, g)d)
        int m_chapter;    ///< Chapter number from the REACLIB database, defining the reaction structure.
        std::string m_sourceLabel; ///< Source label for the rate data, indicating the origin of the rate coefficients (e.g., "wc12w", "st08").
        REACLIBRateCoefficientSet m_rateCoefficients;
    };

    class REACLIBReactionSet final : public ReactionSet {
    public:
        explicit REACLIBReactionSet(std::vector<REACLIBReaction>);
        std::unordered_set<std::string> peNames() const;
        friend std::ostream& operator<<(std::ostream& os, const REACLIBReactionSet& set);
    };

    class REACLIBLogicalReaction final : public Reaction {
    public:
        explicit REACLIBLogicalReaction(const std::vector<REACLIBReaction> &reactions);
        explicit REACLIBLogicalReaction(const REACLIBReaction &reaction);
        void add_reaction(const REACLIBReaction& reaction);

        [[nodiscard]] std::unique_ptr<Reaction> clone() const override;

        [[nodiscard]] std::string_view peName() const override { return m_id; };
        [[nodiscard]] size_t size() const { return m_rates.size(); }

        [[nodiscard]] std::vector<std::string> sources() const { return m_sources; }

        [[nodiscard]] double calculate_rate(const double T9) const override;

        [[nodiscard]] CppAD::AD<double> calculate_rate(const CppAD::AD<double> T9) const override;

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

        [[nodiscard]] int chapter() const { return m_chapter; }

        auto begin() { return m_rates.begin(); }
        auto begin() const { return m_rates.cbegin(); }
        auto end() { return m_rates.end(); }
        auto end() const { return m_rates.cend(); }

    private:
        int m_chapter;
        std::vector<std::string> m_sources;
        std::vector<REACLIBRateCoefficientSet> m_rates;

    };

    class REACLIBLogicalReactionSet final : public ReactionSet {
    public:
        REACLIBLogicalReactionSet() = delete;
        explicit REACLIBLogicalReactionSet(const REACLIBReactionSet& reactionSet);

        [[nodiscard]] std::unordered_set<std::string> peNames() const;
    private:
        std::unordered_set<std::string> m_peNames;
    };

}