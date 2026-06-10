/**
 * @file chains.h
 * @brief Concrete implementations of ReactionChainPolicy for key stellar reaction chains.
 *
 * This file provides concrete policies for fundamental nuclear reaction chains, such as the
 * Proton-Proton chain and the CNO cycle. These policies inherit from ReactionChainPolicy (see
 * `policy_abstract.h`) and provide a pre-defined set of reactions.
 *
 * They are typically used by higher-level NetworkPolicy implementations (e.g., `MainSequencePolicy`
 * in `stellar_policy.h`) to compose a complete set of required reactions for a particular
 * stellar environment.
 *
 * @note Constructors for these policies may throw `gridfire::exceptions::MissingBaseReactionError`
 *       if a required reaction is not found in the underlying REACLIB database. This usually
 *       indicates an issue with the GridFire installation or the bundled reaction data.
 */
#pragma once

#include "gridfire/policy/policy_abstract.h"
#include "gridfire/policy/policy_logical.h"
#include "gridfire/reaction/reaction.h"

#include <memory>
#include <optional>

namespace gridfire::policy {
    /**
     * @brief Base class for reaction chain policies that are active only within specific temperature ranges.
     *
     * Such chains may only operate effectively within certain temperature regimes, reflecting
     * the physical conditions required for the reactions to proceed. This class allows defining
     * such temperature-dependent behavior. This is one of the locations where domain specific knowledge is allowed
     * within GridFire.
     */
    class TemperatureDependentChainPolicy : public ReactionChainPolicy {
    public:
        /**
         * @brief Construct a new Temperature Dependent Chain Policy object
         * @param reactionIDs Vector of reaction IDs that comprise the chain
         *
         * @note The chain is considered active at all temperatures if no min/max T9 are provided.
         */
        explicit TemperatureDependentChainPolicy(const std::vector<std::string>& reactionIDs);

        /**
         * @brief Construct a new Temperature Dependent Chain Policy object
         * @param reactionIDs Vector of reaction IDs that comprise the chain
         * @param minT9 Minimum temperature (in T9) for which the chain is active
         *
         * @note The chain is considered active at all temperatures above minT9.
         */
        explicit TemperatureDependentChainPolicy(const std::vector<std::string>& reactionIDs, std::optional<double> minT9);

        /**
         * @brief Construct a new Temperature Dependent Chain Policy object
         * @param reactionIDs Vector of reaction IDs that comprise the chain
         * @param minT9 Minimum temperature (in T9) for which the chain is active
         * @param maxT9 Maximum temperature (in T9) for which the chain is active
         */
        explicit TemperatureDependentChainPolicy(const std::vector<std::string>& reactionIDs, std::optional<double> minT9, std::optional<double> maxT9);

        /**
         * @brief Virtual destructor
         */
        ~TemperatureDependentChainPolicy() override = default;

        /**
         * @brief Get the reactions that comprise the chain
         * @return The set of reactions comprising the chain
         */
        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override;

        /**
         * @brief Check if a reaction or reaction ID is part of the chain
         * @param id The reaction ID to check
         * @return True if the reaction ID is part of the chain, false otherwise
         */
        [[nodiscard]] bool contains(const std::string& id) const override;

        /**
         * @brief Check if a reaction is part of the chain
         * @param reaction The reaction to check
         * @return True if the reaction is part of the chain, false otherwise
         */
        [[nodiscard]] bool contains(const reaction::Reaction& reaction) const override;

        /**
         * @brief Generate the hash for the reaction chain policy
         * @param seed Seed value for the hash
         * @return Unique hash representing the reaction chain policy. XXHash64 is used internally for speed and collision resistance.
         */
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;

        /**
         * @brief Equality operator
         * @param other The other ReactionChainPolicy to compare against
         * @return True if the two policies are equal, false otherwise
         */
        [[nodiscard]] bool operator==(const ReactionChainPolicy& other) const override;

        /**
         * @brief Inequality operator
         * @param other The other ReactionChainPolicy to compare against
         * @return True if the two policies are not equal, false otherwise
         */
        [[nodiscard]] bool operator!=(const ReactionChainPolicy& other) const override;

        /**
         * @brief Check if the reaction chain is active at a given temperature
         * @param T9 Temperature in billions of Kelvin (T9)
         * @return True if the chain is active at the given temperature, false otherwise
         */
        [[nodiscard]] bool is_active(double T9) const;

    protected:
        /**
         * @brief Struct to hold the active temperature range for the reaction chain
         */
        struct ActiveTempRange {
            std::optional<double> minT9; ///< Minimum temperature (in T9) for which the chain is active
            std::optional<double> maxT9; ///< Maximum temperature (in T9) for which the chain is active
        };

        ActiveTempRange m_tempRange; ///< Active temperature range for the reaction chain
        std::vector<std::string> m_reactionIDs; ///< Vector of reaction IDs that comprise the chain
        reaction::ReactionSet m_reactions; ///< Set of reactions that comprise the chain
    };

    /**
     * @brief Proton-Proton I Chain Policy
     *
     * This class implements the Proton-Proton I chain of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - p(p,e+)d
     *  - d(p,g)he3
     *  - he3(he3,2p)he4
     */
    class ProtonProtonIChainPolicy final: public TemperatureDependentChainPolicy {
    public:
        ProtonProtonIChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Proton-Proton II Chain Policy
     *
     * This class implements the Proton-Proton II chain of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - p(p,e+)d
     *  - d(p,g)he3
     *  - he4(he3,g)be7
     *  - be7(e-,)li7
     *  - li7(p,a)he4
     */
    class ProtonProtonIIChainPolicy final: public TemperatureDependentChainPolicy {
    public:
        ProtonProtonIIChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Proton-Proton III Chain Policy
     *
     * This class implements the Proton-Proton III chain of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - p(p,e+)d
     *  - d(p,g)he3
     *  - he4(he3,g)be7
     *  - be7(p,g)b8
     *  - b8(,e+ a)he4
     */
    class ProtonProtonIIIChainPolicy final: public TemperatureDependentChainPolicy {
    public:
        ProtonProtonIIIChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Proton-Proton Chain Policy
     *
     * This class implements the overall Proton-Proton chain of nuclear reactions, combining the
     * Proton-Proton I, II, and III chains. Enforcing this chain in the policy will ensure that all the
     * Proton-Proton reactions are included in the network.
     *
     * @see ProtonProtonIChainPolicy
     * @see ProtonProtonIIChainPolicy
     * @see ProtonProtonIIIChainPolicy
     */
    class ProtonProtonChainPolicy final : public MultiReactionChainPolicy {
    public:
        ProtonProtonChainPolicy();
        [[nodiscard]] std::string name() const override;
    private:
        std::vector<std::unique_ptr<ReactionChainPolicy>> m_chain_policies;
    };

    /**
     * @brief CNO I Chain Policy
     *
     * This class implements the CNO I cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - c12(p,g)n13
     *  - n13(,e+)c13
     *  - c13(p,g)n14
     *  - n14(p,g)o15
     *  - o15(,e+)n15
     *  - n15(p,a)c12
     */
    class CNOIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief CNO II Chain Policy
     *
     * This class implements the CNO II cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - n15(p,g)o16
     *  - o16(p,g)f17
     *  - f17(p,g)ne18
     *  - ne18(,e+)f18
     *  - f18(p,a)o15
     *  - o15(,e+)n15
     */
    class CNOIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief CNO III Chain Policy
     *
     * This class implements the CNO III cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - o17(p,g)f18
     *  - f18(,e+)o18
     *  - o18(p,a)n15
     *  - n15(p,g)o16
     *  - o16(p,g)f17
     *  - f17(,e+)o17
     */
    class CNOIIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief CNO IV Chain Policy
     *
     * This class implements the CNO IV cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.001 (or 1e6K). This chain includes reactions:
     *  - o18(p,g)f19
     *  - f19(p,a)o16
     *  - o16(p,g)f17
     *  - f17(,e+)o17
     *  - o17(p,g)f18
     *  - f18(,e+)o18
     */
    class CNOIVChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIVChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief CNO Chain Policy
     *
     * This class implements the overall CNO cycle of nuclear reactions, combining the
     * CNO I, II, III, and IV chains. Enforcing this chain in the policy will ensure that all the
     * CNO reactions are included in the network.
     *
     * @see CNOIChainPolicy
     * @see CNOIIChainPolicy
     * @see CNOIIIChainPolicy
     * @see CNOIVChainPolicy
     */
    class CNOChainPolicy final : public MultiReactionChainPolicy {
    public:
        CNOChainPolicy();
        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Hot CNO I Chain Policy
     *
     * This class implements the Hot CNO I cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.1 (or 1e8K). This chain includes reactions:
     *  - c12(p,g)n13
     *  - n13(p,g)o14
     *  - o14(,e+)n14
     *  - n14(p,g)o15
     *  - o15(,e+)n15
     *  - n15(p,a)c12
     */
    class HotCNOIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        HotCNOIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Hot CNO II Chain Policy
     *
     * This class implements the Hot CNO II cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.1 (or 1e8K). This chain includes reactions:
     *  - n15(p,g)o16
     *  - o16(p,g)f17
     *  - f17(p,g)ne18
     *  - ne18(,e+)f18
     *  - f18(p,a)o15
     *  - o15(,e+)n15
     */
    class HotCNOIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        HotCNOIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Hot CNO III Chain Policy
     *
     * This class implements the Hot CNO III cycle of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.1 (or 1e8K). This chain includes reactions:
     *  - f18(p,g)ne19
     *  - ne19(,e+)f19
     *  - f19(p,a)o16
     *  - o16(p,g)f17
     *  - f17(p,g)ne18
     *  - ne18(,e+)f18
     */
    class HotCNOIIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        HotCNOIIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Hot CNO Chain Policy
     *
     * This class implements the overall Hot CNO cycle of nuclear reactions, combining the
     * Hot CNO I, II, and III chains. Enforcing this chain in the policy will ensure that all the
     * Hot CNO reactions are included in the network.
     *
     * @see HotCNOIChainPolicy
     * @see HotCNOIIChainPolicy
     * @see HotCNOIIIChainPolicy
     */
    class HotCNOChainPolicy final : public MultiReactionChainPolicy {
    public:
        HotCNOChainPolicy();
        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Triple-Alpha Chain Policy
     *
     * This class implements the Triple-Alpha process of nuclear reactions. This chain's minimum temperature is
     * set to T9=0.01 (or 1e7K). This chain includes reactions:
     *  - he4(he4,a)be8
     *  - be8(he4,g)c12
     */
    class TripleAlphaChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        TripleAlphaChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    /**
     * @brief Main Sequence Reaction Chain Policy
     *
     * This class implements the main sequence reaction chains, combining the
     * Proton-Proton chain and the CNO cycle. Enforcing this chain in the policy will ensure that all the
     * primary reactions for main sequence stars are included in the network.
     *
     * @see ProtonProtonChainPolicy
     * @see CNOChainPolicy
     */
    class MainSequenceReactionChainPolicy final : public MultiReactionChainPolicy {
    public:
        MainSequenceReactionChainPolicy();

        [[nodiscard]] std::string name() const override;

    };
}
