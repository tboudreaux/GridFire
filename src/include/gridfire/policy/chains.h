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

namespace gridfire::policy {

    class TemperatureDependentChainPolicy : public ReactionChainPolicy {
    public:
        explicit TemperatureDependentChainPolicy(const std::vector<std::string>& reactionIDs);
        explicit TemperatureDependentChainPolicy(const std::vector<std::string>& reactionIDs, std::optional<double> minT9);
        explicit TemperatureDependentChainPolicy(const std::vector<std::string>& reactionIDs, std::optional<double> minT9, std::optional<double> maxT9);
        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override;
        [[nodiscard]] bool contains(const std::string& id) const override;
        [[nodiscard]] bool contains(const reaction::Reaction& reaction) const override;
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;
        [[nodiscard]] bool operator==(const ReactionChainPolicy& other) const override;
        [[nodiscard]] bool operator!=(const ReactionChainPolicy& other) const override;

        [[nodiscard]] bool is_active(double T9) const;

    protected:
        struct ActiveTempRange {
            std::optional<double> minT9;
            std::optional<double> maxT9;
        };

        ActiveTempRange m_tempRange;
        std::vector<std::string> m_reactionIDs;
        reaction::ReactionSet m_reactions;
    };

    class ProtonProtonIChainPolicy final: public TemperatureDependentChainPolicy {
    public:
        ProtonProtonIChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class ProtonProtonIIChainPolicy final: public TemperatureDependentChainPolicy {
    public:
        ProtonProtonIIChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class ProtonProtonIIIChainPolicy final: public TemperatureDependentChainPolicy {
    public:
        ProtonProtonIIIChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class ProtonProtonChainPolicy final : public MultiReactionChainPolicy {
    public:
        ProtonProtonChainPolicy();
        [[nodiscard]] std::string name() const override;
    private:
        std::vector<std::unique_ptr<ReactionChainPolicy>> m_chain_policies;
    };

    class CNOIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class CNOIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class CNOIIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class CNOIVChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        CNOIVChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };

    class CNOChainPolicy final : public MultiReactionChainPolicy {
    public:
        CNOChainPolicy();
        [[nodiscard]] std::string name() const override;
    };

    class HotCNOIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        HotCNOIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
    };

    class HotCNOIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        HotCNOIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
    };

    class HotCNOIIIChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        HotCNOIIIChainPolicy();
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
    };

    class HotCNOChainPolicy final : public MultiReactionChainPolicy {
    public:
        HotCNOChainPolicy();
        [[nodiscard]] std::string name() const override;
    };

    class TripleAlphaChainPolicy final : public TemperatureDependentChainPolicy {
    public:
        TripleAlphaChainPolicy();

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        [[nodiscard]] std::string name() const override;
    };


    class MainSequenceReactionChainPolicy final : public MultiReactionChainPolicy {
    public:
        MainSequenceReactionChainPolicy();

        [[nodiscard]] std::string name() const override;

    };
}
