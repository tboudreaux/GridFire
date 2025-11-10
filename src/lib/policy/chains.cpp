#include "gridfire/policy/policy_abstract.h"
#include "gridfire/policy/policy_logical.h"
#include "gridfire/policy/chains.h"

#include "gridfire/exceptions/error_policy.h"

#include "gridfire/utils/hashing.h"
#include "gridfire/reaction/reaclib.h"

#include "xxhash64.h"


namespace gridfire::policy {
    TemperatureDependentChainPolicy::TemperatureDependentChainPolicy(
        const std::vector<std::string>& reactionIDs
    ) : TemperatureDependentChainPolicy(reactionIDs, std::nullopt, std::nullopt) {}

    TemperatureDependentChainPolicy::TemperatureDependentChainPolicy(
        const std::vector<std::string>& reactionIDs,
        const std::optional<double> minT9
    ) : TemperatureDependentChainPolicy(reactionIDs, minT9, std::nullopt) {}

    TemperatureDependentChainPolicy::TemperatureDependentChainPolicy(
        const std::vector<std::string>& reactionIDs,
        const std::optional<double> minT9,
        const std::optional<double> maxT9
    ) :
    m_reactionIDs(reactionIDs)
    {
        m_tempRange.minT9 = minT9;
        m_tempRange.maxT9 = maxT9;

        const auto& all_reaclib_reactions = reaclib::get_all_reaclib_reactions();

        for (const auto& reactionID : m_reactionIDs) {
            auto reaction = all_reaclib_reactions.get(reactionID);
            if (!reaction) {
                throw exceptions::MissingBaseReactionError("The Underlying REACLIB reaction set is missing the reaction " + std::string(reactionID) + " needed for a chain. This indicates that there is an issue with the GridFire binary you are using. Please try to recompile and if that fails please report this issue to the developers.");
            }
            m_reactions.add_reaction(reaction.value()->clone());
        }
    }

    const reaction::ReactionSet& TemperatureDependentChainPolicy::get_reactions() const {
        return m_reactions;
    }

    bool TemperatureDependentChainPolicy::contains(const std::string &id) const {
        return m_reactions.contains(id);
    }

    bool TemperatureDependentChainPolicy::contains(const reaction::Reaction &reaction) const {
        return m_reactions.contains(reaction);
    }

    uint64_t TemperatureDependentChainPolicy::hash(const uint64_t seed) const {
        std::vector<size_t> hashes;
        for (const auto& reaction : m_reactions) {
            hashes.push_back(reaction->hash(seed));
        }
        return XXHash64::hash(hashes.data(), hashes.size(), seed);
    }

    bool TemperatureDependentChainPolicy::operator==(const ReactionChainPolicy &other) const {
        return this->hash(0) == other.hash(0);
    }

    bool TemperatureDependentChainPolicy::operator!=(const ReactionChainPolicy &other) const {
        return this->hash(0) != other.hash(0);
    }

    bool TemperatureDependentChainPolicy::is_active(const double T9) const {
        return (!m_tempRange.minT9.has_value() || T9 >= m_tempRange.minT9.value()) &&
               (!m_tempRange.maxT9.has_value() || T9 <= m_tempRange.maxT9.value());
    }

    /**
     * Specific Implementations *
     **/

    ProtonProtonIChainPolicy::ProtonProtonIChainPolicy()
    : TemperatureDependentChainPolicy({
        // Proton-Proton I
        "p(p,e+)d",
        "d(p,g)he3",
        "he3(he3,2p)he4",
        },0.001) {}

    ProtonProtonIIChainPolicy::ProtonProtonIIChainPolicy()
    : TemperatureDependentChainPolicy({
        "p(p,e+)d",
        "d(p,g)he3",
        "he4(he3,g)be7",
        "be7(e-,)li7",
        "li7(p,a)he4",
        }, 0.001) {}

    ProtonProtonIIIChainPolicy::ProtonProtonIIIChainPolicy()
    : TemperatureDependentChainPolicy({
        "p(p,e+)d",
        "d(p,g)he3",
        "he4(he3,g)be7",
        "be7(p,g)b8",
        "b8(,e+ a)he4"
        }, 0.001) {}

    std::unique_ptr<ReactionChainPolicy> ProtonProtonIChainPolicy::clone() const {
        return std::make_unique<ProtonProtonIChainPolicy>(*this);
    }

    std::string ProtonProtonIChainPolicy::name() const {
        return "ProtonProtonIChainPolicy";
    }


    std::unique_ptr<ReactionChainPolicy> ProtonProtonIIChainPolicy::clone() const {
        return std::make_unique<ProtonProtonIIChainPolicy>(*this);
    }

    std::string ProtonProtonIIChainPolicy::name() const {
        return "ProtonProtonIIChainPolicy";
    }


    std::unique_ptr<ReactionChainPolicy> ProtonProtonIIIChainPolicy::clone() const {
        return std::make_unique<ProtonProtonIIIChainPolicy>(*this);
    }

    std::string ProtonProtonIIIChainPolicy::name() const {
        return "ProtonProtonIIIChainPolicy";
    }

    ProtonProtonChainPolicy::ProtonProtonChainPolicy() :
    MultiReactionChainPolicy(
        []() {
            std::vector<std::unique_ptr<ReactionChainPolicy>> chain_policies;
            chain_policies.push_back(std::make_unique<ProtonProtonIChainPolicy>());
            chain_policies.push_back(std::make_unique<ProtonProtonIIChainPolicy>());
            chain_policies.push_back(std::make_unique<ProtonProtonIIIChainPolicy>());
            return chain_policies;
        }()
    ){}

    std::string ProtonProtonChainPolicy::name() const {
        return "ProtonProtonChainPolicy";
    }

    CNOIChainPolicy::CNOIChainPolicy()
    : TemperatureDependentChainPolicy({
        "c12(p,g)n13",
        "n13(,e+)c13",
        "c13(p,g)n14",
        "n14(p,g)o15",
        "o15(,e+)n15",
        "n15(p,a)c12",
        }, 0.001){}

    CNOIIChainPolicy::CNOIIChainPolicy()
    : TemperatureDependentChainPolicy({
        "n15(p,g)o16",
        "o16(p,g)f17",
        "f17(,e+)o17",
        "o17(p,a)n14",
        "n14(p,g)o15",
        "o15(,e+)n15",
        }, 0.001){}

    CNOIIIChainPolicy::CNOIIIChainPolicy()
    : TemperatureDependentChainPolicy({
        "o17(p,g)f18",
        "f18(,e+)o18",
        "o18(p,a)n15",
        "n15(p,g)o16",
        "o16(p,g)f17",
        "f17(,e+)o17",
        }, 0.001){}

    CNOIVChainPolicy::CNOIVChainPolicy()
    : TemperatureDependentChainPolicy({
        "o18(p,g)f19",
        "f19(p,a)o16",
        "o16(p,g)f17",
        "f17(,e+)o17",
        "o17(p,g)f18",
        "f18(,e+)o18"
        }, 0.001){}


    std::unique_ptr<ReactionChainPolicy> CNOIChainPolicy::clone() const {
        return std::make_unique<CNOIChainPolicy>(*this);
    }

    std::string CNOIChainPolicy::name() const {
        return "CNOIChainPolicy";
    }

    std::unique_ptr<ReactionChainPolicy> CNOIIChainPolicy::clone() const {
        return std::make_unique<CNOIIChainPolicy>(*this);
    }

    std::string CNOIIChainPolicy::name() const {
        return "CNOIIChainPolicy";
    }

    std::unique_ptr<ReactionChainPolicy> CNOIIIChainPolicy::clone() const {
        return std::make_unique<CNOIIIChainPolicy>(*this);
    }

    std::string CNOIIIChainPolicy::name() const {
        return "CNOIIIChainPolicy";
    }

    std::unique_ptr<ReactionChainPolicy> CNOIVChainPolicy::clone() const {
        return std::make_unique<CNOIVChainPolicy>(*this);
    }

    std::string CNOIVChainPolicy::name() const {
        return "CNOIVChainPolicy";
    }

    CNOChainPolicy::CNOChainPolicy() :
    MultiReactionChainPolicy(
        []() {
            std::vector<std::unique_ptr<ReactionChainPolicy>> chain_policies;
            chain_policies.push_back(std::make_unique<CNOIChainPolicy>());
            chain_policies.push_back(std::make_unique<CNOIIChainPolicy>());
            chain_policies.push_back(std::make_unique<CNOIIIChainPolicy>());
            chain_policies.push_back(std::make_unique<CNOIVChainPolicy>());

            return chain_policies;
        }()
    ){}

    std::string CNOChainPolicy::name() const {
        return "CNOChainPolicy";
    }

    HotCNOIChainPolicy::HotCNOIChainPolicy()
    : TemperatureDependentChainPolicy({
        "c12(p,g)n13",
        "n13(p,g)o14",
        "o14(,e+)n14",
        "n14(p,g)o15",
        "o15(,e+)n15",
        "n15(p,a)c12",
        }, 0.1) {}

    HotCNOIIChainPolicy::HotCNOIIChainPolicy()
    : TemperatureDependentChainPolicy({
        "n15(p,g)o16",
        "o16(p,g)f17",
        "f17(p,g)ne18",
        "ne18(,e+)f18",
        "f18(p,a)o15",
        "o15(,e+)n15",
        }, 0.1) {}

    HotCNOIIIChainPolicy::HotCNOIIIChainPolicy()
    : TemperatureDependentChainPolicy({
        "f18(p,g)ne19",
        "ne19(,e+)f19",
        "f19(p,a)o16",
        "o16(p,g)f17",
        "f17(p,g)ne18",
        "ne18(,e+)f18"
        }, 0.1) {}

    std::unique_ptr<ReactionChainPolicy> HotCNOIChainPolicy::clone() const {
        return std::make_unique<HotCNOIChainPolicy>(*this);
    }

    std::string HotCNOIChainPolicy::name() const {
        return "HotCNOIChainPolicy";
    }

    std::unique_ptr<ReactionChainPolicy> HotCNOIIChainPolicy::clone() const {
        return std::make_unique<HotCNOIIChainPolicy>(*this);
    }

    std::string HotCNOIIChainPolicy::name() const {
        return "HotCNOIIChainPolicy";
    }

    std::unique_ptr<ReactionChainPolicy> HotCNOIIIChainPolicy::clone() const {
        return std::make_unique<HotCNOIIIChainPolicy>(*this);
    }

    std::string HotCNOIIIChainPolicy::name() const {
        return "HotCNOIIIChainPolicy";
    }

    HotCNOChainPolicy::HotCNOChainPolicy()
    : MultiReactionChainPolicy(
        []() {
            std::vector<std::unique_ptr<ReactionChainPolicy>> chain_policies;
            chain_policies.push_back(std::make_unique<HotCNOIChainPolicy>());
            chain_policies.push_back(std::make_unique<HotCNOIIChainPolicy>());
            chain_policies.push_back(std::make_unique<HotCNOIIIChainPolicy>());

            return chain_policies;
        }()
    ){}

    std::string HotCNOChainPolicy::name() const {
        return "HotCNOChainPolicy";
    }

    TripleAlphaChainPolicy::TripleAlphaChainPolicy()
    : TemperatureDependentChainPolicy({
        "he4(he4,a)be8",
        "be8(he4,g)c12"
        }, 0.01) {}

    std::unique_ptr<ReactionChainPolicy> TripleAlphaChainPolicy::clone() const {
        return std::make_unique<TripleAlphaChainPolicy>(*this);
    }

    std::string TripleAlphaChainPolicy::name() const {
        return "TripleAlphaChainPolicy";
    }

    MainSequenceReactionChainPolicy::MainSequenceReactionChainPolicy()
    : MultiReactionChainPolicy(
        []() {
            std::vector<std::unique_ptr<ReactionChainPolicy>> chain_policies;
            chain_policies.push_back(std::make_unique<ProtonProtonChainPolicy>());
            chain_policies.push_back(std::make_unique<CNOChainPolicy>());

            return chain_policies;
        }()
    ){}

    std::string MainSequenceReactionChainPolicy::name() const {
        return "MainSequenceReactionChainPolicy";
    }
}
