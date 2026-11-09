#include "gridfire/utils/provenance.h"

#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/weak/weak.h"
#include "gridfire/utils/config.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_view_abstract.h"

#include "xxhash64.h"

#include <algorithm>
#include <format>
#include <string>
#include <vector>

namespace gridfire {
    namespace {
        /// The GraphEngine at the bottom of a stack of views (null if there is none)
        const engine::GraphEngine* find_graph_engine(const engine::DynamicEngine* engine) {
            while (engine != nullptr) {
                if (const auto* graph = dynamic_cast<const engine::GraphEngine*>(engine)) {
                    return graph;
                }
                const auto* view = dynamic_cast<const engine::EngineView<engine::DynamicEngine>*>(engine);
                engine = view != nullptr ? &view->getBaseEngine() : nullptr;
            }
            return nullptr;
        }

        std::string_view screening_name(const screening::ScreeningType type) {
            switch (type) {
                case screening::ScreeningType::BARE: return "BARE (no screening)";
                case screening::ScreeningType::WEAK: return "WEAK (Salpeter 1954)";
                case screening::ScreeningType::CD09: return "CD09 (Chugunov & DeWitt 2009)";
                case screening::ScreeningType::BLENDED: return "BLENDED (Salpeter 1954 -> Chugunov & DeWitt 2009)";
            }
            return "UNKNOWN";
        }

        std::string_view jacobian_name(const engine::JacobianMethod method) {
            switch (method) {
                case engine::JacobianMethod::ANALYTIC: return "ANALYTIC";
                case engine::JacobianMethod::AUTOMATIC_DIFFERENTIATION: return "AUTOMATIC_DIFFERENTIATION (CppAD)";
            }
            return "UNKNOWN";
        }
    }

    std::string version_string() {
        return version::toString();
    }

    uint64_t network_hash(const reaction::ReactionSet& reactions) {
        std::vector<std::string> ids;
        ids.reserve(reactions.size());
        for (const auto& reaction : reactions) {
            ids.emplace_back(reaction->id());
        }
        std::ranges::sort(ids);
        XXHash64 hasher(0);
        for (const auto& id : ids) {
            hasher.add(id.data(), id.size());
            hasher.add("\n", 1);
        }
        return hasher.hash();
    }

    std::string describe_provenance(
        const config::GridFireConfig& config,
        const engine::DynamicEngine* engine,
        engine::scratch::StateBlob* ctx,
        const std::string_view policy_name
    ) {
        std::string out;
        const auto line = [&out](const std::string_view key, const std::string& value) {
            out += std::format("{}: {}\n", key, value);
        };

        line("gridfire_version", version_string());
        if (!policy_name.empty()) {
            line("policy", std::string(policy_name));
        }
        if (engine != nullptr && ctx != nullptr) {
            const auto& reactions = engine->getNetworkReactions(*ctx);
            line("network", std::format("{} species, {} reactions, hash {:016x}",
                engine->getNetworkSpecies(*ctx).size(), reactions.size(), network_hash(reactions)));
        }
        line("strong_rates", std::format("{}, data hash {:016x}", reaclib::DATA_SOURCE, reaclib::data_hash()));
        line("weak_rates", std::format("{}, data hash {:016x}", rates::weak::DATA_SOURCE, rates::weak::data_hash()));
        std::string multipliers = "none";
        if (const auto* graph = find_graph_engine(engine)) {
            const auto factors = graph->getRateMultipliers();
            for (size_t i = 0; i < factors.size(); ++i) {
                multipliers = (i ? multipliers + ", " : std::string()) + std::format("{}={}", factors[i].first, factors[i].second);
            }
        }
        line("rate_multipliers", multipliers);

        const auto& eng = config.engine;
        std::string screening(screening_name(eng.screening));
        if (eng.screening == screening::ScreeningType::BLENDED) {
            screening += std::format(", handover for weak-screening exponent {} to {}",
                eng.screeningBlend.weakExponent, eng.screeningBlend.strongExponent);
        }
        line("screening", screening);
        line("jacobian", std::string(jacobian_name(eng.jacobian)));
        line("reverse_reactions", eng.useReverseReactions ? "on (detailed balance, two-body strong reactions)" : "off");

        std::string passive = eng.passiveSpecies.empty() ? "none" : "";
        for (size_t i = 0; i < eng.passiveSpecies.size(); ++i) {
            passive += (i ? ", " : "") + eng.passiveSpecies[i];
        }
        line("passive_species", passive);
        line("energy_channels", eng.energyChannels.splitHydrogenPPCNO ? "hydrogen split into pp and CNO" : "hydrogen as one channel");

        const auto& point = config.solver.pointSolver;
        line("solver", std::format("CVODE, abs_tol {}, rel_tol {}, clamp_negative_iterates {}, max_step_failure_restarts {}",
            point.absTol, point.relTol, point.clampNegativeIterates, point.maxStepFailureRestarts));
        return out;
    }
}
