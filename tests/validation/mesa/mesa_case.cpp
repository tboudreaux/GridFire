#include "mesa_case.h"

#include "gridfire/partition/composite/partition_composite.h"

#include "fourdst/composition/utils/utils.h"

#include "nlohmann/json.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <stdexcept>

namespace gridfire::test::mesa {

    std::vector<double> BurnCase::step_durations() const {
        std::vector<double> durations;
        durations.reserve(num_times);
        for (int i = 0; i < num_times; ++i) {
            const double t = min_log_time + i * (max_log_time - min_log_time) / (num_times - 1);
            durations.push_back(std::pow(10.0, t));
        }
        return durations;
    }

    std::filesystem::path validation_dir() {
        return GF_MESA_VALIDATION_DIR;
    }

    std::vector<BurnCase> load_cases(const std::filesystem::path& cases_file) {
        std::ifstream in(cases_file);
        if (!in) {
            throw std::runtime_error("Cannot open MESA validation cases file " + cases_file.string());
        }
        const nlohmann::json doc = nlohmann::json::parse(in);

        std::vector<BurnCase> cases;
        for (const auto& c : doc.at("cases")) {
            BurnCase burn_case;
            burn_case.name = c.at("name").get<std::string>();
            burn_case.logT = c.at("logT").get<double>();
            burn_case.logRho = c.at("logRho").get<double>();
            burn_case.min_log_time = c.at("min_log_time").get<double>();
            burn_case.max_log_time = c.at("max_log_time").get<double>();
            burn_case.num_times = c.at("num_times").get<int>();
            for (const auto& [species, X] : c.at("composition").items()) {
                burn_case.mass_fractions[species] = X.get<double>();
            }
            cases.push_back(std::move(burn_case));
        }
        return cases;
    }

    CaseNetwork build_network(const BurnCase& burn_case, const std::vector<std::string>& excluded) {
        std::vector<std::string> symbols;
        std::vector<double> X;
        for (const auto& [species, x] : burn_case.mass_fractions) {
            symbols.push_back(species);
            X.push_back(x);
        }

        CaseNetwork network;
        network.initial_composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);

        // Same construction as MainSequencePolicy, but without its main-sequence-specific seed checks so that
        // non-hydrogen cases (e.g. helium burning) can be built too.
        network.partition_function = std::make_unique<partition::CompositePartitionFunction>(
            std::vector{partition::BasePartitionType::RauscherThielemann, partition::BasePartitionType::GroundState}
        );
        network.engine = std::make_unique<engine::GraphEngine>(
            network.initial_composition,
            *network.partition_function,
            engine::NetworkBuildDepth::ThirdOrder,
            engine::NetworkConstructionFlags::DEFAULT
        );
        if (!excluded.empty()) {
            auto full_context = network.engine->constructStateBlob(nullptr);
            reaction::ReactionSet kept;
            size_t removed = 0;
            for (const auto& reaction : network.engine->getNetworkReactions(*full_context)) {
                if (std::ranges::contains(excluded, std::string(reaction->id()))) {
                    ++removed;
                } else {
                    kept.add_reaction(*reaction);
                }
            }
            if (removed != excluded.size()) {
                throw std::runtime_error(std::format("{}: {} excluded reactions listed but only {} are in the network",
                                                     burn_case.name, excluded.size(), removed));
            }
            network.engine = std::make_unique<engine::GraphEngine>(kept, *network.partition_function);
        }

        network.context = network.engine->constructStateBlob(nullptr);
        return network;
    }

    std::vector<std::string> read_excluded_reactions(const std::filesystem::path& path) {
        std::vector<std::string> ids;
        std::ifstream in(path);
        for (std::string line; std::getline(in, line);) {
            const auto hash = line.find('#');
            std::string id = line.substr(0, hash);
            id.erase(0, id.find_first_not_of(" \t"));
            id.erase(id.find_last_not_of(" \t") + 1);
            if (!id.empty()) ids.push_back(id);
        }
        return ids;
    }

    std::vector<std::string> sorted_reaction_ids(const CaseNetwork& network) {
        std::vector<std::string> ids;
        for (const auto& reaction : network.engine->getNetworkReactions(*network.context)) {
            ids.emplace_back(reaction->id());
        }
        std::ranges::sort(ids);
        return ids;
    }
}
