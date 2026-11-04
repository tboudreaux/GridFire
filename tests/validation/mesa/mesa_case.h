/**
 * @file mesa_case.h
 * @brief Burn cases shared by the MESA/bbq reference generator and the GridFire comparison test.
 *
 * Both sides build each case's network through build_network(), so the network exported to MESA is exactly the one
 * GridFire integrates. The sorted reaction list is stored with the reference data; the test refuses to compare against
 * references produced for a different network.
 */
#pragma once

#include "gridfire/gridfire.h"
#include "gridfire/engine/engine_graph.h"

#include "fourdst/composition/composition.h"

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace gridfire::test::mesa {

    struct BurnCase {
        std::string name;
        double logT;
        double logRho;
        double min_log_time;
        double max_log_time;
        int num_times;
        std::map<std::string, double> mass_fractions; ///< Initial composition by GridFire species name

        /// @brief Duration of each step, matching bbq's hydrostatic grid: 10^(min + i (max - min) / (n - 1)).
        [[nodiscard]] std::vector<double> step_durations() const;
    };

    /// @brief Directory holding cases.json and reference/ (set at build time).
    std::filesystem::path validation_dir();

    std::vector<BurnCase> load_cases(const std::filesystem::path& cases_file);

    /// @brief A case's GraphEngine (same construction as MainSequencePolicy) plus its initial state.
    struct CaseNetwork {
        std::unique_ptr<partition::PartitionFunction> partition_function;
        std::unique_ptr<engine::GraphEngine> engine;
        std::unique_ptr<engine::scratch::StateBlob> context;
        fourdst::composition::Composition initial_composition;
    };

    /**
     * @brief Build a case's network, optionally without some reactions.
     *
     * @param excluded GridFire reaction ids to remove (reactions MESA cannot evaluate; recorded with the references).
     */
    CaseNetwork build_network(const BurnCase& burn_case, const std::vector<std::string>& excluded = {});

    /// @brief Reaction ids listed in an excluded_reactions.txt file ("<id>  # reason" per line); empty if absent.
    std::vector<std::string> read_excluded_reactions(const std::filesystem::path& path);

    /// @brief Reaction ids of the network, sorted; the identity of the network stored alongside references.
    std::vector<std::string> sorted_reaction_ids(const CaseNetwork& network);
}
