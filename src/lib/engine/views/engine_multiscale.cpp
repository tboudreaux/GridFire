#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/exceptions/error_engine.h"
#include "gridfire/engine/procedures/priming.h"

#include <stdexcept>
#include <vector>
#include <ranges>
#include <unordered_map>
#include <unordered_set>

#include <queue>

#include <algorithm>

#include "quill/LogMacros.h"
#include "quill/Logger.h"

namespace {
    using namespace fourdst::atomic;
    //TODO: Replace all calls to this function with composition.getMolarAbundanceVector() so that
    //      we don't have to keep this function around. (Cant do this right now because there is not a
    //      guarantee that this function will return the same ordering as the canonical vector representation)
    std::vector<double> packCompositionToVector(
        const fourdst::composition::Composition& composition,
        const gridfire::DynamicEngine& engine
    ) {
        std::vector<double> Y(engine.getNetworkSpecies().size(), 0.0);
        const auto& allSpecies = engine.getNetworkSpecies();
        for (size_t i = 0; i < allSpecies.size(); ++i) {
            const auto& species = allSpecies[i];
            if (!composition.contains(species)) {
                Y[i] = 0.0; // Species not in the composition, set to zero
            } else {
                Y[i] = composition.getMolarAbundance(species);
            }
        }
        return Y;
    }

    template <class T>
    void hash_combine(std::size_t& seed, const T& v) {
        std::hash<T> hashed;
        seed ^= hashed(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    std::vector<std::vector<Species>> findConnectedComponentsBFS(
        const std::unordered_map<Species, std::vector<Species>>& graph,
        const std::vector<Species>& nodes
    ) {
        std::vector<std::vector<Species>> components;
        std::unordered_set<Species> visited;

        for (const Species& start_node : nodes) {
            if (!visited.contains(start_node)) {
                std::vector<Species> current_component;
                std::queue<Species> q;

                q.push(start_node);
                visited.insert(start_node);

                while (!q.empty()) {
                    Species u = q.front();
                    q.pop();
                    current_component.push_back(u);

                    if (graph.contains(u)) {
                        for (const auto& v : graph.at(u)) {
                            if (!visited.contains(v)) {
                                visited.insert(v);
                                q.push(v);
                            }
                        }
                    }
                }
                components.push_back(current_component);
            }
        }
        return components;
    }

    struct SpeciesSetIntersection {
        const Species species;
        std::size_t count;
    };

    std::expected<SpeciesSetIntersection, std::string> get_intersection_info (
        const std::unordered_set<Species>& setA,
        const std::unordered_set<Species>& setB
    ) {
        // Iterate over the smaller of the two
        auto* outerSet = &setA;
        auto* innerSet = &setB;
        if (setA.size() > setB.size()) {
            outerSet = &setB;
            innerSet = &setA;
        }

        std::size_t matchCount = 0;
        const Species* firstMatch = nullptr;

        for (const Species& sp : *outerSet) {
            if (innerSet->contains(sp)) {
                if (matchCount == 0) {
                    firstMatch = &sp;
                }
                ++matchCount;
                if (matchCount > 1) {
                    break;
                }
            }
        }
        if (!firstMatch) {
            // No matches found
            return std::unexpected{"Intersection is empty"};
        }
        if (matchCount == 0) {
            // No matches found
            return std::unexpected{"No intersection found"};
        }

        // Return the first match and the count of matches
        return SpeciesSetIntersection{*firstMatch, matchCount};

    }

    bool has_distinct_reactant_and_product_species (
        const std::unordered_set<Species>& poolSpecies,
        const std::unordered_set<Species>& reactants,
        const std::unordered_set<Species>& products
    ) {
        const auto reactant_result = get_intersection_info(poolSpecies, reactants);
        if (!reactant_result) {
            return false; // No reactants found
        }
        const auto [reactantSample, reactantCount] = reactant_result.value();

        const auto product_result = get_intersection_info(poolSpecies, products);
        if (!product_result) {
            return false; // No products found
        }

        const auto [productSample, productCount] = product_result.value();

        // If either side has ≥2 distinct matches, we can always pick
        // one from each that differ.
        if (reactantCount > 1 || productCount > 1) {
            return true;
        }

        // Exactly one match on each side → they must differ
        return reactantSample != productSample;
    }

    const std::unordered_map<Eigen::LevenbergMarquardtSpace::Status, std::string> lm_status_map = {
        {Eigen::LevenbergMarquardtSpace::Status::NotStarted, "NotStarted"},
        {Eigen::LevenbergMarquardtSpace::Status::Running, "Running"},
        {Eigen::LevenbergMarquardtSpace::Status::ImproperInputParameters, "ImproperInputParameters"},
        {Eigen::LevenbergMarquardtSpace::Status::RelativeReductionTooSmall, "RelativeReductionTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::RelativeErrorTooSmall, "RelativeErrorTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::RelativeErrorAndReductionTooSmall, "RelativeErrorAndReductionTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::CosinusTooSmall, "CosinusTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::TooManyFunctionEvaluation, "TooManyFunctionEvaluation"},
        {Eigen::LevenbergMarquardtSpace::Status::FtolTooSmall, "FtolTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::XtolTooSmall, "XtolTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::GtolTooSmall, "GtolTooSmall"},
        {Eigen::LevenbergMarquardtSpace::Status::UserAsked, "UserAsked"}
    };

}

namespace gridfire {
    using fourdst::atomic::Species;

    MultiscalePartitioningEngineView::MultiscalePartitioningEngineView(
        DynamicEngine& baseEngine
    ) : m_baseEngine(baseEngine) {}

    const std::vector<Species> & MultiscalePartitioningEngineView::getNetworkSpecies() const {
        return m_baseEngine.getNetworkSpecies();
    }

    std::expected<StepDerivatives<double>, expectations::StaleEngineError> MultiscalePartitioningEngineView::calculateRHSAndEnergy(
        const fourdst::composition::Composition& comp,
        const double T9,
        const double rho
    ) const {
        // Check the cache to see if the network needs to be repartitioned. Note that the QSECacheKey manages binning of T9, rho, and Y_full to ensure that small changes (which would likely not result in a repartitioning) do not trigger a cache miss.
        const auto result = m_baseEngine.calculateRHSAndEnergy(comp, T9, rho);
        if (!result) {
            return std::unexpected{result.error()};
        }
        auto deriv = result.value();

        for (const auto& species : m_algebraic_species) {
            deriv.dydt[species] = 0.0; // Fix the algebraic species to the equilibrium abundances we calculate.
        }
        return deriv;
    }

    EnergyDerivatives MultiscalePartitioningEngineView::calculateEpsDerivatives(
        const fourdst::composition::Composition& comp,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateEpsDerivatives(comp, T9, rho);
    }

    void MultiscalePartitioningEngineView::generateJacobianMatrix(
        const fourdst::composition::Composition& comp,
        const double T9,
        const double rho
    ) const {
        // TODO: Add sparsity pattern to this to prevent base engine from doing unnecessary work.
        m_baseEngine.generateJacobianMatrix(comp, T9, rho);
    }

    double MultiscalePartitioningEngineView::getJacobianMatrixEntry(
        const Species& rowSpecies,
        const Species& colSpecies
    ) const {
        // Check if the species we are differentiating with respect to is algebraic or dynamic. If it is algebraic we can reduce the work significantly...
        if (std::ranges::contains(m_algebraic_species, colSpecies)) {
            return 0.0;
        }
        if (std::ranges::contains(m_algebraic_species, rowSpecies)) {
            return 0.0;
        }
        return m_baseEngine.getJacobianMatrixEntry(rowSpecies, colSpecies);
    }

    void MultiscalePartitioningEngineView::generateStoichiometryMatrix() {
        m_baseEngine.generateStoichiometryMatrix();
    }

    int MultiscalePartitioningEngineView::getStoichiometryMatrixEntry(
        const Species& species,
        const reaction::Reaction& reaction
    ) const {
        return m_baseEngine.getStoichiometryMatrixEntry(species, reaction);
    }

    double MultiscalePartitioningEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        // Fix the algebraic species to the equilibrium abundances we calculate.
        fourdst::composition::Composition comp_mutable = comp;
        for (const auto& species : m_algebraic_species) {
            // TODO: Check this conversion to mass fraction (also consider adding the ability to set molar abundance directly)
            const double Yi = m_algebraic_abundances.at(species);
            comp_mutable.setMassFraction(species, Yi * species.a() / (rho * 1e-3)); // Convert Yi (mol/g) to Xi (mass fraction)
        }
        if (!comp_mutable.finalize()) {
            LOG_ERROR(m_logger, "Failed to finalize composition after setting algebraic species abundances.");
            m_logger->flush_log();
            throw std::runtime_error("Failed to finalize composition after setting algebraic species abundances.");
        }
        return m_baseEngine.calculateMolarReactionFlow(reaction, comp_mutable, T9, rho);
    }

    const reaction::ReactionSet & MultiscalePartitioningEngineView::getNetworkReactions() const {
        return m_baseEngine.getNetworkReactions();
    }

    void MultiscalePartitioningEngineView::setNetworkReactions(const reaction::ReactionSet &reactions) {
        LOG_CRITICAL(m_logger, "setNetworkReactions is not supported in MultiscalePartitioningEngineView. Did you mean to call this on the base engine?");
        throw exceptions::UnableToSetNetworkReactionsError("setNetworkReactions is not supported in MultiscalePartitioningEngineView. Did you mean to call this on the base engine?");
    }

    std::expected<std::unordered_map<Species, double>, expectations::StaleEngineError> MultiscalePartitioningEngineView::getSpeciesTimescales(
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        const auto result  = m_baseEngine.getSpeciesTimescales(comp, T9, rho);
        if (!result) {
            return std::unexpected{result.error()};
        }
        std::unordered_map<Species, double> speciesTimescales = result.value();
        for (const auto& algebraicSpecies : m_algebraic_species) {
            speciesTimescales[algebraicSpecies] = std::numeric_limits<double>::infinity(); // Algebraic species have infinite timescales.
        }
        return speciesTimescales;
    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError>
    MultiscalePartitioningEngineView::getSpeciesDestructionTimescales(
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        const auto result = m_baseEngine.getSpeciesDestructionTimescales(comp, T9, rho);
        if (!result) {
            return std::unexpected{result.error()};
        }
        std::unordered_map<Species, double> speciesDestructionTimescales = result.value();
        for (const auto& algebraicSpecies : m_algebraic_species) {
            speciesDestructionTimescales[algebraicSpecies] = std::numeric_limits<double>::infinity(); // Algebraic species have infinite destruction timescales.
        }
        return speciesDestructionTimescales;
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::update(const NetIn &netIn) {
        const fourdst::composition::Composition baseUpdatedComposition = m_baseEngine.update(netIn);

        NetIn baseUpdatedNetIn = netIn;
        baseUpdatedNetIn.composition = baseUpdatedComposition;
        const fourdst::composition::Composition equilibratedComposition = equilibrateNetwork(baseUpdatedNetIn);
        std::unordered_map<Species, double> algebraicAbundances;
        for (const auto& species : m_algebraic_species) {
            algebraicAbundances[species] = equilibratedComposition.getMolarAbundance(species);
        }

        m_algebraic_abundances = std::move(algebraicAbundances);

        return equilibratedComposition;
    }

    bool MultiscalePartitioningEngineView::isStale(const NetIn &netIn) {
        const auto key = QSECacheKey(
            netIn.temperature,
            netIn.density,
            packCompositionToVector(netIn.composition, m_baseEngine)
        );
        if (m_qse_abundance_cache.contains(key)) {
            return m_baseEngine.isStale(netIn); // The cache hit indicates the engine is not stale for the given conditions.
        }
        return true;
    }

    void MultiscalePartitioningEngineView::setScreeningModel(
        const screening::ScreeningType model
    ) {
        m_baseEngine.setScreeningModel(model);
    }

    screening::ScreeningType MultiscalePartitioningEngineView::getScreeningModel() const {
        return m_baseEngine.getScreeningModel();
    }

    const DynamicEngine & MultiscalePartitioningEngineView::getBaseEngine() const {
        return m_baseEngine;
    }

    std::vector<std::vector<Species>> MultiscalePartitioningEngineView::analyzeTimescalePoolConnectivity(
        const std::vector<std::vector<Species>> &timescale_pools,
        const fourdst::composition::Composition &comp,
        double T9,
        double rho
    ) const {
        std::vector<std::vector<Species>> final_connected_pools;

        for (const auto& pool : timescale_pools) {
            if (pool.empty()) {
                continue; // Skip empty pools
            }

            // For each timescale pool, we need to analyze connectivity.
            auto connectivity_graph = buildConnectivityGraph(pool);
            auto components = findConnectedComponentsBFS(connectivity_graph, pool);
            final_connected_pools.insert(final_connected_pools.end(), components.begin(), components.end());
        }
        return final_connected_pools;

    }

    void MultiscalePartitioningEngineView::partitionNetwork(
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) {
        // --- Step 0. Clear previous state ---
        LOG_TRACE_L1(m_logger, "Partitioning network...");
        LOG_TRACE_L1(m_logger, "Clearing previous state...");
        m_qse_groups.clear();
        m_dynamic_species.clear();
        m_algebraic_species.clear();

        // --- Step 1. Identify distinct timescale regions ---
        LOG_TRACE_L1(m_logger, "Identifying fast reactions...");
        const std::vector<std::vector<Species>> timescale_pools = partitionByTimescale(comp, T9, rho);
        LOG_TRACE_L1(m_logger, "Found {} timescale pools.", timescale_pools.size());

        // --- Step 2. Select the mean slowest pool as the base dynamical group ---
        LOG_TRACE_L1(m_logger, "Identifying mean slowest pool...");
        const size_t mean_slowest_pool_index = identifyMeanSlowestPool(timescale_pools, comp, T9, rho);
        LOG_TRACE_L1(m_logger, "Mean slowest pool index: {}", mean_slowest_pool_index);

        // --- Step 3. Push the slowest pool into the dynamic species list ---
        for (const auto& slowSpecies : timescale_pools[mean_slowest_pool_index]) {
            m_dynamic_species.push_back(slowSpecies);
        }

        // --- Step 4. Pack Candidate QSE Groups ---
        std::vector<std::vector<Species>> candidate_pools;
        for (size_t i = 0; i < timescale_pools.size(); ++i) {
            if (i == mean_slowest_pool_index) continue; // Skip the slowest pool
            LOG_TRACE_L1(m_logger, "Group {} with {} species identified for potential QSE.", i, timescale_pools[i].size());
            candidate_pools.push_back(timescale_pools[i]);
        }

        LOG_TRACE_L1(m_logger, "Preforming connectivity analysis on timescale pools...");
        const std::vector<std::vector<Species>> connected_pools = analyzeTimescalePoolConnectivity(candidate_pools, comp, T9, rho);
        LOG_TRACE_L1(m_logger, "Found {} connected pools (compared to {} timescale pools) for QSE analysis.", connected_pools.size(), timescale_pools.size());


        // --- Step 5. Identify potential seed species for each candidate pool ---
        LOG_TRACE_L1(m_logger, "Identifying potential seed species for candidate pools...");
        const std::vector<QSEGroup> candidate_groups = constructCandidateGroups(connected_pools, comp, T9, rho);
        LOG_TRACE_L1(m_logger, "Found {} candidate QSE groups for further analysis", candidate_groups.size());
        LOG_TRACE_L2(
            m_logger,
            "{}",
             [&]() -> std::string {
                 std::stringstream ss;
                 int j = 0;
                 for (const auto& group : candidate_groups) {
                     ss << "CandidateQSEGroup(Algebraic: {";
                     int i = 0;
                     for (const auto& species : group.algebraic_species) {
                         ss << species.name();
                         if (i < group.algebraic_species.size() - 1) {
                             ss << ", ";
                         }
                     }
                     ss << "}, Seed: {";
                     i = 0;
                     for (const auto& species : group.seed_species) {
                         ss << species.name();
                         if (i < group.seed_species.size() - 1) {
                             ss << ", ";
                         }
                         i++;
                     }
                     ss << "})";
                     if (j < candidate_groups.size() - 1) {
                         ss << ", ";
                     }
                     j++;
                 }
                 return ss.str();
             }()
        );

        LOG_TRACE_L1(m_logger, "Validating candidate groups with flux analysis...");
        const auto [validated_groups, invalidate_groups] = validateGroupsWithFluxAnalysis(candidate_groups, comp, T9, rho);
        LOG_TRACE_L1(
            m_logger,
            "Validated {} group(s) QSE groups. {}",
            validated_groups.size(),
            [&]() -> std::string {
                std::stringstream ss;
                int count = 0;
                for (const auto& group : validated_groups) {
                    ss << "Group " << count + 1;
                    if (group.is_in_equilibrium) {
                        ss << " is in equilibrium";
                    } else {
                        ss << " is not in equilibrium";
                    }
                    if (count < validated_groups.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                return ss.str();
            }()
        );

        // Push the invalidated groups' species into the dynamic set
        for (const auto& group : invalidate_groups) {
            for (const auto& species : group.algebraic_species) {
                m_dynamic_species.push_back(species);
            }
        }

        m_qse_groups = validated_groups;
        LOG_TRACE_L1(m_logger, "Identified {} QSE groups.", m_qse_groups.size());

        for (const auto& group : m_qse_groups) {
            // Add algebraic species to the algebraic set
            for (const auto& species : group.algebraic_species) {
                if (std::ranges::find(m_algebraic_species, species) == m_algebraic_species.end()) {
                    m_algebraic_species.push_back(species);
                }
            }
        }

        LOG_INFO(
            m_logger,
            "Partitioning complete. Found {} dynamic species, {} algebraic (QSE) species ({}) spread over {} QSE group{}.",
            m_dynamic_species.size(),
            m_algebraic_species.size(),
            [&]() -> std::string {
                std::stringstream ss;
                size_t count = 0;
                for (const auto& species : m_algebraic_species) {
                    ss << species.name();
                    if (m_algebraic_species.size() > 1 && count < m_algebraic_species.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                return ss.str();
            }(),
            m_qse_groups.size(),
            m_qse_groups.size() == 1 ? "" : "s"
        );
    }

    void MultiscalePartitioningEngineView::partitionNetwork(
        const NetIn &netIn
    ) {
        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        partitionNetwork(netIn.composition, T9, rho);
    }

    void MultiscalePartitioningEngineView::exportToDot(
        const std::string &filename,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file for writing: {}", filename);
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }

        const auto& all_species = m_baseEngine.getNetworkSpecies();
        const auto& all_reactions = m_baseEngine.getNetworkReactions();

        // --- 1. Pre-computation and Categorization ---

        // Categorize species into algebraic, seed, and core dynamic
        std::unordered_set<Species> algebraic_species;
        std::unordered_set<Species> seed_species;
        for (const auto& group : m_qse_groups) {
            if (group.is_in_equilibrium) {
                algebraic_species.insert(group.algebraic_species.begin(), group.algebraic_species.end());
                seed_species.insert(group.seed_species.begin(), group.seed_species.end());
            }
        }

        // Calculate reaction flows and find min/max for logarithmic scaling of transparency
        std::vector<double> reaction_flows;
        reaction_flows.reserve(all_reactions.size());
        double min_log_flow = std::numeric_limits<double>::max();
        double max_log_flow = std::numeric_limits<double>::lowest();

        for (const auto& reaction : all_reactions) {
            double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(*reaction, comp, T9, rho));
            reaction_flows.push_back(flow);
            if (flow > 1e-99) { // Avoid log(0)
                double log_flow = std::log10(flow);
                min_log_flow = std::min(min_log_flow, log_flow);
                max_log_flow = std::max(max_log_flow, log_flow);
            }
        }
        const double log_flow_range = (max_log_flow > min_log_flow) ? (max_log_flow - min_log_flow) : 1.0;

        // --- 2. Write DOT file content ---

        dotFile << "digraph PartitionedNetwork {\n";
        dotFile << "    graph [rankdir=TB, splines=true, overlap=false, bgcolor=\"#f8fafc\", label=\"Multiscale Partitioned Network View\", fontname=\"Helvetica\", fontsize=16, labeljust=l];\n";
        dotFile << "    node [shape=circle, style=filled, fontname=\"Helvetica\", width=0.8, fixedsize=true];\n";
        dotFile << "    edge [fontname=\"Helvetica\", fontsize=10];\n\n";

        // --- Node Definitions ---
        // Define all species nodes first, so they can be referenced by clusters and ranks later.
        dotFile << "    // --- Species Nodes Definitions ---\n";
        std::map<int, std::vector<std::string>> species_by_mass;
        for (const auto & species : all_species) {
            std::string fillcolor = "#f1f5f9"; // Default: Other/Uninvolved

            // Determine color based on category. A species can be a seed and also in the core dynamic group.
            // The more specific category (algebraic, then seed) takes precedence.
            if (algebraic_species.contains(species)) {
                fillcolor = "#e0f2fe"; // Light Blue: Algebraic (in QSE)
            } else if (seed_species.contains(species)) {
                fillcolor = "#a7f3d0"; // Light Green: Seed (Dynamic, feeds a QSE group)
            } else if (std::ranges::contains(m_dynamic_species, species)) {
                fillcolor = "#dcfce7"; // Pale Green: Core Dynamic
            }
            dotFile << "    \"" << species.name() << "\" [label=\"" << species.name() << "\", fillcolor=\"" << fillcolor << "\"];\n";

            // Group species by mass number for ranked layout.
            // If species.a() returns incorrect values (e.g., 0 for many species), they will be grouped together here.
            species_by_mass[species.a()].emplace_back(species.name());
        }
        dotFile << "\n";

        // --- Layout and Ranking ---
        // Enforce a top-down layout based on mass number.
        dotFile << "    // --- Layout using Ranks ---\n";
        for (const auto &species_list: species_by_mass | std::views::values) {
            dotFile << "    { rank=same; ";
            for (const auto& name : species_list) {
                dotFile << "\"" << name << "\"; ";
            }
            dotFile << "}\n";
        }
        dotFile << "\n";

        // Chain by mass to get top down ordering
        dotFile << "    // --- Chain by Mass ---\n";
        for (const auto& [mass, species_list] : species_by_mass) {
            // Find the next largest mass in the species list
            int minLargestMass = std::numeric_limits<int>::max();
            for (const auto &next_mass: species_by_mass | std::views::keys) {
                if (next_mass > mass && next_mass < minLargestMass) {
                    minLargestMass = next_mass;
                }
            }
            if (minLargestMass != std::numeric_limits<int>::max()) {
                // Connect the current mass to the next largest mass
                dotFile << "    \"" << species_list[0] << "\" -> \"" << species_by_mass[minLargestMass][0] << "\" [style=invis];\n";
            }
        }

        // --- QSE Group Clusters ---
        // Draw a prominent box around the algebraic species of each valid QSE group.
        dotFile << "    // --- QSE Group Clusters ---\n";
        int group_counter = 0;
        for (const auto& group : m_qse_groups) {
            if (!group.is_in_equilibrium || group.algebraic_species.empty()) {
                continue;
            }
            dotFile << "    subgraph cluster_qse_" << group_counter++ << " {\n";
            dotFile << "        label = \"QSE Group " << group_counter << "\";\n";
            dotFile << "        style = \"filled,rounded\";\n";
            dotFile << "        color = \"#38bdf8\";\n"; // A bright, visible blue for the border
            dotFile << "        penwidth = 2.0;\n";      // Thicker border
            dotFile << "        bgcolor = \"#f0f9ff80\";\n"; // Light blue fill with transparency
            dotFile << "        subgraph cluster_seed_" << group_counter << " {\n";
            dotFile << "            label = \"Seed Species\";\n";
            dotFile << "            style = \"filled,rounded\";\n";
            dotFile << "            color = \"#a7f3d0\";\n"; // Light green for seed species
            dotFile << "            penwidth = 1.5;\n"; // Thinner border for seed cluster
            std::vector<std::string> seed_node_ids;
            seed_node_ids.reserve(group.seed_species.size());
            for (const auto& species : group.seed_species) {
                std::stringstream ss;
                ss << "node_" << group_counter << "_seed_" << species.name();
                dotFile << "            " << ss.str() << "    [label=\"" << species.name() << "\"];\n";
                seed_node_ids.push_back(ss.str());
            }
            for (size_t i = 0; i < seed_node_ids.size() - 1; ++i) {
                dotFile << "            " << seed_node_ids[i] << " -> " << seed_node_ids[i + 1] << " [style=invis];\n";
            }
            dotFile << "        }\n";
            dotFile << "        subgraph cluster_algebraic_" << group_counter << " {\n";
            dotFile << "            label = \"Algebraic Species\";\n";
            dotFile << "            style = \"filled,rounded\";\n";
            dotFile << "            color = \"#e0f2fe\";\n"; // Light blue for algebraic species
            dotFile << "            penwidth = 1.5;\n"; // Thinner border for algebraic cluster
            std::vector<std::string> algebraic_node_ids;
            algebraic_node_ids.reserve(group.algebraic_species.size());
            for (const Species& species : group.algebraic_species) {
                std::stringstream ss;
                ss << "node_" << group_counter << "_algebraic_" << species.name();
                dotFile << "            " << ss.str() << "    [label=\"" << species.name() << "\"];\n";
                algebraic_node_ids.push_back(ss.str());
            }
            // Make invisible edges between algebraic indices to keep them in top-down order
            for (size_t i = 0; i < algebraic_node_ids.size() - 1; ++i) {
                dotFile << "            " << algebraic_node_ids[i] << " -> " << algebraic_node_ids[i + 1] << " [style=invis];\n";
            }
            dotFile << "        }\n";
            dotFile << "    }\n";
        }
        dotFile << "\n";


        // --- Legend ---
        // Add a legend to explain colors and conventions.
        dotFile << "    // --- Legend ---\n";
        dotFile << "    subgraph cluster_legend {\n";
        dotFile << "        rank = sink"; // Try to push the legend to the bottom
        dotFile << "        label = \"Legend\";\n";
        dotFile << "        bgcolor = \"#ffffff\";\n";
        dotFile << "        color = \"#e2e8f0\";\n";
        dotFile << "        node [shape=box, style=filled, fontname=\"Helvetica\"];\n";
        dotFile << "        key_core [label=\"Core Dynamic\", fillcolor=\"#dcfce7\"];\n";
        dotFile << "        key_seed [label=\"Seed (Dynamic)\", fillcolor=\"#a7f3d0\"];\n";
        dotFile << "        key_qse [label=\"Algebraic (QSE)\", fillcolor=\"#e0f2fe\"];\n";
        dotFile << "        key_other [label=\"Other\", fillcolor=\"#f1f5f9\"];\n";
        dotFile << "        key_info [label=\"Edge Opacity ~ log(Reaction Flow)\", shape=plaintext];\n";
        dotFile << "        ";// Use invisible edges to stack legend items vertically
        dotFile << "        key_core -> key_seed -> key_qse -> key_other -> key_info [style=invis];\n";
        dotFile << "    }\n\n";

        // --- Reaction Edges ---
        // Draw edges with transparency scaled by the log of the molar reaction flow.
        dotFile << "    // --- Reaction Edges ---\n";
        for (size_t i = 0; i < all_reactions.size(); ++i) {
            const auto& reaction = all_reactions[i];
            const double flow = reaction_flows[i];

            if (flow < 1e-99) continue; // Don't draw edges for negligible flows

            double log_flow_val = std::log10(flow);
            double norm_alpha = (log_flow_val - min_log_flow) / log_flow_range;
            int alpha_val = 0x30 + static_cast<int>(norm_alpha * (0xFF - 0x30)); // Scale from ~20% to 100% opacity
            alpha_val = std::clamp(alpha_val, 0x00, 0xFF);

            std::stringstream alpha_hex;
            alpha_hex << std::setw(2) << std::setfill('0') << std::hex << alpha_val;
            std::string edge_color = "#475569" + alpha_hex.str();

            std::string reactionNodeId = "reaction_" + std::string(reaction.id());
            dotFile << "    \"" << reactionNodeId << "\" [shape=point, fillcolor=black, width=0.05, height=0.05];\n";
            for (const auto& reactant : reaction.reactants()) {
                dotFile << "    \"" << reactant.name() << "\" -> \"" << reactionNodeId << "\" [color=\"" << edge_color << "\", arrowhead=none];\n";
            }
            for (const auto& product : reaction.products()) {
                dotFile << "    \"" << reactionNodeId << "\" -> \"" << product.name() << "\" [color=\"" << edge_color << "\"];\n";
            }
            dotFile << "\n";
        }

        dotFile << "}\n";
        dotFile.close();
    }


    std::vector<double> MultiscalePartitioningEngineView::mapNetInToMolarAbundanceVector(const NetIn &netIn) const {
        std::vector<double> Y(m_dynamic_species.size(), 0.0); // Initialize with zeros
        for (const auto& [symbol, entry] : netIn.composition) {
            Y[getSpeciesIndex(entry.isotope())] = netIn.composition.getMolarAbundance(symbol); // Map species to their molar abundance
        }
        return Y; // Return the vector of molar abundances
    }

    std::vector<Species> MultiscalePartitioningEngineView::getFastSpecies() const {
        const auto& all_species = m_baseEngine.getNetworkSpecies();
        std::vector<Species> fast_species;
        fast_species.reserve(all_species.size() - m_dynamic_species.size());
        for (const auto& species : all_species) {
            auto it = std::ranges::find(m_dynamic_species, species);
            if (it == m_dynamic_species.end()) {
                fast_species.push_back(species);
            }
        }
        return fast_species;
    }

    const std::vector<Species> & MultiscalePartitioningEngineView::getDynamicSpecies() const {
        return m_dynamic_species;
    }

    PrimingReport MultiscalePartitioningEngineView::primeEngine(const NetIn &netIn) {
        return m_baseEngine.primeEngine(netIn);
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::equilibrateNetwork(
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) {
        partitionNetwork(comp, T9, rho);
        fourdst::composition::Composition qseComposition = solveQSEAbundances(comp, T9, rho);

        for (const auto &symbol: qseComposition | std::views::keys) {
            const double speciesMassFraction = qseComposition.getMassFraction(symbol);
            if (speciesMassFraction < 0.0 && std::abs(speciesMassFraction) < 1e-20) {
                qseComposition.setMassFraction(symbol, 0.0); // Avoid negative mass fractions
            }
        }

        bool didFinalize = qseComposition.finalize(true);
        if (!didFinalize) {
            LOG_ERROR(m_logger, "Failed to finalize composition after solving QSE abundances.");
            m_logger->flush_log();
            throw std::runtime_error("Failed to finalize composition after solving QSE abundances.");
        }

        return qseComposition;
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::equilibrateNetwork(
        const NetIn &netIn
    ) {
        const PrimingReport primingReport = m_baseEngine.primeEngine(netIn);

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        return equilibrateNetwork(primingReport.primedComposition, T9, rho);
    }

    size_t MultiscalePartitioningEngineView::getSpeciesIndex(const Species &species) const {
        return m_baseEngine.getSpeciesIndex(species);
    }


    std::vector<std::vector<Species>> MultiscalePartitioningEngineView::partitionByTimescale(
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        LOG_TRACE_L1(m_logger, "Partitioning by timescale...");
        const auto destructionTimescale= m_baseEngine.getSpeciesDestructionTimescales(comp, T9, rho);
        const auto netTimescale = m_baseEngine.getSpeciesTimescales(comp, T9, rho);

        if (!destructionTimescale) {
            LOG_ERROR(m_logger, "Failed to get species destruction timescales due to stale engine state");
            m_logger->flush_log();
            throw exceptions::StaleEngineError("Failed to get species destruction timescales due to stale engine state");
        }
        if (!netTimescale) {
            LOG_ERROR(m_logger, "Failed to get net species timescales due to stale engine state");
            m_logger->flush_log();
            throw exceptions::StaleEngineError("Failed to get net species timescales due to stale engine state");
        }
        const std::unordered_map<Species, double>& destruction_timescales = destructionTimescale.value();
        const std::unordered_map<Species, double>& net_timescales = netTimescale.value();


        for (const auto& [species, destruction_timescale] : destruction_timescales) {
            LOG_TRACE_L3(m_logger, "For {} destruction timescale is {} s", species.name(), destruction_timescale);
        }

        const auto& all_species = m_baseEngine.getNetworkSpecies();

        std::vector<std::pair<double, Species>> sorted_destruction_timescales;
        for (const auto & species : all_species) {
            double destruction_timescale = destruction_timescales.at(species);
            double net_timescale = net_timescales.at(species);
            if (std::isfinite(destruction_timescale) && destruction_timescale > 0) {
                LOG_TRACE_L3(m_logger, "Species {} has finite destruction timescale: destruction: {} s, net: {} s", species.name(), destruction_timescale, net_timescale);
                sorted_destruction_timescales.emplace_back(destruction_timescale, species);
            } else {
                LOG_TRACE_L3(m_logger, "Species {} has infinite or negative destruction timescale: destruction: {} s, net: {} s", species.name(), destruction_timescale, net_timescale);
            }
        }

        std::ranges::sort(
            sorted_destruction_timescales,
            [](const auto& a, const auto& b)
            {
                return a.first > b.first;
            }
        );

        std::vector<std::vector<Species>> final_pools;
        if (sorted_destruction_timescales.empty()) {
            return final_pools;
        }

        constexpr double ABSOLUTE_QSE_TIMESCALE_THRESHOLD = 3.156e7; // Absolute threshold for QSE timescale (1 yr)
        constexpr double MIN_GAP_THRESHOLD = 2.0; // Require a 2 order of magnitude gap
        constexpr double MIN_MOLAR_ABUNDANCE_THRESHOLD = 1e-10; // Minimum abundance threshold to consider a species for QSE. Any species above this will always be considered dynamic.

        LOG_TRACE_L1(m_logger, "Found {} species with finite timescales.", sorted_destruction_timescales.size());
        LOG_TRACE_L1(m_logger, "Absolute QSE timescale threshold: {} seconds ({} years).",
            ABSOLUTE_QSE_TIMESCALE_THRESHOLD, ABSOLUTE_QSE_TIMESCALE_THRESHOLD / 3.156e7);
        LOG_TRACE_L1(m_logger, "Minimum gap threshold: {} orders of magnitude.", MIN_GAP_THRESHOLD);
        LOG_TRACE_L1(m_logger, "Minimum molar abundance threshold: {}.", MIN_MOLAR_ABUNDANCE_THRESHOLD);

        std::vector<Species> dynamic_pool_species;
        std::vector<std::pair<double, Species>> fast_candidates;

        // 1. First Pass: Absolute Timescale Cutoff
        for (const auto& [destruction_timescale, species] : sorted_destruction_timescales) {
            if (species == n_1) {
                LOG_TRACE_L3(m_logger, "Skipping neutron (n) from timescale analysis. Neutrons are always considered dynamic due to their extremely high connectivity.");
                dynamic_pool_species.push_back(species);
                continue;
            }
            if (destruction_timescale > ABSOLUTE_QSE_TIMESCALE_THRESHOLD) {
                LOG_TRACE_L3(m_logger, "Species {} with timescale {} is considered dynamic (slower than qse timescale threshold).",
                    species.name(), destruction_timescale);
                dynamic_pool_species.push_back(species);
            } else {
                const double Yi = comp.getMolarAbundance(species);
                if (Yi > MIN_MOLAR_ABUNDANCE_THRESHOLD) {
                    LOG_TRACE_L3(m_logger, "Species {} with abundance {} is considered dynamic (above minimum abundance threshold of {}).",
                        species.name(), Yi, MIN_MOLAR_ABUNDANCE_THRESHOLD);
                    dynamic_pool_species.push_back(species);
                    continue;
                }
                LOG_TRACE_L3(m_logger, "Species {} with timescale {} and molar abundance {} is a candidate fast species (faster than qse timescale threshold and less than the molar abundance threshold).",
                    species.name(), destruction_timescale, Yi);
                fast_candidates.emplace_back(destruction_timescale, species);
            }
        }

        if (!dynamic_pool_species.empty()) {
            LOG_TRACE_L1(m_logger, "Found {} dynamic species (slower than QSE timescale threshold).", dynamic_pool_species.size());
            final_pools.push_back(dynamic_pool_species);
        }

        if (fast_candidates.empty()) {
            LOG_TRACE_L1(m_logger, "No fast candidates found.");
            return final_pools;
        }

        // 2. Second Pass: Gap Detection on the remaining "fast" candidates
        std::vector<size_t> split_points;
        for (size_t i = 0; i < fast_candidates.size() - 1; ++i) {
            const double t1 = fast_candidates[i].first;
            const double t2 = fast_candidates[i+1].first;
            if (std::log10(t1) - std::log10(t2) > MIN_GAP_THRESHOLD) {
                LOG_TRACE_L3(m_logger, "Detected gap between species {} (timescale {:0.2E}) and {} (timescale {:0.2E}).",
                    fast_candidates[i].second.name(), t1,
                    fast_candidates[i+1].second.name(), t2);
                split_points.push_back(i + 1);
            }
        }

        size_t last_split = 0;
        for (const size_t split : split_points) {
            std::vector<Species> pool;
            for (size_t i = last_split; i < split; ++i) {
                pool.push_back(fast_candidates[i].second);
            }
            final_pools.push_back(pool);
            last_split = split;
        }

        std::vector<Species> final_fast_pool;
        for (size_t i = last_split; i < fast_candidates.size(); ++i) {
            final_fast_pool.push_back(fast_candidates[i].second);
        }
        final_pools.push_back(final_fast_pool);

        LOG_TRACE_L1(m_logger, "Final partitioned pools: {}",
            [&]() -> std::string {
                std::stringstream ss;
                int oc = 0;
                for (const auto& pool : final_pools) {
                    ss << "[";
                    int ic = 0;
                    for (const auto& species : pool) {
                        ss << species.name();
                        if (ic < pool.size() - 1) {
                            ss << ", ";
                        }
                        ic++;
                    }
                    ss << "]";
                    if (oc < final_pools.size() - 1) {
                        ss << ", ";
                    }
                    oc++;
                }
                return ss.str();
            }());
        return final_pools;

    }

    std::pair<std::vector<MultiscalePartitioningEngineView::QSEGroup>, std::vector<MultiscalePartitioningEngineView::
    QSEGroup>>
    MultiscalePartitioningEngineView::validateGroupsWithFluxAnalysis(
        const std::vector<QSEGroup> &candidate_groups,
        const fourdst::composition::Composition &comp,
        const double T9, const double rho
    ) const {
        std::vector<QSEGroup> validated_groups;
        std::vector<QSEGroup> invalidated_groups;
        validated_groups.reserve(candidate_groups.size());
        for (auto& group : candidate_groups) {
            constexpr double FLUX_RATIO_THRESHOLD = 5;
            constexpr double LOG_FLOW_RATIO_THRESHOLD = 2;

            const std::unordered_set<Species> algebraic_group_members(
                group.algebraic_species.begin(),
                group.algebraic_species.end()
            );

            const std::unordered_set<Species> seed_group_members(
                group.seed_species.begin(),
                group.seed_species.end()
            );

            // Values for measuring the flux coupling vs leakage
            double coupling_flux = 0.0;
            double leakage_flux = 0.0;

            // Values for validating if the group could physically be in equilibrium
            double creationFlux = 0.0;
            double destructionFlux = 0.0;

            for (const auto& reaction: m_baseEngine.getNetworkReactions()) {
                const double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(*reaction, comp, T9, rho));
                if (flow == 0.0) {
                    continue; // Skip reactions with zero flow
                }
                bool has_internal_algebraic_reactant = false;

                for (const auto& reactant : reaction->reactants()) {
                    if (algebraic_group_members.contains(reactant)) {
                        has_internal_algebraic_reactant = true;
                        LOG_TRACE_L3(m_logger, "Adjusting destruction flux (+= {} mol g^-1 s^-1) for QSEGroup due to reactant {} from reaction {}",
                            flow, reactant.name(), reaction->id());
                        destructionFlux += flow;
                    }

                }

                bool has_internal_algebraic_product = false;

                for (const auto& product : reaction->products()) {
                    if (algebraic_group_members.contains(product)) {
                        has_internal_algebraic_product = true;
                        LOG_TRACE_L3(m_logger, "Adjusting creation flux (+= {} mol g^-1 s^-1) for QSEGroup due to product {} from reaction {}",
                            flow, product.name(), reaction->id());
                        creationFlux += flow;
                    }
                }

                if (!has_internal_algebraic_product && !has_internal_algebraic_reactant) {
                    LOG_TRACE_L3(m_logger, "{}: Skipping reaction {} as it has no internal algebraic species in reactants or products.", group.toString(), reaction->id());
                    continue;
                }

                double algebraic_participants = 0;
                double seed_participants = 0;
                double external_participants = 0;

                std::unordered_set<Species> participants;
                for(const auto& p : reaction->reactants()) participants.insert(p);
                for(const auto& p : reaction->products()) participants.insert(p);

                for (const auto& species : participants) {
                    if (algebraic_group_members.contains(species)) {
                        LOG_TRACE_L3(m_logger, "{}: Species {} is an algebraic participant in reaction {}.", group.toString(), species.name(), reaction->id());
                        algebraic_participants++;
                    } else if (seed_group_members.contains(species)) {
                        LOG_TRACE_L3(m_logger, "{}: Species {} is a seed participant in reaction {}.", group.toString(), species.name(), reaction->id());
                        seed_participants++;
                    } else {
                        LOG_TRACE_L3(m_logger, "{}: Species {} is an external participant in reaction {}.", group.toString(), species.name(), reaction->id());
                        external_participants++;
                    }
                }

                const double total_participants = algebraic_participants + seed_participants + external_participants;

                if (total_participants == 0) {
                    LOG_CRITICAL(m_logger, "Some catastrophic error has occurred. Reaction {} has no participants.", reaction->id());
                    throw std::runtime_error("Some catastrophic error has occurred. Reaction " + std::string(reaction->id()) + " has no participants.");
                }

                const double leakage_fraction = external_participants / total_participants;
                const double coupling_fraction = (algebraic_participants + seed_participants) / total_participants;

                leakage_flux += flow * leakage_fraction;
                coupling_flux += flow * coupling_fraction;


            }

            bool group_is_coupled = (coupling_flux / leakage_flux) > FLUX_RATIO_THRESHOLD;
            bool group_is_balanced = std::abs(std::log(creationFlux) - std::log(destructionFlux)) < LOG_FLOW_RATIO_THRESHOLD;

            if (group_is_coupled) {
                LOG_TRACE_L1(
                    m_logger,
                    "Group containing {} is in equilibrium due to high coupling flux and balanced creation and destruction rate: <coupling: leakage flux = {}, coupling flux = {}, ratio = {} (Threshold: {})>, <creation: creation flux = {}, destruction flux = {}, ratio = {} order of mag (Threshold: {} order of mag)>",
                    [&]() -> std::string {
                        std::stringstream ss;
                        int count = 0;
                        for (const auto& species: group.algebraic_species) {
                            ss << species.name();
                            if (count < group.algebraic_species.size() - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }(),
                    leakage_flux,
                    coupling_flux,
                    coupling_flux / leakage_flux,
                    FLUX_RATIO_THRESHOLD,
                    std::log10(creationFlux),
                    std::log10(destructionFlux),
                    std::abs(std::log10(creationFlux) - std::log10(destructionFlux)),
                    LOG_FLOW_RATIO_THRESHOLD
                );
                validated_groups.emplace_back(group);
                validated_groups.back().is_in_equilibrium = true;
            } else {
                LOG_TRACE_L1(
                    m_logger,
                    "Group containing {} is NOT in equilibrium: <coupling: leakage flux = {}, coupling flux = {}, ratio = {} (Threshold: {})>, <creation: creation flux = {}, destruction flux = {}, ratio = {} order of mag (Threshold: {} order of mag)>",
                    [&]() -> std::string {
                        std::stringstream ss;
                        int count = 0;
                        for (const auto& species : group.algebraic_species) {
                            ss << species.name();
                            if (count < group.algebraic_species.size() - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }(),
                    leakage_flux,
                    coupling_flux,
                    coupling_flux / leakage_flux,
                    FLUX_RATIO_THRESHOLD,
                    std::log10(creationFlux),
                    std::log10(destructionFlux),
                    std::abs(std::log10(creationFlux) -  std::log10(destructionFlux)),
                    LOG_FLOW_RATIO_THRESHOLD
                );
                invalidated_groups.emplace_back(group);
                invalidated_groups.back().is_in_equilibrium = false;
            }
        }
        return {validated_groups, invalidated_groups};
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::solveQSEAbundances(
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) {
        LOG_TRACE_L1(m_logger, "Solving for QSE abundances...");

        // Sort by timescale to solve fastest QSE groups first (these can seed slower groups)
        std::ranges::sort(m_qse_groups, [](const QSEGroup& a, const QSEGroup& b) {
            return a.mean_timescale < b.mean_timescale;
        });
        fourdst::composition::Composition outputComposition = comp;

        for (const auto&[is_in_equilibrium, algebraic_species, seed_species, mean_timescale] : m_qse_groups) {
            if (!is_in_equilibrium || (algebraic_species.empty() && seed_species.empty())) {
                continue;
            }
            fourdst::composition::Composition normalized_composition = comp;
            for (const auto& species: algebraic_species) {
                if (!normalized_composition.contains(species)) {
                    normalized_composition.registerSpecies(species);
                    normalized_composition.setMassFraction(species, 0.0);
                }
            }
            for (const auto& species: seed_species) {
                if (!normalized_composition.contains(species)) {
                    normalized_composition.registerSpecies(species);
                    normalized_composition.setMassFraction(species, 0.0);
                }
            }
            bool normCompFinalizedOkay = normalized_composition.finalize(true);
            if (!normCompFinalizedOkay) {
                LOG_ERROR(m_logger, "Failed to finalize composition before QSE solve.");
                throw std::runtime_error("Failed to finalize composition before QSE solve.");
            }

            Eigen::VectorXd Y_scale(algebraic_species.size());
            Eigen::VectorXd v_initial(algebraic_species.size());
            long i = 0;
            std::unordered_map<Species, size_t> species_to_index_map;
            for (const auto& species : algebraic_species) {
                constexpr double abundance_floor = 1.0e-15;
                const double initial_abundance = normalized_composition.getMolarAbundance(species);
                Y_scale(i) = std::max(initial_abundance, abundance_floor);
                v_initial(i) = std::asinh(initial_abundance / Y_scale(i)); // Scale the initial abundances using asinh
                species_to_index_map.emplace(species, i);
                i++;
            }

            EigenFunctor functor(*this, algebraic_species, normalized_composition, T9, rho, Y_scale, species_to_index_map);
            Eigen::LevenbergMarquardt lm(functor);
            lm.parameters.ftol = 1.0e-10;
            lm.parameters.xtol = 1.0e-10;

            LOG_TRACE_L1(m_logger, "Minimizing functor...");
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(v_initial);

            if (status <= 0 || status >= 4) {
                std::stringstream msg;
                msg << "While working on QSE group with algebraic species: ";
                int count = 0;
                for (const auto& species: algebraic_species) {
                    msg << species;
                    if (count < algebraic_species.size() - 1) {
                        msg << ", ";
                    }
                    count++;
                }
                msg << " the QSE solver failed to converge with status: " << lm_status_map.at(status);
                msg << ". This likely indicates that the QSE groups were not properly partitioned";
                msg << " (for example if the algebraic set here contains species which should be dynamic like He-4 or H-1).";
                LOG_ERROR(m_logger, "{}", msg.str());
                throw std::runtime_error(msg.str());
            }
            LOG_TRACE_L1(m_logger, "QSE Group minimization succeeded with status: {}", lm_status_map.at(status));
            Eigen::VectorXd Y_final_qse = Y_scale.array() * v_initial.array().sinh(); // Convert back to physical abundances using asinh scaling
            i = 0;
            for (const auto& species: algebraic_species) {
                LOG_TRACE_L1(
                    m_logger,
                    "During QSE solving species {} started with a molar abundance of {} and ended with an abundance of {}.",
                    species.name(),
                    normalized_composition.getMolarAbundance(species),
                    Y_final_qse(i)
                );
                //TODO: Check this conversion
                double Xi = Y_final_qse(i) * species.mass(); // Convert from molar abundance to mass fraction
                if (!outputComposition.hasSpecies(species)) {
                    outputComposition.registerSpecies(species);
                }
                outputComposition.setMassFraction(species, Xi);
                i++;
            }
        }
        bool didFinalize = outputComposition.finalize(false);
        if (!didFinalize) {
            LOG_ERROR(m_logger, "Failed to finalize composition after solving QSE abundances.");
            m_logger->flush_log();
            throw std::runtime_error("Failed to finalize composition after solving QSE abundances.");
        }
        return outputComposition;
    }

    size_t MultiscalePartitioningEngineView::identifyMeanSlowestPool(
        const std::vector<std::vector<Species>> &pools,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        const auto& result = m_baseEngine.getSpeciesDestructionTimescales(comp, T9, rho);
        if (!result) {
            LOG_ERROR(m_logger, "Failed to get species timescales due to stale engine state");
            m_logger->flush_log();
            throw exceptions::StaleEngineError("Failed to get species timescales due to stale engine state");
        }
        const std::unordered_map<Species, double> all_timescales = result.value();

        size_t slowest_pool_index = 0; // Default to the first pool if no valid pool is found
        double slowest_mean_timescale = std::numeric_limits<double>::min();
        size_t count = 0;
        for (const auto& pool : pools) {
            double mean_timescale = 0.0;
            for (const auto& species : pool) {
                const double timescale = all_timescales.at(species);
                mean_timescale += timescale;
            }
            mean_timescale = mean_timescale / static_cast<double>(pool.size());
            if (std::isinf(mean_timescale)) {
                LOG_CRITICAL(m_logger, "Encountered infinite mean timescale for pool {} with species: {}",
                    count, [&]() -> std::string {
                        std::stringstream ss;
                        size_t iCount = 0;
                        for (const auto& species : pool) {
                            ss << species.name() << ": " << all_timescales.at(species);
                            if (iCount < pool.size() - 1) {
                                ss << ", ";
                            }
                            iCount++;
                        }
                        return ss.str();
                    }()
                );
                m_logger->flush_log();
                throw std::logic_error("Encountered infinite mean destruction timescale for a pool while identifying the mean slowest pool set, indicating a potential issue with species timescales. Check log file for more details on specific pool composition...");
            }
            if (mean_timescale > slowest_mean_timescale) {
                slowest_mean_timescale = mean_timescale;
                slowest_pool_index = &pool - &pools[0]; // Get the index of the pool
            }
        }
        return slowest_pool_index;
    }

    std::unordered_map<Species, std::vector<Species>> MultiscalePartitioningEngineView::buildConnectivityGraph(
        const std::vector<Species> &species_pool
    ) const {
        std::unordered_map<Species, std::vector<Species>> connectivity_graph;
        const std::set<Species> pool_set(species_pool.begin(), species_pool.end());

        const std::unordered_set<Species> pool_species = [&]() -> std::unordered_set<Species> {
            std::unordered_set<Species> result;
            for (const auto& species : species_pool) {
                result.insert(species);
            }
            return result;
        }();

        std::map<size_t, std::vector<reaction::LogicalReaclibReaction*>> speciesReactionMap;
        std::vector<const reaction::LogicalReaclibReaction*> candidate_reactions;

        for (const auto& reaction : m_baseEngine.getNetworkReactions()) {
            const std::vector<Species> &reactants = reaction->reactants();
            const std::vector<Species> &products = reaction->products();

            std::unordered_set<Species> reactant_set(reactants.begin(), reactants.end());
            std::unordered_set<Species> product_set(products.begin(), products.end());

            // Only consider reactions where at least one distinct reactant and product are in the pool
            if (has_distinct_reactant_and_product_species(pool_species, reactant_set, product_set)) {
                std::set<Species> involvedSet;
                involvedSet.insert(reactants.begin(), reactants.end());
                involvedSet.insert(products.begin(), products.end());

                std::vector<Species> intersection;
                intersection.reserve(involvedSet.size());

                for (const auto& s : pool_species) { // Find intersection with pool species
                    if (involvedSet.contains(s)) {
                        intersection.push_back(s);
                    }
                }

                // Add clique
                for (const auto& u : intersection) {
                    for (const auto& v : intersection) {
                        if (u != v) { // Avoid self-loops
                            connectivity_graph[u].push_back(v);
                        }
                    }
                }
            }
        }


        return connectivity_graph;
    }

    std::vector<MultiscalePartitioningEngineView::QSEGroup> MultiscalePartitioningEngineView::constructCandidateGroups(
        const std::vector<std::vector<Species>> &candidate_pools,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        const auto& all_reactions = m_baseEngine.getNetworkReactions();
        const auto& result = m_baseEngine.getSpeciesDestructionTimescales(comp, T9, rho);
        if (!result) {
            LOG_ERROR(m_logger, "Failed to get species timescales due to stale engine state");
            m_logger->flush_log();
            throw exceptions::StaleEngineError("Failed to get species timescales due to stale engine state");
        }
        const std::unordered_map<Species, double> destruction_timescales = result.value();

        std::vector<QSEGroup> candidate_groups;
        for (const auto& pool : candidate_pools) {
            if (pool.empty()) continue; // Skip empty pools

            // For each pool first identify all topological bridge connections
            std::vector<std::pair<const reaction::Reaction*, double>> bridge_reactions;
            for (const auto& ash: pool) {
                for (const auto& reaction : all_reactions) {
                    if (reaction->contains(ash)) {
                        // Check to make sure there is at least one reactant that is not in the pool
                        // This lets seed nuclei bring mass into the QSE group.
                        bool has_external_reactant = false;
                        for (const auto& reactant : reaction->reactants()) {
                            if (std::ranges::find(pool, reactant) == pool.end()) {
                                has_external_reactant = true;
                                LOG_TRACE_L3(m_logger, "Found external reactant {} in reaction {} for species {}.", reactant.name(), reaction->id(), ash.name());
                                break; // Found an external reactant, no need to check further
                            }
                        }
                        if (has_external_reactant) {
                            double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(*reaction, comp, T9, rho));
                            LOG_TRACE_L3(m_logger, "Found bridge reaction {} with flow {} for species {}.", reaction->id(), flow, ash.name());
                            bridge_reactions.emplace_back(reaction.get(), flow);
                        }
                    }
                }
            }
            std::ranges::sort(
                bridge_reactions,
                [](const auto& a, const auto& b) {
                   return a.second > b.second; // Sort by flow in descending order
                });

            constexpr double MIN_GAP_THRESHOLD = 1; // Minimum logarithmic molar flow gap threshold for bridge reactions
            std::vector<size_t> split_points;
            for (size_t i = 0; i < bridge_reactions.size() - 1; ++i) {
                const double f1 = bridge_reactions[i].second;
                const double f2 = bridge_reactions[i + 1].second;
                if (std::log10(f1) - std::log10(f2) > MIN_GAP_THRESHOLD) {
                    LOG_TRACE_L3(m_logger, "Detected gap between bridge reactions with flows {} and {}.", f1, f2);
                    split_points.push_back(i + 1);
                }
            }

            if (split_points.empty()) { // If no split points were found, we consider the whole set of bridge reactions as one group.
                split_points.push_back(bridge_reactions.size() - 1);
            }

            std::vector<Species> seed_species;
            for (auto &reaction: bridge_reactions | std::views::keys) {
                for (const auto& fuel : reaction->reactants()) {
                    // Only add the fuel if it is not already in the pool
                    if (std::ranges::find(pool, fuel) == pool.end()) {
                        seed_species.push_back(fuel);
                    }
                }
            }
            std::set<Species> pool_species(pool.begin(), pool.end());
            for (const auto& species : seed_species) {
                pool_species.insert(species);
            }
            const std::set<Species> poolSet(pool.begin(), pool.end());
            const std::set<Species> seedSet(seed_species.begin(), seed_species.end());

            double mean_timescale = 0.0;
            for (const auto& species : poolSet) {
                if (destruction_timescales.contains(species)) {
                    mean_timescale += std::min(destruction_timescales.at(species), species.halfLife()); // Use the minimum of destruction timescale and half-life
                } else {
                    mean_timescale += species.halfLife();
                }
            }
            mean_timescale /= static_cast<double>(poolSet.size());
            QSEGroup qse_group(false, poolSet, seedSet, mean_timescale);
            candidate_groups.push_back(qse_group);
        }
        return candidate_groups;
    }


    int MultiscalePartitioningEngineView::EigenFunctor::operator()(const InputType &v_qse, OutputType &f_qse) const {
        fourdst::composition::Composition comp_trial = m_initial_comp;
        Eigen::VectorXd y_qse = m_Y_scale.array() * v_qse.array().sinh(); // Convert to physical abundances using asinh scaling

        for (const auto& species: m_qse_solve_species) {
            if (!comp_trial.hasSymbol(std::string(species.name()))) {
                comp_trial.registerSpecies(species);
            }
            const double molarAbundance = y_qse[static_cast<long>(m_qse_solve_species_index_map.at(species))];
            double massFraction = molarAbundance * species.mass();
            if (massFraction < 0 && std::abs(massFraction) < 1e-20) { // if there is a larger negative mass fraction, let the composition module throw an error
                massFraction = 0.0; // Avoid setting minuscule negative mass fractions due to numerical noise
            }
            comp_trial.setMassFraction(species, massFraction);
        }

        const bool didFinalize = comp_trial.finalize(false);
        if (!didFinalize) {
            std::string msg = std::format("Failed to finalize composition (comp_trial) in {} at line {}", __FILE__, __LINE__);
            throw std::runtime_error(msg);
        }

        const auto result = m_view->getBaseEngine().calculateRHSAndEnergy(comp_trial, m_T9, m_rho);
        if (!result) {
            throw exceptions::StaleEngineError("Failed to calculate RHS and energy due to stale engine state");
        }
        const auto&[dydt, nuclearEnergyGenerationRate] = result.value();
        f_qse.resize(static_cast<long>(m_qse_solve_species.size()));
        long i = 0;
        // TODO: make sure that just counting up i is a valid approach, this is a possible place an indexing bug may have crept in
        for (const auto& species: m_qse_solve_species) {
            f_qse(i) = dydt.at(species);
            i++;
        }

        return 0; // Success
    }

    int MultiscalePartitioningEngineView::EigenFunctor::df(const InputType &v_qse, JacobianType &J_qse) const {
        fourdst::composition::Composition comp_trial = m_initial_comp;
        Eigen::VectorXd y_qse = m_Y_scale.array() * v_qse.array().sinh(); // Convert to physical abundances using asinh scaling

        for (const auto& species: m_qse_solve_species) {
            if (!comp_trial.hasSymbol(std::string(species.name()))) {
                comp_trial.registerSpecies(species);
            }
            const double molarAbundance = y_qse[static_cast<long>(m_qse_solve_species_index_map.at(species))];
            const double massFraction = molarAbundance * species.mass();
            comp_trial.setMassFraction(species, massFraction);
        }

        const bool didFinalize = comp_trial.finalize(false);
        if (!didFinalize) {
            const std::string msg = std::format("Failed to finalize composition (comp_trial) in {} at line {}", __FILE__, __LINE__);
            throw std::runtime_error(msg);
        }

        m_view->getBaseEngine().generateJacobianMatrix(comp_trial, m_T9, m_rho);

        const long N = static_cast<long>(m_qse_solve_species.size());
        J_qse.resize(N, N);
        long rowID = 0;
        for (const auto& rowSpecies : m_qse_solve_species) {
            long colID = 0;
            for (const auto& colSpecies: m_qse_solve_species) {
                J_qse(rowID, colID) = m_view->getBaseEngine().getJacobianMatrixEntry(
                    rowSpecies,
                    colSpecies
                );
                colID += 1;
                LOG_TRACE_L3(m_view->m_logger, "Jacobian[{}, {}] (d(dY({}))/dY({})) = {}", rowID, colID - 1, rowSpecies.name(), colSpecies.name(), J_qse(rowID, colID - 1));
            }
            rowID += 1;
        }

        // Chain rule for asinh scaling:
        for (long j = 0; j < J_qse.cols(); ++j) {
            const double dY_dv = m_Y_scale(j) * std::cosh(v_qse(j));
            J_qse.col(j) *= dY_dv; // Scale the column by the derivative of the asinh scaling
        }
        return 0; // Success
    }


    QSECacheKey::QSECacheKey(
        const double T9,
        const double rho,
        const std::vector<double> &Y
    ) :
    m_T9(T9),
    m_rho(rho),
    m_Y(Y) {
        m_hash = hash();
    }

    size_t QSECacheKey::hash() const {
        std::size_t seed = 0;

        hash_combine(seed, m_Y.size());

        hash_combine(seed, bin(m_T9, m_cacheConfig.T9_tol));
        hash_combine(seed, bin(m_rho, m_cacheConfig.rho_tol));

        double negThresh = 1e-10; // Threshold for considering a value as negative.
        for (double Yi : m_Y) {
            if (Yi < 0.0 && std::abs(Yi) < negThresh) {
                Yi = 0.0; // Avoid negative abundances
            } else if (Yi < 0.0 && std::abs(Yi) >= negThresh) {
                throw std::invalid_argument("Yi should be positive for valid hashing (expected Yi > 0, received " + std::to_string(Yi) + ")");
            }
            hash_combine(seed, bin(Yi, m_cacheConfig.Yi_tol));
        }

        return seed;

    }

    long QSECacheKey::bin(const double value, const double tol) {
        return static_cast<long>(std::floor(value / tol));
    }

    bool QSECacheKey::operator==(const QSECacheKey &other) const {
        return m_hash == other.m_hash;
    }

    bool MultiscalePartitioningEngineView::QSEGroup::operator==(const QSEGroup &other) const {
       return mean_timescale == other.mean_timescale;
    }

    bool MultiscalePartitioningEngineView::QSEGroup::operator<(const QSEGroup &other) const {
        return mean_timescale < other.mean_timescale;
    }

    bool MultiscalePartitioningEngineView::QSEGroup::operator>(const QSEGroup &other) const {
        return mean_timescale > other.mean_timescale;
    }

    bool MultiscalePartitioningEngineView::QSEGroup::operator!=(const QSEGroup &other) const {
        return !(*this == other);
    }

    std::string MultiscalePartitioningEngineView::QSEGroup::toString() const {
        std::stringstream ss;
        ss << "QSEGroup(Algebraic: [";
        size_t count = 0;
        for (const auto& species  : algebraic_species) {
            ss << species.name();
            if (count < algebraic_species.size() - 1) {
                ss << ", ";
            }
            count++;
        }
        ss << "], Seed: [";
        count = 0;
        for (const auto& species : seed_species) {
            ss << species.name();
            if (count < seed_species.size() - 1) {
                ss << ", ";
            }
            count++;
        }
        ss << "], Mean Timescale: " << mean_timescale << ", Is In Equilibrium: " << (is_in_equilibrium ? "True" : "False") << ")";
        return ss.str();

    }

    void MultiscalePartitioningEngineView::CacheStats::hit(const operators op) {
        if (op == operators::All) {
            throw std::invalid_argument("Cannot use 'ALL' as an operator for a hit");
        }

        m_hit ++;
        m_operatorHits[op]++;
    }
    void MultiscalePartitioningEngineView::CacheStats::miss(const operators op) {
        if (op == operators::All) {
            throw std::invalid_argument("Cannot use 'ALL' as an operator for a miss");
        }

        m_miss ++;
        m_operatorMisses[op]++;
    }

    size_t MultiscalePartitioningEngineView::CacheStats::hits(const operators op) const {
        if (op == operators::All) {
            return m_hit;
        }
        return m_operatorHits.at(op);
    }

    size_t MultiscalePartitioningEngineView::CacheStats::misses(const operators op) const {
        if (op == operators::All) {
            return m_miss;
        }
        return m_operatorMisses.at(op);
    }

}
