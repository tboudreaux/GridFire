#include "gridfire/engine/views/engine_multiscale.h"

#include <numeric>

#include "fourdst/logging/logging.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "quill/LogMacros.h"
#include "quill/Logger.h"

namespace gridfire {
    static int s_operator_parens_called = 0;
    using fourdst::atomic::Species;

    MultiscalePartitioningEngineView::MultiscalePartitioningEngineView(
        GraphEngine& baseEngine
    ) : m_baseEngine(baseEngine) {}

    const std::vector<Species> & MultiscalePartitioningEngineView::getNetworkSpecies() const {
        return m_dynamic_species;
    }

    StepDerivatives<double> MultiscalePartitioningEngineView::calculateRHSAndEnergy(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateRHSAndEnergy(Y, T9, rho);
    }

    void MultiscalePartitioningEngineView::generateJacobianMatrix(
        const std::vector<double> &Y_dynamic,
        const double T9,
        const double rho
    ) {
        std::vector<double> Y_full(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (size_t i = 0; i < m_dynamic_species_indices.size(); ++i) {
            Y_full[m_dynamic_species_indices[i]] = Y_dynamic[i];
        }
        // solveQSEAbundances(Y_full, T9, rho);
        m_baseEngine.generateJacobianMatrix(Y_dynamic, T9, rho);
    }

    double MultiscalePartitioningEngineView::getJacobianMatrixEntry(
        const int i,
        const int j
    ) const {
        return m_baseEngine.getJacobianMatrixEntry(i, j);
    }

    void MultiscalePartitioningEngineView::generateStoichiometryMatrix() {
        m_baseEngine.generateStoichiometryMatrix();
    }

    int MultiscalePartitioningEngineView::getStoichiometryMatrixEntry(
        const int speciesIndex,
        const int reactionIndex
    ) const {
        return m_baseEngine.getStoichiometryMatrixEntry(speciesIndex, reactionIndex);
    }

    double MultiscalePartitioningEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho);
    }

    const reaction::LogicalReactionSet & MultiscalePartitioningEngineView::getNetworkReactions() const {
        return m_baseEngine.getNetworkReactions();
    }

    std::unordered_map<Species, double> MultiscalePartitioningEngineView::getSpeciesTimescales(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.getSpeciesTimescales(Y, T9, rho);
    }

    void MultiscalePartitioningEngineView::update(const NetIn &netIn) {
        return; // call partition eventually
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

    void MultiscalePartitioningEngineView::partitionNetwork(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) {
        // --- Step 0. Clear previous state ---
        LOG_TRACE_L1(m_logger, "Partitioning network...");
        LOG_TRACE_L1(m_logger, "Clearing previous state...");
        m_qse_groups.clear();
        m_dynamic_species.clear();
        m_dynamic_species_indices.clear();
        m_connectivity_graph.clear();
        m_algebraic_species.clear();
        m_algebraic_species_indices.clear();

        // --- Step 1. Identify distinct timescale regions ---
        LOG_TRACE_L1(m_logger, "Identifying fast reactions...");
        const std::vector<std::vector<size_t>> timescale_pools = partitionByTimescale(Y, T9, rho);
        LOG_TRACE_L1(m_logger, "Found {} timescale pools.", timescale_pools.size());

        // --- Step 2. Select the mean slowest pool as the base dynamical group ---
        LOG_TRACE_L1(m_logger, "Identifying mean slowest pool...");
        const size_t mean_slowest_pool_index = identifyMeanSlowestPool(timescale_pools, Y, T9, rho);
        LOG_TRACE_L1(m_logger, "Mean slowest pool index: {}", mean_slowest_pool_index);

        // --- Step 3. Push the slowest pool into the dynamic species list ---
        m_dynamic_species_indices = timescale_pools[mean_slowest_pool_index];
        for (const auto& index : m_dynamic_species_indices) {
            m_dynamic_species.push_back(m_baseEngine.getNetworkSpecies()[index]);
        }

        // --- Step 4. Pack Candidate QSE Groups ---
        std::vector<std::vector<size_t>> candidate_pools;
        for (size_t i = 0; i < timescale_pools.size(); ++i) {
            if (i == mean_slowest_pool_index) continue; // Skip the slowest pool
            LOG_TRACE_L1(m_logger, "Group {} with {} species identified for potential QSE.", i, timescale_pools[i].size());
            candidate_pools.push_back(timescale_pools[i]);
        }

        // --- Step 5. Identify potential seed species for each candidate pool ---
        LOG_TRACE_L1(m_logger, "Identifying potential seed species for candidate pools...");
        const std::vector<QSEGroup> candidate_groups = constructCandidateGroups(candidate_pools, Y, T9, rho);
        LOG_TRACE_L1(
            m_logger,
            "Found {} candidate QSE groups for further analysis. {}",
            candidate_groups.size(),
            [&]() -> std::string {
                std::stringstream ss;
                int count = 0;
                for (const auto& group : candidate_groups) {
                    ss << "CandidateGroup(seeds: [";
                    int icount = 0;
                    for (const auto& seed_index : group.seed_indices) {
                        ss << m_baseEngine.getNetworkSpecies()[seed_index].name();
                        if (icount < group.seed_indices.size() - 1) {
                            ss << ", ";
                        }
                        icount++;
                    }
                    ss << "], qse species: [";
                    icount = 0;
                    for (const auto& qse_index : group.algebraic_indices) {
                        ss << m_baseEngine.getNetworkSpecies()[qse_index].name();
                        if (icount < group.algebraic_indices.size() - 1) {
                            ss << ", ";
                        }
                        icount++;
                    }
                    ss << "])";
                    if (count < candidate_groups.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                return ss.str();
            }()
        );

        LOG_TRACE_L1(m_logger, "Validating candidate groups with flux analysis...");
        const std::vector<QSEGroup> validated_groups = validateGroupsWithFluxAnalysis(candidate_groups, Y, T9, rho);
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

        m_qse_groups = std::move(validated_groups);
        LOG_TRACE_L1(m_logger, "Identified {} QSE groups.", m_qse_groups.size());

        for (const auto& group : m_qse_groups) {
            // Add algebraic species to the algebraic set
            for (const auto& index : group.algebraic_indices) {
                if (std::ranges::find(m_algebraic_species_indices, index) == m_algebraic_species_indices.end()) {
                    m_algebraic_species.push_back(m_baseEngine.getNetworkSpecies()[index]);
                    m_algebraic_species_indices.push_back(index);

                }
            }
        }

    }

    void MultiscalePartitioningEngineView::partitionNetwork(
        const NetIn &netIn
    ) {
        std::vector<double> Y(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (const auto& [symbol, entry]: netIn.composition ) {
            if (!m_baseEngine.involvesSpecies(entry.isotope())) {
                LOG_WARNING(m_logger, "Species {} is not part of the network. This is likely due to the base network not being deep enough. For rare species this may not be an issue. Skipping...", entry.isotope().name());
                continue; // Skip species not in the network
            }
            const size_t species_index = m_baseEngine.getSpeciesIndex(entry.isotope());
            Y[species_index] = netIn.composition.getMolarAbundance(symbol);
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        partitionNetwork(Y, T9, rho);
    }

    void MultiscalePartitioningEngineView::exportToDot(
        const std::string &filename,
        const std::vector<double>& Y,
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
        std::unordered_set<size_t> algebraic_indices;
        std::unordered_set<size_t> seed_indices;
        for (const auto& group : m_qse_groups) {
            if (group.is_in_equilibrium) {
                algebraic_indices.insert(group.algebraic_indices.begin(), group.algebraic_indices.end());
                seed_indices.insert(group.seed_indices.begin(), group.seed_indices.end());
            }
        }

        // Calculate reaction flows and find min/max for logarithmic scaling of transparency
        std::vector<double> reaction_flows;
        reaction_flows.reserve(all_reactions.size());
        double min_log_flow = std::numeric_limits<double>::max();
        double max_log_flow = std::numeric_limits<double>::lowest();

        for (const auto& reaction : all_reactions) {
            double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho));
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
        for (size_t i = 0; i < all_species.size(); ++i) {
            const auto& species = all_species[i];
            std::string fillcolor = "#f1f5f9"; // Default: Other/Uninvolved

            // Determine color based on category. A species can be a seed and also in the core dynamic group.
            // The more specific category (algebraic, then seed) takes precedence.
            if (algebraic_indices.count(i)) {
                fillcolor = "#e0f2fe"; // Light Blue: Algebraic (in QSE)
            } else if (seed_indices.count(i)) {
                fillcolor = "#a7f3d0"; // Light Green: Seed (Dynamic, feeds a QSE group)
            } else if (std::ranges::contains(m_dynamic_species_indices, i)) {
                fillcolor = "#dcfce7"; // Pale Green: Core Dynamic
            }
            dotFile << "    \"" << species.name() << "\" [label=\"" << species.name() << "\", fillcolor=\"" << fillcolor << "\"];\n";

            // Group species by mass number for ranked layout.
            // If species.a() returns incorrect values (e.g., 0 for many species), they will be grouped together here.
            species_by_mass[species.a()].push_back(std::string(species.name()));
        }
        dotFile << "\n";

        // --- Layout and Ranking ---
        // Enforce a top-down layout based on mass number.
        dotFile << "    // --- Layout using Ranks ---\n";
        for (const auto& [mass, species_list] : species_by_mass) {
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
            for (const auto& [next_mass, _] : species_by_mass) {
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
            if (!group.is_in_equilibrium || group.algebraic_indices.empty()) {
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
            seed_node_ids.reserve(group.seed_indices.size());
            for (const size_t species_idx : group.seed_indices) {
                std::stringstream ss;
                ss << "node_" << group_counter << "_seed_" << species_idx;
                dotFile << "            " << ss.str() << "    [label=\"" << all_species[species_idx].name() << "\"];\n";
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
            algebraic_node_ids.reserve(group.algebraic_indices.size());
            for (const size_t species_idx : group.algebraic_indices) {
                std::stringstream ss;
                ss << "node_" << group_counter << "_algebraic_" << species_idx;
                dotFile << "            " << ss.str() << "    [label=\"" << all_species[species_idx].name() << "\"];\n";
                algebraic_node_ids.push_back(ss.str());
            }
            // Make invisible edges between algebraic indicies to keep them in top down order
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
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) {
        partitionNetwork(Y, T9, rho);
        const std::vector<double> Y_equilibrated = solveQSEAbundances(Y, T9, rho);
        fourdst::composition::Composition composition;

        std::vector<std::string> symbols;
        symbols.reserve(m_baseEngine.getNetworkSpecies().size());
        for (const auto& species : m_baseEngine.getNetworkSpecies()) {
            symbols.emplace_back(species.name());
        }
        composition.registerSymbol(symbols);

        std::vector<double> X;
        X.reserve(Y_equilibrated.size());
        for (size_t i = 0; i < Y_equilibrated.size(); ++i) {
            const double molarMass = m_baseEngine.getNetworkSpecies()[i].mass();
            X.push_back(Y_equilibrated[i] * molarMass); // Convert from molar abundance to mass fraction
        }

        for (size_t i = 0; i < Y_equilibrated.size(); ++i) {
            const auto& species = m_baseEngine.getNetworkSpecies()[i];
            if (X[i] < 0.0 && std::abs(X[i]) < 1e-20) {
                composition.setMassFraction(std::string(species.name()), 0.0); // Avoid negative mass fractions
            } else {
                composition.setMassFraction(std::string(species.name()), X[i]);
            }
        }

        composition.finalize(true);
        return composition;
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::equilibrateNetwork(
        const NetIn &netIn
    ) {
        std::vector<double> Y(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (const auto& [symbol, entry]: netIn.composition ) {
            if (!m_baseEngine.involvesSpecies(entry.isotope())) {
                LOG_WARNING(m_logger, "Species {} is not part of the network. This is likely due to the base network not being deep enough. For rare species this may not be an issue. Skipping...", entry.isotope().name());
                continue; // Skip species not in the network
            }
            const size_t species_index = m_baseEngine.getSpeciesIndex(entry.isotope());
            Y[species_index] = netIn.composition.getMolarAbundance(symbol);
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        return equilibrateNetwork(Y, T9, rho);
    }

    int MultiscalePartitioningEngineView::getSpeciesIndex(const fourdst::atomic::Species &species) const {
        return m_baseEngine.getSpeciesIndex(species);
    }

    std::vector<std::vector<size_t>> MultiscalePartitioningEngineView::partitionByTimescale(
        const std::vector<double>& Y_full,
        const double T9,
        const double rho
    ) const {
        LOG_TRACE_L1(m_logger, "Partitioning by timescale...");
        const auto all_timescales = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);
        const auto& all_species = m_baseEngine.getNetworkSpecies();

        std::vector<std::pair<double, size_t>> sorted_timescales;
        for (size_t i = 0; i < all_species.size(); ++i) {
            double timescale = all_timescales.at(all_species[i]);
            if (std::isfinite(timescale) && timescale > 0) {
                sorted_timescales.push_back({timescale, i});
            }
        }

        std::ranges::sort(
            sorted_timescales,
            [](const auto& a, const auto& b)
            {
                return a.first > b.first;
            }
        );

        std::vector<std::vector<size_t>> final_pools;
        if (sorted_timescales.empty()) {
            return final_pools;
        }

        constexpr double ABSOLUTE_QSE_TIMESCALE_THRESHOLD = 3.156e7; // Absolute threshold for QSE timescale (1 yr)
        constexpr double MIN_GAP_THRESHOLD = 3.0; // Require a 3 order of magnitude gap

        LOG_TRACE_L1(m_logger, "Found {} species with finite timescales.", sorted_timescales.size());
        LOG_TRACE_L1(m_logger, "Absolute QSE timescale threshold: {} seconds ({} years).",
            ABSOLUTE_QSE_TIMESCALE_THRESHOLD, ABSOLUTE_QSE_TIMESCALE_THRESHOLD / 3.156e7);
        LOG_TRACE_L1(m_logger, "Minimum gap threshold: {} orders of magnitude.", MIN_GAP_THRESHOLD);

        std::vector<size_t> dynamic_pool_indices;
        std::vector<std::pair<double, size_t>> fast_candidates;

        // 1. First Pass: Absolute Timescale Cutoff
        for (const auto& ts_pair : sorted_timescales) {
            if (ts_pair.first > ABSOLUTE_QSE_TIMESCALE_THRESHOLD) {
                LOG_TRACE_L1(m_logger, "Species {} with timescale {} is considered dynamic (slower than qse timescale threshold).",
                    all_species[ts_pair.second].name(), ts_pair.first);
                dynamic_pool_indices.push_back(ts_pair.second);
            } else {
                LOG_TRACE_L1(m_logger, "Species {} with timescale {} is a candidate fast species (faster than qse timescale threshold).",
                    all_species[ts_pair.second].name(), ts_pair.first);
                fast_candidates.push_back(ts_pair);
            }
        }

        if (!dynamic_pool_indices.empty()) {
            LOG_TRACE_L1(m_logger, "Found {} dynamic species (slower than QSE timescale threshold).", dynamic_pool_indices.size());
            final_pools.push_back(dynamic_pool_indices);
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
                LOG_TRACE_L1(m_logger, "Detected gap between species {} (timescale {:0.2E}) and {} (timescale {:0.2E}).",
                    all_species[fast_candidates[i].second].name(), t1,
                    all_species[fast_candidates[i+1].second].name(), t2);
                split_points.push_back(i + 1);
            }
        }

        size_t last_split = 0;
        for (const size_t split : split_points) {
            std::vector<size_t> pool;
            for (size_t i = last_split; i < split; ++i) {
                pool.push_back(fast_candidates[i].second);
            }
            final_pools.push_back(pool);
            last_split = split;
        }

        std::vector<size_t> final_fast_pool;
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
                    for (const auto& idx : pool) {
                        ss << all_species[idx].name();
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

    std::unordered_map<size_t, std::vector<size_t>> MultiscalePartitioningEngineView::buildConnectivityGraph(
        const std::unordered_set<size_t> &fast_reaction_indices
    ) const {
        const auto& all_reactions = m_baseEngine.getNetworkReactions();
        std::unordered_map<size_t, std::vector<size_t>> connectivity;
        for (const size_t reaction_idx : fast_reaction_indices) {
            const auto& reaction = all_reactions[reaction_idx];
            const auto& reactants = reaction.reactants();
            const auto& products = reaction.products();

            // For each fast reaction, create edges between all reactants and all products.
            // This represents that nucleons can flow quickly between these species.
            for (const auto& reactant : reactants) {
                const size_t reactant_idx = m_baseEngine.getSpeciesIndex(reactant);
                for (const auto& product : products) {
                    const size_t product_idx = m_baseEngine.getSpeciesIndex(product);

                    // Add a two-way edge to the adjacency list.
                    connectivity[reactant_idx].push_back(product_idx);
                    connectivity[product_idx].push_back(reactant_idx);
                }
            }
        }
        return connectivity;
    }

    // std::vector<MultiscalePartitioningEngineView::QSEGroup> MultiscalePartitioningEngineView::refineCandidateGroup(
    //     const std::vector<size_t> &candidate_group,
    //     const std::vector<double> &Y,
    //     const double T9,
    //     const double rho
    // ) const {
    //     LOG_TRACE_L1(m_logger, "Refining candidate group with {} species.", candidate_group.size());
    //
    //     // --- Step 1. Graph Construction ---
    //     std::unordered_map<size_t, std::vector<std::pair<size_t, double>>> weighted_graph;
    //     std::vector<std::tuple<size_t, size_t, double>> edge_list;
    //     std::unordered_map<size_t, double> node_strength;
    //
    //     for (const auto& node_idx : candidate_group) {
    //         node_strength[node_idx] = 0.0;
    //     }
    //
    //     std::unordered_set<size_t> candidate_set(candidate_group.begin(), candidate_group.end());
    //
    //     for (const auto& reaction : m_baseEngine.getNetworkReactions()) {
    //         bool is_internal = true;
    //         std::vector<size_t> involved_indices;
    //
    //         // Check if all reactants and products are in the candidate group.
    //         std::vector<Species> participating_species = reaction.reactants();
    //         participating_species.insert(participating_species.end(), reaction.products().begin(), reaction.products().end());
    //
    //         for (const auto& sp : participating_species) {
    //             size_t idx = m_baseEngine.getSpeciesIndex(sp);
    //             if (! candidate_set.contains(idx)) {
    //                 is_internal = false;
    //                 break;
    //             }
    //         }
    //         if (!is_internal) {
    //             continue; // Skip reactions that are not fully internal to the candidate group.
    //         }
    //
    //         const double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho));
    //         if (flow < 1e-99) continue;
    //
    //         for (const auto& reactant : reaction.reactants()) {
    //             for (const auto& product : reaction.products()) {
    //                 size_t u = m_baseEngine.getSpeciesIndex(reactant);
    //                 size_t v = m_baseEngine.getSpeciesIndex(product);
    //
    //                 if (u == v) {
    //                     continue; // Skip self-loops
    //                 }
    //
    //                 weighted_graph[u].push_back(std::make_tuple(v, flow));
    //                 weighted_graph[v].push_back(std::make_tuple(u, flow));
    //
    //                 node_strength[u] += flow;
    //                 node_strength[v] += flow;
    //
    //                 if (u < v) {
    //                     edge_list.emplace_back(u, v, flow);
    //                 }
    //             }
    //         }
    //     }
    //
    //     // --- Phase 2. Prune weak bridges ---
    //     auto pruned_graph = weighted_graph;
    //     constexpr double beta = 0.01; // Weak link threshold: flux < 1% of the node strength
    //
    //     for (const auto& edge : edge_list) {
    //         size_t u = std::get<0>(edge);
    //         size_t v = std::get<1>(edge);
    //         double flow = std::get<2>(edge);
    //
    //         if (flow < (beta * node_strength[u]) && flow < (beta * node_strength[v])) {
    //             LOG_TRACE_L1(m_logger, "Pruning weak link ({}, {}) with flow {} below threshold {}",
    //                          u, v, flow, beta * std::min(node_strength[u], node_strength[v]));
    //
    //             auto& u_adj = pruned_graph.at(u);
    //             std::erase_if(
    //                 u_adj,
    //                 [v](const auto& p) {
    //                     return p.first == v;
    //                 }
    //             );
    //             auto& v_adj = pruned_graph.at(v);
    //             std::erase_if(
    //                 v_adj,
    //                 [u](const auto& p) {
    //                     return p.first == u;
    //                 }
    //             );
    //         }
    //     }
    //
    //     // --- Phase 3. Find final connected components ---
    //     auto final_components = findConnectedComponentsBFS(pruned_graph, candidate_group);
    //
    //     std::vector<QSEGroup> refined_groups;
    //     for (const auto& component : final_components) {
    //         if (!component.empty()) {
    //             std::set<size_t> tempSet(component.begin(), component.end());
    //             refined_groups.push_back({tempSet, false}); // Initialize with is_in_equilibrium = false
    //         }
    //     }
    //
    //     return refined_groups;
    // }


    std::vector<MultiscalePartitioningEngineView::QSEGroup>
    MultiscalePartitioningEngineView::validateGroupsWithFluxAnalysis(
        const std::vector<QSEGroup> &candidate_groups,
        const std::vector<double> &Y,
        const double T9, const double rho
    ) const {
        constexpr double FLUX_RATIO_THRESHOLD = 100;
        std::vector<QSEGroup> validated_groups = candidate_groups;
        for (auto& group : validated_groups) {
            double internal_flux = 0.0;
            double external_flux = 0.0;

            const std::unordered_set<size_t> group_members(
                group.species_indices.begin(),
                group.species_indices.end()
            );

            for (const auto& reaction: m_baseEngine.getNetworkReactions()) {
                const double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho));
                if (flow == 0.0) {
                    continue; // Skip reactions with zero flow
                }
                bool has_internal_reactant = false;
                bool has_external_reactant = false;

                for (const auto& reactant : reaction.reactants()) {
                    if (group_members.contains(m_baseEngine.getSpeciesIndex(reactant))) {
                        has_internal_reactant = true;
                    } else {
                        has_external_reactant = true;
                    }
                }

                bool has_internal_product = false;
                bool has_external_product = false;

                for (const auto& product : reaction.products()) {
                    if (group_members.contains(m_baseEngine.getSpeciesIndex(product))) {
                        has_internal_product = true;
                    } else {
                        has_external_product = true;
                    }
                }

                // Classify the reaction based on its participants
                if ((has_internal_reactant && has_internal_product) && !(has_external_reactant || has_external_product)) {
                    internal_flux += flow; // Internal flux within the group
                } else if ((has_internal_reactant || has_internal_product) && (has_external_reactant || has_external_product)) {
                    external_flux += flow; // External flux to/from the group
                }
                // otherwise the reaction is fully decoupled from the QSE group and can be ignored.
            }
            if (internal_flux > FLUX_RATIO_THRESHOLD * external_flux) {
                group.is_in_equilibrium = true; // This group is in equilibrium if internal flux is significantly larger than external flux.
            } else {
                group.is_in_equilibrium = false;
            }
        }
        return validated_groups;
    }

    std::vector<double> MultiscalePartitioningEngineView::solveQSEAbundances(
        const std::vector<double> &Y_full,
        const double T9,
        const double rho
    ) {
        LOG_TRACE_L1(m_logger, "Solving for QSE abundances...");
        auto Y_out = Y_full;
        auto species_timescales = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);

        for (const auto& qse_group : m_qse_groups) {
            if (!qse_group.is_in_equilibrium || qse_group.species_indices.empty()) {
                LOG_TRACE_L1(
                    m_logger,
                    "Skipping QSE group with {} species ({} algebraic, {} seeds) as it is not in equilibrium.",
                    qse_group.species_indices.size(),
                    qse_group.algebraic_indices.size(),
                    qse_group.seed_indices.size()
                );
                continue;
            }

            LOG_TRACE_L1(
                m_logger,
                "Solving for QSE abundances for group with {} species ([{}] algebraic, [{}] seeds).",
                qse_group.species_indices.size(),
                [&]() -> std::string {
                    std::stringstream ss;
                    int count = 0;
                    for (const auto& idx : qse_group.algebraic_indices) {
                        ss << m_baseEngine.getNetworkSpecies()[idx].name();
                        if (count < qse_group.algebraic_indices.size() - 1) {
                            ss << ", ";
                        }
                        count++;
                    }
                    return ss.str();
                }(),
                [&]() -> std::string {
                    std::stringstream ss;
                    int count = 0;
                    for (const auto& idx : qse_group.seed_indices) {
                        ss << m_baseEngine.getNetworkSpecies()[idx].name();
                        if (count < qse_group.seed_indices.size() - 1) {
                            ss << ", ";
                        }
                        count++;
                    }
                    return ss.str();
                }()
            );

            std::vector<size_t> qse_solve_indices;
            std::vector<size_t> seed_indices;

            seed_indices.reserve(qse_group.species_indices.size());
            qse_solve_indices.reserve(qse_group.species_indices.size());

            for (size_t seed_idx : qse_group.seed_indices) {
                seed_indices.emplace_back(seed_idx);
            }

            for (size_t algebraic_idx : qse_group.algebraic_indices) {
                qse_solve_indices.emplace_back(algebraic_idx);
            }

            if (qse_solve_indices.empty()) continue;

            Eigen::VectorXd Y_scale(qse_solve_indices.size());
            Eigen::VectorXd v_initial(qse_solve_indices.size());
            for (size_t i = 0; i < qse_solve_indices.size(); ++i) {
                constexpr double abundance_floor = 1.0e-15;
                const double initial_abundance = Y_full[qse_solve_indices[i]];
                Y_scale(i) = std::max(initial_abundance, abundance_floor);
                v_initial(i) = std::asinh(initial_abundance / Y_scale(i)); // Scale the initial abundances using asinh
            }

            EigenFunctor functor(*this, qse_solve_indices, Y_full, T9, rho, Y_scale);
            Eigen::LevenbergMarquardt lm(functor);
            lm.parameters.ftol = 1.0e-10;
            lm.parameters.xtol = 1.0e-10;

            LOG_TRACE_L1(m_logger, "Minimizing functor...");
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(v_initial);

            if (status <= 0 || status >= 4) {
                std::stringstream msg;
                msg << "QSE solver failed with status: " << status << " for group with seed ";
                if (seed_indices.size() == 1) {
                    msg << "nucleus " << m_baseEngine.getNetworkSpecies()[seed_indices[0]].name();
                } else {
                    msg << "nuclei ";
                    // TODO: Refactor nice list printing into utility function somewhere
                    int counter = 0;
                    for (const auto& seed_idx : seed_indices) {
                        msg << m_baseEngine.getNetworkSpecies()[seed_idx].name();
                        if (counter < seed_indices.size() - 2) {
                            msg << ", ";
                        } else if (counter == seed_indices.size() - 2) {
                            if (seed_indices.size() < 2) {
                                msg << " and ";
                            } else {
                                msg << ", and ";
                            }
                        }
                        ++counter;
                    }
                }
                throw std::runtime_error(msg.str());
            }
            LOG_TRACE_L1(m_logger, "Minimization succeeded!");
            Eigen::VectorXd Y_final_qse = Y_scale.array() * v_initial.array().sinh(); // Convert back to physical abundances using asinh scaling
            for (size_t i = 0; i < qse_solve_indices.size(); ++i) {
                LOG_TRACE_L1(
                    m_logger,
                    "Species {} (index {}) started with abundance {} and ended with {}.",
                    m_baseEngine.getNetworkSpecies()[qse_solve_indices[i]].name(),
                    qse_solve_indices[i],
                    Y_full[qse_solve_indices[i]],
                    Y_final_qse(i)
                );
                Y_out[qse_solve_indices[i]] = Y_final_qse(i);
            }
        }
        return Y_out;
    }

    size_t MultiscalePartitioningEngineView::identifyMeanSlowestPool(
        const std::vector<std::vector<size_t>> &pools,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        const auto& all_timescales = m_baseEngine.getSpeciesTimescales(Y, T9, rho);
        const auto& all_species = m_baseEngine.getNetworkSpecies();
        size_t slowest_pool_index = 0; // Default to the first pool if no valid pool is found
        double slowest_mean_timescale = std::numeric_limits<double>::min();
        for (const auto& pool : pools) {
            double mean_timescale = 0.0;
            for (const auto& species_idx : pool) {
                double timescale = all_timescales.at(all_species[species_idx]);
                mean_timescale += timescale;
            }
            mean_timescale = mean_timescale / pool.size();
            if (mean_timescale > slowest_mean_timescale) {
                slowest_mean_timescale = mean_timescale;
                slowest_pool_index = &pool - &pools[0]; // Get the index of the pool
            }
        }
        return slowest_pool_index;
    }

    std::vector<MultiscalePartitioningEngineView::QSEGroup> MultiscalePartitioningEngineView::constructCandidateGroups(
        const std::vector<std::vector<size_t>> &timescale_pools,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        const auto& all_species = m_baseEngine.getNetworkSpecies();
        const auto& all_reactions = m_baseEngine.getNetworkReactions();

        std::vector<QSEGroup> candidate_groups;
        for (const auto& pool : timescale_pools) {
            if (pool.empty()) continue; // Skip empty pools
            // For each pool first identify all topological bridge connections
            std::vector<std::pair<reaction::LogicalReaction, double>> bridge_reactions;
            for (const auto& species_idx : pool) {
                Species ash = all_species[species_idx];
                for (const auto& reaction : all_reactions) {
                    if (reaction.contains(ash)) {
                        // Check to make sure there is at least one reactant that is not in the pool
                        // This lets seed nuclei bring mass into the QSE group.
                        bool has_external_reactant = false;
                        for (const auto& reactant : reaction.reactants()) {
                            if (std::ranges::find(pool, m_baseEngine.getSpeciesIndex(reactant)) == pool.end()) {
                                has_external_reactant = true;
                                LOG_TRACE_L2(m_logger, "Found external reactant {} in reaction {} for species {}.",
                                             reactant.name(), reaction.id(), ash.name());
                                break; // Found an external reactant, no need to check further
                            }
                        }
                        if (has_external_reactant) {
                            double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho));
                            LOG_TRACE_L2(m_logger, "Found bridge reaction {} with flow {} for species {}.",
                                         reaction.id(), flow, ash.name());
                            bridge_reactions.push_back({reaction, flow});
                        }
                    }
                }
            }
            std::ranges::sort(
                bridge_reactions,
                [](const auto& a, const auto& b) {
                   return a.second > b.second; // Sort by flow in descending order
                });
            LOG_TRACE_L1(
                m_logger,
                "Found a total of {} bridge reactions for the pool with qse species: {}.",
                bridge_reactions.size(),
                [&]() -> std::string {
                    std::stringstream ss;
                    int count = 0;
                    for (const auto& idx : pool) {
                        ss << all_species[idx].name();
                        if (count < pool.size() - 1) {
                            ss << ", ";
                        }
                        count++;
                    }
                    return ss.str();
                }()
            );

            constexpr double MIN_GAP_THRESHOLD = 1; // Minimum logarithmic molar flow gap threshold for bridge reactions
            std::vector<size_t> split_points;
            for (size_t i = 0; i < bridge_reactions.size() - 1; ++i) {
                const double f1 = bridge_reactions[i].second;
                const double f2 = bridge_reactions[i + 1].second;
                if (std::log10(f1) - std::log10(f2) > MIN_GAP_THRESHOLD) {
                    LOG_TRACE_L2(m_logger, "Detected gap between bridge reactions with flows {} and {}.",
                                 f1, f2);
                    split_points.push_back(i + 1);
                }
            }
            LOG_TRACE_L1(m_logger, "Found {} candidate molar flow split points in bridge reactions.", split_points.size() == 0 ? 1 : split_points.size());
            LOG_TRACE_L1(m_logger, "Split point summary: {}",
                [&]() -> std::string {
                    std::stringstream ss;
                    int count = 0;
                    // Print out which reactions are in each molar reaction flow pool
                    int lastSplitPoint = 0;
                    for (const auto& split_point : split_points) {
                        ss << "Reactions in molar flow pool " << count + 1 << ": [";
                        int icount = 0;
                        for (size_t i = lastSplitPoint; i < split_point; ++i) {
                            ss << bridge_reactions[i].first.id();
                            if (icount < (split_point - lastSplitPoint) - 1) {
                                ss << ", ";
                            }
                            icount++;
                        }
                        ss << "]";
                        if (count < split_points.size() - 1) {
                            ss << ", ";
                        }
                        count++;
                        lastSplitPoint = split_point;
                    }
                    return ss.str();
                }()
            );

            if (split_points.empty()) { // If no split points were found, we consider the whole set of bridge reactions as one group.
                split_points.push_back(bridge_reactions.size() - 1);
            }

            std::vector<size_t> seed_indices;
            for (size_t i = 0; i < bridge_reactions.size(); ++i) {
                for (const auto& fuel : bridge_reactions[i].first.reactants()) {
                    size_t fuel_idx = m_baseEngine.getSpeciesIndex(fuel);
                    // Only add the fuel if it is not already in the pool
                    if (std::ranges::find(pool, fuel_idx) == pool.end()) {
                        seed_indices.push_back(fuel_idx);
                    }
                }
            }
            std::set<size_t> all_indices(pool.begin(), pool.end());
            for (const auto& seed_idx : seed_indices) {
                all_indices.insert(seed_idx);
            }
            const std::set<size_t> poolSet(pool.begin(), pool.end());
            const std::set<size_t> seedSet(seed_indices.begin(), seed_indices.end());
            QSEGroup qse_group(all_indices, false, poolSet, seedSet);
            candidate_groups.push_back(qse_group);
        }
        return candidate_groups;
    }

    int MultiscalePartitioningEngineView::EigenFunctor::operator()(const InputType &v_qse, OutputType &f_qse) const {
        s_operator_parens_called++;
        std::vector<double> y_trial = m_Y_full_initial;
        Eigen::VectorXd y_qse = m_Y_scale.array() * v_qse.array().sinh(); // Convert to physical abundances using asinh scaling

        for (size_t i = 0; i < m_qse_solve_indices.size(); ++i) {
            y_trial[m_qse_solve_indices[i]] = y_qse(i);
        }

        const auto [dydt, energy] = m_view->calculateRHSAndEnergy(y_trial, m_T9, m_rho);
        f_qse.resize(m_qse_solve_indices.size());
        for (size_t i = 0; i < m_qse_solve_indices.size(); ++i) {
            f_qse(i) = dydt[m_qse_solve_indices[i]];
        }

        return 0; // Success
    }

    int MultiscalePartitioningEngineView::EigenFunctor::df(const InputType &v_qse, JacobianType &J_qse) const {
        std::vector<double> y_trial = m_Y_full_initial;
        Eigen::VectorXd y_qse = m_Y_scale.array() * v_qse.array().sinh(); // Convert to physical abundances using asinh scaling

        for (size_t i = 0; i < m_qse_solve_indices.size(); ++i) {
            y_trial[m_qse_solve_indices[i]] = y_qse(i);
        }

        // TODO: Think about if the jacobian matrix should be mutable so that generateJacobianMatrix can be const
        m_view->generateJacobianMatrix(y_trial, m_T9, m_rho);

        // TODO: Think very carefully about the indices here.
        J_qse.resize(m_qse_solve_indices.size(), m_qse_solve_indices.size());
        for (size_t i = 0; i < m_qse_solve_indices.size(); ++i) {
            for (size_t j = 0; j < m_qse_solve_indices.size(); ++j) {
                J_qse(i, j) = m_view->getJacobianMatrixEntry(
                    m_qse_solve_indices[i],
                    m_qse_solve_indices[j]
                );
            }
        }

        // Chain rule for asinh scaling:
        for (long j = 0; j < J_qse.cols(); ++j) {
            const double dY_dv = m_Y_scale(j) * std::cosh(v_qse(j));
            J_qse.col(j) *= dY_dv; // Scale the column by the derivative of the asinh scaling
        }
        return 0; // Success
    }
}
