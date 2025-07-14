#include "gridfire/engine/views/engine_multiscale.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>

namespace gridfire {
    static int s_operator_parens_called = 0;
    using fourdst::atomic::Species;

    MultiscalePartitioningEngineView::MultiscalePartitioningEngineView(
        GraphEngine& baseEngine
    ) : m_baseEngine(baseEngine) {}

    const std::vector<Species> & MultiscalePartitioningEngineView::getNetworkSpecies() const {
        return m_baseEngine.getNetworkSpecies();
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
        const double rho,
        const double dt_control
    ) {
        // --- Step 0. Clear previous state ---
        m_qse_groups.clear();
        m_dynamic_species.clear();
        m_dynamic_species_indices.clear();
        m_connectivity_graph.clear();

        // --- Step 1. Identify fast reactions ---
        const std::unordered_set<size_t> fast_reaction_indices = identifyFastReactions(Y, T9, rho, dt_control);

        // --- Step 2. Build connectivity graph based on fast reactions ---
        buildConnectivityGraph(fast_reaction_indices);

        // --- Step 3. Find connected components (candidate QSE clusters) ---
        findConnectedComponents();

        // --- Step 4. Validate clusters via flux analysis ---
        validateGroupsWithFluxAnalysis(Y, T9, rho);

        // --- Step 5. Identify dynamic species ---
        const auto [fst, snd] = identifyDynamicSpecies(Y, m_qse_groups, T9, rho);
        m_dynamic_species = fst;
        m_dynamic_species_indices = snd;

    }

    void MultiscalePartitioningEngineView::partitionNetwork(
        const NetIn &netIn,
        const double dt_control
    ) {
        std::vector<double> Y(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (const auto& [symbol, entry]: netIn.composition ) {
            if (!m_baseEngine.involvesSpecies(entry.isotope())) {
                LOG_ERROR(m_logger, "Species {} is not part of the network. Exiting...", entry.isotope().name());
                throw std::runtime_error("Species " + std::string(entry.isotope().name()) + " is not part of the network.");
            }
            const size_t species_index = m_baseEngine.getSpeciesIndex(entry.isotope());
            Y[species_index] = netIn.composition.getMolarAbundance(symbol);
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        partitionNetwork(Y, T9, rho, dt_control);
    }

    void MultiscalePartitioningEngineView::exportToDot(const std::string &filename) const {
        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file for writing: {}", filename);
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }

        dotFile << "digraph PartitionedNetwork {\n";
        dotFile << "    graph [rankdir=LR, splines=true, overlap=false, bgcolor=\"#fdfdfd\", label=\"Multiscale Partitioned Network View\"];\n";
        dotFile << "    node [shape=circle, style=filled, fillcolor=\"#eafaf1\", fontname=\"Helvetica\"];\n";
        dotFile << "    edge [fontname=\"Helvetica\", fontsize=10, color=\"#5d6d7e\"];\n\n";

        // 1. Define all species from the base network as nodes
        const auto& all_species = m_baseEngine.getNetworkSpecies();
        dotFile << "    // --- Species Nodes ---\n";
        for (const auto& species : all_species) {
            dotFile << "    \"" << species.name() << "\" [label=\"" << species.name() << "\"];\n";
        }
        dotFile << "\n";

        // 2. Define all reactions and their connections from the base network
        dotFile << "    // --- Reaction Edges ---\n";
        for (const auto& reaction : m_baseEngine.getNetworkReactions()) {
            std::string reactionNodeId = "reaction_" + std::string(reaction.id());
            dotFile << "    \"" << reactionNodeId << "\" [shape=point, fillcolor=black, width=0.08, height=0.08];\n";
            for (const auto& reactant : reaction.reactants()) {
                dotFile << "    \"" << reactant.name() << "\" -> \"" << reactionNodeId << "\";\n";
            }
            for (const auto& product : reaction.products()) {
                dotFile << "    \"" << reactionNodeId << "\" -> \"" << product.name() << "\";\n";
            }
            dotFile << "\n";
        }

        // 3. Draw clusters around the identified QSE groups
        dotFile << "    // --- QSE Group Clusters ---\n";
        int group_counter = 0;
        for (const auto& group : m_qse_groups) {
            dotFile << "    subgraph cluster_" << group_counter << " {\n";
            dotFile << "        style = \"filled,rounded\";\n";
            if (group.is_in_equilibrium) {
                dotFile << "        label = \"Valid QSE Group" << group_counter + 1 << "\";\n";
                dotFile << "        color = \"#fdebd0\";\n";
            } else {
                dotFile << "        label = \"Non valid QSE Group" << group_counter + 1 << "\";\n";
                dotFile << "        color = \"#d5f5e3\";\n"; // Different color for non-equilibrium groups
            }

            for (const size_t species_idx : group.species_indices) {
                dotFile << "        \"" << all_species[species_idx].name() << "\";\n";
            }

            dotFile << "    }\n";
            group_counter++;
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
        const double rho,
        const double dt_control
    ) {
        partitionNetwork(Y, T9, rho, dt_control);
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
            composition.setMassFraction(std::string(species.name()), X[i]);
        }

        composition.finalize(true);
        return composition;
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::equilibrateNetwork(
        const NetIn &netIn,
        const double dt_control
    ) {
        std::vector<double> Y(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (const auto& [symbol, entry]: netIn.composition ) {
            if (!m_baseEngine.involvesSpecies(entry.isotope())) {
                LOG_ERROR(m_logger, "Species {} is not part of the network. Exiting...", entry.isotope().name());
                throw std::runtime_error("Species " + std::string(entry.isotope().name()) + " is not part of the network.");
            }
            const size_t species_index = m_baseEngine.getSpeciesIndex(entry.isotope());
            Y[species_index] = netIn.composition.getMolarAbundance(symbol);
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        return equilibrateNetwork(Y, T9, rho, dt_control);
    }

    int MultiscalePartitioningEngineView::getSpeciesIndex(const fourdst::atomic::Species &species) const {
        return m_baseEngine.getSpeciesIndex(species);
    }

    std::unordered_set<size_t> MultiscalePartitioningEngineView::identifyFastReactions(
        const std::vector<double> &Y_full,
        const double T9,
        const double rho,
        const double dt_control
    ) const {
        std::unordered_set<size_t> fast_reaction_indices;
        std::unordered_set<size_t> fast_species_indices;

        const double timescale_threshold = 0.01 * dt_control;
        constexpr int FAST_BOOTSTRAPPING_ITERATIONS = 2;

        const auto& all_reactions = m_baseEngine.getNetworkReactions();

        // Loop a few times to refine the fast species and reactions.
        for (int iteration = 0; iteration < FAST_BOOTSTRAPPING_ITERATIONS; ++iteration) {
            // --- Step A: Identify fast species based on their individual reaction timescales ---
            for (size_t i = 0; i < all_reactions.size(); ++i) {
                const auto& reaction = all_reactions[i];

                // Calculate the forward molar flow, which represents the rate of destruction.
                // Note: We use the base engine's direct calculation method.
                const double forward_flow = m_baseEngine.calculateMolarReactionFlow(reaction, Y_full, T9, rho);

                if (forward_flow == 0.0) continue;

                // Determine the timescale for each reactant against this specific reaction.
                for (const auto& reactant : reaction.reactants()) {
                    // ~ Timescale = Abundance / Destruction Rate

                    // TODO: Should this call base engine getSpeciesIndex or this view's getSpeciesIndex?
                    const size_t reactant_idx = m_baseEngine.getSpeciesIndex(reactant);
                    const double abundance = Y_full[reactant_idx];

                    if (abundance > 0) {
                        const double timescale = abundance / forward_flow;
                        if (timescale < timescale_threshold) {
                            fast_species_indices.insert(reactant_idx);
                        }
                    }
                }
            }

            // --- Step B: Identify fast reactions based on reactants ---
            for (size_t i = 0; i < all_reactions.size(); ++i) {
                const auto& reaction = all_reactions[i];
                bool all_reactants_are_fast = true;
                if (reaction.reactants().empty()) {
                    all_reactants_are_fast = false;
                }
                for (const auto& reactant : reaction.reactants()) {
                    if (!fast_species_indices.contains(m_baseEngine.getSpeciesIndex(reactant))) {
                        all_reactants_are_fast = false;
                        break; // No need to check further if one reactant is not fast
                    }
                }
                if (all_reactants_are_fast) {
                    fast_reaction_indices.insert(i);
                }
            }

            // --- Step C: Propagate "fast" status to products of the fast reactions (this handles things like n-1) ---
            for (size_t reaction_idx : fast_reaction_indices) {
                for (const auto& product : all_reactions[reaction_idx].products()) {
                    fast_species_indices.insert(m_baseEngine.getSpeciesIndex(product));
                }
            }
        }
        return fast_reaction_indices;
    }

    void MultiscalePartitioningEngineView::buildConnectivityGraph(
        const std::unordered_set<size_t> &fast_reaction_indices
    ) {
        const auto& all_reactions = m_baseEngine.getNetworkReactions();
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
                    m_connectivity_graph[reactant_idx].push_back(product_idx);
                    m_connectivity_graph[product_idx].push_back(reactant_idx);
                }
            }
        }


    }

    void MultiscalePartitioningEngineView::findConnectedComponents() {
        const size_t num_total_species = m_baseEngine.getNetworkSpecies().size();
        std::vector<bool> visited(num_total_species, false);

        for (size_t i = 0; i < num_total_species; ++i) {
            if (!visited[i]) {
                // Start a new BFS traversal from this unvisited node.
                // This will find one complete connected component.
                QSEGroup new_group;
                std::queue<size_t> q;

                q.push(i);
                visited[i] = true;

                while (!q.empty()) {
                    const size_t current_species_idx = q.front();
                    q.pop();
                    new_group.species_indices.push_back(current_species_idx);

                    // Check if the node exists in the graph (it might be an isolated species)
                    if (m_connectivity_graph.contains(current_species_idx)) {
                        for (const size_t neighbor_idx : m_connectivity_graph.at(current_species_idx)) {
                            if (!visited[neighbor_idx]) {
                                visited[neighbor_idx] = true;
                                q.push(neighbor_idx);
                            }
                        }
                    }
                }

                // A "group" must contain more than one species to be a cluster.
                // Isolated species are treated as dynamic by default.
                if (new_group.species_indices.size() > 1) {
                    m_qse_groups.push_back(new_group);
                }
            }
        }

    }

    void MultiscalePartitioningEngineView::validateGroupsWithFluxAnalysis(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) {
        constexpr double FLUX_RATIO_THRESHOLD = 100;
        for (auto& group : m_qse_groups) {
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
    }

    std::pair<std::vector<Species>, std::vector<size_t>> MultiscalePartitioningEngineView::identifyDynamicSpecies(
        const std::vector<double>& Y,
        const std::vector<QSEGroup>& qse_groups,
        const double T9,
        const double rho
    ) const {
        // This set will contain the indices of species that are truly in qse (as opposed to those which might be in a qse group but are not in equilibrium).
        std::unordered_set<size_t> algebraic_qse_species_indices;
        auto species_timescales = m_baseEngine.getSpeciesTimescales(Y, T9, rho);

        for (auto& group : qse_groups) {
            if (!group.is_in_equilibrium || group.species_indices.empty()) continue;

            for (size_t species_idx : group.species_indices) {
                const double timescale = species_timescales.at(m_baseEngine.getNetworkSpecies()[species_idx]);

                if (timescale < 1) { // If the timescale is less than 1 second, we consider it algebraic.
                    algebraic_qse_species_indices.insert(species_idx);
                }
            }
        }

        std::vector<fourdst::atomic::Species> dynamic_species;
        std::unordered_set<size_t> qse_species_indices_set;
        std::vector<size_t> dynamic_species_indices;
        const auto& all_base_species = m_baseEngine.getNetworkSpecies();

        for (size_t i = 0; i < all_base_species.size(); ++i) {
            if (!algebraic_qse_species_indices.contains(i)) {
                dynamic_species.push_back(all_base_species[i]);
                dynamic_species_indices.push_back(i);
            }
        }
        return {dynamic_species, dynamic_species_indices};
    }

    std::vector<double> MultiscalePartitioningEngineView::solveQSEAbundances(
        const std::vector<double> &Y_full,
        const double T9,
        const double rho
    ) {
        auto Y_out = Y_full;
        auto species_timescales = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);

        for (const auto& qse_group : m_qse_groups) {
            if (!qse_group.is_in_equilibrium || qse_group.species_indices.empty()) continue;
            std::vector<size_t> qse_solve_indices;
            std::vector<size_t> seed_indices;

            for (size_t species_idx : qse_group.species_indices) {
                const double timescale = species_timescales.at(m_baseEngine.getNetworkSpecies()[species_idx]);
                if (timescale < 1.0) { // If the timescale is less than 1 second, we consider it algebraic.
                    qse_solve_indices.push_back(species_idx);
                } else {
                    seed_indices.push_back(species_idx);
                }
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
            Eigen::VectorXd Y_final_qse = Y_scale.array() * v_initial.array().sinh(); // Convert back to physical abundances using asinh scaling
            for (size_t i = 0; i < qse_solve_indices.size(); ++i) {
                Y_out[qse_solve_indices[i]] = Y_final_qse(i);
            }
        }
        return Y_out;
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
