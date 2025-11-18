#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/exceptions/error_engine.h"
#include "gridfire/engine/procedures/priming.h"

#include <stdexcept>
#include <vector>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

#include <queue>

#include <algorithm>

#include "fourdst/atomic/species.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"

#include  "xxhash64.h"
#include "fourdst/composition/utils/composition_hash.h"

namespace {
    constexpr double MIN_ABS_NORM_ALWAYS_CONVERGED = 1e-30;

    using namespace fourdst::atomic;


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

    std::vector<std::vector<Species>> findConnectedComponentsBFS(
        const std::unordered_map<Species, std::set<Species>>& graph,
        const std::vector<Species>& nodes
    ) {
        std::unordered_map<Species, std::vector<Species>> adjList;
        for (const auto& [u, neighbors] : graph) {
            adjList[u] = std::vector<Species>(neighbors.begin(), neighbors.end());
        }
        return findConnectedComponentsBFS(adjList, nodes);
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
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        LOG_TRACE_L2(m_logger, "Calculating RHS and Energy in MultiscalePartitioningEngineView at T9 = {}, rho = {}.", T9, rho);
        LOG_TRACE_L2(m_logger, "Input composition is {}", [&comp]() -> std::string {
            std::stringstream ss;
            size_t i = 0;
            for (const auto& [species, abundance] : comp) {
                ss << species.name() << ": " << abundance;
                if (i < comp.size() - 1) {
                    ss << ", ";
                }
                i++;
            }
            return ss.str();
        }());
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        LOG_TRACE_L2(m_logger, "Equilibrated composition prior to calling base engine is {}", [&qseComposition, &comp]() -> std::string {
            std::stringstream ss;
            size_t i = 0;
            for (const auto& [species, abundance] : qseComposition) {
                ss << species.name() << ": " << abundance;
                if (comp.contains(species)) {
                    ss << " (input: " << comp.getMolarAbundance(species) << ")";
                }
                if (i < qseComposition.size() - 1) {
                    ss << ", ";
                }
                i++;
            }
            return ss.str();
        }());

        const auto result = m_baseEngine.calculateRHSAndEnergy(qseComposition, T9, rho);
        LOG_TRACE_L2(m_logger, "Base engine calculation of RHS and Energy complete.");

        if (!result) {
            LOG_TRACE_L2(m_logger, "Base engine returned stale error during RHS and Energy calculation.");
            return std::unexpected{result.error()};
        }

        auto deriv = result.value();

        LOG_TRACE_L2(m_logger, "Zeroing out algebraic species derivatives.");
        for (const auto& species : m_algebraic_species) {
            deriv.dydt[species] = 0.0; // Fix the algebraic species to the equilibrium abundances we calculate.
        }
        LOG_TRACE_L2(m_logger, "Done Zeroing out algebraic species derivatives.");
        return deriv;
    }

    EnergyDerivatives MultiscalePartitioningEngineView::calculateEpsDerivatives(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        return m_baseEngine.calculateEpsDerivatives(qseComposition, T9, rho);
    }

    NetworkJacobian MultiscalePartitioningEngineView::generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        return m_baseEngine.generateJacobianMatrix(qseComposition, T9, rho, m_dynamic_species);
    }

    NetworkJacobian MultiscalePartitioningEngineView::generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::vector<Species> &activeSpecies
    ) const {
        bool activeSpeciesIsSubset = true;
        for (const auto& species : activeSpecies) {
            if (!involvesSpecies(species)) activeSpeciesIsSubset = false;
        }
        if (!activeSpeciesIsSubset) {
            std::string msg = std::format(
                "Active species set contains species ({}) not present in network partition. Cannot generate jacobian matrix due to this.",
                [&]() -> std::string {
                    std::stringstream ss;
                    for (const auto& species : activeSpecies) {
                        if (!this->involvesSpecies(species)) {
                            ss << species << " ";
                        }
                    }
                    return ss.str();
                }()
            );
            LOG_CRITICAL(m_logger, "{}", msg);
            throw std::runtime_error(msg);
        }

        std::vector<Species> dynamicActiveSpeciesIntersection;
        for (const auto& species : activeSpecies) {
            if (involvesSpeciesInDynamic(species)) {
                dynamicActiveSpeciesIntersection.push_back(species);
            }
        }

        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);

        return m_baseEngine.generateJacobianMatrix(qseComposition, T9, rho, dynamicActiveSpeciesIntersection);
    }

    NetworkJacobian MultiscalePartitioningEngineView::generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const SparsityPattern &sparsityPattern
    ) const {
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        return m_baseEngine.generateJacobianMatrix(qseComposition, T9, rho, sparsityPattern);
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
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);

        return m_baseEngine.calculateMolarReactionFlow(reaction, qseComposition, T9, rho);
    }

    const reaction::ReactionSet & MultiscalePartitioningEngineView::getNetworkReactions() const {
        return m_baseEngine.getNetworkReactions();
    }

    void MultiscalePartitioningEngineView::setNetworkReactions(const reaction::ReactionSet &reactions) {
        LOG_CRITICAL(m_logger, "setNetworkReactions is not supported in MultiscalePartitioningEngineView. Did you mean to call this on the base engine?");
        throw exceptions::UnableToSetNetworkReactionsError("setNetworkReactions is not supported in MultiscalePartitioningEngineView. Did you mean to call this on the base engine?");
    }

    std::expected<std::unordered_map<Species, double>, expectations::StaleEngineError> MultiscalePartitioningEngineView::getSpeciesTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        const auto result  = m_baseEngine.getSpeciesTimescales(qseComposition, T9, rho);
        if (!result) {
            return std::unexpected{result.error()};
        }
        std::unordered_map<Species, double> speciesTimescales = result.value();
        for (const auto& algebraicSpecies : m_algebraic_species) {
            speciesTimescales[algebraicSpecies] = std::numeric_limits<double>::infinity(); // Algebraic species have infinite timescales.
        }
        return speciesTimescales;
    }

    std::expected<std::unordered_map<Species, double>, expectations::StaleEngineError>
    MultiscalePartitioningEngineView::getSpeciesDestructionTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        const auto result = m_baseEngine.getSpeciesDestructionTimescales(qseComposition, T9, rho);
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
        fourdst::composition::Composition equilibratedComposition = partitionNetwork(baseUpdatedNetIn);

        return equilibratedComposition;
    }

    bool MultiscalePartitioningEngineView::isStale(const NetIn &netIn) {
        return m_baseEngine.isStale(netIn);
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
        const std::vector<std::vector<Species>> &timescale_pools
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

    std::vector<MultiscalePartitioningEngineView::QSEGroup> MultiscalePartitioningEngineView::pruneValidatedGroups(
        const std::vector<QSEGroup> &groups,
        const std::vector<reaction::ReactionSet> &groupReactions,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        const auto result = m_baseEngine.getSpeciesTimescales(comp, T9, rho);
        if (!result) {
            throw std::runtime_error("Base engine returned stale error during pruneValidatedGroups timescale retrieval.");
        }
        std::unordered_map<Species, double> speciesTimescales = result.value();
        std::vector<QSEGroup> newGroups;
        for (const auto &[group, reactions] : std::views::zip(groups, groupReactions)) {
            std::unordered_map<size_t, double> reactionFluxes;
            std::unordered_map<size_t, const reaction::Reaction&> reactionLookup;
            reactionFluxes.reserve(reactions.size());

            double mean_molar_abundance = 0;
            for (const auto& species : group.algebraic_species) {
                mean_molar_abundance += comp.getMolarAbundance(species);
            }
            mean_molar_abundance /= group.algebraic_species.size();
            { // Safety Valve to ensure valid log scaling
                if (mean_molar_abundance <= 0) {
                    LOG_CRITICAL(m_logger, "Non-positive mean molar abundance {} calculated for QSE group during pruning analysis.", mean_molar_abundance);
                    throw std::runtime_error("Non-positive mean molar abundance calculated during pruneValidatedGroups flux analysis.");
                }
            }

            for (const auto& reaction : reactions) {
                const double flux = m_baseEngine.calculateMolarReactionFlow(*reaction, comp, T9, rho);
                size_t hash = reaction->hash(0);
                if (reactionFluxes.contains(hash)) {
                    throw std::runtime_error("Duplicate reaction hash found during pruneValidatedGroups flux analysis.");
                }

                { // Safety Valve to ensure valid log scaling
                    if (flux <= 0) {
                        LOG_CRITICAL(m_logger, "Non-positive flux {} calculated for reaction {} during pruning analysis.", flux, reaction->id());
                        throw std::runtime_error("Non-positive flux calculated during pruneValidatedGroups flux analysis.");
                    }
                }
                double lAbundanceNormalizedFlux = std::log(flux/mean_molar_abundance);
                reactionFluxes.emplace(hash, lAbundanceNormalizedFlux);

                assert(!std::isnan(lAbundanceNormalizedFlux) && !std::isinf(lAbundanceNormalizedFlux) && "Invalid log abundance normalized flux calculated during pruneValidatedGroups flux analysis.");

                reactionLookup.emplace(hash, *reaction);
            }

            std::vector<size_t> sorted_reactions_based_on_flow;
            for (const auto &hash: reactionFluxes | std::views::keys) {
                sorted_reactions_based_on_flow.push_back(hash);
            }

            std::ranges::sort(sorted_reactions_based_on_flow, [&reactionFluxes](const size_t a, const size_t b) {
                return std::abs(reactionFluxes.at(a)) < std::abs(reactionFluxes.at(b));
            });

            std::unordered_map<size_t, double> pruned_reaction_fluxes;
            for (const auto& [hash, normalizedFlux] : reactionFluxes) {
                if (normalizedFlux > -30) { // TODO: replace -30 with some more physically motivated value
                    pruned_reaction_fluxes.emplace(hash, normalizedFlux);
                    LOG_TRACE_L2(m_logger, "Retaining reaction {} with log(mean abundance normalized flux) {} during pruning.", reactionLookup.at(hash).id(), normalizedFlux);
                } else {
                    LOG_TRACE_L2(m_logger, "Pruning reaction {} with log(mean abundance normalized flux) {} during pruning.", reactionLookup.at(hash).id(), normalizedFlux);
                }
            }

            std::set<Species> reachableAlgebraicSpecies;
            std::set<Species> reachableSeedSpecies;

            std::unordered_map<Species, std::set<Species>> connectivity_graph;

            for (const auto& reactionHash : pruned_reaction_fluxes | std::views::keys) {
                const auto& reaction = reactionLookup.at(reactionHash);

                for (const auto& reactant : reaction.reactants()) {
                    if (group.algebraic_species.contains(reactant)) {
                        reachableAlgebraicSpecies.insert(reactant);
                    } else if (group.seed_species.contains(reactant)) {
                        reachableSeedSpecies.insert(reactant);
                    }

                    if (!connectivity_graph.contains(reactant)) {
                        connectivity_graph.emplace(reactant, std::set<Species>{});
                    }
                    for (const auto& product : reaction.products()) {
                        connectivity_graph.at(reactant).insert(product);
                    }
                }

                for (const auto& product : reaction.products()) {
                    if (group.algebraic_species.contains(product)) {
                        reachableAlgebraicSpecies.insert(product);
                    } else if (group.seed_species.contains(product)) {
                        reachableSeedSpecies.insert(product);
                    }
                }

            }

            LOG_TRACE_L2(
                m_logger,
                "{}",
                [&group, &reachableAlgebraicSpecies, &reachableSeedSpecies]() -> std::string {
                    std::stringstream ss;
                    ss << "Pruned QSE Group. Group Started with Algebraic Species: {";
                    int i = 0;
                    for (const auto& species : group.algebraic_species) {
                        ss << species.name();
                        if (i < group.algebraic_species.size() - 1) {
                            ss << ", ";
                        }
                        i++;
                    }
                    ss << "} and Seed Species: {";
                    i = 0;
                    for (const auto& species : group.seed_species) {
                        ss << species.name();
                        if (i < group.seed_species.size() - 1) {
                            ss << ", ";
                        }
                        i++;
                    }
                    ss << "}. After pruning, reachable Algebraic Species: {";
                    i = 0;
                    for (const auto& species : reachableAlgebraicSpecies) {
                        ss << species.name();
                        if (i < reachableAlgebraicSpecies.size() - 1) {
                            ss << ", ";
                        }
                        i++;
                    }
                    ss << "} and reachable Seed Species: {";
                    i = 0;
                    for (const auto& species : reachableSeedSpecies) {
                        ss << species.name();
                        if (i < reachableSeedSpecies.size() - 1) {
                            ss << ", ";
                        }
                        i++;
                    }
                    ss << "}.";
                    return ss.str();
                }()
            );

            std::vector<std::vector<Species>> connected_components = findConnectedComponentsBFS(
                connectivity_graph,
                std::vector<Species>(
                    reachableAlgebraicSpecies.begin(),
                    reachableAlgebraicSpecies.end()
                )
            );

            for (const auto& subgraph: connected_components) {
                QSEGroup g;
                for (const auto& species: subgraph) {
                    if (reachableAlgebraicSpecies.contains(species)) {
                        g.algebraic_species.insert(species);
                    } else if (reachableSeedSpecies.contains(species)) {
                        g.seed_species.insert(species);
                    }
                }
                if (!g.seed_species.empty() && !g.algebraic_species.empty()) {
                    double meanTimescale = 0;
                    for (const auto &species : g.algebraic_species) {
                        meanTimescale += speciesTimescales.at(species);
                    }
                    meanTimescale /= g.algebraic_species.size();
                    g.mean_timescale = meanTimescale;
                    newGroups.push_back(g);
                }
            }
        }
        return newGroups;
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::partitionNetwork(
        const NetIn &netIn
    ) {
        // --- Step 0. Prime the network ---
        const PrimingReport primingReport = m_baseEngine.primeEngine(netIn);
        const fourdst::composition::Composition& comp = primingReport.primedComposition;
        const double T9 = netIn.temperature / 1e9;
        const double rho = netIn.density;

        // --- Step 0.5 Clear previous state ---
        LOG_TRACE_L1(m_logger, "Partitioning network...");
        LOG_TRACE_L1(m_logger, "Clearing previous state...");
        m_qse_groups.clear();
        m_dynamic_species.clear();
        m_algebraic_species.clear();
        m_composition_cache.clear(); // We need to clear the cache now cause the same comp, temp, and density may result in a different value

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
        const std::vector<std::vector<Species>> connected_pools = analyzeTimescalePoolConnectivity(candidate_pools);
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
        const auto [validated_groups, invalidate_groups, validated_group_reactions] = validateGroupsWithFluxAnalysis(candidate_groups, comp, T9, rho);
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

        LOG_TRACE_L1(m_logger, "Pruning groups based on log abundance-normalized flux analysis...");
        const std::vector<QSEGroup> prunedGroups = pruneValidatedGroups(validated_groups, validated_group_reactions, comp, T9, rho);
        LOG_TRACE_L1(m_logger, "After Pruning remaining groups are: {}", [&]() -> std::string {
            std::stringstream ss;
            int j = 0;
            for (const auto& group : prunedGroups) {
                ss << "PrunedQSEGroup(Algebraic: {";
                int i = 0;
                for (const auto& species : group.algebraic_species) {
                    ss << species.name();
                    if (i < group.algebraic_species.size() - 1) {
                        ss << ", ";
                    }
                    i++;
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
                if (j < prunedGroups.size() - 1) {
                    ss << ", ";
                }
                j++;
            }
            return ss.str();
        }());
        LOG_TRACE_L1(m_logger, "Re-validating pruned groups with flux analysis...");
        auto [pruned_validated_groups, _, __] = validateGroupsWithFluxAnalysis(prunedGroups, comp, T9, rho);
        LOG_TRACE_L1(
            m_logger,
            "After re-validation, {} QSE groups remain. ({})",
            pruned_validated_groups.size(),
            [&pruned_validated_groups]()->std::string {
                std::stringstream ss;
                size_t count = 0;
                for (const auto& group : pruned_validated_groups) {
                    ss << group.toString();
                    if (pruned_validated_groups.size() > 1 && count < pruned_validated_groups.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                return ss.str();
            }()
        );

        m_qse_groups = pruned_validated_groups;
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

        // Sort the QSE groups by mean timescale so that fastest groups get equilibrated first (as these may feed slower groups)
        std::ranges::sort(m_qse_groups, [](const QSEGroup& a, const QSEGroup& b) {
            return a.mean_timescale < b.mean_timescale;
        });

        for (const auto& species : m_baseEngine.getNetworkSpecies()) {
            const bool involvesAlgebraic = involvesSpeciesInQSE(species);
            if (std::ranges::find(m_dynamic_species, species) == m_dynamic_species.end() && !involvesAlgebraic) {
                // Species is classed as neither dynamic nor algebraic at end of partitioning → add to dynamic set. This ensures that all species are classified.
                m_dynamic_species.push_back(species);
            }
        }
        return getNormalizedEquilibratedComposition(comp, T9, rho);
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

        const fourdst::composition::Composition qseComposition = getNormalizedEquilibratedComposition(comp, T9, rho);
        // Calculate reaction flows and find min/max for logarithmic scaling of transparency
        std::vector<double> reaction_flows;
        reaction_flows.reserve(all_reactions.size());
        double min_log_flow = std::numeric_limits<double>::max();
        double max_log_flow = std::numeric_limits<double>::lowest();

        for (const auto& reaction : all_reactions) {
            double flow = std::abs(m_baseEngine.calculateMolarReactionFlow(*reaction, qseComposition, T9, rho));
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
        std::vector<double> Y(netIn.composition.size(), 0.0); // Initialize with zeros
        for (const auto& [sp, y] : netIn.composition) {
            Y[getSpeciesIndex(sp)] = y; // Map species to their molar abundance
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


    bool MultiscalePartitioningEngineView::involvesSpecies(
        const Species &species
    ) const {
        if (involvesSpeciesInQSE(species)) return true; // check this first since the vector is likely to be smaller so short circuit cost is less on average
        if (involvesSpeciesInDynamic(species)) return true;
        return false;
    }

    bool MultiscalePartitioningEngineView::involvesSpeciesInQSE(
        const Species &species
    ) const {
        return std::ranges::find(m_algebraic_species, species) != m_algebraic_species.end();
    }

    bool MultiscalePartitioningEngineView::involvesSpeciesInDynamic(
        const Species &species
    ) const {
        return std::ranges::find(m_dynamic_species, species) != m_dynamic_species.end();
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::getNormalizedEquilibratedComposition(
        const fourdst::composition::CompositionAbstract& comp,
        const double T9,
        const double rho
    ) const {
        // Caching mechanism to avoid redundant QSE solves
        const std::array<uint64_t, 3> hashes = {
            fourdst::composition::utils::CompositionHash::hash_exact(comp),
            std::hash<double>()(T9),
            std::hash<double>()(rho)
        };

        const uint64_t composite_hash = XXHash64::hash(hashes.begin(), sizeof(uint64_t) * 3, 0);
        if (m_composition_cache.contains(composite_hash)) {
            LOG_TRACE_L3(m_logger, "Cache Hit in Multiscale Partitioning Engine View for composition at T9 = {}, rho = {}.", T9, rho);
            return m_composition_cache.at(composite_hash);
        }
        LOG_TRACE_L3(m_logger, "Cache Miss in Multiscale Partitioning Engine View for composition at T9 = {}, rho = {}. Solving QSE abundances...", T9, rho);

        // Only solve if the composition and thermodynamic conditions have not been cached yet
        fourdst::composition::Composition qseComposition = solveQSEAbundances(comp, T9, rho);

        for (const auto &[sp, y]: qseComposition) {
            if (y < 0.0 && std::abs(y) < 1e-20) {
                qseComposition.setMolarAbundance(sp, 0.0); // normalize small negative abundances to zero
            }
        }
        m_composition_cache[composite_hash] = qseComposition;

        return qseComposition;
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::collectComposition(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const fourdst::composition::Composition result = m_baseEngine.collectComposition(comp, T9, rho);
        fourdst::composition::Composition qseComposition = solveQSEAbundances(result, T9, rho);

        return qseComposition;
    }

    SpeciesStatus MultiscalePartitioningEngineView::getSpeciesStatus(const Species &species) const {
        const SpeciesStatus status = m_baseEngine.getSpeciesStatus(species);
        if (status == SpeciesStatus::ACTIVE && involvesSpeciesInQSE(species)) {
            return SpeciesStatus::EQUILIBRIUM;
        }
        return status;
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
        [[maybe_unused]] const std::unordered_map<Species, double>& net_timescales = netTimescale.value();

        LOG_TRACE_L3(
            m_logger,
            "{}",
            [&]() -> std::string {
                std::stringstream ss;
                for (const auto& [species, destruction_timescale] : destruction_timescales) {
                    ss << std::format("For {} destruction timescale is {}s\n", species.name(), destruction_timescale);
                }
                return ss.str();
            }()
        );

        const auto& all_species = m_baseEngine.getNetworkSpecies();

        std::vector<std::pair<double, Species>> sorted_destruction_timescales;
        for (const auto & species : all_species) {
            double destruction_timescale = destruction_timescales.at(species);
            if (std::isfinite(destruction_timescale) && destruction_timescale > 0) {
                LOG_TRACE_L2(m_logger, "Species {} has finite destruction timescale: destruction: {} s, net: {} s", species.name(), destruction_timescale, net_timescales.at(species));
                sorted_destruction_timescales.emplace_back(destruction_timescale, species);
            } else {
                LOG_TRACE_L2(m_logger, "Species {} has infinite or negative destruction timescale: destruction: {} s, net: {} s", species.name(), destruction_timescale, net_timescales.at(species));
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
        constexpr double MAX_MOLAR_ABUNDANCE_THRESHOLD = 1e-10; // Maximum molar abundance which a fast species is allowed to have (anything more abundant is always considered dynamic)
        constexpr double MIN_MOLAR_ABUNDANCE_THRESHOLD = 1e-50; // Minimum molar abundance to consider a species at all (anything less abundance will be classed as dynamic but with the intent that some latter view will deal with it)

        LOG_TRACE_L2(m_logger, "Found {} species with finite timescales.", sorted_destruction_timescales.size());
        LOG_TRACE_L2(m_logger, "Absolute QSE timescale threshold: {} seconds ({} years).",
            ABSOLUTE_QSE_TIMESCALE_THRESHOLD, ABSOLUTE_QSE_TIMESCALE_THRESHOLD / 3.156e7);
        LOG_TRACE_L2(m_logger, "Minimum gap threshold: {} orders of magnitude.", MIN_GAP_THRESHOLD);
        LOG_TRACE_L2(m_logger, "Maximum molar abundance threshold for fast species consideration : {}.", MAX_MOLAR_ABUNDANCE_THRESHOLD);
        LOG_TRACE_L2(m_logger, "Minimum molar abundance threshold for species consideration : {}.", MIN_MOLAR_ABUNDANCE_THRESHOLD);

        std::vector<Species> dynamic_pool_species;
        std::vector<std::pair<double, Species>> fast_candidates;

        // 1. First Pass: Absolute Timescale Cutoff
        for (const auto& [destruction_timescale, species] : sorted_destruction_timescales) {
            if (species == n_1) {
                LOG_TRACE_L2(m_logger, "Skipping neutron (n) from timescale analysis. Neutrons are always considered dynamic due to their extremely high connectivity.");
                dynamic_pool_species.push_back(species);

                continue;
            }
            if (destruction_timescale > ABSOLUTE_QSE_TIMESCALE_THRESHOLD) {
                LOG_TRACE_L2(m_logger, "Species {} with timescale {} is considered dynamic (slower than qse timescale threshold).",
                    species.name(), destruction_timescale);
                dynamic_pool_species.push_back(species);
            } else {
                const double Yi = comp.getMolarAbundance(species);
                if (Yi > MAX_MOLAR_ABUNDANCE_THRESHOLD) {
                    LOG_TRACE_L2(m_logger, "Species {} with abundance {} is considered dynamic (above minimum abundance threshold of {}).",
                        species.name(), Yi, MAX_MOLAR_ABUNDANCE_THRESHOLD);
                    dynamic_pool_species.push_back(species);
                    continue;
                }
                if (Yi < MIN_MOLAR_ABUNDANCE_THRESHOLD) {
                    LOG_TRACE_L2(m_logger, "Species {} with abundance {} is considered dynamic (below minimum abundance threshold of {}). Likely another network view (such as adaptive engine view) will be needed to deal with this species",
                        species.name(), Yi, MIN_MOLAR_ABUNDANCE_THRESHOLD);
                    dynamic_pool_species.push_back(species);
                    continue;
                }
                LOG_TRACE_L2(m_logger, "Species {} with timescale {} and molar abundance {} is a candidate fast species (faster than qse timescale threshold and less than the molar abundance threshold).",
                    species.name(), destruction_timescale, Yi);
                fast_candidates.emplace_back(destruction_timescale, species);
            }
        }

        if (!dynamic_pool_species.empty()) {
            LOG_TRACE_L2(m_logger, "Found {} dynamic species (slower than QSE timescale threshold).", dynamic_pool_species.size());
            final_pools.push_back(dynamic_pool_species);
        }

        if (fast_candidates.empty()) {
            LOG_TRACE_L2(m_logger, "No fast candidates found.");
            return final_pools;
        }

        // 2. Second Pass: Gap Detection on the remaining "fast" candidates
        std::vector<size_t> split_points;
        for (size_t i = 0; i < fast_candidates.size() - 1; ++i) {
            const double t1 = fast_candidates[i].first;
            const double t2 = fast_candidates[i+1].first;
            if (std::log10(t1) - std::log10(t2) > MIN_GAP_THRESHOLD) {
                LOG_TRACE_L2(m_logger, "Detected gap between species {} (timescale {:0.2E}) and {} (timescale {:0.2E}).",
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

        LOG_TRACE_L2(m_logger, "Final partitioned pools: {}",
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

        LOG_TRACE_L2(
            m_logger,
            "Species Timescales: {}",
            [&]() -> std::string {
                std::stringstream ss;
                size_t poolID = 0;
                for (const auto& pool : final_pools) {
                    ss << "Pool #" << poolID << " [";
                    int ic = 0;
                    for (const auto& species : pool) {
                        ss << species << ": " << destruction_timescales.at(species);
                        if (ic < pool.size() - 1) {
                            ss << ", ";
                        }
                        ic++;
                    }
                    ss << "]";
                    poolID++;
                }
                return ss.str();
            }()
            );
        return final_pools;

    }

    MultiscalePartitioningEngineView::FluxValidationResult MultiscalePartitioningEngineView::validateGroupsWithFluxAnalysis(
        const std::vector<QSEGroup> &candidate_groups,
        const fourdst::composition::Composition &comp,
        const double T9, const double rho
    ) const {
        std::vector<QSEGroup> validated_groups;
        std::vector<reaction::ReactionSet> group_reactions;
        std::vector<QSEGroup> invalidated_groups;
        validated_groups.reserve(candidate_groups.size());
        group_reactions.reserve(candidate_groups.size());
        for (auto& group : candidate_groups) {
            reaction::ReactionSet group_reaction_set;
            constexpr double FLUX_RATIO_THRESHOLD = 5;

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
                    }

                }

                bool has_internal_algebraic_product = false;

                for (const auto& product : reaction->products()) {
                    if (algebraic_group_members.contains(product)) {
                        has_internal_algebraic_product = true;
                        LOG_TRACE_L3(m_logger, "Adjusting creation flux (+= {} mol g^-1 s^-1) for QSEGroup due to product {} from reaction {}",
                            flow, product.name(), reaction->id());
                    }
                }

                if (!has_internal_algebraic_product && !has_internal_algebraic_reactant) {
                    LOG_TRACE_L3(m_logger, "{}: Skipping reaction {} as it has no internal algebraic species in reactants or products.", group.toString(), reaction->id());
                    continue;
                }

                group_reaction_set.add_reaction(reaction->clone());

                LOG_TRACE_L2(
                    m_logger,
                    "Reaction {} (Coupling {}) has flow {} mol g^-1 s^-1 contributing to QSEGroup {}",
                    reaction->id(),
                    [&group, &reaction]() -> std::string {
                        std::ostringstream ss;
                        if (group.algebraic_species.empty()) {
                            ss << "N/A (no algebraic species)";
                        } else {
                            // Make a string of all the group species coupled from the reaction in the form of
                            // "A, B -> C, D"
                            int count = 0;
                            for (const auto& species : group.algebraic_species) {
                                if (std::ranges::find(reaction->reactants(), species) != reaction->reactants().end()) {
                                    ss << species.name();
                                    if (count < group.algebraic_species.size() - 1) {
                                        ss << ", ";
                                    }
                                    count++;
                                }
                            }
                            ss << " -> ";
                            count = 0;
                            for (const auto& species : group.algebraic_species) {
                                if (std::ranges::find(reaction->products(), species) != reaction->products().end()) {
                                    ss << species.name();
                                    if (count < group.algebraic_species.size() - 1) {
                                        ss << ", ";
                                    }
                                    count++;
                                }
                            }
                        }
                        return ss.str();

                    }(),
                    flow,
                    group.toString()
                );

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

            if (coupling_flux / leakage_flux > FLUX_RATIO_THRESHOLD) {
                LOG_TRACE_L1(
                    m_logger,
                    "Group containing {} is in equilibrium due to high coupling flux and balanced creation and destruction rate: <coupling: leakage flux = {}, coupling flux = {}, ratio = {} (Threshold: {})>",
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
                    FLUX_RATIO_THRESHOLD
                );
                validated_groups.emplace_back(group);
                validated_groups.back().is_in_equilibrium = true;
                group_reactions.emplace_back(group_reaction_set);
            } else {
                LOG_TRACE_L1(
                    m_logger,
                    "Group containing {} is NOT in equilibrium: <coupling: leakage flux = {}, coupling flux = {}, ratio = {} (Threshold: {})>",
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
                    FLUX_RATIO_THRESHOLD
                );
                invalidated_groups.emplace_back(group);
                invalidated_groups.back().is_in_equilibrium = false;
            }
        }

        LOG_TRACE_L1(m_logger, "Validated {} QSE groups and invalidated {} QSE groups after flux analysis.", validated_groups.size(), invalidated_groups.size());
        return {validated_groups, invalidated_groups, group_reactions};
    }

    fourdst::composition::Composition MultiscalePartitioningEngineView::solveQSEAbundances(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        LOG_TRACE_L1(m_logger, "Solving for QSE abundances...");
        LOG_TRACE_L2(m_logger, "Composition before QSE solving: {}", [&comp]() -> std::string {
            std::stringstream ss;
            size_t i = 0;
            for (const auto& [sp, y] : comp) {
                ss << std::format("{}: {}", sp.name(), y);
                if (i < comp.size() - 1) {
                    ss << ", ";
                }
                i++;
            }
            return ss.str();
        }());

        fourdst::composition::Composition outputComposition(comp);

        for (const auto&[is_in_equilibrium, algebraic_species, seed_species, mean_timescale] : m_qse_groups) {
            LOG_TRACE_L2(m_logger, "Working on QSE group with algebraic species: {}",
                [&]() -> std::string {
                    std::stringstream ss;
                    size_t count = 0;
                    for (const auto& species: algebraic_species) {
                        ss << species.name();
                        if (count < algebraic_species.size() - 1) {
                            ss << ", ";
                        }
                        count++;
                    }
                    return ss.str();
                }());
            if (!is_in_equilibrium || (algebraic_species.empty() && seed_species.empty())) {
                continue;
            }

            Eigen::VectorXd Y_scale(algebraic_species.size());
            Eigen::VectorXd v_initial(algebraic_species.size());
            long i = 0;
            std::unordered_map<Species, size_t> species_to_index_map;
            for (const auto& species : algebraic_species) {
                constexpr double abundance_floor = 1.0e-100;
                const double initial_abundance = comp.getMolarAbundance(species);
                const double Y = std::max(initial_abundance, abundance_floor);
                v_initial(i) = std::log(Y);
                species_to_index_map.emplace(species, i);
                LOG_TRACE_L2(m_logger, "For species {} initial molar abundance is {}, log scaled to {}. Species placed at index {}.", species.name(), Y, v_initial(i), i);
                i++;
            }

            LOG_TRACE_L2(m_logger, "Setting up Eigen Levenberg-Marquardt solver for QSE group...");
            EigenFunctor functor(*this, algebraic_species, comp, T9, rho, Y_scale, species_to_index_map);
            Eigen::LevenbergMarquardt lm(functor);
            lm.parameters.ftol = 1.0e-10;
            lm.parameters.xtol = 1.0e-10;

            LOG_TRACE_L2(m_logger, "Minimizing functor...");
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(v_initial);
            LOG_TRACE_L2(m_logger, "Minimizing functor status: {}", lm_status_map.at(status));

            if (status <= 0 || status > 4) {
                std::stringstream msg;
                msg << "While working on QSE group with algebraic species: ";
                size_t count = 0;
                for (const auto& species: algebraic_species) {
                    msg << species;
                    if (count < algebraic_species.size() - 1) {
                        msg << ", ";
                    }
                    count++;
                }
                msg << " the QSE solver failed to converge with status: " << lm_status_map.at(status) << " (No. " << status << ").";
                LOG_ERROR(m_logger, "{}", msg.str());
                throw std::runtime_error(msg.str());
            }
            LOG_TRACE_L1(m_logger, "QSE Group minimization succeeded with status: {}", lm_status_map.at(status));
            Eigen::VectorXd Y_final_qse = v_initial.array().exp(); // Convert back to physical abundances using exponential scaling
            i = 0;
            for (const auto& species: algebraic_species) {
                LOG_TRACE_L1(
                    m_logger,
                    "During QSE solving species {} started with a molar abundance of {} and ended with an abundance of {}.",
                    species.name(),
                    comp.getMolarAbundance(species),
                    Y_final_qse(i)
                );
                // double Xi = Y_final_qse(i) * species.mass(); // Convert from molar abundance to mass fraction
                if (!outputComposition.contains (species)) {
                    outputComposition.registerSpecies(species);
                }
                outputComposition.setMolarAbundance(species, Y_final_qse(i));
                i++;
            }
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


     //////////////////////////////////////
    /// Eigen Functor Member Functions ///
   /////////////////////////////////////


    int MultiscalePartitioningEngineView::EigenFunctor::operator()(const InputType &v_qse, OutputType &f_qse) const {
        fourdst::composition::Composition comp_trial(m_initial_comp.getRegisteredSpecies());
        for (const auto& [sp, y] : m_initial_comp) {
            comp_trial.setMolarAbundance(sp, y);
        }
        Eigen::VectorXd y_qse = v_qse.array().exp(); // Convert to physical abundances using exponential scaling

        for (const auto& species: m_qse_solve_species) {
            auto index = static_cast<long>(m_qse_solve_species_index_map.at(species));
            comp_trial.setMolarAbundance(species, y_qse[index]);
        }

        const auto result = m_view.getBaseEngine().calculateRHSAndEnergy(comp_trial, m_T9, m_rho);
        if (!result) {
            throw exceptions::StaleEngineError("Failed to calculate RHS and energy due to stale engine state");
        }
        const auto&[dydt, nuclearEnergyGenerationRate, _] = result.value();
        f_qse.resize(static_cast<long>(m_qse_solve_species.size()));
        long i = 0;
        // TODO: make sure that just counting up i is a valid approach, this is a possible place an indexing bug may have crept in
        for (const auto& species: m_qse_solve_species) {
            const double dydt_i = dydt.at(species);
            f_qse(i) = dydt_i/y_qse(i); // We square the residuals to improve numerical stability in the solver
            i++;
        }
        LOG_TRACE_L2(
            m_view.m_logger,
            "Functor evaluation at T9 = {}, rho = {}, y_qse (v_qse) = <{}> complete. ||f|| = {}",
            m_T9,
            m_rho,
            [&]() -> std::string {
                std::stringstream ss;
                for (long j = 0; j < y_qse.size(); ++j) {
                    ss << y_qse(j);
                    ss << "(" << v_qse(j) << ")";
                    if (j < y_qse.size() - 1) {
                        ss << ", ";
                    }
                }
                return ss.str();
            }(),
            f_qse.norm()
        );
        LOG_TRACE_L3(
            m_view.m_logger,
            "{}",
            [&]() -> std::string {
                std::stringstream ss;
                const std::vector species(m_qse_solve_species.begin(), m_qse_solve_species.end());
                for (long j = 0; j < f_qse.size(); ++j) {
                    ss << "Residual for species " << species.at(j).name() << " f(" << j << ") = " << f_qse(j) << "\n";
                }
                return ss.str();
            }()
        );
        return 0;
    }

    int MultiscalePartitioningEngineView::EigenFunctor::df(const InputType &v_qse, JacobianType &J_qse) const {
        fourdst::composition::Composition comp_trial(m_initial_comp.getRegisteredSpecies());
        for (const auto& [sp, y] : m_initial_comp) {
            comp_trial.setMolarAbundance(sp, y);
        }
        Eigen::VectorXd y_qse = v_qse.array().exp(); // Convert to physical abundances using exponential scaling

        for (const auto& species: m_qse_solve_species) {
            const double molarAbundance = y_qse[static_cast<long>(m_qse_solve_species_index_map.at(species))];
            comp_trial.setMolarAbundance(species, molarAbundance);
        }

        std::vector<Species> qse_species_vector(m_qse_solve_species.begin(), m_qse_solve_species.end());
        NetworkJacobian jac = m_view.getBaseEngine().generateJacobianMatrix(comp_trial, m_T9, m_rho, qse_species_vector);
        const auto result = m_view.getBaseEngine().calculateRHSAndEnergy(comp_trial, m_T9, m_rho);
        if (!result) {
            throw exceptions::StaleEngineError("Failed to calculate RHS and energy due to stale engine state");
        }
        const auto&[dydt, nuclearEnergyGenerationRate, _] = result.value();

        const long N = static_cast<long>(m_qse_solve_species.size());
        J_qse.resize(N, N);
        long rowID = 0;
        for (const auto& rowSpecies : m_qse_solve_species) {
            long colID = 0;
            for (const auto& colSpecies: m_qse_solve_species) {
                J_qse(rowID, colID) = jac(rowSpecies, colSpecies);
                colID += 1;
                LOG_TRACE_L3(m_view.m_logger, "Jacobian[{}, {}] (d(dY({}))/dY({})) = {}", rowID, colID - 1, rowSpecies.name(), colSpecies.name(), J_qse(rowID, colID - 1));
            }
            rowID += 1;
        }

        for (long i = 0; i < J_qse.rows(); ++i) {
            for (long j = 0; j < J_qse.cols(); ++j) {
                double on_diag_correction = 0.0;
                if (i == j) {
                    auto rowSpecies = *(std::next(m_qse_solve_species.begin(), i));
                    const double Fi = dydt.at(rowSpecies);
                    on_diag_correction = Fi / y_qse(i);
                }
                J_qse(i, j) = y_qse(j) * (J_qse(i, j) - on_diag_correction) / y_qse(i); // Apply chain rule J'(i,j) = y_j * (J(i,j) - δ_ij(F_i/Y_i)) / Y_i
            }
        }

        return 0; // Success
    }



      /////////////////////////////////
     /// QSEGroup Member Functions ///
    ////////////////////////////////



    bool MultiscalePartitioningEngineView::QSEGroup::operator==(const QSEGroup &other) const {
       return mean_timescale == other.mean_timescale;
    }

    void MultiscalePartitioningEngineView::QSEGroup::removeSpecies(const Species &species) {
        if (algebraic_species.contains(species)) {
            algebraic_species.erase(species);
        }
        if (seed_species.contains(species)) {
            seed_species.erase(species);
        }
    }

    void MultiscalePartitioningEngineView::QSEGroup::addSpeciesToAlgebraic(const Species &species) {
        if (seed_species.contains(species)) {
            const std::string msg = std::format("Cannot add species {} to algebraic set as it is already in the seed set.", species.name());
            throw std::invalid_argument(msg);
        }
        if (!algebraic_species.contains(species)) {
            algebraic_species.insert(species);
        }
    }

    void MultiscalePartitioningEngineView::QSEGroup::addSpeciesToSeed(const Species &species) {
        if (algebraic_species.contains(species)) {
            const std::string msg = std::format("Cannot add species {} to seed set as it is already in the algebraic set.", species.name());
            throw std::invalid_argument(msg);
        }
        if (!seed_species.contains(species)) {
            seed_species.insert(species);
        }
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
}
