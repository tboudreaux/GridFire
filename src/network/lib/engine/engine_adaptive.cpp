#include "gridfire/engine/engine_adaptive.h"

#include <ranges>
#include <queue>

#include "gridfire/network.h"

#include "quill/LogMacros.h"
#include "quill/Logger.h"

namespace gridfire {
    using fourdst::atomic::Species;
    AdaptiveEngineView::AdaptiveEngineView(
        DynamicEngine &baseEngine
    ) :
    m_baseEngine(baseEngine),
    m_activeSpecies(baseEngine.getNetworkSpecies()),
    m_activeReactions(baseEngine.getNetworkReactions()),
    m_speciesIndexMap(constructSpeciesIndexMap()),
    m_reactionIndexMap(constructReactionIndexMap())
    {
    }

    std::vector<size_t> AdaptiveEngineView::constructSpeciesIndexMap() const {
        LOG_TRACE_L1(m_logger, "Constructing species index map for adaptive engine view...");
        std::unordered_map<Species, size_t> fullSpeciesReverseMap;
        const auto& fullSpeciesList = m_baseEngine.getNetworkSpecies();

        fullSpeciesReverseMap.reserve(fullSpeciesList.size());

        for (size_t i = 0; i < fullSpeciesList.size(); ++i) {
            fullSpeciesReverseMap[fullSpeciesList[i]] = i;
        }

        std::vector<size_t> speciesIndexMap;
        speciesIndexMap.reserve(m_activeSpecies.size());

        for (const auto& active_species : m_activeSpecies) {
            auto it = fullSpeciesReverseMap.find(active_species);
            if (it != fullSpeciesReverseMap.end()) {
                speciesIndexMap.push_back(it->second);
            } else {
                LOG_ERROR(m_logger, "Species '{}' not found in full species map.", active_species.name());
                m_logger -> flush_log();
                throw std::runtime_error("Species not found in full species map: " + std::string(active_species.name()));
            }
        }
        LOG_TRACE_L1(m_logger, "Species index map constructed with {} entries.", speciesIndexMap.size());
        return speciesIndexMap;

    }

    std::vector<size_t> AdaptiveEngineView::constructReactionIndexMap() const {
        LOG_TRACE_L1(m_logger, "Constructing reaction index map for adaptive engine view...");

        // --- Step 1: Create a reverse map using the reaction's unique ID as the key. ---
        std::unordered_map<std::string_view, size_t> fullReactionReverseMap;
        const auto& fullReactionSet = m_baseEngine.getNetworkReactions();
        fullReactionReverseMap.reserve(fullReactionSet.size());

        for (size_t i_full = 0; i_full < fullReactionSet.size(); ++i_full) {
            fullReactionReverseMap[fullReactionSet[i_full].id()] = i_full;
        }

        // --- Step 2: Build the final index map using the active reaction set. ---
        std::vector<size_t> reactionIndexMap;
        reactionIndexMap.reserve(m_activeReactions.size());

        for (const auto& active_reaction_ptr : m_activeReactions) {
            auto it = fullReactionReverseMap.find(active_reaction_ptr.id());

            if (it != fullReactionReverseMap.end()) {
                reactionIndexMap.push_back(it->second);
            } else {
                LOG_ERROR(m_logger, "Active reaction '{}' not found in base engine during reaction index map construction.", active_reaction_ptr.id());
                m_logger->flush_log();
                throw std::runtime_error("Mismatch between active reactions and base engine.");
            }
        }

        LOG_TRACE_L1(m_logger, "Reaction index map constructed with {} entries.", reactionIndexMap.size());
        return reactionIndexMap;
    }

    void AdaptiveEngineView::update(const NetIn& netIn) {
        LOG_TRACE_L1(m_logger, "Updating adaptive engine view...");

        const auto& fullSpeciesList = m_baseEngine.getNetworkSpecies();
        std::vector<double>Y_full;
        Y_full.reserve(fullSpeciesList.size());

        for (const auto& species : fullSpeciesList) {
            if (netIn.composition.contains(species)) {
                Y_full.push_back(netIn.composition.getMolarAbundance(std::string(species.name())));
            } else {
                LOG_TRACE_L2(m_logger, "Species '{}' not found in composition. Setting abundance to 0.0.", species.name());
                Y_full.push_back(0.0);
            }
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3


        std::vector<ReactionFlow> reactionFlows;
        const auto& fullReactionSet = m_baseEngine.getNetworkReactions();
        reactionFlows.reserve(fullReactionSet.size());

        for (const auto& reactionPtr : fullReactionSet) {
            const double flow = m_baseEngine.calculateMolarReactionFlow(reactionPtr, Y_full, T9, rho);
            reactionFlows.push_back({&reactionPtr, flow});
        }

        double max_flow = 0.0;
        double min_flow = std::numeric_limits<double>::max();
        double flowSum = 0.0;
        for (const auto&[reactionPtr, flowRate] : reactionFlows) {
            if (flowRate > max_flow) {
                max_flow = flowRate;
            } else if (flowRate < min_flow) {
                min_flow = flowRate;
            }
            flowSum += flowRate;
            LOG_TRACE_L2(m_logger, "Reaction '{}' has flow rate: {:0.3E} [mol/s]", reactionPtr->id(), flowRate);
        }
        flowSum /= reactionFlows.size();

        LOG_DEBUG(m_logger, "Maximum reaction flow rate in adaptive engine view: {:0.3E} [mol/s]", max_flow);
        LOG_DEBUG(m_logger, "Minimum reaction flow rate in adaptive engine view: {:0.3E} [mol/s]", min_flow);
        LOG_DEBUG(m_logger, "Average reaction flow rate in adaptive engine view: {:0.3E} [mol/s]", flowSum);

        const double relative_culling_threshold = m_config.get<double>("gridfire:AdaptiveEngineView:RelativeCullingThreshold", 1e-75);

        double absoluteCullingThreshold = relative_culling_threshold * max_flow;
        LOG_DEBUG(m_logger, "Relative culling threshold: {:0.3E} ({})", relative_culling_threshold, absoluteCullingThreshold);

        // --- Reaction Culling ---
        LOG_TRACE_L1(m_logger, "Culling reactions based on reaction flow rates...");
        std::vector<const reaction::Reaction*> flowCulledReactions;
        for (const auto&[reactionPtr, flowRate] : reactionFlows) {
            if (flowRate > absoluteCullingThreshold) {
                LOG_TRACE_L2(m_logger, "Maintaining reaction '{}' with relative (abs) flow rate: {:0.3E} ({:0.3E} [mol/s])", reactionPtr->id(), flowRate/max_flow, flowRate);
                flowCulledReactions.push_back(reactionPtr);
            }
        }
        LOG_DEBUG(m_logger, "Selected {} (total: {}, culled: {}) reactions based on flow rates.", flowCulledReactions.size(), fullReactionSet.size(), fullReactionSet.size() - flowCulledReactions.size());

        // --- Connectivity Analysis ---
        std::queue<Species> species_to_visit;
        std::unordered_set<Species> reachable_species;

        constexpr double ABUNDANCE_FLOOR = 1e-12; // Abundance floor for a species to be considered part of the initial fuel
        for (const auto& species : fullSpeciesList) {
            if (netIn.composition.contains(species) && netIn.composition.getMassFraction(std::string(species.name())) > ABUNDANCE_FLOOR) {
                species_to_visit.push(species);
                reachable_species.insert(species);
                LOG_TRACE_L2(m_logger, "Species '{}' is part of the initial fuel.", species.name());
            }
        }
        std::unordered_map<Species, std::vector<const reaction::Reaction*>> reactant_to_reactions_map;
        for (const auto* reaction_ptr : flowCulledReactions) {
            for (const auto& reactant : reaction_ptr->reactants()) {
                reactant_to_reactions_map[reactant].push_back(reaction_ptr);
            }
        }

        while (!species_to_visit.empty()) {
            Species currentSpecies = species_to_visit.front();
            species_to_visit.pop();

            auto it = reactant_to_reactions_map.find(currentSpecies);
            if (it == reactant_to_reactions_map.end()) {
                continue; // The species does not initiate any further reactions
            }

            const auto& reactions = it->second;
            for (const auto* reaction_ptr : reactions) {
                for (const auto& product : reaction_ptr->products()) {
                    if (!reachable_species.contains(product)) {
                        reachable_species.insert(product);
                        species_to_visit.push(product);
                        LOG_TRACE_L2(m_logger, "Species '{}' is reachable via reaction '{}'.", product.name(), reaction_ptr->id());
                    }
                }
            }
        }
        LOG_DEBUG(m_logger, "Reachable species count: {}", reachable_species.size());

        m_activeSpecies.assign(reachable_species.begin(), reachable_species.end());
        std::ranges::sort(m_activeSpecies,
                          [](const Species &a, const Species &b) { return a.mass() < b.mass(); });

        m_activeReactions.clear();
        for (const auto* reaction_ptr : flowCulledReactions) {
            bool all_reactants_present = true;
            for (const auto& reactant : reaction_ptr->reactants()) {
                if (!reachable_species.contains(reactant)) {
                    all_reactants_present = false;
                    break;
                }
            }

            if (all_reactants_present) {
                m_activeReactions.add_reaction(*reaction_ptr);
                LOG_TRACE_L2(m_logger, "Maintaining reaction '{}' with all reactants present.", reaction_ptr->id());
            } else {
                LOG_TRACE_L1(m_logger, "Culling reaction '{}' due to missing reactants.", reaction_ptr->id());
            }
        }
        LOG_DEBUG(m_logger, "Active reactions count: {} (total: {}, culled: {}, culled due to connectivity: {})", m_activeReactions.size(),
                  fullReactionSet.size(), fullReactionSet.size() - m_activeReactions.size(), flowCulledReactions.size() - m_activeReactions.size());

        m_speciesIndexMap = constructSpeciesIndexMap();
        m_reactionIndexMap = constructReactionIndexMap();

        m_isStale = false;
    }

    const std::vector<Species> & AdaptiveEngineView::getNetworkSpecies() const {
        return m_activeSpecies;
    }

    StepDerivatives<double> AdaptiveEngineView::calculateRHSAndEnergy(
        const std::vector<double> &Y_culled,
        const double T9,
        const double rho
    ) const {
        validateState();

        const auto Y_full = mapCulledToFull(Y_culled);

        const auto [dydt, nuclearEnergyGenerationRate] = m_baseEngine.calculateRHSAndEnergy(Y_full, T9, rho);

        StepDerivatives<double> culledResults;
        culledResults.nuclearEnergyGenerationRate = nuclearEnergyGenerationRate;
        culledResults.dydt = mapFullToCulled(dydt);

        return culledResults;
    }

    void AdaptiveEngineView::generateJacobianMatrix(
        const std::vector<double> &Y_culled,
        const double T9,
        const double rho
    ) {
        validateState();
        const auto Y_full = mapCulledToFull(Y_culled);

        m_baseEngine.generateJacobianMatrix(Y_full, T9, rho);
    }

    double AdaptiveEngineView::getJacobianMatrixEntry(
        const int i_culled,
        const int j_culled
    ) const {
        validateState();
        const size_t i_full = mapCulledToFullSpeciesIndex(i_culled);
        const size_t j_full = mapCulledToFullSpeciesIndex(j_culled);

        return m_baseEngine.getJacobianMatrixEntry(i_full, j_full);
    }

    void AdaptiveEngineView::generateStoichiometryMatrix() {
        validateState();
        m_baseEngine.generateStoichiometryMatrix();
    }

    int AdaptiveEngineView::getStoichiometryMatrixEntry(
        const int speciesIndex_culled,
        const int reactionIndex_culled
    ) const {
        validateState();
        const size_t speciesIndex_full = mapCulledToFullSpeciesIndex(speciesIndex_culled);
        const size_t reactionIndex_full = mapCulledToFullReactionIndex(reactionIndex_culled);
        return m_baseEngine.getStoichiometryMatrixEntry(speciesIndex_full, reactionIndex_full);
    }

    double AdaptiveEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<double> &Y_culled,
        const double T9,
        const double rho
    ) const {
        validateState();
        if (!m_activeReactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the adaptive engine view.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }
        const auto Y = mapCulledToFull(Y_culled);

        return m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho);
    }

    const reaction::LogicalReactionSet & AdaptiveEngineView::getNetworkReactions() const {
        return m_activeReactions;
    }

    std::unordered_map<Species, double> AdaptiveEngineView::getSpeciesTimescales(
        const std::vector<double> &Y_culled,
        const double T9,
        const double rho
    ) const {
        validateState();
        const auto Y_full = mapCulledToFull(Y_culled);
        const auto fullTimescales = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);

        std::unordered_map<Species, double> culledTimescales;
        culledTimescales.reserve(m_activeSpecies.size());
        for (const auto& active_species : m_activeSpecies) {
            if (fullTimescales.contains(active_species)) {
                culledTimescales[active_species] = fullTimescales.at(active_species);
            }
        }
        return culledTimescales;

    }

    std::vector<double> AdaptiveEngineView::mapCulledToFull(const std::vector<double>& culled) const {
        std::vector<double> full(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (size_t i_culled = 0; i_culled < culled.size(); ++i_culled) {
            const size_t i_full = m_speciesIndexMap[i_culled];
            full[i_full] += culled[i_culled];
        }
        return full;
    }

    std::vector<double> AdaptiveEngineView::mapFullToCulled(const std::vector<double>& full) const {
        std::vector<double> culled(m_activeSpecies.size(), 0.0);
        for (size_t i_culled = 0; i_culled < m_activeSpecies.size(); ++i_culled) {
            const size_t i_full = m_speciesIndexMap[i_culled];
            culled[i_culled] = full[i_full];
        }
        return culled;
    }

    size_t AdaptiveEngineView::mapCulledToFullSpeciesIndex(size_t culledSpeciesIndex) const {
        if (culledSpeciesIndex < 0 || culledSpeciesIndex >= static_cast<int>(m_speciesIndexMap.size())) {
            LOG_ERROR(m_logger, "Culled index {} is out of bounds for species index map of size {}.", culledSpeciesIndex, m_speciesIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Culled index " + std::to_string(culledSpeciesIndex) + " is out of bounds for species index map of size " + std::to_string(m_speciesIndexMap.size()) + ".");
        }
        return m_speciesIndexMap[culledSpeciesIndex];
    }

    size_t AdaptiveEngineView::mapCulledToFullReactionIndex(size_t culledReactionIndex) const {
        if (culledReactionIndex < 0 || culledReactionIndex >= static_cast<int>(m_reactionIndexMap.size())) {
            LOG_ERROR(m_logger, "Culled index {} is out of bounds for reaction index map of size {}.", culledReactionIndex, m_reactionIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Culled index " + std::to_string(culledReactionIndex) + " is out of bounds for reaction index map of size " + std::to_string(m_reactionIndexMap.size()) + ".");
        }
        return m_reactionIndexMap[culledReactionIndex];
    }

    void AdaptiveEngineView::validateState() const {
        if (m_isStale) {
            LOG_ERROR(m_logger, "AdaptiveEngineView is stale. Please call update() before calculating RHS and energy.");
            m_logger->flush_log();
            throw std::runtime_error("AdaptiveEngineView is stale. Please call update() before calculating RHS and energy.");
        }
    }
}

