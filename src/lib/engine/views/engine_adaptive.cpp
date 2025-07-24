#include "gridfire/engine/views/engine_adaptive.h"

#include <ranges>
#include <queue>
#include <algorithm>


#include "gridfire/network.h"
#include "gridfire/exceptions/error_engine.h"

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

    fourdst::composition::Composition AdaptiveEngineView::update(const NetIn &netIn) {
        fourdst::composition::Composition baseUpdatedComposition = m_baseEngine.update(netIn);
        NetIn updatedNetIn = netIn;

        // for (const auto &entry: netIn.composition | std::views::values) {
        //     if (baseUpdatedComposition.contains(entry.isotope())) {
        //         updatedNetIn.composition.setMassFraction(entry.isotope(), baseUpdatedComposition.getMassFraction(entry.isotope()));
        //     }
        // }
        updatedNetIn.composition = baseUpdatedComposition;

        updatedNetIn.composition.finalize(false);

        LOG_TRACE_L1(m_logger, "Updating AdaptiveEngineView with new network input...");

        std::vector<double> Y_Full;
        std::vector<ReactionFlow> allFlows = calculateAllReactionFlows(updatedNetIn, Y_Full);

        double maxFlow = 0.0;

        for (const auto&[reactionPtr, flowRate]: allFlows) {
            if (flowRate > maxFlow) {
                maxFlow = flowRate;
            }
        }
        LOG_DEBUG(m_logger, "Maximum reaction flow rate in adaptive engine view: {:0.3E} [mol/s]", maxFlow);

        const std::unordered_set<Species> reachableSpecies = findReachableSpecies(updatedNetIn);
        LOG_DEBUG(m_logger, "Found {} reachable species in adaptive engine view.", reachableSpecies.size());

        const std::vector<const reaction::LogicalReaction*> finalReactions = cullReactionsByFlow(allFlows, reachableSpecies, Y_Full, maxFlow);

        finalizeActiveSet(finalReactions);

        auto [rescuedReactions, rescuedSpecies] = rescueEdgeSpeciesDestructionChannel(Y_Full, netIn.temperature/1e9, netIn.density, m_activeSpecies, m_activeReactions);

        for (const auto& reactionPtr : rescuedReactions) {
            m_activeReactions.add_reaction(*reactionPtr);
        }

        for (const auto& species : rescuedSpecies) {
            if (!std::ranges::contains(m_activeSpecies, species)) {
                m_activeSpecies.push_back(species);
            }
        }

        m_speciesIndexMap = constructSpeciesIndexMap();
        m_reactionIndexMap = constructReactionIndexMap();

        m_isStale = false;

        LOG_INFO(m_logger, "AdaptiveEngineView updated successfully with {} active species and {} active reactions.", m_activeSpecies.size(), m_activeReactions.size());

        return updatedNetIn.composition;
    }

    bool AdaptiveEngineView::isStale(const NetIn &netIn) {
        return m_isStale || m_baseEngine.isStale(netIn);
    }

    const std::vector<Species> & AdaptiveEngineView::getNetworkSpecies() const {
        return m_activeSpecies;
    }

    std::expected<StepDerivatives<double>, expectations::StaleEngineError> AdaptiveEngineView::calculateRHSAndEnergy(
        const std::vector<double> &Y_culled,
        const double T9,
        const double rho
    ) const {
        validateState();

        const auto Y_full = mapCulledToFull(Y_culled);

        auto result = m_baseEngine.calculateRHSAndEnergy(Y_full, T9, rho);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const auto [dydt, nuclearEnergyGenerationRate] = result.value();
        StepDerivatives<double> culledResults;
        culledResults.nuclearEnergyGenerationRate = nuclearEnergyGenerationRate;
        culledResults.dydt = mapFullToCulled(dydt);

        return culledResults;
    }

    void AdaptiveEngineView::generateJacobianMatrix(
        const std::vector<double> &Y_dynamic,
        const double T9,
        const double rho
    ) const {
        validateState();
        const auto Y_full = mapCulledToFull(Y_dynamic);

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

    void AdaptiveEngineView::setNetworkReactions(const reaction::LogicalReactionSet &reactions) {
        LOG_CRITICAL(m_logger, "AdaptiveEngineView does not support setting network reactions directly. Use update() with NetIn instead. Perhaps you meant to call this on the base engine?");
        throw exceptions::UnableToSetNetworkReactionsError("AdaptiveEngineView does not support setting network reactions directly. Use update() with NetIn instead. Perhaps you meant to call this on the base engine?");
    }

    std::expected<std::unordered_map<Species, double>, expectations::StaleEngineError> AdaptiveEngineView::getSpeciesTimescales(
        const std::vector<double> &Y_culled,
        const double T9,
        const double rho
    ) const {
        validateState();
        const auto Y_full = mapCulledToFull(Y_culled);
        const auto result = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const std::unordered_map<Species, double> fullTimescales = result.value();


        std::unordered_map<Species, double> culledTimescales;
        culledTimescales.reserve(m_activeSpecies.size());
        for (const auto& active_species : m_activeSpecies) {
            if (fullTimescales.contains(active_species)) {
                culledTimescales[active_species] = fullTimescales.at(active_species);
            }
        }
        return culledTimescales;

    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError>
    AdaptiveEngineView::getSpeciesDestructionTimescales(
        const std::vector<double> &Y,
        double T9,
        double rho
    ) const {
        validateState();

        const auto Y_full = mapCulledToFull(Y);
        const auto result = m_baseEngine.getSpeciesDestructionTimescales(Y_full, T9, rho);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const std::unordered_map<Species, double> destructionTimescales = result.value();

        std::unordered_map<Species, double> culledTimescales;
        culledTimescales.reserve(m_activeSpecies.size());
        for (const auto& active_species : m_activeSpecies) {
            if (destructionTimescales.contains(active_species)) {
                culledTimescales[active_species] = destructionTimescales.at(active_species);
            }
        }
        return culledTimescales;
    }

    void AdaptiveEngineView::setScreeningModel(const screening::ScreeningType model) {
        m_baseEngine.setScreeningModel(model);
    }

    screening::ScreeningType AdaptiveEngineView::getScreeningModel() const {
        return m_baseEngine.getScreeningModel();
    }

    std::vector<double> AdaptiveEngineView::mapNetInToMolarAbundanceVector(const NetIn &netIn) const {
        std::vector<double> Y(m_activeSpecies.size(), 0.0); // Initialize with zeros
        for (const auto& [symbol, entry] : netIn.composition) {
            Y[getSpeciesIndex(entry.isotope())] = netIn.composition.getMolarAbundance(symbol); // Map species to their molar abundance
        }
        return Y; // Return the vector of molar abundances
    }

    PrimingReport AdaptiveEngineView::primeEngine(const NetIn &netIn) {
        return m_baseEngine.primeEngine(netIn);
    }

    int AdaptiveEngineView::getSpeciesIndex(const fourdst::atomic::Species &species) const {
        const auto it = std::ranges::find(m_activeSpecies, species);
        if (it != m_activeSpecies.end()) {
            return static_cast<int>(std::distance(m_activeSpecies.begin(), it));
        } else {
            LOG_ERROR(m_logger, "Species '{}' not found in active species list.", species.name());
            m_logger->flush_log();
            throw std::runtime_error("Species not found in active species list: " + std::string(species.name()));
        }
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
        if (culledSpeciesIndex < 0 || culledSpeciesIndex >= m_speciesIndexMap.size()) {
            LOG_ERROR(m_logger, "Culled index {} is out of bounds for species index map of size {}.", culledSpeciesIndex, m_speciesIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Culled index " + std::to_string(culledSpeciesIndex) + " is out of bounds for species index map of size " + std::to_string(m_speciesIndexMap.size()) + ".");
        }
        return m_speciesIndexMap[culledSpeciesIndex];
    }

    size_t AdaptiveEngineView::mapCulledToFullReactionIndex(size_t culledReactionIndex) const {
        if (culledReactionIndex < 0 || culledReactionIndex >= m_reactionIndexMap.size()) {
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

    // TODO: Change this to use a return value instead of an output parameter.
    std::vector<AdaptiveEngineView::ReactionFlow> AdaptiveEngineView::calculateAllReactionFlows(
        const NetIn &netIn,
        std::vector<double> &out_Y_Full
    ) const {
        const auto& fullSpeciesList = m_baseEngine.getNetworkSpecies();
        out_Y_Full.clear();
        out_Y_Full.reserve(fullSpeciesList.size());

        for (const auto& species: fullSpeciesList) {
            if (netIn.composition.contains(species)) {
                out_Y_Full.push_back(netIn.composition.getMolarAbundance(std::string(species.name())));
            } else {
                LOG_TRACE_L2(m_logger, "Species '{}' not found in composition. Setting abundance to 0.0.", species.name());
                out_Y_Full.push_back(0.0);
            }
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        std::vector<ReactionFlow> reactionFlows;
        const auto& fullReactionSet = m_baseEngine.getNetworkReactions();
        reactionFlows.reserve(fullReactionSet.size());
        for (const auto& reaction : fullReactionSet) {
            const double flow = m_baseEngine.calculateMolarReactionFlow(reaction, out_Y_Full, T9, rho);
            reactionFlows.push_back({&reaction, flow});
            LOG_TRACE_L1(m_logger, "Reaction '{}' has flow rate: {:0.3E} [mol/s/g]", reaction.id(), flow);
        }
        return reactionFlows;
    }

    std::unordered_set<Species> AdaptiveEngineView::findReachableSpecies(
        const NetIn &netIn
    ) const {
        std::unordered_set<Species> reachable;
        std::queue<Species> to_vist;

        constexpr double ABUNDANCE_FLOOR = 1e-12; // Abundance floor for a species to be considered part of the initial fuel
        for (const auto& species: m_baseEngine.getNetworkSpecies()) {
            if (netIn.composition.contains(species) && netIn.composition.getMassFraction(std::string(species.name())) > ABUNDANCE_FLOOR) {
                if (!reachable.contains(species)) {
                    to_vist.push(species);
                    reachable.insert(species);
                    LOG_TRACE_L2(m_logger, "Network Connectivity Analysis: Species '{}' is part of the initial fuel.", species.name());
                }
            }
        }

        bool new_species_found_in_pass = true;
        while (new_species_found_in_pass) {
            new_species_found_in_pass = false;
            for (const auto& reaction: m_baseEngine.getNetworkReactions()) {
                bool all_reactants_reachable = true;
                for (const auto& reactant: reaction.reactants()) {
                    if (!reachable.contains(reactant)) {
                        all_reactants_reachable = false;
                        break;
                    }
                }
                if (all_reactants_reachable) {
                    for (const auto& product: reaction.products()) {
                        if (!reachable.contains(product)) {
                            reachable.insert(product);
                            new_species_found_in_pass = true;
                            LOG_TRACE_L2(m_logger, "Network Connectivity Analysis: Species '{}' is reachable via reaction '{}'.", product.name(), reaction.id());
                        }
                    }
                }
            }
        }

        return reachable;
    }

    std::vector<const reaction::LogicalReaction *> AdaptiveEngineView::cullReactionsByFlow(
        const std::vector<ReactionFlow> &allFlows,
        const std::unordered_set<fourdst::atomic::Species> &reachableSpecies,
        const std::vector<double> &Y_full,
        const double maxFlow
    ) const {
        LOG_TRACE_L1(m_logger, "Culling reactions based on flow rates...");
        const auto relative_culling_threshold = m_config.get<double>("gridfire:AdaptiveEngineView:RelativeCullingThreshold", 1e-75);
        double absoluteCullingThreshold = relative_culling_threshold * maxFlow;
        LOG_DEBUG(m_logger, "Relative culling threshold: {:0.3E} ({})", relative_culling_threshold, absoluteCullingThreshold);
        std::vector<const reaction::LogicalReaction*> culledReactions;
        for (const auto& [reactionPtr, flowRate]: allFlows) {
            bool keepReaction = false;
            if (flowRate > absoluteCullingThreshold) {
                LOG_TRACE_L2(m_logger, "Maintaining reaction '{}' with relative (abs) flow rate: {:0.3E} ({:0.3E} [mol/s])", reactionPtr->id(), flowRate/maxFlow, flowRate);
                keepReaction = true;
            } else {
                bool zero_flow_due_to_reachable_reactants = false;
                if (flowRate < 1e-99 && flowRate > 0.0) {
                    for (const auto& reactant: reactionPtr->reactants()) {
                        const auto it = std::ranges::find(m_baseEngine.getNetworkSpecies(), reactant);
                        const size_t index = std::distance(m_baseEngine.getNetworkSpecies().begin(), it);
                        if (Y_full[index] < 1e-99 && reachableSpecies.contains(reactant)) {
                            LOG_TRACE_L1(m_logger, "Maintaining reaction '{}' with low flow ({:0.3E} [mol/s/g]) due to reachable reactant '{}'.", reactionPtr->id(), flowRate, reactant.name());
                            zero_flow_due_to_reachable_reactants = true;
                            break;
                        }
                    }
                }
                if (zero_flow_due_to_reachable_reactants) {
                    keepReaction = true;
                }
            }
            if (keepReaction) {
                culledReactions.push_back(reactionPtr);
            } else {
                LOG_TRACE_L1(m_logger, "Culling reaction '{}' due to low flow rate or lack of connectivity.", reactionPtr->id());
            }
        }
        LOG_DEBUG(m_logger, "Selected {} (total: {}, culled: {}) reactions based on flow rates.", culledReactions.size(), allFlows.size(), allFlows.size() - culledReactions.size());
        return culledReactions;
    }

    AdaptiveEngineView::RescueSet AdaptiveEngineView::rescueEdgeSpeciesDestructionChannel(
        const std::vector<double> &Y_full,
        const double T9,
        const double rho,
        const std::vector<Species> &activeSpecies,
        const reaction::LogicalReactionSet &activeReactions
    ) const {
        const auto result = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);
        if (!result) {
            LOG_ERROR(m_logger, "Failed to get species timescales due to stale engine state.");
            throw exceptions::StaleEngineError("Failed to get species timescales");
        }
        std::unordered_map<Species, double> timescales = result.value();
        std::set<Species> onlyProducedSpecies;
        for (const auto& reaction : activeReactions) {
            const std::vector<Species> products = reaction.products();
            onlyProducedSpecies.insert(products.begin(), products.end());
        }

        // Remove species that are consumed by any one of the active reactions.
        std::erase_if(
            onlyProducedSpecies,
            [&](const Species &species) {
                for (const auto& reaction : activeReactions) {
                    if (reaction.contains_reactant(species)) {
                        return true; // If any active reaction consumes the species then erase it from the set.
                    }
                }
                return false;
            }
        );

        // Remove species that have a non-zero timescale (these are expected to be +inf as they should be the equilibrium species if using with a MultiscalePartitioningEngineView)
        std::erase_if(
            onlyProducedSpecies,
            [&](const Species &species) {
                return std::isinf(timescales.at(species));
            }
        );

        LOG_TRACE_L1(
            m_logger,
            "Found {} species {}that are produced but not consumed in any reaction and do not have an inf timescale (those are expected to be equilibrium species and do not contribute to the stiffness of the network).",
            onlyProducedSpecies.size(),
            [&]() -> std::string {
                std::ostringstream ss;
                if (onlyProducedSpecies.empty()) {
                    return "";
                }
                int count = 0;
                ss << "(";
                for (const auto& species : onlyProducedSpecies) {
                    ss << species.name();
                    if (count < onlyProducedSpecies.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                ss << ") ";
                return ss.str();
            }()
        );

        std::unordered_map<Species, const reaction::LogicalReaction*> reactionsToRescue;
        for (const auto& species : onlyProducedSpecies) {
            double maxSpeciesConsumptionRate = 0.0;
            for (const auto& reaction : m_baseEngine.getNetworkReactions()) {
                const bool speciesToCheckIsConsumed = reaction.contains_reactant(species);
                if (!speciesToCheckIsConsumed) {
                    continue; // If the species is not consumed by this reaction, skip it.
                }
                bool allOtherReactantsAreAvailable = true;
                for (const auto& reactant : reaction.reactants()) {
                    const bool reactantIsAvailable = std::ranges::contains(activeSpecies, reactant);
                    if (!reactantIsAvailable && reactant != species) {
                        allOtherReactantsAreAvailable = false;
                    }
                }
                if (allOtherReactantsAreAvailable && speciesToCheckIsConsumed) {
                    double rate = reaction.calculate_rate(T9);
                    if (rate > maxSpeciesConsumptionRate) {
                        maxSpeciesConsumptionRate = rate;
                        reactionsToRescue[species] = &reaction;
                    }
                }
            }
        }
        LOG_TRACE_L1(
            m_logger,
            "Rescuing {} {}reactions",
            reactionsToRescue.size(),
            [&]() -> std::string {
                std::ostringstream ss;
                if (reactionsToRescue.empty()) {
                    return "";
                }
                int count = 0;
                ss << "(";
                for (const auto &reaction : reactionsToRescue | std::views::values) {
                    ss << reaction->id();
                    if (count < reactionsToRescue.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                ss << ") ";
                return ss.str();
            }()
        );

        LOG_TRACE_L1(
            m_logger,
            "Timescale adjustments due to reaction rescue: {}",
            [&]() -> std::string {
                std::stringstream ss;
                if (reactionsToRescue.empty()) {
                    return "No reactions rescued...";
                }
                int count = 0;
                for (const auto& [species, reaction] : reactionsToRescue) {
                    ss << "(Species: " << species.name() << " started with a timescale of " << timescales.at(species);
                    ss << ", rescued by reaction: " << reaction->id();
                    ss << " whose product timescales are --- [";
                    int icount = 0;
                    for (const auto& product : reaction->products()) {
                        ss << product.name() << ": " << timescales.at(product);
                        if (icount < reaction->products().size() - 1) {
                            ss << ", ";
                        }
                        icount++;
                    }
                    ss << "])";

                    if (count < reactionsToRescue.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }

                return ss.str();
            }()
        );

        RescueSet rescueSet;
        std::unordered_set<const reaction::LogicalReaction*> newReactions;
        std::unordered_set<Species> newSpecies;

        for (const auto &reactionPtr: reactionsToRescue | std::views::values) {
            newReactions.insert(reactionPtr);
            for (const auto& product : reactionPtr->products()) {
                newSpecies.insert(product);
            }
        }
        return {std::move(newReactions), std::move(newSpecies)};
    }

    void AdaptiveEngineView::finalizeActiveSet(
        const std::vector<const reaction::LogicalReaction *> &finalReactions
    ) {
        std::unordered_set<Species>finalSpeciesSet;
        m_activeReactions.clear();
        for (const auto* reactionPtr: finalReactions) {
            m_activeReactions.add_reaction(*reactionPtr);
            for (const auto& reactant : reactionPtr->reactants()) {
                if (!finalSpeciesSet.contains(reactant)) {
                    LOG_TRACE_L1(m_logger, "Adding reactant '{}' to active species set through reaction {}.", reactant.name(), reactionPtr->id());
                } else {
                    LOG_TRACE_L1(m_logger, "Reactant '{}' already in active species set through another reaction.", reactant.name());
                }
                finalSpeciesSet.insert(reactant);
            }
            for (const auto& product : reactionPtr->products()) {
                if (!finalSpeciesSet.contains(product)) {
                    LOG_TRACE_L1(m_logger, "Adding product '{}' to active species set through reaction {}.", product.name(), reactionPtr->id());
                } else {
                    LOG_TRACE_L1(m_logger, "Product '{}' already in active species set through another reaction.", product.name());
                }
                finalSpeciesSet.insert(product);
            }
        }
        m_activeSpecies.clear();
        m_activeSpecies = std::vector<Species>(finalSpeciesSet.begin(), finalSpeciesSet.end());
        std::ranges::sort(
            m_activeSpecies,
            [](const Species &a, const Species &b) { return a.mass() < b.mass(); }
        );
    }
}

