#include "gridfire/engine/engine_culled.h"

#include <ranges>

#include "gridfire/network.h"

#include "quill/LogMacros.h"

namespace gridfire {
    using fourdst::atomic::Species;
    AdaptiveEngineView::AdaptiveEngineView(
        DynamicEngine &baseEngine
    ) :
    m_baseEngine(baseEngine),
    m_activeSpecies(baseEngine.getNetworkSpecies()),
    m_activeReactions(baseEngine.getNetworkReactions()),
    m_speciesIndexMap(constructSpeciesIndexMap())
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
                throw std::runtime_error("Species not found in full species map: " + std::string(active_species.name()));
            }
        }
        LOG_TRACE_L1(m_logger, "Successfully constructed species index map with {} entries.", speciesIndexMap.size());
        return speciesIndexMap;

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
                LOG_DEBUG(m_logger, "Species '{}' not found in composition. Setting abundance to 0.0.", species.name());
                Y_full.push_back(0.0);
            }
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        m_isStale = false;

        std::vector<ReactionFlow> reactionFlows;
        const auto& fullReactionSet = m_baseEngine.getNetworkReactions();
        reactionFlows.reserve(fullReactionSet.size());

        for (const auto& reactionPtr : fullReactionSet) {
            const double flow = m_baseEngine.calculateMolarReactionFlow(*reactionPtr, Y_full, T9, rho);
            reactionFlows.push_back({reactionPtr.get(), flow});
        }

        double max_flow = 0.0;
        for (const auto&[reactionPtr, flowRate] : reactionFlows) {
            if (flowRate > max_flow) {
                max_flow = flowRate;
            }
        }

        LOG_DEBUG(m_logger, "Maximum reaction flow rate in adaptive engine view: {:0.3E} [mol/s]", max_flow);
    }

    const std::vector<Species> & AdaptiveEngineView::getNetworkSpecies() const {
        return m_activeSpecies;
    }

    StepDerivatives<double> AdaptiveEngineView::calculateRHSAndEnergy(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateRHSAndEnergy(Y, T9, rho);
    }

    void AdaptiveEngineView::generateJacobianMatrix(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) {
        m_baseEngine.generateJacobianMatrix(Y, T9, rho);
    }

    double AdaptiveEngineView::getJacobianMatrixEntry(
        const int i,
        const int j
    ) const {
        return m_baseEngine.getJacobianMatrixEntry(i, j);
    }

    void AdaptiveEngineView::generateStoichiometryMatrix() {
        m_baseEngine.generateStoichiometryMatrix();
    }

    int AdaptiveEngineView::getStoichiometryMatrixEntry(
        const int speciesIndex,
        const int reactionIndex
    ) const {
        return m_baseEngine.getStoichiometryMatrixEntry(speciesIndex, reactionIndex);
    }

    double AdaptiveEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateMolarReactionFlow(reaction, Y, T9, rho);
    }

    const reaction::REACLIBLogicalReactionSet & AdaptiveEngineView::getNetworkReactions() const {
        return m_activeReactions;
    }

    std::unordered_map<fourdst::atomic::Species, double> AdaptiveEngineView::getSpeciesTimescales(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        auto timescales = m_baseEngine.getSpeciesTimescales(Y, T9, rho);
        for (const auto &species: timescales | std::views::keys) {
            // remove species that are not in the active species list
            if (std::ranges::find(m_activeSpecies, species) == m_activeSpecies.end()) {
                timescales.erase(species);
            }
        }
        return timescales;
    }
}

