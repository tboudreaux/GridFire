#include "gridfire/engine/views/engine_defined.h"

#include <ranges>

#include "quill/LogMacros.h"

#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace gridfire {
    using fourdst::atomic::Species;

    DefinedEngineView::DefinedEngineView(const std::vector<std::string>& peNames, DynamicEngine& baseEngine) :
    m_baseEngine(baseEngine) {
        collect(peNames);
    }

    const DynamicEngine & DefinedEngineView::getBaseEngine() const {
        return m_baseEngine;
    }

    const std::vector<Species> & DefinedEngineView::getNetworkSpecies() const {
        return m_activeSpecies;
    }

    std::expected<StepDerivatives<double>, expectations::StaleEngineError> DefinedEngineView::calculateRHSAndEnergy(
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_defined);
        const auto result = m_baseEngine.calculateRHSAndEnergy(Y_full, T9, rho);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const auto [dydt, nuclearEnergyGenerationRate] = result.value();
        StepDerivatives<double> definedResults;
        definedResults.nuclearEnergyGenerationRate = nuclearEnergyGenerationRate;
        definedResults.dydt = mapFullToView(dydt);
        return definedResults;
    }

    void DefinedEngineView::generateJacobianMatrix(
        const std::vector<double> &Y_dynamic,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_dynamic);
        m_baseEngine.generateJacobianMatrix(Y_full, T9, rho);
    }

    double DefinedEngineView::getJacobianMatrixEntry(
        const int i_defined,
        const int j_defined
    ) const {
        validateNetworkState();

        const size_t i_full = mapViewToFullSpeciesIndex(i_defined);
        const size_t j_full = mapViewToFullSpeciesIndex(j_defined);

        return m_baseEngine.getJacobianMatrixEntry(i_full, j_full);
    }

    void DefinedEngineView::generateStoichiometryMatrix() {
        validateNetworkState();

        m_baseEngine.generateStoichiometryMatrix();
    }

    int DefinedEngineView::getStoichiometryMatrixEntry(
        const int speciesIndex_defined,
        const int reactionIndex_defined
    ) const {
        validateNetworkState();

        const size_t i_full = mapViewToFullSpeciesIndex(speciesIndex_defined);
        const size_t j_full = mapViewToFullReactionIndex(reactionIndex_defined);
        return m_baseEngine.getStoichiometryMatrixEntry(i_full, j_full);
    }

    double DefinedEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        if (!m_activeReactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the DefinedEngineView.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }
        const auto Y_full = mapViewToFull(Y_defined);
        return m_baseEngine.calculateMolarReactionFlow(reaction, Y_full, T9, rho);
    }

    const reaction::LogicalReactionSet & DefinedEngineView::getNetworkReactions() const {
        validateNetworkState();

        return m_activeReactions;
    }

    void DefinedEngineView::setNetworkReactions(const reaction::LogicalReactionSet &reactions) {
        std::vector<std::string> peNames;
        for (const auto& reaction : reactions) {
            peNames.push_back(std::string(reaction.id()));
        }
        collect(peNames);
    }

    std::expected<std::unordered_map<Species, double>, expectations::StaleEngineError> DefinedEngineView::getSpeciesTimescales(
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_defined);
        const auto result = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);
        if (!result) {
            return std::unexpected{result.error()};
        }
        const auto& fullTimescales = result.value();

        std::unordered_map<Species, double> definedTimescales;
        for (const auto& active_species : m_activeSpecies) {
            if (fullTimescales.contains(active_species)) {
                definedTimescales[active_species] = fullTimescales.at(active_species);
            }
        }
        return definedTimescales;
    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError>
    DefinedEngineView::getSpeciesDestructionTimescales(
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_defined);
        const auto result = m_baseEngine.getSpeciesDestructionTimescales(Y_full, T9, rho);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const auto& destructionTimescales = result.value();

        std::unordered_map<Species, double> definedTimescales;
        for (const auto& active_species : m_activeSpecies) {
            if (destructionTimescales.contains(active_species)) {
                definedTimescales[active_species] = destructionTimescales.at(active_species);
            }
        }
        return definedTimescales;
    }

    fourdst::composition::Composition DefinedEngineView::update(const NetIn &netIn) {
        return m_baseEngine.update(netIn);
    }

    bool DefinedEngineView::isStale(const NetIn &netIn) {
        return m_baseEngine.isStale(netIn);
    }


    void DefinedEngineView::setScreeningModel(const screening::ScreeningType model) {
        m_baseEngine.setScreeningModel(model);
    }

    screening::ScreeningType DefinedEngineView::getScreeningModel() const {
        return m_baseEngine.getScreeningModel();
    }

    int DefinedEngineView::getSpeciesIndex(const Species &species) const {
        validateNetworkState();

        const auto it = std::ranges::find(m_activeSpecies, species);
        if (it != m_activeSpecies.end()) {
            return static_cast<int>(std::distance(m_activeSpecies.begin(), it));
        } else {
            LOG_ERROR(m_logger, "Species '{}' not found in active species list.", species.name());
            m_logger->flush_log();
            throw std::runtime_error("Species not found in active species list: " + std::string(species.name()));
        }
    }

    std::vector<double> DefinedEngineView::mapNetInToMolarAbundanceVector(const NetIn &netIn) const {
        std::vector<double> Y(m_activeSpecies.size(), 0.0); // Initialize with zeros
        for (const auto& [symbol, entry] : netIn.composition) {
            auto it = std::ranges::find(m_activeSpecies, entry.isotope());
            if (it != m_activeSpecies.end()) {
                Y[getSpeciesIndex(entry.isotope())] = netIn.composition.getMolarAbundance(symbol); // Map species to their molar abundance
            }
        }
        return Y; // Return the vector of molar abundances
    }

    PrimingReport DefinedEngineView::primeEngine(const NetIn &netIn) {
        return m_baseEngine.primeEngine(netIn);
    }

    std::vector<size_t> DefinedEngineView::constructSpeciesIndexMap() const {
        LOG_TRACE_L3(m_logger, "Constructing species index map for DefinedEngineView...");
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
        LOG_TRACE_L3(m_logger, "Species index map constructed with {} entries.", speciesIndexMap.size());
        return speciesIndexMap;

    }

    std::vector<size_t> DefinedEngineView::constructReactionIndexMap() const {
        LOG_TRACE_L3(m_logger, "Constructing reaction index map for DefinedEngineView...");

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

        LOG_TRACE_L3(m_logger, "Reaction index map constructed with {} entries.", reactionIndexMap.size());
        return reactionIndexMap;
    }

    std::vector<double> DefinedEngineView::mapViewToFull(const std::vector<double>& culled) const {
        std::vector<double> full(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (size_t i_culled = 0; i_culled < culled.size(); ++i_culled) {
            const size_t i_full = m_speciesIndexMap[i_culled];
            full[i_full] += culled[i_culled];
        }
        return full;
    }

    std::vector<double> DefinedEngineView::mapFullToView(const std::vector<double>& full) const {
        std::vector<double> culled(m_activeSpecies.size(), 0.0);
        for (size_t i_culled = 0; i_culled < m_activeSpecies.size(); ++i_culled) {
            const size_t i_full = m_speciesIndexMap[i_culled];
            culled[i_culled] = full[i_full];
        }
        return culled;
    }

    size_t DefinedEngineView::mapViewToFullSpeciesIndex(size_t culledSpeciesIndex) const {
        if (culledSpeciesIndex < 0 || culledSpeciesIndex >= m_speciesIndexMap.size()) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for species index map of size {}.", culledSpeciesIndex, m_speciesIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledSpeciesIndex) + " is out of bounds for species index map of size " + std::to_string(m_speciesIndexMap.size()) + ".");
        }
        return m_speciesIndexMap[culledSpeciesIndex];
    }

    size_t DefinedEngineView::mapViewToFullReactionIndex(size_t culledReactionIndex) const {
        if (culledReactionIndex < 0 || culledReactionIndex >= m_reactionIndexMap.size()) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for reaction index map of size {}.", culledReactionIndex, m_reactionIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledReactionIndex) + " is out of bounds for reaction index map of size " + std::to_string(m_reactionIndexMap.size()) + ".");
        }
        return m_reactionIndexMap[culledReactionIndex];
    }

    void DefinedEngineView::validateNetworkState() const {
        if (m_isStale) {
            LOG_ERROR(m_logger, "DefinedEngineView is stale. Please call update() with a valid NetIn object.");
            m_logger->flush_log();
            throw std::runtime_error("DefinedEngineView is stale. Please call update() with a valid NetIn object.");
        }
    }

    void DefinedEngineView::collect(const std::vector<std::string> &peNames) {
        std::unordered_set<Species> seenSpecies;

        const auto& fullNetworkReactionSet = m_baseEngine.getNetworkReactions();
        for (const auto& peName : peNames) {
            if (!fullNetworkReactionSet.contains(peName)) {
                LOG_ERROR(m_logger, "Reaction with name '{}' not found in the base engine's network reactions. Aborting...", peName);
                m_logger->flush_log();
                throw std::runtime_error("Reaction with name '" + std::string(peName) + "' not found in the base engine's network reactions.");
            }
            auto reaction = fullNetworkReactionSet[peName];
            for (const auto& reactant : reaction.reactants()) {
                if (!seenSpecies.contains(reactant)) {
                    seenSpecies.insert(reactant);
                    m_activeSpecies.push_back(reactant);
                }
            }
            for (const auto& product : reaction.products()) {
                if (!seenSpecies.contains(product)) {
                    seenSpecies.insert(product);
                    m_activeSpecies.push_back(product);
                }
            }
            m_activeReactions.add_reaction(reaction);
        }
        LOG_TRACE_L3(m_logger, "DefinedEngineView built with {} active species and {} active reactions.", m_activeSpecies.size(), m_activeReactions.size());
        LOG_TRACE_L3(m_logger, "Active species: {}", [this]() -> std::string {
            std::string result;
            for (const auto& species : m_activeSpecies) {
                result += std::string(species.name()) + ", ";
            }
            if (!result.empty()) {
                result.pop_back(); // Remove last space
                result.pop_back(); // Remove last comma
            }
            return result;
        }());
        LOG_TRACE_L3(m_logger, "Active reactions: {}", [this]() -> std::string {
            std::string result;
            for (const auto& reaction : m_activeReactions) {
                result += std::string(reaction.id()) + ", ";
            }
            if (!result.empty()) {
                result.pop_back(); // Remove last space
                result.pop_back(); // Remove last comma
            }
            return result;
        }());
        m_speciesIndexMap = constructSpeciesIndexMap();
        m_reactionIndexMap = constructReactionIndexMap();
        m_isStale = false;
    }


    ////////////////////////////////////////////
    /// FileDefinedEngineView Implementation ///
    /////////////////////////////////////////////

    FileDefinedEngineView::FileDefinedEngineView(
        DynamicEngine &baseEngine,
        const std::string &fileName,
        const io::NetworkFileParser &parser
    ):
    DefinedEngineView(parser.parse(fileName), baseEngine),
    m_fileName(fileName),
    m_parser(parser) {}
}
