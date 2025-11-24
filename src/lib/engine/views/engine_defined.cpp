#include "gridfire/engine/views/engine_defined.h"
#include "gridfire/engine/engine_graph.h"

#include "fourdst/atomic/species.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/composition/decorators/composition_masked.h"

#include "quill/LogMacros.h"

#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "fourdst/composition/exceptions/exceptions_composition.h"

namespace gridfire::engine {
    using fourdst::atomic::Species;

    DefinedEngineView::DefinedEngineView(
        const std::vector<std::string>& peNames,
        GraphEngine& baseEngine
    ) :
    m_baseEngine(baseEngine) {
        collect(peNames);
    }

    const DynamicEngine & DefinedEngineView::getBaseEngine() const {
        return m_baseEngine;
    }

    const std::vector<Species> & DefinedEngineView::getNetworkSpecies() const {
        if (m_activeSpeciesVectorCache.has_value()) {
            return m_activeSpeciesVectorCache.value();
        }
        m_activeSpeciesVectorCache = std::vector<Species>(m_activeSpecies.begin(), m_activeSpecies.end());
        return m_activeSpeciesVectorCache.value();
    }

    std::expected<StepDerivatives<double>, EngineStatus> DefinedEngineView::calculateRHSAndEnergy(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);
        const auto result = m_baseEngine.calculateRHSAndEnergy(masked, T9, rho, m_activeReactions);

        if (!result) {
            return std::unexpected{result.error()};
        }

        return result.value();
    }

    EnergyDerivatives DefinedEngineView::calculateEpsDerivatives(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);

        return m_baseEngine.calculateEpsDerivatives(masked, T9, rho, m_activeReactions);
    }

    NetworkJacobian DefinedEngineView::generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();
        if (!m_activeSpeciesVectorCache.has_value()) {
            m_activeSpeciesVectorCache = std::vector<Species>(m_activeSpecies.begin(), m_activeSpecies.end());
        }
        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);
        return m_baseEngine.generateJacobianMatrix(masked, T9, rho, m_activeSpeciesVectorCache.value());
    }

    NetworkJacobian DefinedEngineView::generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::vector<fourdst::atomic::Species> &activeSpecies
    ) const {
        validateNetworkState();

        const std::set<fourdst::atomic::Species> activeSpeciesSet(
            activeSpecies.begin(),
            activeSpecies.end()
        );

        const fourdst::composition::MaskedComposition masked(comp, activeSpeciesSet);
        return m_baseEngine.generateJacobianMatrix(masked, T9, rho, activeSpecies);
    }

    NetworkJacobian DefinedEngineView::generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const SparsityPattern &sparsityPattern
    ) const {
        validateNetworkState();
        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);
        return m_baseEngine.generateJacobianMatrix(masked, T9, rho, sparsityPattern);
    }

    void DefinedEngineView::generateStoichiometryMatrix() {
        validateNetworkState();

        m_baseEngine.generateStoichiometryMatrix();
    }

    int DefinedEngineView::getStoichiometryMatrixEntry(
        const Species& species,
        const reaction::Reaction& reaction
    ) const {
        validateNetworkState();

        if (!m_activeSpecies.contains(species)) {
            LOG_ERROR(m_logger, "Species '{}' is not part of the active species in the DefinedEngineView.", species.name());
            m_logger -> flush_log();
            throw std::runtime_error("Species not found in active species: " + std::string(species.name()));
        }

        if (!m_activeReactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the DefinedEngineView.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }

        return m_baseEngine.getStoichiometryMatrixEntry(species, reaction);
    }

    double DefinedEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        if (!m_activeReactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the DefinedEngineView.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }

        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);
        return m_baseEngine.calculateMolarReactionFlow(reaction, masked, T9, rho);
    }

    const reaction::ReactionSet & DefinedEngineView::getNetworkReactions() const {
        validateNetworkState();

        return m_activeReactions;
    }

    void DefinedEngineView::setNetworkReactions(const reaction::ReactionSet &reactions) {
        std::vector<std::string> peNames;
        for (const auto& reaction : reactions) {
            peNames.emplace_back(reaction->id());
        }
        collect(peNames);
        m_activeSpeciesVectorCache = std::nullopt; // Invalidate species vector cache
    }

    std::expected<std::unordered_map<Species, double>, EngineStatus> DefinedEngineView::getSpeciesTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();
        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);

        const auto result = m_baseEngine.getSpeciesTimescales(masked, T9, rho, m_activeReactions);
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

    std::expected<std::unordered_map<Species, double>, EngineStatus> DefinedEngineView::getSpeciesDestructionTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();
        const fourdst::composition::MaskedComposition masked(comp, m_activeSpecies);

        const auto result = m_baseEngine.getSpeciesDestructionTimescales(masked, T9, rho, m_activeReactions);

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

    size_t DefinedEngineView::getSpeciesIndex(const Species &species) const {
        // TODO: We are working to phase out all of these methods, its probably broken but it also should no longer be used and will be removed soon
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
        for (const auto& [sp, y] : netIn.composition) {
            auto it = std::ranges::find(m_activeSpecies, sp);
            if (it != m_activeSpecies.end()) {
                Y[getSpeciesIndex(sp)] = y; // Map species to their molar abundance
            }
        }
        return Y; // Return the vector of molar abundances
    }

    PrimingReport DefinedEngineView::primeEngine(const NetIn &netIn) {
        return m_baseEngine.primeEngine(netIn);
    }

    fourdst::composition::Composition DefinedEngineView::collectComposition(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        fourdst::composition::Composition result = m_baseEngine.collectComposition(comp, T9, rho);

        for (const auto& species : m_activeSpecies) {
            if (!result.contains(species)) {
                result.registerSpecies(species);
            }
        }
        return result;
    }

    SpeciesStatus DefinedEngineView::getSpeciesStatus(const Species &species) const {
        const SpeciesStatus status = m_baseEngine.getSpeciesStatus(species);
        if (status == SpeciesStatus::ACTIVE && !m_activeSpecies.contains(species)) {
            return SpeciesStatus::INACTIVE_FLOW;
        }
        return status;
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
            auto it = fullReactionReverseMap.find(active_reaction_ptr->id());

            if (it != fullReactionReverseMap.end()) {
                reactionIndexMap.push_back(it->second);
            } else {
                LOG_ERROR(m_logger, "Active reaction '{}' not found in base engine during reaction index map construction.", active_reaction_ptr->id());
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
        if (culledSpeciesIndex >= m_speciesIndexMap.size()) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for species index map of size {}.", culledSpeciesIndex, m_speciesIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledSpeciesIndex) + " is out of bounds for species index map of size " + std::to_string(m_speciesIndexMap.size()) + ".");
        }
        return m_speciesIndexMap[culledSpeciesIndex];
    }

    size_t DefinedEngineView::mapViewToFullReactionIndex(size_t culledReactionIndex) const {
        if (culledReactionIndex >= m_reactionIndexMap.size()) {
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
            const reaction::Reaction* reaction = &fullNetworkReactionSet[peName];
            for (const auto& reactant : reaction->reactants()) {
                if (!seenSpecies.contains(reactant)) {
                    seenSpecies.insert(reactant);
                    m_activeSpecies.emplace(reactant);
                }
            }
            for (const auto& product : reaction->products()) {
                if (!seenSpecies.contains(product)) {
                    seenSpecies.insert(product);
                    m_activeSpecies.emplace(product);
                }
            }
            m_activeReactions.add_reaction(*reaction);
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
                result += std::string(reaction->id()) + ", ";
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
        GraphEngine &baseEngine,
        const std::string &fileName,
        const io::NetworkFileParser &parser
    ):
    DefinedEngineView(parser.parse(fileName), baseEngine),
    m_fileName(fileName),
    m_parser(parser) {}
}
