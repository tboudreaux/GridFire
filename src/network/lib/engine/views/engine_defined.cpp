#include "gridfire/engine/views/engine_defined.h"

#include "quill/LogMacros.h"

namespace gridfire {
    using fourdst::atomic::Species;

    FileDefinedEngineView::FileDefinedEngineView(
        DynamicEngine &baseEngine,
        const std::string &fileName,
        const io::NetworkFileParser &parser
    ):
    m_baseEngine(baseEngine),
    m_fileName(fileName),
    m_parser(parser),
    m_activeSpecies(baseEngine.getNetworkSpecies()),
    m_activeReactions(baseEngine.getNetworkReactions()) {
        buildFromFile(fileName);
    }

    const DynamicEngine & FileDefinedEngineView::getBaseEngine() const {
        return m_baseEngine;
    }

    const std::vector<Species> & FileDefinedEngineView::getNetworkSpecies() const {
        return m_activeSpecies;
    }

    StepDerivatives<double> FileDefinedEngineView::calculateRHSAndEnergy(
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_defined);
        const auto [dydt, nuclearEnergyGenerationRate] = m_baseEngine.calculateRHSAndEnergy(Y_full, T9, rho);

        StepDerivatives<double> definedResults;
        definedResults.nuclearEnergyGenerationRate = nuclearEnergyGenerationRate;
        definedResults.dydt = mapFullToView(dydt);
        return definedResults;
    }

    void FileDefinedEngineView::generateJacobianMatrix(
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_defined);
        m_baseEngine.generateJacobianMatrix(Y_full, T9, rho);
    }

    double FileDefinedEngineView::getJacobianMatrixEntry(
        const int i_defined,
        const int j_defined
    ) const {
        validateNetworkState();

        const size_t i_full = mapViewToFullSpeciesIndex(i_defined);
        const size_t j_full = mapViewToFullSpeciesIndex(j_defined);

        return m_baseEngine.getJacobianMatrixEntry(i_full, j_full);
    }

    void FileDefinedEngineView::generateStoichiometryMatrix() {
        validateNetworkState();

        m_baseEngine.generateStoichiometryMatrix();
    }

    int FileDefinedEngineView::getStoichiometryMatrixEntry(
        const int speciesIndex_defined,
        const int reactionIndex_defined
    ) const {
        validateNetworkState();

        const size_t i_full = mapViewToFullSpeciesIndex(speciesIndex_defined);
        const size_t j_full = mapViewToFullReactionIndex(reactionIndex_defined);
        return m_baseEngine.getStoichiometryMatrixEntry(i_full, j_full);
    }

    double FileDefinedEngineView::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        if (!m_activeReactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the file defined engine view.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }
        const auto Y_full = mapViewToFull(Y_defined);
        return m_baseEngine.calculateMolarReactionFlow(reaction, Y_full, T9, rho);
    }

    const reaction::LogicalReactionSet & FileDefinedEngineView::getNetworkReactions() const {
        validateNetworkState();

        return m_activeReactions;
    }

    std::unordered_map<Species, double> FileDefinedEngineView::getSpeciesTimescales(
        const std::vector<double> &Y_defined,
        const double T9,
        const double rho
    ) const {
        validateNetworkState();

        const auto Y_full = mapViewToFull(Y_defined);
        const auto fullTimescales = m_baseEngine.getSpeciesTimescales(Y_full, T9, rho);

        std::unordered_map<Species, double> definedTimescales;
        for (const auto& active_species : m_activeSpecies) {
            if (fullTimescales.contains(active_species)) {
                definedTimescales[active_species] = fullTimescales.at(active_species);
            }
        }
        return definedTimescales;
    }

    void FileDefinedEngineView::update(const NetIn &netIn) {
        if (m_isStale) {
            buildFromFile(m_fileName);
        }
    }

    void FileDefinedEngineView::setNetworkFile(const std::string &fileName) {
        m_isStale = true;
        m_fileName = fileName;
        LOG_DEBUG(m_logger, "File '{}' set to '{}'. FileDefinedNetworkView is now stale! You MUST call update() before you use it!", m_fileName, fileName);
    }

    void FileDefinedEngineView::setScreeningModel(const screening::ScreeningType model) {
        m_baseEngine.setScreeningModel(model);
    }

    screening::ScreeningType FileDefinedEngineView::getScreeningModel() const {
        return m_baseEngine.getScreeningModel();
    }

    std::vector<size_t> FileDefinedEngineView::constructSpeciesIndexMap() const {
        LOG_TRACE_L1(m_logger, "Constructing species index map for file defined engine view...");
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

    std::vector<size_t> FileDefinedEngineView::constructReactionIndexMap() const {
        LOG_TRACE_L1(m_logger, "Constructing reaction index map for file defined engine view...");

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

    void FileDefinedEngineView::buildFromFile(const std::string &fileName) {
        LOG_TRACE_L1(m_logger, "Building file defined engine view from {}...", fileName);
        auto [reactionPENames] = m_parser.parse(fileName);

        m_activeReactions.clear();
        m_activeSpecies.clear();

        std::unordered_set<Species> seenSpecies;

        const auto& fullNetworkReactionSet = m_baseEngine.getNetworkReactions();
        for (const auto& peName : reactionPENames) {
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
        LOG_TRACE_L1(m_logger, "File defined engine view built with {} active species and {} active reactions.", m_activeSpecies.size(), m_activeReactions.size());
        LOG_DEBUG(m_logger, "Active species: {}", [this]() -> std::string {
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
        LOG_DEBUG(m_logger, "Active reactions: {}", [this]() -> std::string {
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

    std::vector<double> FileDefinedEngineView::mapViewToFull(const std::vector<double>& culled) const {
        std::vector<double> full(m_baseEngine.getNetworkSpecies().size(), 0.0);
        for (size_t i_culled = 0; i_culled < culled.size(); ++i_culled) {
            const size_t i_full = m_speciesIndexMap[i_culled];
            full[i_full] += culled[i_culled];
        }
        return full;
    }

    std::vector<double> FileDefinedEngineView::mapFullToView(const std::vector<double>& full) const {
        std::vector<double> culled(m_activeSpecies.size(), 0.0);
        for (size_t i_culled = 0; i_culled < m_activeSpecies.size(); ++i_culled) {
            const size_t i_full = m_speciesIndexMap[i_culled];
            culled[i_culled] = full[i_full];
        }
        return culled;
    }

    size_t FileDefinedEngineView::mapViewToFullSpeciesIndex(size_t culledSpeciesIndex) const {
        if (culledSpeciesIndex < 0 || culledSpeciesIndex >= static_cast<int>(m_speciesIndexMap.size())) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for species index map of size {}.", culledSpeciesIndex, m_speciesIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledSpeciesIndex) + " is out of bounds for species index map of size " + std::to_string(m_speciesIndexMap.size()) + ".");
        }
        return m_speciesIndexMap[culledSpeciesIndex];
    }

    size_t FileDefinedEngineView::mapViewToFullReactionIndex(size_t culledReactionIndex) const {
        if (culledReactionIndex < 0 || culledReactionIndex >= static_cast<int>(m_reactionIndexMap.size())) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for reaction index map of size {}.", culledReactionIndex, m_reactionIndexMap.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledReactionIndex) + " is out of bounds for reaction index map of size " + std::to_string(m_reactionIndexMap.size()) + ".");
        }
        return m_reactionIndexMap[culledReactionIndex];
    }

    void FileDefinedEngineView::validateNetworkState() const {
        if (m_isStale) {
            LOG_ERROR(m_logger, "File defined engine view is stale. Please call update() with a valid NetIn object.");
            m_logger->flush_log();
            throw std::runtime_error("File defined engine view is stale. Please call update() with a valid NetIn object.");
        }
    }
}
