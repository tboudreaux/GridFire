#pragma once

#include "engine_view_abstract.h"
#include "../engine_abstract.h"

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "quill/Logger.h"

#include <string>

namespace gridfire{
    class FileDefinedEngineView final: public DynamicEngine, public EngineView<DynamicEngine> {
    public:
        explicit FileDefinedEngineView(DynamicEngine& baseEngine, const std::string& fileName);

        // --- EngineView Interface ---
        const DynamicEngine& getBaseEngine() const override;

        // --- Engine Interface ---
        const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

        // --- DynamicEngine Interface ---
        StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;
        void generateJacobianMatrix(
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) override;
        double getJacobianMatrixEntry(
            const int i_defined,
            const int j_defined
        ) const override;
        void generateStoichiometryMatrix() override;
        int getStoichiometryMatrixEntry(
            const int speciesIndex_defined,
            const int reactionIndex_defined
        ) const override;
        double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;
        const reaction::LogicalReactionSet& getNetworkReactions() const override;
        std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        DynamicEngine& m_baseEngine;

        std::vector<fourdst::atomic::Species> m_activeSpecies; ///< Active species in the defined engine.
        reaction::LogicalReactionSet m_activeReactions; ///< Active reactions in the defined engine.

        std::vector<size_t> m_speciesIndexMap; ///< Maps indices of active species to indices in the full network.
        std::vector<size_t> m_reactionIndexMap; ///< Maps indices of active reactions to indices in the full network.
    private:
        void buildFromFile(const std::string& fileName);

        /**
         * @brief Constructs the species index map.
         *
         * @return A vector mapping culled species indices to full species indices.
         *
         * This method creates a map from the indices of the active species to the indices of the
         * corresponding species in the full network.
         *
         * @see AdaptiveEngineView::update()
         */
        std::vector<size_t> constructSpeciesIndexMap() const;

        /**
         * @brief Constructs the reaction index map.
         *
         * @return A vector mapping culled reaction indices to full reaction indices.
         *
         * This method creates a map from the indices of the active reactions to the indices of the
         * corresponding reactions in the full network.
         *
         * @see AdaptiveEngineView::update()
         */
        std::vector<size_t> constructReactionIndexMap() const;

        /**
         * @brief Maps a vector of culled abundances to a vector of full abundances.
         *
         * @param culled A vector of abundances for the active species.
         * @return A vector of abundances for the full network, with the abundances of the active
         *         species copied from the culled vector.
         */
        std::vector<double> mapCulledToFull(const std::vector<double>& culled) const;

        /**
         * @brief Maps a vector of full abundances to a vector of culled abundances.
         *
         * @param full A vector of abundances for the full network.
         * @return A vector of abundances for the active species, with the abundances of the active
         *         species copied from the full vector.
         */
        std::vector<double> mapFullToCulled(const std::vector<double>& full) const;

        /**
         * @brief Maps a culled species index to a full species index.
         *
         * @param culledSpeciesIndex The index of the species in the culled species list.
         * @return The index of the corresponding species in the full network.
         *
         * @throws std::out_of_range If the culled index is out of bounds for the species index map.
         */
        size_t mapCulledToFullSpeciesIndex(size_t culledSpeciesIndex) const;

        /**
         * @brief Maps a culled reaction index to a full reaction index.
         *
         * @param culledReactionIndex The index of the reaction in the culled reaction list.
         * @return The index of the corresponding reaction in the full network.
         *
         * @throws std::out_of_range If the culled index is out of bounds for the reaction index map.
         */
        size_t mapCulledToFullReactionIndex(size_t culledReactionIndex) const;

    };
}