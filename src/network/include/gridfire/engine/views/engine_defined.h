#pragma once

#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/io/network_file.h"
#include "gridfire/network.h"

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "quill/Logger.h"

#include <string>

namespace gridfire{
    class DefinedEngineView : public DynamicEngine, public EngineView<DynamicEngine> {
    public:
        DefinedEngineView(const std::vector<std::string>& peNames, DynamicEngine& baseEngine);
        const DynamicEngine& getBaseEngine() const override;

        // --- Engine Interface ---
        /**
         * @brief Gets the list of active species in the network defined by the file.
         * @return A const reference to the vector of active species.
         */
        const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

        // --- DynamicEngine Interface ---
        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation for the active species.
         *
         * @param Y_defined A vector of abundances for the active species.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @return A StepDerivatives struct containing the derivatives of the active species and the
         *         nuclear energy generation rate.
         *
         * @throws std::runtime_error If the view is stale (i.e., `update()` has not been called after `setNetworkFile()`).
         */
        std::expected<StepDerivatives<double>, expectations::StaleEngineError> calculateRHSAndEnergy(
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;
        /**
         * @brief Generates the Jacobian matrix for the active species.
         *
         * @param Y_dynamic A vector of abundances for the active species.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         *
         * @throws std::runtime_error If the view is stale.
         */
        void generateJacobianMatrix(
            const std::vector<double>& Y_dynamic,
            const double T9,
            const double rho
        ) const override;
        /**
         * @brief Gets an entry from the Jacobian matrix for the active species.
         *
         * @param i_defined The row index (species index) in the defined matrix.
         * @param j_defined The column index (species index) in the defined matrix.
         * @return The value of the Jacobian matrix at (i_defined, j_defined).
         *
         * @throws std::runtime_error If the view is stale.
         * @throws std::out_of_range If an index is out of bounds.
         */
        double getJacobianMatrixEntry(
            const int i_defined,
            const int j_defined
        ) const override;
        /**
         * @brief Generates the stoichiometry matrix for the active reactions and species.
         *
         * @throws std::runtime_error If the view is stale.
         */
        void generateStoichiometryMatrix() override;
        /**
         * @brief Gets an entry from the stoichiometry matrix for the active species and reactions.
         *
         * @param speciesIndex_defined The index of the species in the defined species list.
         * @param reactionIndex_defined The index of the reaction in the defined reaction list.
         * @return The stoichiometric coefficient for the given species and reaction.
         *
         * @throws std::runtime_error If the view is stale.
         * @throws std::out_of_range If an index is out of bounds.
         */
        int getStoichiometryMatrixEntry(
            const int speciesIndex_defined,
            const int reactionIndex_defined
        ) const override;
        /**
         * @brief Calculates the molar reaction flow for a given reaction in the active network.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param Y_defined Vector of current abundances for the active species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * @throws std::runtime_error If the view is stale or if the reaction is not in the active set.
         */
        double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;
        /**
         * @brief Gets the set of active logical reactions in the network.
         *
         * @return Reference to the LogicalReactionSet containing all active reactions.
         *
         * @throws std::runtime_error If the view is stale.
         */
        const reaction::LogicalReactionSet& getNetworkReactions() const override;

        void setNetworkReactions(const reaction::LogicalReactionSet& reactions) override;
        /**
         * @brief Computes timescales for all active species in the network.
         *
         * @param Y_defined Vector of current abundances for the active species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;

        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const std::vector<double>& Y_defined,
            const double T9,
            const double rho
        ) const override;

        /**
         * @brief Updates the engine view if it is marked as stale.
         *
         * This method checks if the view is stale (e.g., after `setNetworkFile` was called).
         * If it is, it rebuilds the active network from the currently set file.
         * The `netIn` parameter is not used by this implementation but is required by the interface.
         *
         * @param netIn The current network input (unused).
         *
         * @post If the view was stale, it is rebuilt and is no longer stale.
         */
        fourdst::composition::Composition update(const NetIn &netIn) override;

        bool isStale(const NetIn& netIn) override;

        /**
         * @brief Sets the screening model for the base engine.
         *
         * @param model The screening model to set.
         */
        void setScreeningModel(screening::ScreeningType model) override;

        /**
         * @brief Gets the screening model from the base engine.
         *
         * @return The current screening model type.
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        [[nodiscard]] int getSpeciesIndex(const fourdst::atomic::Species &species) const override;

        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

        [[nodiscard]] PrimingReport primeEngine(const NetIn &netIn) override;
    protected:
        bool m_isStale = true;
        DynamicEngine& m_baseEngine;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log"); ///< Logger instance for trace and debug information.
        ///< Active species in the defined engine.
        std::vector<fourdst::atomic::Species> m_activeSpecies;
        ///< Active reactions in the defined engine.
        reaction::LogicalReactionSet m_activeReactions;

        ///< Maps indices of active species to indices in the full network.
        std::vector<size_t> m_speciesIndexMap;
        ///< Maps indices of active reactions to indices in the full network.
        std::vector<size_t> m_reactionIndexMap;
    private:
        /**
         * @brief Constructs the species index map.
         *
         * @return A vector mapping defined species indices to full species indices.
         *
         * This method creates a map from the indices of the active species to the indices of the
         * corresponding species in the full network.
         *
         * @throws std::runtime_error If an active species is not found in the base engine's species list.
         */
        std::vector<size_t> constructSpeciesIndexMap() const;

        /**
         * @brief Constructs the reaction index map.
         *
         * @return A vector mapping defined reaction indices to full reaction indices.
         *
         * This method creates a map from the indices of the active reactions to the indices of the
         * corresponding reactions in the full network.
         *
         * @throws std::runtime_error If an active reaction is not found in the base engine's reaction list.
         */
        std::vector<size_t> constructReactionIndexMap() const;

        /**
         * @brief Maps a vector of culled abundances to a vector of full abundances.
         *
         * @param defined A vector of abundances for the active species.
         * @return A vector of abundances for the full network, with the abundances of the active
         *         species copied from the defined vector.
         */
        std::vector<double> mapViewToFull(const std::vector<double>& defined) const;

        /**
         * @brief Maps a vector of full abundances to a vector of culled abundances.
         *
         * @param full A vector of abundances for the full network.
         * @return A vector of abundances for the active species, with the abundances of the active
         *         species copied from the full vector.
         */
        std::vector<double> mapFullToView(const std::vector<double>& full) const;

        /**
         * @brief Maps a culled species index to a full species index.
         *
         * @param definedSpeciesIndex The index of the species in the defined species list.
         * @return The index of the corresponding species in the full network.
         *
         * @throws std::out_of_range If the defined index is out of bounds for the species index map.
         */
        size_t mapViewToFullSpeciesIndex(size_t definedSpeciesIndex) const;

        /**
         * @brief Maps a culled reaction index to a full reaction index.
         *
         * @param definedReactionIndex The index of the reaction in the defined reaction list.
         * @return The index of the corresponding reaction in the full network.
         *
         * @throws std::out_of_range If the defined index is out of bounds for the reaction index map.
         */
        size_t mapViewToFullReactionIndex(size_t definedReactionIndex) const;

        void validateNetworkState() const;

        void collect(const std::vector<std::string>& peNames);

    };

    class FileDefinedEngineView final: public DefinedEngineView {
    public:
        explicit FileDefinedEngineView(
            DynamicEngine& baseEngine,
            const std::string& fileName,
            const io::NetworkFileParser& parser
        );
        std::string getNetworkFile() const { return m_fileName; }
        const io::NetworkFileParser& getParser() const { return m_parser; }
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        std::string m_fileName;
        ///< Parser for the network file.
        const io::NetworkFileParser& m_parser;
    };
}