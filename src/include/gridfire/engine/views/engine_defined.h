#pragma once

#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/io/network_file.h"
#include "gridfire/types/types.h"

#include "gridfire/config/config.h"

#include "gridfire/engine/scratchpads/blob.h"

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "quill/Logger.h"

#include <string>

namespace gridfire::engine {
    class DefinedEngineView : public DynamicEngine, public EngineView<DynamicEngine> {
    public:
        DefinedEngineView(const std::vector<std::string>& peNames, GraphEngine& baseEngine);

        /** @brief Get the base engine associated with this defined engine view.
         * @return A const reference to the base DynamicEngine.
         */
        [[nodiscard]] const DynamicEngine& getBaseEngine() const override;

        // --- Engine Interface ---
        /**
         * @brief Gets the list of active species in the network defined by the file.
         * @return A const reference to the vector of active species.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getNetworkSpecies(
            scratch::StateBlob& ctx
        ) const override;

        // --- DynamicEngine Interface ---
        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation for the active species.
         *
         * @param comp A Composition object containing the current composition of the system
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param trust
         * @return A StepDerivatives struct containing the derivatives of the active species and the
         *         nuclear energy generation rate.
         *
         * @throws std::runtime_error If the view is stale (i.e., `update()` has not been called after `setNetworkFile()`).
         */
        [[nodiscard]] std::expected<StepDerivatives<double>, EngineStatus> calculateRHSAndEnergy(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            bool trust
        ) const override;

        [[nodiscard]] EnergyDerivatives calculateEpsDerivatives(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the active species.
         *
         * @param comp A Composition object containing the current composition of the system
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the active species.
         *
         * @param comp A Composition object containing the current composition of the system
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param activeSpecies The vector of active species to include in the Jacobian.
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const std::vector<fourdst::atomic::Species> &activeSpecies
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for a given sparsity pattern
         *
         * @param comp A Composition object containing the current composition of the system
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param sparsityPattern The sparsity pattern to use for the Jacobian matrix.
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const SparsityPattern &sparsityPattern
        ) const override;

        /**
         * @brief Calculates the molar reaction flow for a given reaction in the active network.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param comp A Composition object containing the current composition of the system
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * @throws std::runtime_error If the view is stale or if the reaction is not in the active set.
         */
        [[nodiscard]] double calculateMolarReactionFlow(
            scratch::StateBlob& ctx,
            const reaction::Reaction& reaction,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;
        /**
         * @brief Gets the set of active logical reactions in the network.
         *
         * @return Reference to the LogicalReactionSet containing all active reactions.
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] const reaction::ReactionSet& getNetworkReactions(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Computes timescales for all active species in the network.
         *
         * @param comp A Composition object containing the current composition of the system
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus>getSpeciesTimescales(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes destruction timescales for all active species in the network.
         *
         * @param comp A Composition object containing the current composition of the system
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their destruction timescales (s).
         *
         * @throws std::runtime_error If the view is stale.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesDestructionTimescales(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
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
        fourdst::composition::Composition project(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const override;

        /**
         * @brief Gets the screening model from the base engine.
         *
         * @return The current screening model type.
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel(
            scratch::StateBlob& ctx
        ) const override;

        /** @brief Maps a species from the full network to its index in the defined active network.
         *
         * @param species The species to map.
         * @return The index of the species in the active network.
         *
         * @throws std::runtime_error If the species is not in the active set.
         */
        [[nodiscard]] size_t getSpeciesIndex(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        ) const override;

        /**
         * @brief Prime the engine view for calculations. This will delegate to the base engine.
         * @param netIn The current network input.
         * @return The PrimingReport from the base engine.
         */
        [[nodiscard]] PrimingReport primeEngine(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const override;

        /**
         * @brief Collects a Composition object from the base engine.
         *
         * @param comp The full Composition object.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @return A composition object representing the state of the engine stack and the current view.
         */
        [[nodiscard]] fourdst::composition::Composition collectComposition(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the status of a species in the active network.
         *
         * @param species The species for which to get the status.
         * @return The SpeciesStatus indicating if the species is active, inactive, or not present.
         */
        [[nodiscard]] SpeciesStatus getSpeciesStatus(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        ) const override;

        [[nodiscard]] std::optional<StepDerivatives<double>>getMostRecentRHSCalculation(
            scratch::StateBlob &ctx
        ) const override;
    protected:
        bool m_isStale = true;
        GraphEngine& m_baseEngine;
    private:
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log"); ///< Logger instance for trace and debug information.
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
        [[nodiscard]] std::vector<size_t> constructSpeciesIndexMap(
            scratch::StateBlob& ctx
        ) const;

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
        [[nodiscard]] std::vector<size_t> constructReactionIndexMap(
            scratch::StateBlob& ctx
        ) const;

        /**
         * @brief Maps a vector of culled abundances to a vector of full abundances.
         *
         * @param defined A vector of abundances for the active species.
         * @return A vector of abundances for the full network, with the abundances of the active
         *         species copied from the defined vector.
         */
        [[nodiscard]] std::vector<double> mapViewToFull(
            scratch::StateBlob& ctx,
            const std::vector<double>& defined
        ) const;

        /**
         * @brief Maps a vector of full abundances to a vector of culled abundances.
         *
         * @param full A vector of abundances for the full network.
         * @return A vector of abundances for the active species, with the abundances of the active
         *         species copied from the full vector.
         */
        [[nodiscard]] static std::vector<double> mapFullToView(
            scratch::StateBlob& ctx,
            const std::vector<double>& full
        );

        /**
         * @brief Maps a culled species index to a full species index.
         *
         * @param definedSpeciesIndex The index of the species in the defined species list.
         * @return The index of the corresponding species in the full network.
         *
         * @throws std::out_of_range If the defined index is out of bounds for the species index map.
         */
        [[nodiscard]] size_t mapViewToFullSpeciesIndex(
            scratch::StateBlob& ctx,
            size_t definedSpeciesIndex
        ) const;

        /**
         * @brief Maps a culled reaction index to a full reaction index.
         *
         * @param definedReactionIndex The index of the reaction in the defined reaction list.
         * @return The index of the corresponding reaction in the full network.
         *
         * @throws std::out_of_range If the defined index is out of bounds for the reaction index map.
         */
        [[nodiscard]] size_t mapViewToFullReactionIndex(
            scratch::StateBlob& ctx,
            size_t definedReactionIndex
        ) const;

        void collect(
            scratch::StateBlob& ctx,
            const std::vector<std::string>& peNames
        ) const;

    };

    class FileDefinedEngineView final: public DefinedEngineView {
    public:
        explicit FileDefinedEngineView(
            GraphEngine& baseEngine,
            const std::string& fileName,
            const io::NetworkFileParser& parser
        );
        [[nodiscard]] std::string getNetworkFile() const { return m_fileName; }
        [[nodiscard]] const io::NetworkFileParser& getParser() const { return m_parser; }



    private:
        Config<config::GridFireConfig> m_config;

        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        std::string m_fileName;
        ///< Parser for the network file.
        const io::NetworkFileParser& m_parser;
    };
}
