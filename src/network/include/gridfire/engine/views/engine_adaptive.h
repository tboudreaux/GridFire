#pragma once
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/network.h"

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "gridfire/engine/procedures/priming.h"
#include "gridfire/engine/procedures/construction.h"

#include "quill/Logger.h"

namespace gridfire {
    /**
     * @class AdaptiveEngineView
     * @brief An engine view that dynamically adapts the reaction network based on runtime conditions.
     *
     * This class implements an EngineView that dynamically culls species and reactions from the
     * full reaction network based on their reaction flow rates and connectivity. This allows for
     * efficient simulation of reaction networks by focusing computational effort on the most
     * important species and reactions.
     *
     * The AdaptiveEngineView maintains a subset of "active" species and reactions, and maps
     * between the full network indices and the active subset indices. This allows the base engine
     * to operate on the full network data, while the AdaptiveEngineView provides a reduced view
     * for external clients.
     *
     * The adaptation process is driven by the `update()` method, which performs the following steps:
     *   1. **Reaction Flow Calculation:** Calculates the molar reaction flow rate for each reaction
     *      in the full network based on the current temperature, density, and composition.
     *   2. **Reaction Culling:** Culls reactions with flow rates below a threshold, determined by
     *      a relative culling threshold multiplied by the maximum flow rate.
     *   3. **Connectivity Analysis:** Performs a connectivity analysis to identify species that are
     *      reachable from the initial fuel species through the culled reaction network.
     *   4. **Species Culling:** Culls species that are not reachable from the initial fuel.
     *   5. **Index Map Construction:** Constructs index maps to map between the full network indices
     *      and the active subset indices for species and reactions.
     *
     * @implements DynamicEngine
     * @implements EngineView<DynamicEngine>
     *
     * @see engine_abstract.h
     * @see engine_view_abstract.h
     * @see AdaptiveEngineView::update()
     */
    class AdaptiveEngineView final : public DynamicEngine, public EngineView<DynamicEngine> {
    public:
        /**
         * @brief Constructs an AdaptiveEngineView.
         *
         * @param baseEngine The underlying DynamicEngine to which this view delegates calculations.
         *
         * Initializes the active species and reactions to the full network, and constructs the
         * initial index maps.
         */
        explicit AdaptiveEngineView(DynamicEngine& baseEngine);

        /**
         * @brief Updates the active species and reactions based on the current conditions.
         *
         * @param netIn The current network input, containing temperature, density, and composition.
         *
         * This method performs the reaction flow calculation, reaction culling, connectivity analysis,
         * and index map construction steps described above.
         *
         * The culling thresholds are read from the configuration using the following keys:
         *   - `gridfire:AdaptiveEngineView:RelativeCullingThreshold` (default: 1e-75)
         *
         * @throws std::runtime_error If there is a mismatch between the active reactions and the base engine.
         * @post The active species and reactions are updated, and the index maps are reconstructed.
         * @see AdaptiveEngineView
         * @see AdaptiveEngineView::constructSpeciesIndexMap()
         * @see AdaptiveEngineView::constructReactionIndexMap()
         */
        fourdst::composition::Composition update(const NetIn &netIn) override;

        bool isStale(const NetIn& netIn) override;

        /**
         * @brief Gets the list of active species in the network.
         * @return A const reference to the vector of active species.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation for the active species.
         *
         * @param Y_culled A vector of abundances for the active species.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @return A StepDerivatives struct containing the derivatives of the active species and the
         *         nuclear energy generation rate.
         *
         * This method maps the culled abundances to the full network abundances, calls the base engine
         * to calculate the RHS and energy generation, and then maps the full network derivatives back
         * to the culled derivatives.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @see AdaptiveEngineView::update()
         */
        [[nodiscard]] std::expected<StepDerivatives<double>, expectations::StaleEngineError> calculateRHSAndEnergy(
            const std::vector<double> &Y_culled,
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
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to generate the Jacobian matrix.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @see AdaptiveEngineView::update()
         */
        void generateJacobianMatrix(
            const std::vector<double> &Y_dynamic,
            const double T9,
            const double rho
        ) const override;

        /**
         * @brief Gets an entry from the Jacobian matrix for the active species.
         *
         * @param i_culled The row index (species index) in the culled matrix.
         * @param j_culled The column index (species index) in the culled matrix.
         * @return The value of the Jacobian matrix at (i_culled, j_culled).
         *
         * This method maps the culled indices to the full network indices and calls the base engine
         * to get the Jacobian matrix entry.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @throws std::out_of_range If the culled index is out of bounds for the species index map.
         * @see AdaptiveEngineView::update()
         */
        [[nodiscard]] double getJacobianMatrixEntry(
            const int i_culled,
            const int j_culled
        ) const override;

        /**
         * @brief Generates the stoichiometry matrix for the active reactions and species.
         *
         * This method calls the base engine to generate the stoichiometry matrix.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @note The stoichiometry matrix generated by the base engine is assumed to be consistent with
         *       the active species and reactions in this view.
         */
        void generateStoichiometryMatrix() override;

        /**
         * @brief Gets an entry from the stoichiometry matrix for the active species and reactions.
         *
         * @param speciesIndex_culled The index of the species in the culled species list.
         * @param reactionIndex_culled The index of the reaction in the culled reaction list.
         * @return The stoichiometric coefficient for the given species and reaction.
         *
         * This method maps the culled indices to the full network indices and calls the base engine
         * to get the stoichiometry matrix entry.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @throws std::out_of_range If the culled index is out of bounds for the species or reaction index map.
         * @see AdaptiveEngineView::update()
         */
        [[nodiscard]] int getStoichiometryMatrixEntry(
            const int speciesIndex_culled,
            const int reactionIndex_culled
        ) const override;

        /**
         * @brief Calculates the molar reaction flow for a given reaction in the active network.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param Y_culled Vector of current abundances for the active species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to calculate the molar reaction flow.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @throws std::runtime_error If the reaction is not part of the active reactions in the adaptive engine view.
         */
        [[nodiscard]] double calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<double> &Y_culled,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the set of active logical reactions in the network.
         *
         * @return Reference to the LogicalReactionSet containing all active reactions.
         */
        [[nodiscard]] const reaction::LogicalReactionSet& getNetworkReactions() const override;

        void setNetworkReactions(const reaction::LogicalReactionSet& reactions) override;

        /**
         * @brief Computes timescales for all active species in the network.
         *
         * @param Y_culled Vector of current abundances for the active species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to compute the species timescales.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const std::vector<double> &Y_culled,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the base engine.
         * @return A const reference to the base engine.
         */
        [[nodiscard]] const DynamicEngine& getBaseEngine() const override { return m_baseEngine; }

        /**
         * @brief Sets the screening model for the base engine.
         *
         * This method delegates the call to the base engine to set the electron screening model.
         *
         * @param model The electron screening model to set.
         *
         * @par Usage Example:
         * @code
         * AdaptiveEngineView engineView(...);
         * engineView.setScreeningModel(screening::ScreeningType::WEAK);
         * @endcode
         *
         * @post The screening model of the base engine is updated.
         */
        void setScreeningModel(screening::ScreeningType model) override;

        /**
         * @brief Gets the screening model from the base engine.
         *
         * This method delegates the call to the base engine to get the screening model.
         *
         * @return The current screening model type.
         *
         * @par Usage Example:
         * @code
         * AdaptiveEngineView engineView(...);
         * screening::ScreeningType model = engineView.getScreeningModel();
         * @endcode
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        [[nodiscard]] int getSpeciesIndex(const fourdst::atomic::Species &species) const override;

        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

        [[nodiscard]] PrimingReport primeEngine(const NetIn &netIn) override;
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;
        /** @brief A reference to the singleton Config instance, used for retrieving configuration parameters. */
        Config& m_config = Config::getInstance();
        /** @brief A pointer to the logger instance, used for logging messages. */
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        /** @brief The underlying engine to which this view delegates calculations. */
        DynamicEngine& m_baseEngine;

        /** @brief The set of species that are currently active in the network. */
        std::vector<fourdst::atomic::Species> m_activeSpecies;
        /** @brief The set of reactions that are currently active in the network. */
        reaction::LogicalReactionSet m_activeReactions;

        /** @brief A map from the indices of the active species to the indices of the corresponding species in the full network. */
        std::vector<size_t> m_speciesIndexMap;
        /** @brief A map from the indices of the active reactions to the indices of the corresponding reactions in the full network. */
        std::vector<size_t> m_reactionIndexMap;

        /** @brief A flag indicating whether the view is stale and needs to be updated. */
        bool m_isStale = true;

    private:
        /**
         * @brief A struct to hold a reaction and its flow rate.
         */
        struct ReactionFlow {
            const reaction::LogicalReaction* reactionPtr;
            double flowRate;
        };
    private:
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
        [[nodiscard]] std::vector<size_t> constructSpeciesIndexMap() const;

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
        [[nodiscard]] std::vector<size_t> constructReactionIndexMap() const;

        /**
         * @brief Maps a vector of culled abundances to a vector of full abundances.
         *
         * @param culled A vector of abundances for the active species.
         * @return A vector of abundances for the full network, with the abundances of the active
         *         species copied from the culled vector.
         */
        [[nodiscard]] std::vector<double> mapCulledToFull(const std::vector<double>& culled) const;

        /**
         * @brief Maps a vector of full abundances to a vector of culled abundances.
         *
         * @param full A vector of abundances for the full network.
         * @return A vector of abundances for the active species, with the abundances of the active
         *         species copied from the full vector.
         */
        [[nodiscard]] std::vector<double> mapFullToCulled(const std::vector<double>& full) const;

        /**
         * @brief Maps a culled species index to a full species index.
         *
         * @param culledSpeciesIndex The index of the species in the culled species list.
         * @return The index of the corresponding species in the full network.
         *
         * @throws std::out_of_range If the culled index is out of bounds for the species index map.
         */
        [[nodiscard]] size_t mapCulledToFullSpeciesIndex(size_t culledSpeciesIndex) const;

        /**
         * @brief Maps a culled reaction index to a full reaction index.
         *
         * @param culledReactionIndex The index of the reaction in the culled reaction list.
         * @return The index of the corresponding reaction in the full network.
         *
         * @throws std::out_of_range If the culled index is out of bounds for the reaction index map.
         */
        [[nodiscard]] size_t mapCulledToFullReactionIndex(size_t culledReactionIndex) const;

        /**
         * @brief Validates that the AdaptiveEngineView is not stale.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         */
        void validateState() const;

        /**
         * @brief Calculates the molar reaction flow rate for all reactions in the full network.
         *
         * This method iterates through all reactions in the base engine's network and calculates
         * their molar flow rates based on the provided network input conditions (temperature, density,
         * and composition). It also constructs a vector of molar abundances for all species in the
         * full network.
         *
         * @param netIn The current network input, containing temperature, density, and composition.
         * @param out_Y_Full A vector that will be populated with the molar abundances of all species in the full network.
         * @return A vector of ReactionFlow structs, each containing a pointer to a reaction and its calculated flow rate.
         *
         * @par Algorithm:
         * 1. Clears and reserves space in `out_Y_Full`.
         * 2. Iterates through all species in the base engine's network.
         * 3. For each species, it retrieves the molar abundance from `netIn.composition`. If the species is not found, its abundance is set to 0.0.
         * 4. Converts the temperature from Kelvin to T9.
         * 5. Iterates through all reactions in the base engine's network.
         * 6. For each reaction, it calls the base engine's `calculateMolarReactionFlow` to get the flow rate.
         * 7. Stores the reaction pointer and its flow rate in a `ReactionFlow` struct and adds it to the returned vector.
         */
        std::vector<ReactionFlow> calculateAllReactionFlows(
            const NetIn& netIn,
            std::vector<double>& out_Y_Full
        ) const;
        /**
         * @brief Finds all species that are reachable from the initial fuel through the reaction network.
         *
         * This method performs a connectivity analysis to identify all species that can be produced
         * starting from the initial fuel species. A species is considered part of the initial fuel if its
         * mass fraction is above a certain threshold (`ABUNDANCE_FLOOR`).
         *
         * @param netIn The current network input, containing the initial composition.
         * @return An unordered set of all reachable species.
         *
         * @par Algorithm:
         * 1. Initializes a set `reachable` and a queue `to_visit` with the initial fuel species.
         * 2. Iteratively processes the reaction network until no new species can be reached.
         * 3. In each pass, it iterates through all reactions in the base engine's network.
         * 4. If all reactants of a reaction are in the `reachable` set, all products of that reaction are added to the `reachable` set.
         * 5. The process continues until a full pass over all reactions does not add any new species to the `reachable` set.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> findReachableSpecies(
            const NetIn& netIn
        ) const;
        /**
         * @brief Culls reactions from the network based on their flow rates.
         *
         * This method filters the list of all reactions, keeping only those with a flow rate
         * above an absolute culling threshold. The threshold is calculated by multiplying the
         * maximum flow rate by a relative culling threshold read from the configuration.
         *
         * @param allFlows A vector of all reactions and their flow rates.
         * @param reachableSpecies A set of all species reachable from the initial fuel.
         * @param Y_full A vector of molar abundances for all species in the full network.
         * @param maxFlow The maximum reaction flow rate in the network.
         * @return A vector of pointers to the reactions that have been kept after culling.
         *
         * @par Algorithm:
         * 1. Retrieves the `RelativeCullingThreshold` from the configuration.
         * 2. Calculates the `absoluteCullingThreshold` by multiplying `maxFlow` with the relative threshold.
         * 3. Iterates through `allFlows`.
         * 4. A reaction is kept if its `flowRate` is greater than the `absoluteCullingThreshold`.
         * 5. The pointers to the kept reactions are stored in a vector and returned.
         */
        [[nodiscard]] std::vector<const reaction::LogicalReaction*> cullReactionsByFlow(
            const std::vector<ReactionFlow>& allFlows,
            const std::unordered_set<fourdst::atomic::Species>& reachableSpecies,
            const std::vector<double>& Y_full,
            double maxFlow
        ) const;

        typedef std::pair<std::unordered_set<const reaction::LogicalReaction*>, std::unordered_set<fourdst::atomic::Species>> RescueSet;
        [[nodiscard]] RescueSet rescueEdgeSpeciesDestructionChannel(
            const std::vector<double>& Y_full,
            const double T9,
            const double rho,
            const std::vector<fourdst::atomic::Species>& activeSpecies,
            const reaction::LogicalReactionSet& activeReactions
        ) const;
        /**
         * @brief Finalizes the set of active species and reactions.
         *
         * This method takes the final list of culled reactions and populates the
         * `m_activeReactions` and `m_activeSpecies` members. The active species are
         * determined by collecting all reactants and products from the final reactions.
         * The active species list is then sorted by mass.
         *
         * @param finalReactions A vector of pointers to the reactions to be included in the active set.
         *
         * @post
         * - `m_activeReactions` is cleared and populated with the reactions from `finalReactions`.
         * - `m_activeSpecies` is cleared and populated with all unique species present in `finalReactions`.
         * - `m_activeSpecies` is sorted by atomic mass.
         */
        void finalizeActiveSet(
            const std::vector<const reaction::LogicalReaction*>& finalReactions
        );
    };
}
