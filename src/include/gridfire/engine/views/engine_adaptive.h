#pragma once
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/types/types.h"
#include "gridfire/config/config.h"

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"

#include "gridfire/engine/scratchpads/blob.h"

#include "quill/Logger.h"

#include <expected>
#include <optional>

namespace gridfire::engine {
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
         * @param ctx The scratchpad context for storing thread-local data.
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
        fourdst::composition::Composition project(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const override;

        /**
         * @brief Gets the list of active species in the network.
         * @return A const reference to the vector of active species.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getNetworkSpecies(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation for the active species.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp The current composition of the system.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param trust
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
        [[nodiscard]] std::expected<StepDerivatives<double>, EngineStatus> calculateRHSAndEnergy(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            bool trust
        ) const override;


        /**
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp The current composition of the system.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @return A struct containing the derivatives of the energy generation rate with respect to temperature and density.
         */
        [[nodiscard]] EnergyDerivatives calculateEpsDerivatives(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the active species.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp The current composition of the system.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to generate the Jacobian matrix.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @see AdaptiveEngineView::update()
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for some set of active species such that that set is a subset of the active species in the view.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp The current composition of the system.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param activeSpecies The list of active species for which to generate the Jacobian.
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to generate the Jacobian matrix.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @see AdaptiveEngineView::update()
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const std::vector<fourdst::atomic::Species> &activeSpecies
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the active species with a given sparsity pattern.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp The current composition of the system.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param sparsityPattern The sparsity pattern to use for the Jacobian matrix.
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to generate the Jacobian matrix.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @see AdaptiveEngineView::update()
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
         * @param ctx The scratchpad context for storing thread-local data.
         * @param reaction The reaction for which to calculate the flow.
         * @param comp Composition object containing current abundances.
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
            scratch::StateBlob& ctx,
            const reaction::Reaction &reaction,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the set of active logical reactions in the network.
         *
         * @return Reference to the LogicalReactionSet containing all active reactions.
         */
        [[nodiscard]] const reaction::ReactionSet& getNetworkReactions(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Gets the set of inactive logical reactions in the network.
         *
         * @return ReactionSet containing all inactive reactions.
         *
         * This method returns the set of reactions that have been culled from the active
         * network based on the adaptation criteria.
         */
        [[nodiscard]] reaction::ReactionSet getInactiveNetworkReactions(
            scratch::StateBlob &ctx
        ) const override;

        [[nodiscard]] double getInactiveReactionMolarReactionFlow(
            scratch::StateBlob& ctx,
            const reaction::Reaction &reaction,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes timescales for all active species in the network.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to compute the species timescales.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesTimescales(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes destruction timescales for all active species in the network.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their destruction timescales (s).
         *
         * This method maps the culled abundances to the full network abundances and calls the base engine
         * to compute the species destruction timescales.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesDestructionTimescales(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the base engine.
         * @return A const reference to the base engine.
         */
        [[nodiscard]] const DynamicEngine& getBaseEngine() const override {
            return m_baseEngine;
        }

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
        [[nodiscard]] screening::ScreeningType getScreeningModel(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Gets the index of a species in the active species list.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param species The species for which to get the index.
         * @return The index of the species in the active species list.
         *
         * @throws std::runtime_error If the AdaptiveEngineView is stale (i.e., `update()` has not been called).
         * @throws std::out_of_range If the species is not part of the active species in the adaptive engine view.
         */
        [[nodiscard]] size_t getSpeciesIndex(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        ) const override;

        /**
         * @brief Primes the engine with the given network input.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param netIn The current network input, containing temperature, density, and composition.
         * @return A PrimingReport indicating the result of the priming operation.
         *
         * This method delegates the priming operation to the base engine.
         */
        [[nodiscard]] PrimingReport primeEngine(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const override;

        /**
         * @brief Collect the composition of the base engine, ensure all active species are registered, and pass
         * the composition back to the caller.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param comp The current composition of the system.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         *
         * @note This function ensures that the state of both the base engine and the adaptive view are synchronized in the
         * result back to the caller
         */
        [[nodiscard]] fourdst::composition::Composition collectComposition(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the status of a species in the network.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param species The species for which to get the status.
         * @return The SpeciesStatus indicating the status of the species.
         *
         * This method delegates the call to the base engine to get the species status. If the base engine says that
         * the species is active, but it is not in the active species list of this view, the status is returned as
         * INACTIVE_FLOW.
         */
        [[nodiscard]] SpeciesStatus getSpeciesStatus(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        ) const override;

        [[nodiscard]] std::optional<StepDerivatives<double>>getMostRecentRHSCalculation(
            scratch::StateBlob &ctx
        ) const override;

        [[nodiscard]] std::unique_ptr<scratch::StateBlob> constructStateBlob(const scratch::StateBlob *blob) const override;
    private:
        using LogManager = fourdst::logging::LogManager;

        fourdst::config::Config<config::GridFireConfig> m_config;

        /** @brief A pointer to the logger instance, used for logging messages. */
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        /** @brief The underlying engine to which this view delegates calculations. */
        DynamicEngine& m_baseEngine;

    private:
        /**
         * @brief A struct to hold a reaction and its flow rate.
         */
        struct ReactionFlow {
            const reaction::Reaction* reactionPtr;
            double flowRate;
        };
    private:
        /**
         * @brief Calculates the molar reaction flow rate for all reactions in the full network.
         *
         * This method iterates through all reactions in the base engine's network and calculates
         * their molar flow rates based on the provided network input conditions (temperature, density,
         * and composition). It also constructs a vector of molar abundances for all species in the
         * full network.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param netIn The current network input, containing temperature, density, and composition.
         * @return A pair with the first element a vector of ReactionFlow structs, each containing a pointer to a
         * reaction and its calculated flow rate and the second being a composition object where species which were not
         * present in netIn but are present in the definition of the base engine are registered but have 0 mass fraction.
         *
         * @par Algorithm:
         * 1. Iterates through all species in the base engine's network.
         * 2. For each species, it retrieves the molar abundance from `netIn.composition`. If the species is not found, its abundance is set to 0.0.
         * 3. Converts the temperature from Kelvin to T9.
         * 4. Iterates through all reactions in the base engine's network.
         * 5. For each reaction, it calls the base engine's `calculateMolarReactionFlow` to get the flow rate.
         * 6. Stores the reaction pointer and its flow rate in a `ReactionFlow` struct and adds it to the returned vector.
         */
        [[nodiscard]] std::pair<std::vector<ReactionFlow>, fourdst::composition::Composition> calculateAllReactionFlows(
            scratch::StateBlob& ctx,
            const NetIn& netIn
        ) const;
        /**
         * @brief Finds all species that are reachable from the initial fuel through the reaction network.
         *
         * This method performs a connectivity analysis to identify all species that can be produced
         * starting from the initial fuel species. A species is considered part of the initial fuel if its
         * mass fraction is above a certain threshold (`ABUNDANCE_FLOOR`).
         *
         * @param ctx The scratchpad context for storing thread-local data.
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
            scratch::StateBlob& ctx,
            const NetIn& netIn
        ) const;
        /**
         * @brief Culls reactions from the network based on their flow rates.
         *
         * This method filters the list of all reactions, keeping only those with a flow rate
         * above an absolute culling threshold. The threshold is calculated by multiplying the
         * maximum flow rate by a relative culling threshold read from the configuration.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param allFlows A vector of all reactions and their flow rates.
         * @param reachableSpecies A set of all species reachable from the initial fuel.
         * @param comp The current composition of the system.
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
        [[nodiscard]] std::vector<const reaction::Reaction*> cullReactionsByFlow(
            scratch::StateBlob& ctx,
            const std::vector<ReactionFlow>& allFlows,
            const std::unordered_set<fourdst::atomic::Species>& reachableSpecies,
            const fourdst::composition::Composition& comp,
            double maxFlow
        ) const;

        typedef std::pair<std::unordered_set<const reaction::Reaction*>, std::unordered_set<fourdst::atomic::Species>> RescueSet;
        [[nodiscard]] RescueSet rescueEdgeSpeciesDestructionChannel(
            scratch::StateBlob& ctx,
            const fourdst::composition::Composition& comp,
            double T9,
            double rho
        ) const;
        /**
         * @brief Finalizes the set of active species and reactions.
         *
         * This method takes the final list of culled reactions and populates the
         * `m_activeReactions` and `m_activeSpecies` members. The active species are
         * determined by collecting all reactants and products from the final reactions.
         * The active species list is then sorted by mass.
         *
         * @param ctx The scratchpad context for storing thread-local data.
         * @param finalReactions A vector of pointers to the reactions to be included in the active set.
         *
         * @post
         * - `m_activeReactions` is cleared and populated with the reactions from `finalReactions`.
         * - `m_activeSpecies` is cleared and populated with all unique species present in `finalReactions`.
         * - `m_activeSpecies` is sorted by atomic mass.
         */
        void finalizeActiveSet(
            scratch::StateBlob& ctx,
            const std::vector<const reaction::Reaction*>& finalReactions
        ) const;
    };
}
