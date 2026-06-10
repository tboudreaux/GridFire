#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_graph.h"

#include "gridfire/engine/scratchpads/blob.h"

#include "sundials/sundials_linearsolver.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

#include <expected>
#include <optional>

namespace gridfire::engine {
    /**
     * @class MultiscalePartitioningEngineView
     * @brief An engine view that partitions the reaction network into multiple groups based on timescales.
     *
     * @par Purpose
     * This class is designed to accelerate the integration of stiff nuclear reaction networks.
     * It identifies species that react on very short timescales ("fast" species) and treats them
     * as being in Quasi-Steady-State Equilibrium (QSE). Their abundances are solved for algebraically,
     * removing their stiff differential equations from the system. The remaining "slow" or "dynamic"
     * species are integrated normally. This significantly improves the stability and performance of
     * the solver.
     *
     * @par How
     * The core logic resides in the `partitionNetwork()` and `equilibrateNetwork()` methods.
     * The partitioning process involves:
     *   1.  **Timescale Analysis:** Using `getSpeciesDestructionTimescales` from the base engine,
     *       all species are sorted by their characteristic timescales.
     *   2.  **Gap Detection:** The sorted list of timescales is scanned for large gaps (e.g., several
     *       orders of magnitude) to create distinct "timescale pools".
     *   3.  **Connectivity Analysis:** Each pool is analyzed for internal reaction connectivity to
     *       form cohesive groups.
     *   4.  **Flux Validation:** Candidate QSE groups are validated by comparing the total reaction
     *       flux *within* the group to the flux *leaving* the group. A high internal-to-external
     *       flux ratio indicates a valid QSE group.
     *   5.  **QSE Solve:** For valid QSE groups, `solveQSEAbundances` uses a Levenberg-Marquardt
     *       nonlinear solver (`Eigen::LevenbergMarquardt`) to find the equilibrium abundances of the
     *       "algebraic" species, holding the "seed" species constant.
     *
     * All calculations are cached using `QSECacheKey` to avoid re-partitioning and re-solving for
     * similar thermodynamic conditions.
     *
     * @par Usage Example:
     * @code
     * // 1. Create a base engine (e.g., GraphEngine)
     * gridfire::GraphEngine baseEngine(composition);
     *
     * // 2. Wrap it with the MultiscalePartitioningEngineView
     * gridfire::MultiscalePartitioningEngineView multiscaleEngine(baseEngine);
     *
     * // 3. Before integration, update the view to partition the network
     * //    and find the initial equilibrium state.
     * NetIn initialConditions = { .composition = composition, .temperature = 1e8, .density = 1e3 };
     * fourdst::composition::Composition equilibratedComp = multiscaleEngine.update(initialConditions);
     *
     * // 4. Use the multiscaleEngine for integration. It will use the cached QSE solution.
     * //    The integrator will call calculateRHSAndEnergy, etc. on the multiscaleEngine.
     * auto Y_initial = multiscaleEngine.mapNetInToMolarAbundanceVector({equilibratedComp, ...});
     * auto derivatives = multiscaleEngine.calculateRHSAndEnergy(Y_initial, T9, rho);
     * @endcode
     *
     * @implements DynamicEngine
     * @implements EngineView<DynamicEngine>
     */
    class MultiscalePartitioningEngineView final: public DynamicEngine, public EngineView<DynamicEngine> {
        /**
         * @brief Type alias for a QSE partition.
         *
         * A QSE partition is a tuple containing the fast species, their indices,
         * the slow species, and their indices.
         */
        typedef std::tuple<std::vector<fourdst::atomic::Species>, std::vector<size_t>, std::vector<fourdst::atomic::Species>, std::vector<size_t>> QSEPartition;
    public:
        /**
         * @brief Constructs a MultiscalePartitioningEngineView.
         *
         * @param baseEngine The underlying GraphEngine to which this view delegates calculations.
         *                   It must be a `GraphEngine` and not a more general `DynamicEngine`
         *                   because this view relies on its specific implementation details.
         */
        explicit MultiscalePartitioningEngineView(DynamicEngine& baseEngine);

        /**
         * @brief Gets the list of species in the network.
         * @return A const reference to the vector of `Species` objects representing all species
         *         in the underlying base engine. This view does not alter the species list itself,
         *         only how their abundances are evolved.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species> & getNetworkSpecies(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param trust
         * @return A `std::expected` containing `StepDerivatives<double>` on success, or a
         *         `StaleEngineError` if the engine's QSE cache does not contain a solution
         *         for the given state.
         *
         * @par Purpose
         * To compute the time derivatives for the ODE solver. This implementation
         * modifies the derivatives from the base engine to enforce the QSE condition.
         *
         * @par How
         * It first performs a lookup in the QSE abundance cache (`m_qse_abundance_cache`).
         * If a cache hit occurs, it calls the base engine's `calculateRHSAndEnergy`. It then
         * manually sets the time derivatives (`dydt`) of all identified algebraic species to zero,
         * effectively removing their differential equations from the system being solved.
         *
         * @pre The engine must have been updated via `update()` or `equilibrateNetwork()` for the
         *      current thermodynamic conditions, so that a valid entry exists in the QSE cache.
         * @post The returned derivatives will have `dydt=0` for all algebraic species.
         *
         * @throws StaleEngineError If the QSE cache does not contain an entry for the given
         *         (T9, rho, Y_full). This indicates `update()` was not called recently enough.
         */
        [[nodiscard]] std::expected<StepDerivatives<double>, engine::EngineStatus> calculateRHSAndEnergy(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            bool trust
        ) const override;

        /**
         * @brief Calculates the energy generation rate derivatives with respect to abundances.
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @return The energy generation rate derivatives (dEps/dT and dEps/drho).
         */
        [[nodiscard]] EnergyDerivatives calculateEpsDerivatives(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * @par Purpose
         * To compute the Jacobian matrix required by implicit ODE solvers.
         *
         * @par How
         * It first performs a QSE cache lookup. On a hit, it delegates the full Jacobian
         * calculation to the base engine. While this view could theoretically return a
         * modified, sparser Jacobian reflecting the QSE constraints, the current implementation
         * returns the full Jacobian from the base engine. The solver is expected to handle the
         * algebraic constraints (e.g., via `dydt=0` from `calculateRHSAndEnergy`).
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @post The base engine's internal Jacobian is updated.
         *
         * @throws exceptions::StaleEngineError If the QSE cache misses, as it cannot proceed
         *         without a valid partition.
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for a subset of active species.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param activeSpecies The subset of species to include in the Jacobian.
         *
         * @par Purpose
         * To compute a reduced Jacobian matrix for implicit solvers that only
         * consider a subset of species.
         *
         * @par How
         * Similar to the full Jacobian generation, it first checks the QSE cache.
         * On a hit, it calls the base engine's `generateJacobianMatrix` with
         * the specified active species. The returned Jacobian still reflects
         * the full network, but only for the active species subset.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @post The base engine's internal Jacobian is updated for the active species.
         *
         * @throws exceptions::StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const std::vector<fourdst::atomic::Species> &activeSpecies
        ) const override;

        /**
         * @brief Generates the Jacobian matrix using a sparsity pattern.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param sparsityPattern The sparsity pattern to use for the Jacobian.
         *
         * @par Purpose
         * To compute the Jacobian matrix while leveraging a known sparsity pattern
         * for efficiency. This is effectively a lower level version of the active species method.
         *
         * @par How
         * It first checks the QSE cache. On a hit, it delegates to the base engine's
         * `generateJacobianMatrix` method with the provided sparsity pattern.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @post The base engine's internal Jacobian is updated according to the sparsity pattern.
         *
         * @throws exceptions::StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const SparsityPattern &sparsityPattern
        ) const override;

        /**
         * @brief Calculates the molar reaction flow for a given reaction.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param reaction The reaction for which to calculate the flow.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * @par Purpose
         * To compute the net rate of a single reaction.
         *
         * @par How
         * It first checks the QSE cache. On a hit, it retrieves the cached equilibrium
         *      abundances for the algebraic species. It creates a mutable copy of `Y_full`,
         *      overwrites the algebraic species abundances with the cached equilibrium values,
         *      and then calls the base engine's `calculateMolarReactionFlow` with this modified
         *      abundance vector.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @throws StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] double calculateMolarReactionFlow(
            scratch::StateBlob& ctx,
            const reaction::Reaction &reaction,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the set of logical reactions in the network.
         *
         * @return A const reference to the `LogicalReactionSet` from the base engine,
         *         containing all reactions in the full network.
         */
        [[nodiscard]] const reaction::ReactionSet & getNetworkReactions(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Computes timescales for all species in the network.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A `std::expected` containing a map from `Species` to their characteristic
         *         timescales (s) on success, or a `StaleEngineError` on failure.
         *
         * @par Purpose
         * To get the characteristic timescale `Y / (dY/dt)` for each species.
         *
         * @par How
         * It delegates the calculation to the base engine. For any species identified
         *      as algebraic (in QSE), it manually sets their timescale to 0.0 to signify
         *      that they equilibrate instantaneously on the timescale of the solver.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @throws StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesTimescales(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes destruction timescales for all species in the network.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A `std::expected` containing a map from `Species` to their characteristic
         *         destruction timescales (s) on success, or a `StaleEngineError` on failure.
         *
         * @par Purpose
         * To get the timescale for species destruction, which is used as the primary
         *          metric for network partitioning.
         *
         * @par How
         * It delegates the calculation to the base engine. For any species identified
         *      as algebraic (in QSE), it manually sets their timescale to 0.0.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @throws StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesDestructionTimescales(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Updates the internal state of the engine, performing partitioning and QSE equilibration.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param netIn A struct containing the current network input: temperature, density, and composition.
         * @return The new composition after QSE species have been brought to equilibrium.
         *
         * @par Purpose
         * This is the main entry point for preparing the multiscale engine for use. It
         * triggers the network partitioning and solves for the initial QSE abundances, caching the result.
         *
         * @how
         * 1.  It first checks the QSE cache. If a valid entry already exists for the input state,
         *     it returns the input composition, as no work is needed.
         * 2.  If the cache misses, it calls `equilibrateNetwork()`.
         * 3.  `equilibrateNetwork()` in turn calls `partitionNetwork()` to define the dynamic and
         *     algebraic species sets.
         * 4.  It then calls `solveQSEAbundances()` to compute the equilibrium abundances.
         * 5.  The resulting equilibrium abundances for the algebraic species are stored in the
         *     `m_qse_abundance_cache`.
         * 6.  A new `fourdst::composition::Composition` object reflecting the equilibrated state
         *     is created and returned.
         *
         * @pre The `netIn` struct should contain a valid physical state.
         * @post The engine is partitioned (`m_dynamic_species`, `m_algebraic_species`, etc. are populated).
         *       The `m_qse_abundance_cache` is populated with the QSE solution for the given state.
         *       The returned composition reflects the new equilibrium.
         */
        fourdst::composition::Composition project(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const override;


        /**
         * @brief Gets the current electron screening model.
         *
         * @return The currently active screening model type.
         *
         * @par How
         * This method delegates directly to the base engine's `getScreeningModel()`.
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel(
            scratch::StateBlob& ctx
        ) const override;

        /**
         * @brief Gets the base engine.
         *
         * @return A const reference to the base engine.
         */
        [[nodiscard]] const DynamicEngine & getBaseEngine() const override;



        /**
         * @brief Partitions the network based on timescales from a `NetIn` struct.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param netIn A struct containing the current network input.
         *
         * @par Purpose
         * A convenience overload for `partitionNetwork`.
         *
         * @par How
         * It unpacks the `netIn` struct into `Y`, `T9`, and `rho` and then calls the
         *      primary `partitionNetwork` method.
         */
        fourdst::composition::Composition partitionNetwork(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const;

        /**
         * @brief Exports the network to a DOT file for visualization.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param filename The name of the DOT file to create.
         * @param comp Composition object
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * @par Purpose
         * To visualize the partitioned network graph.
         *
         * @par How
         * This method delegates the DOT file export to the base engine. It does not
         *      currently add any partitioning information to the output graph.
         */
        void exportToDot(
            scratch::StateBlob &ctx,
            const std::string& filename,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Gets the index of a species in the full network.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param species The species to get the index of.
         * @return The index of the species in the base engine's network.
         *
         * @par How
         * This method delegates directly to the base engine's `getSpeciesIndex()`.
         */
        [[nodiscard]] size_t getSpeciesIndex(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        ) const override;

        /**
         * @brief Primes the engine with a specific species.
         *
         * @param netIn A struct containing the current network input.
         * @return A `PrimingReport` struct containing information about the priming process.
         *
         * @par Purpose
         * To prepare the network for ignition or specific pathway studies.
         *
         * @par How
         * This method delegates directly to the base engine's `primeEngine()`. The
         *      multiscale view does not currently interact with the priming process.
         */
        [[nodiscard]] PrimingReport primeEngine(
            scratch::StateBlob& ctx,
            const NetIn &netIn
        ) const override;

        /**
         * @brief Gets the fast species in the network.
         *
         * @return A vector of species identified as "fast" or "algebraic" by the partitioning.
         *
         * @par Purpose
         * To allow external queries of the partitioning results.
         *
         * @par How
         * It returns a copy of the `m_algebraic_species` member vector.
         *
         * @pre `partitionNetwork()` must have been called.
         */
        [[nodiscard]] std::vector<fourdst::atomic::Species> getFastSpecies(
            scratch::StateBlob& ctx
        ) const;

        /**
         * @brief Gets the dynamic species in the network.
         *
         * @return A const reference to the vector of species identified as "dynamic" or "slow".
         *
         * @par Purpose
         * To allow external queries of the partitioning results.
         *
         * @par How
         * It returns a const reference to the `m_dynamic_species` member vector.
         *
         * @pre `partitionNetwork()` must have been called.
         */
        [[nodiscard]] static const std::vector<fourdst::atomic::Species>& getDynamicSpecies(
            scratch::StateBlob& ctx
        );

        /**
         * @brief Checks if a species is involved in the partitioned network.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param species The species to check.
         * @return `true` if the species is in either the dynamic or algebraic sets, `false` otherwise.
         *
         * @par Purpose
         * To allow external queries about species involvement in the partitioned network.
         *
         * @par How
         * It checks for membership in both `m_dynamic_species` and `m_algebraic_species`.
         *
         * @pre `partitionNetwork()` must have been called.
         */
        static bool involvesSpecies(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        );

        /**
         * @brief Check if a species is involved in the QSE (algebraic) set.
         * @param ctx The scratch data for thread-local storage.
         * @param species The species to check.
         * @return Boolean indicating if the species is in the algebraic set.
         */
        static bool involvesSpeciesInQSE(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        );

        /**
         * @brief Check if a species is involved in the dynamic set.
         * @param ctx The scratch data for thread-local storage.
         * @param species The species to check.
         * @return Boolean indicating if the species is in the dynamic set.
         */
        static bool involvesSpeciesInDynamic(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        );

        /**
         * @brief Gets a normalized composition with QSE species equilibrated.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp The input composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param trust
         * @return A new `Composition` object with algebraic species set to their equilibrium values.
         *
         * @par Purpose
         * To provide a way to get the equilibrated composition without modifying the internal state.
         *
         * @par How
         * It calls `solveQSEAbundances()` to compute the equilibrium abundances for the algebraic species,
         * then constructs a new `Composition` object reflecting these values.
         *
         * @pre The engine must have a valid QSE partition for the given state.
         * @throws StaleEngineError If the QSE cache misses.
         */
        fourdst::composition::Composition getNormalizedEquilibratedComposition(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract& comp,
            double T9,
            double rho,
            bool trust
        ) const;

        /**
         * @brief Collect the composition from this and sub engines.
         * @details This method operates by injecting the current equilibrium abundances for algebraic species into
         * the composition object so that they can be bubbled up to the caller.
         * @param ctx The scratch data for thread-local storage.
         * @param comp Input Composition
         * @param T9
         * @param rho
         * @return New composition which is comp + any edits from lower levels + the equilibrium abundances of all algebraic species.
         * @throws BadCollectionError: if there is a species in the algebraic species set which does not show up in the reported composition from the base engine.:w
         */
        fourdst::composition::Composition collectComposition(
            scratch::StateBlob& ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the status of a species in the network.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param species The species to query.
         * @return The `SpeciesStatus` indicating if the species is dynamic, algebraic, or not involved.
         *
         * @par Purpose
         * To allow external queries about the role of a species in the partitioned network.
         *
         * @par How
         * It checks for membership in `m_dynamic_species` and `m_algebraic_species` to determine
         * the appropriate status.
         *
         * @pre `partitionNetwork()` must have been called.
         */
        SpeciesStatus getSpeciesStatus(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species &species
        ) const override;

        [[nodiscard]] std::optional<StepDerivatives<double>>getMostRecentRHSCalculation(
            scratch::StateBlob &
        ) const override;

        [[nodiscard]] std::unique_ptr<scratch::StateBlob> constructStateBlob(const scratch::StateBlob *blob) const override;
    public:
        /**
         * @brief Struct representing a QSE group.
         *
         * @par Purpose
         * A container to hold all information about a set of species that are potentially
         * in quasi-steady-state equilibrium with each other.
         */
        struct QSEGroup {
            bool is_in_equilibrium = false;                       ///< Flag set by flux analysis.
            std::set<fourdst::atomic::Species> algebraic_species; ///< Algebraic species in this group.
            std::set<fourdst::atomic::Species> seed_species;      ///< Dynamic species in this group.
            double mean_timescale;                                ///< Mean timescale of the group.

            /**
             * @brief Less-than operator for QSEGroup, used for sorting.
             * @param other The other QSEGroup to compare to.
             * @return True if this group's mean timescale is less than the other's.
             */
            bool operator<(const QSEGroup& other) const;
            /**
             * @brief Greater-than operator for QSEGroup.
             * @param other The other QSEGroup to compare to.
             * @return True if this group's mean timescale is greater than the other's.
             */
            bool operator>(const QSEGroup& other) const;
            /**
             * @brief Equality operator for QSEGroup.
             * @param other The other QSEGroup to compare to.
             * @return True if the sets of species indices are identical.
             */
            bool operator==(const QSEGroup& other) const;
            /**
             * @brief Inequality operator for QSEGroup.
             * @param other The other QSEGroup to compare to.
             * @return True if the sets of species indices are not identical.
             */
            bool operator!=(const QSEGroup& other) const;

            [[nodiscard]] [[maybe_unused]] std::string toString(bool verbose) const;

            friend std::ostream& operator<<(std::ostream& os, const QSEGroup& group) {
                os << group.toString(false);
                return os;
            }

            [[nodiscard]] bool contains(const fourdst::atomic::Species& species) const;
            [[nodiscard]] bool containsAlgebraic(const fourdst::atomic::Species &species) const;
            [[nodiscard]] bool containsSeed(const fourdst::atomic::Species &species) const;
        };

        class QSESolver {
        public:
            QSESolver(
                const std::vector<fourdst::atomic::Species>& species,
                const DynamicEngine& engine,
                SUNContext sun_ctx
            );

            QSESolver(
                const QSESolver& other
            ) = delete;
            QSESolver& operator=(
                const QSESolver& other
            ) = delete;

            ~QSESolver();

            fourdst::composition::Composition solve(
                scratch::StateBlob& ctx,
                const fourdst::composition::Composition& comp,
                double T9,
                double rho
            ) const;

            size_t solves() const;

            void log_diagnostics(const QSEGroup &group, const fourdst::composition::Composition &comp) const;

            std::unique_ptr<QSESolver> clone() const;
            std::unique_ptr<QSESolver> clone(SUNContext sun_ctx) const;
        private:

            static int sys_func(
                N_Vector y,
                N_Vector f,
                void *user_data
            );
            static int sys_jac(
                N_Vector y,
                N_Vector fy,
                SUNMatrix J,
                void *user_data,
                N_Vector tmp1,
                N_Vector tmp2
            );

            struct UserData {
                const DynamicEngine& engine;
                double T9;
                double rho;
                fourdst::composition::Composition& comp;
                const std::unordered_map<fourdst::atomic::Species, size_t>& qse_solve_species_index_map;
                const std::vector<fourdst::atomic::Species>& qse_solve_species;
                const QSESolver& instance;
                scratch::StateBlob& ctx;
            };

        private:
            // Cache members
            mutable size_t m_solves = 0;
            mutable bool m_has_jacobian = false;

            // Solver members
            size_t m_N;
            const DynamicEngine& m_engine;
            std::vector<fourdst::atomic::Species> m_species;
            std::unordered_map<fourdst::atomic::Species, size_t> m_speciesMap;

            // --- SUNDIALS / KINSOL Persistent resources
            SUNContext m_sun_ctx   = nullptr;
            void* m_kinsol_mem     = nullptr;

            N_Vector m_Y           = nullptr;
            N_Vector m_scale       = nullptr;
            N_Vector m_f_scale     = nullptr;
            N_Vector m_constraints = nullptr;
            N_Vector m_func_tmpl   = nullptr;

            SUNMatrix m_J          = nullptr;
            SUNLinearSolver m_LS   = nullptr;

            static quill::Logger* getLogger() {
                return LogManager::getInstance().getLogger("log");
            }


        };

        struct FluxValidationResult {
            std::vector<QSEGroup> valid_groups;
            std::vector<QSEGroup> invalid_groups;
            std::vector<reaction::ReactionSet> validatedGroupReactions;
        };

    private:
        /**
         * @brief Logger instance for logging messages.
         */
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        /**
         * @brief The base engine to which this view delegates calculations.
         */
        DynamicEngine& m_baseEngine;

    private:
        /**
         * @brief Partitions the network by timescale.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp Vector of current molar abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of vectors of species indices, where each inner vector represents a
         *         timescale pool.
         *
         * @par Purpose
         * To group species into "pools" based on their destruction timescales.
         *
         * @par How
         * It retrieves all species destruction timescales from the base engine, sorts them,
         *      and then iterates through the sorted list, creating a new pool whenever it detects
         *      a gap between consecutive timescales that is larger than a predefined threshold
         *      (e.g., a factor of 100).
         */
        std::vector<std::vector<fourdst::atomic::Species>> partitionByTimescale(
            scratch::StateBlob& ctx,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        std::pair<bool, reaction::ReactionSet> group_is_a_qse_cluster(
            scratch::StateBlob& ctx,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho,
            const QSEGroup &group
        ) const;

        bool group_is_a_qse_pipeline(
            scratch::StateBlob& ctx,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho,
            const QSEGroup &group
        ) const;

        /**
         * @brief Validates candidate QSE groups using flux analysis.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param candidate_groups A vector of candidate QSE groups.
         * @param comp Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of validated QSE groups that meet the flux criteria.
         *
         * @par Purpose
         * To ensure that a candidate QSE group is truly in equilibrium by checking that
         *          the reaction fluxes *within* the group are much larger than the fluxes
         *          *leaving* the group.
         *
         * @par How
         * For each candidate group, it calculates the sum of all internal reaction fluxes and
         *      the sum of all external (bridge) reaction fluxes. If the ratio of internal to external
         *      flux exceeds a configurable threshold, the group is considered valid and is added
         *      to the returned vector.
         */
        FluxValidationResult validateGroupsWithFluxAnalysis(
            scratch::StateBlob& ctx,
            const std::vector<QSEGroup> &candidate_groups,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Solves for the QSE abundances of the algebraic species in a given state.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param comp Vector of current molar abundances for all species in the base engine.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of molar abundances for the algebraic species.
         *
         * @par Purpose
         * To find the equilibrium abundances of the algebraic species that satisfy
         *          the QSE conditions.
         *
         * @par How
         * It uses the Levenberg-Marquardt algorithm via Eigen's `LevenbergMarquardt` class.
         *      The problem is defined by the `EigenFunctor` which computes the residuals and
         *      Jacobian for the QSE equations.
         *
         * @pre The input state (Y_full, T9, rho) must be a valid physical state.
         * @post The algebraic species in the QSE cache are updated with the new equilibrium abundances.
         */
        auto solveQSEAbundances(
            scratch::StateBlob &ctx,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const -> fourdst::composition::Composition;

        /**
         * @brief Identifies the pool with the slowest mean timescale.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param pools A vector of vectors of species indices, where each inner vector represents a
         *              timescale pool.
         * @param comp Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return The index of the pool with the largest (slowest) mean destruction timescale.
         *
         * @par Purpose
         * To identify the core set of dynamic species that will not be part of any QSE group.
         *
         * @par How
         * It calculates the geometric mean of the destruction timescales for all species in each
         *      pool and returns the index of the pool with the maximum mean timescale.
         */
        size_t identifyMeanSlowestPool(
            scratch::StateBlob& ctx,
            const std::vector<std::vector<fourdst::atomic::Species>>& pools,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Builds a connectivity graph from a species pool.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param species_pool A vector of species indices representing a species pool.
         * @param comp
         * @param T9
         * @param rho
         * @return An unordered map representing the adjacency list of the connectivity graph.
         *
         * @par Purpose
         * To find reaction connections within a specific group of species.
         *
         * @par How
         * It iterates through all reactions in the base engine. If a reaction involves
         *      at least two distinct species from the input `species_pool` (one as a reactant
         *      and one as a product), it adds edges between all reactants and products from
         *      that reaction that are also in the pool.
         */
        std::unordered_map<fourdst::atomic::Species, std::vector<fourdst::atomic::Species>> buildConnectivityGraph(
            scratch::StateBlob& ctx,
            const std::vector<fourdst::atomic::Species>& species_pool,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Constructs candidate QSE groups from connected timescale pools.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param candidate_pools A vector of vectors of species indices, where each inner vector
         *                        represents a connected pool of species with similar fast timescales.
         * @param comp Vector of current molar abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of `QSEGroup` structs, ready for flux validation.
         *
         * @par How
         * For each input pool, it identifies "bridge" reactions that connect the pool to
         * species outside the pool. The reactants of these bridge reactions that are *not* in the
         * pool are identified as "seed" species. The original pool members are the "algebraic"
         * species. It then bundles the seed and algebraic species into a `QSEGroup` struct.
         *
         * @pre The `candidate_pools` should be connected components from `analyzeTimescalePoolConnectivity`.
         * @post A list of candidate `QSEGroup` objects is returned.
         */
        std::vector<QSEGroup> constructCandidateGroups(
            scratch::StateBlob& ctx,
            const std::vector<std::vector<fourdst::atomic::Species>>& candidate_pools,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Analyzes the connectivity of timescale pools.
         *
         * @param ctx The scratch data for thread-local storage.
         * @param timescale_pools A vector of vectors of species indices, where each inner vector
         *                        represents a timescale pool.
         * @param comp
         * @param T9
         * @param rho
         * @return A vector of vectors of species indices, where each inner vector represents a
         *         single connected component.
         *
         * @par Purpose
         * To merge timescale pools that are strongly connected by reactions, forming
         *          cohesive groups for QSE analysis.
         *
         * @par How
         * For each pool, it builds a reaction connectivity graph using `buildConnectivityGraph`.
         *      It then finds the connected components within that graph using a Breadth-First Search (BFS).
         *      The resulting components from all pools are collected and returned.
         */
        std::vector<std::vector<fourdst::atomic::Species>> analyzeTimescalePoolConnectivity(
            scratch::StateBlob& ctx,
            const std::vector<std::vector<fourdst::atomic::Species>> &timescale_pools,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        std::vector<QSEGroup> pruneValidatedGroups(
            scratch::StateBlob& ctx,
            const std::vector<QSEGroup> &groups,
            const std::vector<reaction::ReactionSet> &groupReactions,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        static std::vector<QSEGroup> merge_coupled_groups(
            scratch::StateBlob& ctx,
            const std::vector<QSEGroup> &groups,
            const std::vector<reaction::ReactionSet> &groupReactions
        );
    };
}

