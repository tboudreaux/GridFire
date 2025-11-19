#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_graph.h"
#include "sundials/sundials_linearsolver.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

#include "unsupported/Eigen/NonLinearOptimization"

namespace gridfire {
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

        ~MultiscalePartitioningEngineView() override;

        /**
         * @brief Gets the list of species in the network.
         * @return A const reference to the vector of `Species` objects representing all species
         *         in the underlying base engine. This view does not alter the species list itself,
         *         only how their abundances are evolved.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species> & getNetworkSpecies() const override;

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation.
         *
         * @param comp The current composition.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
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
        [[nodiscard]] std::expected<StepDerivatives<double>, expectations::StaleEngineError> calculateRHSAndEnergy(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] EnergyDerivatives calculateEpsDerivatives(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state.
         *
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
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for a subset of active species.
         *
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
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const std::vector<fourdst::atomic::Species> &activeSpecies
        ) const override;

        /**
         * @brief Generates the Jacobian matrix using a sparsity pattern.
         *
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
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const SparsityPattern &sparsityPattern
        ) const override;

        /**
         * @brief Generates the stoichiometry matrix for the network.
         *
         * @par Purpose
         * To prepare the stoichiometry matrix for later queries.
         *
         * @par How
         * This method delegates directly to the base engine's `generateStoichiometryMatrix()`.
         *      The stoichiometry is based on the full, unpartitioned network.
         */
        void generateStoichiometryMatrix() override;

        /**
         * @brief Gets an entry from the stoichiometry matrix.
         *
         * @param species Species to look up stoichiometry for.
         * @param reaction Reaction to find.
         * @return Stoichiometric coefficient for the species in the reaction.
         *
         * @par Purpose
         * To query the stoichiometric relationship between a species and a reaction.
         *
         * @par How
         * This method delegates directly to the base engine's `getStoichiometryMatrixEntry()`.
         *
         * @pre `generateStoichiometryMatrix()` must have been called.
         */
        [[nodiscard]] int getStoichiometryMatrixEntry(
            const fourdst::atomic::Species& species,
            const reaction::Reaction& reaction
        ) const override;

        /**
         * @brief Calculates the molar reaction flow for a given reaction.
         *
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
        [[nodiscard]] const reaction::ReactionSet & getNetworkReactions() const override;

        /**
         * @brief Sets the set of logical reactions in the network.
         *
         * @param reactions The set of logical reactions to use.
         *
         * @par Purpose
         * To modify the reaction network.
         *
         * @par How
         * This operation is not supported by the `MultiscalePartitioningEngineView` as it
         *      would invalidate the partitioning logic. It logs a critical error and throws an
         *      exception. Network modifications should be done on the base engine before it is
         *      wrapped by this view.
         *
         * @throws exceptions::UnableToSetNetworkReactionsError Always.
         */
        void setNetworkReactions(
            const reaction::ReactionSet &reactions
        ) override;

        /**
         * @brief Computes timescales for all species in the network.
         *
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
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes destruction timescales for all species in the network.
         *
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
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Updates the internal state of the engine, performing partitioning and QSE equilibration.
         *
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
        fourdst::composition::Composition update(
            const NetIn &netIn
        ) override;

        /**
         * @brief Checks if the engine's internal state is stale relative to the provided conditions.
         *
         * @param netIn A struct containing the current network input.
         * @return `true` if the engine is stale, `false` otherwise.
         *
         * @par Purpose
         * To determine if `update()` needs to be called.
         *
         * @par How
         * It creates a `QSECacheKey` from the `netIn` data and checks for its
         * existence in the `m_qse_abundance_cache`. A cache miss indicates the engine is
         * stale because it does not have a valid QSE partition for the current conditions.
         * It also queries the base engine's `isStale()` method.
         */
        bool isStale(const NetIn& netIn) override;

        /**
         * @brief Sets the electron screening model.
         *
         * @param model The type of screening model to use for reaction rate calculations.
         *
         * @par How
         * This method delegates directly to the base engine's `setScreeningModel()`.
         */
        void setScreeningModel(
            screening::ScreeningType model
        ) override;

        /**
         * @brief Gets the current electron screening model.
         *
         * @return The currently active screening model type.
         *
         * @par How
         * This method delegates directly to the base engine's `getScreeningModel()`.
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        /**
         * @brief Gets the base engine.
         *
         * @return A const reference to the base engine.
         */
        const DynamicEngine & getBaseEngine() const override;



        /**
         * @brief Partitions the network based on timescales from a `NetIn` struct.
         *
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
            const NetIn &netIn
        );

        /**
         * @brief Exports the network to a DOT file for visualization.
         *
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
            const std::string& filename,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Gets the index of a species in the full network.
         *
         * @param species The species to get the index of.
         * @return The index of the species in the base engine's network.
         *
         * @par How
         * This method delegates directly to the base engine's `getSpeciesIndex()`.
         */
        [[nodiscard]] size_t getSpeciesIndex(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Maps a `NetIn` struct to a molar abundance vector for the full network.
         *
         * @param netIn A struct containing the current network input.
         * @return A vector of molar abundances corresponding to the species order in the base engine.
         *
         * @par How
         * This method delegates directly to the base engine's `mapNetInToMolarAbundanceVector()`.
         */
        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

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
        [[nodiscard]] PrimingReport primeEngine(const NetIn &netIn) override;

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
        [[nodiscard]] std::vector<fourdst::atomic::Species> getFastSpecies() const;
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
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getDynamicSpecies() const;


        bool involvesSpecies(const fourdst::atomic::Species &species) const;

        bool involvesSpeciesInQSE(const fourdst::atomic::Species &species) const;

        bool involvesSpeciesInDynamic(const fourdst::atomic::Species &species) const;

        fourdst::composition::Composition getNormalizedEquilibratedComposition(const fourdst::composition::CompositionAbstract& comp, double T9, double rho) const;

        /**
         * @brief Collect the composition from this and sub engines.
         * @details This method operates by injecting the current equilibrium abundances for algebraic species into
         * the composition object so that they can be bubbled up to the caller.
         * @param comp Input Composition
         * @param T9
         * @param rho
         * @return New composition which is comp + any edits from lower levels + the equilibrium abundances of all algebraic species.
         * @throws BadCollectionError: if there is a species in the algebraic species set which does not show up in the reported composition from the base engine.:w
         */
        fourdst::composition::Composition collectComposition(const fourdst::composition::CompositionAbstract &comp, double T9, double rho) const override;

        SpeciesStatus getSpeciesStatus(const fourdst::atomic::Species &species) const override;


    private:
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

            // DEBUG METHODS.
            // THESE SHOULD NOT BE USED IN PRODUCTION CODE.
            [[deprecated("Use for debug only")]] void removeSpecies(const fourdst::atomic::Species& species);

            [[deprecated("Use for debug only")]] void addSpeciesToAlgebraic(const fourdst::atomic::Species& species);

            [[deprecated("Use for debug only")]] void addSpeciesToSeed(const fourdst::atomic::Species& species);


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

            [[nodiscard]] [[maybe_unused]] std::string toString() const;
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
                const fourdst::composition::Composition& comp,
                double T9,
                double rho
            ) const;

            size_t solves() const;

            void log_diagnostics() const;
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
            };

        private:
            mutable size_t m_solves = 0;
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
        /**
         * @brief The list of identified equilibrium groups.
         */
        std::vector<QSEGroup> m_qse_groups;

        /**
         * @brief A set of solvers, one for each QSE group
         */
        std::vector<std::unique_ptr<QSESolver>> m_qse_solvers;
        /**
         * @brief The simplified set of species presented to the solver (the "slow" species).
         */
        std::vector<fourdst::atomic::Species> m_dynamic_species;
        /**
         * @brief Species that are treated as algebraic (in QSE) in the QSE groups.
         */
        std::vector<fourdst::atomic::Species> m_algebraic_species;

        /**
         * @brief Map from species to their calculated abundances in the QSE state.
         */
        std::unordered_map<fourdst::atomic::Species, double> m_algebraic_abundances;

        /**
         * @brief Indices of all species considered active in the current partition (dynamic + algebraic).
         */
        std::vector<size_t> m_activeSpeciesIndices;
        /**
         * @brief Indices of all reactions involving only active species.
         */
        std::vector<size_t> m_activeReactionIndices;

        mutable std::unordered_map<uint64_t, fourdst::composition::Composition> m_composition_cache;

        SUNContext m_sun_ctx = nullptr; // TODO: initialize and safely destroy this

    private:
        /**
         * @brief Partitions the network by timescale.
         *
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
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Validates candidate QSE groups using flux analysis.
         *
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
            const std::vector<QSEGroup> &candidate_groups,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Solves for the QSE abundances of the algebraic species in a given state.
         *
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
        fourdst::composition::Composition solveQSEAbundances(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Identifies the pool with the slowest mean timescale.
         *
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
            const std::vector<std::vector<fourdst::atomic::Species>>& pools,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Builds a connectivity graph from a species pool.
         *
         * @param species_pool A vector of species indices representing a species pool.
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
            const std::vector<fourdst::atomic::Species>& species_pool
        ) const;

        /**
         * @brief Constructs candidate QSE groups from connected timescale pools.
         *
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
            const std::vector<std::vector<fourdst::atomic::Species>>& candidate_pools,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;

        /**
         * @brief Analyzes the connectivity of timescale pools.
         *
         * @param timescale_pools A vector of vectors of species indices, where each inner vector
         *                        represents a timescale pool.
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
            const std::vector<std::vector<fourdst::atomic::Species>> &timescale_pools
        ) const;

        std::vector<QSEGroup> pruneValidatedGroups(
            const std::vector<QSEGroup> &groups,
            const std::vector<reaction::ReactionSet> &groupReactions,
            const fourdst::composition::Composition &comp,
            double T9,
            double rho
        ) const;
    };
}

