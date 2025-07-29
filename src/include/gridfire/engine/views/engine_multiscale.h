#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_graph.h"

#include "unsupported/Eigen/NonLinearOptimization"

namespace gridfire {
    /**
     * @brief Configuration struct for the QSE cache.
     *
     * @purpose This struct defines the tolerances used to determine if a QSE cache key
     * is considered a hit. It allows for tuning the sensitivity of the cache.
     *
     * @how It works by providing binning widths for temperature, density, and abundances.
     * When a `QSECacheKey` is created, it uses these tolerances to discretize the
     * continuous physical values into bins. If two sets of conditions fall into the
     * same bins, they will produce the same hash and be considered a cache hit.
     *
     * @par Usage Example:
     * Although not typically set by the user directly, the `QSECacheKey` uses this
     * internally. A smaller tolerance (e.g., `T9_tol = 1e-4`) makes the cache more
     * sensitive, leading to more frequent re-partitions, while a larger tolerance
     * (`T9_tol = 1e-2`) makes it less sensitive.
     */
    struct QSECacheConfig {
        double T9_tol; ///< Absolute tolerance to produce the same hash for T9.
        double rho_tol; ///< Absolute tolerance to produce the same hash for rho.
        double Yi_tol; ///< Absolute tolerance to produce the same hash for species abundances.
    };

    /**
     * @brief Key struct for the QSE abundance cache.
     *
     * @purpose This struct is used as the key for the QSE abundance cache (`m_qse_abundance_cache`)
     * within the `MultiscalePartitioningEngineView`. Its primary goal is to avoid
     * expensive re-partitioning and QSE solves for thermodynamic conditions that are
     * "close enough" to previously computed ones.
     *
     * @how It works by storing the temperature (`m_T9`), density (`m_rho`), and species
     * abundances (`m_Y`). A pre-computed hash is generated in the constructor by
     * calling the `hash()` method. This method discretizes the continuous physical
     * values into bins using the tolerances defined in `QSECacheConfig`. The `operator==`
     * simply compares the pre-computed hash values for fast lookups in the `std::unordered_map`.
     */
    struct QSECacheKey {
        double m_T9;
        double m_rho;
        std::vector<double> m_Y; ///< Note that the ordering of Y must match the dynamic species indices in the view.

        std::size_t m_hash = 0; ///< Precomputed hash value for this key.

        // TODO: We should probably sort out how to adjust these from absolute to relative tolerances.
        QSECacheConfig m_cacheConfig = {
            1e-3, // Default tolerance for T9
            1e-1, // Default tolerance for rho
            1e-3  // Default tolerance for species abundances
        };

        /**
         * @brief Constructs a QSECacheKey.
         *
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param Y Species molar abundances.
         *
         * @post The `m_hash` member is computed and stored.
         */
        QSECacheKey(
            const double T9,
            const double rho,
            const std::vector<double>& Y
        );

        /**
         * @brief Computes the hash value for this key.
         *
         * @return The computed hash value.
         *
         * @how This method combines the hashes of the binned temperature, density, and
         * each species abundance. The `bin()` static method is used for discretization.
         */
        size_t hash() const;

        /**
         * @brief Converts a value to a discrete bin based on a tolerance.
         * @param value The value to bin.
         * @param tol The tolerance (bin width) to use for binning.
         * @return The bin number as a long integer.
         *
         * @how The algorithm is `floor(value / tol)`.
         */
        static long bin(double value, double tol);

        /**
         * @brief Equality operator for QSECacheKey.
         * @param other The other QSECacheKey to compare to.
         * @return True if the pre-computed hashes are equal, false otherwise.
         */
        bool operator==(const QSECacheKey& other) const;

    };
}

// Needs to be in this order (splitting gridfire namespace up) to avoid some issues with forward declarations and the () operator.
namespace std {
    template <>
    struct hash<gridfire::QSECacheKey> {
        /**
         * @brief Computes the hash of a QSECacheKey for use in `std::unordered_map`.
         * @param key The QSECacheKey to hash.
         * @return The pre-computed hash value of the key.
         */
        size_t operator()(const gridfire::QSECacheKey& key) const noexcept {
            // The hash is pre-computed, so we just return it.
            return key.m_hash;
        }
    };
} // namespace std

namespace gridfire {
    /**
     * @class MultiscalePartitioningEngineView
     * @brief An engine view that partitions the reaction network into multiple groups based on timescales.
     *
     * @purpose This class is designed to accelerate the integration of stiff nuclear reaction networks.
     * It identifies species that react on very short timescales ("fast" species) and treats them
     * as being in Quasi-Steady-State Equilibrium (QSE). Their abundances are solved for algebraically,
     * removing their stiff differential equations from the system. The remaining "slow" or "dynamic"
     * species are integrated normally. This significantly improves the stability and performance of
     * the solver.
     *
     * @how The core logic resides in the `partitionNetwork()` and `equilibrateNetwork()` methods.
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
        explicit MultiscalePartitioningEngineView(GraphEngine& baseEngine);

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
         * @param Y_full Vector of current molar abundances for all species in the base engine.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A `std::expected` containing `StepDerivatives<double>` on success, or a
         *         `StaleEngineError` if the engine's QSE cache does not contain a solution
         *         for the given state.
         *
         * @purpose To compute the time derivatives for the ODE solver. This implementation
         * modifies the derivatives from the base engine to enforce the QSE condition.
         *
         * @how It first performs a lookup in the QSE abundance cache (`m_qse_abundance_cache`).
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
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state.
         *
         * @param Y_full Vector of current molar abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * @purpose To compute the Jacobian matrix required by implicit ODE solvers.
         *
         * @how It first performs a QSE cache lookup. On a hit, it delegates the full Jacobian
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
        void generateJacobianMatrix(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets an entry from the previously generated Jacobian matrix.
         *
         * @param i_full Row index (species index) in the full network.
         * @param j_full Column index (species index) in the full network.
         * @return Value of the Jacobian matrix at (i_full, j_full).
         *
         * @purpose To provide Jacobian entries to an implicit solver.
         *
         * @how This method directly delegates to the base engine's `getJacobianMatrixEntry`.
         *      It does not currently modify the Jacobian to reflect the QSE algebraic constraints,
         *      as these are handled by setting `dY/dt = 0` in `calculateRHSAndEnergy`.
         *
         * @pre `generateJacobianMatrix()` must have been called for the current state.
         */
        [[nodiscard]] double getJacobianMatrixEntry(
            int i_full,
            int j_full
        ) const override;

        /**
         * @brief Generates the stoichiometry matrix for the network.
         *
         * @purpose To prepare the stoichiometry matrix for later queries.
         *
         * @how This method delegates directly to the base engine's `generateStoichiometryMatrix()`.
         *      The stoichiometry is based on the full, unpartitioned network.
         */
        void generateStoichiometryMatrix() override;

        /**
         * @brief Gets an entry from the stoichiometry matrix.
         *
         * @param speciesIndex Index of the species in the full network.
         * @param reactionIndex Index of the reaction in the full network.
         * @return Stoichiometric coefficient for the species in the reaction.
         *
         * @purpose To query the stoichiometric relationship between a species and a reaction.
         *
         * @how This method delegates directly to the base engine's `getStoichiometryMatrixEntry()`.
         *
         * @pre `generateStoichiometryMatrix()` must have been called.
         */
        [[nodiscard]] int getStoichiometryMatrixEntry(
            int speciesIndex,
            int reactionIndex
        ) const override;

        /**
         * @brief Calculates the molar reaction flow for a given reaction.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param Y_full Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * @purpose To compute the net rate of a single reaction.
         *
         * @how It first checks the QSE cache. On a hit, it retrieves the cached equilibrium
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
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the set of logical reactions in the network.
         *
         * @return A const reference to the `LogicalReactionSet` from the base engine,
         *         containing all reactions in the full network.
         */
        [[nodiscard]] const reaction::LogicalReactionSet & getNetworkReactions() const override;

        /**
         * @brief Sets the set of logical reactions in the network.
         *
         * @param reactions The set of logical reactions to use.
         *
         * @purpose To modify the reaction network.
         *
         * @how This operation is not supported by the `MultiscalePartitioningEngineView` as it
         *      would invalidate the partitioning logic. It logs a critical error and throws an
         *      exception. Network modifications should be done on the base engine before it is
         *      wrapped by this view.
         *
         * @throws exceptions::UnableToSetNetworkReactionsError Always.
         */
        void setNetworkReactions(
            const reaction::LogicalReactionSet &reactions
        ) override;

        /**
         * @brief Computes timescales for all species in the network.
         *
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A `std::expected` containing a map from `Species` to their characteristic
         *         timescales (s) on success, or a `StaleEngineError` on failure.
         *
         * @purpose To get the characteristic timescale `Y / (dY/dt)` for each species.
         *
         * @how It delegates the calculation to the base engine. For any species identified
         *      as algebraic (in QSE), it manually sets their timescale to 0.0 to signify
         *      that they equilibrate instantaneously on the timescale of the solver.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @throws StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes destruction timescales for all species in the network.
         *
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A `std::expected` containing a map from `Species` to their characteristic
         *         destruction timescales (s) on success, or a `StaleEngineError` on failure.
         *
         * @purpose To get the timescale for species destruction, which is used as the primary
         *          metric for network partitioning.
         *
         * @how It delegates the calculation to the base engine. For any species identified
         *      as algebraic (in QSE), it manually sets their timescale to 0.0.
         *
         * @pre The engine must have a valid QSE cache entry for the given state.
         * @throws StaleEngineError If the QSE cache misses.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Updates the internal state of the engine, performing partitioning and QSE equilibration.
         *
         * @param netIn A struct containing the current network input: temperature, density, and composition.
         * @return The new composition after QSE species have been brought to equilibrium.
         *
         * @purpose This is the main entry point for preparing the multiscale engine for use. It
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
         * @purpose To determine if `update()` needs to be called.
         *
         * @how It creates a `QSECacheKey` from the `netIn` data and checks for its
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
         * @how This method delegates directly to the base engine's `setScreeningModel()`.
         */
        void setScreeningModel(
            screening::ScreeningType model
        ) override;

        /**
         * @brief Gets the current electron screening model.
         *
         * @return The currently active screening model type.
         *
         * @how This method delegates directly to the base engine's `getScreeningModel()`.
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        /**
         * @brief Gets the base engine.
         *
         * @return A const reference to the base engine.
         */
        const DynamicEngine & getBaseEngine() const override;

        /**
         * @brief Analyzes the connectivity of timescale pools.
         *
         * @param timescale_pools A vector of vectors of species indices, where each inner vector
         *                        represents a timescale pool.
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of vectors of species indices, where each inner vector represents a
         *         single connected component.
         *
         * @purpose To merge timescale pools that are strongly connected by reactions, forming
         *          cohesive groups for QSE analysis.
         *
         * @how For each pool, it builds a reaction connectivity graph using `buildConnectivityGraph`.
         *      It then finds the connected components within that graph using a Breadth-First Search (BFS).
         *      The resulting components from all pools are collected and returned.
         */
        std::vector<std::vector<size_t>> analyzeTimescalePoolConnectivity(
            const std::vector<std::vector<size_t>> &timescale_pools,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const;

        /**
         * @brief Partitions the network into dynamic and algebraic (QSE) groups based on timescales.
         *
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * @purpose To perform the core partitioning logic that identifies which species are "fast"
         * (and can be treated algebraically) and which are "slow" (and must be integrated dynamically).
         *
         * @how
         * 1.  **`partitionByTimescale`**: Gets species destruction timescales from the base engine,
         *     sorts them, and looks for large gaps to create timescale "pools".
         * 2.  **`identifyMeanSlowestPool`**: The pool with the slowest average timescale is designated
         *     as the core set of dynamic species.
         * 3.  **`analyzeTimescalePoolConnectivity`**: The other (faster) pools are analyzed for
         *     reaction connectivity to form cohesive groups.
         * 4.  **`constructCandidateGroups`**: These connected groups are processed to identify "seed"
         *     species (dynamic species that feed the group) and "algebraic" species (the rest).
         * 5.  **`validateGroupsWithFluxAnalysis`**: The groups are validated by ensuring their internal
         *     reaction flux is much larger than the flux connecting them to the outside network.
         *
         * @pre The input state (Y, T9, rho) must be a valid physical state.
         * @post The internal member variables `m_qse_groups`, `m_dynamic_species`, and
         *       `m_algebraic_species` (and their index maps) are populated with the results of the
         *       partitioning.
         */
        void partitionNetwork(
            const std::vector<double>& Y,
            double T9,
            double rho
        );

        /**
         * @brief Partitions the network based on timescales from a `NetIn` struct.
         *
         * @param netIn A struct containing the current network input.
         *
         * @purpose A convenience overload for `partitionNetwork`.
         *
         * @how It unpacks the `netIn` struct into `Y`, `T9`, and `rho` and then calls the
         *      primary `partitionNetwork` method.
         */
        void partitionNetwork(
            const NetIn& netIn
        );

        /**
         * @brief Exports the network to a DOT file for visualization.
         *
         * @param filename The name of the DOT file to create.
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * @purpose To visualize the partitioned network graph.
         *
         * @how This method delegates the DOT file export to the base engine. It does not
         *      currently add any partitioning information to the output graph.
         */
        void exportToDot(
            const std::string& filename,
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const;

        /**
         * @brief Gets the index of a species in the full network.
         *
         * @param species The species to get the index of.
         * @return The index of the species in the base engine's network.
         *
         * @how This method delegates directly to the base engine's `getSpeciesIndex()`.
         */
        [[nodiscard]] int getSpeciesIndex(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Maps a `NetIn` struct to a molar abundance vector for the full network.
         *
         * @param netIn A struct containing the current network input.
         * @return A vector of molar abundances corresponding to the species order in the base engine.
         *
         * @how This method delegates directly to the base engine's `mapNetInToMolarAbundanceVector()`.
         */
        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

        /**
         * @brief Primes the engine with a specific species.
         *
         * @param netIn A struct containing the current network input.
         * @return A `PrimingReport` struct containing information about the priming process.
         *
         * @purpose To prepare the network for ignition or specific pathway studies.
         *
         * @how This method delegates directly to the base engine's `primeEngine()`. The
         *      multiscale view does not currently interact with the priming process.
         */
        [[nodiscard]] PrimingReport primeEngine(const NetIn &netIn) override;

        /**
         * @brief Gets the fast species in the network.
         *
         * @return A vector of species identified as "fast" or "algebraic" by the partitioning.
         *
         * @purpose To allow external queries of the partitioning results.
         *
         * @how It returns a copy of the `m_algebraic_species` member vector.
         *
         * @pre `partitionNetwork()` must have been called.
         */
        [[nodiscard]] std::vector<fourdst::atomic::Species> getFastSpecies() const;
        /**
         * @brief Gets the dynamic species in the network.
         *
         * @return A const reference to the vector of species identified as "dynamic" or "slow".
         *
         * @purpose To allow external queries of the partitioning results.
         *
         * @how It returns a const reference to the `m_dynamic_species` member vector.
         *
         * @pre `partitionNetwork()` must have been called.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getDynamicSpecies() const;

        /**
         * @brief Equilibrates the network by partitioning and solving for QSE abundances.
         *
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A new composition object with the equilibrated abundances.
         *
         * @purpose A convenience method to run the full QSE analysis and get an equilibrated
         * composition object as a result.
         *
         * @how It first calls `partitionNetwork()` with the given state to define the QSE groups.
         * Then, it calls `solveQSEAbundances()` to compute the new equilibrium abundances for the
         * algebraic species. Finally, it packs the resulting full abundance vector into a new
         * `fourdst::composition::Composition` object and returns it.
         *
         * @pre The input state (Y, T9, rho) must be a valid physical state.
         * @post The engine's internal partition is updated. A new composition object is returned.
         */
        fourdst::composition::Composition equilibrateNetwork(
            const std::vector<double> &Y,
            double T9,
            double rho
        );

        /**
         * @brief Equilibrates the network using QSE from a `NetIn` struct.
         *
         * @param netIn A struct containing the current network input.
         * @return The equilibrated composition.
         *
         * @purpose A convenience overload for `equilibrateNetwork`.
         *
         * @how It unpacks the `netIn` struct into `Y`, `T9`, and `rho` and then calls the
         *      primary `equilibrateNetwork` method.
         */
        fourdst::composition::Composition equilibrateNetwork(
            const NetIn &netIn
        );


    private:
        /**
         * @brief Struct representing a QSE group.
         *
         * @purpose A container to hold all information about a set of species that are potentially
         * in quasi-steady-state equilibrium with each other.
         */
        struct QSEGroup {
            std::set<size_t> species_indices; ///< Indices of all species in this group.
            bool is_in_equilibrium = false;      ///< Flag set by flux analysis.
            std::set<size_t> algebraic_indices; ///< Indices of algebraic species in this group.
            std::set<size_t> seed_indices; ///< Indices of dynamic species in this group.
            double mean_timescale; ///< Mean timescale of the group.

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
        };

        /**
         * @brief Functor for solving QSE abundances using Eigen's nonlinear optimization.
         *
         * @purpose This struct provides the objective function (`operator()`) and its Jacobian
         * (`df`) to Eigen's Levenberg-Marquardt solver. The goal is to find the abundances
         * of algebraic species that make their time derivatives (`dY/dt`) equal to zero.
         *
         * @how
         * - **`operator()`**: Takes a vector `v_qse` (scaled abundances of algebraic species) as input.
         *   It constructs a full trial abundance vector `y_trial`, calls the base engine's
         *   `calculateRHSAndEnergy`, and returns the `dY/dt` values for the algebraic species.
         *   The solver attempts to drive this return vector to zero.
         * - **`df`**: Computes the Jacobian of the objective function. It calls the base engine's
         *   `generateJacobianMatrix` and extracts the sub-matrix corresponding to the algebraic
         *   species. It applies the chain rule to account for the `asinh` scaling used on the
         *   abundances.
         *
         * The abundances are scaled using `asinh` to handle the large dynamic range and ensure positivity.
         */
        struct EigenFunctor {
            using InputType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
            using OutputType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
            using JacobianType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
            enum {
                InputsAtCompileTime = Eigen::Dynamic,
                ValuesAtCompileTime = Eigen::Dynamic
            };

            /**
             * @brief Pointer to the MultiscalePartitioningEngineView instance.
             */
            MultiscalePartitioningEngineView* m_view;
            /**
             * @brief Indices of the species to solve for in the QSE group.
             */
            const std::vector<size_t>& m_qse_solve_indices;
            /**
             * @brief Initial abundances of all species in the full network.
             */
            const std::vector<double>& m_Y_full_initial;
            /**
             * @brief Temperature in units of 10^9 K.
             */
            const double m_T9;
            /**
             * @brief Density in g/cm^3.
             */
            const double m_rho;
            /**
             * @brief Scaling factors for the species abundances, used to improve solver stability.
             */
            const Eigen::VectorXd& m_Y_scale;

            /**
             * @brief Constructs an EigenFunctor.
             *
             * @param view The MultiscalePartitioningEngineView instance.
             * @param qse_solve_indices Indices of the species to solve for in the QSE group.
             * @param Y_full_initial Initial abundances of all species.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             * @param Y_scale Scaling factors for the species abundances.
             */
            EigenFunctor(
                MultiscalePartitioningEngineView& view,
                const std::vector<size_t>& qse_solve_indices,
                const std::vector<double>& Y_full_initial,
                const double T9,
                const double rho,
                const Eigen::VectorXd& Y_scale
            ) :
            m_view(&view),
            m_qse_solve_indices(qse_solve_indices),
            m_Y_full_initial(Y_full_initial),
            m_T9(T9),
            m_rho(rho),
            m_Y_scale(Y_scale) {}

            /**
             * @brief Gets the number of output values from the functor (size of the residual vector).
             * @return The number of algebraic species being solved.
             */
            [[nodiscard]] int values() const { return m_qse_solve_indices.size(); }
            /**
             * @brief Gets the number of input values to the functor (size of the variable vector).
             * @return The number of algebraic species being solved.
             */
            [[nodiscard]] int inputs() const { return m_qse_solve_indices.size(); }

            /**
             * @brief Evaluates the functor's residual vector `f_qse = dY_alg/dt`.
             * @param v_qse The input vector of scaled algebraic abundances.
             * @param f_qse The output residual vector.
             * @return 0 on success.
             */
            int operator()(const InputType& v_qse, OutputType& f_qse) const;
            /**
             * @brief Evaluates the Jacobian of the functor, `J_qse = d(f_qse)/d(v_qse)`.
             * @param v_qse The input vector of scaled algebraic abundances.
             * @param J_qse The output Jacobian matrix.
             * @return 0 on success.
             */
            int df(const InputType& v_qse, JacobianType& J_qse) const;
        };


        /**
         * @brief Struct for tracking cache statistics.
         *
         * @purpose A simple utility to monitor the performance of the QSE cache by counting
         * hits and misses for various engine operations.
         */
        struct CacheStats {
            enum class operators {
                CalculateRHSAndEnergy,
                GenerateJacobianMatrix,
                CalculateMolarReactionFlow,
                GetSpeciesTimescales,
                GetSpeciesDestructionTimescales,
                Other,
                All
            };

            /**
             * @brief Map from operators to their string names for logging.
             */
            std::map<operators, std::string> operatorsNameMap = {
                {operators::CalculateRHSAndEnergy, "calculateRHSAndEnergy"},
                {operators::GenerateJacobianMatrix, "generateJacobianMatrix"},
                {operators::CalculateMolarReactionFlow, "calculateMolarReactionFlow"},
                {operators::GetSpeciesTimescales, "getSpeciesTimescales"},
                {operators::GetSpeciesDestructionTimescales, "getSpeciesDestructionTimescales"},
                {operators::Other, "other"}
            };

            /**
             * @brief Total number of cache hits.
             */
            size_t m_hit = 0;
            /**
             * @brief Total number of cache misses.
             */
            size_t m_miss = 0;

            /**
             * @brief Map from operators to the number of cache hits for that operator.
             */
            std::map<operators, size_t> m_operatorHits = {
                {operators::CalculateRHSAndEnergy, 0},
                {operators::GenerateJacobianMatrix, 0},
                {operators::CalculateMolarReactionFlow, 0},
                {operators::GetSpeciesTimescales, 0},
                {operators::GetSpeciesDestructionTimescales, 0},
                {operators::Other, 0}
            };

            /**
             * @brief Map from operators to the number of cache misses for that operator.
             */
            std::map<operators, size_t> m_operatorMisses = {
                {operators::CalculateRHSAndEnergy, 0},
                {operators::GenerateJacobianMatrix, 0},
                {operators::CalculateMolarReactionFlow, 0},
                {operators::GetSpeciesTimescales, 0},
                {operators::GetSpeciesDestructionTimescales, 0},
                {operators::Other, 0}
            };

            /**
             * @brief Increments the hit counter for a given operator.
             * @param op The operator that resulted in a cache hit.
             * @throws std::invalid_argument if `op` is `All`.
             */
            void hit(const operators op=operators::Other);
            /**
             * @brief Increments the miss counter for a given operator.
             * @param op The operator that resulted in a cache miss.
             * @throws std::invalid_argument if `op` is `All`.
             */
            void miss(const operators op=operators::Other);

            /**
             * @brief Gets the number of hits for a specific operator or all operators.
             * @param op The operator to get the number of hits for. Defaults to `All`.
             * @return The number of hits.
             */
            [[nodiscard]] size_t hits(const operators op=operators::All) const;
            /**
             * @brief Gets the number of misses for a specific operator or all operators.
             * @param op The operator to get the number of misses for. Defaults to `All`.
             * @return The number of misses.
             */
            [[nodiscard]] size_t misses(const operators op=operators::All) const;
        };


    private:
        /**
         * @brief Logger instance for logging messages.
         */
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        /**
         * @brief The base engine to which this view delegates calculations.
         */
        GraphEngine& m_baseEngine;
        /**
         * @brief The list of identified equilibrium groups.
         */
        std::vector<QSEGroup> m_qse_groups;
        /**
         * @brief The simplified set of species presented to the solver (the "slow" species).
         */
        std::vector<fourdst::atomic::Species> m_dynamic_species;
        /**
         * @brief Indices mapping the dynamic species back to the base engine's full species list.
         */
        std::vector<size_t> m_dynamic_species_indices;
        /**
         * @brief Species that are treated as algebraic (in QSE) in the QSE groups.
         */
        std::vector<fourdst::atomic::Species> m_algebraic_species;
        /**
         * @brief Indices of algebraic species in the full network.
         */
        std::vector<size_t> m_algebraic_species_indices;

        /**
         * @brief Indices of all species considered active in the current partition (dynamic + algebraic).
         */
        std::vector<size_t> m_activeSpeciesIndices;
        /**
         * @brief Indices of all reactions involving only active species.
         */
        std::vector<size_t> m_activeReactionIndices;

        // TODO: Enhance the hashing for the cache to consider not just T and rho but also the current abundance in some careful way that automatically ignores small changes (i.e. network should only be repartitioned sometimes)
        /**
         * @brief Cache for QSE abundances based on T9, rho, and Y.
         *
         * @purpose This is the core of the caching mechanism. It stores the results of QSE solves
         * to avoid re-computation. The key is a `QSECacheKey` which hashes the thermodynamic
         * state, and the value is the vector of solved molar abundances for the algebraic species.
         */
        mutable std::unordered_map<QSECacheKey, std::vector<double>> m_qse_abundance_cache;
        /**
         * @brief Statistics for the QSE abundance cache.
         */
        mutable CacheStats m_cacheStats;




    private:
        /**
         * @brief Partitions the network by timescale.
         *
         * @param Y_full Vector of current molar abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of vectors of species indices, where each inner vector represents a
         *         timescale pool.
         *
         * @purpose To group species into "pools" based on their destruction timescales.
         *
         * @how It retrieves all species destruction timescales from the base engine, sorts them,
         *      and then iterates through the sorted list, creating a new pool whenever it detects
         *      a gap between consecutive timescales that is larger than a predefined threshold
         *      (e.g., a factor of 100).
         */
        std::vector<std::vector<size_t>> partitionByTimescale(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const;

        /**
         * @brief Builds a connectivity graph from a set of fast reaction indices.
         *
         * @param fast_reaction_indices A set of indices for reactions considered "fast".
         * @return An unordered map representing the adjacency list of the connectivity graph,
         *         where keys are species indices and values are vectors of connected species indices.
         *
         * @purpose To represent the reaction pathways among a subset of reactions.
         *
         * @how It iterates through the specified fast reactions. For each reaction, it creates
         *      a two-way edge in the graph between every reactant and every product, signifying
         *      that mass can flow between them.
         */
        std::unordered_map<size_t, std::vector<size_t>> buildConnectivityGraph(
            const std::unordered_set<size_t> &fast_reaction_indices
        ) const;

        /**
         * @brief Validates candidate QSE groups using flux analysis.
         *
         * @param candidate_groups A vector of candidate QSE groups.
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of validated QSE groups that meet the flux criteria.
         *
         * @purpose To ensure that a candidate QSE group is truly in equilibrium by checking that
         *          the reaction fluxes *within* the group are much larger than the fluxes
         *          *leaving* the group.
         *
         * @how For each candidate group, it calculates the sum of all internal reaction fluxes and
         *      the sum of all external (bridge) reaction fluxes. If the ratio of internal to external
         *      flux exceeds a configurable threshold, the group is considered valid and is added
         *      to the returned vector.
         */
        std::vector<QSEGroup> validateGroupsWithFluxAnalysis(
            const std::vector<QSEGroup> &candidate_groups,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const;

        /**
         * @brief Solves for the QSE abundances of the algebraic species in a given state.
         *
         * @param Y_full Vector of current molar abundances for all species in the base engine.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of molar abundances for the algebraic species.
         *
         * @purpose To find the equilibrium abundances of the algebraic species that satisfy
         *          the QSE conditions.
         *
         * @how It uses the Levenberg-Marquardt algorithm via Eigen's `LevenbergMarquardt` class.
         *      The problem is defined by the `EigenFunctor` which computes the residuals and
         *      Jacobian for the QSE equations.
         *
         * @pre The input state (Y_full, T9, rho) must be a valid physical state.
         * @post The algebraic species in the QSE cache are updated with the new equilibrium abundances.
         */
        std::vector<double> solveQSEAbundances(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        );

        /**
         * @brief Identifies the pool with the slowest mean timescale.
         *
         * @param pools A vector of vectors of species indices, where each inner vector represents a
         *              timescale pool.
         * @param Y Vector of current molar abundances for the full network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return The index of the pool with the largest (slowest) mean destruction timescale.
         *
         * @purpose To identify the core set of dynamic species that will not be part of any QSE group.
         *
         * @how It calculates the geometric mean of the destruction timescales for all species in each
         *      pool and returns the index of the pool with the maximum mean timescale.
         */
        size_t identifyMeanSlowestPool(
            const std::vector<std::vector<size_t>>& pools,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const;

        /**
         * @brief Builds a connectivity graph from a species pool.
         *
         * @param species_pool A vector of species indices representing a species pool.
         * @return An unordered map representing the adjacency list of the connectivity graph.
         *
         * @purpose To find reaction connections within a specific group of species.
         *
         * @how It iterates through all reactions in the base engine. If a reaction involves
         *      at least two distinct species from the input `species_pool` (one as a reactant
         *      and one as a product), it adds edges between all reactants and products from
         *      that reaction that are also in the pool.
         */
        std::unordered_map<size_t, std::vector<size_t>> buildConnectivityGraph(
            const std::vector<size_t>& species_pool
        ) const;

        /**
         * @brief Constructs candidate QSE groups from connected timescale pools.
         *
         * @param candidate_pools A vector of vectors of species indices, where each inner vector
         *                        represents a connected pool of species with similar fast timescales.
         * @param Y Vector of current molar abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A vector of `QSEGroup` structs, ready for flux validation.
         *
         * @how For each input pool, it identifies "bridge" reactions that connect the pool to
         * species outside the pool. The reactants of these bridge reactions that are *not* in the
         * pool are identified as "seed" species. The original pool members are the "algebraic"
         * species. It then bundles the seed and algebraic species into a `QSEGroup` struct.
         *
         * @pre The `candidate_pools` should be connected components from `analyzeTimescalePoolConnectivity`.
         * @post A list of candidate `QSEGroup` objects is returned.
         */
        std::vector<QSEGroup> constructCandidateGroups(
            const std::vector<std::vector<size_t>>& candidate_pools,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const;
    };
}

