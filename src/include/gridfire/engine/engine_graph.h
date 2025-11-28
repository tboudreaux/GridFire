#pragma once

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/composition/composition.h"
#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "gridfire/types/types.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/partition/partition_abstract.h"
#include "gridfire/engine/procedures/construction.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <ranges>
#include <functional>


#include "cppad/cppad.hpp"
#include "cppad/utility/sparse_rc.hpp"
#include "cppad/speed/sparse_jac_fun.hpp"

#include "gridfire/reaction/weak/weak_interpolator.h"
#include "gridfire/reaction/weak/weak_rate_library.h"

// PERF: The function getNetReactionStoichiometry returns a map of species to their stoichiometric coefficients for a given reaction.
//       this makes extra copies of the species, which is not ideal and could be optimized further.
//       Even more relevant is the member m_reactionIDMap which makes copies of a REACLIBReaction for each reaction ID.
//       REACLIBReactions are quite large data structures, so this could be a performance bottleneck.
namespace gridfire::engine {
    /**
     * @brief Alias for CppAD AD type for double precision.
     *
     * This alias simplifies the use of the CppAD automatic differentiation type.
     */
    typedef CppAD::AD<double> ADDouble;

    using fourdst::config::Config;
    using fourdst::logging::LogManager;
    using fourdst::constant::Constants;

    /**
     * @brief Minimum density threshold below which reactions are ignored.
     *
     * Reactions are not calculated if the density falls below this threshold.
     * This helps to improve performance by avoiding unnecessary calculations
     * in very low-density regimes.
     */
    static constexpr double MIN_DENSITY_THRESHOLD = 1e-18;

    /**
     * @brief Minimum abundance threshold below which species are ignored.
     *
     * Species with abundances below this threshold are treated as zero in
     * reaction rate calculations. This helps to improve performance by
     * avoiding unnecessary calculations for trace species.
     */
    static constexpr double MIN_ABUNDANCE_THRESHOLD = 1e-18;

    /**
     * @brief Minimum value for Jacobian matrix entries.
     *
     * Jacobian matrix entries with absolute values below this threshold are
     * treated as zero to maintain sparsity and improve performance.
     */
    static constexpr double MIN_JACOBIAN_THRESHOLD = 1e-24;




    /**
     * @class GraphEngine
     * @brief A reaction network engine that uses a graph-based representation.
     *
     * The GraphEngine class implements the DynamicEngine interface using a
     * graph-based representation of the reaction network. It uses sparse
     * matrices for efficient storage and computation of the stoichiometry
     * and Jacobian matrices. Automatic differentiation (AD) is used to
     * calculate the Jacobian matrix.
     *
     * The engine supports:
     *   - Calculation of the right-hand side (dY/dt) and energy generation rate.
     *   - Generation and access to the Jacobian matrix.
     *   - Generation and access to the stoichiometry matrix.
     *   - Calculation of molar reaction flows.
     *   - Access to the set of logical reactions in the network.
     *   - Computation of timescales for each species.
     *   - Exporting the network to DOT and CSV formats for visualization and analysis.
     *
     * @implements DynamicEngine
     *
     * @see engine_abstract.h
     */
    class GraphEngine final : public DynamicEngine{
    public:
        /**
         * @brief Constructs a GraphEngine from a composition.
         *
         * @param composition The composition of the material.
         *
         * This constructor builds the reaction network from the given composition
         * using the `build_reaclib_nuclear_network` function.
         *
         * @see build_reaclib_nuclear_network
         */
        explicit GraphEngine(
            const fourdst::composition::Composition &composition,
            BuildDepthType = NetworkBuildDepth::Full
        );

        explicit GraphEngine(
            const fourdst::composition::Composition &composition,
            const partition::PartitionFunction& partitionFunction,
            BuildDepthType buildDepth = NetworkBuildDepth::Full
        );

        explicit GraphEngine(
            const fourdst::composition::Composition &composition,
            const partition::PartitionFunction& partitionFunction,
            BuildDepthType buildDepth,
            NetworkConstructionFlags reactionTypes
        );

        /**
         * @brief Constructs a GraphEngine from a set of reactions.
         *
         * @param reactions The set of reactions to use in the network.
         *
         * This constructor uses the given set of reactions to construct the
         * reaction network.
         */
        explicit GraphEngine(const reaction::ReactionSet &reactions);

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation rate.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<double> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate for the current state.
         *
         * @see StepDerivatives
         */
        [[nodiscard]] std::expected<StepDerivatives<double>, engine::EngineStatus> calculateRHSAndEnergy(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation rate for a subset of reactions.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param activeReactions The set of reactions to include in the calculation.
         * @return StepDerivatives<double> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate considering only the specified
         * subset of reactions. This allows for flexible calculations with
         * different reaction sets without modifying the engine's internal state.
         *
         * @see StepDerivatives
         */
        [[nodiscard]] std::expected<StepDerivatives<double>, EngineStatus> calculateRHSAndEnergy(
            const fourdst::composition::CompositionAbstract& comp,
            double T9,
            double rho,
            const reaction::ReactionSet &activeReactions
        ) const;

        /**
         * @brief Calculates the derivatives of the energy generation rate with respect to temperature and density
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return EnergyDerivatives struct containing the derivatives.
         *
         * This method computes the partial derivatives of the specific nuclear
         * energy generation rate with respect to temperature (∂ε/∂T) and
         * density (∂ε/∂ρ)
         *
         * @see EnergyDerivatives
         */
        [[nodiscard]] EnergyDerivatives calculateEpsDerivatives(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Calculates the derivatives of the energy generation rate with respect to temperature and density for a subset of reactions
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param activeReactions The set of reactions to include in the calculation.
         * @return EnergyDerivatives struct containing the derivatives.
         *
         * This method computes the partial derivatives of the specific nuclear
         * energy generation rate with respect to temperature (∂ε/∂T) and
         * density (∂ε/∂ρ) considering
         * only the specified subset of reactions. This allows for flexible
         * calculations with different reaction sets without modifying the
         * engine's internal state.
         *
         * @see EnergyDerivatives
         */
        [[nodiscard]] EnergyDerivatives calculateEpsDerivatives(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const reaction::ReactionSet &activeReactions
        ) const;

        /**
         * @brief Generates the Jacobian matrix for the current state.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * This method computes and stores the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state using automatic differentiation. The matrix can
         * then be accessed via `getJacobianMatrixEntry()`.
         *
         * @see getJacobianMatrixEntry()
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state with a specified set of active species.
         *        generally this will be much faster than the full matrix generation. Here we use forward mode
         *        to generate the Jacobian only for the active species.
         * @param comp The Composition object containing current abundances.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The density in g/cm^3.
         * @param activeSpecies A vector of Species objects representing the active species.
         *
         * @see getJacobianMatrixEntry()
         * @see generateJacobianMatrix()
         */
        [[nodiscard]] NetworkJacobian generateJacobianMatrix(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const std::vector<fourdst::atomic::Species>& activeSpecies
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state with a specified sparsity pattern.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param sparsityPattern The sparsity pattern to use for the Jacobian matrix.
         *
         * This method computes and stores the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state using automatic differentiation, taking into
         * account the provided sparsity pattern. The matrix can then be accessed
         * via `getJacobianMatrixEntry()`.
         *
         * @see getJacobianMatrixEntry()
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
         * This method computes and stores the stoichiometry matrix,
         * which encodes the net change of each species in each reaction.
         */
        void generateStoichiometryMatrix() override;

        /**
         * @brief Calculates the molar reaction flow for a given reaction.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the net rate at which the given reaction proceeds
         * under the current state.
         *
         */
        [[nodiscard]] double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Gets the list of species in the network.
         * @return Vector of Species objects representing all network species.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

        /**
         * @brief Gets the set of logical reactions in the network.
         * @return Reference to the LogicalReactionSet containing all reactions.
         */
        [[nodiscard]] const reaction::ReactionSet& getNetworkReactions() const override;

        /**
         * @brief Sets the reactions for the network.
         *
         * @param reactions The set of reactions to use in the network.
         *
         * This method replaces the current set of reactions in the network
         * with the provided set. It marks the engine as stale, requiring
         * regeneration of matrices and recalculation of rates.
         */
        void setNetworkReactions(const reaction::ReactionSet& reactions) override;

        /**
         * @brief Gets the net stoichiometry for a given reaction.
         *
         * @param reaction The reaction for which to get the stoichiometry.
         * @return Map of species to their stoichiometric coefficients.
         */
        [[nodiscard]] static std::unordered_map<fourdst::atomic::Species, int> getNetReactionStoichiometry(
            const reaction::Reaction& reaction
        );

        /**
         * @brief Gets an entry from the stoichiometry matrix.
         *
         * @param species Species to look up stoichiometry for.
         * @param reaction Reaction to find.
         * @return Stoichiometric coefficient for the species in the reaction.
         *
         * The stoichiometry matrix must have been generated by `generateStoichiometryMatrix()`.
         *
         * @see generateStoichiometryMatrix()
         */
        [[nodiscard]] int getStoichiometryMatrixEntry(
            const fourdst::atomic::Species& species,
            const reaction::Reaction& reaction
        ) const override;

        /**
         * @brief Computes timescales for all species in the network.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method estimates the timescale for abundance change of each species,
         * which can be used for timestep control or diagnostics.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, engine::EngineStatus>
        getSpeciesTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes timescales for all species in the network considering a subset of reactions.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param activeReactions The set of reactions to include in the calculation.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method estimates the timescale for abundance change of each species,
         * considering only the specified subset of reactions. This allows for flexible
         * calculations with different reaction sets without modifying the engine's internal state.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const reaction::ReactionSet &activeReactions
        ) const;


        /**
         * @brief Computes destruction timescales for all species in the network.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their destruction timescales (s).
         *
         * This method estimates the destruction timescale for each species,
         * which can be useful for understanding reaction flows and equilibrium states.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, engine::EngineStatus>
        getSpeciesDestructionTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Computes destruction timescales for all species in the network considering a subset of reactions.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param activeReactions The set of reactions to include in the calculation.
         * @return Map from Species to their destruction timescales (s).
         *
         * This method estimates the destruction timescale for each species,
         * considering only the specified subset of reactions. This allows for flexible
         * calculations with different reaction sets without modifying the engine's internal state.
         */
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesDestructionTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const reaction::ReactionSet &activeReactions
        ) const;

        /**
         * @brief Updates the state of the network and the composition to be usable for the current network.
         *
         * @details For graph engine all this does is ensure that the returned composition has all the species in the network registered.
         * if a species was already in the composition is will keep its abundance, otherwise it will be added with zero abundance.
         *
         * @param netIn The input netIn to use, this includes the composition, temperature, and density
         *
         * @return The updated composition that includes all species in the network.
         */
        fourdst::composition::Composition update(
            const NetIn &netIn
        ) override;


        /**
         * @brief Checks if the engine view is stale and needs to be updated.
         *
         * @param netIn The current network input (unused).
         * @return True if the view is stale, false otherwise.
         *
         * @deprecated This method is deprecated and will be removed in future versions.
         *             Stale states are returned as part of the results of methods that
         *             require the ability to report them.
         */
        [[deprecated]] bool isStale(
            const NetIn &netIn
        ) override;

        /**
         * @brief Checks if a given species is involved in the network.
         *
         * @param species The species to check.
         * @return True if the species is involved in the network, false otherwise.
         */
        [[nodiscard]] bool involvesSpecies(
            const fourdst::atomic::Species& species
        ) const;

        /**
         * @brief Exports the network to a DOT file for visualization.
         *
         * @param filename The name of the DOT file to create.
         *
         * This method generates a DOT file that can be used to visualize the
         * reaction network as a graph. The DOT file can be converted to a
         * graphical image using Graphviz.
         *
         * @throws std::runtime_error If the file cannot be opened for writing.
         *
         * Example usage:
         * @code
         * engine.exportToDot("network.dot");
         * @endcode
         */
        void exportToDot(
            const std::string& filename
        ) const;

        /**
         * @brief Exports the network to a CSV file for analysis.
         *
         * @param filename The name of the CSV file to create.
         *
         * This method generates a CSV file containing information about the
         * reactions in the network, including the reactants, products, Q-value,
         * and reaction rate coefficients.
         *
         * @throws std::runtime_error If the file cannot be opened for writing.
         *
         * Example usage:
         * @code
         * engine.exportToCSV("network.csv");
         * @endcode
         */
        void exportToCSV(
            const std::string& filename
        ) const;

        /**
         * @brief Sets the electron screening model for reaction rate calculations.
         *
         * @param model The type of screening model to use.
         *
         * This method allows changing the screening model at runtime. Screening corrections
         * account for the electrostatic shielding of nuclei by electrons, which affects
         * reaction rates in dense stellar plasmas.
         */
        void setScreeningModel(
            screening::ScreeningType model
        ) override;

        /**
         * @brief Gets the current electron screening model.
         *
         * @return The currently active screening model type.
         *
         * Example usage:
         * @code
         * screening::ScreeningType currentModel = engine.getScreeningModel();
         * @endcode
         */
        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        /**
         * @brief Sets whether to precompute reaction rates.
         *
         * @param precompute True to enable precomputation, false to disable.
         *
         * This method allows enabling or disabling precomputation of reaction rates
         * for performance optimization. When enabled, reaction rates are computed
         * once and stored for later use.
         *
         * @post If precomputation is enabled, reaction rates will be precomputed and cached.
         *       If disabled, reaction rates will be computed on-the-fly as needed.
         */
        void setPrecomputation(
            bool precompute
        );

        /**
         * @brief Checks if precomputation of reaction rates is enabled.
         *
         * @return True if precomputation is enabled, false otherwise.
         *
         * This method allows checking the current state of precomputation for
         * reaction rates in the engine.
         */
        [[nodiscard]] bool isPrecomputationEnabled() const;

        /**
         * @brief Gets the partition function used for reaction rate calculations.
         *
         * @return Reference to the PartitionFunction object.
         *
         * This method provides access to the partition function used in the engine,
         * which is essential for calculating thermodynamic properties and reaction rates.
         */
        [[nodiscard]] const partition::PartitionFunction& getPartitionFunction() const;

        /**
         * @brief Calculates the reverse rate for a given reaction.
         *
         * @param reaction The reaction for which to calculate the reverse rate.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho
         * @param comp Composition object containing current abundances.
         * @return Reverse rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the reverse rate based on the forward rate and
         * thermodynamic properties of the reaction.
         */
        [[nodiscard]] double calculateReverseRate(
            const reaction::Reaction &reaction,
            double T9,
            double rho,
            const fourdst::composition::CompositionAbstract &comp
        ) const;

        /**
         * @brief Calculates the reverse rate for a two-body reaction.
         *
         * @param reaction The reaction for which to calculate the reverse rate.
         * @param T9 Temperature in units of 10^9 K.
         * @param forwardRate The forward rate of the reaction.
         * @param expFactor Exponential factor for the reaction.
         * @return Reverse rate for the two-body reaction (e.g., mol/g/s).
         *
         * This method computes the reverse rate using the forward rate and
         * thermodynamic properties of the reaction.
         */
        [[nodiscard]] double calculateReverseRateTwoBody(
            const reaction::Reaction &reaction,
            double T9,
            double forwardRate,
            double expFactor
        ) const;

        /**
         * @brief Calculates the derivative of the reverse rate for a two-body reaction with respect to temperature.
         *
         * @param reaction The reaction for which to calculate the derivative.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param comp Composition object containing current abundances.
         * @param reverseRate The reverse rate of the reaction.
         * @return Derivative of the reverse rate with respect to temperature.
         *
         * This method computes the derivative of the reverse rate using automatic differentiation.
         */
        [[nodiscard]] double calculateReverseRateTwoBodyDerivative(
            const reaction::Reaction &reaction,
            double T9,
            double rho,
            const fourdst::composition::Composition& comp,
            double reverseRate
        ) const;

        /**
         * @brief Checks if reverse reactions are enabled.
         *
         * @return True if reverse reactions are enabled, false otherwise.
         *
         * This method allows checking whether the engine is configured to use
         * reverse reactions in its calculations.
         */
        [[nodiscard]] bool isUsingReverseReactions() const;

        /**
         * @brief Sets whether to use reverse reactions in the engine.
         *
         * @param useReverse True to enable reverse reactions, false to disable.
         *
         * This method allows enabling or disabling reverse reactions in the engine.
         * If disabled, only forward reactions will be considered in calculations.
         *
         * @post If reverse reactions are enabled, the engine will consider both
         *       forward and reverse reactions in its calculations. If disabled,
         *       only forward reactions will be considered.
         */
        void setUseReverseReactions(
            bool useReverse
        );

        /**
         * @brief Gets the index of a species in the network.
         *
         * @param species The species for which to get the index.
         * @return Index of the species in the network, or -1 if not found.
         *
         * This method returns the index of the given species in the network's
         * species vector. If the species is not found, it returns -1.
         */
        [[nodiscard]] size_t getSpeciesIndex(
            const fourdst::atomic::Species &species
        ) const override;

        /**
         * @brief Maps the NetIn object to a vector of molar abundances.
         *
         * @param netIn The NetIn object containing the input conditions.
         * @return Vector of molar abundances corresponding to the species in the network.
         *
         * This method converts the NetIn object into a vector of molar abundances
         * for each species in the network, which can be used for further calculations.
         */
        [[deprecated]] [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(
            const NetIn &netIn
        ) const override;

        /**
         * @brief Prepares the engine for calculations with initial conditions.
         *
         * @param netIn The input conditions for the network.
         * @return PrimingReport containing information about the priming process.
         *
         * This method initializes the engine with the provided input conditions,
         * setting up reactions, species, and precomputing necessary data.
         */
        [[nodiscard]] PrimingReport primeEngine(
            const NetIn &netIn
        ) override;

        /**
         * @brief Gets the depth of the network.
         *
         * @return The build depth of the network.
         *
         * This method returns the current build depth of the reaction network,
         * which indicates how many levels of reactions are included in the network.
         */
        [[nodiscard]] BuildDepthType getDepth() const override;

        /**
         * @brief Rebuilds the reaction network based on a new composition.
         *
         * @param comp The new composition to use for rebuilding the network.
         * @param depth The build depth to use for the network.
         *
         * This method rebuilds the reaction network using the provided composition
         * and build depth. It updates all internal data structures accordingly.
         */
        void rebuild(
            const fourdst::composition::CompositionAbstract &comp,
            BuildDepthType depth
        ) override;

        /**
         * @brief This will return the input comp with the molar abundances of any species not registered in that but
         * registered in the engine active species set to 0.0.
         * @note Effectively this method does not change input composition; rather it ensures that all species which
         * can be tracked by an instance of GraphEngine are registered in the composition object.
         * @note If a species is in the input comp but not in the network
         * @param comp Input Composition
         * @param T9
         * @param rho
         * @param T9
         * @param rho
         * @return A new composition where all members of the active species set are registered. And any members not in comp
         * have a molar abundance set to 0.
         * @throws BadCollectionError If the input composition contains species not present in the network species set
         */
        fourdst::composition::Composition collectComposition(const fourdst::composition::CompositionAbstract &comp, double T9, double rho) const override;

        /**
         * @brief Gets the status of a species in the network.
         *
         * @param species The species for which to get the status.
         * @return SpeciesStatus indicating the status of the species.
         *
         * This method checks if the given species is part of the network and
         * returns its status (e.g., Active, Inactive, NotFound).
         */
        [[nodiscard]]
        SpeciesStatus getSpeciesStatus(const fourdst::atomic::Species &species) const override;


    private:
        struct PrecomputedReaction {
            // Forward cacheing
            size_t reaction_index{};
            reaction::ReactionType reaction_type{};
            uint64_t reaction_hash{};
            std::vector<size_t> unique_reactant_indices{};
            std::vector<int> reactant_powers{};
            double symmetry_factor{};
            std::vector<size_t> affected_species_indices{};
            std::vector<int> stoichiometric_coefficients{};

            // Reverse cacheing
            std::vector<size_t> unique_product_indices{}; ///< Unique product indices for reverse reactions.
            std::vector<int> product_powers{}; ///< Powers of each unique product in the reverse reaction.
            double reverse_symmetry_factor{}; ///< Symmetry factor for reverse reactions.
        };

        struct constants {
            const double u = Constants::getInstance().get("u").value; ///< Atomic mass unit in g.
            const double Na = Constants::getInstance().get("N_a").value; ///< Avogadro's number.
            const double c = Constants::getInstance().get("c").value; ///< Speed of light in cm/s.
            const double kB = Constants::getInstance().get("kB").value; ///< Boltzmann constant in erg/K.
            const double MeV_to_erg = Constants::getInstance().get("MeV_to_erg").value; ///< Conversion factor from MeV to erg.
        };

        enum class JacobianMatrixState {
            UNINITIALIZED,
            STALE,
            READY_DENSE,
            READY_SPARSE
        };

        std::unordered_map<JacobianMatrixState, std::string> m_jacobianMatrixStateNameMap = {
            {JacobianMatrixState::UNINITIALIZED, "Uninitialized"},
            {JacobianMatrixState::STALE, "Stale"},
            {JacobianMatrixState::READY_DENSE, "Ready (dense)"},
            {JacobianMatrixState::READY_SPARSE, "Ready (sparse)"},
        };
    private:
        class AtomicReverseRate final : public CppAD::atomic_base<double> {
        public:
            AtomicReverseRate(
                const reaction::Reaction& reaction,
                const GraphEngine& engine
            ):
            atomic_base<double>("AtomicReverseRate"),
            m_reaction(reaction),
            m_engine(engine) {}

            bool forward(
                size_t p,
                size_t q,
                const CppAD::vector<bool>& vx,
                CppAD::vector<bool>& vy,
                const CppAD::vector<double>& tx,
                CppAD::vector<double>& ty
            ) override;
            bool reverse(
                size_t q,
                const CppAD::vector<double>& tx,
                const CppAD::vector<double>& ty,
                CppAD::vector<double>& px,
                const CppAD::vector<double>& py
            ) override;
            bool for_sparse_jac(
                size_t q,
                const CppAD::vector<std::set<size_t>>&r,
                CppAD::vector<std::set<size_t>>& s
            ) override;
            bool rev_sparse_jac(
                size_t q,
                const CppAD::vector<std::set<size_t>>&rt,
                CppAD::vector<std::set<size_t>>& st
            ) override;
            bool for_sparse_jac(
                size_t q,
                const CppAD::vector<bool> &r,
                CppAD::vector<bool> &s,
                const CppAD::vector<double> &x
            ) override;
            bool rev_sparse_jac(
                size_t q,
                const CppAD::vector<bool> &rt,
                CppAD::vector<bool> &st,
                const CppAD::vector<double> &x
            ) override;

        private:
            const reaction::Reaction& m_reaction;
            const GraphEngine& m_engine;
        };
    private:
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        constants m_constants;

        rates::weak::WeakRateInterpolator m_weakRateInterpolator; ///< Interpolator for weak reaction rates.

        reaction::ReactionSet m_reactions; ///< Set of REACLIB reactions in the network.
        std::unordered_map<std::string_view, reaction::Reaction*> m_reactionIDMap; ///< Map from reaction ID to REACLIBReaction. //PERF: This makes copies of REACLIBReaction and could be a performance bottleneck.

        std::vector<fourdst::atomic::Species> m_networkSpecies; ///< Vector of unique species in the network.
        std::unordered_map<std::string_view, fourdst::atomic::Species> m_networkSpeciesMap; ///< Map from species name to Species object.
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap; ///< Map from species to their index in the stoichiometry matrix.
        std::unordered_map<size_t, fourdst::atomic::Species> m_indexToSpeciesMap; ///< Map from index to species in the stoichiometry matrix.


        mutable CppAD::ADFun<double> m_rhsADFun; ///< CppAD function for the right-hand side of the ODE.
        mutable CppAD::ADFun<double> m_epsADFun; ///< CppAD function for the energy generation rate.
        mutable CppAD::sparse_jac_work m_jac_work; ///< Work object for sparse Jacobian calculations.

        bool m_has_been_primed = false; ///< Flag indicating if the engine has been primed.

        CppAD::sparse_rc<std::vector<size_t>> m_full_jacobian_sparsity_pattern; ///< Full sparsity pattern for the Jacobian matrix.
        std::set<std::pair<size_t, size_t>> m_full_sparsity_set; ///< For quick lookups of the base sparsity pattern

        std::vector<std::unique_ptr<AtomicReverseRate>> m_atomicReverseRates;

        screening::ScreeningType m_screeningType = screening::ScreeningType::BARE; ///< Screening type for the reaction network. Default to no screening.
        std::unique_ptr<screening::ScreeningModel> m_screeningModel = screening::selectScreeningModel(m_screeningType);

        bool m_usePrecomputation = true; ///< Flag to enable or disable using precomputed reactions for efficiency. Mathematically, this should not change the results. Generally end users should not need to change this.
        bool m_useReverseReactions = true; ///< Flag to enable or disable reverse reactions. If false, only forward reactions are considered.

        BuildDepthType m_depth;

        std::vector<PrecomputedReaction> m_precomputedReactions; ///< Precomputed reactions for efficiency.
        std::unordered_map<uint64_t, size_t> m_precomputedReactionIndexMap; ///< Set of hashed precomputed reactions for quick lookup.
        std::unique_ptr<partition::PartitionFunction> m_partitionFunction; ///< Partition function for the network.

    private:
        /**
         * @brief Synchronizes the internal maps.
         *
         * This method synchronizes the internal maps used by the engine,
         * including the species map, reaction ID map, and species-to-index map.
         * It also generates the stoichiometry matrix and records the AD tape.
         */
        void syncInternalMaps();

        /**
         * @brief Collects the unique species in the network.
         *
         * This method collects the unique species in the network from the
         * reactants and products of all reactions.
         */
        void collectNetworkSpecies();

        /**
         * @brief Populates the reaction ID map.
         *
         * This method populates the reaction ID map, which maps reaction IDs
         * to REACLIBReaction objects.
         */
        void populateReactionIDMap();

        /**
         * @brief Populates the species-to-index map.
         *
         * This method populates the species-to-index map, which maps species
         * to their index in the stoichiometry matrix.
         */
        void populateSpeciesToIndexMap();


        /**
         * @brief Records the AD tape for the right-hand side of the ODE.
         *
         * This method records the AD tape for the right-hand side of the ODE,
         * which is used to calculate the Jacobian matrix using automatic
         * differentiation.
         *
         * @throws std::runtime_error If there are no species in the network.
         */
        void recordADTape() const;

        void collectAtomicReverseRateAtomicBases();

        void precomputeNetwork();


        /**
         * @brief Validates mass and charge conservation across all reactions.
         *
         * @return True if all reactions conserve mass and charge, false otherwise.
         *
         * This method checks that all reactions in the network conserve mass
         * and charge. If any reaction does not conserve mass or charge, an
         * error message is logged and false is returned.
         */
        [[nodiscard]] bool validateConservation() const;


        [[nodiscard]] StepDerivatives<double> calculateAllDerivativesUsingPrecomputation(
            const fourdst::composition::CompositionAbstract &comp,
            const std::vector<double> &bare_rates,
            const std::vector<double> &bare_reverse_rates,
            double T9,
            double rho, const reaction::ReactionSet &activeReactions
        ) const;

        /**
         * @brief Calculates the molar reaction flow for a given reaction.
         *
         * @tparam T The numeric type to use for the calculation.
         * @param reaction The reaction for which to calculate the flow.
         * @param Y Vector of molar abundances for all species in the network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param Ye
         * @param mue
         * @param speciesIDLookup
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the net rate at which the given reaction proceeds
         * under the current state.
         */
        template <IsArithmeticOrAD T>
        T calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<T>& Y,
            T T9,
            T rho,
            T Ye,
            T mue,
            const std::function<std::optional<size_t>(const fourdst::atomic::Species &)>&speciesIDLookup
        ) const;

        template<IsArithmeticOrAD T>
        T calculateReverseMolarReactionFlow(
            T T9,
            T rho,
            std::vector<T> screeningFactors,
            const std::vector<T>& Y,
            size_t reactionIndex,
            const reaction::Reaction &reaction
        ) const;

        /**
         * @brief Calculates all derivatives (dY/dt) and the energy generation rate.
         *
         * @tparam T The numeric type to use for the calculation.
         * @param Y_in Vector of molar abundances for all species in the network.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param Ye
         * @param mue
         * @param speciesLookup
         * @param reactionLookup
         * @return StepDerivatives<T> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate for the current state.
         */
        template<IsArithmeticOrAD T>
        [[nodiscard]] StepDerivatives<T> calculateAllDerivatives(
            const std::vector<T>& Y_in,
            T T9,
            T rho,
            T Ye,
            T mue,
            std::function<std::optional<size_t>(const fourdst::atomic::Species &)> speciesLookup, const std::function<bool(const
                reaction::Reaction &)>& reactionLookup
        ) const;
    };

    template <IsArithmeticOrAD T>
    T GraphEngine::calculateReverseMolarReactionFlow(
        T T9,
        T rho,
        std::vector<T> screeningFactors,
        const std::vector<T>& Y,
        size_t reactionIndex,
        const reaction::Reaction &reaction
    ) const {
        if (!m_useReverseReactions) {
            return static_cast<T>(0.0); // If reverse reactions are not used, return zero
        }
        T reverseMolarFlow = static_cast<T>(0.0);

        if (reaction.qValue() != 0.0) {
            T reverseRateConstant = static_cast<T>(0.0);
            if constexpr (std::is_same_v<T, ADDouble>) { // Check if T is an AD type at compile time
                const auto& atomic_func_ptr = m_atomicReverseRates[reactionIndex];
                if (atomic_func_ptr != nullptr) {
                    // A. Instantiate the atomic operator for the specific reaction
                    // B. Marshal the input vector
                    std::vector<T> ax = { T9 };

                    std::vector<T> ay(1);
                    (*atomic_func_ptr)(ax, ay);
                    reverseRateConstant = static_cast<T>(ay[0]);
                } else {
                    return reverseMolarFlow; // If no atomic function is available, return zero
                }
            } else { // The case where T is of type double
                // A,B If not calling with an AD type, calculate the reverse rate directly
                std::vector<std::string> symbols;
                symbols.reserve(m_networkSpecies.size());
                for (const auto& species : m_networkSpecies) {
                    symbols.emplace_back(species.name());
                }
                std::vector<double> X;
                X.reserve(m_networkSpecies.size());
                for (const auto& species: m_networkSpecies) {
                    double Xi = species.mass() * Y[m_speciesToIndexMap.at(species)];
                    X.push_back(Xi);
                }
                fourdst::composition::Composition comp(symbols, X);
                reverseRateConstant = calculateReverseRate(reaction, T9, rho, comp);
            }

            // C. Get product multiplicities
            std::unordered_map<fourdst::atomic::Species, int> productCounts;
            for (const auto& product : reaction.products()) {
                productCounts[product]++;
            }

            // D. Calculate the symmetry factor
            T reverseSymmetryFactor = static_cast<T>(1.0);
            for (const auto &count: productCounts | std::views::values) {
                reverseSymmetryFactor /= static_cast<T>(std::tgamma(static_cast<double>(count + 1))); // Gamma function for factorial
            }

            // E. Calculate the abundance term
            T productAbundanceTerm = static_cast<T>(1.0);
            for (const auto& [species, count] : productCounts) {
                const size_t speciesIndex = m_speciesToIndexMap.at(species);
                productAbundanceTerm *= CppAD::pow(Y[speciesIndex], count);
            }

            // F. Determine the power for the density term
            const size_t num_products = reaction.products().size();
            const T rho_power = CppAD::pow(rho, static_cast<T>(num_products > 1 ? num_products - 1 : 0)); // Density raised to the power of (N-1) for N products

            // G. Assemble the reverse molar flow rate
            reverseMolarFlow = screeningFactors[reactionIndex] *
                               reverseRateConstant *
                               productAbundanceTerm *
                               reverseSymmetryFactor *
                               rho_power;
        }
        return reverseMolarFlow;

    }

    template<IsArithmeticOrAD T>
    StepDerivatives<T> GraphEngine::calculateAllDerivatives(
        const std::vector<T>& Y_in,
        const T T9,
        const T rho,
        const T Ye,
        const T mue,
        const std::function<std::optional<size_t>(const fourdst::atomic::Species &)> speciesLookup,
        const std::function<bool(const reaction::Reaction &)>& reactionLookup
    ) const {
        std::vector<T> screeningFactors = m_screeningModel->calculateScreeningFactors(
            m_reactions,
            m_networkSpecies,
            Y_in,
            T9,
            rho
        );

        // --- Setup output derivatives structure ---
        // We use a vector internally since indexed lookups are much cheeper, O(1)
        std::vector<T> dydt_vec;
        dydt_vec.resize(m_reactions.size(), static_cast<T>(0.0));

        // --- AD Pre-setup (flags to control conditionals in an AD safe / branch aware manner) ---
        // ----- Constants for AD safe calculations ---
        const T zero = static_cast<T>(0.0);
        const T one = static_cast<T>(1.0);

        // ----- Initialize variables for molar concentration product and thresholds ---
        // Note: the logic here is that we use CppAD::CondExprLt to test thresholds and if they are less we set the flag
        //       to zero so that the final returned reaction flow is 0. This is as opposed to standard if statements
        //       which create branches that break the AD tape.
        const T rho_threshold = static_cast<T>(MIN_DENSITY_THRESHOLD);

        // --- Check if the density is below the threshold where we ignore reactions ---
        T threshold_flag = CppAD::CondExpLt(rho, rho_threshold, zero, one); // If rho < threshold, set flag to 0

        std::vector<T> Y = Y_in;
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            // We use CppAD::CondExpLt to handle AD taping and prevent branching
            // Note that while this is syntactically more complex this is equivalent to
            // if (Y[i] < 0) {Y[i] = 0;}
            // The issue is that this would introduce a branch which would require the auto diff tape to be re-recorded
            // each timestep, which is very inefficient.
            Y[i] = CppAD::CondExpLt(Y[i], zero, zero, Y[i]); // Ensure no negative abundances
        }

        const T u = static_cast<T>(m_constants.u); // Atomic mass unit in grams
        const T N_A = static_cast<T>(m_constants.Na); // Avogadro's number in mol^-1
        const T c = static_cast<T>(m_constants.c); // Speed of light in cm/s

        // --- SINGLE LOOP OVER ALL REACTIONS ---
        StepDerivatives<T> result{};
        for (size_t reactionIndex = 0; reactionIndex < m_reactions.size(); ++reactionIndex) {
            bool skipReaction = false;
            const auto& reaction = m_reactions[reactionIndex];

            if (!reactionLookup(reaction)) {
                continue; // Skip this reaction if not in the "active" reaction set
            }
            for (const auto& reactant : reaction.reactant_species()) {
                if (!speciesLookup(reactant).has_value()) {
                    skipReaction = true;
                    break;
                }
            }
            if (skipReaction) {
                continue; // Skip this reaction if any reactant is not present
            }

            // 1. Calculate forward reaction rate
            const T forwardMolarReactionFlow = screeningFactors[reactionIndex] *
                calculateMolarReactionFlow<T>(
                    reaction,
                    Y,
                    T9,
                    rho,
                    Ye,
                    mue,
                    speciesLookup
                );

            // 2. Calculate reverse reaction rate
            T reverseMolarFlow = static_cast<T>(0.0);
            // Do not calculate reverse flow for weak reactions since photodisintegration does not apply
            if ((reaction.type() == reaction::ReactionType::LOGICAL_REACLIB || reaction.type() == reaction::ReactionType::REACLIB) && m_useReverseReactions) {
                reverseMolarFlow = calculateReverseMolarReactionFlow<T>(
                    T9,
                    rho,
                    screeningFactors,
                    Y,
                    reactionIndex,
                    reaction
                );
            }


            const T molarReactionFlow = forwardMolarReactionFlow - reverseMolarFlow; // Net molar reaction flow

            // 3. Use the rate to update all relevant species derivatives (dY/dt)
            for (size_t speciesIdx = 0; speciesIdx < m_networkSpecies.size(); ++speciesIdx) {
                const auto& species = m_networkSpecies[speciesIdx];
                const T nu_ij = static_cast<T>(reaction.stoichiometry(species));
                const T dydt_increment = threshold_flag * molarReactionFlow * nu_ij;
                dydt_vec[speciesIdx] += dydt_increment;
                result.reactionContributions[species][std::string(reaction.id())] = dydt_increment;
            }

        }

        T massProductionRate = static_cast<T>(0.0); // [mol][s^-1]
        for (const auto& [species, deriv] : std::views::zip(m_networkSpecies, dydt_vec)) {
            massProductionRate += deriv * species.mass() * u;
            result.dydt[species] = deriv; // [mol][s^-1][g^-1]
        }

        result.nuclearEnergyGenerationRate = -massProductionRate * N_A * c * c; // [cm^2][s^-3] = [erg][s^-1][g^-1]

        return result;
    }


    template <IsArithmeticOrAD T>
    T GraphEngine::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<T>& Y,
        const T T9,
        const T rho,
        const T Ye,
        const T mue,
        const std::function<std::optional<size_t>(const fourdst::atomic::Species &)>& speciesIDLookup
    ) const {
        // --- Pre-setup (flags to control conditionals in an AD safe / branch aware manner) ---
        // ----- Constants for AD safe calculations ---
        const T zero = static_cast<T>(0.0);

        // --- Calculate the molar reaction rate (in units of [s^-1][cm^3(N-1)][mol^(1-N)] for N reactants) ---
        const T k_reaction = reaction.calculate_rate(T9, rho, Ye, mue, Y, m_indexToSpeciesMap);

        // --- Cound the number of each reactant species to account for species multiplicity ---
        std::unordered_map<fourdst::atomic::Species, int> reactant_counts;
        reactant_counts.reserve(reaction.reactants().size());
        for (const auto& reactant : reaction.reactants()) {
            reactant_counts[reactant] = reaction.countReactantOccurrences(reactant);
        }
        const int totalReactants = static_cast<int>(reaction.reactants().size());

        // --- Accumulator for the molar concentration ---
        auto molar_concentration_product = static_cast<T>(1.0);

        // --- Loop through each unique reactant species and calculate the molar concentration for that species then multiply that into the accumulator ---
        for (const auto& [species, count] : reactant_counts) {
            // --- Resolve species to molar abundance ---
            const std::optional<size_t> species_index = speciesIDLookup(species);
            if (!species_index.has_value()) {
                return static_cast<T>(0.0); // If any reactant is not present, the reaction cannot proceed
            }
            const T Yi = Y[species_index.value()];

            // --- If count is > 1 , we need to raise the molar concentration to the power of count since there are really count bodies in that reaction ---
            molar_concentration_product *= CppAD::pow(Yi, static_cast<T>(count)); // ni^count

            // --- Apply factorial correction for identical reactions ---
            if (count > 1) {
                molar_concentration_product /= static_cast<T>(std::tgamma(static_cast<double>(count + 1))); // Gamma function for factorial
            }

        }
        // --- Final reaction flow calculation [mol][s^-1][g^-1] ---
        // Note: If the threshold flag ever gets set to zero this will return zero.
        //       This will result basically in multiple branches being written to the AD tape, which will make
        //       the tape more expensive to record, but it will also mean that we only need to record it once for
        //       the entire network.
        const T densityTerm = CppAD::pow(rho, totalReactants > 1 ? static_cast<T>(totalReactants - 1) : zero); // Density raised to the power of (N-1) for N reactants
        return molar_concentration_product * k_reaction * densityTerm;
    }

};
