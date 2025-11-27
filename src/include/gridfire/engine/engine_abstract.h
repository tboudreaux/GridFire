#pragma once

#include "gridfire/reaction/reaction.h"
#include "gridfire/types/types.h"
#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"

#include "gridfire/engine/types/reporting.h"
#include "gridfire/engine/types/building.h"
#include "gridfire/engine/types/jacobian.h"

#include "fourdst/composition/composition_abstract.h"

#include <vector>
#include <unordered_map>
#include <utility>
#include <expected>

/**
 * @file engine_abstract.h
 * @brief Abstract interfaces for reaction network engines in GridFire.
 *
 * This header defines the abstract base classes and concepts for implementing
 * reaction network solvers in the GridFire framework. It provides the contract
 * for calculating right-hand sides, energy generation, Jacobians, stoichiometry,
 * and other core operations required for time integration of nuclear reaction networks.
 *
 * @author
 * Emily M. Boudreaux
 */

namespace gridfire::engine {


    /**
     * @brief Structure holding derivatives and energy generation for a network step.
     *
     * @tparam T Numeric type (double or CppAD::AD<double>).
     *
     * This struct is used to return both the time derivatives of all species abundances
     * and the specific nuclear energy generation rate for a single network evaluation.
     *
     * Example usage:
     * @code
     * StepDerivatives<double> result = engine.calculateRHSAndEnergy(Y, T9, rho);
     * for (double dydt_i : result.dydt) {
     *     // Use derivative
     * }
     * double energyRate = result.nuclearEnergyGenerationRate;
     * @endcode
     */
    template <IsArithmeticOrAD T>
    struct StepDerivatives {
        std::map<fourdst::atomic::Species, T> dydt{}; ///< Derivatives of abundances (dY/dt for each species).
        T nuclearEnergyGenerationRate = T(0.0); ///< Specific energy generation rate (e.g., erg/g/s).
        std::map<fourdst::atomic::Species, std::unordered_map<std::string, T>> reactionContributions{};
        T neutrinoEnergyLossRate = T(0.0); // (erg/g/s)
        T totalNeutrinoFlux = T(0.0);      // (neutrinos/g/s)


        StepDerivatives() : dydt(), nuclearEnergyGenerationRate(T(0.0)) {}
    };

    using SparsityPattern = std::vector<std::pair<size_t, size_t>>; ///< Type alias for sparsity pattern representation.

    /**
     * @brief Structure holding derivatives of energy generation rate with respect to T and rho.
     *
     * This struct encapsulates the partial derivatives of the specific nuclear energy
     * generation rate with respect to temperature and density.
     */
    struct EnergyDerivatives {
        double dEps_dT = 0.0;  ///< Partial derivative of energy generation rate with respect to temperature.
        double dEps_dRho = 0.0;///< Partial derivative of energy generation rate with respect to density.

        friend std::ostream& operator<<(std::ostream& os, const EnergyDerivatives& ed) {
            os << "<dε/dT: " << ed.dEps_dT << ", dε/dρ: " << ed.dEps_dRho << ">";
            return os;
        }
    };

    /**
     * @brief Enumeration of possible engine statuses.
     *
     * This enum defines the various states an engine can be in after performing
     * calculations, such as being up-to-date (OKAY), needing an update (STALE),
     * or encountering an error (ERROR).
     */
    enum class EngineStatus {
        OKAY,
        STALE,
        ERROR,
        COUNT
    };

    /**
     * @brief Convert EngineStatus enum to string representation.
     *
     * @param status The EngineStatus value to convert.
     * @return A string_view representing the name of the EngineStatus.
     */
    constexpr std::string_view EngineStatus_to_string(const EngineStatus status) {
        constexpr std::array<std::string_view, static_cast<size_t>(EngineStatus::COUNT)> names = {
            "OKAY",
            "STALE",
            "ERROR"
        };
        return names[static_cast<size_t>(status)];
    }

    /**
     * @brief Abstract base class for a reaction network engine.
     *
     * This class defines the minimal interface for a reaction network engine,
     * which is responsible for evaluating the right-hand side (dY/dt) and
     * energy generation for a given set of abundances, temperature, and density.
     *
     * Intended usage: Derive from this class to implement a concrete engine
     * for a specific network or integration method.
     *
     * Example:
     * @code
     * class MyEngine : public gridfire::Engine {
     *     // Implement required methods...
     * };
     * @endcode
     */
    class Engine {
    public:
        /**
         * @brief Virtual destructor.
         */
        virtual ~Engine() = default;

        /**
         * @brief Get the list of species in the network.
         * @return Vector of Species objects representing all network species.
         */
        [[nodiscard]] virtual const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const = 0;

        /**
         * @brief Calculate the right-hand side (dY/dt) and energy generation.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return expected<StepDerivatives<double>> containing either dY/dt and energy generation rate or a stale engine
         * error indicating that the engine must be updated
         *
         * This function must be implemented by derived classes to compute the
         * time derivatives of all species and the specific nuclear energy generation
         * rate for the current state.
         */
        [[nodiscard]] virtual std::expected<StepDerivatives<double>, EngineStatus> calculateRHSAndEnergy(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;
    };

    /**
     * @brief Abstract class for engines supporting Jacobian and stoichiometry operations.
     *
     * Extends Engine with additional methods for:
     *  - Generating and accessing the Jacobian matrix (for implicit solvers).
     *  - Generating and accessing the stoichiometry matrix.
     *  - Calculating molar reaction flows for individual reactions.
     *  - Accessing the set of logical reactions in the network.
     *  - Computing timescales for each species.
     *
     * Intended usage: Derive from this class to implement engines that support
     * advanced solver features such as implicit integration, sensitivity analysis,
     * QSE (Quasi-Steady-State Equilibrium) handling, and more. Generally this will be the main engine type
     */
    class DynamicEngine : public Engine {
    public:
        /**
         * @brief Generate the Jacobian matrix for the current state.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * This method must compute and store the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state. The matrix can then be accessed via getJacobianMatrixEntry().
         */
        [[nodiscard]] virtual NetworkJacobian generateJacobianMatrix(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Generate the Jacobian matrix for the current state using a subset of active species.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param activeSpecies The set of species to include in the Jacobian calculation.
         *
         * This method must compute and store the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state, considering only the specified subset of active species.
         * The matrix can then be accessed via getJacobianMatrixEntry().
         */
        [[nodiscard]] virtual NetworkJacobian generateJacobianMatrix(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const std::vector<fourdst::atomic::Species>& activeSpecies
        ) const = 0;


        /**
         * @brief Generate the Jacobian matrix for the current state with a specified sparsity pattern.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param sparsityPattern The sparsity pattern to use for the Jacobian matrix.
         *
         * This method must compute and store the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state using automatic differentiation, taking into
         * account the provided sparsity pattern. The matrix can then be accessed
         * via `getJacobianMatrixEntry()`.
         *
         * @see getJacobianMatrixEntry()
         */
        [[nodiscard]] virtual NetworkJacobian generateJacobianMatrix(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho,
            const SparsityPattern& sparsityPattern
        ) const = 0;


        /**
         * @brief Generate the stoichiometry matrix for the network.
         *
         * This method must compute and store the stoichiometry matrix,
         * which encodes the net change of each species in each reaction.
         */
        virtual void generateStoichiometryMatrix() = 0;

        /**
         * @brief Get an entry from the stoichiometry matrix.
         *
         * @param species species to look up stoichiometry for.
         * @param reaction reaction to find
         * @return Stoichiometric coefficient for the species in the reaction.
         *
         * The stoichiometry matrix must have been generated by generateStoichiometryMatrix().
         */
        [[nodiscard]] virtual int getStoichiometryMatrixEntry(
            const fourdst::atomic::Species& species,
            const reaction::Reaction& reaction
        ) const = 0;

        /**
         * @brief Calculate the molar reaction flow for a given reaction.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the net rate at which the given reaction proceeds
         * under the current state.
         */
        [[nodiscard]] virtual double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Calculate the derivatives of the energy generation rate with respect to T and rho.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return EnergyDerivatives containing dEps/dT and dEps/dRho.
         *
         * This method computes the partial derivatives of the specific nuclear energy
         * generation rate with respect to temperature and density for the current state.
         */
        [[nodiscard]] virtual EnergyDerivatives calculateEpsDerivatives(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Get the set of logical reactions in the network.
         *
         * @return Reference to the LogicalReactionSet containing all reactions.
         */
        [[nodiscard]] virtual const reaction::ReactionSet& getNetworkReactions() const = 0;

        /**
         * @brief Set the reactions for the network.
         *
         * @param reactions The set of reactions to use in the network.
         *
         * This method replaces the current set of reactions in the network
         * with the provided set. It marks the engine as stale, requiring
         * regeneration of matrices and recalculation of rates.
         */
        virtual void setNetworkReactions(const reaction::ReactionSet& reactions) = 0;

        /**
         * @brief Compute timescales for all species in the network.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method estimates the timescale for abundance change of each species,
         * which can be used for timestep control, diagnostics, and reaction network culling.
         */
        [[nodiscard]] virtual std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Compute destruction timescales for all species in the network.
         *
         * @param comp Composition object containing current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their destruction timescales (s).
         *
         * This method estimates the destruction timescale for each species,
         * which can be useful for understanding reaction flows and equilibrium states.
         */
        [[nodiscard]] virtual std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> getSpeciesDestructionTimescales(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Update the internal state of the engine.
         *
         * @param netIn A struct containing the current network input, such as
         *              temperature, density, and composition.
         *
         * This method is intended to be implemented by derived classes to update
         * their internal state based on the provided network conditions. For example,
         * an adaptive engine might use this to re-evaluate which reactions and species
         * are active. For other engines that do not support manually updating, this
         * method might do nothing.
         *
         * @par Usage Example:
         * @code
         * NetIn input = { ... };
         * myEngine.update(input);
         * @endcode
         *
         * @post The internal state of the engine is updated to reflect the new conditions.
         */
        virtual fourdst::composition::Composition update(const NetIn &netIn) = 0;

        /**
         * @brief Check if the engine's internal state is stale.
         *
         * @param netIn A struct containing the current network input, such as
         *              temperature, density, and composition.
         * @return True if the engine's state is stale and needs to be updated; false otherwise.
         *
         * This method allows derived classes to determine if their internal state
         * is out-of-date with respect to the provided network conditions. If the engine
         * is stale, it may require a call to `update()` before performing calculations.
         *
         * @par Usage Example:
         * @code
         * NetIn input = { ... };
         * if (myEngine.isStale(input)) {
         *     // Update the engine before proceeding
         * }
         * @endcode
         */
        [[nodiscard]]
        virtual bool isStale(const NetIn& netIn) = 0;

        /**
         * @brief Set the electron screening model.
         *
         * @param model The type of screening model to use for reaction rate calculations.
         *
         * This method allows changing the screening model at runtime. Screening corrections
         * account for the electrostatic shielding of nuclei by electrons, which affects
         * reaction rates in dense stellar plasmas.
         *
         * @par Usage Example:
         * @code
         * myEngine.setScreeningModel(screening::ScreeningType::WEAK);
         * @endcode
         *
         * @post The engine will use the specified screening model for subsequent rate calculations.
         */
        virtual void setScreeningModel(screening::ScreeningType model) = 0;

        /**
         * @brief Get the current electron screening model.
         *
         * @return The currently active screening model type.
         *
         * @par Usage Example:
         * @code
         * screening::ScreeningType currentModel = myEngine.getScreeningModel();
         * @endcode
         */
        [[nodiscard]] virtual screening::ScreeningType getScreeningModel() const = 0;

        /**
         * @brief Get the index of a species in the network.
         *
         * @param species The species to look up.
         *
         * This method allows querying the index of a specific species in the
         * engine's internal representation. It is useful for accessing species
         * data efficiently.
         */
        [[nodiscard]] virtual size_t getSpeciesIndex(const fourdst::atomic::Species &species) const = 0;

        /**
         * @brief Map a NetIn object to a vector of molar abundances.
         *
         * @param netIn The input conditions for the network.
         * @return A vector of molar abundances corresponding to the species in the network.
         *
         * This method converts the input conditions into a vector of molar abundances,
         * which can be used for further calculations or diagnostics.
         */
        [[nodiscard]] virtual std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const = 0;

        /**
         * @brief Prime the engine with initial conditions.
         *
         * @param netIn The input conditions for the network.
         * @return PrimingReport containing information about the priming process.
         *
         * This method is used to prepare the engine for calculations by setting up
         * initial conditions, reactions, and species. It may involve compiling reaction
         * rates, initializing internal data structures, and performing any necessary
         * pre-computation.
         */
        [[nodiscard]] virtual PrimingReport primeEngine(const NetIn &netIn) = 0;

        /**
         * @brief Get the depth of the network.
         *
         * @return The depth of the network, which may indicate the level of detail or
         *         complexity in the reaction network.
         *
         * This method is intended to provide information about the network's structure,
         * such as how many layers of reactions or species are present. It can be useful
         * for diagnostics and understanding the network's complexity.
         */
        [[nodiscard]] virtual BuildDepthType getDepth() const {
            throw std::logic_error("Network depth not supported by this engine.");
        }

        /**
         * @brief Rebuild the network with a specified depth.
         *
         * @param comp The composition to rebuild the network with.
         * @param depth The desired depth of the network.
         *
         * This method is intended to allow dynamic adjustment of the network's depth,
         * which may involve adding or removing species and reactions based on the
         * specified depth. However, not all engines support this operation.
         */
        virtual void rebuild(const fourdst::composition::CompositionAbstract &comp, BuildDepthType depth) {
            throw std::logic_error("Setting network depth not supported by this engine.");
            // ReSharper disable once CppDFAUnreachableCode
        }

        /**
         * @brief Recursively collect composition from current engine and any sub engines if they exist.
         * @details If species i is defined in comp and in any sub engine or self composition then the molar abundance of
         * species i in the returned composition will be that defined in comp. If there are species defined in sub engine
         * compositions which are not defined in comp then their molar abundances will be based on the reported values
         * from each sub engine.
         * @note It is up to each engine to decide how to handle filling in the return composition.
         * @note These methods return an unfinalized composition which must then be finalized by the caller
         * @param comp Input composition to "normalize".
         * @param T9
         * @param rho
         * @return An updated composition which is a superset of comp. This may contain species which were culled, for
         * example, by either QSE partitioning or reaction flow rate culling
         */
        virtual fourdst::composition::Composition collectComposition(
            const fourdst::composition::CompositionAbstract &comp,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Get the status of a species in the network.
         *
         * @param species The species to check.
         * @return SpeciesStatus indicating whether the species is active, inactive, or culled.
         *
         * This method allows querying the current status of a specific species
         * within the engine's network.
         */
        [[nodiscard]] virtual SpeciesStatus getSpeciesStatus(const fourdst::atomic::Species& species) const = 0;

    };
}