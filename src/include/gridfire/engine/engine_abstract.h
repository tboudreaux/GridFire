#pragma once

#include "gridfire/reaction/reaction.h"
#include "gridfire/network.h"
#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"

#include "gridfire/engine/types/reporting.h"
#include "gridfire/engine/types/building.h"

#include "gridfire/expectations/expected_engine.h"

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

namespace gridfire {

    /**
     * @brief Concept for types allowed in engine calculations.
     *
     * This concept restricts template parameters to either double or CppAD::AD<double>,
     * enabling both standard and automatic differentiation types.
     */
    template<typename T>
    concept IsArithmeticOrAD = std::is_same_v<T, double> || std::is_same_v<T, CppAD::AD<double>>;

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
        std::vector<T> dydt; ///< Derivatives of abundances (dY/dt for each species).
        T nuclearEnergyGenerationRate = T(0.0); ///< Specific energy generation rate (e.g., erg/g/s).
    };

    using SparsityPattern = std::vector<std::pair<size_t, size_t>>;

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
         * @param Y Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<double> containing dY/dt and energy generation rate.
         *
         * This function must be implemented by derived classes to compute the
         * time derivatives of all species and the specific nuclear energy generation
         * rate for the current state.
         */
        [[nodiscard]] virtual std::expected<StepDerivatives<double>, expectations::StaleEngineError> calculateRHSAndEnergy(
            const std::vector<double>& Y,
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
     * QSE (Quasi-Steady-State Equilibrium) handling, and more.
     */
    class DynamicEngine : public Engine {
    public:
        /**
         * @brief Generate the Jacobian matrix for the current state.
         *
         * @param Y_dynamic Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * This method must compute and store the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state. The matrix can then be accessed via getJacobianMatrixEntry().
         */
        virtual void generateJacobianMatrix(
            const std::vector<double>& Y_dynamic,
            double T9,
            double rho
        ) const = 0;

        virtual void generateJacobianMatrix(
            const std::vector<double>& Y_dynamic,
            double T9,
            double rho,
            const SparsityPattern& sparsityPattern
        ) const {
            throw std::logic_error("Sparsity pattern not supported by this engine.");
        }

        /**
         * @brief Get an entry from the previously generated Jacobian matrix.
         *
         * @param i Row index (species index).
         * @param j Column index (species index).
         * @return Value of the Jacobian matrix at (i, j).
         *
         * The Jacobian must have been generated by generateJacobianMatrix() before calling this.
         */
        [[nodiscard]] virtual double getJacobianMatrixEntry(
            int i,
            int j
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
         * @param speciesIndex Index of the species.
         * @param reactionIndex Index of the reaction.
         * @return Stoichiometric coefficient for the species in the reaction.
         *
         * The stoichiometry matrix must have been generated by generateStoichiometryMatrix().
         */
        [[nodiscard]] virtual int getStoichiometryMatrixEntry(
            int speciesIndex,
            int reactionIndex
        ) const = 0;

        /**
         * @brief Calculate the molar reaction flow for a given reaction.
         *
         * @param reaction The reaction for which to calculate the flow.
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the net rate at which the given reaction proceeds
         * under the current state.
         */
        [[nodiscard]] virtual double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const = 0;

        /**
         * @brief Get the set of logical reactions in the network.
         *
         * @return Reference to the LogicalReactionSet containing all reactions.
         */
        [[nodiscard]] virtual const reaction::LogicalReactionSet& getNetworkReactions() const = 0;

        virtual void setNetworkReactions(const reaction::LogicalReactionSet& reactions) = 0;

        /**
         * @brief Compute timescales for all species in the network.
         *
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method estimates the timescale for abundance change of each species,
         * which can be used for timestep control, diagnostics, and reaction network culling.
         */
        [[nodiscard]] virtual std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const = 0;

        [[nodiscard]] virtual std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const std::vector<double>& Y,
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

        [[nodiscard]] virtual int getSpeciesIndex(const fourdst::atomic::Species &species) const = 0;

        [[nodiscard]] virtual std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const = 0;

        [[nodiscard]] virtual PrimingReport primeEngine(const NetIn &netIn) = 0;

        [[nodiscard]] virtual BuildDepthType getDepth() const {
            throw std::logic_error("Network depth not supported by this engine.");
        }

        virtual void rebuild(const fourdst::composition::Composition& comp, BuildDepthType depth) {
            throw std::logic_error("Setting network depth not supported by this engine.");
        }

    };
}