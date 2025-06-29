#pragma once

#include "gridfire/reaction/reaction.h"

#include <vector>
#include <unordered_map>

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
        virtual const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const = 0;

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
        virtual StepDerivatives<double> calculateRHSAndEnergy(
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
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * This method must compute and store the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state. The matrix can then be accessed via getJacobianMatrixEntry().
         */
        virtual void generateJacobianMatrix(
            const std::vector<double>& Y,
            double T9, double rho
        ) = 0;

        /**
         * @brief Get an entry from the previously generated Jacobian matrix.
         *
         * @param i Row index (species index).
         * @param j Column index (species index).
         * @return Value of the Jacobian matrix at (i, j).
         *
         * The Jacobian must have been generated by generateJacobianMatrix() before calling this.
         */
        virtual double getJacobianMatrixEntry(
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
        virtual int getStoichiometryMatrixEntry(
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
        virtual double calculateMolarReactionFlow(
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
        virtual const reaction::LogicalReactionSet& getNetworkReactions() const = 0;

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
        virtual std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const = 0;
    };
}