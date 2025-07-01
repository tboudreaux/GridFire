#pragma once

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/composition.h"
#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "gridfire/network.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "cppad/cppad.hpp"

// PERF: The function getNetReactionStoichiometry returns a map of species to their stoichiometric coefficients for a given reaction.
//       this makes extra copies of the species, which is not ideal and could be optimized further.
//       Even more relevant is the member m_reactionIDMap which makes copies of a REACLIBReaction for each reaction ID.
//       REACLIBReactions are quite large data structures, so this could be a performance bottleneck.

namespace gridfire {
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
        explicit GraphEngine(const fourdst::composition::Composition &composition);

        /**
         * @brief Constructs a GraphEngine from a set of reactions.
         *
         * @param reactions The set of reactions to use in the network.
         *
         * This constructor uses the given set of reactions to construct the
         * reaction network.
         */
        explicit GraphEngine(reaction::LogicalReactionSet reactions);

        /**
         * @brief Calculates the right-hand side (dY/dt) and energy generation rate.
         *
         * @param Y Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<double> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate for the current state.
         *
         * @see StepDerivatives
         */
        [[nodiscard]] StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state.
         *
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         *
         * This method computes and stores the Jacobian matrix (∂(dY/dt)_i/∂Y_j)
         * for the current state using automatic differentiation. The matrix can
         * then be accessed via `getJacobianMatrixEntry()`.
         *
         * @see getJacobianMatrixEntry()
         */
        void generateJacobianMatrix(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) override;

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
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the net rate at which the given reaction proceeds
         * under the current state.
         */
        [[nodiscard]] double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const std::vector<double>&Y,
            const double T9,
            const double rho
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
        [[nodiscard]] const reaction::LogicalReactionSet& getNetworkReactions() const override;

        /**
         * @brief Gets an entry from the previously generated Jacobian matrix.
         *
         * @param i Row index (species index).
         * @param j Column index (species index).
         * @return Value of the Jacobian matrix at (i, j).
         *
         * The Jacobian must have been generated by `generateJacobianMatrix()` before calling this.
         *
         * @see generateJacobianMatrix()
         */
        [[nodiscard]] double getJacobianMatrixEntry(
            const int i,
            const int j
        ) const override;

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
         * @param speciesIndex Index of the species.
         * @param reactionIndex Index of the reaction.
         * @return Stoichiometric coefficient for the species in the reaction.
         *
         * The stoichiometry matrix must have been generated by `generateStoichiometryMatrix()`.
         *
         * @see generateStoichiometryMatrix()
         */
        [[nodiscard]] int getStoichiometryMatrixEntry(
            const int speciesIndex,
            const int reactionIndex
        ) const override;

        /**
         * @brief Computes timescales for all species in the network.
         *
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Map from Species to their characteristic timescales (s).
         *
         * This method estimates the timescale for abundance change of each species,
         * which can be used for timestep control or diagnostics.
         */
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        void update(const NetIn& netIn) override;

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

        void setScreeningModel(screening::ScreeningType) override;

        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;


    private:
        reaction::LogicalReactionSet m_reactions; ///< Set of REACLIB reactions in the network.
        std::unordered_map<std::string_view, reaction::Reaction*> m_reactionIDMap; ///< Map from reaction ID to REACLIBReaction. //PERF: This makes copies of REACLIBReaction and could be a performance bottleneck.

        std::vector<fourdst::atomic::Species> m_networkSpecies; ///< Vector of unique species in the network.
        std::unordered_map<std::string_view, fourdst::atomic::Species> m_networkSpeciesMap; ///< Map from species name to Species object.
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap; ///< Map from species to their index in the stoichiometry matrix.

        boost::numeric::ublas::compressed_matrix<int> m_stoichiometryMatrix; ///< Stoichiometry matrix (species x reactions).
        boost::numeric::ublas::compressed_matrix<double> m_jacobianMatrix; ///< Jacobian matrix (species x species).

        CppAD::ADFun<double> m_rhsADFun; ///< CppAD function for the right-hand side of the ODE.

        screening::ScreeningType m_screeningType = screening::ScreeningType::BARE; ///< Screening type for the reaction network. Default to no screening.
        std::unique_ptr<screening::ScreeningModel> m_screeningModel = screening::selectScreeningModel(m_screeningType);

        Config& m_config = Config::getInstance();
        Constants& m_constants = Constants::getInstance(); ///< Access to physical constants.
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

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
         * @brief Reserves space for the Jacobian matrix.
         *
         * This method reserves space for the Jacobian matrix, which is used
         * to store the partial derivatives of the right-hand side of the ODE
         * with respect to the species abundances.
         */
        void reserveJacobianMatrix();

        /**
         * @brief Records the AD tape for the right-hand side of the ODE.
         *
         * This method records the AD tape for the right-hand side of the ODE,
         * which is used to calculate the Jacobian matrix using automatic
         * differentiation.
         *
         * @throws std::runtime_error If there are no species in the network.
         */
        void recordADTape();

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

        /**
         * @brief Validates the composition against the current reaction set.
         *
         * @param composition The composition to validate.
         * @param culling The culling threshold to use.
         * @param T9 The temperature to use.
         *
         * This method validates the composition against the current reaction set.
         * If the composition is not compatible with the reaction set, the
         * reaction set is rebuilt from the composition.
         */
        void validateComposition(
            const fourdst::composition::Composition &composition,
            double culling,
            double T9
        );

        /**
         * @brief Calculates the molar reaction flow for a given reaction.
         *
         * @tparam T The numeric type to use for the calculation.
         * @param reaction The reaction for which to calculate the flow.
         * @param Y Vector of current abundances.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Molar flow rate for the reaction (e.g., mol/g/s).
         *
         * This method computes the net rate at which the given reaction proceeds
         * under the current state.
         */
        template <IsArithmeticOrAD T>
        T calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<T> &Y,
            const T T9,
            const T rho
        ) const;

        /**
         * @brief Calculates all derivatives (dY/dt) and the energy generation rate.
         *
         * @tparam T The numeric type to use for the calculation.
         * @param Y_in Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<T> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate for the current state.
         */
        template<IsArithmeticOrAD T>
        [[nodiscard]] StepDerivatives<T> calculateAllDerivatives(
            const std::vector<T> &Y_in,
            T T9,
            T rho
        ) const;

        /**
         * @brief Calculates all derivatives (dY/dt) and the energy generation rate (double precision).
         *
         * @param Y_in Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<double> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate for the current state using
         * double precision arithmetic.
         */
        [[nodiscard]] StepDerivatives<double> calculateAllDerivatives(
            const std::vector<double>& Y_in,
            const double T9,
            const double rho
        ) const;

        /**
         * @brief Calculates all derivatives (dY/dt) and the energy generation rate (automatic differentiation).
         *
         * @param Y_in Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<ADDouble> containing dY/dt and energy generation rate.
         *
         * This method calculates the time derivatives of all species and the
         * specific nuclear energy generation rate for the current state using
         * automatic differentiation.
         */
        [[nodiscard]] StepDerivatives<ADDouble> calculateAllDerivatives(
            const std::vector<ADDouble>& Y_in,
            const ADDouble &T9,
            const ADDouble &rho
        ) const;
    };


    template<IsArithmeticOrAD T>
    StepDerivatives<T> GraphEngine::calculateAllDerivatives(
        const std::vector<T> &Y_in, T T9, T rho) const {
        std::vector<T> screeningFactors = m_screeningModel->calculateScreeningFactors(
            m_reactions,
            m_networkSpecies,
            Y_in,
            T9,
            rho
        );

        // --- Setup output derivatives structure ---
        StepDerivatives<T> result;
        result.dydt.resize(m_networkSpecies.size(), static_cast<T>(0.0));

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

        const T u = static_cast<T>(m_constants.get("u").value); // Atomic mass unit in grams
        const T N_A = static_cast<T>(m_constants.get("N_a").value); // Avogadro's number in mol^-1
        const T c = static_cast<T>(m_constants.get("c").value); // Speed of light in cm/s

        // --- SINGLE LOOP OVER ALL REACTIONS ---
        for (size_t reactionIndex = 0; reactionIndex < m_reactions.size(); ++reactionIndex) {
            const auto& reaction = m_reactions[reactionIndex];

            // 1. Calculate reaction rate
            const T molarReactionFlow = screeningFactors[reactionIndex] * calculateMolarReactionFlow<T>(reaction, Y, T9, rho);

            // 2. Use the rate to update all relevant species derivatives (dY/dt)
            for (size_t speciesIndex = 0; speciesIndex < m_networkSpecies.size(); ++speciesIndex) {
                const T nu_ij = static_cast<T>(m_stoichiometryMatrix(speciesIndex, reactionIndex));
                result.dydt[speciesIndex] += threshold_flag * nu_ij * molarReactionFlow / rho;
            }
        }

        T massProductionRate = static_cast<T>(0.0); // [mol][s^-1]
        for (const auto& [species, index] : m_speciesToIndexMap) {
            massProductionRate += result.dydt[index] * species.mass() * u;
        }

        result.nuclearEnergyGenerationRate = -massProductionRate * N_A * c * c; // [cm^2][s^-3] = [erg][s^-1][g^-1]

        return result;
    }


    template <IsArithmeticOrAD T>
    T GraphEngine::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<T> &Y,
        const T T9,
        const T rho
    ) const {

        // --- Pre-setup (flags to control conditionals in an AD safe / branch aware manner) ---
        // ----- Constants for AD safe calculations ---
        const T zero = static_cast<T>(0.0);
        const T one = static_cast<T>(1.0);

        // ----- Initialize variables for molar concentration product and thresholds ---
        // Note: the logic here is that we use CppAD::CondExprLt to test thresholds and if they are less we set the flag
        //       to zero so that the final returned reaction flow is 0. This is as opposed to standard if statements
        //       which create branches that break the AD tape.
        const T Y_threshold = static_cast<T>(MIN_ABUNDANCE_THRESHOLD);
        T threshold_flag = one;

        // --- Calculate the molar reaction rate (in units of [s^-1][cm^3(N-1)][mol^(1-N)] for N reactants) ---
        const T k_reaction = reaction.calculate_rate(T9);

        // --- Cound the number of each reactant species to account for species multiplicity ---
        std::unordered_map<std::string, int> reactant_counts;
        reactant_counts.reserve(reaction.reactants().size());
        for (const auto& reactant : reaction.reactants()) {
            reactant_counts[std::string(reactant.name())]++;
        }

        // --- Accumulator for the molar concentration ---
        auto molar_concentration_product = static_cast<T>(1.0);

        // --- Loop through each unique reactant species and calculate the molar concentration for that species then multiply that into the accumulator ---
        for (const auto& [species_name, count] : reactant_counts) {
            // --- Resolve species to molar abundance ---
            // PERF: Could probably optimize out this lookup
            const auto species_it = m_speciesToIndexMap.find(m_networkSpeciesMap.at(species_name));
            const size_t species_index = species_it->second;
            const T Yi = Y[species_index];

            // --- Check if the species abundance is below the threshold where we ignore reactions ---
            threshold_flag *= CppAD::CondExpLt(Yi, Y_threshold, zero, one);

            // --- Convert from molar abundance to molar concentration ---
            T molar_concentration = Yi * rho;

            // --- If count is > 1 , we need to raise the molar concentration to the power of count since there are really count bodies in that reaction ---
            molar_concentration_product *= CppAD::pow(molar_concentration, static_cast<T>(count)); // ni^count

            // --- Apply factorial correction for identical reactions ---
            if (count > 1) {
                molar_concentration_product /= static_cast<T>(std::tgamma(static_cast<double>(count + 1))); // Gamma function for factorial
            }
        }
        // --- Final reaction flow calculation [mol][s^-1][cm^-3] ---
        // Note: If the threshold flag ever gets set to zero this will return zero.
        //       This will result basically in multiple branches being written to the AD tape, which will make
        //       the tape more expensive to record, but it will also mean that we only need to record it once for
        //       the entire network.
        return molar_concentration_product * k_reaction * threshold_flag;
    }
};