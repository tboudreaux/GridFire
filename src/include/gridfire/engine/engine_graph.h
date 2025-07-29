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
#include "gridfire/partition/partition_abstract.h"
#include "gridfire/engine/procedures/construction.h"

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <ranges>

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "cppad/cppad.hpp"
#include "cppad/utility/sparse_rc.hpp"
#include "cppad/speed/sparse_jac_fun.hpp"

#include "procedures/priming.h"


#include "quill/LogMacros.h"

// PERF: The function getNetReactionStoichiometry returns a map of species to their stoichiometric coefficients for a given reaction.
//       this makes extra copies of the species, which is not ideal and could be optimized further.
//       Even more relevant is the member m_reactionIDMap which makes copies of a REACLIBReaction for each reaction ID.
//       REACLIBReactions are quite large data structures, so this could be a performance bottleneck.
// static bool isF17 = false;
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
        explicit GraphEngine(
            const fourdst::composition::Composition &composition,
            const BuildDepthType = NetworkBuildDepth::Full
        );

        explicit GraphEngine(
            const fourdst::composition::Composition &composition,
            const partition::PartitionFunction& partitionFunction,
            const BuildDepthType buildDepth = NetworkBuildDepth::Full
        );

        /**
         * @brief Constructs a GraphEngine from a set of reactions.
         *
         * @param reactions The set of reactions to use in the network.
         *
         * This constructor uses the given set of reactions to construct the
         * reaction network.
         */
        explicit GraphEngine(const reaction::LogicalReactionSet &reactions);

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
        [[nodiscard]] std::expected<StepDerivatives<double>, expectations::StaleEngineError> calculateRHSAndEnergy(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const override;

        /**
         * @brief Generates the Jacobian matrix for the current state.
         *
         * @param Y_dynamic Vector of current abundances.
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
            const std::vector<double>& Y_dynamic,
            const double T9,
            const double rho
        ) const override;

        void generateJacobianMatrix(
            const std::vector<double> &Y_dynamic,
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

        void setNetworkReactions(const reaction::LogicalReactionSet& reactions) override;

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
        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        fourdst::composition::Composition update(const NetIn &netIn) override;

        bool isStale(const NetIn &netIn) override;

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

        void setPrecomputation(bool precompute);

        [[nodiscard]] bool isPrecomputationEnabled() const;

        [[nodiscard]] const partition::PartitionFunction& getPartitionFunction() const;

        [[nodiscard]] double calculateReverseRate(
            const reaction::Reaction &reaction,
            double T9
        ) const;

        [[nodiscard]] double calculateReverseRateTwoBody(
            const reaction::Reaction &reaction,
            const double T9,
            const double forwardRate,
            const double expFactor
        ) const;

        [[nodiscard]] double calculateReverseRateTwoBodyDerivative(
            const reaction::Reaction &reaction,
            const double T9,
            const double reverseRate
        ) const;

        [[nodiscard]] bool isUsingReverseReactions() const;

        void setUseReverseReactions(bool useReverse);

        [[nodiscard]] int getSpeciesIndex(
            const fourdst::atomic::Species& species
        ) const override;

        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

        [[nodiscard]] PrimingReport primeEngine(const NetIn &netIn) override;

        [[nodiscard]] BuildDepthType getDepth() const override;

        void rebuild(const fourdst::composition::Composition& comp, const BuildDepthType depth) override;


    private:
        struct PrecomputedReaction {
            // Forward cacheing
            size_t reaction_index;
            std::vector<size_t> unique_reactant_indices;
            std::vector<int> reactant_powers;
            double symmetry_factor;
            std::vector<size_t> affected_species_indices;
            std::vector<int> stoichiometric_coefficients;

            // Reverse cacheing
            std::vector<size_t> unique_product_indices; ///< Unique product indices for reverse reactions.
            std::vector<int> product_powers; ///< Powers of each unique product in the reverse reaction.
            double reverse_symmetry_factor; ///< Symmetry factor for reverse reactions.
        };

        struct constants {
            const double u = Constants::getInstance().get("u").value; ///< Atomic mass unit in g.
            const double Na = Constants::getInstance().get("N_a").value; ///< Avogadro's number.
            const double c = Constants::getInstance().get("c").value; ///< Speed of light in cm/s.
            const double kB = Constants::getInstance().get("kB").value; ///< Boltzmann constant in erg/K.
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

        private:
            const reaction::Reaction& m_reaction;
            const GraphEngine& m_engine;
        };
    private:
        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

        constants m_constants;

        reaction::LogicalReactionSet m_reactions; ///< Set of REACLIB reactions in the network.
        std::unordered_map<std::string_view, reaction::Reaction*> m_reactionIDMap; ///< Map from reaction ID to REACLIBReaction. //PERF: This makes copies of REACLIBReaction and could be a performance bottleneck.

        std::vector<fourdst::atomic::Species> m_networkSpecies; ///< Vector of unique species in the network.
        std::unordered_map<std::string_view, fourdst::atomic::Species> m_networkSpeciesMap; ///< Map from species name to Species object.
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap; ///< Map from species to their index in the stoichiometry matrix.

        boost::numeric::ublas::compressed_matrix<int> m_stoichiometryMatrix; ///< Stoichiometry matrix (species x reactions).

        mutable boost::numeric::ublas::compressed_matrix<double> m_jacobianMatrix; ///< Jacobian matrix (species x species).
        mutable CppAD::ADFun<double> m_rhsADFun; ///< CppAD function for the right-hand side of the ODE.
        mutable CppAD::sparse_jac_work m_jac_work; ///< Work object for sparse Jacobian calculations.
        CppAD::sparse_rc<std::vector<size_t>> m_full_jacobian_sparsity_pattern; ///< Full sparsity pattern for the Jacobian matrix.

        std::vector<std::unique_ptr<AtomicReverseRate>> m_atomicReverseRates;

        screening::ScreeningType m_screeningType = screening::ScreeningType::BARE; ///< Screening type for the reaction network. Default to no screening.
        std::unique_ptr<screening::ScreeningModel> m_screeningModel = screening::selectScreeningModel(m_screeningType);

        bool m_usePrecomputation = true; ///< Flag to enable or disable using precomputed reactions for efficiency. Mathematically, this should not change the results. Generally end users should not need to change this.

        bool m_useReverseReactions = true; ///< Flag to enable or disable reverse reactions. If false, only forward reactions are considered.

        BuildDepthType m_depth;

        std::vector<PrecomputedReaction> m_precomputedReactions; ///< Precomputed reactions for efficiency.
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
         * @brief Reserves space for the Jacobian matrix.
         *
         * This method reserves space for the Jacobian matrix, which is used
         * to store the partial derivatives of the right-hand side of the ODE
         * with respect to the species abundances.
         */
        void reserveJacobianMatrix() const;

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
            const std::vector<double> &Y_in,
            const std::vector<double>& bare_rates,
            const std::vector<double> &bare_reverse_rates,
            double T9, double rho
        ) const;

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

        template<IsArithmeticOrAD T>
        T calculateReverseMolarReactionFlow(
            T T9,
            T rho,
            std::vector<T> screeningFactors,
            std::vector<T> Y,
            size_t reactionIndex,
            const reaction::LogicalReaction &reaction
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



    template <IsArithmeticOrAD T>
    T GraphEngine::calculateReverseMolarReactionFlow(
        T T9,
        T rho,
        std::vector<T> screeningFactors,
        std::vector<T> Y,
        size_t reactionIndex,
        const reaction::LogicalReaction &reaction
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
            } else {
                // A,B If not calling with an AD type, calculate the reverse rate directly
                reverseRateConstant = calculateReverseRate(reaction, T9);
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
                const unsigned long speciesIndex = m_speciesToIndexMap.at(species);
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

        const T u = static_cast<T>(m_constants.u); // Atomic mass unit in grams
        const T N_A = static_cast<T>(m_constants.Na); // Avogadro's number in mol^-1
        const T c = static_cast<T>(m_constants.c); // Speed of light in cm/s

        // --- SINGLE LOOP OVER ALL REACTIONS ---
        for (size_t reactionIndex = 0; reactionIndex < m_reactions.size(); ++reactionIndex) {
            const auto& reaction = m_reactions[reactionIndex];

            // 1. Calculate forward reaction rate
            const T forwardMolarReactionFlow = screeningFactors[reactionIndex] *
                calculateMolarReactionFlow<T>(reaction, Y, T9, rho);

            // 2. Calculate reverse reaction rate
            T reverseMolarFlow = calculateReverseMolarReactionFlow<T>(
                T9,
                rho,
                screeningFactors,
                Y,
                reactionIndex,
                reaction
            );

            const T molarReactionFlow = forwardMolarReactionFlow - reverseMolarFlow; // Net molar reaction flow

            // 3. Use the rate to update all relevant species derivatives (dY/dt)
            for (size_t speciesIndex = 0; speciesIndex < m_networkSpecies.size(); ++speciesIndex) {
                const T nu_ij = static_cast<T>(m_stoichiometryMatrix(speciesIndex, reactionIndex));
                result.dydt[speciesIndex] += threshold_flag * nu_ij * molarReactionFlow;
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

        // --- Calculate the molar reaction rate (in units of [s^-1][cm^3(N-1)][mol^(1-N)] for N reactants) ---
        const T k_reaction = reaction.calculate_rate(T9);

        // --- Cound the number of each reactant species to account for species multiplicity ---
        std::unordered_map<std::string, int> reactant_counts;
        reactant_counts.reserve(reaction.reactants().size());
        for (const auto& reactant : reaction.reactants()) {
            reactant_counts[std::string(reactant.name())]++;
        }
        const int totalReactants = static_cast<int>(reaction.reactants().size());

        // --- Accumulator for the molar concentration ---
        auto molar_concentration_product = static_cast<T>(1.0);

        // --- Loop through each unique reactant species and calculate the molar concentration for that species then multiply that into the accumulator ---
        for (const auto& [species_name, count] : reactant_counts) {
            // --- Resolve species to molar abundance ---
            // PERF: Could probably optimize out this lookup
            const auto species_it = m_speciesToIndexMap.find(m_networkSpeciesMap.at(species_name));
            const size_t species_index = species_it->second;
            const T Yi = Y[species_index];

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