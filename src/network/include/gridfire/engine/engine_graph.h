#pragma once

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/composition.h"
#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "gridfire/network.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"

#include <string>
#include <unordered_map>
#include <vector>

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "cppad/cppad.hpp"

// PERF: The function getNetReactionStoichiometry returns a map of species to their stoichiometric coefficients for a given reaction.
//       this makes extra copies of the species, which is not ideal and could be optimized further.
//       Even more relevant is the member m_reactionIDMap which makes copies of a REACLIBReaction for each reaction ID.
//       REACLIBReactions are quite large data structures, so this could be a performance bottleneck.

namespace gridfire {
    typedef CppAD::AD<double> ADDouble; ///< Alias for CppAD AD type for double precision.

    using fourdst::config::Config;
    using fourdst::logging::LogManager;
    using fourdst::constant::Constants;

    static constexpr double MIN_DENSITY_THRESHOLD = 1e-18;
    static constexpr double MIN_ABUNDANCE_THRESHOLD = 1e-18;
    static constexpr double MIN_JACOBIAN_THRESHOLD = 1e-24;

    class GraphEngine final : public DynamicEngine{
    public:
        explicit GraphEngine(const fourdst::composition::Composition &composition);
        explicit GraphEngine(reaction::REACLIBLogicalReactionSet reactions);

        StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const override;

        void generateJacobianMatrix(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) override;

        void generateStoichiometryMatrix() override;

        double calculateMolarReactionFlow(
            const reaction::Reaction& reaction,
            const std::vector<double>&Y,
            const double T9,
            const double rho
        ) const override;

        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;
        [[nodiscard]] const reaction::REACLIBLogicalReactionSet& getNetworkReactions() const override;
        [[nodiscard]] double getJacobianMatrixEntry(
            const int i,
            const int j
        ) const override;
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, int> getNetReactionStoichiometry(
            const reaction::Reaction& reaction
        ) const;
        [[nodiscard]] int getStoichiometryMatrixEntry(
            const int speciesIndex,
            const int reactionIndex
        ) const override;
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] bool involvesSpecies(
            const fourdst::atomic::Species& species
        ) const;

        void exportToDot(
            const std::string& filename
        ) const;
        void exportToCSV(
            const std::string& filename
        ) const;


    private:
        reaction::REACLIBLogicalReactionSet m_reactions; ///< Set of REACLIB reactions in the network.
        std::unordered_map<std::string_view, reaction::Reaction*> m_reactionIDMap; ///< Map from reaction ID to REACLIBReaction. //PERF: This makes copies of REACLIBReaction and could be a performance bottleneck.

        std::vector<fourdst::atomic::Species> m_networkSpecies; ///< Vector of unique species in the network.
        std::unordered_map<std::string_view, fourdst::atomic::Species> m_networkSpeciesMap; ///< Map from species name to Species object.
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap; ///< Map from species to their index in the stoichiometry matrix.

        boost::numeric::ublas::compressed_matrix<int> m_stoichiometryMatrix; ///< Stoichiometry matrix (species x reactions).
        boost::numeric::ublas::compressed_matrix<double> m_jacobianMatrix; ///< Jacobian matrix (species x species).

        CppAD::ADFun<double> m_rhsADFun; ///< CppAD function for the right-hand side of the ODE.

        Config& m_config = Config::getInstance();
        Constants& m_constants = Constants::getInstance(); ///< Access to physical constants.
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");

    private:
        void syncInternalMaps();
        void collectNetworkSpecies();
        void populateReactionIDMap();
        void populateSpeciesToIndexMap();
        void reserveJacobianMatrix();
        void recordADTape();

        [[nodiscard]] bool validateConservation() const;
        void validateComposition(
            const fourdst::composition::Composition &composition,
            double culling,
            double T9
        );

        template <IsArithmeticOrAD T>
        T calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<T> &Y,
            const T T9,
            const T rho
        ) const;

        template<IsArithmeticOrAD T>
        StepDerivatives<T> calculateAllDerivatives(
            const std::vector<T> &Y_in,
            T T9,
            T rho
        ) const;

        StepDerivatives<double> calculateAllDerivatives(
            const std::vector<double>& Y_in,
            const double T9,
            const double rho
        ) const;

        StepDerivatives<ADDouble> calculateAllDerivatives(
            const std::vector<ADDouble>& Y_in,
            const ADDouble T9,
            const ADDouble rho
        ) const;
    };


    template<IsArithmeticOrAD T>
    StepDerivatives<T> GraphEngine::calculateAllDerivatives(
        const std::vector<T> &Y_in, T T9, T rho) const {

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
            const T molarReactionFlow = calculateMolarReactionFlow<T>(reaction, Y, T9, rho);

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