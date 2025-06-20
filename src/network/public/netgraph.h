#pragma once

#include "atomicSpecies.h"
#include "composition.h"
#include "network.h"
#include "reaclib.h"

#include <string>
#include <unordered_map>
#include <vector>

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "cppad/cppad.hpp"

#include "quill/LogMacros.h"


// PERF: The function getNetReactionStoichiometry returns a map of species to their stoichiometric coefficients for a given reaction.
//       this makes extra copies of the species, which is not ideal and could be optimized further.
//       Even more relevant is the member m_reactionIDMap which makes copies of a REACLIBReaction for each reaction ID.
//       REACLIBReactions are quite large data structures, so this could be a performance bottleneck.

namespace serif::network {
    // Define a concept to check if a type is one of our allowed scalar types
    template<typename T>
    concept IsArithmeticOrAD = std::is_same_v<T, double> || std::is_same_v<T, CppAD::AD<double>>;

    static constexpr double MIN_DENSITY_THRESHOLD = 1e-18;
    static constexpr double MIN_ABUNDANCE_THRESHOLD = 1e-18;
    static constexpr double MIN_JACOBIAN_THRESHOLD = 1e-24;

    class GraphNetwork final : public Network {
    public:
        explicit GraphNetwork(const serif::composition::Composition &composition);
        explicit GraphNetwork(const serif::composition::Composition &composition,
                              const double cullingThreshold, const double T9);

        NetOut evaluate(const NetIn &netIn) override;

        [[nodiscard]] const std::vector<serif::atomic::Species>& getNetworkSpecies() const;
        [[nodiscard]] const reaclib::REACLIBReactionSet& getNetworkReactions() const;
        [[nodiscard]] std::unordered_map<serif::atomic::Species, int> getNetReactionStoichiometry(
            const reaclib::REACLIBReaction &reaction) const;

        [[nodiscard]] bool involvesSpecies(const serif::atomic::Species& species) const;

        [[deprecated("not implimented")]] std::vector<std::vector<std::string>> detectCycles() const;
        [[deprecated("not implimented")]] std::vector<std::string> topologicalSortReactions() const;
    private:
        reaclib::REACLIBReactionSet m_reactions; ///< Set of REACLIB reactions for the network.

        std::vector<serif::atomic::Species> m_networkSpecies; ///< The species in the network.
        std::unordered_map<std::string_view, serif::atomic::Species> m_networkSpeciesMap;
        std::unordered_map<std::string_view, const reaclib::REACLIBReaction> m_reactionIDMap; ///< Map of reaction IDs to REACLIB reactions.

        boost::numeric::ublas::compressed_matrix<int> m_stoichiometryMatrix; ///< Stoichiometry matrix for the network.
        boost::numeric::ublas::compressed_matrix<double> m_jacobianMatrix; ///< Jacobian matrix for the network.
        std::unordered_map<serif::atomic::Species, size_t> m_speciesToIndexMap; ///< Map of species to their index in the stoichiometry matrix.

        CppAD::ADFun<double> m_rhsADFun; ///< AD tape

        struct ODETerm {
            GraphNetwork& m_network; ///< Reference to the network
            const double m_T9;
            const double m_rho;
            const size_t m_numSpecies;

            ODETerm(GraphNetwork& network, const double T9, double rho) :
             m_network(network), m_T9(T9), m_rho(rho), m_numSpecies(m_network.m_networkSpecies.size()) {}

            void operator()(const boost::numeric::ublas::vector<double>&Y, boost::numeric::ublas::vector<double>& dYdt, double /* t */) {
                const std::vector<double> y(Y.begin(), Y.begin() + m_numSpecies);

                StepDerivatives<double> derivatives = m_network.calculateAllDerivatives<double>(y, m_T9, m_rho);

                dYdt.resize(m_numSpecies + 1);
                std::copy(derivatives.dydt.begin(), derivatives.dydt.end(), dYdt.begin());
                dYdt(m_numSpecies) = derivatives.specificEnergyRate; // Last element is the specific energy rate
            }
        };

        struct JacobianTerm {
            GraphNetwork& m_network;
            const double m_T9;
            const double m_rho;
            const size_t m_numSpecies;

            JacobianTerm(GraphNetwork& network, const double T9, const double rho) :
                m_network(network), m_T9(T9), m_rho(rho), m_numSpecies(m_network.m_networkSpecies.size()) {}

            void operator()(const boost::numeric::ublas::vector<double>& Y, boost::numeric::ublas::matrix<double>& J, double /* t */, boost::numeric::ublas::vector<double>& /* dfdt */) {
                const std::vector<double> y_species(Y.begin(), Y.begin() + m_numSpecies);

                m_network.generateJacobianMatrix(y_species, m_T9, m_rho);

                J.resize(m_numSpecies + 1, m_numSpecies + 1);
                J.clear(); // Zero out the matrix

                // PERF: Could probably be optimized
                for (auto it1 = m_network.m_jacobianMatrix.begin1(); it1 != m_network.m_jacobianMatrix.end1(); ++it1) {
                    for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
                        J(it2.index1(), it2.index2()) = *it2;
                    }
                }
            }
        };

        template <IsArithmeticOrAD T>
        struct StepDerivatives {
            std::vector<T> dydt;
            T specificEnergyRate = T(0.0);
        };


    private:
        void syncInternalMaps();
        void collectNetworkSpecies();
        void populateReactionIDMap();
        void populateSpeciesToIndexMap();
        void reserveJacobianMatrix();
        void recordADTape();

        // --- Validation methods ---
        bool validateConservation() const;

        void validateComposition(const serif::composition::Composition &composition, double culling, double T9);

        // --- Simple Derivative Calculations ---
        template <IsArithmeticOrAD T>
        StepDerivatives<T> calculateAllDerivatives(const std::vector<T>& Y, T T9, T rho) const;

        // --- Generate the system matrices ---
        void generateStoichiometryMatrix();
        void generateJacobianMatrix(const std::vector<double>& Y, const double T9, const double rho);

        template <IsArithmeticOrAD GeneralScalarType>
        std::vector<GeneralScalarType> calculateRHS(const std::vector<GeneralScalarType> &Y, const GeneralScalarType T9,
                                                    const GeneralScalarType rho) const;
        template <IsArithmeticOrAD GeneralScalarType>
        GeneralScalarType calculateReactionRate(const reaclib::REACLIBReaction &reaction,
                                                const std::vector<GeneralScalarType> &Y, const GeneralScalarType T9,
                                                const GeneralScalarType rho) const;

        void detectStiff(const NetIn &netIn, double T9, double numSpecies, const boost::numeric::ublas::vector<double>& Y);
    };


    template<IsArithmeticOrAD T>
    typename GraphNetwork::template StepDerivatives<T> GraphNetwork::calculateAllDerivatives(
        const std::vector<T> &Y, T T9, T rho) const {
        StepDerivatives<T> result;
        result.dydt.resize(m_networkSpecies.size(), static_cast<T>(0.0));

        if (rho < MIN_DENSITY_THRESHOLD) {
            return result; // Return zero for everything if density is too low
        }

        const T u = static_cast<T>(m_constants.get("u").value); // Atomic mass unit in grams
        const T MeV_to_erg = static_cast<T>(m_constants.get("MeV_to_erg").value);

        T volumetricEnergyRate = static_cast<T>(0.0); // Accumulator for erg / (cm^3 * s)

        // --- SINGLE LOOP OVER ALL REACTIONS ---
        for (size_t reactionIndex = 0; reactionIndex < m_reactions.size(); ++reactionIndex) {
            const auto& reaction = m_reactions[reactionIndex];

            // 1. Calculate reaction rate
            const T reactionRate = calculateReactionRate(reaction, Y, T9, rho);

            // 2. Use the rate to update all relevant species derivatives (dY/dt)
            for (size_t speciesIndex = 0; speciesIndex < m_networkSpecies.size(); ++speciesIndex) {
                const T nu_ij = static_cast<T>(m_stoichiometryMatrix(speciesIndex, reactionIndex));
                if (nu_ij != 0) {
                    const T speciesAtomicMassAMU = static_cast<T>(m_networkSpecies[speciesIndex].mass());
                    const T speciesAtomicMassGrams = speciesAtomicMassAMU * u;
                    result.dydt[speciesIndex] += (nu_ij * reactionRate * speciesAtomicMassGrams) / rho;
                }
            }

            // 3. Use the same rate to update the energy generation rate
            const T q_value_ergs = static_cast<T>(reaction.qValue()) * MeV_to_erg;
            volumetricEnergyRate += reactionRate * q_value_ergs;
        }

        result.specificEnergyRate = volumetricEnergyRate / rho;

        return result;
    }


    template <IsArithmeticOrAD GeneralScalarType>
    GeneralScalarType GraphNetwork::calculateReactionRate(const reaclib::REACLIBReaction &reaction, const std::vector<GeneralScalarType> &Y,
        const GeneralScalarType T9, const GeneralScalarType rho) const {
        const auto &constants  = serif::constant::Constants::getInstance();

        const auto u = constants.get("u"); // Atomic mass unit in g/mol
        const GeneralScalarType uValue = static_cast<GeneralScalarType>(u.value); // Convert to double for calculations
        const GeneralScalarType k_reaction = reaction.calculate_rate(T9);

        GeneralScalarType reactant_product_or_particle_densities = static_cast<GeneralScalarType>(1.0);

        std::unordered_map<std::string, int> reactant_counts;
        reactant_counts.reserve(reaction.reactants().size());
        for (const auto& reactant : reaction.reactants()) {
            reactant_counts[std::string(reactant.name())]++;
        }

        const GeneralScalarType minAbundanceThreshold = static_cast<GeneralScalarType>(MIN_ABUNDANCE_THRESHOLD);
        const GeneralScalarType minDensityThreshold = static_cast<GeneralScalarType>(MIN_DENSITY_THRESHOLD);

        if (rho < minDensityThreshold) {
            // LOG_INFO(m_logger, "Density is below the minimum threshold ({} g/cm^3), returning zero reaction rate for reaction '{}'.",
            //          CppAD::Value(rho), reaction.id()); // CppAD::Value causes a compilation error here, not clear why...
            return static_cast<GeneralScalarType>(0.0); // If density is below a threshold, return zero rate.
        }

        for (const auto& [species_name, count] : reactant_counts) {
            auto species_it = m_speciesToIndexMap.find(m_networkSpeciesMap.at(species_name));
            if (species_it == m_speciesToIndexMap.end()) {
                LOG_ERROR(m_logger, "Reactant species '{}' not found in species to index map for reaction '{}'.",
                         species_name, reaction.id());
                throw std::runtime_error("Reactant species not found in species to index map: " + species_name);
            }
            const size_t species_index = species_it->second;
            const GeneralScalarType Yi = Y[species_index];

            if (Yi < minAbundanceThreshold) {
                return static_cast<GeneralScalarType>(0.0); // If any reactant is below a threshold, return zero rate.
            }

            GeneralScalarType atomicMassAMU = static_cast<GeneralScalarType>(m_networkSpecies[species_index].mass());

            // Convert to number density
            GeneralScalarType ni;
            const GeneralScalarType denominator = atomicMassAMU * uValue;
            if (denominator > minDensityThreshold) {
                ni = (Yi * rho) / (denominator);
            } else {
                ni = static_cast<GeneralScalarType>(0.0);
            }

            reactant_product_or_particle_densities *= ni;

            // Apply factorial correction for identical reactions
            if (count > 1) {
                reactant_product_or_particle_densities /= static_cast<GeneralScalarType>(std::tgamma(static_cast<double>(count + 1))); // Gamma function for factorial
            }
        }
        const GeneralScalarType Na = static_cast<GeneralScalarType>(constants.get("N_a").value); // Avogadro's number in mol^-1
        const int numReactants = static_cast<int>(reaction.reactants().size());
        GeneralScalarType molarCorrectionFactor = static_cast<GeneralScalarType>(1.0); // No correction needed for single reactant reactions
        if (numReactants > 1) {
            molarCorrectionFactor = CppAD::pow(Na, static_cast<GeneralScalarType>(reaction.reactants().size() - 1));
        }
        return (reactant_product_or_particle_densities * k_reaction) / molarCorrectionFactor; // reaction rate in per volume per time (particles/cm^3/s)
    }

    template <IsArithmeticOrAD GeneralScalarType>
    std::vector<GeneralScalarType> GraphNetwork::calculateRHS(
        const std::vector<GeneralScalarType>& Y,
        const GeneralScalarType T9,
        const GeneralScalarType rho) const {

        auto derivatives = calculateAllDerivatives<GeneralScalarType>(Y, T9, rho);
        return derivatives.dydt;
    }

};