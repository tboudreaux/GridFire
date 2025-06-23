#pragma once

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/composition.h"

#include "gridfire/network.h"
#include "gridfire/reaclib.h"

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

/**
 * @file netgraph.h
 * @brief Defines the GraphNetwork class for representing and evolving nuclear reaction networks as a heavily connected graph.
 *
 * This file provides the GraphNetwork class, which implements a graph-based nuclear reaction network
 * using REACLIB reactions and supports both stiff and non-stiff ODE integration. The network is constructed
 * from a composition and can be queried for species, reactions, and stoichiometry.
 *
 * This is a general nuclear reaction network; however, it does not currently manage reverse reactions, weak reactions,
 * or other reactions which become much more relevant for more extreme astrophysical sources.
 *
 * @see gridfire::GraphNetwork
 */

namespace gridfire {

    /**
     * @brief Concept to check if a type is either double or CppAD::AD<double>.
     *
     * Used to enable templated functions for both arithmetic and automatic differentiation types.
     */
    template<typename T>
    concept IsArithmeticOrAD = std::is_same_v<T, double> || std::is_same_v<T, CppAD::AD<double>>;

    /// Minimum density threshold (g/cm^3) below which reactions are ignored.
    static constexpr double MIN_DENSITY_THRESHOLD = 1e-18;
    /// Minimum abundance threshold below which reactions are ignored.
    static constexpr double MIN_ABUNDANCE_THRESHOLD = 1e-18;
    /// Minimum Jacobian entry threshold for sparsity.
    static constexpr double MIN_JACOBIAN_THRESHOLD = 1e-24;

    /**
     * @brief Graph-based nuclear reaction network using REACLIB reactions.
     *
     * GraphNetwork constructs a reaction network from a given composition, builds the associated
     * stoichiometry and Jacobian matrices, and provides methods for evaluating the network's evolution
     * using ODE solvers. It supports both stiff and non-stiff integration and can be queried for
     * network species, reactions, and stoichiometry.
     *
     * @note Currently does not handle reverse reactions, weak reactions, or other complex reaction types. These will be added in future versions.
     *
     * Example usage:
     * @code
     * serif::composition::Composition comp = ...;
     * serif::network::GraphNetwork net(comp);
     * serif::network::NetIn input;
     * input.composition = comp;
     * input.temperature = 1.5e9;
     * input.density = 1e6;
     * input.tMax = 1.0;
     * input.dt0 = 1e-3;
     * serif::network::NetOut result = net.evaluate(input);
     * @endcode
     */
    class GraphNetwork final : public Network {
    public:
        /**
         * @brief Construct a GraphNetwork from a composition.
         * @param composition The composition specifying the initial isotopic abundances.
         */
        explicit GraphNetwork(const fourdst::composition::Composition &composition);

        /**
         * @brief Construct a GraphNetwork from a composition with culling and temperature.
         * @param composition The composition specifying the initial isotopic abundances.
         * @param cullingThreshold specific reaction rate threshold to cull reactions below.
         * @param T9 Temperature in units of 10^9 K where culling is evaluated at.
         *
         * @see serif::network::build_reaclib_nuclear_network
         */
        explicit GraphNetwork(const fourdst::composition::Composition &composition,
                              const double cullingThreshold, const double T9);

        explicit GraphNetwork(const reaclib::REACLIBReactionSet& reactions);

        /**
         * @brief Evolve the network for the given input conditions.
         *
         * This is the primary entry point for users of GridFire. This function integrates the network ODEs using either a stiff or non-stiff solver,
         * depending on the detected stiffness of the system.
         *
         * @param netIn The input structure containing composition, temperature, density, and integration parameters.
         * @return NetOut The output structure containing the final composition, number of steps, and energy.
         *
         * Example:
         * @code
         * serif::network::NetIn input;
         * // ... set up input ...
         * serif::network::NetOut result = net.evaluate(input);
         * @endcode
         */
        NetOut evaluate(const NetIn &netIn) override;

        /**
         * @brief Get the vector of unique species in the network.
         * @return Reference to the vector of species.
         *
         * Example:
         * @code
         * const auto& species = net.getNetworkSpecies();
         * for (const auto& sp : species) {
         *     std::cout << sp.name() << std::endl;
         * }
         * @endcode
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const;

        /**
         * @brief Get the set of REACLIB reactions in the network.
         * @return Reference to the set of reactions.
         */
        [[nodiscard]] const reaclib::REACLIBReactionSet& getNetworkReactions() const;

        /**
         * @brief Get the net stoichiometric coefficients for all species in a reaction.
         *
         * Returns a map from species to their net stoichiometric coefficient (products minus reactants).
         *
         * @param reaction The REACLIB reaction to analyze.
         * @return Map from species to stoichiometric coefficient.
         *
         * @throws std::runtime_error If a species in the reaction is not found in the network.
         *
         * Example:
         * @code
         * for (const auto& reaction : net.getNetworkReactions()) {
         *     auto stoichiometry = net.getNetReactionStoichiometry(reaction);
         *     // ...
         * }
         * @endcode
         */
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, int> getNetReactionStoichiometry(
            const reaclib::REACLIBReaction &reaction) const;

        /**
         * @brief Check if a species is present in the network.
         * @param species The species to check.
         * @return True if the species is present, false otherwise.
         *
         * Example:
         * @code
         * if (net.involvesSpecies(mySpecies)) { ... }
         * @endcode
         */
        [[nodiscard]] bool involvesSpecies(const fourdst::atomic::Species& species) const;

        /**
         * @brief Detect cycles in the reaction network (not implemented).
         * @return Vector of cycles, each represented as a vector of reaction IDs.
         *
         * @note This is not yet implemented; however, it will allow for easy detection of things like the CNO cycle.
         * @deprecated Not implemented.
         */
        [[deprecated("not implimented")]] [[nodiscard]] std::vector<std::vector<std::string>> detectCycles() const = delete;

        /**
         * @brief Perform a topological sort of the reactions (not implemented).
         * @return Vector of reaction IDs in topologically sorted order.
         * @deprecated Not implemented.
         */
        [[deprecated("not implimented")]] [[nodiscard]] std::vector<std::string> topologicalSortReactions() const = delete;

        /**
         * @brief Export the network to a DOT file for visualization.
         * @param filename The name of the output DOT file.
         */
        void exportToDot(const std::string& filename) const;

    private:
        reaclib::REACLIBReactionSet m_reactions; ///< Set of REACLIB reactions in the network.
        std::unordered_map<std::string_view, const reaclib::REACLIBReaction> m_reactionIDMap; ///< Map from reaction ID to REACLIBReaction. //PERF: This makes copies of REACLIBReaction and could be a performance bottleneck.

        std::vector<fourdst::atomic::Species> m_networkSpecies; ///< Vector of unique species in the network.
        std::unordered_map<std::string_view, fourdst::atomic::Species> m_networkSpeciesMap; ///< Map from species name to Species object.
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap; ///< Map from species to their index in the stoichiometry matrix.

        boost::numeric::ublas::compressed_matrix<int> m_stoichiometryMatrix; ///< Stoichiometry matrix (species x reactions).
        boost::numeric::ublas::compressed_matrix<double> m_jacobianMatrix; ///< Jacobian matrix (species x species).

        CppAD::ADFun<double> m_rhsADFun; ///< CppAD function for the right-hand side of the ODE.

        /**
         * @brief Functor for ODE right-hand side evaluation.
         *
         * Used by ODE solvers to compute dY/dt and energy generation rate. This is the only functor used in the non-NSE case.
         */
        struct ODETerm {
            GraphNetwork& m_network; ///< Reference to the network
            const double m_T9;       ///< Temperature in 10^9 K
            const double m_rho;      ///< Density in g/cm^3
            const size_t m_numSpecies; ///< Number of species

            /**
             * @brief Construct an ODETerm functor.
             * @param network Reference to the GraphNetwork.
             * @param T9 Temperature in 10^9 K.
             * @param rho Density in g/cm^3.
             */
            ODETerm(GraphNetwork& network, const double T9, double rho) :
             m_network(network), m_T9(T9), m_rho(rho), m_numSpecies(m_network.m_networkSpecies.size()) {}

            /**
             * @brief Compute dY/dt and energy rate for the ODE solver.
             * @param Y Input vector of abundances.
             * @param dYdt Output vector for derivatives (resized).
             *
             * @note this functor does not need auto differentiation to it called the <double> version of calculateAllDerivatives.
             */
            void operator()(const boost::numeric::ublas::vector<double>&Y, boost::numeric::ublas::vector<double>& dYdt, double /* t */) const {
                const std::vector<double> y(Y.begin(), m_numSpecies + Y.begin());

                StepDerivatives<double> derivatives = m_network.calculateAllDerivatives<double>(y, m_T9, m_rho);

                dYdt.resize(m_numSpecies + 1);
                std::ranges::copy(derivatives.dydt, dYdt.begin());
                dYdt(m_numSpecies) = derivatives.specificEnergyRate; // Last element is the specific energy rate
            }
        };

        /**
         * @brief Functor for Jacobian evaluation for stiff ODE solvers. (used in the NSE case)
         */
        struct JacobianTerm {
            GraphNetwork& m_network; ///< Reference to the network
            const double m_T9;       ///< Temperature in 10^9 K
            const double m_rho;      ///< Density in g/cm^3
            const size_t m_numSpecies; ///< Number of species

            /**
             * @brief Construct a JacobianTerm functor.
             * @param network Reference to the GraphNetwork.
             * @param T9 Temperature in 10^9 K.
             * @param rho Density in g/cm^3.
             */
            JacobianTerm(GraphNetwork& network, const double T9, const double rho) :
                m_network(network), m_T9(T9), m_rho(rho), m_numSpecies(m_network.m_networkSpecies.size()) {}

            /**
             * @brief Compute the Jacobian matrix for the ODE solver.
             * @param Y Input vector of abundances.
             * @param J Output matrix for the Jacobian (resized).
             */
            void operator()(const boost::numeric::ublas::vector<double>& Y, boost::numeric::ublas::matrix<double>& J, double /* t */, boost::numeric::ublas::vector<double>& /* dfdt */) const {
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

        /**
         * @brief Struct holding derivatives for a single ODE step.
         * @tparam T Scalar type (double or CppAD::AD<double>).
         */
        template <IsArithmeticOrAD T>
        struct StepDerivatives {
            std::vector<T> dydt; ///< Derivatives of abundances.
            T specificEnergyRate = T(0.0); ///< Specific energy generation rate.
        };

    private:
        /**
         * @brief Synchronize all internal maps and matrices after network changes.
         *
         * Called after the reaction set is updated to rebuild all derived data structures.
         *
         * @note This method must be called any time the network topology changes.
         */
        void syncInternalMaps();

        /**
         * @brief Collect all unique species from the reaction set.
         *
         * Populates m_networkSpecies and m_networkSpeciesMap.
         * @throws std::runtime_error If a species is not found in the global atomic species database.
         *
         * @note Called by syncInternalMaps() to ensure the species list is up-to-date.
         */
        void collectNetworkSpecies();

        /**
         * @brief Populate the reaction ID map from the reaction set.
         *
         * @note Called by syncInternalMaps() to ensure the reaction ID map is up-to-date.
         */
        void populateReactionIDMap();

        /**
         * @brief Populate the species-to-index map for matrix construction.
         *
         * @note Called by syncInternalMaps() to ensure the species-to-index map is up-to-date.
         */
        void populateSpeciesToIndexMap();

        /**
         * @brief Reserve and resize the Jacobian matrix.
         *
         * @note Called by syncInternalMaps() to ensure the Jacobian matrix is ready for use.
         */
        void reserveJacobianMatrix();

        /**
         * @brief Record the CppAD tape for automatic differentiation.
         * @throws std::runtime_error If there are no species in the network.
         *
         * @note Called by syncInternalMaps() to ensure the AD tape is recorded for the current network state.
         */
        void recordADTape();

        // --- Validation methods ---

        /**
         * @brief Validate mass and charge conservation for all reactions.
         * @return True if all reactions conserve mass and charge, false otherwise.
         *
         * @note Logs errors for any violations.
         */
        [[nodiscard]] bool validateConservation() const;

        /**
         * @brief Validate and update the network composition if needed.
         *
         * If the composition or culling/temperature has changed, rebuilds the reaction set and synchronizes maps.
         * This is primarily used to enable caching in dynamic network situations where the composition, temperature, and density
         * may change but the relevant reaction set remains equivalent.
         *
         * @param composition The composition to validate.
         * @param culling Culling threshold.
         * @param T9 Temperature in 10^9 K.
         */
        void validateComposition(const fourdst::composition::Composition &composition, double culling, double T9);

        // --- Simple Derivative Calculations ---

        /**
         * @brief Calculate all derivatives (dY/dt and energy rate) for the current state.
         * @tparam T Scalar type (double or CppAD::AD<double>).
         * @param Y Vector of abundances.
         * @param T9 Temperature in 10^9 K.
         * @param rho Density in g/cm^3.
         * @return StepDerivatives<T> containing dY/dt and energy rate.
         */
        template <IsArithmeticOrAD T>
        StepDerivatives<T> calculateAllDerivatives(const std::vector<T>& Y, T T9, T rho) const;

        // --- Generate the system matrices ---

        /**
         * @brief Generate the stoichiometry matrix for the network.
         *
         * Populates m_stoichiometryMatrix based on the current reactions and species.
         * @throws std::runtime_error If a species is not found in the species-to-index map.
         */
        void generateStoichiometryMatrix();

        /**
         * @brief Generate the Jacobian matrix for the network.
         *
         * Populates m_jacobianMatrix using automatic differentiation via the precomputed CppAD tape.
         *
         * @param Y Vector of abundances.
         * @param T9 Temperature in 10^9 K.
         * @param rho Density in g/cm^3.
         */
        void generateJacobianMatrix(const std::vector<double>& Y, double T9, const double rho);

        /**
         * @brief Calculate the right-hand side (dY/dt) for the ODE system.
         * @tparam GeneralScalarType Scalar type (double or CppAD::AD<double>).
         * @param Y Vector of abundances.
         * @param T9 Temperature in 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Vector of dY/dt values.
         */
        template <IsArithmeticOrAD GeneralScalarType>
        std::vector<GeneralScalarType> calculateRHS(const std::vector<GeneralScalarType> &Y, const GeneralScalarType T9,
                                                    GeneralScalarType rho) const;

        /**
         * @brief Calculate the reaction rate for a given reaction in dimensions of particles per cm^3 per second.
         * @tparam GeneralScalarType Scalar type (double or CppAD::AD<double>).
         * @param reaction The REACLIB reaction.
         * @param Y Vector of abundances.
         * @param T9 Temperature in 10^9 K.
         * @param rho Density in g/cm^3.
         * @return Reaction rate.
         * @throws std::runtime_error If a reactant species is not found in the species-to-index map.
         *
         * @note reaclib uses molar units that vary depending on the number of reactants, It's pretty straightforward
         * to convert these into particle based units. Specifically, we just need to divide the final result by
         * Avogadro's number raised to the number of reactants - 1;
         */
        template <IsArithmeticOrAD GeneralScalarType>
        GeneralScalarType calculateReactionRate(const reaclib::REACLIBReaction &reaction,
                                                const std::vector<GeneralScalarType> &Y, const GeneralScalarType T9,
                                                const GeneralScalarType rho) const;

        /**
         * @brief Detect if the network is stiff and select the appropriate ODE solver.
         *
         * Heuristically determines stiffness based on the ratio of timescales. The stiffness heuristic is as follows:
         *
         * 1. For each species, compute the timescale as |Y_i / (dY_i/dt)|, where Y_i is the abundance and dY_i/dt is its time derivative.
         * 2. Find the minimum and maximum timescales across all species.
         * 3. Compute the stiffness ratio as (maximum timescale) / (minimum timescale).
         * 4. If the stiffness ratio exceeds a threshold (default: 1e6), the system is considered stiff and a stiff ODE solver is used.
         *    Otherwise, a non-stiff ODE solver is used.
         *
         * This heuristic is based on the observation that stiff systems have widely varying timescales, which can cause explicit
         * solvers to become unstable or inefficient. The threshold can be tuned based on the characteristics of the network.
         *
         * @param netIn The input structure.
         * @param T9 Temperature in 10^9 K.
         * @param numSpecies Number of species.
         * @param Y Vector of abundances.
         */
        void detectStiff(const NetIn &netIn, double T9, unsigned long numSpecies, const boost::numeric::ublas::vector<double>& Y);

    };


    template<IsArithmeticOrAD T>
    GraphNetwork::StepDerivatives<T> GraphNetwork::calculateAllDerivatives(
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
        const auto &constants  = fourdst::constant::Constants::getInstance();

        const auto u = constants.get("u"); // Atomic mass unit in g/mol
        const auto uValue = static_cast<GeneralScalarType>(u.value); // Convert to double for calculations
        const GeneralScalarType k_reaction = reaction.calculate_rate(T9);

        auto reactant_product_or_particle_densities = static_cast<GeneralScalarType>(1.0);

        std::unordered_map<std::string, int> reactant_counts;
        reactant_counts.reserve(reaction.reactants().size());
        for (const auto& reactant : reaction.reactants()) {
            reactant_counts[std::string(reactant.name())]++;
        }

        const auto minAbundanceThreshold = static_cast<GeneralScalarType>(MIN_ABUNDANCE_THRESHOLD);
        const auto minDensityThreshold = static_cast<GeneralScalarType>(MIN_DENSITY_THRESHOLD);

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

            auto atomicMassAMU = static_cast<GeneralScalarType>(m_networkSpecies[species_index].mass());

            // Convert to number density
            GeneralScalarType ni;
            const GeneralScalarType denominator = atomicMassAMU * uValue;
            assert (denominator > 0.0);
            ni = (Yi * rho) / (denominator);

            reactant_product_or_particle_densities *= ni;

            // Apply factorial correction for identical reactions
            if (count > 1) {
                reactant_product_or_particle_densities /= static_cast<GeneralScalarType>(std::tgamma(static_cast<double>(count + 1))); // Gamma function for factorial
            }
        }
        const GeneralScalarType Na = static_cast<GeneralScalarType>(constants.get("N_a").value); // Avogadro's number in mol^-1
        const int numReactants = static_cast<int>(reaction.reactants().size());
        auto molarCorrectionFactor = static_cast<GeneralScalarType>(1.0); // No correction needed for single reactant reactions
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