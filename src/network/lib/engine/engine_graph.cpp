#include "gridfire/engine/engine_graph.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/network.h"
#include "gridfire/screening/screening_types.h"

#include "fourdst/composition/species.h"
#include "fourdst/composition/atomicSpecies.h"

#include "quill/LogMacros.h"

#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fstream>
#include <ranges>

#include <boost/numeric/odeint.hpp>


namespace gridfire {
    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition
    ):
     m_reactions(build_reaclib_nuclear_network(composition, false)),
     m_partitionFunction(std::make_unique<partition::GroundStatePartitionFunction>()){
        syncInternalMaps();
        precomputeNetwork();
    }

    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition,
        const partition::PartitionFunction& partitionFunction) :
    m_reactions(build_reaclib_nuclear_network(composition, false)),
    m_partitionFunction(partitionFunction.clone())  // Clone the partition function to ensure ownership
    {
        syncInternalMaps();
        precomputeNetwork();
    }

    GraphEngine::GraphEngine(
        const reaction::LogicalReactionSet &reactions
    ) :
    m_reactions(reactions) {
        syncInternalMaps();
        precomputeNetwork();
    }

    StepDerivatives<double> GraphEngine::calculateRHSAndEnergy(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        if (m_usePrecomputation) {
            std::vector<double> bare_rates;
            bare_rates.reserve(m_reactions.size());
            for (const auto& reaction: m_reactions) {
                bare_rates.push_back(reaction.calculate_rate(T9));
            }

            // --- The public facing interface can always use the precomputed version since taping is done internally ---
            return calculateAllDerivativesUsingPrecomputation(Y, bare_rates, T9, rho);
        } else {
            return calculateAllDerivatives<double>(Y, T9, rho);
        }
    }


    void GraphEngine::syncInternalMaps() {
        collectNetworkSpecies();
        populateReactionIDMap();
        populateSpeciesToIndexMap();
        generateStoichiometryMatrix();
        reserveJacobianMatrix();
        recordADTape();
    }

    // --- Network Graph Construction Methods ---
    void GraphEngine::collectNetworkSpecies() {
        m_networkSpecies.clear();
        m_networkSpeciesMap.clear();

        std::set<std::string_view> uniqueSpeciesNames;

        for (const auto& reaction: m_reactions) {
            for (const auto& reactant: reaction.reactants()) {
                uniqueSpeciesNames.insert(reactant.name());
            }
            for (const auto& product: reaction.products()) {
                uniqueSpeciesNames.insert(product.name());
            }
        }

        for (const auto& name: uniqueSpeciesNames) {
            auto it = fourdst::atomic::species.find(std::string(name));
            if (it != fourdst::atomic::species.end()) {
                m_networkSpecies.push_back(it->second);
                m_networkSpeciesMap.insert({name, it->second});
            } else {
                LOG_ERROR(m_logger, "Species '{}' not found in global atomic species database.", name);
                m_logger->flush_log();
                throw std::runtime_error("Species not found in global atomic species database: " + std::string(name));
            }
        }

    }

    void GraphEngine::populateReactionIDMap() {
        LOG_TRACE_L1(m_logger, "Populating reaction ID map for REACLIB graph network (serif::network::GraphNetwork)...");
        m_reactionIDMap.clear();
        for (auto& reaction: m_reactions) {
            m_reactionIDMap.emplace(reaction.id(), &reaction);
        }
        LOG_TRACE_L1(m_logger, "Populated {} reactions in the reaction ID map.", m_reactionIDMap.size());
    }

    void GraphEngine::populateSpeciesToIndexMap() {
        m_speciesToIndexMap.clear();
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            m_speciesToIndexMap.insert({m_networkSpecies[i], i});
        }
    }

    void GraphEngine::reserveJacobianMatrix() {
        // The implementation of this function (and others) constrains this nuclear network to a constant temperature and density during
        // each evaluation.
        size_t numSpecies = m_networkSpecies.size();
        m_jacobianMatrix.clear();
        m_jacobianMatrix.resize(numSpecies, numSpecies, false); // Sparse matrix, no initial values
        LOG_TRACE_L2(m_logger, "Jacobian matrix resized to {} rows and {} columns.",
                 m_jacobianMatrix.size1(), m_jacobianMatrix.size2());
    }

    // --- Basic Accessors and Queries ---
    const std::vector<fourdst::atomic::Species>& GraphEngine::getNetworkSpecies() const {
        // Returns a constant reference to the vector of unique species in the network.
        LOG_TRACE_L3(m_logger, "Providing access to network species vector. Size: {}.", m_networkSpecies.size());
        return m_networkSpecies;
    }

    const reaction::LogicalReactionSet& GraphEngine::getNetworkReactions() const {
        // Returns a constant reference to the set of reactions in the network.
        LOG_TRACE_L3(m_logger, "Providing access to network reactions set. Size: {}.", m_reactions.size());
        return m_reactions;
    }

    bool GraphEngine::involvesSpecies(const fourdst::atomic::Species& species) const {
        // Checks if a given species is present in the network's species map for efficient lookup.
        const bool found = m_networkSpeciesMap.contains(species.name());
        LOG_DEBUG(m_logger, "Checking if species '{}' is involved in the network: {}.", species.name(), found ? "Yes" : "No");
        return found;
    }

    // --- Validation Methods ---
    bool GraphEngine::validateConservation() const {
        LOG_TRACE_L1(m_logger, "Validating mass (A) and charge (Z) conservation across all reactions in the network.");

        for (const auto& reaction : m_reactions) {
            uint64_t totalReactantA = 0;
            uint64_t totalReactantZ = 0;
            uint64_t totalProductA = 0;
            uint64_t totalProductZ = 0;

            // Calculate total A and Z for reactants
            for (const auto& reactant : reaction.reactants()) {
                auto it = m_networkSpeciesMap.find(reactant.name());
                if (it != m_networkSpeciesMap.end()) {
                    totalReactantA += it->second.a();
                    totalReactantZ += it->second.z();
                } else {
                    // This scenario indicates a severe data integrity issue:
                    // a reactant is part of a reaction but not in the network's species map.
                    LOG_ERROR(m_logger, "CRITICAL ERROR: Reactant species '{}' in reaction '{}' not found in network species map during conservation validation.",
                             reactant.name(), reaction.id());
                    return false;
                }
            }

            // Calculate total A and Z for products
            for (const auto& product : reaction.products()) {
                auto it = m_networkSpeciesMap.find(product.name());
                if (it != m_networkSpeciesMap.end()) {
                    totalProductA += it->second.a();
                    totalProductZ += it->second.z();
                } else {
                    // Similar critical error for product species
                    LOG_ERROR(m_logger, "CRITICAL ERROR: Product species '{}' in reaction '{}' not found in network species map during conservation validation.",
                             product.name(), reaction.id());
                    return false;
                }
            }

            // Compare totals for conservation
            if (totalReactantA != totalProductA) {
                LOG_ERROR(m_logger, "Mass number (A) not conserved for reaction '{}': Reactants A={} vs Products A={}.",
                         reaction.id(), totalReactantA, totalProductA);
                return false;
            }
            if (totalReactantZ != totalProductZ) {
                LOG_ERROR(m_logger, "Atomic number (Z) not conserved for reaction '{}': Reactants Z={} vs Products Z={}.",
                         reaction.id(), totalReactantZ, totalProductZ);
                return false;
            }
        }

        LOG_TRACE_L1(m_logger, "Mass (A) and charge (Z) conservation validated successfully for all reactions.");
        return true; // All reactions passed the conservation check
    }

    void GraphEngine::validateComposition(const fourdst::composition::Composition &composition, double culling, double T9) {
        // Check if the requested network has already been cached.
        // PERF: Rebuilding this should be pretty fast but it may be a good point of optimization in the future.
        const reaction::LogicalReactionSet validationReactionSet = build_reaclib_nuclear_network(composition, false);
        // TODO: need some more robust method here to
        //       A. Build the basic network from the composition's species with non zero mass fractions.
        //       B. rebuild a new composition from the reaction set's reactants + products (with the mass fractions from the things that are only products set to 0)
        //       C. Rebuild the reaction set from the new composition
        //       D. Cull reactions where all reactants have mass fractions below the culling threshold.
        //       E. Be careful about maintaining caching through all of this


        // This allows for dynamic network modification while retaining caching for networks which are very similar.
        if (validationReactionSet != m_reactions) {
            LOG_DEBUG(m_logger, "Reaction set not cached. Rebuilding the reaction set for T9={} and culling={}.", T9, culling);
            m_reactions = validationReactionSet;
            syncInternalMaps(); // Re-sync internal maps after updating reactions. Note this will also retrace the AD tape.
        }
    }

    double GraphEngine::calculateReverseRate(
        const reaction::Reaction &reaction,
        const double T9,
        const double expFactor
    ) const {
        double reverseRate = 0.0;
        const double forwardRate = reaction.calculate_rate(T9);

        if (reaction.reactants().size() == 2 && reaction.products().size() == 2) {
            reverseRate = calculateReverseRateTwoBody(reaction, T9, forwardRate, expFactor);
        } else {
            LOG_WARNING(m_logger, "Reverse rate calculation for reactions with more than two reactants or products is not implemented.");
        }
        return reverseRate;
    }

    double GraphEngine::calculateReverseRateTwoBody(
        const reaction::Reaction &reaction,
        const double T9,
        const double forwardRate,
        const double expFactor
    ) const {
        std::vector<double> reactantPartitionFunctions;
        std::vector<double> productPartitionFunctions;

        reactantPartitionFunctions.reserve(reaction.reactants().size());
        productPartitionFunctions.reserve(reaction.products().size());

        std::unordered_map<fourdst::atomic::Species, int> reactantMultiplicity;
        std::unordered_map<fourdst::atomic::Species, int> productMultiplicity;

        reactantMultiplicity.reserve(reaction.reactants().size());
        productMultiplicity.reserve(reaction.products().size());

        for (const auto& reactant : reaction.reactants()) {
            reactantMultiplicity[reactant] += 1;
        }
        for (const auto& product : reaction.products()) {
            productMultiplicity[product] += 1;
        }
        double reactantSymmetryFactor = 1.0;
        double productSymmetryFactor = 1.0;
        for (const auto& count : reactantMultiplicity | std::views::values) {
            reactantSymmetryFactor *= std::tgamma(count + 1);
        }
        for (const auto& count : productMultiplicity | std::views::values) {
            productSymmetryFactor *= std::tgamma(count + 1);
        }
        const double symmetryFactor = reactantSymmetryFactor / productSymmetryFactor;

        // Accumulate mass terms
        auto mass_op = [](double acc, const auto& species) { return acc * species.a(); };
        const double massNumerator = std::accumulate(
            reaction.reactants().begin(),
            reaction.reactants().end(),
            1.0,
            mass_op
        );
        const double massDenominator = std::accumulate(
            reaction.products().begin(),
            reaction.products().end(),
            1.0,
            mass_op
        );

        // Accumulate partition functions
        auto pf_op = [&](double acc, const auto& species) { return acc * m_partitionFunction->evaluate(species.z(), species.a(), T9); };
        const double partitionFunctionNumerator = std::accumulate(
            reaction.reactants().begin(),
            reaction.reactants().end(),
            1.0,
            pf_op
        );
        const double partitionFunctionDenominator = std::accumulate(
            reaction.products().begin(),
            reaction.products().end(),
            1.0,
            pf_op
        );

        const double CT = std::pow(massNumerator/massDenominator, 1.5) * (partitionFunctionNumerator/partitionFunctionDenominator);

        return forwardRate * symmetryFactor * CT * expFactor;

    }

    double GraphEngine::calculateReverseRateTwoBodyDerivative(
        const reaction::Reaction &reaction,
        const double T9,
        const double reverseRate
    ) const {
        const double d_log_kFwd = reaction.calculate_forward_rate_log_derivative(T9);

        auto log_deriv_pf_op = [&](double acc, const auto& species) {
            const double g = m_partitionFunction->evaluate(species.z(), species.a(), T9);
            const double dg_dT = m_partitionFunction->evaluateDerivative(species.z(), species.a(), T9);
            return (g == 0.0) ? acc : acc + (dg_dT / g);
        };

        const double reactant_log_derivative_sum = std::accumulate(
            reaction.reactants().begin(),
            reaction.reactants().end(),
            0.0,
            log_deriv_pf_op
        );

        const double product_log_derivative_sum = std::accumulate(
            reaction.products().begin(),
            reaction.products().end(),
            0.0,
            log_deriv_pf_op
        );

        const double d_log_C = reactant_log_derivative_sum - product_log_derivative_sum;

        const double d_log_exp = reaction.qValue() / (m_constants.kB * T9 * T9);

        const double log_total_derivative = d_log_kFwd + d_log_C + d_log_exp;

        return reverseRate * log_total_derivative; // Return the derivative of the reverse rate with respect to T9

    }

    StepDerivatives<double> GraphEngine::calculateAllDerivativesUsingPrecomputation(
        const std::vector<double> &Y_in,
        const std::vector<double> &bare_rates,
        const double T9,
        const double rho
    ) const {
        // --- Calculate screening factors ---
        const std::vector<double> screeningFactors = m_screeningModel->calculateScreeningFactors(
            m_reactions,
            m_networkSpecies,
            Y_in,
            T9,
            rho
        );

        // --- Optimized loop ---
        std::vector<double> molarReactionFlows;
        molarReactionFlows.reserve(m_precomputedReactions.size());

        for (const auto& precomp : m_precomputedReactions) {
            double abundanceProduct = 1.0;
            bool below_threshold = false;
            for (size_t i = 0; i < precomp.unique_reactant_indices.size(); ++i) {
                const size_t reactantIndex = precomp.unique_reactant_indices[i];
                const int power = precomp.reactant_powers[i];
                const double abundance = Y_in[reactantIndex];
                if (abundance < MIN_ABUNDANCE_THRESHOLD) {
                    below_threshold = true;
                    break;
                }

                abundanceProduct *= std::pow(Y_in[reactantIndex], power);
            }
            if (below_threshold) {
                molarReactionFlows.push_back(0.0);
                continue; // Skip this reaction if any reactant is below the abundance threshold
            }

            const double bare_rate = bare_rates[precomp.reaction_index];
            const double screeningFactor = screeningFactors[precomp.reaction_index];
            const size_t numReactants = m_reactions[precomp.reaction_index].reactants().size();

            const double molarReactionFlow =
                    screeningFactor *
                    bare_rate *
                    precomp.symmetry_factor *
                    abundanceProduct *
                    std::pow(rho, numReactants);

            molarReactionFlows.push_back(molarReactionFlow);
        }

        // --- Assemble molar abundance derivatives ---
        StepDerivatives<double> result;
        result.dydt.assign(m_networkSpecies.size(), 0.0); // Initialize derivatives to zero
        for (size_t j = 0; j < m_precomputedReactions.size(); ++j) {
            const auto& precomp = m_precomputedReactions[j];
            const double R_j = molarReactionFlows[j];

            for (size_t i = 0; i < precomp.affected_species_indices.size(); ++i) {
                const size_t speciesIndex = precomp.affected_species_indices[i];
                const int stoichiometricCoefficient = precomp.stoichiometric_coefficients[i];

                // Update the derivative for this species
                result.dydt[speciesIndex] += static_cast<double>(stoichiometricCoefficient) * R_j / rho;
            }
        }

        // --- Calculate the nuclear energy generation rate ---
        double massProductionRate = 0.0; // [mol][s^-1]
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            const auto& species = m_networkSpecies[i];
            massProductionRate += result.dydt[i] * species.mass() * m_constants.u;
        }
        result.nuclearEnergyGenerationRate = -massProductionRate * m_constants.Na * m_constants.c * m_constants.c; // [erg][s^-1][g^-1]
        return result;

    }

    // --- Generate Stoichiometry Matrix ---
    void GraphEngine::generateStoichiometryMatrix() {
        LOG_TRACE_L1(m_logger, "Generating stoichiometry matrix...");

        // Task 1: Set dimensions and initialize the matrix
        size_t numSpecies = m_networkSpecies.size();
        size_t numReactions = m_reactions.size();
        m_stoichiometryMatrix.resize(numSpecies, numReactions, false);

        LOG_TRACE_L1(m_logger, "Stoichiometry matrix initialized with dimensions: {} rows (species) x {} columns (reactions).",
                 numSpecies, numReactions);

        // Task 2: Populate the stoichiometry matrix
        // Iterate through all reactions, assign them a column index, and fill in their stoichiometric coefficients.
        size_t reactionColumnIndex = 0;
        for (const auto& reaction : m_reactions) {
            // Get the net stoichiometry for the current reaction
            std::unordered_map<fourdst::atomic::Species, int> netStoichiometry = reaction.stoichiometry();

            // Iterate through the species and their coefficients in the stoichiometry map
            for (const auto& [species, coefficient] : netStoichiometry) {
                // Find the row index for this species
                auto it = m_speciesToIndexMap.find(species);
                if (it != m_speciesToIndexMap.end()) {
                    const size_t speciesRowIndex = it->second;
                    // Set the matrix element. Boost.uBLAS handles sparse insertion.
                    m_stoichiometryMatrix(speciesRowIndex, reactionColumnIndex) = coefficient;
                } else {
                    // This scenario should ideally not happen if m_networkSpeciesMap and m_speciesToIndexMap are correctly synced
                    LOG_ERROR(m_logger, "CRITICAL ERROR: Species '{}' from reaction '{}' stoichiometry not found in species to index map.",
                             species.name(), reaction.id());
                    m_logger -> flush_log();
                    throw std::runtime_error("Species not found in species to index map: " + std::string(species.name()));
                }
            }
            reactionColumnIndex++; // Move to the next column for the next reaction
        }

        LOG_TRACE_L1(m_logger, "Stoichiometry matrix population complete. Number of non-zero elements: {}.",
                 m_stoichiometryMatrix.nnz()); // Assuming nnz() exists for compressed_matrix
    }

    StepDerivatives<double> GraphEngine::calculateAllDerivatives(
        const std::vector<double> &Y_in,
        const double T9,
        const double rho
    ) const {
        return calculateAllDerivatives<double>(Y_in, T9, rho);
    }

    StepDerivatives<ADDouble> GraphEngine::calculateAllDerivatives(
        const std::vector<ADDouble> &Y_in,
        const ADDouble &T9,
        const ADDouble &rho
    ) const {
        return calculateAllDerivatives<ADDouble>(Y_in, T9, rho);
    }

    void GraphEngine::setScreeningModel(const screening::ScreeningType model) {
        m_screeningModel = screening::selectScreeningModel(model);
        m_screeningType = model;
    }

    screening::ScreeningType GraphEngine::getScreeningModel() const {
        return m_screeningType;
    }

    void GraphEngine::setPrecomputation(const bool precompute) {
        m_usePrecomputation = precompute;
    }

    bool GraphEngine::isPrecomputationEnabled() const {
        return m_usePrecomputation;
    }

    const partition::PartitionFunction & GraphEngine::getPartitionFunction() const {
        return *m_partitionFunction;
    }

    double GraphEngine::calculateMolarReactionFlow(
        const reaction::Reaction &reaction,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return calculateMolarReactionFlow<double>(reaction, Y, T9, rho);
    }

    void GraphEngine::generateJacobianMatrix(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) {

        LOG_TRACE_L1(m_logger, "Generating jacobian matrix for T9={}, rho={}..", T9, rho);
        const size_t numSpecies = m_networkSpecies.size();

        // 1. Pack the input variables into a vector for CppAD
        std::vector<double> adInput(numSpecies + 2, 0.0); // +2 for T9 and rho
        for (size_t i = 0; i < numSpecies; ++i) {
            adInput[i] = Y[i];
        }
        adInput[numSpecies]     = T9;  // T9
        adInput[numSpecies + 1] = rho; // rho

        // 2. Calculate the full jacobian
        const std::vector<double> dotY = m_rhsADFun.Jacobian(adInput);

        // 3. Pack jacobian vector into sparse matrix
        m_jacobianMatrix.clear();
        for (size_t i = 0; i < numSpecies; ++i) {
            for (size_t j = 0; j < numSpecies; ++j) {
                const double value = dotY[i * (numSpecies + 2) + j];
                if (std::abs(value) > MIN_JACOBIAN_THRESHOLD) {
                    m_jacobianMatrix(i, j) = value;
                }
            }
        }
        LOG_TRACE_L1(m_logger, "Jacobian matrix generated with dimensions: {} rows x {} columns.", m_jacobianMatrix.size1(), m_jacobianMatrix.size2());
    }

    double GraphEngine::getJacobianMatrixEntry(const int i, const int j) const {
        return m_jacobianMatrix(i, j);
    }

    std::unordered_map<fourdst::atomic::Species, int> GraphEngine::getNetReactionStoichiometry(
        const reaction::Reaction &reaction
    ) {
        return reaction.stoichiometry();
    }

    int GraphEngine::getStoichiometryMatrixEntry(
        const int speciesIndex,
        const int reactionIndex
    ) const {
        return m_stoichiometryMatrix(speciesIndex, reactionIndex);
    }

    void GraphEngine::exportToDot(const std::string &filename) const {
        LOG_TRACE_L1(m_logger, "Exporting network graph to DOT file: {}", filename);

        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file for writing: {}", filename);
            m_logger->flush_log();
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }

        dotFile << "digraph NuclearReactionNetwork {\n";
        dotFile << "    graph [rankdir=LR, splines=true, overlap=false, bgcolor=\"#f0f0f0\"];\n";
        dotFile << "    node [shape=circle, style=filled, fillcolor=\"#a7c7e7\", fontname=\"Helvetica\"];\n";
        dotFile << "    edge [fontname=\"Helvetica\", fontsize=10];\n\n";

        // 1. Define all species as nodes
        dotFile << "    // --- Species Nodes ---\n";
        for (const auto& species : m_networkSpecies) {
            dotFile << "    \"" << species.name() << "\" [label=\"" << species.name() << "\"];\n";
        }
        dotFile << "\n";

        // 2. Define all reactions as intermediate nodes and connect them
        dotFile << "    // --- Reaction Edges ---\n";
        for (const auto& reaction : m_reactions) {
            // Create a unique ID for the reaction node
            std::string reactionNodeId = "reaction_" + std::string(reaction.id());

            // Define the reaction node (small, black dot)
            dotFile << "    \"" << reactionNodeId << "\" [shape=point, fillcolor=black, width=0.1, height=0.1, label=\"\"];\n";

            // Draw edges from reactants to the reaction node
            for (const auto& reactant : reaction.reactants()) {
                dotFile << "    \"" << reactant.name() << "\" -> \"" << reactionNodeId << "\";\n";
            }

            // Draw edges from the reaction node to products
            for (const auto& product : reaction.products()) {
                dotFile << "    \"" << reactionNodeId << "\" -> \"" << product.name() << "\" [label=\"" << reaction.qValue() << " MeV\"];\n";
            }
            dotFile << "\n";
        }

        dotFile << "}\n";
        dotFile.close();
        LOG_TRACE_L1(m_logger, "Successfully exported network to {}", filename);
    }

    void GraphEngine::exportToCSV(const std::string &filename) const {
        LOG_TRACE_L1(m_logger, "Exporting network graph to CSV file: {}", filename);

        std::ofstream csvFile(filename, std::ios::out | std::ios::trunc);
        if (!csvFile.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file for writing: {}", filename);
            m_logger->flush_log();
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        csvFile << "Reaction;Reactants;Products;Q-value;sources;rates\n";
        for (const auto& reaction : m_reactions) {
            // Dynamic cast to REACLIBReaction to access specific properties
            csvFile << reaction.id() << ";";
            // Reactants
            int count = 0;
            for (const auto& reactant : reaction.reactants()) {
                csvFile << reactant.name();
                if (++count < reaction.reactants().size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";";
            count = 0;
            for (const auto& product : reaction.products()) {
                csvFile << product.name();
                if (++count < reaction.products().size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";" << reaction.qValue() << ";";
            // Reaction coefficients
            auto sources = reaction.sources();
            count = 0;
            for (const auto& source : sources) {
                csvFile << source;
                if (++count < sources.size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";";
            // Reaction coefficients
            count = 0;
            for (const auto& rates : reaction) {
                csvFile << rates;
                if (++count < reaction.size()) {
                    csvFile << ",";
                }
            }
            csvFile << "\n";
        }
        csvFile.close();
        LOG_TRACE_L1(m_logger, "Successfully exported network graph to {}", filename);
    }

    std::unordered_map<fourdst::atomic::Species, double> GraphEngine::getSpeciesTimescales(const std::vector<double> &Y, const double T9,
        const double rho) const {
        auto [dydt, _] = calculateAllDerivatives<double>(Y, T9, rho);
        std::unordered_map<fourdst::atomic::Species, double> speciesTimescales;
        speciesTimescales.reserve(m_networkSpecies.size());
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            double timescale = std::numeric_limits<double>::infinity();
            const auto species = m_networkSpecies[i];
            if (std::abs(dydt[i]) > 0.0) {
                timescale = std::abs(Y[i] / dydt[i]);
            }
            speciesTimescales.emplace(species, timescale);
        }
        return speciesTimescales;
    }

    void GraphEngine::update(const NetIn &netIn) {
        // No-op for GraphEngine, as it does not support manually triggering updates.
    }

    void GraphEngine::recordADTape() {
        LOG_TRACE_L1(m_logger, "Recording AD tape for the RHS calculation...");

        // Task 1: Set dimensions and initialize the matrix
        const size_t numSpecies = m_networkSpecies.size();
        if (numSpecies == 0) {
            LOG_ERROR(m_logger, "Cannot record AD tape: No species in the network.");
            m_logger->flush_log();
            throw std::runtime_error("Cannot record AD tape: No species in the network.");
        }
        const size_t numADInputs = numSpecies + 2; // Note here that by not letting T9 and rho be independent variables, we are constraining the network to a constant temperature and density during each evaluation.

        // --- CppAD Tape Recording ---
        // 1. Declare independent variable (adY)
        //    We also initialize the dummy variable for tape recording (these tell CppAD what the derivative chain looks like).
        //    Their numeric values are irrelevant except for in so far as they avoid numerical instabilities.

        // Distribute total mass fraction uniformly between species in the dummy variable space
        const auto uniformMassFraction = static_cast<CppAD::AD<double>>(1.0 / static_cast<double>(numSpecies));
        std::vector<CppAD::AD<double>> adInput(numADInputs, uniformMassFraction);
        adInput[numSpecies]     = 1.0; // Dummy T9
        adInput[numSpecies + 1] = 1.0; // Dummy rho

        // 3. Declare independent variables (what CppAD will differentiate wrt.)
        //    This also beings the tape recording process.
        CppAD::Independent(adInput);

        std::vector<CppAD::AD<double>> adY(numSpecies);
        for(size_t i = 0; i < numSpecies; ++i) {
            adY[i] = adInput[i];
        }
        const CppAD::AD<double> adT9  = adInput[numSpecies];
        const CppAD::AD<double> adRho = adInput[numSpecies + 1];


        // 5. Call the actual templated function
        // We let T9 and rho be constant, so we pass them as fixed values.
        auto [dydt, nuclearEnergyGenerationRate] = calculateAllDerivatives<CppAD::AD<double>>(adY, adT9, adRho);

        m_rhsADFun.Dependent(adInput, dydt);

        LOG_TRACE_L1(m_logger, "AD tape recorded successfully for the RHS calculation. Number of independent variables: {}.",
                 adInput.size());
    }

    void GraphEngine::precomputeNetwork() {
        LOG_TRACE_L1(m_logger, "Pre-computing constant components of GraphNetwork state...");

        // --- Reverse map for fast species lookups ---
        std::unordered_map<fourdst::atomic::Species, size_t> speciesIndexMap;
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            speciesIndexMap[m_networkSpecies[i]] = i;
        }

        m_precomputedReactions.clear();
        m_precomputedReactions.reserve(m_reactions.size());

        for (size_t i = 0; i < m_reactions.size(); ++i) {
            const auto& reaction = m_reactions[i];
            PrecomputedReaction precomp;
            precomp.reaction_index = i;

            // --- Precompute reactant information ---
            // Count occurrences for each reactant to determine powers and symmetry
            std::unordered_map<size_t, int> reactantCounts;
            for (const auto& reactant: reaction.reactants()) {
                size_t reactantIndex = speciesIndexMap.at(reactant);
                reactantCounts[reactantIndex]++;
            }

            double symmetryDenominator = 1.0;
            for (const auto& [index, count] : reactantCounts) {
                precomp.unique_reactant_indices.push_back(index);
                precomp.reactant_powers.push_back(count);

                symmetryDenominator *= 1.0/std::tgamma(count + 1);
            }

            precomp.symmetry_factor = symmetryDenominator;

            // --- Precompute stoichiometry information ---
            const auto stoichiometryMap = reaction.stoichiometry();
            precomp.affected_species_indices.reserve(stoichiometryMap.size());
            precomp.stoichiometric_coefficients.reserve(stoichiometryMap.size());

            for (const auto& [species, coeff] : stoichiometryMap) {
                precomp.affected_species_indices.push_back(speciesIndexMap.at(species));
                precomp.stoichiometric_coefficients.push_back(coeff);
            }

            m_precomputedReactions.push_back(std::move(precomp));
        }
    }

    bool AtomicReverseRate::forward(
        const size_t p,
        const size_t q,
        const CppAD::vector<bool> &vx,
        CppAD::vector<bool> &vy,
        const CppAD::vector<double> &tx,
        CppAD::vector<double> &ty
    ) {
        const double T9 = tx[0];
        const double expFactor = std::exp(-m_reaction.qValue() / (m_kB * T9 * 1e9)); // Convert MeV to erg
        if (p == 0) {
            // --- Zeroth order forward sweep ---
            const auto k_rev = m_engine.calculateReverseRate(m_reaction, T9, expFactor);
        }
    }

}
