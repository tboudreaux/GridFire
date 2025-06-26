#include "gridfire/engine/engine_graph.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/network.h"

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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/odeint.hpp>


namespace gridfire {
    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition
    ):
     m_reactions(build_reaclib_nuclear_network(composition, false)) {
        syncInternalMaps();
    }

    GraphEngine::GraphEngine(reaction::REACLIBLogicalReactionSet reactions) :
        m_reactions(std::move(reactions)) {
            syncInternalMaps();
        }

    StepDerivatives<double> GraphEngine::calculateRHSAndEnergy(
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return calculateAllDerivatives<double>(Y, T9, rho);
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
            for (const auto& reactant: reaction->reactants()) {
                uniqueSpeciesNames.insert(reactant.name());
            }
            for (const auto& product: reaction->products()) {
                uniqueSpeciesNames.insert(product.name());
            }
        }

        for (const auto& name: uniqueSpeciesNames) {
            auto it = fourdst::atomic::species.find(name);
            if (it != fourdst::atomic::species.end()) {
                m_networkSpecies.push_back(it->second);
                m_networkSpeciesMap.insert({name, it->second});
            } else {
                LOG_ERROR(m_logger, "Species '{}' not found in global atomic species database.", name);
                throw std::runtime_error("Species not found in global atomic species database: " + std::string(name));
            }
        }

    }

    void GraphEngine::populateReactionIDMap() {
        LOG_TRACE_L1(m_logger, "Populating reaction ID map for REACLIB graph network (serif::network::GraphNetwork)...");
        m_reactionIDMap.clear();
        for (const auto& reaction: m_reactions) {
            m_reactionIDMap.emplace(reaction->id(), reaction.get());
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
        LOG_DEBUG(m_logger, "Providing access to network species vector. Size: {}.", m_networkSpecies.size());
        return m_networkSpecies;
    }

    const reaction::REACLIBLogicalReactionSet& GraphEngine::getNetworkReactions() const {
        // Returns a constant reference to the set of reactions in the network.
        LOG_DEBUG(m_logger, "Providing access to network reactions set. Size: {}.", m_reactions.size());
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
            for (const auto& reactant : reaction->reactants()) {
                auto it = m_networkSpeciesMap.find(reactant.name());
                if (it != m_networkSpeciesMap.end()) {
                    totalReactantA += it->second.a();
                    totalReactantZ += it->second.z();
                } else {
                    // This scenario indicates a severe data integrity issue:
                    // a reactant is part of a reaction but not in the network's species map.
                    LOG_ERROR(m_logger, "CRITICAL ERROR: Reactant species '{}' in reaction '{}' not found in network species map during conservation validation.",
                             reactant.name(), reaction->id());
                    return false;
                }
            }

            // Calculate total A and Z for products
            for (const auto& product : reaction->products()) {
                auto it = m_networkSpeciesMap.find(product.name());
                if (it != m_networkSpeciesMap.end()) {
                    totalProductA += it->second.a();
                    totalProductZ += it->second.z();
                } else {
                    // Similar critical error for product species
                    LOG_ERROR(m_logger, "CRITICAL ERROR: Product species '{}' in reaction '{}' not found in network species map during conservation validation.",
                             product.name(), reaction->id());
                    return false;
                }
            }

            // Compare totals for conservation
            if (totalReactantA != totalProductA) {
                LOG_ERROR(m_logger, "Mass number (A) not conserved for reaction '{}': Reactants A={} vs Products A={}.",
                         reaction->id(), totalReactantA, totalProductA);
                return false;
            }
            if (totalReactantZ != totalProductZ) {
                LOG_ERROR(m_logger, "Atomic number (Z) not conserved for reaction '{}': Reactants Z={} vs Products Z={}.",
                         reaction->id(), totalReactantZ, totalProductZ);
                return false;
            }
        }

        LOG_TRACE_L1(m_logger, "Mass (A) and charge (Z) conservation validated successfully for all reactions.");
        return true; // All reactions passed the conservation check
    }

    void GraphEngine::validateComposition(const fourdst::composition::Composition &composition, double culling, double T9) {
        // Check if the requested network has already been cached.
        // PERF: Rebuilding this should be pretty fast but it may be a good point of optimization in the future.
        const reaction::REACLIBLogicalReactionSet validationReactionSet = build_reaclib_nuclear_network(composition, false);
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
            std::unordered_map<fourdst::atomic::Species, int> netStoichiometry = reaction->stoichiometry();

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
                             species.name(), reaction->id());
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
        const ADDouble T9,
        const ADDouble rho
    ) const {
        return calculateAllDerivatives<ADDouble>(Y_in, T9, rho);
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
        LOG_DEBUG(m_logger, "Jacobian matrix generated with dimensions: {} rows x {} columns.", m_jacobianMatrix.size1(), m_jacobianMatrix.size2());
    }

    double GraphEngine::getJacobianMatrixEntry(const int i, const int j) const {
        return m_jacobianMatrix(i, j);
    }

    std::unordered_map<fourdst::atomic::Species, int> GraphEngine::getNetReactionStoichiometry(
        const reaction::Reaction &reaction
    ) const {
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
            std::string reactionNodeId = "reaction_" + std::string(reaction->id());

            // Define the reaction node (small, black dot)
            dotFile << "    \"" << reactionNodeId << "\" [shape=point, fillcolor=black, width=0.1, height=0.1, label=\"\"];\n";

            // Draw edges from reactants to the reaction node
            for (const auto& reactant : reaction->reactants()) {
                dotFile << "    \"" << reactant.name() << "\" -> \"" << reactionNodeId << "\";\n";
            }

            // Draw edges from the reaction node to products
            for (const auto& product : reaction->products()) {
                dotFile << "    \"" << reactionNodeId << "\" -> \"" << product.name() << "\" [label=\"" << reaction->qValue() << " MeV\"];\n";
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
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        csvFile << "Reaction;Reactants;Products;Q-value;sources;rates\n";
        for (const auto& reaction : m_reactions) {
            // Dynamic cast to REACLIBReaction to access specific properties
            csvFile << reaction->id() << ";";
            // Reactants
            int count = 0;
            for (const auto& reactant : reaction->reactants()) {
                csvFile << reactant.name();
                if (++count < reaction->reactants().size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";";
            count = 0;
            for (const auto& product : reaction->products()) {
                csvFile << product.name();
                if (++count < reaction->products().size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";" << reaction->qValue() << ";";
            // Reaction coefficients
            auto* reaclibReaction = dynamic_cast<const reaction::REACLIBLogicalReaction*>(reaction.get());
            if (!reaclibReaction) {
                LOG_ERROR(m_logger, "Failed to cast Reaction to REACLIBLogicalReaction in GraphNetwork::exportToCSV().");
                throw std::runtime_error("Failed to cast Reaction to REACLIBLogicalReaction in GraphNetwork::exportToCSV(). This should not happen, please check your reaction setup.");
            }
            auto sources = reaclibReaction->sources();
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
            for (const auto& rates : *reaclibReaction) {
                csvFile << rates;
                if (++count < reaclibReaction->size()) {
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

    void GraphEngine::recordADTape() {
        LOG_TRACE_L1(m_logger, "Recording AD tape for the RHS calculation...");

        // Task 1: Set dimensions and initialize the matrix
        const size_t numSpecies = m_networkSpecies.size();
        if (numSpecies == 0) {
            LOG_ERROR(m_logger, "Cannot record AD tape: No species in the network.");
            throw std::runtime_error("Cannot record AD tape: No species in the network.");
        }
        const size_t numADInputs = numSpecies + 2; // Note here that by not letting T9 and rho be independent variables, we are constraining the network to a constant temperature and density during each evaluation.

        // --- CppAD Tape Recording ---
        // 1. Declare independent variable (adY)
        //    We also initialize the dummy variable for tape recording (these tell CppAD what the derivative chain looks like).
        //    Their numeric values are irrelevant except for in so far as they avoid numerical instabilities.

        // Distribute total mass fraction uniformly between species in the dummy variable space
        const auto uniformMassFraction = static_cast<CppAD::AD<double>>(1.0 / numSpecies);
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
}
