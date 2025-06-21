#include "netgraph.h"
#include "atomicSpecies.h"
#include "const.h"
#include "network.h"
#include "reaclib.h"
#include "species.h"

#include "quill/LogMacros.h"

#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/odeint.hpp>


namespace gridfire {
    GraphNetwork::GraphNetwork(
        const fourdst::composition::Composition &composition
    ):
     Network(REACLIB),
     m_reactions(build_reaclib_nuclear_network(composition)) {
        syncInternalMaps();
    }

    GraphNetwork::GraphNetwork(
        const fourdst::composition::Composition &composition,
        const double cullingThreshold,
        const double T9
    ):
     Network(REACLIB),
     m_reactions(build_reaclib_nuclear_network(composition, cullingThreshold, T9)) {
        syncInternalMaps();
    }

    void GraphNetwork::syncInternalMaps() {
        collectNetworkSpecies();
        populateReactionIDMap();
        populateSpeciesToIndexMap();
        reserveJacobianMatrix();
        recordADTape();
    }

    // --- Network Graph Construction Methods ---
    void GraphNetwork::collectNetworkSpecies() {
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

    void GraphNetwork::populateReactionIDMap() {
        LOG_INFO(m_logger, "Populating reaction ID map for REACLIB graph network (serif::network::GraphNetwork)...");
        m_reactionIDMap.clear();
        for (const auto& reaction: m_reactions) {
            m_reactionIDMap.insert({reaction.id(), reaction});
        }
        LOG_INFO(m_logger, "Populated {} reactions in the reaction ID map.", m_reactionIDMap.size());
    }

    void GraphNetwork::populateSpeciesToIndexMap() {
        m_speciesToIndexMap.clear();
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            m_speciesToIndexMap.insert({m_networkSpecies[i], i});
        }
    }

    void GraphNetwork::reserveJacobianMatrix() {
        // The implementation of this function (and others) constrains this nuclear network to a constant temperature and density during
        // each evaluation.
        size_t numSpecies = m_networkSpecies.size();
        m_jacobianMatrix.clear();
        m_jacobianMatrix.resize(numSpecies, numSpecies, false); // Sparse matrix, no initial values
        LOG_INFO(m_logger, "Jacobian matrix resized to {} rows and {} columns.",
                 m_jacobianMatrix.size1(), m_jacobianMatrix.size2());
    }

    // --- Basic Accessors and Queries ---
    const std::vector<fourdst::atomic::Species>& GraphNetwork::getNetworkSpecies() const {
        // Returns a constant reference to the vector of unique species in the network.
        LOG_DEBUG(m_logger, "Providing access to network species vector. Size: {}.", m_networkSpecies.size());
        return m_networkSpecies;
    }

    const reaclib::REACLIBReactionSet& GraphNetwork::getNetworkReactions() const {
        // Returns a constant reference to the set of reactions in the network.
        LOG_DEBUG(m_logger, "Providing access to network reactions set. Size: {}.", m_reactions.size());
        return m_reactions;
    }

    bool GraphNetwork::involvesSpecies(const fourdst::atomic::Species& species) const {
        // Checks if a given species is present in the network's species map for efficient lookup.
        const bool found = m_networkSpeciesMap.contains(species.name());
        LOG_DEBUG(m_logger, "Checking if species '{}' is involved in the network: {}.", species.name(), found ? "Yes" : "No");
        return found;
    }

    std::unordered_map<fourdst::atomic::Species, int> GraphNetwork::getNetReactionStoichiometry(const reaclib::REACLIBReaction& reaction) const {
        // Calculates the net stoichiometric coefficients for species in a given reaction.
        std::unordered_map<fourdst::atomic::Species, int> stoichiometry;

        // Iterate through reactants, decrementing their counts
        for (const auto& reactant : reaction.reactants()) {
            auto it = m_networkSpeciesMap.find(reactant.name());
            if (it != m_networkSpeciesMap.end()) {
                stoichiometry[it->second]--; // Copy Species by value (PERF: Future performance improvements by using pointers or references (std::reference_wrapper<const ...>) or something like that)
            } else {
                LOG_WARNING(m_logger, "Reactant species '{}' in reaction '{}' not found in network species map during stoichiometry calculation.",
                         reactant.name(), reaction.id());
            }
        }

        // Iterate through products, incrementing their counts
        for (const auto& product : reaction.products()) {
            auto it = m_networkSpeciesMap.find(product.name());
            if (it != m_networkSpeciesMap.end()) {
                stoichiometry[it->second]++; // Copy Species by value (PERF: Future performance improvements by using pointers or references (std::reference_wrapper<const ...>) or something like that)
            } else {
                LOG_WARNING(m_logger, "Product species '{}' in reaction '{}' not found in network species map during stoichiometry calculation.",
                         product.name(), reaction.id());
            }
        }
        LOG_DEBUG(m_logger, "Calculated net stoichiometry for reaction '{}'. Total unique species in stoichiometry: {}.", reaction.id(), stoichiometry.size());
        return stoichiometry;
    }

    // --- Validation Methods ---
    bool GraphNetwork::validateConservation() const {
        LOG_INFO(m_logger, "Validating mass (A) and charge (Z) conservation across all reactions in the network.");

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

        LOG_INFO(m_logger, "Mass (A) and charge (Z) conservation validated successfully for all reactions.");
        return true; // All reactions passed the conservation check
    }

    void GraphNetwork::validateComposition(const fourdst::composition::Composition &composition, double culling, double T9) {

        // Check if the requested network has already been cached.
        // PERF: Rebuilding this should be pretty fast but it may be a good point of optimization in the future.
        const reaclib::REACLIBReactionSet validationReactionSet = build_reaclib_nuclear_network(composition, culling, T9);
        // TODO: need some more robust method here to
        //       A. Build the basic network from the composition's species with non zero mass fractions.
        //       B. rebuild a new composition from the reaction set's reactants + products (with the mass fractions from the things that are only products set to 0)
        //       C. Rebuild the reaction set from the new composition
        //       D. Cull reactions where all reactants have mass fractions below the culling threshold.
        //       E. Be careful about maintaining caching through all of this


        // This allows for dynamic network modification while retaining caching for networks which are very similar.
        if (validationReactionSet != m_reactions) {
            LOG_INFO(m_logger, "Reaction set not cached. Rebuilding the reaction set for T9={} and culling={}.", T9, culling);
            m_reactions = validationReactionSet;
            syncInternalMaps(); // Re-sync internal maps after updating reactions. Note this will also retrace the AD tape.
        }
    }

    // --- Generate Stoichiometry Matrix ---
    void GraphNetwork::generateStoichiometryMatrix() {
        LOG_INFO(m_logger, "Generating stoichiometry matrix...");

        // Task 1: Set dimensions and initialize the matrix
        size_t numSpecies = m_networkSpecies.size();
        size_t numReactions = m_reactions.size();
        m_stoichiometryMatrix.resize(numSpecies, numReactions, false);

        LOG_INFO(m_logger, "Stoichiometry matrix initialized with dimensions: {} rows (species) x {} columns (reactions).",
                 numSpecies, numReactions);

        // Task 2: Populate the stoichiometry matrix
        // Iterate through all reactions, assign them a column index, and fill in their stoichiometric coefficients.
        size_t reactionColumnIndex = 0;
        for (const auto& reaction : m_reactions) {
            // Get the net stoichiometry for the current reaction
            std::unordered_map<fourdst::atomic::Species, int> netStoichiometry = getNetReactionStoichiometry(reaction);

            // Iterate through the species and their coefficients in the stoichiometry map
            for (const auto& pair : netStoichiometry) {
                const fourdst::atomic::Species& species = pair.first; // The Species object
                const int coefficient = pair.second;                // The stoichiometric coefficient

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
                    throw std::runtime_error("Species not found in species to index map: " + std::string(species.name()));
                }
            }
            reactionColumnIndex++; // Move to the next column for the next reaction
        }

        LOG_INFO(m_logger, "Stoichiometry matrix population complete. Number of non-zero elements: {}.",
                 m_stoichiometryMatrix.nnz()); // Assuming nnz() exists for compressed_matrix
    }

    void GraphNetwork::generateJacobianMatrix(const std::vector<double> &Y, const double T9,
        const double rho) {

        LOG_INFO(m_logger, "Generating jacobian matrix for T9={}, rho={}..", T9, rho);
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
        LOG_INFO(m_logger, "Jacobian matrix generated with dimensions: {} rows x {} columns.", m_jacobianMatrix.size1(), m_jacobianMatrix.size2());
    }

    void GraphNetwork::detectStiff(const NetIn &netIn, const double T9, const double numSpecies, const boost::numeric::ublas::vector<double>& Y) {
        // --- Heuristic for automatic stiffness detection ---
        const std::vector<double> initial_y_stl(Y.begin(), Y.begin() + numSpecies); // Copy only the species abundances, not the specific energy rate
        const auto derivatives = calculateAllDerivatives<double>(initial_y_stl, T9, netIn.density);
        const std::vector<double>& initial_dotY = derivatives.dydt;

        double min_timescale = std::numeric_limits<double>::max();
        double max_timescale = 0.0;

        for (size_t i = 0; i < numSpecies; ++i) {
            if (std::abs(initial_dotY[i]) > 0) {
                const double timescale = std::abs(Y(i) / initial_dotY[i]);
                if (timescale > max_timescale) {max_timescale = timescale;}
                if (timescale < min_timescale) {min_timescale = timescale;}
            }
        }

        const double stiffnessRatio = max_timescale / min_timescale;

        // TODO: Pull this out into a configuration option or a more sophisticated heuristic.
        constexpr double stiffnessThreshold = 1.0e6; // This is a heuristic threshold, can be tuned based on network characteristics.

        LOG_INFO(m_logger, "Stiffness ratio is {} (max timescale: {}, min timescale: {}).", stiffnessRatio, max_timescale, min_timescale);
        if (stiffnessRatio > stiffnessThreshold) {
            LOG_INFO(m_logger, "Network is detected to be stiff. Using stiff ODE solver.");
            m_stiff = true;
        } else {
            LOG_INFO(m_logger, "Network is detected to be non-stiff. Using non-stiff ODE solver.");
            m_stiff = false;
        }
    }

    void GraphNetwork::exportToDot(const std::string &filename) const {
        LOG_INFO(m_logger, "Exporting network graph to DOT file: {}", filename);

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
        LOG_INFO(m_logger, "Successfully exported network to {}", filename);
    }

    NetOut GraphNetwork::evaluate(const NetIn &netIn) {
        namespace ublas = boost::numeric::ublas;
        namespace odeint = boost::numeric::odeint;

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        validateComposition(netIn.composition, netIn.culling, T9);

        const double numSpecies = m_networkSpecies.size();
        constexpr double abs_tol = 1.0e-8;
        constexpr double rel_tol = 1.0e-8;

        int stepCount = 0;

        // TODO: Pull these out into configuration options

        ODETerm rhs_functor(*this, T9, netIn.density);


        ublas::vector<double> Y(numSpecies + 1);
        for (size_t i = 0; i < numSpecies; ++i) {
            const auto& species = m_networkSpecies[i];
            // Get the mass fraction for this specific species from the input object
            Y(i) = netIn.composition.getMassFraction(std::string(species.name()));
        }
        Y(numSpecies) = 0; // initial specific energy rate, will be updated later

        detectStiff(netIn, T9, numSpecies, Y);

        if (m_stiff) {
            JacobianTerm jacobian_functor(*this, T9, netIn.density);
            LOG_INFO(m_logger, "Making use of stiff ODE solver for network evaluation.");
            auto stepper = odeint::make_controlled<odeint::rosenbrock4<double>>(abs_tol, rel_tol);
            stepCount = odeint::integrate_adaptive(
                stepper,
                std::make_pair(rhs_functor, jacobian_functor),
                Y,
                0.0, // Start time
                netIn.tMax,
                netIn.dt0
                );

        } else {
            LOG_INFO(m_logger, "Making use of ODE solver (non-stiff) for network evaluation.");
            using state_type = ublas::vector<double>;
            auto stepper = odeint::make_controlled<odeint::runge_kutta_dopri5<state_type>>(abs_tol, rel_tol);
            stepCount = odeint::integrate_adaptive(
                stepper,
                rhs_functor,
                Y,
                0.0, // Start time
                netIn.tMax,
                netIn.dt0
            );


        }

        double sumY = 0.0;
        for (int i = 0; i < numSpecies; ++i) { sumY += Y(i); }
        for (int i = 0; i < numSpecies; ++i) { Y(i) /= sumY; }

        // --- Marshall output variables ---
        // PERF: Im sure this step could be tuned to avoid so many copies, that is a job for another day
        std::vector<std::string> speciesNames;
        speciesNames.reserve(numSpecies);
        for (const auto& species : m_networkSpecies) {
            speciesNames.push_back(std::string(species.name()));
        }

        std::vector<double> finalAbundances(Y.begin(), Y.begin() + numSpecies);
        fourdst::composition::Composition outputComposition(speciesNames, finalAbundances);
        outputComposition.finalize(true);

        NetOut netOut;
        netOut.composition = outputComposition;
        netOut.num_steps = stepCount;
        netOut.energy = Y(numSpecies); // The last element in Y is the specific energy rate

        return netOut;

    }

    void GraphNetwork::recordADTape() {
        LOG_INFO(m_logger, "Recording AD tape for the RHS calculation...");

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
        auto derivatives = calculateAllDerivatives<CppAD::AD<double>>(adY, adT9, adRho);

        m_rhsADFun.Dependent(adInput, derivatives.dydt);

        LOG_INFO(m_logger, "AD tape recorded successfully for the RHS calculation. Number of independent variables: {}.",
                 adInput.size());
    }
}
