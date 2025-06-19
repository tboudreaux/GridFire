#include "netgraph.h"
#include "atomicSpecies.h"
#include "reaclib.h"
#include "network.h"
#include "species.h"
#include "const.h"

#include "quill/LogMacros.h"

#include <set>
#include <unordered_map>
#include <string>
#include <vector>
#include <stdexcept>
#include <string_view>
#include <cstdint>
#include <iostream>

#include "const.h"


namespace serif::network {
    GraphNetwork::GraphNetwork(
        const serif::composition::Composition &composition
    ):
     Network(REACLIB),
     m_reactions(build_reaclib_nuclear_network(composition)),
     m_initialComposition(composition),
     m_currentComposition(composition),
     m_cullingThreshold(0),
     m_T9(0) {
        syncInternalMaps();
    }

    GraphNetwork::GraphNetwork(
        const serif::composition::Composition &composition,
        const double cullingThreshold,
        const double T9
    ):
     Network(REACLIB),
     m_reactions(build_reaclib_nuclear_network(composition, cullingThreshold, T9)),
     m_initialComposition(composition),
     m_currentComposition(composition),
     m_cullingThreshold(cullingThreshold),
     m_T9(T9) {
        syncInternalMaps();
    }

    void GraphNetwork::syncInternalMaps() {
        collectNetworkSpecies();
        populateReactionIDMap();
        populateSpeciesToIndexMap();
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
            auto it = serif::atomic::species.find(name);
            if (it != serif::atomic::species.end()) {
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

    void GraphNetwork::buildNetworkGraph() {
        LOG_INFO(m_logger, "Building network graph...");
        m_reactions = build_reaclib_nuclear_network(m_initialComposition, m_cullingThreshold, m_T9);
        syncInternalMaps();
        LOG_INFO(m_logger, "Network graph built with {} reactions and {} species.", m_reactions.size(), m_networkSpecies.size());
    }
        // --- Basic Accessors and Queries ---

    const std::vector<serif::atomic::Species>& GraphNetwork::getNetworkSpecies() const {
        // Returns a constant reference to the vector of unique species in the network.
        LOG_DEBUG(m_logger, "Providing access to network species vector. Size: {}.", m_networkSpecies.size());
        return m_networkSpecies;
    }

    const reaclib::REACLIBReactionSet& GraphNetwork::getNetworkReactions() const {
        // Returns a constant reference to the set of reactions in the network.
        LOG_DEBUG(m_logger, "Providing access to network reactions set. Size: {}.", m_reactions.size());
        return m_reactions;
    }

    bool GraphNetwork::involvesSpecies(const serif::atomic::Species& species) const {
        // Checks if a given species is present in the network's species map for efficient lookup.
        const bool found = m_networkSpeciesMap.contains(species.name());
        LOG_DEBUG(m_logger, "Checking if species '{}' is involved in the network: {}.", species.name(), found ? "Yes" : "No");
        return found;
    }

    std::unordered_map<serif::atomic::Species, int> GraphNetwork::getNetReactionStoichiometry(const reaclib::REACLIBReaction& reaction) const {
        // Calculates the net stoichiometric coefficients for species in a given reaction.
        std::unordered_map<serif::atomic::Species, int> stoichiometry;

        // Iterate through reactants, decrementing their counts
        for (const auto& reactant : reaction.reactants()) {
            auto it = m_networkSpeciesMap.find(reactant.name());
            if (it != m_networkSpeciesMap.end()) {
                stoichiometry[it->second]--; // Copy Species by value (PERF: Future performance improvements by using pointers or references)
            } else {
                LOG_WARNING(m_logger, "Reactant species '{}' in reaction '{}' not found in network species map during stoichiometry calculation.",
                         reactant.name(), reaction.id());
            }
        }

        // Iterate through products, incrementing their counts
        for (const auto& product : reaction.products()) {
            auto it = m_networkSpeciesMap.find(product.name());
            if (it != m_networkSpeciesMap.end()) {
                stoichiometry[it->second]++; // Copy Species by value (PERF: Future performance improvements by using pointers or references)
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
            std::unordered_map<serif::atomic::Species, int> netStoichiometry = getNetReactionStoichiometry(reaction);

            // Iterate through the species and their coefficients in the stoichiometry map
            for (const auto& pair : netStoichiometry) {
                const serif::atomic::Species& species = pair.first; // The Species object
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

    void GraphNetwork::generateJacobianMatrix() {
        LOG_INFO(m_logger, "Generating jacobian matrix...");

        // Task 1: Set dimensions and initialize the matrix
        size_t numSpecies = m_networkSpecies.size();
        m_stoichiometryMatrix.resize(numSpecies, numSpecies, false);

        LOG_INFO(m_logger, "Jacobian matrix initialized with dimensions: {} rows (species) x {} columns (species).",
                 numSpecies, numSpecies);
    }

    std::vector<double> GraphNetwork::calculateRHS(const std::vector<double>& Y, const double T9, const double rho) const {
        std::vector<double> dotY(m_networkSpecies.size(), 0.0);
        const size_t numReactions = m_reactions.size();
        const size_t numSpecies = m_networkSpecies.size();

        if (m_stoichiometryMatrix.size1() != numSpecies || m_stoichiometryMatrix.size2() != numReactions) {
            LOG_ERROR(m_logger, "Stoichiometry matrix dimensions do not match network species and reactions sizes.");
            throw std::runtime_error("Stoichiometry matrix dimensions mismatch.");
        }

        for (size_t reactionIndex = 0; reactionIndex < numReactions; ++reactionIndex) {
            const auto& currentReaction = m_reactions[reactionIndex];
            const double reactionRatePerVolumePerTime = calculateReactionRate(currentReaction, Y, T9, rho);

            // dY_i/dt = ∑ υ_ij * R_j * mass_in_grams
            // TODO: make sure the unit conversion is correct after calculate reaction rate.

            for (size_t speciesIndex = 0; speciesIndex < numSpecies; ++speciesIndex) {
                const int nu_ij = m_stoichiometryMatrix(speciesIndex, reactionIndex);
                if (nu_ij != 0) {
                    const double speciesAtomicMassAMU = m_networkSpecies[speciesIndex].mass();
                    dotY[speciesIndex] += (nu_ij * reactionRatePerVolumePerTime * speciesAtomicMassAMU)/rho;
                }
            }
        }
        return dotY;
    }

    double GraphNetwork::calculateReactionRate(const reaclib::REACLIBReaction &reaction, const std::vector<double> &Y,
        const double T9, const double rho) const {
        const auto &constants  = serif::constant::Constants::getInstance();

        const auto u = constants.get("u"); // Atomic mass unit in g/mol
        const double k_reaction = reaction.calculate_rate(T9); // PERF: Consider precomputing all of these and putting them into an O(1) lookup table.

        double reactant_product = 1.0;

        std::unordered_map<std::string, int> reactant_counts;
        reactant_counts.reserve(reaction.reactants().size());
        for (const auto& reactant : reaction.reactants()) {
            reactant_counts[std::string(reactant.name())]++;
        }

        for (const auto& [species_name, count] : reactant_counts) {
            constexpr double minThreshold = 1e-18;
            auto species_it = m_speciesToIndexMap.find(m_networkSpeciesMap.at(species_name));
            if (species_it == m_speciesToIndexMap.end()) {
                LOG_ERROR(m_logger, "Reactant species '{}' not found in species to index map for reaction '{}'.",
                         species_name, reaction.id());
                throw std::runtime_error("Reactant species not found in species to index map: " + species_name);
            }
            const size_t species_index = species_it->second;
            const double Yi = Y[species_index];
            double Ai = m_networkSpecies[species_index].a();

            if (Yi < minThreshold) {
                return 0.0; // If any reactant is below a threshold, return zero rate.
            }

            double atomicMassAMU = m_networkSpecies[species_index].mass();

            // Convert to number density
            double ni;
            const double denominator = atomicMassAMU * u.value;
            if (denominator > minThreshold) {
                ni = (Yi * rho) / (denominator);
            } else {
                ni = 0.0;
            }

            reactant_product *= ni;

            // Apply factorial correction for identical reactions
            if (count > 1) {
                reactant_product /= static_cast<double>(std::tgamma(count + 1)); // Gamma function for factorial
            }
        }
        const double Na = constants.get("N_a").value; // Avogadro's number in mol^-1
        const double molarCorrectionFactor = std::pow(Na, reaction.reactants().size() - 1);
        return (reactant_product * k_reaction) / molarCorrectionFactor; // reaction rate in per volume per time (particles/cm^3/s)
    }

    NetOut GraphNetwork::evaluate(const NetIn &netIn) {
        return Network::evaluate(netIn);
    }
}
