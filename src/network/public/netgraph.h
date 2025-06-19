#pragma once

#include "network.h"
#include "reaclib.h"
#include "atomicSpecies.h"
#include "composition.h"

#include <unordered_map>
#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix_sparse.hpp>


// PERF: The function getNetReactionStoichiometry returns a map of species to their stoichiometric coefficients for a given reaction.
//       this makes extra copies of the species, which is not ideal and could be optimized further.
//       Even more relevant is the member m_reactionIDMap which makes copies of a REACLIBReaction for each reaction ID.
//       REACLIBReactions are quite large data structures, so this could be a performance bottleneck.

namespace serif::network {
    template <typename GeneralScalarType>
    class GraphNetwork final : public Network {
    public:
        explicit GraphNetwork(const serif::composition::Composition &composition);
        explicit GraphNetwork(const serif::composition::Composition &composition, const double cullingThreshold, const double T9);

        NetOut evaluate(const NetIn &netIn) override;

        [[nodiscard]] const std::vector<serif::atomic::Species>& getNetworkSpecies() const;
        [[nodiscard]] const reaclib::REACLIBReactionSet& getNetworkReactions() const;
        [[nodiscard]] std::unordered_map<serif::atomic::Species, int> getNetReactionStoichiometry(const reaclib::REACLIBReaction& reaction) const;

        [[nodiscard]] bool involvesSpecies(const serif::atomic::Species& species) const;

        std::vector<std::vector<std::string>> detectCycles() const;

        std::vector<std::string> topologicalSortReactions() const;
    private:
        reaclib::REACLIBReactionSet m_reactions; ///< Set of REACLIB reactions for the network.
        serif::composition::Composition m_initialComposition;
        serif::composition::Composition m_currentComposition;
        GeneralScalarType m_cullingThreshold; ///< Threshold for culling reactions.
        GeneralScalarType m_T9; ///< Temperature in T9 units to evaluate culling.

        std::vector<serif::atomic::Species> m_networkSpecies; ///< The species in the network.
        std::unordered_map<std::string_view, serif::atomic::Species> m_networkSpeciesMap;
        std::unordered_map<std::string_view, const reaclib::REACLIBReaction> m_reactionIDMap; ///< Map of reaction IDs to REACLIB reactions.

        boost::numeric::ublas::compressed_matrix<int> m_stoichiometryMatrix; ///< Stoichiometry matrix for the network.
        boost::numeric::ublas::compressed_matrix<GeneralScalarType> m_jacobianMatrix; ///< Jacobian matrix for the network.
        std::unordered_map<serif::atomic::Species, size_t> m_speciesToIndexMap; ///< Map of species to their index in the stoichiometry matrix.


    private:
        void syncInternalMaps();
        void collectNetworkSpecies();
        void populateReactionIDMap();
        void populateSpeciesToIndexMap();
        void buildNetworkGraph();
        bool validateConservation() const;

        // --- Generate the system matrices ---
        void generateStoichiometryMatrix();
        void generateJacobianMatrix();

        std::vector<GeneralScalarType> calculateRHS(const std::vector<GeneralScalarType>& Y, const GeneralScalarType T9, const GeneralScalarType rho) const;
        GeneralScalarType calculateReactionRate(const reaclib::REACLIBReaction &reaction, const std::vector<GeneralScalarType>& Y, const GeneralScalarType T9, const GeneralScalarType rho) const;

        void pruneNetworkGraph();
    };
};