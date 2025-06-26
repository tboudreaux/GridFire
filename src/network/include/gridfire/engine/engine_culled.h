#pragma once
#include "gridfire/engine/engine_abstract.h"

#include "fourdst/composition/atomicSpecies.h"

namespace gridfire {
    class CulledEngine final : public DynamicEngine {
    public:
        explicit CulledEngine(DynamicEngine& baseEngine);

        const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;
        StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        void generateJacobianMatrix(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) override;
        double getJacobianMatrixEntry(
            int i,
            int j
        ) const override;

        void generateStoichiometryMatrix() override;
        int getStoichiometryMatrixEntry(
            int speciesIndex,
            int reactionIndex
        ) const override;

        double calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        const reaction::REACLIBLogicalReactionSet& getNetworkReactions() const override;
        std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;
    private:
        DynamicEngine& m_baseEngine;
        std::unordered_map<fourdst::atomic::Species, size_t> m_fullSpeciesToIndexMap;
        std::vector<fourdst::atomic::Species> m_culledSpecies;
    private:
        std::unordered_map<fourdst::atomic::Species, size_t> populatedSpeciesToIndexMap;
        std::vector<fourdst::atomic::Species> determineCullableSpecies;
    };
}
