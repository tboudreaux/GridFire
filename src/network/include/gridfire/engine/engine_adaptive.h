#pragma once
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_view_abstract.h"

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/logging/logging.h"

#include "quill/Logger.h"

namespace gridfire {
    class AdaptiveEngineView final : public DynamicEngine, public EngineView<DynamicEngine> {
    public:
        explicit AdaptiveEngineView(DynamicEngine& baseEngine);

        void update(const NetIn& netIn);

        const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;
        StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double> &Y,
            const double T9,
            const double rho
        ) const override;

        void generateJacobianMatrix(
            const std::vector<double> &Y,
            const double T9,
            const double rho
        ) override;
        double getJacobianMatrixEntry(
            const int i,
            const int j
        ) const override;

        void generateStoichiometryMatrix() override;
        int getStoichiometryMatrixEntry(
            const int speciesIndex,
            const int reactionIndex
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

        const DynamicEngine& getBaseEngine() const override { return m_baseEngine; }
    private:
        using Config = fourdst::config::Config;
        using LogManager = fourdst::logging::LogManager;

        DynamicEngine& m_baseEngine;

        std::vector<fourdst::atomic::Species> m_activeSpecies;
        reaction::REACLIBLogicalReactionSet m_activeReactions;

        std::vector<size_t> m_speciesIndexMap;

        bool m_isStale = true;

        Config& m_config = Config::getInstance();
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
    private:
        std::vector<size_t> constructSpeciesIndexMap() const;
    private:
        struct ReactionFlow {
            const reaction::Reaction* reactionPtr;
            double flowRate;
        };
    };
}
