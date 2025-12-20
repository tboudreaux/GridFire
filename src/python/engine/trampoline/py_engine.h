#pragma once

#include "gridfire/engine/engine.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <vector>
#include <expected>


class PyEngine final : public gridfire::engine::Engine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies(
        gridfire::engine::scratch::StateBlob& ctx
    ) const override;

    std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> calculateRHSAndEnergy(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho,
        bool trust
    ) const override;
private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;
};

class PyDynamicEngine final : public gridfire::engine::DynamicEngine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies(
        gridfire::engine::scratch::StateBlob& ctx
    ) const override;

    std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> calculateRHSAndEnergy(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho,
        bool trust
    ) const override;

    gridfire::engine::NetworkJacobian generateJacobianMatrix(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract& comp,
        double T9,
        double rho
    ) const override;

    gridfire::engine::NetworkJacobian generateJacobianMatrix(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho,
        const std::vector<fourdst::atomic::Species> &activeSpecies
    ) const override;

    gridfire::engine::NetworkJacobian generateJacobianMatrix(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract& comp,
        double T9,
        double rho,
        const gridfire::engine::SparsityPattern &sparsityPattern
    ) const override;

    double calculateMolarReactionFlow(
        gridfire::engine::scratch::StateBlob& ctx,
        const gridfire::reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    const gridfire::reaction::ReactionSet& getNetworkReactions(
        gridfire::engine::scratch::StateBlob& ctx
    ) const override;

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> getSpeciesTimescales(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> getSpeciesDestructionTimescales(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    fourdst::composition::Composition project(
        gridfire::engine::scratch::StateBlob& ctx,
        const gridfire::NetIn &netIn
    ) const override;

    gridfire::screening::ScreeningType getScreeningModel(
        gridfire::engine::scratch::StateBlob& ctx
    ) const override;

    size_t getSpeciesIndex(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::atomic::Species &species
    ) const override;

    gridfire::engine::PrimingReport primeEngine(
        gridfire::engine::scratch::StateBlob& ctx,
        const gridfire::NetIn &netIn
    ) const override;

    [[nodiscard]] gridfire::engine::EnergyDerivatives calculateEpsDerivatives(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    fourdst::composition::Composition collectComposition(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    gridfire::engine::SpeciesStatus getSpeciesStatus(
        gridfire::engine::scratch::StateBlob& ctx,
        const fourdst::atomic::Species &species
    ) const override;

    std::optional<gridfire::engine::StepDerivatives<double>> getMostRecentRHSCalculation(
        gridfire::engine::scratch::StateBlob &ctx
    ) const override;

private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;
};

class PyEngineView final : public gridfire::engine::EngineView<gridfire::engine::Engine> {
    [[nodiscard]] const gridfire::engine::Engine& getBaseEngine() const override;
};

class PyDynamicEngineView final : public gridfire::engine::EngineView<gridfire::engine::DynamicEngine> {
    [[nodiscard]] const gridfire::engine::DynamicEngine& getBaseEngine() const override;
};