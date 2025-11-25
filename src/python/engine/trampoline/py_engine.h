#pragma once

#include "gridfire/engine/engine.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <vector>
#include <expected>


class PyEngine final : public gridfire::engine::Engine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

    std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> calculateRHSAndEnergy(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;
private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;
};

class PyDynamicEngine final : public gridfire::engine::DynamicEngine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

    std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> calculateRHSAndEnergy(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    gridfire::engine::NetworkJacobian generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract& comp,
        double T9,
        double rho
    ) const override;

    gridfire::engine::NetworkJacobian generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho,
        const std::vector<fourdst::atomic::Species> &activeSpecies
    ) const override;

    gridfire::engine::NetworkJacobian generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract& comp,
        double T9,
        double rho,
        const gridfire::engine::SparsityPattern &sparsityPattern
    ) const override;

    void generateStoichiometryMatrix() override;

    int getStoichiometryMatrixEntry(
        const fourdst::atomic::Species& species,
        const gridfire::reaction::Reaction& reaction
    ) const override;

    double calculateMolarReactionFlow(
        const gridfire::reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    const gridfire::reaction::ReactionSet& getNetworkReactions() const override;

    void setNetworkReactions(
        const gridfire::reaction::ReactionSet& reactions
    ) override;

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> getSpeciesTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> getSpeciesDestructionTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    fourdst::composition::Composition update(
        const gridfire::NetIn &netIn
    ) override;

    bool isStale(
        const gridfire::NetIn &netIn
    ) override;

    void setScreeningModel(
        gridfire::screening::ScreeningType model
    ) override;

    gridfire::screening::ScreeningType getScreeningModel() const override;

    size_t getSpeciesIndex(
        const fourdst::atomic::Species &species
    ) const override;

    std::vector<double> mapNetInToMolarAbundanceVector(
        const gridfire::NetIn &netIn
    ) const override;

    gridfire::engine::PrimingReport primeEngine(
        const gridfire::NetIn &netIn
    ) override;

    gridfire::engine::BuildDepthType getDepth() const override {
        throw std::logic_error("Network depth not supported by this engine.");
    }
    void rebuild(
        const fourdst::composition::CompositionAbstract &comp,
        gridfire::engine::BuildDepthType depth
    ) override {
        throw std::logic_error("Setting network depth not supported by this engine.");
    }

    [[nodiscard]] gridfire::engine::EnergyDerivatives calculateEpsDerivatives(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    fourdst::composition::Composition collectComposition(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    gridfire::engine::SpeciesStatus getSpeciesStatus(
        const fourdst::atomic::Species &species
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