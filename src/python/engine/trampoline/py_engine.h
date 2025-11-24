#pragma once

#include "gridfire/engine/engine.h"
#include "gridfire/expectations/expected_engine.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <vector>
#include <expected>


class PyEngine final : public gridfire::Engine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

    std::expected<gridfire::StepDerivatives<double>,gridfire::expectations::StaleEngineError> calculateRHSAndEnergy(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;
private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;
};

class PyDynamicEngine final : public gridfire::DynamicEngine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;

    std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError> calculateRHSAndEnergy(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    void generateJacobianMatrix(
        const fourdst::composition::Composition& comp,
        double T9,
        double rho
    ) const override;

    void generateJacobianMatrix(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho,
        const std::vector<fourdst::atomic::Species> &activeSpecies
    ) const override;

    void generateJacobianMatrix(
        const fourdst::composition::Composition& comp,
        double T9,
        double rho,
        const gridfire::SparsityPattern &sparsityPattern
    ) const override;

    double getJacobianMatrixEntry(
        const fourdst::atomic::Species& rowSpecies,
        const fourdst::atomic::Species& colSpecies
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

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> getSpeciesTimescales(
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const override;

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> getSpeciesDestructionTimescales(
        const fourdst::composition::Composition &comp,
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

    gridfire::PrimingReport primeEngine(
        const gridfire::NetIn &netIn
    ) override;

    gridfire::BuildDepthType getDepth() const override {
        throw std::logic_error("Network depth not supported by this engine.");
    }
    void rebuild(
        const fourdst::composition::CompositionAbstract &comp,
        gridfire::BuildDepthType depth
    ) override {
        throw std::logic_error("Setting network depth not supported by this engine.");
    }

    [[nodiscard]] gridfire::EnergyDerivatives calculateEpsDerivatives(
        const fourdst::composition::Composition &comp,
        double T9,
        double rho
    ) const override;

    fourdst::composition::Composition collectComposition(
        fourdst::composition::Composition &comp
    ) const override;

private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;
};

class PyEngineView final : public gridfire::EngineView<gridfire::Engine> {
    [[nodiscard]] const gridfire::Engine& getBaseEngine() const override;
};

class PyDynamicEngineView final : public gridfire::EngineView<gridfire::DynamicEngine> {
    [[nodiscard]] const gridfire::DynamicEngine& getBaseEngine() const override;
};