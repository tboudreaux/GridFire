#pragma once

#include "gridfire/engine/engine.h"
#include "gridfire/expectations/expected_engine.h"

#include "fourdst/composition/atomicSpecies.h"

#include <vector>
#include <expected>


class PyEngine final : public gridfire::Engine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;
    std::expected<gridfire::StepDerivatives<double>,gridfire::expectations::StaleEngineError> calculateRHSAndEnergy(const std::vector<double> &Y, double T9, double rho) const override;
private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;
};

class PyDynamicEngine final : public gridfire::DynamicEngine {
public:
    const std::vector<fourdst::atomic::Species>& getNetworkSpecies() const override;
    std::expected<gridfire::StepDerivatives<double>,gridfire::expectations::StaleEngineError> calculateRHSAndEnergy(const std::vector<double> &Y, double T9, double rho) const override;
    void generateJacobianMatrix(const std::vector<double> &Y_dynamic, double T9, double rho) const override;
    void generateJacobianMatrix(const std::vector<double> &Y_dynamic, double T9, double rho, const gridfire::SparsityPattern &sparsityPattern) const override;
    double getJacobianMatrixEntry(int i, int j) const override;
    void generateStoichiometryMatrix() override;
    int getStoichiometryMatrixEntry(int speciesIndex, int reactionIndex) const override;
    double calculateMolarReactionFlow(const gridfire::reaction::Reaction &reaction, const std::vector<double> &Y, double T9, double rho) const override;
    const gridfire::reaction::LogicalReactionSet& getNetworkReactions() const override;
    void setNetworkReactions(const gridfire::reaction::LogicalReactionSet& reactions) override;
    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> getSpeciesTimescales(const std::vector<double> &Y, double T9, double rho) const override;
    std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> getSpeciesDestructionTimescales(const std::vector<double> &Y, double T9, double rho) const override;
    fourdst::composition::Composition update(const gridfire::NetIn &netIn) override;
    bool isStale(const gridfire::NetIn &netIn) override;
    void setScreeningModel(gridfire::screening::ScreeningType model) override;
    gridfire::screening::ScreeningType getScreeningModel() const override;
    int getSpeciesIndex(const fourdst::atomic::Species &species) const override;
    std::vector<double> mapNetInToMolarAbundanceVector(const gridfire::NetIn &netIn) const override;
    gridfire::PrimingReport primeEngine(const gridfire::NetIn &netIn) override;
    gridfire::BuildDepthType getDepth() const override {
        throw std::logic_error("Network depth not supported by this engine.");
    }
    void rebuild(const fourdst::composition::Composition& comp, gridfire::BuildDepthType depth) override {
        throw std::logic_error("Setting network depth not supported by this engine.");
    }
private:
    mutable std::vector<fourdst::atomic::Species> m_species_cache;


};

class PyEngineView final : public gridfire::EngineView<gridfire::Engine> {
    const gridfire::Engine& getBaseEngine() const override;
};

class PyDynamicEngineView final : public gridfire::EngineView<gridfire::DynamicEngine> {
    const gridfire::DynamicEngine& getBaseEngine() const override;
};