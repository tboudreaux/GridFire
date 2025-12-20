#include "py_engine.h"

#include "gridfire/engine/engine.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> // Needed for std::function

#include <expected>
#include <unordered_map>
#include <vector>


namespace py = pybind11;

const std::vector<fourdst::atomic::Species>& PyEngine::getNetworkSpecies(
    gridfire::engine::scratch::StateBlob& ctx
) const {
    PYBIND11_OVERRIDE_PURE(
        const std::vector<fourdst::atomic::Species>&,
        gridfire::engine::Engine,
        getNetworkSpecies,
        ctx
    );
}

std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> PyEngine::calculateRHSAndEnergy(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho,
    bool trust
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus>),
        gridfire::engine::Engine,
        calculateRHSAndEnergy,
        ctx, comp, T9, rho, trust
    );
}

///////////////////////////////////////
/// PyDynamicEngine Implementation ///
/////////////////////////////////////

const std::vector<fourdst::atomic::Species>& PyDynamicEngine::getNetworkSpecies(
    gridfire::engine::scratch::StateBlob& ctx
) const {
        PYBIND11_OVERRIDE_PURE(
        const std::vector<fourdst::atomic::Species>&,
        gridfire::engine::DynamicEngine,
        getNetworkSpecies,
        ctx
    );
}


std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> PyDynamicEngine::calculateRHSAndEnergy(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho,
    bool trust
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus>),
        gridfire::engine::DynamicEngine,
        calculateRHSAndEnergy,
        ctx, comp, T9, rho, trust
    );
}

gridfire::engine::NetworkJacobian PyDynamicEngine::generateJacobianMatrix(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract& comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::NetworkJacobian,
        gridfire::engine::DynamicEngine,
        generateJacobianMatrix,
        comp,
        T9,
        rho
    );
}

gridfire::engine::NetworkJacobian PyDynamicEngine::generateJacobianMatrix(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    const double T9,
    const double rho,
    const std::vector<fourdst::atomic::Species> &activeSpecies
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::NetworkJacobian,
        gridfire::engine::DynamicEngine,
        generateJacobianMatrix,
        ctx,
        comp,
        T9,
        rho,
        activeSpecies
    );
}

gridfire::engine::NetworkJacobian PyDynamicEngine::generateJacobianMatrix(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho,
    const gridfire::engine::SparsityPattern &sparsityPattern
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::NetworkJacobian,
        gridfire::engine::DynamicEngine,
        generateJacobianMatrix,
        ctx,
        comp,
        T9,
        rho,
        sparsityPattern
    );
}

double PyDynamicEngine::calculateMolarReactionFlow(
    gridfire::engine::scratch::StateBlob& ctx,
    const gridfire::reaction::Reaction &reaction,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::engine::DynamicEngine,
        calculateMolarReactionFlow,
        ctx,
        reaction,
        comp,
        T9,
        rho
    );
}

const gridfire::reaction::ReactionSet& PyDynamicEngine::getNetworkReactions(
    gridfire::engine::scratch::StateBlob& ctx
) const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::reaction::ReactionSet&,
        gridfire::engine::DynamicEngine,
        getNetworkReactions,
        ctx
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> PyDynamicEngine::getSpeciesTimescales(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus>),
        gridfire::engine::DynamicEngine,
        getSpeciesTimescales,
        ctx,
        comp,
        T9,
        rho
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> PyDynamicEngine::getSpeciesDestructionTimescales(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus>),
        gridfire::engine::DynamicEngine,
        getSpeciesDestructionTimescales,
        ctx, comp, T9, rho
    );
}

fourdst::composition::Composition PyDynamicEngine::project(
    gridfire::engine::scratch::StateBlob& ctx,
    const gridfire::NetIn &netIn
) const {
    PYBIND11_OVERRIDE_PURE(
        fourdst::composition::Composition,
        gridfire::engine::DynamicEngine,
        project,
        ctx,
        netIn
    );
}

gridfire::screening::ScreeningType PyDynamicEngine::getScreeningModel(
    gridfire::engine::scratch::StateBlob& ctx
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::screening::ScreeningType,
        gridfire::engine::DynamicEngine,
        getScreeningModel,
        ctx
    );
}

size_t PyDynamicEngine::getSpeciesIndex(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::atomic::Species &species
) const {
    PYBIND11_OVERRIDE_PURE(
        int,
        gridfire::engine::DynamicEngine,
        getSpeciesIndex,
        ctx,
        species
    );
}

gridfire::engine::PrimingReport PyDynamicEngine::primeEngine(
    gridfire::engine::scratch::StateBlob& ctx,
    const gridfire::NetIn &netIn
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::PrimingReport,
        gridfire::engine::DynamicEngine,
        primeEngine,
        ctx,
        netIn
    );
}

gridfire::engine::EnergyDerivatives PyDynamicEngine::calculateEpsDerivatives(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    const double T9,
    const double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::EnergyDerivatives,
        gridfire::engine::DynamicEngine,
        calculateEpsDerivatives,
        ctx,
        comp,
        T9,
        rho
    );
}

fourdst::composition::Composition PyDynamicEngine::collectComposition(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::composition::CompositionAbstract &comp,
    const double T9,
    const double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        fourdst::composition::Composition,
        gridfire::engine::DynamicEngine,
        collectComposition,
        ctx,
        comp,
        T9,
        rho
    );
}

gridfire::engine::SpeciesStatus PyDynamicEngine::getSpeciesStatus(
    gridfire::engine::scratch::StateBlob& ctx,
    const fourdst::atomic::Species &species
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::SpeciesStatus,
        gridfire::engine::DynamicEngine,
        getSpeciesStatus,
        ctx,
        species
    );
}

std::optional<gridfire::engine::StepDerivatives<double>> PyDynamicEngine::getMostRecentRHSCalculation(
    gridfire::engine::scratch::StateBlob &ctx
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::optional<gridfire::engine::StepDerivatives<double>>),
        gridfire::engine::DynamicEngine,
        getMostRecentRHSCalculation,
        ctx
    );
}

const gridfire::engine::Engine& PyEngineView::getBaseEngine() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::engine::Engine&,
        gridfire::engine::EngineView<gridfire::engine::Engine>,
        getBaseEngine
    );
}

const gridfire::engine::DynamicEngine& PyDynamicEngineView::getBaseEngine() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::engine::DynamicEngine&,
        gridfire::engine::EngineView<gridfire::engine::DynamicEngine>,
        getBaseEngine
    );
}


