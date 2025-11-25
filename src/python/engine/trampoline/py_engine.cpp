#include "py_engine.h"

#include "gridfire/engine/engine.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> // Needed for std::function

#include <expected>
#include <unordered_map>
#include <vector>


namespace py = pybind11;

const std::vector<fourdst::atomic::Species>& PyEngine::getNetworkSpecies() const {
    /*
     * Acquire the GIL (Global Interpreter Lock) for thread safety
     * with the Python interpreter.
     */
    py::gil_scoped_acquire gil;

    /*
     * get_override() looks for a Python method that overrides this C++ one.
     */

    if (const py::function override = py::get_override(this, "getNetworkSpecies")) {
        const py::object result = override();
        m_species_cache = result.cast<std::vector<fourdst::atomic::Species>>();
        return m_species_cache;
    }

    py::pybind11_fail("Tried to call pure virtual function \"DynamicEngine::getNetworkSpecies\"");
}

std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> PyEngine::calculateRHSAndEnergy(
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus>),
        gridfire::engine::Engine,
        calculateRHSAndEnergy,
        comp, T9, rho
    );
}

///////////////////////////////////////
/// PyDynamicEngine Implementation ///
/////////////////////////////////////

const std::vector<fourdst::atomic::Species>& PyDynamicEngine::getNetworkSpecies() const {
    /*
     * Acquire the GIL (Global Interpreter Lock) for thread safety
     * with the Python interpreter.
     */
    py::gil_scoped_acquire gil;

    /*
     * get_override() looks for a Python method that overrides this C++ one.
     */

    if (const py::function override = py::get_override(this, "getNetworkSpecies")) {
        const py::object result = override();
        m_species_cache = result.cast<std::vector<fourdst::atomic::Species>>();
        return m_species_cache;
    }

    py::pybind11_fail("Tried to call pure virtual function \"DynamicEngine::getNetworkSpecies\"");
}


std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus> PyDynamicEngine::calculateRHSAndEnergy(
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::engine::StepDerivatives<double>, gridfire::engine::EngineStatus>),
        gridfire::engine::DynamicEngine,
        calculateRHSAndEnergy,
        comp, T9, rho
    );
}

gridfire::engine::NetworkJacobian PyDynamicEngine::generateJacobianMatrix(
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
    const fourdst::composition::CompositionAbstract &comp,
    const double T9,
    const double rho,
    const std::vector<fourdst::atomic::Species> &activeSpecies
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::NetworkJacobian,
        gridfire::engine::DynamicEngine,
        generateJacobianMatrix,
        comp,
        T9,
        rho,
        activeSpecies
    );
}

gridfire::engine::NetworkJacobian PyDynamicEngine::generateJacobianMatrix(
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho,
    const gridfire::engine::SparsityPattern &sparsityPattern
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::NetworkJacobian,
        gridfire::engine::DynamicEngine,
        generateJacobianMatrix,
        comp,
        T9,
        rho,
        sparsityPattern
    );
}

void PyDynamicEngine::generateStoichiometryMatrix() {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::engine::DynamicEngine,
        generateStoichiometryMatrix
    );
}

int PyDynamicEngine::getStoichiometryMatrixEntry(
    const fourdst::atomic::Species& species,
    const gridfire::reaction::Reaction& reaction
) const {
    PYBIND11_OVERRIDE_PURE(
        int,
        gridfire::engine::DynamicEngine,
        getStoichiometryMatrixEntry,
        species,
        reaction
    );
}

double PyDynamicEngine::calculateMolarReactionFlow(
    const gridfire::reaction::Reaction &reaction,
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::engine::DynamicEngine,
        calculateMolarReactionFlow,
        reaction,
        comp,
        T9,
        rho
    );
}

const gridfire::reaction::ReactionSet& PyDynamicEngine::getNetworkReactions() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::reaction::ReactionSet&,
        gridfire::engine::DynamicEngine,
        getNetworkReactions
    );
}

void PyDynamicEngine::setNetworkReactions(const gridfire::reaction::ReactionSet& reactions) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::engine::DynamicEngine,
        setNetworkReactions,
        reactions
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> PyDynamicEngine::getSpeciesTimescales(
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus>),
        gridfire::engine::DynamicEngine,
        getSpeciesTimescales,
        comp,
        T9,
        rho
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus> PyDynamicEngine::getSpeciesDestructionTimescales(
    const fourdst::composition::CompositionAbstract &comp,
    double T9,
    double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::engine::EngineStatus>),
        gridfire::engine::DynamicEngine,
        getSpeciesDestructionTimescales,
        comp, T9, rho
    );
}

fourdst::composition::Composition PyDynamicEngine::update(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        fourdst::composition::Composition,
        gridfire::engine::DynamicEngine,
        update,
        netIn
    );
}

bool PyDynamicEngine::isStale(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::engine::DynamicEngine,
        isStale,
        netIn
    );
}

void PyDynamicEngine::setScreeningModel(gridfire::screening::ScreeningType model) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::engine::DynamicEngine,
        setScreeningModel,
        model
    );
}

gridfire::screening::ScreeningType PyDynamicEngine::getScreeningModel() const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::screening::ScreeningType,
        gridfire::engine::DynamicEngine,
        getScreeningModel
    );
}

size_t PyDynamicEngine::getSpeciesIndex(const fourdst::atomic::Species &species) const {
    PYBIND11_OVERRIDE_PURE(
        int,
        gridfire::engine::DynamicEngine,
        getSpeciesIndex,
        species
    );
}

std::vector<double> PyDynamicEngine::mapNetInToMolarAbundanceVector(const gridfire::NetIn &netIn) const {
    PYBIND11_OVERRIDE_PURE(
        std::vector<double>,
        gridfire::engine::DynamicEngine,
        mapNetInToMolarAbundanceVector,
        netIn
    );
}

gridfire::engine::PrimingReport PyDynamicEngine::primeEngine(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::PrimingReport,
        gridfire::engine::DynamicEngine,
        primeEngine,
        netIn
    );
}

gridfire::engine::EnergyDerivatives PyDynamicEngine::calculateEpsDerivatives(
    const fourdst::composition::CompositionAbstract &comp,
    const double T9,
    const double rho) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::EnergyDerivatives,
        gridfire::engine::DynamicEngine,
        calculateEpsDerivatives,
        comp,
        T9,
        rho
    );
}

fourdst::composition::Composition PyDynamicEngine::collectComposition(
    const fourdst::composition::CompositionAbstract &comp,
    const double T9,
    const double rho
) const {
    PYBIND11_OVERRIDE_PURE(
        fourdst::composition::Composition,
        gridfire::engine::DynamicEngine,
        collectComposition,
        comp,
        T9,
        rho
    );
}

gridfire::engine::SpeciesStatus PyDynamicEngine::getSpeciesStatus(const fourdst::atomic::Species &species) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::SpeciesStatus,
        gridfire::engine::DynamicEngine,
        getSpeciesStatus,
        species
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


