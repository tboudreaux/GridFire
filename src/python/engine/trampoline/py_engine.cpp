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
    py::function override = py::get_override(this, "getNetworkSpecies");

    if (override) {
        py::object result = override();
        m_species_cache = result.cast<std::vector<fourdst::atomic::Species>>();
        return m_species_cache;
    }

    py::pybind11_fail("Tried to call pure virtual function \"DynamicEngine::getNetworkSpecies\"");
}

std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError> PyEngine::calculateRHSAndEnergy(const std::vector<double> &Y, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError>),
        gridfire::Engine,
        calculateRHSAndEnergy,
        Y, T9, rho
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
    py::function override = py::get_override(this, "getNetworkSpecies");

    if (override) {
        py::object result = override();
        m_species_cache = result.cast<std::vector<fourdst::atomic::Species>>();
        return m_species_cache;
    }

    py::pybind11_fail("Tried to call pure virtual function \"DynamicEngine::getNetworkSpecies\"");
}
std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError> PyDynamicEngine::calculateRHSAndEnergy(const std::vector<double> &Y, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError>),
        gridfire::Engine,
        calculateRHSAndEnergy,
        Y, T9, rho
    );
}

void PyDynamicEngine::generateJacobianMatrix(const std::vector<double> &Y_dynamic, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateJacobianMatrix,
        Y_dynamic, T9, rho
    );
}

void PyDynamicEngine::generateJacobianMatrix(const std::vector<double> &Y_dynamic, double T9, double rho, const gridfire::SparsityPattern &sparsityPattern) const {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateJacobianMatrix,
        Y_dynamic, T9, rho, sparsityPattern
    );
}

double PyDynamicEngine::getJacobianMatrixEntry(int i, int j) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::DynamicEngine,
        getJacobianMatrixEntry,
        i, j
    );
}

void PyDynamicEngine::generateStoichiometryMatrix() {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateStoichiometryMatrix
    );
}

int PyDynamicEngine::getStoichiometryMatrixEntry(int speciesIndex, int reactionIndex) const {
    PYBIND11_OVERRIDE_PURE(
        int,
        gridfire::DynamicEngine,
        getStoichiometryMatrixEntry,
        speciesIndex, reactionIndex
    );
}

double PyDynamicEngine::calculateMolarReactionFlow(const gridfire::reaction::Reaction &reaction, const std::vector<double> &Y, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::DynamicEngine,
        calculateMolarReactionFlow,
        reaction, Y, T9, rho
    );
}

const gridfire::reaction::LogicalReactionSet& PyDynamicEngine::getNetworkReactions() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::reaction::LogicalReactionSet&,
        gridfire::DynamicEngine,
        getNetworkReactions
    );
}

void PyDynamicEngine::setNetworkReactions(const gridfire::reaction::LogicalReactionSet& reactions) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        setNetworkReactions,
        reactions
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> PyDynamicEngine::getSpeciesTimescales(const std::vector<double> &Y, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError>),
        gridfire::DynamicEngine,
        getSpeciesTimescales,
        Y, T9, rho
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> PyDynamicEngine::getSpeciesDestructionTimescales(const std::vector<double> &Y, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError>),
        gridfire::DynamicEngine,
        getSpeciesDestructionTimescales,
        Y, T9, rho
    );
}

fourdst::composition::Composition PyDynamicEngine::update(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        fourdst::composition::Composition,
        gridfire::DynamicEngine,
        update,
        netIn
    );
}

bool PyDynamicEngine::isStale(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::DynamicEngine,
        isStale,
        netIn
    );
}

void PyDynamicEngine::setScreeningModel(gridfire::screening::ScreeningType model) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        setScreeningModel,
        model
    );
}

gridfire::screening::ScreeningType PyDynamicEngine::getScreeningModel() const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::screening::ScreeningType,
        gridfire::DynamicEngine,
        getScreeningModel
    );
}

int PyDynamicEngine::getSpeciesIndex(const fourdst::atomic::Species &species) const {
    PYBIND11_OVERRIDE_PURE(
        int,
        gridfire::DynamicEngine,
        getSpeciesIndex,
        species
    );
}

std::vector<double> PyDynamicEngine::mapNetInToMolarAbundanceVector(const gridfire::NetIn &netIn) const {
    PYBIND11_OVERRIDE_PURE(
        std::vector<double>,
        gridfire::DynamicEngine,
        mapNetInToMolarAbundanceVector,
        netIn
    );
}

gridfire::PrimingReport PyDynamicEngine::primeEngine(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        gridfire::PrimingReport,
        gridfire::DynamicEngine,
        primeEngine,
        netIn
    );
}

const gridfire::Engine& PyEngineView::getBaseEngine() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::Engine&,
        gridfire::EngineView<gridfire::Engine>,
        getBaseEngine
    );
}

const gridfire::DynamicEngine& PyDynamicEngineView::getBaseEngine() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::DynamicEngine&,
        gridfire::EngineView<gridfire::DynamicEngine>,
        getBaseEngine
    );
}


