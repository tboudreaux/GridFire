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

std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError> PyEngine::calculateRHSAndEnergy(const fourdst::composition::Composition &comp, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError>),
        gridfire::Engine,
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
std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError> PyDynamicEngine::calculateRHSAndEnergy(const fourdst::composition::Composition &comp, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<gridfire::StepDerivatives<double>, gridfire::expectations::StaleEngineError>),
        gridfire::Engine,
        calculateRHSAndEnergy,
        comp, T9, rho
    );
}

void PyDynamicEngine::generateJacobianMatrix(const fourdst::composition::Composition& comp, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateJacobianMatrix,
        comp,
        T9,
        rho
    );
}

void PyDynamicEngine::generateJacobianMatrix(
    const fourdst::composition::Composition &comp,
    const double T9,
    const double rho,
    const std::vector<fourdst::atomic::Species> &activeSpecies
) const {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateJacobianMatrix,
        comp,
        T9,
        rho,
        activeSpecies
    );
}

void PyDynamicEngine::generateJacobianMatrix(const fourdst::composition::Composition &comp, double T9, double rho, const gridfire::SparsityPattern &sparsityPattern) const {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateJacobianMatrix,
        comp,
        T9,
        rho,
        sparsityPattern
    );
}

double PyDynamicEngine::getJacobianMatrixEntry(const fourdst::atomic::Species& rowSpecies, const fourdst::atomic::Species& colSpecies) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::DynamicEngine,
        getJacobianMatrixEntry,
        rowSpecies,
        colSpecies
    );
}

void PyDynamicEngine::generateStoichiometryMatrix() {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        generateStoichiometryMatrix
    );
}

int PyDynamicEngine::getStoichiometryMatrixEntry(const fourdst::atomic::Species& species, const gridfire::reaction::Reaction& reaction) const {
    PYBIND11_OVERRIDE_PURE(
        int,
        gridfire::DynamicEngine,
        getStoichiometryMatrixEntry,
        species,
        reaction
    );
}

double PyDynamicEngine::calculateMolarReactionFlow(const gridfire::reaction::Reaction &reaction, const fourdst::composition::Composition &comp, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::DynamicEngine,
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
        gridfire::DynamicEngine,
        getNetworkReactions
    );
}

void PyDynamicEngine::setNetworkReactions(const gridfire::reaction::ReactionSet& reactions) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::DynamicEngine,
        setNetworkReactions,
        reactions
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> PyDynamicEngine::getSpeciesTimescales(const fourdst::composition::Composition &comp, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError>),
        gridfire::DynamicEngine,
        getSpeciesTimescales,
        comp,
        T9,
        rho
    );
}

std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError> PyDynamicEngine::getSpeciesDestructionTimescales(const fourdst::composition::Composition &comp, double T9, double rho) const {
    PYBIND11_OVERRIDE_PURE(
        PYBIND11_TYPE(std::expected<std::unordered_map<fourdst::atomic::Species, double>, gridfire::expectations::StaleEngineError>),
        gridfire::DynamicEngine,
        getSpeciesDestructionTimescales,
        comp, T9, rho
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

size_t PyDynamicEngine::getSpeciesIndex(const fourdst::atomic::Species &species) const {
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

gridfire::EnergyDerivatives PyDynamicEngine::calculateEpsDerivatives(
    const fourdst::composition::Composition &comp,
    const double T9,
    const double rho) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::EnergyDerivatives,
        gridfire::DynamicEngine,
        calculateEpsDerivatives,
        comp,
        T9,
        rho
    );
}

fourdst::composition::Composition PyDynamicEngine::collectComposition(
    fourdst::composition::Composition &comp
) const {
    PYBIND11_OVERRIDE_PURE(
        fourdst::composition::Composition,
        gridfire::DynamicEngine,
        collectComposition,
        comp
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


