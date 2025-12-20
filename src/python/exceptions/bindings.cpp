#include <pybind11/pybind11.h>

#include "bindings.h"

#include "gridfire/exceptions/error_scratchpad.h"

namespace py = pybind11;

#include "gridfire/exceptions/exceptions.h"

void register_exception_bindings(const py::module &m) {
    py::register_exception<gridfire::exceptions::GridFireError>(m, "GridFireError");

    py::register_exception<gridfire::exceptions::DebugException>(m, "DebugException", m.attr("GridFireError"));

    py::register_exception<gridfire::exceptions::EngineError>(m, "EngineError", m.attr("GridFireError"));

    py::register_exception<gridfire::exceptions::FailedToPartitionEngineError>(m, "FailedToPartitionEngineError", m.attr("EngineError"));
    py::register_exception<gridfire::exceptions::NetworkResizedError>(m, "NetworkResizedError", m.attr("EngineError"));
    py::register_exception<gridfire::exceptions::UnableToSetNetworkReactionsError>(m, "UnableToSetNetworkReactionsError", m.attr("EngineError"));
    py::register_exception<gridfire::exceptions::BadCollectionError>(m, "BadCollectionError", m.attr("EngineError"));
    py::register_exception<gridfire::exceptions::InvalidQSESolutionError>(m, "InvalidQSESolutionError", m.attr("EngineError"));
    py::register_exception<gridfire::exceptions::BadRHSEngineError>(m, "BadRHSEngineError", m.attr("EngineError"));

    py::register_exception<gridfire::exceptions::JacobianError>(m, "JacobianError", m.attr("EngineError"));

    py::register_exception<gridfire::exceptions::StaleJacobianError>(m, "StaleJacobianError", m.attr("JacobianError"));
    py::register_exception<gridfire::exceptions::UninitializedJacobianError>(m, "UninitializedJacobianError", m.attr("JacobianError"));
    py::register_exception<gridfire::exceptions::UnknownJacobianError>(m, "UnknownJacobianError", m.attr("JacobianError"));

    py::register_exception<gridfire::exceptions::UtilityError>(m, "UtilityError", m.attr("GridFireError"));
    py::register_exception<gridfire::exceptions::HashingError>(m, "HashingError", m.attr("UtilityError"));

    py::register_exception<gridfire::exceptions::PolicyError>(m, "PolicyError", m.attr("GridFireError"));
    py::register_exception<gridfire::exceptions::MissingBaseReactionError>(m, "MissingBaseReactionError", m.attr("PolicyError"));
    py::register_exception<gridfire::exceptions::MissingSeedSpeciesError>(m, "MissingSeedSpeciesError", m.attr("PolicyError"));
    py::register_exception<gridfire::exceptions::MissingKeyReactionError>(m, "MissingKeyReactionError", m.attr("PolicyError"));

    py::register_exception<gridfire::exceptions::ReactionError>(m, "ReactionError", m.attr("GridFireError"));
    py::register_exception<gridfire::exceptions::ReactionParsingError>(m, "ReactionParsingError", m.attr("ReactionError"));

    py::register_exception<gridfire::exceptions::SolverError>(m, "SolverError", m.attr("GridFireError"));
    py::register_exception<gridfire::exceptions::SingularJacobianError>(m, "SingularJacobianError", m.attr("SolverError"));
    py::register_exception<gridfire::exceptions::IllConditionedJacobianError>(m, "IllConditionedJacobianError", m.attr("SolverError"));
    py::register_exception<gridfire::exceptions::SUNDIALSError>(m, "SUNDIALSError", m.attr("SolverError"));
    py::register_exception<gridfire::exceptions::CVODESolverFailureError>(m, "CVODESolverFailureError", m.attr("SUNDIALSError"));
    py::register_exception<gridfire::exceptions::KINSolSolverFailureError>(m, "KINSolSolverFailureError", m.attr("SUNDIALSError"));

    py::register_exception<gridfire::exceptions::ScratchPadError>(m, "ScratchPadError", m.attr("GridFireError"));

}
