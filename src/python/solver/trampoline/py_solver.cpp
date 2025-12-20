#include "gridfire/solver/solver.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>
#include <tuple>
#include <string>
#include <any>

#include "py_solver.h"


namespace py = pybind11;

gridfire::NetOut PySingleZoneDynamicNetworkSolver::evaluate(
    gridfire::solver::SolverContextBase &solver_ctx,
    const gridfire::NetIn &netIn
) const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::NetOut,
        gridfire::solver::SingleZoneDynamicNetworkSolver,
        evaluate,
        solver_ctx,
        netIn
    );
}

std::vector<gridfire::NetOut> PyMultiZoneDynamicNetworkSolver::evaluate(
    gridfire::solver::SolverContextBase &solver_ctx,
    const std::vector<gridfire::NetIn> &netIns
) const {
    PYBIND11_OVERRIDE_PURE(
        std::vector<gridfire::NetOut>,
        gridfire::solver::MultiZoneDynamicNetworkSolver,
        evaluate,
        solver_ctx,
        netIns
    );
}

std::vector<std::tuple<std::string, std::string>> PyTimestepContextBase::describe() const {
    using ReturnType = std::vector<std::tuple<std::string, std::string>>;
    PYBIND11_OVERRIDE_PURE(
        ReturnType,
        gridfire::solver::TimestepContextBase,
        describe
    );
}

void PySolverContextBase::init() {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::solver::SolverContextBase,
        init
    );
}

void PySolverContextBase::set_stdout_logging(bool enable) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::solver::SolverContextBase,
        set_stdout_logging,
        enable
    );
}

void PySolverContextBase::set_detailed_logging(bool enable) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::solver::SolverContextBase,
        set_detailed_logging,
        enable
    );
}