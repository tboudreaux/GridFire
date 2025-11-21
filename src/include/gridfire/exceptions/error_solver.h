#pragma once

#include "gridfire/exceptions/gridfire_exception.h"

namespace gridfire::exceptions {
    class SolverError : GridFireError {
        using GridFireError::GridFireError;
    };

    class SingularJacobianError final : public SolverError {
        using SolverError::SolverError;
    };

    class IllConditionedJacobianError final : public SolverError {
        using SolverError::SolverError;
    };

    class SUNDIALSError : public SolverError {
        using SolverError::SolverError;
    };

    class CVODESolverFailureError final : public SUNDIALSError {
        using SUNDIALSError::SUNDIALSError;
    };

    class KINSolSolverFailureError final : public SUNDIALSError {
        using SUNDIALSError::SUNDIALSError;
    };

}