#pragma once

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    /**
     * @class SolverError
     * @brief Base class for all exceptions related to solver operations.
     *
     * This exception is the parent for more specific solver-related errors. Catching this
     * type will catch any exception originating from the solver system.
     */
    class SolverError : public GridFireError {
        using GridFireError::GridFireError;
    };

    /**
     * @class SingularJacobianError
     * @brief Exception thrown when the Jacobian matrix is singular.
     *
     * This exception is thrown when the solver encounters a singular Jacobian matrix,
     * which prevents it from proceeding with the solution of the system.
     */
    class SingularJacobianError final : public SolverError {
        using SolverError::SolverError;
    };

    /**
     * @class IllConditionedJacobianError
     * @brief Exception thrown when the Jacobian matrix is ill-conditioned.
     *
     * This exception is thrown when the solver detects that the Jacobian matrix
     * is ill-conditioned, which may lead to inaccurate or unstable solutions.
     */
    class IllConditionedJacobianError final : public SolverError {
        using SolverError::SolverError;
    };

    /**
     * @class SUNDIALSError
     * @brief Base class for exceptions related to SUNDIALS solver operations.
     *
     * This class serves as the base for all exceptions specific to
     * SUNDIALS solver errors in the GridFire library.
     */
    class SUNDIALSError : public SolverError {
        using SolverError::SolverError;
    };

    /**
     * @class CVODESolverFailureError
     * @brief Exception for failures in the CVODE solver.
     *
     * This exception is thrown when the CVODE solver from the SUNDIALS suite
     * encounters a failure during its operation.
     */
    class CVODESolverFailureError final : public SUNDIALSError {
        using SUNDIALSError::SUNDIALSError;
    };

    /**
     * @class KINSolSolverFailureError
     * @brief Exception for failures in the KINSOL solver.
     *
     * This exception is thrown when the KINSOL solver from the SUNDIALS suite
     * encounters a failure during its operation.
     */
    class KINSolSolverFailureError final : public SUNDIALSError {
        using SUNDIALSError::SUNDIALSError;
    };

}