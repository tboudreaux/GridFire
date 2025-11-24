#pragma once

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    /**
     * @brief Base class for engine-related exceptions.
     *
     * This class serves as the base for all exceptions specific to the
     * reaction network engine in the GridFire library. It extends the
     * GridFireError class and allows for custom error messages related
     * to engine operations.
     */
    class EngineError : public GridFireError {
        using GridFireError::GridFireError;
    };

    /**
     * @brief Exception for failures in partitioning the engine.
     *
     * This exception is thrown when the engine fails to partition
     * the reaction network as required for certain calculations or
     * optimizations.
     */
    class FailedToPartitionEngineError final : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Exception for errors during network resizing.
     *
     * This exception is thrown when the engine encounters an error
     * while attempting to resize the reaction network, such as when
     * adding or removing species or reactions.
     */
    class NetworkResizedError final : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Exception for failures in setting network reactions.
     *
     * This exception is thrown when the engine fails to properly
     * set or initialize the reactions in the reaction network.
     */
    class UnableToSetNetworkReactionsError final : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Exception for invalid composition collection in the engine.
     *
     * This exception is thrown when the engine encounters an
     * invalid state while trying to collect the composition from the entire engine stack
     */
    class BadCollectionError final : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Exception for invalid QSE solution in the engine.
     *
     * This exception is thrown when the engine computes an
     * invalid solution while performing Quasi-Statistical Equilibrium (QSE) calculations.
     */
    class InvalidQSESolutionError final : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Exception for errors in calculating the right-hand side (RHS).
     *
     * This exception is thrown when the engine encounters an error
     * while calculating the right-hand side of the reaction network equations.
     */
    class BadRHSEngineError final : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Base class for Jacobian-related exceptions.
     *
     * This class serves as the base for all exceptions specific to
     * Jacobian matrix operations in the reaction network engine.
     */
    class JacobianError : public EngineError {
        using EngineError::EngineError;
    };

    /**
     * @brief Exception for stale Jacobian matrix access.
     *
     * This exception is thrown when an attempt is made to access
     * a Jacobian matrix that is stale and needs to be regenerated.
     */
    class StaleJacobianError final : public JacobianError {
        using JacobianError::JacobianError;
    };

    /**
     * @brief Exception for uninitialized Jacobian matrix access.
     *
     * This exception is thrown when an attempt is made to access
     * a Jacobian matrix that has not been initialized.
     */
    class UninitializedJacobianError final: public JacobianError {
        using JacobianError::JacobianError;
    };

    /**
     * @brief Exception for unknown Jacobian matrix access.
     *
     * This exception is thrown when an attempt is made to access
     * a Jacobian matrix that is unknown or not recognized by the engine.
     */
    class UnknownJacobianError final : public JacobianError {
        using JacobianError::JacobianError;
    };


}