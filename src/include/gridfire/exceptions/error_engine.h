#pragma once

#include <exception>
#include <string>
#include <utility>
#include <vector>

#include "gridfire/exceptions/gridfire_exception.h"

namespace gridfire::exceptions {
    class EngineError : public GridFireError {
        using GridFireError::GridFireError;
    };

    class FailedToPartitionEngineError final : public EngineError {
        using EngineError::EngineError;
    };

    class NetworkResizedError final : public EngineError {
        using EngineError::EngineError;
    };

    class UnableToSetNetworkReactionsError final : public EngineError {
        using EngineError::EngineError;
    };

    class BadCollectionError final : public EngineError {
        using EngineError::EngineError;
    };

    class InvalidQSESolutionError final : public EngineError {
        using EngineError::EngineError;
    };

    class BadRHSEngineError final : public EngineError {
        using EngineError::EngineError;
    };

    class JacobianError : public EngineError {
        using EngineError::EngineError;
    };

    class StaleJacobianError final : public JacobianError {
        using JacobianError::JacobianError;
    };

    class UninitializedJacobianError final: public JacobianError {
        using JacobianError::JacobianError;
    };

    class UnknownJacobianError final : public JacobianError {
        using JacobianError::JacobianError;
    };


}