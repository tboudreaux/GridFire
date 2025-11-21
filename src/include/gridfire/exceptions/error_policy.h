/**
 * @file error_policy.h
 * @brief Defines custom exception types related to network policy construction and verification.
 *
 * This file contains a hierarchy of exception classes that are thrown by `NetworkPolicy`
 * implementations (see `stellar_policy.h`) and `ReactionChainPolicy` implementations
 * (see `chains.h`) when errors occur.
 *
 * The base class `PolicyError` inherits from `std::exception`, and specific error
 * conditions are represented by derived classes.
 */
#pragma once

#include <exception>
#include <string>
#include <utility>

#include "gridfire/exceptions/gridfire_exception.h"

namespace gridfire::exceptions {
    /**
     * @class PolicyError
     * @brief Base class for all exceptions related to network policy operations.
     *
     * This exception is the parent for more specific policy-related errors. Catching this
     * type will catch any exception originating from the policy system.
     */
    class PolicyError : public GridFireError {
        using GridFireError::GridFireError;
    };

    /**
     * @class MissingBaseReactionError
     * @brief Exception thrown when a required reaction is missing from the underlying database (e.g., REACLIB).
     *
     * This typically occurs during the construction of a `ReactionChainPolicy` if the
     * reaction library used by GridFire does not contain a reaction specified by the policy.
     */
    class MissingBaseReactionError final : public PolicyError {
        using PolicyError::PolicyError;
    };

    /**
     * @class MissingSeedSpeciesError
     * @brief Exception thrown when a required seed species is not found in the initial composition.
     *
     * This error occurs when a `NetworkPolicy` is initialized with a composition that lacks
     * one or more of the essential species needed to construct the network.
     */
    class MissingSeedSpeciesError final : public PolicyError {
        using PolicyError::PolicyError;
    };

    /**
     * @class MissingKeyReactionError
     * @brief Exception thrown when a constructed network fails to include a key reaction required by the policy.
     *
     * This error is typically thrown from the `construct()` method of a `NetworkPolicy` after
     * the network has been built but fails the final verification step.
     */
    class MissingKeyReactionError final : public PolicyError {
        using PolicyError::PolicyError;
    };
}