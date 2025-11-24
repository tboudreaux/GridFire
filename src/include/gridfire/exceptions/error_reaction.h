#pragma once

#include <string>
#include <format>

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    /**
     * @class ReactionError
     * @brief Base class for all exceptions related to reaction operations.
     *
     * This exception is the parent for more specific reaction-related errors. Catching this
     * type will catch any exception originating from reaction handling.
     */
    class ReactionError : public GridFireError {
    private:
        std::string m_message;
        std::string m_reactionID;
    public:
        ReactionError(const std::string& msg, const std::string& reactionId) : GridFireError(msg) {
            m_reactionID = reactionId;
            m_message = std::format("Reaction {}: {}", reactionId, msg);
        }

        const char* what() const noexcept override {
            return m_message.c_str();
        }
    };

    class ReactionParsingError final : public ReactionError {
        using ReactionError::ReactionError;
    };
}