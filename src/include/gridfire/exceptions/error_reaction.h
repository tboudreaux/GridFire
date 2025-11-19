#pragma once

#include <exception>
#include <string>
#include <format>

namespace gridfire::exceptions {
    class ReactionError : public std::exception {
    private:
        std::string m_message;
        std::string m_reactionID;
    public:
        ReactionError(const std::string& msg, const std::string& reactionId) {
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