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
        ReactionError(const std::string& msg, const std::string& reactionId): m_message(msg), m_reactionID(reactionId) {}
        const char* what() const noexcept override {
            return std::format("Reaction {}: {}", m_reactionID, m_message).c_str();
        }
    };

    class ReactionParsingError final : public ReactionError {
        using ReactionError::ReactionError;
    };
}