#pragma once

#include <exception>
#include <string>

namespace gridfire::exceptions {
    class PolicyError : std::exception {
    public:
        explicit PolicyError(const std::string& msg) : m_message(msg) {};

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

    class MissingBaseReactionError final : public PolicyError {
    public:
        explicit MissingBaseReactionError(const std::string& msg) : PolicyError(msg) {};
    };

    class MissingSeedSpeciesError final : public PolicyError {
    public:
        explicit MissingSeedSpeciesError(const std::string& msg) : PolicyError(msg) {};
    };

    class MissingKeyReactionError final : public PolicyError {
    public:
        explicit MissingKeyReactionError(const std::string& msg) : PolicyError(msg) {}
    };
}