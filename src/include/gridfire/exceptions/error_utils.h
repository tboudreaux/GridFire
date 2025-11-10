#pragma once

#include <exception>
#include <string>
#include <utility>

namespace gridfire::exceptions {
    class UtilityError : public std::exception {
    public:
        explicit UtilityError(std::string  message) : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

    class HashingError final : public UtilityError {
    public:

        explicit HashingError(const std::string &message) : UtilityError(message) {}

    };
}
