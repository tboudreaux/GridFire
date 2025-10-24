#pragma once

#include <exception>
#include <string>
#include <utility>

namespace gridfire::exceptions {
    class UtilityError : public std::exception {};

    class HashingError final : public UtilityError {
    public:
        explicit HashingError() = default;

        explicit HashingError(std::string message)
            : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };
}
