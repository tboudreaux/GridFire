#pragma once

#include <exception>
#include <string>
#include <iostream>

namespace gridfire::exceptions {
    class EngineError : public std::exception {};

    class StaleEngineError final : public EngineError {
    public:
        explicit StaleEngineError(const std::string& message)
            : m_message(message) {}

        const char* what() const noexcept override{
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

    class FailedToPartitionEngineError final : public EngineError {
    public:
        explicit FailedToPartitionEngineError(const std::string& message)
            : m_message(message) {}

        const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

}