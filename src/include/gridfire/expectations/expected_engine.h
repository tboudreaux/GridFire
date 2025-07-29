#pragma once

#include <string>
#include <iostream>

namespace gridfire::expectations {
    enum class EngineErrorTypes {
        FAILURE,
        INDEX,
        STALE
    };

    enum class StaleEngineErrorTypes {
        SYSTEM_RESIZED
    };

    // TODO: rename this to EngineExpectation or something similar
    struct EngineError {
        std::string m_message;
        const EngineErrorTypes type = EngineErrorTypes::FAILURE;

        explicit EngineError(const std::string &message, const EngineErrorTypes type)
            : m_message(message), type(type) {}

        virtual ~EngineError() = default;

        friend std::ostream& operator<<(std::ostream& os, const EngineError& e) {
            os << "EngineError: " << e.m_message;
            return os;
        }
    };

    struct EngineIndexError : EngineError {
        int m_index;

       explicit EngineIndexError(const int index)
            : EngineError("Index error occurred", EngineErrorTypes::INDEX), m_index(index) {}
        friend std::ostream& operator<<(std::ostream& os, const EngineIndexError& e) {
            os << "EngineIndexError: " << e.m_message << " at index " << e.m_index;
            return os;
        }
    };

    struct StaleEngineError : EngineError {
        StaleEngineErrorTypes staleType;

        explicit StaleEngineError(const StaleEngineErrorTypes sType)
            : EngineError("Stale engine error occurred", EngineErrorTypes::STALE), staleType(sType) {}

        explicit operator std::string() const {
            switch (staleType) {
                case (StaleEngineErrorTypes::SYSTEM_RESIZED):
                    return "StaleEngineError: System resized, please update the engine.";
                default:
                    return "StaleEngineError: Unknown stale error type.";
            }
        }
    };
}