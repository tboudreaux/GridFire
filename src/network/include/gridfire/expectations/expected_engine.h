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

    struct EngineError {
        std::string m_message;
        EngineErrorTypes type = EngineErrorTypes::FAILURE;
        friend std::ostream& operator<<(std::ostream& os, const EngineError& e) {
            os << "EngineError: " << e.m_message;
            return os;
        }
    };

    struct EngineIndexError : EngineError {
        int m_index;
        EngineErrorTypes type = EngineErrorTypes::INDEX;
        friend std::ostream& operator<<(std::ostream& os, const EngineIndexError& e) {
            os << "EngineIndexError: " << e.m_message << " at index " << e.m_index;
            return os;
        }
    };

    struct StaleEngineError : EngineError {
        EngineErrorTypes type = EngineErrorTypes::STALE;
        StaleEngineErrorTypes staleType;

        explicit StaleEngineError(StaleEngineErrorTypes staleType) : staleType(staleType) {}

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