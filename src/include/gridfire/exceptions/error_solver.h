#pragma once

#include <exception>
#include <string>

namespace gridfire::exceptions {
    class SolverError : public std::exception {
    public:
        SolverError(std::string msg) : m_msg(std::move(msg)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_msg.c_str();
        }
    private:
        std::string m_msg;
    };

    class CVODESolverFailureError final : public SolverError {
        using SolverError::SolverError;
    };

    class SingularJacobianError final : public SolverError {
        using SolverError::SolverError;
    };
}