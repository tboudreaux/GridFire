#pragma once

#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace gridfire::exceptions {
    class EngineError : public std::exception {};

    class StaleEngineTrigger final : public EngineError {
    public:
        struct state {
            double m_T9;
            double m_rho;
            std::vector<double> m_Y;
            double m_t;
            int m_total_steps;
            double m_eps_nuc;
        };
        explicit StaleEngineTrigger(state s)
            :  m_state(std::move(s)) {}

        [[nodiscard]] const char* what() const noexcept override{
            return "Engine reports stale state. This means that the caller should trigger a update of the engine state before continuing with the integration. If you as an end user are seeing this error, it is likely a bug in the code that should be reported. Please provide the input parameters and the context in which this error occurred. Thank you for your help!";
        }

        [[nodiscard]] state getState() const {
            return m_state;
        }

        [[nodiscard]] size_t numSpecies() const {
            return m_state.m_Y.size();
        }

        [[nodiscard]] size_t totalSteps() const {
            return m_state.m_total_steps;
        }

        [[nodiscard]] double energy() const {
            return m_state.m_eps_nuc;
        }

        [[nodiscard]] double getMolarAbundance(const size_t index) const {
            if (index > m_state.m_Y.size() - 1) {
                throw std::out_of_range("Index out of bounds for molar abundance vector.");
            }
            return m_state.m_Y[index];
        }

        [[nodiscard]] double temperature() const {
            return m_state.m_T9 * 1e9; // Convert T9 back to Kelvin
        }

        [[nodiscard]] double density() const {
            return m_state.m_rho;
        }
    private:
        state m_state;

    };

    class StaleEngineError final : public EngineError {
    public:
        explicit StaleEngineError(std::string  message)
            : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }

    private:
        std::string m_message;
    };

    class FailedToPartitionEngineError final : public EngineError {
    public:
        explicit FailedToPartitionEngineError(std::string  message)
            : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

    class NetworkResizedError final : public EngineError {
    public:
        explicit NetworkResizedError(std::string  message)
            : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

    class UnableToSetNetworkReactionsError final : public EngineError {
    public:
        explicit UnableToSetNetworkReactionsError(std::string  message)
            : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }

    private:
        std::string m_message;
    };

}