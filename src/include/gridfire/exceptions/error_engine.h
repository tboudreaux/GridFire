#pragma once

#include <exception>
#include <string>
#include <iostream>

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
        explicit StaleEngineTrigger(const state &s)
            :  m_state(s) {}

        const char* what() const noexcept override{
            return "Engine reports stale state. This means that the caller should trigger a update of the engine state before continuing with the integration. If you as an end user are seeing this error, it is likely a bug in the code that should be reported. Please provide the input parameters and the context in which this error occurred. Thank you for your help!";
        }

        state getState() const {
            return m_state;
        }

        size_t numSpecies() const {
            return m_state.m_Y.size();
        }

        size_t totalSteps() const {
            return m_state.m_total_steps;
        }

        double energy() const {
            return m_state.m_eps_nuc;
        }

        double getMolarAbundance(const size_t index) const {
            if (index > m_state.m_Y.size() - 1) {
                throw std::out_of_range("Index out of bounds for molar abundance vector.");
            }
            return m_state.m_Y[index];
        }

        double temperature() const {
            return m_state.m_T9 * 1e9; // Convert T9 back to Kelvin
        }

        double density() const {
            return m_state.m_rho;
        }
    private:
        state m_state;

    };

    class StaleEngineError final : public EngineError {
    public:
        explicit StaleEngineError(const std::string& message)
            : m_message(message) {}

        const char* what() const noexcept override {
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

    class NetworkResizedError final : public EngineError {
    public:
        explicit NetworkResizedError(const std::string& message)
            : m_message(message) {}

        const char* what() const noexcept override {
            return m_message.c_str();
        }
    private:
        std::string m_message;
    };

    class UnableToSetNetworkReactionsError final : public EngineError {
    public:
        explicit UnableToSetNetworkReactionsError(const std::string& message)
            : m_message(message) {}

        const char* what() const noexcept override {
            return m_message.c_str();
        }

    private:
        std::string m_message;
    };

}