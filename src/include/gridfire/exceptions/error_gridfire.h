#pragma once
#include <exception>

#include <string>

namespace gridfire::exceptions {
    class GridFireError : public std::exception {
    public:
        explicit GridFireError(std::string msg) : m_msg(std::move(msg)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_msg.c_str();
        }
    private:
        std::string m_msg;
    };
}
