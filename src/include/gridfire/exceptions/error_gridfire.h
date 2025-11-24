#pragma once
#include <exception>

#include <string>

namespace gridfire::exceptions {
    /**
     * @brief Base class for GridFire exceptions.
     *
     * This class serves as the base for all exceptions specific to the GridFire library.
     * It extends the standard std::exception class and allows for custom error messages.
     */
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
