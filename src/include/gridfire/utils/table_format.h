#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <memory>


namespace gridfire::utils {
    class ColumnBase {
    public:
        virtual ~ColumnBase() = default;
        // Gets the string representation of the data at a given row
        virtual std::string getCellData(size_t rowIndex) const = 0;
        // Gets the header text for the column
        virtual std::string getHeader() const = 0;
        // Gets the number of data rows in the column
        virtual size_t getRowCount() const = 0;
    };

    template<typename T>
    class Column final : public ColumnBase {
    public:
        Column(const std::string& header, const std::vector<T>& data)
            : m_header(header), m_data(data) {}

        std::string getCellData(size_t rowIndex) const override {
            std::stringstream ss;
            if (rowIndex < m_data.size()) {
                ss << m_data[rowIndex];
            }
            return ss.str();
        }

        std::string getHeader() const override {
            return m_header;
        }

        size_t getRowCount() const override {
            return m_data.size();
        }
    private:
        std::string m_header;
        std::vector<T> m_data;

    };

    inline std::string format_table(const std::string& tableName, const std::vector<std::unique_ptr<ColumnBase>>& columns) {
        // --- 1. Handle Empty Table ---
        if (columns.empty()) {
            return tableName + "\n(Table has no columns)\n";
        }

        // --- 2. Determine dimensions and calculate column widths ---
        size_t num_cols = columns.size();
        size_t num_rows = 0;
        for(const auto& col : columns) {
            num_rows = std::max(num_rows, col->getRowCount());
        }

        std::vector<size_t> col_widths(num_cols);
        for (size_t j = 0; j < num_cols; ++j) {
            col_widths[j] = columns[j]->getHeader().length();
            for (size_t i = 0; i < num_rows; ++i) {
                col_widths[j] = std::max(col_widths[j], columns[j]->getCellData(i).length());
            }
        }

        // --- 3. Build the table string using stringstream ---
        std::stringstream table_ss;

        // --- Table Title ---
        size_t total_width = std::accumulate(col_widths.begin(), col_widths.end(), 0) + (num_cols * 3) + 1;
        size_t title_padding = (total_width > tableName.length()) ? (total_width - tableName.length()) / 2 : 0;
        table_ss << std::string(title_padding, ' ') << tableName << "\n";

        // --- Helper to draw horizontal border ---
        auto draw_border = [&]() {
            table_ss << "+";
            for (size_t width : col_widths) {
                table_ss << std::string(width + 2, '-'); // +2 for padding
                table_ss << "+";
            }
            table_ss << "\n";
        };

        // --- Draw Top Border ---
        draw_border();

        // --- Draw Header Row ---
        table_ss << "|";
        for (size_t j = 0; j < num_cols; ++j) {
            table_ss << " " << std::left << std::setw(col_widths[j]) << columns[j]->getHeader() << " |";
        }
        table_ss << "\n";

        // --- Draw Separator ---
        draw_border();

        // --- Draw Data Rows ---
        for (size_t i = 0; i < num_rows; ++i) {
            table_ss << "|";
            for (size_t j = 0; j < num_cols; ++j) {
                table_ss << " " << std::left << std::setw(col_widths[j]) << columns[j]->getCellData(i) << " |";
            }
            table_ss << "\n";
        }

        // --- Draw Bottom Border ---
        draw_border();

        return table_ss.str();
    }



}